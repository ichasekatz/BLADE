"""Step 3: ORB screening relaxations with exactly the settings BLADE used for its SQS energies.

Settings (from examples/full_framework.toml + materialsframework.tools.sqs2tdb): ORBCalculator
model "orb-v3-conservative-inf-omat", precision "float32-high", FIRE optimiser, full cell relaxation
(FrechetCellFilter), fmax = 0.001 eV/A ([tdb_fit] fmax), steps = 1000 ([mlip_kwargs] steps).

Jobs:
  check     single-point energies of BLADE's relaxed SQS CONTCARs (x = 0, 0.25, 0.5, 0.75, 1) -> confirms the
            installed ORB reproduces BLADE's stored energies, so BLADE energies can be reused as references.
  compete   relax the MP competitor structures (inputs/competitors) -> consistent ORB hull for mu_B limits.
  structs   relax every structure in inputs/structure_manifest.csv.

Resumable: a structure with results/orb/<id>/result.json is skipped. MLIP screening only; no VASP, no phonons.
Usage: python s03_orb_relax.py [check|compete|structs|all] [--limit N]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

BLADE = Path(__file__).resolve().parents[3]
WF = Path(__file__).resolve().parents[1]
RES = WF / "results" / "orb"
SETTINGS = {"model": "orb-v3-conservative-inf-omat", "precision": "float32-high", "device": "cuda",
            "steps": 1000, "fmax": 0.001, "optimizer": "FIRE", "relax_cell": True}


def calculator():
    from materialsframework.calculators.orb import ORBCalculator

    calc = ORBCalculator(model=SETTINGS["model"], precision=SETTINGS["precision"], device=SETTINGS["device"],
                         steps=SETTINGS["steps"], optimizer=SETTINGS["optimizer"], relax_cell=SETTINGS["relax_cell"])
    calc.fmax = SETTINGS["fmax"]
    return calc


def run_one(calc, sid: str, poscar: Path, relax: bool = True) -> dict:
    from ase.io import read, write

    out = RES / sid
    if (out / "result.json").exists():
        return json.loads((out / "result.json").read_text())
    out.mkdir(parents=True, exist_ok=True)
    atoms = read(poscar)
    t0 = time.time()
    if relax:
        res = calc.relax(atoms)
        final = res["final_structure"]
        n_steps = len(getattr(res["trajectory"], "forces", [])) - 1
        from pymatgen.io.ase import AseAtomsAdaptor

        write(out / "CONTCAR", AseAtomsAdaptor.get_atoms(final), format="vasp", direct=True)
        fmax_final = float(np.abs(np.asarray(res["forces"])).max()) if "forces" in res else float("nan")
    else:
        res = calc.calculate(atoms)
        n_steps, fmax_final = 0, float(np.abs(np.asarray(res["forces"])).max())
    energy = float(res["energy"])
    info = {"structure_id": sid, "poscar": str(poscar), "relaxed": relax, "energy_eV": energy, "n_atoms": len(atoms),
            "energy_eV_atom": energy / len(atoms), "n_steps": n_steps, "max_force_component_eV_A": fmax_final,
            "converged": bool(n_steps < SETTINGS["steps"]) if relax else None,
            "wall_s": round(time.time() - t0, 1), "settings": SETTINGS}
    (out / "energy").write_text(f"{energy}\n")
    (out / "result.json").write_text(json.dumps(info, indent=1))
    return info


def jobs_check() -> list[tuple[str, Path, bool]]:
    base = BLADE / "Files" / "system_structures"
    out = []
    for sysdir, m in (("NbZrB", "Nb"), ("TaZrB", "Ta")):
        for d in sorted((base / sysdir / "blade").iterdir()):
            out.append((f"check_{sysdir}_{d.name.split('_sqs')[0]}", d / "CONTCAR", False))
    return out


def jobs_compete() -> list[tuple[str, Path, bool]]:
    man = pd.read_csv(WF / "inputs" / "competitors" / "manifest.csv")
    return [(f"comp_{r.material_id}_{r.formula}", WF / r.poscar, True) for r in man.itertuples()]


def jobs_structs() -> list[tuple[str, Path, bool]]:
    man = pd.read_csv(WF / "inputs" / "structure_manifest.csv")
    order = {"Y": -1, "C": 0, "L": 1, "S": 2, "V": 3}  # configurations first (H2 answers fastest)
    man = man.assign(o=man.family.map(order)).sort_values(["o", "n_atoms", "structure_id"])
    return [(r.structure_id, WF / r.poscar, True) for r in man.itertuples()]


def collect() -> pd.DataFrame:
    rows = [json.loads(p.read_text()) for p in RES.glob("*/result.json")]
    df = pd.DataFrame(rows).drop(columns=["settings"], errors="ignore").sort_values("structure_id")
    df.to_csv(WF / "results" / "s03_orb_energies.csv", index=False)
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("what", nargs="?", default="all", choices=["check", "compete", "structs", "all", "collect"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--reverse", action="store_true", help="process the job list backwards (second GPU worker)")
    args = ap.parse_args()
    if args.what != "collect":
        jobs = {"check": jobs_check, "compete": jobs_compete, "structs": jobs_structs}
        todo = [j for k in (["check", "compete", "structs"] if args.what == "all" else [args.what]) for j in jobs[k]()]
        if args.reverse:
            todo = todo[::-1]
        if args.limit:
            todo = todo[: args.limit]
        calc = calculator()
        for i, (sid, path, relax) in enumerate(todo, 1):
            if (RES / sid / "result.json").exists() or (RES / sid / ".lock").exists():
                continue
            (RES / sid).mkdir(parents=True, exist_ok=True)
            (RES / sid / ".lock").touch()
            try:
                info = run_one(calc, sid, path, relax)
            finally:
                (RES / sid / ".lock").unlink(missing_ok=True)
            print(f"[{i}/{len(todo)}] {sid}: {info['energy_eV_atom']:.6f} eV/atom, steps={info['n_steps']}, "
                  f"conv={info['converged']}, {info['wall_s']} s", flush=True)
    collect()


if __name__ == "__main__":
    main()
