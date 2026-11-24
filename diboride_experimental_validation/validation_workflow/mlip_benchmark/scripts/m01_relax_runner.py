"""Model-agnostic relaxation runner, executed INSIDE each isolated MLIP environment (needs only ase, numpy, pandas).

Protocol identical to the ORB screen (materialsframework BaseCalculator.relax as configured by BLADE):
  ASE FIRE optimiser, FrechetCellFilter (full cell + positions), fmax = 0.001 eV/A, max 1000 steps,
  starting from the same geometry file for every model. Energies are total potential energies (eV) of each model;
  they are only ever compared through formation / reaction energies built within ONE model.

Models (checkpoint names recorded in each result.json):
  mace-mpa    mace_mp(model="medium-mpa-0")      MPtrj + sAlex (PBE/PBE+U), MACE 0.3.x
  mace-omat   mace_mp(model="medium-omat-0")     OMat24 (PBE/PBE+U) - same training corpus family as ORB-v3-omat
  mattersim   MatterSimCalculator("MatterSim-v1.0.0-5M.pth")
  sevennet    SevenNetCalculator("7net-mf-ompa", modal="mpa")

Usage: <env>/bin/python m01_relax_runner.py MODEL [--limit N] [--only role]
Resumable: skips structures with results/<MODEL>/<id>/result.json.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

MB = Path(__file__).resolve().parents[1]
FMAX, STEPS = 0.001, 1000


def make_calc(model: str):
    if model == "orb":  # BLADE .venv; same model + precision as materialsframework ORBCalculator used by BLADE
        from materialsframework.calculators.orb import ORBCalculator

        mf = ORBCalculator(model="orb-v3-conservative-inf-omat", device="cuda", precision="float32-high")
        return mf.calculator, "orb-v3-conservative-inf-omat (float32-high, materialsframework ORBCalculator)"
    if model == "mace-mpa":
        from mace.calculators import mace_mp

        return mace_mp(model="medium-mpa-0", device="cuda", default_dtype="float64"), "mace-torch medium-mpa-0 (float64)"
    if model == "mace-omat":
        from mace.calculators import mace_mp

        return mace_mp(model="medium-omat-0", device="cuda", default_dtype="float64"), "mace-torch medium-omat-0 (float64)"
    if model == "mattersim":
        from mattersim.forcefield import MatterSimCalculator

        return MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device="cuda"), "MatterSim-v1.0.0-5M"
    if model == "sevennet":
        from sevenn.calculator import SevenNetCalculator

        return SevenNetCalculator(model="7net-mf-ompa", modal="mpa", device="cuda"), "SevenNet 7net-mf-ompa (modal mpa)"
    raise ValueError(model)


def relax(atoms, calc) -> tuple[object, int, float]:
    from ase.filters import FrechetCellFilter
    from ase.optimize import FIRE

    atoms.calc = calc
    filt = FrechetCellFilter(atoms)
    opt = FIRE(filt, logfile=None)
    steps = 0
    for _ in opt.irun(fmax=FMAX, steps=STEPS):
        steps += 1
    fmax = float(np.sqrt((filt.get_forces() ** 2).sum(axis=1).max()))
    return atoms, steps - 1, fmax


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="")
    ap.add_argument("--manifest", default="pilot_manifest.csv", help="pilot_manifest.csv (Stage 1) or stage2_manifest.csv")
    args = ap.parse_args()
    from ase.io import read, write

    man = pd.read_csv(MB / args.manifest)
    if args.only:
        man = man[man.role == args.only]
    man = man.sort_values("n_atoms")  # small references first (quick failure detection)
    if args.limit:
        man = man.head(args.limit)
    calc, ckpt = make_calc(args.model)
    out_root = MB / "results" / args.model
    for i, r in enumerate(man.itertuples(), 1):
        out = out_root / r.structure_id
        if (out / "result.json").exists():
            continue
        out.mkdir(parents=True, exist_ok=True)
        atoms0 = read(r.start_geometry, format="vasp")
        atoms = atoms0.copy()
        t0 = time.time()
        try:
            atoms, steps, fmax = relax(atoms, calc)
            e = float(atoms.get_potential_energy())
            status = "ok"
        except Exception as exc:  # record, never hide, a failure
            e, steps, fmax, status = float("nan"), -1, float("nan"), f"error: {exc!r}"[:300]
        disp = atoms.get_positions() - atoms0.get_positions() if status == "ok" else None
        info = {"structure_id": r.structure_id, "model": args.model, "checkpoint": ckpt, "status": status, "energy_eV": e,
                "n_atoms": len(atoms), "energy_eV_atom": e / len(atoms), "n_steps": steps, "fmax_final_eV_A": fmax,
                "converged": bool(status == "ok" and steps < STEPS),
                "volume_A3": float(atoms.get_volume()), "volume_start_A3": float(atoms0.get_volume()),
                "cell_lengths": [float(v) for v in atoms.cell.lengths()], "cell_angles": [float(v) for v in atoms.cell.angles()],
                "max_atomic_displacement_A": float(np.linalg.norm(disp, axis=1).max()) if disp is not None else None,
                "wall_s": round(time.time() - t0, 1), "fmax_target": FMAX, "max_steps": STEPS, "optimizer": "FIRE+FrechetCellFilter"}
        if status == "ok":
            write(out / "CONTCAR", atoms, format="vasp", direct=True)
        (out / "result.json").write_text(json.dumps(info, indent=1))
        print(f"[{i}/{len(man)}] {args.model} {r.structure_id}: {status} E={e:.5f} steps={steps} fmax={fmax:.4f} {info['wall_s']}s", flush=True)


if __name__ == "__main__":
    main()
