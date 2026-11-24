"""BLADE binary-boride correction (commit 774db35) reproduced at the competitor level, in a separate directory.

The fix makes BLADE's database stage also query single-metal + B chemical systems from Materials Project. Its
effect on the stability stage is that the stable binary borides of every metal become fixed-phase competitors of
every system containing that metal. Instead of regenerating BLADE's whole database (which would rewrite Files/),
this script reproduces exactly that effect:
  1. Materials Project query for chemsys M-B, M in BLADE's 13 primary metals, is_stable=True (BLADE default
     stable_only = True). MB2 entries are skipped (BLADE supplies its own SQS diboride end members, as with
     skip_blade_formulas), as are pure elements (already present).
  2. Relaxation with BLADE's own database relaxation: materialsframework Sqs2tdb._calculate, ORB
     orb-v3-conservative-inf-omat, fmax 1e-3, 1000 steps, cuda -> CONTCAR + energy (total eV) in the same format.
  3. Copies of Files/system_structures, system_structures_4metal, system_structures_5metal under
     bfix_rerun/ with the relevant binary borides added to each system's ORB/ competitor folder.
Nothing under Files/ is modified. The MP key is read as BLADE does and never printed.
Outputs: bfix_rerun/MP_binary_borides/<M-B>/<mpid_formula>/, bfix_rerun/system_structures*/, bfix_rerun/added_competitors.csv
"""

from __future__ import annotations

import os
import shutil
import sys
import tomllib
from pathlib import Path

import pandas as pd

BLADE = Path("/home/ichasekatz/BLADE/BLADE")
OUT = BLADE / "bfix_rerun"
METALS = ["Cr", "Cu", "Fe", "Ge", "Hf", "Mn", "Mo", "Nb", "Ni", "Ta", "Ti", "Y", "Zr"]


def api_key() -> str:
    if os.environ.get("MP_API_KEY"):
        return os.environ["MP_API_KEY"]
    db = tomllib.loads((BLADE / "examples" / "full_framework.toml").read_text()).get("database", {})
    if str(db.get("api_key", "")).strip():
        return str(db["api_key"]).strip()
    name = str(db.get("api_key_env", "MP_API_KEY")).strip()
    return os.environ.get(name, name)


def fetch() -> pd.DataFrame:
    from mp_api.client import MPRester
    from pymatgen.core import Composition
    from pymatgen.io.vasp import Poscar

    rows = []
    with MPRester(api_key()) as mpr:
        for m in METALS:
            chemsys = "-".join(sorted(["B", m]))
            docs = mpr.materials.summary.search(chemsys=chemsys, is_stable=True,
                                                fields=["material_id", "formula_pretty", "structure", "energy_above_hull",
                                                        "formation_energy_per_atom"])
            for d in docs:
                comp = Composition(d.formula_pretty)
                if comp.reduced_composition == Composition(f"{m}B2").reduced_composition:
                    continue  # any MB2 polymorph is skipped, as BLADE's skip_blade_formulas does (formula match)
                folder = OUT / "MP_binary_borides" / chemsys / f"{d.material_id}_{d.formula_pretty}"
                folder.mkdir(parents=True, exist_ok=True)
                if not (folder / "POSCAR").exists():
                    Poscar(d.structure).write_file(folder / "POSCAR")
                rows.append({"metal": m, "chemsys": chemsys, "material_id": str(d.material_id), "formula": d.formula_pretty,
                             "nsites": len(d.structure), "mp_formation_energy_eV_atom": d.formation_energy_per_atom, "dir": str(folder)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "MP_binary_borides" / "manifest.csv", index=False)
    return df


def relax(df: pd.DataFrame) -> None:
    sys.path.insert(0, str(BLADE / "src"))
    from materialsframework.calculators.registry import get_calculator
    from materialsframework.tools.sqs2tdb import Sqs2tdb

    calc = get_calculator("orb", steps=1000, device="cuda")
    solver = Sqs2tdb(fmax=1e-3, verbose=False, track_trajectory=False, calculator=calc)
    for r in df.itertuples():
        d = Path(r.dir)
        if (d / "energy").exists():
            continue
        try:
            solver._calculate(d, relax=True)
            print(f"relaxed {r.formula} ({r.material_id})", flush=True)
        except Exception as exc:  # recorded, not hidden
            print(f"FAILED {r.formula} ({r.material_id}): {exc}", flush=True)


def augment(df: pd.DataFrame) -> None:
    import re

    added = []
    for tree in ("system_structures", "system_structures_4metal", "system_structures_5metal"):
        src, dst = BLADE / "Files" / tree, OUT / tree
        if not dst.exists():
            shutil.copytree(src, dst, symlinks=False)
        for sysdir in sorted(p for p in dst.iterdir() if p.is_dir()):
            metals = set(re.findall(r"[A-Z][a-z]?", sysdir.name[:-1]))
            orb = sysdir / "ORB"
            orb.mkdir(exist_ok=True)
            for r in df[df.metal.isin(metals)].itertuples():
                d = Path(r.dir)
                if not (d / "energy").exists():
                    continue
                target = orb / f"{r.formula}_{r.material_id}_{r.formula}"
                if not target.exists():
                    target.mkdir()
                    for f in ("CONTCAR", "energy"):
                        shutil.copy(d / f, target / f)
                added.append({"tree": tree, "system": sysdir.name, "added": f"{r.formula} ({r.material_id})"})
    pd.DataFrame(added).to_csv(OUT / "added_competitors.csv", index=False)
    print(pd.DataFrame(added).groupby("tree").size())


def main() -> None:
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    man = OUT / "MP_binary_borides" / "manifest.csv"
    df = fetch() if which in ("fetch", "all") or not man.exists() else pd.read_csv(man)
    print(df[["metal", "formula", "material_id", "nsites"]].to_string(index=False))
    if which in ("relax", "all"):
        relax(df)
    if which in ("augment", "all"):
        augment(df)


if __name__ == "__main__":
    main()
