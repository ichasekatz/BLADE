"""Step 1: fetch the Zr-B, Nb-B, Ta-B, Zr-Nb-B and Zr-Ta-B phases that the BLADE competitor sets lack.

The step-0 audit found that BLADE's competitor sets for NbZrB / TaZrB contain only the elements, the MB2 SQS
end members and the ordered ZrMB4 (mp-1215209 / mp-1215211). They have no metal-rich borides (MB, M3B4, M5B6,
M3B2, M2B) and no B-rich ones (ZrB12). These phases bound the boron chemical potential, so vacancy formation
energies and B-deficient decomposition cannot be evaluated without them.

Downloads every MP entry in those chemical systems with energy_above_hull <= EHULL_MAX (MP's own hull),
writes POSCARs to inputs/competitors/<mp-id>_<formula>/POSCAR and a manifest. Nothing in BLADE's Files/ is
touched. The MP key is read like examples/full_framework.py does (env MP_API_KEY, else the toml) and never
printed.
"""

from __future__ import annotations

import os
import tomllib
from pathlib import Path

import pandas as pd

BLADE = Path(__file__).resolve().parents[3]
WF = Path(__file__).resolve().parents[1]
CHEMSYS = ["Zr-B", "Nb-B", "Ta-B", "Zr-Nb-B", "Zr-Ta-B", "B", "Zr", "Nb", "Ta"]
EHULL_MAX = 0.10  # eV/atom; metastable phases kept for screening, flagged by their MP e_hull


def api_key() -> str:
    if os.environ.get("MP_API_KEY"):
        return os.environ["MP_API_KEY"]
    db = tomllib.loads((BLADE / "examples" / "full_framework.toml").read_text()).get("database", {})
    if str(db.get("api_key", "")).strip():
        return str(db["api_key"]).strip()
    name = str(db.get("api_key_env", "MP_API_KEY")).strip()
    return os.environ.get(name, name)  # older inputs store the literal key in api_key_env


def main() -> None:
    from mp_api.client import MPRester

    out = WF / "inputs" / "competitors"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    with MPRester(api_key()) as mpr:
        docs = mpr.materials.summary.search(
            chemsys=CHEMSYS, energy_above_hull=(0, EHULL_MAX),
            fields=["material_id", "formula_pretty", "structure", "energy_above_hull", "formation_energy_per_atom",
                    "symmetry", "nsites", "theoretical"])
    for d in docs:
        mid, formula = str(d.material_id), d.formula_pretty
        folder = out / f"{mid}_{formula}"
        folder.mkdir(exist_ok=True)
        d.structure.to(filename=str(folder / "POSCAR"), fmt="poscar")
        rows.append({"material_id": mid, "formula": formula, "spacegroup": d.symmetry.symbol, "nsites": d.nsites,
                     "mp_e_above_hull_eV_atom": d.energy_above_hull, "mp_formation_energy_eV_atom": d.formation_energy_per_atom,
                     "theoretical": d.theoretical, "poscar": str(folder.relative_to(WF) / "POSCAR")})
    df = pd.DataFrame(rows).sort_values(["formula", "mp_e_above_hull_eV_atom"])
    df.to_csv(out / "manifest.csv", index=False)
    print(df.drop(columns="poscar").to_string(index=False))


if __name__ == "__main__":
    main()
