"""Collect CONTCAR + energy files for each formula in system_energy_dH Excel files.

For each system (e.g. CrHfB), gathers:
  - GRACE-relaxed MP structures from MaterialsProject_Comps_GRACE/
  - BLADE boride structures from Files/Comps/

For GRACE: keeps lowest-energy CONTCAR per formula.
For BLADE: keeps all SQS levels (different compositions within the system).

Output layout:
    Files/system_structures/<system>/grace/<formula>_<entry_id>/{CONTCAR,energy}
    Files/system_structures/<system>/blade/<sqs_folder>/{CONTCAR,energy}
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition, Structure

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
files_dir = path1 / "Files"

system_dH_dir    = files_dir / "system_energy_dH"
parent_excel_dir = files_dir / "parent_formula_excel_files"
grace_base       = files_dir / "MaterialsProject_Comps_GRACE"
blade_comp_dir   = files_dir / "Comps"
output_dir       = files_dir / "system_structures"

output_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Build stable material_id set from MP summary
# ------------------------------------------------------------------
_mp_summary = pd.read_excel(grace_base / "materials_project_summary.xlsx")
_mp_summary["energy_above_hull"] = pd.to_numeric(_mp_summary["energy_above_hull"], errors="coerce")
_mp_summary["is_stable"] = _mp_summary["is_stable"].apply(
    lambda v: str(v).strip().lower() in {"true", "1", "yes"}
)
stable_ids: set[str] = set(
    _mp_summary[_mp_summary["is_stable"] | (_mp_summary["energy_above_hull"] == 0)]["material_id"]
    .astype(str)
    .str.strip()
)
print(f"Stable MP entries: {len(stable_ids)}")


def normalize(formula: str) -> str:
    try:
        return Composition(str(formula)).reduced_formula
    except Exception:
        return str(formula).strip()


def formula_elements(formula: str) -> frozenset[str]:
    try:
        return frozenset(Composition(str(formula)).get_el_amt_dict().keys())
    except Exception:
        return frozenset()


def read_energy(path: Path) -> float | None:
    try:
        return float(path.read_text().strip().split()[0])
    except Exception:
        return None


# ------------------------------------------------------------------
# Build GRACE index: reduced_formula → lowest-energy entry
# ------------------------------------------------------------------
print("Indexing GRACE CONTCARs...")
# GRACE structure: energy file contains total energy; divide by n_atoms for epa
grace_index: dict[str, dict] = {}   # formula → {entry_id, contcar, energy_path, epa}

for contcar in sorted(grace_base.rglob("CONTCAR")):
    energy_path = contcar.parent / "energy"
    entry_id = contcar.parent.name          # e.g. mp-12345_BHf
    mp_id    = entry_id.split("_")[0]       # e.g. mp-12345
    if mp_id not in stable_ids:
        continue
    try:
        s = Structure.from_file(str(contcar))
        f = s.composition.reduced_formula
        n_atoms = len(s)
    except Exception:
        continue
    total_e = read_energy(energy_path)
    epa = (total_e / n_atoms) if total_e is not None else None

    if f not in grace_index or (epa is not None and (grace_index[f]["epa"] is None or epa < grace_index[f]["epa"])):
        grace_index[f] = {
            "entry_id": entry_id,
            "contcar":  contcar,
            "energy_path": energy_path if energy_path.exists() else None,
            "epa": epa,
        }

print(f"  {len(grace_index)} unique formulas in GRACE")

# ------------------------------------------------------------------
# Build POSCAR fallback index for formulas missing from GRACE
# ------------------------------------------------------------------
poscar_base = files_dir / "MaterialsProject_Comps"
poscar_fallback: dict[str, dict] = {}
if poscar_base.exists():
    for poscar in sorted(poscar_base.rglob("POSCAR")):
        entry_id = poscar.parent.name
        mp_id    = entry_id.split("_")[0]
        if mp_id not in stable_ids:
            continue
        try:
            s = Structure.from_file(str(poscar))
            f = s.composition.reduced_formula
        except Exception:
            continue
        if f not in poscar_fallback:
            poscar_fallback[f] = {"entry_id": entry_id, "poscar": poscar, "energy_path": None, "epa": None}
    print(f"  {len(poscar_fallback)} unique formulas in POSCAR fallback")

# ------------------------------------------------------------------
# Build BLADE index: reduced_formula → list of SQS entries
# ------------------------------------------------------------------
print("Indexing BLADE CONTCARs...")
blade_index: dict[str, list[dict]] = {}   # formula → [{entry_id, contcar, energy_path}]

for contcar in sorted(blade_comp_dir.rglob("CONTCAR")):
    energy_path = contcar.parent / "energy"
    try:
        s = Structure.from_file(str(contcar))
        f = s.composition.reduced_formula
    except Exception:
        continue
    entry_id = contcar.parent.name   # sqs_lev=... folder name
    blade_index.setdefault(f, []).append({
        "entry_id": entry_id,
        "contcar":  contcar,
        "energy_path": energy_path if energy_path.exists() else None,
        "comp_dir": contcar.parent.parent.parent.name,  # e.g. CrHf
    })

print(f"  {len(blade_index)} unique formulas in BLADE")


# ------------------------------------------------------------------
# Helper: copy CONTCAR + energy into dest dir
# ------------------------------------------------------------------
def copy_files(contcar: Path, energy_path: Path | None, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(contcar, dest / "CONTCAR")
    if energy_path and energy_path.exists():
        shutil.copy2(energy_path, dest / "energy")


# ------------------------------------------------------------------
# Process each system
# ------------------------------------------------------------------
for xlsx_path in sorted(p for p in system_dH_dir.glob("*.xlsx") if not p.name.startswith("~$")):
    system_name = xlsx_path.stem
    print(f"\nSystem: {system_name}")

    df = pd.read_excel(xlsx_path)
    if "formula" not in df.columns:
        print(f"  SKIP: no formula column")
        continue

    try:
        sys_els = frozenset(Composition(system_name).get_el_amt_dict().keys())
    except Exception:
        sys_els = frozenset()

    sys_out   = output_dir / system_name
    grace_out = sys_out / "grace"
    blade_out = sys_out / "blade"

    # GRACE structures: formulas listed in the system Excel file
    n_grace = n_grace_missing = 0
    for raw in df["formula"].dropna().unique():
        f = normalize(str(raw))
        hit = grace_index.get(f)
        if hit:
            dest = grace_out / f"{f}_{hit['entry_id']}"
            copy_files(hit["contcar"], hit["energy_path"], dest)
            n_grace += 1
        else:
            fb = poscar_fallback.get(f)
            if fb:
                dest = grace_out / f"{f}_{fb['entry_id']}"
                copy_files(fb["poscar"], None, dest)
                n_grace += 1
            else:
                n_grace_missing += 1

    # BLADE borides: all SQS entries whose (elements - O) ⊆ system elements
    n_blade = 0
    blade_seen: set[Path] = set()
    for f, entries in blade_index.items():
        f_els = formula_elements(f)
        if not f_els:
            continue
        if not (f_els - {"O"}).issubset(sys_els):
            continue
        for hit in entries:
            if hit["contcar"] in blade_seen:
                continue
            # Only use structures from a comp_dir whose elements ⊆ sys_els
            # (prevents CrV-system endmembers like VB2 appearing in MoVB)
            comp_els = formula_elements(hit["comp_dir"])
            if comp_els and not comp_els.issubset(sys_els):
                continue
            blade_seen.add(hit["contcar"])
            dest = blade_out / f"{f}_{hit['entry_id']}"
            copy_files(hit["contcar"], hit["energy_path"], dest)
            n_blade += 1

    # TDB files: comp dirs whose metal elements ⊆ sys_els
    tdb_out = sys_out / "tdb"
    n_tdb = 0
    for comp_dir in blade_comp_dir.iterdir():
        if not comp_dir.is_dir():
            continue
        comp_els = formula_elements(comp_dir.name)
        if not comp_els or not comp_els.issubset(sys_els):
            continue
        for tdb in comp_dir.glob("*.tdb"):
            tdb_out.mkdir(parents=True, exist_ok=True)
            shutil.copy2(tdb, tdb_out / tdb.name)
            n_tdb += 1

    # system_energy_dH Excel for this system
    _sys_xlsx = system_dH_dir / f"{system_name}.xlsx"
    if _sys_xlsx.exists():
        shutil.copy2(_sys_xlsx, sys_out / f"{system_name}_system_energy.xlsx")

    # parent_formula_excel_files: all source subfolders, files whose formula elements ⊆ sys_els
    n_parent = 0
    if parent_excel_dir.exists():
        for src_subdir in parent_excel_dir.iterdir():
            if not src_subdir.is_dir():
                continue
            dest_sub = sys_out / "parent_excel" / src_subdir.name
            for pxlsx in src_subdir.glob("*.xlsx"):
                if pxlsx.name.startswith("~$"):
                    continue
                try:
                    f_els = formula_elements(pxlsx.stem)
                except Exception:
                    f_els = frozenset()
                if not f_els or not (f_els - {"O"}).issubset(sys_els):
                    continue
                dest_sub.mkdir(parents=True, exist_ok=True)
                shutil.copy2(pxlsx, dest_sub / pxlsx.name)
                n_parent += 1

    print(f"  GRACE: {n_grace} copied, {n_grace_missing} not found in index")
    print(f"  BLADE: {n_blade} SQS entries copied, {n_tdb} TDB files")
    print(f"  Parent Excel: {n_parent} files across source subfolders")
    print(f"  → {sys_out}")

print(f"\nDone. Output: {output_dir.resolve()}")
