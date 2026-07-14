"""Oxide database built for phase systems.

Reads blade_generated_data.xlsx (from 01_compositions.py) and Materials Project
reference data to build energy tables, per-parent Excel sheets, and system energy summaries.

Run after:
    02_poscars.py   (download MP reference structures)
    03_energy.py    (relax MP structures with MLIP)

Output files (written to BLADE/Files/):
    all_energy_data.xlsx
    blade_oxide_data.xlsx
    system_energy_dH/<system>.xlsx
    parent_formula_excel_files/<source>/<formula>.xlsx
    pure_element_energies.xlsx
"""

from pathlib import Path
from blade.oxidation.database import OxideDatabase

path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
files_dir = path1 / "Files"

mlip_sources = [
    ("MP+ORB", "MaterialsProject_Comps_ORB"),
    # ("MP+GRACE", "MaterialsProject_Comps_GRACE"),
]

fallback_refs = {
    "C":  -9.200, "Cr": -9.632, "Hf": -9.956, "Mo": -10.850,
    "Nb": -10.094, "O": -4.949, "Ta": -11.853, "Ti":  -7.897,
    "V":  -9.080,  "W": -12.960, "Zr":  -8.547, "B": -6.68,
}

db = OxideDatabase(
    files_dir=files_dir,
    mlip_sources=mlip_sources,
    fallback_refs=fallback_refs,
    poscar_folder="MaterialsProject_Comps_POSCARs",
    fixed_elements=frozenset({"B"}),
)

# Scan BLADE CONTCARs → blade_generated_data.xlsx (run before db.run())
db.scan_blade_data(
    blade_comp_dir=files_dir / "Comps",
    mlip_ref_label="MP+ORB",
    mlip_ref_folder=files_dir / "MaterialsProject_Comps_ORB",
    oxygen_element="O",
)

db.run()

# Collect CONTCAR + energy + TDB files per system into Files/system_structures/
db.collect_structures(
    blade_comp_dirs=[
        files_dir / "Comps",          # add additional Comps folders here
        # files_dir / "Comps_scraps",
        # files_dir / "Comps_MAX",
    ],
    mlip_ref_folder=files_dir / "MaterialsProject_Comps_ORB",
    oxygen_element="O",
)
