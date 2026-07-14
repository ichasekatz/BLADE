"""Example: extracting volume, lattice parameters, and energy values with BLADEData.

This script shows how to scan a composition directory for POSCAR and energy files and
inspect the resulting DataFrame. The DataFrame is suitable for exporting to
CSV, Excel, or feeding into downstream plotting and filtering workflows.

Prerequisites:
    - numpy, pandas must be installed.
    - POSCAR and energy files must exist at the paths used below.

    python examples/structures/05_data_analysis.py
"""

from pathlib import Path

import pandas as pd

from blade.analysis.blade_data import BLADEData

# ------------------------------------------------------------------
# Paths — adjust to your environment
# ------------------------------------------------------------------
path1 = Path("/Users/chasekatz/Desktop/School/Research/BLADE")

# ------------------------------------------------------------------
# Example 1: scan a single composition directory
# ------------------------------------------------------------------
print("=" * 60)
print("Example 1: Single composition (CrHf)")
print("=" * 60)

data = BLADEData()
comp_dir = path1 / "Files" / "Comps" / "CrHf"

if comp_dir.exists():
    df = data.scan_poscars(comp_dir)
    print(f"\nFound {len(df)} POSCAR(s)")
    if not df.empty:
        cols = ["phase_folder", "sqs_level", "natoms", "volume_per_atom_A3", "energy_eV", "energy_per_atom_eV"]
        print(df[[c for c in cols if c in df.columns]].to_string())
else:
    print(f"Directory not found: {comp_dir} — skipping.")

# ------------------------------------------------------------------
# Example 2: scan all composition directories and combine results
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 2: All compositions in Comps/")
print("=" * 60)

comps_dir = path1 / "Files" / "Comps"
all_frames = []

if comps_dir.exists():
    for comp_dir in sorted(p for p in comps_dir.iterdir() if p.is_dir()):
        df_comp = data.scan_poscars(comp_dir)
        if not df_comp.empty:
            all_frames.append(df_comp)

if all_frames:
    combined = pd.concat(all_frames, ignore_index=True)
    print(f"\nTotal POSCARs found: {len(combined)}")

    print("\nVolume per atom (Å³) by phase:")
    print(combined.groupby("phase_folder")["volume_per_atom_A3"].describe().round(3).to_string())

    if "energy_per_atom_eV" in combined.columns:
        has_energy = combined["energy_per_atom_eV"].notna()
        print(f"\nPOSCARs with energy files: {has_energy.sum()} / {len(combined)}")
        if has_energy.any():
            print("\nEnergy per atom (eV) by phase:")
            print(combined[has_energy].groupby("phase_folder")["energy_per_atom_eV"].describe().round(4).to_string())

    # Export to CSV
    out_csv = path1 / "data_summary.csv"
    combined.to_csv(out_csv, index=False)
    print(f"\nSaved summary -> {out_csv}")
else:
    print("No POSCARs found in any composition directory.")
