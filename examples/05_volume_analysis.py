"""Example: extracting volume and lattice parameters with BLADEVolume.

This script shows how to scan a composition directory for POSCAR files and
inspect the resulting DataFrame. The DataFrame is suitable for exporting to
CSV, Excel, or feeding into downstream plotting and filtering workflows.

Prerequisites:
    - numpy, pandas must be installed.
    - POSCAR files must exist at the paths used below.

    python examples/05_volume_analysis.py
"""

from pathlib import Path

import pandas as pd

from blade.analysis.blade_volume import BLADEVolume

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

vol = BLADEVolume()
comp_dir = path1 / "Files" / "Comps" / "CrHf"

if comp_dir.exists():
    df = vol.scan_poscars(comp_dir)
    print(f"\nFound {len(df)} POSCAR(s)")
    if not df.empty:
        print(df[["phase_folder", "sqs_level", "natoms", "volume_per_atom_A3"]].to_string())
else:
    print(f"Directory not found: {comp_dir} — skipping.")

# ------------------------------------------------------------------
# Example 2: scan all composition directories and combine results
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 2: All compositions in path1")
print("=" * 60)

all_frames = []
for comp_dir in sorted(p for p in path1.iterdir() if p.is_dir()):
    df_comp = vol.scan_poscars(comp_dir)
    if not df_comp.empty:
        all_frames.append(df_comp)

if all_frames:
    combined = pd.concat(all_frames, ignore_index=True)
    print(f"\nTotal POSCARs found: {len(combined)}")
    print("\nVolume per atom (Å³) by phase:")
    print(
        combined.groupby("phase_folder")["volume_per_atom_A3"]
        .describe()
        .round(3)
        .to_string()
    )

    # Export to CSV
    out_csv = path1 / "volume_summary.csv"
    combined.to_csv(out_csv, index=False)
    print(f"\nSaved summary -> {out_csv}")
else:
    print("No POSCARs found in any composition directory.")
