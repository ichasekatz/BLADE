"""Calculate MLIP energies for Materials Project structures (Section 2b).

Relaxes all POSCARs under MaterialsProject_Comps/ using a configurable MLIP
calculator and writes an energy file per structure directory.

Prerequisites:
    materialsframework (ichasekatz fork)

Input:  BLADE/Files/MaterialsProject_Comps/**/**/POSCAR
Output: BLADE/Files/MaterialsProject_Comps/**/**/energy  (one per structure)
"""

from __future__ import annotations

from pathlib import Path

from materialsframework.calculators.registry import get_calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"

files_dir = path1 / "Files"
files_dir.mkdir(parents=True, exist_ok=True)

mp_comps_dir = files_dir / "MaterialsProject_Comps_ORB"


# ==============================================================================
# Calculate MLIP energies for MP structures using GRACECalculator
# ==============================================================================


def count_atoms_from_poscar(poscar_path: Path) -> int:
    lines = poscar_path.read_text().splitlines()

    tokens = lines[5].split()

    if all(tok.isdigit() for tok in tokens):
        counts = [int(tok) for tok in tokens]
    else:
        counts = [int(tok) for tok in lines[6].split()]

    return sum(counts)


def write_energy_per_atom(subdir: Path) -> None:
    energy_path = subdir / "energy"
    poscar_path = subdir / "POSCAR"
    energy_per_atom_path = subdir / "energy_per_atom"

    if not energy_path.exists():
        return

    energy_text = energy_path.read_text().strip()
    if not energy_text:
        return

    energy = float(energy_text.split()[0])
    n_atoms = count_atoms_from_poscar(poscar_path)

    energy_per_atom = energy / n_atoms
    energy_per_atom_path.write_text(f"{energy_per_atom:.16f}\n")


# ------------------------------------------------------------------
# MLIP calculator — change mlip to swap potentials (see MaterialsFramework registry)
# ------------------------------------------------------------------
mlip = "orb"  # e.g. "grace", "mace", "uma", "chgnet", "orb"
mlip_kwargs = {"steps": 1000, "device": "cpu"}  # kwargs forwarded to the calculator constructor

skip_existing_mp_energies = True

calc = get_calculator(mlip, **mlip_kwargs)
track_trajectory = False  # set False to skip relaxation_live.xyz
sqs = Sqs2tdb(fmax=1e-3, verbose=True, track_trajectory=track_trajectory, calculator=calc)

poscar_paths = sorted(mp_comps_dir.rglob("POSCAR"))
print(f"Found {len(poscar_paths)} MP POSCARs — calculating MLIP energies")

for poscar_path in poscar_paths:
    subdir = poscar_path.parent
    if skip_existing_mp_energies and (subdir / "energy").exists():
        write_energy_per_atom(subdir)
        continue
    print(f"  Relaxing: {subdir.relative_to(mp_comps_dir)}")
    try:
        sqs._calculate(subdir, relax=True)
        write_energy_per_atom(subdir)
    except Exception as e:
        print(f"  Failed {subdir.name}: {e}")
