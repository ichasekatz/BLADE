"""Example: running the SQS-to-TDB fitting pipeline with BladeTDBGen.

This script shows how to configure and run the CALPHAD database fitting
step. BladeTDBGen wraps MaterialsFramework's Sqs2tdb workflow across many
chemical compositions, running one fit per composition in sequence.

Key difference from the original BLADE: construction is side-effect-free.
The fitting loop only executes when gen.fit() is called explicitly.

Prerequisites:
    - MaterialsFramework must be installed with the GRACE calculator.
    - SQS structures must already exist (run 02_sqs_generation.py first).
    - A CUDA-capable GPU is strongly recommended.

    python examples/03_tdb_generation.py
"""

from pathlib import Path

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_tdb_gen import BladeTDBGen

# ------------------------------------------------------------------
# Paths — adjust to your environment
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# ------------------------------------------------------------------
# Phase list (must match what was used for SQS generation)
# ------------------------------------------------------------------
phase_list = [
    {"generator_name": "PHASE", "lattice": "PHASE1", "supercell_size": (4, 3, 2)},
]

# Fixed-sublattice species (B sublattice in HEDB1)
sqsgen_levels2 = [
    {"element": "B", "compositions": "1.0", "letter": "b", "count": "2"},
]

# ------------------------------------------------------------------
# Compositions
# ------------------------------------------------------------------
composer = BladeCompositions(
    primary_elements=["Cr", "Hf"],
    secondary_elements=[],
    system_size=2,
    primary_min=2,
    primary_max=2,
    secondary_min=0,
    secondary_max=0,
    allow_lower_order=False,
)
composition_list = composer.generate_compositions()
print(f"Fitting TDBs for {len(composition_list)} composition(s): {composition_list}")

# ------------------------------------------------------------------
# Configure and run the fitting pipeline
# ------------------------------------------------------------------
# Step 1: construct the object — no computation happens here
gen = BladeTDBGen(
    phases=phase_list,
    liquid=False,
    paths=paths,
    composition_list=composition_list,
    level=5,
    sqsgen_levels2=sqsgen_levels2,
    skip_existing=True,   # skip compositions that already have a .tdb file
)

# Step 2: run the fitting loop explicitly
gen.fit()

print("\nTDB generation complete.")
print(f"Output TDB files are located in: {path1}")
