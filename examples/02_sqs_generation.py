"""Example: generating SQS inputs and running ATAT mcsqs with BladeSQS.

This script shows how to configure and run the SQS generation step for a
single phase prototype. BladeSQS writes the ATAT input files (rndstr.skel,
sqsgen.in), runs sqs2tdb -mk to create sqsdb_lev=* directories, then
executes corrdump and parallel mcsqs in each directory.

Prerequisites:
    - ATAT must be installed and on the PATH (corrdump, mcsqs, sqs2tdb).
    - The output directory (path1) must exist.

    python examples/02_sqs_generation.py
"""

from pathlib import Path

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_sqsgen import BladeSQS

# ------------------------------------------------------------------
# Paths — adjust to your environment
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# ------------------------------------------------------------------
# Phase prototype: HEDB1 (hexagonal, a-sublattice mixing)
# ------------------------------------------------------------------
phases = {
    "PHASE1": {
        "a": 1,
        "b": 1,
        "c": 1,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.333333 0.666667 0.500000 B\n"
            "0.666667 0.333333 0.500000 B\n"
        ),
    }
}

phase_list = [{"generator_name": "PHASE", "lattice": "PHASE1", "supercell_size": (4, 3, 2)}]

# ------------------------------------------------------------------
# SQS composition levels
# ------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]], "letter": ["a"]},
    {"level": 1, "compositions": [[0.5, 0.5]], "letter": ["a"]},
    {"level": 2, "compositions": [[0.75, 0.25]], "letter": ["a"]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]], "letter": ["a"]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]], "letter": ["a"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]], "letter": ["a"]},
    {"level": 6, "compositions": [[0.75, 0.125, 0.125]], "letter": ["a"]},
]

# ------------------------------------------------------------------
# mcsqs runtime parameters
# ------------------------------------------------------------------
mcsqs_params = {
    "use_time": True,
    "time": 30,           # seconds per sqsdb directory
    "2": 5,               # use the 5th-nearest neighbor as the pair cutoff
    "3": 4,               # 4th-nearest for triplets
    "4": 3,               # 3rd-nearest for quadruplets
    "wr": 20,
    "wn": 0.75,
    "wd": 1,
    "parallel_runs": 10,
}

# ------------------------------------------------------------------
# Generate compositions and run SQS
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
unique_len_comps = composer.get_systems()

print(f"Compositions: {composition_list}")
print(f"System sizes: {unique_len_comps}")

level = 5

for specific_phase in phase_list:
    for len_comp in unique_len_comps:
        sqs_gen = BladeSQS(
            phases_dict=phases[specific_phase["lattice"]],
            sqsgen_levels=sqsgen_levels,
            level=level,
            len_comp=len_comp,
            skip_existing_sqs=True,   # skip if bestcorr.out already exists
        )
        params = mcsqs_params | {"super_cell_size": specific_phase["supercell_size"]}
        sqs_gen.sqs_gen(phase=specific_phase, paths=paths, iter=1_000_000, params=params)

# ------------------------------------------------------------------
# Preview only (no ATAT required) — inspect the generated file text
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Preview: rndstr.skel and sqsgen.in for len_comp=3, level=5")
print("=" * 60)

preview = BladeSQS(
    phases_dict=phases["PHASE1"],
    sqsgen_levels=sqsgen_levels,
    level=5,
    len_comp=3,
)
sqsgen_text, rndstr_text = preview.sqs_struct()
print("--- rndstr.skel ---")
print(rndstr_text)
print("--- sqsgen.in ---")
print(sqsgen_text)
