"""Example: generating SQS inputs and running ATAT mcsqs with BladeSQS.

This script shows how to configure and run the SQS generation step across
multiple phase prototypes. BladeSQS writes the ATAT input files (rndstr.skel,
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
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# ------------------------------------------------------------------
# Run flags
# ------------------------------------------------------------------
level = 2
skip_existing_sqs = False

sqsgen_in: dict | None = None
# sqsgen_in = {"ALLOY3": "level=0\ta=1\tb=1\nlevel=1\ta=0.5,0.5\tb=1\n"}

# ------------------------------------------------------------------
# Elements and composition constraints
# ------------------------------------------------------------------
primary_elements = ["Hf", "Cr", "Zr", "Ti"]
secondary_elements: list[str] = []

primary_min = 2
primary_max = 2
secondary_min = 0
secondary_max = 0

# ------------------------------------------------------------------
# FCC lattice constants (Å, conventional cubic cell)
# ------------------------------------------------------------------
_FCC_A = {
    "Cr": 3.52,
    "Hf": 4.11,
    "Mo": 3.96,
    "Nb": 4.30,
    "Ta": 4.29,
    "Ti": 4.10,
    "V": 3.80,
    "W": 4.01,
    "Zr": 4.15,
}
_active_f = [el for el in primary_elements if el in _FCC_A]
_avg_a_fcc = sum(_FCC_A[el] for el in _active_f) / len(_active_f)
print(f"FCC lattice estimate:      a={_avg_a_fcc:.4f} Å  (avg of {_active_f})")

# ------------------------------------------------------------------
# Phase prototypes
# ------------------------------------------------------------------
phases: dict[str, dict] = {
    # Phase with 4 variable a sites (e.g. Cr, Hf, Zr, Ti)
    "ALLOY1": {
        "a": _avg_a_fcc,
        "b": _avg_a_fcc,
        "c": _avg_a_fcc,
        "alpha": 90,
        "beta": 90,
        "gamma": 90,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.000000 0.500000 0.500000 a\n"
            "0.500000 0.000000 0.500000 a\n"
            "0.500000 0.500000 0.000000 a\n"
        ),
    },
    # Phase with 2 variable a sites (e.g. Cr, Hf) and 2 variable b sites (e.g. Zr, Ti)
    "ALLOY2": {
        "a": _avg_a_fcc,
        "b": _avg_a_fcc,
        "c": _avg_a_fcc,
        "alpha": 90,
        "beta": 90,
        "gamma": 90,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.000000 0.500000 0.500000 a\n"
            "0.500000 0.000000 0.500000 b\n"
            "0.500000 0.500000 0.000000 b\n"
        ),
    },
    # Phase with 2 variable a sites (e.g. Cr, Hf, Zr, Ti) and 2 fixed b sites (e.g. Zr)
    "ALLOY3": {
        "a": _avg_a_fcc,
        "b": _avg_a_fcc,
        "c": _avg_a_fcc,
        "alpha": 90,
        "beta": 90,
        "gamma": 90,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.000000 0.500000 0.500000 a\n"
            "0.500000 0.000000 0.500000 Zr\n"
            "0.500000 0.500000 0.000000 Zr\n"
        ),
    },
}

phase_list = [
    # supercell_size controls n_atoms = unit_cell_sites × product(supercell_size)
    # (4,3,2) → 72 atoms   (4,4,3) → 144   (6,6,2) → 216   (8,6,2) → 288
    {"generator_name": "ALLOY1", "lattice": "ALLOY1_1", "supercell_size": (2, 2, 2)},
    {"generator_name": "ALLOY2", "lattice": "ALLOY2_1", "supercell_size": (2, 2, 2)},
    {"generator_name": "ALLOY3", "lattice": "ALLOY3_1", "supercell_size": (2, 2, 2)},
]

# ------------------------------------------------------------------
# mcsqs run parameters
# ------------------------------------------------------------------
mcsqs_params = {
    "time": 30,  # seconds per sqsdb directory
    "cutoff_mode": "nn",  # "nn" = NN shell index (decimals OK), "angstrom" = direct Å
    "2": 5,  # pair cutoff  = 5th-nearest shell
    "3": 4,  # triplet cutoff
    "4": 3,  # quadruplet cutoff
    "wr": 20,
    "wn": 0.75,
    "wd": 1,
    "parallel_runs": 10,
}

# ------------------------------------------------------------------
# SQS composition levels
# ------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]], "letter": ["a", "b"]},
    {"level": 1, "compositions": [[0.5, 0.5]], "letter": ["a", "b"]},
    {"level": 2, "compositions": [[0.75, 0.25]], "letter": ["a", "b"]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]], "letter": ["a", "b"]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]], "letter": ["a", "b"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]], "letter": ["a", "b"]},
    {"level": 6, "compositions": [[0.75, 0.125, 0.125]], "letter": ["a", "b"]},
]

# ------------------------------------------------------------------
# 1. Generate compositions
# ------------------------------------------------------------------
composer = BladeCompositions(
    primary_elements=primary_elements,
    secondary_elements=secondary_elements,
    primary_min=primary_min,
    primary_max=primary_max,
    secondary_min=secondary_min,
    secondary_max=secondary_max,
)
composition_list = composer.generate_compositions()
unique_len_comps = composer.get_systems()

print(f"Compositions: {composition_list}")
print(f"System sizes: {unique_len_comps}")

# ------------------------------------------------------------------
# 2. Generate SQS structures
# ------------------------------------------------------------------
for specific_phase in phase_list:
    for len_comp in unique_len_comps:
        lattice = specific_phase["lattice"]
        sqs_gen = BladeSQS(
            phases_dict=phases[lattice],
            sqsgen_levels=sqsgen_levels,
            level=level,
            len_comp=len_comp,
            skip_existing_sqs=skip_existing_sqs,
            sqsgen_in=sqsgen_in.get(lattice) if sqsgen_in else None,
        )
        params = mcsqs_params | {"super_cell_size": specific_phase["supercell_size"]}
        sqs_gen.sqs_gen(phase=specific_phase, paths=paths, params=params)

# ------------------------------------------------------------------
# Preview only (no ATAT required) — inspect the generated file text
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Preview: rndstr.skel and sqsgen.in for len_comp=2, level=2")
print("=" * 60)

preview = BladeSQS(
    phases_dict=phases["ALLOY1_1"],
    sqsgen_levels=sqsgen_levels,
    level=level,
    len_comp=2,
)
sqsgen_text, rndstr_text = preview.sqs_struct()
print("--- rndstr.skel ---")
print(rndstr_text)
print("--- sqsgen.in ---")
print(sqsgen_text)
