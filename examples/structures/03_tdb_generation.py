"""Example: running the SQS-to-TDB fitting pipeline with BladeTDBGen.

This script shows how to configure and run the CALPHAD database fitting
step across multiple phase prototypes. BladeTDBGen wraps MaterialsFramework's
Sqs2tdb workflow across many chemical compositions, running one fit per
composition in sequence.

Key difference from the original BLADE: construction is side-effect-free.
The fitting loop only executes when gen.fit() is called explicitly.

Prerequisites:
    - SQS structures must already exist (run 02_sqs_generation.py first).
    - A CUDA-capable GPU is strongly recommended.

    python examples/03_tdb_generation.py
"""

from pathlib import Path

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_tdb_gen import BladeTDBGen

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# ------------------------------------------------------------------
# MLIP calculator — change mlip to swap potentials (see MaterialsFramework registry)
# ------------------------------------------------------------------
mlip        = "orb"                            # e.g. "grace", "mace", "uma", "chgnet", "orb"
mlip_kwargs = {"steps": 1000, "device": "cpu"} # kwargs forwarded to the calculator constructor

# ------------------------------------------------------------------
# Run flags
# ------------------------------------------------------------------
level = 2
skip_existing = False

tdb_params = {
    "fmax": 1e-4,
    "verbose": True,
    "calculator": mlip,
    "calculator_kwargs": mlip_kwargs,
    "t_min": 298.15,
    "t_max": 10000.0,
    "sro": False,
    "bv": 1e-3,
    "phonon": False,
    "open_calphad": False,
    "track_trajectory": True,   # set False to skip relaxation_live.xyz
    "terms": None,
}

terms_in: dict | None = None
terms_in = {
    "ALLOY1": "1,0:1,0\n2,2:1,0\n",
    "ALLOY2": "1,0:1,0\n2,2:2,2\n",
    "ALLOY3": "1,0:1,0\n2,2:1,0\n"
}

mult_in: dict | None = None
mult_in = {
    "ALLOY1": "a=1",
    "ALLOY2": "a=1\tb=1",
}

# sublattice_map: dict | None = None
sublattice_map = {"ALLOY2": {"a": ["Cr", "Hf"], "b": ["Zr", "Ti"]}}

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
    "Cr": 3.52, "Hf": 4.11, "Mo": 3.96, "Nb": 4.30,
    "Ta": 4.29, "Ti": 4.10, "V":  3.80, "W":  4.01, "Zr": 4.15,
}
_active_f = [el for el in primary_elements if el in _FCC_A]
_avg_a_fcc = sum(_FCC_A[el] for el in _active_f) / len(_active_f)
print(f"FCC lattice estimate:      a={_avg_a_fcc:.4f} Å  (avg of {_active_f})")

# ------------------------------------------------------------------
# Phase prototypes and run list (must match what was used for SQS generation)
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

liquid = False

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
print(f"Fitting TDBs for {len(composition_list)} composition(s): {composition_list}")

# ------------------------------------------------------------------
# 2. Fit TDB databases
# ------------------------------------------------------------------
comps_dir = path1 / "Files" / "Comps"

gen = BladeTDBGen(
    phases=phase_list,
    phases_dict=phases,
    liquid=liquid,
    paths=paths,
    composition_list=composition_list,
    level=level,
    skip_existing=skip_existing,
    output_dir=comps_dir,
    terms_in=terms_in,
    mult_in=mult_in,
    sublattice_map=sublattice_map,
    tdb_params=tdb_params,
)
gen.fit()

print("\nTDB generation complete.")
print(f"Output TDB files are located in: {comps_dir}")
