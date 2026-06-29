"""Example: running the SQS-to-TDB fitting pipeline with BladeTDBGen.

This script shows how to configure and run the CALPHAD database fitting
step across multiple phase prototypes. BladeTDBGen wraps MaterialsFramework's
Sqs2tdb workflow across many chemical compositions, running one fit per
composition in sequence.

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
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# ------------------------------------------------------------------
# Run flags
# ------------------------------------------------------------------
level = 5
skip_existing = False

tdb_params = {
    "fmax": 1e-4,
    "verbose": True,
    "calculator": "cuda",
    "t_min": 298.15,
    "t_max": 10000.0,
    "sro": False,
    "bv": 1e-3,
    "phonon": False,
    "open_calphad": False,
    "terms": None,
}

terms_in: dict | None = None
terms_in = {"HEDB1": "1,0:1,0\n2,2:1,0\n"}

mult_in: dict | None = None
mult_in = {
    "FCC1": "a=1",
    "FCC2": "a=1\tb=1",
}

# sublattice_map: dict | None = None
sublattice_map = {"FCC2": {"a": ["Cr", "Hf"], "b": ["Zr", "Ti"]}}

# ------------------------------------------------------------------
# Elements and composition constraints
# ------------------------------------------------------------------
primary_elements = ["Hf", "Cr"]
secondary_elements: list[str] = []

primary_min = 2
primary_max = 2
secondary_min = 0
secondary_max = 0

# ------------------------------------------------------------------
# AlB2-type diboride lattice constants (Å, from DFT literature)
# ------------------------------------------------------------------
_DIBORIDE_A = {
    "Cr": 2.969, "Hf": 3.141, "Mo": 3.053, "Nb": 3.086,
    "Ta": 3.098, "Ti": 3.028, "V":  2.998, "W":  3.020, "Zr": 3.169,
}
_DIBORIDE_C = {
    "Cr": 3.066, "Hf": 3.470, "Mo": 3.169, "Nb": 3.269,
    "Ta": 3.227, "Ti": 3.228, "V":  3.057, "W":  3.137, "Zr": 3.530,
}
_active_d = [el for el in primary_elements if el in _DIBORIDE_A]
_avg_a_diboride = sum(_DIBORIDE_A[el] for el in _active_d) / len(_active_d)
_avg_c_diboride = sum(_DIBORIDE_C[el] for el in _active_d) / len(_active_d)
print(f"Diboride lattice estimate: a={_avg_a_diboride:.4f} Å  c={_avg_c_diboride:.4f} Å  (avg of {_active_d})")

# FCC lattice constants (Å, conventional cubic cell)
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
    "HEDB1": {
        "a": _avg_a_diboride,
        "b": _avg_a_diboride,
        "c": _avg_c_diboride,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.333333 0.666667 0.500000 B\n"
            "0.666667 0.333333 0.500000 B\n"
        ),
    },
    "FCC1": {
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
    "FCC2": {
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
}

phase_list = [
    {"generator_name": "HEDB", "lattice": "HEDB1", "supercell_size": (4, 3, 2)},
    {"generator_name": "FCC",  "lattice": "FCC1",  "supercell_size": (2, 2, 2)},
    {"generator_name": "FCC",  "lattice": "FCC2",  "supercell_size": (2, 2, 2)},
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
