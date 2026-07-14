"""BLADE workflow: mcsqs SQS generation for Cr-Hf-Zr FCC ternary alloy.

Single-sublattice FCC phase. mcsqs generates optimal SQS structures
for each composition level. BladeTDBGen handles MLIP relaxation and
CALPHAD fitting.

Run on HPC:
    nohup python -u examples/structures/tdb_gen_alloy_mcsqs.py > tdb_gen_alloy_mcsqs.log 2>&1 &
    tail -f tdb_gen_alloy_mcsqs.log
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from pycalphad import Database

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_sqsgen import BladeSQS
from blade.tools.blade_tdb_gen import BladeTDBGen
from blade.analysis.blade_visual import BladeVisualizer

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
level    = 6
sqs_iter = 1_000_000
run_sqs  = True
skip_existing_sqs = False

run_tdb  = True
skip_existing_tdb = False

# ------------------------------------------------------------------
# MLIP calculator
# ------------------------------------------------------------------
mlip        = "orb"
mlip_kwargs = {"steps": 1000, "device": "cpu"}

tdb_params = {
    "fmax": 1e-4,
    "verbose": True,
    "calculator": mlip,
    "calculator_kwargs": mlip_kwargs,
    "t_min": 298.15,
    "t_max": 6000.0,
    "sro": False,
    "bv": 1e-3,
    "phonon": False,
    "open_calphad": False,
    "terms": None,
}

terms_in: dict | None = {"ALLOY1": "1,0:1,0\n2,2:1,0\n3,0:1,0\n"}
mult_in:  dict | None = None
sublattice_map: dict | None = None
sqsgen_in: dict | None = None

# ------------------------------------------------------------------
# Elements
# ------------------------------------------------------------------
primary_elements   = ["Cr", "Hf", "Zr"]
secondary_elements: list[str] = []
primary_min = 3
primary_max = 3
secondary_min = 0
secondary_max = 0

_FCC_A = {"Cr": 3.52, "Hf": 4.11, "Zr": 4.15}
_avg_a = sum(_FCC_A[el] for el in primary_elements) / len(primary_elements)
print(f"FCC lattice estimate: a={_avg_a:.4f} Å (avg of {primary_elements})")

# ------------------------------------------------------------------
# Phase prototype  (FCC primitive cell, single variable site)
# a=b=c, alpha=beta=gamma=60° is the standard ATAT FCC primitive
# ------------------------------------------------------------------
phases: dict[str, dict] = {
    "ALLOY1": {
        "a": _avg_a, 
        "b": _avg_a, 
        "c": _avg_a,
        "alpha": 60, 
        "beta": 60, 
        "gamma": 60,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": "0.000000 0.000000 0.000000 a\n",
    }
}

phase_list = [
    # 3×2×2 primitive = 12 M-atoms; accommodates 1/3 fractions exactly
    {"generator_name": "ALLOY", "lattice": "ALLOY1", "supercell_size": (3, 2, 2)},
]

# ------------------------------------------------------------------
# mcsqs parameters
# Cutoff values are nearest-neighbor shell indices (BladeSQS resolves
# them to actual Å using BladeCutoff). "3": 0 disables triplets.
# ------------------------------------------------------------------
mcsqs_params = {
    "use_time":    True,
    "time":        60,          # seconds per sqsdb directory
    "2":           4,           # pair cutoff  = 4th-NN shell
    "3":           3,           # triplet cutoff = 3rd-NN shell
    "4":           0,           # quadruplet — disabled
    "wr":          20,
    "wn":          0.75,
    "wd":          1,
    "parallel_runs": 10,
}

# ------------------------------------------------------------------
# SQS composition levels
#
# Each composition can carry an optional "elements" list specifying exactly
# which elements occupy that site. When set, dir_filter in BladeTDBGen
# deletes all ATAT-generated permutation dirs that don't match, keeping
# only the intended element assignment. Omit "elements" to allow all perms.
#
#   Level 0: pure Cr endmember
#   Level 1: binary 50/50 CrHf
#   Level 2: CrHfZr 0.5/0.25/0.25 + CrHf 0.5
#   Level 3: CrHfZr 0.5/0.25/0.25 + binary 0.75/0.25 (CrHf 0.25 + Cr 0.25)
#   Level 4: CrHf 0.5/0.5 and CrHf 0.75/0.25
#   Level 5: CrHf 0.5/0.5 + CrHfZr 0.5/0.25/0.25
#   Level 6: binary 0.5/0.5 (via perms) + equal-thirds ternary
# ------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0, 0.0]], "letter": ["a"], "elements": [["Cr"]]},
    {"level": 1, "compositions": [[0.5, 0.5, 0.0]], "letter": ["a"], "elements": [["Cr", "Hf"]]},
    {"level": 2, "compositions": [[0.5, 0.25, 0.25], [0.5, 0.5, 0.0]], "letter": ["a"], "elements": [["Cr", "Hf", "Zr"], ["Cr", "Hf"]]},
    {"level": 3, "compositions": [[0.5, 0.25, 0.25], [0.75, 0.25, 0.0]], "letter": ["a"], "elements": [["Cr", "Hf", "Zr"], ["Cr", "Hf"]]},
    {"level": 4, "compositions": [[0.5, 0.5, 0.0], [0.75, 0.25, 0.0]], "letter": ["a"], "elements": [["Cr", "Hf"], ["Cr", "Hf"]]},
    {"level": 5, "compositions": [[0.5, 0.5, 0.0], [0.5, 0.25, 0.25]], "letter": ["a"], "elements": [["Cr", "Hf"], ["Cr", "Hf", "Zr"]]},
    {"level": 6, "compositions": [[0.5, 0.5, 0.0], [1/3, 1/3, 1/3]], "letter": ["a"], "elements": [["Cr", "Hf", "Zr"], ["Cr", "Hf", "Zr"]]},
]

# Auto-build composition_elements from sqsgen_levels "elements" keys.
# Maps (level, sorted-desc-frac-str) → list of elements for dir_filter.
def _frac_key(fracs: list[float]) -> str:
    s = sorted(fracs, reverse=True)
    while len(s) > 1 and round(s[-1], 8) == 0.0:
        s.pop()
    return ",".join(str(round(f, 6)) for f in s)

composition_elements: dict[tuple, list[str]] = {}
for _lvl in sqsgen_levels:
    if "elements" not in _lvl:
        continue
    for _comp, _els in zip(_lvl["compositions"], _lvl["elements"]):
        composition_elements[(_lvl["level"], _frac_key(_comp))] = _els

liquid    = False
comps_dir = path1 / "Files" / "Comps_alloy"

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

composition_list  = composer.generate_compositions()
unique_len_comps  = composer.get_systems()

print(f"Compositions ({len(composition_list)} total): {composition_list}")
print(f"System sizes: {unique_len_comps}")

# ------------------------------------------------------------------
# 2. Generate SQS structures via mcsqs
# ------------------------------------------------------------------
if run_sqs:
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
            sqs_gen.sqs_gen(phase=specific_phase, paths=paths, iter=sqs_iter, params=params)

# ------------------------------------------------------------------
# 3. Fit TDB databases
# ------------------------------------------------------------------
if run_tdb:
    gen = BladeTDBGen(
        phases=phase_list,
        phases_dict=phases,
        liquid=liquid,
        paths=paths,
        composition_list=composition_list,
        level=level,
        skip_existing=skip_existing_tdb,
        output_dir=comps_dir,
        terms_in=terms_in,
        mult_in=mult_in,
        sublattice_map=sublattice_map,
        tdb_params=tdb_params,
    )
    gen.composition_elements = composition_elements
    gen.fit()

# ------------------------------------------------------------------
# 4. Plot Gibbs energy and phase diagrams
# ------------------------------------------------------------------
_coords = phases[phase_list[0]["lattice"]]["coords"]
_labels = [ln.strip().split()[-1] for ln in _coords.strip().splitlines() if ln.strip()]
_fixed  = [l for l in _labels if not (len(l) == 1 and l.islower())]
remove_elements = set(_fixed)
fixed_species   = {el: count / len(_labels) for el, count in Counter(_fixed).items()}

filt_comp_list = [
    [el for el in comp if el not in remove_elements]
    for comp in composition_list
]

viz = BladeVisualizer()

for comp, comp_filt in zip(composition_list, filt_comp_list):
    comp_name  = "".join(comp_filt)
    comp_dir   = comps_dir / comp_name
    if not comp_dir.exists():
        continue
    phase_name = f"{phase_list[0]['lattice']}_{len(comp_filt)}"
    for tdb_file in comp_dir.glob("*.tdb"):
        tdb = Database(str(tdb_file))
        viz.plot_gibbs_energy(
            tdb=tdb, metals=comp_filt, phase=phase_name,
            fixed_species=fixed_species,
            temperatures=[300, 1000, 2000, 3000],
            output_path=comp_dir / f"{comp_name}_Gibbs_Energy.png",
        )
