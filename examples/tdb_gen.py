import os
os.environ["OMP_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"

import torch
torch.set_num_threads(8)
torch.set_num_interop_threads(2)

import itertools
from pathlib import Path

import matplotlib.pyplot as plt
from pycalphad import Database, binplot
from pycalphad import variables as v

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_sqs import BladeSQS
from blade.tools.blade_tdb_gen import BladeTDBGen
from blade.analysis.blade_visual import BLADEVisualizer


# Define phases, pathways, and SQS generation settings
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path2 = path0 / "BLADE/"

level = 5
sqs_iter = 1000000
sqs = True
fit_tdb = True
skip_existing_tdb = False

# Define elements and composition settings
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
rare_earths = ["Sc", "Lu", "Er", "Sm", "Ho", "Yb", "Tm", "La", "Y", "Dy", "Gd", "Nd", "Pr", "Eu", "Tb"]
#transition_metals = ["Hf", "Cr"]
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "Nb", "Mo", "Ni", "Mn", "Ge", "Cu", "Fe", "Y"]
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
transition_metals = ["Ti", "Cr", "W"]
#transition_metals = ["Ni", "Re"]
system_size = 2
tm_element_range = [2, 2]
re_element_range = [0, 0]
allow_lower_order = False

# Define phases
phase_list = [
    {"generator_name": "FCC", "lattice": "FCC_A1", "supercell_size": (2, 3, 4)},
    {"generator_name": "HCP", "lattice": "HCP_A3", "supercell_size": (4, 4, 3)},
    {"generator_name": "HEDB", "lattice": "HEDB1", "supercell_size": (4, 4, 2)},
]

phase_list = [
    #{"generator_name": "FCC",  "lattice": "FCC_A1", "supercell_size": (3, 3, 4)},  # 144 atoms
    #{"generator_name": "HCP",  "lattice": "HCP_A3", "supercell_size": (6, 4, 3)},  # 144 atoms
    {"generator_name": "HEDB", "lattice": "HEDB1",  "supercell_size": (4, 4, 3)},  # 144 atoms
]

liquid = False

shell_weights = {
    1: 1.0,
    2: 1.0,
    3: 1.0,
    4: 1.0,
}

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

phases = {}

phases["HEDB1"] = {
    "a": 4.58,
    "b": 4.58,
    "c": 5.04,
    "alpha": 90,
    "beta": 90,
    "gamma": 120,
    "vectors":"""
1 0 0
0 1 0
0 0 1
""",
    "coords": """
0.000000 0.000000 0.000000  a
0.333333 0.666667 0.500000  B
0.666667 0.333333 0.500000  B
""",
}

phases["BCC_A2"] = {
    "a": 1,
    "b": 1,
    "c": 1,
    "alpha": 90,
    "beta": 90,
    "gamma": 90,
    "vectors": """
1 0 0
0 1 0
0 0 1
""",
    "coords": """
0.000000 0.000000 0.000000 a
0.500000 0.500000 0.500000 a
""",
}

phases["FCC_A1"] = {
    "a": 3.818376618407357,
    "b": 3.818376618407357,
    "c": 3.818376618407357,
    "alpha": 90,
    "beta": 90,
    "gamma": 90,
    "vectors": """
1 0 0
0 1 0
0 0 1
""",
    "coords": """
0.000000 0.000000 0.000000 a
0.500000 0.500000 0.000000 a
0.500000 0.000000 0.500000 a
0.000000 0.500000 0.500000 a
""",
}

phases["HCP_A3"] = {
    "a": 2.7,
    "b": 2.7,
    "c": 4.409081537009721,
    "alpha": 90,
    "beta": 90,
    "gamma": 120,
    "vectors": """
1 0 0
0 1 0
0 0 1
""",
    "coords": """
0.333333 0.666667 0.250000 a
0.666667 0.333333 0.750000 a
""",
}

# Specify SQS composition levels
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0]]},
    {"level": 1, "compositions": [[0.5, 0.5]]},
    {"level": 2, "compositions": [[0.75, 0.25]]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]]},
    {"level": 6, "compositions": [[0.75, 0.125, 0.125]]},
]

paths = [path0, path2]

# Generate compositions
compositions = BladeCompositions(
    transition_metals,
    rare_earths,
    system_size,
    tm_min=tm_element_range[0],
    tm_max=tm_element_range[1],
    re_min=re_element_range[0],
    re_max=re_element_range[1],
    allow_lower_order=allow_lower_order,
)

composition_list = compositions.generate_compositions()

unique_len_comps = compositions.get_systems()

print("Compositions: ", composition_list)
print("Total # compositions: ", len(composition_list))
print("Unique length compositions: ", unique_len_comps)

# Generate SQS structures for every composition system in each phase
if sqs:
    for specific_phase in phase_list:
        sqs_gen = BladeSQS(phases[specific_phase["lattice"]], sqsgen_levels, level)
        for len_comp in unique_len_comps:
            sqs_gen.sqs_gen(len_comp, specific_phase, path2, sqs_iter, shell_weights)

if fit_tdb:
    tdb_gen = BladeTDBGen(
        phases=phases,
        liquid=liquid,
        paths=paths,
        composition_list=composition_list,
        level=level,
        skip_existing=skip_existing_tdb,
    )

    tdb_gen.run_all_compositions(
        sqsgen_levels=sqsgen_levels,
        phase_dicts=phase_list,
        default_supercell_size=(2, 2, 2),
        params=tdb_params
    )

PHASE_DIAGRAM_SYSTEM_SIZE = 3

# Phase diagram plotting function
def plot(tdb, elements, phases, file, comp):
    if len(elements) >= PHASE_DIAGRAM_SYSTEM_SIZE:
        return
    fig = plt.figure(figsize=(9, 7))
    axes = fig.gca()

    phase_list = []
    for phase in phases:
        phase_list += [f"{phase}"]
    if liquid:
        phase_list.append("LIQUID")
    print("phase_list:", phase_list)
    print("db phases:", sorted(tdb.phases.keys()))

    binplot(
        tdb,
        elements,
        phase_list,
        {v.X(elements[0]): (0, 1, 0.02), v.T: (1, 6000, 10), v.P: 101325, v.N: 1},
        plot_kwargs={"ax": axes},
    )
    plt.tight_layout()
    plt.savefig(f"{file}_Phase_Diagram.png", dpi=300)


# Generate TDB files and plot phase diagrams for every composition
for comp in composition_list:
    file = Path(path2) / "".join(comp)
    os.chdir(file)
    elements = [el.upper() for el in comp]
    file_names = ["_".join(p) for p in itertools.permutations(elements)]
    for files in file_names:
        if Path(f"{files}.tdb").is_file():
            tdb = Database(f"{files}.tdb")
            plot(tdb, elements, phases, files, comp)

# Combine generated phase diagrams into one image
viz = BLADEVisualizer()
pngs = []
for comp in composition_list:
    comp_dir = Path(path2) / "".join(comp)
    pngs.extend(comp_dir.glob("*_Phase_Diagram.png"))
if pngs:
    viz.phase_diagram(pngs, save=Path(f"{path2}{''.join('Combined_Phase_Diagrams.png')}"))
    print("Combined phase diagrams saved as Combined_Phase_Diagrams.png")


# Combine generated POSCAR files into one image for each phase in every composition
for comp in composition_list:
    comp_dir = Path(path2) / "".join(comp)
    if not comp_dir.exists():
        continue
    for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
        poscars = sorted(phase_dir.glob("sqs_lev=*/POSCAR"))
        if not poscars:
            continue
        out = comp_dir / f"Combined_POSCARs_{''.join(comp)}_{phase_dir.name}.png"
        viz.poscar(poscars, save=out)
        print(f"Saved combined POSCARs → {out}")
