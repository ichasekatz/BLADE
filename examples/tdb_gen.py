import os
from pathlib import Path
import itertools

import matplotlib.pyplot as plt
from pycalphad import Database, binplot
from pycalphad import variables as v

from blade.blade_compositions import BladeCompositions
from blade.blade_sqs import BladeSQS
from blade.blade_tdb_gen import BladeTDBGen
from blade.blade_visual import BLADEVisualizer

phases = ["HEDB1"]
liquid = True
liquid = False
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "PhaseForge/PhaseForge/atat/data/sqsdb/"
path2 = path0 / "BLADE/BLADE/"
level = 6
time = 30

# Specify elements and system size (Total # elements)
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
rare_earths = ["Sc", "Y", "La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
system_size = 2
tm_element_range = [2, 2]
re_element_range = [0, 0]
allow_lower_order = True

a = 1
b = 1
c = 1.63299
alpha = 90
beta = 90
gamma = 120

rndstr = """
0.166667 0.833333 0.5       B
0.833333 0.166667 0.5       B
0.250000 0.250000 0.5       B
0.750000 0.750000 0.5       B
0.600000 0.300000 0.5       B
0.300000 0.600000 0.5       B
0.000000 0.000000 1.000000  a
0.333333 0.666667 1.000000  a
0.666667 0.333333 1.000000  a
0.500000 0.000000 1.000000  a
0.000000 0.500000 1.000000  a
0.500000 0.500000 1.000000  a
0.250000 0.750000 1.000000  a
"""

sqsgen_levels = [
    """level=0         a=1""",
    """level=1         a=0.5,0.5""",
    """level=2         a=0.75,0.25""",
    """level=3         a=0.33333,0.33333,0.33333""",
    """level=4         a=0.5,0.25,0.25""",
    """level=5         a=0.875,0.125\nlevel=5         a=0.625,0.375""",
    """level=6         a=0.75,0.125,0.125""",
]

paths = [path0, path1, path2]
composition_settings = [transition_metals, rare_earths, system_size,
                        tm_element_range, re_element_range, allow_lower_order]


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


for phase in phases:
    sqs_gen = BladeSQS(a, b, c, alpha, beta, gamma, rndstr, sqsgen_levels, level)
    sqs_gen.sqs_gen(unique_len_comps, phase, path1, time)

BladeTDBGen(
    phases,
    liquid,
    paths,
    composition_list,
    level,
)

fraction = "0.75, 0.125, 0.125"
fractions = [float(x) for x in fraction.split(",")]
sqs_gen = BladeSQS(a, b, c, alpha, beta, gamma, rndstr, sqsgen_levels, level)
N_atoms, counts = sqs_gen.supercell_size(fractions)
print(N_atoms, counts)

PHASE_DIAGRAM_SYSTEM_SIZE = 3

def plot(tdb, elements, phases, file, comp):
    if len(elements) >= PHASE_DIAGRAM_SYSTEM_SIZE:
        return
    fig = plt.figure(figsize=(9, 7))
    axes = fig.gca()

    phase_list = []
    for phase in phases:
        phase_list += [f"{phase}_{len(comp)}"]
    if liquid:
        phase_list.append("LIQUID")

    # Compute the phase diagram and plot it
    binplot(
        tdb,
        elements,
        phase_list,
        {v.X(elements[0]): (0, 1, 0.02), v.T: (1, 4000, 10), v.P: 101325, v.N: 1},
        plot_kwargs={"ax": axes},
    )
    plt.tight_layout()
    plt.savefig(f"{file}_Phase_Diagram.png", dpi=300)

for comp in composition_list:
    file = f"{path2}{''.join(comp)}"
    os.chdir(file)
    elements = [el.upper() for el in comp]
    file_names = ["_".join(p) for p in itertools.permutations(elements)]
    for files in file_names:
        if Path(f"{files}.tdb").is_file():
            tdb = Database(f"{files}.tdb")
            plot(tdb, elements, phases, files, comp)

viz = BLADEVisualizer()
pngs = []
for comp in composition_list:
    comp_dir = Path(f"{path2}{''.join(comp)}")
    pngs.extend(comp_dir.glob("*_Phase_Diagram.png"))
if pngs:
    viz.phase_diagram(
        pngs,
        save = Path(f"{path2}{''.join('Combined_Phase_Diagrams.png')}")
    )
    print("Combined phase diagrams saved as Combined_Phase_Diagrams.png")



for comp in composition_list:
    comp_dir = Path(f"{path2}{''.join(comp)}")
    if not comp_dir.exists():
        continue
    for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
        poscars = sorted(phase_dir.glob("sqs_lev=*/POSCAR"))
        if not poscars:
            continue
        out = comp_dir / "Combined_POSCARs.png"
        viz.poscar(
            poscars,
            save=out
        )
        print(f"Saved combined POSCARs → {out}")
