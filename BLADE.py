from materialsframework.calculators import GraceCalculator as Calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb
from blade.blade_compositions import BladeCompositions
from blade.blade_sqs import BladeSQS
from blade.blade_phase_diagram import BladePhaseDiagram
import matplotlib.pyplot as plt
from pycalphad import Database, binplot, ternplot
from pycalphad import variables as v, binplot
import pycalphad.variables as v
import os

phase = 'HEDB1'
liquid = True
liquid = False
path0 = '/Users/chasekatz/Desktop/School/Research'
path1 = os.path.join(path0, 'PhaseForge/PhaseForge/atat/data/sqsdb/')
path2 = os.path.join(path0, 'MF_tests/p10/')
level = 4
time = 30

# Specify elements and system size (Total # elements)
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
rare_earths = ["Sc", "Y", "La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W", "Sc"]
system_size = 3
tm_element_range = [3, 3]
re_element_range = [0, 0]
allow_lower_order = True
#allow_lower_order = False

comp_gen = BladeCompositions(transition_metals, rare_earths, system_size, 
                                 tm_min=tm_element_range[0], tm_max=tm_element_range[1], 
                                 re_min=re_element_range[0], re_max=re_element_range[1], 
                                 allow_lower_order=allow_lower_order)

compositions = comp_gen.generate_compositions()

unique_len_comps = comp_gen.get_systems()

print("Total # compositions: ", compositions)
print("Unique length compositions: ", unique_len_comps)


a = 1
b = 1
c = 1.63299
alpha = 90
beta = 90
gamma = 120

rndstr = """
0.000000 0.000000 0.000000  a
0.333333 0.666667 0.000000  a
0.666667 0.333333 0.000000  a
0.500000 0.000000 0.000000  a
0.000000 0.500000 0.000000  a
0.500000 0.500000 0.000000  a
0.250000 0.750000 0.000000  a
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
"""level=6         a=0.75,0.125,0.125"""]

os.chdir(path0)

sqs_gen = BladeSQS(a, b, c, alpha, beta, gamma, rndstr, sqsgen_levels, level)

struct, sqs_struct = sqs_gen.sqs_struct()

fraction = "0.75, 0.125, 0.125"
fractions = [float(x) for x in fraction.split(",")]
N_atoms, counts = sqs_gen.supercell_size(fractions)
print(N_atoms, counts)

sqs_gen.sqs_gen(unique_len_comps, phase, path1, time)

def sqsfit_func(comp, phases, level):
    # Initialize the wrapper
    #sqs = Sqs2tdb(fmax=0.001, verbose=True, calculator=Calculator(device="cpu"), md_timestep=100, md_temperature=2000)
    sqs = Sqs2tdb(fmax=0.01, verbose=True, calculator=Calculator(device="cpu"))

    # Fit SQS
    sqs.fit(
        species=comp,
        lattices=phases,
        level=level,
    )

def plot(tdb, elements, phases, file):
    if len(elements) == 3:
        return
    fig = plt.figure(figsize=(9,7))
    axes = fig.gca()

    # Compute the phase diagram and plot it
    binplot(
        tdb, elements, phases,
        {v.X(elements[0]): (0, 1, 0.02), v.T: (1, 4000, 10), v.P: 101325, v.N: 1},
        plot_kwargs={'ax': axes}
    )
    plt.tight_layout()
    plt.savefig(f'{file}_Phase_Diagram.png', dpi=300)

phasediagram = BladePhaseDiagram(liquid, path2, phase)

for comp in compositions:
    file, phases = phasediagram.directory(comp)
    os.chdir(file)
    sqsfit_func(comp, phases, level)

    file_names, elements = phasediagram.file_names(comp)
    # Load database and choose the phases that will be considered

    for files in file_names:
        if os.path.isfile(f'{files}.tdb'):
            tdb = Database(f'{files}.tdb')
            #plot(tdb, elements, phases, files)

    # Change back to main folder
    os.chdir(path0)