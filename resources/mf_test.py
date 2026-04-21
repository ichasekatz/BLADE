from materialsframework.calculators import GraceCalculator as Calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb

# # Initialize the wrapper
# sqs = Sqs2tdb(fmax=0.001, verbose=True, calculator=Calculator(device="cuda"))

# # Fit SQS for Ni-Re system on FCC_A1 and HCP_A3 lattices
# sqs.fit(species=["Ni", "Re"], lattices=["FCC_A1", "HCP_A3"], level=5)

import matplotlib.pyplot as plt
from pycalphad import Database, binplot
import pycalphad.variables as v

# Load database and choose the phases that will be considered
#db_alzn = sqs.database
#db_alzn = Database("/Users/chasekatz/Desktop/School/Research/BLADE/NI_RE.tdb")
db_alzn = Database("/Users/chasekatz/Desktop/School/Research/BLADE/Tests/NiRe_Before/BLADENiRe/NI_RE.tdb")
my_phases_alzn = ["FCC_A1", "HCP_A3"]

# Create a matplotlib Figure object and get the active Axes
fig = plt.figure(figsize=(9, 7))
axes = fig.gca()

# Compute the phase diagram and plot it on the existing axes using the `plot_kwargs={'ax': axes}` keyword argument
binplot(
    db_alzn,
    ["NI", "RE"],
    my_phases_alzn,
    {v.X("NI"): (0, 1, 0.02), v.T: (1, 6000, 10), v.P: 101325, v.N: 1},
    plot_kwargs={"ax": axes},
)
plt.tight_layout()
plt.savefig(
    "/Users/chasekatz/Desktop/School/Research/BLADE//Ni-Re_Phase_Diagram.png", dpi=300
)


fig = plt.figure(figsize=(9, 5))
axes = fig.gca()

binplot(
    db_alzn,
    ["NI", "RE"],
    my_phases_alzn,
    {v.X("NI"): (0, 1, 0.02), v.T: (1, 6000, 10), v.P: 101325, v.N: 1},
    plot_kwargs={"ax": axes},
)

# Increase all font sizes
axes.set_xlabel("X(NI)", fontsize=18)
axes.set_ylabel("Temperature (K)", fontsize=18)
axes.set_title("Ni-Re Phase Diagram", fontsize=20)

axes.tick_params(axis="both", which="major", labelsize=16)

# Make legend text larger if there is a legend
legend = axes.get_legend()
if legend is not None:
    legend.set_title(legend.get_title().get_text(), prop={"size": 16})
    for text in legend.get_texts():
        text.set_fontsize(14)

plt.tight_layout()
plt.savefig(
    "/Users/chasekatz/Desktop/School/Research/BLADE/Ni-Re_Phase_Diagram.png",
    dpi=600,
)
