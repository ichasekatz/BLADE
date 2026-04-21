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
phases = ["HEDB2"]
liquid = False
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "PhaseForge/PhaseForge/atat/data/sqsdb/"
path2 = path0 / "BLADE/BLADEO/"
level = 6
sqs_iter = 1000000

# Define elements and composition settings
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
rare_earths = ["Sc", "Y", "La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
transition_metals = ["Cr"]
system_size = 1
tm_element_range = [1, 1]
re_element_range = [0, 0]
allow_lower_order = True

# Define phases
phases = {}

phases["Cr2O3_CORUNDUM"] = {
    "a": 1,
    "b": 1,
    "c": 2.741347,
    "alpha": 90,
    "beta": 90,
    "gamma": 120,
    "coords": """
0.000000 0.000000 0.347500  a
0.666667 0.333333 0.680833  a
0.333333 0.666667 0.014167  a
0.000000 0.000000 0.847500  a
0.666667 0.333333 0.180833  a
0.333333 0.666667 0.514167  a
0.000000 0.000000 0.152500  a
0.666667 0.333333 0.485833  a
0.333333 0.666667 0.819167  a
0.000000 0.000000 0.652500  a
0.666667 0.333333 0.985833  a
0.333333 0.666667 0.319167  a
0.306000 0.000000 0.250000  O
0.972667 0.333333 0.583333  O
0.639333 0.666667 0.916667  O
0.306000 0.306000 0.750000  O
0.972667 0.639333 0.083333  O
0.639333 0.972667 0.416667  O
0.000000 0.306000 0.750000  O
0.666667 0.639333 0.083333  O
0.333333 0.972667 0.416667  O
0.694000 0.694000 0.250000  O
0.360667 0.027333 0.583333  O
0.027333 0.360667 0.916667  O
0.694000 0.000000 0.250000  O
0.360667 0.333333 0.583333  O
0.027333 0.666667 0.916667  O
0.000000 0.694000 0.250000  O
0.666667 0.027333 0.583333  O
0.333333 0.360667 0.916667  O
""",
}

phases["CrO_ROCKSALT"] = {
    "a": 1,
    "b": 1,
    "c": 1,
    "alpha": 60,
    "beta": 60,
    "gamma": 60,
    "coords": """
0.000000 0.000000 0.000000  a
0.500000 0.500000 0.500000  O
""",
}


phases["CrO2_RUTILE_META"] = {
    "a": 1,
    "b": 1,
    "c": 0.644,
    "alpha": 90,
    "beta": 90,
    "gamma": 90,
    "coords": """
0.000000 0.000000 0.000000  a
0.500000 0.500000 0.500000  a
0.305000 0.305000 0.000000  O
0.695000 0.695000 0.000000  O
0.805000 0.195000 0.500000  O
0.195000 0.805000 0.500000  O
""",
}
