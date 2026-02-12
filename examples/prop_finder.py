import itertools
import os
from pathlib import Path
import pandas as pd

from pycalphad import Database, binplot
from pycalphad import variables as v

from blade.blade_compositions import BladeCompositions
from blade.blade_volume import BLADEVolume

# Define phases, pathways, and SQS generation settings
phases = ["HEDB1"]
liquid = False
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "PhaseForge/PhaseForge/atat/data/sqsdb/"
path2 = path0 / "BLADE/BLADE/"
level = 6
sqs_iter = 1000000

# Define elements and composition settings
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
rare_earths = ["Sc", "Y", "La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
transition_metals = ["Hf", "Mo", "Cr"]
system_size = 2
tm_element_range = [2, 2]
re_element_range = [0, 0]
allow_lower_order = True

# Define phases
phases = {}

phases["HEDB1"] = {
    "a": 1,
    "b": 1,
    "c": 1.63299,
    "alpha": 90,
    "beta": 90,
    "gamma": 120,
    "coords": """
0.000000 0.000000 0.000000  a
0.333333 0.666667 0.500000  B
0.666667 0.333333 0.500000  B
""",
}

# Specify SQS composition levels
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

# Optimize SQS structures, generate TDB files, and plot phase diagrams for every composition
for comp in composition_list:
    file = f"{path2}{''.join(comp)}"
    os.chdir(file)
    elements = [el.upper() for el in comp]
    file_names = ["_".join(p) for p in itertools.permutations(elements)]
    for files in file_names:
        if Path(f"{files}.tdb").is_file():
            tdb = Database(f"{files}.tdb")


path0 = Path("/Users/chasekatz/Desktop/School/Research")
path2 = path0 / "BLADE"
vol = BLADEVolume()
all_dfs = []
for comp in composition_list:
    comp_tag = "".join(comp)
    comp_dir = path2 / f"BLADE{comp_tag}"
    if not comp_dir.exists():
        print(f"Skipping (not found): {comp_dir}")
        continue
    df = vol.scan_poscars(comp_dir)
    if df.empty:
        print(f"No readable POSCARs found under: {comp_dir}")
        continue
    all_dfs.append(df)

if not all_dfs:
    raise RuntimeError("No readable POSCAR files found for any composition in composition_list.")

df_all = pd.concat(all_dfs, ignore_index=True)
df_all = df_all.sort_values(["composition_folder", "phase_folder", "sqs_level", "poscar_path"]).reset_index(drop=True)

out_csv = path2 / "poscar_phase_volumes.csv"
df_all.to_csv(out_csv, index=False)
print(f"Wrote: {out_csv}")
