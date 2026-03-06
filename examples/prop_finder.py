import itertools
import os
from pathlib import Path
import time
import pandas as pd
from itertools import combinations
import argparse
import re

from pycalphad import Database, binplot
from pycalphad import variables as v

from blade.blade_compositions import BladeCompositions
from blade.blade_volume import BLADEVolume

# Define phases, pathways, and SQS generation settings
phases = ["HEDB1"]
liquid = False
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "PhaseForge/PhaseForge/atat/data/sqsdb/"
path2 = path0 / "BLADE/TDBs/BLADE/"
level = 6
sqs_iter = 2000
use_time = [30, False]

# Define elements and composition settings
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "Nb", "Mo"]
rare_earths = ["Sc", "Y", "La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
rare_earths = ["Sc", "Lu", "Er", "Sm", "Ho", "Yb", "Tm", "La", "Y", "Dy", "Gd", "Nd", "Pr", "Eu", "Tb"]
rare_earths = ["Sc"]

["Zr", "Hf", "Ta", "Cr", "Ti", "Nb", "Mo"]

# transition_metals = ["Ni", "Re"]
system_size = 3
tm_element_range = [1, 3]
re_element_range = [0, 1]
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

# phases["FCC1"] = {
#     "a": 1,
#     "b": 1,
#     "c": 1,
#     "alpha": 60,
#     "beta": 60,
#     "gamma": 60,
#     "coords": """
# 0.000000 0.000000 0.000000  a
# """,
# }

# phases["HCP1"] = {
#     "a": 1,
#     "b": 1,
#     "c": 1.63299,
#     "alpha": 90,
#     "beta": 90,
#     "gamma": 120,
#     "coords": """
# 0.000000 0.000000 0.000000  a
# 0.666667 0.333333 0.500000  a
# """,
# }

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
tdb_paths = []
for comp in composition_list:
    file = f"{path2}{''.join(comp)}"
    if not Path(file).exists():
        print(f"Directory not found, skipping: {file}")
        continue
    os.chdir(file)
    elements = [el.upper() for el in comp]
    file_names = ["_".join(p) for p in itertools.permutations(elements)]
    for files in file_names:
        if Path(f"{files}.tdb").is_file():
            tdb = Database(f"{files}.tdb")
            tdb_paths.append(f"{file}/{files}.tdb")


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
df_all = df_all.sort_values(["composition_folder", "phase_folder", "sqs_level"]).reset_index(drop=True)

out_csv = path2 / "poscar_phase_volumes.csv"
df_all.to_csv(out_csv, index=False)
print(f"Wrote: {out_csv}")

import csv
import numpy as np
from pathlib import Path
from pycalphad import Database, equilibrium, variables as v


def compute_row(tdb_path, T=2000.0, P=101325.0, xB=0.5):
    db = Database(str(tdb_path))
    comps = sorted([str(e) for e in db.elements if str(e).upper() != "VA"])
    phases = list(db.phases.keys())

    row = {"file": str(tdb_path), "T_K": T, "P_Pa": P, "elements": ",".join(comps)}

    conds = {v.T: T, v.P: P}

    # ---- NEW PART (only change) ----
    if len(comps) == 2:
        A, B = comps[0], comps[1]
        conds[v.X(B)] = float(xB)
        row[f"x_{A}"] = 1.0 - float(xB)
        row[f"x_{B}"] = float(xB)
    # --------------------------------

    eq = equilibrium(db, comps, phases, conds, output="GM")

    try:
        row["G_system"] = float(np.nanmin(eq["GM"].values))
    except Exception:
        row["G_system"] = np.nan

    if "MU" in eq:
        for c in comps:
            try:
                row[f"mu_{c}"] = float(np.nanmin(eq["MU"].sel(component=c).values))
            except Exception:
                row[f"mu_{c}"] = np.nan
    else:
        for c in comps:
            row[f"mu_{c}"] = np.nan

    try:
        phase_axis = eq.coords["Phase"].values
        npv = eq["NP"].values
        present = []
        fracs = []
        for i, ph in enumerate(phase_axis):
            val = float(np.nanmax(npv[..., i]))
            if val > 1e-8:
                present.append(str(ph))
                fracs.append(val)
        row["eq_phases"] = ",".join(present)
        row["eq_phase_fractions"] = ",".join([f"{f:.6g}" for f in fracs])
    except Exception:
        row["eq_phases"] = ""
        row["eq_phase_fractions"] = ""

    return row

def write_csv(rows, out_csv):
    outp = Path(out_csv)
    outp.parent.mkdir(parents=True, exist_ok=True)
    cols = sorted({k for r in rows for k in r.keys()})
    with outp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

# ---- EDIT THESE ----
T_list = [1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
out_csv = path2 / "tdb_summary.csv"
# -------------------

x_grid = np.linspace(0.0, 1.0, 21)
rows = []
for p in tdb_paths:
    for T in T_list:
        for xB in x_grid:
            try:
                rows.append(compute_row(p, T=T, P=101325.0, xB=xB))
            except Exception:
                rows.append({"file": str(p), "T_K": T, "P_Pa": 101325.0, "xB": xB, "error": "failed"})

write_csv(rows, out_csv)



df = pd.read_csv(out_csv)

# Make sure Gibbs energy exists
df = df.dropna(subset=["G_system"])

results = []

# ---- rank for each temperature ----
temps = sorted(df["T_K"].unique())
xcols = [c for c in df.columns if c.startswith("x_")]
for T in temps:
    print("\n==============================")
    print(f"TOP 10 MOST STABLE @ {T} K")
    print("==============================")

    sub = df[df["T_K"] == T].sort_values("G_system").head(10)

    for _, row in sub.iterrows():

        # only keep X columns that have a real value
        comps = []
        for c in xcols:
            val = row[c]
            if not pd.isna(val) and abs(val) > 1e-12:
                comps.append(f"{c}={val:.3f}")

        comp_string = ", ".join(comps)

        print(f"{comp_string}   |   G = {row['G_system']:.2f} J/mol")