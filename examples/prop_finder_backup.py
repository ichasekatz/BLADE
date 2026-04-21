import itertools
import os
from pathlib import Path
import time
import pandas as pd
from itertools import combinations
import argparse
import re
import matplotlib.pyplot as plt

from pycalphad import Database, binplot
from pycalphad import variables as v

import csv
import numpy as np
from pycalphad import Database, equilibrium, variables as v

from blade.tools.blade_compositions import BladeCompositions
from blade.analysis.blade_visual import BLADEVisualizer
from blade.analysis.blade_volume import BLADEVolume


# Define phases, pathways, and SQS generation settings
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path2 = path0 / "BLADE/"

# Define elements and composition settings
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "Nb", "Mo"]
rare_earths = ["Sc", "Y", "La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
rare_earths = ["Sc", "Lu", "Er", "Sm", "Ho", "Yb", "Tm", "La", "Y", "Dy", "Gd", "Nd", "Pr", "Eu", "Tb"]
rare_earths = ["Sc"]
transition_metals = ["Ni", "Re"]

system_size = 2
tm_element_range = [2, 2]
re_element_range = [0, 0]
allow_lower_order = False

liquid = False

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

tdb_paths = []
for comp in composition_list:
    file = Path(path2) / "".join(comp)
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


# Optimize SQS structures, generate TDB files, and plot phase diagrams for every composition
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






vol = BLADEVolume()
all_dfs = []
for comp in composition_list:
    comp_dir = Path(path2) / "".join(comp)
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





import matplotlib.pyplot as plt

data = {
1000: [(0.10,-929189.53),(0.15,-929131.51),(0.05,-928972.60),(0.20,-928873.42),(0.25,-928449.97)],
1500: [(0.20,-971214.97),(0.25,-971160.40),(0.15,-971038.66),(0.30,-970904.33),(0.10,-970579.41),(0.35,-970464.93),(0.40,-969853.80),(0.05,-969724.89),(0.45,-969078.16),(0.50,-968141.99)],
2000: [(0.30,-1019190.35),(0.25,-1019085.15),(0.35,-1019062.55),(0.20,-1018722.74),(0.40,-1018717.21),(0.45,-1018163.98),(0.15,-1018063.92),(0.50,-1017408.14),(0.10,-1017039.30),(0.55,-1016451.42)],
2500: [(0.35,-1071936.72),(0.40,-1071908.58),(0.30,-1071701.51),(0.45,-1071629.15),(0.25,-1071183.65),(0.50,-1071105.06),(0.20,-1070352.85),(0.55,-1070338.43),(0.60,-1069327.13),(0.15,-1069160.12)],
3000: [(0.45,-1128929.53),(0.40,-1128890.62),(0.50,-1128681.72),(0.35,-1128557.01),(0.55,-1128149.73),(0.30,-1127914.27),(0.60,-1127331.01),(0.25,-1126939.20),(0.65,-1126217.60),(0.20,-1125595.47)]
}

# plt.figure(figsize=(8,6))

# for T, points in data.items():
#     x=[p[0] for p in points]
#     g=[p[1] for p in points]
#     plt.scatter(x,g,label=f'{T} K')

# plt.xlabel('Hf Composition')
# plt.ylabel('Gibbs Free Energy (J/mol)')
# plt.title('Composition Stability vs Temperature')
# plt.legend()
# plt.tight_layout()
# plt.savefig("/Users/chasekatz/Desktop/School/Research/BLADE/G_vs_x.png", dpi=300)
# plt.show()

plt.figure(figsize=(8,6))

markers = ['o','s','^','D','P']

for i,(T, points) in enumerate(data.items()):
    points = sorted(points, key=lambda p: p[0])

    # convert HF → TA
    x = [1 - p[0] for p in points]   # <-- change here
    g = [p[1] for p in points]

    plt.scatter(x, g, s=70, marker=markers[i], label=f'{T} K')

plt.xlabel('Ta Composition')
plt.ylabel('Gibbs Free Energy (J/mol)')
plt.title('Composition Stability vs Temperature')
plt.xlim(0, 1)
plt.legend()
plt.tight_layout()
plt.savefig("/Users/chasekatz/Desktop/School/Research/BLADE/G_vs_xTa.png", dpi=300)
plt.show()