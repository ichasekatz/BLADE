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
import itertools
import csv
import numpy as np
from pycalphad import Database, equilibrium, variables as v

from blade.tools.blade_compositions import BladeCompositions
from blade.analysis.blade_visual import BLADEVisualizer
from blade.analysis.blade_volume import BLADEVolume

from pymatgen.io.vasp import Poscar
from ase.io import read

from materialsframework.calculators import GraceCalculator
from materialsframework.analysis.cte import CTEAnalyzer
from materialsframework.analysis.phono3py import Phono3pyAnalyzer
from materialsframework.analysis.formation_energy import FormationEnergyAnalyzer


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
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
transition_metals = ["Ti", "Cr", "W"]

system_size = 2
tm_element_range = [2, 2]
re_element_range = [0, 0]
allow_lower_order = True

liquid = False

# Define phases
phases = {}

phases["HEDB1"] = {
    "a": 4.58,
    "b": 4.58,
    "c": 5.04,
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





# PHASE_DIAGRAM_SYSTEM_SIZE = 3

# # Phase diagram plotting function
# def plot(tdb, elements, phases, file, comp):
#     if len(elements) >= PHASE_DIAGRAM_SYSTEM_SIZE:
#         return
#     fig = plt.figure(figsize=(9, 7))
#     axes = fig.gca()

#     phase_list = []
#     for phase in phases:
#         phase_list += [f"{phase}"]
#     if liquid:
#         phase_list.append("LIQUID")
#     print("phase_list:", phase_list)
#     print("db phases:", sorted(tdb.phases.keys()))

#     binplot(
#         tdb,
#         elements,
#         phase_list,
#         {v.X(elements[0]): (0, 1, 0.02), v.T: (1, 6000, 10), v.P: 101325, v.N: 1},
#         plot_kwargs={"ax": axes},
#     )
#     plt.tight_layout()
#     plt.savefig(f"{file}_Phase_Diagram.png", dpi=300)


# # Optimize SQS structures, generate TDB files, and plot phase diagrams for every composition
# for comp in composition_list:
#     file = Path(path2) / "".join(comp)
#     os.chdir(file)
#     elements = [el.upper() for el in comp]
#     file_names = ["_".join(p) for p in itertools.permutations(elements)]
#     for files in file_names:
#         if Path(f"{files}.tdb").is_file():
#             tdb = Database(f"{files}.tdb")
#             plot(tdb, elements, phases, files, comp)

# # Combine generated phase diagrams into one image
# viz = BLADEVisualizer()
# pngs = []
# for comp in composition_list:
#     comp_dir = Path(path2) / "".join(comp)
#     pngs.extend(comp_dir.glob("*_Phase_Diagram.png"))
# if pngs:
#     viz.phase_diagram(pngs, save=Path(f"{path2}{''.join('Combined_Phase_Diagrams.png')}"))
#     print("Combined phase diagrams saved as Combined_Phase_Diagrams.png")


# # Combine generated POSCAR files into one image for each phase in every composition
# for comp in composition_list:
#     comp_dir = Path(path2) / "".join(comp)
#     if not comp_dir.exists():
#         continue
#     for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
#         poscars = sorted(phase_dir.glob("sqs_lev=*/POSCAR"))
#         if not poscars:
#             continue
#         out = comp_dir / f"Combined_POSCARs_{''.join(comp)}_{phase_dir.name}.png"
#         viz.poscar(poscars, save=out)
#         print(f"Saved combined POSCARs → {out}")






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

from pathlib import Path
import csv
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pycalphad import Database, equilibrium
from pycalphad import variables as v


def compute_row(tdb_path, T=2000.0, P=101325.0, x=None):
    db = Database(str(tdb_path))
    comps = sorted([str(e) for e in db.elements if str(e).upper() != "VA"])
    phases = list(db.phases.keys())

    row = {
        "file": str(tdb_path),
        "T_K": T,
        "P_Pa": P,
        "elements": ",".join(comps),
    }

    conds = {v.T: T, v.P: P}

    if x is not None:
        x_full = {c: float(x.get(c, 0.0)) for c in comps}
        total = sum(x_full.values())
        if not np.isclose(total, 1.0, atol=1e-8):
            raise ValueError(f"Composition must sum to 1. Got {total} for {x_full}")

        for c in comps:
            row[f"x_{c}"] = x_full[c]

        # set N-1 independent composition variables
        for c in comps[:-1]:
            conds[v.X(c)] = x_full[c]

    # ask for NP too so eq phases can be extracted
    eq = equilibrium(db, comps, phases, conds, output=["GM", "MU", "NP"])

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
        npv = np.asarray(eq["NP"].values)

        present = []
        fracs = []
        for i, ph in enumerate(phase_axis):
            val = float(np.nanmax(npv[..., i]))
            if np.isfinite(val) and val > 1e-8:
                present.append(str(ph))
                fracs.append(val)

        row["eq_phases"] = ",".join(present)
        row["eq_phase_fractions"] = ",".join(f"{f:.6g}" for f in fracs)
    except Exception as e:
        row["eq_phases"] = ""
        row["eq_phase_fractions"] = ""
        row["phase_error"] = str(e)

    return row


def write_csv(rows, out_csv):
    outp = Path(out_csv)
    outp.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        raise RuntimeError(f"No rows were generated for {outp}")

    cols = sorted({k for r in rows for k in r.keys()})
    if not cols:
        raise RuntimeError(f"No columns found for {outp}")

    with outp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def composition_grid(elements, n_div=20):
    n = len(elements)
    for counts in itertools.product(range(n_div + 1), repeat=n):
        if sum(counts) != n_div:
            continue
        yield {el: c / n_div for el, c in zip(elements, counts)}


def metal_grid_with_fixed_boron(comps, boron_fraction=2/3, n_div=20):
    boron_names = [c for c in comps if c.upper() == "B"]
    if len(boron_names) != 1:
        raise ValueError(f"Expected exactly one boron component in {comps}, found {boron_names}")

    boron_name = boron_names[0]
    metals = [c for c in comps if c.upper() != "B"]

    if len(metals) == 0:
        raise ValueError("No non-B elements found after fixing boron.")

    metal_total = 1.0 - float(boron_fraction)
    if metal_total < 0:
        raise ValueError("boron_fraction must be <= 1")

    for counts in itertools.product(range(n_div + 1), repeat=len(metals)):
        if sum(counts) != n_div:
            continue

        x = {boron_name: float(boron_fraction)}
        for m, c in zip(metals, counts):
            x[m] = metal_total * c / n_div

        yield x


# -----------------------------
# use a NEW variable name here
# -----------------------------
tdb_summary_csv = path2 / "tdb_summary.csv"

T_list = [1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
boron_fraction = 2.0 / 3.0
n_div = 20

rows = []

if not tdb_paths:
    raise RuntimeError("tdb_paths is empty. No TDB files were found to evaluate.")

for p in tdb_paths:
    db = Database(str(p))
    comps = sorted([str(e) for e in db.elements if str(e).upper() != "VA"])

    has_boron = any(c.upper() == "B" for c in comps)

    # IMPORTANT: make this a list, not a generator
    if has_boron:
        x_list = list(
            metal_grid_with_fixed_boron(
                comps,
                boron_fraction=boron_fraction,
                n_div=n_div,
            )
        )
    else:
        x_list = list(
            composition_grid(
                comps,
                n_div=n_div,
            )
        )

    if not x_list:
        print(f"No compositions generated for {p}")
        continue

    for T in T_list:
        for x in x_list:
            try:
                rows.append(compute_row(p, T=T, P=101325.0, x=x))
            except Exception as e:
                fail_row = {
                    "file": str(p),
                    "T_K": T,
                    "P_Pa": 101325.0,
                    "error": str(e),
                }
                for c in comps:
                    fail_row[f"x_{c}"] = x.get(c, np.nan)
                rows.append(fail_row)

write_csv(rows, tdb_summary_csv)
print(f"Wrote: {tdb_summary_csv}")

# ----------------------------
# read and rank most stable rows
# ----------------------------
df = pd.read_csv(tdb_summary_csv)

if df.empty:
    raise RuntimeError(f"{tdb_summary_csv} was created but contains no data.")

df = df.dropna(subset=["G_system"])

if df.empty:
    raise RuntimeError("No valid Gibbs energy rows were found after dropping NaNs.")

temps = sorted(df["T_K"].unique())
xcols = [c for c in df.columns if c.startswith("x_")]

for T in temps:
    print("\n==============================")
    print(f"TOP 10 MOST STABLE @ {T} K")
    print("==============================")

    sub = df[df["T_K"] == T].sort_values("G_system").head(10)

    for _, row in sub.iterrows():
        comps_present = []
        for c in xcols:
            val = row[c]
            if not pd.isna(val) and abs(val) > 1e-12:
                comps_present.append(f"{c}={val:.3f}")

        comp_string = ", ".join(comps_present)
        print(f"{comp_string}   |   G = {row['G_system']:.2f}")

# ----------------------------
# plot best compositions
# ----------------------------
data = {}
temps = sorted(df["T_K"].unique())

for T in temps:
    sub = df[df["T_K"] == T].sort_values("G_system").head(10)

    xcols = [c for c in df.columns if c.startswith("x_")]
    if len(xcols) == 2:
        comp_col = sorted(xcols)[1]
    else:
        comp_col = xcols[0]

    data[T] = [(row[comp_col], row["G_system"]) for _, row in sub.iterrows()]

plt.figure(figsize=(8, 6))
markers = ['o', 's', '^', 'D', 'P']

for i, (T, points) in enumerate(data.items()):
    points = sorted(points, key=lambda p: p[0])
    x = [p[0] for p in points]
    g = [p[1] for p in points]
    plt.scatter(x, g, s=70, marker=markers[i % len(markers)], label=f"{T:.0f} K")

plt.xlabel(comp_col.replace("x_", "") + " Composition")
plt.ylabel("Gibbs Free Energy (J/mol)")
plt.title("Composition Stability vs Temperature")
plt.xlim(0, 1)
plt.legend()
plt.tight_layout()
plt.savefig("/Users/chasekatz/Desktop/School/Research/BLADE/G_vs_xTa.png", dpi=300)
plt.show()



base_root = Path("/Users/chasekatz/Desktop/School/Research/BLADE")

cte_out_path = base_root / "cte_results.csv"
tc_out_path = base_root / "thermal_conductivity_results.csv"
fe_out_path = base_root / "formation_energy_results.csv"

# ----------------------------
# CTE settings
# ----------------------------
cte_temperatures = [400, 2000]
cte_steps = 20

# ----------------------------
# thermal conductivity settings
# ----------------------------
tc_is_relaxed = False
tc_distance = 0.01
tc_supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
tc_primitive_matrix = None
tc_phonon_supercell_matrix = None
tc_mesh = [10, 10, 10]
tc_is_lbte = False
tc_is_isotope = False
tc_conductivity_type = None
tc_boundary_mfp = None
tc_t_min = 300
tc_t_max = 1500
tc_t_step = 300
tc_log_level = 0

# ----------------------------
# formation energy settings
# ----------------------------
fe_is_relaxed = False

# ----------------------------
# analyzers
# ----------------------------
calc = GraceCalculator()
cte_analyzer = CTEAnalyzer(calculator=calc)
phono_analyzer = Phono3pyAnalyzer(calculator=calc)
fe_analyzer = FormationEnergyAnalyzer(calculator=calc)


# ----------------------------
# helpers
# ----------------------------
def ensure_parent(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

def csv_exists_and_has_content(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0

def append_row(path: Path, row_dict: dict, header_order: list[str]):
    ensure_parent(path)
    exists = csv_exists_and_has_content(path)
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header_order)
        if not exists:
            writer.writeheader()
        writer.writerow({k: row_dict.get(k, "") for k in header_order})

def extract_tc_data(tc):
    temps = None
    for attr in ["temperatures", "_temperatures"]:
        if hasattr(tc, attr):
            temps = np.array(getattr(tc, attr), dtype=float)
            break

    kappa = None
    for attr in ["kappa", "_kappa", "mode_kappa", "_mode_kappa", "kappa_RTA", "kappa_LBTE"]:
        if hasattr(tc, attr):
            kappa = getattr(tc, attr)
            break

    result = {
        "tc_method": "LBTE" if tc_is_lbte else "RTA",
        "tc_temperatures": "",
        "kappa_shape": "",
        "kappa_xx_WmK": "",
        "kappa_yy_WmK": "",
        "kappa_zz_WmK": "",
        "kappa_iso_WmK": "",
        "kappa_iso_first_WmK": "",
        "kappa_raw": "",
    }

    if temps is not None:
        result["tc_temperatures"] = ";".join(map(str, temps.tolist()))

    if kappa is None:
        result["kappa_raw"] = "Could not parse thermal conductivity object"
        return result

    kappa = np.array(kappa)
    result["kappa_shape"] = str(kappa.shape)

    if kappa.ndim == 2 and kappa.shape[1] >= 3:
        kxx = kappa[:, 0]
        kyy = kappa[:, 1]
        kzz = kappa[:, 2]
        kiso = (kxx + kyy + kzz) / 3.0

        result["kappa_xx_WmK"] = ";".join(map(str, kxx.tolist()))
        result["kappa_yy_WmK"] = ";".join(map(str, kyy.tolist()))
        result["kappa_zz_WmK"] = ";".join(map(str, kzz.tolist()))
        result["kappa_iso_WmK"] = ";".join(map(str, kiso.tolist()))
        result["kappa_iso_first_WmK"] = float(kiso[0]) if len(kiso) > 0 else ""

    elif kappa.ndim == 3 and kappa.shape[1:] == (3, 3):
        kxx = kappa[:, 0, 0]
        kyy = kappa[:, 1, 1]
        kzz = kappa[:, 2, 2]
        kiso = (kxx + kyy + kzz) / 3.0

        result["kappa_xx_WmK"] = ";".join(map(str, kxx.tolist()))
        result["kappa_yy_WmK"] = ";".join(map(str, kyy.tolist()))
        result["kappa_zz_WmK"] = ";".join(map(str, kzz.tolist()))
        result["kappa_iso_WmK"] = ";".join(map(str, kiso.tolist()))
        result["kappa_iso_first_WmK"] = float(kiso[0]) if len(kiso) > 0 else ""
    else:
        result["kappa_raw"] = str(kappa)

    return result

# dynamic per-temperature volume columns
cte_volume_columns = [f"volume_at_{int(T)}K" for T in cte_temperatures]

cte_header = [
    "system",
    "phase",
    "composition_folder",
    "volumetric_cte_per_k",
    "volumetric_cte_ppm_per_k",
    "cte_reference_temperature",
    "cte_reference_volume",
    "cte_temperatures",
    "cte_volumes",
    *cte_volume_columns,
    "cte_error",
]

tc_header = [
    "system",
    "phase",
    "composition_folder",
    "tc_method",
    "tc_temperatures",
    "kappa_shape",
    "kappa_xx_WmK",
    "kappa_yy_WmK",
    "kappa_zz_WmK",
    "kappa_iso_WmK",
    "kappa_iso_first_WmK",
    "kappa_raw",
    "tc_error",
]

fe_header = [
    "system",
    "phase",
    "composition_folder",
    "formation_energy_eV_per_atom",
    "fe_error",
]

# ----------------------------
# main loop
# ----------------------------
for composition in composition_list:
    system_name = "".join(el.title() for el in composition)
    system_dir = base_root / system_name

    if not system_dir.exists():
        print(f"Missing system folder: {system_dir}")
        continue

    for phase_dir in sorted(system_dir.iterdir()):
        if not phase_dir.is_dir():
            continue

        phase_name = phase_dir.name

        for comp_dir in sorted(phase_dir.iterdir()):
            if not comp_dir.is_dir():
                continue

            poscar_path = comp_dir / "POSCAR"
            if not poscar_path.exists():
                print(f"Missing POSCAR: {poscar_path}")
                continue

            try:
                poscar = Poscar.from_file(str(poscar_path))
                fixed_path = comp_dir / "POSCAR_fixed"
                poscar.write_file(str(fixed_path))
                struct = read(str(fixed_path))
            except Exception as e:
                print(f"Skipping {poscar_path} due to POSCAR/read error: {e}")
                continue

            # ---------------- CTE ----------------
            cte_row = {
                "system": system_name,
                "phase": phase_name,
                "composition_folder": comp_dir.name,
                "volumetric_cte_per_k": "",
                "volumetric_cte_ppm_per_k": "",
                "cte_reference_temperature": "",
                "cte_reference_volume": "",
                "cte_temperatures": "",
                "cte_volumes": "",
                "cte_error": "",
            }
            for col in cte_volume_columns:
                cte_row[col] = ""

            try:
                cte_res = cte_analyzer.calculate(
                    structure=struct,
                    temperatures=cte_temperatures,
                    steps=cte_steps,
                )
                cte = cte_res["cte"]

                cte_row["volumetric_cte_per_k"] = cte["volumetric_per_k"]
                cte_row["volumetric_cte_ppm_per_k"] = cte["volumetric_ppm_per_k"]
                cte_row["cte_reference_temperature"] = cte["reference_temperature"]
                cte_row["cte_reference_volume"] = cte["reference_volume"]
                cte_row["cte_temperatures"] = ";".join(map(str, cte_res["temperatures"]))
                cte_row["cte_volumes"] = ";".join(map(str, cte_res["volumes"]))

                # add explicit columns like volume_at_300K
                for T, V in zip(cte_res["temperatures"], cte_res["volumes"]):
                    key = f"volume_at_{int(float(T))}K"
                    if key in cte_row:
                        cte_row[key] = V

            except Exception as e:
                cte_row["cte_error"] = str(e)

            append_row(cte_out_path, cte_row, cte_header)
            print(f"Wrote CTE row for {system_name} / {phase_name} / {comp_dir.name}")

            # ---------------- thermal conductivity ----------------
            # tc_row = {
            #     "system": system_name,
            #     "phase": phase_name,
            #     "composition_folder": comp_dir.name,
            #     "tc_method": "",
            #     "tc_temperatures": "",
            #     "kappa_shape": "",
            #     "kappa_xx_WmK": "",
            #     "kappa_yy_WmK": "",
            #     "kappa_zz_WmK": "",
            #     "kappa_iso_WmK": "",
            #     "kappa_iso_first_WmK": "",
            #     "kappa_raw": "",
            #     "tc_error": "",
            # }

            # try:
            #     phono_res = phono_analyzer.calculate(
            #         structure=struct,
            #         is_relaxed=tc_is_relaxed,
            #         distance=tc_distance,
            #         supercell_matrix=tc_supercell_matrix,
            #         primitive_matrix=tc_primitive_matrix,
            #         phonon_supercell_matrix=tc_phonon_supercell_matrix,
            #         mesh=tc_mesh,
            #         is_lbte=tc_is_lbte,
            #         is_isotope=tc_is_isotope,
            #         conductivity_type=tc_conductivity_type,
            #         boundary_mfp=tc_boundary_mfp,
            #         t_min=tc_t_min,
            #         t_max=tc_t_max,
            #         t_step=tc_t_step,
            #         log_level=tc_log_level,
            #     )

            #     tc = phono_res["thermal_conductivity"]
            #     tc_row.update(extract_tc_data(tc))

            # except Exception as e:
            #     tc_row["tc_error"] = str(e)

            # append_row(tc_out_path, tc_row, tc_header)
            # print(f"Wrote TC row for {system_name} / {phase_name} / {comp_dir.name}")

            # ---------------- formation energy ----------------
            fe_row = {
                "system": system_name,
                "phase": phase_name,
                "composition_folder": comp_dir.name,
                "formation_energy_eV_per_atom": "",
                "fe_error": "",
            }

            try:
                fe_res = fe_analyzer.calculate(
                    structure=struct,
                    is_relaxed=fe_is_relaxed,
                )
                fe_row["formation_energy_eV_per_atom"] = fe_res["formation_energy"]

            except Exception as e:
                fe_row["fe_error"] = str(e)

            append_row(fe_out_path, fe_row, fe_header)
            print(f"Wrote FE row for {system_name} / {phase_name} / {comp_dir.name}")

print(f"Done. CTE CSV: {cte_out_path}")
print(f"Done. TC CSV: {tc_out_path}")
print(f"Done. FE CSV: {fe_out_path}")