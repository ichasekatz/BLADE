"""Benchmark every MLIP against published DFT / experimental values found by the literature search.

References (validation_workflow/literature/*.csv; values copied from text actually read, access level noted there):
  R1  Ektarawong et al., J. Phys. Mater. 2023 (GGA-PAW VASP; bcc Ta + alpha-B references = our B-rich convention):
      dilute B vacancy in TaB2, dE = -0.389 eV per vacancy (full-text snippet).
  R2  same: TaB2-x with 0.167 <~ x <~ 0.25 on the Ta-B convex hull at 0 K (cluster expansion, 850 configurations).
  R3  same: DFT TaB2 a = 3.099 A, c = 3.326 A; c shrinks ~2.8 % at TaB1.75 (ordered), a changes < 0.3 %.
  R4  Dahlqvist & Rosen, Materialia 2022 (PBE): vacancies destabilise ZrB2 (group 4); B vacancies stabilise NbB2 and
      TaB2 (ordered B vacancies stable up to ~10 % (NbB2) and ~17 % (TaB2) of B sites); B vacancies reduce c in NbB2.
Model quantities: Stage 1 m02 tables + Stage 2 relaxed Y cells (TaB2 k = 8 spread is TaB1.833, x = 0.167).
Output: mlip_benchmark/results/m08_literature_benchmark.csv
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

MB = Path(__file__).resolve().parents[1]
WF = MB.parent
OUT = MB / "results"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
R1_TAB2_VB = -0.389


def cell(model: str, sid: str) -> list[float] | None:
    from ase.io import read

    if model == "orb":
        man = pd.concat([pd.read_csv(MB / "pilot_manifest.csv"), pd.read_csv(MB / "stage2_manifest.csv")]).set_index("structure_id")
        p = WF / "results" / "orb" / sid / "CONTCAR"
        if not p.exists():
            p = Path(man.loc[sid, "start_geometry"])
    else:
        p = MB / "results" / model / sid / "CONTCAR"
    return [float(v) for v in read(p, format="vasp").cell.lengths()] if p.exists() else None


def main() -> None:
    em = pd.read_csv(OUT / "m02_end_member_benchmarks.csv").set_index(["model", "end_member"])
    z = pd.read_csv(OUT / "m06_zeroK_Bdeficient_stability.csv")
    rows = []
    for m in MODELS:
        tab = em.loc[(m, "TaB2")]
        r = {"model": m, "TaB2_VB_dilute_eV": tab.E_f_VB_Brich_eV, "R1_DFT_eV": R1_TAB2_VB, "error_vs_R1_eV": tab.E_f_VB_Brich_eV - R1_TAB2_VB,
             "sign_TaB2_VB_matches_DFT": bool(tab.E_f_VB_Brich_eV < 0),
             "ZrB2_VB_positive (R4)": bool(em.loc[(m, "ZrB2")].E_f_VB_Brich_eV > 0),
             "a_TaB2_A": tab.a_A, "c_TaB2_A": tab.c_A, "a_err_vs_R3_pct": 100 * (tab.a_A / 3.099 - 1), "c_err_vs_R3_pct": 100 * (tab.c_A / 3.326 - 1)}
        for em_, sid in (("TaB2", "Y_TaB2_k8_spread"), ("NbB2", "Y_NbB2_k8_spread")):
            c0 = cell(m, f"host_{em_}")
            c8 = cell(m, sid)
            r[f"{em_}_c_change_pct_at_x0.167"] = 100 * (c8[2] / c0[2] - 1) if c0 and c8 else np.nan
            r[f"{em_}_a_change_pct_at_x0.167"] = 100 * (c8[0] / c0[0] - 1) if c0 and c8 else np.nan
        zt = z[(z.model == m) & (z.end_member == "TaB2") & z.cell.str.startswith("Y_TaB2_k8")]
        r["TaB2-x (x=0.167) below own Ta-B hull (R2 says yes)"] = bool((zt.e_above_own_hull_meV_atom < 0).any()) if len(zt) else None
        r["TaB2-x min e_hull_meV_atom"] = float(zt.e_above_own_hull_meV_atom.min()) if len(zt) else np.nan
        zn = z[(z.model == m) & (z.end_member == "NbB2") & z.cell.str.startswith("Y_NbB2")]
        r["NbB2-x below own hull (R4: ordered B vacancies stable up to ~10 %)"] = bool((zn.e_above_own_hull_meV_atom < 0).any()) if len(zn) else None
        rows.append(r)
    df = pd.DataFrame(rows)
    df["R3_note"] = "DFT c contraction ~2.8 % at x = 0.25 (ordered); linear scaling gives ~1.9 % at x = 0.167 (assumption)"
    df.to_csv(OUT / "m08_literature_benchmark.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.round(3).T.to_string())


if __name__ == "__main__":
    main()
