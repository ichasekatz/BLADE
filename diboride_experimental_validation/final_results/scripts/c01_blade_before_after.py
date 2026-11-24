"""Before/after comparison of BLADE stability and Pareto results for the binary-boride competitor correction.

before: Files/Stability_blade (2-3 metal), Files/Stability_blade_4metal, Files/Stability_blade_5metal, and
        Files/Pareto_2to4metal, Files/Pareto_2to5metal (original, read only)
after : bfix_rerun/Stability_23, Stability_4, Stability_5r (+ Stability_5s when finished), and the Pareto rerun
        bfix_rerun/Pareto_2to4metal (same command and settings as the original: --max-rank 1000, 200-2500 K)
Per (system, composition, T): dG_hull, ss_fraction, is_single_phase_strict, is_unique_homogeneous, dominant phase.
Per composition (Pareto): pareto_rank, strict_frac, mean_dG_hull, front membership.
Outputs: final_results/supplementary_data/blade_before_after_*.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BLADE = Path("/home/ichasekatz/BLADE/BLADE")
FR = BLADE / "diboride_experimental_validation" / "final_results"
SD = FR / "supplementary_data"
PAIRS = [("2-3 metal", BLADE / "Files" / "Stability_blade", BLADE / "bfix_rerun" / "Stability_23"),
         ("4 metal", BLADE / "Files" / "Stability_blade_4metal", BLADE / "bfix_rerun" / "Stability_4"),
         ("5 metal (ranked subset)", BLADE / "Files" / "Stability_blade_5metal", BLADE / "bfix_rerun" / "Stability_5r"),
         ("5 metal (stable-but-unranked rest)", BLADE / "Files" / "Stability_blade_5metal", BLADE / "bfix_rerun" / "Stability_5s")]
COLS = ["system", "comp_label", "T", "dG_hull", "ss_fraction", "is_single_phase_strict", "is_unique_homogeneous", "dominant_phase",
        "n_competitors"]
GROUP4 = {"Nb", "Ta", "Mo", "Cr", "Mn", "Fe", "Ni", "Ti", "Y"}  # metals whose binary borides were added


def compare_stability() -> pd.DataFrame:
    summ = []
    per_system = []
    for label, before, after in PAIRS:
        if not after.exists():
            summ.append({"set": label, "status": "not run / not finished"})
            continue
        files = sorted(after.glob("*_stability.parquet"))
        for f in files:
            b = before / f.name
            if not b.exists():
                continue
            A = pq.read_table(f, columns=COLS).to_pandas()
            B = pq.read_table(b, columns=COLS).to_pandas()
            m = B.merge(A, on=["system", "comp_label", "T"], suffixes=("_before", "_after"))
            per_system.append({
                "set": label, "system": f.name.split("_")[0], "n_rows": len(m),
                "competitors_before": int(m.n_competitors_before.max()), "competitors_after": int(m.n_competitors_after.max()),
                "max_increase_dG_hull_eV_fu": float((m.dG_hull_after - m.dG_hull_before).max()),
                "mean_increase_dG_hull_eV_fu": float((m.dG_hull_after - m.dG_hull_before).mean()),
                "rows_strict_lost": int((m.is_single_phase_strict_before & ~m.is_single_phase_strict_after).sum()),
                "rows_strict_gained": int((~m.is_single_phase_strict_before & m.is_single_phase_strict_after).sum()),
                "rows_unique_lost": int((m.is_unique_homogeneous_before & ~m.is_unique_homogeneous_after).sum()),
                "strict_rows_before": int(m.is_single_phase_strict_before.sum()), "strict_rows_after": int(m.is_single_phase_strict_after.sum()),
                "mean_ss_fraction_change": float((m.ss_fraction_after - m.ss_fraction_before).mean()),
                "rows_dominant_phase_changed": int((m.dominant_phase_before != m.dominant_phase_after).sum()),
            })
        print(label, "systems compared:", sum(1 for p in per_system if p["set"] == label), flush=True)
    ps = pd.DataFrame(per_system)
    ps.to_csv(SD / "blade_before_after_per_system.csv", index=False)
    for label, g in ps.groupby("set"):
        summ.append({"set": label, "status": "compared", "systems": len(g), "rows": int(g.n_rows.sum()),
                     "systems_with_new_competitors": int((g.competitors_after > g.competitors_before).sum()),
                     "systems_with_any_dG_change_>1meV": int((g.max_increase_dG_hull_eV_fu > 1e-3).sum()),
                     "strict_rows_before": int(g.strict_rows_before.sum()), "strict_rows_after": int(g.strict_rows_after.sum()),
                     "strict_rows_lost": int(g.rows_strict_lost.sum()), "strict_rows_gained": int(g.rows_strict_gained.sum()),
                     "unique_rows_lost": int(g.rows_unique_lost.sum())})
    s = pd.DataFrame(summ)
    s.to_csv(SD / "blade_before_after_summary.csv", index=False)
    return s


def compare_pareto(before: Path, after: Path, tag: str) -> pd.DataFrame | None:
    if not (after / "pareto_metrics.parquet").exists():
        return None
    cols = ["system", "comp_label", "n_metals", "pareto_rank", "strict_frac", "mean_dG_hull", "global_candidate_rank"]
    B = pq.read_table(before / "pareto_metrics.parquet", columns=cols).to_pandas()
    A = pq.read_table(after / "pareto_metrics.parquet", columns=cols).to_pandas()
    m = B.merge(A, on=["system", "comp_label", "n_metals"], suffixes=("_before", "_after"), how="outer")
    m["front_before"] = (m.pareto_rank_before == 1).fillna(False).astype(bool)
    m["front_after"] = (m.pareto_rank_after == 1).fillna(False).astype(bool)
    m.to_csv(SD / f"blade_pareto_before_after_{tag}.csv.gz", index=False)
    out = []
    for n, g in m.groupby("n_metals"):
        out.append({"pareto": tag, "n_metals": n, "front_before": int(g.front_before.sum()), "front_after": int(g.front_after.sum()),
                    "left_front": int((g.front_before & ~g.front_after).sum()), "joined_front": int((~g.front_before & g.front_after).sum()),
                    "ranked_before": int(g.pareto_rank_before.notna().sum()), "ranked_after": int(g.pareto_rank_after.notna().sum()),
                    "strict_frac_mean_change_ranked": float((g.strict_frac_after - g.strict_frac_before)[g.pareto_rank_before.notna()].mean())})
    r = pd.DataFrame(out)
    r.to_csv(SD / f"blade_pareto_front_change_{tag}.csv", index=False)
    # which systems lost front members
    lost = m[m.front_before & ~m.front_after].groupby("system").size().sort_values(ascending=False)
    lost.to_csv(SD / f"blade_pareto_front_losses_by_system_{tag}.csv", header=["compositions_left_front"])
    return r


def main() -> None:
    pd.set_option("display.width", 250)
    print(compare_stability().to_string(index=False))
    for before, after, tag in ((BLADE / "Files" / "Pareto_2to4metal", BLADE / "bfix_rerun" / "Pareto_2to4metal", "2to4"),
                               (BLADE / "Files" / "Pareto_2to5metal", BLADE / "bfix_rerun" / "Pareto_2to5metal", "2to5")):
        r = compare_pareto(before, after, tag)
        if r is not None:
            print(r.to_string(index=False))


if __name__ == "__main__":
    main()
