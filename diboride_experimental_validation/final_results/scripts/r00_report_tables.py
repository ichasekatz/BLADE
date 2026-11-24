"""Markdown tables for FINAL_REPORT.md, generated from the machine-readable results (no hand transcription).

Writes supplementary_data/report_tables.json: {table_name: markdown}. FINAL_REPORT.md contains markers
<!--TABLE:name--> ... <!--/TABLE:name--> and <!--CAPTION:figname--> ... <!--/CAPTION:figname--> that r00 fills in place
(python r00_report_tables.py fill).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
SD = FR / "supplementary_data"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
LAB = {"orb": "ORB", "mace-mpa": "MACE-MPA", "mace-omat": "MACE-OMAT", "mattersim": "MatterSim", "sevennet": "SevenNet"}
SYS = ["TiNb", "TiTa", "ZrNb", "ZrTa", "HfNb", "HfTa"]
SYSLAB = {s: f"{s[:2]}B₂–{s[2:]}B₂" for s in SYS}


def md(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        out.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(out)


def f(v, nd=2, sign=False):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "–"
    return f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"


def tables() -> dict[str, str]:
    T = {}
    em = pd.read_csv(SD / "x01_endmember_vacancy.csv").drop_duplicates(["end_member", "model"]).pivot(index="model", columns="end_member", values="E_v_Brich_eV")
    p = pd.read_csv(SD / "x01_params.csv")
    p = p[p.status == "ok"]
    z = pd.read_csv(SD / "z01_ce_fits.csv").set_index(["model", "end_member"])
    rows = []
    for m in MODELS:
        w = p[p.model == m]
        rows.append({"Model": LAB[m],
                     "E_v TiB₂ / ZrB₂ / HfB₂ (eV)": f"{em.loc[m, 'TiB2']:.2f} / {em.loc[m, 'ZrB2']:.2f} / {em.loc[m, 'HfB2']:.2f}",
                     "E_v NbB₂ (eV)": f(em.loc[m, "NbB2"]), "E_v TaB₂ (eV) [DFT −0.39]": f(em.loc[m, "TaB2"]),
                     "slope (eV per Group-V nn), 6 systems": f"{w.slope_eV_per_nn.min():.2f} to {w.slope_eV_per_nn.max():.2f}",
                     "W NbB₂ / TaB₂ (eV/fu)": f"{w[w.system.str.endswith('Nb')].W_MB2_eV_fu.iloc[0]:.2f} / {w[w.system.str.endswith('Ta')].W_MB2_eV_fu.iloc[0]:.2f}",
                     "lowest ordered TaB₂₋ₓ / NbB₂₋ₓ vs MB₂+B (meV/atom)": f"{1000 * z.loc[(m, 'TaB2'), 'min_dE_eV_fu'] / 3:+.1f} / {1000 * z.loc[(m, 'NbB2'), 'min_dE_eV_fu'] / 3:+.1f}",
                     "pair-CE LOO-CV TaB₂ / NbB₂ (meV/vac)": f"{1000 * z.loc[(m, 'TaB2'), 'loo_cv_eV_per_vacancy']:.0f} / {1000 * z.loc[(m, 'NbB2'), 'loo_cv_eV_per_vacancy']:.0f}"})
    T["models"] = md(pd.DataFrame(rows))

    s = pd.read_csv(FR / "SIX_SYSTEM_VALIDATION.csv")
    s["System"] = [SYSLAB[x] for x in SYS]
    T["stoich"] = md(s[["System", "stoich_dHmix_x0.5_meV_atom_DFT (Tang / Wen)", "stoich_dHmix_x0.5_meV_atom_MLIP_range",
                        "stoich_prediction (x=0.5, B/M=2, own, 1500-3000 K)"]].rename(columns={
        "stoich_dHmix_x0.5_meV_atom_DFT (Tang / Wen)": "ΔH_mix(0.5), DFT Tang LDA / Wen PBE (meV/atom)",
        "stoich_dHmix_x0.5_meV_atom_MLIP_range": "ΔH_mix(0.5), five MLIPs (meV/atom)",
        "stoich_prediction (x=0.5, B/M=2, own, 1500-3000 K)": "Assemblage at x = 0.5, B/M = 2.00, 1500–3000 K (own energetics)"}))

    k = pd.read_csv(SD / "x01_critical.csv")
    k = k[(k.x == 0.5) & (k.T_K == 2000.0)]
    short = {"mp-10255_Nb3B4": "Nb₃B₄", "mp-1102394_Nb5B6": "Nb₅B₆", "mp-2580_NbB": "NbB", "mp-1097_TaB": "TaB", "mp-10142_Ta3B4": "Ta₃B₄",
             "mp-28629_Ta5B6": "Ta₅B₆", "mp-1080021_Nb2B3": "Nb₂B₃"}
    rows = []
    for sid in SYS:
        for scen in ["own"] + (["DFT-anch"] if sid.endswith("Ta") else []):
            for m in MODELS:
                r = k[(k.system == sid) & (k.model == m) & (k.scenario == scen)]
                if r.empty:
                    continue
                r = r.iloc[0]
                prod = (f"AlB₂ (x={r.second_x_at_5:.2f}, B/M={r.second_BM_at_5:.2f})" if r.class_at_5 == "two AlB2"
                        else short.get(str(r.second_phase_at_5), str(r.second_phase_at_5)))
                d30 = "not reached ≤ 20 %" if r.delta_30 != r.delta_30 else f"{100 * r.delta_30:.1f}"
                rows.append({"System": SYSLAB[sid], "Scenario": scen, "Model": LAB[m], "δ_onset (%)": f"{100 * r.delta_onset:.1f}",
                             "δ₅ (%)": f"{100 * r.delta_5:.1f}", "δ₃₀ (%)": d30, "Product at δ₅": prod,
                             "Majority phase x at δ₅": f"{r.main_x_at_5:.2f}"})
    T["equilibria"] = md(pd.DataFrame(rows))

    r = pd.read_csv(SD / "y02_delta_to_match_ROI.csv")
    r = r[r.T_K.isin([2000.0, 2500.0])]
    rows = []
    for (sid, scen), g in r.groupby(["system", "scenario"], sort=False):
        row = {"System": SYSLAB[sid], "Scenario": scen, "ROI fraction (Ness Fig. 5)": f"{g.ROI.iloc[0]:.2f}"}
        for m in MODELS:
            q = g[g.model == m].set_index("T_K").delta_match
            row[LAB[m]] = " / ".join("n.r." if (t not in q or q[t] != q[t]) else f"{100 * q[t]:.1f}" for t in (2000.0, 2500.0))
        rows.append(row)
    order = {s: i for i, s in enumerate(SYSLAB.values())}
    T["roi"] = md(pd.DataFrame(rows).sort_values(["System", "Scenario"], key=lambda c: c.map(order) if c.name == "System" else c))

    y = pd.read_csv(SD / "y01_diffraction_six.csv")
    rows = []
    for sid in SYS:
        for scen in ["own"] + (["DFT-anch"] if sid.endswith("Ta") else []):
            g = y[(y.system == sid) & (y.scenario == scen)]
            two = g[g.class_at_delta5 == "two AlB2"]
            bor = g[g.class_at_delta5 == "metal-rich boride"]
            rows.append({"System": SYSLAB[sid], "Scenario": scen,
                         "two-AlB₂ models": len(two),
                         "Δ2θ(001) (°)": f"{two.split_001_deg.min():.2f}–{two.split_001_deg.max():.2f}" if len(two) else "–",
                         "Δ2θ(002) (°)": f"{two.split_002_deg.min():.2f}–{two.split_002_deg.max():.2f}" if len(two) else "–",
                         "Δ2θ(100) (°)": f"{two.split_100_deg.min():+.2f} to {two.split_100_deg.max():+.2f}" if len(two) else "–",
                         "Δc/c of B-deficient phase (%)": f"{two.c_change_phase2_vs_1_pct.min():+.1f} to {two.c_change_phase2_vs_1_pct.max():+.1f}" if len(two) else "–",
                         "boride models: strongest peak 2θ (°), distance to nearest AlB₂ peak (°)":
                             "; ".join(f"{short.get(b.boride, b.boride)} {b.boride_strongest_2th:.1f}, {b.boride_min_distance_to_AlB2_deg:.1f}" for b in bor.itertuples()) or "–"})
    T["xrd"] = md(pd.DataFrame(rows))

    fl = pd.read_csv(SD / "x02_flip_thresholds.csv")
    rows = []
    for (sid, scen), g in fl.groupby(["system", "scenario"], sort=False):
        rows.append({"System": SYSLAB[sid], "Scenario": scen,
                     "Reference product (2000 K, B/M 1.95), by model": "; ".join(f"{LAB[r.model]}: {r.reference_class}" for r in g.itertuples()),
                     "V_B(MB₂) shift that flips it (eV)": f"{g.d_m_flip.min():+.2f} to {g.d_m_flip.max():+.2f}",
                     "boride shift that flips it (eV/atom)": (f"{g.comp_flip.abs().min():.3f}–{g.comp_flip.abs().max():.3f} ({g.comp_flip.notna().sum()}/{len(g)} models)"
                                                              if g.comp_flip.notna().any() else "no flip within ±0.10")})
    T["flips"] = md(pd.DataFrame(rows))

    b = pd.read_csv(SD / "blade_before_after_summary.csv")
    b = b[b.status == "compared"]
    T["blade"] = md(pd.DataFrame({"Set": b.set, "Systems": b.systems.astype(int),
                                  "Systems with dG_hull change > 1 meV/fu": b["systems_with_any_dG_change_>1meV"].astype(int),
                                  "Strict single-phase rows before → after": [f"{int(x):,} → {int(y):,}" for x, y in zip(b.strict_rows_before, b.strict_rows_after)],
                                  "Lost (%)": [f"{100 * x / y:.0f}" for x, y in zip(b.strict_rows_lost, b.strict_rows_before)]}))
    pf = pd.read_csv(SD / "blade_pareto_front_change_2to5.csv")
    T["pareto"] = md(pd.DataFrame({"Metals": pf.n_metals, "Front before → after": [f"{a} → {c}" for a, c in zip(pf.front_before, pf.front_after)],
                                   "Left the front": pf.left_front, "Joined": pf.joined_front,
                                   "Ranked compositions before → after": [f"{a:,} → {c:,}" for a, c in zip(pf.ranked_before, pf.ranked_after)]}))
    v = pd.read_csv(SD / "v01_mass_balance.csv")
    g = pd.read_csv(SD / "v01_grid_convergence.csv")
    T["validation_text"] = (f"{len(v)} equilibria: maximum mass-balance error {max(v.x_err.max(), v.BM_err.max(), v.metal_sum_err.max()):.1e}; "
                            f"{int((g.class_default != g.class_fine).sum())} of {len(g)} states change class on a 2× finer x and y grid "
                            f"(all at a second-phase fraction < 0.1), maximum change in the second-phase fraction {g.abs_frac_diff.max():.3f}")
    return T


SUB = str.maketrans("0123456789+-−x", "₀₁₂₃₄₅₆₇₈₉₊₋₋ₓ")


def plain(t: str) -> str:
    """Convert matplotlib-style $_{..}$ / $_x$ subscripts in captions to Unicode (DOCX renders stray $ as math)."""
    def sub(m):
        body = m.group(1)
        return body.translate(SUB) if all(c in "0123456789+-−x" for c in body) else "_" + body
    t = re.sub(r"\$_\{([^}]*)\}\$", sub, t)
    t = re.sub(r"\$_([^$])\$", sub, t)
    return t.replace("$", "")


def captions() -> dict[str, str]:
    txt = plain((FR / "figures" / "CAPTIONS.md").read_text())
    out = {}
    for blk in txt.split("\n## ")[1:]:
        name, body = blk.split("\n", 1)
        out[name.strip()] = body.strip().split("\n\n")[0].strip()
    return out


def fill() -> None:
    T, C = tables(), captions()
    p = FR / "FINAL_REPORT.md"
    s = p.read_text()
    for k, v in T.items():
        s = re.sub(rf"<!--TABLE:{k}-->.*?<!--/TABLE:{k}-->", lambda _m, v=v, k=k: f"<!--TABLE:{k}-->\n{v}\n<!--/TABLE:{k}-->", s, flags=re.S)
    for k, v in C.items():
        s = re.sub(rf"<!--CAPTION:{k}-->.*?<!--/CAPTION:{k}-->", lambda _m, v=v, k=k: f"<!--CAPTION:{k}-->{v}<!--/CAPTION:{k}-->", s, flags=re.S)
    p.write_text(s)
    left = re.findall(r"«[^»]*»", s)
    print("filled; remaining placeholders:", left)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "fill":
        fill()
    else:
        json.dump(tables(), open(SD / "report_tables.json", "w"), indent=1)
        print("\n\n".join(f"### {k}\n{v}" for k, v in tables().items()))
