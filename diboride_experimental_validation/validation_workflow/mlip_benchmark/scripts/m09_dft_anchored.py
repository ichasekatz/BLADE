"""DFT-anchored scenario: each model's own Stage 3 free energy with ONLY its TaB2 B-vacancy energy replaced by
the published DFT value (Ektarawong et al., J. Phys. Mater. 2023: -0.389 eV, bcc Ta + alpha-B references, i.e. our
B-rich convention). Everything else (mixing, ZrB2 vacancy, environment slope, W, competitors) stays the model's own.
The local-environment energies shift consistently (s11 `perturbed`: E_v(n) += d_m n/6).
No NbB2 DFT number is available (Dahlqvist & Rosen 2022 report ordered B vacancies stable in NbB2 up to ~10 % of
B sites but no single-vacancy value), so ZrNb is reported with each model's own energies only.

Reports phase assemblage, compositions (x_M), B/M of each AlB2 phase and metal fractions at x = 0.5, plus
composition-temperature maps for the anchored ZrTa case.
Output: mlip_benchmark/results/m09_dft_anchored_phases.csv, figures/fig_m09_dft_anchored_ZrTa.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

MB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MB / "scripts"))
sys.path.insert(0, str(MB.parent / "scripts"))
from m06_stage3_thermo import MODELS, build_system, classify, energies  # noqa: E402
from s11_vacancy_sensitivity import perturbed  # noqa: E402

DFT_TAB2_VB = -0.389


def describe(ph: pd.DataFrame) -> dict:
    al = ph[ph.phase == "AlB2"].sort_values("x_M")
    oth = ph[ph.phase != "AlB2"]
    return {"AlB2_x": ";".join(f"{v:.2f}" for v in al.x_M), "AlB2_B_per_M": ";".join(f"{2 * (1 - v):.2f}" for v in al.y_Va),
            "AlB2_metal_frac": ";".join(f"{v:.2f}" for v in al.metal_fraction), "other": ";".join(oth.phase),
            "other_metal_frac": ";".join(f"{v:.3f}" for v in oth.metal_fraction)}


def main() -> None:
    e, E = energies()
    natoms = dict(zip(e.structure_id, e.n_atoms))
    em = pd.read_csv(MB / "results" / "m02_end_member_benchmarks.csv").set_index(["model", "end_member"])
    rows = []
    for sid in ("ZrTa", "ZrNb"):
        for m in MODELS:
            base = build_system(m, sid, E, natoms)
            scen = {"own": base}
            if sid == "ZrTa":
                d = DFT_TAB2_VB - em.loc[(m, "TaB2"), "E_f_VB_Brich_eV"]
                scen["DFT-anchored TaB2 V_B"] = perturbed(base, d_m=float(d))
                scen["DFT-anchored TaB2 V_B"].compounds = base.compounds
            for sname, s in scen.items():
                for t in (1500.0, 2000.0, 2500.0, 3000.0):
                    for bm in (2.0, 1.98, 1.95, 1.90, 1.80):
                        ph = s.equilibrium(0.5, bm, t)
                        rows.append({"system": sid, "model": m, "scenario": sname, "T_K": t, "B_per_metal": bm,
                                     "class": classify(ph, base.compounds), **describe(ph)})
    df = pd.DataFrame(rows)
    df.to_csv(MB / "results" / "m09_dft_anchored_phases.csv", index=False)
    pd.set_option("display.width", 250)
    a = df[(df.system == "ZrTa")]
    print(a.pivot_table(index=["scenario", "model"], columns=["B_per_metal"], values="class", aggfunc=lambda v: "/".join(sorted(set(v)))).to_string())
    print(df[(df.T_K == 2000.0) & (df.B_per_metal == 1.95)][["system", "model", "scenario", "class", "AlB2_x", "AlB2_B_per_M",
                                                             "AlB2_metal_frac", "other", "other_metal_frac"]].to_string(index=False))
    plot(e, E, natoms, em)


def plot(e, E, natoms, em) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cls_col = {"HOMOGENEOUS": "#4daf4a", "AlB2_GAP": "#e41a1c", "AlB2+BORIDE": "#377eb8", "AlB2+BORIDE+AlB2_GAP": "#984ea3"}
    fig, axes = plt.subplots(2, len(MODELS), figsize=(3.3 * len(MODELS), 6.4), sharex=True, sharey=True)
    for j, m in enumerate(MODELS):
        base = build_system(m, "ZrTa", E, natoms)
        anch = perturbed(base, d_m=float(DFT_TAB2_VB - em.loc[(m, "TaB2"), "E_f_VB_Brich_eV"]))
        anch.compounds = base.compounds
        for i, (lab, s) in enumerate((("own energies", base), ("TaB2 V_B = DFT", anch))):
            ax = axes[i, j]
            for t in np.arange(1250.0, 3001.0, 250.0):
                for x in np.round(np.arange(0.1, 0.91, 0.1), 2):
                    c = classify(s.equilibrium(x, 1.95, t), base.compounds)
                    ax.scatter(x, t, marker="s", s=70, color=cls_col.get(c, "0.5"))
            ax.set_title(f"{m}\n{lab}", fontsize=8)
            if j == 0:
                ax.set_ylabel("T (K)")
            if i == 1:
                ax.set_xlabel("x (Ta)")
    from matplotlib.lines import Line2D

    fig.legend(handles=[Line2D([], [], marker="s", ls="", color=c, label=k) for k, c in cls_col.items()], loc="lower center", ncol=4,
               fontsize=8, frameon=False)
    fig.suptitle("ZrB$_2$–TaB$_2$ at B/M = 1.95: each model, before and after anchoring the TaB$_2$ B-vacancy energy to DFT", fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    (MB / "figures").mkdir(exist_ok=True)
    fig.savefig(MB / "figures" / "fig_m09_dft_anchored_ZrTa.png", dpi=170)


if __name__ == "__main__":
    main()
