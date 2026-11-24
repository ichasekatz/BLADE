"""Figures for the quantitative boron-loss analysis (reads p01 / p02 outputs only).

fig_p07_deltaT_maps_<system>_<scenario>.png  assemblage class over (T, delta) at x = 0.5 for each model,
                                              with 5 % and 30 % second-phase contours; hatched = vacancy model
                                              extrapolated (y > 0.125)
fig_p07_thresholds_x0.5.png                   delta_onset, delta_5, delta_30 vs T, all models and scenarios
fig_p07_fraction_vs_delta_2000K.png           second-phase metal fraction vs delta at 2000 K, x = 0.5, with the
                                              Ness lamellar light-phase area-fraction band (p01). Area fraction in
                                              lamellar ROIs is NOT a whole-sample phase fraction (see report)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
FIG = HERE / "figures"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
CLS = {"homogeneous": "#4daf4a", "two AlB2": "#e41a1c", "metal-rich boride": "#377eb8", "vacancy-ordered AlB2 cell": "#ff7f00",
       "ordered ZrMB4": "#984ea3"}


def cls_of(r) -> str:
    if r.n_phases <= 1 or r.second_identity == "none":
        return "homogeneous"
    return {"AlB2 (second solid solution)": "two AlB2"}.get(r.second_identity, r.second_identity)


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    raw = pd.read_csv(HERE / "results" / "p02_delta_scan_raw.csv")
    raw["cls"] = [cls_of(r) for r in raw.itertuples()]
    crit = pd.read_csv(HERE / "results" / "p02_critical_deficiency.csv")
    FIG.mkdir(exist_ok=True)
    for (sid, scen), g in raw[raw.x == 0.5].groupby(["system", "scenario"]):
        fig, axes = plt.subplots(1, len(MODELS), figsize=(3.4 * len(MODELS), 3.8), sharey=True)
        for ax, m in zip(axes, MODELS):
            h = g[g.model == m]
            if h.empty:
                ax.set_axis_off()
                continue
            ax.scatter(h.T_K, 100 * h.delta, c=[CLS.get(c, "0.5") for c in h.cls], marker="s", s=55)
            ext = h[h.validity == "extrapolated"]
            ax.scatter(ext.T_K, 100 * ext.delta, marker="x", color="k", s=18, lw=0.8)
            piv = h.pivot_table(index="delta", columns="T_K", values="second_metal_frac")
            if piv.shape[0] > 1 and piv.shape[1] > 1:
                cs = ax.contour(piv.columns, 100 * piv.index, piv.values, levels=[0.05, 0.30], colors=["k", "k"], linestyles=["--", "-"], linewidths=1)
                ax.clabel(cs, fmt={0.05: "5 %", 0.30: "30 %"}, fontsize=7)
            ax.set_title(m, fontsize=10)
            ax.set_xlabel("T (K)")
        axes[0].set_ylabel("boron deficiency δ (% of B)")
        h = [Line2D([], [], marker="s", ls="", color=c, label=k) for k, c in CLS.items()]
        h += [Line2D([], [], marker="x", ls="", color="k", label="vacancy model extrapolated (y > 0.125)"),
              Line2D([], [], ls="--", color="k", label="5 % second phase"), Line2D([], [], ls="-", color="k", label="30 % second phase")]
        fig.legend(handles=h, loc="lower center", ncol=4, fontsize=8, frameon=False)
        label = {"ZrNb": "ZrB$_2$–NbB$_2$", "ZrTa": "ZrB$_2$–TaB$_2$"}[sid]
        fig.suptitle(f"{label}, x = 0.5: equilibrium assemblage vs temperature and boron loss ({scen})", fontsize=11)
        fig.tight_layout(rect=(0, 0.12, 1, 0.94))
        fig.savefig(FIG / f"fig_p07_deltaT_maps_{sid}_{scen.replace(' ', '_').replace('(', '').replace(')', '').replace('<=', 'le')}.png", dpi=170)
        plt.close(fig)

    # thresholds vs T
    c5 = crit[crit.x == 0.5]
    scens = sorted(c5.scenario.unique())
    fig, axes = plt.subplots(2, len(scens) // 1 if len(scens) <= 4 else 4, figsize=(16, 7.5), squeeze=False)
    for j, scen in enumerate(scens[:4]):
        for i, sid in enumerate(("ZrNb", "ZrTa")):
            ax = axes[i, j]
            h = c5[(c5.scenario == scen) & (c5.system == sid)]
            if h.empty:
                ax.set_axis_off()
                continue
            for k, m in enumerate(MODELS):
                q = h[h.model == m].sort_values("T_K")
                ax.plot(q.T_K, 100 * q.delta_onset, ":", color=f"C{k}")
                ax.plot(q.T_K, 100 * q.delta_5, "--", color=f"C{k}")
                ax.plot(q.T_K, 100 * q.delta_30, "-", color=f"C{k}", label=m)
            ax.set_title(f"{sid} | {scen}", fontsize=9)
            ax.set_ylim(-0.5, 21)
            ax.set_xlabel("T (K)")
            ax.set_ylabel("δ (% of B)")
    axes[0, 0].legend(fontsize=7)
    fig.legend(handles=[Line2D([], [], ls=":", color="k", label="onset"), Line2D([], [], ls="--", color="k", label="5 % second phase"),
                        Line2D([], [], ls="-", color="k", label="30 % second phase")], loc="lower center", ncol=3, frameon=False)
    fig.suptitle("Boron deficiency needed at x = 0.5 (missing line = not reached within δ ≤ 20 %)", fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(FIG / "fig_p07_thresholds_x0.5.png", dpi=170)
    plt.close(fig)

    # fraction vs delta at 2000 K with Ness band
    p01 = pd.read_csv(HERE / "results" / "p01_ness_lamellae.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, sid, sysname in zip(axes, ("ZrNb", "ZrTa"), ("ZrB2-NbB2", "ZrB2-TaB2")):
        g = raw[(raw.system == sid) & (raw.x == 0.5) & (raw.T_K == 2000.0)]
        for k, m in enumerate(MODELS):
            for scen, ls in (("own (y<=0.125)", "-"), ("DFT-anchored (y<=0.125)", "--")):
                q = g[(g.model == m) & (g.scenario == scen)].sort_values("delta")
                if len(q):
                    ax.plot(100 * q.delta, q.second_metal_frac, ls, color=f"C{k}", label=f"{m} {scen.split(' ')[0]}" if sid == "ZrTa" or ls == "-" else None)
        e = p01[p01.system == sysname]
        lo, hi = e[["A_light_minus10", "A_light_plus10"]].min().min(), e[["A_light_minus10", "A_light_plus10"]].max().max()
        ax.axhspan(min(lo, 1 - hi), max(hi, 1 - lo) if False else hi, color="0.85", zorder=0,
                   label="Ness light-lamella area fraction\n(inside lamellar ROIs only)")
        ax.set_title(f"{sid}, x = 0.5, 2000 K (vacancy fraction capped at y ≤ 0.125)", fontsize=9)
        ax.set_xlabel("boron deficiency δ (% of B)")
    axes[0].set_ylabel("second-phase fraction (metal basis)")
    axes[1].legend(fontsize=6, ncol=2)
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(FIG / "fig_p07_fraction_vs_delta_2000K.png", dpi=170)


if __name__ == "__main__":
    main()
