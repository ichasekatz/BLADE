"""Publication figure set (16 figures) for the final report. Reads only CSV outputs and stored structures.

Visual conventions (all figures):
  experiment          black squares / grey bands         published DFT          black stars / black dashed lines
  MLIP (own)          solid, colour per model            DFT-anchored           dashed, same colour
  illustrative        grey dotted
Model colours: ORB blue, MACE-MPA orange, MACE-OMAT green, MatterSim red, SevenNet purple.
Each figure: figures/figNN_<name>.png (300 dpi) and .pdf; source data paths and captions in figures/CAPTIONS.md.
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
VAL = FR.parent
SD = FR / "supplementary_data"
FIG = FR / "figures"
sys.path.insert(0, str(FR / "scripts"))

MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
LABEL = {"orb": "ORB", "mace-mpa": "MACE-MPA", "mace-omat": "MACE-OMAT", "mattersim": "MatterSim", "sevennet": "SevenNet"}
COL = {"orb": "#1f77b4", "mace-mpa": "#ff7f0e", "mace-omat": "#2ca02c", "mattersim": "#d62728", "sevennet": "#9467bd"}
SYSTEMS = ["TiNb", "TiTa", "ZrNb", "ZrTa", "HfNb", "HfTa"]
SYSLAB = {s: f"{s[:2]}B$_2$–{s[2:]}B$_2$" for s in SYSTEMS}
CLS_COL = {"homogeneous": "#4daf4a", "two AlB2": "#e41a1c", "metal-rich boride": "#377eb8", "vacancy-ordered MB2-x": "#ff7f00",
           "ordered metal (ZrMB4-type)": "#984ea3"}
CAPTIONS: dict[str, str] = {}


def setup():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 7, "figure.dpi": 100,
                         "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
    return plt


def save(fig, name: str, caption: str, sources: list[str]):
    FIG.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight")
    CAPTIONS[name] = caption + "\n\n*Source data:* " + ", ".join(f"`{s}`" for s in sources)


def read(name):
    p = SD / name
    return pd.read_csv(p) if p.exists() else None


# ---------------------------------------------------------------------------------------------------------------
def fig01(plt):
    from PIL import Image

    s11 = pd.read_csv(VAL / "simulation_comparison" / "discrepancy" / "s11_otani_xrd_composition_spread.csv")
    p01 = pd.read_csv(VAL / "validation_workflow" / "phase_fraction_analysis" / "results" / "p01_ness_lamellae.csv")
    img = np.asarray(Image.open(VAL / "validation_workflow" / "phase_fraction_analysis" / "ness_images" / "main" / "fig-004.png").convert("L"))
    fig = plt.figure(figsize=(12, 6.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.28])
    ax = fig.add_subplot(gs[0, :2])
    ax.imshow(img, cmap="gray")
    from matplotlib.patches import Rectangle

    for r in p01.itertuples():
        ax.add_patch(Rectangle((r.x0, r.y0), r.x1 - r.x0, r.y1 - r.y0, fill=False, ec="yellow", lw=1.2))
    ax.set_axis_off()
    ax.set_title("(a) Ness et al. 2026 Fig. 5 (SE micrographs, 50 mol% melts) with lamellar ROIs used here", loc="left")
    ax = fig.add_subplot(gs[0, 2])
    for (system, cooling), g in s11.groupby(["system", "cooling"]):
        st = {"ZrB2-NbB2": ("s", "k"), "ZrB2-TaB2": ("D", "0.35"), "ZrB2-TiB2": ("^", "0.6"), "ZrB2-HfB2": ("o", "0.8")}[system]
        ax.plot(g.x_solute, g.dx_equiv_linear, marker=st[0], color=st[1], ls="-" if cooling.startswith("slow") else ":", ms=4,
                label=f"{system.replace('B2', 'B$_2$')} ({'slow' if cooling.startswith('slow') else 'quenched'})")
    ax.set_xlabel("x in Zr$_{1-x}$M$_x$B$_2$")
    ax.set_ylabel("equivalent composition spread Δx")
    ax.set_title("(b) Otani 2009: (10-11) XRD broadening\nexpressed as composition spread (exp.)", loc="left")
    ax.legend(fontsize=6)
    ax = fig.add_subplot(gs[1, :])
    ax.set_axis_off()
    boxes = [("Experiments\nNess 2026 lamellae; Otani, Eckert,\nVorotilo, Wen", 0.02, "#dddddd"),
             ("Stoichiometric theory\nDFT (Tang, Wen, user) +\nBLADE/ORB SQS: ΔH$_{mix}$ < 0", 0.215, "#dddddd"),
             ("5 MLIPs, 1 protocol\nvacancies, environments,\nordering, competitors", 0.41, "#cfe2f3"),
             ("Nonstoichiometric G(x,δ,T)\nCEF (M,M')(B,Va)$_2$ + borides\nclosed (δ) / open (μ$_B$)", 0.605, "#cfe2f3"),
             ("Predictions vs experiment\nfractions, compositions,\nXRD; BLADE implications", 0.80, "#fde9d9")]
    for txt, x, c in boxes:
        ax.add_patch(Rectangle((x, 0.02), 0.17, 0.96, transform=ax.transAxes, fc=c, ec="k"))
        ax.text(x + 0.085, 0.5, txt, transform=ax.transAxes, ha="center", va="center", fontsize=8)
    for x in (0.19, 0.385, 0.58, 0.775):
        ax.annotate("", xy=(x + 0.025, 0.5), xytext=(x, 0.5), xycoords="axes fraction", arrowprops=dict(arrowstyle="->"))
    ax.set_title("(c) Workflow", loc="left")
    fig.tight_layout()
    save(fig, "fig01_evidence_workflow", "Experimental evidence and computational workflow. (a) The six 50 mol% melt-pool micrographs of Ness et al. "
         "(2026, Fig. 5; secondary-electron contrast, published JPEG) with the lamella-only regions of interest analysed in Section 11. (b) Otani et al. "
         "(2009) peak broadening of the (10-11) reflection converted into an equivalent composition spread (Cu Kα1 assumed); ZrB2–TiB2 (known gap) "
         "and ZrB2–HfB2 (complete solubility) are controls. (c) Workflow.", ["validation_workflow/phase_fraction_analysis/results/p01_ness_lamellae.csv",
                                                                         "simulation_comparison/discrepancy/s11_otani_xrd_composition_spread.csv"])


def fig02(plt):
    mix = read("x01_mix_points.csv")
    s8 = pd.read_csv(VAL / "simulation_comparison" / "discrepancy" / "s8_six_system_validation_table.csv").set_index("system_id")
    user = pd.read_csv(VAL / "raw_data" / "digitized_figures" / "user_dft_ZrTaB2_mixing_enthalpy.csv")
    from x01_six_system_thermo import TANG

    fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), sharex=True)
    xs = np.linspace(0, 1, 101)
    for ax, sid in zip(axes.ravel(), SYSTEMS):
        L0, L1 = TANG[sid]
        ax.plot(xs, xs * (1 - xs) * (L0 + L1 * (2 * xs - 1)) / 3 / 96.485 * 1000, "k--", lw=1.2, label="DFT LDA (Tang 2025, RK)")
        ax.plot(0.5, s8.loc[sid, "dHmix_0.5_Wen_PBE_meV_per_atom"], "k*", ms=9, label="DFT PBE (Wen 2020)")
        if sid == "ZrTa":
            xu = user.columns[[c.lower().startswith("x") for c in user.columns]][0]
            yu = [c for c in user.columns if "meV" in c or "enthalpy" in c.lower()][0]
            ax.plot(user[xu], user[yu], "k^", mfc="none", ms=5, label="DFT PBE (user)")
        if mix is not None:
            for m in MODELS:
                g = mix[(mix.system == sid) & (mix.model == m)].sort_values("x")
                if len(g):
                    ax.plot(np.r_[0, g.x, 1], np.r_[0, g.dH_mix_meV_atom, 0], "-o", color=COL[m], ms=3, lw=1, label=LABEL[m])
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_title(SYSLAB[sid])
        ax.set_xlabel("x (Group-V fraction)")
        ax.set_ylabel("ΔH$_{mix}$ (meV/atom)")
    axes[0, 0].legend(fontsize=6)
    fig.tight_layout()
    save(fig, "fig02_mixing_enthalpy", "Stoichiometric mixing enthalpy of the six AlB2-type pseudobinaries (SQS, 0 K). Lines with circles: five MLIPs, each "
         "relative to its own end members (x = 0.25/0.5/0.75 SQS). Dashed: DFT-LDA Redlich–Kister fits (Tang et al. 2025, Table 2). Star: DFT-PBE "
         "(Wen et al. 2020). Triangles: user DFT-PBE (ZrB2–TaB2, digitized). Negative values at every composition in every DFT source and in all MLIPs "
         "imply complete miscibility of the stoichiometric solid.", ["supplementary_data/x01_mix_points.csv", "raw_data/digitized_figures/user_dft_ZrTaB2_mixing_enthalpy.csv"])


def fig03(plt):
    e = read("x01_endmember_vacancy.csv")
    ems = ["TiB2", "ZrB2", "HfB2", "NbB2", "TaB2"]
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    w = 0.15
    for i, m in enumerate(MODELS):
        g = e[e.model == m].drop_duplicates("end_member").set_index("end_member").reindex(ems)
        ax.bar(np.arange(5) + (i - 2) * w, g.E_v_Brich_eV, w, color=COL[m], label=LABEL[m])
    ax.plot(4, -0.389, "k*", ms=12, label="DFT (Ektarawong 2023)")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(range(5), [m.replace("B2", "B$_2$") for m in ems])
    ax.set_ylabel("E$_f$(V$_B$), dilute, B-rich limit (eV)")
    ax.legend(ncol=3, fontsize=7)
    ax.text(0.0, 1.02, "Group IV hosts", transform=ax.transAxes, fontsize=8)
    ax.text(0.62, 1.02, "Group V", transform=ax.transAxes, fontsize=8)
    save(fig, "fig03_vacancy_formation_energies", "Dilute boron-vacancy formation energy in the pure diborides at the boron-rich limit (μ_B = α-B), each model with its own "
         "references (144-atom cells, one vacancy). Star: published DFT for TaB2 (Ektarawong et al. 2023, GGA, same reference convention). "
         "All models agree that vacancies are costly in the Group IV diborides (2.6–3.3 eV) and much cheaper in NbB2 and TaB2; only ORB reproduces the "
         "negative DFT value for TaB2.", ["supplementary_data/x01_endmember_vacancy.csv"])


def fig04(plt):
    env = read("x01_vacancy_env.csv")
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), sharey=True)
    for ax, sid in zip(axes.ravel(), SYSTEMS):
        for m in MODELS:
            g = env[(env.system == sid) & (env.model == m)].sort_values("n_solute")
            if len(g):
                ax.plot(g.n_solute, g.E_v_Brich_eV, "-o", color=COL[m], ms=3.5, label=LABEL[m])
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_title(f"{SYSLAB[sid]}, x = 0.5 SQS")
        ax.set_xlabel(f"{sid[2:]} among 6 metal neighbours of the B site")
    for ax in axes[:, 0]:
        ax.set_ylabel("E$_f$(V$_B$), B-rich (eV)")
    axes[0, 0].legend(fontsize=6)
    fig.tight_layout()
    save(fig, "fig04_vacancy_vs_coordination", "Boron-vacancy formation energy versus local metal coordination in the equimolar SQS of all six systems (one vacancy per "
         "144-atom cell, all environments n = 0–6 present). The monotonic decrease with the number of Group-V neighbours is the vacancy–metal "
         "coupling that drives boron-loss-induced segregation.", ["supplementary_data/x01_vacancy_env.csv"])


def fig05(plt):
    vv = read("x01_vacancy_interactions.csv")
    ems = ["TiB2", "ZrB2", "HfB2", "NbB2", "TaB2"]
    arr = ["k2_spread", "k4_spread", "k8_spread", "k4_cluster", "k8_plane"]
    fig, axes = plt.subplots(1, 5, figsize=(13, 3.6), sharey=True)
    for ax, em in zip(axes, ems):
        for i, m in enumerate(MODELS):
            g = vv[(vv.end_member == em) & (vv.model == m)].set_index("arrangement").reindex(arr)
            ax.plot(np.arange(5) + (i - 2) * 0.08, g.interaction_eV_per_vacancy, "o", color=COL[m], ms=4, label=LABEL[m])
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(range(5), ["2 spread", "4 spread", "8 spread", "4 cluster", "8 plane"], rotation=45)
        ax.set_title(em.replace("B2", "B$_2$"))
    axes[0].set_ylabel("interaction energy (eV per vacancy)")
    axes[0].legend(fontsize=6)
    fig.tight_layout()
    save(fig, "fig05_vacancy_interactions", "Boron-vacancy–vacancy interaction energy per vacancy, [E(k vacancies) − E(host) − k·E1]/k, for spread, clustered and "
         "single-plane arrangements in 144-atom cells. Clustering is favoured in the Group IV diborides and disfavoured in NbB2 and TaB2 for every "
         "model; spread vacancies interact weakly.", ["supplementary_data/x01_vacancy_interactions.csv"])


def fig06(plt):
    o = read("x01_ordering.csv")
    h = read("z01_ordered_hull.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    ax = axes[0]
    cfg = o[o.config.isin(["qsqs1", "clay1", "clay2"])].pivot_table(index=["system", "model"], columns="config", values="dH_mix_meV_atom").reset_index()
    for k, sid in enumerate(SYSTEMS):
        for i, m in enumerate(MODELS):
            r = cfg[(cfg.system == sid) & (cfg.model == m)]
            if len(r):
                r = r.iloc[0]
                ax.plot(k + (i - 2) * 0.1, r.clay1 - r.qsqs1, "D", color=COL[m], ms=4)
                ax.plot(k + (i - 2) * 0.1, r.clay2 - r.qsqs1, "v", color=COL[m], ms=4, mfc="none")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(6), [SYSLAB[s] for s in SYSTEMS], rotation=40, fontsize=7)
    ax.set_ylabel("E(config) − E(random), meV/atom")
    ax.set_title("(a) metal ordering (◆ alternating c-layers)\nvs segregation (▽ AABB slab), x = 0.5")
    for j, em in enumerate(("TaB2", "NbB2")):
        ax = axes[1 + j]
        if h is None:
            ax.set_axis_off()
            continue
        for m in MODELS:
            g = h[(h.model == m) & (h.end_member == em)]
            if g.empty:
                continue
            ax.plot(g.B_per_M, g.dE_eV_fu * 1000 / 3, ".", color=COL[m], ms=3, alpha=0.5)
            hull = g[g.on_0K_hull].sort_values("B_per_M")
            low = g.groupby("B_per_M").dE_eV_fu.min().sort_index(ascending=False)
            ax.plot(np.r_[2.0, low.index], np.r_[0, low.values * 1000 / 3], "--", color=COL[m], lw=0.9, label=LABEL[m])
            if len(hull):
                ax.plot(np.r_[2.0, hull.B_per_M], np.r_[0, hull.dE_eV_fu * 1000 / 3], "-", color=COL[m], lw=2.2)
        ax.axhline(0, color="k", lw=0.6)
        if em == "TaB2":
            ax.axvspan(1.75, 1.833, color="0.85", zorder=0, label="DFT: TaB$_{2-x}$ on hull (Ektarawong)")
        ax.set_xlabel("B/M of ordered cell")
        ax.set_ylabel("ΔE vs MB$_2$ + α-B (meV/atom)")
        ax.set_title(f"({'b' if j == 0 else 'c'}) enumerated ordered {em.replace('B2', 'B$_2$')}$_{{-x}}$ (223 cells)")
        ax.invert_xaxis()
        ax.legend(fontsize=6)
    fig.tight_layout()
    save(fig, "fig06_ordering_vs_clustering", "Ordering versus clustering. (a) Energy of the metal-ordered alternating-c-layer arrangement (ZrMB4-type) and of a "
         "segregated AABB slab relative to a same-cell quasi-random arrangement at x = 0.5. (b, c) Formation energy of all 223 symmetry-distinct "
         "ordered B-vacancy arrangements (2–8 formula units, 1–2 vacancies) relative to stoichiometric MB2 plus α-B; dashed lines: each model's lowest-energy cell at each B/M; thick solid line: "
         "stable ordered phases (0 K convex hull below MB2 + B; ORB TaB2 only). Grey band: B/M range where DFT cluster expansion places ordered TaB2−x on the hull.", ["supplementary_data/x01_ordering.csv", "supplementary_data/z01_ordered_hull.csv"])


def fig07(plt):
    from x01_six_system_thermo import build, energies

    E, nat = energies()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (model, sid, scen) in zip(axes, (("orb", "ZrTa", "DFT-anch"), ("mace-mpa", "ZrNb", "own"))):
        s = build(model, sid, E, nat, scen)
        if s is None:
            ax.set_axis_off()
            continue
        T = 2000.0
        xs = np.linspace(0.001, 0.999, 241)
        ys = np.linspace(0, 0.125, 121)
        X, Y = np.meshgrid(xs, ys)
        G = s.alb2.gibbs(X, Y, T)
        # tangent plane at the homogeneous state of the delta = 5 % composition (x = 0.5, y = 0.025)
        g0 = s.alb2.gibbs(np.array([0.5]), np.array([0.025]), T)[0]
        gx = (s.alb2.gibbs(np.array([0.501]), np.array([0.025]), T)[0] - s.alb2.gibbs(np.array([0.499]), np.array([0.025]), T)[0]) / 0.002
        gy = (s.alb2.gibbs(np.array([0.5]), np.array([0.026]), T)[0] - s.alb2.gibbs(np.array([0.5]), np.array([0.024]), T)[0]) / 0.002
        D = 1000 * (G - (g0 + gx * (X - 0.5) + gy * (Y - 0.025))) / 3
        lv = np.linspace(-np.nanmax(np.abs(D)), np.nanmax(np.abs(D)), 31)
        cs = ax.contourf(X, 2 * (1 - Y), D, levels=lv, cmap="RdBu_r")
        ax.contour(X, 2 * (1 - Y), D, levels=[0], colors="k", linewidths=0.8)
        ph = s.equilibrium(0.5, 1.9, T)
        al = ph[ph.phase == "AlB2"]
        ax.plot(al.x_M, 2 * (1 - al.y_Va), "w*", ms=11, mec="k", label="equilibrium phases at δ = 5 %")
        ax.plot([0.5], [1.95], "ko", mfc="none", label="homogeneous reference (x = 0.5, B/M 1.95)")
        fig.colorbar(cs, ax=ax, label="G − tangent plane (meV/atom)")
        ax.set_xlabel("x (Group-V fraction)")
        ax.set_ylabel("B/M of the AlB$_2$ phase")
        ax.set_title(f"{SYSLAB[sid]}, {LABEL[model]} ({scen}), {T:.0f} K")
        ax.legend(fontsize=6, loc="lower left")
    fig.tight_layout()
    save(fig, "fig07_gibbs_surface", "Gibbs energy of the (M,M')(B,Va)2 AlB2 phase at 2000 K, relative to the plane tangent at the homogeneous state x = 0.5, "
         "B/M = 1.95. Blue regions lie below the tangent plane: the homogeneous boron-deficient solution can lower its energy by splitting into a "
         "Zr-rich near-stoichiometric and a Group-V-rich boron-deficient AlB2 phase (white stars: equilibrium phases from the full Gibbs minimisation "
         "with competing borides at overall B/M = 1.90). The vacancy fraction is restricted to y ≤ 0.125 (B/M ≥ 1.75).", ["scripts/x01_six_system_thermo.py (computed on the fly)"])


def fig08(plt):
    c = read("x01_closed_scan.csv")
    for scen, tag in (("own", "own energetics"), ("DFT-anch", "DFT-anchored TaB$_2$ V$_B$")):
        sub = c[(c.scenario == scen) & (c.x == 0.5)]
        systems = [s for s in SYSTEMS if s in set(sub.system)]
        if not systems:
            continue
        fig, axes = plt.subplots(len(systems), 5, figsize=(13, 2.2 * len(systems)), sharex=True, sharey=True, squeeze=False)
        for i, sid in enumerate(systems):
            for j, m in enumerate(MODELS):
                ax = axes[i, j]
                g = sub[(sub.system == sid) & (sub.model == m)]
                if len(g):
                    ts, ds = np.sort(g.T_K.unique()), np.sort(g.delta.unique())
                    te = np.r_[ts[0] - (ts[1] - ts[0]) / 2 if len(ts) > 1 else ts[0] - 125, (ts[1:] + ts[:-1]) / 2, ts[-1] + ((ts[-1] - ts[-2]) / 2 if len(ts) > 1 else 125)]
                    de = np.r_[0, (ds[1:] + ds[:-1]) / 2, ds[-1] + (ds[-1] - ds[-2]) / 2]
                    keys = list(CLS_COL)
                    cl = g.pivot_table(index="delta", columns="T_K", values="class", aggfunc="first").reindex(index=ds, columns=ts)
                    idx = cl.apply(lambda col: col.map(lambda k: keys.index(k) if k in keys else np.nan)).values.astype(float)
                    from matplotlib.colors import BoundaryNorm, ListedColormap
                    ax.pcolormesh(te, 100 * de, idx, cmap=ListedColormap([CLS_COL[k] for k in keys]),
                                  norm=BoundaryNorm(np.arange(len(keys) + 1) - 0.5, len(keys)), shading="flat")
                    for lev, ls in ((0.05, "--"), (0.30, "-")):
                        fr = g.pivot_table(index="delta", columns="T_K", values="second_metal_frac").reindex(index=ds, columns=ts)
                        thr = [100 * fr.index[fr[t] >= lev][0] if (fr[t] >= lev).any() else np.nan for t in ts]
                        ax.plot(ts, thr, ls, color="k", lw=1.1, marker="o", ms=2.5)
                    ax.set_xlim(te[0], te[-1])
                    ax.set_ylim(0, 100 * de[-1])
                if i == 0:
                    ax.set_title(LABEL[m])
                if j == 0:
                    ax.set_ylabel(f"{SYSLAB[sid]}\nδ (% B lost)")
                if i == len(systems) - 1:
                    ax.set_xlabel("T (K)")
        from matplotlib.lines import Line2D

        fig.legend(handles=[Line2D([], [], marker="s", ls="", color=v, label=k) for k, v in CLS_COL.items()] +
                   [Line2D([], [], ls="--", color="k", label="5 % second phase"), Line2D([], [], ls="-", color="k", label="30 % second phase")],
                   loc="lower center", ncol=4, frameon=False)
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        save(fig, f"fig08_phase_maps_{scen.replace('-', '')}", f"Equilibrium assemblage at x = 0.5 versus temperature and boron loss δ (closed system, full Zr/M/B mass balance; {tag}; "
             "AlB2 vacancy fraction ≤ 0.125). Colours: phase assemblage (each cell = one computed (T, δ) point); dashed and solid lines: smallest δ at which 5 % and 30 % of the metal is in the second phase (missing points: not reached by δ = 20 %).",
             ["supplementary_data/x01_closed_scan.csv"])


def group_v_fraction(g: pd.DataFrame) -> np.ndarray:
    """Metal fraction in the Group-V-rich product phase(s): a metal-rich boride, or the higher-x AlB2 phase whichever
    of the two AlB2 phases is larger; 0 for a single phase. Unlike the 'minority' fraction this does not jump when the
    Group-V-rich phase becomes the majority."""
    two = (g["class"] == "two AlB2") & (g.main_x > g.second_x)
    out = np.where(two, 1 - g.second_metal_frac, g.second_metal_frac)
    return np.asarray(out, float)


def fig09_10(plt):
    c = read("x01_closed_scan.csv")
    g0 = c[(c.x == 0.5) & (c.T_K == 2000.0)]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), sharex=True, sharey=True)
    for ax, sid in zip(axes.ravel(), SYSTEMS):
        for m in MODELS:
            for scen, ls in (("own", "-"), ("DFT-anch", "--")):
                g = g0[(g0.system == sid) & (g0.model == m) & (g0.scenario == scen)].sort_values("delta")
                if len(g):
                    ax.plot(100 * g.delta, group_v_fraction(g), ls, color=COL[m], lw=1.2, marker=".", ms=3,
                            label=f"{LABEL[m]}{' (DFT-anch)' if scen != 'own' else ''}")
        ax.set_title(SYSLAB[sid])
        ax.set_xlabel("δ (% of B lost)")
        ax.set_ylabel("metal fraction in Group-V-rich\nproduct phase(s)")
    axes[0, 1].legend(fontsize=5.5, ncol=2)
    fig.tight_layout()
    save(fig, "fig09_phase_fractions", "Fraction of the metal in the Group-V-rich product phase(s) versus boron loss at x = 0.5 and 2000 K for every system "
         "and model (solid: own energetics; dashed: TaB2 vacancy energy anchored to DFT). Product = the higher-x, B-deficient AlB2 phase or a "
         "metal-rich boride. Discontinuities are real changes of assemblage: e.g. at δ = 12.5 % the whole alloy reaches the imposed vacancy cap "
         "(B/M = 1.75) as one phase, and beyond it a metal-rich boride takes up further boron loss. Non-zero values at δ = 0 (DFT-anchored "
         "MACE-OMAT/MatterSim/SevenNet) reflect stoichiometric TaB2 itself being unstable toward TaB2−x + B with the DFT vacancy energy.", ["supplementary_data/x01_closed_scan.csv"])
    g1 = g0[(g0.delta == 0.05) & (g0.scenario.isin(["own", "DFT-anch"]))]
    fig, ax = plt.subplots(figsize=(9, 4))
    for k, sid in enumerate(SYSTEMS):
        for i, m in enumerate(MODELS):
            for scen, mk in (("own", "o"), ("DFT-anch", "s")):
                r = g1[(g1.system == sid) & (g1.model == m) & (g1.scenario == scen)]
                if r.empty:
                    continue
                r = r.iloc[0]
                xpos = k + (i - 2) * 0.1 + (0.04 if scen != "own" else 0)
                ax.plot([xpos, xpos], [r.main_x, r.second_x], "-", color=COL[m], lw=1)
                ax.plot(xpos, r.main_x, mk, color=COL[m], ms=4)
                fill = "full" if r["class"] == "two AlB2" else "none"
                ax.plot(xpos, r.second_x, mk, color=COL[m], ms=6, fillstyle=fill)
                if r["class"] == "two AlB2" and r.second_B_per_M == r.second_B_per_M:
                    ax.text(xpos, r.second_x + 0.03, f"{r.second_B_per_M:.2f}", fontsize=5, ha="center", rotation=90)
    ax.set_xticks(range(6), [SYSLAB[s] for s in SYSTEMS], fontsize=8)
    ax.set_ylabel("Group-V fraction x of each phase")
    ax.set_title("Phase compositions at δ = 5 %, x = 0.5, 2000 K (filled = second AlB$_2$ phase, numbers = its B/M; open = metal-rich boride; ■ DFT-anchored)")
    fig.tight_layout()
    save(fig, "fig10_phase_compositions", "Compositions of the coexisting phases at δ = 5 %, x = 0.5, 2000 K. Lower marker: majority (Zr/Ti/Hf-rich) AlB2 phase; upper marker: "
         "second phase (filled: B-deficient AlB2 solid solution, with its B/M; open: metal-rich boride, x = 1).", ["supplementary_data/x01_closed_scan.csv"])


def fig11(plt):
    f = read("x02_flip_thresholds.csv")
    p = read("x01_params.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    ax = axes[0]
    if f is not None:
        tab = f.copy()
        tab["row"] = tab.system + " " + tab.scenario
        piv = tab.pivot_table(index="row", columns="model", values="d_m_flip").reindex(columns=MODELS)
        im = ax.imshow(piv.values, cmap="coolwarm", vmin=-1.5, vmax=1.5, aspect="auto")
        for (i, j), v in np.ndenumerate(piv.values):
            ax.text(j, i, "none" if np.isnan(v) else f"{v:+.2f}", ha="center", va="center", fontsize=6)
        ax.set_xticks(range(5), [LABEL[m] for m in MODELS], rotation=30)
        ax.set_yticks(range(len(piv)), piv.index, fontsize=7)
        fig.colorbar(im, ax=ax, label="Group-V V$_B$ shift that flips the assemblage (eV)")
        ax.set_title("(a) flip threshold at x = 0.5, 2000 K, B/M = 1.95")
    ax = axes[1]
    if p is not None:
        q = p[p.status == "ok"]
        cats = [("E$_v$(NbB$_2$)\n(eV)", "Nb", "dEv_MB2_eV"), ("E$_v$(TaB$_2$)\n(eV)", "Ta", "dEv_MB2_eV"),
                ("W(NbB$_2$)\n(eV/fu)", "Nb", "W_MB2_eV_fu"), ("W(TaB$_2$)\n(eV/fu)", "Ta", "W_MB2_eV_fu"),
                ("coupling slope ×10\n(eV per Group-V nn,\nmean of systems)", None, "slope_eV_per_nn")]
        for k, (lab, m_el, col) in enumerate(cats):
            for i, m in enumerate(MODELS):
                qm = q[q.model == m]
                if m_el:
                    qm = qm[qm.system.str.endswith(m_el)]
                    val = qm[col].iloc[0] if len(qm) else np.nan
                else:
                    val = 10 * qm[col].mean()
                ax.plot(k + (i - 2) * 0.08, val, "o", color=COL[m], ms=6, label=LABEL[m] if k == 0 else None)
        ax.plot(1, -0.389, "k*", ms=12, label="DFT TaB$_2$ (Ektarawong 2023)")
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_xticks(range(len(cats)), [c[0] for c in cats], fontsize=7)
        ax.set_ylabel("value")
        ax.legend(fontsize=6)
        ax.set_title("(b) controlling terms, per model")
    fig.tight_layout()
    save(fig, "fig11_sensitivity", "Model disagreement and sensitivity. (a) Shift of the Group-V boron-vacancy energy needed to change the predicted phase "
         "assemblage at the reference condition, per system, model and scenario ('none': no change within ±1.5 eV). (b) The controlling terms for each MLIP: "
         "Group-V vacancy energies, vacancy–vacancy interaction W and the environment coupling slope (×10).", ["supplementary_data/x02_flip_thresholds.csv", "supplementary_data/x01_params.csv"])


def fig12(plt):
    p01 = pd.read_csv(VAL / "validation_workflow" / "phase_fraction_analysis" / "results" / "p01_ness_lamellae.csv")
    c = read("x01_closed_scan.csv")
    name = {"TiB2-NbB2": "TiNb", "TiB2-TaB2": "TiTa", "ZrB2-NbB2": "ZrNb", "ZrB2-TaB2": "ZrTa", "HfB2-NbB2": "HfNb", "HfB2-TaB2": "HfTa"}
    fig, ax = plt.subplots(figsize=(9, 4.2))
    for k, sid in enumerate(SYSTEMS):
        e = p01[p01.system.map(name) == sid]
        for r in e.itertuples():
            lo, hi = min(r.A_light_minus10, r.A_light_plus10), max(r.A_light_minus10, r.A_light_plus10)
            v, lo, hi = [min(t, 1 - t) for t in (r.A_light_otsu, lo, hi)]
            ax.errorbar(k - 0.25, v, yerr=[[max(v - min(lo, hi), 0)], [max(max(lo, hi) - v, 0)]], fmt="ks", ms=5, capsize=2)
        for j, (d, off) in enumerate(((0.02, -0.05), (0.04, 0.1), (0.06, 0.25))):
            for scen, mk in (("own", "o"), ("DFT-anch", "s")):
                g = c[(c.system == sid) & (c.x == 0.5) & (c.T_K == 2000.0) & (c.scenario == scen) & np.isclose(c.delta, d)]
                if len(g):
                    ax.vlines(k + off, g.second_metal_frac.min(), g.second_metal_frac.max(), color="C%d" % (j + 1), lw=3 if scen == "own" else 1.5,
                              alpha=0.6 if scen == "own" else 1)
    from matplotlib.lines import Line2D

    ax.legend(handles=[Line2D([], [], marker="s", ls="", color="k", label="Ness ROI minority area fraction (exp.)"),
                       Line2D([], [], color="C1", lw=3, label="predicted, δ = 2 % (range of 5 MLIPs)"),
                       Line2D([], [], color="C2", lw=3, label="predicted, δ = 4 %"), Line2D([], [], color="C3", lw=3, label="predicted, δ = 6 %"),
                       Line2D([], [], color="k", lw=1.5, label="thin: DFT-anchored")], fontsize=6)
    ax.set_xticks(range(6), [SYSLAB[s] for s in SYSTEMS])
    ax.set_ylabel("minority-phase fraction")
    ax.set_title("Experimental lamellar-region fractions vs predicted equilibrium second-phase fractions (x = 0.5, 2000 K)")
    fig.tight_layout()
    save(fig, "fig12_exp_vs_pred_fractions", "Experimental minority-lamella area fractions inside lamella-only regions of Ness's Figure 5 (black; ±10 % threshold "
         "sensitivity) compared with the predicted equilibrium minority-phase fraction after 2, 4 and 6 % boron loss (bars: range over the five MLIPs; "
         "thick: own energetics; thin: DFT-anchored TaB2). Region fractions are not whole-sample phase fractions; the comparison tests plausibility only.",
         ["validation_workflow/phase_fraction_analysis/results/p01_ness_lamellae.csv", "supplementary_data/x01_closed_scan.csv"])


def fig13(plt):
    pk = pd.read_csv(VAL / "validation_workflow" / "phase_fraction_analysis" / "results" / "p05_alb2_peak_positions.csv")
    br = pd.read_csv(VAL / "validation_workflow" / "phase_fraction_analysis" / "results" / "p05_boride_reflections.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, sid, m_bor in zip(axes, ("ZrTa", "ZrNb"), (["TaB", "Ta3B4", "Ta5B6"], ["Nb3B4", "Nb5B6", "NbB"])):
        rows = []
        h = pk[(pk.system == sid) & (pk.model == "orb") & (pk.wavelength == "CuKa1")]
        hom = h[h.case.str.startswith("homogeneous")]
        rows.append(("homogeneous (Zr$_{0.5}$M$_{0.5}$)B$_2$", hom.two_theta_phase1.values, "k"))
        two = pk[(pk.system == sid) & (pk.wavelength == "CuKa1") & pk.case.str.startswith("two")]
        if len(two):
            t = two[two.model == two.model.iloc[0]]
            rows.append((f"two AlB$_2$ ({LABEL[t.model.iloc[0]]}, {t.scenario.iloc[0].split(' ')[0]}): Zr-rich", t.two_theta_phase1.values, "C3"))
            rows.append(("two AlB$_2$: Group-V-rich, B-deficient", t.two_theta_phase2.values, "C3"))
        for b in m_bor:
            g = br[(br.model == "orb") & (br.phase.str.endswith(b))]
            rows.append((f"{b} (metal-rich boride)", g.two_theta_CuKa1.values, "C0"))
        for i, (lab, pos, c) in enumerate(rows):
            ax.vlines(pos, i, i + 0.8, color=c, lw=1.4)
            ax.text(80.5, i + 0.3, lab, fontsize=6.5, va="center")
        ax.set_xlim(20, 80)
        ax.set_yticks([])
        ax.axvspan(33.8, 35.0, color="C0", alpha=0.08)
        ax.set_xlabel("2θ (deg, Cu Kα$_1$ 1.5406 Å)")
        ax.set_title(f"{SYSLAB[sid]}: predicted reflection positions (prediction, no experimental pattern)")
    fig.tight_layout()
    save(fig, "fig13_diffraction", "Predicted powder-diffraction reflection positions (Cu Kα1) for the competing assemblages, from ORB-relaxed structures and "
         "the lattice model of the B-deficient phase. Two AlB2 phases split the (00l) family by 1–3.6°; metal-rich borides add strong reflections at "
         "2θ ≈ 33.8–35.0° (shaded). No experimental diffraction of the Ness melts exists; these are predictions.",
         ["validation_workflow/phase_fraction_analysis/results/p05_alb2_peak_positions.csv", "validation_workflow/phase_fraction_analysis/results/p05_boride_reflections.csv"])


def fig14(plt):
    v = read("SIX_SYSTEM_VALIDATION_numeric.csv")
    if v is None:
        raise FileNotFoundError("SIX_SYSTEM_VALIDATION_numeric.csv")
    cols = [c for c in v.columns if c != "system"]
    fig, ax = plt.subplots(figsize=(12, 3.6))
    M = v[cols].to_numpy(float)
    im = ax.imshow(M, cmap="RdYlGn_r", vmin=0, vmax=1, aspect="auto")
    for (i, j), val in np.ndenumerate(M):
        ax.text(j, i, "n/a" if np.isnan(val) else f"{val:.2f}", ha="center", va="center", fontsize=7)
    ax.set_xticks(range(len(cols)), cols, rotation=25, ha="right", fontsize=7)
    ax.set_yticks(range(len(v)), [SYSLAB[s] for s in v.system])
    fig.colorbar(im, ax=ax, label="normalised indicator (see table)")
    ax.set_title("Six-system validation summary")
    fig.tight_layout()
    save(fig, "fig14_six_system_summary", "Six-system validation summary (normalised indicators; definitions in SIX_SYSTEM_VALIDATION.csv).",
         ["SIX_SYSTEM_VALIDATION.csv"])


def fig15(plt):
    hc = read("blade_MB2_hull_check_vs_MP.csv")
    s = read("blade_before_after_summary.csv")
    fr = [read(f) for f in ("blade_pareto_front_change_2to4.csv", "blade_pareto_front_change_2to5.csv")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    ax = axes[0]
    hc = hc.dropna(subset=["MP_AlB2_ehull_eV_atom"])
    ax.plot(hc.MP_AlB2_ehull_eV_atom * 1000, hc.ORB_dG_hull_eV_atom * 1000, "ko")
    for r in hc.itertuples():
        if max(r.MP_AlB2_ehull_eV_atom, r.ORB_dG_hull_eV_atom) * 1000 > 10:  # on-hull MB2 are summarised in one label
            ax.annotate(f"{r.metal}B$_2$", (r.MP_AlB2_ehull_eV_atom * 1000, r.ORB_dG_hull_eV_atom * 1000), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.annotate("Hf, Nb, Ta, Ti, Y, Zr B$_2$ (on hull)", (0, 0), fontsize=6.5, xytext=(8, 14), textcoords="offset points")
    lim = max(hc.MP_AlB2_ehull_eV_atom.max(), hc.ORB_dG_hull_eV_atom.max()) * 1000 * 1.1
    ax.plot([0, lim], [0, lim], "k:", lw=0.8)
    ax.set_xlabel("MP (GGA) E$_{hull}$ of AlB$_2$-type MB$_2$ (meV/atom)")
    ax.set_ylabel("BLADE/ORB after fix (meV/atom)")
    ax.set_title("(a) diboride hull distance with binary borides")
    ax = axes[1]
    s = s[s.status == "compared"]
    ax.bar(np.arange(len(s)) - 0.2, s.strict_rows_before, 0.4, color="0.6", label="original")
    ax.bar(np.arange(len(s)) + 0.2, s.strict_rows_after, 0.4, color="k", label="with binary borides")
    ax.set_xticks(range(len(s)), s.set, rotation=20, fontsize=7)
    ax.set_ylabel("strict single-phase (composition, T) rows")
    ax.set_title("(b) strict single-phase stability")
    ax.legend(fontsize=7)
    ax = axes[2]
    fr = pd.concat([f for f in fr if f is not None]) if any(f is not None for f in fr) else None
    if fr is not None and len(fr):
        g = fr.groupby("n_metals")[["front_before", "front_after"]].max()
        ax.bar(g.index - 0.2, g.front_before, 0.4, color="0.6", label="original front")
        ax.bar(g.index + 0.2, g.front_after, 0.4, color="k", label="corrected front")
        ax.set_xlabel("number of metals")
        ax.set_ylabel("primary Pareto front size")
        ax.legend(fontsize=7)
    ax.set_title("(c) Pareto-front membership")
    fig.tight_layout()
    save(fig, "fig15_blade_before_after", "Effect of the binary-boride competitor correction (commit 774db35) on BLADE. (a) Hull distance of the AlB2-type "
         "diborides after adding the stable Materials Project binary borides (ORB, BLADE settings) versus Materials Project GGA. (b) Number of strictly "
         "single-phase (composition, temperature) points before and after. (c) Primary Pareto-front size by number of metals.",
         ["supplementary_data/blade_MB2_hull_check_vs_MP.csv", "supplementary_data/blade_before_after_summary.csv", "supplementary_data/blade_pareto_front_change_*.csv"])


def fig16(plt):
    k = read("x01_critical.csv")
    g = k[(k.x == 0.5) & (k.T_K == 2000.0)]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), sharey=False)
    for ax, col, lab in zip(axes, ("delta_onset", "delta_5", "delta_30"), ("δ$_c$ (onset of instability)", "δ$_5$ (≥5 % second phase)", "δ$_{30}$ (≥30 % second phase)")):
        for k_, sid in enumerate(SYSTEMS):
            for scen, off, mk in (("own", -0.12, "o"), ("DFT-anch", 0.12, "s")):
                q = g[(g.system == sid) & (g.scenario == scen)]
                for i_m, m in enumerate(MODELS):
                    r = q[q.model == m]
                    if len(r):
                        ax.plot(k_ + off + (i_m - 2) * 0.035, 100 * r[col].iloc[0] if r[col].iloc[0] == r[col].iloc[0] else 21, mk, color=COL[m], ms=4,
                                mfc=COL[m] if scen == "own" else "none", label=(LABEL[m] if (k_ == 0 and scen == "own") else None))
        ax.axhline(20.5, color="0.7", ls=":")
        ax.text(5.3, 20.8, "not reached ≤ 20 %", fontsize=6, ha="right")
        ax.set_xticks(range(6), [SYSLAB[s] for s in SYSTEMS], rotation=30, fontsize=7)
        ax.set_ylabel("boron loss (% of B)")
        ax.set_title(lab)
    from matplotlib.lines import Line2D
    axes[0].legend(handles=axes[0].get_legend_handles_labels()[0] + [Line2D([], [], marker="s", ls="", mfc="none", color="k", label="DFT-anchored TaB$_2$")],
                   fontsize=6, loc="upper left")
    fig.tight_layout()
    save(fig, "fig16_vacancy_aware_metrics", "Proposed vacancy-aware stability metrics for the six systems at x = 0.5, 2000 K: critical boron loss for destabilisation "
         "(δc), for a detectable (5 %) and for a substantial (30 %) second phase. Circles: own MLIP energetics; open squares: DFT-anchored TaB2. Points at "
         "the top: threshold not reached within 20 % boron loss.", ["supplementary_data/x01_critical.csv"])


def main() -> None:
    plt = setup()
    done = []
    for f in (fig01, fig02, fig03, fig04, fig05, fig06, fig07, fig08, fig09_10, fig11, fig12, fig13, fig14, fig15, fig16):
        try:
            f(plt)
            done.append(f.__name__)
        except Exception as exc:
            print(f"{f.__name__} skipped: {exc!r}")
            traceback.print_exc(limit=1)
        plt.close("all")
    with open(FIG / "CAPTIONS.md", "w") as fh:
        fh.write("# Figure captions\n\n")
        for k, v in CAPTIONS.items():
            fh.write(f"## {k}\n\n{v}\n\n")
    print("figures done:", done)


if __name__ == "__main__":
    main()
