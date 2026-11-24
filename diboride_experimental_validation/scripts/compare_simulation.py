"""Compare BLADE stability predictions with experiments for six Group IV-V diboride pseudobinaries.

Reads (never writes) the BLADE stability results and SQS structures, and writes:

  simulation_comparison/matched_compositions.csv   experiment -> simulated composition mapping
  simulation_comparison/validation_results.csv     predictions at x = 0.5 for every temperature
  simulation_comparison/model_boundaries.csv       AlB2-only binodal / spinodal vs T (RK model)
  simulation_comparison/model_check.csv            reconstructed model vs stored G_ss
  figures/*.png, *.pdf

Two kinds of prediction are kept apart throughout:

* FULL EQUILIBRIUM (stored BLADE results): convex hull over the AlB2 solid solution,
  its SQS structures and competing phases (elements, other borides) from ORB/ and
  parent_excel/MP_ORB; 0.05 composition grid; 200-2500 K.
* AlB2-ONLY MODEL (reconstructed here): the Redlich-Kister (order 3) + ideal-mixing
  G_ss(x, T) used by BLADE, on a fine grid. Its binodal/spinodal describe phase
  separation *within* the AlB2 structure only. Temperatures above 2500 K are model
  extrapolations; the model contains no liquid, so it cannot describe melting/eutectics.

Usage (from the BLADE repository root):
    .venv/bin/python diboride_experimental_validation/scripts/compare_simulation.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
VAL = HERE.parent
REPO = VAL.parent
sys.path.insert(0, str(REPO / "src"))

from blade.analysis._stability_tables import build_phase_data, detect_metals  # noqa: E402
from blade.oxidation.framework._lp import solve_grand_lp  # noqa: E402

# ── Configuration ─────────────────────────────────────────────────────────────
SIM_DIR = REPO / "Files" / "Stability_blade"          # BLADE stability output (2-3 metals)
STRUCT_DIR = REPO / "Files" / "system_structures"     # SQS + competitor inputs
COMPETITORS = ["ORB", "parent_excel/MP_ORB"]           # as used for SIM_DIR (settings_sig v4)
RK_ORDER = 3
X_EXP = 0.50                                           # nominal Ness et al. composition
MODEL_DX = 0.001                                       # fine grid for the AlB2-only model
MODEL_T = np.arange(200.0, 4000.0 + 1e-9, 25.0)        # model temperatures (K); >2500 K extrapolated
GAP_TOL = 1e-6                                          # eV/f.u. above the AlB2 hull = in gap
KB = 8.617333262145e-5                                  # eV/K

# system_id: (simulation tag, Group IV metal, Group V metal)
SYSTEMS = {
    "TiNb": ("NbTiB", "Ti", "Nb"),
    "TiTa": ("TaTiB", "Ti", "Ta"),
    "ZrNb": ("NbZrB", "Zr", "Nb"),
    "ZrTa": ("TaZrB", "Zr", "Ta"),
    "HfNb": ("HfNbB", "Hf", "Nb"),
    "HfTa": ("HfTaB", "Hf", "Ta"),
}
ORDER = list(SYSTEMS)

OUT_CMP = VAL / "simulation_comparison"
OUT_FIG = VAL / "figures"


# ── Helpers ───────────────────────────────────────────────────────────────────
def classify(row: pd.Series) -> str:
    """Equilibrium assemblage class from the stored LP flags."""
    if row["is_unique_homogeneous"]:
        return "homogeneous_single_phase"
    if row["n_active_fixed"] == 0 and row["n_active_ss"] >= 2:
        return "two_or_more_AlB2_solid_solutions"
    if row["n_active_fixed"] >= 1:
        return "decomposes_incl_non_AlB2_phases"
    return "other"


def lower_hull_gap(x: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Height of g above its lower convex envelope on a sorted 1-D grid."""
    hull = [0]
    for i in range(1, len(x)):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            # keep b only if it lies strictly below the chord a -> i
            if (g[b] - g[a]) * (x[i] - x[a]) <= (g[i] - g[a]) * (x[b] - x[a]):
                break
            hull.pop()
        hull.append(i)
    env = np.interp(x, x[hull], g[hull])
    return g - env


def intervals(mask: np.ndarray, x: np.ndarray) -> list[tuple[float, float]]:
    """Contiguous True intervals of mask as (x_start, x_end)."""
    out, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if (not m or i == len(mask) - 1) and start is not None:
            end = i if m else i - 1
            out.append((float(x[start]), float(x[end])))
            start = None
    return out


def load_sim(tag: str) -> pd.DataFrame:
    return pd.read_parquet(SIM_DIR / f"{tag}_stability.parquet")


def model_gss(tag: str, xv_metal: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fine-grid AlB2-only G_ss model: returns (x_V, H0, mix_shape), sorted by x_V."""
    bdir = STRUCT_DIR / tag / "blade"
    metals = detect_metals(bdir, "B")
    pdat = build_phase_data(bdir, metals, "B", 2.0, rk_order=RK_ORDER, y_step=MODEL_DX)
    xv = pdat["y_grid"][:, metals.index(xv_metal)]
    order = np.argsort(xv)
    return xv[order], pdat["phase_H0"][order], pdat["phase_mix_shape"][order]


def assemblage_at(tag: str, x_v: float, xv_metal: str, temps: np.ndarray) -> list[dict]:
    """Recompute the full-equilibrium LP at x_V to list active phases (same settings as SIM_DIR)."""
    sdir = STRUCT_DIR / tag
    metals = detect_metals(sdir / "blade", "B")
    comp_dirs = [sdir / c for c in COMPETITORS if (sdir / c).is_dir()]
    pdat = build_phase_data(sdir / "blade", metals, "B", 2.0, rk_order=RK_ORDER, y_step=0.05, competitors_dir=comp_dirs)
    y = np.array([x_v if m == xv_metal else 1.0 - x_v for m in metals])
    nf = pdat["n_fixed"]
    rows = []
    for t in temps:
        grand = np.concatenate([pdat["fixed_energy_formula"], pdat["phase_H0"] + KB * t * pdat["phase_mix_shape"]])
        amt, obj, ok = solve_grand_lp(pdat["A_eq"], np.concatenate([y, [2.0]]), grand)
        if not ok:
            rows.append({"T_K": t, "assemblage": "LP failed"})
            continue
        parts = []
        for k in np.flatnonzero(amt > 1e-8):
            if k >= nf:
                xs = pdat["y_grid"][k - nf][metals.index(xv_metal)]
                parts.append(f"AlB2-SS(x={xs:.2f}):{amt[k]:.3f}")
            else:
                parts.append(f"{pdat['phase_ids'][k]}:{amt[k]:.3f}")
        rows.append({"T_K": t, "assemblage": "; ".join(parts), "G_hull_eV_fu": obj})
    return rows


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    warnings.filterwarnings("ignore")
    OUT_CMP.mkdir(exist_ok=True)
    OUT_FIG.mkdir(exist_ok=True)
    exp = pd.read_csv(VAL / "processed_data" / "experimental_observations.csv")
    sol = pd.read_csv(VAL / "processed_data" / "solubility_limits.csv")
    for col in ("x_solute_lower_bound", "x_solute_upper_bound", "temperature_K"):
        sol[col] = pd.to_numeric(sol[col], errors="coerce")  # 'missing' -> NaN

    matched, results, bounds, checks, curves = [], [], [], [], {}
    sims, models = {}, {}
    for sid, (tag, m4, m5) in SYSTEMS.items():
        sim = load_sim(tag)
        sims[sid] = sim
        xcol = f"x_{m5}"
        grid = np.sort(sim[xcol].unique())
        j = int(np.argmin(np.abs(grid - X_EXP)))
        x_sim = float(grid[j])
        matched.append({
            "system_id": sid, "simulation_system": tag, "simulation_file": str((SIM_DIR / f"{tag}_stability.parquet").relative_to(REPO)),
            "composition_variable": f"{xcol} (mole fraction of {m5}B2 on metal sublattice)",
            "x_experiment_nominal": X_EXP, "x_simulation": x_sim, "exact_match": bool(abs(x_sim - X_EXP) < 1e-9),
            "grid_spacing": float(np.round(np.diff(grid).min(), 6)),
            "T_grid_K": f"{sim['T'].min():.0f}-{sim['T'].max():.0f} step {np.diff(np.sort(sim['T'].unique())).min():.0f}",
            "settings_sig": sim["settings_sig"].iloc[0],
        })

        # Stored full-equilibrium predictions at x = 0.5; chemical Hessian from stored G_ss (central FD, h = 0.05).
        temps = np.sort(sim["T"].unique())
        asm = {r["T_K"]: r for r in assemblage_at(tag, x_sim, m5, temps)}
        for t in temps:
            st = sim[sim["T"] == t].set_index(sim[sim["T"] == t][xcol].round(6))
            r = st.loc[round(x_sim, 6)]
            h = 0.05
            try:
                d2 = (st.loc[round(x_sim + h, 6), "G_ss"] - 2 * r["G_ss"] + st.loc[round(x_sim - h, 6), "G_ss"]) / h**2
            except KeyError:
                d2 = np.nan
            results.append({
                "system_id": sid, "T_K": t, "x_groupV": x_sim,
                "dG_hull_eV_per_fu": r["dG_hull"], "ss_fraction": r["ss_fraction"],
                "is_parent_structure_retained": bool(r["is_parent_structure_retained"]),
                "is_single_phase_strict": bool(r["is_single_phase_strict"]),
                "is_unique_homogeneous": bool(r["is_unique_homogeneous"]),
                "n_active_AlB2_phases": int(r["n_active_ss"]), "n_active_other_phases": int(r["n_active_fixed"]),
                "predicted_class": classify(r),
                "hessian_d2Gss_dx2_eV_per_fu": d2,
                "lp_assemblage_recomputed": asm[t]["assemblage"],
                "lp_recompute_matches_stored": bool(abs(asm[t].get("G_hull_eV_fu", np.nan) - r["G_lp"]) < 1e-6),
                "source": "BLADE stability (full equilibrium, stored)",
            })

        # AlB2-only fine model; verify it reproduces the stored G_ss on the 0.05 grid.
        xv, h0, mix = model_gss(tag, m5)
        models[sid] = (xv, h0, mix)
        for t in (300.0, 1500.0, 2500.0):
            st = sim[sim["T"] == t]
            g_model = np.interp(st[xcol], xv, h0 + KB * t * mix)
            checks.append({"system_id": sid, "T_K": t, "max_abs_diff_eV_per_fu": float(np.max(np.abs(g_model - st["G_ss"])))})
        for t in MODEL_T:
            g = h0 + KB * t * mix
            gap = lower_hull_gap(xv, g) > GAP_TOL
            d2 = np.gradient(np.gradient(g, xv), xv)
            spin = d2 < 0
            for kind, mask in (("binodal_AlB2_only", gap), ("spinodal_AlB2_only", spin)):
                for a, b in intervals(mask, xv):
                    bounds.append({"system_id": sid, "T_K": t, "boundary": kind, "x_low": a, "x_high": b,
                                   "extrapolated_beyond_sim_grid": bool(t > 2500), "source": "AlB2-only RK model (reconstructed)"})
        curves[sid] = (xv, h0, mix)

    pd.DataFrame(matched).to_csv(OUT_CMP / "matched_compositions.csv", index=False)
    res = pd.DataFrame(results)
    res.to_csv(OUT_CMP / "validation_results.csv", index=False)
    bnd = pd.DataFrame(bounds)
    bnd.to_csv(OUT_CMP / "model_boundaries.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT_CMP / "model_check.csv", index=False)

    # Per-system summary used by the report.
    summ = []
    for sid in ORDER:
        r = res[res.system_id == sid]
        b = bnd[(bnd.system_id == sid)]
        gap_half = b[(b.boundary == "binodal_AlB2_only") & (b.x_low <= X_EXP) & (b.x_high >= X_EXP)]
        spin_half = b[(b.boundary == "spinodal_AlB2_only") & (b.x_low <= X_EXP) & (b.x_high >= X_EXP)]
        tc_b = b[b.boundary == "binodal_AlB2_only"].T_K.max() if (b.boundary == "binodal_AlB2_only").any() else np.nan
        tc_s = b[b.boundary == "spinodal_AlB2_only"].T_K.max() if (b.boundary == "spinodal_AlB2_only").any() else np.nan
        summ.append({
            "system_id": sid,
            "experiment": "phase separated (lamellae after solidification from laser melt)",
            "classes_at_x0.5_200-2500K": "; ".join(f"{c}:{n}" for c, n in r.predicted_class.value_counts().items()),
            "homogeneous_T_range_K": (f"{r[r.is_unique_homogeneous].T_K.min():.0f}-{r[r.is_unique_homogeneous].T_K.max():.0f}"
                                      if r.is_unique_homogeneous.any() else "none"),
            "class_at_2500K": r[r.T_K == r.T_K.max()].predicted_class.iloc[0],
            "dG_hull_at_2500K_eV_fu": float(r[r.T_K == r.T_K.max()].dG_hull_eV_per_fu.iloc[0]),
            "hessian_negative_T_count": int((r.hessian_d2Gss_dx2_eV_per_fu < 0).sum()),
            "AlB2_gap_covers_x0.5_up_to_K": float(gap_half.T_K.max()) if len(gap_half) else np.nan,
            "AlB2_spinodal_covers_x0.5_up_to_K": float(spin_half.T_K.max()) if len(spin_half) else np.nan,
            "AlB2_gap_Tc_K_(max T with gap; 4000 K = still open)": tc_b,
            "AlB2_spinodal_Tc_K": tc_s,
        })
    pd.DataFrame(summ).to_csv(OUT_CMP / "system_summary.csv", index=False)

    make_figures(res, bnd, sims, curves, exp, sol)
    print("wrote", OUT_CMP, "and", OUT_FIG)


# ── Figures ───────────────────────────────────────────────────────────────────
def make_figures(res, bnd, sims, curves, exp, sol) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "savefig.dpi": 200, "figure.dpi": 100})
    labels = {"TiNb": "TiB$_2$–NbB$_2$", "TiTa": "TiB$_2$–TaB$_2$", "ZrNb": "ZrB$_2$–NbB$_2$",
              "ZrTa": "ZrB$_2$–TaB$_2$", "HfNb": "HfB$_2$–NbB$_2$", "HfTa": "HfB$_2$–TaB$_2$"}
    # End-member melting points (K): Thorpe2026 (S06) where available; NbB2 from Eckert 1966 Fig. 42 label (S07).
    tm = {"Ti": 3313, "Zr": 3404, "Hf": 3529, "Ta": 3284, "Nb": 3301}
    tm_src = {"Ti": "Thorpe", "Zr": "Thorpe", "Hf": "Thorpe", "Ta": "Thorpe", "Nb": "Eckert"}

    def save(fig, name):
        for ext in ("png", "pdf"):
            fig.savefig(OUT_FIG / f"{name}.{ext}", bbox_inches="tight")
        plt.close(fig)

    def otani(ax, sid, ymax):
        """Experimental solubility constraints in ZrB2 from processed_data/solubility_limits.csv.

        Otani 2009 floating-zone brackets (records L01/L02; growth temperature not stated, near
        solidus) are drawn as a vertical band; Vorotilo 2022 (L09, hot pressed at 2173 K) as a
        horizontal bar at its known temperature.
        """
        if not sid.startswith("Zr"):
            return
        m5 = SYSTEMS[sid][2]
        sub = sol[(sol.host_diboride == "ZrB2") & (sol.solute_diboride == f"{m5}B2")]
        fz = sub[sub.record_id.isin(["L01", "L02"])]
        for _, r in fz.iterrows():
            ax.axvspan(r.x_solute_lower_bound, r.x_solute_upper_bound, color="tab:green", alpha=0.35, lw=0)
            ax.annotate(f"Otani FZ: {r.x_solute_lower_bound:.3g}–{r.x_solute_upper_bound:.3g}\n(near solidus, T not stated)",
                        (r.x_solute_upper_bound, ymax), xytext=(3, -45), textcoords="offset points", fontsize=6.5, color="tab:green")
        for _, r in sub[sub.record_id == "L09"].iterrows():
            ax.plot([r.x_solute_lower_bound, r.x_solute_upper_bound], [r.temperature_K] * 2, color="tab:purple", lw=3)
            ax.annotate("Vorotilo: max Ta in SS\n(hot press 2173 K)", (r.x_solute_upper_bound, r.temperature_K),
                        xytext=(4, 2), textcoords="offset points", fontsize=6.5, color="tab:purple")

    # 1. Classification: experiment vs prediction at x = 0.5 across T.
    cls = ["homogeneous_single_phase", "two_or_more_AlB2_solid_solutions", "decomposes_incl_non_AlB2_phases", "other"]
    colors = ["#4daf4a", "#e41a1c", "#984ea3", "#999999"]
    temps = np.sort(res.T_K.unique())
    mat = np.array([[cls.index(res[(res.system_id == s) & (res.T_K == t)].predicted_class.iloc[0]) for t in temps] for s in ORDER])
    fig, (ax, axe) = plt.subplots(1, 2, figsize=(10, 3.2), gridspec_kw={"width_ratios": [12, 1]})
    ax.imshow(mat, aspect="auto", cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-0.5, 4.5), 4), interpolation="nearest")
    ax.set_yticks(range(len(ORDER)), [labels[s] for s in ORDER])
    ax.set_xticks(range(0, len(temps), 2), [f"{t:.0f}" for t in temps[::2]], rotation=45)
    ax.set_xlabel("Temperature (K)")
    ax.set_title("Predicted equilibrium at x = 0.50 (BLADE full equilibrium)")
    axe.imshow(np.full((len(ORDER), 1), 1), aspect="auto", cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-0.5, 4.5), 4))
    axe.set_xticks([0], ["Expt."])
    axe.set_yticks([])
    axe.set_title("Ness 2026", fontsize=8)
    fig.legend(handles=[Patch(color=c, label=l.replace("_", " ")) for c, l in zip(colors[:3], cls[:3])],
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.22), frameon=False)
    fig.text(0.5, -0.27, "Experiment: lamellar phase separation after cooling from a laser-induced melt (temperature of "
             "formation not reported in accessible text).", ha="center", fontsize=7)
    save(fig, "fig1_classification_x050")

    # 2. dG_hull contours (stored, full equilibrium) with experimental markers.
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharey=True)
    for ax, sid in zip(axes.flat, ORDER):
        sim = sims[sid]
        m5 = SYSTEMS[sid][2]
        p = sim.pivot_table(index="T", columns=f"x_{m5}", values="dG_hull")
        # LP round-off gives |ΔG_hull| ~ 1e-15 eV on the hull; clip below 1e-9 eV to exactly 0 so that
        # on-hull points are not left unfilled (they were white in earlier versions of this figure).
        vals = np.where(np.abs(p.values) < 1e-9, 0.0, p.values) * 1000
        cs = ax.contourf(p.columns, p.index, vals, levels=[0, 1e-3, 0.5, 1, 2, 3, 5], cmap="viridis", extend="max")
        strict = sim.pivot_table(index="T", columns=f"x_{m5}", values="is_single_phase_strict").astype(float)
        ax.contour(strict.columns, strict.index, strict.values, levels=[0.5], colors="white", linewidths=1.2)
        ax.axvline(X_EXP, color="red", ls="--", lw=1.2)
        otani(ax, sid, p.index.max())
        ax.set_title(labels[sid])
        ax.set_xlabel(f"x ({m5}B$_2$)")
        fig.colorbar(cs, ax=ax, label="ΔG above hull (meV/f.u.)")
    axes[0, 0].set_ylabel("T (K)")
    axes[1, 0].set_ylabel("T (K)")
    fig.legend(handles=[Line2D([], [], color="white", lw=1.2, label="strict single-phase boundary (model)"),
                        Line2D([], [], color="red", ls="--", label="Ness 2026 nominal x=0.5: phase separated (T not reported)"),
                        Patch(color="tab:green", alpha=0.35, label="Otani 2009 FZ solubility bracket (near solidus)"),
                        Line2D([], [], color="tab:purple", lw=3, label="Vorotilo 2022 max Ta in SS (2173 K)")],
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.04), frameon=False, facecolor="0.8")
    fig.suptitle("BLADE full-equilibrium ΔG above convex hull (incl. competing phases); darkest = 0 (on hull, |ΔG| < 1e-9 eV clipped)", y=1.0, fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, "fig2_dG_hull_contours")

    # 3. Chemical Hessian d2G_ss/dx2 (AlB2-only model) with zero contour = spinodal.
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharey=True)
    for ax, sid in zip(axes.flat, ORDER):
        xv, h0, mix = curves[sid]
        tt = MODEL_T
        g = h0[None, :] + KB * tt[:, None] * mix[None, :]
        d2 = np.gradient(np.gradient(g, xv, axis=1), xv, axis=1)
        sl = slice(5, -5)  # avoid the ln(x) singularity at the endpoints
        lim = np.nanpercentile(np.abs(d2[:, sl]), 95)
        cs = ax.contourf(xv[sl], tt, d2[:, sl], levels=np.linspace(-lim, lim, 21), cmap="RdBu", extend="both")
        ax.contour(xv[sl], tt, d2[:, sl], levels=[0], colors="k", linewidths=1.2)
        ax.axhline(2500, color="0.3", ls="-.", lw=0.8)
        ax.axvline(X_EXP, color="red", ls="--", lw=1.2)
        ax.set_title(labels[sid])
        ax.set_xlabel(f"x ({SYSTEMS[sid][2]}B$_2$)")
        fig.colorbar(cs, ax=ax, label="d²G$_{ss}$/dx² (eV/f.u.)")
    axes[0, 0].set_ylabel("T (K)")
    axes[1, 0].set_ylabel("T (K)")
    fig.suptitle("AlB$_2$-only chemical Hessian (RK model); black = spinodal (d²G/dx²=0); above dash-dot line = extrapolated beyond simulated T grid", y=1.0, fontsize=9)
    fig.tight_layout()
    save(fig, "fig3_hessian_contours")

    # 4/5. Binodal + spinodal (AlB2-only model) with experimental overlays.
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharey=True)
    for ax, sid in zip(axes.flat, ORDER):
        b = bnd[bnd.system_id == sid]
        for kind, col, lab in (("binodal_AlB2_only", "tab:blue", "binodal (AlB$_2$ miscibility gap)"),
                               ("spinodal_AlB2_only", "tab:orange", "spinodal")):
            s = b[b.boundary == kind]
            ax.plot(s.x_low, s.T_K, ".", ms=1.6, color=col, label=lab)
            ax.plot(s.x_high, s.T_K, ".", ms=1.6, color=col)
        m4, m5 = SYSTEMS[sid][1], SYSTEMS[sid][2]
        tms = [tm[m] for m in (m4, m5) if m in tm]
        if tms:
            ax.axhline(min(tms), color="k", lw=0.8)
            low = min((m for m in (m4, m5)), key=lambda m: tm[m])
            ax.text(0.01, min(tms) + 30, f"lower end-member T$_m$ ({low}B$_2$, {tm_src[low]})", fontsize=6)
        ax.axhspan(2500, 4000, color="0.92", zorder=0)
        if not (b.boundary == "binodal_AlB2_only").any() or b[b.boundary == "binodal_AlB2_only"].T_K.max() < 300:
            ax.text(0.5, 1400, "no AlB$_2$ miscibility gap\npredicted (200–4000 K)", ha="center", va="center",
                    fontsize=8, color="0.25", bbox=dict(fc="white", ec="0.7"))
        ax.axvline(X_EXP, color="red", ls="--", lw=1.2)
        otani(ax, sid, 3900)
        ax.set_xlim(0, 1)
        ax.set_ylim(200, 4000)
        ax.set_title(labels[sid])
        ax.set_xlabel(f"x ({m5}B$_2$)")
    axes[0, 0].set_ylabel("T (K)")
    axes[1, 0].set_ylabel("T (K)")
    axes[0, 0].legend(loc="lower right", fontsize=6.5, markerscale=5)
    fig.suptitle("AlB$_2$-only miscibility gap from the RK model (no liquid; shaded = beyond simulated 2500 K)\n"
                 "red dashed: Ness 2026 x=0.5 phase separated after solidification; green band: Otani 2009 FZ bracket; purple: Vorotilo 2022 (2173 K)", y=1.02, fontsize=9)
    fig.tight_layout()
    save(fig, "fig4_binodal_spinodal_with_experiments")

    # 6. G vs x with common tangents (AlB2-only model).
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5))
    for ax, sid in zip(axes.flat, ORDER):
        xv, h0, mix = curves[sid]
        for t, col in ((1000.0, "tab:blue"), (2000.0, "tab:orange"), (2500.0, "tab:red")):
            g = h0 + KB * t * mix
            gmix = (g - (g[0] + (g[-1] - g[0]) * xv)) * 1000  # meV/f.u. relative to end members
            gap = lower_hull_gap(xv, gmix / 1000) > GAP_TOL
            ax.plot(xv, gmix, color=col, lw=1.2, label=f"{t:.0f} K")
            for a, bb in intervals(gap, xv):
                ia, ib = np.searchsorted(xv, a) - 1, np.searchsorted(xv, bb) + 1
                ia, ib = max(ia, 0), min(ib, len(xv) - 1)
                ax.plot([xv[ia], xv[ib]], [gmix[ia], gmix[ib]], color=col, ls="--", lw=1)
                ax.plot([xv[ia], xv[ib]], [gmix[ia], gmix[ib]], "o", color=col, ms=3)
        ax.axvline(X_EXP, color="red", ls=":", lw=1)
        ax.set_title(labels[sid])
        ax.set_xlabel(f"x ({SYSTEMS[sid][2]}B$_2$)")
        ax.set_ylabel("G$_{mix}$ (meV/f.u.)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("AlB$_2$ solid-solution G$_{mix}$(x) (RK model); dashed = common tangent (equilibrium AlB$_2$ + AlB$_2$ tie line)", y=1.0, fontsize=9)
    fig.tight_layout()
    save(fig, "fig5_Gx_common_tangent")


if __name__ == "__main__":
    main()
