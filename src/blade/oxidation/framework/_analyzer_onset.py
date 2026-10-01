"""Onset-AUC diagram helpers for SystemAnalyzer."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from .utils import make_animation, system_key


def _ensure_imports() -> None:
    """Add the framework directory and its python subdirectory to sys.path."""
    root = Path(__file__).parent
    for p in [str(root), str(root / "python")]:
        if p not in sys.path:
            sys.path.insert(0, p)


def _composition_grid_for_onset(config, n_metals: int):
    """Build a simplex composition grid for the onset AUC scan.

    The step size is chosen per the config and may be doubled if the grid
    exceeds ``config.onset_max_comp_points``.

    Args:
        config: Config instance providing ``onset_comp_step_binary``,
            ``onset_comp_step_ternary``, ``onset_comp_step_higher``, and
            ``onset_max_comp_points``.
        n_metals: Number of metal components.

    Returns:
        Tuple of (comp_grid, step) where *comp_grid* is a list of composition
        arrays and *step* is the final step size used.
    """
    _ensure_imports()
    from thermodynamics import simplex_grid_nd

    cfg = config
    step = (
        cfg.onset_comp_step_binary
        if n_metals == 2
        else cfg.onset_comp_step_ternary
        if n_metals == 3
        else cfg.onset_comp_step_higher
    )
    grid = simplex_grid_nd(n_metals, step)
    while len(grid) > cfg.onset_max_comp_points and step < 0.5:
        step *= 2
        grid = simplex_grid_nd(n_metals, step)
    return grid, step


def _run_onset_auc_diagrams(config, sys_cfg, pd_data) -> None:
    """Compute onset μO and AUC metrics for all compositions and temperatures.

    For each temperature in ``config.onset_auc_T_values`` this function either
    reads a valid cache CSV or runs the grand-LP for every composition on the
    simplex grid, computing the oxidation-onset μO and phase-fraction AUC
    integrals.  Saves per-temperature onset diagrams and an animation.

    Args:
        config: Config instance (controls all onset parameters, paths, and
            plot flags).
        sys_cfg: SystemConfig for the current system.
        pd_data: Phase-data dict produced by ``load_phases()``.
    """
    import pandas as pd

    _ensure_imports()
    from thermodynamics import KB

    from ._analyzer_cache import _onset_cache_matches
    from ._analyzer_metrics import _load_reusable_scalar_states, _solve_metrics, _state_key

    cfg = config
    metals = sys_cfg.metals
    n_metals = len(metals)
    comp_grid, comp_step = _composition_grid_for_onset(cfg, n_metals)
    mu_values = np.asarray(cfg.onset_auc_mu_O, dtype=float)
    sys_name = system_key(metals, cfg.phase_element)
    out_root = cfg.figures_dir / sys_name / "onset_auc"
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"--- Onset diagrams ({len(comp_grid)} comps, step={comp_step:g}) ---")
    onset_frames = []
    all_summary = []
    for T in cfg.onset_auc_T_values:
        csv_path = cfg.tables_dir / f"{sys_cfg.tag}_onset_auc_T{int(T)}.csv"
        if _onset_cache_matches(cfg, csv_path, comp_grid, comp_step, mu_values, T, metals):
            df = pd.read_csv(csv_path)
        else:
            if not cfg.run_calculations:
                raise RuntimeError(f"plot-only mode requires a current onset cache: {csv_path}")
            phase_G = pd_data["phase_H0"] + KB * float(T) * pd_data["phase_mix_shape"]
            reusable = _load_reusable_scalar_states(cfg, sys_cfg, T)
            reused_states = 0
            solved_states = 0
            rows = []
            for comp in comp_grid:
                phase_fraction_curve, parent_fraction_curve, oxide_fraction_curve, sampled_mu = [], [], [], []
                onset_mu = np.nan
                first_infeasible_mu = np.nan
                for mu_o in mu_values:
                    cached = reusable.get(_state_key(comp, mu_o))
                    if cached is not None:
                        phase_fraction, parent_fraction, oxide_fraction, absorbed_o = cached
                        ok = not np.isnan(phase_fraction)
                        reused_states += 1
                    else:
                        r = _solve_metrics(cfg, pd_data, comp, float(T), float(mu_o), phase_G)
                        phase_fraction = r["phase_fraction"] if r["ok"] else np.nan
                        parent_fraction = r["parent_phase_fraction"] if r["ok"] else np.nan
                        oxide_fraction = r["oxide_phase_fraction"] if r["ok"] else np.nan
                        absorbed_o = r["absorbed_O"] if r["ok"] else np.nan
                        ok = r["ok"]
                        solved_states += 1
                    if not ok:
                        first_infeasible_mu = float(mu_o)
                        break
                    phase_fraction_curve.append(phase_fraction)
                    parent_fraction_curve.append(parent_fraction)
                    oxide_fraction_curve.append(oxide_fraction)
                    sampled_mu.append(float(mu_o))
                    if np.isnan(onset_mu) and absorbed_o > cfg.onset_threshold:
                        onset_mu = float(mu_o)
                curve = np.asarray(phase_fraction_curve, dtype=float)
                parent_curve = np.asarray(parent_fraction_curve, dtype=float)
                oxide_curve = np.asarray(oxide_fraction_curve, dtype=float)
                sampled_mu_array = np.asarray(sampled_mu, dtype=float)
                if len(curve) == 0:
                    auc = np.nan
                    parent_auc = np.nan
                    oxide_auc = np.nan
                    feasible_mu_max = np.nan
                else:
                    auc = 0.0 if len(curve) == 1 else float(np.trapezoid(curve, sampled_mu_array))
                    parent_auc = 0.0 if len(parent_curve) == 1 else float(np.trapezoid(parent_curve, sampled_mu_array))
                    oxide_auc = 0.0 if len(oxide_curve) == 1 else float(np.trapezoid(oxide_curve, sampled_mu_array))
                    feasible_mu_max = float(sampled_mu_array[-1])
                rows.append(
                    {
                        "T_K": float(T),
                        **{f"x_{m}": float(comp[i]) for i, m in enumerate(metals)},
                        "onset_muO_eV": onset_mu,
                        "phase_fraction_auc_eV": auc,
                        "parent_phase_fraction_auc_eV": parent_auc,
                        "oxide_phase_fraction_auc_eV": oxide_auc,
                        "phase_element": cfg.phase_element or "",
                        "feasible_muO_max_eV": feasible_mu_max,
                        "first_infeasible_muO_eV": first_infeasible_mu,
                        "feasible_muO_points": len(sampled_mu),
                        "comp_step": float(comp_step),
                        "muO_min": float(mu_values[0]),
                        "muO_max": float(mu_values[-1]),
                        "muO_step": float(mu_values[1] - mu_values[0]) if len(mu_values) > 1 else np.nan,
                    }
                )
            df = pd.DataFrame(rows)
            df.to_csv(csv_path, index=False)
            print(f"    T={int(T)} K: reused {reused_states} states, solved {solved_states} new states")
        all_summary.append(df)

        if not cfg.run_plots:
            continue

        t_dir = out_root / f"T{int(T)}"
        t_dir.mkdir(parents=True, exist_ok=True)
        if n_metals == 2:
            onset_png = _plot_onset_auc_binary(df, sys_cfg, int(T), t_dir)
        elif n_metals == 3:
            onset_png = _plot_onset_auc_ternary(df, sys_cfg, int(T), t_dir)
        else:
            onset_png = None
            (t_dir / "README.txt").write_text(
                f"{''.join(metals)}{cfg.phase_element or ''} has {n_metals} metals.\n"
                "n-ary onset written to CSV; no static simplex plot above ternary.\n"
            )
        if onset_png:
            onset_frames.append(onset_png)

    if all_summary:
        import pandas as pd

        pd.concat(all_summary, ignore_index=True).to_csv(
            cfg.tables_dir / f"{sys_cfg.tag}_onset_auc_all_temperatures.csv",
            index=False,
        )
    if cfg.run_animations:
        make_animation(
            onset_frames,
            out_root / "onset_diagram.gif",
            out_root / "onset_diagram.mp4",
            cfg.anim_fps,
            cfg.mp4_crf,
            cfg.mp4_preset,
        )


def _plot_onset_auc_binary(df, sys_cfg, T, out_dir):
    """Plot onset μO vs composition for a binary system and save as PNG.

    Args:
        df: DataFrame with columns ``x_{metals[0]}`` and ``onset_muO_eV``.
        sys_cfg: SystemConfig providing ``metals`` and ``phase_label``.
        T: Temperature in K (displayed in the title).
        out_dir: Directory where the PNG is saved.

    Returns:
        Path to the saved ``onset_diagram.png``.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    m0 = sys_cfg.metals[0]
    x = df[f"x_{m0}"].values
    order = np.argsort(x)
    x = x[order]
    onset = df["onset_muO_eV"].values[order]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, onset, "o-", lw=1.8, ms=4)
    ax.set_xlabel(f"x_{m0}")
    ax.set_ylabel(r"onset $\mu_O$ (eV/O)")
    ax.set_title(f"{sys_cfg.phase_label} onset | T={T} K")
    ax.grid(True, alpha=0.35)
    onset_png = out_dir / "onset_diagram.png"
    fig.savefig(onset_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return onset_png


def _ternary_frame(ax, metals):
    """Draw the triangular frame and corner labels for a ternary diagram.

    Args:
        ax: Matplotlib axes on which to draw.
        metals: Sequence of three metal element symbols (corners: left, right,
            top).
    """
    sqrt3_2 = np.sqrt(3.0) / 2.0
    tri = np.array([[0.5, sqrt3_2], [0.0, 0.0], [1.0, 0.0], [0.5, sqrt3_2]])
    ax.plot(tri[:, 0], tri[:, 1], "k-", lw=1.5)
    off = 0.06
    ax.text(-off, -off * 0.7, metals[0], ha="right", va="top", fontsize=12, fontweight="bold")
    ax.text(0.5, sqrt3_2 + off, metals[2], ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.text(1 + off, -off * 0.7, metals[1], ha="left", va="top", fontsize=12, fontweight="bold")
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, sqrt3_2 + 0.14)
    ax.set_aspect("equal")
    ax.axis("off")


def _plot_ternary_scalar(df, sys_cfg, T, out_dir, value_col, title, color_label, filename, cmap="viridis"):
    """Render a scalar field on a ternary simplex and save as PNG.

    Args:
        df: DataFrame with columns ``x_{metals[1]}``, ``x_{metals[2]}``, and
            *value_col*.
        sys_cfg: SystemConfig providing ``metals`` and ``phase_label``.
        T: Temperature in K (displayed in the title).
        out_dir: Directory where the PNG is saved.
        value_col: Name of the DataFrame column to colour-map.
        title: String appended to the plot title after ``phase_label``.
        color_label: Label for the colourbar.
        filename: Output filename (relative to *out_dir*).
        cmap: Matplotlib colourmap name.

    Returns:
        Path to the saved PNG.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri_mod

    metals = sys_cfg.metals
    x2 = df[f"x_{metals[1]}"].values
    x3 = df[f"x_{metals[2]}"].values
    vals = df[value_col].astype(float).values
    xcart = x2 + 0.5 * x3
    ycart = (np.sqrt(3.0) / 2.0) * x3
    triang = tri_mod.Triangulation(xcart, ycart)
    tri_vals = vals[triang.triangles]
    triang.set_mask(np.any(np.isnan(tri_vals), axis=1))
    fig, ax = plt.subplots(figsize=(8.5, 7.8))
    _ternary_frame(ax, metals)
    valid = ~np.isnan(vals)
    if valid.sum() >= 3:
        v_min, v_max = float(np.nanmin(vals)), float(np.nanmax(vals))
        if v_max > v_min:
            pc = ax.tricontourf(triang, vals, levels=30, cmap=cmap, vmin=v_min, vmax=v_max)
        else:
            pc = ax.tripcolor(triang, vals, shading="flat", cmap=cmap)
    else:
        pc = ax.tripcolor(triang, vals, shading="flat", cmap=cmap)
    fig.colorbar(pc, ax=ax, fraction=0.04, pad=0.03, label=color_label)
    ax.set_title(f"{sys_cfg.phase_label} {title} | T={T} K", fontsize=12)
    out = out_dir / filename
    fig.savefig(out, dpi=190, bbox_inches="tight")
    plt.close(fig)
    return out


def _plot_onset_auc_ternary(df, sys_cfg, T, out_dir):
    """Plot the oxidation-onset μO field on a ternary simplex and save as PNG.

    Args:
        df: DataFrame produced by the onset AUC calculation loop.
        sys_cfg: SystemConfig providing ``metals`` and ``phase_label``.
        T: Temperature in K.
        out_dir: Directory where the PNG is saved.

    Returns:
        Path to the saved ``onset_diagram.png``.
    """
    return _plot_ternary_scalar(
        df,
        sys_cfg,
        T,
        out_dir,
        "onset_muO_eV",
        "oxidation onset",
        r"onset $\mu_O$ (eV/O)",
        "onset_diagram.png",
        cmap="RdYlGn",
    )
