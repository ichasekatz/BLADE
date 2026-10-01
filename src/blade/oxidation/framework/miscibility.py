"""Miscibility gap plots for ternary CALPHAD systems."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .utils import phase_formula, system_key

if TYPE_CHECKING:
    from .config import Config


def plot_miscibility_gaps(cfg: Config) -> None:
    """Compute and plot miscibility gaps for all ternary TDB systems.

    Args:
        cfg: Framework configuration object.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri_mod
    from matplotlib import cm

    try:
        from pycalphad import Database, calculate
        from scipy.interpolate import LinearNDInterpolator
        from scipy.spatial import ConvexHull
    except ImportError as e:
        print(f"  plot_miscibility_gaps: skipped ({e})")
        return

    comps_dir = cfg.comps_dir
    out_dir = cfg.figures_dir / "miscibility_gaps"
    out_dir.mkdir(parents=True, exist_ok=True)
    t_start, t_end, t_step = cfg.miscibility_t_start, cfg.miscibility_t_end, cfg.miscibility_t_step
    mf = cfg.miscibility_metal_fraction
    n_grid = cfg.miscibility_n_grid
    gap_thr = cfg.miscibility_gap_threshold
    fixed_elements = {cfg.phase_element} if cfg.phase_element else set()
    pressure = 101325.0
    sqrt3_2 = np.sqrt(3) / 2

    def _cart(x1, x2, x3):
        return x2 + 0.5 * x3, sqrt3_2 * x3

    def _ternary_grid_norm(n):
        pts = []
        for i in range(n + 1):
            for j in range(n + 1 - i):
                x1, x2 = i / n, j / n
                pts.append((x1, x2, max(0.0, 1.0 - x1 - x2)))
        return np.array(pts)

    def _frame(ax, labels):
        v = np.array([[0.5, sqrt3_2], [0, 0], [1, 0], [0.5, sqrt3_2]])
        ax.plot(v[:, 0], v[:, 1], "k-", lw=1.5)
        off = 0.08
        ax.text(-off, -off * 0.7, labels[0], ha="right", va="top", fontsize=10, fontweight="bold")
        ax.text(0.5, sqrt3_2 + off, labels[2], ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.text(1 + off, -off * 0.7, labels[1], ha="left", va="top", fontsize=10, fontweight="bold")
        ax.set_aspect("equal")
        ax.axis("off")

    def _gmix(tdb, all_comps, metals, phases, T, grid_norm, x_fixed):
        fixed_els = [e for e in tdb.elements if e.upper() in {f.upper() for f in fixed_elements} and e != "VA"]
        x1, x2, x3 = grid_norm[:, 0], grid_norm[:, 1], grid_norm[:, 2]
        pts = {metals[0]: x1 * mf, metals[1]: x2 * mf, metals[2]: x3 * mf}
        for el in fixed_els:
            pts[el] = np.full(len(x1), x_fixed)
        res = calculate(tdb, all_comps, phases[0], T=float(T), P=pressure, points=pts)
        gm = res.GM.values.ravel()[: len(x1)]
        g_ends = []
        for m in metals:
            ep = {mm: np.array([1e-10]) for mm in metals}
            ep[m] = np.array([mf - 2e-10])
            for el in fixed_els:
                ep[el] = np.array([x_fixed])
            r = calculate(tdb, all_comps, phases[0], T=float(T), P=pressure, points=ep)
            g_ends.append(float(r.GM.values.ravel()[0]))
        return np.array([gm[i] - sum(grid_norm[i, j] * g_ends[j] for j in range(3)) for i in range(len(x1))])

    def _two_phase(xc, yc, gmix):
        pts3d = np.column_stack([xc, yc, gmix])
        try:
            hull = ConvexHull(pts3d)
        except Exception:
            return np.zeros(len(gmix), bool)
        lower = set()
        for simplex, eq in zip(hull.simplices, hull.equations, strict=False):
            if eq[2] < 0:
                lower.update(simplex)
        if len(lower) < 3:
            return np.zeros(len(gmix), bool)
        lidx = np.array(list(lower))
        interp = LinearNDInterpolator(pts3d[lidx, :2], pts3d[lidx, 2])
        gc = interp(np.column_stack([xc, yc]))
        nan_m = np.isnan(gc)
        gc[nan_m] = gmix[nan_m]
        return (gmix - gc) > gap_thr

    def _draw_gaps(ax, triang, gap_masks, temps, cmap_t, norm_t):
        any_gap = False
        for T, mask in zip(temps, gap_masks, strict=False):
            if mask is not None and mask.any():
                any_gap = True
                c = cmap_t(norm_t(T))
                ax.tricontourf(triang, mask.astype(float), levels=[0.5, 1.5], colors=[c], alpha=0.40)
                ax.tricontour(triang, mask.astype(float), levels=[0.5], colors=[c], linewidths=1.2)
        return any_gap

    tdb_files = sorted(comps_dir.rglob("*.tdb")) if comps_dir.exists() else []
    ternary_tdb = []
    for tp in tdb_files:
        if " 2" in tp.name:
            continue
        try:
            tdb = Database(str(tp))
        except Exception:
            continue
        metals = sorted(e for e in tdb.elements if e.upper() not in {f.upper() for f in fixed_elements} and e != "VA")
        if len(metals) == 3:
            ternary_tdb.append((tp, tdb, metals))
    if not ternary_tdb:
        print("  plot_miscibility_gaps: no ternary TDB files found")
        return
    print(f"  [miscibility gaps] {len(ternary_tdb)} ternary TDB files")

    grid_norm = _ternary_grid_norm(n_grid)
    x1n, x2n, x3n = grid_norm[:, 0], grid_norm[:, 1], grid_norm[:, 2]
    xc, yc = _cart(x1n, x2n, x3n)
    x_fixed = 1.0 - mf
    triang = tri_mod.Triangulation(xc, yc)
    temps = np.arange(t_start, t_end + t_step, t_step)
    norm_t = plt.Normalize(vmin=t_start, vmax=t_end)
    cmap_t = cm.plasma

    for _, tdb, metals in ternary_tdb:
        all_comps = list(tdb.elements) + ["VA"]
        phases = list(tdb.phases.keys())
        sys_name = system_key(metals, cfg.phase_element)
        formula = phase_formula(metals, cfg.phase_element, cfg.phase_element_stoichiometry)
        gap_masks = []
        for T in temps:
            try:
                gm = _gmix(tdb, all_comps, metals, phases, T, grid_norm, x_fixed)
                gap_masks.append(None if np.isnan(gm).all() else _two_phase(xc, yc, gm))
            except Exception:
                gap_masks.append(None)
        fig_g, ax_g = plt.subplots(figsize=(7, 6.5))
        _frame(ax_g, [m.title() for m in metals])
        any_gap = _draw_gaps(ax_g, triang, gap_masks, temps, cmap_t, norm_t)
        sm = cm.ScalarMappable(cmap=cmap_t, norm=norm_t)
        sm.set_array([])
        fig_g.colorbar(sm, ax=ax_g, pad=0.02, shrink=0.7, label="T (K)")
        note = "shaded = miscibility gap" if any_gap else "single phase"
        ax_g.set_title(f"{formula}  {t_start}-{t_end} K\n{note}", fontsize=11, fontweight="bold")
        fig_g.tight_layout()
        fname = f"{sys_name}_miscibility"
        fig_g.savefig(out_dir / f"{fname}.png", dpi=180, bbox_inches="tight")
        plt.close(fig_g)
