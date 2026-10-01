"""Combined three-panel figures: assemblage map, onset resistance, miscibility gap."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .utils import phase_formula, system_key, system_tag

if TYPE_CHECKING:
    from .config import Config


def compile_combined_figures(cfg: Config) -> None:
    """Build per-system three-panel combined figures and an overview grid.

    Args:
        cfg: Framework configuration object.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri_mod
    import pandas as pd
    from matplotlib import cm

    try:
        from pycalphad import Database, calculate
        from scipy.interpolate import LinearNDInterpolator
        from scipy.spatial import ConvexHull
    except ImportError as e:
        print(f"  compile_combined_figures: skipped ({e})")
        return

    comps_dir = cfg.comps_dir
    sqrt3_2 = np.sqrt(3) / 2
    t_start, t_end, t_step = cfg.miscibility_t_start, cfg.miscibility_t_end, cfg.miscibility_t_step
    mf = cfg.miscibility_metal_fraction
    n_grid = cfg.miscibility_n_grid
    gap_thr = cfg.miscibility_gap_threshold
    fixed_elements = {cfg.phase_element} if cfg.phase_element else set()
    pressure = 101325.0
    norm_t = plt.Normalize(vmin=t_start, vmax=t_end)
    cmap_t = cm.plasma

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
        ax.text(0.5, sqrt3_2 + off, labels[0], ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.text(-off, -off * 0.7, labels[1], ha="right", va="top", fontsize=10, fontweight="bold")
        ax.text(1 + off, -off * 0.7, labels[2], ha="left", va="top", fontsize=10, fontweight="bold")
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

    def _draw_gaps(ax, triang, gap_masks, temps):
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
        print("  compile_combined_figures: no ternary TDB files")
        return
    print(f"  [combined figures] {len(ternary_tdb)} ternary systems")

    grid_norm = _ternary_grid_norm(n_grid)
    x1n, x2n, x3n = grid_norm[:, 0], grid_norm[:, 1], grid_norm[:, 2]
    xc, yc = _cart(x1n, x2n, x3n)
    x_fixed = 1.0 - mf
    triang = tri_mod.Triangulation(xc, yc)
    temps = np.arange(t_start, t_end + t_step, t_step)

    for _, tdb, metals in ternary_tdb:
        all_comps = list(tdb.elements) + ["VA"]
        phases = list(tdb.phases.keys())
        tag = system_tag(metals, cfg.phase_element)
        sys_name = system_key(metals, cfg.phase_element)
        formula = phase_formula(metals, cfg.phase_element, cfg.phase_element_stoichiometry)
        label = "-".join(m.title() for m in metals)
        gap_masks = []
        for T in temps:
            try:
                gm = _gmix(tdb, all_comps, metals, phases, T, grid_norm, x_fixed)
                gap_masks.append(None if np.isnan(gm).all() else _two_phase(xc, yc, gm))
            except Exception:
                gap_masks.append(None)
        any_gap = any(m is not None and m.any() for m in gap_masks)

        fig_c, axes = plt.subplots(1, 3, figsize=(20, 8), gridspec_kw={"wspace": 0.08})

        def _panel_label(ax, txt, fontsize=10):
            ax.text(0.5, -0.08, txt, ha="center", va="top", transform=ax.transAxes, fontsize=fontsize, style="italic")

        for candidate in [
            cfg.figures_dir / sys_name / "muO_x_phase_map" / "assemblage_region_map.png",
            cfg.figures_dir / f"{sys_name}_muO_x_phase_map" / "assemblage_region_map.png",
        ]:
            if candidate.exists():
                axes[0].imshow(plt.imread(candidate))
                axes[0].axis("off")
                _panel_label(axes[0], "μO–x assemblage map")
                break
        else:
            axes[0].axis("off")
            axes[0].text(
                0.5,
                0.5,
                f"{formula}\nassemblage map\nnot yet generated",
                ha="center",
                va="center",
                transform=axes[0].transAxes,
                fontsize=9,
                color="grey",
            )
            _panel_label(axes[0], "μO–x assemblage map")

        onset_csv = cfg.tables_dir / sys_name / f"{tag}_ternary_3d_onset.csv"
        axes[1].set_aspect("equal")
        axes[1].axis("off")
        axes[1].set_xlim(-0.05, 1.05)
        axes[1].set_ylim(-0.12, sqrt3_2 + 0.18)
        _frame(axes[1], [m.title() for m in metals])
        if onset_csv.exists():
            df = pd.read_csv(onset_csv)
            xcols = [c for c in df.columns if c.startswith("x_") and len(c) > 2]
            if len(xcols) >= 3:
                o1 = df[xcols[0]].values
                o2 = df[xcols[1]].values
                o3 = df[xcols[2]].values
                resist = (
                    df["resistance_eV"].values
                    if "resistance_eV" in df.columns
                    else (df["onset_muO_eV"].values - df["onset_muO_eV"].min())
                )
                resist_f = np.where(np.isnan(resist), np.nanmax(resist), resist)
                xco, yco = _cart(o1, o2, o3)
                from matplotlib.patches import PathPatch
                from matplotlib.path import Path as MplPath
                from matplotlib.tri import Triangulation

                tri = Triangulation(xco, yco)
                hb = axes[1].tripcolor(tri, resist_f, shading="flat", cmap="RdYlGn", edgecolors="black", linewidth=0.4)
                hb.set_clip_path(
                    PathPatch(MplPath([(0.0, 0.0), (1.0, 0.0), (0.5, sqrt3_2), (0.0, 0.0)]), transform=axes[1].transData)
                )
                fig_c.colorbar(hb, ax=axes[1], fraction=0.03, pad=0.02, label="resistance (eV)")
        _panel_label(axes[1], "Oxidation resistance (onset μO − baseline)")

        axes[2].set_aspect("equal")
        axes[2].axis("off")
        axes[2].set_xlim(-0.05, 1.05)
        axes[2].set_ylim(-0.12, sqrt3_2 + 0.18)
        _frame(axes[2], [m.title() for m in metals])
        _draw_gaps(axes[2], triang, gap_masks, temps)
        sm2 = cm.ScalarMappable(cmap=cmap_t, norm=norm_t)
        sm2.set_array([])
        fig_c.colorbar(sm2, ax=axes[2], fraction=0.03, pad=0.02, label="T (K)")
        _panel_label(axes[2], f"Miscibility gap ({'found' if any_gap else 'none'})")
        fig_c.suptitle(formula, fontsize=13, fontweight="bold", y=1.01)
        (cfg.figures_dir / f"{sys_name}").mkdir(parents=True, exist_ok=True)
        fig_c.savefig(cfg.figures_dir / f"{sys_name}" / f"{sys_name}_combined.png", dpi=180, bbox_inches="tight")
        plt.close(fig_c)
        print(f"    {label}: saved")

    combined_pngs = []
    for _, _, metals in ternary_tdb:
        sn = system_key(metals, cfg.phase_element)
        combined_pngs.append((sn, cfg.figures_dir / sn / f"{sn}_combined.png"))
    n = len(combined_pngs)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols
    fig_all, axes_all = plt.subplots(
        nrows, ncols, figsize=(ncols * 7.0, nrows * 2.8), gridspec_kw={"wspace": 0.02, "hspace": 0.15}
    )
    axes_all = np.array(axes_all).flatten()
    for idx, (sn, p) in enumerate(combined_pngs):
        ax = axes_all[idx]
        try:
            ax.imshow(plt.imread(str(p.resolve())))
        except Exception:
            ax.text(
                0.5,
                0.5,
                f"{sn}\n(not yet\ngenerated)",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=7,
                color="grey",
            )
        ax.axis("off")
        ax.set_title(sn, fontsize=6, pad=1)
    for ax in axes_all[len(combined_pngs) :]:
        ax.axis("off")
    fig_all.suptitle("All ternary — assemblage map | oxidation resistance | miscibility gap", fontsize=11, y=1.005)
    fig_all.savefig(cfg.figures_dir / "compiled_combined_grid.png", dpi=150, bbox_inches="tight")
    plt.close(fig_all)
    print(f"  compiled_combined_grid.png  ({n} systems)")

    # Standalone miscibility gap grid (small, gap only)
    ncols2 = min(n, 8)
    nrows2 = (n + ncols2 - 1) // ncols2
    fig_g2, axes_g2 = plt.subplots(
        nrows2, ncols2, figsize=(ncols2 * 3.0, nrows2 * 3.0), gridspec_kw={"wspace": 0.05, "hspace": 0.15}
    )
    axes_g2 = np.array(axes_g2).flatten()
    for idx, (_, tdb, metals) in enumerate(ternary_tdb):
        ax = axes_g2[idx]
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, sqrt3_2 + 0.15)
        _frame(ax, [m.title() for m in metals])
        all_comps = list(tdb.elements) + ["VA"]
        phases = list(tdb.phases.keys())
        for T in temps:
            try:
                gm = _gmix(tdb, all_comps, metals, phases, T, grid_norm, x_fixed)
                if np.isnan(gm).all():
                    continue
                mask = _two_phase(xc, yc, gm)
                if mask.any():
                    ax.tricontourf(triang, mask.astype(float), levels=[0.5, 1.5], colors=[cmap_t(norm_t(T))], alpha=0.40)
            except Exception:
                pass
        ax.set_title("-".join(m.title() for m in metals), fontsize=6, pad=1)
    for ax in axes_g2[len(ternary_tdb) :]:
        ax.axis("off")
    sm = cm.ScalarMappable(cmap=cmap_t, norm=norm_t)
    sm.set_array([])
    fig_g2.colorbar(
        sm,
        ax=axes_g2[: len(ternary_tdb)],
        orientation="vertical",
        fraction=0.012,
        pad=0.02,
        aspect=50,
        label="T (K)",
    )
    fig_g2.suptitle("Ternary phase miscibility gaps — all systems", fontsize=10, y=1.005)
    fig_g2.savefig(cfg.figures_dir / "compiled_miscibility_gaps.png", dpi=180, bbox_inches="tight")
    plt.close(fig_g2)
    print(f"  compiled_miscibility_gaps.png  ({n} systems)")
