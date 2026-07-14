"""Compiler — post-processing: combined figures, onset plots, miscibility gaps."""

from __future__ import annotations

import numpy as np

from .config import Config
from .utils import phase_formula, system_key, system_tag


class Compiler:
    """Compile post-processing figures after BatchRunner finishes.

    Usage::

        from .compiler import Compiler; from .config import Config
        cfg = Config(...)
        comp = Compiler(cfg)
        comp.compile(tags_info)   # tags_info from BatchRunner.run()
    """

    def __init__(self, config: Config):
        self.config = config

    def compile(self, tags_info: list) -> None:
        """Run all compilation steps."""
        self.compile_onset_lines(tags_info)
        self.compile_ternary_onset_grid()
        self.plot_miscibility_gaps()
        self.compile_combined_figures()
        print("Post-processing done.")

    # ---------------------------------------------------------------- onset lines

    def compile_onset_lines(self, tags_info: list) -> None:
        """Copy assemblage maps + build combined onset line plots."""
        import shutil

        import matplotlib
        import matplotlib.lines as mlines
        import matplotlib.pyplot as plt
        import pandas as pd
        import plotly.graph_objects as go

        matplotlib.use("Agg")

        cfg = self.config
        cfg.figures_dir.mkdir(parents=True, exist_ok=True)

        compilations = {
            "compiled_muO_x_assemblage": ("muO_x_phase_map", "assemblage_region_map.png"),
            "compiled_fixed_phase_assemblage": ("fixed_phase_muO_x_phase_map", "assemblage_region_map.png"),
            "compiled_oxidation_onset": ("fixed_phase_muO_x_phase_map", "oxidation_onset_comparison.png"),
        }
        for folder_name, (src_subdir, src_file) in compilations.items():
            out_dir = cfg.figures_dir / folder_name
            out_dir.mkdir(parents=True, exist_ok=True)
            n_copied = 0
            for tag, m1, m2, sys_name in tags_info:
                src = cfg.figures_dir / sys_name / src_subdir / src_file
                if src.exists():
                    shutil.copy(src, out_dir / f"{sys_name}_{src_file}")
                    n_copied += 1
            print(f"  {folder_name}/  ({n_copied} files)")

        # Combined onset (fixed=dashed, flexible=solid, color=M2)
        all_m2s = sorted({m2 for _, m1, m2, sn in tags_info})
        m2_cmap = plt.get_cmap("tab10", max(len(all_m2s), 1))
        m2_color = {m2: m2_cmap(i) for i, m2 in enumerate(all_m2s)}
        fig, ax = plt.subplots(figsize=(16, 9))
        legend_handles = {}
        for tag, m1, m2, sys_name in tags_info:
            color = m2_color.get(m2, "grey")
            for csv_suf, y_col, ls in [
                ("_fixed_phase_muO_x_phase_map_onset_vs_flexible.csv", "onset_muO_flexible_phase_eV", "-"),
                ("_fixed_phase_muO_x_phase_map_oxidation_onset.csv", "onset_muO_fixed_phase_eV", "--"),
            ]:
                f = cfg.tables_dir / sys_name / f"{tag}{csv_suf}"
                if not f.exists():
                    continue
                df = pd.read_csv(f)
                x_col = f"x_{m1}"
                if x_col not in df.columns or y_col not in df.columns:
                    continue
                mask = df[y_col].notna()
                if not mask.any():
                    continue
                (line,) = ax.plot(
                    df.loc[mask, x_col],
                    df.loc[mask, y_col],
                    lw=1.6,
                    color=color,
                    linestyle=ls,
                    label=sys_name if ls == "-" else "_",
                )
                if ls == "-":
                    legend_handles[sys_name] = line
        solid_patch = mlines.Line2D([], [], color="k", lw=1.6, linestyle="-", label="flexible phase")
        dash_patch = mlines.Line2D([], [], color="k", lw=1.6, linestyle="--", label="fixed phase")
        m2_handles = [mlines.Line2D([], [], color=m2_color[m], lw=3, label=f"M2={m}") for m in all_m2s]
        ncol = max(1, (len(legend_handles) + len(m2_handles) + 2) // 25 + 1)
        ax.legend(
            handles=list(legend_handles.values()) + [solid_patch, dash_patch] + m2_handles,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            fontsize=8,
            ncol=ncol,
            framealpha=0.9,
            borderaxespad=0,
            title="System  |  — flexible  -- fixed  |  M2 color",
        )
        ax.set_xlabel("initial $x_{M1}$ fraction", fontsize=13)
        ax.set_ylabel(r"oxidation onset $\mu_O$ (eV per O atom)", fontsize=13)
        ax.set_title("Oxidation onset — all systems\n(solid=flexible, dashed=fixed, color=M2)", fontsize=13)
        ax.set_xlim(0.0, 1.0)
        ax.grid(True, alpha=0.4)
        fig.tight_layout()
        fig.savefig(cfg.figures_dir / "combined_onset.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("  combined_onset.png")

        # Separate fixed / flexible
        cmap_sys = plt.get_cmap("tab20", max(len(tags_info), 1))
        for mode, y_col, csv_suf, title, fname in [
            (
                "fixed",
                "onset_muO_fixed_phase_eV",
                "_fixed_phase_muO_x_phase_map_oxidation_onset.csv",
                "Oxidation onset — fixed-composition phase (all systems)",
                "combined_onset_fixed",
            ),
            (
                "flexible",
                "onset_muO_flexible_phase_eV",
                "_fixed_phase_muO_x_phase_map_onset_vs_flexible.csv",
                "Oxidation onset — flexible-composition phase (all systems)",
                "combined_onset_flexible",
            ),
        ]:
            fig2, ax2 = plt.subplots(figsize=(15, 8))
            pfig = go.Figure()
            n_plotted = 0
            for i, (tag, m1, m2, sys_name) in enumerate(tags_info):
                f = cfg.tables_dir / sys_name / f"{tag}{csv_suf}"
                if not f.exists():
                    continue
                df = pd.read_csv(f)
                x_col = f"x_{m1}"
                if x_col not in df.columns or y_col not in df.columns:
                    continue
                mask = df[y_col].notna()
                if not mask.any():
                    continue
                xs, ys = df.loc[mask, x_col].values, df.loc[mask, y_col].values
                c_mpl = cmap_sys(i % 20)
                c_hex = "#{:02x}{:02x}{:02x}".format(int(c_mpl[0] * 255), int(c_mpl[1] * 255), int(c_mpl[2] * 255))
                ax2.plot(xs, ys, lw=1.5, color=c_mpl, label=sys_name)
                pfig.add_trace(
                    go.Scatter(
                        x=xs,
                        y=ys,
                        mode="lines",
                        name=sys_name,
                        line=dict(color=c_hex, width=2),
                        hovertemplate=f"<b>{sys_name}</b><br>x_M1({m1})=%{{x:.3f}}<br>onset μO=%{{y:.3f}} eV/O<extra></extra>",
                    )
                )
                n_plotted += 1
            ax2.set_xlabel("initial x_M1 fraction", fontsize=13)
            ax2.set_ylabel("onset μO (eV per O atom)", fontsize=13)
            ax2.set_title(title, fontsize=14)
            ax2.set_xlim(0.0, 1.0)
            ax2.grid(True, alpha=0.4)
            ncol = max(1, n_plotted // 25 + 1)
            ax2.legend(
                loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, ncol=ncol, framealpha=0.9, borderaxespad=0
            )
            fig2.tight_layout()
            fig2.savefig(cfg.figures_dir / f"{fname}.png", dpi=150, bbox_inches="tight")
            plt.close(fig2)
            pfig.update_layout(
                title=title,
                width=1100,
                height=700,
                xaxis=dict(title="initial x_M1 fraction", range=[0, 1]),
                yaxis=dict(title="onset μO (eV per O atom)"),
                hovermode="closest",
                legend=dict(x=1.01, y=1, font=dict(size=9)),
            )
            pfig.write_html(str(cfg.figures_dir / f"{fname}.html"))
            print(f"  {fname}.png + .html  ({n_plotted} systems)")

    # ---------------------------------------------------------------- ternary onset grid

    def compile_ternary_onset_grid(self) -> None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        import pandas as pd
        import plotly.graph_objects as go

        cfg = self.config
        onset_files = sorted(cfg.tables_dir.glob("*/*_ternary_3d_onset.csv"))
        if not onset_files:
            print("  compiled_ternary_onset_grid: no ternary onset CSVs found")
            return

        systems = []
        for f in onset_files:
            df = pd.read_csv(f)
            x_cols = [c for c in df.columns if c.startswith("x_") and len(c) > 2]
            if len(x_cols) < 3:
                continue
            m1, m2, m3 = [c[2:] for c in x_cols[:3]]
            resist = (
                df["resistance_eV"].values
                if "resistance_eV" in df.columns
                else df["onset_muO_eV"].values
                - (df["onset_muO_eV"].min() if df["onset_muO_eV"].notna().any() else -10.0)
            )
            suffix = f"–{cfg.phase_element}" if cfg.phase_element else ""
            systems.append(
                {
                    "label": f"{m1}–{m2}–{m3}{suffix}",
                    "x1": df[x_cols[0]].values,
                    "x2": df[x_cols[1]].values,
                    "x3": df[x_cols[2]].values,
                    "m1": m1,
                    "m2": m2,
                    "m3": m3,
                    "resist": resist,
                    "onset": df["onset_muO_eV"].values if "onset_muO_eV" in df.columns else resist,
                }
            )
        if not systems:
            print("  compiled_ternary_onset_grid: no valid data")
            return

        all_resist = np.concatenate([s["resist"][~np.isnan(s["resist"])] for s in systems])
        vmin, vmax = float(np.nanmin(all_resist)), float(np.nanmax(all_resist))
        if vmax <= vmin:
            vmax = vmin + 1.0
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap_name = "RdYlGn"
        n = len(systems)
        ncols = min(n, 6)
        nrows = (n + ncols - 1) // ncols
        sqrt3_2 = np.sqrt(3) / 2
        tri_verts = np.array([[0.5, sqrt3_2], [0.0, 0.0], [1.0, 0.0], [0.5, sqrt3_2]])
        from matplotlib.gridspec import GridSpec

        fig = plt.figure(figsize=(ncols * 3.0 + 0.7, nrows * 3.4))
        gs = GridSpec(nrows, ncols + 1, figure=fig, width_ratios=[1] * ncols + [0.12], wspace=0.05, hspace=0.25)
        axes = [fig.add_subplot(gs[r, c]) for r in range(nrows) for c in range(ncols)]
        cbar_ax = fig.add_subplot(gs[:, ncols])
        for idx, sys in enumerate(systems):
            ax = axes[idx]
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set_xlim(-0.12, 1.12)
            ax.set_ylim(-0.18, sqrt3_2 + 0.18)
            ax.plot(tri_verts[:, 0], tri_verts[:, 1], "k-", lw=0.8)
            xc = sys["x2"] + 0.5 * sys["x3"]
            yc = sqrt3_2 * sys["x3"]
            r = np.where(np.isnan(sys["resist"]), vmax, sys["resist"])
            from matplotlib.patches import PathPatch
            from matplotlib.path import Path as MplPath
            from matplotlib.tri import Triangulation

            tri = Triangulation(xc, yc)
            hb = ax.tripcolor(tri, r, shading="flat", cmap=cmap_name, norm=norm, edgecolors="black", linewidth=0.3)
            hb.set_clip_path(
                PathPatch(MplPath([(0.0, 0.0), (1.0, 0.0), (0.5, sqrt3_2), (0.0, 0.0)]), transform=ax.transData)
            )
            off, fs = 0.06, 6
            ax.text(0.5, sqrt3_2 + off, sys["m1"], ha="center", va="bottom", fontsize=fs, fontweight="bold")
            ax.text(-off, -off * 0.5, sys["m2"], ha="right", va="top", fontsize=fs, fontweight="bold")
            ax.text(1 + off, -off * 0.5, sys["m3"], ha="left", va="top", fontsize=fs, fontweight="bold")
            ax.text(0.5, -0.13, sys["label"], ha="center", va="top", fontsize=5.5, style="italic")
        for ax in axes[len(systems) :]:
            ax.axis("off")
        sm = cm.ScalarMappable(norm=norm, cmap=cmap_name)
        sm.set_array([])
        fig.colorbar(sm, cax=cbar_ax, label="resistance (eV)", aspect=30)
        fig.suptitle("All ternary systems — oxidation resistance", fontsize=9, y=1.005)
        fig.savefig(cfg.figures_dir / "compiled_ternary_onset_grid.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"  compiled_ternary_onset_grid.png  ({n} systems)")

        try:
            from plotly.subplots import make_subplots

            pfig = make_subplots(
                rows=nrows,
                cols=ncols,
                specs=[[{"type": "ternary"}] * ncols for _ in range(nrows)],
                subplot_titles=[s["label"] for s in systems] + [""] * (nrows * ncols - n),
                horizontal_spacing=0.04,
                vertical_spacing=0.06,
            )
            for idx, sys in enumerate(systems):
                row = idx // ncols + 1
                col = idx % ncols + 1
                r = np.where(np.isnan(sys["resist"]), vmax, sys["resist"])
                pfig.add_trace(
                    go.Scatterternary(
                        a=sys["x1"],
                        b=sys["x2"],
                        c=sys["x3"],
                        mode="markers",
                        marker=dict(
                            size=6,
                            color=r,
                            colorscale="RdYlGn",
                            cmin=vmin,
                            cmax=vmax,
                            showscale=(idx == 0),
                            colorbar=dict(title="resistance (eV)", len=0.4, y=0.8),
                        ),
                        showlegend=False,
                    ),
                    row=row,
                    col=col,
                )
                tern_key = "" if idx == 0 else str(idx + 1)
                pfig.update_layout(
                    **{
                        f"ternary{tern_key}": dict(
                            aaxis=dict(title=sys["m1"], tickfont=dict(size=7)),
                            baxis=dict(title=sys["m2"], tickfont=dict(size=7)),
                            caxis=dict(title=sys["m3"], tickfont=dict(size=7)),
                        )
                    }
                )
            pfig.update_layout(
                title="All ternary — oxidation resistance", height=max(400, nrows * 300), width=max(600, ncols * 250)
            )
            pfig.write_html(str(cfg.figures_dir / "compiled_ternary_onset_grid.html"))
            print("  compiled_ternary_onset_grid.html")
        except Exception as e:
            print(f"  plotly ternary grid skipped: {e}")

        try:
            all_onset = np.concatenate([s["onset"][~np.isnan(s["onset"])] for s in systems])
            o_min, o_max = float(np.nanmin(all_onset)), float(np.nanmax(all_onset))
            p3d = go.Figure()
            for idx, sys in enumerate(systems):
                on_mask = ~np.isnan(sys["onset"])
                if on_mask.sum() < 3:
                    continue
                x1 = sys["x1"][on_mask]
                x2 = sys["x2"][on_mask]
                z = sys["onset"][on_mask]
                p3d.add_trace(
                    go.Mesh3d(
                        x=x1,
                        y=x2,
                        z=z,
                        intensity=z,
                        colorscale="RdYlGn",
                        cmin=o_min,
                        cmax=o_max,
                        reversescale=True,
                        showscale=(idx == 0),
                        opacity=0.75,
                        name=sys["label"],
                        colorbar=dict(title="onset μO<br>(eV/O)", len=0.5),
                        hovertemplate=(
                            f'{sys["label"]}<br>x_{sys["m1"]}=%{{x:.2f}}<br>'
                            f'x_{sys["m2"]}=%{{y:.2f}}<br>onset μO=%{{z:.3f}}<extra></extra>'
                        ),
                    )
                )
            p3d.update_layout(
                title="All ternary — onset μO",
                scene=dict(
                    xaxis_title="x_M1",
                    yaxis_title="x_M2",
                    zaxis_title="onset μO (eV/O)",
                    xaxis=dict(range=[0, 1]),
                    yaxis=dict(range=[0, 1]),
                ),
                width=1000,
                height=800,
                legend=dict(x=1.05, y=1.0),
            )
            p3d.write_html(str(cfg.figures_dir / "compiled_ternary_onset_3d.html"))
            print("  compiled_ternary_onset_3d.html")
        except Exception as e:
            print(f"  plotly 3D skipped: {e}")

    # ---------------------------------------------------------------- miscibility gaps

    def plot_miscibility_gaps(self) -> None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
        import matplotlib.tri as tri_mod

        try:
            from pycalphad import Database, calculate
            from scipy.interpolate import LinearNDInterpolator
            from scipy.spatial import ConvexHull
        except ImportError as e:
            print(f"  plot_miscibility_gaps: skipped ({e})")
            return

        cfg = self.config
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
            for simplex, eq in zip(hull.simplices, hull.equations):
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
            for T, mask in zip(temps, gap_masks):
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
            metals = sorted(
                e for e in tdb.elements if e.upper() not in {f.upper() for f in fixed_elements} and e != "VA"
            )
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

        for tp, tdb, metals in ternary_tdb:
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

    # ---------------------------------------------------------------- combined figures

    def compile_combined_figures(self) -> None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.cm as cm
        import matplotlib.pyplot as plt
        import matplotlib.tri as tri_mod
        import pandas as pd

        try:
            from pycalphad import Database, calculate
            from scipy.interpolate import LinearNDInterpolator
            from scipy.spatial import ConvexHull
        except ImportError as e:
            print(f"  compile_combined_figures: skipped ({e})")
            return

        cfg = self.config
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
            for simplex, eq in zip(hull.simplices, hull.equations):
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
            for T, mask in zip(temps, gap_masks):
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
            metals = sorted(
                e for e in tdb.elements if e.upper() not in {f.upper() for f in fixed_elements} and e != "VA"
            )
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

        for tp, tdb, metals in ternary_tdb:
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
                ax.text(
                    0.5, -0.08, txt, ha="center", va="top", transform=ax.transAxes, fontsize=fontsize, style="italic"
                )

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
                    hb = axes[1].tripcolor(
                        tri, resist_f, shading="flat", cmap="RdYlGn", edgecolors="black", linewidth=0.4
                    )
                    hb.set_clip_path(
                        PathPatch(
                            MplPath([(0.0, 0.0), (1.0, 0.0), (0.5, sqrt3_2), (0.0, 0.0)]), transform=axes[1].transData
                        )
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
        for idx, (tp, tdb, metals) in enumerate(ternary_tdb):
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
                        ax.tricontourf(
                            triang, mask.astype(float), levels=[0.5, 1.5], colors=[cmap_t(norm_t(T))], alpha=0.40
                        )
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
