"""Onset line and ternary onset grid compilation plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .config import Config


def compile_onset_lines(cfg: Config, tags_info: list) -> None:
    """Copy assemblage maps + build combined onset line plots.

    Args:
        cfg: Framework configuration object.
        tags_info: List of (tag, m1, m2, sys_name) tuples from BatchRunner.
    """
    import shutil

    import matplotlib as mpl
    import matplotlib.lines as mlines
    import matplotlib.pyplot as plt
    import pandas as pd
    import plotly.graph_objects as go

    mpl.use("Agg")

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
        for _tag, _m1, _m2, sys_name in tags_info:
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
    for _mode, y_col, csv_suf, title, fname in [
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
        for i, (tag, m1, _m2, sys_name) in enumerate(tags_info):
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
            c_hex = f"#{int(c_mpl[0] * 255):02x}{int(c_mpl[1] * 255):02x}{int(c_mpl[2] * 255):02x}"
            ax2.plot(xs, ys, lw=1.5, color=c_mpl, label=sys_name)
            pfig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="lines",
                    name=sys_name,
                    line={"color": c_hex, "width": 2},
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
        ax2.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, ncol=ncol, framealpha=0.9, borderaxespad=0)
        fig2.tight_layout()
        fig2.savefig(cfg.figures_dir / f"{fname}.png", dpi=150, bbox_inches="tight")
        plt.close(fig2)
        pfig.update_layout(
            title=title,
            width=1100,
            height=700,
            xaxis={"title": "initial x_M1 fraction", "range": [0, 1]},
            yaxis={"title": "onset μO (eV per O atom)"},
            hovermode="closest",
            legend={"x": 1.01, "y": 1, "font": {"size": 9}},
        )
        pfig.write_html(str(cfg.figures_dir / f"{fname}.html"))
        print(f"  {fname}.png + .html  ({n_plotted} systems)")


def compile_ternary_onset_grid(cfg: Config) -> None:
    """Build a grid of ternary onset plots across all systems.

    Args:
        cfg: Framework configuration object.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import pandas as pd
    import plotly.graph_objects as go
    from matplotlib import cm

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
            else df["onset_muO_eV"].values - (df["onset_muO_eV"].min() if df["onset_muO_eV"].notna().any() else -10.0)
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
        hb.set_clip_path(PathPatch(MplPath([(0.0, 0.0), (1.0, 0.0), (0.5, sqrt3_2), (0.0, 0.0)]), transform=ax.transData))
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
                    marker={
                        "size": 6,
                        "color": r,
                        "colorscale": "RdYlGn",
                        "cmin": vmin,
                        "cmax": vmax,
                        "showscale": (idx == 0),
                        "colorbar": {"title": "resistance (eV)", "len": 0.4, "y": 0.8},
                    },
                    showlegend=False,
                ),
                row=row,
                col=col,
            )
            tern_key = "" if idx == 0 else str(idx + 1)
            pfig.update_layout(
                **{
                    f"ternary{tern_key}": {
                        "aaxis": {"title": sys["m1"], "tickfont": {"size": 7}},
                        "baxis": {"title": sys["m2"], "tickfont": {"size": 7}},
                        "caxis": {"title": sys["m3"], "tickfont": {"size": 7}},
                    }
                }
            )
        pfig.update_layout(title="All ternary — oxidation resistance", height=max(400, nrows * 300), width=max(600, ncols * 250))
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
                    colorbar={"title": "onset μO<br>(eV/O)", "len": 0.5},
                    hovertemplate=(
                        f"{sys['label']}<br>x_{sys['m1']}=%{{x:.2f}}<br>"
                        f"x_{sys['m2']}=%{{y:.2f}}<br>onset μO=%{{z:.3f}}<extra></extra>"
                    ),
                )
            )
        p3d.update_layout(
            title="All ternary — onset μO",
            scene={
                "xaxis_title": "x_M1",
                "yaxis_title": "x_M2",
                "zaxis_title": "onset μO (eV/O)",
                "xaxis": {"range": [0, 1]},
                "yaxis": {"range": [0, 1]},
            },
            width=1000,
            height=800,
            legend={"x": 1.05, "y": 1.0},
        )
        p3d.write_html(str(cfg.figures_dir / "compiled_ternary_onset_3d.html"))
        print("  compiled_ternary_onset_3d.html")
    except Exception as e:
        print(f"  plotly 3D skipped: {e}")
