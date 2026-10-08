"""Plotting helpers for SingleCompositionAnalyzer.

Contains _plot_muO_T — the μO vs T assemblage map renderer —
extracted from single.py to keep that module slim.
"""

from __future__ import annotations

import numpy as np


def _plot_muO_T(
    self,
    mu_vals,
    T_vals,
    region_grid,
    rlabels,
    region_details,
    cell_texts,
    out_dir,
    sys_str,
    M1,
    M2,
    boundary_lw,
    region_label_fontsize,
    phase_ids,
    phase_kinds,
) -> None:
    """Render and save the μO vs T phase-assemblage map.

    Args:
        self: SingleCompositionAnalyzer instance.
        mu_vals: 1-D array of μO grid values (eV/O).
        T_vals: 1-D array of temperature grid values (K).
        region_grid: 2-D integer array of region IDs, shape (n_T, n_mu).
        rlabels: List of assemblage-label strings, indexed by region_id - 1.
        region_details: Dict mapping region_id → detail dict (phase ranges, etc.).
        cell_texts: 2-D object array of per-cell phase-fraction strings.
        out_dir: Path to directory where outputs are written.
        sys_str: System label string for the figure title.
        M1: First metal symbol.
        M2: Second metal symbol (empty string if binary with one metal).
        boundary_lw: Line width for region-boundary lines.
        region_label_fontsize: Font size for in-map region labels.
        phase_ids: Sequence of phase identifier strings/integers.
        phase_kinds: Sequence of phase kind strings ('phase', 'fixed', …).
    """
    import hashlib as _hl

    import matplotlib as mpl
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    from .thermodynamics import (
        add_region_annotation,
        format_phase_detail_line,
        grid_edges,
        region_annotation_text,
        separate_region_annotations,
        write_region_map_html,
    )

    mpl.use("Agg")
    out_dir.mkdir(parents=True, exist_ok=True)
    _cmap = plt.get_cmap("tab20")
    # Deterministic hash color — same assemblage always same color
    colors = {r + 1: _cmap(int(_hl.md5(lbl.encode()).hexdigest(), 16) % 20 / 20) for r, lbl in enumerate(rlabels)}
    mu_edges = grid_edges(mu_vals)
    T_edges = grid_edges(T_vals)

    if self.region_label_mode == "id":
        fig, (ax, ax_leg) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05})
        ax_leg.axis("off")
    else:
        fig, ax = plt.subplots(figsize=(13, 8))
        ax_leg = None

    for iT in range(len(T_vals)):
        for imu in range(len(mu_vals)):
            rid = region_grid[iT, imu]
            ax.add_patch(
                plt.Rectangle(
                    (mu_edges[imu], T_edges[iT]),
                    mu_edges[imu + 1] - mu_edges[imu],
                    T_edges[iT + 1] - T_edges[iT],
                    color=colors.get(rid, (0.85, 0.85, 0.85, 1.0)),
                    linewidth=0,
                )
            )
    for iT in range(len(T_vals)):
        for imu in range(len(mu_vals) - 1):
            if region_grid[iT, imu] != region_grid[iT, imu + 1]:
                ax.plot([mu_edges[imu + 1]] * 2, [T_edges[iT], T_edges[iT + 1]], "k-", lw=boundary_lw)
    for iT in range(len(T_vals) - 1):
        for imu in range(len(mu_vals)):
            if region_grid[iT, imu] != region_grid[iT + 1, imu]:
                ax.plot([mu_edges[imu], mu_edges[imu + 1]], [T_edges[iT + 1]] * 2, "k-", lw=boundary_lw)
    annotation_artists = []
    for rid, lbl in enumerate(rlabels, start=1):
        mask = region_grid == rid
        if not np.any(mask):
            continue
        rows, cols = np.where(mask)
        r = max(0, min(len(T_vals) - 1, int(np.round(np.median(rows)))))
        c = max(0, min(len(mu_vals) - 1, int(np.round(np.median(cols)))))
        annotation = region_annotation_text(rid, lbl, region_details, self.region_label_mode)
        region_color = colors[rid]

        lightness = 0.5  # 0.0 = unchanged, larger = lighter
        annotation_color = (
            region_color[0] + (1.0 - region_color[0]) * lightness,
            region_color[1] + (1.0 - region_color[1]) * lightness,
            region_color[2] + (1.0 - region_color[2]) * lightness,
            region_color[3],
        )

        annotation_artists.append(
            add_region_annotation(
                ax,
                mu_vals[c],
                T_vals[r],
                annotation,
                9,
                annotation_color,
            )
        )
    ax.set_xlim(mu_vals[0], mu_vals[-1])
    ax.set_ylim(T_vals[0], T_vals[-1])
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    ax.set_xlabel(r"$\mu_O$ (eV per O atom)", fontsize=20)
    ax.set_ylabel("T (K)", fontsize=20)
    ax.set_title(f"{sys_str} | {self._x_label}", fontsize=24)
    separate_region_annotations(ax, annotation_artists)

    def _fmt(detail):
        if not detail:
            return ""
        pr = detail.get("phase_ranges") or []
        return format_phase_detail_line(pr) if pr else (detail.get("exact_label") or "")

    handles = []
    for rid, lbl in enumerate(rlabels, start=1):
        if not np.any(region_grid == rid):
            continue
        detail = _fmt(region_details.get(rid) if region_details else None)
        handles.append(mpatches.Patch(color=colors[rid], label=f"{rid}: {detail if detail else lbl}"))
    leg_fs = max(5, min(8, int(110 / max(len(handles), 1))))
    if ax_leg is not None:
        ax_leg.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            fontsize=leg_fs,
            framealpha=0.9,
            handlelength=1.2,
            borderaxespad=0,
            title="Region key",
            title_fontsize=leg_fs + 5,
        )

    png = out_dir / "muO_T_map.png"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    try:
        write_region_map_html(
            mu_vals,
            T_vals,
            region_grid,
            rlabels,
            out_dir / "muO_T_map.html",
            title=f"{sys_str} | {self._x_label}",
            x_label="mu_O (eV/O)",
            y_label="T (K)",
            region_details=region_details,
            cell_text_grid=cell_texts,
            region_label_mode=self.region_label_mode,
            region_label_fontsize=region_label_fontsize,
        )
    except Exception as e:
        print(f"  [html skipped] {e}")
    plt.close(fig)
    print(f"  Figure → {png}")
