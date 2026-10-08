"""Region-map and side-by-side diagram helpers for SystemAnalyzer."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from ._analyzer_cache import _composition_from_slice
from .utils import system_key


def _ensure_imports() -> None:
    """Add the framework directory and its python subdirectory to sys.path."""
    root = Path(__file__).parent
    for p in [str(root), str(root / "python")]:
        if p not in sys.path:
            sys.path.insert(0, p)


def _find_external_diagram(config, metals, T=None) -> Path | None:
    """Locate an external phase diagram image for the given metal system.

    Searches ``config.ternary_diagrams_root / system_key(metals)`` for
    known filename patterns, preferring animated GIFs, then per-temperature
    PNGs, then generic phase-map PNGs.

    Args:
        config: Config instance providing ``ternary_diagrams_root``.
        metals: Sequence of metal element symbols.
        T: Optional temperature in K; used to select per-temperature images.

    Returns:
        Path to the best matching image, or None if none is found.
    """
    base = config.ternary_diagrams_root / system_key(metals)
    if not base.exists():
        return None
    # Prefer clean GIF
    clean_gif = base / f"{system_key(metals)}_phase_evolution_clean.gif"
    if clean_gif.exists():
        return clean_gif
    if T is not None:
        for p in [
            base / "per_temperature" / f"{int(round(T))}K.png",
            base / "per_temperature" / f"{int(round(T))}K_gibbs.png",
        ]:
            if p.exists():
                return p
    for p in [
        base / f"{system_key(metals)}_phase_map.png",
        base / f"{system_key(metals)}_phase_diagram.png",
        base / f"{system_key(metals)}_ternary.png",
    ]:
        if p.exists():
            return p
    # Case-insensitive fallback: prefer files containing "phase" and "diagram"
    pngs = sorted(base.glob("*.png"))
    for p in pngs:
        n = p.name.lower()
        if "phase" in n and "diagram" in n:
            return p
    for p in pngs:
        if "phase" in p.name.lower():
            return p
    return pngs[0] if pngs else None


def _read_diagram_image(diag: Path, T=None, t_start: float = 200.0, gif_t_step: float = 10.0):
    """Read an image file (PNG or GIF) as a NumPy RGB array.

    For animated GIFs the frame closest to temperature *T* is selected.

    Args:
        diag: Path to a PNG or GIF image.
        T: Optional temperature in K used to select a GIF frame.
        t_start: Temperature corresponding to frame 0 of the GIF (K).
        gif_t_step: Temperature step between consecutive GIF frames (K).

    Returns:
        NumPy uint8 array of shape (H, W, 3).
    """
    if diag.suffix.lower() == ".gif":
        from PIL import Image as _PILImg

        im = _PILImg.open(diag)
        frame_idx = max(0, int(round((float(T) - t_start) / gif_t_step))) if T is not None else getattr(im, "n_frames", 1) - 1
        im.seek(min(frame_idx, getattr(im, "n_frames", 1) - 1))
        return np.array(im.convert("RGB"))
    import matplotlib.pyplot as plt

    return plt.imread(diag)


def _plot_region_map_png(
    config,
    mu_values,
    x_values,
    region_grid,
    region_labels,
    sys_cfg,
    T,
    spec,
    out_dir,
    region_details=None,
    cell_text_grid=None,
    y_label=None,
    y_axis_label=None,
) -> Path:
    """Render and save a coloured assemblage-region map as PNG (and HTML).

    Args:
        config: Config instance providing plot styling options.
        mu_values: 1-D array of μO values (x-axis of the map).
        x_values: 1-D array of composition or temperature values (y-axis).
        region_grid: 2-D integer array (len(x_values) × len(mu_values)) of
            region IDs.
        region_labels: List of region-label strings indexed by region ID − 1.
        sys_cfg: SystemConfig for axis labels and titles.
        T: Temperature in K (displayed in the title), or None.
        spec: Slice specification dict with keys ``axis`` and ``label``.
        out_dir: Directory where output files are written (created if absent).
        region_details: Optional dict mapping region ID → detail dict.
        cell_text_grid: Optional 2-D object array of per-cell annotation strings.
        y_label: Override for the y-axis label shown on the plot.
        y_axis_label: Unused; reserved for future use.

    Returns:
        Path to the saved ``assemblage_region_map.png``.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import hashlib as _hl

    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    _ensure_imports()
    from .thermodynamics import (
        add_region_annotation,
        format_phase_detail_line,
        grid_edges,
        region_annotation_text,
        separate_region_annotations,
        write_region_map_html,
    )

    cfg = config
    out_dir.mkdir(parents=True, exist_ok=True)
    cmap = plt.get_cmap("tab20")

    def _lc(lbl):
        return cmap(int(_hl.md5(lbl.encode()).hexdigest(), 16) % 20 / 20)

    colors = {r + 1: _lc(lbl) for r, lbl in enumerate(region_labels)}
    xe = grid_edges(x_values)
    me = grid_edges(mu_values)

    if cfg.region_label_mode == "id":
        fig, (ax, ax_leg) = plt.subplots(1, 2, figsize=(17, 8), gridspec_kw={"width_ratios": [3.2, 1.2], "wspace": 0.05})
        ax_leg.axis("off")
    else:
        fig, ax = plt.subplots(figsize=(13, 8))
        ax_leg = None
    for ix in range(len(x_values)):
        for imu in range(len(mu_values)):
            rid = int(region_grid[ix, imu])
            ax.add_patch(
                plt.Rectangle(
                    (me[imu], xe[ix]),
                    me[imu + 1] - me[imu],
                    xe[ix + 1] - xe[ix],
                    color=colors.get(rid, (0.9, 0.9, 0.9, 1.0)),
                    linewidth=0,
                )
            )
    for ix in range(len(x_values)):
        for imu in range(len(mu_values) - 1):
            if region_grid[ix, imu] != region_grid[ix, imu + 1]:
                ax.plot([me[imu + 1]] * 2, [xe[ix], xe[ix + 1]], "k-", lw=cfg.boundary_lw)
    for ix in range(len(x_values) - 1):
        for imu in range(len(mu_values)):
            if region_grid[ix, imu] != region_grid[ix + 1, imu]:
                ax.plot([me[imu], me[imu + 1]], [xe[ix + 1]] * 2, "k-", lw=cfg.boundary_lw)
    annotation_artists = []
    for rid, _ in enumerate(region_labels, start=1):
        mask = region_grid == rid
        if not np.any(mask):
            continue
        rows, cols = np.where(mask)
        r = max(0, min(len(x_values) - 1, int(np.round(np.median(rows)))))
        c = max(0, min(len(mu_values) - 1, int(np.round(np.median(cols)))))
        annotation = region_annotation_text(
            rid,
            region_labels[rid - 1],
            region_details,
            cfg.region_label_mode,
        )
        annotation_artists.append(
            add_region_annotation(
                ax,
                mu_values[c],
                x_values[r],
                annotation,
                cfg.region_label_fontsize,
                colors[rid],
            )
        )
    axis_metal = sys_cfg.metals[spec["axis"]] if sys_cfg else ""
    resolved_y_label = y_label if y_label else f"initial x_{axis_metal}"
    ax.set_xlim(float(mu_values[0]), float(mu_values[-1]))
    ax.set_ylim(float(x_values[0]), float(x_values[-1]))
    ax.set_xlabel(r"$\mu_O$ (eV per O atom)", fontsize=13)
    ax.set_ylabel(resolved_y_label, fontsize=13)
    title_T = f" | T={T} K" if T is not None else ""
    ax.set_title(f"{sys_cfg.phase_label} | {spec['label']}{title_T}", fontsize=13)
    separate_region_annotations(ax, annotation_artists)

    def _format_detail(detail):
        if not detail:
            return ""
        phase_ranges = detail.get("phase_ranges") or []
        if phase_ranges:
            return format_phase_detail_line(phase_ranges)
        exact = detail.get("exact_label")
        return exact or ""

    handles = []
    for rid, lbl in enumerate(region_labels, start=1):
        if not np.any(region_grid == rid):
            continue
        detail = _format_detail(region_details.get(rid) if region_details else None)
        label = f"{rid}: {detail if detail else lbl}"
        handles.append(mpatches.Patch(color=colors[rid], label=label))
    fs = max(5, min(8, int(110 / max(len(handles), 1))))
    if ax_leg is not None:
        ax_leg.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            fontsize=fs,
            framealpha=0.9,
            handlelength=1.2,
            borderaxespad=0,
            title="Region key",
            title_fontsize=fs + 1,
        )
    out = out_dir / "assemblage_region_map.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    try:
        write_region_map_html(
            mu_values,
            x_values,
            region_grid,
            region_labels,
            out_dir / "assemblage_region_map.html",
            title=f"{sys_cfg.phase_label} | {spec['label']}{title_T}",
            x_label="mu_O (eV per O atom)",
            y_label=resolved_y_label,
            region_details=region_details,
            cell_text_grid=cell_text_grid,
            region_label_mode=cfg.region_label_mode,
            region_label_fontsize=cfg.region_label_fontsize,
        )
    except Exception as e:
        print(f"[html skipped] {e}")
    plt.close(fig)
    return out


def _plot_side_by_side(config, src_map, metals, T, spec, out_dir, composition=None) -> Path | None:
    """Place the region map beside an external miscibility diagram and save.

    Args:
        config: Config instance providing ``make_side_by_side_miscibility``,
            ``slice_x_values``, ``phase_diagram_t_start``,
            ``phase_diagram_t_end``, and ``phase_element``.
        src_map: Path to the assembled region-map PNG.
        metals: Sequence of metal element symbols.
        T: Temperature in K, or None.
        spec: Slice specification dict with keys ``axis``, ``remainder``,
            and ``label``.
        out_dir: Directory where the output PNG is saved.
        composition: Optional composition array used to draw a vertical line on
            binary diagrams.

    Returns:
        Path to the saved ``assemblage_with_miscibility.png``, or None when the
        side-by-side is disabled or the source map does not exist.
    """
    if len(metals) < 2 or not config.make_side_by_side_miscibility:
        return None
    diag = _find_external_diagram(config, metals, T)
    if not src_map.exists():
        return None
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), gridspec_kw={"wspace": 0.04})
    axes[0].imshow(plt.imread(src_map))
    axes[0].axis("off")
    if diag is not None:
        is_gif = Path(diag).suffix.lower() == ".gif"
        if len(metals) == 3:
            # Ternary: show diagram with composition scan-line overlay
            sqrt3_2 = np.sqrt(3.0) / 2.0
            axes[1].imshow(_read_diagram_image(diag, T), extent=(0, 1, 0, sqrt3_2), origin="upper")
            axes[1].set_xlim(0, 1)
            axes[1].set_ylim(0, sqrt3_2)
            axes[1].axis("off")
            x0 = float(config.slice_x_values[0])
            x1 = float(config.slice_x_values[-1])
            c0 = _composition_from_slice(x0, 3, spec["axis"], spec["remainder"])
            c1 = _composition_from_slice(x1, 3, spec["axis"], spec["remainder"])
            xa, ya = c0[1] + 0.5 * c0[2], sqrt3_2 * c0[2]
            xb, yb = c1[1] + 0.5 * c1[2], sqrt3_2 * c1[2]
            axes[1].plot([xa, xb], [ya, yb], color="black", lw=5, alpha=0.85)
            axes[1].plot([xa, xb], [ya, yb], color="white", lw=2.5, alpha=0.95)
            off = 0.05
            axes[1].text(0.50, sqrt3_2 + off, metals[2], ha="center", va="bottom", fontsize=12, fontweight="bold")
            axes[1].text(-off, -off * 0.4, metals[0], ha="right", va="top", fontsize=12, fontweight="bold")
            axes[1].text(1.0 + off, -off * 0.4, metals[1], ha="left", va="top", fontsize=12, fontweight="bold")
        else:
            # Binary (or higher-order): show diagram with T line + element labels.
            # BLADEvisualizer: x-axis = x_{M2} (0=pure M1 left, 1=pure M2 right)
            # y-axis: Temperature (K), bottom=T_min, top=T_max
            # data x range ≈ transAxes [0.12, 0.92], y range ≈ [0.12, 0.88]
            axes[1].imshow(_read_diagram_image(diag, T), aspect="auto")
            axes[1].axis("off")
            if len(metals) == 2:
                M1, M2 = metals[0], metals[1]
                # Element labels at x-axis ends (x_M2: left=0=M1, right=1=M2)
                x_left, x_right = 0.12, 0.92  # data x bounds in transAxes
                y_xlabel = 0.08  # just below data area
                for x_pos, ha, label in [
                    (x_left, "center", M1),
                    (x_right, "center", M2),
                ]:
                    axes[1].text(
                        x_pos,
                        y_xlabel,
                        label,
                        transform=axes[1].transAxes,
                        ha=ha,
                        va="top",
                        fontsize=13,
                        fontweight="bold",
                        bbox={"fc": "white", "ec": "none", "alpha": 0.7, "pad": 1},
                    )
                if not is_gif and T is not None:
                    t0 = float(getattr(config, "phase_diagram_t_start", 300))
                    t1 = float(getattr(config, "phase_diagram_t_end", 4500))
                    y_data_lo, y_data_hi = 0.12, 0.88
                    if t1 > t0:
                        raw = min(1.0, max(0.0, (float(T) - t0) / (t1 - t0)))
                        y_line = y_data_lo + raw * (y_data_hi - y_data_lo)
                        axes[1].plot(
                            [x_left, x_right],
                            [y_line, y_line],
                            transform=axes[1].transAxes,
                            color="black",
                            lw=2.8,
                            solid_capstyle="round",
                        )
                        axes[1].plot(
                            [x_left, x_right],
                            [y_line, y_line],
                            transform=axes[1].transAxes,
                            color="white",
                            lw=1.2,
                            solid_capstyle="round",
                        )
                if composition is not None:
                    comp = np.asarray(composition, dtype=float)
                    if comp.shape == (2,):
                        x_second = min(1.0, max(0.0, float(comp[1])))
                        x_line = x_left + x_second * (x_right - x_left)
                        y_data_lo, y_data_hi = 0.12, 0.88
                        axes[1].plot(
                            [x_line, x_line],
                            [y_data_lo, y_data_hi],
                            transform=axes[1].transAxes,
                            color="black",
                            lw=3.2,
                            solid_capstyle="round",
                        )
                        axes[1].plot(
                            [x_line, x_line],
                            [y_data_lo, y_data_hi],
                            transform=axes[1].transAxes,
                            color="white",
                            lw=1.4,
                            solid_capstyle="round",
                        )
        # Temperature stamp on all diagram types
        axes[1].text(
            0.50,
            0.02,
            (f"T = {int(T)} K  |  " if T is not None else "") + spec["label"],
            transform=axes[1].transAxes,
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            bbox={"fc": "white", "ec": "none", "alpha": 0.8, "pad": 2},
        )
    else:
        axes[1].axis("off")
        axes[1].text(
            0.5,
            0.5,
            "No external miscibility diagram found",
            transform=axes[1].transAxes,
            ha="center",
            va="center",
            fontsize=11,
        )
    fig.suptitle(f"{''.join(metals)}{config.phase_element or ''} | {spec['label']}", fontsize=13)
    out = out_dir / "assemblage_with_miscibility.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out
