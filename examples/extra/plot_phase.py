"""Plot N-component phase diagrams with shaded two-phase regions from BLADE TDB files.

For each TDB, computes G_mix on a composition grid, finds the lower
convex hull to identify two-phase (miscibility gap) regions, shades them, and
overlays results for multiple temperatures — all on one plot.

Works for binary, ternary, or higher-order systems with multi-sublattice phases
(e.g. diborides with fixed B sublattice).

Usage:
    pixi run python examples/extra/plot_phase.py
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/private/tmp") / "blade_matplotlib"))

import imageio.v2 as imageio
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from pycalphad import Database, calculate
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull
from blade.analysis.blade_visual import BladeVisualizer

viz = BladeVisualizer()

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"

comps_dir  = path1 / "Files" / "Comps"
output_dir = path1 / "Files" / "Phase_Diagrams"
output_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------
t_start = 200       # K
t_end   = 2000      # K
t_step  = 200       # K increment per overlay/per-temperature plot

# True = plot only from existing CSV data; False = use CSV data and calculate missing temperatures.
skip = False

# GIF/movie settings
make_gif   = True   # save animated GIF of diagram evolving over temperature
gif_t_step = 10  # K per GIF frame; set smaller for smoother GIFs
gif_fps    = 20     # frames per second

# Fixed sublattice elements (excluded from composition axes)
fixed_elements: set[str] = {"B"}

# Metal sublattice total fraction: 1/3 for MB2, 1.0 for pure alloy
metal_fraction = 1 / 3

# Composition grid density — more points = sharper boundaries but slower
n_grid = 40

# Minimum G above convex hull (J/mol) to classify as two-phase
two_phase_threshold = 5.0

pressure = 101325.0

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def get_variable_elements(tdb: Database, fixed: set[str]) -> list[str]:
    return sorted(
        el for el in tdb.elements
        if el.upper() not in {f.upper() for f in fixed} and el != "VA"
    )


def generate_polygon_vertices(n: int, radius: float = 1.0, center=(0, 0)) -> np.ndarray:
    """Generate vertices of a regular n-gon.

    Matches the notebook's generate_regular_pentagon logic:
      n=4: uses π/n start so edges are horizontal/vertical (square, not diamond)
      other n: uses π/2 start (apex at top)
    """
    if n == 4:
        angles = np.linspace(np.pi / n, 2 * np.pi + np.pi / n, n, endpoint=False)
    else:
        angles = np.linspace(np.pi / 2, 2 * np.pi + np.pi / 2, n, endpoint=False)
    return np.array([[center[0] + radius * np.cos(a),
                      center[1] + radius * np.sin(a)] for a in angles])


def barycentric_to_cartesian(barycentric_coords: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """Project N-D barycentric coords to 2D via dot product (from notebook)."""
    return np.dot(barycentric_coords, vertices)


def simplex_grid(n_metals: int, n: int, x_total: float) -> np.ndarray:
    """Grid of compositions summing to x_total on N-metal simplex.

    Uses recursive generation — O(C(n+n_metals-1, n_metals-1)) iterations
    instead of O((n+1)^n_metals). Scales to any n_metals.
    """
    step = x_total / n
    pts: list[tuple] = []

    def _recurse(remaining: int, dims_left: int, current: list) -> None:
        if dims_left == 1:
            pts.append(tuple(current + [remaining * step]))
            return
        for i in range(remaining + 1):
            _recurse(remaining - i, dims_left - 1, current + [i * step])

    _recurse(n, n_metals, [])
    return np.array(pts)


def compositions_to_xy(grid: np.ndarray, x_total: float) -> tuple[np.ndarray, np.ndarray]:
    """Project N-metal composition grid to 2D via affine polygon projection."""
    n_metals = grid.shape[1]
    vertices = generate_polygon_vertices(n_metals)
    # Normalize to barycentric (fractions summing to 1)
    bary = grid / x_total
    xy = barycentric_to_cartesian(bary, vertices)
    return xy[:, 0], xy[:, 1]


def draw_polygon_frame(ax, labels: list[str]) -> np.ndarray:
    """Draw polygon frame: gainsboro fill, black outline, labels outside each vertex.

    Matches notebook approach — filled Polygon patch behind data, labels at
    radius_adjustment=0.2 outside every vertex (no special-casing for top).
    """
    from matplotlib.patches import Polygon as MplPolygon
    n = len(labels)
    vertices = generate_polygon_vertices(n)

    # Gray filled polygon behind data (same as notebook Polygon patch)
    poly_patch = MplPolygon(vertices, closed=True,
                            edgecolor='k', facecolor='none',
                            linewidth=2, zorder=0)
    ax.add_patch(poly_patch)

    # Labels outside each vertex — radius_adjustment=0.2 (matches notebook)
    for j, label in enumerate(labels):
        theta = np.arctan2(vertices[j, 1], vertices[j, 0])
        radius_adjustment = 0.2
        nx = vertices[j, 0] + radius_adjustment * np.cos(theta)
        ny = vertices[j, 1] + radius_adjustment * np.sin(theta)
        ax.text(nx, ny, label, fontsize=12, ha='center', va='center',
                fontweight='bold')

    ax.set_aspect('equal')
    ax.axis('off')
    # Axis limits: accommodate labels + room for title above
    pad = 0.35
    ax.set_xlim(vertices[:, 0].min() - pad, vertices[:, 0].max() + pad)
    ax.set_ylim(vertices[:, 1].min() - pad, vertices[:, 1].max() + pad + 0.2)
    return vertices


def temperature_values(step: int) -> np.ndarray:
    temperatures = np.arange(t_start, t_end + 1, step, dtype=int)
    if temperatures.size == 0 or temperatures[-1] != t_end:
        temperatures = np.append(temperatures, t_end)
    return temperatures


def _gmix_cache_path(data_dir: Path) -> Path:
    return data_dir / "gmix_data.csv"


def _load_cached_gmix(data_dir: Path, grid, x_cart, y_cart) -> dict[int, np.ndarray]:
    cache_path = _gmix_cache_path(data_dir)
    if not cache_path.exists():
        return {}

    data = np.genfromtxt(cache_path, delimiter=",", names=True)
    if data.size == 0:
        return {}
    if data.shape == ():
        data = data.reshape(1)

    n_metals = grid.shape[1]
    x_cols = [f"x{i+1}" for i in range(n_metals)]
    required_cols = {"temperature_K", "x_cart", "y_cart", "g_mix"} | set(x_cols)
    cols = set(data.dtype.names or [])
    missing_cols = required_cols - cols
    if missing_cols:
        raise ValueError(f"missing columns: {sorted(missing_cols)}")

    cached: dict[int, np.ndarray] = {}
    for raw_T in np.unique(data["temperature_K"]):
        rows = data[np.isclose(data["temperature_K"], raw_T)]
        T = int(round(float(raw_T)))

        if len(rows) != len(grid):
            raise ValueError(f"{T} K has {len(rows)} rows, expected {len(grid)}")

        cached_grid = np.column_stack([rows[c] for c in x_cols])
        cached_xy = np.column_stack([rows["x_cart"], rows["y_cart"]])
        current_xy = np.column_stack([x_cart, y_cart])
        if not np.allclose(cached_grid, grid, rtol=0.0, atol=1e-12):
            raise ValueError(f"{T} K grid does not match current n_grid/metal_fraction settings")
        if not np.allclose(cached_xy, current_xy, rtol=0.0, atol=1e-12):
            raise ValueError(f"{T} K plot coordinates do not match current settings")

        cached[T] = np.asarray(rows["g_mix"], dtype=float)

    return cached


def _write_cached_gmix(
    data_dir: Path,
    cached: dict[int, np.ndarray],
    grid,
    x_cart,
    y_cart,
) -> None:
    if not cached:
        return

    data_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for T in sorted(cached):
        g_mix = cached[T]
        if len(g_mix) != len(grid):
            continue
        n_metals = grid.shape[1]
        cols = [np.full(len(grid), int(T), dtype=int)]
        cols += [grid[:, i] for i in range(n_metals)]
        cols += [x_cart, y_cart, g_mix]
        rows.append(np.column_stack(cols))

    if not rows:
        return

    n_metals = grid.shape[1]
    x_headers = ",".join(f"x{i+1}" for i in range(n_metals))
    header = f"temperature_K,{x_headers},x_cart,y_cart,g_mix"
    fmt = ["%d"] + ["%.16g"] * (n_metals + 3)
    np.savetxt(_gmix_cache_path(data_dir), np.vstack(rows), delimiter=",", header=header, comments="", fmt=fmt)


def _gm_values_by_temperature(result, n_points: int) -> tuple[np.ndarray, np.ndarray]:
    gm = result.GM
    indexer = {dim: 0 for dim in gm.dims if dim not in {"T", "points"}}
    if indexer:
        gm = gm.isel(indexer)
    gm = gm.transpose("T", "points")
    temperatures = np.asarray(gm.coords["T"].values, dtype=float)
    values = np.asarray(gm.values, dtype=float)[:, :n_points]
    return temperatures, values


def compute_gmix(tdb, all_comps, metals, phases, T, grid, x_fixed_val) -> np.ndarray:
    """Compute G_mix relative to pure endmembers on full grid via calculate. Works for N metals."""
    fixed_els = [el for el in tdb.elements
                 if el.upper() in {f.upper() for f in fixed_elements} and el != "VA"]
    n_metals = len(metals)
    n_pts = len(grid)

    pts = {metals[i]: grid[:, i] for i in range(n_metals)}
    for el in fixed_els:
        pts[el] = np.full(n_pts, x_fixed_val)

    result = calculate(tdb, all_comps, phases[0], T=float(T), P=pressure, points=pts)
    gm = result.GM.values.ravel()[:n_pts]

    # Pure endmember GMs — one per metal, all others at epsilon
    g_ends = []
    for k, m in enumerate(metals):
        end_pts = {metals[i]: np.array([1e-10]) for i in range(n_metals)}
        end_pts[m] = np.array([metal_fraction - (n_metals - 1) * 1e-10])
        for el in fixed_els:
            end_pts[el] = np.array([x_fixed_val])
        r = calculate(tdb, all_comps, phases[0], T=float(T), P=pressure, points=end_pts)
        g_ends.append(float(r.GM.values.ravel()[0]))

    # G_mix = G - sum(xᵢ/x_total * Gᵢ_endmember)
    g_mix = np.array([
        gm[i] - sum(grid[i, j] / metal_fraction * g_ends[j] for j in range(n_metals))
        for i in range(n_pts)
    ])
    return g_mix


def compute_gmix_many(tdb, all_comps, metals, phases, temperatures, grid, x_fixed_val) -> dict[int, np.ndarray]:
    """Compute G_mix for several temperatures in one pycalphad call when possible."""
    temperatures = np.asarray(temperatures, dtype=float)
    if temperatures.size == 0:
        return {}

    fixed_els = [el for el in tdb.elements
                 if el.upper() in {f.upper() for f in fixed_elements} and el != "VA"]
    n_metals = len(metals)
    n_pts = len(grid)

    pts = {metals[i]: grid[:, i] for i in range(n_metals)}
    for el in fixed_els:
        pts[el] = np.full(n_pts, x_fixed_val)

    try:
        result = calculate(tdb, all_comps, phases[0], T=temperatures, P=pressure, points=pts)
        calc_temperatures, gm = _gm_values_by_temperature(result, len(grid))

        g_ends = []
        for m in metals:
            end_pts = {metals[i]: np.array([1e-10]) for i in range(n_metals)}
            end_pts[m] = np.array([metal_fraction - (n_metals - 1) * 1e-10])
            for el in fixed_els:
                end_pts[el] = np.array([x_fixed_val])
            r = calculate(tdb, all_comps, phases[0], T=temperatures, P=pressure, points=end_pts)
            _, end_gm = _gm_values_by_temperature(r, 1)
            g_ends.append(end_gm[:, 0])

        weights = grid / metal_fraction
        endmember_gm = np.column_stack(g_ends)
        baseline = endmember_gm @ weights.T
        g_mix = gm - baseline
        return {int(round(T)): g_mix[i] for i, T in enumerate(calc_temperatures)}
    except Exception as e:
        print(f"    Batch calculation failed ({e}); falling back to per-temperature calculations.")

    computed = {}
    for T in temperatures:
        try:
            computed[int(round(float(T)))] = compute_gmix(tdb, all_comps, metals, phases, T, grid, x_fixed_val)
        except Exception as e:
            print(f"    {int(T)} K calculation failed: {e}")
    return computed


def load_or_compute_gmix(
    tdb,
    all_comps,
    metals,
    phases,
    temperatures,
    grid,
    x_cart,
    y_cart,
    x_fixed_val,
    data_dir: Path,
) -> dict[int, np.ndarray]:
    requested = [int(T) for T in temperatures]
    cache_path = _gmix_cache_path(data_dir)

    try:
        cached = _load_cached_gmix(data_dir, grid, x_cart, y_cart)
    except ValueError as e:
        cached = {}
        print(f"    Ignoring incompatible CSV cache ({cache_path.name}): {e}")
        if skip:
            print("    skip=True, so no recalculation will be attempted.")

    gmix_by_t = {T: cached[T] for T in requested if T in cached}
    missing = [T for T in requested if T not in cached]

    if missing and skip:
        missing_str = ", ".join(f"{T}K" for T in missing)
        print(f"    Missing cached CSV data for: {missing_str}")
        return gmix_by_t

    if missing:
        print(f"    Calculating {len(missing)} temperature(s); saving CSV cache to {cache_path}")
        computed = compute_gmix_many(tdb, all_comps, metals, phases, missing, grid, x_fixed_val)
        cached.update(computed)
        gmix_by_t.update({T: computed[T] for T in missing if T in computed})
        _write_cached_gmix(data_dir, cached, grid, x_cart, y_cart)
    elif gmix_by_t:
        print(f"    Loaded CSV cache: {cache_path}")

    return gmix_by_t


def find_two_phase_mask(x_cart, y_cart, g_mix, threshold) -> np.ndarray:
    """Points above lower convex hull by > threshold are in two-phase region."""
    pts_3d = np.column_stack([x_cart, y_cart, g_mix])
    pts_2d = np.column_stack([x_cart, y_cart])

    try:
        hull = ConvexHull(pts_3d)
    except Exception:
        return np.zeros(len(g_mix), dtype=bool)

    # Collect lower hull vertices (facet normals with negative z-component)
    lower_verts = set()
    for simplex, eq in zip(hull.simplices, hull.equations):
        if eq[2] < 0:
            lower_verts.update(simplex)

    if len(lower_verts) < 3:
        return np.zeros(len(g_mix), dtype=bool)

    lower_idx = np.array(list(lower_verts))
    interp = LinearNDInterpolator(pts_3d[lower_idx, :2], pts_3d[lower_idx, 2])
    g_convex = interp(pts_2d)

    # Fill NaN at hull boundary with the G_mix value (not above hull)
    nan_mask = np.isnan(g_convex)
    g_convex[nan_mask] = g_mix[nan_mask]

    gap = g_mix - g_convex
    return gap > threshold


def _mask_outside_simplex(triang, x_cart, y_cart, x_total: float, n_metals: int):
    """Mask Delaunay triangles whose centroid lies outside the simplex polygon.

    Uses matplotlib Path.contains_points — works for any n_metals without the
    underdetermined barycentric issue that occurs for n≥4 with lstsq.
    """
    from matplotlib.path import Path
    cx = x_cart[triang.triangles].mean(axis=1)
    cy = y_cart[triang.triangles].mean(axis=1)
    vertices = generate_polygon_vertices(n_metals)
    closed = np.vstack([vertices, vertices[0]])
    poly_path = Path(closed)
    inside = poly_path.contains_points(np.column_stack([cx, cy]), radius=1e-6)
    triang.set_mask(~inside)
    return triang


def _clip_to_simplex(ax, n_metals: int) -> None:
    """Clip all ax collections and lines to the simplex polygon boundary."""
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    verts = generate_polygon_vertices(n_metals)
    closed = np.vstack([verts, verts[0]])
    clip_path = Path(closed)
    clip_patch = PathPatch(clip_path, transform=ax.transData)
    for coll in ax.collections:
        coll.set_clip_path(clip_patch)
    for line in ax.lines:
        line.set_clip_path(clip_patch)


def _plot_single_phase(ax, triang, g_mix, x_cart, y_cart, metals, T, color=None):
    """Plot one temperature's two-phase shading onto ax. Returns True if two-phase found."""
    two_phase = find_two_phase_mask(x_cart, y_cart, g_mix, two_phase_threshold)
    c = color or "steelblue"
    if two_phase.any():
        shade_vals = two_phase.astype(float)
        ax.tricontourf(triang, shade_vals, levels=[0.5, 1.5], colors=[c], alpha=0.35)
        ax.tricontour(triang, shade_vals, levels=[0.5], colors=[c], linewidths=1.2, alpha=0.9)
        _clip_to_simplex(ax, len(metals))
    return bool(two_phase.any())


def plot_phase_overlay(
    tdb,
    metals,
    phases,
    output_path,
    per_t_dir: Path | None = None,
    data_dir: Path | None = None,
):
    all_comps = list(tdb.elements) + ["VA"]
    temperatures = temperature_values(t_step)
    cmap = cm.plasma
    norm = plt.Normalize(vmin=t_start, vmax=t_end)

    phase_grid = simplex_grid(len(metals), n_grid, metal_fraction)
    x_cart, y_cart = compositions_to_xy(phase_grid, metal_fraction)
    x_fixed_val = 1.0 - metal_fraction
    triang = tri.Triangulation(x_cart, y_cart)
    # Mask spurious Delaunay triangles that fall outside the simplex.
    # These form at boundary edges and create artifacts in the phase map.
    triang = _mask_outside_simplex(triang, x_cart, y_cart, metal_fraction, len(metals))
    label_str = '-'.join(m.title() for m in metals)
    data_dir = data_dir or output_path.parent / "gmix_data"
    gmix_by_t = load_or_compute_gmix(
        tdb, all_comps, metals, phases, temperatures, phase_grid, x_cart, y_cart, x_fixed_val, data_dir
    )

    fig_ov, ax_ov = plt.subplots(figsize=(9, 8))
    draw_polygon_frame(ax_ov, [m.title() for m in metals])
    any_two_phase = False
    plotted_any = False

    for T in temperatures:
        T = int(T)
        g_mix = gmix_by_t.get(T)
        if g_mix is None:
            print(f"    {T} K: no CSV data")
            continue

        color = cmap(norm(T))
        try:
            if np.isnan(g_mix).all():
                print(f"    {int(T)} K: all NaN")
                continue

            plotted_any = True
            found = _plot_single_phase(ax_ov, triang, g_mix, x_cart, y_cart, metals, T, color)
            if found:
                any_two_phase = True
                print(f"    {int(T)} K: two-phase ({(g_mix > 0).sum()} pts above hull)")
            else:
                print(f"    {int(T)} K: single phase")

            # Individual per-temperature plot
            if per_t_dir is not None:
                fig_t, ax_t = plt.subplots(figsize=(7, 6))
                draw_polygon_frame(ax_t, [m.title() for m in metals])
                _plot_single_phase(ax_t, triang, g_mix, x_cart, y_cart, metals, T, color=color)
                note = "two-phase (shaded)" if found else "single phase"
                ax_t.set_title(f"{label_str}  {int(T)} K  —  {note}", fontsize=11, fontweight="bold")
                fig_t.tight_layout()
                fig_t.savefig(per_t_dir / f"{int(T)}K.png", dpi=150, bbox_inches="tight")
                plt.close(fig_t)

        except Exception as e:
            print(f"    {int(T)} K failed: {e}")

    if not plotted_any:
        plt.close(fig_ov)
        print("  No phase map data plotted.")
        return False

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig_ov.colorbar(sm, ax=ax_ov, pad=0.02, shrink=0.7)
    cbar.set_label("Temperature (K)", fontsize=11)
    note = "shaded = two-phase" if any_two_phase else "single phase across all T"
    ax_ov.set_title(
        f"{label_str}  {int(t_start)}-{int(t_end)} K  (every {int(t_step)} K)\n{note}",
        fontsize=12, fontweight="bold",
    )
    fig_ov.tight_layout()
    fig_ov.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig_ov)
    print(f"  Saved overlay: {output_path.name}")
    return True


def make_phase_gif(tdb, metals, phases, gif_path: Path, data_dir: Path | None = None) -> None:
    """Render one frame per gif_t_step K and save as animated GIF.

    Works for binary, ternary, or higher-order systems.
    """
    import tempfile
    all_comps = list(tdb.elements) + ["VA"]
    temperatures = temperature_values(gif_t_step)
    x_fixed_val = 1.0 - metal_fraction
    label_str = '-'.join(m.title() for m in metals)

    phase_grid = simplex_grid(len(metals), n_grid, metal_fraction)
    x_cart, y_cart = compositions_to_xy(phase_grid, metal_fraction)
    triang = tri.Triangulation(x_cart, y_cart)
    triang = _mask_outside_simplex(triang, x_cart, y_cart, metal_fraction, len(metals))

    frames = []
    clean_frames = []
    tmpdir = tempfile.mkdtemp()
    print(f"  Generating {len(temperatures)} GIF frames ({int(gif_t_step)} K steps)...")

    cmap = cm.plasma
    norm = plt.Normalize(vmin=t_start, vmax=t_end)
    data_dir = data_dir or gif_path.parent / "gmix_data"
    gmix_by_t = load_or_compute_gmix(
        tdb, all_comps, metals, phases, temperatures, phase_grid, x_cart, y_cart, x_fixed_val, data_dir
    )

    for idx, T in enumerate(temperatures):
        T = int(T)
        g_mix = gmix_by_t.get(T)
        color = cmap(norm(T))
        has_data = g_mix is not None and not np.isnan(g_mix).all()

        if not has_data:
            # No G_mix data — still emit a clean frame (empty triangle) for grid sync
            fig_c, ax_c = plt.subplots(figsize=(4, 3.5))
            ax_c.set_facecolor("white"); fig_c.patch.set_facecolor("white")
            _v = generate_polygon_vertices(len(metals))
            _cv = np.vstack([_v, _v[0]])
            ax_c.plot(_cv[:, 0], _cv[:, 1], "k-", lw=1.2)
            _p = 0.15
            ax_c.set_xlim(_v[:, 0].min()-_p, _v[:, 0].max()+_p)
            ax_c.set_ylim(_v[:, 1].min()-_p, _v[:, 1].max()+_p)
            ax_c.set_aspect("equal"); ax_c.axis("off")
            fig_c.tight_layout(pad=0)
            clean_path = os.path.join(tmpdir, f"clean_{idx:04d}.png")
            fig_c.savefig(clean_path, dpi=80)
            plt.close(fig_c)
            clean_frames.append(clean_path)
            continue

        two_phase = find_two_phase_mask(x_cart, y_cart, g_mix, two_phase_threshold)


        # Labeled frame (with title + colorbar)
        fig, axes = plt.subplots(1, 2, figsize=(7, 5.5),
                                 gridspec_kw={"width_ratios": [6, 0.4]})
        ax, cax = axes
        draw_polygon_frame(ax, [m.title() for m in metals])
        if two_phase.any():
            shade_vals = two_phase.astype(float)
            ax.tricontourf(triang, shade_vals, levels=[0.5, 1.5],
                           colors=[color], alpha=0.45)
            ax.tricontour(triang, shade_vals, levels=[0.5],
                          colors=['black'], linewidths=1.5)
            _clip_to_simplex(ax, len(metals))
            note = "two-phase"
        else:
            note = "single phase"
        ax.set_title(f"{label_str}  {int(T)} K  —  {note}", fontsize=11, fontweight="bold")
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label("T (K)", fontsize=9)
        cax.axhline(norm(T), color="white", linewidth=2.5)
        cax.axhline(norm(T), color="black", linewidth=1.0)
        fig.tight_layout()
        frame_path = os.path.join(tmpdir, f"frame_{idx:04d}.png")
        fig.savefig(frame_path, dpi=100)
        plt.close(fig)
        frames.append(frame_path)

        # Clean frame — use same coords as data (generate_polygon_vertices)
        fig_c, ax_c = plt.subplots(figsize=(4, 3.5))
        ax_c.set_facecolor("white")
        fig_c.patch.set_facecolor("white")
        _verts = generate_polygon_vertices(len(metals))
        _closed = np.vstack([_verts, _verts[0]])
        ax_c.plot(_closed[:, 0], _closed[:, 1], "k-", lw=1.2)
        if two_phase.any():
            ax_c.tricontourf(triang, shade_vals, levels=[0.5, 1.5],
                             colors=[color], alpha=0.45)
            ax_c.tricontour(triang, shade_vals, levels=[0.5],
                            colors=['black'], linewidths=1.0)
            _clip_to_simplex(ax_c, len(metals))
        _pad = 0.15
        ax_c.set_xlim(_verts[:, 0].min() - _pad, _verts[:, 0].max() + _pad)
        ax_c.set_ylim(_verts[:, 1].min() - _pad, _verts[:, 1].max() + _pad)
        ax_c.set_aspect("equal")
        ax_c.axis("off")
        fig_c.tight_layout(pad=0)
        clean_path = os.path.join(tmpdir, f"clean_{idx:04d}.png")
        fig_c.savefig(clean_path, dpi=80)
        plt.close(fig_c)
        clean_frames.append(clean_path)

    # Save clean GIF first — always, even if no labeled frames exist
    if clean_frames:
        from PIL import Image as _PILImg
        h0, w0 = imageio.imread(clean_frames[0]).shape[:2]
        clean_imgs = []
        for cf in clean_frames:
            img = imageio.imread(cf)
            if img.shape[:2] != (h0, w0):
                img = np.array(_PILImg.fromarray(img).resize((w0, h0)))
            clean_imgs.append(img)
        clean_gif_path = gif_path.with_name(gif_path.stem + "_clean.gif")
        imageio.mimsave(clean_gif_path, clean_imgs, fps=gif_fps, loop=0)
        print(f"  Saved clean GIF: {clean_gif_path.name} ({len(clean_imgs)} frames, {clean_gif_path.stat().st_size / 1e3:.0f} KB)")

    if not frames:
        print("  No labeled frames — skipping labeled GIF and MP4")
        return

    imgs = [imageio.imread(f) for f in frames]
    imageio.mimsave(gif_path, imgs, fps=gif_fps, loop=0)

    # MP4 via ffmpeg — high CRF (smaller file, minor quality loss)
    mp4_path = gif_path.with_suffix(".mp4")
    frame0 = imgs[0]
    h, w = frame0.shape[:2]
    mp4_cmd = [
        "ffmpeg", "-y", "-f", "rawvideo", "-vcodec", "rawvideo",
        "-s", f"{w}x{h}", "-pix_fmt", "rgb24", "-r", str(gif_fps), "-i", "-",
        "-vcodec", "libx264", "-crf", "28", "-pix_fmt", "yuv420p",
        "-vf", f"scale={w // 2 * 2}:{h // 2 * 2},setpts=PTS/1.5",  # 3x speed
        str(mp4_path),
    ]
    mp4_proc = subprocess.Popen(mp4_cmd, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for img in imgs:
        if img.shape[2] == 4:
            img = img[:, :, :3]
        mp4_proc.stdin.write(img.tobytes())
    mp4_proc.stdin.close()
    mp4_proc.wait()
    print(f"  Saved MP4: {mp4_path.name} ({mp4_path.stat().st_size / 1e3:.0f} KB)")

    for f in frames + clean_frames:
        try: os.remove(f)
        except: pass
    os.rmdir(tmpdir)
    print(f"  Saved GIF: {gif_path.name} ({len(imgs)} frames @ {gif_fps} fps)")


def _plot_binary_phase_diagram(tdb, metals, phase, fixed_sp, temperatures, out_path):
    """Custom T-x binary phase diagram using G_mix convex hull to detect miscibility gaps."""
    from pycalphad import equilibrium, variables as v

    m1, m2 = metals[0].upper(), metals[1].upper()
    comps = [m1, m2] + [el.upper() for el in fixed_sp]
    x_fixed = 1.0 - metal_fraction  # exact
    metal_frac = metal_fraction
    n_pts = 150
    x_site = np.linspace(1e-4, 1 - 1e-4, n_pts)
    x_m2 = np.clip(x_site * metal_frac, 1e-5, metal_frac - 1e-5)

    fig, ax = plt.subplots(figsize=(10, 7))
    two_phase_regions = []  # (T, x_lo, x_hi)

    for T in sorted(temperatures):
        # Fix non-metal elements as scalars, sweep metal_2 as array
        conds = {v.N: 1, v.P: pressure, v.T: float(T),
                 v.X(m2): x_m2}
        for el, frac in fixed_sp.items():  # use exact fraction
            conds[v.X(el.upper())] = x_fixed / len(fixed_sp)
        result = equilibrium(tdb, comps, [phase], conds, output="GM")
        gm = np.squeeze(result.GM.values)
        if np.isnan(gm).all():
            continue
        g0, g1 = gm[0], gm[-1]
        dgmix = gm - ((1 - x_site) * g0 + x_site * g1)

        # Convex hull of (x, dgmix) — points above hull = two-phase region
        from scipy.spatial import ConvexHull as _CH
        try:
            pts = np.column_stack([x_site, dgmix])
            hull = _CH(pts)
            lower_verts = sorted(set(
                v_idx for simplex, eq in zip(hull.simplices, hull.equations)
                if eq[1] < 0 for v_idx in simplex
            ))
            if len(lower_verts) >= 2:
                hull_x = x_site[lower_verts]
                hull_g = dgmix[lower_verts]
                from scipy.interpolate import interp1d
                hull_interp = interp1d(hull_x, hull_g, kind='linear',
                                       bounds_error=False, fill_value='extrapolate')
                gap = dgmix - hull_interp(x_site)
                in_gap = gap > two_phase_threshold
                # find contiguous two-phase spans
                changes = np.diff(in_gap.astype(int))
                starts = np.where(changes == 1)[0] + 1
                ends   = np.where(changes == -1)[0] + 1
                if in_gap[0]:  starts = np.concatenate([[0], starts])
                if in_gap[-1]: ends   = np.concatenate([ends, [n_pts]])
                for s, e in zip(starts, ends):
                    two_phase_regions.append((T, x_site[s], x_site[e-1]))
        except Exception:
            pass

    # Draw two-phase regions as horizontal bars
    cmap_fn = cm.plasma
    norm_T = plt.Normalize(vmin=min(temperatures), vmax=max(temperatures))
    for T, xlo, xhi in two_phase_regions:
        ax.fill_betweenx([T - (temperatures[1]-temperatures[0])/2,
                          T + (temperatures[1]-temperatures[0])/2],
                         xlo, xhi, color=cmap_fn(norm_T(T)), alpha=0.6, linewidth=0)

    # Draw phase boundary lines
    if two_phase_regions:
        T_vals = sorted(set(r[0] for r in two_phase_regions))
        xlo_vals = [r[1] for r in sorted(two_phase_regions)]
        xhi_vals = [r[2] for r in sorted(two_phase_regions)]
        ax.plot([r[1] for r in sorted(two_phase_regions)],
                [r[0] for r in sorted(two_phase_regions)], 'k-', lw=1.5)
        ax.plot([r[2] for r in sorted(two_phase_regions)],
                [r[0] for r in sorted(two_phase_regions)], 'k-', lw=1.5)

    sm = plt.cm.ScalarMappable(cmap=cmap_fn, norm=norm_T)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Temperature (K)")
    ax.set_xlabel(rf"Metal-site fraction $x_{{{m2}}}$", fontsize=14)
    ax.set_ylabel("Temperature (K)", fontsize=14)
    ax.set_title(rf"$({m1},{m2})$ Phase Diagram — Miscibility Gap", fontsize=15)
    ax.set_xlim(0, 1); ax.set_ylim(min(temperatures), max(temperatures))
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_binary_gmix(tdb, metals, phase, fixed_sp, temperatures, out_path):
    """Custom binary G_mix plot — includes VA in components to avoid NaN."""
    from pycalphad import equilibrium, variables as v

    m1, m2 = metals[0].upper(), metals[1].upper()
    comps = [m1, m2] + [el.upper() for el in fixed_sp]
    x_fixed = 1.0 - metal_fraction  # exact
    metal_frac = metal_fraction
    x_site = np.linspace(0.0, 1.0, 201)
    x_m2 = np.clip(x_site * metal_frac, 1e-5, metal_frac - 1e-5)

    cmap_fn = cm.plasma
    norm_T = plt.Normalize(vmin=min(temperatures), vmax=max(temperatures))
    fig, ax = plt.subplots(figsize=(10, 7))

    for T in temperatures:
        conds = {v.N: 1, v.P: pressure, v.T: float(T), v.X(m2): x_m2}
        for el, frac in fixed_sp.items():  # use exact fraction
            conds[v.X(el.upper())] = x_fixed / len(fixed_sp)
        result = equilibrium(tdb, comps, [phase], conds, output="GM")
        gm = np.squeeze(result.GM.values)
        if np.isnan(gm).all():
            continue
        g0, g1 = gm[0], gm[-1]
        dgmix = gm - ((1 - x_site) * g0 + x_site * g1)
        color = cmap_fn(norm_T(T))
        ax.plot(x_site, dgmix, linewidth=2, color=color, label=f"{int(T)} K")

    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    sm = plt.cm.ScalarMappable(cmap=cmap_fn, norm=norm_T)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Temperature (K)")
    ax.set_xlabel(rf"Metal-site fraction $x_{{{m2}}}$", fontsize=14)
    ax.set_ylabel(r"$\Delta G_{\rm mix}$ (J/mol)", fontsize=14)
    ax.set_title(rf"$\Delta G_{{\rm mix}}$ — $({m1},{m2})$ Phase", fontsize=15)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------
# Only take top-level TDB per composition folder (not nested lattice-level TDBs)
tdb_files = sorted(p for p in comps_dir.rglob("*.tdb")
                   if p.parent.parent == comps_dir)
if not tdb_files:
    print(f"No TDB files found under {comps_dir}")
else:
    print(f"Found {len(tdb_files)} TDB files\n")

for tdb_path in tdb_files:
    comp_name = tdb_path.parent.name
    print(f"Processing: {comp_name} / {tdb_path.name}")

    try:
        tdb = Database(str(tdb_path))
    except Exception as e:
        print(f"  Failed to load TDB: {e}")
        continue

    metals = get_variable_elements(tdb, fixed_elements)
    phases = list(tdb.phases.keys())
    print(f"  Elements: {metals}  |  Phases: {phases}")

    comp_out = output_dir / comp_name
    comp_out.mkdir(parents=True, exist_ok=True)

    # Use exact fraction — rounding breaks pycalphad mass balance
    _fixed_sp = {el.title(): (1.0 - metal_fraction) / max(len(fixed_elements), 1)
                 for el in fixed_elements} if fixed_elements else {}
    metals_title = [m.title() for m in metals]
    name_str = "".join(metals_title)
    T_list = [int(T) for T in temperature_values(t_step)]

    if len(metals) == 2:
        ge_path  = comp_out / f"{name_str}_Gibbs_Energy.png"
        gmx_path = comp_out / f"{name_str}_Gibbs_Mixing.png"
        pd_path  = comp_out / f"{name_str}_Phase_Diagram.png"

        if not (skip and ge_path.exists()):
            viz.plot_gibbs_energy(tdb=tdb, metals=metals_title, phase=phases[0],
                fixed_species=_fixed_sp, temperatures=T_list, output_path=ge_path)
            print(f"  Gibbs energy: {ge_path.name}")

        if not (skip and gmx_path.exists()):
            viz.plot_gibbs_mixing(tdb=tdb, metals=metals_title, phase=phases[0],
                fixed_species=_fixed_sp, temperatures=T_list, output_path=gmx_path)
            print(f"  Gibbs mixing: {gmx_path.name}")

        if not (skip and pd_path.exists()):
            viz.plot_binary_phase_diagram(tdb=tdb, metals=metals_title, phases=phases,
                fixed_species=_fixed_sp, temperature_range=(300, 4500, 50),
                output_path=pd_path)
            print(f"  Phase diagram: {pd_path.name}")

    elif len(metals) >= 3:
        # Works for N≥3 metals via affine polygon projection (from notebook)
        n_label = f"{len(metals)}-component"
        per_t_dir = comp_out / "per_temperature"
        per_t_dir.mkdir(exist_ok=True)
        out_path = comp_out / f"{name_str}_phase_map.png"
        data_dir = comp_out / "gmix_data" / tdb_path.stem
        plotted = plot_phase_overlay(tdb, metals, phases, out_path, per_t_dir=per_t_dir, data_dir=data_dir)
        if plotted:
            print(f"  {n_label} phase map: {out_path.name}")
        if make_gif:
            make_phase_gif(tdb, metals, phases,
                           comp_out / f"{name_str}_phase_evolution.gif",
                           data_dir=data_dir)
    else:
        print(f"  Skipping: {len(metals)} variable elements (need ≥2)")

print("\nDone.")
