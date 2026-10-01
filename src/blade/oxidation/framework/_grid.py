"""Simplex composition grid and plotting edge helpers."""

from __future__ import annotations

import numpy as np


def simplex_grid_nd(n_metals: int, step: float = 0.01) -> np.ndarray:
    """Composition grid on the (n_metals-1)-simplex.

    Args:
        n_metals: Number of metal components (must be >= 2).
        step: Grid spacing along each simplex edge.

    Returns:
        Array of shape (n_pts, n_metals) where each row sums to 1.

    Raises:
        ValueError: If n_metals < 2 or step is invalid.
    """
    if n_metals < 2:
        raise ValueError(f"simplex_grid_nd needs at least 2 metals; got {n_metals}")
    if n_metals == 2:
        y0 = np.arange(0.0, 1.0 + step / 2.0, step)
        if abs(y0[-1] - 1.0) > 1e-12:
            y0 = np.append(y0, 1.0)
        return np.column_stack([y0, 1.0 - y0])
    elif n_metals == 3:
        pts = []
        n = int(round(1.0 / step))
        for i in range(n + 1):
            y0 = i / n
            for j in range(n - i + 1):
                y1 = j / n
                y2 = 1.0 - y0 - y1
                if y2 < -1e-10:
                    break
                pts.append([y0, y1, max(0.0, y2)])
        return np.array(pts)

    n = int(round(1.0 / step))
    if n <= 0:
        raise ValueError(f"Invalid simplex step: {step}")

    pts = []

    def _fill(prefix, remaining_units, remaining_dims):
        if remaining_dims == 1:
            pts.append(prefix + [remaining_units / n])
            return
        for units in range(remaining_units + 1):
            _fill(prefix + [units / n], remaining_units - units, remaining_dims - 1)

    _fill([], n, n_metals)
    return np.array(pts)


def grid_edges(values):
    """Compute cell-edge positions for a pcolormesh grid.

    Args:
        values: 1-D array of cell-centre coordinates.

    Returns:
        Array of length ``len(values) + 1`` containing the cell edges.
    """
    v = np.asarray(values, dtype=float)
    if len(v) == 1:
        return np.array([v[0] - 0.5, v[0] + 0.5])
    mids = 0.5 * (v[:-1] + v[1:])
    return np.concatenate([[v[0] - 0.5 * (v[1] - v[0])], mids, [v[-1] + 0.5 * (v[-1] - v[-2])]])
