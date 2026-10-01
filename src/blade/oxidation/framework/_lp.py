"""Grand-potential linear-programming solver for phase equilibria."""

from __future__ import annotations

import numpy as np
from scipy.optimize import linprog


def solve_grand_lp_batch(A_eq, b_eq, grand_stack):
    """Solve an objective stack without per-batch thread-pool overhead.

    Each row still uses the exact same cold LP solve as :func:`solve_grand_lp`,
    preserving deterministic phase amounts for degenerate optima.

    Args:
        A_eq: Equality constraint matrix shared across all objectives.
        b_eq: Equality constraint right-hand side vector.
        grand_stack: Iterable of grand-potential objective vectors.

    Returns:
        List of (amounts, obj, ok) tuples, one per objective in grand_stack.
    """
    return [solve_grand_lp(A_eq, b_eq, grand) for grand in grand_stack]


def solve_grand_lp(A_eq, b_eq, grand):
    """Minimize grand @ n s.t. A_eq @ n = b_eq, n >= 0.

    Args:
        A_eq: Equality constraint matrix of shape (n_constraints, n_phases).
        b_eq: Right-hand side vector of shape (n_constraints,).
        grand: Grand-potential objective vector of shape (n_phases,).

    Returns:
        Tuple of (amounts, obj, ok) where amounts is the phase-amount vector,
        obj is the minimized objective value, and ok is True on success.
    """
    bounds = [(0.0, None)] * len(grand)
    opts = {"disp": False}
    # revised simplex is ~1.5× faster for small dense LPs; fall back to HiGHS if removed
    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = linprog(grand, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="revised simplex", options=opts)
    except Exception:
        result = linprog(grand, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs", options=opts)
    if result.status == 0:
        x = np.maximum(result.x, 0.0)
        x[x < 1e-10] = 0.0
        return x, float(result.fun), True
    return np.full(len(grand), np.nan), np.nan, False
