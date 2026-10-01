"""N-component Muggianu mixing enthalpy and ideal mixing entropy."""

from __future__ import annotations

import numpy as np


def muggianu_energy_nd(y_mat: np.ndarray, h_endpoints: np.ndarray, binary_coeffs: dict, ternary_coeff: float = 0.0) -> np.ndarray:
    """Muggianu mixing enthalpy for an N-metal phase family (eV/formula unit).

    Args:
        y_mat: Composition array of shape (n_pts, n_metals); rows sum to 1.
        h_endpoints: Pure-endpoint energies of shape (n_metals,).
        binary_coeffs: Redlich-Kister coefficients keyed by (i, j) pair with i < j.
        ternary_coeff: Single L^{012} ternary interaction parameter.

    Returns:
        Enthalpy array of shape (n_pts,) in eV per formula unit.
    """
    H = y_mat @ h_endpoints
    for (i, j), L in binary_coeffs.items():
        yi, yj = y_mat[:, i], y_mat[:, j]
        z, z_pow, poly = yi - yj, np.ones(len(y_mat)), np.zeros(len(y_mat))
        for c in L:
            poly += c * z_pow
            z_pow *= z
        H += yi * yj * poly
    if y_mat.shape[1] >= 3 and ternary_coeff != 0.0:
        H += ternary_coeff * y_mat[:, 0] * y_mat[:, 1] * y_mat[:, 2]
    return H


def ideal_mixing_nd(y_mat: np.ndarray) -> np.ndarray:
    """N-component ideal mixing entropy shape: sum_i y_i * ln(y_i).

    Args:
        y_mat: Composition array of shape (n_pts, n_metals); rows sum to 1.

    Returns:
        Entropy shape array of shape (n_pts,). Multiply by k_B * T for energy.
    """
    s = np.zeros(len(y_mat))
    for j in range(y_mat.shape[1]):
        y = y_mat[:, j]
        m = y > 0.0
        s[m] += y[m] * np.log(y[m])
    return s
