"""Oxygen thermochemistry: Shomate equations and mu_O / log10(pO2) conversions."""

from __future__ import annotations

import numpy as np

R = 8.31446261815324  # J / (mol K)
EV_KJ_PER_MOL = 96.48533212331002  # eV -> kJ/mol conversion


def _oxygen_shomate_coeffs(T):
    """Return NIST Shomate coefficients for O2 at temperature T (K).

    Args:
        T: Temperature in Kelvin (valid range 100–6000 K).

    Returns:
        Tuple of eight Shomate coefficients (A, B, C, D, E, F, G, Hc).

    Raises:
        ValueError: If T is outside the valid range.
    """
    if 100 <= T < 700:
        return 31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0
    elif 700 <= T < 2000:
        return 30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0
    elif 2000 <= T <= 6000:
        return 20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185, 0.0
    else:
        raise ValueError(f"Shomate O2 valid 100-6000 K; got T={T}")


def _oxygen_delta_g0(T):
    """Standard Gibbs free energy change of O2 formation at T relative to 298.15 K.

    Args:
        T: Temperature in Kelvin.

    Returns:
        Delta G0 in kJ/mol.
    """
    A, B, C, D, E, F, G, Hc = _oxygen_shomate_coeffs(T)
    t = T / 1000.0
    Hr = A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - Hc
    S = A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G
    T0 = 298.15
    A0, B0, C0, D0, E0, _, G0, _ = _oxygen_shomate_coeffs(T0)
    t0 = T0 / 1000.0
    S0 = A0 * np.log(t0) + B0 * t0 + C0 * t0**2 / 2 + D0 * t0**3 / 3 - E0 / (2 * t0**2) + G0
    return Hr - T * S / 1000.0 + T0 * S0 / 1000.0


def mu_o_from_log10_po2(T, log10_po2, mu_o_offset=0.0):
    """Convert log10(pO2) to oxygen chemical potential mu_O (eV per O atom).

    Args:
        T: Temperature in Kelvin.
        log10_po2: Base-10 logarithm of the oxygen partial pressure (in atm).
        mu_o_offset: Optional rigid shift applied to the result (eV).

    Returns:
        Oxygen chemical potential in eV per O atom.
    """
    return (_oxygen_delta_g0(T) + R * T * np.log(10.0) * log10_po2 / 1000.0) / (2.0 * EV_KJ_PER_MOL) + mu_o_offset


def log10_po2_from_mu_o(T, mu_o, mu_o_offset=0.0):
    """Convert oxygen chemical potential mu_O to log10(pO2).

    Args:
        T: Temperature in Kelvin.
        mu_o: Oxygen chemical potential in eV per O atom.
        mu_o_offset: Optional rigid shift that was applied to mu_o (eV).

    Returns:
        Base-10 logarithm of the oxygen partial pressure (in atm).
    """
    return (2.0 * (mu_o - mu_o_offset) * EV_KJ_PER_MOL - _oxygen_delta_g0(T)) / (R * T * np.log(10.0) / 1000.0)
