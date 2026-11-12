"""Bartel G(T) temperature correction for fixed reference oxide phases.

Applies the SISSO-learned descriptor from Bartel et al. (2018) [1]_ to add
finite-temperature Gibbs corrections to 0 K DFT formation enthalpies.

The descriptor for the vibrational / configurational correction is:

    G_delta(V, m, T) = (-2.48e-4 * ln(V) - 8.94e-5 * m/V) * T
                       + 0.181 * ln(T) - 0.882          [eV/atom]

where V is the volume per atom [Å³/atom] and m is a reduced-mass descriptor
computed from atomic masses and stoichiometry.

The corrected Gibbs energy of a compound with n atoms per formula unit is:

    G(T) = H_f(0 K) + G_delta(V, m, T)   [eV/atom]

so the per-formula-unit energy shift is simply G_delta * n_atoms.

Requires a ``vol_per_atom`` column in the fixed-phase table CSV.

.. [1] Bartel et al., Nature Communications 9, 4168 (2018).
       https://doi.org/10.1038/s41467-018-06682-4
"""

from __future__ import annotations

import math
from itertools import combinations

import numpy as np

# ---------------------------------------------------------------------------
# IUPAC atomic masses [amu] — from Bartel et al. supplementary data
# ---------------------------------------------------------------------------
_MASSES: dict[str, float] = {
    "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811,
    "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797,
    "Na": 22.98976928, "Mg": 24.305, "Al": 26.9815386, "Si": 28.0855,
    "P": 30.973762, "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.0983,
    "Ca": 40.078, "Sc": 44.955912, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
    "Mn": 54.938045, "Fe": 55.845, "Co": 58.933195, "Ni": 58.6934,
    "Cu": 63.546, "Zn": 65.409, "Ga": 69.723, "Ge": 72.64, "As": 74.9216,
    "Se": 78.96, "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62,
    "Y": 88.90585, "Zr": 91.224, "Nb": 92.90638, "Mo": 95.94, "Tc": 98.0,
    "Ru": 101.07, "Rh": 102.9055, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411,
    "In": 114.818, "Sn": 118.71, "Sb": 121.76, "Te": 127.6, "I": 126.90447,
    "Xe": 131.293, "Cs": 132.9054519, "Ba": 137.327, "La": 138.90547,
    "Ce": 140.116, "Pr": 140.90765, "Nd": 144.242, "Pm": 145.0, "Sm": 150.36,
    "Eu": 151.964, "Gd": 157.25, "Tb": 158.92535, "Dy": 162.5,
    "Ho": 164.93032, "Er": 167.259, "Tm": 168.93421, "Yb": 173.04,
    "Lu": 174.967, "Hf": 178.49, "Ta": 180.94788, "W": 183.84, "Re": 186.207,
    "Os": 190.23, "Ir": 192.217, "Pt": 195.084, "Au": 196.966569,
    "Hg": 200.59, "Tl": 204.3833, "Pb": 207.2, "Bi": 208.9804,
    "Th": 232.03806, "Pa": 231.03588, "U": 238.02891,
}


def _reduced_mass(counts: dict[str, int | float]) -> float:
    """Compute the reduced-mass descriptor m for the Bartel SISSO model.

    Args:
        counts: Mapping of element symbol to stoichiometric count per formula unit.

    Returns:
        Reduced mass m [amu] for the SISSO descriptor, or NaN if any element
        mass is unknown.
    """
    names = list(counts.keys())
    nums = [float(counts[el]) for el in names]
    n_els = len(names)
    n_atoms = sum(nums)
    denom = (n_els - 1) * n_atoms
    if denom <= 0:
        return float("nan")

    masses = []
    for el in names:
        m = _MASSES.get(el, float("nan"))
        if math.isnan(m):
            return float("nan")
        masses.append(m)

    pair_red_sum = 0.0
    for i, j in combinations(range(n_els), 2):
        pair_coeff = nums[i] + nums[j]
        pair_prod = masses[i] * masses[j]
        pair_sum = masses[i] + masses[j]
        pair_red_sum += pair_coeff * pair_prod / pair_sum

    return pair_red_sum / denom


def _gd_sisso(vol_per_atom: float, m: float, T: float) -> float:
    """Evaluate the Bartel SISSO G_delta descriptor.

    Args:
        vol_per_atom: Volume per atom [Å³/atom].
        m: Reduced-mass descriptor [amu].
        T: Temperature [K].

    Returns:
        G_delta(V, m, T) [eV/atom].
    """
    V = vol_per_atom
    return (-2.48e-4 * math.log(V) - 8.94e-5 * m / V) * T + 0.181 * math.log(T) - 0.882


def bartel_fixed_energy_corrections(pd_data: dict, T: float) -> np.ndarray:
    """Return per-formula-unit Bartel G(T) corrections for the fixed oxide phases.

    Requires ``pd_data`` to contain:

    - ``"fixed_vol_per_atom"`` — 1-D array [Å³/atom], shape (n_fixed,).
    - ``"fixed_formula_counts"`` — list of ``{element: count}`` dicts, length n_fixed.

    If a phase has an unknown volume, missing masses, or is a pure element
    (n_elements == 1), the correction for that phase is set to zero.

    Args:
        pd_data: Phase-data dict produced by :func:`_phase_data._load_nd`.
        T: Temperature [K].

    Returns:
        Array of shape (n_fixed,) with energy corrections [eV per formula unit].

    Raises:
        KeyError: If ``vol_per_atom`` was not loaded (Bartel data not available).
    """
    n_fixed = pd_data["n_fixed"]
    vol = pd_data.get("fixed_vol_per_atom")
    formula_counts = pd_data.get("fixed_formula_counts")

    if vol is None or formula_counts is None:
        raise KeyError(
            "bartel_corrections=True requires a 'vol_per_atom' column in the "
            "fixed-phase table CSV. Regenerate tables with vol_per_atom populated "
            "or set bartel_corrections=False."
        )

    corrections = np.zeros(n_fixed, dtype=float)

    for i in range(n_fixed):
        counts = formula_counts[i]
        n_atoms = sum(counts.values())

        if len(counts) <= 1 or n_atoms == 0:
            continue

        v = float(vol[i])
        if not math.isfinite(v) or v <= 0.0:
            continue

        m = _reduced_mass(counts)
        if not math.isfinite(m):
            continue

        corrections[i] = _gd_sisso(v, m, T) * n_atoms

    return corrections
