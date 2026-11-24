"""Multi-sublattice utilities for BLADE post-TDB stages.

Provides sublattice-aware dir-name parsing, composition grids, and formula-unit
helpers used by the stability and phase-diagram stages for two-sublattice phases
like MAX (M₂AX: sublattice a = M-site mult 2, sublattice b = A-site mult 1,
sublattice c = C fixed mult 1).

Single-sublattice paths (borides, alloys) are unaffected: callers that do not
pass ``sublattice_metals`` / ``sublattice_mults`` continue to use the existing
code.
"""

from __future__ import annotations

import re
from itertools import product as iterproduct
from pathlib import Path

import numpy as np

# Pattern: matches "<letter>_<Element>=<value>" tokens in an sqs_lev dir name.
_SUBL_PAT = re.compile(r"([a-z])_([A-Z][a-z]?)=([0-9]*\.?[0-9]+)")


def parse_subl_fracs(dirname: str) -> dict[str, dict[str, float]]:
    """Parse per-sublattice fractional compositions from an ``sqs_lev=*`` dir name.

    Args:
        dirname: Directory basename, e.g.
            ``"sqs_lev=1_a_Ti=0.5,a_Nb=0.5,b_Al=1.0,c_C=1.0"``.

    Returns:
        Nested dict ``{sublattice_letter: {element: fraction}}``, e.g.
        ``{"a": {"Ti": 0.5, "Nb": 0.5}, "b": {"Al": 1.0}, "c": {"C": 1.0}}``.
        Empty dict when no ``<letter>_<El>=<x>`` tokens are found.
    """
    result: dict[str, dict[str, float]] = {}
    for letter, element, value in _SUBL_PAT.findall(dirname):
        result.setdefault(letter, {})[element] = float(value)
    return result


def product_simplex_grid(
    subl_metals: list[list[str]],
    subl_mults: list[int],
    y_step: float = 0.05,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Build a composition grid as the Cartesian product of per-sublattice simplices.

    For a two-sublattice MAX phase with ``subl_metals = [["Ti","Nb"], ["Al","Si"]]``
    and ``subl_mults = [2, 1]``:

    - Sublattice a (M-site, mult 2): simplex in Ti–Nb space
    - Sublattice b (A-site, mult 1): simplex in Al–Si space
    - Overall composition: ``x_Mi = 2/(2+1) * y_a_i``, ``x_Aj = 1/(2+1) * y_b_j``

    The resulting ``y_overall`` can be used directly in the LP element-balance
    constraint (``b_eq = [y_overall..., phase_stoich]``).

    Args:
        subl_metals: Ordered metal symbols per sublattice,
            e.g. ``[["Ti","Nb"], ["Al","Si"]]``.
        subl_mults: Multiplicity of each sublattice,
            e.g. ``[2, 1]`` for MAX M-site (mult 2) and A-site (mult 1).
        y_step: Simplex grid spacing (default 0.05).  Coarser than the
            single-sublattice default (0.02) because the product grows as n_a × n_b.

    Returns:
        ``(y_overall, subl_grids)`` where:

        - ``y_overall``: ``(n_points, sum(len(m) for m in subl_metals))`` float
          array of overall metal-fraction vectors (each row sums to 1).
        - ``subl_grids``: list of per-sublattice simplex arrays, each
          ``(n_subl_pts, n_metals_on_subl)``.

    Raises:
        ValueError: If ``subl_metals`` and ``subl_mults`` have different lengths,
            or any sublattice is empty.
    """
    if len(subl_metals) != len(subl_mults):
        raise ValueError(
            f"subl_metals and subl_mults must have the same length, "
            f"got {len(subl_metals)} vs {len(subl_mults)}"
        )
    for i, metals in enumerate(subl_metals):
        if not metals:
            raise ValueError(f"sublattice {i} has no metals")

    from blade.oxidation.framework._grid import simplex_grid_nd

    total_mult = sum(subl_mults)
    # simplex_grid_nd requires n >= 2; a single-element sublattice is a trivial 1×1 grid.
    subl_grids = [
        simplex_grid_nd(len(metals), y_step) if len(metals) >= 2 else np.ones((1, 1), dtype=float)
        for metals in subl_metals
    ]

    # Build product: all combinations of one point from each sublattice simplex.
    # Weight each sublattice by its multiplicity share.
    rows = []
    for pts in iterproduct(*subl_grids):
        row = np.concatenate([pts[i] * subl_mults[i] / total_mult for i in range(len(subl_metals))])
        rows.append(row)
    y_overall = np.array(rows, dtype=float)
    return y_overall, subl_grids


def phase_stoich_from_mults(subl_mults: list[int], fixed_mult: int) -> float:
    """Return the fixed-element stoichiometry per total metal atom.

    For MAX with mult a=2, b=1, c=1 (C fixed):
    ``phase_stoich = 1 / (2 + 1) = 0.333...``

    This is the ``phase_stoich`` argument expected by :func:`build_phase_data` and
    :class:`~blade.analysis.blade_stability.BladeStability`.

    Args:
        subl_mults: Multiplicities of variable (metal-bearing) sublattices.
        fixed_mult: Multiplicity of the fixed (non-metal) sublattice.

    Returns:
        ``fixed_mult / sum(subl_mults)``
    """
    total_metal = sum(subl_mults)
    if total_metal <= 0:
        raise ValueError("subl_mults must contain at least one positive multiplicity")
    return fixed_mult / total_metal


def detect_sublattice_metals(
    lattice_dir: Path,
    phase_element: str = "B",
) -> tuple[dict[str, list[str]], bool]:
    """Scan sqs_lev=* dirs and return per-sublattice metal lists.

    Args:
        lattice_dir: Blade structure directory containing ``sqs_lev=*`` subdirs.
        phase_element: Fixed element to exclude (default ``"B"``).

    Returns:
        ``(subl_metals, is_multi_sublattice)`` where ``subl_metals`` maps sublattice
        letter to ordered list of metal symbols, and ``is_multi_sublattice`` is True
        when more than one letter has variable metals.
    """
    _fixed = {phase_element} if phase_element else set()
    subl: dict[str, list[str]] = {}
    for d in sorted(Path(lattice_dir).iterdir()):
        if not d.is_dir() or not d.name.startswith("sqs_lev="):
            continue
        fracs = parse_subl_fracs(d.name)
        for letter, el_map in fracs.items():
            for sym in el_map:
                if sym in _fixed:
                    continue
                subl.setdefault(letter, [])
                if sym not in subl[letter]:
                    subl[letter].append(sym)
    is_multi = len(subl) > 1
    return subl, is_multi


def cross_endmember_fracs(
    subl_metals: list[list[str]],
    subl_mults: list[int],
) -> list[dict[str, float]]:
    """Return one overall-composition dict per cross-product endmember.

    For MAX with M=[Ti,Nb] and A=[Al,Si] the four endmembers are:
    Ti₂AlC, Ti₂SiC, Nb₂AlC, Nb₂SiC.

    Args:
        subl_metals: Per-sublattice metal lists.
        subl_mults: Per-sublattice multiplicities.

    Returns:
        List of ``{element: fraction}`` dicts (one per cross-product endpoint).
        Each dict has exactly one entry per sublattice (the pure-metal endpoint
        on that sublattice) with fraction = mult / sum(mults).
    """
    total_mult = sum(subl_mults)
    endpoints = []
    # One pure-element choice per sublattice, Cartesian product.
    for chosen in iterproduct(*subl_metals):
        d: dict[str, float] = {}
        for el, mult in zip(chosen, subl_mults, strict=True):
            d[el] = mult / total_mult
        endpoints.append(d)
    return endpoints
