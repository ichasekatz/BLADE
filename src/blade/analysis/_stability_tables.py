"""Parse BLADE Comps lattice directories and fit the Redlich-Kister interaction model."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from blade.oxidation.framework._grid import simplex_grid_nd
from blade.oxidation.framework._mixing import ideal_mixing_nd, muggianu_energy_nd


class MissingEndpointError(ValueError):
    """Raised when a pure-metal endpoint composition is absent from the SQS data."""


def build_phase_data(
    lattice_dir: Path,
    metals: list[str],
    phase_element: str,
    phase_stoich: float,
    rk_order: int = 3,
    y_step: float = 0.02,
) -> dict:
    r"""Build the full phase data dictionary for stability LP from a BLADE Comps lattice dir.

    Reads GRACE-relaxed SQS structures from ``sqs_lev=*`` subdirectories, fits the
    Redlich-Kister/Muggianu interaction model over all structures simultaneously, and
    returns a dict suitable for the grand-potential LP with no oxygen chemical potential.

    The same SQS structures serve two roles:
    - **Fixed-phase competitors**: discrete line compounds in the LP (integer stoichs).
    - **RK training set**: endpoints and mixed compositions for fitting H(y).

    Args:
        lattice_dir: Path to ``<Comps>/<System>/<LatticeKey>/`` containing ``sqs_lev=*``
            subdirs, each holding ``energy`` (scalar eV) and ``CONTCAR`` files.
        metals: Ordered list of metal element symbols (e.g. ``["Cr", "Cu"]``).
        phase_element: Non-metal sublattice element (e.g. ``"B"``).
        phase_stoich: Stoichiometry of ``phase_element`` per formula unit (e.g. ``2.0``
            for diborides AlB₂-type). Units: atoms per 1 metal atom.
        rk_order: Redlich-Kister polynomial order (default 3, i.e. L0–L3).
        y_step: Composition simplex grid spacing for the continuous SS phase
            (default 0.02). Coarser grids run faster; finer grids improve dG_hull
            resolution near tie-line vertices.

    Returns:
        Dict with keys: ``phase_ids``, ``phase_kinds``, ``phase_y_nd``,
        ``phase_metal_stoich``, ``phase_H0``, ``phase_mix_shape``, ``A_eq``,
        ``n_fixed``, ``n_metals``, ``metals``, ``phase_element``,
        ``phase_element_stoichiometry``, ``fixed_energy_formula``, ``y_grid``,
        ``h_endpoints``, ``binary_coeffs``, ``ternary_coeff``.

    Raises:
        MissingEndpointError: If any pure-metal endpoint is absent from the SQS data.
        RuntimeError: If no valid SQS structures are found in ``lattice_dir``.
    """
    n = len(metals)
    elements = metals + [phase_element, "O"]
    metal_pats = [re.compile(rf"a_{m}=([0-9.]+)") for m in metals]
    y_cols = [f"y_{m}_metal_site" for m in metals]

    rows = _scan_sqs_dirs(lattice_dir, metals, elements, metal_pats, y_cols, n)
    if not rows:
        raise RuntimeError(f"No valid SQS structures found in {lattice_dir}")
    rows.sort(key=lambda r: [r[y_cols[k]] for k in range(n)])

    h_ep = _endpoint_energies(rows, metals, y_cols, n, lattice_dir)
    binary_coeffs, ternary_coeff = _fit_rk(rows, h_ep, pairs=_pairs(n), n=n, rk_order=rk_order)

    fixed_rows = _build_fixed_rows(rows, lattice_dir, metals, elements, phase_element)
    n_fixed = len(fixed_rows)

    fmetal = np.array([fr["metals"] for fr in fixed_rows])
    fphase_el = np.array([fr["phase_el"] for fr in fixed_rows])
    fo = np.array([fr["O"] for fr in fixed_rows])
    fatoms = fmetal.sum(axis=1) + fphase_el + fo
    fef = np.array([fr["energy_eV_per_atom"] for fr in fixed_rows]) * fatoms

    y_grid = simplex_grid_nd(n, y_step)
    n_pts = len(y_grid)
    phase_h0 = muggianu_energy_nd(y_grid, h_ep, binary_coeffs, ternary_coeff)
    phase_mix = ideal_mixing_nd(y_grid)

    pmetal = np.vstack([fmetal, y_grid])
    pphase_el = np.concatenate([fphase_el, phase_stoich * np.ones(n_pts)])
    a_eq = np.vstack([pmetal.T, pphase_el[np.newaxis, :]])

    fixed_ids = np.array([fr["phase_id"] for fr in fixed_rows])
    ss_ids = _ss_phase_ids(y_grid, metals, n, n_pts)
    phase_ids = np.concatenate([fixed_ids, ss_ids])
    phase_kinds = np.concatenate([np.full(n_fixed, "fixed"), np.full(n_pts, "ss")])
    phase_y_nd = np.vstack([np.full((n_fixed, n), np.nan), y_grid])

    return {
        "phase_ids": phase_ids,
        "phase_kinds": phase_kinds,
        "phase_y_nd": phase_y_nd,
        "phase_metal_stoich": pmetal,
        "phase_H0": phase_h0,
        "phase_mix_shape": phase_mix,
        "fixed_energy_formula": fef,
        "A_eq": a_eq,
        "n_fixed": n_fixed,
        "n_metals": n,
        "metals": metals,
        "phase_element": phase_element,
        "phase_element_stoichiometry": phase_stoich,
        "y_grid": y_grid,
        "h_endpoints": h_ep,
        "binary_coeffs": binary_coeffs,
        "ternary_coeff": ternary_coeff,
    }


def detect_metals(lattice_dir: Path) -> list[str]:
    """Infer ordered metal list from ``a_<Metal>=x`` patterns in SQS dir names.

    Args:
        lattice_dir: Path to a BLADE Comps lattice directory.

    Returns:
        Ordered list of metal element symbols (by first appearance in sorted dir names).
    """
    pat = re.compile(r"a_([A-Z][a-z]?)=")
    seen: list[str] = []
    for d in sorted(lattice_dir.iterdir()):
        if not d.is_dir() or not d.name.startswith("sqs_lev="):
            continue
        for sym in pat.findall(d.name):
            if sym not in seen:
                seen.append(sym)
    return seen


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _pairs(n: int) -> list[tuple[int, int]]:
    """Return all ordered (i, j) pairs with i < j for n elements."""
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def _read_scalar(path: Path) -> float:
    """Read the first float from a plain-text file."""
    return float(path.read_text().strip().split()[0])


def _read_contcar_counts(path: Path, elements: list[str]) -> dict[str, int]:
    """Parse element counts from VASP CONTCAR/POSCAR lines 6–7.

    Args:
        path: Path to CONTCAR file.
        elements: Element symbols to collect; others are ignored.

    Returns:
        Dict mapping element symbol → count (0 for elements not present).
    """
    lines = path.read_text().splitlines()
    if len(lines) < 7:
        raise ValueError(f"CONTCAR too short: {path}")
    syms = lines[5].split()
    counts = list(map(int, lines[6].split()))
    result: dict[str, int] = dict.fromkeys(elements, 0)
    for s, c in zip(syms, counts, strict=False):
        if s in result:
            result[s] += c
    return result


def _gcd(a: int, b: int) -> int:
    """Greatest common divisor via Euclidean algorithm."""
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _reduce_counts(counts: dict[str, int]) -> dict[str, int]:
    """Divide all atom counts by their GCD to obtain the reduced formula.

    Args:
        counts: Dict of element → atom count.

    Returns:
        Counts divided by GCD (or unchanged if GCD ≤ 1).
    """
    vals = [v for v in counts.values() if v > 0]
    if not vals:
        return counts
    g = vals[0]
    for v in vals[1:]:
        g = _gcd(g, v)
    return {k: v // g if g > 1 else v for k, v in counts.items()}


def _scan_sqs_dirs(
    lattice_dir: Path,
    metals: list[str],
    elements: list[str],
    metal_pats: list[re.Pattern],
    y_cols: list[str],
    n: int,
) -> list[dict]:
    """Walk sqs_lev= subdirs and collect composition + energy rows.

    Args:
        lattice_dir: Lattice directory containing sqs_lev= subdirs.
        metals: Ordered metal symbols.
        elements: All elements including phase_element and O.
        metal_pats: Compiled regex patterns for a_<Metal>=x per metal.
        y_cols: Column names for metal-site fractions.
        n: Number of metals.

    Returns:
        List of row dicts with phase_id, y_<Metal>_metal_site, energy_eV_per_formula.
    """
    rows: list[dict] = []
    for d in lattice_dir.iterdir():
        if not d.is_dir() or not d.name.startswith("sqs_lev="):
            continue
        ep, cp = d / "energy", d / "CONTCAR"
        if not ep.exists() or not cp.exists():
            continue
        fracs = []
        for pat in metal_pats:
            m = pat.search(d.name)
            fracs.append(float(m.group(1)) if m else 0.0)
        total_metal = sum(fracs)
        if total_metal <= 0:
            continue
        y = [c / total_metal for c in fracs]
        counts = _read_contcar_counts(cp, elements)
        n_fu = sum(counts[m] for m in metals)
        if n_fu <= 0:
            continue
        e_total = _read_scalar(ep)
        rows.append({"phase_id": d.name, **{y_cols[k]: y[k] for k in range(n)}, "energy_eV_per_formula": e_total / n_fu})
    return rows


def _endpoint_energies(rows: list[dict], metals: list[str], y_cols: list[str], n: int, lattice_dir: Path) -> np.ndarray:
    """Extract pure-endpoint energies from the row list.

    Args:
        rows: Sorted SQS rows from :func:`_scan_sqs_dirs`.
        metals: Ordered metal symbols.
        y_cols: Column names for metal-site fractions.
        n: Number of metals.
        lattice_dir: Used in the error message.

    Returns:
        Array of endpoint energies in eV/formula, shape (n,).

    Raises:
        MissingEndpointError: If any endpoint is absent.
    """
    h_ep = np.zeros(n)
    for i, m in enumerate(metals):
        eye = [1.0 if k == i else 0.0 for k in range(n)]
        eps = [r["energy_eV_per_formula"] for r in rows if all(abs(r[y_cols[k]] - eye[k]) < 1e-10 for k in range(n))]
        if not eps:
            raise MissingEndpointError(f"Missing pure endpoint y_{m}=1 in {lattice_dir}")
        h_ep[i] = eps[0]
    return h_ep


def _fit_rk(
    rows: list[dict],
    h_ep: np.ndarray,
    pairs: list[tuple[int, int]],
    n: int,
    rk_order: int,
) -> tuple[dict[tuple[int, int], np.ndarray], float]:
    """Global simultaneous least-squares fit of Redlich-Kister/Muggianu parameters.

    Fits all binary L parameters and an optional ternary L parameter in a single
    lstsq over all non-endpoint SQS structures. Endpoints are excluded because
    their excess enthalpy is zero by construction.

    References:
        Muggianu, Y. M. et al. (1975). J. Chim. Phys., 72, 83–88.

    Args:
        rows: SQS rows; pure-endpoint rows (any y_k ≥ 1) are automatically skipped.
        h_ep: Endpoint energies in eV/formula, shape (n,).
        pairs: List of (i, j) index pairs with i < j.
        n: Number of metals.
        rk_order: Polynomial order (number of L parameters per pair = rk_order + 1).

    Returns:
        Tuple of (binary_coeffs, ternary_coeff) where binary_coeffs maps (i, j)
        to an array of L0, L1, ... coefficients in eV/formula.
    """
    # rebuild actual column names from rows (metal order comes from the row keys)
    actual_y_cols = [k for k in rows[0] if k.startswith("y_") and "metal_site" in k]
    actual_y_cols.sort()

    n_params = len(pairs) * (rk_order + 1) + (1 if n == 3 else 0)
    basis_rows: list[list[float]] = []
    targets: list[float] = []

    for r in rows:
        y_vec = np.array([r[c] for c in actual_y_cols])
        if np.any(y_vec >= 1.0 - 1e-10):
            continue
        h_exc = r["energy_eV_per_formula"] - float(np.dot(y_vec, h_ep))
        basis: list[float] = []
        for i, j in pairs:
            yi, yj = y_vec[i], y_vec[j]
            z, zpow = yi - yj, 1.0
            for _ in range(rk_order + 1):
                basis.append(yi * yj * zpow)
                zpow *= z
        if n == 3:
            basis.append(y_vec[0] * y_vec[1] * y_vec[2])
        basis_rows.append(basis)
        targets.append(h_exc)

    if basis_rows:
        coeffs_all, _, _, _ = np.linalg.lstsq(np.array(basis_rows), np.array(targets), rcond=None)
    else:
        coeffs_all = np.zeros(n_params)

    binary_coeffs: dict[tuple[int, int], np.ndarray] = {}
    idx = 0
    for pair in pairs:
        binary_coeffs[pair] = coeffs_all[idx : idx + rk_order + 1]
        idx += rk_order + 1
    ternary_coeff = float(coeffs_all[idx]) if n == 3 and idx < len(coeffs_all) else 0.0
    return binary_coeffs, ternary_coeff


def _build_fixed_rows(
    rows: list[dict],
    lattice_dir: Path,
    metals: list[str],
    elements: list[str],
    phase_element: str,
) -> list[dict]:
    """Build fixed-phase competitor rows with reduced-formula stoichs and energy/atom.

    Args:
        rows: SQS rows from :func:`_scan_sqs_dirs`.
        lattice_dir: Lattice directory (used to locate CONTCAR/energy files).
        metals: Ordered metal symbols.
        elements: All elements including phase_element and O.
        phase_element: Non-metal sublattice element.

    Returns:
        List of dicts with keys: phase_id, metals (ndarray), phase_el, O,
        energy_eV_per_atom.
    """
    fixed: list[dict] = []
    for r in rows:
        d = lattice_dir / r["phase_id"]
        cp, ep = d / "CONTCAR", d / "energy"
        counts = _read_contcar_counts(cp, elements)
        n_atoms = sum(counts.values())
        if n_atoms <= 0:
            continue
        reduced = _reduce_counts(counts)
        e_total = _read_scalar(ep)
        fixed.append(
            {
                "phase_id": r["phase_id"],
                "metals": np.array([float(reduced[m]) for m in metals]),
                "phase_el": float(reduced[phase_element]),
                "O": float(reduced.get("O", 0)),
                "energy_eV_per_atom": e_total / n_atoms,
            }
        )
    return fixed


def _ss_phase_ids(y_grid: np.ndarray, metals: list[str], n: int, n_pts: int) -> np.ndarray:
    """Generate string phase IDs for the continuous SS phase grid.

    Args:
        y_grid: Composition grid, shape (n_pts, n).
        metals: Ordered metal symbols.
        n: Number of metals.
        n_pts: Number of grid points.

    Returns:
        Array of ID strings, shape (n_pts,).
    """
    if n == 2:
        return np.array([f"SS_y={y[0]:.4f}" for y in y_grid])
    return np.array(["SS_" + "_".join(f"{metals[k]}={y_grid[r, k]:.4f}" for k in range(n - 1)) for r in range(n_pts)])
