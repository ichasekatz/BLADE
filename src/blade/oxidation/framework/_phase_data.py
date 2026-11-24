"""Phase data loading for binary and N-metal oxidation systems.

The 2-sublattice (MAX-phase) path is auto-detected from the interaction
coefficient CSV: when ``pair == "endmember"`` rows are present, :func:`_load_nd`
dispatches to :func:`_load_nd_2subl` which uses a Compound Energy Formalism (CEF)
product-simplex grid instead of the single-sublattice Redlich-Kister model.
Single-sublattice callers (borides, alloys) are unaffected.
"""

from __future__ import annotations

from itertools import product as _iterproduct

import numpy as np
import pandas as pd

from ._grid import simplex_grid_nd
from ._mixing import ideal_mixing_nd, muggianu_energy_nd


def load_phase_data(
    phase_table_file,
    phase_points_file,
    rk_coeff_file,
    metals_or_M1,
    M2=None,
    y_column=None,
    phase_label=None,
    y_step=0.01,
    phase_element=None,
    phase_element_stoichiometry=0.0,
    exclude_oxide_phases: bool = False,
):
    """Load phase data for binary or N-metal system.

    The flexible phase's non-metal element and stoichiometry are explicit inputs.

    Args:
        phase_table_file: Path to the fixed-phase CSV table.
        phase_points_file: Path to the phase grid points CSV.
        rk_coeff_file: Path to the Redlich-Kister coefficients CSV.
        metals_or_M1: List of metal symbols (N-metal) or first metal symbol (binary).
        M2: Second metal symbol when using the binary (M1, M2) calling convention.
        y_column: Unused; retained for backward compatibility.
        phase_label: Override the default phase-family label string.
        y_step: Simplex grid spacing (default 0.01).
        phase_element: Non-oxygen non-metal sublattice element (e.g. ``"B"``).
        phase_element_stoichiometry: Stoichiometry of phase_element per formula unit.

    Returns:
        Dict with keys: phase_ids, phase_kinds, phase_y, phase_y_nd,
        phase_metal_stoich, phase_element_stoich, phase_O, fixed_energy_formula,
        phase_H0, phase_mix_shape, y_grid, h_endpoints, binary_coeffs,
        ternary_coeff, A_eq, n_fixed, n_metals, metals, phase_label,
        phase_element, phase_element_stoichiometry, phase_suffix.
    """
    if isinstance(metals_or_M1, (list, tuple)):
        metals = list(metals_or_M1)
        if phase_label is None:
            phase_label = "".join(metals)
    else:
        metals = [metals_or_M1, M2]
        if phase_label is None:
            phase_label = f"{metals_or_M1}{M2}"
    return _load_nd(
        phase_table_file,
        phase_points_file,
        rk_coeff_file,
        metals,
        phase_label,
        y_step,
        phase_element,
        phase_element_stoichiometry,
        exclude_oxide_phases=exclude_oxide_phases,
    )


def _load_nd(
    phase_table_file,
    phase_points_file,
    rk_coeff_file,
    metals,
    phase_label,
    y_step,
    phase_element,
    phase_element_stoichiometry,
    exclude_oxide_phases: bool = False,
):
    """Internal N-metal phase data loader.

    Args:
        phase_table_file: Path to the fixed-phase CSV table.
        phase_points_file: Path to the phase grid points CSV.
        rk_coeff_file: Path to the Redlich-Kister coefficients CSV.
        metals: List of metal element symbols.
        phase_label: Phase-family label string.
        y_step: Simplex grid spacing.
        phase_element: Non-oxygen non-metal sublattice element or None.
        phase_element_stoichiometry: Stoichiometry of phase_element per formula unit.

    Returns:
        Dict containing all phase data arrays and metadata required by the LP solver.

    Raises:
        ValueError: If a pure-endpoint row is missing from phase_points_file.
    """
    n = len(metals)
    pt = pd.read_csv(phase_table_file)
    if exclude_oxide_phases and "O" in pt.columns:
        pt = pt[pt["O"].astype(float) == 0].reset_index(drop=True)
    bt = pd.read_csv(phase_points_file)
    rk = pd.read_csv(rk_coeff_file)

    # 2-sublattice CEF path: auto-detected from "endmember" pair rows.
    if "pair" in rk.columns and (rk["pair"] == "endmember").any():
        return _load_nd_2subl(
            pt, rk, metals, phase_label, y_step,
            phase_element, phase_element_stoichiometry,
        )

    # Fixed phases
    fids = pt["phase_id"].astype(str).values
    fmetal = np.column_stack([pt[m].astype(float).values for m in metals])
    phase_element_stoichiometry = float(phase_element_stoichiometry) if phase_element else 0.0
    fphase_element = (
        pt[phase_element].astype(float).values if phase_element and phase_element in pt.columns else np.zeros(len(pt))
    )
    fO = pt["O"].astype(float).values
    fatoms = fmetal.sum(axis=1) + fphase_element + fO
    fef = pt["energy_eV_per_atom"].astype(float).values * fatoms
    n_fixed = len(fids)

    # Optional Bartel data: vol_per_atom column and formula counts
    fixed_vol_per_atom = None
    if "vol_per_atom" in pt.columns:
        fixed_vol_per_atom = pt["vol_per_atom"].astype(float).values

    fixed_formula_counts: list[dict[str, float]] = []
    all_els = metals + ([phase_element] if phase_element else []) + ["O"]
    all_counts = [fmetal[:, j] for j in range(len(metals))]
    all_counts += [fphase_element] if phase_element else []
    all_counts += [fO]
    for i in range(n_fixed):
        counts = {el: float(all_counts[k][i]) for k, el in enumerate(all_els) if float(all_counts[k][i]) > 0}
        fixed_formula_counts.append(counts)

    # Simplex grid
    y_grid = simplex_grid_nd(n, y_step)  # (n_pts, n)
    n_pts = len(y_grid)

    # Endpoint energies from the phase grid CSV
    y_cols = [f"y_{m}_metal_site" for m in metals]
    if all(c in bt.columns for c in y_cols):
        by = bt[y_cols].astype(float).values
    else:
        # Legacy binary CSV: only y_M1_metal_site present; derive y_M2 = 1 - y_M1
        old = [c for c in bt.columns if c.startswith("y_") and "metal" in c]
        y0 = bt[old[0]].astype(float).values
        by = np.column_stack([y0, 1.0 - y0])
    be = bt["energy_eV_per_formula"].astype(float).values

    h_ep = np.zeros(n)
    for i, m in enumerate(metals):
        eye = np.zeros(n)
        eye[i] = 1.0
        mask = np.all(np.abs(by - eye) < 1e-10, axis=1)
        if not np.any(mask):
            raise ValueError(f"Missing pure endpoint y_{m}=1 in {phase_points_file}")
        h_ep[i] = float(be[mask][0])

    # Interaction coefficients
    bc: dict = {}
    tc = 0.0
    if "pair" in rk.columns:
        for pstr, grp in rk.groupby("pair"):
            nd = pstr.count("-")
            if nd == 1:
                a, b = pstr.split("-")
                if a in metals and b in metals:
                    i, j = metals.index(a), metals.index(b)
                    if i > j:
                        i, j = j, i
                    bc[(i, j)] = grp.sort_values("term")["value_eV_per_formula"].values
            elif nd >= 2:
                tc = float(grp["value_eV_per_formula"].iloc[0])
    else:
        bc[(0, 1)] = rk["value_eV_per_formula"].values

    phase_H0 = muggianu_energy_nd(y_grid, h_ep, bc, tc)
    phase_mix = ideal_mixing_nd(y_grid)

    # Phase IDs for the phase grid
    if n == 2:
        phase_ids = np.array([f"{phase_label}_y={y[0]:.4f}" for y in y_grid])
    else:
        phase_ids = np.array(
            [phase_label + "_" + "_".join(f"{metals[k]}={y_grid[r, k]:.4f}" for k in range(n - 1)) for r in range(n_pts)]
        )

    phase_ids = np.concatenate([fids, phase_ids])
    phase_kinds = np.concatenate([np.full(n_fixed, "fixed"), np.full(n_pts, "phase")])
    pmetal = np.vstack([fmetal, y_grid])
    pphase_element = np.concatenate(
        [
            fphase_element,
            phase_element_stoichiometry * np.ones(n_pts),
        ]
    )
    pO = np.concatenate([fO, np.zeros(n_pts)])
    A_eq = np.vstack([pmetal.T, pphase_element[np.newaxis, :]])
    phase_y = np.concatenate([np.full(n_fixed, np.nan), y_grid[:, 0]])
    nan_f = np.full((n_fixed, n), np.nan)
    phase_y_nd = np.vstack([nan_f, y_grid])

    suffix = phase_label[len("".join(metals)) :] if phase_label.startswith("".join(metals)) else ""

    return {
        "phase_ids": phase_ids,
        "phase_kinds": phase_kinds,
        "phase_y": phase_y,
        "phase_y_nd": phase_y_nd,
        "phase_metal_stoich": pmetal,
        "phase_element_stoich": pphase_element,
        "phase_O": pO,
        "fixed_energy_formula": fef,
        "phase_H0": phase_H0,
        "phase_mix_shape": phase_mix,
        "y_grid": y_grid,
        "h_endpoints": h_ep,
        "binary_coeffs": bc,
        "ternary_coeff": tc,
        "A_eq": A_eq,
        "n_fixed": n_fixed,
        "n_metals": n,
        "metals": metals,
        "phase_label": phase_label,
        "phase_element": phase_element,
        "phase_element_stoichiometry": phase_element_stoichiometry,
        "phase_suffix": suffix,
        "fixed_vol_per_atom": fixed_vol_per_atom,
        "fixed_formula_counts": fixed_formula_counts,
    }


# ---------------------------------------------------------------------------
# 2-sublattice (CEF / MAX) helpers — local copies avoid circular imports with
# blade.analysis._sublattice which itself imports from this package's _grid.
# ---------------------------------------------------------------------------


def _product_simplex_grid_2subl(
    subl_metals: list[list[str]],
    subl_mults: list[int],
    y_step: float,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Cartesian product of per-sublattice simplex grids, multiplicity-weighted."""
    total = sum(subl_mults)
    grids = [
        simplex_grid_nd(len(metals), y_step) if len(metals) >= 2 else np.ones((1, 1), dtype=float)
        for metals in subl_metals
    ]
    rows = []
    for pts in _iterproduct(*grids):
        row = np.concatenate([pts[i] * subl_mults[i] / total for i in range(len(subl_metals))])
        rows.append(row)
    return np.array(rows, dtype=float), grids


def _cef_enthalpy_2subl(y_a_grid: np.ndarray, y_b_grid: np.ndarray, h_ref: np.ndarray) -> np.ndarray:
    """CEF enthalpy over the product grid: H0[k] = Σ_ij y_a_i * y_b_j * H_ref[i,j]."""
    return (y_a_grid @ h_ref @ y_b_grid.T).ravel()


def _ideal_mixing_2subl(y_a_grid: np.ndarray, y_b_grid: np.ndarray, subl_mults: list[int]) -> np.ndarray:
    """2-sublattice configurational entropy shape over the product grid."""
    total = sum(subl_mults)
    s_a = ideal_mixing_nd(y_a_grid)  # (n_a_pts,)
    s_b = ideal_mixing_nd(y_b_grid)  # (n_b_pts,)
    s_mix = (subl_mults[0] * s_a[:, np.newaxis] + subl_mults[1] * s_b[np.newaxis, :]) / total
    return s_mix.ravel()


def _load_nd_2subl(
    pt: "pd.DataFrame",
    rk: "pd.DataFrame",
    metals: list[str],
    phase_label: str,
    y_step: float,
    phase_element: str | None,
    phase_element_stoichiometry: float,
) -> dict:
    """2-sublattice path for :func:`_load_nd`.

    Reads cross-product endmember energies from the ``endmember`` rows in *rk*
    and sublattice multiplicities from the ``_subl_mults`` metadata row.
    Returns a dict with the same keys as :func:`_load_nd` so the LP solver
    (:mod:`_compute`) is unaffected.

    Args:
        pt: Already-loaded phase-table DataFrame (fixed phases / oxides).
        rk: Already-loaded interaction-coefficient DataFrame.
        metals: Ordered list of all metal symbols across all sublattices.
        phase_label: Phase-family label string.
        y_step: Composition grid spacing.
        phase_element: Fixed non-metal element (e.g. ``"C"``), or None.
        phase_element_stoichiometry: Atoms of *phase_element* per metal atom.

    Returns:
        Phase-data dict identical in structure to :func:`_load_nd`.
    """
    n = len(metals)
    phase_element_stoichiometry = float(phase_element_stoichiometry) if phase_element else 0.0

    # ---- Parse sublattice mults from metadata row ---------------------------
    mults_rows = rk[rk["pair"] == "_subl_mults"]
    if mults_rows.empty:
        raise ValueError("2-subl rk CSV is missing the _subl_mults metadata row.")
    mults_str = str(mults_rows["term"].iloc[0])  # e.g. "a=2;b=1"
    subl_letters: list[str] = []
    subl_mults: list[int] = []
    for part in mults_str.split(";"):
        letter, val = part.split("=")
        subl_letters.append(letter.strip())
        subl_mults.append(int(val.strip()))

    # ---- Parse endmember energies into h_ref[i,j] ---------------------------
    em_rows = rk[rk["pair"] == "endmember"]
    # Determine subl_metals by collecting symbols from term col ("Ti_Al" → [Ti], [Al])
    subl_metals_sets: list[set[str]] = [set() for _ in subl_letters]
    for _, row in em_rows.iterrows():
        parts = str(row["term"]).split("_")
        for idx, sym in enumerate(parts):
            if idx < len(subl_metals_sets):
                subl_metals_sets[idx].add(sym)
    subl_metals = [sorted(s) for s in subl_metals_sets]

    n_a, n_b = len(subl_metals[0]), len(subl_metals[1])
    h_ref = np.full((n_a, n_b), np.nan)
    for _, row in em_rows.iterrows():
        parts = str(row["term"]).split("_")
        m_a, m_b = parts[0], parts[1]
        i = subl_metals[0].index(m_a)
        j = subl_metals[1].index(m_b)
        h_ref[i, j] = float(row["value_eV_per_formula"])

    # ---- Product simplex grid -----------------------------------------------
    y_overall, (y_a_grid, y_b_grid) = _product_simplex_grid_2subl(subl_metals, subl_mults, y_step)
    n_pts = len(y_overall)

    phase_H0 = _cef_enthalpy_2subl(y_a_grid, y_b_grid, h_ref)
    phase_mix = _ideal_mixing_2subl(y_a_grid, y_b_grid, subl_mults)

    # ---- Fixed phases from phase_table CSV ----------------------------------
    fids = pt["phase_id"].astype(str).values
    fmetal = np.column_stack([pt[m].astype(float).values for m in metals])
    fphase_element = (
        pt[phase_element].astype(float).values
        if phase_element and phase_element in pt.columns
        else np.zeros(len(pt))
    )
    fO = pt["O"].astype(float).values
    fatoms = fmetal.sum(axis=1) + fphase_element + fO
    fef = pt["energy_eV_per_atom"].astype(float).values * fatoms
    n_fixed = len(fids)

    fixed_vol_per_atom = pt["vol_per_atom"].astype(float).values if "vol_per_atom" in pt.columns else None

    all_els = metals + ([phase_element] if phase_element else []) + ["O"]
    all_counts = [fmetal[:, j] for j in range(len(metals))]
    all_counts += [fphase_element] if phase_element else []
    all_counts += [fO]
    fixed_formula_counts: list[dict] = []
    for i in range(n_fixed):
        counts = {el: float(all_counts[k][i]) for k, el in enumerate(all_els) if float(all_counts[k][i]) > 0}
        fixed_formula_counts.append(counts)

    # ---- Combined arrays for LP solver --------------------------------------
    pmetal = np.vstack([fmetal, y_overall])
    pphase_element_arr = np.concatenate([fphase_element, phase_element_stoichiometry * np.ones(n_pts)])
    pO = np.concatenate([fO, np.zeros(n_pts)])
    A_eq = np.vstack([pmetal.T, pphase_element_arr[np.newaxis, :]])

    # Phase IDs for grid points
    phase_ids_grid = np.array(
        [
            phase_label + "_" + "_".join(f"{metals[k]}={y_overall[r, k]:.4f}" for k in range(n - 1))
            for r in range(n_pts)
        ]
    )
    phase_ids = np.concatenate([fids, phase_ids_grid])
    phase_kinds = np.concatenate([np.full(n_fixed, "fixed"), np.full(n_pts, "phase")])

    nan_f = np.full((n_fixed, n), np.nan)
    phase_y_nd = np.vstack([nan_f, y_overall])
    phase_y = np.concatenate([np.full(n_fixed, np.nan), y_overall[:, 0]])

    suffix = phase_label[len("".join(metals)) :] if phase_label.startswith("".join(metals)) else ""

    # h_endpoints stub — not meaningful for 2-subl but kept for key compatibility
    h_ep = h_ref.ravel()[: n] if h_ref.size >= n else np.zeros(n)

    return {
        "phase_ids": phase_ids,
        "phase_kinds": phase_kinds,
        "phase_y": phase_y,
        "phase_y_nd": phase_y_nd,
        "phase_metal_stoich": pmetal,
        "phase_element_stoich": pphase_element_arr,
        "phase_O": pO,
        "fixed_energy_formula": fef,
        "phase_H0": phase_H0,
        "phase_mix_shape": phase_mix,
        "y_grid": y_overall,
        "h_endpoints": h_ep,
        "binary_coeffs": {},
        "ternary_coeff": 0.0,
        "A_eq": A_eq,
        "n_fixed": n_fixed,
        "n_metals": n,
        "metals": metals,
        "phase_label": phase_label,
        "phase_element": phase_element,
        "phase_element_stoichiometry": phase_element_stoichiometry,
        "phase_suffix": suffix,
        "fixed_vol_per_atom": fixed_vol_per_atom,
        "fixed_formula_counts": fixed_formula_counts,
        # 2-subl extras for diagnostics
        "h_ref_2subl": h_ref,
        "subl_metals": subl_metals,
        "subl_mults": subl_mults,
    }
