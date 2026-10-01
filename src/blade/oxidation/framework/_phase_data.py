"""Phase data loading for binary and N-metal oxidation systems."""

from __future__ import annotations

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
    bt = pd.read_csv(phase_points_file)
    rk = pd.read_csv(rk_coeff_file)

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
    }
