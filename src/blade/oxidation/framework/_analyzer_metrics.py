"""LP-solver and scalar-metric helpers for SystemAnalyzer."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def _ensure_imports() -> None:
    """Add the framework directory and its python subdirectory to sys.path."""
    root = Path(__file__).parent
    for p in [str(root), str(root / "python")]:
        if p not in sys.path:
            sys.path.insert(0, p)


def _solve_metrics(config, pd_data, comp, T, mu_o, phase_G):
    """Solve the grand-potential LP and return per-state scalar metrics.

    Args:
        config: Config instance providing ``active_threshold`` and
            ``include_0p01_to_0p05_components``.
        pd_data: Phase-data dict produced by ``load_phases()``.
        comp: Metal composition vector (mole fractions, length = n_metals).
        T: Temperature in K (unused directly; caller pre-computes *phase_G*).
        mu_o: Chemical potential of oxygen (eV per O atom).
        phase_G: Pre-computed Gibbs energies for the flexible phase SQS
            structures (length = n_SQS).

    Returns:
        Dict with keys: ``ok``, ``amounts``, ``fracs``, ``omega``,
        ``absorbed_O``, ``phase_fraction``, ``parent_phase_fraction``,
        ``oxide_phase_fraction``, ``exact_label``, ``family_label``.
        All numeric fields are ``np.nan`` when the LP is infeasible.
    """
    _ensure_imports()
    from .thermodynamics import KB, build_assemblage_labels, solve_grand_lp  # noqa: F401

    phase_O = pd_data["phase_O"]
    n_fixed = pd_data["n_fixed"]
    fixed_ef = pd_data["fixed_energy_formula"]
    grand = np.concatenate([fixed_ef - phase_O[:n_fixed] * mu_o, phase_G])
    b_eq = np.concatenate([comp, [pd_data.get("phase_element_stoichiometry", 0.0)]])
    amounts, omega, ok = solve_grand_lp(pd_data["A_eq"], b_eq, grand)
    n = len(pd_data["phase_ids"])
    threshold = config.active_threshold
    if not ok:
        return {
            "ok": False,
            "amounts": np.full(n, np.nan),
            "fracs": np.full(n, np.nan),
            "omega": np.nan,
            "absorbed_O": np.nan,
            "phase_fraction": np.nan,
            "parent_phase_fraction": np.nan,
            "oxide_phase_fraction": np.nan,
            "exact_label": "no feasible assemblage",
            "family_label": "no feasible assemblage",
        }
    total = float(np.nansum(amounts))
    raw_fracs = amounts / total if total > 0 else amounts.copy()
    fracs = raw_fracs.copy()
    fracs[np.abs(fracs) < threshold] = 0.0
    phase_mask = pd_data["phase_kinds"] == "phase"
    phase_y_nd = np.asarray(pd_data["phase_y_nd"], dtype=float)
    # Preserve phase identity by element set, not by its initial ratios.
    parent_mask = phase_mask & np.all(phase_y_nd > 0.01, axis=1)
    oxide_mask = np.asarray(phase_O, dtype=float) > 1e-12
    exact_label, family_label = build_assemblage_labels(
        pd_data["phase_ids"],
        pd_data["phase_kinds"],
        fracs,
        threshold,
        pd_data.get("phase_label", ""),
        metals=pd_data.get("metals", []),
        phase_suffix=pd_data.get("phase_suffix", ""),
        family_threshold=(0.01 if config.include_0p01_to_0p05_components else 0.05),
        family_values=amounts,
        family_threshold_inclusive=config.include_0p01_to_0p05_components,
    )
    return {
        "ok": True,
        "amounts": amounts,
        "fracs": fracs,
        "omega": float(omega),
        "absorbed_O": float(np.nansum(amounts * phase_O)),
        "phase_fraction": float(np.nansum(fracs[phase_mask])),
        "parent_phase_fraction": float(np.nansum(raw_fracs[parent_mask])),
        "oxide_phase_fraction": float(np.nansum(raw_fracs[oxide_mask])),
        "exact_label": exact_label,
        "family_label": family_label,
    }


def _state_key(comp, mu_o):
    """Return a stable coordinate key for an already-calculated equilibrium state.

    Args:
        comp: Composition array or sequence of floats.
        mu_o: Oxygen chemical potential (eV per O atom).

    Returns:
        Tuple of rounded floats suitable for use as a dict key.
    """
    return tuple(round(float(x), 10) for x in comp) + (round(float(mu_o), 10),)


def _load_reusable_scalar_states(config, sys_cfg, T):
    """Read onset inputs from existing full-grid calculation tables.

    Scans all composition-slice CSVs and (for binary systems) the muO–x phase
    map cache to build a dict mapping (comp…, mu_o) → scalar fractions.

    Args:
        config: Config instance providing ``tables_dir``.
        sys_cfg: SystemConfig for the current system (provides ``tag`` and
            ``metals``).
        T: Temperature in K; used to locate per-temperature CSV files.

    Returns:
        Dict mapping state keys (from :func:`_state_key`) to
        ``(phase_fraction, parent_phase_fraction, oxide_phase_fraction,
        absorbed_O_atoms)`` tuples.
    """
    import pandas as pd

    metals = sys_cfg.metals
    reusable = {}
    paths = sorted(config.tables_dir.glob(f"{sys_cfg.tag}_slice_*_T{int(T)}.csv"))
    if len(metals) == 2:
        paths.append(config.tables_dir / f"{sys_cfg.tag}_muO_x_phase_map_T{int(T)}_cache.csv")

    required = {
        "muO_eV_per_O",
        "phase_fraction",
        "parent_phase_fraction",
        "oxide_phase_fraction",
        "absorbed_O_atoms",
    }
    for path in paths:
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if not required.issubset(df.columns):
            continue

        metal_cols = [f"x_{m}" for m in metals]
        has_full_composition = all(c in df.columns for c in metal_cols)
        if not has_full_composition and len(metals) != 2:
            continue

        for values in df.to_dict("records"):
            if has_full_composition:
                comp = [values[c] for c in metal_cols]
            else:
                x0 = float(values[f"x_{metals[0]}"])
                comp = [x0, 1.0 - x0]
            reusable[_state_key(comp, values["muO_eV_per_O"])] = (
                float(values["phase_fraction"]),
                float(values["parent_phase_fraction"]),
                float(values["oxide_phase_fraction"]),
                float(values["absorbed_O_atoms"]),
            )
    return reusable
