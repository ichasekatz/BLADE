"""Cache-validation helpers for SystemAnalyzer."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from .utils import csv_has_rows


def _ensure_imports() -> None:
    """Add the framework directory and its python subdirectory to sys.path."""
    root = Path(__file__).parent
    for p in [str(root), str(root / "python")]:
        if p not in sys.path:
            sys.path.insert(0, p)


def _composition_from_slice(x_axis: float, n_metals: int, axis: int, remainder_weights):
    """Compute a composition vector from a 1-D slice parameterisation.

    Args:
        x_axis: Mole fraction assigned to the slice axis component.
        n_metals: Total number of metal components.
        axis: Index of the slice-axis component.
        remainder_weights: Sequence of weights distributing (1 - x_axis) among
            the non-axis components.  Weights are normalised internally.

    Returns:
        numpy array of shape (n_metals,) summing to 1.
    """
    comp = np.zeros(n_metals, dtype=float)
    comp[axis] = x_axis
    rest = max(0.0, 1.0 - x_axis)
    others = [i for i in range(n_metals) if i != axis]
    weights = np.asarray(remainder_weights, dtype=float)
    weights = weights / weights.sum() if weights.sum() > 0 else np.ones(len(others)) / len(others)
    for i, w in zip(others, weights, strict=False):
        comp[i] = rest * w
    return comp


def _cache_columns_match(csv_path, expected_columns) -> bool:
    """Return True only when cached coordinates exactly match the requested grid.

    Args:
        csv_path: Path to the CSV cache file.
        expected_columns: Mapping of column name → 1-D array of expected values.

    Returns:
        True if all columns match within atol=1e-10; False on any mismatch or
        I/O error.
    """
    import pandas as pd

    try:
        columns = list(expected_columns)
        df = pd.read_csv(csv_path, usecols=columns)
        if len(df) != len(next(iter(expected_columns.values()))):
            return False
        return all(
            np.allclose(
                df[column].to_numpy(dtype=float),
                np.asarray(expected, dtype=float),
                atol=1e-10,
                rtol=0.0,
                equal_nan=True,
            )
            for column, expected in expected_columns.items()
        )
    except (OSError, KeyError, TypeError, ValueError):
        return False


def _component_cache_matches(config, csv_path) -> bool:
    """Return True when the cache was built with the current component threshold.

    Args:
        config: Config instance providing ``include_0p01_to_0p05_components``.
        csv_path: Path to the CSV cache file.

    Returns:
        True if the stored ``component_presence_threshold`` column matches
        config; False otherwise.
    """
    try:
        import pandas as pd

        values = pd.read_csv(csv_path, usecols=["component_presence_threshold"])["component_presence_threshold"].to_numpy(
            dtype=float
        )
        expected = 0.01 if config.include_0p01_to_0p05_components else 0.05
        return len(values) > 0 and np.allclose(values, expected)
    except Exception:
        return False


def _slice_cache_matches(config, csv_path, T, spec, metals, x_values, mu_values) -> bool:
    """Return True when the μO–x slice CSV is valid and matches the current grid.

    Args:
        config: Config instance (used for ``skip_if_analysis_exists`` and
            ``include_0p01_to_0p05_components``).
        csv_path: Path to the CSV cache file.
        T: Temperature in K.
        spec: Slice specification dict with keys ``axis`` and ``remainder``.
        metals: Sequence of metal element symbols.
        x_values: 1-D array of composition axis values.
        mu_values: 1-D array of μO values (eV per O atom).

    Returns:
        True when the CSV exists, has the correct row count, and all coordinate
        columns match the requested grid.
    """
    if not (
        config.skip_if_analysis_exists
        and csv_has_rows(csv_path, len(x_values) * len(mu_values))
        and _component_cache_matches(config, csv_path)
    ):
        return False
    compositions = np.asarray(
        [_composition_from_slice(float(x), len(metals), spec["axis"], spec["remainder"]) for x in x_values],
        dtype=float,
    )
    expected = {
        "T_K": np.full(len(x_values) * len(mu_values), float(T)),
        "muO_eV_per_O": np.tile(mu_values, len(x_values)),
        f"x_{metals[spec['axis']]}_axis": np.repeat(x_values, len(mu_values)),
    }
    expected.update({f"x_{metal}": np.repeat(compositions[:, i], len(mu_values)) for i, metal in enumerate(metals)})
    return _cache_columns_match(csv_path, expected)


def _slice_muT_cache_matches(config, csv_path, comp, axis_m, x, metals, T_values, mu_values) -> bool:
    """Return True when the μO–T slice CSV is valid and matches the current grid.

    Args:
        config: Config instance.
        csv_path: Path to the CSV cache file.
        comp: Fixed composition array for this slice point.
        axis_m: Metal symbol for the slice axis.
        x: Composition axis value for this slice point.
        metals: Sequence of metal element symbols.
        T_values: 1-D array of temperature values (K).
        mu_values: 1-D array of μO values (eV per O atom).

    Returns:
        True when the CSV exists, has the correct row count, coordinate columns
        match, and the required result columns are present.
    """
    if not (
        config.skip_if_analysis_exists
        and csv_has_rows(csv_path, len(T_values) * len(mu_values))
        and _component_cache_matches(config, csv_path)
    ):
        return False
    n_states = len(T_values) * len(mu_values)
    expected = {
        "T_K": np.repeat(T_values, len(mu_values)),
        "muO_eV_per_O": np.tile(mu_values, len(T_values)),
        f"x_{axis_m}_axis": np.full(n_states, float(x)),
    }
    expected.update({f"x_{metal}": np.full(n_states, float(comp[i])) for i, metal in enumerate(metals)})
    if not _cache_columns_match(csv_path, expected):
        return False
    try:
        import pandas as pd

        columns = set(pd.read_csv(csv_path, nrows=0).columns)
        return {
            "phase_fraction",
            "parent_phase_fraction",
            "oxide_phase_fraction",
            "assemblage_exact",
            "phase_fraction_summary",
        }.issubset(columns)
    except (OSError, ValueError):
        return False


def _onset_cache_matches(config, csv_path, comp_grid, comp_step, mu_values, T, metals) -> bool:
    """Return True when the onset AUC CSV is valid and matches the current grid.

    Args:
        config: Config instance (used for ``skip_if_analysis_exists``).
        csv_path: Path to the CSV cache file.
        comp_grid: List of composition arrays spanning the simplex.
        comp_step: Composition grid step size used to generate *comp_grid*.
        mu_values: 1-D array of μO values (eV per O atom).
        T: Temperature in K.
        metals: Sequence of metal element symbols.

    Returns:
        True when the CSV exists, has the correct row count, coordinate columns
        match, and the required onset result columns are present.
    """
    if not (config.skip_if_analysis_exists and csv_has_rows(csv_path, len(comp_grid))):
        return False
    compositions = np.asarray(comp_grid, dtype=float)
    n_rows = len(compositions)
    mu_step = float(mu_values[1] - mu_values[0]) if len(mu_values) > 1 else np.nan
    expected = {
        "T_K": np.full(n_rows, float(T)),
        "comp_step": np.full(n_rows, float(comp_step)),
        "muO_min": np.full(n_rows, float(mu_values[0])),
        "muO_max": np.full(n_rows, float(mu_values[-1])),
        "muO_step": np.full(n_rows, mu_step),
    }
    expected.update({f"x_{metal}": compositions[:, i] for i, metal in enumerate(metals)})
    if not _cache_columns_match(csv_path, expected):
        return False
    try:
        import pandas as pd

        columns = set(pd.read_csv(csv_path, nrows=0).columns)
        return {
            "feasible_muO_max_eV",
            "first_infeasible_muO_eV",
            "feasible_muO_points",
            "parent_phase_fraction_auc_eV",
            "oxide_phase_fraction_auc_eV",
        }.issubset(columns)
    except (OSError, ValueError):
        return False
