"""Pareto analysis of AlB₂-type diboride single-phase stability across compositions.

Enabled from the stability stage with ``[stability] run_pareto = true`` (see
:class:`blade.analysis.blade_stability.BladeStability`), or called directly via
:func:`run_pareto_analysis`.

Loads stability_raw.parquet and stability_summary.parquet from ``db_dir``,
computes per-composition metrics, assigns three distinct Pareto frontiers,
writes an Excel workbook, static and interactive plots, temperature profiles,
and prints a terminal summary of the top candidates.

Three Pareto frontiers (all computed on the same eligible composition set):

  1. Primary (``pareto_rank``)      — minimize ``mean_dG_hull``, maximize ``strict_frac``
  2. Alternative (``pareto_rank_alt``) — minimize ``mean_dG_hull``, maximize ``ss_auc``
  3. Supplemental 3-obj (``pareto_rank_3obj``) — minimize ``mean_dG_hull``,
     maximize ``strict_frac``, maximize ``ss_auc``

Additional reported metrics (not Pareto objectives):

- ``unique_frac``       — fraction of T steps where ``is_unique_homogeneous``
- ``max_dG_hull``       — maximum ΔG above hull across T steps
- ``mean_G_ss_formation`` — mean formation Gibbs free energy of the solid solution
- ``T_min_retained``    — lowest T where the AlB₂ parent structure is retained
- ``strict_T_span``     — longest continuous temperature window (K) with strict stability

Usage::

    from blade.analysis.stability_pareto import run_pareto_analysis
    run_pareto_analysis(db_dir="Files/Stability", out_dir="Files/Stability/pareto_analysis",
                        t_min=300, t_max=2000, t_step=100)

Output written to ``out_dir``:
  Diboride_Synthesizability_Ranking.xlsx
  pareto_metrics.parquet
  pareto_2D.png / .pdf
  pareto_2D_interactive.html
  pareto_3D_interactive.html
  temperature_stability.png / .pdf
  analysis_summary.txt
  temperature_profiles/    (per-composition detail plots, top N)
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# ── CONFIG ────────────────────────────────────────────────────────────────────

db_dir  = Path("Files/Stability")
out_dir = Path("pareto_analysis")

# Temperature grid of the stability database ([T_min, T_max] at T_step) and the
# sub-range used for all evaluations: metrics, Pareto ranking, shortlist, plots.
# Narrow the analysis range (e.g. 1500–2500) for a synthesis-focused analysis
# without penalizing high-T candidates for low-T instability.
# Set both with run_pareto_analysis(); the values below are the CLI defaults.
_DB_T_MIN, _DB_T_MAX, _DB_T_STEP = 200.0, 2500.0, 100.0
analysis_T_min: float = 200.0
analysis_T_max: float = 2500.0

# Keep aliases so any external code that referenced t_norm_min/t_norm_max still works.
t_norm_min = analysis_T_min
t_norm_max = analysis_T_max

# Minimum fraction of T steps with valid dG data to include a composition.
# 1.0 requires complete temperature coverage; lower values allow partial data.
min_valid_frac: float = 1.0

# Interpolate missing interior T steps for otherwise-incomplete compositions.
# Only applies to interior gaps bounded by valid T on both sides; never extrapolates.
# Interpolated compositions are flagged and reported separately.
interpolate_incomplete: bool = True
max_interpolation_gap_K: float = 500.0   # skip interpolation if gap width exceeds this
discontinuity_threshold: float = 0.05    # flag if dG jumps > 50 meV across interpolated gap

# Numerical tolerance for Pareto dominance comparisons.
# 1e-6 (1 μeV) is below any physically meaningful energy difference and avoids
# treating LP numerical zeros (~1e-15 eV) as distinct from intentional near-zero values.
pareto_eps: float = 1e-6

# Number of top compositions for terminal table, Excel Top Candidates, T profiles.
top_n: int = 20

# Maximum compositions per system to include in Top Candidates supplemental rows.
top_per_system: int = 3

# Scientific eligibility filters — applied before Pareto ranking.
# These are not mathematical pruning rules; they define the composition set
# eligible for the primary candidate ranking.
# rank_filter_min_strict_frac = 0.01 excludes compositions with zero strict
# stability (strict_frac = 0 at every temperature), which are not viable
# synthesis candidates regardless of their hull distance.
# Excluded compositions are retained in the Excel workbook for transparency.
rank_filter_max_dG_hull:     float | None = None
rank_filter_min_strict_frac: float | None = 0.01   # require at least one strictly stable T step

# Use the fast exact front algorithms (sorted O(n log n) where applicable, else chunked).
use_fast_pareto: bool = True

# Rows per Excel sheet; longer sheets are truncated (Excel's hard limit is 1,048,576).
# The full data is always in pareto_metrics.parquet.
excel_max_rows: int = 500_000

# Stop peeling after this many fronts (None = rank every composition). Compositions
# beyond the cap get pareto_rank = <NA>. Useful for millions of 4-5 metal compositions.
max_pareto_rank: int | None = None

# Print progress during Pareto ranking: logs after each front is peeled.
verbose_ranking: bool = True

# Experimental shortlist: select this many compositions per system from rank-1
# compositions that have passed all thermodynamic eligibility filters.  Ordered
# by candidate_rank (entropy-first) — NOT by thermodynamic stability rank.
shortlist_per_system: int = 3

# ── Experimental candidate ranking ───────────────────────────────────────────
# candidate_rank / global_candidate_rank are ENTROPY-PRIORITIZED SELECTION
# RANKINGS, not measures of thermodynamic stability.  They apply only to
# compositions that have already passed all stability eligibility filters
# (pareto_rank=1, strict_frac=1, unique_frac=1, dG≈0, not interpolated,
# curvature_flagged=False).  Among those thermodynamically equivalent candidates,
# higher normalized configurational entropy is preferred.
#
# Within-system (candidate_rank): s_config_norm↓ → Hessian bin↓ → FSM bin↓.
# Global (global_candidate_rank): s_config_norm↓ → FSM bin↓.
#   Hessian excluded from global sort: eigenvalue magnitudes are not comparable
#   across systems with different composition-simplex dimension (n_metals).
#
# Bin widths define the tie threshold below which entropy preference applies:
#   candidate_hessian_tol = 1 meV — below experimental reproducibility
#   candidate_fsm_tol     = 0.1 meV — below LP precision
# Entropy preference is NOT applied within a bin — it IS the primary sort key.
candidate_hessian_tol: float = 1e-3   # eV — Hessian tiebreaker bin width
candidate_fsm_tol:     float = 1e-4   # eV — FSM tiebreaker bin width

# ── Supplementary stability metrics ──────────────────────────────────────────
# Finite-separation decomposition margin (ss_fsm_* columns):
#   AlB₂ SS-only metric — competing assemblage restricted to AlB₂ SS phases at
#   Chebyshev distance ≥ fsm_min_sep.  Does NOT include non-AlB₂ competing
#   phases.  Positive value does NOT imply equilibrium stability — always
#   interpret alongside is_unique_homogeneous / pareto_rank.
#
#   Near-endmember compositions return NaN when min_sep makes the LP
#   geometrically infeasible (insufficient compositions on one side of target).
#   This is not a stability signal; equimolar interior compositions are unaffected.
#
#   Recommended values: 0.05 (one grid step, maximises coverage),
#   0.10 (excludes two nearest neighbours, default), 0.15 (conservative).
#   Run at multiple values to check ordering robustness before final selection.
fsm_min_sep: float = 0.10
# Compute finite-separation margin and curvature for rank-1 only (fast) or all ranked.
# Rank-1 only recommended: ~30 seconds vs ~10 minutes for all ranked.
supp_metrics_rank1_only: bool = True

# Primary Pareto objectives — always used for pareto_rank and system_pareto_rank.
# ss_auc is always computed and used in the alternative and 3-objective frontiers.
_PRIMARY_OBJECTIVES: list[tuple[str, str]] = [
    ("mean_dG_hull", "min"),
    ("strict_frac",  "max"),
]

# ── END CONFIG ────────────────────────────────────────────────────────────────

# Expected temperature grid for the selected analysis range (100 K increments).
# Determined by analysis_T_min / analysis_T_max. Used to detect T steps missing
# from the parquet entirely and to normalize all fraction / AUC metrics.
_EXPECTED_T_GRID: list[float] = [
    float(t) for t in np.arange(analysis_T_min, analysis_T_max + _DB_T_STEP / 2, _DB_T_STEP)
]
_N_EXPECTED_T: int = len(_EXPECTED_T_GRID)
_EXPECTED_T_SET: set[float] = set(_EXPECTED_T_GRID)

# mean_dG_hull_binned is added by _add_pareto_columns; rounds to pareto_eps so
# compositions within LP noise sort identically on the energy criterion.
_SECONDARY_SORT = ["pareto_rank", "strict_frac", "mean_dG_hull_binned", "unique_frac"]
_SECONDARY_ASC  = [True, False, True, False]

# Human-readable label for the analysis temperature range — used in titles and summaries.
_T_RANGE_LABEL: str = f"{int(analysis_T_min)}–{int(analysis_T_max)} K"

_REQUIRED_RAW_COLS: set[str] = {
    "system", "metals", "composition", "comp_label", "T",
    "ss_fraction", "dG_above_hull",
    "is_parent_structure_retained",
    "is_single_phase_strict",
    "is_unique_homogeneous",
    "G_ss_formation", "G_formation",
}


# ── Data loading ──────────────────────────────────────────────────────────────

def _validate_schema(raw: pd.DataFrame) -> None:
    """Raise if the raw parquet is missing new-schema columns.

    Args:
        raw: Raw stability DataFrame.

    Raises:
        ValueError: If any required column is absent.
    """
    missing = _REQUIRED_RAW_COLS - set(raw.columns)
    if missing:
        raise ValueError(
            f"stability_raw.parquet is missing columns: {sorted(missing)}.\n"
            "Rerun build_stability_db.py to regenerate with the updated schema."
        )


def _load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and validate both stability parquets.

    Returns:
        Tuple of (raw, summary) DataFrames.

    Raises:
        FileNotFoundError: If either parquet is absent.
        ValueError: If raw is missing required new-schema columns.
    """
    raw_path     = db_dir / "stability_raw.parquet"
    summary_path = db_dir / "stability_summary.parquet"
    for p in (raw_path, summary_path):
        if not p.exists():
            raise FileNotFoundError(f"Not found: {p}\nRun build_stability_db.py first.")
    raw     = pd.read_parquet(raw_path)
    summary = pd.read_parquet(summary_path)
    _validate_schema(raw)
    print(
        f"Loaded {len(raw):,} raw rows  |  {len(summary):,} summary rows  |  "
        f"{raw['system'].nunique()} systems"
    )
    return raw, summary


# ── Interpolation ─────────────────────────────────────────────────────────────

def _interpolate_group(
    grp_s: pd.DataFrame,
    valid: pd.DataFrame,
) -> tuple[pd.DataFrame, int, list[str]]:
    """Linearly interpolate missing interior T steps for continuous quantities.

    Interpolates ``dG_above_hull``, ``ss_fraction``, and ``G_ss_formation`` at
    missing T steps that are bounded by valid data on both sides.  Does not
    extrapolate at leading or trailing T edges.  Boolean flags
    (``is_single_phase_strict``, ``is_unique_homogeneous``,
    ``is_parent_structure_retained``) are left NaN at interpolated steps.

    Args:
        grp_s: Full group sorted by T (may contain NaN rows).
        valid: Subset of grp_s with non-NaN dG_above_hull and ss_fraction.

    Returns:
        Tuple of (augmented_valid, n_interpolated, warnings) where
        augmented_valid contains original valid rows plus interpolated rows,
        n_interpolated is the count of new rows, and warnings lists any
        discontinuity warnings.
    """
    # Use the full expected grid so rows absent from the parquet are also interpolated.
    all_T   = np.array(_EXPECTED_T_GRID, dtype=float)
    valid_T = valid["T"].to_numpy(dtype=float)
    missing = np.array(sorted(_EXPECTED_T_SET - set(valid_T.tolist())))

    if len(missing) == 0:
        return valid, 0, []

    t_min, t_max = valid_T.min(), valid_T.max()
    interior = missing[(missing > t_min) & (missing < t_max)]
    if len(interior) == 0:
        return valid, 0, []

    # Group consecutive interior T values into runs; filter runs exceeding max gap.
    def _split_runs(ts: list[float]) -> list[list[float]]:
        runs: list[list[float]] = []
        cur: list[float] = []
        for t in ts:
            if not cur or t - cur[-1] <= _DB_T_STEP + 1e-3:
                cur.append(t)
            else:
                runs.append(cur)
                cur = [t]
        if cur:
            runs.append(cur)
        return runs

    eligible: list[float] = []
    for run in _split_runs(interior.tolist()):
        gap_K = run[-1] - run[0] + _DB_T_STEP
        if gap_K <= max_interpolation_gap_K:
            eligible.extend(run)

    if not eligible:
        return valid, 0, []

    interp_rows: list[dict] = []
    warns: list[str] = []
    template = valid.iloc[0].to_dict()

    for t in eligible:
        lo_mask = valid_T < t
        hi_mask = valid_T > t
        if not lo_mask.any() or not hi_mask.any():
            continue
        t_lo = float(valid_T[lo_mask][-1])
        t_hi = float(valid_T[hi_mask][0])
        alpha = (t - t_lo) / (t_hi - t_lo)

        row_lo = valid[valid["T"] == t_lo].iloc[0]
        row_hi = valid[valid["T"] == t_hi].iloc[0]

        dg_lo, dg_hi = float(row_lo["dG_above_hull"]), float(row_hi["dG_above_hull"])
        dg_interp    = dg_lo + alpha * (dg_hi - dg_lo)
        ss_interp    = float(np.clip(
            float(row_lo["ss_fraction"]) + alpha * (float(row_hi["ss_fraction"]) - float(row_lo["ss_fraction"])),
            0.0, 1.0,
        ))

        actual_jump = abs(dg_hi - dg_lo)
        if actual_jump > discontinuity_threshold:
            warns.append(
                f"T={t_lo:.0f}→{t_hi:.0f} K: dG jump {actual_jump:.3f} eV — "
                "possible phase transition; interpolation may be unreliable"
            )

        gsf = float("nan")
        if ("G_ss_formation" in row_lo.index
                and pd.notna(row_lo["G_ss_formation"])
                and pd.notna(row_hi["G_ss_formation"])):
            gsf = float(row_lo["G_ss_formation"]) + alpha * (
                float(row_hi["G_ss_formation"]) - float(row_lo["G_ss_formation"])
            )

        new_row = {
            **template,
            "T":                           t,
            "dG_above_hull":               dg_interp,
            "ss_fraction":                 ss_interp,
            "G_ss_formation":              gsf,
            "is_single_phase_strict":      float("nan"),
            "is_unique_homogeneous":       float("nan"),
            "is_parent_structure_retained": float("nan"),
            "_interpolated":               True,
        }
        interp_rows.append(new_row)

    if not interp_rows:
        return valid, 0, []

    augmented = pd.concat(
        [valid, pd.DataFrame(interp_rows)], ignore_index=True
    ).sort_values("T").reset_index(drop=True)
    return augmented, len(interp_rows), warns


# ── Metric computation ────────────────────────────────────────────────────────

def _compute_pareto_metrics(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute per-composition metrics from the raw stability table.

    For each (system, composition) computes:

    - ``mean_dG_hull``       mean dG_above_hull over valid T (eV/f.u.) — minimize.
    - ``max_dG_hull``        maximum dG_above_hull over valid T.
    - ``ss_auc``             trapz(ss_fraction, T) / (T_max − T_min) — maximize.
    - ``ss_auc_K``           unnormalized AUC in K.
    - ``strict_frac``        fraction of T steps where is_single_phase_strict — maximize.
    - ``unique_frac``        fraction of T steps where is_unique_homogeneous.
    - ``strict_T_span``      longest continuous T window (K) with is_single_phase_strict.
    - ``mean_G_ss_formation`` mean G_ss_formation over valid T.
    - ``T_min_retained``     lowest T with is_parent_structure_retained.

    Compositions with fewer than ``min_valid_frac`` of T steps with valid dG
    data are excluded and returned separately as incomplete.

    Args:
        raw: Raw stability DataFrame (must pass ``_validate_schema``).

    Returns:
        Tuple of (metrics, incomplete) DataFrames.
    """
    t_range       = t_norm_max - t_norm_min
    records:        list[dict]        = []
    incomplete:     list[dict]        = []
    interp_rows_all: list[pd.DataFrame] = []   # interpolated rows for plotting

    groups  = list(raw.groupby(["system", "composition"], sort=False))
    n_total = len(groups)
    print(f"Computing metrics for {n_total:,} (system, composition) pairs …")

    for k, ((system, comp), grp) in enumerate(groups):
        if k % 10_000 == 0 and k > 0:
            print(f"  {k:>8,} / {n_total:,} …", flush=True)

        grp_s   = grp.sort_values("T")
        valid   = grp_s.dropna(subset=["dG_above_hull", "ss_fraction"])
        n_steps = len(grp_s)
        n_valid = len(valid)
        row0    = grp_s.iloc[0]
        metals  = str(row0["metals"])
        n_m     = row0.get("n_metals", len(metals.split(",")))

        # Validate against expected grid — detect rows absent from parquet entirely.
        observed_T_set = set(grp_s["T"].tolist())
        grid_missing_T = sorted(_EXPECTED_T_SET - observed_T_set)      # absent rows
        data_missing_T = sorted(                                        # NaN rows
            observed_T_set - set(valid["T"].tolist())
        )
        all_missing_T  = sorted(_EXPECTED_T_SET - set(valid["T"].tolist()))
        n_grid_missing = len(grid_missing_T)

        base = {
            "system":          system,
            "metals":          metals,
            "n_metals":        int(n_m),
            "composition":     comp,
            "comp_label":      row0.get("comp_label", comp),
            "n_valid_T":       n_valid,
            "n_total_T":       n_steps,
            "n_grid_missing":  n_grid_missing,
        }
        x_cols = [c for c in grp_s.columns if c.startswith("x_")]
        for xc in x_cols:
            base[xc] = float(grp_s[xc].iloc[0]) if not grp_s[xc].isna().all() else float("nan")

        if n_valid / _N_EXPECTED_T < min_valid_frac:
            missing_T = all_missing_T
            base["missing_T_steps"]       = missing_T
            base["grid_missing_T_steps"]  = grid_missing_T
            base["data_missing_T_steps"]  = data_missing_T

            if interpolate_incomplete and n_valid > 0:
                aug, n_interp, interp_warns = _interpolate_group(grp_s, valid)
                if len(aug) / _N_EXPECTED_T >= min_valid_frac:
                    for w in interp_warns:
                        print(f"  [interp-warn] {system}  {comp}  {w}")
                    base["n_interpolated_T"] = n_interp
                    base["interpolated"]     = True
                    base["interp_warnings"]  = "; ".join(interp_warns) if interp_warns else ""
                    # Collect interpolated rows for plot use.
                    if "_interpolated" in aug.columns:
                        new_rows = aug[aug["_interpolated"].fillna(False) == True].copy()
                        if not new_rows.empty:
                            interp_rows_all.append(new_rows)
                    valid   = aug
                    n_valid = len(aug)
                    # fall through to metric computation below
                else:
                    incomplete.append(base)
                    desc = (
                        f"grid-absent={grid_missing_T}, data-invalid={data_missing_T}"
                        if grid_missing_T else f"missing T={missing_T}"
                    )
                    print(
                        f"  [incomplete] {system}  {comp}  {desc}  "
                        f"({n_valid}/{_N_EXPECTED_T} valid)"
                    )
                    continue
            else:
                incomplete.append(base)
                desc = (
                    f"grid-absent={grid_missing_T}, data-invalid={data_missing_T}"
                    if grid_missing_T else f"missing T={missing_T}"
                )
                print(
                    f"  [incomplete] {system}  {comp}  {desc}  "
                    f"({n_valid}/{_N_EXPECTED_T} valid)"
                )
                continue

        # Observed-only subset for boolean metrics (excludes interpolated rows).
        if "_interpolated" in valid.columns:
            obs = valid[valid["_interpolated"].fillna(False) != True]
        else:
            obs = valid
        n_obs     = len(obs)
        obs_valid = n_obs if n_obs > 0 else 1   # guard against zero division

        t_vals  = valid["T"].to_numpy(dtype=float)
        ss_vals = valid["ss_fraction"].fillna(0.0).to_numpy(dtype=float)
        ss_auc_K = float(np.trapezoid(ss_vals, t_vals)) if len(t_vals) >= 2 else 0.0

        # Longest continuous temperature window where is_single_phase_strict is True.
        # Uses observed rows only (NaN at interpolated steps → fillna(False) = conservative).
        # Gaps in the T grid (missing steps) break the continuous interval.
        strict_flags = obs["is_single_phase_strict"].fillna(False).to_numpy(dtype=bool)
        obs_T_vals   = obs["T"].to_numpy(dtype=float)
        strict_T_span = 0.0
        if strict_flags.any():
            run_start_T: float | None = None
            prev_obs_T:  float | None = None
            for flag, t in zip(strict_flags, obs_T_vals):
                # Gap > 100 K in observed T values means missing T steps break the interval.
                if prev_obs_T is not None and t - prev_obs_T > _DB_T_STEP + 1e-3:
                    run_start_T = None
                if flag:
                    if run_start_T is None:
                        run_start_T = t
                    span = t - run_start_T
                    strict_T_span = max(strict_T_span, span)
                else:
                    run_start_T = None
                prev_obs_T = t

        obs_strict = int(obs["is_single_phase_strict"].fillna(False).sum())
        obs_unique = int(obs["is_unique_homogeneous"].fillna(False).sum())

        records.append({
            **base,
            "n_interpolated_T":    base.get("n_interpolated_T", 0),
            "interpolated":        base.get("interpolated", False),
            "interp_warnings":     base.get("interp_warnings", ""),
            "mean_dG_hull":        float(valid["dG_above_hull"].mean()),
            "max_dG_hull":         float(valid["dG_above_hull"].max()),
            "ss_auc":              ss_auc_K / t_range,
            "ss_auc_K":            ss_auc_K,
            "strict_frac":         float(obs_strict) / obs_valid,
            "strict_frac_lb":      float(obs_strict) / _N_EXPECTED_T,
            "unique_frac":         float(obs_unique) / obs_valid,
            "unique_frac_lb":      float(obs_unique) / _N_EXPECTED_T,
            "strict_T_span":       strict_T_span,
            "mean_G_ss_formation": (
                float(valid["G_ss_formation"].mean())
                if "G_ss_formation" in valid.columns else float("nan")
            ),
            "T_min_retained": (
                float(obs.loc[obs["is_parent_structure_retained"].fillna(False), "T"].min())
                if obs["is_parent_structure_retained"].fillna(False).any() else float("nan")
            ),
        })

    metrics    = pd.DataFrame(records)
    incompl_df = pd.DataFrame(incomplete)
    interp_df  = (
        pd.concat(interp_rows_all, ignore_index=True) if interp_rows_all
        else pd.DataFrame()
    )
    n_interp   = int(metrics["interpolated"].fillna(False).sum()) if not metrics.empty else 0
    print(
        f"  → {len(metrics):,} with sufficient T coverage  "
        f"({n_interp} recovered via interpolation)  |  "
        f"{len(incompl_df):,} excluded (< {min_valid_frac:.0%} valid T steps on "
        f"{_N_EXPECTED_T}-step expected grid)"
    )
    return metrics, incompl_df, interp_df


def _nanmean_rows(m: np.ndarray) -> np.ndarray:
    """Row means skipping NaN (pandas ``mean`` semantics); NaN for all-NaN rows."""
    count = np.sum(~np.isnan(m), axis=1)
    total = np.nansum(m, axis=1)
    return np.where(count > 0, total / np.maximum(count, 1), np.nan)


def _compute_pareto_metrics_fast(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Vectorized :func:`_compute_pareto_metrics` with identical outputs.

    Compositions observed at every expected temperature with valid ``dG_above_hull``
    and ``ss_fraction`` are reshaped to (n_compositions, n_T) arrays and processed at
    once; any other composition (gaps, NaN, interpolation) goes through the original
    per-composition code, so interpolation and incompleteness handling are unchanged.

    Args:
        raw: Raw stability DataFrame (must pass ``_validate_schema``).

    Returns:
        Tuple of (metrics, incomplete, interpolated_rows) as the original function.
    """
    n_t = _N_EXPECTED_T
    grid = np.asarray(_EXPECTED_T_GRID, dtype=float)
    key = raw["system"].astype(str) + "\x1f" + raw["composition"].astype(str)
    codes, _ = pd.factorize(key, sort=False)  # first-appearance order, like groupby(sort=False)
    n_groups = int(codes.max()) + 1 if len(codes) else 0

    t_idx = np.searchsorted(grid, raw["T"].to_numpy(dtype=float))
    on_grid = (t_idx < n_t) & np.isclose(grid[np.minimum(t_idx, n_t - 1)], raw["T"].to_numpy(dtype=float))
    valid_row = on_grid & raw["dG_above_hull"].notna().to_numpy() & raw["ss_fraction"].notna().to_numpy()
    n_rows = np.bincount(codes, minlength=n_groups)
    n_valid_rows = np.bincount(codes, weights=valid_row, minlength=n_groups)
    complete = (n_rows == n_t) & (n_valid_rows == n_t)
    # Exactly one valid row per (composition, T) slot.
    slot = codes.astype(np.int64) * n_t + np.minimum(t_idx, n_t - 1)
    dup = np.bincount(slot[complete[codes]], minlength=n_groups * n_t) if complete.any() else np.zeros(0)
    if len(dup) and (dup > 1).any():
        complete[np.unique((np.flatnonzero(dup > 1)) // n_t)] = False

    sel = complete[codes]
    sub = raw[sel]
    order = np.lexsort((t_idx[sel], codes[sel]))
    sub = sub.iloc[order]
    grp_codes = codes[sel][order][::n_t]

    def mat(col: str, fill=None) -> np.ndarray:
        v = sub[col]
        if fill is not None:
            v = v.fillna(fill)
        return v.to_numpy().reshape(-1, n_t)

    dg = mat("dG_above_hull").astype(float)
    ss = mat("ss_fraction", 0.0).astype(float)
    strict = mat("is_single_phase_strict", False).astype(bool)
    unique = mat("is_unique_homogeneous", False).astype(bool)
    retained = mat("is_parent_structure_retained", False).astype(bool)

    ss_auc_k = np.trapezoid(ss, grid, axis=1) if n_t >= 2 else np.zeros(len(ss))
    run = np.zeros(len(strict))
    longest = np.zeros(len(strict))
    for j in range(n_t):  # longest consecutive run of strict-stable T steps
        run = np.where(strict[:, j], run + 1, 0)
        longest = np.maximum(longest, run)
    strict_span = np.where(longest > 0, (longest - 1) * _DB_T_STEP, 0.0)
    first_ret = np.argmax(retained, axis=1)
    t_min_ret = np.where(retained.any(axis=1), grid[first_ret], np.nan)

    first = sub.iloc[::n_t]
    metals = first["metals"].astype(str)
    n_m = first["n_metals"] if "n_metals" in first.columns else metals.str.count(",") + 1
    fast = pd.DataFrame({
        "system":          first["system"].to_numpy(),
        "metals":          metals.to_numpy(),
        "n_metals":        n_m.astype(int).to_numpy(),
        "composition":     first["composition"].to_numpy(),
        "comp_label":      (first["comp_label"] if "comp_label" in first.columns else first["composition"]).to_numpy(),
        "n_valid_T":       n_t,
        "n_total_T":       n_t,
        "n_grid_missing":  0,
    })
    for xc in [c for c in raw.columns if c.startswith("x_")]:
        fast[xc] = first[xc].astype(float).to_numpy()
    fast = fast.assign(
        n_interpolated_T=0,
        interpolated=False,
        interp_warnings="",
        mean_dG_hull=dg.mean(axis=1),
        max_dG_hull=dg.max(axis=1),
        ss_auc=ss_auc_k / (t_norm_max - t_norm_min),
        ss_auc_K=ss_auc_k,
        strict_frac=strict.sum(axis=1) / n_t,
        strict_frac_lb=strict.sum(axis=1) / n_t,
        unique_frac=unique.sum(axis=1) / n_t,
        unique_frac_lb=unique.sum(axis=1) / n_t,
        strict_T_span=strict_span,
        mean_G_ss_formation=(_nanmean_rows(mat("G_ss_formation").astype(float))
                             if "G_ss_formation" in sub.columns else np.nan),
        T_min_retained=t_min_ret,
    )
    fast["_order"] = grp_codes

    rest = raw[~sel]
    if len(rest):
        print(f"  {int((~complete).sum()):,} compositions with gaps → per-composition path")
        slow, incompl_df, interp_df = _compute_pareto_metrics(rest)
        if len(slow):
            slow_codes = pd.Series(np.arange(n_groups), index=pd.Index(key.iloc[np.unique(codes, return_index=True)[1]]))
            slow["_order"] = (slow["system"].astype(str) + "\x1f" + slow["composition"].astype(str)).map(slow_codes).to_numpy()
    else:
        slow, incompl_df, interp_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    metrics = pd.concat([fast, slow], ignore_index=True) if len(slow) else fast
    metrics = metrics.sort_values("_order", kind="stable").drop(columns="_order").reset_index(drop=True)
    print(
        f"  → {len(metrics):,} with sufficient T coverage ({len(fast):,} vectorized)  |  "
        f"{len(incompl_df):,} excluded (< {min_valid_frac:.0%} valid T steps on {_N_EXPECTED_T}-step grid)"
    )
    return metrics, incompl_df, interp_df


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_metrics(raw: pd.DataFrame, metrics: pd.DataFrame) -> list[str]:
    """Run consistency checks and return a list of warning strings.

    Checks:
    - No duplicate (system, composition, T) entries in raw.
    - dG_above_hull ≥ 0 within tolerance.
    - ss_fraction ∈ [0, 1].
    - ss_auc ∈ [0, 1+eps] for complete compositions.
    - All compositions use the same number of T steps.
    - strict_frac ≤ 1.

    Args:
        raw: Raw stability DataFrame.
        metrics: Computed metrics DataFrame.

    Returns:
        List of human-readable warning strings (empty = all clear).
    """
    issues: list[str] = []

    # Duplicate T entries.
    n_dup = int(raw.duplicated(subset=["system", "composition", "T"]).sum())
    if n_dup > 0:
        issues.append(
            f"{n_dup} duplicate (system, composition, T) rows in raw data. "
            "This may skew AUC and fraction metrics."
        )

    neg_dG = raw["dG_above_hull"].dropna()
    n_neg  = int((neg_dG < -pareto_eps).sum())
    if n_neg > 0:
        worst = float(neg_dG.min())
        issues.append(
            f"{n_neg} rows have dG_above_hull < 0 (worst {worst:.2e} eV). "
            "May indicate LP infeasibility at those T steps."
        )

    bad_ss = raw["ss_fraction"].dropna()
    n_bad  = int(((bad_ss < -pareto_eps) | (bad_ss > 1 + pareto_eps)).sum())
    if n_bad > 0:
        issues.append(f"{n_bad} rows have ss_fraction outside [0, 1].")

    if "ss_auc" in metrics.columns:
        bad_auc = metrics["ss_auc"].dropna()
        n_bad_auc = int(((bad_auc < -pareto_eps) | (bad_auc > 1 + pareto_eps)).sum())
        if n_bad_auc > 0:
            issues.append(
                f"{n_bad_auc} compositions have ss_auc outside [0, 1]. "
                "Check AUC normalization — T range may not span the full data grid."
            )

    t_counts = raw.groupby(["system", "composition"])["T"].count()
    if t_counts.nunique() > 1:
        issues.append(
            f"Unequal T step counts across compositions "
            f"(min {t_counts.min()}, max {t_counts.max()}). "
            "strict_frac denominators differ; direct comparison is valid but AUC "
            "may be skewed for compositions with fewer T steps."
        )

    if "strict_frac" in metrics.columns:
        if (metrics["strict_frac"] > 1 + pareto_eps).any():
            issues.append("strict_frac > 1 detected; check is_single_phase_strict column.")

    # Validate that the raw data covers the expected temperature grid.
    actual_T_in_db = set(float(t) for t in raw["T"].unique())
    grid_absent = sorted(_EXPECTED_T_SET - actual_T_in_db)
    if grid_absent:
        issues.append(
            f"{len(grid_absent)} expected T steps absent from raw data for all compositions: "
            f"{grid_absent}.  These steps cannot be computed or interpolated."
        )
    extra_T = sorted(actual_T_in_db - _EXPECTED_T_SET)
    if extra_T:
        issues.append(
            f"{len(extra_T)} T steps in raw data not on the expected grid: {extra_T}.  "
            "Check t_norm_min / t_norm_max / step size in CONFIG."
        )

    return issues


# ── Pareto ranking ────────────────────────────────────────────────────────────

def _pareto_front_mask_fast(data: np.ndarray) -> np.ndarray:
    """Chunked vectorized Pareto front mask.

    Processes rows in chunks to bound peak memory at O(CHUNK × n × k) instead
    of O(n² × k).  Produces identical results to ``_pareto_front_mask_slow``.

    Args:
        data: Array (n, k) — all objectives in 'higher is better' space.

    Returns:
        Boolean array of length n; True = not dominated.
    """
    n, k = data.shape
    chunk_size = 2000
    dominated  = np.zeros(n, dtype=bool)

    for i_start in range(0, n, chunk_size):
        i_end = min(i_start + chunk_size, n)
        chunk = data[i_start:i_end]          # (c, k)
        c     = i_end - i_start

        all_ge = np.ones((c, n), dtype=bool)
        any_gt = np.zeros((c, n), dtype=bool)
        for d in range(k):
            # diff_d[i, j] = data[j, d] - chunk[i, d]
            diff_d = data[None, :, d] - chunk[:, None, d]   # (c, n)
            all_ge &= diff_d >= -pareto_eps
            any_gt |= diff_d >  pareto_eps

        dom = all_ge & any_gt
        for ci in range(c):
            dom[ci, i_start + ci] = False   # never self-dominated

        dominated[i_start:i_end] = np.any(dom, axis=1)

    return ~dominated


def _pareto_front_mask_slow(data: np.ndarray) -> np.ndarray:
    """Cull-algorithm Pareto front mask (original implementation).

    Args:
        data: Array (n, k) — all objectives in 'higher is better' space.

    Returns:
        Boolean array of length n; True = not dominated.
    """
    n         = len(data)
    is_pareto = np.ones(n, dtype=bool)

    for i in range(n):
        if not is_pareto[i]:
            continue
        cur_idx = np.where(is_pareto)[0]
        cur     = data[cur_idx]
        diffs   = cur - data[i]
        dom     = np.all(diffs >= -pareto_eps, axis=1) & np.any(diffs > pareto_eps, axis=1)
        self_p  = int(np.searchsorted(cur_idx, i))
        dom[self_p] = False
        if dom.any():
            is_pareto[i] = False

    return is_pareto


def _range_max_table(values: np.ndarray) -> list[np.ndarray]:
    """Sparse table for O(1) range-maximum queries over ``values``."""
    table = [values]
    span = 1
    while 2 * span <= len(values):
        prev = table[-1]
        table.append(np.maximum(prev[:-span], prev[span:]))
        span *= 2
    return table


def _range_max(table: list[np.ndarray], lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Maximum of the tabled values over half-open ranges [lo, hi); -inf if empty."""
    out = np.full(len(lo), -np.inf)
    ok = hi > lo
    if not ok.any():
        return out
    length = hi[ok] - lo[ok]
    level = np.floor(np.log2(length)).astype(int)
    res = np.empty(len(length))
    for lv in np.unique(level):
        sel = level == lv
        t = table[lv]
        res[sel] = np.maximum(t[lo[ok][sel]], t[hi[ok][sel] - (1 << lv)])
    out[ok] = res
    return out


def _dominated_2d(a: np.ndarray, c: np.ndarray, eps: float) -> np.ndarray:
    """Exact eps-dominance for two continuous objectives (higher is better).

    j dominates i  ⇔  a_j ≥ a_i − ε, c_j ≥ c_i − ε, and (a_j > a_i + ε or c_j > c_i + ε),
    split into (A) a_j > a_i + ε with c_j ≥ c_i − ε, answered with a suffix maximum of
    c over a-sorted points, and (B) |a_j − a_i| ≤ ε with c_j > c_i + ε, answered with a
    range maximum over the ε-window in a.
    """
    n = len(a)
    if n == 0:
        return np.zeros(0, dtype=bool)
    order = np.argsort(a, kind="stable")
    sa, sc = a[order], c[order]
    suffix = np.maximum.accumulate(sc[::-1])[::-1]
    suffix = np.append(suffix, -np.inf)
    k = np.searchsorted(sa, a + eps, side="right")
    dom = suffix[k] >= c - eps
    lo = np.searchsorted(sa, a - eps, side="left")
    hi = np.searchsorted(sa, a + eps, side="right")
    dom |= _range_max(_range_max_table(sc), lo, hi) > c + eps
    return dom


def _discrete_column(data: np.ndarray, eps: float) -> int | None:
    """Index of an objective whose distinct values are spaced by more than 2ε (e.g.
    strict_frac = k/24), or None. Such a column compares exactly under the ε rule."""
    for d in range(data.shape[1]):
        u = np.unique(data[:, d])
        if len(u) <= 512 and (len(u) == 1 or np.min(np.diff(u)) > 2 * eps):
            return d
    return None


def _pareto_front_mask_sorted(data: np.ndarray) -> np.ndarray | None:
    """Exact ε-dominance front mask in O(n log n) for the 2- and 3-objective cases.

    Gives the same result as :func:`_pareto_front_mask_fast` (same ε rule) for:
    two objectives (any values), or three objectives where one has well-separated
    discrete values (``strict_frac``). Returns None when neither applies.

    Args:
        data: Array (n, k) — all objectives in 'higher is better' space.

    Returns:
        Boolean front mask, or None if the fast path does not apply.
    """
    n, k = data.shape
    eps = pareto_eps
    if k == 2:
        return ~_dominated_2d(data[:, 0], data[:, 1], eps)
    if k != 3:
        return None
    d_disc = _discrete_column(data, eps)
    if d_disc is None:
        return None
    a, c = (data[:, d] for d in range(3) if d != d_disc)
    b = data[:, d_disc]
    dominated = np.zeros(n, dtype=bool)
    levels = np.unique(b)[::-1]  # best (highest) first
    higher_a = np.empty(0)
    higher_c = np.empty(0)
    for lvl in levels:
        idx = np.flatnonzero(b == lvl)
        # Same b: ordinary 2-D ε-dominance on (a, c).
        dominated[idx] = _dominated_2d(a[idx], c[idx], eps)
        # Higher b is a strict improvement; dominated if any higher point is not worse
        # by more than ε in both a and c.
        if len(higher_a):
            order = np.argsort(higher_a, kind="stable")
            sa = higher_a[order]
            suffix = np.append(np.maximum.accumulate(higher_c[order][::-1])[::-1], -np.inf)
            kk = np.searchsorted(sa, a[idx] - eps, side="left")
            dominated[idx] |= suffix[kk] >= c[idx] - eps
        higher_a = np.concatenate([higher_a, a[idx]])
        higher_c = np.concatenate([higher_c, c[idx]])
    return ~dominated


def _pareto_front_mask(data: np.ndarray) -> np.ndarray:
    """Dispatch to fast or slow Pareto front implementation.

    Uses the exact O(n log n) sorted algorithm when it applies (2 objectives, or 3 with
    one discrete objective), else the chunked O(n²) mask.

    Args:
        data: Array (n, k) — all objectives in 'higher is better' space.

    Returns:
        Boolean array of length n; True = not dominated.
    """
    if use_fast_pareto:
        mask = _pareto_front_mask_sorted(data)
        return mask if mask is not None else _pareto_front_mask_fast(data)
    return _pareto_front_mask_slow(data)


def _assign_pareto_ranks(
    df: pd.DataFrame,
    objectives: list[tuple[str, str]],
) -> np.ndarray:
    """Assign Pareto ranks by iterative front peeling (NSGA-II style).

    Args:
        df: DataFrame containing objective columns.
        objectives: List of (column, direction) where direction is ``"min"``
            or ``"max"``.

    Returns:
        Array of shape (n,) with Pareto rank (1 = Pareto front); NaN beyond
        ``max_pareto_rank`` when a cap is set.
    """
    data = np.column_stack([
        df[col].to_numpy(dtype=float) * (1.0 if d == "max" else -1.0)
        for col, d in objectives
    ])
    n         = len(data)
    ranks     = np.full(n, np.nan) if max_pareto_rank else np.zeros(n, dtype=int)
    remaining = np.ones(n, dtype=bool)
    rank      = 1
    ranked_so_far = 0
    while remaining.any() and (max_pareto_rank is None or rank <= max_pareto_rank):
        sub   = np.where(remaining)[0]
        front = _pareto_front_mask(data[sub])
        n_front = int(front.sum())
        ranks[sub[front]] = rank
        remaining[sub[front]] = False
        ranked_so_far += n_front
        if verbose_ranking:
            print(
                f"    rank {rank}: {n_front} compositions  "
                f"({ranked_so_far}/{n} ranked, {remaining.sum()} remaining)",
                flush=True,
            )
        rank += 1
    return ranks


def _add_pareto_columns(metrics: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Compute three Pareto frontiers on the same eligible composition set.

    All three frontiers use identical eligibility criteria so their results
    are directly comparable.

    Columns added:

    - ``pareto_rank``        — primary (mean_dG_hull ↓, strict_frac ↑).
    - ``pareto_rank_alt``    — alternative (mean_dG_hull ↓, ss_auc ↑).
    - ``pareto_rank_3obj``   — 3-objective (mean_dG_hull ↓, strict_frac ↑, ss_auc ↑).
    - ``system_pareto_rank`` — per-system primary rank.

    Args:
        metrics: DataFrame from ``_compute_pareto_metrics``.

    Returns:
        Tuple of (metrics_with_ranks, comparison_dict).
    """
    metrics = metrics.copy()

    # Round mean_dG_hull to pareto_eps for secondary sorting so that compositions
    # with energy differences smaller than the Pareto tolerance sort identically.
    # This avoids imposing a thermodynamic ranking based on LP numerical noise.
    metrics["mean_dG_hull_binned"] = (
        (metrics["mean_dG_hull"] / pareto_eps).round() * pareto_eps
    )

    # All three objective columns must be non-NaN to be eligible.
    all_obj_cols = ["mean_dG_hull", "strict_frac", "ss_auc"]
    valid_m = metrics[all_obj_cols].notna().all(axis=1)

    n_invalid = int((~valid_m).sum())
    if n_invalid > 0:
        warnings.warn(
            f"{n_invalid} compositions excluded from all Pareto rankings (NaN in objectives)."
        )

    # Optional physical screening thresholds — not Pareto criteria.
    rank_m = valid_m.copy()
    if rank_filter_max_dG_hull is not None:
        prev   = int(rank_m.sum())
        rank_m &= metrics["mean_dG_hull"] <= rank_filter_max_dG_hull
        n_drop  = prev - int(rank_m.sum())
        if n_drop > 0:
            print(
                f"Optional screen: removed {n_drop:,} compositions with "
                f"mean_dG_hull > {rank_filter_max_dG_hull} eV/f.u."
            )
    n_excluded_strict_frac = 0
    if rank_filter_min_strict_frac is not None:
        prev   = int(rank_m.sum())
        rank_m &= metrics["strict_frac"] >= rank_filter_min_strict_frac
        n_excluded_strict_frac = prev - int(rank_m.sum())
        if n_excluded_strict_frac > 0:
            print(
                f"Eligibility filter: excluded {n_excluded_strict_frac:,} compositions with "
                f"strict_frac < {rank_filter_min_strict_frac} (zero strict stability; "
                "retained in 'Zero-Stability' Excel sheet)"
            )

    valid = metrics[rank_m].copy().reset_index(drop=True)
    n_v   = len(valid)

    _primary   = [("mean_dG_hull", "min"), ("strict_frac", "max")]
    _alt       = [("mean_dG_hull", "min"), ("ss_auc",      "max")]
    _three_obj = [("mean_dG_hull", "min"), ("strict_frac", "max"), ("ss_auc", "max")]

    print(f"Ranking {n_v:,} compositions across 3 Pareto frontiers …")
    r_primary  = _assign_pareto_ranks(valid, _primary)
    print("  Primary (dG+strict) done.")
    r_alt      = _assign_pareto_ranks(valid, _alt)
    print("  Alternative (dG+ss_auc) done.")
    r_3obj     = _assign_pareto_ranks(valid, _three_obj)
    print("  3-objective done.")

    print("  Per-system primary ranks …")
    sys_r = np.full(n_v, np.nan)
    for _, sg in valid.groupby("system", sort=False):
        idx = sg.index.to_numpy()
        sys_r[idx] = _assign_pareto_ranks(sg, _primary)

    valid_idx = metrics[rank_m].index
    metrics.loc[valid_idx, "pareto_rank"]        = r_primary
    metrics.loc[valid_idx, "pareto_rank_alt"]    = r_alt
    metrics.loc[valid_idx, "pareto_rank_3obj"]   = r_3obj
    metrics.loc[valid_idx, "system_pareto_rank"] = sys_r

    for col in ("pareto_rank", "pareto_rank_alt", "pareto_rank_3obj", "system_pareto_rank"):
        if col in metrics.columns:
            metrics[col] = metrics[col].astype("Int64")

    # Flag compositions with high strict stability but low homogeneous-equilibrium fraction.
    # These warrant extra scrutiny: the LP may favour phase-separated AlB₂ assemblages.
    sf = metrics["strict_frac"] if "strict_frac" in metrics.columns else pd.Series(np.nan, index=metrics.index)
    uf = metrics["unique_frac"] if "unique_frac" in metrics.columns else pd.Series(np.nan, index=metrics.index)
    flag = sf.notna() & uf.notna() & (sf >= 0.5) & (uf < 0.5 * sf)
    metrics["strict_unique_flag"] = np.where(flag, "!", "")

    # Frontier overlap summary.
    front_p = set(np.where(r_primary == 1)[0])
    front_a = set(np.where(r_alt     == 1)[0])
    front_3 = set(np.where(r_3obj    == 1)[0])

    comparison = {
        "n_primary":               len(front_p),
        "n_alt":                   len(front_a),
        "n_3obj":                  len(front_3),
        "n_primary_and_alt":       len(front_p & front_a),
        "n_only_primary":          len(front_p - front_a),
        "n_only_alt":              len(front_a - front_p),
        "n_3obj_new":              len(front_3 - (front_p | front_a)),
        "n_eligible":              n_v,
        "n_excluded_strict_frac":  n_excluded_strict_frac,
        "description": (
            f"Primary front (dG+strict): {len(front_p)}  |  "
            f"Alt front (dG+ss_auc): {len(front_a)}  |  "
            f"3-obj front: {len(front_3)}  |  "
            f"In both primary+alt: {len(front_p & front_a)}  |  "
            f"Primary-only: {len(front_p - front_a)}  |  "
            f"Alt-only: {len(front_a - front_p)}"
        ),
    }
    print(comparison["description"])
    return metrics, comparison


# ── Epsilon sensitivity ───────────────────────────────────────────────────────

def _rank1_eps_sensitivity(metrics: pd.DataFrame) -> dict:
    """Validate rank-1 composition set across multiple pareto_eps tolerances.

    Runs the primary Pareto ranking with eps = 0 (exact), 1e-9, 1e-6 (current),
    and 1e-3 on the same eligible set used for the main ranking.  Reports how many
    compositions are on rank-1 at each tolerance.

    A large jump between eps=0 and eps=1e-6 means the front is sensitive to
    numerical noise.  A large jump between eps=1e-6 and eps=1e-3 means physical
    gradients are being masked at the larger tolerance.

    Args:
        metrics: Full ranked metrics DataFrame (must include ``pareto_rank``).

    Returns:
        Dict mapping eps → rank-1 count, plus a ``description`` string.
    """
    eligible = metrics.dropna(subset=["pareto_rank"])
    if eligible.empty:
        return {"description": "No ranked compositions."}

    data_max = np.column_stack([
        -eligible["mean_dG_hull"].to_numpy(dtype=float),
        eligible["strict_frac"].to_numpy(dtype=float),
    ])
    n = len(data_max)

    # Map eps → (count, frozenset of integer positions in eligible that are rank-1).
    eps_to_n:    dict[float, int]            = {}
    eps_to_set:  dict[float, frozenset[int]] = {}

    for eps in [0.0, 1e-9, 1e-6, 1e-3]:
        CHUNK     = 2000
        dominated = np.zeros(n, dtype=bool)
        for i_start in range(0, n, CHUNK):
            i_end = min(i_start + CHUNK, n)
            chunk = data_max[i_start:i_end]
            c     = i_end - i_start
            all_ge = np.ones((c, n), dtype=bool)
            any_gt = np.zeros((c, n), dtype=bool)
            for d in range(2):
                diff_d = data_max[None, :, d] - chunk[:, None, d]
                all_ge &= diff_d >= -eps
                any_gt |= diff_d >  eps
            dom = all_ge & any_gt
            for ci in range(c):
                dom[ci, i_start + ci] = False
            dominated[i_start:i_end] = np.any(dom, axis=1)
        front = frozenset(np.where(~dominated)[0].tolist())
        eps_to_n[eps]   = len(front)
        eps_to_set[eps] = front

    # Membership comparison: how similar are the fronts?
    ref_set = eps_to_set[pareto_eps]   # current eps=1e-6 as reference

    def _jaccard(a: frozenset, b: frozenset) -> float:
        u = len(a | b)
        return len(a & b) / u if u else 1.0

    membership: dict[float, dict] = {}
    for eps in [0.0, 1e-9, 1e-6, 1e-3]:
        s = eps_to_set[eps]
        membership[eps] = {
            "n":          len(s),
            "jaccard_vs_current": round(_jaccard(s, ref_set), 4),
            "only_in_this":  len(s - ref_set),
            "only_in_curr":  len(ref_set - s),
        }

    # Stability assessment.
    n_exact  = eps_to_n[0.0]
    n_curr   = eps_to_n[pareto_eps]
    n_coarse = eps_to_n[1e-3]
    j_exact  = membership[0.0]["jaccard_vs_current"]
    if j_exact < 0.9:
        note = (
            f"Front at eps=0 (n={n_exact}) shares only {j_exact:.0%} membership with "
            f"eps=1e-6 front (n={n_curr}): LP noise inflating rank-1 set."
        )
    elif n_coarse - n_curr > 0.1 * n_curr:
        note = (
            f"Front stable (Jaccard={j_exact:.0%}) from exact to eps=1e-6 "
            f"({n_exact}→{n_curr}).  "
            f"Grows to {n_coarse} at eps=1e-3: physical gradients present at ~1 meV scale."
        )
    else:
        note = (
            f"Front stable across all tolerances (Jaccard={j_exact:.0%}, "
            f"{n_exact}→{n_curr}→{n_coarse}).  Rank-1 set is robust."
        )

    result: dict = {f"eps_{eps:.0e}".replace("e-0", "e-"): v for eps, v in eps_to_n.items()}
    result["membership"]  = membership
    result["description"] = (
        f"Rank-1 count at eps=0: {n_exact}  |  eps=1e-9: {eps_to_n[1e-9]}  |  "
        f"eps=1e-6 (current): {n_curr}  |  eps=1e-3: {n_coarse}.  {note}"
    )
    return result


# ── Sensitivity comparison ────────────────────────────────────────────────────

def _sensitivity_comparison(metrics: pd.DataFrame) -> dict:
    """Compare Pareto front membership with and without interpolated compositions.

    Re-ranks the complete-only subset and compares actual front membership
    (which specific compositions enter/leave the front) rather than just sizes.

    Args:
        metrics: Full ranked metrics DataFrame (must include ``pareto_rank``).

    Returns:
        Dict with comparison statistics and a human-readable ``description``.
    """
    if "interpolated" not in metrics.columns or not metrics["interpolated"].fillna(False).any():
        return {"n_interpolated_total": 0, "description": "No interpolated compositions."}

    ranked        = metrics.dropna(subset=["pareto_rank"])
    interp_mask   = ranked["interpolated"].fillna(False)
    complete_only = ranked[~interp_mask].copy()
    interp_comps  = ranked[interp_mask].copy()

    # Interpolated comps on the full front.
    interp_on_full_front = interp_comps[interp_comps["pareto_rank"] == 1]
    n_interp_on_front    = len(interp_on_full_front)

    # Re-rank without interpolated compositions to get complete-only front.
    _primary     = [("mean_dG_hull", "min"), ("strict_frac", "max")]
    r_complete   = _assign_pareto_ranks(complete_only, _primary)
    complete_front_idx  = set(complete_only.index[r_complete == 1])
    full_front_idx      = set(ranked.index[ranked["pareto_rank"] == 1])
    noninterp_on_full   = full_front_idx - set(interp_comps.index)

    # Compositions displaced: on complete-only front but NOT on full front
    # (a newly included interpolated comp dominated them).
    displaced_idx = complete_front_idx - noninterp_on_full

    # Compositions gained: on full front but not on complete-only front
    # (shouldn't happen — adding more compositions can only dominate existing ones).
    gained_idx = noninterp_on_full - complete_front_idx

    displaced_rows = complete_only.loc[list(displaced_idx), ["system", "comp_label", "strict_frac", "mean_dG_hull"]] if displaced_idx else pd.DataFrame()
    interp_on_front_rows = interp_on_full_front[["system", "comp_label", "strict_frac", "mean_dG_hull"]]

    changed = n_interp_on_front > 0 or len(displaced_idx) > 0

    result = {
        "n_interpolated_total":      len(interp_comps),
        "n_interpolated_on_front":   n_interp_on_front,
        "n_front_complete_only":     len(complete_front_idx),
        "n_front_full":              len(full_front_idx),
        "n_displaced_by_interp":     len(displaced_idx),
        "n_gained_by_interp":        len(gained_idx),
        "interp_changed_front":      changed,
        "displaced_compositions":    displaced_rows.to_dict("records"),
        "interpolated_on_front":     interp_on_front_rows.to_dict("records"),
        "description": (
            f"{len(interp_comps)} interpolated compositions in ranking.  "
            f"{n_interp_on_front} interpolated comps on primary Pareto front.  "
            f"{len(displaced_idx)} non-interpolated comps displaced from front by "
            f"interpolated comps.  "
            f"Front size: {len(complete_front_idx)} (complete-only) vs "
            f"{len(full_front_idx)} (with interpolated).  "
            + (
                "Front membership CHANGED by interpolation."
                if changed
                else "Primary front membership unchanged by inclusion of interpolated compositions."
            )
        ),
    }
    return result


# ── Supplementary stability metrics ──────────────────────────────────────────

def _build_gs_lookup(raw: pd.DataFrame) -> dict[str, dict[float, dict[tuple, float]]]:
    """Build {system → {T → {comp_tuple → G_ss}}} lookup from raw parquet.

    Used by curvature and finite-separation margin computations.
    Only includes rows where G_ss is not NaN.
    """
    if "G_ss" not in raw.columns:
        return {}
    valid = raw[raw["G_ss"].notna()][["system", "T", "composition", "G_ss"]].copy()
    valid["_comp"] = valid["composition"].apply(
        lambda s: tuple(round(float(x), 4) for x in s.split(","))
    )
    lookup: dict = {}
    for sys_name, sys_df in valid.groupby("system", sort=False):
        lookup[sys_name] = {}
        for t_val, t_df in sys_df.groupby("T", sort=False):
            lookup[sys_name][float(t_val)] = dict(zip(t_df["_comp"], t_df["G_ss"]))
    return lookup


def _simplex_hessian_min_eig(
    comp: tuple,
    t_gs: dict,
    h: float = 0.05,
) -> float:
    """Minimum eigenvalue of the G_ss Hessian projected onto the composition simplex.

    Uses central finite differences on the existing 0.05 composition grid.
    Returns NaN if any required neighbor composition is absent.

    **Sign interpretation (coordinate-invariant for the sign itself):**

    - Positive: G_ss is locally convex along this direction — the homogeneous SS
      is a local energy minimum within the G_ss landscape (above the G_ss spinodal).
    - Negative: G_ss is locally concave — the pure G_ss landscape has a spinodal
      at this composition.  However, **a negative eigenvalue does NOT imply that
      the LP will find a competing assemblage**: if the neighbouring SS compositions
      (x ± h) are themselves dominated by fixed competing phases (e.g. binary
      borides), the lever-rule mixture involving them costs more energy than the
      target, so the equilibrium LP still selects the target as the global minimum.
      Approximately 3–4% of rank-1 compositions (``is_unique_homogeneous = True``
      at every T) show negative hessian eigenvalues — their stability is maintained
      by fixed-phase dominance of nearby SS compositions, not by G_ss convexity.
      Always cross-reference with ``is_unique_homogeneous`` and ``unique_frac`` for
      actual equilibrium stability.

    **Magnitude interpretation (coordinate-dependent):** the parametrisation
    (x₁, …, x_{n-1}) with xₙ = 1 − Σxᵢ fixes a specific coordinate system on
    the simplex.  For binary systems (n=2) the Hessian is 1×1 — d²G/dx² — and
    the eigenvalue is determined by the parametrisation up to a sign-preserving
    scale.  For ternary and higher (n≥3) the Hessian matrix and therefore
    eigenvalue magnitudes depend on which n−1 components are chosen as the free
    variables.  **Do not compare eigenvalue magnitudes across systems with
    different numbers of components.**  Within a single chemical system the sign
    is a reliable indicator of G_ss curvature; magnitude ordering is valid only
    within the same system.

    Finite-difference error is O(h²) ≈ O(2.5e-3).

    Args:
        comp: Target composition tuple (x_1, ..., x_n) summing to 1.
        t_gs: {comp_tuple → G_ss} for all compositions in this system at this T.
        h: Grid spacing (0.05 for standard BLADE grid).

    Returns:
        Minimum eigenvalue in eV / f.u. (per composition unit²).
    """
    n = len(comp)
    if n <= 1:
        return float("nan")

    g0 = t_gs.get(comp)
    if g0 is None:
        return float("nan")

    ndim = n - 1  # independent variables on the simplex

    def _neighbor(deltas: list) -> tuple:
        new = list(comp)
        for i, d in enumerate(deltas):
            new[i] = round(new[i] + d, 4)
        new[-1] = round(1.0 - sum(new[:-1]), 4)
        return tuple(new)

    hess = np.zeros((ndim, ndim))
    for i in range(ndim):
        dp = [h if j == i else 0.0 for j in range(ndim)]
        dm = [-h if j == i else 0.0 for j in range(ndim)]
        gp = t_gs.get(_neighbor(dp))
        gm = t_gs.get(_neighbor(dm))
        if gp is None or gm is None:
            return float("nan")
        hess[i, i] = (gp - 2.0 * g0 + gm) / (h * h)

    for i in range(ndim):
        for j in range(i + 1, ndim):
            corners = {}
            for si, sj, key in [(1, 1, "pp"), (1, -1, "pm"), (-1, 1, "mp"), (-1, -1, "mm")]:
                d = [h * si if k == i else h * sj if k == j else 0.0 for k in range(ndim)]
                g = t_gs.get(_neighbor(d))
                if g is None:
                    return float("nan")
                corners[key] = g
            hess[i, j] = hess[j, i] = (
                corners["pp"] - corners["pm"] - corners["mp"] + corners["mm"]
            ) / (4.0 * h * h)

    try:
        return float(np.linalg.eigvalsh(hess).min())
    except np.linalg.LinAlgError:
        return float("nan")


def _finite_sep_margin(
    comp: tuple,
    t_gs: dict,
    min_sep: float = 0.10,
) -> float:
    """AlB₂ SS-only decomposition margin: ΔG = G_competing_restricted − G_ss(target).

    Measures resistance to separating into AlB₂ solid-solution phases whose
    Chebyshev distance from the target exceeds ``min_sep``.  The competing
    assemblage satisfies mass balance and partition-of-unity constraints.

    **Scope — AlB₂ SS phases only.**  Non-AlB₂ competing phases (intermetallics,
    borides with different crystal structure, etc.) are NOT included because their
    individual Gibbs energies are not stored in the parquet.  For rank-1
    compositions where ``n_active_fixed = 0`` the original LP has already
    established that no non-SS phase is thermodynamically competitive at
    equilibrium; this metric then quantifies additional compositional robustness
    within the SS free-energy landscape.

    **A positive value does NOT imply equilibrium stability.**  If ``min_sep``
    excludes the actual lowest-energy decomposition pathway (e.g. the composition
    is at an interior of a miscibility gap whose binodal compositions lie within
    ``min_sep``), the constrained LP will return a positive value even when
    ``is_unique_homogeneous = False`` and ``dG_above_hull > 0``.  Always
    interpret this metric alongside ``is_unique_homogeneous``.

    NaN is returned when:
    - G_ss is unavailable for the target composition.
    - Fewer than ``n_metals`` eligible competing compositions exist (LP
      underdetermined); this occurs for near-endmember compositions at large
      ``min_sep`` and is a geometric infeasibility, not a stability signal.

    Args:
        comp: Target composition tuple (metal fractions summing to 1).
        t_gs: {comp_tuple → G_ss} mapping for this system at this T;
              values must be finite (no NaN).
        min_sep: Minimum Chebyshev distance to any competing composition.
                 Typical values: 0.05 (one grid step), 0.10, 0.15.

    Returns:
        ΔG in eV/f.u.  Positive = target SS favoured over constrained assemblage.
    """
    from scipy.optimize import linprog

    g0 = t_gs.get(comp)
    if g0 is None or not np.isfinite(g0):
        return float("nan")

    n_metals = len(comp)

    # Exclude compositions within min_sep AND those with non-finite G_ss.
    eligible = [
        (c, g)
        for c, g in t_gs.items()
        if np.isfinite(g)
        and max(abs(c[i] - comp[i]) for i in range(n_metals)) >= min_sep - 1e-9
    ]

    if len(eligible) < n_metals:
        return float("nan")

    n_e = len(eligible)
    c_obj = np.array([g for _, g in eligible])

    A_eq = np.zeros((n_metals + 1, n_e))
    b_eq = np.zeros(n_metals + 1)
    for i in range(n_metals):
        A_eq[i] = [c[i] for c, _ in eligible]
        b_eq[i] = comp[i]
    A_eq[n_metals] = 1.0
    b_eq[n_metals] = 1.0

    result = linprog(
        c_obj,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=[(0.0, None)] * n_e,
        method="highs",
        options={"disp": False},
    )

    if not result.success:
        return float("nan")

    return float(result.fun - g0)


def _compute_supplementary_metrics(
    raw: pd.DataFrame,
    metrics: pd.DataFrame,
) -> pd.DataFrame:
    """Compute curvature and finite-separation margin supplementary stability metrics.

    For each composition (rank-1 only if ``supp_metrics_rank1_only``), calculates
    per-temperature:

    - Minimum Hessian eigenvalue of G_ss on the composition simplex (curvature).
    - Finite-separation decomposition margin at ``fsm_min_sep`` (ΔG_decomp).

    Aggregates over ``analysis_T_min`` to ``analysis_T_max`` to produce summary columns.

    Args:
        raw: Full raw stability DataFrame.
        metrics: Per-composition metrics DataFrame with ``pareto_rank`` column.

    Returns:
        ``metrics`` copy with supplementary columns added.
    """
    # Initialise all supplementary columns up front so downstream code (Excel
    # export, parquet write) always sees them, even when computation is skipped.
    _supp_float_cols = [
        "hessian_min_eig_mean", "hessian_min_eig_min", "hessian_all_stable_frac",
        "ss_fsm_mean", "ss_fsm_min",
    ]
    out = metrics.copy()
    for col in _supp_float_cols:
        out[col] = float("nan")
    out["n_valid_T_hessian"] = 0  # integer count; 0 = no valid T steps
    out["n_valid_T_fsm"] = 0
    # curvature_flagged: True when hessian was computed at ≥1 T step but all eigenvalues
    # are negative (G_ss locally concave at every T).  These compositions are kept in
    # Pareto rankings unchanged but excluded from the experimental shortlist.
    out["curvature_flagged"] = False

    if "G_ss" not in raw.columns:
        print("  G_ss column absent — supplementary metrics skipped (columns set to NaN/0).")
        return out

    print(
        f"\nComputing supplementary metrics "
        f"({'rank-1 only' if supp_metrics_rank1_only else 'all ranked'}, "
        f"T range: {_T_RANGE_LABEL}) …"
    )

    raw_t = raw[raw["T"].between(analysis_T_min, analysis_T_max, inclusive="both")]
    gs_lookup = _build_gs_lookup(raw_t)

    if supp_metrics_rank1_only:
        target_comps = metrics[metrics["pareto_rank"] == 1][
            ["system", "composition", "n_metals"]
        ].copy()
    else:
        target_comps = metrics.dropna(subset=["pareto_rank"])[
            ["system", "composition", "n_metals"]
        ].copy()

    n_target = len(target_comps)
    t_grid = [t for t in _EXPECTED_T_GRID if analysis_T_min <= t <= analysis_T_max]
    n_t_expected = len(t_grid)
    print(f"  {n_target:,} compositions × {n_t_expected} T steps …")

    results: dict[tuple, dict] = {}

    for k, (_, row) in enumerate(target_comps.iterrows()):
        sys = row["system"]
        comp_str = row["composition"]
        comp = tuple(round(float(x), 4) for x in comp_str.split(","))

        sys_gs = gs_lookup.get(sys, {})

        eigs: list[float] = []
        margins: list[float] = []

        for t in t_grid:
            t_gs = sys_gs.get(t, {})
            if not t_gs:
                continue
            eig = _simplex_hessian_min_eig(comp, t_gs)
            if not np.isnan(eig):
                eigs.append(eig)
            fsm = _finite_sep_margin(comp, t_gs, min_sep=fsm_min_sep)
            if not np.isnan(fsm):
                margins.append(fsm)

        results[(sys, comp_str)] = {
            "hessian_min_eig_mean": float(np.mean(eigs)) if eigs else float("nan"),
            "hessian_min_eig_min":  float(np.min(eigs))  if eigs else float("nan"),
            "hessian_all_stable_frac": (
                float(sum(e > 0 for e in eigs) / len(eigs)) if eigs else float("nan")
            ),
            "n_valid_T_hessian": len(eigs),
            "ss_fsm_mean": float(np.mean(margins))  if margins else float("nan"),
            "ss_fsm_min":  float(np.min(margins))   if margins else float("nan"),
            "n_valid_T_fsm": len(margins),
        }

        if (k + 1) % 100 == 0:
            print(f"  {k + 1:>5} / {n_target} …", flush=True)

    print(f"  Done — {len(results)} compositions processed.")

    # Coverage diagnostics — make incomplete results visible.
    n_zero_hessian = sum(1 for v in results.values() if v["n_valid_T_hessian"] == 0)
    n_zero_fsm     = sum(1 for v in results.values() if v["n_valid_T_fsm"] == 0)
    n_partial_fsm  = sum(
        1 for v in results.values()
        if 0 < v["n_valid_T_fsm"] < n_t_expected
    )
    if n_zero_hessian:
        print(
            f"  Warning: {n_zero_hessian}/{n_target} compositions have zero valid "
            f"curvature T steps (all neighbours absent at grid boundary)."
        )
    if n_zero_fsm:
        print(
            f"  Warning: {n_zero_fsm}/{n_target} compositions have zero valid FSM "
            f"T steps (geometric infeasibility at all T for min_sep={fsm_min_sep:.2f} — "
            f"likely near-endmember compositions)."
        )
    if n_partial_fsm:
        print(
            f"  Info: {n_partial_fsm}/{n_target} compositions have partial FSM "
            f"coverage (< {n_t_expected} valid T steps). "
            f"Check n_valid_T_fsm before comparing ss_fsm_mean across compositions."
        )

    for (sys, comp_str), vals in results.items():
        mask = (out["system"] == sys) & (out["composition"] == comp_str)
        for col, val in vals.items():
            out.loc[mask, col] = val

    # Flag compositions where hessian was computed but is negative at every T step.
    # Globally stable (rank-1) but G_ss locally concave — stability maintained by
    # fixed-phase dominance of neighbouring SS compositions, not G_ss convexity.
    out["curvature_flagged"] = (
        (out["n_valid_T_hessian"] > 0) &
        (out["hessian_all_stable_frac"].fillna(1.0) == 0.0)
    )
    n_flagged = int(out["curvature_flagged"].sum())
    if n_flagged:
        print(
            f"  Curvature flag: {n_flagged} compositions have hessian_min_eig < 0 at "
            f"all valid T steps — marked curvature_flagged=True, excluded from shortlist."
        )

    return out


# ── Terminal output ───────────────────────────────────────────────────────────

def _print_top_n(metrics: pd.DataFrame, n: int = top_n) -> None:
    """Print two terminal tables: Pareto-ranked and candidate-ranked top-N.

    Table 1 — Thermodynamic (unchanged): top-N by primary Pareto rank, then
    strict_frac ↓, mean_dG_hull ↑, unique_frac ↓.

    Table 2 — Experimental candidates: top-N by global_candidate_rank
    (entropy-prioritized, eligibility-filtered).  Absent when no eligible
    candidates exist (global_candidate_rank column missing or all NaN).

    Args:
        metrics: DataFrame with ``pareto_rank`` column.
        n: Number of compositions to display per table.
    """
    # ── Table 1: Pareto ──────────────────────────────────────────────────────
    ranked = (
        metrics.dropna(subset=["pareto_rank"])
        .sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)
        .head(n)
        .reset_index(drop=True)
    )
    hdr = (
        f"{'#':>3}  {'System':<14} {'Composition':<30} "
        f"{'P-Rank':>7}  {'A-Rank':>7}  {'3-Rank':>7}  "
        f"{'dG_hull':>9}  {'strict':>7}  {'unique':>7}  "
        f"{'T_span(K)':>9}  {'!':>2}"
    )
    print(f"\nTop {n} compositions by primary Pareto rank (thermodynamic):")
    print(hdr)
    print("─" * len(hdr))
    for i, row in ranked.iterrows():
        label   = str(row.get("comp_label", row["composition"]))[:30]
        p_rank  = int(row["pareto_rank"])
        a_rank  = int(row["pareto_rank_alt"])  if pd.notna(row.get("pareto_rank_alt"))  else -1
        t3_rank = int(row["pareto_rank_3obj"]) if pd.notna(row.get("pareto_rank_3obj")) else -1
        flag    = str(row.get("strict_unique_flag", ""))
        t_span  = row.get("strict_T_span", float("nan"))
        interp  = "[i]" if row.get("interpolated", False) else "   "
        print(
            f"{i+1:>3}  {row['system']:<14} {label:<30} "
            f"{p_rank:>7}  {a_rank:>7}  {t3_rank:>7}  "
            f"{row['mean_dG_hull']:>9.4f}  {row['strict_frac']:>7.4f}  "
            f"{row['unique_frac']:>7.4f}  {t_span:>9.0f}  {flag:>2}  {interp}"
        )
    print()

    # ── Table 2: Experimental candidates (entropy-prioritized) ───────────────
    if "global_candidate_rank" not in metrics.columns:
        return
    cand = (
        metrics[metrics["global_candidate_rank"].notna()]
        .sort_values("global_candidate_rank")
        .head(n)
        .reset_index(drop=True)
    )
    if cand.empty:
        return
    hdr2 = (
        f"{'#':>3}  {'System':<14} {'Composition':<30} "
        f"{'G-Cand':>7}  {'s_norm':>7}  "
        f"{'FSM':>9}  {'Hessian':>9}  {'P-Rank':>7}"
    )
    print(f"Top {n} experimental candidates by global_candidate_rank")
    print("  (entropy-first; eligible rank-1 only — NOT a thermodynamic ranking):")
    print(hdr2)
    print("─" * len(hdr2))
    for i, row in cand.iterrows():
        label  = str(row.get("comp_label", row["composition"]))[:30]
        g_rank = int(row["global_candidate_rank"])
        s_norm = row.get("s_config_norm", float("nan"))
        fsm    = row.get("ss_fsm_mean",   float("nan"))
        hess   = row.get("hessian_min_eig_mean", float("nan"))
        p_rank = int(row["pareto_rank"]) if pd.notna(row.get("pareto_rank")) else -1
        print(
            f"{i+1:>3}  {row['system']:<14} {label:<30} "
            f"{g_rank:>7}  {s_norm:>7.4f}  "
            f"{fsm:>9.4f}  {hess:>9.4f}  {p_rank:>7}"
        )
    print()


# ── Analysis summary text ─────────────────────────────────────────────────────

def _write_summary_txt(
    metrics: pd.DataFrame,
    comparison: dict,
    issues: list[str],
    path: Path,
    sensitivity: dict | None = None,
    eps_sensitivity: dict | None = None,
    n_incomplete: int = 0,
) -> None:
    """Write a plain-text analysis summary file.

    Args:
        metrics: Metrics DataFrame with Pareto rank columns.
        comparison: Dict from ``_add_pareto_columns``.
        issues: Validation issue strings from ``_validate_metrics``.
        path: Output path.
        sensitivity: Optional dict from ``_sensitivity_comparison``.
        eps_sensitivity: Optional dict from ``_rank1_eps_sensitivity``.
        n_incomplete: Count of compositions excluded for insufficient T coverage.
    """
    ranked = (
        metrics.dropna(subset=["pareto_rank"])
        .sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)
    )
    front  = ranked[ranked["pareto_rank"] == 1]
    top20  = ranked.head(20)

    n_total_metrics   = len(metrics)
    n_interpolated    = int(metrics["interpolated"].fillna(False).sum())
    n_excl_sf         = int(comparison.get("n_excluded_strict_frac", 0))
    n_ranked          = len(ranked)

    # Entropy-prioritized candidate top-20.
    cand_top20: pd.DataFrame = pd.DataFrame()
    if "global_candidate_rank" in metrics.columns:
        cand_top20 = (
            metrics[metrics["global_candidate_rank"].notna()]
            .sort_values("global_candidate_rank")
            .head(20)
        )

    lines: list[str] = [
        "=" * 72,
        "  AlB₂-type Diboride Stability Pareto Analysis",
        f"  Temperature range: {_T_RANGE_LABEL}  ({_N_EXPECTED_T} steps, {_DB_T_STEP:g} K increments)",
        "=" * 72,
        "",
        "Composition accounting:",
        f"  Loaded (system, composition) pairs         : {n_total_metrics + n_incomplete:,}",
        f"  Excluded — insufficient T coverage         : {n_incomplete:,}",
        f"    └─ recovered via interpolation           : {n_interpolated:,}",
        f"  With sufficient T coverage (ranked + excl) : {n_total_metrics:,}",
        f"  Excluded — zero strict stability filter    : {n_excl_sf:,}",
        f"             (strict_frac < {rank_filter_min_strict_frac}; in 'Zero-Stability' sheet)",
        f"  Eligible for Pareto ranking                : {comparison['n_eligible']:,}",
        f"  Ranked (have pareto_rank)                  : {n_ranked:,}",
        "",
        "Pareto frontier sizes:",
        f"  Primary   (dG+strict)  : {comparison['n_primary']}",
        f"  Alternative (dG+ss_auc): {comparison['n_alt']}",
        f"  3-objective            : {comparison['n_3obj']}",
        f"  In both primary + alt  : {comparison['n_primary_and_alt']}",
        f"  Primary-only           : {comparison['n_only_primary']}",
        f"  Alt-only               : {comparison['n_only_alt']}",
        "",
        "Metric ranges (ranked compositions only):",
        f"  mean_dG_hull  : {ranked['mean_dG_hull'].min():.4f} – "
          f"{ranked['mean_dG_hull'].max():.4f} eV/f.u.",
        f"  strict_frac   : {ranked['strict_frac'].min():.4f} – "
          f"{ranked['strict_frac'].max():.4f}",
        f"  ss_auc        : {ranked['ss_auc'].min():.4f} – "
          f"{ranked['ss_auc'].max():.4f}  (normalized)",
        f"  unique_frac   : {ranked['unique_frac'].min():.4f} – "
          f"{ranked['unique_frac'].max():.4f}",
        "",
        "Primary front (rank = 1) best values:",
        f"  lowest mean_dG_hull  : {front['mean_dG_hull'].min():.4f} eV/f.u.  "
          f"({front.loc[front['mean_dG_hull'].idxmin(), 'system']} "
          f"{front.loc[front['mean_dG_hull'].idxmin(), 'comp_label']})",
        f"  highest strict_frac  : {front['strict_frac'].max():.4f}  "
          f"({front.loc[front['strict_frac'].idxmax(), 'system']} "
          f"{front.loc[front['strict_frac'].idxmax(), 'comp_label']})",
        f"  highest ss_auc       : {front['ss_auc'].max():.4f}  "
          f"({front.loc[front['ss_auc'].idxmax(), 'system']} "
          f"{front.loc[front['ss_auc'].idxmax(), 'comp_label']})",
        "",
        "Top 20 compositions (primary rank, then strict_frac↓, dG_hull↑, unique_frac↓):",
        f"  {'#':>2}  {'System':<14}  {'Composition':<30}  "
          f"{'P':>4}  {'dG':>8}  {'strict':>7}  {'unique':>7}  {'T_span':>7}  {'!':>2}",
    ]
    for i, (_, row) in enumerate(top20.iterrows(), 1):
        t_span = row.get("strict_T_span", float("nan"))
        flag   = str(row.get("strict_unique_flag", ""))
        lines.append(
            f"  {i:>2}. {row['system']:<14}  "
            f"{str(row.get('comp_label',''))[:30]:<30}  "
            f"{int(row['pareto_rank']):>4}  "
            f"{row['mean_dG_hull']:>8.4f}  "
            f"{row['strict_frac']:>7.4f}  "
            f"{row['unique_frac']:>7.4f}  "
            f"{t_span:>7.0f}  {flag:>2}"
        )

    # Entropy-prioritized experimental candidate top-20.
    if not cand_top20.empty:
        n_elig = int(metrics["global_candidate_rank"].notna().sum()) if "global_candidate_rank" in metrics.columns else 0
        lines += [
            "",
            f"Top 20 experimental candidates by global_candidate_rank ({n_elig} eligible):",
            "  (entropy-first within eligibility-filtered rank-1 — NOT a thermodynamic ranking)",
            f"  {'#':>2}  {'System':<14}  {'Composition':<30}  "
            f"{'G-Cand':>7}  {'s_norm':>7}  {'FSM':>8}  {'P-Rank':>7}",
        ]
        for i, (_, row) in enumerate(cand_top20.iterrows(), 1):
            g_rank = int(row["global_candidate_rank"])
            s_norm = row.get("s_config_norm", float("nan"))
            fsm    = row.get("ss_fsm_mean",   float("nan"))
            p_rank = int(row["pareto_rank"]) if pd.notna(row.get("pareto_rank")) else -1
            lines.append(
                f"  {i:>2}. {row['system']:<14}  "
                f"{str(row.get('comp_label',''))[:30]:<30}  "
                f"{g_rank:>7}  {s_norm:>7.4f}  {fsm:>8.4f}  {p_rank:>7}"
            )

    # Sensitivity analysis section.
    if sensitivity and sensitivity.get("n_interpolated_total", 0) > 0:
        n_disp = sensitivity.get("n_displaced_by_interp", 0)
        n_on   = sensitivity.get("n_interpolated_on_front", 0)
        lines += [
            "",
            "Sensitivity Analysis (interpolated compositions):",
            f"  Interpolated compositions            : {sensitivity['n_interpolated_total']}",
            f"  Interpolated comps on Pareto front   : {n_on}",
            f"  Non-interp comps displaced from front: {n_disp}",
            f"  Front size (complete-only)           : {sensitivity['n_front_complete_only']}",
            f"  Front size (with interpolated)       : {sensitivity['n_front_full']}",
            f"  Conclusion: {sensitivity['description']}",
        ]
        if n_on > 0:
            lines.append("  Interpolated comps on front:")
            for r in sensitivity.get("interpolated_on_front", []):
                lines.append(f"    {r.get('system','')} {r.get('comp_label','')}")
        if n_disp > 0:
            lines.append("  Displaced non-interpolated comps:")
            for r in sensitivity.get("displaced_compositions", []):
                lines.append(f"    {r.get('system','')} {r.get('comp_label','')}")
        lines += [
            "  Note: Interpolated dG_above_hull and ss_fraction values are linearly",
            "  estimated from neighboring calculated T steps.  Boolean stability flags",
            "  (is_single_phase_strict, is_unique_homogeneous) are not interpolated.",
            "  strict_frac_lb / unique_frac_lb give conservative lower bounds treating",
            "  all missing steps (including interpolated gaps) as not-stable.",
            "  Interpolated compositions are flagged in 'Overall Rankings' sheet.",
            "  Interpolated compositions are EXCLUDED from the high-confidence shortlist.",
        ]

    if eps_sensitivity and "description" in eps_sensitivity:
        mem = eps_sensitivity.get("membership", {})
        lines += [
            "",
            "Numerical Tolerance Sensitivity (rank-1 count and membership vs pareto_eps):",
            f"  eps=0 (exact)   : {eps_sensitivity.get('eps_0', 'n/a')} compositions",
            f"  eps=1e-9        : {eps_sensitivity.get('eps_1e-9', 'n/a')} compositions",
            f"  eps=1e-6 (used) : {eps_sensitivity.get('eps_1e-6', 'n/a')} compositions",
            f"  eps=1e-3        : {eps_sensitivity.get('eps_1e-3', 'n/a')} compositions",
        ]
        for eps_val in [0.0, 1e-9, 1e-3]:
            m = mem.get(eps_val, {})
            j = m.get("jaccard_vs_current", float("nan"))
            only_this = m.get("only_in_this", 0)
            only_curr = m.get("only_in_curr", 0)
            lines.append(
                f"  Jaccard(eps={eps_val:.0e} vs current): {j:.4f}  "
                f"(+{only_this} unique to this eps, -{only_curr} present only in current)"
            )
        lines += [
            f"  Assessment: {eps_sensitivity['description']}",
            "  A large jump at eps=1e-6 vs eps=0 means LP noise is inflating the front.",
            "  A large jump at eps=1e-3 vs eps=1e-6 means physical gradients are present.",
        ]

    lines += ["", "Validation:"]
    if issues:
        for iss in issues:
            lines.append("  [WARN] " + iss)
    else:
        lines.append("  All checks passed.")

    lines += [
        "",
        "Methodology:",
        "  Primary Pareto objectives (pareto_rank):",
        "    1. mean_dG_hull = mean(dG_above_hull over valid T)  [eV/f.u.]  — minimize",
        "       Mean energy above the LP-computed convex hull across all temperature steps.",
        "       Lower values indicate the target composition is closer to thermodynamic",
        "       equilibrium with respect to competing phases.",
        "    2. strict_frac = sum(is_single_phase_strict) / n_valid_T          — maximize",
        "       Fraction of T steps where ss_fraction >= 0.95 AND dG_above_hull <= 1 meV.",
        "       Requires the target composition to lie on the hull as a single phase.",
        "       Excludes phase-separated two-AlB₂ assemblages.",
        "",
        "  Alternative Pareto frontier (pareto_rank_alt):",
        "    1. mean_dG_hull  — minimize (same as above)",
        "    2. ss_auc = trapz(ss_fraction, T) / (T_max - T_min)             — maximize",
        "       Normalized AlB₂ parent-structure retention AUC.  Counts all AlB₂-type",
        "       phase weight including phase-separated assemblages.  Less strict than",
        "       strict_frac but captures near-miss compositions.",
        "",
        "  Supplemental 3-objective frontier (pareto_rank_3obj):",
        "    Minimizes mean_dG_hull, maximizes strict_frac AND ss_auc simultaneously.",
        "",
        "  Reported metrics (not Pareto objectives):",
        "    unique_frac  = sum(is_unique_homogeneous) / n_valid_T",
        "    strict_T_span = longest continuous T window (K) with is_single_phase_strict",
        "    max_dG_hull, mean_G_ss_formation, T_min_retained",
        "",
        "  Secondary ordering within each Pareto rank:",
        "    strict_frac ↓  →  mean_dG_hull ↑  →  unique_frac ↓",
        "",
        "  Stability flag '!' marks compositions with strict_frac >= 0.5 but",
        "  unique_frac < 0.5 × strict_frac.  The LP-assigned homogeneous equilibrium",
        "  fraction is substantially lower than the strict-stability fraction,",
        "  suggesting the equilibrium assemblage may contain multiple AlB₂ phases.",
        "",
        f"  T range          : {int(t_norm_min)} – {int(t_norm_max)} K",
        f"  min_valid_frac   : {min_valid_frac}",
        f"  pareto_eps       : {pareto_eps}",
        "  ss_threshold     : 0.95  (defined in build_stability_db.py)",
        "  dG_hull_tol      : 0.001 eV/f.u.  (defined in build_stability_db.py)",
        "",
        "Limitations:",
        "  Thermodynamic stability is relative to competing phases included in the upstream",
        "  thermodynamic database only.  Stability claims do not extend to phases not",
        "  included in build_stability_db.py.",
        "  ss_fraction counts all AlB₂-type phases, including phase-separated assemblages.",
        "  is_unique_homogeneous describes the LP solution at a sampled grid point, not a",
        "    mathematically unique equilibrium over the full composition space.",
        "  Thermodynamic stability does not guarantee experimental synthesizability.",
        "=" * 72,
    ]

    path.write_text("\n".join(lines))
    print(f"Summary → {path}")


# ── Excel workbook ────────────────────────────────────────────────────────────

def _autoformat_sheet(ws, df: pd.DataFrame) -> None:
    """Apply auto-width, frozen header, and column filters to a worksheet.

    Args:
        ws: openpyxl Worksheet.
        df: DataFrame that was written to the sheet.
    """
    from openpyxl.utils import get_column_letter

    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions

    for col_idx, col_name in enumerate(df.columns, start=1):
        col_letter = get_column_letter(col_idx)
        max_len    = max(len(str(col_name)), 8)
        for cell in ws[col_letter][1:]:
            try:
                max_len = max(max_len, len(str(cell.value)) if cell.value is not None else 0)
            except TypeError:
                pass
        ws.column_dimensions[col_letter].width = min(max_len + 2, 40)


def _apply_rank_colorscale(ws, df: pd.DataFrame) -> None:
    """Apply a green→red color scale on the pareto_rank column if present.

    Args:
        ws: openpyxl Worksheet.
        df: DataFrame written to the sheet.
    """
    from openpyxl.formatting.rule import ColorScaleRule
    from openpyxl.utils import get_column_letter

    if "pareto_rank" not in df.columns:
        return
    rank_col   = df.columns.get_loc("pareto_rank") + 1
    col_letter = get_column_letter(rank_col)
    n_rows     = len(df)
    ws.conditional_formatting.add(
        f"{col_letter}2:{col_letter}{n_rows + 1}",
        ColorScaleRule(
            start_type="min", start_color="63BE7B",
            mid_type="percentile", mid_value=50, mid_color="FFEB84",
            end_type="max", end_color="F8696B",
        ),
    )


def _assign_candidate_ranks(metrics: pd.DataFrame) -> pd.DataFrame:
    """Add entropy and candidate-rank columns to the metrics DataFrame.

    New columns (all rows):
    - ``s_config``      J/mol/K = −R·Σ xᵢ·ln(xᵢ)
    - ``s_config_norm`` ∈ [0,1] = −Σ xᵢ·ln(xᵢ) / ln(n_metals).  Equals 1 for
      perfectly equimolar, 0 for a pure endmember.  Comparable across systems
      with different numbers of metals.

    New columns (eligible rank-1 only; NaN for ineligible):
    - ``candidate_rank``        within-system 1-based rank
    - ``global_candidate_rank`` cross-system 1-based rank

    **Eligibility** (same as experimental shortlist):
    pareto_rank=1, strict_frac=1, unique_frac=1, |mean_dG_hull|<1 μeV,
    not interpolated, curvature_flagged=False.

    **Within-system ordering** (``candidate_rank``):
    1. ``s_config_norm`` (desc) — PRIMARY; equimolarity distinguishes candidates
       that satisfy the same thermodynamic eligibility requirements.
    2. ``hessian_min_eig_mean`` binned to ``candidate_hessian_tol`` (desc) —
       tiebreaker; Hessian magnitude is coordinate-dependent but comparable
       within one system (all compositions share the same n_metals).
    3. ``ss_fsm_mean`` binned to ``candidate_fsm_tol`` (desc) — tiebreaker.
    NaN supplementary values rank after all valid values within each key.

    **Global ordering** (``global_candidate_rank``):
    Hessian eigenvalue magnitudes are NOT comparable across systems with
    different n_metals (composition simplex dimension differs); they are
    excluded from the global sort.
    1. ``s_config_norm`` (desc) — PRIMARY; a 50/50 binary ranks above a 90/10
       binary regardless of their Hessian or FSM values.
    2. ``ss_fsm_mean`` binned to ``candidate_fsm_tol`` (desc) — tiebreaker
       only; does not override the entropy preference.
    NaN FSM ranks after all valid values within the same s_config_norm.

    Both rankings are compositional/robustness-based only and do **not**
    alter or re-interpret the thermodynamic Pareto rank.

    Args:
        metrics: Per-composition metrics DataFrame after Pareto ranking and
            supplementary-metric computation.

    Returns:
        New DataFrame with the four new columns added.
    """
    _R = 8.314472  # J/(mol·K)
    out = metrics.copy()

    # ── Configurational entropy for every row ────────────────────────────────
    x_cols = sorted(c for c in out.columns if c.startswith("x_"))
    n_metals_col = out["n_metals"] if "n_metals" in out.columns else None

    def _entropy(row: pd.Series) -> tuple[float, float]:
        fracs = [row[c] for c in x_cols if pd.notna(row[c]) and row[c] > 1e-9]
        if not fracs:
            return float("nan"), float("nan")
        s = -sum(f * np.log(f) for f in fracs)   # s > 0
        n = len(fracs) if n_metals_col is None else int(row["n_metals"]) if pd.notna(row.get("n_metals")) else len(fracs)
        s_norm = s / np.log(n) if n > 1 else 0.0
        return float(_R * s), float(s_norm)

    if x_cols:
        # Vectorized _entropy: same result for every row.
        X = out[x_cols].to_numpy(dtype=float)
        present = np.isfinite(X) & (X > 1e-9)
        with np.errstate(divide="ignore", invalid="ignore"):
            s = -np.where(present, X * np.log(np.where(present, X, 1.0)), 0.0).sum(axis=1)
        n_present = present.sum(axis=1)
        if n_metals_col is not None:
            nm = out["n_metals"].to_numpy(dtype=float)
            n = np.where(np.isfinite(nm), nm, n_present)
        else:
            n = n_present.astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            s_norm = np.where(n > 1, s / np.log(np.maximum(n, 2)), 0.0)
        out["s_config"]      = np.where(n_present > 0, _R * s, np.nan)
        out["s_config_norm"] = np.where(n_present > 0, s_norm, np.nan)
    else:
        out["s_config"]      = float("nan")
        out["s_config_norm"] = float("nan")

    # ── Eligibility filter ───────────────────────────────────────────────────
    ranked = out.dropna(subset=["pareto_rank"])
    _curvature_ok = (
        ~ranked["curvature_flagged"].fillna(False)
        if "curvature_flagged" in ranked.columns
        else pd.Series(True, index=ranked.index)
    )
    eligible_mask = (
        (ranked["pareto_rank"] == 1) &
        (ranked["strict_frac"] == 1.0) &
        (ranked["unique_frac"] == 1.0) &
        (ranked["mean_dG_hull"].abs() < 1e-6) &
        (~ranked["interpolated"].fillna(False)) &
        _curvature_ok
    )
    eligible_idx = ranked.index[eligible_mask]

    out["candidate_rank"]        = pd.array([pd.NA] * len(out), dtype=pd.Int64Dtype())
    out["global_candidate_rank"] = pd.array([pd.NA] * len(out), dtype=pd.Int64Dtype())

    if len(eligible_idx) == 0:
        return out

    elig = out.loc[eligible_idx].copy()

    # Binned columns for near-tie detection.
    def _bin_col(src_col: str, tol: float) -> pd.Series:
        v = elig[src_col] if src_col in elig.columns else pd.Series(float("nan"), index=elig.index)
        return (v / tol).round() * tol

    elig["_h_bin"] = _bin_col("hessian_min_eig_mean", candidate_hessian_tol)
    elig["_f_bin"] = _bin_col("ss_fsm_mean",           candidate_fsm_tol)

    # ── Within-system candidate_rank ─────────────────────────────────────────
    # s_config_norm (primary) → hessian bin → FSM bin
    within_sort = ["s_config_norm", "_h_bin", "_f_bin"]
    within_asc  = [False,           False,    False]
    candidate_ranks: dict[int, int] = {}
    for _sys, grp in elig.groupby("system", sort=False):
        sorted_grp = grp.sort_values(within_sort, ascending=within_asc, na_position="last")
        for rank_i, idx in enumerate(sorted_grp.index, start=1):
            candidate_ranks[idx] = rank_i

    # ── Global global_candidate_rank ─────────────────────────────────────────
    # s_config_norm (primary) → FSM bin tiebreaker; no Hessian (dimension-dependent)
    global_sort = ["s_config_norm", "_f_bin"]
    global_asc  = [False,           False]
    global_sorted = elig.sort_values(global_sort, ascending=global_asc, na_position="last")
    global_ranks: dict[int, int] = {
        idx: rank_i
        for rank_i, idx in enumerate(global_sorted.index, start=1)
    }

    for idx in eligible_idx:
        out.at[idx, "candidate_rank"]        = candidate_ranks[idx]
        out.at[idx, "global_candidate_rank"] = global_ranks[idx]

    # ── Validation print ─────────────────────────────────────────────────────
    if verbose_ranking:
        top_global = (
            elig
            .assign(_gr=elig.index.map(global_ranks))
            .sort_values("_gr")
            .head(8)
        )
        print("\n  Candidate ranking validation — top 8 globally:")
        for _, row in top_global.iterrows():
            print(
                f"    global#{int(row['_gr']):<4} {row.get('system','?'):<20} "
                f"comp={row.get('composition','?'):<20} "
                f"s_norm={row.get('s_config_norm', float('nan')):.3f}  "
                f"fsm_bin={row.get('_f_bin', float('nan')):.4f}  "
                f"hess_bin={row.get('_h_bin', float('nan')):.4f}"
            )

    return out


def _build_experimental_shortlist(metrics: pd.DataFrame) -> pd.DataFrame:
    """Select the top-``shortlist_per_system`` candidates per system from the
    high-confidence rank-1 eligible set, using pre-computed ``candidate_rank``.

    Call ``_assign_candidate_ranks`` before this function so that ``metrics``
    already carries ``candidate_rank``, ``s_config``, and ``s_config_norm``.

    **Eligibility criteria** (applied by ``_assign_candidate_ranks`` before ranking):
    - ``pareto_rank == 1``
    - ``strict_frac == 1.0``    — strictly stable at every observed T step
    - ``unique_frac == 1.0``    — uniquely homogeneous at every observed T step
    - ``|mean_dG_hull| < 1e-6`` — ΔG within LP numerical noise
    - ``interpolated == False`` — all T steps fully observed
    - ``curvature_flagged == False`` (or absent)

    Eligibility is applied BEFORE entropy-based ordering so that higher entropy
    can never allow an ineligible composition to outrank an eligible one.

    Ordering within each system (by ``candidate_rank``, ascending):
    ``s_config_norm``↓ → Hessian bin↓ → FSM bin↓ (NaN last for Hessian/FSM).
    This is an ENTROPY-PRIORITIZED SELECTION RANKING, not a thermodynamic rank.

    Args:
        metrics: Per-composition metrics DataFrame after ``_assign_candidate_ranks``.

    Returns:
        Shortlist DataFrame (top ``shortlist_per_system`` per system).
    """
    if "candidate_rank" not in metrics.columns:
        return pd.DataFrame()

    eligible = metrics[metrics["candidate_rank"].notna()].copy()
    if eligible.empty:
        return pd.DataFrame()

    ranked = metrics.dropna(subset=["pareto_rank"])
    n_rank1        = int((ranked["pareto_rank"] == 1).sum())
    n_hc           = len(eligible)
    n_curv_flagged = int(
        ((ranked["pareto_rank"] == 1) &
         ranked.get("curvature_flagged", pd.Series(False, index=ranked.index)).fillna(False)
        ).sum()
    ) if "curvature_flagged" in ranked.columns else 0

    eligible["shortlist_group"] = "high_confidence"
    eligible["shortlist_note"] = (
        f"High-confidence shortlist: verified pareto_rank=1, strict_frac=1.0, "
        f"unique_frac=1.0 (is_unique_homogeneous at ALL observed T steps — "
        f"primary equilibrium stability screen), mean_dG_hull <1 μeV (LP noise), "
        f"all T steps fully observed (not interpolated), curvature_flagged=False.  "
        f"{n_hc} of {n_rank1} rank-1 compositions eligible.  "
        f"Excluded: {n_rank1 - n_hc} total — "
        f"unique_frac/dG/interpolation filters and {n_curv_flagged} with "
        f"curvature_flagged=True (G_ss locally concave at all T).  "
        "Within-system ordering by candidate_rank: s_config_norm↓ → Hessian bin↓ → FSM bin↓ "
        "(NaN last for Hessian/FSM).  Entropy is the primary criterion — higher normalized "
        "configurational entropy preferred among candidates passing the same stability filters.  "
        "Hessian/FSM are tiebreakers only; bin widths: Hessian 1 meV, FSM 0.1 meV.  "
        "Missing Hessian/FSM values marked 'unavailable', not unstable."
    )

    shortlist = (
        eligible
        .sort_values(["n_metals", "system", "candidate_rank"], ascending=[False, True, True])
        .groupby("system", sort=False)
        .head(shortlist_per_system)
        .sort_values(["n_metals", "system", "candidate_rank"], ascending=[False, True, True])
        .reset_index(drop=True)
    )
    return shortlist


def _build_excel(
    metrics: pd.DataFrame,
    raw: pd.DataFrame,
    summary: pd.DataFrame,
    comparison: dict,
    incompl_df: pd.DataFrame,
    path: Path,
) -> None:
    """Write a 7-sheet Excel workbook of Pareto analysis results.

    Sheets: Overall Rankings, Pareto Front, Top Candidates, Temperature Data,
    System Summary, Methodology, Excluded Compositions.

    Args:
        metrics: Per-composition metrics with Pareto rank columns.
        raw: Raw stability DataFrame.
        summary: Summary stability DataFrame.
        comparison: Dict from ``_add_pareto_columns``.
        incompl_df: DataFrame of excluded compositions.
        path: Output .xlsx path.
    """
    ranked_all = (
        metrics.dropna(subset=["pareto_rank"])
        .sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)
        .reset_index(drop=True)
    )

    # Top Candidates: scan global_candidate_rank ascending, enforce top_per_system
    # cap while building the global top_n block; then fill per-system slots to the
    # same cap.  Falls back to Pareto sort if global_candidate_rank is absent.
    _has_cand_rank = "global_candidate_rank" in metrics.columns
    if _has_cand_rank:
        _cand_pool = metrics[metrics["global_candidate_rank"].notna()].copy()
        _pool_sorted = _cand_pool.sort_values("global_candidate_rank")
    else:
        _pool_sorted = ranked_all.copy()

    seen: set[tuple] = set()
    sys_counts: dict[str, int] = {}     # total selected per system
    top_overall_rows: list[pd.Series] = []

    for _, row in _pool_sorted.iterrows():
        if len(top_overall_rows) >= top_n:
            break
        sys_ = row["system"]
        if sys_counts.get(sys_, 0) >= top_per_system:
            continue                    # system already at cap — skip
        key = (sys_, row["composition"])
        top_overall_rows.append(row)
        seen.add(key)
        sys_counts[sys_] = sys_counts.get(sys_, 0) + 1

    top_overall = pd.DataFrame(top_overall_rows)

    sys_extra: list[pd.Series] = []
    # Fill remaining per-system slots (systems with fewer than top_per_system
    # selections that weren't reached by the global top_n scan).
    for _sys, sg in (
        (_cand_pool if _has_cand_rank else ranked_all).groupby("system", sort=False)
    ):
        already = sys_counts.get(_sys, 0)
        slots   = top_per_system - already
        if slots <= 0:
            continue
        sg_sorted = (
            sg.sort_values("candidate_rank")
            if _has_cand_rank and "candidate_rank" in sg.columns
            else sg.sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)
        )
        added = 0
        for _, row in sg_sorted.iterrows():
            if added >= slots:
                break
            key = (row["system"], row["composition"])
            if key not in seen:
                sys_extra.append(row)
                seen.add(key)
                added += 1

    _extra_df = pd.DataFrame(sys_extra)
    top_candidates = pd.concat(
        [top_overall, _extra_df] if not _extra_df.empty else [top_overall],
        ignore_index=True,
    )
    # Final sort: global candidates first, extras after (by candidate_rank within system).
    if _has_cand_rank:
        top_candidates = top_candidates.sort_values(
            ["global_candidate_rank", "system", "candidate_rank"],
            ascending=[True, True, True],
            na_position="last",
        )
    else:
        top_candidates = top_candidates.sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)

    # Pareto Front sheet: primary rank = 1 with all three frontier flags.
    pareto_front = ranked_all[ranked_all["pareto_rank"] == 1].copy()
    pareto_front["in_primary_front"] = True
    pareto_front["in_alt_front"]     = pareto_front["pareto_rank_alt"]  == 1
    pareto_front["in_3obj_front"]    = pareto_front["pareto_rank_3obj"] == 1

    # Temperature Data: top_n compositions only.
    top_keys  = top_candidates[["system", "composition"]].assign(_k=True)
    temp_data = (
        raw.merge(top_keys, on=["system", "composition"], how="inner")
        .drop(columns=["_k"])
        .sort_values(["system", "composition", "T"])
    )

    # System Summary.
    sys_sum = (
        metrics.dropna(subset=["pareto_rank"])
        .groupby("system", sort=False)
        .agg(
            n_compositions          = ("composition",      "count"),
            n_primary_front         = ("pareto_rank",      lambda x: int((x == 1).sum())),
            n_alt_front             = ("pareto_rank_alt",  lambda x: int((x == 1).sum())),
            n_3obj_front            = ("pareto_rank_3obj", lambda x: int((x == 1).sum())),
            best_primary_rank       = ("pareto_rank",      "min"),
            mean_dG_hull_best       = ("mean_dG_hull",     "min"),
            strict_frac_best        = ("strict_frac",      "max"),
            ss_auc_best             = ("ss_auc",           "max"),
            unique_frac_best        = ("unique_frac",      "max"),
            strict_T_span_best      = ("strict_T_span",    "max"),
        )
        .reset_index()
        .sort_values("best_primary_rank")
    )

    # Methodology sheet.
    method_df = pd.DataFrame({
        "Parameter": [
            "Script",
            "Primary frontier objectives",
            "Alternative frontier objectives",
            "3-objective frontier objectives",
            "Secondary sort within rank",
            "strict_unique_flag condition",
            "Analysis temperature range (K)",
            "min_valid_frac",
            "pareto_eps",
            "ss_threshold (from build_stability_db.py)",
            "dG_hull_tol (from build_stability_db.py)",
            "Data source",
            "Frontier comparison",
            "Interpolation",
            "max_interpolation_gap_K",
            "discontinuity_threshold",
            "Candidate ranking — eligibility",
            "Candidate ranking — within-system (candidate_rank)",
            "Candidate ranking — global (global_candidate_rank)",
            "Candidate ranking — note",
        ],
        "Value": [
            "analyze_stability_pareto.py",
            "minimize mean_dG_hull | maximize strict_frac",
            "minimize mean_dG_hull | maximize ss_auc",
            "minimize mean_dG_hull | maximize strict_frac | maximize ss_auc",
            "strict_frac↓, mean_dG_hull↑, unique_frac↓",
            "strict_frac >= 0.5 AND unique_frac < 0.5 × strict_frac",
            _T_RANGE_LABEL,
            str(min_valid_frac),
            str(pareto_eps),
            "0.95",
            "0.001 eV/f.u.",
            str((db_dir / "stability_raw.parquet").resolve()),
            comparison["description"],
            (
                "Linear interpolation of dG_above_hull, ss_fraction, G_ss_formation "
                "at interior missing T steps.  Boolean flags NOT interpolated.  "
                "strict_frac_lb / unique_frac_lb are conservative lower bounds."
                if interpolate_incomplete else "Disabled"
            ),
            str(max_interpolation_gap_K),
            str(discontinuity_threshold),
            (
                "pareto_rank=1, strict_frac=1.0, unique_frac=1.0, |mean_dG_hull|<1 μeV, "
                "not interpolated, curvature_flagged=False.  "
                "Eligibility applied BEFORE entropy ordering."
            ),
            (
                "s_config_norm↓ (primary) → Hessian bin↓ → FSM bin↓ (NaN last).  "
                f"Hessian bin width: {candidate_hessian_tol*1e3:.0f} meV; "
                f"FSM bin width: {candidate_fsm_tol*1e3:.1f} meV."
            ),
            (
                "s_config_norm↓ (primary) → FSM bin↓ (NaN last).  "
                "Hessian excluded: eigenvalue magnitudes are not comparable across "
                "systems with different n_metals."
            ),
            (
                "candidate_rank and global_candidate_rank are ENTROPY-PRIORITIZED "
                "SELECTION RANKINGS, not thermodynamic stability ranks.  "
                "They apply only to compositions already passing all stability filters.  "
                "s_config_norm = -Σ xᵢ·ln(xᵢ) / ln(n_metals) ∈ [0,1]; 1 = equimolar."
            ),
        ],
    })

    shortlist = _build_experimental_shortlist(metrics)

    # Supplementary Metrics sheet: rank-1 compositions with curvature + FSM data.
    # ss_fsm_* columns measure AlB₂ SS-only decomposition margin; they do NOT
    # include non-AlB₂ competing phases and do NOT imply equilibrium stability.
    # n_valid_T_* report how many temperature steps had valid calculations.
    supp_cols_excel = [
        "system", "composition", "comp_label", "n_metals",
        "pareto_rank", "candidate_rank", "global_candidate_rank",
        "s_config", "s_config_norm",
        "strict_frac", "unique_frac", "mean_dG_hull",
        "hessian_min_eig_mean", "hessian_min_eig_min", "hessian_all_stable_frac",
        "n_valid_T_hessian", "curvature_flagged",
        "ss_fsm_mean", "ss_fsm_min", "n_valid_T_fsm",
    ]
    supp_df = metrics[
        [c for c in supp_cols_excel if c in metrics.columns]
    ].dropna(
        subset=["hessian_min_eig_mean", "ss_fsm_mean"], how="all"
    ).sort_values(
        ["pareto_rank", "ss_fsm_mean"], ascending=[True, False]
    ).reset_index(drop=True)

    # Candidate Rankings sheet: ALL eligible rank-1 compositions with candidate_rank.
    # This is the full within-system ordering, of which the shortlist is the top-N subset.
    cand_cols = [
        "global_candidate_rank", "system", "composition", "comp_label", "n_metals",
        "candidate_rank", "s_config_norm", "s_config",
        "hessian_min_eig_mean", "hessian_min_eig_min", "hessian_all_stable_frac",
        "n_valid_T_hessian", "curvature_flagged",
        "ss_fsm_mean", "ss_fsm_min", "n_valid_T_fsm",
        "pareto_rank", "strict_frac", "unique_frac", "mean_dG_hull",
    ]
    cand_df = pd.DataFrame()
    if "global_candidate_rank" in metrics.columns:
        cand_df = (
            metrics[metrics["global_candidate_rank"].notna()]
            [[c for c in cand_cols if c in metrics.columns]]
            .sort_values("global_candidate_rank")  # interleaved cross-system by global rank
            .reset_index(drop=True)
        )

    # Compositions excluded by the strict-stability eligibility filter.
    # They are present in the metrics table but have no pareto_rank.
    zero_strict = pd.DataFrame()
    if rank_filter_min_strict_frac is not None:
        _zero_mask = (
            metrics["strict_frac"].notna() &
            (metrics["strict_frac"] < rank_filter_min_strict_frac) &
            metrics["pareto_rank"].isna()
        )
        zero_strict = (
            metrics[_zero_mask]
            .sort_values(["system", "mean_dG_hull"])
            .reset_index(drop=True)
        )

    # Write to a local temp file first, then move to the final path.
    # Avoids iCloud / network-filesystem write timeouts: openpyxl flushes a
    # large ZIP archive in one shot; writing to /tmp is local, then one fast
    # rename delivers it to the iCloud-synced destination.
    import shutil as _shutil
    import tempfile as _tempfile
    _tmp_fd, _tmp_path_str = _tempfile.mkstemp(suffix=".xlsx", dir="/tmp")
    import os as _os; _os.close(_tmp_fd)
    _excel_path = Path(_tmp_path_str)

    try:
        with pd.ExcelWriter(_excel_path, engine="openpyxl") as writer:
            sheets: list[tuple[str, pd.DataFrame]] = [
                ("Overall Rankings",      ranked_all),
                ("Pareto Front",          pareto_front),
                ("Top Candidates",        top_candidates),
                ("Experimental Shortlist", shortlist),
                ("Candidate Rankings",    cand_df),
                ("Supplementary Metrics", supp_df),
                ("Temperature Data",      temp_data),
                ("System Summary",        sys_sum),
                ("Methodology",           method_df),
                ("Zero-Stability",        zero_strict),
                ("Excluded Compositions", incompl_df),
            ]
            for sheet_name, df in sheets:
                if len(df) > excel_max_rows:
                    print(
                        f"  Excel '{sheet_name}': {len(df):,} rows truncated to {excel_max_rows:,} "
                        "(full data in pareto_metrics.parquet)"
                    )
            sheets = [(name, df.head(excel_max_rows)) for name, df in sheets]
            for sheet_name, df in sheets:
                df.to_excel(writer, sheet_name=sheet_name, index=False)

            wb = writer.book
            for sheet_name, df in sheets:
                ws = wb[sheet_name]
                _autoformat_sheet(ws, df)
                if sheet_name in ("Overall Rankings", "Pareto Front", "Top Candidates"):
                    _apply_rank_colorscale(ws, df)

        _shutil.move(str(_excel_path), str(path))
    except Exception:
        _excel_path.unlink(missing_ok=True)
        raise

    print(f"Excel workbook → {path}")


# ── Plots ─────────────────────────────────────────────────────────────────────

def _plot_2d_static(metrics: pd.DataFrame, plot_dir: Path) -> None:
    """Save the primary 2D Pareto scatter as PNG and PDF.

    X = strict_frac (primary objective, maximize),
    Y = mean_dG_hull (primary objective, minimize),
    color = unique_frac.
    Primary Pareto front (rank 1) circled in crimson.
    Top 15 compositions labelled.

    Args:
        metrics: DataFrame with Pareto rank and objective columns.
        plot_dir: Output directory.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ranked = metrics.dropna(
        subset=["pareto_rank", "mean_dG_hull", "strict_frac", "unique_frac"]
    ).sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)
    if ranked.empty:
        print("No valid data for 2D static plot; skipping.")
        return

    fig, ax = plt.subplots(figsize=(10, 7.5))

    sc = ax.scatter(
        ranked["strict_frac"],
        ranked["mean_dG_hull"],
        c=ranked["unique_frac"],
        cmap="plasma",
        s=20,
        alpha=0.55,
        linewidths=0,
        vmin=0,
        vmax=1,
        zorder=3,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("unique_frac (homogeneous equilibrium)", fontsize=10)

    front = ranked[ranked["pareto_rank"] == 1]
    ax.scatter(
        front["strict_frac"],
        front["mean_dG_hull"],
        s=80, c="none",
        edgecolors="crimson", linewidths=1.6,
        zorder=6, label="Primary Pareto front (rank 1)",
    )

    label_top = ranked.head(15)
    for _, row in label_top.iterrows():
        ax.annotate(
            row["system"],
            xy=(row["strict_frac"], row["mean_dG_hull"]),
            xytext=(5, 3), textcoords="offset points",
            fontsize=6.5, color="#222222", alpha=0.85,
            zorder=7,
        )

    ax.set_xlabel("Strict single-phase fraction", fontsize=12)
    ax.set_ylabel("Mean ΔG above hull (eV / f.u.)", fontsize=12)
    ax.set_title(
        f"AlB₂ diboride stability — primary Pareto front (strict_frac vs mean_dG_hull)  [{_T_RANGE_LABEL}]",
        fontsize=12,
    )
    ax.legend(fontsize=9)
    plt.tight_layout()

    for ext in ("png", "pdf"):
        p = plot_dir / f"pareto_2D.{ext}"
        fig.savefig(p, dpi=150)
        print(f"2D static → {p}")
    plt.close(fig)


def _plot_interactive_2d(metrics: pd.DataFrame, plot_dir: Path) -> None:
    """Save a self-contained 2D interactive Pareto scatter as HTML.

    X = strict_frac, Y = mean_dG_hull, color = unique_frac.
    All three Pareto frontiers shown as toggleable overlays.

    Args:
        metrics: DataFrame with Pareto rank and objective columns.
        plot_dir: Output directory.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("plotly not installed; skipping 2D interactive plot.")
        return

    ranked = metrics.dropna(
        subset=["pareto_rank", "mean_dG_hull", "strict_frac", "unique_frac"]
    ).copy().sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)
    if ranked.empty:
        return

    def _hover(row: pd.Series) -> str:
        metals  = str(row.get("metals", ""))
        x_parts = [
            f"  x_{m}={row[f'x_{m}']:.3f}"
            for m in metals.split(",")
            if f"x_{m}" in row.index and pd.notna(row.get(f"x_{m}"))
        ]
        fracs   = "<br>".join(x_parts) if x_parts else ""
        t_span  = row.get("strict_T_span", float("nan"))
        p_rank  = int(row["pareto_rank"])
        a_rank  = int(row["pareto_rank_alt"])  if pd.notna(row.get("pareto_rank_alt"))  else -1
        t3_rank = int(row["pareto_rank_3obj"]) if pd.notna(row.get("pareto_rank_3obj")) else -1
        flag    = str(row.get("strict_unique_flag", ""))
        return (
            f"<b>{row['system']}</b><br>"
            f"{row.get('comp_label', row['composition'])}<br>"
            f"{fracs}<br>"
            f"mean_dG_hull: {row['mean_dG_hull']:.4f} eV/f.u.<br>"
            f"strict_frac: {row['strict_frac']:.4f}<br>"
            f"ss_auc: {row.get('ss_auc', float('nan')):.4f}<br>"
            f"unique_frac: {row.get('unique_frac', float('nan')):.4f}<br>"
            f"strict_T_span: {t_span:.0f} K<br>"
            f"Primary rank: {p_rank}<br>"
            f"Alt rank: {a_rank}<br>"
            f"3-obj rank: {t3_rank}<br>"
            f"System rank: {int(row['system_pareto_rank']) if pd.notna(row.get('system_pareto_rank')) else -1}"
            + ("<br><b>! high strict / low unique</b>" if flag == "!" else "")
        )

    hover_text = [_hover(r) for _, r in ranked.iterrows()]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=ranked["strict_frac"].tolist(),
        y=ranked["mean_dG_hull"].tolist(),
        mode="markers",
        marker={
            "color":     ranked["unique_frac"].tolist(),
            "colorscale": "Plasma",
            "cmin": 0, "cmax": 1,
            "size": 6,
            "opacity": 0.7,
            "colorbar": {"title": "unique_frac"},
        },
        text=hover_text,
        hoverinfo="text",
        name="All compositions",
    ))

    # Frontier overlays.
    front_p  = ranked[ranked["pareto_rank"] == 1]
    front_a  = ranked[ranked["pareto_rank_alt"] == 1]
    front_3  = ranked[ranked["pareto_rank_3obj"] == 1]

    def _safe_hover_idx(subset: pd.DataFrame) -> list[str]:
        return [hover_text[ranked.index.get_loc(i)] for i in subset.index]

    fig.add_trace(go.Scatter(
        x=front_p["strict_frac"].tolist(),
        y=front_p["mean_dG_hull"].tolist(),
        mode="markers",
        marker={"size": 12, "color": "rgba(0,0,0,0)",
                "line": {"color": "crimson", "width": 2}},
        text=_safe_hover_idx(front_p),
        hoverinfo="text",
        name="Primary front (dG+strict)",
    ))
    fig.add_trace(go.Scatter(
        x=front_a["strict_frac"].tolist(),
        y=front_a["mean_dG_hull"].tolist(),
        mode="markers",
        marker={"size": 12, "color": "rgba(0,0,0,0)",
                "line": {"color": "royalblue", "width": 2}},
        text=_safe_hover_idx(front_a),
        hoverinfo="text",
        name="Alt front (dG+ss_auc)",
        visible="legendonly",
    ))
    fig.add_trace(go.Scatter(
        x=front_3["strict_frac"].tolist(),
        y=front_3["mean_dG_hull"].tolist(),
        mode="markers",
        marker={"size": 12, "color": "rgba(0,0,0,0)",
                "line": {"color": "darkorange", "width": 2}},
        text=_safe_hover_idx(front_3),
        hoverinfo="text",
        name="3-obj front",
        visible="legendonly",
    ))

    fig.update_layout(
        title="AlB₂ diboride stability — 2D Pareto (strict_frac vs mean_dG_hull)",
        xaxis_title="Strict single-phase fraction",
        yaxis_title="Mean ΔG above hull (eV / f.u.)",
        template="plotly_white",
        legend={"itemsizing": "constant"},
        height=650,
    )

    p = plot_dir / "pareto_2D_interactive.html"
    fig.write_html(str(p), include_plotlyjs=True)
    print(f"2D interactive → {p}")


def _plot_interactive_3d(metrics: pd.DataFrame, plot_dir: Path) -> None:
    """Save a self-contained 3D interactive Pareto scatter as HTML.

    X = strict_frac, Y = mean_dG_hull, Z = ss_auc.

    Args:
        metrics: DataFrame with Pareto rank and objective columns.
        plot_dir: Output directory.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("plotly not installed; skipping 3D interactive plot.")
        return

    ranked = metrics.dropna(
        subset=["pareto_rank", "mean_dG_hull", "ss_auc", "strict_frac"]
    ).copy().sort_values("pareto_rank")
    if ranked.empty:
        return

    front = ranked[ranked["pareto_rank"] == 1]
    rest  = ranked[ranked["pareto_rank"] != 1]

    def _hover3(row: pd.Series) -> str:
        p_rank  = int(row["pareto_rank"])
        a_rank  = int(row["pareto_rank_alt"])  if pd.notna(row.get("pareto_rank_alt"))  else -1
        t3_rank = int(row["pareto_rank_3obj"]) if pd.notna(row.get("pareto_rank_3obj")) else -1
        return (
            f"<b>{row['system']}</b><br>"
            f"{row.get('comp_label', row['composition'])}<br>"
            f"mean_dG_hull: {row['mean_dG_hull']:.4f} eV<br>"
            f"strict_frac: {row['strict_frac']:.4f}<br>"
            f"ss_auc: {row['ss_auc']:.4f}<br>"
            f"unique_frac: {row.get('unique_frac', float('nan')):.4f}<br>"
            f"Primary rank: {p_rank} | Alt: {a_rank} | 3-obj: {t3_rank}"
        )

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=rest["strict_frac"].tolist(),
        y=rest["mean_dG_hull"].tolist(),
        z=rest["ss_auc"].tolist(),
        mode="markers",
        marker={"size": 3, "color": rest["strict_frac"].tolist(),
                "colorscale": "Viridis", "opacity": 0.4, "cmin": 0, "cmax": 1},
        text=[_hover3(r) for _, r in rest.iterrows()],
        hoverinfo="text",
        name="All",
    ))
    fig.add_trace(go.Scatter3d(
        x=front["strict_frac"].tolist(),
        y=front["mean_dG_hull"].tolist(),
        z=front["ss_auc"].tolist(),
        mode="markers",
        marker={"size": 7, "color": "crimson", "opacity": 0.9,
                "line": {"color": "white", "width": 1}},
        text=[_hover3(r) for _, r in front.iterrows()],
        hoverinfo="text",
        name="Primary Pareto front",
    ))
    fig.update_layout(
        title="AlB₂ diboride — 3-objective Pareto space (strict_frac, dG_hull, ss_auc)",
        scene={
            "xaxis_title": "Strict single-phase fraction",
            "yaxis_title": "Mean ΔG above hull (eV)",
            "zaxis_title": "AlB₂ retention AUC",
        },
        template="plotly_white",
        height=700,
    )

    p = plot_dir / "pareto_3D_interactive.html"
    fig.write_html(str(p), include_plotlyjs=True)
    print(f"3D interactive → {p}")


def _bool_to_float_nan(series: pd.Series) -> np.ndarray:
    """Convert boolean column to float, preserving NaN as np.nan.

    ``fillna(False).astype(float)`` silently treats unconfirmed steps as 0
    (not stable), which is misleading for interpolated or missing T steps.
    This function leaves NaN as NaN so matplotlib renders a gap instead.
    """
    return np.where(series.isna(), np.nan, series.astype(float).values)


def _pad_with_nan_gaps(grp: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Pad a per-composition T-sorted DataFrame with NaN rows at missing T steps.

    Entirely absent T steps (rows not in the raw parquet at all) would otherwise
    be silently skipped, connecting the line across a gap.  This function inserts
    a NaN row at each missing expected T step so matplotlib renders a visual gap.

    Args:
        grp: DataFrame for one composition, sorted by T.
        cols: Columns to carry forward as NaN in gap rows.

    Returns:
        Padded DataFrame with NaN rows inserted at missing T steps.
    """
    observed = set(grp["T"].tolist())
    missing  = sorted(_EXPECTED_T_SET - observed)
    if not missing:
        return grp
    gap_rows = [{
        **{c: float("nan") for c in cols},
        "T": t,
    } for t in missing]
    return (
        pd.concat([grp, pd.DataFrame(gap_rows)], ignore_index=True)
        .sort_values("T")
        .reset_index(drop=True)
    )


def _plot_temperature_stability(
    raw: pd.DataFrame,
    metrics: pd.DataFrame,
    plot_dir: Path,
    n: int = min(top_n, 8),
    interp_rows: pd.DataFrame | None = None,
) -> None:
    """Save combined temperature stability figures as PNG and PDF.

    4-panel figure (ss_fraction, dG_above_hull, is_single_phase_strict,
    is_unique_homogeneous) for up to N top compositions.  Also saves individual
    per-composition profiles to temperature_profiles/.

    Entirely absent T steps (not in raw, not NaN) are shown as explicit NaN
    gaps by padding the per-composition data with NaN rows at missing T steps.
    Interpolated continuous data (from ``interp_rows``) is overlaid with dashed
    markers to distinguish it from calculated data.
    NaN boolean flags remain as gaps rather than being treated as False.
    Interpolated compositions are labelled ``[i]`` in the legend.

    Args:
        raw: Raw stability DataFrame (calculated data only).
        metrics: Metrics DataFrame with ``pareto_rank`` column.
        plot_dir: Output directory.
        n: Number of top compositions in the combined figure.
        interp_rows: Optional DataFrame of interpolated rows from
            ``_compute_pareto_metrics``.  Plotted with dashed style.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib as _mpl
    import matplotlib.pyplot as plt
    get_cmap = (
        _mpl.colormaps.get_cmap
        if hasattr(_mpl.colormaps, "get_cmap")
        else _mpl.cm.get_cmap
    )

    top = (
        metrics.dropna(subset=["pareto_rank"])
        .sort_values(_SECONDARY_SORT, ascending=_SECONDARY_ASC)
        .head(top_n)
    )
    if top.empty:
        return

    top_keys = top[["system", "composition"]].assign(_k=True)
    raw_top  = (
        raw.merge(top_keys, on=["system", "composition"], how="inner")
        .drop(columns=["_k"])
    )

    # Build a lookup of interpolated rows per (system, composition) for overlay.
    _plot_cols = ["T", "ss_fraction", "dG_above_hull",
                  "is_single_phase_strict", "is_unique_homogeneous"]
    interp_lookup: dict[tuple[str, str], pd.DataFrame] = {}
    if interp_rows is not None and not interp_rows.empty:
        for (sys_, comp_), ig in interp_rows.groupby(["system", "composition"], sort=False):
            interp_lookup[(str(sys_), str(comp_))] = ig.sort_values("T")

    n_show = min(n, len(top))
    cmap   = get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(n_show)]

    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    ax_ss, ax_dg, ax_st, ax_un = axes

    for ci, (_, row) in enumerate(list(top.iterrows())[:n_show]):
        grp = raw_top[
            (raw_top["system"] == row["system"]) &
            (raw_top["composition"] == row["composition"])
        ].sort_values("T")
        if grp.empty:
            continue
        is_interp = bool(row.get("interpolated", False))
        label = (
            f"{row['system']} [{int(row['pareto_rank'])}]"
            + (" [i]" if is_interp else "")
        )
        col = colors[ci]

        # Insert NaN at entirely absent T steps so the line shows a visible gap.
        grp_padded = _pad_with_nan_gaps(grp, _plot_cols)
        t = grp_padded["T"].to_numpy(dtype=float)

        # Continuous series: NaN → gap (not 0).
        ax_ss.plot(t, grp_padded["ss_fraction"].to_numpy(dtype=float),
                   color=col, lw=1.3, label=label)
        ax_dg.plot(t, grp_padded["dG_above_hull"].to_numpy(dtype=float),
                   color=col, lw=1.3)
        # Boolean flags: NaN → gap (not False).
        ax_st.plot(t, _bool_to_float_nan(grp_padded["is_single_phase_strict"]),
                   color=col, lw=1.3)
        ax_un.plot(t, _bool_to_float_nan(grp_padded["is_unique_homogeneous"]),
                   color=col, lw=1.3)

        # Overlay interpolated data points with dashed markers.
        key = (str(row["system"]), str(row["composition"]))
        if key in interp_lookup:
            ig = interp_lookup[key]
            t_i = ig["T"].to_numpy(dtype=float)
            ax_ss.plot(t_i, ig["ss_fraction"].to_numpy(dtype=float),
                       color=col, lw=0, marker="x", ms=5, alpha=0.7)
            ax_dg.plot(t_i, ig["dG_above_hull"].to_numpy(dtype=float),
                       color=col, lw=0, marker="x", ms=5, alpha=0.7)

    ax_ss.axhline(0.95, color="gray", ls="--", lw=0.8, label="ss_threshold=0.95")
    ax_ss.set_ylabel("ss_fraction")
    ax_ss.set_ylim(-0.05, 1.05)
    ax_ss.legend(fontsize=7, ncol=2, loc="lower left")

    ax_dg.axhline(0.001, color="gray", ls="--", lw=0.8, label="dG_hull_tol=1 meV")
    ax_dg.set_ylabel("ΔG above hull (eV)")
    ax_dg.legend(fontsize=7, loc="upper right")

    ax_st.set_ylabel("is_single_phase_strict")
    ax_st.set_ylim(-0.05, 1.15)
    ax_st.text(0.01, 1.08, "gap = unconfirmed / interpolated T step",
               transform=ax_st.transAxes, fontsize=6.5, color="gray")

    ax_un.set_ylabel("is_unique_homogeneous")
    ax_un.set_ylim(-0.05, 1.15)
    ax_un.set_xlabel("Temperature (K)")
    ax_un.text(0.01, 1.08, "gap = unconfirmed / interpolated T step",
               transform=ax_un.transAxes, fontsize=6.5, color="gray")

    fig.suptitle(
        f"Temperature stability profiles — top AlB₂ diboride compositions  "
        f"([i] = interpolated)  [{_T_RANGE_LABEL}]",
        fontsize=12, y=0.999,
    )
    plt.tight_layout()

    for ext in ("png", "pdf"):
        p = plot_dir / f"temperature_stability.{ext}"
        fig.savefig(p, dpi=150)
        print(f"Temperature stability → {p}")
    plt.close(fig)

    # Individual per-composition profiles.
    prof_dir = plot_dir / "temperature_profiles"
    prof_dir.mkdir(parents=True, exist_ok=True)

    for _, row in top.iterrows():
        grp = raw_top[
            (raw_top["system"] == row["system"]) &
            (raw_top["composition"] == row["composition"])
        ].sort_values("T")
        if grp.empty:
            continue

        is_interp = bool(row.get("interpolated", False))
        fig2, ax2s = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

        # Insert NaN at entirely absent T steps.
        grp_padded2 = _pad_with_nan_gaps(
            grp, ["ss_fraction", "dG_above_hull",
                  "is_parent_structure_retained",
                  "is_single_phase_strict", "is_unique_homogeneous"]
        )
        t = grp_padded2["T"].to_numpy(dtype=float)

        # Continuous: NaN → gap.
        ax2s[0].plot(t, grp_padded2["ss_fraction"].to_numpy(dtype=float),
                     color="steelblue", lw=1.5, label="calculated")
        ax2s[0].axhline(0.95, color="gray", ls="--", lw=0.8, label="ss_threshold=0.95")
        # Overlay interpolated data with dashed markers.
        key2 = (str(row["system"]), str(row["composition"]))
        if key2 in interp_lookup:
            ig2 = interp_lookup[key2]
            t_i2 = ig2["T"].to_numpy(dtype=float)
            ax2s[0].plot(t_i2, ig2["ss_fraction"].to_numpy(dtype=float),
                         color="steelblue", lw=0, marker="x", ms=6,
                         alpha=0.8, label="interpolated")
        ax2s[0].set_ylabel("ss_fraction")
        ax2s[0].set_ylim(-0.05, 1.05)
        ax2s[0].legend(fontsize=8)

        ax2s[1].plot(t, grp_padded2["dG_above_hull"].to_numpy(dtype=float),
                     color="darkorange", lw=1.5, label="calculated")
        if key2 in interp_lookup:
            ig2 = interp_lookup[key2]
            ax2s[1].plot(ig2["T"].to_numpy(dtype=float),
                         ig2["dG_above_hull"].to_numpy(dtype=float),
                         color="darkorange", lw=0, marker="x", ms=6,
                         alpha=0.8, label="interpolated")
        ax2s[1].axhline(0.001, color="gray", ls="--", lw=0.8, label="dG_hull_tol=1 meV")
        ax2s[1].set_ylabel("ΔG above hull (eV)")
        ax2s[1].legend(fontsize=8)

        # Boolean flags: NaN → gap (unconfirmed steps are not shown as 0).
        retained = _bool_to_float_nan(grp_padded2["is_parent_structure_retained"])
        strict   = _bool_to_float_nan(grp_padded2["is_single_phase_strict"])
        unique   = _bool_to_float_nan(grp_padded2["is_unique_homogeneous"])
        ax2s[2].fill_between(t, np.where(np.isnan(retained), 0, retained),
                             alpha=0.25, color="steelblue",  label="retained")
        ax2s[2].fill_between(t, np.where(np.isnan(strict),   0, strict),
                             alpha=0.40, color="darkorange", label="strict")
        ax2s[2].fill_between(t, np.where(np.isnan(unique),   0, unique),
                             alpha=0.55, color="crimson",    label="unique")
        # Overlay NaN regions as hatched grey to mark unconfirmed steps.
        nan_mask_strict = np.isnan(strict)
        if nan_mask_strict.any():
            ax2s[2].fill_between(t, np.where(nan_mask_strict, 1, 0),
                                 alpha=0.15, color="gray", hatch="//",
                                 label="unconfirmed")
        ax2s[2].set_ylabel("Stability flag")
        ax2s[2].set_xlabel("Temperature (K)")
        ax2s[2].set_ylim(-0.05, 1.15)
        ax2s[2].legend(fontsize=8)

        interp_note = " [interpolated — provisional stability flags]" if is_interp else ""
        fig2.suptitle(
            f"{row['system']} — {row.get('comp_label', row['composition'])}"
            f"  [rank {int(row['pareto_rank'])}]{interp_note}",
            fontsize=10, y=0.998,
        )
        plt.tight_layout()
        safe = f"{row['system']}_{row['composition'].replace(',', '_')}.png"
        fig2.savefig(prof_dir / safe, dpi=120)
        plt.close(fig2)

    print(f"Individual profiles → {prof_dir}/  ({len(top)} files)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    """Run the full Pareto analysis pipeline."""
    print(
        "Pareto frontiers to compute:\n"
        "  Primary   : minimize mean_dG_hull | maximize strict_frac\n"
        "  Alternative: minimize mean_dG_hull | maximize ss_auc\n"
        "  3-objective: minimize mean_dG_hull | maximize strict_frac | maximize ss_auc"
    )

    out_dir.mkdir(parents=True, exist_ok=True)

    raw, summary = _load_data()

    # Restrict analysis to the configured temperature range.
    if analysis_T_min != _DB_T_MIN or analysis_T_max != _DB_T_MAX:
        n_before = len(raw)
        raw = raw[
            raw["T"].between(analysis_T_min, analysis_T_max, inclusive="both")
        ].copy()
        print(
            f"Temperature filter ({_T_RANGE_LABEL}): "
            f"{n_before:,} → {len(raw):,} rows retained."
        )

    metrics, incompl_df, interp_df = _compute_pareto_metrics_fast(raw)

    if metrics.empty:
        raise RuntimeError("No compositions with sufficient T coverage.")

    issues = _validate_metrics(raw, metrics)
    _finish(
        metrics, incompl_df, interp_df, issues, raw, summary,
        supplementary=lambda m: _compute_supplementary_metrics(raw, m),
    )  # fmt: skip


def _finish(
    metrics: pd.DataFrame,
    incompl_df: pd.DataFrame,
    interp_df: pd.DataFrame,
    issues: list[str],
    raw: pd.DataFrame | None,
    summary: pd.DataFrame,
    supplementary,
    raw_loader=None,
) -> None:
    """Rank, add supplementary metrics, and write every output (shared by both drivers).

    Args:
        metrics: Per-composition metrics.
        incompl_df: Compositions excluded for incomplete T coverage.
        interp_df: Interpolated raw rows (for plots).
        issues: Validation warnings.
        raw: Raw rows needed for outputs (top compositions at least).
        summary: Summary table (may be empty).
        supplementary: Callable ``metrics -> metrics`` adding curvature / FSM columns.
        raw_loader: When ``raw`` is None, callable ``metrics -> raw`` returning the raw
            rows of the compositions shown in the workbook and plots (after ranking).
    """
    if issues:
        for iss in issues:
            warnings.warn(iss)
    else:
        print("Validation: all checks passed.")

    # Save metrics parquet before ranking so it survives if ranking is interrupted.
    mp = out_dir / "pareto_metrics.parquet"
    metrics.to_parquet(mp, index=False)
    print(f"Metrics parquet (pre-rank) → {mp}")

    metrics, comparison = _add_pareto_columns(metrics)

    # Supplementary stability metrics (curvature + finite-separation margin).
    metrics = supplementary(metrics)

    # Within-system candidate ranking (s_config + candidate_rank).
    # Must run after supplementary metrics so hessian/FSM columns exist.
    metrics = _assign_candidate_ranks(metrics)
    n_eligible = int(metrics["candidate_rank"].notna().sum())
    print(f"  Candidate ranking: {n_eligible} eligible rank-1 compositions ranked within their systems.")

    # Overwrite with final ranked metrics + supplementary + candidate_rank columns.
    metrics.to_parquet(mp, index=False)
    print(f"Metrics parquet (ranked) → {mp}")

    if raw is None:
        raw = raw_loader(metrics) if raw_loader is not None else pd.DataFrame()

    _print_top_n(metrics, n=top_n)

    sensitivity = _sensitivity_comparison(metrics)
    if sensitivity.get("n_interpolated_total", 0) > 0:
        print(f"\nInterpolation sensitivity: {sensitivity['description']}")

    print("\nEpsilon sensitivity (rank-1 count across pareto_eps values):")
    eps_sens = _rank1_eps_sensitivity(metrics)
    print(f"  {eps_sens['description']}")

    _write_summary_txt(
        metrics, comparison, issues, out_dir / "analysis_summary.txt",
        sensitivity=sensitivity,
        eps_sensitivity=eps_sens,
        n_incomplete=len(incompl_df),
    )

    _build_excel(
        metrics, raw, summary, comparison, incompl_df,
        out_dir / "Diboride_Synthesizability_Ranking.xlsx",
    )

    _plot_2d_static(metrics, out_dir)
    _plot_interactive_2d(metrics, out_dir)
    _plot_interactive_3d(metrics, out_dir)
    _plot_temperature_stability(raw, metrics, out_dir, interp_rows=interp_df)

    front_n = int((metrics["pareto_rank"] == 1).sum())
    ranked  = metrics.dropna(subset=["pareto_rank"])
    hc_n = int(
        (
            (ranked["pareto_rank"] == 1) &
            (ranked["strict_frac"] == 1.0) &
            (ranked["unique_frac"] == 1.0) &
            (ranked["mean_dG_hull"].abs() < 1e-6) &
            (~ranked["interpolated"].fillna(False))
        ).sum()
    )
    print(
        f"\nDone.  Primary Pareto front: {front_n} compositions  |  "
        f"High-confidence homogeneous rank-1 (observed-only): {hc_n}  |  "
        f"Output → {out_dir.resolve()}/"
    )


def _configure_temperature_grid(t_min: float, t_max: float, t_step: float, a_min: float, a_max: float) -> None:
    """Set the database grid and analysis range, and recompute derived module state.

    Raises:
        ValueError: If the analysis range lies outside the grid or is not aligned to it.
    """
    global _DB_T_MIN, _DB_T_MAX, _DB_T_STEP, analysis_T_min, analysis_T_max, t_norm_min, t_norm_max
    global _EXPECTED_T_GRID, _N_EXPECTED_T, _EXPECTED_T_SET, _T_RANGE_LABEL
    if t_step <= 0 or t_max <= t_min:
        raise ValueError(f"Invalid temperature grid: t_min={t_min}, t_max={t_max}, t_step={t_step}")
    if not (t_min <= a_min < a_max <= t_max):
        raise ValueError(f"Analysis range {a_min}–{a_max} K must lie within the stability grid {t_min}–{t_max} K.")
    for name, val in (("analysis_t_min", a_min), ("analysis_t_max", a_max)):
        if abs(((val - t_min) / t_step) - round((val - t_min) / t_step)) > 1e-9:
            raise ValueError(f"{name}={val} is not on the {t_step:g} K grid starting at {t_min:g} K.")
    _DB_T_MIN, _DB_T_MAX, _DB_T_STEP = float(t_min), float(t_max), float(t_step)
    analysis_T_min, analysis_T_max = float(a_min), float(a_max)
    t_norm_min, t_norm_max = analysis_T_min, analysis_T_max
    _EXPECTED_T_GRID = [round(float(t), 6) for t in np.arange(a_min, a_max + t_step / 2, t_step)]
    _N_EXPECTED_T = len(_EXPECTED_T_GRID)
    _EXPECTED_T_SET = set(_EXPECTED_T_GRID)
    _T_RANGE_LABEL = f"{int(analysis_T_min)}–{int(analysis_T_max)} K"


_OPTION_NAMES = {
    "min_valid_frac", "interpolate_incomplete", "max_interpolation_gap_K", "discontinuity_threshold",
    "pareto_eps", "top_n", "top_per_system", "rank_filter_max_dG_hull", "rank_filter_min_strict_frac",
    "use_fast_pareto", "max_pareto_rank", "excel_max_rows", "verbose_ranking", "shortlist_per_system", "candidate_hessian_tol",
    "candidate_fsm_tol", "fsm_min_sep", "supp_metrics_rank1_only",
}  # fmt: skip


def run_pareto_analysis(
    db_dir: str | Path,
    out_dir: str | Path,
    t_min: float,
    t_max: float,
    t_step: float,
    analysis_t_min: float | None = None,
    analysis_t_max: float | None = None,
    **options,
) -> Path:
    """Run the full Pareto analysis on a stability database.

    Args:
        db_dir: Directory holding ``stability_raw.parquet`` and ``stability_summary.parquet``.
        out_dir: Output directory (created).
        t_min: Lowest temperature of the stability grid (K).
        t_max: Highest temperature of the stability grid (K).
        t_step: Temperature step of the stability grid (K).
        analysis_t_min: Lower bound of the analysed range (default ``t_min``).
        analysis_t_max: Upper bound of the analysed range (default ``t_max``).
        **options: Overrides for module settings, e.g. ``top_n``, ``pareto_eps``,
            ``rank_filter_min_strict_frac``, ``fsm_min_sep`` (see the CONFIG block).

    Returns:
        The output directory.

    Raises:
        ValueError: On unknown options or an invalid temperature range.
    """
    unknown = sorted(set(options) - _OPTION_NAMES)
    if unknown:
        raise ValueError(f"Unknown Pareto option(s) {unknown}; valid: {sorted(_OPTION_NAMES)}")
    _configure_temperature_grid(
        t_min, t_max, t_step,
        t_min if analysis_t_min is None else analysis_t_min,
        t_max if analysis_t_max is None else analysis_t_max,
    )  # fmt: skip
    module = globals()
    module.update(options)
    module["db_dir"] = Path(db_dir)
    module["out_dir"] = Path(out_dir)
    main()
    return Path(out_dir)


def _read_stability_file(path: Path) -> pd.DataFrame:
    """Read one per-system ``<system>_stability.parquet`` in the raw-database schema."""
    raw = pd.read_parquet(path).rename(columns={"dG_hull": "dG_above_hull", "G_lp": "G_total"})
    if analysis_T_min != _DB_T_MIN or analysis_T_max != _DB_T_MAX:
        raw = raw[raw["T"].between(analysis_T_min, analysis_T_max, inclusive="both")]
    return raw


def run_pareto_analysis_multi(
    stability_dirs: list[str | Path],
    out_dir: str | Path,
    t_min: float,
    t_max: float,
    t_step: float,
    analysis_t_min: float | None = None,
    analysis_t_max: float | None = None,
    **options,
) -> Path:
    """Pareto analysis over the per-system results of one or more stability runs.

    Scales to millions of compositions (e.g. 2-5 metal systems together): each
    ``<system>_stability.parquet`` is read and reduced to per-composition metrics on
    its own, so the raw rows are never all in memory. Ranking, supplementary metrics
    (computed system by system), the Excel workbook and plots then run on the combined
    metrics exactly as in :func:`run_pareto_analysis`. Large Excel sheets are capped at
    ``excel_max_rows``; ``max_pareto_rank`` limits how many fronts are peeled.

    Args:
        stability_dirs: Output folders of :class:`BladeStability` runs (all on the same
            temperature grid and without oxides, i.e. one row per composition and T).
        out_dir: Output directory (created).
        t_min: Lowest temperature of the stability grid (K).
        t_max: Highest temperature of the stability grid (K).
        t_step: Temperature step of the stability grid (K).
        analysis_t_min: Lower bound of the analysed range (default ``t_min``).
        analysis_t_max: Upper bound of the analysed range (default ``t_max``).
        **options: Same overrides as :func:`run_pareto_analysis`.

    Returns:
        The output directory.
    """
    import contextlib
    import io

    unknown = sorted(set(options) - _OPTION_NAMES)
    if unknown:
        raise ValueError(f"Unknown Pareto option(s) {unknown}; valid: {sorted(_OPTION_NAMES)}")
    _configure_temperature_grid(
        t_min, t_max, t_step,
        t_min if analysis_t_min is None else analysis_t_min,
        t_max if analysis_t_max is None else analysis_t_max,
    )  # fmt: skip
    module = globals()
    module.update(options)
    module["db_dir"] = Path(stability_dirs[0])
    module["out_dir"] = Path(out_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files: dict[str, Path] = {}
    for d in stability_dirs:
        for f in sorted(Path(d).glob("*_stability.parquet")):
            system = f.name[: -len("_stability.parquet")]
            if system in files:
                raise ValueError(f"System {system} appears in more than one stability folder")
            files[system] = f
    print(f"Pareto (multi): {len(files):,} systems from {len(stability_dirs)} folder(s)")

    parts, incompl, interp, issues = [], [], [], []
    for k, (system, f) in enumerate(files.items(), 1):
        raw = _read_stability_file(f)
        if "mu_O" in raw.columns and raw["mu_O"].notna().any():
            raise ValueError(f"{f} includes oxides (mu_O); run the Pareto analysis per mu_O instead")
        _validate_schema(raw)
        with contextlib.redirect_stdout(io.StringIO()):
            m, inc, itp = _compute_pareto_metrics_fast(raw)
            issues += [f"{system}: {msg}" for msg in _validate_metrics(raw, m)]
        parts.append(m)
        incompl.append(inc)
        interp.append(itp)
        if k % 200 == 0 or k == len(files):
            print(f"  metrics: {k:,}/{len(files):,} systems", flush=True)
    metrics = pd.concat(parts, ignore_index=True)
    incompl_df = pd.concat([x for x in incompl if len(x)], ignore_index=True) if any(len(x) for x in incompl) else pd.DataFrame()
    interp_df = pd.concat([x for x in interp if len(x)], ignore_index=True) if any(len(x) for x in interp) else pd.DataFrame()
    if metrics.empty:
        raise RuntimeError("No compositions with sufficient T coverage.")
    print(f"  → {len(metrics):,} compositions with metrics  |  {len(incompl_df):,} excluded")

    def supplementary(m_all: pd.DataFrame) -> pd.DataFrame:
        """Curvature / FSM per system, loading only systems that have rank-1 compositions."""
        target = m_all["pareto_rank"] == 1 if supp_metrics_rank1_only else m_all["pareto_rank"].notna()
        systems = m_all.loc[target, "system"].unique()
        print(f"Computing supplementary metrics for {int(target.sum()):,} compositions in {len(systems):,} systems …")
        groups = m_all.groupby("system", sort=False).indices
        # Defaults for every row (NaN / 0 / False), exactly as the single-file path sets
        # them for compositions it does not compute: raw without G_ss returns defaults only.
        with contextlib.redirect_stdout(io.StringIO()):
            out = _compute_supplementary_metrics(pd.DataFrame(columns=["system", "T"]), m_all)
        new_cols = [c for c in out.columns if c not in m_all.columns]
        for k, system in enumerate(systems, 1):
            idx = groups[system]
            with contextlib.redirect_stdout(io.StringIO()):
                res = _compute_supplementary_metrics(_read_stability_file(files[system]), m_all.iloc[idx])
            out.loc[res.index, new_cols] = res[new_cols].to_numpy()
            if k % 200 == 0 or k == len(systems):
                print(f"  supplementary: {k:,}/{len(systems):,} systems", flush=True)
        return out

    def raw_for_outputs(m_all: pd.DataFrame) -> pd.DataFrame:
        """Raw rows of the compositions shown in the workbook and plots."""
        keep = m_all["pareto_rank"].le(2).fillna(False)
        if "candidate_rank" in m_all.columns:
            keep |= m_all["candidate_rank"].notna()
        keys = m_all.loc[keep, ["system", "composition"]]
        frames = []
        for system, comps in keys.groupby("system")["composition"]:
            r = _read_stability_file(files[system])
            frames.append(r[r["composition"].isin(set(comps))])
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    _finish(
        metrics, incompl_df, interp_df, issues, raw=None, summary=pd.DataFrame(),
        supplementary=supplementary, raw_loader=raw_for_outputs,
    )  # fmt: skip
    return out_dir


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Pareto analysis of BLADE stability results")
    parser.add_argument("dirs", nargs="+", help="BladeStability output folders (per-system parquet files)")
    parser.add_argument("--out", required=True, help="output folder")
    parser.add_argument("--t-range", nargs=3, type=float, default=[200.0, 2500.0, 100.0], metavar=("MIN", "MAX", "STEP"))
    parser.add_argument("--max-rank", type=int, default=None, help="peel at most this many fronts")
    a = parser.parse_args()
    run_pareto_analysis_multi(a.dirs, a.out, *a.t_range, max_pareto_rank=a.max_rank)
