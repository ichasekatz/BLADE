"""Pure solid-state thermodynamic stability analysis for BLADE diboride systems.

Computes equilibrium phase fractions and ΔG from the convex hull over the
Redlich-Kister SS phase + fixed SQS-structure competitors, with no oxygen
chemical potential. Outputs per-system parquet files and a combined summary.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from blade.oxidation.framework._lp import solve_grand_lp

from ._stability_tables import MissingEndpointError

gas_constant = 8.314472  # J / (mol·K)
kb_ev_per_k = 8.617333262145e-5  # eV / K


class StabilityConfigError(ValueError):
    """Raised when a required configuration value is missing or invalid."""


@dataclass(frozen=True)
class StabilityThresholds:
    """Criteria for the per-(composition, T) single-phase flags.

    Attributes:
        parent_threshold: Minimum ``ss_fraction`` (metal fraction in AlB₂-type
            phases) for ``is_parent_structure_retained``.
        strict_dg_tol: Maximum ``dG_hull`` (eV/f.u.) for ``is_single_phase_strict``.
        min_metal_frac: Every metal must reach this site fraction for the
            composition to count as a mixed alloy (excludes trivial endpoints).
        lp_tol: Phase amounts at or below this are treated as inactive.
    """

    parent_threshold: float = 0.95
    strict_dg_tol: float = 1e-3
    min_metal_frac: float = 0.05
    lp_tol: float = 1e-6


class BladeStability:
    """Pure solid-state phase stability scanner for BLADE system_structures outputs.

    Reads relaxed SQS structures from each system's ``blade/`` subdirectory,
    fits the Redlich-Kister/Muggianu interaction model, and runs the grand-potential
    LP (without oxygen) over a composition simplex grid × temperature range.

    Additional competing phases from xlsx files in a configurable subdir within
    each system directory are included as fixed-phase LP competitors when present,
    expanding the convex hull beyond the AlB₂-type SS phase alone.

    For each (composition, T) point the LP finds the minimum free energy over all
    competing phases (fixed SQS structures + continuous SS grid). The SS-phase
    stability is reported as:

        dG_hull = G_ss(y, T) − G_lp(y, T)   [eV per formula unit]

    where G_lp ≤ G_ss always, so dG_hull ≥ 0; dG_hull ≈ 0 (≤ 1e-6 eV) means the SS
    phase is on the hull (stable) at that composition and temperature.

    Oxides: with ``exclude_oxides=True`` (default) every O-bearing competitor is
    dropped and the scan is purely solid-state. With ``exclude_oxides=False`` oxide
    competitors are loaded and the system is treated as open to oxygen at each
    chemical potential in ``mu_o_values`` (absolute eV per O atom, same convention as
    the oxidation stage): fixed-phase energies become ``G_j − n_O,j · μ_O`` and the
    output gains one block of rows per μ_O. A closed, O-free system cannot form
    oxides at all, so including oxides without a μ_O is rejected.

    Args:
        system_structures_root: Path to the ``system_structures/`` directory produced
            by BLADE (each subdirectory is one composition system).
        output_dir: Directory for stability results (per-system parquet + combined).
        phase_element: Non-metal sublattice element (e.g. ``"B"``).
        phase_stoich: Stoichiometry per formula unit (e.g. ``2.0`` for AlB₂-type).
        competitors_subdir: Subpath (or list of subpaths) within each system dir
            holding competing phases — structure subdirs with ``CONTCAR`` + ``energy``
            (e.g. ``"ORB"``) and/or xlsx tables (e.g. ``"parent_excel/MP_ORB"``).
            Merged with the lowest energy per formula kept. Missing dirs are skipped;
            None loads no competitors.
        exclude_oxides: Drop O-bearing competitors (default True).
        mu_o_values: Oxygen chemical potentials (eV per O atom) — required when
            ``exclude_oxides`` is False, ignored otherwise.
        allowed_elements: Only process systems whose elements are all in this set.
        t_min: Minimum temperature in K (default 300).
        t_max: Maximum temperature in K (default 2000).
        t_step: Temperature step in K (default 100).
        composition_step: Simplex grid spacing (default 0.02).
        skip_existing: Skip systems whose parquet already exists (default True).
        rk_order: Redlich-Kister polynomial order (default 3).
        n_workers: Number of parallel worker processes (default: CPU count).
        thresholds: Criteria for the single-phase flags (see :class:`StabilityThresholds`).
        run_pareto: After the scan, write ``stability_raw.parquet`` /
            ``stability_summary.parquet`` and run the Pareto analysis
            (:mod:`blade.analysis.stability_pareto`). With oxides included this runs
            once per μ_O value, in ``mu_O=<value>/`` subfolders.
        pareto_options: Overrides for the Pareto analysis, e.g. ``output_folder``,
            ``analysis_t_min``, ``analysis_t_max``, ``top_n``,
            ``rank_filter_min_strict_frac`` (see ``run_pareto_analysis``).
    """

    def __init__(
        self,
        system_structures_root: str | Path,
        output_dir: str | Path,
        phase_element: str = "B",
        phase_stoich: float = 2.0,
        competitors_subdir: str | list[str] | None = None,
        exclude_oxides: bool = True,
        mu_o_values: list[float] | None = None,
        allowed_elements: set[str] | None = None,
        t_min: float = 300.0,
        t_max: float = 2000.0,
        t_step: float = 100.0,
        composition_step: float = 0.02,
        skip_existing: bool = True,
        rk_order: int = 3,
        n_workers: int | None = None,
        thresholds: StabilityThresholds | None = None,
        run_pareto: bool = False,
        pareto_options: dict | None = None,
    ) -> None:
        """Initialize with file paths and stability scan parameters."""
        self._system_structures_root = Path(system_structures_root)
        self._output_dir = Path(output_dir)
        self.phase_element = phase_element
        self.phase_stoich = phase_stoich
        if isinstance(competitors_subdir, str):
            competitors_subdir = [competitors_subdir] if competitors_subdir else None
        self._competitors_subdir = tuple(competitors_subdir) if competitors_subdir else None
        self._exclude_oxides = bool(exclude_oxides)
        if self._exclude_oxides:
            self._mu_o_values: tuple[float, ...] = ()
        else:
            if not mu_o_values:
                raise StabilityConfigError(
                    "exclude_oxides = false requires mu_o_values (oxygen chemical potential, eV per O atom): "
                    "without an oxygen reservoir no oxide can form, so including oxides would change nothing."
                )
            self._mu_o_values = tuple(float(m) for m in mu_o_values)
        self._allowed_elements = allowed_elements
        self.t_min = t_min
        self.t_max = t_max
        self.t_step = t_step
        self.composition_step = composition_step
        self.skip_existing = skip_existing
        self.rk_order = rk_order
        self.n_workers = n_workers or mp.cpu_count()
        self.thresholds = thresholds or StabilityThresholds()
        self.run_pareto = bool(run_pareto)
        self.pareto_options = dict(pareto_options or {})

    @cached_property
    def output_dir(self) -> Path:
        """Absolute path to the Stability output directory (created on access)."""
        self._output_dir.mkdir(parents=True, exist_ok=True)
        return self._output_dir

    @cached_property
    def t_values(self) -> np.ndarray:
        """Temperature grid in K from t_min to t_max (inclusive) at t_step intervals."""
        return np.arange(self.t_min, self.t_max + 1e-9, self.t_step)

    def run(self) -> Path:
        """Run stability analysis for all systems found in system_structures_root.

        Returns:
            Path to the combined ``stability_results.parquet`` file.
        """
        system_dirs = self._find_system_dirs()
        print(f"BladeStability: found {len(system_dirs)} systems in {self._system_structures_root}")

        args = [
            (
                sdir,
                bdir,
                self._competitors_subdir,
                self._exclude_oxides,
                self._mu_o_values,
                self.phase_element,
                self.phase_stoich,
                self.t_values,
                self.composition_step,
                self.rk_order,
                self.output_dir,
                self.skip_existing,
                self.thresholds,
            )
            for sdir, bdir in system_dirs
        ]

        oxide_note = f"oxides INCLUDED at mu_O = {list(self._mu_o_values)} eV/O" if self._mu_o_values else "oxides EXCLUDED"
        print(f"BladeStability: competitors_subdir={self._competitors_subdir!r}, {oxide_note}")

        all_frames: list[pd.DataFrame] = []
        # spawn (not fork) so workers import numpy fresh with single-threaded BLAS;
        # forked workers inherit a multi-threaded BLAS pool and oversubscribe the CPU.
        # The env vars must be in place before the children start (they import numpy
        # while unpickling the worker function, before any initializer runs).
        ctx = mp.get_context("spawn")
        thread_vars = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        saved_env = {v: os.environ.get(v) for v in thread_vars}
        os.environ.update(dict.fromkeys(thread_vars, "1"))
        try:
            with ctx.Pool(processes=self.n_workers, initializer=_worker_init) as pool:
                for result, err in pool.imap_unordered(_stability_worker, args):
                    if err:
                        print(f"  SKIPPED: {err}")
                    elif result is not None:
                        all_frames.append(result)
        finally:
            for var, val in saved_env.items():
                if val is None:
                    os.environ.pop(var, None)
                else:
                    os.environ[var] = val

        if not all_frames:
            print("No stability results produced.")
            combined = pd.DataFrame()
        else:
            combined = pd.concat(all_frames, ignore_index=True)
            combined = combined.sort_values(["system", "composition", "T"]).reset_index(drop=True)

        out = self.output_dir / "stability_results.parquet"
        combined.to_parquet(out, index=False)
        print(f"Combined results → {out}  ({len(combined)} rows, {combined['system'].nunique() if len(combined) else 0} systems)")
        if self.run_pareto and len(combined):
            self._run_pareto(combined)
        return out

    def _run_pareto(self, combined: pd.DataFrame) -> list[Path]:
        """Write the Pareto input database(s) and run the Pareto analysis.

        Args:
            combined: Combined per-(system, composition, T[, μ_O]) results.

        Returns:
            Output directories of the Pareto analyses (one per μ_O when oxides are
            included, otherwise one).
        """
        from blade.analysis.stability_pareto import run_pareto_analysis

        opts = dict(self.pareto_options)
        folder = opts.pop("output_folder", "pareto_analysis")
        if self._mu_o_values:
            blocks = [(self.output_dir / f"mu_O={mu:g}", combined[np.isclose(combined["mu_O"], mu)]) for mu in self._mu_o_values]
        else:
            blocks = [(self.output_dir, combined)]

        out_dirs: list[Path] = []
        for db_dir, frame in blocks:
            db_dir.mkdir(parents=True, exist_ok=True)
            write_pareto_database(frame, db_dir)
            out_dir = db_dir / folder
            run_pareto_analysis(
                db_dir=db_dir, out_dir=out_dir, t_min=self.t_min, t_max=self.t_max, t_step=self.t_step, **opts
            )
            out_dirs.append(out_dir)
        return out_dirs

    def _find_system_dirs(self) -> list[tuple[Path, Path]]:
        """Return (system_dir, blade_dir) pairs for all valid system_structures systems.

        A valid system has a ``blade/`` subdirectory containing ≥ 3 structure
        subdirs with both ``energy`` and ``CONTCAR`` present, and all its
        metal elements must be in ``allowed_elements`` when that set is provided.

        Returns:
            Sorted list of (system_dir, blade_dir) tuples.
        """
        result: list[tuple[Path, Path]] = []
        for sdir in sorted(self._system_structures_root.iterdir()):
            if not sdir.is_dir():
                continue
            if self._allowed_elements is not None and not _system_in_pool(sdir.name, self._allowed_elements):
                continue
            blade_dir = sdir / "blade"
            if blade_dir.is_dir() and _count_valid_sqs(blade_dir) >= 3:
                result.append((sdir, blade_dir))
        return result


# ---------------------------------------------------------------------------
# LP grid runner (pure stability, no oxygen)
# ---------------------------------------------------------------------------


def _run_stability_grid(
    phase_data: dict,
    t_values: np.ndarray,
    system_tag: str,
    mu_o_values: tuple[float, ...] = (),
    thresholds: StabilityThresholds | None = None,
) -> pd.DataFrame:
    """Run the stability LP over a composition simplex grid × temperature range.

    For each composition y_i and temperature T, solves the LP:

        min  Σ_j G_j(T) · n_j
        s.t. A_eq @ n = b_eq(y_i),  n ≥ 0

    where G_j(T) = energy_per_formula − n_O,j · μ_O for fixed phases and
    G_j(T) = H_RK(y_j) + k_B T · Σ y_k ln(y_k) for the SS phase grid. Oxygen is not
    mass-balanced (open system); with no μ_O values every fixed phase must be O-free.

    Args:
        phase_data: Output of :func:`_stability_tables.build_phase_data`.
        t_values: Temperature array in K.
        system_tag: System identifier string for the ``system`` column.
        mu_o_values: Oxygen chemical potentials (eV per O atom); empty for the
            purely solid-state (oxides excluded) scan.
        thresholds: Criteria for the parent-retention / strict / unique flags.

    Returns:
        DataFrame with columns: system, metals, composition, comp_label, x_<Metal>,
        T, mu_O, oxides_included, dG_hull (= dG_above_hull), G_ss, G_lp, stable,
        ss_fraction, n_active_ss, n_active_fixed, is_parent_structure_retained,
        is_single_phase_strict, is_unique_homogeneous, G_ss_formation, G_formation,
        s_config, s_config_norm, dominant_phase.

        ``ss_fraction`` is the fraction of metal atoms in AlB₂-type phases (the SS
        grid plus the SQS structures). Formation energies are relative to the
        pure-metal parent endpoints, ``G_ref = Σ y_k · H_endpoint_k``.
    """
    n_fixed = phase_data["n_fixed"]
    y_grid = phase_data["y_grid"]
    metals = phase_data["metals"]
    phase_stoich = phase_data["phase_element_stoichiometry"]
    phase_h0 = phase_data["phase_H0"]
    phase_mix = phase_data["phase_mix_shape"]
    fef = phase_data["fixed_energy_formula"]
    a_eq = phase_data["A_eq"]
    phase_ids = phase_data["phase_ids"]
    fixed_o = phase_data["fixed_O"]
    n = len(metals)
    thr = thresholds or StabilityThresholds()
    metals_str = ",".join(metals)
    h_ep = np.asarray(phase_data["h_endpoints"], dtype=float)
    # Per-phase metal content and AlB₂-parent mask (fixed SQS + SS grid).
    metal_content = a_eq[:n, :].sum(axis=0)
    is_parent = np.concatenate([np.asarray(phase_data["fixed_is_parent"], dtype=bool), np.ones(len(y_grid), dtype=bool)])
    parent_y = np.vstack([np.asarray(phase_data["fixed_y"], dtype=float), y_grid])

    oxides_included = len(mu_o_values) > 0
    if not oxides_included and np.any(fixed_o > 0):
        raise ValueError("O-bearing fixed phases present but no mu_O given (oxides must be excluded)")
    mu_list = list(mu_o_values) if oxides_included else [np.nan]

    rows: list[dict] = []
    for t, mu_o in ((t, mu) for t in t_values for mu in mu_list):
        g_ss_grid = phase_h0 + kb_ev_per_k * float(t) * phase_mix
        fixed_grand = fef - fixed_o * mu_o if oxides_included else fef
        grand = np.concatenate([fixed_grand, g_ss_grid])

        for i, y_i in enumerate(y_grid):
            b_eq = np.concatenate([y_i, [phase_stoich]])
            amounts, g_lp, ok = solve_grand_lp(a_eq, b_eq, grand)
            if not ok:
                # The revised-simplex solver occasionally fails numerically; the LP is
                # always feasible (the SS phase at y_i itself), so retry with HiGHS
                # rather than leave a gap in the (composition, T) grid.
                amounts, g_lp, ok = _solve_lp_highs(a_eq, b_eq, grand)
            if not ok:
                continue

            g_ss_i = grand[n_fixed + i]
            dg = g_ss_i - g_lp

            fracs = y_i[y_i > 1e-9]
            s_cfg = -gas_constant * float(np.sum(fracs * np.log(fracs)))
            s_norm = s_cfg / (-gas_constant * np.log(1.0 / n)) if n > 1 else 0.0

            dom_idx = int(np.nanargmax(amounts)) if amounts is not None else -1
            dom_phase = str(phase_ids[dom_idx]) if dom_idx >= 0 else ""

            metal_amt = amounts * metal_content
            total_metal = float(metal_amt.sum())
            ss_frac = float(metal_amt[is_parent].sum() / total_metal) if total_metal > 0 else float("nan")
            active = amounts > thr.lp_tol
            n_active_ss = int((active & is_parent).sum())
            n_active_fixed = int((active & ~is_parent).sum())
            unique = False
            if n_active_fixed == 0 and n_active_ss == 1:
                (k_act,) = np.flatnonzero(active)
                unique = bool(np.allclose(parent_y[k_act], y_i, atol=1e-6))
            all_present = bool(np.all(y_i >= thr.min_metal_frac - 1e-12))
            retained = all_present and ss_frac >= thr.parent_threshold
            g_ref = float(y_i @ h_ep)

            comp_str = ",".join(f"{y_i[k]:.4f}" for k in range(n))
            row: dict = {
                "system": system_tag,
                "metals": metals_str,
                "composition": comp_str,
                "comp_label": " ".join(f"x_{m}={y_i[k]:.3f}" for k, m in enumerate(metals)),
                "T": float(t),
                "mu_O": float(mu_o),
                "oxides_included": oxides_included,
                "dG_hull": float(dg),
                "G_ss": float(g_ss_i),
                "G_lp": float(g_lp),
                "stable": bool(dg <= 1e-6),
                "ss_fraction": ss_frac,
                "n_active_ss": n_active_ss,
                "n_active_fixed": n_active_fixed,
                "is_parent_structure_retained": retained,
                "is_single_phase_strict": bool(retained and dg <= thr.strict_dg_tol),
                "is_unique_homogeneous": bool(all_present and unique),
                "G_ss_formation": float(g_ss_i - g_ref),
                "G_formation": float(g_lp - g_ref),
                "s_config": float(s_cfg),
                "s_config_norm": float(s_norm),
                "dominant_phase": dom_phase,
            }
            for k, m in enumerate(metals):
                row[f"x_{m}"] = float(y_i[k])
            rows.append(row)

    return pd.DataFrame(rows)


def _solve_lp_highs(a_eq: np.ndarray, b_eq: np.ndarray, grand: np.ndarray) -> tuple[np.ndarray, float, bool]:
    """Solve the grand-potential LP with HiGHS (fallback for revised-simplex failures)."""
    from scipy.optimize import linprog

    res = linprog(grand, A_eq=a_eq, b_eq=b_eq, bounds=[(0.0, None)] * len(grand), method="highs")
    if res.status != 0:
        return np.full(len(grand), np.nan), float("nan"), False
    x = np.maximum(res.x, 0.0)
    x[x < 1e-10] = 0.0
    return x, float(res.fun), True


# ---------------------------------------------------------------------------
# Multiprocessing helpers
# ---------------------------------------------------------------------------


def _worker_init() -> None:
    """Ensure blade src tree is on sys.path in each spawned worker process.

    Lazy path injection so worker processes (which don't inherit the uv venv
    .pth file on macOS spawn) can still import blade.
    """
    import sys as _sys
    from pathlib import Path as _Path

    _src = _Path(__file__).parent.parent.parent
    if _src.is_dir() and str(_src) not in _sys.path:
        _sys.path.insert(0, str(_src))


def _stability_worker(
    args: tuple,
) -> tuple[pd.DataFrame | None, str | None]:
    """Process one system in a worker process.

    Intended as a ``multiprocessing.Pool.imap_unordered`` target. Imports
    lazily so worker processes don't pay import cost they don't need.

    Args:
        args: Tuple of (system_dir, blade_dir, competitors_subdir, exclude_oxides,
            mu_o_values, phase_element, phase_stoich, t_values, composition_step,
            rk_order, output_dir, skip_existing, thresholds).

    Returns:
        (DataFrame, None) on success, or (None, error_message) on failure.
    """
    from blade.analysis._stability_tables import build_phase_data, detect_metals

    (
        system_dir,
        blade_dir,
        competitors_subdir,
        exclude_oxides,
        mu_o_values,
        phase_element,
        phase_stoich,
        t_values,
        composition_step,
        rk_order,
        output_dir,
        skip_existing,
        thresholds,
    ) = args

    if isinstance(competitors_subdir, str):  # a bare string is one subdir, not characters
        competitors_subdir = (competitors_subdir,) if competitors_subdir else None
    system_tag = system_dir.name
    out_path = output_dir / f"{system_tag}_stability.parquet"
    settings_sig = _settings_signature(
        competitors_subdir, exclude_oxides, mu_o_values, phase_element, phase_stoich, t_values, composition_step, rk_order
    ) + f"|thr={thresholds.parent_threshold:g},{thresholds.strict_dg_tol:g},{thresholds.min_metal_frac:g},{thresholds.lp_tol:g}"

    if skip_existing and out_path.exists():
        try:
            cached = pd.read_parquet(out_path)
            # Reuse only results produced with identical settings — otherwise toggling
            # exclude_oxides (or T range, grid, μ_O...) would silently return stale data.
            if "settings_sig" in cached.columns and (cached["settings_sig"] == settings_sig).all():
                return cached, None
        except Exception:
            pass

    competitors_dir: list[Path] | None = None
    if competitors_subdir:
        competitors_dir = [system_dir / sub for sub in competitors_subdir if (system_dir / sub).is_dir()] or None

    try:
        metals = detect_metals(blade_dir, phase_element)
        if len(metals) < 2:
            return None, f"{system_tag}: fewer than 2 metals detected, skipping"

        phase_data = build_phase_data(
            blade_dir,
            metals,
            phase_element,
            phase_stoich,
            rk_order=rk_order,
            y_step=composition_step,
            competitors_dir=competitors_dir,
            exclude_oxides=exclude_oxides,
        )
    except MissingEndpointError as e:
        return None, f"{system_tag}: {e}"
    except Exception as e:
        return None, f"{system_tag}: phase data build failed — {e}"

    try:
        df = _run_stability_grid(phase_data, t_values, system_tag, mu_o_values, thresholds)
        df["n_metals"] = len(metals)
        df["n_competitors"] = phase_data["n_competitors"]
        df["settings_sig"] = settings_sig
        df.to_parquet(out_path, index=False)
        print(
            f"  [{system_tag}] {len(metals)} metals | {phase_data['n_competitors']} competitors"
            f" | {len(df)} rows → {out_path.name}"
        )
        return df, None
    except Exception as e:
        return None, f"{system_tag}: LP grid failed — {e}"


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def write_pareto_database(frame: pd.DataFrame, db_dir: Path) -> tuple[Path, Path]:
    """Write ``stability_raw.parquet`` and ``stability_summary.parquet`` for Pareto analysis.

    Args:
        frame: Stability rows for a single μ_O (or the oxide-free scan), as produced
            by :func:`_run_stability_grid`.
        db_dir: Destination directory.

    Returns:
        Paths of the raw and summary parquet files.
    """
    raw = frame.rename(columns={"dG_hull": "dG_above_hull", "G_lp": "G_total"}).copy()
    raw_path = db_dir / "stability_raw.parquet"
    raw.to_parquet(raw_path, index=False)

    keys = ["system", "composition"]
    valid = raw.dropna(subset=["ss_fraction"])
    grouped = valid.groupby(keys)
    summary = grouped[["metals", "n_metals", "comp_label"] + [c for c in raw.columns if c.startswith("x_")]].first()
    for flag, tag in (
        ("is_parent_structure_retained", "retained"),
        ("is_single_phase_strict", "strict"),
        ("is_unique_homogeneous", "unique"),
    ):
        summary[f"T_min_{tag}"] = valid["T"].where(valid[flag]).groupby([valid[k] for k in keys]).min()
        summary[f"area_frac_{tag}"] = grouped[flag].mean()
    agg = raw.groupby(keys)
    summary["G_f_min"] = agg["G_formation"].min()
    summary["G_ss_f_min"] = agg["G_ss_formation"].min()
    summary["dG_hull_min"] = agg["dG_above_hull"].min()
    at_tmin = raw.merge(summary["T_min_retained"].rename("T").reset_index(), on=[*keys, "T"])
    summary["G_f_at_Tmin"] = at_tmin.groupby(keys)["G_formation"].first()
    summary_path = db_dir / "stability_summary.parquet"
    summary.reset_index().to_parquet(summary_path, index=False)
    return raw_path, summary_path


def _settings_signature(
    competitors_subdir: tuple[str, ...] | None,
    exclude_oxides: bool,
    mu_o_values: tuple[float, ...],
    phase_element: str,
    phase_stoich: float,
    t_values: np.ndarray,
    composition_step: float,
    rk_order: int,
) -> str:
    """Stable string identifying every setting that changes a system's results."""
    t = np.asarray(t_values, dtype=float)
    return (
        f"v4|comp={','.join(competitors_subdir or ())}|ox={'excl' if exclude_oxides else 'incl'}"
        f"|mu={','.join(f'{m:.6g}' for m in mu_o_values)}|el={phase_element}|st={phase_stoich:.6g}"
        f"|T={t.min():.6g}:{t.max():.6g}:{len(t)}|dy={composition_step:.6g}|rk={rk_order}"
    )


def _system_in_pool(system_name: str, allowed: set[str]) -> bool:
    """Return True if every element symbol in system_name is in allowed.

    Parses element symbols greedily (longest match first) so e.g. ``"CrHfB"``
    yields ``{"Cr", "Hf", "B"}``.
    """
    import re as _re
    # Sort allowed by length descending so two-char symbols match before one-char
    pat = _re.compile("|".join(sorted(allowed, key=len, reverse=True)))
    matched = set(pat.findall(system_name))
    # Reconstruct the name from matches; if leftovers exist, name has unknown elements
    remainder = system_name
    for sym in sorted(matched, key=len, reverse=True):
        remainder = remainder.replace(sym, "", 1)
    # Strip digits (stoichiometry numbers)
    remainder = _re.sub(r"\d", "", remainder)
    return len(remainder) == 0 and len(matched) > 0


def _count_valid_sqs(blade_dir: Path) -> int:
    """Count subdirs with both ``energy`` and ``CONTCAR`` present.

    Accepts both ``sqs_lev=*`` naming (BladeTDBGen output) and
    ``formula_entryid`` naming (database-stage output).

    Args:
        blade_dir: Path to ``<system_structures>/<SystemName>/blade/``.

    Returns:
        Number of valid structure subdirectories.
    """
    return sum(
        1
        for d in blade_dir.iterdir()
        if d.is_dir() and (d / "energy").exists() and (d / "CONTCAR").exists()
    )
