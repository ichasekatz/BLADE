"""Pure solid-state thermodynamic stability analysis for BLADE diboride systems.

Computes equilibrium phase fractions and ΔG from the convex hull over the
Redlich-Kister SS phase + fixed SQS-structure competitors, with no oxygen
chemical potential. Outputs per-system parquet files and a combined summary.
"""

from __future__ import annotations

import multiprocessing as mp
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from blade.oxidation.framework._lp import solve_grand_lp

from ._stability_tables import MissingEndpointError, build_phase_data, detect_metals

gas_constant = 8.314472  # J / (mol·K)
kb_ev_per_k = 8.617333262145e-5  # eV / K


class StabilityConfigError(ValueError):
    """Raised when a required configuration value is missing or invalid."""


class BladeStability:
    """Pure solid-state phase stability scanner for BLADE Comps outputs.

    Reads GRACE-relaxed SQS structures from each system's Comps lattice directory,
    fits the Redlich-Kister/Muggianu interaction model, and runs the grand-potential
    LP (without oxygen) over a composition simplex grid × temperature range.

    For each (composition, T) point the LP finds the minimum free energy over all
    competing phases (fixed SQS structures + continuous SS grid). The SS-phase
    stability is reported as:

        dG_hull = G_ss(y, T) − G_lp(y, T)   [eV per formula unit]

    where G_lp ≤ G_ss always. dG_hull ≤ 0 means the SS phase is thermodynamically
    stable (or tied) at that composition and temperature.

    Args:
        files_dir: Path to the BLADE ``Files/`` directory.
        phase_element: Non-metal sublattice element (e.g. ``"B"``).
        phase_stoich: Stoichiometry per formula unit (e.g. ``2.0`` for AlB₂-type).
        comps_folder: Subfolder of ``files_dir`` holding TDB-fitting outputs.
        output_folder: Subfolder of ``files_dir`` for stability results.
        t_min: Minimum temperature in K (default 300).
        t_max: Maximum temperature in K (default 2000).
        t_step: Temperature step in K (default 100).
        composition_step: Simplex grid spacing (default 0.02).
        skip_existing: Skip systems whose parquet already exists (default True).
        rk_order: Redlich-Kister polynomial order (default 3).
        n_workers: Number of parallel worker processes (default: CPU count).
    """

    def __init__(
        self,
        files_dir: str | Path,
        phase_element: str = "B",
        phase_stoich: float = 2.0,
        comps_folder: str = "Comps",
        output_folder: str = "Stability",
        t_min: float = 300.0,
        t_max: float = 2000.0,
        t_step: float = 100.0,
        composition_step: float = 0.02,
        skip_existing: bool = True,
        rk_order: int = 3,
        n_workers: int | None = None,
    ) -> None:
        """Initialize with file paths and stability scan parameters."""
        self._files_dir = Path(files_dir)
        self.phase_element = phase_element
        self.phase_stoich = phase_stoich
        self._comps_folder = comps_folder
        self._output_folder = output_folder
        self.t_min = t_min
        self.t_max = t_max
        self.t_step = t_step
        self.composition_step = composition_step
        self.skip_existing = skip_existing
        self.rk_order = rk_order
        self.n_workers = n_workers or mp.cpu_count()

    @cached_property
    def comps_dir(self) -> Path:
        """Absolute path to the Comps directory."""
        return self._files_dir / self._comps_folder

    @cached_property
    def output_dir(self) -> Path:
        """Absolute path to the Stability output directory (created on access)."""
        d = self._files_dir / self._output_folder
        d.mkdir(parents=True, exist_ok=True)
        return d

    @cached_property
    def t_values(self) -> np.ndarray:
        """Temperature grid in K from t_min to t_max (inclusive) at t_step intervals."""
        return np.arange(self.t_min, self.t_max + 1e-9, self.t_step)

    def run(self) -> Path:
        """Run stability analysis for all systems found in the Comps directory.

        Returns:
            Path to the combined ``stability_results.parquet`` file.
        """
        system_dirs = self._find_system_dirs()
        print(f"BladeStability: found {len(system_dirs)} systems in {self.comps_dir}")

        args = [
            (
                sdir, ldir, self.phase_element, self.phase_stoich,
                self.t_values, self.composition_step, self.rk_order,
                self.output_dir, self.skip_existing,
            )
            for sdir, ldir in system_dirs
        ]

        all_frames: list[pd.DataFrame] = []
        with mp.Pool(processes=self.n_workers, initializer=_worker_init) as pool:
            for result, err in pool.imap_unordered(_stability_worker, args):
                if err:
                    print(f"  SKIPPED: {err}")
                elif result is not None:
                    all_frames.append(result)

        if not all_frames:
            print("No stability results produced.")
            combined = pd.DataFrame()
        else:
            combined = pd.concat(all_frames, ignore_index=True)
            combined = combined.sort_values(["system", "composition", "T"]).reset_index(drop=True)

        out = self.output_dir / "stability_results.parquet"
        combined.to_parquet(out, index=False)
        print(f"Combined results → {out}  ({len(combined)} rows, {combined['system'].nunique() if len(combined) else 0} systems)")
        return out

    def _find_system_dirs(self) -> list[tuple[Path, Path]]:
        """Return (system_dir, lattice_dir) pairs for all valid Comps systems.

        A valid system has at least one lattice subdir containing ≥ 3 ``sqs_lev=*``
        subdirs with both ``energy`` and ``CONTCAR`` present.

        Returns:
            Sorted list of (system_dir, lattice_dir) tuples.
        """
        result: list[tuple[Path, Path]] = []
        for sdir in sorted(self.comps_dir.iterdir()):
            if not sdir.is_dir():
                continue
            ldir = _find_lattice_dir(sdir)
            if ldir is not None:
                result.append((sdir, ldir))
        return result


# ---------------------------------------------------------------------------
# LP grid runner (pure stability, no oxygen)
# ---------------------------------------------------------------------------


def _run_stability_grid(
    phase_data: dict,
    t_values: np.ndarray,
    system_tag: str,
) -> pd.DataFrame:
    """Run the stability LP over a composition simplex grid × temperature range.

    For each composition y_i and temperature T, solves the LP:

        min  Σ_j G_j(T) · n_j
        s.t. A_eq @ n = b_eq(y_i),  n ≥ 0

    where G_j(T) = energy_per_formula for fixed phases and G_j(T) = H_RK(y_j) +
    k_B T · Σ y_k ln(y_k) for the SS phase grid. No oxide phases or μ_O term.

    Args:
        phase_data: Output of :func:`_stability_tables.build_phase_data`.
        t_values: Temperature array in K.
        system_tag: System identifier string for the ``system`` column.

    Returns:
        DataFrame with columns: system, composition, x_<Metal>, T, dG_hull,
        G_ss, G_lp, stable, s_config, s_config_norm, dominant_phase.
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
    n = len(metals)

    rows: list[dict] = []
    for t in t_values:
        g_ss_grid = phase_h0 + kb_ev_per_k * float(t) * phase_mix
        grand = np.concatenate([fef, g_ss_grid])

        for i, y_i in enumerate(y_grid):
            b_eq = np.concatenate([y_i, [phase_stoich]])
            amounts, g_lp, ok = solve_grand_lp(a_eq, b_eq, grand)
            if not ok:
                continue

            g_ss_i = grand[n_fixed + i]
            dg = g_ss_i - g_lp

            fracs = y_i[y_i > 1e-9]
            s_cfg = -gas_constant * float(np.sum(fracs * np.log(fracs)))
            s_norm = s_cfg / (-gas_constant * np.log(1.0 / n)) if n > 1 else 0.0

            dom_idx = int(np.nanargmax(amounts)) if amounts is not None else -1
            dom_phase = str(phase_ids[dom_idx]) if dom_idx >= 0 else ""

            comp_str = ",".join(f"{y_i[k]:.4f}" for k in range(n))
            row: dict = {
                "system": system_tag,
                "composition": comp_str,
                "T": float(t),
                "dG_hull": float(dg),
                "G_ss": float(g_ss_i),
                "G_lp": float(g_lp),
                "stable": bool(dg <= 1e-6),
                "s_config": float(s_cfg),
                "s_config_norm": float(s_norm),
                "dominant_phase": dom_phase,
            }
            for k, m in enumerate(metals):
                row[f"x_{m}"] = float(y_i[k])
            rows.append(row)

    return pd.DataFrame(rows)


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
        args: Tuple of (system_dir, lattice_dir, phase_element, phase_stoich,
            t_values, composition_step, rk_order, output_dir, skip_existing).

    Returns:
        (DataFrame, None) on success, or (None, error_message) on failure.
    """
    from blade.analysis._stability_tables import build_phase_data, detect_metals

    system_dir, lattice_dir, phase_element, phase_stoich, t_values, composition_step, rk_order, output_dir, skip_existing = args

    system_tag = system_dir.name
    out_path = output_dir / f"{system_tag}_stability.parquet"

    if skip_existing and out_path.exists():
        try:
            return pd.read_parquet(out_path), None
        except Exception:
            pass

    try:
        metals = detect_metals(lattice_dir)
        if len(metals) < 2:
            return None, f"{system_tag}: fewer than 2 metals detected, skipping"

        phase_data = build_phase_data(
            lattice_dir, metals, phase_element, phase_stoich,
            rk_order=rk_order, y_step=composition_step,
        )
    except MissingEndpointError as e:
        return None, f"{system_tag}: {e}"
    except Exception as e:
        return None, f"{system_tag}: phase data build failed — {e}"

    try:
        df = _run_stability_grid(phase_data, t_values, system_tag)
        df["n_metals"] = len(metals)
        df.to_parquet(out_path, index=False)
        print(f"  [{system_tag}] {len(metals)} metals | {len(df)} rows → {out_path.name}")
        return df, None
    except Exception as e:
        return None, f"{system_tag}: LP grid failed — {e}"


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def _find_lattice_dir(system_dir: Path) -> Path | None:
    """Find the first valid lattice subdir within a system Comps directory.

    A valid lattice dir contains at least 3 ``sqs_lev=*`` subdirs each
    with both ``energy`` and ``CONTCAR`` present.

    Args:
        system_dir: Path to ``<Comps>/<SystemName>/``.

    Returns:
        Path to the lattice dir, or None if none qualifies.
    """
    for d in sorted(system_dir.iterdir()):
        if not d.is_dir():
            continue
        valid = sum(
            1
            for sub in d.iterdir()
            if sub.is_dir() and sub.name.startswith("sqs_lev=")
            and (sub / "energy").exists() and (sub / "CONTCAR").exists()
        )
        if valid >= 3:
            return d
    return None
