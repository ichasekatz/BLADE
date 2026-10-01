"""SingleCompositionAnalyzer — all graphs for one specific composition.

Usage in run_simple.py or a notebook::

    from framework import SingleCompositionAnalyzer
    import numpy as np

    ana = SingleCompositionAnalyzer(
        system_dir = Path("/path/to/CrHfB"),
        tables_dir = Path("tables"),
        figures_dir = Path("figures"),
        metals     = ["Cr", "Hf"],          # auto-detected if None
        composition = 0.3,                  # binary: x_M1 scalar
        #composition = [0.3, 0.5, 0.2],    # ternary: [x_M1, x_M2, x_M3]
    )
    ana.run(
        T_values   = np.arange(200, 2001, 10),
        mu_O_values = np.arange(-10, -4, 0.1),
        scan_T     = [700, 1000, 1300],
    )
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._compute import (
    _amounts_from_fractions,
    _conservation_rhs,
    run_muO_T_map,
    run_muO_x_map,
    run_scan,
)
from ._compute import (
    _family_threshold as _family_threshold_fn,
)
from ._plots import _plot_muO_T
from ._prepare import _detect_metals, _ensure_path, load, prepare
from .utils import prepare_system_tables_dir, system_key


class SingleCompositionAnalyzer:
    """Run all standard analyses for one specific composition.

    Parameters
    ----------
    system_dir   : raw GRACE calculation folder (blade/ + grace/ subdirs)
    tables_dir   : root containing one CSV folder per modeled system
    figures_dir  : where figures are saved
    metals       : element list e.g. ['Cr','Hf']. Auto-detected if None.
    composition  : scalar x_M1 for binary; list [x_M1,x_M2,...] for n-ary.
                   Fractions must sum to 1.
    rk_order     : RK polynomial order for table preparation
    y_step       : flexible-phase composition grid step
    mixed_phase_subdir  : subfolder name for SQS/GRACE outputs
    fixed_phases_subdir : subfolder name for reference phase calculations
    """

    # Bind sibling-module functions as methods.
    _ensure_path = _ensure_path
    _detect_metals = _detect_metals
    prepare = prepare
    load = load
    _conservation_rhs = _conservation_rhs
    _amounts_from_fractions = _amounts_from_fractions
    run_scan = run_scan
    run_muO_T_map = run_muO_T_map
    run_muO_x_map = run_muO_x_map
    _plot_muO_T = _plot_muO_T

    @property
    def _family_threshold(self) -> float:
        return _family_threshold_fn(self)

    def __init__(
        self,
        system_dir: str | Path | None = None,
        tables_dir: str | Path = "tables",
        figures_dir: str | Path = "figures",
        metals: list[str] | None = None,
        composition: float | list[float] = 0.5,
        rk_order: int = 3,
        y_step: float = 0.01,
        mixed_phase_subdir: str = "blade",
        fixed_phases_subdir: str = "ORB",
        phase_element: str | None = None,
        phase_element_stoichiometry: float = 0.0,
        include_0p01_to_0p05_components: bool = True,
        region_label_mode: str = "id",
        region_label_fontsize: int = 12,
    ):
        self.system_dir = Path(system_dir) if system_dir else None
        self.tables_root = Path(tables_dir)
        self.figures_dir = Path(figures_dir)
        self.rk_order = rk_order
        self.y_step = y_step
        self.mixed_phase_subdir = mixed_phase_subdir
        self.fixed_phases_subdir = fixed_phases_subdir
        self.phase_element = str(phase_element).strip() if phase_element else None
        self.phase_element_stoichiometry = float(phase_element_stoichiometry) if self.phase_element else 0.0
        self.include_0p01_to_0p05_components = bool(include_0p01_to_0p05_components)
        mode = str(region_label_mode).strip().lower()
        if mode in {"id", "ids", "number", "numbers"}:
            self.region_label_mode = "id"
        elif mode in {"phase", "phases", "full"}:
            self.region_label_mode = "phases"
        else:
            raise ValueError("region_label_mode must be 'id' or 'phases'")
        self.region_label_fontsize = int(region_label_fontsize)
        if self.region_label_fontsize <= 0:
            raise ValueError("region_label_fontsize must be positive")
        if self.phase_element and self.phase_element_stoichiometry <= 0:
            raise ValueError("phase_element_stoichiometry must be positive when phase_element is set")

        self._ensure_path()

        # Resolve metals
        if metals is None:
            metals = self._detect_metals()
        self.metals = metals
        self.n_metals = len(metals)
        self.tables_dir = prepare_system_tables_dir(self.tables_root, self.metals, self.phase_element)

        # Normalise composition
        if np.isscalar(composition):
            x = float(composition)
            if self.n_metals == 2:
                self._comp = np.array([x, 1.0 - x])
            else:
                rest = (1.0 - x) / (self.n_metals - 1)
                self._comp = np.array([x] + [rest] * (self.n_metals - 1))
        else:
            self._comp = np.asarray(composition, dtype=float)
            if abs(self._comp.sum() - 1.0) > 1e-6:
                raise ValueError("Composition fractions must sum to 1.0")

        self._x_label = "  ".join(f"x_{m}={self._comp[i]:.3g}" for i, m in enumerate(self.metals))
        self._x_str = "_".join(f"{m}{self._comp[i]:.2f}".replace(".", "p") for i, m in enumerate(self.metals))

        self._sys_cfg = None
        self._pd_data = None
        self._run_calculations = True
        self._run_plots = True

    # ---------------------------------------------------------------- properties

    @property
    def pd_data(self) -> dict:
        """Phase-data dictionary; loaded on first access."""
        if self._pd_data is None:
            self.load()
        return self._pd_data

    @property
    def sys_cfg(self):
        """System configuration; resolved on first access."""
        if self._sys_cfg is None:
            self.prepare()
        return self._sys_cfg

    @property
    def out_dir(self) -> Path:
        """Per-composition output directory (created on access)."""
        sys_name = system_key(self.metals, self.phase_element)
        d = self.figures_dir / sys_name / "single_composition" / self._x_str
        d.mkdir(parents=True, exist_ok=True)
        return d

    # ---------------------------------------------------------------- run all

    def run(
        self,
        T_values: np.ndarray | None = None,
        mu_O_values: np.ndarray | None = None,
        scan_T: list[float] | None = None,
        x_values: np.ndarray | None = None,
        skip_if_exists: bool = True,
        run_calculations: bool = True,
        run_plots: bool = True,
    ) -> Path:
        """Run all analyses for this composition. Returns output directory.

        When both stage flags are true, all numerical caches are completed
        before a second, cache-only plotting pass starts.

        Args:
            T_values: Temperature grid (K). Defaults to np.arange(200, 2001, 10).
            mu_O_values: μO grid (eV/O). Defaults to np.arange(-10, -4, 0.1).
            scan_T: Temperatures for 1-D scan plots. Defaults to grid midpoint.
            x_values: Composition grid for μO-x map. Defaults to np.arange(0, 1.01, 0.01).
            skip_if_exists: Re-use on-disk caches when available.
            run_calculations: Whether to run LP solves.
            run_plots: Whether to render figures.

        Returns:
            Path to the per-composition output directory.
        """
        import traceback

        T_vals = np.asarray(T_values if T_values is not None else np.arange(200, 2001, 10))
        mu_vals = np.asarray(mu_O_values if mu_O_values is not None else np.arange(-10, -4, 0.1))
        s_T = scan_T if scan_T is not None else [int(T_vals[len(T_vals) // 2])]
        x_vals = np.asarray(x_values if x_values is not None else np.arange(0.0, 1.01, 0.01))

        self.prepare()
        self.load()

        def execute_pass(allow_calculations: bool, allow_plots: bool, use_cache: bool) -> None:
            self._run_calculations = allow_calculations
            self._run_plots = allow_plots
            failures = []
            for label, fn in [
                ("scan", lambda: self.run_scan(s_T, mu_vals, use_cache)),
                ("muO-T map", lambda: self.run_muO_T_map(T_vals, mu_vals, use_cache)),
                ("muO-x map", lambda: [self.run_muO_x_map(T, x_vals, mu_vals, use_cache) for T in s_T]),
            ]:
                try:
                    fn()
                except Exception as e:
                    print(f"  [{label}] FAILED: {e}")
                    traceback.print_exc()
                    failures.append(f"{label}: {e}")
            if failures:
                raise RuntimeError("; ".join(failures))

        if run_calculations and run_plots:
            print("\n=== Single-composition calculation pass ===")
            execute_pass(True, False, skip_if_exists)
            print("\n=== Single-composition plotting pass (cached data only) ===")
            execute_pass(False, True, True)
        elif run_calculations or run_plots:
            execute_pass(run_calculations, run_plots, True if run_plots else skip_if_exists)
        else:
            print("Single-composition oxidation calculations and plots are both disabled.")

        if run_plots:
            print(f"\n  All figures → {self.out_dir}")
        return self.out_dir
