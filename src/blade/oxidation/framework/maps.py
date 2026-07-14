"""MapsRunner — standard analysis class (scan, muO-x, muO-T, fixed, ternary-3D).

run.py is used as an internal library only.  run_simple.py never imports it directly.
All skip checks and per-T caching live here.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np

from .config import Config
from .utils import csv_has_rows, system_key, system_tag


class MapsRunner:
    """Run standard phase-map analyses for one system.

    All configuration comes from a :class:`Config` instance.
    ``run.py`` is an internal implementation detail — callers only see this class.

    Usage::

        runner = MapsRunner(system_dir, metals, config)
        runner.run(sys_cfg, pd_data)
    """

    def __init__(self, system_dir: Path, metals: list[str], config: Config):
        self.system_dir = Path(system_dir)
        self.metals = metals
        self.config = config
        self.tag = system_tag(metals, config.phase_element)

    # ------------------------------------------------------------------ setup

    def _ensure_path(self) -> None:
        root = Path(__file__).parent.parent
        for p in [str(root), str(root / "python")]:
            if p not in sys.path:
                sys.path.insert(0, p)

    def _run_module(self):
        self._ensure_path()
        from . import _run as _r

        return _r

    def _push_globals(self) -> None:
        """Sync Config → run.py module globals so internal functions pick them up."""
        _r = self._run_module()
        cfg = self.config
        m = self.metals

        _r.DATA_ROOT = self.system_dir
        _r.TABLES_DIR = cfg.tables_dir
        _r.FIGURES_DIR = cfg.figures_dir
        _r.M1 = m[0]
        _r.M2 = m[1]
        _r.M3 = m[2] if len(m) >= 3 else None
        _r.RK_ORDER = cfg.rk_order
        _r.PHASE_ELEMENT = cfg.phase_element
        _r.PHASE_ELEMENT_STOICHIOMETRY = cfg.phase_element_stoichiometry
        _r.ALLOW_CALCULATIONS = cfg.run_calculations

        _r.RUN_OXIDATION_SCAN = cfg.run_oxidation_scan
        _r.RUN_MUO_X_MAP = cfg.run_muO_x_map
        _r.RUN_MUO_T_MAP = cfg.run_muO_T_map
        _r.RUN_FIXED_PHASE_MAP = cfg.run_fixed_phase_map
        _r.RUN_TERNARY_3D_MAP = cfg.run_ternary_3d_map

        _r.Y_STEP = cfg.y_step if len(m) < 3 else cfg.y_step_ternary
        _r.ACTIVE_THRESHOLD = cfg.active_threshold
        _r.PLOT_THRESHOLD = cfg.plot_threshold
        _r.INCLUDE_0P01_TO_0P05_COMPONENTS = cfg.include_0p01_to_0p05_components

        _r.SCAN_X = cfg.scan_x
        _r.SCAN_T = cfg.scan_T
        _r.SCAN_MU_O = cfg.scan_mu_O
        _r.SCAN_OXYGEN_MODE = cfg.scan_oxygen_mode
        _r.SCAN_WRITE_ALL_PHASES = cfg.scan_write_all_phases
        _r.SCAN_MAKE_PLOTS = cfg.run_plots and cfg.scan_make_plots
        _r.SCAN_PLOT_ONLY_ACTIVE = cfg.scan_plot_only_active

        _r.MAP_X_T = cfg.map_x_T
        _r.MAP_X_T_VALUES = list(cfg.map_x_T_values)
        _r.ANIM_FPS = cfg.anim_fps
        _r.MAP_X_X_VALUES = cfg.map_x_x_values
        _r.MAP_X_MU_O = cfg.map_x_mu_O
        _r.MAP_X_MAKE_PLOTS = cfg.run_plots and cfg.map_x_make_plots
        _r.MAP_X_WRITE_WIDE = cfg.map_x_write_wide

        _r.MAP_T_X = cfg.map_T_x
        _r.MAP_T_MU_O = cfg.map_T_mu_O
        _r.MAP_T_T_VALUES = cfg.map_T_T_values
        _r.MAP_T_MAKE_PLOTS = cfg.run_plots and cfg.map_T_make_plots
        _r.MAP_T_WRITE_WIDE = cfg.map_T_write_wide

        _r.FIXED_OXIDATION_THRESHOLD = cfg.fixed_oxidation_threshold
        _r.SKIP_IF_ANALYSIS_EXISTS = cfg.skip_if_analysis_exists
        _r.TERNARY_3D_COMP_STEP = cfg.ternary_3d_comp_step
        _r.TERNARY_3D_MU_O = cfg.ternary_3d_mu_O

        _r.PLOT_SCAN_NONPHASE = cfg.run_plots and cfg.plot_scan_nonphase
        _r.PLOT_SCAN_PHASE_Y = cfg.run_plots and cfg.plot_scan_phase_y
        _r.PLOT_MUO_X_REGION = cfg.run_plots and cfg.plot_muO_x_region
        _r.PLOT_MUO_T_REGION = cfg.run_plots and cfg.plot_muO_T_region
        _r.PLOT_MUO_T_SCALARS = cfg.run_plots and cfg.plot_muO_T_scalars
        _r.PLOT_FIXED_REGION = cfg.run_plots and cfg.plot_fixed_region
        _r.PLOT_FIXED_ONSET = cfg.run_plots and cfg.plot_fixed_onset
        _r.PLOT_3D_EQUILIBRIUM = cfg.run_plots and cfg.plot_3d_equilibrium
        _r.PLOT_3D_BOUNDARY = cfg.run_plots and cfg.plot_3d_boundary
        _r.PLOT_3D_ONSET_SURF = cfg.run_plots and cfg.plot_3d_onset_surf
        _r.PLOT_3D_TERNARY_DIAG = cfg.run_plots and cfg.plot_3d_ternary_diag

        _r.REGION_LABEL_MODE = cfg.region_label_mode
        _r.REGION_LABEL_FONTSIZE = cfg.region_label_fontsize
        _r.BOUNDARY_LW = cfg.boundary_lw

    # ------------------------------------------------------------------ skip helpers

    def _per_t_csv(self, analysis: str, T: int | float) -> Path:
        """Stable per-temperature CSV path for an analysis."""
        return self.config.tables_dir / f"{self.tag}_{analysis}_T{int(T)}_cache.csv"

    def _skip(self, csv_path: Path, expected_rows: int, check_component_threshold: bool = False) -> bool:
        if not (self.config.skip_if_analysis_exists and csv_has_rows(csv_path, expected_rows)):
            return False
        if not check_component_threshold:
            return True
        try:
            import pandas as pd

            values = pd.read_csv(csv_path, usecols=["component_presence_threshold"])[
                "component_presence_threshold"
            ].to_numpy(dtype=float)
            expected = 0.01 if self.config.include_0p01_to_0p05_components else 0.05
            return len(values) > 0 and np.allclose(values, expected)
        except Exception:
            return False

    def _muT_cache_matches(self, csv_path: Path, expected_rows: int, composition: float) -> bool:
        """Validate row count, initial composition, and labeling settings."""
        if not self._skip(csv_path, expected_rows, check_component_threshold=True):
            return False
        try:
            import pandas as pd

            column = f"x_{self.metals[0]}"
            df = pd.read_csv(csv_path, usecols=[column, "T_K", "muO_eV_per_O"])
            expected_T = np.repeat(
                np.asarray(self.config.map_T_T_values, dtype=float),
                len(self.config.map_T_mu_O),
            )
            expected_mu = np.tile(
                np.asarray(self.config.map_T_mu_O, dtype=float),
                len(self.config.map_T_T_values),
            )
            return (
                len(df) == expected_rows
                and np.allclose(df[column], float(composition), atol=1e-10, rtol=0.0)
                and np.allclose(df["T_K"], expected_T, atol=1e-10, rtol=0.0)
                and np.allclose(df["muO_eV_per_O"], expected_mu, atol=1e-10, rtol=0.0)
            )
        except Exception:
            return False

    def _scan_cache_matches(self, csv_path: Path) -> bool:
        if not self.config.skip_if_analysis_exists or not csv_path.exists():
            return False
        try:
            import pandas as pd

            cfg = self.config
            column = f"x_{self.metals[0]}"
            df = pd.read_csv(csv_path, usecols=[column, "T_K", "oxygen_axis"])
            temperatures = np.asarray(cfg.scan_T, dtype=float)
            oxygen = np.asarray(cfg.scan_mu_O, dtype=float)
            expected_T = np.repeat(temperatures, len(oxygen))
            expected_oxygen = np.tile(oxygen, len(temperatures))
            return (
                len(df) == len(expected_T)
                and np.allclose(df[column], cfg.scan_x, atol=1e-10, rtol=0.0)
                and np.allclose(df["T_K"], expected_T, atol=1e-10, rtol=0.0)
                and np.allclose(df["oxygen_axis"], expected_oxygen, atol=1e-10, rtol=0.0)
            )
        except Exception:
            return False

    def _muX_cache_matches(self, csv_path: Path, temperature: float) -> bool:
        expected_rows = len(self.config.map_x_x_values) * len(self.config.map_x_mu_O)
        if not self._skip(csv_path, expected_rows, check_component_threshold=True):
            return False
        try:
            import pandas as pd

            column = f"x_{self.metals[0]}"
            df = pd.read_csv(csv_path, usecols=[column, "T_K", "muO_eV_per_O"])
            expected_x = np.repeat(
                np.asarray(self.config.map_x_x_values, dtype=float),
                len(self.config.map_x_mu_O),
            )
            expected_mu = np.tile(
                np.asarray(self.config.map_x_mu_O, dtype=float),
                len(self.config.map_x_x_values),
            )
            return (
                len(df) == expected_rows
                and np.allclose(df[column], expected_x, atol=1e-10, rtol=0.0)
                and np.allclose(df["T_K"], temperature, atol=1e-10, rtol=0.0)
                and np.allclose(df["muO_eV_per_O"], expected_mu, atol=1e-10, rtol=0.0)
            )
        except Exception:
            return False

    # ------------------------------------------------------------------ analyses

    def run_scan(self, sys_cfg, pd_data) -> None:
        if not self.config.run_oxidation_scan:
            return
        _r = self._run_module()
        scan_csv = self.config.tables_dir / f"{self.tag}_scan_summary.csv"
        _r.SKIP_IF_ANALYSIS_EXISTS = self._scan_cache_matches(scan_csv)
        print("--- Oxidation scan ---")
        _r.run_oxidation_scan(sys_cfg, pd_data)  # skip logic inside _run.py
        _r.SKIP_IF_ANALYSIS_EXISTS = self.config.skip_if_analysis_exists

    def run_muO_T_map(self, sys_cfg, pd_data) -> None:
        if not self.config.run_muO_T_map:
            return
        _r = self._run_module()
        cfg = self.config
        M1 = self.metals[0]

        # x-sweep: loop over map_T_x_values if defined, else single map_T_x
        x_vals = getattr(cfg, "map_T_x_values", None)
        if x_vals is None or len(x_vals) == 0:
            x_vals = [cfg.map_T_x]

        import shutil as _sh

        sys_name = system_key(self.metals, cfg.phase_element)
        base_dir = cfg.figures_dir / sys_name / "muO_T_phase_map"
        base_dir.mkdir(parents=True, exist_ok=True)
        frames = []
        summary_csv = cfg.tables_dir / f"{self.tag}_muO_T_phase_map_summary.csv"
        n_expected = len(cfg.map_T_T_values) * len(cfg.map_T_mu_O)

        for x in x_vals:
            _r.MAP_T_X = round(float(x), 6)
            x_str = f"x{round(float(x), 4):.2f}".replace(".", "p")
            x_dir = base_dir / x_str
            per_x_cache = cfg.tables_dir / f"{self.tag}_muO_T_phase_map_{x_str}_cache.csv"

            # Migrate the legacy shared summary when it happens to contain
            # this composition. Other compositions calculate once and then
            # receive their own durable cache.
            if not per_x_cache.exists() and self._muT_cache_matches(summary_csv, n_expected, x):
                _sh.copy(summary_csv, per_x_cache)

            if self._muT_cache_matches(per_x_cache, n_expected, x):
                print(f"\n--- muO-T map (x_{M1}={x:.2f}): " "data complete — skipping LP ---")
                _sh.copy(per_x_cache, summary_csv)
                _r.SKIP_IF_ANALYSIS_EXISTS = True
            else:
                print(f"\n--- muO-T map (x_{M1}={x:.2f}) ---")
                _r.SKIP_IF_ANALYSIS_EXISTS = False
            _r.run_muO_T_map(sys_cfg, pd_data)
            if summary_csv.exists():
                _sh.copy(summary_csv, per_x_cache)
            _r.SKIP_IF_ANALYSIS_EXISTS = cfg.skip_if_analysis_exists

            if cfg.run_plots:
                # Move PNG + HTML from base into the per-composition folder.
                x_dir.mkdir(parents=True, exist_ok=True)
                for fname in ["assemblage_region_map.png", "assemblage_region_map.html"]:
                    src = base_dir / fname
                    dst = x_dir / fname
                    if src.exists():
                        _sh.move(str(src), str(dst))

                frame = x_dir / "assemblage_region_map.png"
                if frame.exists():
                    frames.append(frame)

        # Compile animation in the parent muO_T_phase_map/ folder
        if cfg.run_plots and cfg.run_animations and len(frames) > 1:
            from .utils import make_animation

            make_animation(
                frames,
                base_dir / "muO_T_x_sweep.gif",
                base_dir / "muO_T_x_sweep.mp4",
                fps=cfg.anim_fps,
                mp4_crf=getattr(cfg, "mp4_crf", 18),
                mp4_preset=getattr(cfg, "mp4_preset", "slow"),
            )
            print(f"  muO-T sweep: {len(frames)} frames → muO_T_x_sweep.gif + .mp4")

    def run_muO_x_map(self, sys_cfg, pd_data) -> None:
        """Run muO-x map for each T, skipping temperatures whose cache CSV is complete."""
        if not self.config.run_muO_x_map:
            return
        _r = self._run_module()
        cfg = self.config
        summary_csv = cfg.tables_dir / f"{self.tag}_muO_x_phase_map_summary.csv"

        for T in cfg.map_x_T_values:
            _r.MAP_X_T = T
            per_t = self._per_t_csv("muO_x_phase_map", T)

            if self._muX_cache_matches(per_t, T):
                print(f"\n--- muO-x map (T={T} K): data complete — skipping LP ---")
                # Restore summary CSV for this T so figure generation still works
                shutil.copy(per_t, summary_csv)
                _r.SKIP_IF_ANALYSIS_EXISTS = True  # figures-only mode
                _r.run_muO_x_map(sys_cfg, pd_data)
                _r.SKIP_IF_ANALYSIS_EXISTS = cfg.skip_if_analysis_exists
            else:
                print(f"\n--- muO-x map (T={T} K) ---")
                _r.SKIP_IF_ANALYSIS_EXISTS = False
                _r.run_muO_x_map(sys_cfg, pd_data)
                # Cache per-T result immediately
                if summary_csv.exists():
                    shutil.copy(summary_csv, per_t)
        _r.SKIP_IF_ANALYSIS_EXISTS = cfg.skip_if_analysis_exists

    def run_fixed_phase(self, sys_cfg, pd_data) -> None:
        if not self.config.run_fixed_phase_map:
            return
        _r = self._run_module()
        cfg = self.config
        onset_csv = cfg.tables_dir / f"{self.tag}_fixed_phase_muO_x_phase_map_oxidation_onset.csv"
        n_expected = len(cfg.map_x_x_values)
        if self._skip(onset_csv, n_expected):
            print("  Fixed-composition map: data complete — skipping LP")
            _r.SKIP_IF_ANALYSIS_EXISTS = True
        else:
            print("\n--- Fixed-composition control map ---")
            _r.SKIP_IF_ANALYSIS_EXISTS = False
        _r.run_fixed_phase_map(sys_cfg, pd_data)
        _r.SKIP_IF_ANALYSIS_EXISTS = cfg.skip_if_analysis_exists

    def run_ternary_3d(self, sys_cfg, pd_data) -> None:
        if not self.config.run_ternary_3d_map or len(self.metals) != 3:
            return
        _r = self._run_module()
        cfg = self.config
        for T in cfg.map_x_T_values:
            _r.MAP_X_T = T
            print(f"\n--- Ternary 3D (T={T} K) ---")
            _r.SKIP_IF_ANALYSIS_EXISTS = cfg.skip_if_analysis_exists
            _r.run_ternary_3d_map(sys_cfg, pd_data)

    def compile_animations(self, sys_cfg) -> None:
        if self.config.run_plots and self.config.run_animations and len(self.config.map_x_T_values) > 1:
            _r = self._run_module()
            _r._compile_system_animations(sys_cfg)

    # ------------------------------------------------------------------ main

    def run(self, sys_cfg, pd_data) -> None:
        """Run all enabled analyses with per-T skip checks."""
        self._push_globals()
        self.run_scan(sys_cfg, pd_data)
        self.run_muO_T_map(sys_cfg, pd_data)
        self.run_muO_x_map(sys_cfg, pd_data)
        self.run_fixed_phase(sys_cfg, pd_data)
        self.run_ternary_3d(sys_cfg, pd_data)
        self.compile_animations(sys_cfg)
