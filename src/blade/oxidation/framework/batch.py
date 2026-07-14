"""BatchConfig and BatchRunner — public API for run_simple.py.

All analysis logic lives in framework/ (replicated from nary/ and python/).
BatchConfig  : user-facing settings matching run_simple.py field names.
BatchRunner  : translates BatchConfig → internal Config, calls _batch_internal.BatchRunner.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

_FW = Path(__file__).parent
if str(_FW) not in sys.path:
    sys.path.insert(0, str(_FW))


@dataclass
class BatchConfig:
    # ---- flexible-phase chemistry -------------------------------------------
    phase_element: str | None = None
    phase_element_stoichiometry: float = 0.0

    # ---- paths ---------------------------------------------------------------
    system_structures_root: Path = Path("system_structures")
    tables_dir: Path = Path("tables")
    figures_dir: Path = Path("figures")
    phase_diagrams_root: Path = Path("Phase_Diagrams")

    # ---- subfolder names inside each system directory ------------------------
    mixed_phase_subdir: str = "blade"
    fixed_phases_subdir: str = "ORB"

    # ---- which analyses to run -----------------------------------------------
    run_calculations: bool = True
    run_plots: bool = True
    run_scan: bool = True
    run_muO_x_map: bool = True
    run_onset: bool = True
    run_3d_onset: bool = False
    run_composition_slices: bool = True
    run_composition_slice_muT: bool = True
    run_animations: bool = True
    calculation_first: bool = True

    # ---- skip / filter -------------------------------------------------------
    skip_if_tables_exist: bool = True
    skip_if_analysis_exists: bool = True
    skip_systems: set = field(default_factory=set)
    # Exact system directory names. Empty includes every matching system.
    systems: list[str] = field(default_factory=list)
    # Empty includes all systems; otherwise systems may only contain these metals.
    elements: list[str] = field(default_factory=list)

    # ---- shared grids --------------------------------------------------------
    temperature_values: object = field(default_factory=lambda: np.arange(200, 2001, 200))
    mu_O_values: object = field(default_factory=lambda: np.arange(-10.0, -4.0 + 0.001, 0.10))

    # ---- RK fit / phase grid -------------------------------------------------
    rk_order: int = 3
    phase_grid_step: float = 0.01
    phase_grid_step_ternary: float = 0.05
    active_threshold: float = 1e-9
    include_0p01_to_0p05_components: bool = True

    # ---- region-map annotations ---------------------------------------------
    region_label_mode: str = "id"  # "id" or "phases"
    region_label_fontsize: int = 12

    # ---- composition slice paths --------------------------------------------
    slice_axis: int | str | None = None
    slice_axis_priority: list[str] = field(default_factory=list)
    slice_remainder_ratios: list = field(default_factory=lambda: [0.0, 0.25, 1 / 3, 0.5, 2 / 3, 0.75, 1.0])
    slice_muT_comp_step: float = 0.10

    # ---- 1D scan -------------------------------------------------------------
    scan_x: float = 0.50
    scan_T: list = field(default_factory=lambda: [700])
    scan_mu_O: object = field(default_factory=lambda: np.arange(-10.0, -4.0 + 0.001, 0.05))

    # ---- muO-x map -----------------------------------------------------------
    map_x_x_values: object = field(default_factory=lambda: np.arange(0.0, 1.01, 0.05))
    anim_fps: int = 2

    # ---- onset diagrams ------------------------------------------------------
    onset_comp_step_binary: float = 0.02
    onset_comp_step_ternary: float = 0.025
    onset_comp_step_nary: float = 0.10
    onset_threshold: float = 1e-8

    # ---- phase diagram image T range (for indicator line in side-by-side) ----
    phase_diagram_t_start: int = 300
    phase_diagram_t_end: int = 4500

    def __post_init__(self):
        self.system_structures_root = Path(self.system_structures_root)
        self.tables_dir = Path(self.tables_dir)
        self.figures_dir = Path(self.figures_dir)
        self.phase_diagrams_root = Path(self.phase_diagrams_root)
        self.elements = list(
            dict.fromkeys(str(element).strip().title() for element in self.elements if str(element).strip())
        )
        self.systems = list(dict.fromkeys(str(system).strip() for system in self.systems if str(system).strip()))
        self.slice_muT_comp_step = float(self.slice_muT_comp_step)
        if not 0.0 < self.slice_muT_comp_step <= 1.0:
            raise ValueError("slice_muT_comp_step must be greater than 0 and at most 1")

    def _to_internal_config(self):
        from .config import Config

        return Config(
            phase_element=self.phase_element,
            phase_element_stoichiometry=self.phase_element_stoichiometry,
            system_structures_root=self.system_structures_root,
            blade_subdir=self.mixed_phase_subdir,
            fixed_phases_subdir=self.fixed_phases_subdir,
            tables_dir=self.tables_dir,
            figures_dir=self.figures_dir,
            ternary_diagrams_root=self.phase_diagrams_root,
            run_calculations=self.run_calculations,
            run_plots=self.run_plots,
            run_oxidation_scan=self.run_scan,
            run_muO_x_map=self.run_muO_x_map,
            run_muO_T_map=False,
            run_fixed_phase_map=False,
            run_ternary_3d_map=self.run_3d_onset,
            run_composition_slice_maps=self.run_composition_slices,
            run_composition_slice_muT_maps=self.run_composition_slice_muT,
            run_onset_auc_diagrams=self.run_onset,
            run_animations=self.run_animations,
            skip_if_tables_exist=self.skip_if_tables_exist,
            skip_if_analysis_exists=self.skip_if_analysis_exists,
            skip_systems=self.skip_systems,
            systems=self.systems,
            elements=self.elements,
            rk_order=self.rk_order,
            y_step=self.phase_grid_step,
            y_step_ternary=self.phase_grid_step_ternary,
            active_threshold=self.active_threshold,
            include_0p01_to_0p05_components=self.include_0p01_to_0p05_components,
            region_label_mode=self.region_label_mode,
            region_label_fontsize=self.region_label_fontsize,
            slice_remainder_ratios=self.slice_remainder_ratios,
            slice_axis=self.slice_axis,
            slice_axis_priority=self.slice_axis_priority,
            slice_muT_x_values=np.arange(
                0.0,
                1.0 + self.slice_muT_comp_step / 2.0,
                self.slice_muT_comp_step,
            ),
            temperature_values=self.temperature_values,
            mu_O_values=self.mu_O_values,
            scan_x=self.scan_x,
            scan_T=self.scan_T,
            scan_mu_O=self.scan_mu_O,
            map_x_x_values=self.map_x_x_values,
            anim_fps=self.anim_fps,
            compile_slice_video_grids=self.run_animations,
            onset_comp_step_binary=self.onset_comp_step_binary,
            onset_comp_step_ternary=self.onset_comp_step_ternary,
            onset_comp_step_higher=self.onset_comp_step_nary,
            onset_threshold=self.onset_threshold,
            make_side_by_side_miscibility=True,
            phase_diagram_t_start=self.phase_diagram_t_start,
            phase_diagram_t_end=self.phase_diagram_t_end,
        )


class BatchRunner:
    """Translate BatchConfig → internal Config and run the full batch pipeline."""

    def __init__(self, config: BatchConfig):
        self.config = config

    def run(self) -> list:
        from ._batch_internal import BatchRunner as _InternalRunner

        if not self.config.run_calculations and not self.config.run_plots:
            print("Oxidation calculations and plots are both disabled.")
            return []

        if self.config.run_calculations and self.config.run_plots and self.config.calculation_first:
            print("\n=== Oxidation calculation pass ===")
            calculation_config = replace(
                self.config,
                run_plots=False,
                run_animations=False,
            )
            _InternalRunner(calculation_config._to_internal_config()).run()

            print("\n=== Oxidation plotting pass (cached data only) ===")
            plotting_config = replace(
                self.config,
                run_calculations=False,
                skip_if_tables_exist=True,
                skip_if_analysis_exists=True,
            )
            return _InternalRunner(plotting_config._to_internal_config()).run()

        return _InternalRunner(self.config._to_internal_config()).run()


# (appended — handled in _to_internal_config via config defaults)
