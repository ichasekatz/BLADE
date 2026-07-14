"""All settings for the N-ary phase oxidation analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Config:
    # Element on the non-metal sublattice of the flexible phase.  None models
    # a substitutional phase containing only the configured metals.
    phase_element: str | None = None
    phase_element_stoichiometry: float = 0.0
    # ------------------------------------------------------------------ Paths
    system_structures_root: Path = Path("system_structures")
    blade_subdir: str = "blade"
    fixed_phases_subdir: str = "fixed_phases"
    tables_dir: Path = Path("tables")
    figures_dir: Path = Path("figures")
    # TDB files for pycalphad miscibility gap calculations
    comps_dir: Path = Path("Comps")
    # External phase diagrams to place beside slice maps
    ternary_diagrams_root: Path = Path("Phase_Diagrams")

    # ----------------------------------------------------------- Which to run
    # When both are enabled, the public batch runner completes all calculations
    # before starting a second, cache-only plotting pass.
    run_calculations: bool = True
    run_plots: bool = True
    run_oxidation_scan: bool = True
    run_muO_x_map: bool = True
    run_muO_T_map: bool = False
    run_fixed_phase_map: bool = False
    run_ternary_3d_map: bool = False
    run_composition_slice_maps: bool = True
    run_composition_slice_muT_maps: bool = True
    run_onset_auc_diagrams: bool = True
    run_animations: bool = True

    # ---------------------------------------------------------- Per-plot flags
    plot_scan_nonphase: bool = True
    plot_scan_phase_y: bool = True
    plot_muO_x_region: bool = True
    plot_muO_T_region: bool = True
    plot_muO_T_scalars: bool = True
    plot_fixed_region: bool = False
    plot_fixed_onset: bool = False
    plot_3d_equilibrium: bool = True
    plot_3d_boundary: bool = True
    plot_3d_onset_surf: bool = True
    plot_3d_ternary_diag: bool = True

    # ------------------------------------------------------ Skip / filter
    skip_if_tables_exist: bool = True
    skip_if_analysis_exists: bool = True
    filter_no_miscibility_gap: bool = False
    skip_systems: set = field(default_factory=set)
    systems: list[str] = field(default_factory=list)
    # Empty includes all systems; otherwise systems may only contain these metals.
    elements: list[str] = field(default_factory=list)

    # ---------------------------------------------------------- Grid settings
    rk_order: int = 3
    y_step: float = 0.01
    y_step_ternary: float = 0.05
    active_threshold: float = 1e-9
    plot_threshold: float = 1e-4
    include_0p01_to_0p05_components: bool = True

    # ------------------------------------------------------- Shared grids
    temperature_values: object = field(default_factory=lambda: np.arange(200, 2000 + 1, 200))
    mu_O_values: object = field(default_factory=lambda: np.arange(-10.0, -4.0 + 0.001, 0.10))

    # --------------------------------------------------- Oxidation scan
    scan_x: float = 0.50
    scan_T: list = field(default_factory=lambda: [700])
    scan_mu_O: object = field(default_factory=lambda: np.arange(-10.0, -4.0 + 0.001, 0.05))
    scan_oxygen_mode: str = "muO"
    scan_write_all_phases: bool = False
    scan_make_plots: bool = True
    scan_plot_only_active: bool = True

    # -------------------------------------------------- muO-x phase map
    map_x_T: int = 700
    map_x_x_values: object = field(default_factory=lambda: np.arange(0.0, 1.0 + 0.005, 0.01))
    map_x_make_plots: bool = True
    map_x_write_wide: bool = False
    anim_fps: int = 2

    # -------------------------------------------------- muO-T phase map
    map_T_x: float = 0.50
    map_T_x_values: object = None  # if set, sweeps muO-T over these x values and makes video
    map_T_make_plots: bool = True
    map_T_write_wide: bool = False

    # ----------------------------------------------- Fixed-composition map
    fixed_oxidation_threshold: float = 1e-8

    # ----------------------------------------------- Ternary 3D map
    ternary_3d_comp_step: float = 0.05

    # ---------------------------------------------- Composition slice maps
    # Integer selects by metal order; a symbol selects that element; None uses all.
    slice_axis: int | str | None = None
    # First listed element present in a system overrides slice_axis.
    slice_axis_priority: list[str] = field(default_factory=list)
    slice_remainder_ratios: list = field(default_factory=lambda: [0.0, 0.25, 1 / 3, 0.5, 2 / 3, 0.75, 1.0])
    slice_x_values: object = field(default_factory=lambda: np.arange(0.0, 1.0 + 0.005, 0.01))
    # Fixed compositions sampled by the per-composition muO-T maps only.
    slice_muT_x_values: object = field(default_factory=lambda: np.arange(0.0, 1.0 + 0.05, 0.10))
    slice_comp_bin: float = 0.05

    # ----------------------------------------------- Onset / AUC diagrams
    onset_comp_step_binary: float = 0.02
    onset_comp_step_ternary: float = 0.025
    onset_comp_step_higher: float = 0.10
    onset_max_comp_points: int = 2500
    onset_threshold: float = 1e-8

    # -------------------------------------------- Compilation / video
    make_side_by_side_miscibility: bool = True
    compile_slice_video_grids: bool = True
    max_global_video_grid_items: int = 36
    mp4_crf: int = 18
    mp4_preset: str = "slow"

    # ----------------------------------------------- Plot style
    region_label_mode: str = "id"
    region_label_fontsize: int = 12
    boundary_lw: float = 0.8

    # ----------------------------------------------- Miscibility gap TDB
    miscibility_t_start: int = 300
    miscibility_t_end: int = 3000
    # T range of the external phase diagram images (used to position the indicator line)
    phase_diagram_t_start: int = 300
    phase_diagram_t_end: int = 4500
    miscibility_t_step: int = 500
    miscibility_n_grid: int = 40
    miscibility_gap_threshold: float = 5.0
    miscibility_metal_fraction: float = field(init=False)  # auto

    def __post_init__(self):
        if self.phase_element is not None:
            self.phase_element = str(self.phase_element).strip() or None
        self.phase_element_stoichiometry = float(self.phase_element_stoichiometry)
        if self.phase_element is None:
            self.phase_element_stoichiometry = 0.0
        elif self.phase_element_stoichiometry <= 0:
            raise ValueError("phase_element_stoichiometry must be positive when phase_element is set")
        self.tables_dir = Path(self.tables_dir)
        self.figures_dir = Path(self.figures_dir)
        self.comps_dir = Path(self.comps_dir)
        self.system_structures_root = Path(self.system_structures_root)
        self.ternary_diagrams_root = Path(self.ternary_diagrams_root)
        self.blade_subdir = str(self.blade_subdir)
        self.fixed_phases_subdir = str(self.fixed_phases_subdir)
        self.elements = list(
            dict.fromkeys(str(element).strip().title() for element in self.elements if str(element).strip())
        )
        self.systems = list(dict.fromkeys(str(system).strip() for system in self.systems if str(system).strip()))
        label_mode = str(self.region_label_mode).strip().lower()
        if label_mode in {"id", "ids", "number", "numbers"}:
            self.region_label_mode = "id"
        elif label_mode in {"phase", "phases", "full"}:
            self.region_label_mode = "phases"
        else:
            raise ValueError("region_label_mode must be 'id' or 'phases'")
        self.region_label_fontsize = int(self.region_label_fontsize)
        if self.region_label_fontsize <= 0:
            raise ValueError("region_label_fontsize must be positive")
        self.miscibility_metal_fraction = 1.0 / (1.0 + self.phase_element_stoichiometry) if self.phase_element else 1.0

    # Derived convenience
    @property
    def map_x_T_values(self):
        return list(self.temperature_values)

    muO_T_temperature_values: object = None  # if None, defaults to arange(200,2001,10)

    @property
    def map_T_T_values(self):
        if self.muO_T_temperature_values is not None:
            return self.muO_T_temperature_values
        return np.arange(200, 2001, 10)

    @property
    def map_T_mu_O(self):
        return self.mu_O_values

    @property
    def map_x_mu_O(self):
        return self.mu_O_values

    @property
    def ternary_3d_mu_O(self):
        return self.mu_O_values

    @property
    def slice_mu_O(self):
        return self.mu_O_values

    @property
    def slice_T_values(self):
        return self.temperature_values

    @property
    def onset_auc_T_values(self):
        return self.temperature_values

    @property
    def onset_auc_mu_O(self):
        return self.mu_O_values
