"""Run batch or single-composition oxidation analysis from BLADE."""

from pathlib import Path

import numpy as np

from blade.oxidation import OxidationCalculator, OxidationConfig

# Choose "batch" (run_simple equivalent) or "single" (run_single equivalent).
MODE = "batch"

BLADE_FILES = Path("/Users/chasekatz/Desktop/School/Research/BLADE/Files")

# Shared chemistry and plotting settings. For a metal-only phase, use None/0.0.
config = OxidationConfig(
    files_dir=BLADE_FILES,
    phase_element="B",
    phase_element_stoichiometry=2.0,
    region_label_mode="phases",  # use "id" for numbered regions
    region_label_fontsize=7,
    slice_axis_priority=["Cr", "Hf", "Mo", "Nb", "Ta", "Ti", "V", "W", "Zr"],
)


# ---------------------------------------------------------------------------
# Batch mode: equivalent to the generalized batch oxidation workflow
# ---------------------------------------------------------------------------
BATCH_TEMPERATURES = np.arange(200, 2001, 200)
BATCH_MU_O = np.arange(-10.0, -3.999, 0.1)

BATCH_SETTINGS = {
    # Empty includes every system. Otherwise, systems may only contain these metals.
    "elements": ["Cr", "Hf", "Mo"],
    # Both true performs a complete calculation pass, then a cache-only plot pass.
    "run_calculations": True,
    "run_plots": True,
    "run_scan": True,
    "run_muO_x_map": False,
    "run_onset": True,
    "run_3d_onset": True,
    "run_animations": True,
    "run_composition_slices": True,
    "run_composition_slice_muT": True,
    "skip_if_tables_exist": True,
    "skip_if_analysis_exists": True,
    "temperature_values": BATCH_TEMPERATURES,
    "mu_O_values": BATCH_MU_O,
    "slice_remainder_ratios": [0.0, 0.25, 1 / 3, 0.5, 2 / 3, 0.75, 1.0],
}


# ---------------------------------------------------------------------------
# Single mode: equivalent to the fixed-composition oxidation workflow
# ---------------------------------------------------------------------------
SINGLE_SYSTEM = "CrHfMoB"
SINGLE_METALS = ["Cr", "Hf", "Mo"]
SINGLE_COMPOSITION = [1 / 3, 1 / 3, 1 / 3]
SINGLE_TEMPERATURES = np.arange(200, 2001, 10)
SINGLE_MU_O = np.arange(-12.0, -3.999, 0.1)
SINGLE_SCAN_T = [700, 1000, 1273, 1300]
SINGLE_X_VALUES = np.arange(0.0, 1.001, 0.01)

SINGLE_ANALYZER_SETTINGS = {
    "include_0p01_to_0p05_components": True,
    "rk_order": 3,
    "y_step": 0.01,
}


def run_batch() -> list:
    calculator = OxidationCalculator(config)
    return calculator.run_batch(**BATCH_SETTINGS)


def run_single() -> Path:
    calculator = OxidationCalculator(config)
    analyzer = calculator.single(
        system=SINGLE_SYSTEM,
        metals=SINGLE_METALS,
        composition=SINGLE_COMPOSITION,
        **SINGLE_ANALYZER_SETTINGS,
    )
    return analyzer.run(
        T_values=SINGLE_TEMPERATURES,
        mu_O_values=SINGLE_MU_O,
        scan_T=SINGLE_SCAN_T,
        x_values=SINGLE_X_VALUES,
        skip_if_exists=True,
        run_calculations=True,
        run_plots=True,
    )


if __name__ == "__main__":
    selected_mode = MODE.strip().lower()
    if selected_mode == "batch":
        run_batch()
    elif selected_mode == "single":
        run_single()
    else:
        raise ValueError("MODE must be 'batch' or 'single'")
