"""Demonstrate the oxidation stability analysis framework.

The oxidation framework computes oxygen chemical-potential stability windows
across temperature and composition space for multi-component boride systems.
It produces mu_O–T diagrams, phase-fraction maps, onset-composition curves,
3-D onset surfaces, and animated phase-stability maps.

Thermodynamic equilibrium solver and Gibbs free energy minimization:
    Siya Zhu (siyazhu1@gmail.com)

Adjust files_dir before running.
"""

from __future__ import annotations

from pathlib import Path

from blade.oxidation import OxidationCalculator, OxidationConfig

# ── Configuration ─────────────────────────────────────────────────────────────
blade_root = Path(__file__).parent.parent
files_dir = blade_root / "Files"  # auto-detected in full_framework.py

oxidation_cfg = OxidationConfig(
    files_dir=files_dir,
    phase_element="B",
    phase_element_stoichiometry=2.0,
    mixed_phase_subdir="blade",  # BLADE relaxed structures
    fixed_phases_subdir="ORB",  # MP+ORB reference phases
    region_label_mode="phases",
    region_label_fontsize=7,
    slice_axis_priority=["Cr", "Hf", "Mo", "Nb", "Ta", "Ti", "V", "W", "Zr"],
)

calculator = OxidationCalculator(oxidation_cfg)


def run_batch_example() -> None:
    """Run batch oxidation analysis over a list of systems.

    Produces per-temperature CSVs, onset plots, and 3-D onset surfaces for
    each system in `systems`.

    Returns:
        None. Results are written to files_dir/oxidation/<system>/.
    """
    calculator.run_batch(
        systems=["CrHfB"],
        temperature_values=list(range(250, 2001, 250)),
        mu_O_values=[round(-10.0 + i * 0.1, 1) for i in range(61)],  # -10 to -4 eV
        run_calculations=True,
        run_plots=True,
        run_onset=True,
        run_3d_onset=True,
        run_animations=False,
        run_composition_slices=False,
        run_composition_slice_muT=False,
        run_scan=False,
        run_muO_x_map=False,
        skip_if_tables_exist=True,
        skip_if_analysis_exists=True,
        onset_comp_step_binary=0.02,
        onset_comp_step_ternary=0.025,
        onset_comp_step_nary=0.1,
        onset_threshold=1e-8,
        scan_mu_O=[round(-10.0 + i * 0.05, 2) for i in range(121)],
        scan_x=0.5,
        scan_T=[700, 1273, 2000],
        map_x_x_values=[round(i * 0.01, 2) for i in range(101)],
        slice_remainder_ratios=[0.0, 0.25, 0.5, 0.75, 1.0],
        slice_muT_comp_step=0.1,
    )
    print("Batch complete. Results in:", files_dir / "oxidation")


def run_single_example() -> None:
    """Run oxidation analysis for one specific composition.

    Returns:
        None. mu_O–T maps and scan plots written to files_dir/oxidation/.
    """
    calculator.run_single_composition(
        system="CrHfMoB",
        metals=["Cr", "Hf", "Mo"],
        composition=[1 / 3, 1 / 3, 1 / 3],
        temperatures=list(range(250, 2001, 10)),
        mu_o=[-12.0 + i * 0.1 for i in range(81)],
        scan_temperatures=[700, 1000, 1273],
        x_values=[round(i * 0.01, 2) for i in range(101)],
        run_calculations=True,
        run_plots=True,
        skip_if_exists=True,
        rk_order=3,
        y_step=0.01,
        include_0p01_to_0p05_components=True,
    )
    print("Single composition complete.")


if __name__ == "__main__":
    run_batch_example()
