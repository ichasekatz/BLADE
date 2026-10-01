"""Entry point for running the full BLADE pipeline.

Reads full_framework.toml and executes all pipeline stages in sequence:
BladeCompositions -> BladeSQS -> BladeTDBGen -> phase diagrams -> oxidation.

Each stage is gated by the [stages] table in the TOML. The [elements] table
is the single source of truth for the element pool shared across all stages.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

from blade.analysis.blade_visual import BladeVisualizer
from blade.oxidation import OxidationCalculator, OxidationConfig
from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_sqsgen import BladeSQS
from blade.tools.blade_tdb_gen import BladeTDBGen

if __name__ == "__main__":
    # blade_root auto-detected from this file; only sqsdb_dir must be set explicitly
    blade_root = Path(__file__).parent.parent
    with (Path(__file__).parent / "full_framework.toml").open("rb") as _f:
        cfg = tomllib.load(_f)

    sqsdb_dir = Path(cfg["paths"]["sqsdb_dir"])
    files_dir = blade_root / "Files"
    comps_dir = files_dir / cfg["tdb"].get("comps_folder", "Comps")
    paths = [blade_root.parent, blade_root, sqsdb_dir]
    stages = cfg.get("stages", {})
    tdb_cfg = cfg.get("tdb", {})
    level = tdb_cfg.get("level", 6)

    # [elements] is the single source of truth for the element pool
    _el = cfg["elements"]
    primary_elements = _el["primary"]
    secondary_elements = _el.get("secondary", [])
    # Phase prototype from [phase]
    _ph = cfg["phase"]
    phases_dict = {_ph["key"]: {k: _ph[k] for k in ("a", "b", "c", "alpha", "beta", "gamma", "vectors", "coords")}}
    phase_list = [
        {"generator_name": _ph["generator_name"], "lattice": _ph["key"], "supercell_size": tuple(_ph["supercell_size"])}
    ]

    # --- 1. Compositions -------------------------------------------------
    composer = BladeCompositions(
        primary_elements=primary_elements,
        secondary_elements=secondary_elements,
        primary_min=_el.get("primary_min", 2),
        primary_max=_el.get("primary_max", 3),
        secondary_min=_el.get("secondary_min", 0),
        secondary_max=_el.get("secondary_max", 0),
    )
    composition_list = composer.generate_compositions()
    unique_len_comps = composer.get_systems()
    print(f"Compositions ({len(composition_list)}): {composition_list}")

    # --- 2. SQS generation -----------------------------------------------
    if stages.get("sqs_generation", False):
        sqsgen_levels = cfg.get("tdb_sqs_levels", [])
        mcsqs_params = cfg.get("tdb_sqs", {})
        for sp in phase_list:
            for len_comp in unique_len_comps:
                sqs = BladeSQS(
                    phases_dict=phases_dict[sp["lattice"]],
                    sqsgen_levels=sqsgen_levels,
                    level=level,
                    len_comp=len_comp,
                    skip_existing_sqs=tdb_cfg.get("skip_existing_sqs", True),
                )
                sqs.sqs_gen(phase=sp, paths=paths, params=mcsqs_params | {"super_cell_size": sp["supercell_size"]})

    # --- 3. TDB fitting --------------------------------------------------
    if stages.get("tdb_fitting", False):
        _fit = cfg.get("tdb_fit", {})
        BladeTDBGen(
            phases=phase_list,
            phases_dict=phases_dict,
            liquid=tdb_cfg.get("liquid", False),
            paths=paths,
            composition_list=composition_list,
            level=level,
            skip_existing=tdb_cfg.get("skip_existing_tdb", True),
            refit_existing=tdb_cfg.get("refit_existing_tdb", False),
            output_dir=comps_dir,
            terms_in=cfg.get("tdb_inputs", {}).get("terms_in") or None,
            tdb_params={
                "fmax": _fit.get("fmax", 1e-3),
                "verbose": _fit.get("verbose", True),
                "calculator": tdb_cfg.get("mlip", "orb"),
                "calculator_kwargs": cfg.get("tdb_mlip_kwargs", {}),
                "t_min": _fit.get("t_min", 298.15),
                "t_max": _fit.get("t_max", 10000.0),
                "sro": _fit.get("sro", False),
                "bv": _fit.get("bv", 1e-3),
            },
        ).fit()

    # --- 4. Phase diagrams — generate_plots gates all five plot types ----
    if generate_plots := stages.get("phase_diagrams", False):
        from pycalphad import Database

        viz = BladeVisualizer()
        for comp in (c for c in composition_list if len(c) == 2):
            comp_dir = comps_dir / "".join(comp)
            phase_name = f"{_ph['generator_name']}1_{len(comp)}"
            if not comp_dir.exists():
                continue
            for tdb_file in comp_dir.glob("*.tdb"):
                tdb = Database(str(tdb_file))
                viz.plot_gibbs_energy(
                    tdb=tdb, metals=list(comp), phase=phase_name, output_path=comp_dir / f"{''.join(comp)}_Gibbs_Energy.png"
                )
                viz.plot_gibbs_mixing(
                    tdb=tdb, metals=list(comp), phase=phase_name, output_path=comp_dir / f"{''.join(comp)}_Gibbs_Mixing.png"
                )
                viz.plot_binary_phase_diagram(
                    tdb=tdb, metals=list(comp), phases=[phase_name], output_path=comp_dir / f"{''.join(comp)}_Phase_Diagram.png"
                )
        pngs = list(comps_dir.rglob("*_Phase_Diagram.png"))
        if pngs:
            viz.phase_diagram(pngs, save=comps_dir / "Combined_Phase_Diagrams.png")
        for comp in composition_list:
            comp_dir = comps_dir / "".join(comp)
            if not comp_dir.exists():
                continue
            for pd in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
                contcars = sorted(pd.glob("sqs_lev=*/CONTCAR"))
                if contcars:
                    viz.contcar(contcars, save=comp_dir / f"Combined_CONTCARs_{''.join(comp)}_{pd.name}.png")

    # --- 5. Oxidation analysis -------------------------------------------
    if stages.get("oxidation_graphs", False):
        _ox = cfg.get("oxidation", {})
        _ox_b = cfg.get("oxidation_batch", {})
        _T = _ox_b.get("temperature_range", [250, 2000, 250])
        _mu = _ox_b.get("mu_o_range", [-10.0, -4.0, 0.1])
        OxidationCalculator(
            OxidationConfig(
                files_dir=files_dir,
                phase_element=_ox.get("phase_element", "B"),
                phase_element_stoichiometry=_ox.get("phase_element_stoichiometry", 2.0),
                mixed_phase_subdir=_ox.get("mixed_phase_subdir", "blade"),
                fixed_phases_subdir=_ox.get("fixed_phases_subdir", "ORB"),
                region_label_mode=_ox.get("region_label_mode", "phases"),
                slice_axis_priority=_ox.get("slice_axis_priority", primary_elements),
            )
        ).run_batch(
            systems=_ox_b.get("systems", []),
            temperature_values=list(range(_T[0], _T[1] + 1, _T[2])),
            mu_O_values=[round(_mu[0] + i * _mu[2], 3) for i in range(int((_mu[1] - _mu[0]) / _mu[2]) + 1)],
            run_calculations=_ox_b.get("run_calculations", True),
            run_plots=_ox_b.get("run_plots", True),
            run_onset=_ox_b.get("run_onset", True),
            run_3d_onset=_ox_b.get("run_3d_onset", True),
            skip_if_tables_exist=_ox_b.get("skip_if_tables_exist", True),
            skip_if_analysis_exists=_ox_b.get("skip_if_analysis_exists", True),
        )
