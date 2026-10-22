"""BLADE workflow for MAX system using m_scraps instead of mcsqs for SQS generation.

Mirrors tdb_gen_max_scraps.py exactly, replacing ScrapsSQSGen with MScrapsSQSGen.
m_scraps is a Python package (no tools-directory import needed) and adds support
for SRO targets, higher-order figure correlations, and MC kernel selection.

Output goes to Files/Comps_mscraps/ so mcsqs, SCRAPS, and m_scraps results can
be compared side-by-side without overwriting each other.
"""

from collections import Counter
from pathlib import Path

from pycalphad import Database

from blade.analysis.blade_visual import BladeVisualizer
from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_mscraps_gen import MScrapsSQSGen
from blade.tools.blade_tdb_gen import BladeTDBGen

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# m_scraps repo — binary lives at MSCRAPS_REPO/SCRAPs/scraps
MSCRAPS_REPO = path0 / "m_scraps"
MSCRAPS_BIN = MSCRAPS_REPO / "SCRAPs" / "scraps"

# ------------------------------------------------------------------
# Run flags
# ------------------------------------------------------------------
level = 2
run_sqs = True
skip_existing_sqs = False

run_tdb = True
skip_existing_tdb = False
refit_existing_tdb = False
skip_existing_plots = False
generate_gibbs_energy = True
generate_gibbs_mixing = True
generate_phase_diagram = True
generate_combined_phase_diagram = True
generate_contcar_plots = True

# ------------------------------------------------------------------
# MLIP calculator
# ------------------------------------------------------------------
mlip = "orb"
mlip_kwargs = {"steps": 1000, "device": "cpu"}

tdb_params = {
    "fmax": 1e-4,
    "verbose": True,
    "calculator": mlip,
    "calculator_kwargs": mlip_kwargs,
    "t_min": 298.15,
    "t_max": 10000.0,
    "sro": False,
    "bv": 1e-3,
    "phonon": False,
    "open_calphad": False,
    "track_trajectory": True,
    "terms": None,
}

# ------------------------------------------------------------------
# Optional per-phase model overrides
# ------------------------------------------------------------------
terms_in: dict | None = None
terms_in = {"MAX2": "1,0:1,0:1,0\n" "2,2:2,2:1,0\n" "3,0:2,0:1,0\n"}

mult_in: dict | None = None
mult_in = {"MAX2": "a=2\tb=1\tc=1"}

sublattice_map: dict | None = None
sublattice_map = {
    "MAX2": {
        "a": ["Zr", "Hf", "Ti", "Nb", "V"],
        "b": ["Al", "Sn"],
        "Constant": ["a", "b"],
    }
}

sqsgen_in: dict | None = None
fixed_compositions: dict | None = {
    "a": [0.2, 0.2, 0.2, 0.25, 0.15],
    "b": [0.40, 0.60],
}

run_movie = False

# ------------------------------------------------------------------
# Elements and composition constraints
# ------------------------------------------------------------------
primary_elements = ["Zr", "Hf", "Ti", "Nb", "V"]
secondary_elements: list[str] = []

primary_min = 5
primary_max = 5
secondary_min = 0
secondary_max = 0

# ------------------------------------------------------------------
# 211 MAX-phase (M2AX) lattice constants (Å, from DFT literature)
# ------------------------------------------------------------------
_max_a = {
    "Cr": 2.86, "Hf": 3.28, "Mo": 2.96, "Nb": 3.10,
    "Ta": 3.08, "Ti": 3.04, "V": 2.91, "W": 2.91, "Zr": 3.32,
}
_max_c = {
    "Cr": 12.80, "Hf": 14.36, "Mo": 13.20, "Nb": 13.80,
    "Ta": 14.10, "Ti": 13.60, "V": 13.20, "W": 13.80, "Zr": 14.57,
}
_active_m = [el for el in primary_elements if el in _max_a]
_avg_a = sum(_max_a[el] for el in _active_m) / len(_active_m)
_avg_c = sum(_max_c[el] for el in _active_m) / len(_active_m)
print(f"MAX lattice estimate: a={_avg_a:.4f} Å  c={_avg_c:.4f} Å  (avg of {_active_m})")

# ------------------------------------------------------------------
# Phase prototypes
# ------------------------------------------------------------------
phases: dict[str, dict] = {
    "MAX2": {
        "a": _avg_a,
        "b": _avg_a,
        "c": _avg_c,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.333333 0.666667 0.083333 a\n"
            "0.666667 0.333333 0.916667 a\n"
            "0.666667 0.333333 0.583333 a\n"
            "0.333333 0.666667 0.416667 a\n"
            "0.333333 0.666667 0.750000 b\n"
            "0.666667 0.333333 0.250000 b\n"
            "0.000000 0.000000 0.000000 C\n"
            "0.000000 0.000000 0.500000 C\n"
        ),
    },
}

phase_list = [
    {"generator_name": "MAX", "lattice": "MAX2", "supercell_size": (5, 5, 2)},
]

liquid = False

# ------------------------------------------------------------------
# SQS composition levels
# ------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]], "letter": ["a", "b"]},
    {"level": 1, "compositions": [[0.5, 0.5]], "letter": ["a", "b"]},
    {"level": 2, "compositions": [[0.75, 0.25]], "letter": ["a", "b"]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]], "letter": ["a", "b"]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]], "letter": ["a", "b"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]], "letter": ["a", "b"]},
    {"level": 6, "compositions": [[0.75, 0.125, 0.125]], "letter": ["a", "b"]},
]

# ------------------------------------------------------------------
# m_scraps parameters
# ------------------------------------------------------------------
mscraps_ranks = 10
auto_budget = 3        # 1=basic, 2=thorough, 3=exhaustive
max_shellnum = 4
fix_multibasis_sublattice = True
mc_kernel = "sa"       # "sa" (simulated annealing) or "pt" (parallel tempering)

# SRO targets: list of (species_i, species_j, shell, alpha) — None disables.
# Example: sro_targets = [(0, 1, 1, -0.5)]  # 50/50 species pair, 1st shell, SRO=-0.5
sro_targets: list[tuple[int, int, int, float]] | None = None

# Figure cutoffs: list of cutoff distances in Å for higher-order correlations.
figure_cutoffs: list[float] | None = None

# Figure species: species indices included in figure targets; None = all.
figure_species: list[int] | None = None

figure_objective = "minimize"  # "minimize" or "match"
figure_wd = 1.0
figure_wn = 0.75
figure_wr = 1.0
figure_tol = 1e-4

# ------------------------------------------------------------------
# 1. Generate compositions
# ------------------------------------------------------------------
composer = BladeCompositions(
    primary_elements=primary_elements,
    secondary_elements=secondary_elements,
    primary_min=primary_min,
    primary_max=primary_max,
    secondary_min=secondary_min,
    secondary_max=secondary_max,
)

composition_list = composer.generate_compositions()
unique_len_comps = composer.get_systems()

print(f"Compositions ({len(composition_list)} total): {composition_list}")
print(f"System sizes: {unique_len_comps}")

# ------------------------------------------------------------------
# 2. Generate SQS structures via m_scraps
# ------------------------------------------------------------------
if run_sqs:
    for specific_phase in phase_list:
        for len_comp in unique_len_comps:
            lattice = specific_phase["lattice"]
            sqs_gen = MScrapsSQSGen(
                phases_dict=phases[lattice],
                sqsgen_levels=sqsgen_levels,
                level=level,
                len_comp=len_comp,
                skip_existing_sqs=skip_existing_sqs,
                sublattice_map=sublattice_map,
                sqsgen_in=sqsgen_in.get(lattice) if sqsgen_in else None,
                fixed_compositions=fixed_compositions,
                scraps_bin=MSCRAPS_BIN,
                ranks=mscraps_ranks,
                auto_budget=auto_budget,
                max_shellnum=max_shellnum,
                fix_multibasis_sublattice=fix_multibasis_sublattice,
                sro_targets=sro_targets,
                figure_cutoffs=figure_cutoffs,
                figure_species=figure_species,
                figure_objective=figure_objective,
                figure_wd=figure_wd,
                figure_wn=figure_wn,
                figure_wr=figure_wr,
                figure_tol=figure_tol,
                mc_kernel=mc_kernel,
            )
            params = {"super_cell_size": specific_phase["supercell_size"]}
            sqs_gen.sqs_gen(phase=specific_phase, paths=paths, params=params)

# ------------------------------------------------------------------
# 3. Fit TDB databases
# ------------------------------------------------------------------
comps_dir = path1 / "Files" / "Comps_mscraps"

if run_tdb:
    gen = BladeTDBGen(
        phases=phase_list,
        phases_dict=phases,
        liquid=liquid,
        paths=paths,
        composition_list=composition_list,
        level=level,
        skip_existing=skip_existing_tdb,
        refit_existing=refit_existing_tdb,
        output_dir=comps_dir,
        terms_in=terms_in,
        mult_in=mult_in,
        sublattice_map=sublattice_map,
        tdb_params=tdb_params,
        fixed_compositions=fixed_compositions,
    )
    gen.fit()

# ------------------------------------------------------------------
# 4. Plot Gibbs energy and phase diagrams
# ------------------------------------------------------------------
_coords = phases[phase_list[0]["lattice"]]["coords"]
_labels = [ln.strip().split()[-1] for ln in _coords.strip().splitlines() if ln.strip()]
_fixed = [label for label in _labels if not (len(label) == 1 and label.islower())]
remove_elements = set(_fixed)
fixed_species = {el: count / len(_labels) for el, count in Counter(_fixed).items()}

filt_comp_list = [[el for el in comp if el not in remove_elements] for comp in composition_list]

viz = BladeVisualizer()

for _comp, comp_filt in zip(composition_list, filt_comp_list, strict=False):
    if len(comp_filt) != 2:
        continue
    comp_name = "".join(comp_filt)
    comp_dir = comps_dir / comp_name
    if not comp_dir.exists():
        continue
    phase_name = f"{phase_list[0]['lattice']}_{len(comp_filt)}"
    tdb_phases = [phase_name]
    plot_paths = (
        comp_dir / f"{comp_name}_Gibbs_Energy.png",
        comp_dir / f"{comp_name}_Gibbs_Mixing.png",
        comp_dir / f"{comp_name}_Phase_Diagram.png",
    )
    make_energy = generate_gibbs_energy and not (skip_existing_plots and plot_paths[0].exists())
    make_mixing = generate_gibbs_mixing and not (skip_existing_plots and plot_paths[1].exists())
    make_phase = generate_phase_diagram and not (skip_existing_plots and plot_paths[2].exists())
    if not any((make_energy, make_mixing, make_phase)):
        continue
    for tdb_file in comp_dir.glob("*.tdb"):
        tdb = Database(str(tdb_file))
        if make_energy:
            viz.plot_gibbs_energy(
                tdb=tdb,
                metals=comp_filt,
                phase=phase_name,
                fixed_species=fixed_species,
                temperatures=[300, 1000, 2000, 3000, 4000],
                output_path=plot_paths[0],
            )
        if make_mixing:
            viz.plot_gibbs_mixing(
                tdb=tdb,
                metals=comp_filt,
                phase=phase_name,
                fixed_species=fixed_species,
                temperatures=[300, 1000, 2000, 3000, 4000],
                output_path=plot_paths[1],
            )
        if make_phase:
            viz.plot_binary_phase_diagram(
                tdb=tdb,
                metals=comp_filt,
                phases=tdb_phases,
                fixed_species=fixed_species,
                temperature_range=(300, 4500, 50),
                output_path=plot_paths[2],
            )

# ------------------------------------------------------------------
# 5. Visualize
# ------------------------------------------------------------------
if generate_combined_phase_diagram:
    pngs = []
    for comp_filt in filt_comp_list:
        comp_dir = comps_dir / "".join(comp_filt)
        pngs.extend(comp_dir.glob("*_Phase_Diagram.png"))
    out_pd = comps_dir / "Combined_Phase_Diagrams.png"
    if pngs and not (skip_existing_plots and out_pd.exists()):
        viz.phase_diagram(pngs, save=out_pd)
        print(f"Combined phase diagrams → {out_pd}")

if generate_contcar_plots:
    for comp_filt in filt_comp_list:
        comp_name = "".join(comp_filt)
        comp_dir = comps_dir / comp_name
        if not comp_dir.exists():
            continue
        for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
            contcars = sorted(phase_dir.glob("sqs_lev=*/CONTCAR"))
            if not contcars:
                continue
            out = comp_dir / f"Combined_CONTCARs_{comp_name}_{phase_dir.name}.png"
            if skip_existing_plots and out.exists():
                continue
            viz.contcar(contcars, save=out)
            print(f"Saved combined CONTCARs → {out}")

if run_movie:
    viz.make_combined_relaxation_movie(
        composition_list=filt_comp_list,
        path1=comps_dir,
        traj_name="relaxation_live.xyz",
        fps=10,
    )
