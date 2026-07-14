"""Full BLADE workflow: compositions → SQS → TDB → phase diagrams → visualization.

This is the main driver script showing an example of HEDBs but uses the refactored BLADE API:

  - BladeTDBGen stores config on construction; call gen.fit() to run
  - BladeVisualizer (not BLADEVisualizer)
  - All classes imported from blade.*
"""

from collections import Counter

from pathlib import Path

from pycalphad import Database

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_sqsgen import BladeSQS
from blade.tools.blade_tdb_gen import BladeTDBGen
from blade.analysis.blade_visual import BladeVisualizer

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# ------------------------------------------------------------------
# Run flags
# ------------------------------------------------------------------
level = 10
sqs_iter = 1_000_000
run_sqs = True
skip_existing_sqs = False

run_tdb = True
skip_existing_tdb = False

# ------------------------------------------------------------------
# MLIP calculator — change mlip to swap potentials (see MaterialsFramework registry)
# ------------------------------------------------------------------
mlip        = "orb"                            # e.g. "grace", "mace", "uma", "chgnet", "orb"
mlip_kwargs = {"steps": 1000, "device": "cpu"} # kwargs forwarded to the calculator constructor

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
    "track_trajectory": True,   # set False to skip relaxation_live.xyz
    "terms": None,
}

# Optional per-phase model overrides — set to None to use ATAT defaults.
# Keys are lattice base names (e.g. "CARBIDE1", "FCC1", "BCC1").

terms_in: dict | None = None
terms_in = {"HEDB4": 
    "1,0:1,0\n"
    "2,2:1,0\n"
    "3,0:3,0\n"
    "4,0:4,0\n"
    }

mult_in: dict | None = None
# mult_in = {"HEDB1": "a=1\tb=2"}

sublattice_map: dict | None = None
# sublattice_map = {"HEDB1": {"a": ["Cr", "Hf"]}}

sqsgen_in: dict | None = None
# sqsgen_in = {"HEDB1": "level=0\ta=1\tb=1
# level=1\ta=0.5,0.5\tb=1"}

fixed_compositions: dict | None = None
# fixed_compositions: dict | None = {
#     "a": [0.5, 0.5]}

run_movie = False

# ------------------------------------------------------------------
# Elements and composition constraints
# ------------------------------------------------------------------
primary_elements = ["Zr", "Hf", "Ta", "Cr"]
secondary_elements: list[str] = []

primary_min = 4
primary_max = 4
secondary_min = 0
secondary_max = 0

# ------------------------------------------------------------------
# AlB2-type diboride lattice constants (Å, from DFT literature)
# Used to auto-compute a, b, c from the chosen primary elements.
# ------------------------------------------------------------------
_DIBORIDE_A = {
    "Cr": 2.969, "Hf": 3.141, "Mo": 3.053, "Nb": 3.086,
    "Ta": 3.098, "Ti": 3.028, "V":  2.998, "W":  3.020, "Zr": 3.169,
}
_DIBORIDE_C = {
    "Cr": 3.066, "Hf": 3.470, "Mo": 3.169, "Nb": 3.269,
    "Ta": 3.227, "Ti": 3.228, "V":  3.057, "W":  3.137, "Zr": 3.530,
}
_active_d = [el for el in primary_elements if el in _DIBORIDE_A]
_avg_a = sum(_DIBORIDE_A[el] for el in _active_d) / len(_active_d)
_avg_c = sum(_DIBORIDE_C[el] for el in _active_d) / len(_active_d)
print(f"Diboride lattice estimate: a={_avg_a:.4f} Å  c={_avg_c:.4f} Å  (avg of {_active_d})")

# ------------------------------------------------------------------
# Phase prototypes
# ------------------------------------------------------------------
phases: dict[str, dict] = {
    "HEDB4": {
        "a": _avg_a,
        "b": _avg_a,
        "c": _avg_c,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.333333 0.666667 0.500000 B\n"
            "0.666667 0.333333 0.500000 B\n"
        ),
    },
}

phase_list = [
    # supercell_size controls n_atoms = unit_cell_sites × product(supercell_size)
    # (2,2,2) → 64 atoms   (4,4,2) → 128   (4,4,4) → 256
    {"generator_name": "HEDB", "lattice": "HEDB4",  "supercell_size": (4, 4, 3)},
]

liquid = False

# ------------------------------------------------------------------
# SQS composition levels
# ------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]],                            "letter": ["a"]},
    {"level": 1, "compositions": [[0.5, 0.5]],                             "letter": ["a"]},
    {"level": 2, "compositions": [[0.75, 0.25]],                           "letter": ["a"]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]],            "letter": ["a"]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]],                      "letter": ["a"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]],         "letter": ["a"]},
    {"level": 6, "compositions": [[0.75, 0.125, 0.125]],                   "letter": ["a"]},
    {"level": 7, "compositions": [[0.25, 0.25, 0.25, 0.25]],               "letter": ["a"]},
    {"level": 8, "compositions": [[0.5, 0.16667, 0.16667, 0.16667]],       "letter": ["a"]},
    {"level": 9, "compositions": [[0.625, 0.125, 0.125, 0.125]],           "letter": ["a"]},
    {"level": 10, "compositions": [[0.4, 0.2, 0.2, 0.2]],                  "letter": ["a"]},
]


# ------------------------------------------------------------------
# mcsqs run parameters
# ------------------------------------------------------------------
mcsqs_params = {
    "use_time": True,
    "time": 60,          # seconds per sqsdb directory
    "cutoff_mode": "nn", # "nn" = NN shell index (decimals OK), "angstrom" = direct Å
    "2": 5,              # pair cutoff  = 5th-nearest shell
    "3": 4,              # triplet cutoff
    "4": 3,              # quadruplet cutoff
    "wr": 20,
    "wn": 0.75,
    "wd": 1,
    "parallel_runs": 10,
}

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
# 2. Generate SQS structures
# ------------------------------------------------------------------
if run_sqs:
    for specific_phase in phase_list:
        for len_comp in unique_len_comps:
            lattice = specific_phase["lattice"]
            sqs_gen = BladeSQS(
                phases_dict=phases[lattice],
                sqsgen_levels=sqsgen_levels,
                level=level,
                len_comp=len_comp,
                skip_existing_sqs=skip_existing_sqs,
                sqsgen_in=sqsgen_in.get(lattice) if sqsgen_in else None,
                fixed_compositions=fixed_compositions,
            )
            params = mcsqs_params | {"super_cell_size": specific_phase["supercell_size"]}
            sqs_gen.sqs_gen(phase=specific_phase, paths=paths, iter=sqs_iter, params=params)

# ------------------------------------------------------------------
# 3. Fit TDB databases
# ------------------------------------------------------------------
comps_dir = path1 / "Files" / "Comps"

if run_tdb:
    gen = BladeTDBGen(
        phases=phase_list,
        phases_dict=phases,
        liquid=liquid,
        paths=paths,
        composition_list=composition_list,
        level=level,
        skip_existing=skip_existing_tdb,
        output_dir=comps_dir,
        terms_in=terms_in,
        mult_in=mult_in,
        sublattice_map=sublattice_map,
        tdb_params=tdb_params,
        fixed_compositions=fixed_compositions,
    )
    gen.fit()  # explicit call — no side effects on construction

# ------------------------------------------------------------------
# 4. Plot Gibbs energy and phase diagrams
# ------------------------------------------------------------------
# Derive fixed-sublattice elements from phase coords.
# Lowercase single-letter labels are variable sublattice sites;
# uppercase labels are fixed element symbols (e.g. "B" in diborides).
_coords = phases[phase_list[0]["lattice"]]["coords"]
_labels = [ln.strip().split()[-1] for ln in _coords.strip().splitlines() if ln.strip()]
_fixed = [l for l in _labels if not (len(l) == 1 and l.islower())]
remove_elements = set(_fixed)
fixed_species = {el: count / len(_labels) for el, count in Counter(_fixed).items()}

filt_comp_list = [
    [el for el in comp if el not in remove_elements]
    for comp in composition_list
]

viz = BladeVisualizer()

for comp, comp_filt in zip(composition_list, filt_comp_list):
    comp_name = "".join(comp_filt)
    comp_dir = comps_dir / comp_name
    if not comp_dir.exists():
        continue
    phase_name = f"{phase_list[0]['generator_name']}1_{len(comp_filt)}"
    tdb_phases = [phase_name]
    for tdb_file in comp_dir.glob("*.tdb"):
        tdb = Database(str(tdb_file))
        viz.plot_gibbs_energy(
            tdb=tdb,
            metals=comp_filt,
            phase=phase_name,
            fixed_species=fixed_species,
            temperatures=[300, 1000, 2000, 3000, 4000],
            output_path=comp_dir / f"{comp_name}_Gibbs_Energy.png",
        )
        viz.plot_gibbs_mixing(
            tdb=tdb,
            metals=comp_filt,
            phase=phase_name,
            fixed_species=fixed_species,
            temperatures=[300, 1000, 2000, 3000, 4000],
            output_path=comp_dir / f"{comp_name}_Gibbs_Mixing.png",
        )
        viz.plot_binary_phase_diagram(
            tdb=tdb,
            metals=comp_filt,
            phases=tdb_phases,
            fixed_species=fixed_species,
            temperature_range=(300, 4500, 50),
            output_path=comp_dir / f"{comp_name}_Phase_Diagram.png",
        )

# ------------------------------------------------------------------
# 5. Visualize
# ------------------------------------------------------------------

# Combine phase diagram PNGs
pngs = []
for comp_filt in filt_comp_list:
    comp_dir = comps_dir / "".join(comp_filt)
    pngs.extend(comp_dir.glob("*_Phase_Diagram.png"))

if pngs:
    out_pd = comps_dir / "Combined_Phase_Diagrams.png"
    viz.phase_diagram(pngs, save=out_pd)
    print(f"Combined phase diagrams → {out_pd}")

# Combine CONTCAR structures for each phase in every composition
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
        viz.contcar(contcars, save=out)
        print(f"Saved combined CONTCARs → {out}")

# Relaxation movies (requires trajectory files from TDB fitting)
if run_movie:
    viz.make_combined_relaxation_movie(
        composition_list=filt_comp_list,
        path1=comps_dir,
        traj_name="relaxation_live.xyz",
        fps=10,
    )
