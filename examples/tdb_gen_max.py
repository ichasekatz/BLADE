"""Full BLADE workflow: compositions → SQS → TDB → phase diagrams → visualization.

This is the main driver script showing an example of MAX Phases but uses the refactored BLADE API:

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
level = 2
sqs_iter = 1_000_000
run_sqs = True
skip_existing_sqs = False

run_tdb = True
skip_existing_tdb = False

# ------------------------------------------------------------------
# MLIP calculator — change mlip to swap potentials (see MaterialsFramework registry)
# ------------------------------------------------------------------
mlip        = "grace"                            # e.g. "grace", "mace", "uma", "chgnet", "orb"
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
    "terms": None,
}

# Optional per-phase model overrides — set to None to use ATAT defaults.
# Keys are lattice base names (e.g. "CARBIDE1", "FCC1", "BCC1").

terms_in: dict | None = None
terms_in = {"MAX1": 
    "1,0:1,0:1,0\n"
    "2,0:1,0:1,0\n"
    "1,0:2,0:1,0\n"
    "2,0:2,0:1,0\n"
    }

mult_in: dict | None = None
mult_in = {"MAX1": "a=2\tb=1\tc=1"}

sublattice_map: dict | None = None
sublattice_map = {
    "MAX1": {
        "a": ["Zr", "Hf", "Ti", "Nb", "V"],
        "b": ["Al", "Sn"],
        "Constant": ["a", "b"],   # use exact order from fixed_compositions, no permutations
    }
}

sqsgen_in: dict | None = None
# sqsgen_in = {"MAX1": "level=0\ta=1\tb=1
# level=1\ta=0.5,0.5\tb=1"}

# Fix b sublattice at 0.7 Al + 0.3 Sn for every SQS composition.
# Keys match sublattice labels; values are fractional compositions in the
# same order as the element list in sublattice_map["MAX1"]["b"].
fixed_compositions: dict | None = None
fixed_compositions: dict | None = {
    "a": [0.2, 0.2, 0.2, 0.25, 0.15],
    "b": [0.40, 0.60],
}

# fixed_compositions: dict | None = {
#     "a": [0.2, 0.15, 0.15, 0.4, 0.1],       # Zr=0.2, Hf=0.15, Ti=0.15, Nb=0.4, V=0.1
#     "b": [0.30, 0.70],                      # Sn=0.3, Al=0.7
# }

run_movie = False

# ------------------------------------------------------------------
# Elements and composition constraints
# ------------------------------------------------------------------
primary_elements = ["Zr", "Hf", "Ti", "Nb", "V"]   # must match sublattice_map a-site
secondary_elements: list[str] = []

primary_min = 5   # quinary — all 5 a-site elements in one composition
primary_max = 5
secondary_min = 0
secondary_max = 0

# ------------------------------------------------------------------
# 211 MAX-phase (M2AX) lattice constants (Å, from DFT literature)
# Used to auto-compute a, b, c from the chosen primary elements.
# ------------------------------------------------------------------
_MAX_A = {
    "Cr": 2.86, "Hf": 3.28, "Mo": 2.96, "Nb": 3.10,
    "Ta": 3.08, "Ti": 3.04, "V":  2.91, "W":  2.91, "Zr": 3.32,
}

_MAX_C = {
    "Cr": 12.80, "Hf": 14.36, "Mo": 13.20, "Nb": 13.80,
    "Ta": 14.10, "Ti": 13.60, "V":  13.20, "W":  13.80, "Zr": 14.57,
}
_active_m = [el for el in primary_elements if el in _MAX_A]
_avg_a = sum(_MAX_A[el] for el in _active_m) / len(_active_m)
_avg_c = sum(_MAX_C[el] for el in _active_m) / len(_active_m)
print(f"MAX lattice estimate: a={_avg_a:.4f} Å  c={_avg_c:.4f} Å  (avg of {_active_m})")

# ------------------------------------------------------------------
# Phase prototypes
# ------------------------------------------------------------------
phases: dict[str, dict] = {
    "MAX1": {
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
    {"generator_name": "MAX", "lattice": "MAX1",  "supercell_size": (5, 5, 2)},
]

liquid = False

# ------------------------------------------------------------------
# SQS composition levels
# ------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]],                            "letter": ["a", "b"]},
    {"level": 1, "compositions": [[0.5, 0.5]],                             "letter": ["a", "b"]},
    {"level": 2, "compositions": [[0.75, 0.25]],                           "letter": ["a", "b"]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]],            "letter": ["a", "b"]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]],                      "letter": ["a", "b"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]],         "letter": ["a", "b"]},
    {"level": 6, "compositions": [[0.75, 0.125, 0.125]],                   "letter": ["a", "b"]},
]

# ------------------------------------------------------------------
# mcsqs run parameters
# ------------------------------------------------------------------
mcsqs_params = {
    "use_time": True,
    "time": 30,          # seconds per sqsdb directory
    "cutoff_mode": "nn", # "nn" = NN shell index (decimals OK), "angstrom" = direct Å
    "2": 4.5,            # pair cutoff
    "3": 0,              # triplet cutoff (0 = omit)
    "4": 0,              # quadruplet cutoff (0 = omit)
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