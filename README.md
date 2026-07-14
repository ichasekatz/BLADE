<div align="center">

# BLADE

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
[![Python](https://img.shields.io/badge/python-3.12-brightgreen.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-1.7.0-orange.svg)](pyproject.toml)

**Batch Lattice Analysis and Discovery Engine — configurable structure-to-CALPHAD automation and environmental stability analysis for multicomponent materials.**

[Report a Bug](https://github.com/ichasekatz/BLADE/issues/new?labels=bug) · [Request a Feature](https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement)

</div>

---

## Overview

BLADE automates a structure-to-thermodynamics workflow. Given element pools,
composition constraints, and one or more crystal prototypes, it enumerates
N-component systems, generates SQS supercells, relaxes them with a configurable
interatomic model, fits CALPHAD parameters, and produces `.tdb` databases. Its
environmental-stability layer then combines the modeled phase with configurable
reference compounds to map equilibrium assemblages against composition,
temperature, and chemical potential.

The reusable classes accept configurable elements, crystal prototypes,
system sizes, and variable or fixed sublattice layouts. The included HEDB and
MAX configurations are worked examples illustrating two prototype
descriptions.

BLADE builds on [MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) as a computational backend. **The [ichasekatz fork](https://github.com/ichasekatz/MaterialsFramework) is required** — install it before BLADE.

## What Is Distinctive

BLADE's contribution is workflow integration and generalization rather than a
new CALPHAD formalism or equilibrium solver:

- **Prototype-driven instead of chemistry-driven.** Lattice vectors,
  coordinates, variable sites, fixed sites, and site compositions are input
  data. The same orchestration can therefore be applied to alloys,
  intermetallics, ceramics, and other ordered or disordered prototypes.
- **System-size-aware model fitting.** Base and exact generated phase keys let
  binary, ternary, and higher-order systems use different `terms.in`,
  `mult.in`, and sublattice constraints without separate workflow code.
- **Interchangeable structure generation.** ATAT `mcsqs` and SCRAPS feed the
  same relaxation and TDB-fitting interface.
- **Incremental, artifact-level reuse.** Existing SQS results, relaxed
  structures, TDBs, refits, individual plots, equilibrium tables, and map
  caches can be reused independently. Missing outputs do not force unrelated
  calculations to repeat.
- **One configuration across scales.** A TOML input connects element-system
  enumeration, SQS generation, MLIP relaxation, TDB fitting, phase
  visualization, reference-database construction, and environmental stability
  calculations.
- **Generalized stability visualization.** Binary through N-component phase
  projections, composition slices, exact interactive hover data, onset
  surfaces, and animations are generated from the same cached equilibrium
  representation.

These are implementation-level capabilities. The physical accuracy of a run
still depends on the supplied prototype, configurations, energies, interaction
model, reference phases, and convergence settings.

---

## Modules

### Tools (`blade.tools`)

| Class | Description |
|---|---|
| `BladeCompositions` | Enumerate compositions from primary and secondary element pools |
| `BladeSQS` | Generate SQS structures using ATAT `mcsqs` in parallel |
| `ScrapsSQSGen` | Generate SQS structures using the SCRAPS algorithm |
| `BladeTDBGen` | Relax structures and fit CALPHAD thermodynamic databases |
| `BladeCutoff` | Compute neighbor-shell distances from SQS supercell geometry |

### Analysis (`blade.analysis`)

| Class | Description |
|---|---|
| `BladeVisualizer` | Plot Gibbs energy, mixing energy, and phase diagrams; combine CONTCARs and relaxation movies |
| `BLADEData` | Extract volumes, lattice parameters, and energies from POSCAR/energy files |

### Oxidation (`blade.oxidation`)

| Class | Description |
|---|---|
| `OxideCompositions` | Enumerate oxide-relevant compositions from the primary element pool |
| `OxideDatabase` | Download and relax MP reference structures; build unified energy database; run grand potential minimization |
| `OxidationCalculator` | Run generalized cached grand-potential oxidation calculations using BLADE structures |
| `OxidationVisualizer` | Regenerate oxidation maps from caches with numeric or formatted-phase region labels |

All classes are available via lazy top-level import:

```python
from blade import BladeCompositions, BladeSQS, ScrapsSQSGen, BladeTDBGen
from blade import BladeVisualizer, BLADEData
from blade.oxidation import OxideCompositions, OxideDatabase
from blade.oxidation import OxidationCalculator, OxidationConfig, OxidationVisualizer
```

The generalized oxidation backend is bundled under
`src/blade/oxidation/framework/`; no BOxidation checkout is required. Use
`examples/oxidation/oxidation_framework.py` for calculation-only batch or
single runs, or `examples/full_framework.py` for the complete TOML workflow.

---

## Getting Started

### Prerequisites

- [ichasekatz/MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) — required fork, install before BLADE
- [ATAT](https://axelvandewalle.github.io/www-avdw/atat/) — `mcsqs`, `sqs2tdb`, `corrdump` must be on `$PATH`

### Installation

#### uv

Create the locked Python 3.12 environment and install the complete workflow:

```bash
uv sync --extra workflow
```

For the importable BLADE library without the external CALPHAD, Materials
Project, and default MLIP integrations:

```bash
uv sync
```

Run commands inside the environment without activating it:

```bash
uv run --extra workflow python examples/full_framework.py --check
uv run --extra workflow python examples/full_framework.py --dry-run
uv run ruff check src examples
uv run black --check src examples
uv run isort --check-only src examples
```

The `workflow` extra installs the required `ichasekatz/MaterialsFramework`
fork directly from GitHub together with pycalphad, `mp-api`, and the default
ORB integration. ATAT and SCRAPS remain external executables and must still be
available at the configured paths.

#### Pixi

Install the locked environment with [pixi](https://pixi.sh):

```bash
pixi install
```

When developing MaterialsFramework alongside BLADE, a local editable checkout
can replace the locked Git source:

```bash
git clone https://github.com/ichasekatz/MaterialsFramework.git
uv pip install -e ../MaterialsFramework
# or, in the Pixi environment:
pixi run pip install -e ../MaterialsFramework
```

An editable BLADE install is already created by `uv sync`. For a manual Pixi
editable install:

```bash
pixi run pip install -e .
```

---

## Complete TOML Workflow

[`examples/full_framework.py`](examples/full_framework.py) runs the complete
BLADE workflow from one input file:

1. Generate compositions and SQS structures.
2. Relax structures and fit TDB databases.
3. Plot phase diagrams and optional animation grids.
4. Build the Materials Project/MLIP oxidation database.
5. Run batch or fixed-composition grand-potential oxidation analysis.

The first section of each TOML enables or skips the major workflow parts:

```toml
[stages]
sqs_generation = true
tdb_fitting = true
phase_diagrams = true
oxide_database = true
oxidation_graphs = true
```

These switches override detailed settings in later sections. For example,
`tdb_fitting = false` prevents fitting and its Gibbs/phase/CONTCAR outputs,
while `phase_diagrams = false` skips the separate phase-map, GIF, and grid
stage. Existing SQS or TDB outputs can still be reused by enabling a later
stage while disabling the earlier one.

```bash
export MP_API_KEY="your-materials-project-key"

python examples/full_framework.py --check
python examples/full_framework.py --dry-run
python examples/full_framework.py
```

Pass a different TOML file as the first argument when maintaining several
phase models:

```bash
python examples/full_framework.py examples/full_framework_max.toml --check
python examples/full_framework.py examples/full_framework_max.toml
```

Two complete example inputs are included. They exercise different sublattice
layouts but do not restrict which prototypes can be supplied:

| Input | Purpose |
|---|---|
| [`full_framework.toml`](examples/full_framework.toml) | One-variable-sublattice demonstration with environmental analysis enabled |
| [`full_framework_max.toml`](examples/full_framework_max.toml) | Two-variable-sublattice demonstration with fixed sites |

The second input disables downstream environmental analysis by default because
that example has two independently variable sublattices, while the current
environmental model represents one mixed phase-composition vector. Other
prototypes with that same reduction can use the module directly.

### SQS method

Both TOMLs default to ATAT `mcsqs`:

```toml
[tdb]
sqs_method = "mcsqs" # "mcsqs" or "scraps"
prototype_driver = "standard" # "standard" or "multibasis"
```

`prototype_driver` selects the implementation shape, not a chemistry:
`standard` handles the usual prototype input, while `multibasis` enables the
adapter for independently configured basis/site groups. Phase names and
elements remain arbitrary in both modes.

Setting `sqs_method = "scraps"` selects the SCRAPS structure-generation
backend. Configure its repository under `[paths]` and tuning parameters under
`[tdb_scraps]`:

```toml
[paths]
scraps_repo = "/path/to/SCRAPS/scraps-perpair"

[tdb_scraps]
ranks = 10
auto_budget = 3
max_shellnum = 4
fix_multibasis_sublattice = false
```

Use `fix_multibasis_sublattice = true` whenever a multibasis prototype has a
site that must not be optimized by SCRAPS. SCRAPS and `mcsqs` use the same
phase, composition, fitting, and output settings.

### TDB reuse and refitting

The TDB controls distinguish new calculations, unchanged cached databases,
and parameter-only refits:

```toml
[tdb]
run_tdb = true
skip_existing_tdb = true
refit_existing_tdb = true
skip_existing_plots = true
generate_gibbs_energy = true
generate_gibbs_mixing = true
generate_phase_diagram = true
generate_combined_phase_diagram = true
generate_contcar_plots = true
```

Each `generate_*` switch independently enables its output.
`skip_existing_plots` is evaluated per output, so an existing PNG is retained
while another enabled, missing PNG is generated. It also applies independently
to the combined phase-diagram image and each combined CONTCAR image. Systems
with more or fewer than two variable elements do not produce the three binary
Gibbs/phase summary plots.

The separate `[phase_plots]` stage processes only systems selected by the TDB
element pools and primary/secondary min/max counts, rather than every cached
TDB under `Files/Comps`. To plot only explicit systems, override that automatic
selection:

```toml
[phase_plots]
systems = ["CrHfMo"]
```

The same system filter is applied to the compiled phase-animation grid.

| Existing `.tdb` | `skip_existing_tdb` | `refit_existing_tdb` | Result |
|---|---:|---:|---|
| No | either | either | Full relaxation and TDB fit |
| Yes | `false` | either | Delete and perform the full calculation again |
| Yes | `true` | `false` | Leave the existing TDB unchanged |
| Yes | `true` | `true` | Reuse relaxed structures and rerun only `sqs2tdb -fit/-tdb` |

Refitting rewrites configured `terms.in` and `mult.in`, so changing those
files does not require rerunning the MLIP relaxations.

### Per-size model inputs

`terms.in` and `mult.in` may use a phase base key or a generated phase key.
The generated key takes precedence for that system size:

```toml
[tdb_inputs.terms_in]
PHASE1 = """fallback terms for every size
"""
PHASE1_2 = """binary terms
"""
PHASE1_3 = """ternary terms
"""

[tdb_inputs.mult_in]
PHASE1_2 = "a=1\tb=2"
PHASE1_3 = "a=1\tb=3"
```

For example, a ternary `PHASE1_3` fit uses `PHASE1_3`; a quaternary
`PHASE1_4` fit falls back to `PHASE1` when no exact key exists.

### Lattice averaging and sublattices

When `phase.use_average_lattice = true`, `a`, `b`, and `c` are averaged over
every unique element in `tdb.primary_elements + tdb.secondary_elements`.
Periodic-table reference cells and covalent-diameter estimates provide an
automatic starting geometry for all element symbols. Optional example-family
presets and TOML overrides can replace those estimates. The TOML does not need
one table per element.

```toml
[phase]
use_average_lattice = true
```

When `average_elements` is omitted, every primary and secondary element is
used with the default `a:a`, `b:b`, `c:c` mapping. An optional dictionary can
select elements and map each output axis to a source reference axis:

```toml
average_elements = { elements = ["Cr", "Hf"], a = "a", b = "a", c = "a" }
```

This example averages the Cr/Hf reference `a` values for all three output
parameters. Valid source values are `"a"`, `"b"`, and `"c"`. The optional
`lattice_parameter_family` selects a named preset when one is available;
otherwise periodic-table values are used. Optional
`[phase.element_lattice_parameters.<element>]` tables may override built-in
values but are not required.

Optional ATAT/model inputs are configured through `[tdb_inputs.terms_in]`,
`[tdb_inputs.mult_in]`, `[tdb_inputs.sqsgen_in]`,
`[tdb_inputs.sublattice_map.*]`, and `[tdb_inputs.fixed_compositions]`.
The site form combines element order and fixed composition:

```toml
[tdb_inputs.sites.PHASE1.a]
elements = ["Zr", "Hf", "Ti", "Nb", "V"]
constant = true
composition = [0.20, 0.20, 0.20, 0.25, 0.15]
```

`constant = true` prevents element/fraction permutations. System-specific
overrides use the same syntax under `tdb_system_overrides`; names are order
independent, so `CrHf` and `HfCr` address the same system. Override fractions
are automatically included in SQS generation.

---

## Quick Start

### 1. Generate compositions

```python
from blade.tools.blade_compositions import BladeCompositions

composer = BladeCompositions(
    primary_elements=["Cr", "Hf", "Ta", "Ti", "Mo"],
    secondary_elements=[],
    primary_min=2, primary_max=3,
    secondary_min=0, secondary_max=0,
)
composition_list = composer.generate_compositions()
unique_len_comps  = composer.get_systems()
```

### 2. Generate SQS structures

```python
from blade.tools.blade_sqsgen import BladeSQS

sqs_gen = BladeSQS(
    phases_dict=phases["PHASE1"],
    sqsgen_levels=sqsgen_levels,
    level=5,
    len_comp=2,
    skip_existing_sqs=True,
)
params = mcsqs_params | {"super_cell_size": (4, 3, 2)}
sqs_gen.sqs_gen(phase=phase, paths=paths, params=params)
```

**`mcsqs_params` keys:**

| Key | Description |
|---|---|
| `"2"` | Pair cutoff — nearest-neighbor shell index (resolved to Å by `BladeCutoff`) |
| `"3"` | Triplet cutoff shell index (`0` to disable) |
| `"4"` | Quadruplet cutoff shell index (`0` to disable) |
| `"cutoff_mode"` | `"nn"` (NN shell index) or `"angstrom"` (direct Å value) |
| `"time"` | Seconds per sqsdb directory |
| `"parallel_runs"` | Number of simultaneous `mcsqs` processes |

**`sqsgen_levels` format:**

```python
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]],                    "letter": ["a"]},
    {"level": 1, "compositions": [[0.5, 0.5]],                     "letter": ["a"]},
    {"level": 2, "compositions": [[0.75, 0.25]],                   "letter": ["a"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]], "letter": ["a"]},
]
```

### 2a. Generate SQS with SCRAPS

```python
from blade.tools.blade_scraps_gen import ScrapsSQSGen

sqs_gen = ScrapsSQSGen(
    phases_dict=phases["PHASE1"],
    sqsgen_levels=sqsgen_levels,
    level=2,
    len_comp=5,
    scraps_bin=SCRAPS_BIN,
    scraps_tools=SCRAPS_TOOLS,
    ranks=10,
    auto_budget=3,
    max_shellnum=2,
)
sqs_gen.sqs_gen(phase=specific_phase, paths=paths, params=params)
```

### 3. Fit thermodynamic databases

```python
from blade.tools.blade_tdb_gen import BladeTDBGen

gen = BladeTDBGen(
    phases=phase_list,
    phases_dict=phases,
    liquid=False,
    paths=paths,
    composition_list=composition_list,
    level=5,
    skip_existing=False,
    terms_in={"PHASE1": "1,0:1,0\n2,2:1,0\n"},
    sublattice_map={"PHASE1": {"a": ["Cr", "Hf"], "b": ["Al"]}},
    tdb_params={
        "fmax": 1e-4,
        "calculator": "grace",
        "calculator_kwargs": {"steps": 1000, "device": "cuda"},
        "t_min": 298.15,
        "t_max": 10000.0,
        "sro": False,
        "bv": 1e-3,
    },
)
gen.fit()
```

### 4. Oxidation screening

The bundled oxidation workflow has two preparation layers:

```bash
# Composition, Materials Project, MLIP, and system-structure preparation
python examples/oxidation/database_framework.py

# Batch or one fixed initial composition
python examples/oxidation/oxidation_framework.py
```

The complete TOML runner executes both layers when `database` and `oxidation`
are enabled. Set `oxidation.mode = "batch"` for all detected systems or
`"single"` for the exact system and composition in `[oxidation_single]`.

Direct batch API:

```python
import numpy as np
from blade.oxidation import OxidationCalculator, OxidationConfig

config = OxidationConfig(
    files_dir="/path/to/BLADE/Files",
    phase_element="B",                 # None for a metal-only mixed phase
    phase_element_stoichiometry=2.0,
    mixed_phase_subdir="blade",
    fixed_phases_subdir="ORB",
    region_label_mode="phases",        # or "id"
    slice_axis_priority=["Cr", "Hf", "Mo"],
)

calculator = OxidationCalculator(config)
calculator.run_batch(
    elements=["Cr", "Hf", "Mo"],       # empty list includes all systems
    run_calculations=True,              # calculate/cache every system first
    run_plots=True,                     # then plot strictly from those caches
    run_scan=True,
    run_muO_x_map=True,
    run_onset=True,
    run_3d_onset=True,
    run_composition_slices=True,
    run_composition_slice_muT=True,
    slice_muT_comp_step=0.1,
    run_animations=True,
    skip_if_tables_exist=True,
    skip_if_analysis_exists=True,
    temperature_values=np.arange(200, 2001, 200),
    mu_O_values=np.arange(-10.0, -3.999, 0.1),
    slice_remainder_ratios=[0.0, 0.25, 1/3, 0.5, 2/3, 0.75, 1.0],
)
```

With both execution flags enabled, oxidation runs in two global passes: all
selected systems are calculated first, then all figures and animations are
generated from cached results. Set `run_plots=False` to populate caches only,
or set `run_calculations=False` to redraw existing results without LP solves.
Plot-only mode raises an error identifying any required cache that is missing.
A failed calculation pass reports all failed systems and stops before the
plotting pass begins.

Direct single-composition API:

```python
analyzer = calculator.single(
    system="CrHfMoB",
    metals=["Cr", "Hf", "Mo"],
    composition=[1/3, 1/3, 1/3],
    include_0p01_to_0p05_components=True,
)
analyzer.run(
    T_values=np.arange(200, 2001, 10),
    mu_O_values=np.arange(-12.0, -3.999, 0.1),
    scan_T=[700, 1000, 1273, 1300],
    x_values=np.arange(0.0, 1.001, 0.01),
    run_calculations=True,
    run_plots=True,
    skip_if_exists=True,
)
```

`skip_if_tables_exist` reuses prepared phase and interaction tables.
`skip_if_analysis_exists` reuses compatible map/NPZ data and regenerates plots
without repeating equilibrium solves. Cache compatibility includes the grids,
temperatures, chemical-potential values, chemistry, and composition settings;
changing those settings triggers the required recalculation.

#### Maps and labels

- Oxidation scans include phase-fraction and average-composition plots in a
  separate temperature folder.
- `muO-x` and `muO-T` assemblage maps group regions only when their phase
  components change; internal fraction variation does not create duplicate
  regions.
- Binary and higher-dimensional composition slices use one selected axis plus
  configurable remainder ratios. `slice_axis_priority` chooses the first
  available preferred element and avoids generating all axis permutations.
- Side-by-side slice figures include matching phase diagrams and miscibility
  ranges when those inputs are available.
- Onset/AUC, ternary/higher-dimensional surfaces, GIFs, and MP4s are controlled
  independently by the batch flags. MP4 output requires `ffmpeg`.
- Each onset chemical-potential sweep ends at the first infeasible assemblage.
  Its table records the last feasible and first infeasible values; valid
  no-onset curves use the last feasible value rather than the configured
  `muO_max` endpoint.
- Onset tables also integrate the exact parent multicomponent-phase and oxide
  fractions over the feasible chemical-potential interval. Parent-phase
  composition may change or split into multiple compositions; it retains its
  identity while every configured variable element remains above 0.01.
- `region_label_mode = "id"` uses numbered regions and a key.
  `region_label_mode = "phases"` places formatted phase lines inside each
  region and removes the redundant key.
- Static PNG labels round phase-fraction and phase-composition ranges for
  readability. Interactive HTML maps retain the exact pointwise phase
  fractions and compositions in hover text, including small nonzero values.
- `include_0p01_to_0p05_components` controls whether small composition
  components participate in displayed component identities. Components at or
  below the configured edge tolerance are reduced to the corresponding end
  member consistently for every phase and system size.

### 5. Visualize results

```python
from blade.analysis.blade_visual import BladeVisualizer

viz = BladeVisualizer()
viz.plot_gibbs_energy(
    tdb=tdb, metals=["Cr", "Hf"], phase="PHASE1_2",
    fixed_species={"B": 2/3},
    temperatures=[300, 1000, 2000, 3000, 4000],
    output_path="CrHf/CrHf_Gibbs_Energy.png",
)
viz.plot_gibbs_mixing(
    tdb=tdb, metals=["Cr", "Hf"], phase="PHASE1_2",
    fixed_species={"B": 2/3},
    temperatures=[300, 1000, 2000, 3000, 4000],
    output_path="CrHf/CrHf_Gibbs_Mixing.png",
)
viz.plot_binary_phase_diagram(
    tdb=tdb, metals=["Cr", "Hf"], phases=["PHASE1_2"],
    fixed_species={"B": 2/3},
    temperature_range=(300, 4500, 50),
    output_path="CrHf/CrHf_Phase_Diagram.png",
)
```

---

## Phase Prototype Format

```python
phases = {
    "PHASE1": {
        "a": 3.10, "b": 3.10, "c": 13.80,
        "alpha": 90, "beta": 90, "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.333333 0.666667 0.083333 a\n"   # lowercase = variable sublattice
            "0.333333 0.666667 0.750000 b\n"
            "0.000000 0.000000 0.000000 C\n"   # uppercase = fixed element
        ),
    }
}
```

- Lowercase letter (`a`, `b`) — variable sublattice; mixed occupancy
- Uppercase symbol (`C`, `Al`) — fixed element; never substituted

---

## Example Scripts

Start with the reproducible [examples run guide](examples/README.md). Detailed
step-by-step guides cover [structures and TDB fitting](examples/structures/README.md),
[oxidation](examples/oxidation/README.md), and
[additional visualization tools](examples/extra/README.md). The complete
technical reference is [docs/BLADE_framework_manual.tex](docs/BLADE_framework_manual.tex).

### Step-by-step (`examples/structures/`)

| Script | Description |
|---|---|
| [`01_compositions.py`](examples/structures/01_compositions.py) | Enumerate compositions |
| [`02_sqs_generation.py`](examples/structures/02_sqs_generation.py) | Generate SQS structures with ATAT `mcsqs` |
| [`03_tdb_generation.py`](examples/structures/03_tdb_generation.py) | Relax structures and fit TDB databases |
| [`04_visualization.py`](examples/structures/04_visualization.py) | Plot Gibbs energy and phase diagrams |
| [`05_data_analysis.py`](examples/structures/05_data_analysis.py) | Extract volumes, lattice parameters, and energies from POSCAR files |

### Oxidation screening (`examples/oxidation/`)

| Script | Description |
|---|---|
| [`01_compositions.py`](examples/oxidation/01_compositions.py) | Enumerate oxide-relevant compositions |
| [`02_poscars.py`](examples/oxidation/02_poscars.py) | Download Materials Project structures |
| [`03_energy.py`](examples/oxidation/03_energy.py) | Relax MP structures with MLIP |
| [`04_database.py`](examples/oxidation/04_database.py) | Build unified energy database |
| [`database_framework.py`](examples/oxidation/database_framework.py) | Resumable composition, MP, MLIP, and database workflow |
| [`oxidation_framework.py`](examples/oxidation/oxidation_framework.py) | Generalized batch or single-composition calculations and plots |

### Combined inputs (`examples/`)

| File | Description |
|---|---|
| [`full_framework.py`](examples/full_framework.py) | Validate and run all enabled stages from TOML |
| [`full_framework.toml`](examples/full_framework.toml) | Single-variable-sublattice workflow example |
| [`full_framework_max.toml`](examples/full_framework_max.toml) | Multiple-variable-sublattice workflow example |

### Extra utilities (`examples/extra/`)

| Script | Description |
|---|---|
| [`phase_visualization.py`](examples/extra/phase_visualization.py) | Combined phase-map/GIF generation and synchronized grid-video compilation; controlled by `run_phase_plots` and `run_phase_grid` |
| [`refit_tdb.py`](examples/extra/refit_tdb.py) | Refit existing TDB from relaxed structures with updated `terms.in` / `mult.in` — skips relaxation |
| [`compare_energy.py`](examples/extra/compare_energy.py) | Compare energy per atom between any two folders of relaxed structures; parity and residual plots |
| [`tdb_gen_alloy.py`](examples/extra/tdb_gen_alloy.py) | HPC driver for pure metallic alloy systems (no fixed sublattice) |
| [`tdb_gen_hedb_quaternary.py`](examples/extra/tdb_gen_hedb_quaternary.py) | Legacy quaternary prototype demonstration |

### HPC driver scripts (`examples/structures/`)

| Script | System | SQS method |
|---|---|---|
| [`tdb_gen_hedb.py`](examples/structures/tdb_gen_hedb.py) | Example A: one variable and one fixed sublattice | mcsqs |
| [`tdb_gen_hedb_scraps.py`](examples/structures/tdb_gen_hedb_scraps.py) | Example A: one variable and one fixed sublattice | SCRAPS |
| [`tdb_gen_max.py`](examples/structures/tdb_gen_max.py) | Example B: multiple variable and fixed sublattices | mcsqs |
| [`tdb_gen_max_scraps.py`](examples/structures/tdb_gen_max_scraps.py) | Example B: multiple variable and fixed sublattices | SCRAPS |

---

## Quality Checks

The repository pins its lint and formatting tools in the Pixi environment:

```bash
pixi run check
pixi run python -m compileall -q src examples
pixi run python examples/full_framework.py --check
pixi run python examples/full_framework.py --dry-run

# Equivalent uv checks
uv run ruff check src examples
uv run black --check src examples
uv run isort --check-only src examples
uv run --extra workflow python -m compileall -q src examples
```

`pixi run check` runs Ruff, Black's check mode, and isort's check mode. The
framework checks validate TOML structure, dependencies, stage compatibility,
and configured paths without starting calculations.

## License

Distributed under the GPL-3.0-or-later License.

## Contact

Chase Katz — [ichasekatz@tamu.edu](mailto:ichasekatz@tamu.edu)
