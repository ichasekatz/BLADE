<div align="center">

# BLADE

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
[![Python](https://img.shields.io/badge/python-3.10%2B-brightgreen.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-1.5.0-orange.svg)](pyproject.toml)

**Batch Lattice Analysis and Discovery Engine — automated CALPHAD thermodynamic database generation for multicomponent materials systems.**

[Report a Bug](https://github.com/ichasekatz/BLADE/issues/new?labels=bug) · [Request a Feature](https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement)

</div>

---

## Overview

BLADE automates the full CALPHAD workflow from start to finish: given an element pool and a crystal structure prototype, it enumerates every valid N-component system, generates SQS supercells, relaxes them with an ML interatomic potential, fits CALPHAD parameters, and produces `.tdb` files — without manual intervention at any step.

Inspired by [MaterialsFramework](https://github.com/dogusariturk/MaterialsFramework)'s modular design, BLADE builds on top of it as a computational backend. **BLADE requires the [ichasekatz fork of MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework)** — install it before installing BLADE.

- **Composition generation** — enumerate binary, ternary, or N-component systems from primary and secondary element pools with configurable count bounds
- **SQS generation (mcsqs)** — drive ATAT `mcsqs` in parallel across all compositions; supports per-sublattice element assignment and fixed-composition constraints for multi-sublattice structures
- **SQS generation (SCRAPS)** — generate SQS structures using the SCRAPS algorithm as an alternative to `mcsqs`, with support for multi-basis variable sublattices
- **TDB fitting** — relax structures with any MLIP supported by MaterialsFramework, fit CALPHAD parameters, write `.tdb` files
- **Fit customization** — override ATAT-generated `terms.in`, `mult.in`, and sublattice occupancy maps per phase; filter permutation directories by element assignment per composition level
- **Visualization** — plot Gibbs energy and pseudo-binary phase diagrams, combine structure images, render relaxation movies
- **Volume analysis** — scan POSCAR files and extract lattice parameters and per-atom volumes into a DataFrame

---

## BLADE Capabilities

| Capability | BLADE |
|---|---|
| Structure relaxation | Batched across every composition and SQS level |
| Composition enumeration | `BladeCompositions` traverses the full N-component search space |
| SQS generation (mcsqs) | `BladeSQS` runs parallel `mcsqs`; exact compositions from `sqsgen_levels` |
| SQS generation (SCRAPS) | `ScrapsSQSGen` runs SCRAPS; handles multi-basis sublattices |
| Multi-sublattice customization | `sublattice_map` assigns elements per sublattice; `composition_elements` controls permutations per composition level |
| Fixed-composition constraints | `fixed_compositions` pins sublattice fractions, bypassing mcsqs permutation generation |
| CALPHAD database output | `BladeTDBGen` processes all systems and phases in one call |
| TDB fit accuracy | Per-phase `terms_in`/`mult_in` override ATAT cluster expansion files |
| Fixed sublattice sites | Auto-derived from uppercase element labels in `coords` |
| Phase diagrams | Gibbs energy, Gibbs energy of mixing, and pseudo-binary phase diagrams via pycalphad |
| Visualization | `BladeVisualizer` plots thermodynamic quantities and combines structures and relaxation movies |

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
| `BLADEVolume` | Extract volumes and lattice parameters from POSCAR files |

All classes are available via lazy top-level import:

```python
from blade import BladeCompositions, BladeSQS, ScrapsSQSGen, BladeTDBGen
from blade import BladeVisualizer, BLADEVolume
```

---

## Getting Started

### Prerequisites

- [ichasekatz/MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) — required fork, install before BLADE
- [ATAT](https://axelvandewalle.github.io/www-avdw/atat/) — `mcsqs`, `sqs2tdb`, `corrdump` must be on `$PATH`

### Installation

Install dependencies with [pixi](https://pixi.sh) and pip:

```bash
pixi add git python==3.12 cmake pip psutil pandarallel orb-models pytorch tcsh pymatgen sqsgenerator
pip install matplotlib pillow imageio pycalphad torch ase pynanoflann
pip install tensorflow tensorpotential
```

Install the required MaterialsFramework fork:

```bash
git clone https://github.com/ichasekatz/MaterialsFramework.git
cd MaterialsFramework && pip install -e . && cd ..
```

Install BLADE:

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE && pip install -e .
```

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
sqs_gen.sqs_gen(phase=phase, paths=paths, iter=1_000_000, params=params)
```

**`mcsqs_params` keys:**

| Key | Description |
|---|---|
| `"2"` | Pair cutoff — nearest-neighbor shell index (resolved to Å by `BladeCutoff`) |
| `"3"` | Triplet cutoff shell index (`0` to disable) |
| `"4"` | Quadruplet cutoff shell index (`0` to disable) |
| `"cutoff_mode"` | `"nn"` (default, NN shell index) or `"angstrom"` (direct Å value) |
| `"time"` | Seconds per sqsdb directory |
| `"use_time"` | `True` to use time-based stopping criterion |
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

Each `"letter"` entry names the variable sublattice sites for that level. Fixed sites (uppercase in `coords`) are excluded automatically. Compositions listed are used exactly — no expansion to other fractions at the same denominator.

### 2a. Generate SQS structures with SCRAPS

```python
from blade.tools.blade_scraps_gen import ScrapsSQSGen

sqs_gen = ScrapsSQSGen(
    phases_dict=phases["PHASE1"],
    sqsgen_levels=sqsgen_levels,
    level=2,
    len_comp=5,
    sublattice_map=sublattice_map,
    fixed_compositions=fixed_compositions,
    scraps_bin=SCRAPS_BIN,
    scraps_tools=SCRAPS_TOOLS,
    ranks=10,
    auto_budget=3,
    max_shellnum=2,
    fix_multibasis_sublattice=True,
)
sqs_gen.sqs_gen(phase=specific_phase, paths=paths, iter=0, params=params)
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
    fixed_compositions={"a": [0.5, 0.5]},
    tdb_params={
        "fmax": 1e-4,
        "calculator": "grace",
        "calculator_kwargs": {"steps": 1000, "device": "cuda"},
        "t_min": 298.15,
        "t_max": 10000.0,
        "sro": False,
        "bv": 1e-3,
        "phonon": False,
        "open_calphad": False,
    },
)
gen.fit()
```

### 4. Visualize results

```python
from blade.analysis.blade_visual import BladeVisualizer

viz = BladeVisualizer()
viz.plot_gibbs_energy(
    tdb=tdb, metals=["Cr", "Hf"], phase="PHASE1_2",
    fixed_species={"Al": 1/3},
    temperatures=[300, 1000, 2000, 3000, 4000],
    output_path="CrHf/CrHf_Gibbs_Energy.png",
)
viz.plot_binary_phase_diagram(
    tdb=tdb, metals=["Cr", "Hf"], phases=["PHASE1_2"],
    fixed_species={"Al": 1/3},
    temperature_range=(300, 4500, 50),
    output_path="CrHf/CrHf_Phase_Diagram.png",
)
viz.phase_diagram(pngs, save="Combined_Phase_Diagrams.png")
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

**Label convention in `coords`:**
- Lowercase letter (`a`, `b`, `c`) — variable sublattice; mixed occupancy; driven by `sqsgen_levels`
- Uppercase element symbol (`C`, `Al`, `N`) — fixed element; never substituted

`BladeTDBGen` auto-derives fixed-sublattice definitions from the uppercase labels in `coords`.

---

## Example Scripts

### Step-by-step (`examples/structure/`)

| Script | Description |
|---|---|
| [`01_compositions.py`](examples/structure/01_compositions.py) | Enumerate compositions |
| [`02_sqs_generation.py`](examples/structure/02_sqs_generation.py) | Generate SQS structures with ATAT `mcsqs` |
| [`03_tdb_generation.py`](examples/structure/03_tdb_generation.py) | Relax structures and fit TDB databases |
| [`04_visualization.py`](examples/structure/04_visualization.py) | Plot Gibbs energy and phase diagrams |
| [`05_volume_analysis.py`](examples/structure/05_volume_analysis.py) | Extract per-atom volumes from POSCAR files |

### HPC driver scripts (`examples/structure/`)

| Script | System | SQS method |
|---|---|---|
| [`tdb_gen_hedb.py`](examples/structure/tdb_gen_hedb.py) | Hexagonal two-sublattice phase | mcsqs |
| [`tdb_gen_hedb_scraps.py`](examples/structure/tdb_gen_hedb_scraps.py) | Hexagonal two-sublattice phase | SCRAPS |
| [`tdb_gen_max.py`](examples/structure/tdb_gen_max.py) | 211 MAX phase | mcsqs |
| [`tdb_gen_max_scraps.py`](examples/structure/tdb_gen_max_scraps.py) | 211 MAX phase | SCRAPS |

---

## License

Distributed under the GPL-3.0-or-later License.

## Contact

Chase Katz — [ichasekatz@tamu.edu](mailto:ichasekatz@tamu.edu)
