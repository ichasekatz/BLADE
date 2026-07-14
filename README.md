<div align="center">

# BLADE

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
[![Python](https://img.shields.io/badge/python-3.10%2B-brightgreen.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-1.6.0-orange.svg)](pyproject.toml)

**Batch Lattice Analysis and Discovery Engine — automated CALPHAD thermodynamic database generation and oxidation screening for multicomponent materials systems.**

[Report a Bug](https://github.com/ichasekatz/BLADE/issues/new?labels=bug) · [Request a Feature](https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement)

</div>

---

## Overview

BLADE automates the full CALPHAD workflow: given an element pool and a crystal structure prototype, it enumerates every valid N-component system, generates SQS supercells, relaxes them with an ML interatomic potential, fits CALPHAD parameters, and produces `.tdb` files — without manual intervention at any step. An integrated oxidation screening module extends the pipeline to map phase stability as a function of composition and oxygen chemical potential.

BLADE builds on [MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) as a computational backend. **The [ichasekatz fork](https://github.com/ichasekatz/MaterialsFramework) is required** — install it before BLADE.

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

All classes are available via lazy top-level import:

```python
from blade import BladeCompositions, BladeSQS, ScrapsSQSGen, BladeTDBGen
from blade import BladeVisualizer, BLADEData
from blade.oxidation import OxideCompositions, OxideDatabase
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
| `"cutoff_mode"` | `"nn"` (NN shell index) or `"angstrom"` (direct Å value) |
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

The oxidation workflow is driven by sequential scripts in `examples/oxidation/`. Each script reads output from the previous step.

```bash
# Step 1 — enumerate oxide-relevant compositions and scan BLADE CONTCARs
python examples/oxidation/01_compositions.py

# Step 2 — download MP reference structures
python examples/oxidation/02_poscars.py

# Step 3 — relax MP structures with MLIP
python examples/oxidation/03_energy.py

# Step 4 — build unified energy database and run grand potential minimization
python examples/oxidation/04_database.py

# Step 5 — generate oxidation onset maps and phase assemblage plots
python examples/oxidation/05_oxidation_calc.py
```

The core classes used internally:

```python
from blade.oxidation.compositions import OxideCompositions
from blade.oxidation.database import OxideDatabase

# OxideCompositions enumerates all oxide-relevant compositions
comp = OxideCompositions(
    files_dir=files_dir,
    primary_elements=["Cr", "Hf", "Ta", "Ti", "Mo", "W", "Nb", "V", "Zr"],
    fixed_elements=frozenset({"B"}),
    primary_min=0, primary_max=3,
    include_fixed=True,
    include_added_oxygen=True,
    include_fixed_oxygen=True,
)
comp.run()

# OxideDatabase compiles energies, builds the database, and runs screening
db = OxideDatabase(files_dir=files_dir, ...)
db.scan_blade_data(blade_comp_dir=files_dir / "Comps")
db.run()
```

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
| [`05_oxidation_calc.py`](examples/oxidation/05_oxidation_calc.py) | Grand potential minimization and phase maps |
| [`06_visualization.py`](examples/oxidation/06_visualization.py) | Plot oxidation onset and phase assemblage maps |

### Extra utilities (`examples/extra/`)

| Script | Description |
|---|---|
| [`plot_phase.py`](examples/extra/plot_phase.py) | N-component phase maps with G_mix convex hull and shaded miscibility gaps; animated GIF per system |
| [`compile_phase_grid.py`](examples/extra/compile_phase_grid.py) | Compile per-system GIFs into a synchronized grid video (MP4) grouped by system size |
| [`refit_tdb.py`](examples/extra/refit_tdb.py) | Refit existing TDB from relaxed structures with updated `terms.in` / `mult.in` — skips relaxation |
| [`compare_energy.py`](examples/extra/compare_energy.py) | Compare energy per atom between any two folders of relaxed structures; parity and residual plots |
| [`tdb_gen_alloy.py`](examples/extra/tdb_gen_alloy.py) | HPC driver for pure metallic alloy systems (no fixed sublattice) |
| [`tdb_gen_hedb_quaternary.py`](examples/extra/tdb_gen_hedb_quaternary.py) | HPC driver for quaternary diboride systems |

### HPC driver scripts (`examples/structures/`)

| Script | System | SQS method |
|---|---|---|
| [`tdb_gen_hedb.py`](examples/structures/tdb_gen_hedb.py) | Hexagonal two-sublattice (diboride) | mcsqs |
| [`tdb_gen_hedb_scraps.py`](examples/structures/tdb_gen_hedb_scraps.py) | Hexagonal two-sublattice (diboride) | SCRAPS |
| [`tdb_gen_max.py`](examples/structures/tdb_gen_max.py) | 211 MAX phase | mcsqs |
| [`tdb_gen_max_scraps.py`](examples/structures/tdb_gen_max_scraps.py) | 211 MAX phase | SCRAPS |

---

## License

Distributed under the GPL-3.0-or-later License.

## Contact

Chase Katz — [ichasekatz@tamu.edu](mailto:ichasekatz@tamu.edu)
