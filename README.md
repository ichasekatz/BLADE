<div align="center">

# BLADE

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
[![Python](https://img.shields.io/badge/python-3.10%2B-brightgreen.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-1.4.0-orange.svg)](pyproject.toml)

**Batch Lattice Analysis and Discovery Engine — automated CALPHAD thermodynamic database generation for multicomponent materials systems.**

[Report a Bug](https://github.com/ichasekatz/BLADE/issues/new?labels=bug) · [Request a Feature](https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement)

</div>

---

## Overview

BLADE automates the full CALPHAD workflow from start to finish: given an element pool and a crystal structure prototype, it enumerates every valid N-component system, generates SQS supercells, relaxes them with an ML interatomic potential, fits CALPHAD parameters, and produces `.tdb` files — without manual intervention at any step.

Inspired by [MaterialsFramework](https://github.com/dogusariturk/MaterialsFramework)'s modular design, BLADE builds on top of it as a computational backend. **BLADE requires the [ichasekatz fork of MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework)** — install it before installing BLADE.

- **Composition generation** — enumerate binary, ternary, or N-component systems from primary and secondary element pools with configurable count bounds
- **SQS generation** — drive ATAT `mcsqs` in parallel across all compositions and supercell sizes; supports per-sublattice element assignment for multi-sublattice structures
- **TDB fitting** — relax structures with any MLIP supported by MaterialsFramework, compute formation energies and thermal properties, fit CALPHAD parameters, write `.tdb` files; copies SQS structures from the staging tree into the ATAT sqsdb directory before each fit
- **Fit customization** — override ATAT-generated `terms.in`, `mult.in`, and sublattice occupancy maps per phase prototype to maximize CALPHAD accuracy
- **Visualization** — plot Gibbs energy and pseudo-binary phase diagrams, combine structure images, render relaxation movies
- **Volume analysis** — scan POSCAR files and extract lattice parameters and per-atom volumes into a DataFrame

---

## BLADE Capabilities

Existing tools handle pieces of this workflow in isolation — a calculator relaxes one structure, a fitting tool processes one system, a visualization tool reads one file. BLADE connects all of these into a single automated pipeline driven by a composition list.

| Capability | BLADE |
|---|---|
| Structure relaxation | Batched across every composition and SQS level, instead of a single structure at a time |
| Property calculation | Automated per-SQS-level pipeline, instead of per-structure manual invocation |
| Composition enumeration | `BladeCompositions` traverses the full N-component search space automatically |
| SQS generation | `BladeSQS` runs parallel `mcsqs` across all compositions simultaneously; canonical compositions generated from each level (all sorted-descending fractions sharing the same denominator, no mirrors) |
| Multi-sublattice customization | `sublattice_map` assigns specific elements to each sublattice per phase, enabling ordered intermetallics and structures where different sublattice sites are occupied by different species |
| CALPHAD database output | `BladeTDBGen` processes all systems and phases in one call, instead of one system per run |
| TDB fit accuracy | Per-phase `terms_in` and `mult_in` override ATAT-generated `terms.in` / `mult.in` cluster expansion files; `sublattice_map` controls element-to-sublattice assignment passed to the fitter |
| Fixed sublattice sites | `sqsgen_levels2` holds any species fixed on any sublattice, enabling ceramics and intermetallics |
| Phase diagrams | Gibbs energy, Gibbs energy of mixing, and pseudo-binary phase diagrams via `pycalphad`, instead of a separate workflow |
| Visualization | `BladeVisualizer` plots thermodynamic quantities and combines structures, diagrams, and relaxation movies in one place |

The multi-sublattice and fixed-sublattice support is particularly significant: BLADE is not limited to simple alloys or single-sublattice ceramics. Any crystal structure with multiple distinct sublattice sites — where different sites should be occupied by different elements — can be described through a `sublattice_map` per phase prototype, making BLADE applicable across a wide range of complex intermetallics, ceramics, and ordered compounds.

---

## Modules

### Tools (`blade.tools`)

| Class | Description |
|---|---|
| `BladeCompositions` | Enumerate compositions from primary and secondary element pools |
| `BladeSQS` | Generate SQS structures using ATAT `mcsqs` in parallel |
| `BladeTDBGen` | Relax structures and fit CALPHAD thermodynamic databases |
| `BladeCutoff` | Compute neighbor-shell distances from SQS supercell geometry |

### Analysis (`blade.analysis`)

| Class | Description |
|---|---|
| `BladeVisualizer` | Plot Gibbs energy, mixing energy, and phase diagrams |
| `BLADEVolume` | Extract volumes and lattice parameters from POSCAR files |

```python
from blade import BladeCompositions, BladeSQS, BladeTDBGen
from blade import BladeVisualizer, BLADEVolume
```

------

## Getting Started

### Prerequisites

- [ichasekatz/MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) — required fork, install before BLADE
- [ATAT](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/) (`mcsqs` must be on `$PATH`)

### Installation

Install dependencies with [pixi](https://pixi.sh) and pip:

```bash
pixi add git python==3.12 cmake pip psutil pandarallel orb-models pytorch tcsh pymatgen sqsgenerator
pip install matplotlib pillow imageio pycalphad torch ase pynanoflann
pip install tensorflow
pip install tensorpotential
```

Install the required MaterialsFramework fork:

```bash
git clone https://github.com/ichasekatz/MaterialsFramework.git
cd MaterialsFramework
pip install -e .
cd ..
```

Install BLADE:

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE
pip install -e .
```

---

## Quick Start

### 1. Generate compositions

```python
from blade.tools.blade_compositions import BladeCompositions

composer = BladeCompositions(
    primary_elements=["Hf", "Cr", "Ta", "Ti", "Mo"],
    secondary_elements=[],
    system_size=3,
    primary_min=3, primary_max=3,
    secondary_min=0, secondary_max=0,
    allow_lower_order=False,
)

composition_list = composer.generate_compositions()
print(f"{len(composition_list)} ternary systems: {composition_list[:5]}")
```

### 2. Generate SQS structures

```python
from blade.tools.blade_sqsgen import BladeSQS

sqs_gen = BladeSQS(
    phases_dict=phases["HEDB1"],
    sqsgen_levels=sqsgen_levels,
    level=5,
    len_comp=2,
    skip_existing_sqs=True,
    # sqsgen_in={"HEDB1": "..."},  # verbatim sqsgen.in override; omit to auto-generate
)
sqs_gen.sqs_gen(phase=phase, paths=paths, iter=1_000_000, params=mcsqs_params)
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
    # terms_in={"HEDB1": "1,0:1,0\n2,0:1,0\n"},   # override ATAT terms.in per phase
    # mult_in={"HEDB1": "..."},                     # override ATAT mult.in per phase
    # sublattice_map={"HEDB1": {"a": ["Cr", "Hf"], "b": ["Zr", "Ti"]}},  # per-sublattice elements
)
gen.fit()  # explicit call — no side effects on construction
```

### 4. Visualize results

```python
from blade.analysis.blade_visual import BladeVisualizer

viz = BladeVisualizer()
viz.contcar(contcars, save="CrHf/Combined_CONTCARs.png")
viz.plot_gibbs_energy(tdb, metals=["Cr", "Hf"], phase="HEDB1_2",
                      fixed_species={"B": 2/3}, temperatures=[300, 1000, 2000, 3000, 4000],
                      output_path="CrHf/CrHf_Gibbs_Energy.png")
viz.plot_gibbs_mixing(tdb, metals=["Cr", "Hf"], phase="HEDB1_2",
                      fixed_species={"B": 2/3}, temperatures=[300, 1000, 2000, 3000, 4000],
                      output_path="CrHf/CrHf_Gibbs_Mixing.png")
viz.plot_binary_phase_diagram(tdb, metals=["Cr", "Hf"], phases=["HEDB1_2"],
                              fixed_species={"B": 2/3}, temperature_range=(300, 4500, 50),
                              output_path="CrHf/CrHf_Phase_Diagram.png")
viz.phase_diagram(pngs, save="Combined_Phase_Diagrams.png")
```

---

## Example Scripts

### Core CALPHAD Pipeline

| Script | Description |
|---|---|
| [`01_compositions.py`](examples/01_compositions.py) | Enumerate binary/ternary/N-component compositions |
| [`02_sqs_generation.py`](examples/02_sqs_generation.py) | Generate SQS structures with ATAT `mcsqs` |
| [`03_tdb_generation.py`](examples/03_tdb_generation.py) | Relax structures and fit TDB databases |
| [`04_visualization.py`](examples/04_visualization.py) | Generate Gibbs energy and phase diagram plots, combine structures and movies |
| [`05_volume_analysis.py`](examples/05_volume_analysis.py) | Extract per-atom volumes from POSCAR files |
| [`tdb_gen.py`](examples/tdb_gen.py) | Full end-to-end workflow (HPC driver script) |

---

## License

Distributed under the GPL-3.0-or-later License.

## Contact

Chase Katz — [ichasekatz@tamu.edu](mailto:ichasekatz@tamu.edu)
