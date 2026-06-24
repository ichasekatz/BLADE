<div align="center">

# BLADE

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
[![Python](https://img.shields.io/badge/python-3.10%2B-brightgreen.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-1.3.0-orange.svg)](pyproject.toml)

**Batch Lattice Analysis and Discovery Engine — automated CALPHAD thermodynamic database generation for multicomponent materials systems.**

[Report a Bug](https://github.com/ichasekatz/BLADE/issues/new?labels=bug) · [Request a Feature](https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement)

</div>

---

<div align="center">
<img width="700" alt="BLADE_Framework" src="https://github.com/user-attachments/assets/8a64271d-0427-4e2e-b964-df0af7ff18d0" />
</div>

---

## Overview

BLADE automates the full CALPHAD workflow from start to finish: given an element pool and a crystal structure prototype, it enumerates every valid N-component system, generates SQS supercells, relaxes them with an ML interatomic potential, fits CALPHAD parameters, and produces `.tdb` files.

**BLADE requires the [ichasekatz fork of MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework)** — install it before installing BLADE.

- **Composition generation** — enumerate binary, ternary, or N-component systems from primary and secondary element pools
- **SQS generation** — drive ATAT `mcsqs` in parallel across all compositions; supports per-sublattice element assignment
- **TDB fitting** — relax structures with any MLIP supported by MaterialsFramework, fit CALPHAD parameters, write `.tdb` files
- **Fit customization** — override ATAT-generated `terms.in`, `mult.in`, and sublattice occupancy maps per phase
- **Visualization** — plot Gibbs energy and pseudo-binary phase diagrams, combine structure images, render relaxation movies
- **Volume analysis** — scan POSCAR files and extract lattice parameters and per-atom volumes

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

---

## Getting Started

### Prerequisites

- [ichasekatz/MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) — required fork, install before BLADE
- [ATAT](https://axelvandewalle.github.io/www-avdw/atat/) — `mcsqs`, `sqs2tdb`, `corrdump` must be on `$PATH`

### Installation

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE && pip install -e .
```

---

## Example Scripts (`examples/`)

| Script | Description |
|---|---|
| `01_compositions.py` | Enumerate compositions |
| `02_sqs_generation.py` | Generate SQS structures |
| `03_tdb_generation.py` | Fit TDB databases |
| `04_visualization.py` | Plot Gibbs energy and phase diagrams |
| `05_volume_analysis.py` | Extract per-atom volumes |
| `tdb_gen.py` | Full end-to-end HPC driver |

---

## License

Distributed under the GPL-3.0-or-later License.

## Contact

Chase Katz — [ichasekatz@tamu.edu](mailto:ichasekatz@tamu.edu)
