<div align="center">

# BLADE

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
[![Python](https://img.shields.io/badge/python-3.10%2B-brightgreen.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-1.0.0-orange.svg)](pyproject.toml)

**Batch Lattice Analysis and Discovery Engine — automated CALPHAD thermodynamic database generation for multicomponent materials systems.**

[Report a Bug](https://github.com/ichasekatz/BLADE/issues/new?labels=bug) · [Request a Feature](https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement)

</div>

---

## Overview

BLADE automates CALPHAD database generation for multicomponent materials: given an element pool and a crystal structure prototype, it enumerates N-component systems, generates SQS supercells, and fits thermodynamic parameters.

**BLADE requires the [ichasekatz fork of MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework)** — install it before installing BLADE.

- **Composition generation** — enumerate N-component systems from primary and secondary element pools
- **SQS generation** — drive ATAT `mcsqs` in parallel across all compositions
- **Phase diagram output** — compute and plot Gibbs energy phase diagrams via pycalphad

---

## Modules

| Class | Description |
|---|---|
| `BladeCompositions` | Enumerate compositions from primary and secondary element pools |
| `BladeSQS` | Generate SQS structures using ATAT `mcsqs` |
| `BladePhaseDiagram` | Compute and plot phase diagrams |

---

## Getting Started

### Prerequisites

- [ichasekatz/MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) — required fork
- [ATAT](https://axelvandewalle.github.io/www-avdw/atat/) — `mcsqs` must be on `$PATH`

### Installation

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE && pip install -e .
```

---

## Usage

```python
from src.blade.blade_compositions import BladeCompositions
from src.blade.blade_sqs import BladeSQS

composer = BladeCompositions(
    primary_elements=["Cr", "Hf", "Ta", "Ti", "Mo"],
    secondary_elements=[],
    primary_min=2, primary_max=3,
    secondary_min=0, secondary_max=0,
)
composition_list = composer.generate_compositions()
```

---

## License

Distributed under the GPL-3.0-or-later License.

## Contact

Chase Katz — [ichasekatz@tamu.edu](mailto:ichasekatz@tamu.edu)
