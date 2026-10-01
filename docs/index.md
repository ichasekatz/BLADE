# BLADE

**Automated CALPHAD TDB generation and oxidation-stability analysis for multi-component borides and alloys.**

BLADE (Boride Lattice Automated Discovery Engine) is a Python pipeline for
high-throughput thermodynamic database construction and oxidation analysis
of transition-metal diboride alloys. Given a pool of candidate elements,
BLADE generates SQS structures, relaxes them with a machine-learning
interatomic potential (MLIP), fits CALPHAD Gibbs energy models, and
computes oxygen chemical-potential stability windows across temperature and
composition space.

## Pipeline

| Stage | What it does |
|---|---|
| **BladeCompositions** | Enumerate all N-component element combinations |
| **BladeSQS** | Write ATAT input files, run `mcsqs` in parallel |
| **BladeTDBGen** | MLIP-relax SQS structures, fit CALPHAD TDB |
| **OxidationCalculator** | Build oxide database, compute μ_O–T stability maps |

## Quick Start

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE && pixi install
pixi run python examples/full_framework.py examples/full_framework.toml
```

See [Getting Started](usage/getting-started.md) for the full walkthrough.
