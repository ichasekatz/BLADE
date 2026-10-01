# Getting Started

## Prerequisites

**Python 3.12** via [pixi](https://pixi.sh) or [uv](https://docs.astral.sh/uv/).

**ATAT binaries** — `mcsqs`, `sqs2tdb`, and `corrdump` must be on `$PATH`.
They live in `PhaseForge/atat/bin/`; add that directory to your shell profile:

```bash
export PATH="/path/to/PhaseForge/atat/bin:$PATH"
```

Verify:

```bash
which mcsqs sqs2tdb corrdump
```

**MaterialsFramework** — the `ichasekatz` fork, installed as an editable package
into the BLADE environment. Required for MLIP relaxation and TDB fitting.

## Clone and Install

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE
pixi install          # or: uv sync
```

Install MaterialsFramework into the same environment:

```bash
cd ../PhaseForge/MaterialsFramework
pip install -e .
```

Confirm both packages are visible:

```bash
pixi run python -c "import blade; import materialsframework; print('ok')"
```

## Configure `full_framework.toml`

All settings live in a single TOML file. Copy the bundled example and edit it:

```bash
cp examples/full_framework.toml my_run.toml
```

### [paths]

Three paths are required. `files_dir` defaults to `blade_root/Files` and can be omitted.

```toml
[paths]
blade_root = "/path/to/BLADE"
sqsdb_dir  = "/path/to/PhaseForge/atat/data/sqsdb"
# files_dir defaults to blade_root/Files — override only if needed
```

`blade_root` must point to the BLADE working directory.
`sqsdb_dir` is the ATAT sqsdb base directory that MaterialsFramework reads and writes.

### [tdb] — element pool

The element pool is defined once under `[tdb]` and is the source of truth
for all composition enumeration. Every binary, ternary, … combination
up to `primary_max` elements is generated automatically.

```toml
[tdb]
primary_elements  = ["Cr", "Hf", "Mo"]
secondary_elements = []
primary_min   = 2
primary_max   = 3
secondary_min = 0
secondary_max = 0
```

The `[database]` and `[oxidation_batch]` sections have their own `primary_elements`
fields for independent element scopes (e.g., a broader oxide database covering
more metals than the TDB fit). Override them explicitly or mirror the `[tdb]` pool.

### [stages]

Enable only the stages you need. Each stage is independent — earlier stages
do not need to run in the same invocation as later ones.

```toml
[stages]
sqs_generation   = true   # Step 1: run mcsqs, write ATAT input trees
tdb_fitting      = true   # Step 2: MLIP relax + CALPHAD fit → .tdb files
phase_diagrams   = false  # Step 3: phase diagram plots and GIFs
oxide_database   = false  # Step 4: download MP phases, compute reference energies
oxidation_graphs = false  # Step 5: mu_O–T stability maps and oxidation-onset analysis
```

### [phase] — crystal structure

Define the prototype lattice that SQS structures are built on.
The averaged lattice option computes `a`, `b`, `c` from element-specific
reference values so the supercell scales with your element pool:

```toml
[phase]
key              = "HEDB1"
generator_name   = "HEDB"
a = 3.0624
b = 3.0624
c = 3.2392
alpha = 90
beta  = 90
gamma = 120
supercell_size       = [4, 3, 2]
use_average_lattice  = true       # average a/c over primary_elements
vectors = """1 0 0
0 1 0
0 0 1
"""
coords = """0.000000 0.000000 0.000000 a
0.333333 0.666667 0.500000 B
0.666667 0.333333 0.500000 B
"""
```

### Validate your configuration

Before a long run, check that all paths exist and Python dependencies are present:

```bash
pixi run python examples/full_framework.py my_run.toml --check
```

This exits 0 on success and prints every error found otherwise.

## Run the Full Pipeline

```bash
pixi run python examples/full_framework.py my_run.toml
```

To preview which stages will run without executing them:

```bash
pixi run python examples/full_framework.py my_run.toml --dry-run
```

The pipeline prints a stage header (`=== tdb ===`, `=== oxidation ===`, …)
before each enabled stage so you can track progress.

## Expected Outputs

### After `sqs_generation` + `tdb_fitting`

SQS trees and fitted TDB files are written under `Files/Comps/<system>/`:

```
Files/Comps/
├── CrHf/
│   └── HEDB1_2/
│       ├── sqs_lev=0/          — endmembers (pure Cr, pure Hf)
│       │   ├── CONTCAR
│       │   ├── energy
│       │   ├── force.out
│       │   ├── stress.out
│       │   └── str_relax.out
│       ├── sqs_lev=1/          — 50/50 alloy composition
│       │   └── ...
│       └── CrHf_*.tdb          — fitted CALPHAD database
├── CrMo/
│   └── ...
└── CrHfMo/
    └── ...
```

Each `.tdb` file is a standard CALPHAD database readable by pycalphad or Thermo-Calc.
Gibbs energy, Gibbs mixing, and binary/ternary phase-diagram plots are written
alongside the TDB when `generate_*` flags are enabled under `[tdb]`.

### After `phase_diagrams`

```
Files/Phase_Diagrams/
├── CrHfMo_200K.png
├── CrHfMo_400K.png
├── ...
└── CrHfMo.gif          — animated phase diagram (when make_gif = true)
```

### After `oxidation_graphs`

```
Files/oxidation/<system>/
├── mu_O_T_map.png           — thermodynamic stability map
├── onset_binary.csv         — oxidation-onset compositions (binary tie lines)
├── onset_ternary.csv        — oxidation-onset compositions (ternary)
├── scan_<T>K.png            — mu_O scan at fixed temperatures
└── composition_slices/      — per-composition mu_O–T maps
```

## Next Steps

- [SQS Generation](sqs-generation.md) — tuning `mcsqs` search time, cutoffs, and composition levels
- [TDB Fitting](tdb-fitting.md) — MLIP selection, relaxation parameters, and CALPHAD fit options
- [Oxidation Analysis](oxidation-analysis.md) — batch and single-composition oxidation workflows
