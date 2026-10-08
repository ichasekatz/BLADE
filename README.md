<div align="center">

# BLADE

[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](https://spdx.org/licenses/GPL-3.0-or-later.html)
![Python](https://img.shields.io/badge/python-3.12-blue)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macos-lightgrey)

[![Tests](https://github.com/ichasekatz/BLADE/actions/workflows/tests.yml/badge.svg)](https://github.com/ichasekatz/BLADE/actions/workflows/tests.yml)
[![Lint](https://github.com/ichasekatz/BLADE/actions/workflows/lint.yml/badge.svg)](https://github.com/ichasekatz/BLADE/actions/workflows/lint.yml)

<!-- DOI badge will be added once the Zenodo record is minted -->

**High-throughput Python framework for automated CALPHAD TDB generation, SQS structure creation, MLIP-based structural relaxation, and oxidation-stability analysis for multi-component boride and alloy systems.**

<p>
  <a href="https://github.com/ichasekatz/BLADE/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement">Request a Feature</a> |
  <a href="https://ichasekatz.github.io/BLADE">Documentation</a>
</p>

</div>

---

## Key Features

- Enumerate all valid N-component element combinations from primary and secondary element pools under user-defined size constraints; the composition list feeds every downstream stage without repetition
- Write ATAT input files, run `sqs2tdb -mk` to populate `sqsdb_lev=*` sub-directories, and spawn parallel `mcsqs` instances; cutoff distances are derived automatically from lattice parameters
- Orchestrate [MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework)'s `Sqs2tdb` fitting workflow across many chemical systems; the constructor is side-effect-free and fitting runs only when `.fit()` is called
- Compute oxygen chemical-potential stability windows across temperature and composition space; supports batch mode over many systems, onset-composition mapping, μ_O/T scans, 3-D onset surfaces, composition slice plots, and animated phase-stability maps
- Drive every stage from a single TOML file; the `[elements]` section is the single source of truth for the element pool, and stages are toggled individually so interrupted runs resume without repeating completed work

---

## Pipeline

| Class | Stage | Responsibility |
|---|---|---|
| `BladeCompositions` | Enumeration | Generate all valid element subsets under `primary_min`/`primary_max` and `secondary_min`/`secondary_max` constraints |
| `BladeSQS` | SQS generation | Write ATAT input files, invoke `sqs2tdb -mk` and `corrdump`, run parallel `mcsqs` searches, and stage SQS trees for fitting |
| `BladeTDBGen` | TDB fitting | Copy SQS trees into the sqsdb, call `Sqs2tdb.fit()` per composition via MLIP relaxation, write CALPHAD TDB files and diagnostic plots |
| `OxidationCalculator` | Oxidation analysis | Minimize grand-potential LP at each (T, μ_O) point; produce phase-fraction maps, onset curves, and 3-D stability surfaces |

---

## Installation

BLADE is managed with [uv](https://docs.astral.sh/uv/). Clone the repository and sync the environment — all dependencies, including [MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework), are installed automatically.

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE
uv sync --extra workflow
```

### ATAT binaries

The SQS generation and TDB fitting stages require the ATAT suite (`mcsqs`, `sqs2tdb`, `corrdump`) on `$PATH`. These are not Python packages; obtain them from the [ATAT website](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/) and add the binary directory to your shell profile.

```bash
# confirm ATAT is available before running
which mcsqs sqs2tdb corrdump
```

---

## Quickstart

### Running the full pipeline

Configure your element pool and paths once in `examples/full_framework.toml`, then run:

```bash
uv run --extra workflow python examples/full_framework.py examples/full_framework.toml
```

Enable or disable stages with boolean flags — set a stage to `false` to skip it and resume from where a previous run left off:

```toml
[stages]
sqs_generation   = true
tdb_fitting      = true
phase_diagrams   = false
oxide_database   = false
oxidation_graphs = false
```

### Enumerate compositions

```python
from blade.tools.blade_compositions import BladeCompositions

composer = BladeCompositions(
    primary_elements=["Cr", "Hf", "Mo"],
    secondary_elements=[],
    primary_min=2,
    primary_max=3,
    secondary_min=0,
    secondary_max=0,
)
compositions = composer.generate_compositions()
# [['Cr', 'Hf'], ['Cr', 'Mo'], ['Hf', 'Mo'], ['Cr', 'Hf', 'Mo']]
```

### Generate SQS structures

```python
from blade.tools.blade_sqsgen import BladeSQS

sqs = BladeSQS(
    phases_dict=phases["HEDB1"],
    sqsgen_levels=sqsgen_levels,
    level=6,
    len_comp=2,
    skip_existing_sqs=True,
)
sqs.sqs_gen(phase=phase_list[0], paths=paths, params=mcsqs_params)
```

### Fit CALPHAD TDBs

```python
from blade.tools.blade_tdb_gen import BladeTDBGen

gen = BladeTDBGen(
    phases=phase_list,
    phases_dict=phases,
    liquid=False,
    paths=paths,
    composition_list=compositions,
    level=6,
    skip_existing=True,
    refit_existing=False,
    output_dir=comps_dir,
    tdb_params={"calculator": "orb", "fmax": 1e-4, "t_min": 298.15, "t_max": 10000.0},
)
gen.fit()
```

### Oxidation analysis

```python
from blade.oxidation import OxidationCalculator, OxidationConfig

calculator = OxidationCalculator(OxidationConfig(
    files_dir=blade_root / "Files",
    phase_element="B",
    phase_element_stoichiometry=2.0,
    mixed_phase_subdir="blade",
    fixed_phases_subdir="ORB",
))
calculator.run_batch(
    systems=["CrHfB", "CrMoB"],
    temperature_values=list(range(250, 2001, 250)),
    mu_O_values=[round(-10.0 + i * 0.1, 1) for i in range(61)],
    run_calculations=True,
    run_plots=True,
    run_onset=True,
)
```

See the numbered scripts in [`examples/`](examples/) for complete, runnable demonstrations of each class.

---

## Configuration

All BLADE behavior is controlled by a single TOML file. The annotated reference is [`examples/full_framework.toml`](examples/full_framework.toml).

### Elements — defined once, inherited everywhere

```toml
[elements]
primary       = ["Cr", "Hf", "Mo", "Ta", "Zr"]
secondary     = []
primary_min   = 2
primary_max   = 3
secondary_min = 0
secondary_max = 0
```

Every stage (`[tdb]`, `[database]`, `[oxidation_batch]`) inherits this pool automatically. Override per-stage by adding `primary = [...]` inside the relevant section.

### Paths — only one path required

```toml
[paths]
sqsdb_dir = "/path/to/PhaseForge/atat/data/sqsdb"
```

`blade_root` and `files_dir` are auto-detected from the location of the entry-point script. No other paths need to be set.

### Plot generation — one flag for all plots

```toml
[tdb]
generate_plots = true
```

This fans out to all five plot types (Gibbs energy, Gibbs mixing, phase diagram, combined phase diagram, CONTCAR visualizations). Set to `false` to skip all plots.

---

## Examples

| Script | Demonstrates |
|---|---|
| [`01_compositions.py`](examples/01_compositions.py) | `BladeCompositions`: binary, ternary, mixed pools |
| [`02_sqs_generation.py`](examples/02_sqs_generation.py) | `BladeSQS`: phase prototype, composition levels, mcsqs params |
| [`03_tdb_fitting.py`](examples/03_tdb_fitting.py) | `BladeTDBGen`: tdb_params, terms_in, `.fit()` |
| [`04_oxidation_analysis.py`](examples/04_oxidation_analysis.py) | `OxidationCalculator`: batch mode, single composition |
| [`05_full_pipeline.py`](examples/05_full_pipeline.py) | End-to-end TOML-driven pipeline, all stages |

---

## Running Scripts

Any BLADE script runs with `uv run` — no manual dependency management needed:

```bash
uv run python examples/01_compositions.py

# background job with logging
nohup uv run python -u examples/05_full_pipeline.py > run.log 2>&1 &
```

---

## Tests

Unit tests run offline (no MLIP or ATAT required). Integration tests are gated by `@pytest.mark.integration` and require ATAT binaries and a compatible MLIP on the machine.

```bash
# unit tests only
uv run python -m pytest tests/ -m "not integration"

# all tests (requires ATAT + MLIP)
uv run python -m pytest tests/
```

---

## Credits

- **Chase Katz** — BLADE pipeline, SQS generation, TDB fitting automation, and oxidation analysis framework
- **Siya Zhu** (siyazhu1@gmail.com) — Gibbs free energy minimization and thermodynamic equilibrium solver (`blade/oxidation/framework/`)
- **Doguhan Sariturk** — [MaterialsFramework](https://github.com/dogusariturk/MaterialsFramework) and the underlying CALPHAD tooling (`Sqs2tdb`, `GraceCalculator`) that BLADE builds on

---

## Citing

If you use BLADE in your research, please cite:

> Chase Katz. *BLADE: Boride Lattice Automated Discovery Engine* (Version 1.7.0). Texas A&M University, 2026. https://github.com/ichasekatz/BLADE

```bibtex
@software{katz2026blade,
  author      = {Katz, Chase and Zhu, Siya},
  title       = {{BLADE}: {Boride Lattice Automated Discovery Engine}},
  year        = {2026},
  version     = {1.7.0},
  institution = {Texas A\&M University},
  url         = {https://github.com/ichasekatz/BLADE},
  license     = {GPL-3.0-or-later},
}
```

---

## License

BLADE is distributed under the [GNU General Public License v3.0 or later](LICENSE).
