# BLADE Examples: Reproducible Run Guide

Run every command from the BLADE repository root. Verify the absolute paths
near the top of a direct script before a production run.

## Environment

Choose one environment manager:

```bash
pixi install
pixi run python examples/full_framework.py --check
```

```bash
uv sync --locked --extra workflow
uv run --locked --extra workflow python examples/full_framework.py --check
```

External programs are not installed by `uv`: ATAT supplies `sqs2tdb`,
`corrdump`, and `mcsqs`; SCRAPS is needed only when selected; FFmpeg is
needed for MP4 output. Set `MP_API_KEY` before Materials Project downloads.

## Recommended Complete Run

1. Edit `examples/full_framework.toml`.
2. Set the five booleans in `[stages]`.
3. Set `paths.blade_root`, `paths.files_dir`, and `paths.sqsdb_dir`.
4. Define element pools, system-size limits, phase sites, and SQS levels.
5. Validate without calculations:

```bash
pixi run python examples/full_framework.py --check
pixi run python examples/full_framework.py --dry-run
```

6. Run the enabled workflow:

```bash
export MP_API_KEY="your-key"
pixi run python examples/full_framework.py
```

The positional TOML argument is optional. With no argument,
`examples/full_framework.toml` is used. To select another configuration:

```bash
pixi run python examples/full_framework.py examples/full_framework_max.toml
```

## Stage Order

1. `tdb`: compositions, optional SQS generation, MLIP relaxation, TDB
   fitting/refitting, and selected per-system plots.
2. `phase_visualization`: phase maps/GIFs and synchronized grid videos.
3. `database`: oxidation compositions, reference structures, MLIP energies,
   database tables, and collected system structures.
4. `oxidation`: batch or fixed-composition equilibrium calculations and plots.

Each stage writes durable artifacts under `Files/`. Skip/refit settings
control whether those artifacts are reused.

## Standalone Script Matrix

Every Python file under `examples/` can be launched directly. A direct
script uses in-file settings and does not read TOML, except
`full_framework.py`.

| Script | Prior artifacts or services |
|---|---|
| `structures/01_compositions.py` | None |
| `structures/02_sqs_generation.py` | ATAT and configured SQS database |
| `structures/03_tdb_generation.py` | Existing SQS directories and MLIP |
| `structures/04_visualization.py` | TDB/CONTCAR/trajectory inputs as applicable |
| `structures/05_data_analysis.py` | Structure and energy files |
| `structures/tdb_gen_*.py` | Backend-specific external tools |
| `oxidation/01_compositions.py` | Relevant BLADE `Files/Comps` data |
| `oxidation/02_poscars.py` | Step 01 and Materials Project access |
| `oxidation/03_energy.py` | Step 02 structures and MLIP |
| `oxidation/04_database.py` | Steps 01-03 |
| `oxidation/database_framework.py` | Services required by enabled substeps |
| `oxidation/oxidation_framework.py` | Completed system database |
| `extra/phase_visualization.py` | TDBs for plots; GIFs for grid-only mode |
| Other `extra/*.py` | Inputs documented in each script |

Use `pixi run python <script>` or
`uv run --locked --extra workflow python <script>`. “Runnable” means the
file is an executable example, not that downstream data can be synthesized
when prerequisites are absent. Files under `src/blade/` are importable
library modules, not CLI programs.

## Retracing a Run

1. Preserve the TOML or in-file settings.
2. Record the BLADE Git commit and environment lock file.
3. Keep `terms.in`, `mult.in`, `sqsgen.in`, SQS correlation files,
   relaxed structures, energies, and generated TDBs.
4. Keep oxidation `tables/` and `.npz` caches when plots may need to be
   regenerated without equilibrium solves.
5. Use skip flags only after confirming the expected artifact exists.
6. Use `refit_existing_tdb = true` when model terms changed but energies did
   not.
7. Run `--check` and `--dry-run` before restarting a partial workflow.

See [structures/README.md](structures/README.md),
[oxidation/README.md](oxidation/README.md), and
[extra/README.md](extra/README.md). The complete reference is
`docs/BLADE_framework_manual.tex`.
