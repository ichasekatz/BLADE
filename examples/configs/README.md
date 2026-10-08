# configs/ — Use-case configuration overrides

Each TOML here is a thin override for a specific crystal system. Pass it as a
second argument to `full_framework.py` to merge its settings on top of the base
config:

```bash
uv run --extra workflow python examples/full_framework.py \
    examples/full_framework.toml examples/configs/<system>.toml
```

Keys in the use-case file override matching keys in the base; everything else
inherits from the base.

## Available configs

| File | System | Prototype | Sublattices |
|---|---|---|---|
| `borides_hedb.toml` | AlB₂-type hexagonal borides | HEDB1 | a (metals), B fixed |
| `max_phases.toml` | M₂AX MAX phases | MAX1 | a (M-metals), b (A-metal), C fixed |
| `alloy_bcc.toml` | BCC refractory HEAs | BCC1 | a (metals only) |

## What to set in your base toml

- `[paths] sqsdb_dir` — path to your ATAT sqsdb directory
- `[paths] blade_root` — auto-detected from script location; override if needed
- `[database] api_key_env` — env var holding your Materials Project API key

## Adding a new crystal system

Copy the closest existing config and change:
1. `[phase]` — geometry (`a`, `b`, `c`, angles, `coords`, `supercell_size`)
2. `[elements]` — element pool and composition range
3. `tdb_sqs_levels` — `letter` list matches the variable sites in `coords`
4. `[tdb_inputs.terms_in]` — CE interaction terms per sublattice count
5. `[phase_plots]` — `fixed_elements` and `metal_fraction` for the system
