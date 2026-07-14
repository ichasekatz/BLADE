# Additional Analysis and Visualization

All scripts here run directly with in-file settings. They are optional tools.

## Combined Phase Visualization

`phase_visualization.py` contains phase-map/GIF generation and synchronized
grid-video compilation. There are no separate plot or grid compiler files.

```python
run_phase_plots = True
run_phase_grid = True
```

- Plot only: first true, second false.
- Grid only from cached GIFs: first false, second true.
- Plot then compile: both true.
- Load/check without work: both false.

```bash
pixi run python examples/extra/phase_visualization.py
```

Plotting discovers top-level TDBs, calculates or reloads temperature-specific
free-energy data, writes static plots, and optionally writes GIFs. Grid mode
groups clean GIFs by system size and uses FFmpeg for synchronized MP4 files.
No GIFs is a normal, nonfatal grid result.

`full_framework.py` injects TOML values when it launches this script. Direct
execution does not require TOML.

## Other Utilities

- `compare_energy.py`: energy parity/residual figures and tables.
- `refit_tdb.py`: refit from existing structures/energies without MLIP
  relaxation.
- `tdb_gen_alloy.py`: complete single-sublattice alloy example.
- `tdb_gen_hedb_quaternary.py`: complete higher-order example.

Set each script's path/chemistry block and preserve its source inputs and
outputs so the operation can be independently retraced.
