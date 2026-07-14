# Oxidation Examples

The oxidation workflow is chemistry-general. `fixed_elements`,
`oxygen_element`, and `phase_element` are explicit; boron appears only in
demonstration settings.

## Database Construction

1. Generate systems and scan BLADE data:

   ```bash
   pixi run python examples/oxidation/01_compositions.py
   ```

   This writes the composition workbook and scans relevant BLADE artifacts.

2. Download reference structures:

   ```bash
   export MP_API_KEY="your-key"
   pixi run python examples/oxidation/02_poscars.py
   ```

   This queries Materials Project and writes POSCARs plus a summary workbook.

3. Calculate consistent MLIP energies:

   ```bash
   pixi run python examples/oxidation/03_energy.py
   ```

   Existing energy artifacts can be reused with the script's skip setting.

4. Build and collect the database:

   ```bash
   pixi run python examples/oxidation/04_database.py
   ```

   This combines BLADE/reference energies, creates formation-energy tables,
   and collects structures and TDBs into per-system directories.

## Resumable Database Driver

```bash
export MP_API_KEY="your-key"
pixi run python examples/oxidation/database_framework.py
```

`PipelineConfig` has four independent `run_compositions`,
`run_download`, `run_energies`, and `run_database` switches. Direct
execution uses these in-file settings and does not need TOML.

## Oxidation Calculations and Plots

```bash
pixi run python examples/oxidation/oxidation_framework.py
```

Use `MODE = "batch"` for multiple discovered systems or `MODE = "single"`
for one fixed initial composition. `OxidationConfig` controls paths, fixed
phase species, stoichiometry, labels, and slice-axis priority. Batch and
single arrays define temperature, oxygen chemical potential, scan, and
composition grids. Set `BATCH_SETTINGS["elements"]` to the allowed metal
pool; an empty list includes all systems. Batch muO-T maps are generated per
composition under each slice's `muO_T_plots/`, while GIF/MP4 products remain
at the slice root.

`BATCH_SETTINGS["run_calculations"] = True` and
`BATCH_SETTINGS["run_plots"] = True` use two global passes: all selected
systems are calculated and cached before plotting starts. Disable
`run_plots` for a calculation-only run. Disable `run_calculations` for a
strict cache-only redraw; missing caches are reported instead of being solved.

## Cache and Provenance Rules

- `skip_if_tables_exist` reuses equilibrium tables.
- `skip_if_analysis_exists` reuses derived analysis.
- `run_calculations` and `run_plots` control the two oxidation execution stages.
- HTML hover text uses exact point data; static region keys summarize ranges.
- Keep per-composition `tables/` and `.npz` caches to redraw plots without
  equilibrium solves.
- Grid, database, or solver-setting changes invalidate reuse assumptions.

Retain source TDBs, energy references, phase structures/energies, fixed-species
settings, initial compositions, grids, label settings, package versions, and
the exact script/TOML.
