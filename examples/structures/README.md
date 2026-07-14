# Structure, SQS, and TDB Examples

These scripts demonstrate a general structure-to-CALPHAD workflow. HEDB and
MAX are example prototype names; the BLADE classes accept user-defined phases,
sites, species, and element counts.

## Numbered Workflow

### 1. Enumerate systems

```bash
pixi run python examples/structures/01_compositions.py
```

Set primary/secondary element pools and count limits.
`BladeCompositions.generate_compositions()` returns element lists and
`get_systems()` reports generated system sizes. No external binary or prior
data is required.

### 2. Generate SQS inputs and structures

```bash
pixi run python examples/structures/02_sqs_generation.py
```

Set paths, element constraints, lattice vectors/sites, supercell dimensions,
SQS levels, correlation cutoffs, and `skip_existing_sqs`. The script writes
ATAT inputs and invokes `sqs2tdb`, `corrdump`, and parallel `mcsqs`.
Retain `rndstr.skel`, `sqsgen.in`, `bestcorr.out`, and structures.

### 3. Fit TDBs

```bash
pixi run python examples/structures/03_tdb_generation.py
```

Use the same phase/site definitions and composition levels as step 2. Select
the MLIP, fitting tolerance/range, model terms, site maps, and output folder.
`BladeTDBGen.fit()` performs the work; class construction is side-effect
free. Outputs include relaxed structures, energies, fit inputs, and TDB files.

### 4. Visualize

```bash
pixi run python examples/structures/04_visualization.py
```

This demonstrates combined CONTCAR images, Gibbs and mixing-energy curves,
binary phase diagrams, montages, and relaxation movies. FFmpeg is required for
MP4 output.

### 5. Extract data

```bash
pixi run python examples/structures/05_data_analysis.py
```

`BLADEData` scans structural/energy artifacts into a DataFrame for filtering
and export.

## Complete Prototype Drivers

- `tdb_gen_hedb.py`: one-variable-sublattice ATAT/mcsqs example.
- `tdb_gen_hedb_scraps.py`: corresponding SCRAPS example.
- `tdb_gen_max.py`: multibasis ATAT/mcsqs example.
- `tdb_gen_max_scraps.py`: multibasis SCRAPS example.

Each runs directly using top-of-file settings. `full_framework.py` also uses
these as drivers and injects TOML settings without editing their source.

## Restart and Provenance

- `skip_existing_sqs` reuses SQS results when the expected result exists.
- `skip_existing_tdb` reuses a generated TDB.
- `refit_existing_tdb` reruns fitting from existing energies so changed
  `terms.in` or `mult.in` can be applied without MLIP recalculation.
- `skip_existing_plots` checks each enabled output independently.
- Individual `generate_*` booleans control each plot/CONTCAR product.
- Disable trajectory tracking when movies are not required.

Record element ordering, phase coordinates, lattice parameters, supercell,
SQS backend/version and cutoffs, selected objective, MLIP/version, relaxation
tolerance, model terms, fit range, and fixed-site compositions. A TDB alone is
not complete provenance.
