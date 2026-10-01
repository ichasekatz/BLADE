# SQS Generation

`BladeSQS` handles the full special quasi-random structure (SQS) generation
sub-workflow for a single phase prototype.  It writes the ATAT input files,
runs `sqs2tdb -mk` to expand composition directories, then executes
`corrdump` and parallel `mcsqs` in each one.

## What SQS generation does

BLADE generates SQS structures using ATAT's `mcsqs` program.  The workflow
for each phase proceeds in four steps:

1. **Build input files.**  `BladeSQS` constructs `rndstr.skel` (the random
   structure skeleton) and `sqsgen.in` (the composition level specification)
   from the phase prototype dict and `sqsgen_levels`.

2. **Expand composition directories.**  `sqs2tdb -mk` reads `sqsgen.in` and
   creates one `sqsdb_lev=<level>_<sublattice>=<composition>` sub-directory
   per composition entry.

3. **Run `corrdump` then `mcsqs`.**  In each sub-directory, `corrdump`
   computes cluster correlations using the pair, triplet, and quadruplet
   cutoffs derived from the lattice parameters.  Then `params["parallel_runs"]`
   independent `mcsqs -ip=N` processes run simultaneously.  A `stopsqs`
   sentinel file is written after `params["time"]` seconds to stop them, and
   `mcsqs -best` merges the parallel results into a single `bestsqs.out`.

4. **Summarize.**  `objective_functions.txt` is written in the parent SQS
   directory with the final objective-function value from each sub-directory.

Cutoff distances are derived automatically from the lattice parameters by
`BladeCutoff` using nearest-neighbor shell distances.  `cutoff_mode = "nn"`
(the default) interprets the `"2"`, `"3"`, `"4"` params as NN shell indices
(decimals interpolate between shells).  Set `cutoff_mode = "angstrom"` to
pass distances in Angstroms directly.

## The `phases` dict structure

Each phase prototype is a plain dict with eight required keys:

| Key | Type | Description |
|---|---|---|
| `"a"`, `"b"`, `"c"` | `float` | Lattice parameter lengths in Angstroms |
| `"alpha"`, `"beta"`, `"gamma"` | `float` | Lattice angles in degrees |
| `"vectors"` | `str` | Lattice vector rows as a multi-line string |
| `"coords"` | `str` | Fractional coordinates, one site per line: `"x y z LABEL"` |

The `LABEL` field in `coords` is the key that distinguishes variable from
fixed sublattice sites.  A single lowercase letter (e.g., `a`, `b`) marks a
variable sublattice; an element symbol starting with an uppercase letter
(e.g., `B`, `Zr`) marks a fixed site that is not mixed by `mcsqs`.

```python
phases = {
    # AlB2-type diboride: metal sublattice 'a', fixed boron sites
    "HEDB1": {
        "a": 3.134,
        "b": 3.134,
        "c": 3.379,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.333333 0.666667 0.500000 B\n"
            "0.666667 0.333333 0.500000 B\n"
        ),
    },
    # FCC with two distinct variable sublattices
    "ALLOY2": {
        "a": 4.22,
        "b": 4.22,
        "c": 4.22,
        "alpha": 90,
        "beta": 90,
        "gamma": 90,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.000000 0.500000 0.500000 a\n"
            "0.500000 0.000000 0.500000 b\n"
            "0.500000 0.500000 0.000000 b\n"
        ),
    },
}
```

## Key parameters

### Constructor parameters

| Parameter | Type | Description |
|---|---|---|
| `phases_dict` | `dict` | Phase prototype dict (one entry from `phases`) |
| `sqsgen_levels` | `list[dict]` | Composition seed definitions for `sqsgen.in` |
| `level` | `int` | Highest sqsgen level to include (inclusive) |
| `len_comp` | `int` | Number of elements in the chemical system |
| `skip_existing_sqs` | `bool` | Skip directories that already have `bestcorr.out` |
| `sqsgen_in` | `str \| None` | Verbatim `sqsgen.in` content (bypasses auto-generation) |
| `sublattice_map` | `dict \| None` | Per-sublattice active element lists for multi-sublattice phases |

### `sqsgen_levels`

`sqsgen_levels` is a list of dicts that defines the composition space sampled
by `mcsqs`.  Each entry must contain:

- `"level"` (`int`): Level index.  Higher levels include more off-equimolar
  compositions and require larger supercells to represent faithfully.
- `"letter"` (`list[str]`): Sublattice letter(s) this entry applies to.
- `"compositions"` (`list[list[float]]`): Fractional compositions for that
  sublattice.  Each inner list must sum to 1.  The endmember `[1.0, 0.0]` at
  level 0 is always included.

BLADE auto-expands each non-endmember composition to all canonical
(sorted-descending) compositions sharing the same denominator.  For example,
`[0.75, 0.25]` also generates `[0.5, 0.5]` for a binary system.

```python
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]],                         "letter": ["a"]},
    {"level": 1, "compositions": [[0.5, 0.5]],                         "letter": ["a"]},
    {"level": 2, "compositions": [[0.75, 0.25]],                       "letter": ["a"]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]],        "letter": ["a"]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]],                  "letter": ["a"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]],     "letter": ["a"]},
]
```

### `mcsqs_params`

Passed to `BladeSQS.sqs_gen()` as the `params` argument.  The
`super_cell_size` key is merged in from the phase list entry.

| Key | Type | Description |
|---|---|---|
| `"time"` | `float` | Run duration per `sqsdb_lev=*` directory in seconds |
| `"parallel_runs"` | `int` | Number of simultaneous `mcsqs -ip=N` processes |
| `"cutoff_mode"` | `str` | `"nn"` (default) or `"angstrom"` |
| `"2"` | `float` | Pair cutoff: NN shell index or Angstrom distance |
| `"3"` | `float` | Triplet cutoff (0 to omit) |
| `"4"` | `float` | Quadruplet cutoff (0 to omit) |
| `"wr"`, `"wn"`, `"wd"` | `float` | `mcsqs` weight parameters |
| `"super_cell_size"` | `tuple[int,int,int]` | Supercell repetitions; total atoms = unit cell sites × product |

## Output layout

For each phase and system size, SQS output lands under:

```
BLADE/Files/SQS/<lattice>_<n>/
├── rndstr.skel                 # random structure skeleton (input to sqs2tdb -mk)
├── sqsgen.in                   # composition level specification
├── objective_functions.txt     # summary of objective values across all runs
└── sqsdb_lev=<L>_a=<comp>/     # one directory per composition entry
    ├── rndstr.in               # generated by sqs2tdb -mk
    ├── clusters.out            # cluster definitions from corrdump
    ├── bestcorr.out            # best correlation mismatch (objective function)
    ├── bestsqs.out             # best SQS structure (main output)
    └── objective_history.txt   # objective vs time log from monitoring thread
```

`bestsqs.out` is the SQS structure in ATAT format.  `BladeTDBGen` later
copies each `<lattice>_<n>/` tree into the MaterialsFramework `sqsdb` and
runs MLIP relaxation and CALPHAD fitting on it.

## Code example

```python
from pathlib import Path

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_sqsgen import BladeSQS

path0 = Path("/path/to/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

phases = {
    "HEDB1": {
        "a": 3.134, "b": 3.134, "c": 3.379,
        "alpha": 90, "beta": 90, "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.333333 0.666667 0.500000 B\n"
            "0.666667 0.333333 0.500000 B\n"
        ),
    },
}

phase_list = [
    {"generator_name": "HEDB", "lattice": "HEDB1", "supercell_size": (5, 5, 4)},
]

sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]],    "letter": ["a"]},
    {"level": 1, "compositions": [[0.5, 0.5]],    "letter": ["a"]},
    {"level": 2, "compositions": [[0.75, 0.25]],  "letter": ["a"]},
]

mcsqs_params = {
    "time": 30,
    "cutoff_mode": "nn",
    "2": 5,
    "3": 4,
    "4": 3,
    "wr": 20,
    "wn": 0.75,
    "wd": 1,
    "parallel_runs": 10,
}

composer = BladeCompositions(
    primary_elements=["Cr", "Hf"],
    secondary_elements=[],
    primary_min=2, primary_max=2,
    secondary_min=0, secondary_max=0,
)
composition_list = composer.generate_compositions()
unique_len_comps = composer.get_systems()

for specific_phase in phase_list:
    for len_comp in unique_len_comps:
        lattice = specific_phase["lattice"]
        sqs = BladeSQS(
            phases_dict=phases[lattice],
            sqsgen_levels=sqsgen_levels,
            level=2,
            len_comp=len_comp,
            skip_existing_sqs=True,
        )
        params = mcsqs_params | {"super_cell_size": specific_phase["supercell_size"]}
        sqs.sqs_gen(phase=specific_phase, paths=paths, params=params)
```

To inspect generated file text without running ATAT, call `sqs.sqs_struct()`
instead of `sqs.sqs_gen()`.  It returns `(sqsgen_text, rndstr_text)` and
prints both to stdout.
