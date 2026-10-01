# TDB Fitting

`BladeTDBGen` drives the CALPHAD TDB generation step.  It copies SQS
directories from the BLADE staging tree into the MaterialsFramework `sqsdb`,
then invokes `Sqs2tdb.fit()` for each chemical system in `composition_list`.

## Constructor is side-effect-free

`BladeTDBGen.__init__` stores all configuration but does nothing else — no
files are written, no processes are spawned.  The full fitting loop runs only
when `.fit()` is called explicitly.  This makes it safe to construct a
generator early in a script, inspect its attributes, and defer execution.

```python
gen = BladeTDBGen(...)     # no side effects
gen.fit()                  # runs everything
```

## `tdb_params` keys

`tdb_params` is a dict forwarded to `Sqs2tdb` and `GraceCalculator`.  All
keys are optional; the defaults are shown below.

| Key | Default | Description |
|---|---|---|
| `"fmax"` | `0.005` | Force convergence criterion for MLIP relaxation in eV/Å |
| `"verbose"` | `True` | Print relaxation progress |
| `"calculator"` | `"grace"` | MLIP backend name passed to the MaterialsFramework calculator registry (e.g., `"grace"`, `"mace"`, `"orb"`, `"uma"`, `"chgnet"`).  Legacy device strings `"cuda"` and `"cpu"` are accepted for backwards compatibility and map to GRACE. |
| `"calculator_kwargs"` | `{"steps": 1000, "device": "cuda"}` | Keyword arguments forwarded to the calculator constructor |
| `"t_min"` | `298.15` | Lower temperature bound for CALPHAD polynomial fit in K |
| `"t_max"` | `10000.0` | Upper temperature bound in K |
| `"sro"` | `False` | Include short-range order correction in the CALPHAD model |
| `"bv"` | `5e-3` | Energy bump value applied to degenerate end-members |
| `"phonon"` | `False` | Include phonon contributions for end-member phases |
| `"open_calphad"` | `False` | Write OpenCalphad-compliant `.tdb` output |
| `"track_trajectory"` | `True` | Save `relaxation_live.xyz` trajectory during MLIP relaxation |
| `"terms"` | `None` | Additional interaction terms string appended to `terms.in` |

## CALPHAD model control: `terms_in`, `mult_in`, `sublattice_map`

### `terms_in`

Controls the interaction parameters (excess Gibbs energy terms) written to
`terms.in` between the two `sqs2tdb -fit` calls.  Keys are lattice base names
(e.g., `"HEDB1"`); values are the full `terms.in` file content as a string.
Phases absent from the dict use the ATAT-generated template unchanged.

The format follows the ATAT `terms.in` convention: `<order>,<index>` per
line, where the first integer is the interaction order and the second selects
the Redlich-Kister coefficient.

```python
terms_in = {
    "HEDB1": (
        "1,0:1,0\n"   # end-member terms
        "2,2:1,0\n"   # binary L2 interaction
    ),
}
```

### `mult_in`

Overrides the sublattice multiplicity file (`mult.in`) written between the two
`sqs2tdb -fit` calls.  Keys are lattice base names; values are file content.
Phases absent from the dict keep the ATAT-generated `mult.in` unchanged.

```python
mult_in = {"HEDB1": "a=1\tb=2\n"}
```

### `sublattice_map`

Assigns element symbols to sublattice letters.  This controls which elements
are placed on each sublattice when building the CALPHAD model, and is written
to `species.in` between the two `sqs2tdb -cp` calls.

Structure: outer key is the lattice base name; inner key is the sublattice
letter (or the special key `"Constant"` for sublattices whose composition is
fixed and should be excluded from the mixing cross-product); value is a list
of element symbols.

```python
sublattice_map = {
    "HEDB1": {
        "a": ["Cr", "Hf"],   # metal sublattice — mixing here
    },
}

# Multi-sublattice example with a fixed sublattice
sublattice_map = {
    "CARBIDE1": {
        "a": ["Cr", "Hf"],          # variable sublattice
        "b": ["C"],                  # effectively fixed
        "Constant": ["b"],           # exclude 'b' from composition cross-product
    },
}
```

## Output layout

For each composition in `composition_list`, a directory is created at:

```
BLADE/Files/Comps/<system>/
├── <system>.tdb                        # final CALPHAD TDB (e.g., CrHf.tdb)
└── <PHASE>_<n>/                        # one directory per phase
    └── sqs_lev=<L>_a_<El>=<x>.../     # one directory per SQS structure
        ├── CONTCAR                     # relaxed structure (VASP format)
        ├── energy                      # total energy in eV (single float)
        ├── str_relax.out               # relaxed structure in ATAT format
        ├── force.out                   # forces on all atoms
        └── stress.out                  # stress tensor
```

`<system>` is the concatenation of element symbols in composition order
(e.g., `CrHf` for `["Cr", "Hf"]`).  `skip_existing=True` skips a
composition if any `.tdb` file already exists in its output directory.
`refit_existing=True` combined with `skip_existing=True` rewrites
`terms.in`/`mult.in` and reruns only `sqs2tdb -fit` and `sqs2tdb -tdb`,
reusing the existing relaxed structures.

## Code example

```python
from pathlib import Path

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_tdb_gen import BladeTDBGen

path0 = Path("/path/to/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

phase_list = [
    {"generator_name": "HEDB", "lattice": "HEDB1", "supercell_size": (5, 5, 4)},
]

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

tdb_params = {
    "fmax": 1e-4,
    "verbose": True,
    "calculator": "grace",
    "calculator_kwargs": {"steps": 1000, "device": "cuda"},
    "t_min": 298.15,
    "t_max": 10000.0,
    "sro": False,
    "bv": 1e-3,
    "phonon": False,
    "open_calphad": False,
}

terms_in = {
    "HEDB1": "1,0:1,0\n2,2:1,0\n",
}

composer = BladeCompositions(
    primary_elements=["Cr", "Hf", "Ta"],
    secondary_elements=[],
    primary_min=2, primary_max=3,
    secondary_min=0, secondary_max=0,
)
composition_list = composer.generate_compositions()

gen = BladeTDBGen(
    phases=phase_list,
    phases_dict=phases,
    liquid=False,
    paths=paths,
    composition_list=composition_list,
    level=5,
    skip_existing=True,
    refit_existing=False,
    output_dir=path1 / "Files" / "Comps",
    terms_in=terms_in,
    mult_in=None,
    sublattice_map=None,
    tdb_params=tdb_params,
)
gen.fit()
```

## Per-composition overrides

When different chemical systems require different CALPHAD models, use
`system_overrides` at the script level (see
`examples/structures/tdb_gen_hedb.py`) to pass per-system `terms_in`,
`mult_in`, `sublattice_map`, and `fixed_compositions` dicts.  Each override
is merged with the global settings before constructing a separate `BladeTDBGen`
for that composition.
