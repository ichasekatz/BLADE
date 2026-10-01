# Oxidation Analysis

The oxidation framework computes oxygen chemical-potential stability windows across
temperature and composition space for multi-component boride (and related) systems.
It produces μ_O–T phase diagrams, phase-fraction maps, onset-composition curves,
3-D onset surfaces, composition slice plots, and animated phase-stability maps.

> **Credit**: The thermodynamic equilibrium solver and Gibbs free energy minimization
> were contributed by **Siya Zhu** ([siyazhu1@gmail.com](mailto:siyazhu1@gmail.com)).

---

## Overview

At each (T, μ_O) point the framework solves a grand-potential linear program (LP)
over a discrete simplex grid of competing phases — a flexible mixed phase whose
free energy is described by a Muggianu/RK model, plus fixed line-compound reference
phases from the Materials Project and OQMD/ORB databases.  The LP output gives
equilibrium phase amounts from which phase fractions and oxidation onset are derived.

The full analysis pipeline produces:

| Output | Description |
|---|---|
| **μ_O–x region maps** | Stability regions vs composition at fixed T |
| **μ_O–T region maps** | Stability regions vs temperature at fixed composition |
| **1-D scans** | Phase fractions vs μ_O at one composition and one or more T |
| **Onset diagrams** | Onset μ_O vs composition (binary: line plot; ternary: simplex map) |
| **3-D onset surfaces** | Onset μ_O as a function of both T and composition |
| **Composition slice maps** | μ_O–x or μ_O–T maps swept over fixed remainder ratios |
| **Animations** | Phase boundaries evolving with temperature or composition |

---

## Configuration: `OxidationConfig`

`OxidationConfig` is the single entry point for all path and chemistry settings.
All other classes read from it.

```python
from blade.oxidation import OxidationCalculator, OxidationConfig
from pathlib import Path

cfg = OxidationConfig(
    files_dir=Path("Research/BLADE/Files"),
    phase_element="B",                   # fixed element in the non-metal sublattice
    phase_element_stoichiometry=2.0,     # stoichiometry: M_x B_2 formula
    mixed_phase_subdir="blade",          # subfolder for BLADE/MLIP-relaxed phases
    fixed_phases_subdir="ORB",           # subfolder for MP + ORB reference phases
    region_label_mode="phases",          # "id" (numeric) or "phases" (formula)
    region_label_fontsize=7,
    slice_axis_priority=["Cr", "Hf", "Mo", "Nb", "Ta", "Ti", "V", "W", "Zr"],
)
```

### Key parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `files_dir` | `Path` | required | Root directory for all BLADE oxidation outputs |
| `phase_element` | `str \| None` | `None` | Non-metal sublattice element (e.g. `"B"` for borides) |
| `phase_element_stoichiometry` | `float` | `0.0` | Atoms per formula unit of `phase_element`; must be > 0 when `phase_element` is set |
| `mixed_phase_subdir` | `str` | `"blade"` | Sub-directory name for the BLADE-relaxed (flexible) phase calculations |
| `fixed_phases_subdir` | `str` | `"ORB"` | Sub-directory name for fixed line-compound reference phase calculations |
| `region_label_mode` | `str` | `"phases"` | How stability regions are labelled: `"id"` for numeric IDs, `"phases"` for phase formulae |
| `region_label_fontsize` | `int` | `7` | Font size (pt) for region labels in maps |
| `slice_axis` | `int \| str \| None` | `None` | Composition axis fixed when projecting high-dimensional diagrams to 2-D |
| `slice_axis_priority` | `list[str]` | `[]` | Ordered list of element candidates when `slice_axis` is `None` |

`OxidationConfig` exposes derived paths as properties:

```python
cfg.structures_dir   # files_dir / "system_structures"
cfg.tables_dir       # files_dir / "oxidation" / "tables"
cfg.outputs_dir      # files_dir / "oxidation" / "figures"
```

---

## Batch Mode: `run_batch()`

`OxidationCalculator.run_batch()` processes a list of named systems, running all
enabled analyses for each one.  When both `run_calculations` and `run_plots` are
enabled, the runner completes a full numerical pass for all systems before starting
a second, cache-only plotting pass.

```python
from blade.oxidation import OxidationCalculator, OxidationConfig
import numpy as np
from pathlib import Path

cfg = OxidationConfig(
    files_dir=Path("Research/BLADE/Files"),
    phase_element="B",
    phase_element_stoichiometry=2.0,
    mixed_phase_subdir="blade",
    fixed_phases_subdir="ORB",
    region_label_mode="phases",
    region_label_fontsize=7,
    slice_axis_priority=["Cr", "Hf", "Mo", "Nb", "Ta", "Ti", "V", "W", "Zr"],
)
calculator = OxidationCalculator(cfg)

calculator.run_batch(
    systems=["CrHfB", "CrMoB", "HfMoB"],
    temperature_values=list(range(250, 2001, 250)),        # K
    mu_O_values=[round(-10.0 + i * 0.1, 1) for i in range(61)],  # eV, -10 to -4
    run_calculations=True,
    run_plots=True,
    run_scan=True,
    run_muO_x_map=True,
    run_onset=True,
    run_3d_onset=False,
    run_composition_slices=True,
    run_composition_slice_muT=True,
    run_animations=True,
    skip_if_tables_exist=True,       # re-use cached tables when grid matches
    skip_if_analysis_exists=True,    # re-use cached analysis when grid matches
    onset_threshold=1e-8,            # absorbed-O fraction threshold for onset
    onset_comp_step_binary=0.02,
    onset_comp_step_ternary=0.025,
    onset_comp_step_nary=0.10,
    scan_x=0.5,
    scan_T=[700, 1273, 2000],
    scan_mu_O=[round(-10.0 + i * 0.05, 2) for i in range(121)],
    map_x_x_values=[round(i * 0.01, 2) for i in range(101)],
    slice_remainder_ratios=[0.0, 0.25, 0.5, 0.75, 1.0],
    slice_muT_comp_step=0.10,
)
```

### Key `run_batch()` parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `systems` | `list[str]` | `[]` (all) | Exact system directory names under `structures_dir` to process |
| `temperature_values` | array-like | `arange(200, 2001, 200)` | Temperature grid (K) for all maps |
| `mu_O_values` | array-like | `arange(-10, -4+ε, 0.10)` | μ_O grid (eV/O) for all maps |
| `elements` | `list[str]` | `[]` (all) | Filter: include only systems whose metals are a subset of this list |
| `run_calculations` | `bool` | `True` | Execute LP solves |
| `run_plots` | `bool` | `True` | Render figures |
| `run_scan` | `bool` | `True` | 1-D phase-fraction scan vs μ_O |
| `run_muO_x_map` | `bool` | `True` | μ_O–x region map at fixed T |
| `run_onset` | `bool` | `True` | Onset μ_O diagrams per temperature |
| `run_3d_onset` | `bool` | `False` | 3-D onset surface (T × composition) |
| `run_composition_slices` | `bool` | `True` | μ_O–x slice maps for fixed remainder ratios |
| `run_composition_slice_muT` | `bool` | `True` | Per-composition μ_O–T maps |
| `run_animations` | `bool` | `True` | GIF/MP4 animations from frame sequences |
| `skip_if_tables_exist` | `bool` | `True` | Re-use existing CSV tables when grid parameters match |
| `skip_if_analysis_exists` | `bool` | `True` | Re-use existing analysis figures |
| `onset_threshold` | `float` | `1e-8` | Minimum absorbed-O fraction to declare oxidation onset |

---

## Single-Composition Mode: `run_single_composition()`

`run_single_composition()` produces all three standard plots — a 1-D scan vs μ_O,
a full μ_O–T map, and a μ_O–x map — for exactly one alloy composition.  This is
useful for detailed investigation of a specific point in composition space.

```python
calculator.run_single_composition(
    system="CrHfMoB",
    metals=["Cr", "Hf", "Mo"],
    composition=[1/3, 1/3, 1/3],      # metal fractions; must sum to 1
    temperatures=list(range(250, 2001, 10)),
    mu_o=[-12.0 + i * 0.1 for i in range(81)],
    scan_temperatures=[700, 1000, 1273],
    x_values=[round(i * 0.01, 2) for i in range(101)],
    run_calculations=True,
    run_plots=True,
    skip_if_exists=True,
    rk_order=3,
    y_step=0.01,
    include_0p01_to_0p05_components=True,
)
```

For binary systems `composition` may be a scalar `x_M1`; the complementary metal
fraction is set automatically.  For ternary and higher systems pass a list whose
elements sum to 1.

The output directory for a single-composition run is:

```
oxidation/figures/<SystemName>/single_composition/<M1>x1_<M2>x2_<M3>x3/
```

---

## Outputs

All outputs land under `files_dir / "oxidation"`.

### Tables (`oxidation/tables/`)

Per-temperature CSV files are written for each analysis type:

| File pattern | Content |
|---|---|
| `<tag>_onset_auc_T<T>.csv` | Per-composition onset μ_O, AUC integrals, feasible μ_O range at temperature T |
| `<tag>_onset_auc_all_temperatures.csv` | Concatenated onset table across all temperatures |
| `<tag>_scan_T<T>.csv` | Phase fractions vs μ_O for the 1-D scan at temperature T |

Each onset CSV includes columns for the metal fractions (`x_Cr`, `x_Hf`, ...), `onset_muO_eV`,
`phase_fraction_auc_eV`, `parent_phase_fraction_auc_eV`, `oxide_phase_fraction_auc_eV`,
composition step, and μ_O grid bounds.

### Figures (`oxidation/figures/<SystemName>/`)

```
<SystemName>/
├── onset_auc/
│   ├── T250/
│   │   └── onset_diagram.png       # binary: onset μ_O vs x; ternary: simplex map
│   ├── T500/
│   │   └── onset_diagram.png
│   ├── ...
│   ├── onset_diagram.gif           # animated onset evolution with temperature
│   └── onset_diagram.mp4
├── muO_x_map/
│   └── T<T>/
│       └── assemblage_region_map.png
├── composition_slice_maps/
│   └── axis_<element>_<remainder>/
│       ├── muO_T.gif
│       ├── muO_T.mp4
│       └── muO_T_plots/
│           └── assemblage_region_map_x<val>.png
└── single_composition/
    └── <M1>x1_<M2>x2/
        ├── scan_T<T>.png
        ├── muO_T_map.png
        └── muO_x_map_T<T>.png
```

### 3-D onset surfaces

When `run_3d_onset=True` the framework renders 3-D onset surfaces showing onset μ_O
as a joint function of temperature and composition.  For ternary systems the surface
lives in a 4-D space (T × ternary simplex × μ_O); the figures project slices or
use colour to encode the fourth dimension.

---

## TOML Configuration

When running through `full_framework.py` all oxidation parameters are set in the
`[oxidation_batch]` table of `full_framework.toml`.  The element pool is inherited
from the global `[elements]` section.

```toml
[oxidation_batch]
run_calculations = true
run_plots        = true
run_onset        = true
run_3d_onset     = false
temperature_values = [250, 500, 750, 1000, 1250, 1500, 1750, 2000]
onset_threshold  = 1e-8
onset_comp_step_binary  = 0.02
onset_comp_step_ternary = 0.025
```
