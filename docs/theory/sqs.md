# SQS Method

## What Are Special Quasirandom Structures?

A **Special Quasirandom Structure (SQS)** is a finite periodic supercell designed to
mimic the site-correlation functions of a perfectly disordered (random) alloy at a
given composition. The concept was introduced by Zunger et al. (1990) as a practical
way to compute alloy properties with periodic DFT or MLIP codes, which require a
repeat unit, while still capturing the statistical character of a solid solution.

For a binary A₁₋ₓBₓ random alloy, the key quantities to match are the
**Warren–Cowley pair-correlation functions**:

```
alpha_m = (P_AA(m) − x_A²) / (x_A(1 − x_A))
```

where `P_AA(m)` is the probability that two sites separated by shell `m` are both
occupied by A. In a perfectly random alloy, `alpha_m = 0` for all shells `m >= 1`.
An SQS minimizes `|alpha_m|` for the most important shells — typically the first few
nearest neighbors — within the constraint of a finite supercell.

The total energy of the SQS, referenced to the pure endmembers, gives the **mixing
enthalpy** at composition x:

```
dH_mix(x) = E_SQS(x) − sum_i x_i * E_i°
```

This is the primary input to the CALPHAD Redlich–Kister fit. See
[CALPHAD Method](calphad.md) for how these enthalpies map to interaction parameters.

## ATAT mcsqs Algorithm

BLADE uses the **ATAT** (Alloy Theoretic Automated Toolkit) implementation of SQS
generation, specifically the `mcsqs` program.

`mcsqs` performs a **Monte Carlo search** over the space of supercell shapes and
atomic decorations:

1. A parent lattice and target composition are specified in `rndstr.in` and
   `sqsgen.in`, respectively.
2. An objective function — a weighted sum of squared correlation differences from
   the ideal random alloy — is minimized by randomly swapping atom species on the
   lattice sites.
3. The search runs until either the objective function reaches zero (a perfect SQS
   is found) or a timeout (`stopsqs` sentinel file) is detected.
4. The best structure found is written to `bestsqs.out` and the final objective
   value to `bestcorr.out`.

Multiple independent `mcsqs` processes can run in parallel (each started with a
different random seed via `-ip=N`) and race to converge. BLADE defaults to
`parallel_runs = 10` concurrent instances per composition directory. The winner's
`bestsqs.out` is used; ties resolve to whichever process wrote last.

Pair, triplet, and quadruplet correlation cutoffs are computed automatically from
the lattice parameters by `BladeCutoff`, based on the `cutoff_mode` setting
(`"nn"` for nearest-neighbor multiples or explicit distance values).

## SQS Level Scheme in BLADE

BLADE organizes compositions into a **level hierarchy**. Each level adds finer
composition points, progressively improving the Redlich–Kister fit while keeping
early levels cheap. The level used for a given run is set by `tdb.level` in the
TOML configuration.

| Level | Compositions (binary form)       | Purpose                                      |
|-------|----------------------------------|----------------------------------------------|
| 0     | [1.0, 0.0]                       | Pure endmembers — required by all fits       |
| 1     | [0.5, 0.5]                       | 50:50 binary, constrains L0                  |
| 2     | [0.75, 0.25]                     | Off-center binary, adds L1 resolution        |
| 3     | [1/3, 1/3, 1/3]                  | Equimolar ternary                            |
| 4     | [0.5, 0.25, 0.25]                | Ternary with dominant species                |
| 5     | [0.875, 0.125], [0.625, 0.375]   | Dilute/asymmetric binary corrections         |
| 6     | [0.75, 0.125, 0.125]             | Dilute ternary, highest-order correction     |

Level 0 compositions are pure endmembers — they do not require `mcsqs` runs because
the endmember energy is taken directly from the MLIP relaxation of the ordered
structure. Levels 1 and 2 constrain binary interaction parameters. Levels 3 and 4
introduce ternary terms. Levels 5 and 6 refine the fit in composition regions that
are dilute or asymmetric.

Each composition entry is automatically expanded to all canonical (descending-sorted)
compositions sharing the same denominator. For example, specifying `[0.75, 0.25]` for
a ternary system also generates `[0.75, 0.25, 0.0]`, `[0.5, 0.5, 0.0]`, etc., so
every binary sub-system contained in the ternary is covered.

The sqsdb directory for each composition point is named:

```
sqsdb_lev=<level>_<composition string>/
```

for example `sqsdb_lev=1_0.5_0.5/`.

## Supercell Size: Accuracy vs. Computation Time

The supercell size sets the number of atoms available for the `mcsqs` search and
directly controls the quality of the SQS.

### Why size matters

A larger supercell offers more degrees of freedom to zero out correlation functions
across many neighbor shells. Small supercells (fewer than ~8 atoms) can match only
the first one or two shells and may retain non-negligible correlations at larger
distances, introducing a systematic error in the mixing enthalpy.

### Practical constraints

Each additional atom multiplies the DFT or MLIP relaxation cost roughly linearly
(for MLIP) to cubically (for DFT). For the HEDB1 AlB2 prototype with 3 atoms per
formula unit, the default supercell is:

```
supercell_size = [4, 3, 2]     # 4 x 3 x 2 = 24 unit cells = 72 atoms total
```

This gives 24 metal sites — enough to represent compositions at the 1/24 level (~4%
steps) — while keeping MLIP relaxation times under a few minutes per structure on CPU.

### Recommended trade-off

| Priority            | Supercell atoms | Notes                                          |
|---------------------|-----------------|------------------------------------------------|
| Fast screening      | 12–24           | Binary L0 only; coarse accuracy                |
| Standard production | 48–72           | Covers L0–L1; BLADE default for HEDB1          |
| High accuracy       | 96–144          | Diminishing returns beyond ~120 for diborides  |

For MLIP-based workflows the computational overhead of going from 72 to 96 atoms is
modest; for DFT it can be prohibitive. The BLADE default of 72 atoms reflects a
balance appropriate for GRACE-2L or ORB relaxations used in the HEDB diboride work.
