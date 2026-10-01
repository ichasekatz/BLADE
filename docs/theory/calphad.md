# CALPHAD Method

CALPHAD (CALculation of PHAse Diagrams) is a semi-empirical framework for computing
multi-component thermodynamic phase equilibria from assessed Gibbs energy functions.
Rather than relying on direct enumeration of phases, CALPHAD parameterizes the Gibbs
energy of each phase as a function of temperature and composition, then minimizes the
total Gibbs energy of the system to predict stable phase assemblages.

The assessed parameters for each phase are stored in a **thermodynamic database (TDB)**.
A TDB contains Gibbs energy expressions for all relevant phases — stoichiometric
compounds, solution phases, and the liquid — in a format compatible with CALPHAD
solvers such as OpenCalphad or Thermo-Calc.

BLADE automates the most labor-intensive part of TDB construction: fitting the
composition-dependent interaction parameters for each binary and ternary sub-system
from first-principles energetics, without manual assessment.

## Gibbs Energy Models

### Stoichiometric endmembers

For a pure element or stoichiometric compound, the molar Gibbs energy is expressed as a
polynomial in temperature:

```
G(T) = a + bT + cT·ln(T) + Σ dₙTⁿ
```

The coefficients `a`, `b`, `c`, and `dₙ` are fit to heat capacity data or to
MLIP-computed 0 K enthalpies supplemented by empirical corrections.

### Subregular solution model for excess Gibbs energy

For a substitutional alloy phase — such as the metal sublattice of an AlB₂-type
diboride — the molar Gibbs energy is written as:

```
G = Σᵢ xᵢ Gᵢ° + RT Σᵢ xᵢ ln(xᵢ) + G_xs
```

where `Gᵢ°` is the endmember Gibbs energy, the second term is the ideal mixing
entropy, and `G_xs` is the excess Gibbs energy.

BLADE models `G_xs` with **Redlich–Kister (RK) polynomials**, which form the basis
of the subregular solution model:

```
G_xs = Σᵢ<ⱼ xᵢ xⱼ Σₖ Lᵢⱼ,ₖ (xᵢ − xⱼ)ᵏ
```

where `Lᵢⱼ,ₖ` are the interaction parameters. The order of the series controls the
asymmetry of the excess energy:

- **k = 0** (regular solution): symmetric G_xs, one parameter per binary pair.
- **k = 1** (subregular solution): first-order asymmetry, two parameters per pair.
- **k ≥ 2**: higher-order corrections; rarely needed for diborides.

Ternary interactions add a cross term `xᵢ xⱼ xₖ Lᵢⱼₖ`.

BLADE fits the `Lᵢⱼ,ₖ` coefficients from MLIP-relaxed SQS total energies using
ATAT's `sqs2tdb -fit`.

## How BLADE Generates a TDB

The full pipeline from structure to database proceeds in four stages:

### 1. SQS generation

Special Quasirandom Structures (SQS) are generated for each binary and ternary
composition needed to constrain the RK expansion. Each SQS is a finite supercell
that reproduces the pair-correlation functions of a random alloy at a target
composition (e.g. x = 0.5 for a binary 50:50 mixture). BLADE uses ATAT's `mcsqs`
to search for the best SQS at each composition; see [SQS Method](sqs.md) for
details.

### 2. MLIP relaxation

Each SQS supercell is relaxed with a machine-learned interatomic potential (MLIP)
— currently GRACE-2L or ORB, configurable via `tdb.mlip`. The `GraceCalculator`
(in `MaterialsFramework/calculators/grace.py`) wraps the MLIP and writes the
relaxed geometry and total energy to `CONTCAR`, `energy`, `force.out`, and
`stress.out` inside the corresponding sqsdb sub-directory.

### 3. Energy extraction

The `sqs2tdb -cp` command copies the relaxed structures into the ATAT sqsdb tree.
The mixing enthalpy for each SQS is then:

```
ΔH_mix(x) = E_SQS(x) − Σᵢ xᵢ Eᵢ°
```

where `Eᵢ°` is the MLIP energy of element `i` in its reference structure (or a
database fallback if element data are unavailable).

### 4. Interaction parameter fit

`sqs2tdb -fit` assembles the Redlich–Kister matrix from all available SQS
enthalpies across every binary and ternary sub-system and solves for the
interaction parameters by least squares. The resulting TDB is written by
`sqs2tdb -tdb`. In BLADE, the entire fit–write sequence is driven by
`MaterialsFramework.tools.sqs2tdb.Sqs2tdb.fit()`.

## AlB₂-Type Diboride Prototypes

BLADE targets transition-metal diborides with the **AlB₂ structure** (space group
P6/mmm, Strukturbericht B35). The unit cell contains three atoms per formula unit:

- **1 metal atom** at the origin (Wyckoff 1a: 0, 0, 0)
- **2 boron atoms** at the honeycomb sites (Wyckoff 2d: 1/3, 2/3, 1/2 and 2/3, 1/3, 1/2)

The hexagonal lattice is characterized by two free parameters: the in-plane
nearest-neighbor metal–metal distance `a` and the interlayer spacing `c`. For the
default HEDB1 prototype in BLADE:

```
a = b = 3.0624 Å      (in-plane, averaged over Cr, Hf, Mo)
c = 3.2392 Å           (out-of-plane)
alpha = beta = 90°,  gamma = 120°
```

These are composition-averaged lattice parameters computed from pure-element
endmember relaxations. When `use_average_lattice = true` in the TOML config, BLADE
averages `a` and `c` across all primary elements before generating SQS structures.

The two-sublattice model in ATAT notation is `(M)₁(B)₂`, where the metal sublattice
`M` is substitutionally disordered and the boron sublattice is fixed. This is encoded
in `rndstr.skel`:

```
# hexagonal vectors
3.0624  0.0     0.0
-1.5312 2.6519  0.0
0.0     0.0     3.2392
# metal site — variable (marked with species letter 'a')
0.000  0.000  0.000   a
# boron sites — fixed
0.333  0.667  0.500   B
0.667  0.333  0.500   B
```

Higher-order phases (ternary, quaternary) use the same formalism. The RK model
simply gains additional pairwise and cross-interaction terms.

## Reference State

Energies are referenced to MLIP-relaxed pure elements in their standard structure.
If element relaxation data are unavailable, BLADE falls back to hardcoded values in
`_database_fallback_refs` rather than requiring additional TOML configuration.

## MaterialsFramework Integration

The fitting engine is `MaterialsFramework.tools.sqs2tdb.Sqs2tdb`, located in
`PhaseForge/MaterialsFramework/tools/sqs2tdb.py`. This class wraps the four ATAT
commands (`sqs2tdb -mk`, `-cp`, `-fit`, `-tdb`) and the MLIP relaxation step into a
single `.fit()` call. BLADE's `BladeTDBGen` invokes `Sqs2tdb.fit()` for each
composition, passing the `tdb_params` dict to control device, convergence tolerance
(`fmax`), temperature range, and RK fit options (`sro`, `bv`, `phonon`).
