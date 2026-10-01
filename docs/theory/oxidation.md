# Oxidation Thermodynamics

The equilibrium solver and Gibbs free energy minimization described on this page
were contributed by **Siya Zhu** ([siyazhu1@gmail.com](mailto:siyazhu1@gmail.com)).

---

## Grand-Potential Thermodynamic Model

The central quantity is the grand potential Ω, which replaces the Gibbs free
energy when the oxygen chemical potential μ_O is imposed by the environment
rather than by mass balance.  At each (T, μ_O) point the framework finds the
mixture of phases that minimizes the total grand potential subject to elemental
conservation.

### Phase free energy

Each phase candidate — whether a flexible SQS-derived mixed phase or a fixed
line-compound reference phase — contributes a composition-dependent free energy.

For the flexible mixed phase, whose internal composition vector **y** lives on
the N-metal simplex, the molar Gibbs free energy is:

```
G(y, T) = H_Muggianu(y) + k_B T Σ_i y_i ln(y_i)
```

where:

- **H_Muggianu(y)** is the Muggianu mixing enthalpy (see below).
- The second term is the ideal configurational entropy.  k_B is the Boltzmann
  constant (8.617 × 10⁻⁵ eV K⁻¹).
- The composition coordinates y_i are the metal-sublattice fractions; they sum
  to 1 and are constrained to the simplex.

For fixed line-compound phases the free energy is a single number — the
MLIP-relaxed total energy of the structure — which does not depend on
composition.

### Muggianu mixing enthalpy

The enthalpy is expanded in a Redlich–Kister (RK) polynomial using the Muggianu
extrapolation scheme, which is exact for binary subsystems and generalises
consistently to ternary and higher-order systems:

```
H(y) = Σ_i y_i H_i^0
       + Σ_{i<j} y_i y_j Σ_k L_ij^k (y_i − y_j)^k
       + L^{012} y_0 y_1 y_2          (ternary term, when present)
```

where H_i^0 are the pure-endpoint energies and L_ij^k are the RK interaction
coefficients fitted to SQS training data by `BladeTDBGen`.

### Oxygen chemical potential

The oxygen reservoir is described by its chemical potential μ_O (eV per O atom).
BLADE evaluates μ_O as a function of temperature and oxygen partial pressure
using the NIST Shomate equation for O₂:

```
μ_O(T, pO₂) = [ΔG⁰_O₂(T) + R T ln(pO₂)] / (2 × 96485 J mol⁻¹ eV⁻¹)
```

where ΔG⁰_O₂(T) is the standard Gibbs free energy of O₂ computed from
piecewise Shomate coefficients (valid 100–6000 K).  Internally μ_O is the
primary scan variable; `log10(pO₂)` is a convenience alias derived from the
same equation.

---

## Grand-Potential Linear Program

### Phase-space construction

The solver assembles a set of P candidate phases.  For a system with N metals
and a fixed non-metal element (e.g. B₂):

- The flexible phase family is discretized on a uniform simplex grid with step
  `y_step` (default 0.01), giving one "phase" per grid point.
- Fixed line-compound phases are loaded from the `fixed_phases_subdir`
  directory (MLIP-relaxed MP/ORB reference structures).

The grand potential of phase p at (T, μ_O) is:

```
Ω_p(T, μ_O) = G_p(T) − μ_O × n_O(p)
```

where n_O(p) is the number of oxygen atoms per formula unit of phase p.

### Linear program

At each (T, μ_O) point the equilibrium is found by solving:

```
minimize    Ω · n
subject to  A_eq n = b_eq
            n ≥ 0
```

where:

- **n** is the vector of phase amounts (one entry per candidate phase).
- **Ω** is the vector of grand potentials.
- **A_eq** is the elemental composition matrix: row i encodes the number of
  atoms of element i per formula unit of each phase.
- **b_eq** is the elemental conservation vector derived from the alloy
  composition and the fixed non-metal stoichiometry.

The LP is solved using SciPy's `revised simplex` method (falling back to HiGHS
for degenerate cases).  Phase amounts below 10⁻¹⁰ are zeroed out.

---

## Oxide Stability Criterion

A phase is classified as an oxide if it contains oxygen (n_O > 0).  After the
LP is solved the total oxide phase fraction is:

```
f_oxide = Σ_{p: n_O(p) > 0} n_p / Σ_p n_p
```

The parent (non-oxide) phase fraction is `f_parent = 1 − f_oxide`.

A system is judged to be **oxidizing** at (T, μ_O) when the absorbed-oxygen
metric exceeds the onset threshold:

```
absorbed_O = Σ_{p: n_O(p) > 0} n_p × n_O(p)  >  onset_threshold
```

The default threshold is 10⁻⁸ (configurable via `onset_threshold`).

---

## Phase Fraction Calculation

The raw LP output **n** gives phase amounts in units consistent with the
conservation constraints.  Normalized fractions are:

```
φ_p = n_p / Σ_q n_q
```

These fractions appear in all region maps, scan plots, and CSV tables.

Phases are grouped into families based on the dominant metal composition of the
flexible phase.  A phase is included in a family when its metal-sublattice
fraction exceeds a membership threshold (0.01 when `include_0p01_to_0p05_components`
is enabled, otherwise 0.05).  Contiguous regions in (μ_O, x) or (μ_O, T) space
that share the same phase assemblage are assigned a single region ID.

The area under the phase-fraction curve (AUC) integrated over μ_O is also
reported:

```
AUC = ∫ φ_oxide(μ_O) dμ_O
```

computed by the trapezoidal rule over the μ_O grid.

---

## Onset Definition

The **oxidation onset** at temperature T for a given metal composition is the
smallest μ_O value at which the absorbed-oxygen metric exceeds `onset_threshold`:

```
μ_O^onset(T, x) = inf{ μ_O : absorbed_O(T, μ_O, x) > onset_threshold }
```

A less negative (higher) onset μ_O means the alloy begins oxidizing under
conditions closer to ambient oxygen activity — lower oxidation resistance.
A more negative onset μ_O indicates stability under more oxidizing conditions.

The onset is evaluated on the full composition simplex grid at each temperature
in `temperature_values`.  For binary systems the onset is plotted as a 1-D
curve vs composition; for ternary systems a colour-mapped simplex diagram is
produced; for quaternary and higher systems only the CSV is written.

---

## Gibbs Energy Sources

Phase Gibbs energies come from two sources:

1. **Mixed boride phases** — MLIP-relaxed SQS energies fitted by BLADE's TDB
   generation stage and parameterized by the Muggianu/RK model.  Structures
   live under `mixed_phase_subdir` (default `"blade"`).
2. **Oxide reference phases** — DFT energies from the Materials Project,
   re-relaxed with ORB (or another MLIP) for energy-scale consistency.
   Structures live under `fixed_phases_subdir` (default `"ORB"`).

---

## Numerical Parameters

| Symbol | Config key | Default | Meaning |
|---|---|---|---|
| Δy | `phase_grid_step` | 0.01 | Simplex grid spacing for the flexible phase |
| k | `rk_order` | 3 | Redlich–Kister polynomial order |
| ε_act | `active_threshold` | 10⁻⁹ | Phase amount below which a phase is considered absent |
| ε_onset | `onset_threshold` | 10⁻⁸ | Absorbed-O threshold for onset declaration |
| ε_plot | `plot_threshold` | 10⁻⁴ | Phase fraction below which a phase is omitted from scan plots |
