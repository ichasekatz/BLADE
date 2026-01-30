# BLADE: Boride Learning and Design Engine

High-Entropy Diboride Oxidation Screening Framework

## Overview

This framework is designed to predict and screen **oxidation-relevant properties of high-entropy diborides (HEDBs)** using a combination of **special quasirandom structures (SQS)**, **first-principles energetics**, and **CALPHAD-style fitting**. [1-5]

High-entropy diborides are highly customizable ceramics in which:

* **Boron is held fixed** on the middle sublattice
* **Metal species are substituted** on the top and bottom sublattices
* The parent AlB₂-type prototype contains **14 total atoms per unit cell**

This structural flexibility leads to a vast composition space, which is efficiently explored using automated SQS generation and thermodynamic fitting.

<img width="700" alt="BLADE_Framework" src="https://github.com/user-attachments/assets/8a64271d-0427-4e2e-b964-df0af7ff18d0" />

**Figure 1.** Overview of the current BLADE framework and workflow.

---

## Code Structure and Workflow

---

### `blade_compositions.py`

This module generates **composition permutations** for a given system by specifying:

* Allowed transition-metal elements
* Optional rare-earth elements
* Target system size (number of metal species)
* Minimum/maximum number of transition metals and rare earths

The output is a list of unique metal-sublattice compositions consistent with the specified constraints.

---

### Lattice and SQS Specification

Users define:

* Lattice parameters for the AlB₂-type parent structure
* SQS composition levels (e.g., endmembers, midpoints, off-stoichiometric points)

These inputs control which compositions are explicitly sampled and how well random alloy statistics are approximated.

<img width="700" alt="HEDB_Lattice" src="https://github.com/user-attachments/assets/b3b0bf3b-051a-4d18-82a5-fe79b9f9441e" />

**Figure 2.** Schematic of the HEDB lattice. Boron is held fixed on the central sublattice, and metal atoms are substituted on the outer sublattices. Reproduced from Gild et al. (2016) under the Creative Commons Attribution (CC BY) license. [1]

---

### `blade_sqs.py`

This module generates **SQS supercells** for each composition using ATAT-style SQS construction.

Key features include:

* Automatic calculation of the **minimum supercell size** required for a given composition and SQS level
* Generation of SQS structures consistent with the specified lattice and sublattice disorder
* The `supercell_size` function reports the number of atoms (or unit cells) required for each SQS

Each generated SQS is intended to be relaxed using first-principles calculations prior to fitting.

---

### `blade_tdb_gen.py`

This module relaxes SQS structures, performs first-principles calculations, and organizes results for thermodynamic database construction.

Key features include:
* First-principles relaxation and fitting managed by `Sqs2tdb`, which relaxes the generated SQS structures.
* Produces reusable thermodynamic models suitable for phase-diagram and phase-stability calculations.

The output is a .tdb file that can be used to describe the thermodynamic properties of materials.

---

### Thermodynamic Fitting and Phase Diagrams

After SQS structures are generated and relaxed:

* Each composition is fit using **sqs2tdb**, producing a **CALPHAD-style `.tdb` file**
* For binary metal systems, a **phase diagram** is generated to visualize:

  * Single-phase regions
  * Miscibility gaps
  * Composition–temperature stability trends

These phase diagrams are used to identify metal combinations that are likely to form stable or metastable HEDB solid solutions.

<img width="700" alt="CR_NB_Phase_Diagram" src="https://github.com/user-attachments/assets/9f68e914-e764-49af-91bd-a04f1a9e71c0" />

**Figure 3.** Example phase diagram generated with the BLADE framework showing a miscibility gap in the Cr–Nb–B₂ system.

---

## Outputs

* Relaxed SQS structures (POSCAR format)
* CALPHAD-compatible `.tdb` files
* Binary phase diagrams (where applicable)
* Composition-level miscibility and stability screening data

---

## Intended Use

This framework is intended for **high-throughput screening**, not direct prediction of full oxidation kinetics.
It is particularly useful for:

* Identifying metal combinations with favorable mixing behavior
* Ranking candidate HEDBs prior to more expensive oxidation modeling
* Guiding experimental alloy selection

<img width="700" alt="Multiverse" src="https://github.com/user-attachments/assets/f8de55aa-16ba-427d-baf9-22373655f391" />

**Figure 4.** Overview of the BLADE framework and its integration with external tools and supporting frameworks (the Multi-verse).

---

## Example Usage

An example script demonstrating SQS generation and TDB construction is provided in `examples/tdb_gen.py`.

---

## Future Work

Planned extensions of the framework include:

1. **Explicit oxide modeling**

   * Oxide structures will be introduced in place of HEDB lattices
   * Oxide thermodynamics will be evaluated alongside diboride phases

2. **Volume and lattice-parameter matching**

   * Optimized POSCAR files will be used to compute unit-cell volumes
   * Oxide and diboride lattice parameters will be compared to assess:

     * Volume change upon oxidation
     * Lattice mismatch and strain accommodation

3. **Oxide adhesion and stability metrics**

   * Favorable volume and lattice matching will be used as proxies for oxide adhesion
   * Compositions producing mechanically stable oxide scales will be prioritized

4. **Oxygen diffusivity and oxide-scale growth**

   * Oxygen transport will be evaluated to estimate oxide-scale thickness and growth behavior
   * This will enable ranking of compositions by predicted oxidation resistance

---

## Notes

* Boron is treated as a **fixed sublattice** in the current implementation and is implicit in the thermodynamic fitting.
* The generated `.tdb` files describe **metal-sublattice mixing behavior** within the HEDB prototype and do not represent a full multi-phase thermodynamic database.

---

## Citations

[1] J. Gild et al., “High-entropy metal diborides: A new class of ultrahigh-temperature ceramics,” Scientific Reports, 6, 37946 (2016).

[2] S. Zhu et al., “Machine learning potentials for alloys: A workflow to predict phase diagrams and benchmark accuracy,” arXiv preprint (2025).

[3] H. Tang et al., “Rare-earth compositional screening of high-entropy diborides for improved oxidation resistance,” Journal of the American Ceramic Society, 108, e20123 (2025).

[4] C. Kunselman et al., “Analytical gradient-based optimization of CALPHAD model parameters,” arXiv preprint (2025).

[5] C. Kunselman et al., “Construction and tuning of CALPHAD models using machine-learned interatomic potentials: Pt–W system,” arXiv preprint (2025).
