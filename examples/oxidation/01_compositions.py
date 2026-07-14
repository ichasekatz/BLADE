"""Composition generation and BLADE CONTCAR scanning for oxide screening.

Enumerates N-component compositions from a primary element pool, scans
BLADE-relaxed structures for energies, and writes reference data used by
downstream oxidation analysis scripts.

Adapt this script for any crystal system by adjusting:
  - primary_elements / secondary_elements: element pools
  - fixed_elements: non-metal fixed species on the lattice
      {"B"} for diborides, {"C"} for MAX phases, frozenset() for pure alloys
  - blade_comp_dir: Comps folder from the relevant BLADE TDB run
  - primary_min / primary_max: composition size range

Output files (written to files_dir/):
    composition_list.xlsx      — enumerated compositions
    blade_generated_data.xlsx  — energies and formation enthalpies from BLADE
"""

from pathlib import Path
from blade.oxidation.compositions import OxideCompositions

path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
files_dir = path1 / "Files"

comp = OxideCompositions(
    files_dir=files_dir,

    # Primary (M-site) metals — variable sublattice
    primary_elements=["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"],

    # A-site (or secondary sublattice) elements — do NOT include fixed_elements or O here.
    # C is added automatically from fixed_elements; O is added via oxygen_element below.
    secondary_elements=[],

    # Oxidizing species symbol — added to composition pool when oxygen flags are True
    oxygen_element="O",

    # Primary element count per composition
    primary_min=0,
    primary_max=3,
    # Secondary element count per composition
    secondary_min=0,
    secondary_max=0,

    # Include pure metal/alloy compounds — no oxygen AND no fixed elements (e.g. Zr, ZrHf, ZrHfAl)
    include_no_oxygen=True,

    # Include compounds containing both fixed_elements AND oxygen (e.g. CO2 for MAX, BO2 for diborides)
    # Set True to screen oxy-compounds of the fixed non-metal species
    include_fixed_oxygen=True,

    # Include compositions formed by adding O to a base compound
    include_added_oxygen=True,

    # Include base-phase compounds without oxygen (e.g. the parent MAX/diboride/alloy phase)
    include_fixed=True,

    # Non-metal fixed species on the lattice:
    #   frozenset({"B"})  → diboride (AlB2-type)
    #   frozenset({"C"})  → MAX phase (M2AX with X=C)
    #   frozenset()       → pure metallic alloy
    fixed_elements=frozenset({"B"}),
)
comp.run()
