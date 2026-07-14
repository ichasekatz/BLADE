"""Composition generation and BLADE CONTCAR scanning for oxide screening.

Enumerates N-component compositions from a primary element pool, scans
BLADE-relaxed structures for energies, and writes reference data used by
downstream oxidation analysis scripts.

The values below demonstrate one prototype. Adapt them for any crystal system:
  - primary_elements / secondary_elements: element pools
  - fixed_elements: any species held fixed on the structural sublattice
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
    # Primary variable-sublattice elements for this example
    primary_elements=["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"],
    # Optional secondary variable-sublattice elements. Do not duplicate
    # fixed_elements or oxygen_element here.
    secondary_elements=[],
    # Oxidizing species symbol — added to composition pool when oxygen flags are True
    oxygen_element="O",
    # Primary element count per composition
    primary_min=0,
    primary_max=3,
    # Secondary element count per composition
    secondary_min=0,
    secondary_max=0,
    # Include compounds with neither oxygen nor fixed structural species.
    include_no_oxygen=True,
    # Include compounds containing both fixed_elements and oxygen.
    include_fixed_oxygen=True,
    # Include compositions formed by adding O to a base compound
    include_added_oxygen=True,
    # Include base-phase compounds without oxygen.
    include_fixed=True,
    # Fixed structural species for this demonstration. Use an empty set when
    # the modeled phase has no fixed species.
    fixed_elements=frozenset({"B"}),
)
comp.run()
