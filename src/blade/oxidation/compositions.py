"""OxideCompositions — composition enumeration and BLADE CONTCAR scanning.

Encapsulates the logic from scripts/01_compositions.py:
    - Section 1: enumerate systems, build composition table, save composition_list.xlsx
    - Section 2: scan BLADE CONTCARs, compute dH, save blade_generated_data.xlsx
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pymatgen.core import Composition, Structure

from blade.tools.blade_compositions import BladeCompositions


# ---------------------------------------------------------------------------
# Module-level helper
# ---------------------------------------------------------------------------

def calculate_dH(formula: str, energy_per_atom: float, refs: dict) -> float:
    comp = Composition(formula)
    total_atoms = comp.num_atoms
    ref_total = 0.0
    for el, amount in comp.get_el_amt_dict().items():
        if refs.get(el) is None:
            raise ValueError(f"Missing reference energy for {el}")
        ref_total += amount * refs[el]
    return energy_per_atom - (ref_total / total_atoms)


# ---------------------------------------------------------------------------
# Class
# ---------------------------------------------------------------------------

class OxideCompositions:
    """Enumerate oxide/boride compositions and scan BLADE relaxation data.

    Parameters
    ----------
    files_dir:
        Root output directory; ``composition_list.xlsx`` and
        ``blade_generated_data.xlsx`` are written here.
    primary_elements:
        Transition metals occupying the variable sublattice.
    secondary_elements:
        Fixed/mixed species (typically ``["B", "O"]``).
    mlip_ref_label:
        Human-readable label for the MLIP reference set (used in print output).
    mlip_ref_folder:
        Directory containing MLIP-relaxed single-element CONTCAR + energy files
        used to build element reference energies.
    blade_comp_dir:
        Directory tree of BLADE SQS relaxation outputs (CONTCARs + energy files).
    primary_min:
        Minimum number of primary elements per composition.
    primary_max:
        Maximum number of primary elements per composition.
    include_no_oxygen:
        Include compositions with no O and no B (pure metals).
    include_fixed_oxygen:
        Include compositions containing both B and O.
    include_added_oxygen:
        Include compositions containing O but not B (pure oxides).
    include_fixed:
        Include compositions containing the fixed structural element(s) but not O.
    fixed_elements:
        Set of elements treated as fixed structural species (e.g. ``{"B"}`` for
        diborides, ``{"C"}`` for MAX phases).  Defaults to ``frozenset({"B"})``
        so existing diboride workflows are unchanged.
    """

    def __init__(
        self,
        files_dir: Path,
        primary_elements: list[str],
        secondary_elements: list[str],
        primary_min: int = 0,
        primary_max: int = 2,
        secondary_min: int = 0,
        secondary_max: int = 2,
        include_no_oxygen: bool = True,
        include_fixed_oxygen: bool = True,
        include_added_oxygen: bool = True,
        include_fixed: bool = False,
        fixed_elements: frozenset[str] = frozenset({"B"}),
        oxygen_element: str = "O",
    ) -> None:
        self.files_dir = files_dir
        self.primary_elements = primary_elements
        self.secondary_elements = secondary_elements
        self.primary_min = primary_min
        self.primary_max = primary_max
        self.secondary_min = secondary_min
        self.secondary_max = secondary_max
        self.include_no_oxygen = include_no_oxygen
        self.include_fixed_oxygen = include_fixed_oxygen
        self.include_added_oxygen = include_added_oxygen
        self.include_fixed = include_fixed
        self.fixed_elements = fixed_elements
        self.oxygen_element = oxygen_element

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def run(self) -> None:
        """Generate and save composition_list.xlsx."""
        self.generate_composition_list()

    def generate_composition_list(self) -> None:
        """Section 1 — enumerate systems and save ``composition_list.xlsx``."""
        self.files_dir.mkdir(parents=True, exist_ok=True)

        # Generate primary + user secondary combinations (secondary_min/max control A-site mixing only)
        composer = BladeCompositions(
            primary_elements=self.primary_elements,
            secondary_elements=list(self.secondary_elements),
            primary_min=self.primary_min,
            primary_max=self.primary_max,
            secondary_min=self.secondary_min,
            secondary_max=self.secondary_max,
        )

        base_compositions = composer.generate_compositions()

        # Augment with fixed_elements and oxygen independently of secondary_max
        seen: set[tuple] = {tuple(c) for c in base_compositions}
        composition_list = list(base_compositions)
        oxygen = self.oxygen_element
        fixed = sorted(self.fixed_elements)

        for comp in base_compositions:
            comp_set = set(comp)
            if self.include_fixed:
                for el in fixed:
                    new = tuple(sorted(comp_set | {el}))
                    if new not in seen:
                        seen.add(new)
                        composition_list.append(list(new))
            if self.include_added_oxygen and oxygen not in comp_set:
                new = tuple(sorted(comp_set | {oxygen}))
                if new not in seen:
                    seen.add(new)
                    composition_list.append(list(new))
            if self.include_fixed_oxygen:
                for el in fixed:
                    new = tuple(sorted(comp_set | {el, oxygen}))
                    if new not in seen:
                        seen.add(new)
                        composition_list.append(list(new))
        unique_len_comps = composer.get_systems()

        print(f"Compositions: {composition_list}")
        print(f"Total # compositions: {len(composition_list)}")
        print(f"Unique system sizes: {unique_len_comps}")

        # Add single-element compositions for all primary elements so pure
        # references are present.
        non_primary = set(self.secondary_elements) | self.fixed_elements | {self.oxygen_element}

        for el in self.primary_elements:
            if [el] not in composition_list:
                composition_list.append([el])

        # Build composition table.
        rows: list[dict] = []
        max_len = max(len(comp) for comp in composition_list)

        for comp in composition_list:
            comp_sorted = sorted(comp)
            has_O = self.oxygen_element in comp_sorted
            has_fixed = bool(set(comp_sorted) & self.fixed_elements)
            is_single = len(comp_sorted) == 1

            if is_single:
                pass  # always include pure-element reference rows
            elif has_O and has_fixed:
                if not self.include_fixed_oxygen:
                    continue
            elif has_O and not has_fixed:
                if not self.include_added_oxygen:
                    continue
            elif has_fixed and not has_O:
                if not self.include_fixed:
                    continue
            else:  # no O, no fixed element — pure metals
                if not self.include_no_oxygen:
                    continue

            metals_only = sorted(el for el in comp if el not in non_primary)
            both_only = sorted(el for el in comp if el in non_primary)

            row: dict = {
                "composition": "".join(comp_sorted),
                "metal_composition": "".join(metals_only),
                "both_composition": "".join(both_only),
                "elements": ",".join(comp_sorted),
                "n_elements": len(comp_sorted),
            }
            for i in range(max_len):
                row[f"element_{i+1}"] = comp_sorted[i] if i < len(comp_sorted) else ""
            rows.append(row)

        comp_df = pd.DataFrame(rows)
        comp_df = comp_df.drop_duplicates()
        comp_df = comp_df.sort_values(
            ["metal_composition", "both_composition", "composition"]
        ).reset_index(drop=True)

        comp_df.to_excel(self.files_dir / "composition_list.xlsx", index=False)
        print(f"Saved {(self.files_dir / 'composition_list.xlsx').resolve()}")
