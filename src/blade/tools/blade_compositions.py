"""
This module defines the `BladeCompositions` class for generating material compositions based on specified
transition metals and rare earth elements.

The `BladeCompositions` class constructs compositional subspaces for multicomponent material systems by
enumerating combinations of transition metals and rare-earth elements subject to user-defined constraints
on system size and minimum/maximum element counts. This functionality is primarily used to generate
candidate chemical systems for high-throughput materials discovery and design workflows.
"""
import itertools


class BladeCompositions:
    """
    Generate element combinations for multicomponent materials systems.

    This class constructs all allowed combinations of transition metals and rare-earth elements given constraints on:
      - total system size
      - minimum and maximum number of transition metals
      - minimum and maximum number of rare-earth elements
      - whether lower-order systems (e.g., binaries in a ternary search) are permitted

    The output is a sorted list of element lists representing unique chemical systems.
    """
    def __init__(
        self,
        transition_metals,
        rare_earths,
        system_size,
        tm_min,
        tm_max,
        re_min,
        re_max,
        allow_lower_order,
    ):
        """
        Initializes the `BladeCompositions` object.

        Args:
            transition_metals (list[str]): List of transition metal element symbols (e.g., ["Cr", "Hf", "Ta"]).
            rare_earths (list[str]): List of rare-earth element symbols.
            system_size (int): Target number of elements in the chemical system (e.g., 3 for ternary systems).
            tm_min (int): Minimum number of transition metals allowed in a system.
            tm_max (int): Maximum number of transition metals allowed in a system.
            re_min (int): Minimum number of rare-earth elements allowed in a system.
            re_max (int): Maximum number of rare-earth elements allowed in a system.
            allow_lower_order (bool): If True, include systems with fewer than `system_size` elements.
                If False, only include systems with exactly `system_size` elements.
        """
        self.transition_metals = transition_metals
        self.rare_earths = rare_earths
        self.system_size = system_size
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.re_min = re_min
        self.re_max = re_max
        self.allow_lower_order = allow_lower_order

    def generate_compositions(self):
        """
        Generate all valid chemical compositions.

        This method enumerates combinations of transition metals and rare-earth elements that satisfy
        the constraints specified at initialization. Compositions are alphabetically sorted and filtered by system size.

        Returns:
            list[list[str]]: A sorted list of compositions, where each composition is represented as a list of element symbols.
        """
        tm_elements = list(range(self.tm_min, self.tm_max + 1))
        re_elements = list(range(self.re_min, self.re_max + 1))

        tm_combos = []
        re_combos = []
        combined_comps = []
        compositions = []

        for i in tm_elements:
            if i == 0:
                for j in re_elements:
                    if j == 0:
                        compositions += [[""]]
                    else:
                        combined_comps += [
                            list(c) for c in itertools.combinations(self.rare_earths, j)
                        ]
            if i != 0:
                tm_combos += [list(c) for c in itertools.combinations(self.transition_metals, i)]

        for j in re_elements:
            if j == 0:
                combined_comps += tm_combos
            else:
                re_combos += [list(c) for c in itertools.combinations(self.rare_earths, j)]

        for tm_comp in tm_combos:
            for re_comp in re_combos:
                combined_comps.append(list(tm_comp) + re_comp)

        combined_comps = [c for c in combined_comps if len(c) <= self.system_size]

        for i in combined_comps:
            compositions += [sorted(i)]
        if not self.allow_lower_order:
            compositions = [c for c in compositions if len(c) == self.system_size]
        compositions = sorted(compositions)

        self.compositions = compositions
        return compositions

    def get_systems(self):
        """
        Return the unique system sizes present in the generated compositions.

        Returns:
            set[int]:
                A set of integers corresponding to the number of elements in each unique chemical system.
        """
        compositions = self.compositions
        len_comps = []
        for i in compositions:
            len_comps += [len(i)]
        unique_len_comps = set(len_comps) if len(len_comps) >= 2 else {len_comps[0]}
        return unique_len_comps
