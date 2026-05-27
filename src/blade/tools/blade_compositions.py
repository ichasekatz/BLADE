"""Composition generation for multicomponent materials systems.

This module provides :class:`BladeCompositions`, which enumerates all valid
combinations of primary and secondary element pools subject to
user-defined size and count constraints. The resulting composition list is
the primary input for downstream SQS generation and TDB fitting workflows.

For multicomponent systems, primary elements occupy the variable sublattice and
additional secondary elements provide further mixing.


Example::

    from blade.tools.blade_compositions import BladeCompositions

    composer = BladeCompositions(
        primary_elements=["Cr", "Hf", "Ta", "Ti", "Mo"],
        secondary_elements=[],
        system_size=3,
        primary_min=3, primary_max=3,
        secondary_min=0, secondary_max=0,
        allow_lower_order=False,
    )
    comps = composer.generate_compositions()
    print(comps)  # [['Cr', 'Hf', 'Mo'], ['Cr', 'Hf', 'Ta'], ...]
"""

from __future__ import annotations

import itertools

__author__ = "Chase Katz"


class BladeCompositions:
    """Generate element combinations for multicomponent materials systems.

    Constructs all allowed combinations of primary and secondary element pools
    subject to constraints on total system size and per-group element counts.
    The output feeds directly into :class:`~blade.tools.blade_sqsgen.BladeSQS`
    and :class:`~blade.tools.blade_tdb_gen.BladeTDBGen`.

    Attributes:
        primary_elements (list[str]): Primary element symbols.
        secondary_elements (list[str]): Secondary element symbols.
        system_size (int): Target number of elements per system.
        primary_min (int): Minimum number of primary elements per system.
        primary_max (int): Maximum number of primary elements per system.
        secondary_min (int): Minimum number of secondary elements per system.
        secondary_max (int): Maximum number of secondary elements per system.
        allow_lower_order (bool): Whether to include sub-``system_size`` systems.
        compositions (list[list[str]]): Populated after calling :meth:`generate_compositions`.
    """

    def __init__(
        self,
        primary_elements: list[str],
        secondary_elements: list[str],
        system_size: int,
        primary_min: int,
        primary_max: int,
        secondary_min: int,
        secondary_max: int,
        allow_lower_order: bool,
    ) -> None:
        """Initialize BladeCompositions.

        Args:
            primary_elements (list[str]): Primary element symbols for the
                variable sublattice (e.g., ``["Cr", "Hf", "Ta"]``; transition
                e.g., transition metals).
            secondary_elements (list[str]): Secondary element symbols
                (e.g., ``["Y", "La"]``; e.g., secondary elements).
                Pass an empty list if not needed.
            system_size (int): Target number of elements per chemical system
                (e.g., ``3`` for ternary systems).
            primary_min (int): Minimum number of primary elements in a system.
            primary_max (int): Maximum number of primary elements in a system.
            secondary_min (int): Minimum number of secondary elements in a system.
            secondary_max (int): Maximum number of secondary elements in a system.
            allow_lower_order (bool): If ``True``, include systems with fewer
                than ``system_size`` elements. If ``False``, only include
                systems with exactly ``system_size`` elements.
        """
        self.primary_elements = primary_elements
        self.secondary_elements = secondary_elements
        self.system_size = system_size
        self.primary_min = primary_min
        self.primary_max = primary_max
        self.secondary_min = secondary_min
        self.secondary_max = secondary_max
        self.allow_lower_order = allow_lower_order
        self.compositions: list[list[str]] = []

    def generate_compositions(self) -> list[list[str]]:
        """Generate all valid chemical compositions.

        Enumerates combinations of primary and secondary elements that satisfy
        the constraints set at initialization. Each composition is
        alphabetically sorted. Results are also stored in
        ``self.compositions`` for use by :meth:`get_systems`.

        Returns:
            list[list[str]]: Sorted list of compositions, where each
            composition is a list of element symbols
            (e.g., ``[['Cr', 'Hf', 'Ta'], ['Cr', 'Hf', 'Ti'], ...]``).
        """
        primary_counts = list(range(self.primary_min, self.primary_max + 1))
        secondary_counts = list(range(self.secondary_min, self.secondary_max + 1))

        primary_combos: list[list[str]] = []
        secondary_combos: list[list[str]] = []
        combined_comps: list[list[str]] = []
        compositions: list[list[str]] = []

        for i in primary_counts:
            if i == 0:
                for j in secondary_counts:
                    if j == 0:
                        compositions += [[""]]
                    else:
                        combined_comps += [
                            list(c) for c in itertools.combinations(self.secondary_elements, j)
                        ]
            else:
                primary_combos += [
                    list(c) for c in itertools.combinations(self.primary_elements, i)
                ]

        for j in secondary_counts:
            if j == 0:
                combined_comps += primary_combos
            else:
                secondary_combos += [
                    list(c) for c in itertools.combinations(self.secondary_elements, j)
                ]

        for p_comp in primary_combos:
            for s_comp in secondary_combos:
                combined_comps.append(list(p_comp) + s_comp)

        combined_comps = [c for c in combined_comps if len(c) <= self.system_size]

        for comp in combined_comps:
            compositions.append(sorted(comp))

        if not self.allow_lower_order:
            compositions = [c for c in compositions if len(c) == self.system_size]

        compositions = sorted(compositions)
        self.compositions = compositions
        return compositions

    def get_systems(self) -> set[int]:
        """Return the unique system sizes present in the generated compositions.

        Must be called after :meth:`generate_compositions`.

        Returns:
            set[int]: Set of integers representing the number of elements in
            each unique chemical system (e.g., ``{2, 3}`` for a mixed
            binary/ternary search).

        Raises:
            AttributeError: If called before :meth:`generate_compositions`.
        """
        len_comps = [len(c) for c in self.compositions]
        return set(len_comps) if len(len_comps) >= 2 else {len_comps[0]}
