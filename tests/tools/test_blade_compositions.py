"""Tests for BladeCompositions."""

from __future__ import annotations

import pytest

from blade.tools.blade_compositions import BladeCompositions


class TestGenerateCompositionsBinary:
    """generate_compositions() with a pure-primary binary system."""

    def test_binary_count_is_correct(self) -> None:
        """Three primary elements taken 2 at a time yields exactly 3 compositions."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Ta"],
            secondary_elements=[],
            primary_min=2,
            primary_max=2,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        assert len(result) == 3

    def test_binary_each_entry_is_sorted(self) -> None:
        """Each binary composition is alphabetically sorted."""
        composer = BladeCompositions(
            primary_elements=["Ta", "Cr", "Hf"],
            secondary_elements=[],
            primary_min=2,
            primary_max=2,
            secondary_min=0,
            secondary_max=0,
        )
        for comp in composer.generate_compositions():
            assert comp == sorted(comp)

    def test_binary_all_pairs_present(self) -> None:
        """All three pairwise combinations of three elements are generated."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Ta"],
            secondary_elements=[],
            primary_min=2,
            primary_max=2,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        assert ["Cr", "Hf"] in result
        assert ["Cr", "Ta"] in result
        assert ["Hf", "Ta"] in result


class TestGenerateCompositionsTernary:
    """generate_compositions() with a pure-primary ternary system."""

    def test_ternary_count_is_correct(self) -> None:
        """Five primary elements taken 3 at a time yields 10 compositions."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Mo", "Ta", "Ti"],
            secondary_elements=[],
            primary_min=3,
            primary_max=3,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        assert len(result) == 10

    def test_ternary_all_entries_length_three(self) -> None:
        """Every ternary composition contains exactly three elements."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Mo", "Ta", "Ti"],
            secondary_elements=[],
            primary_min=3,
            primary_max=3,
            secondary_min=0,
            secondary_max=0,
        )
        for comp in composer.generate_compositions():
            assert len(comp) == 3


class TestGenerateCompositionsMixedSize:
    """generate_compositions() spanning multiple system sizes."""

    def test_mixed_contains_both_sizes(self) -> None:
        """Binary-to-ternary range over three elements yields both 2- and 3-element systems."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Ta"],
            secondary_elements=[],
            primary_min=2,
            primary_max=3,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        sizes = {len(c) for c in result}
        assert 2 in sizes
        assert 3 in sizes

    def test_mixed_total_count(self) -> None:
        """Binary-to-ternary over three elements yields 3 + 1 = 4 compositions."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Ta"],
            secondary_elements=[],
            primary_min=2,
            primary_max=3,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        assert len(result) == 4


class TestGenerateCompositionsWithSecondary:
    """generate_compositions() with secondary elements."""

    def test_secondary_elements_included_in_output(self) -> None:
        """Compositions including a secondary element appear when secondary_max >= 1."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf"],
            secondary_elements=["Y"],
            primary_min=2,
            primary_max=2,
            secondary_min=0,
            secondary_max=1,
        )
        result = composer.generate_compositions()
        assert any("Y" in c for c in result)

    def test_without_secondary_subset_present(self) -> None:
        """The pure-primary composition is still present when secondary_min=0."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf"],
            secondary_elements=["Y"],
            primary_min=2,
            primary_max=2,
            secondary_min=0,
            secondary_max=1,
        )
        result = composer.generate_compositions()
        assert ["Cr", "Hf"] in result

    def test_compositions_stored_on_self(self) -> None:
        """generate_compositions() stores the result in self.compositions."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf"],
            secondary_elements=[],
            primary_min=2,
            primary_max=2,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        assert composer.compositions == result


class TestGetSystems:
    """get_systems() returns the correct set of system sizes."""

    def test_pure_binary_returns_singleton_set(self) -> None:
        """A binary-only search returns {2}."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Ta"],
            secondary_elements=[],
            primary_min=2,
            primary_max=2,
            secondary_min=0,
            secondary_max=0,
        )
        composer.generate_compositions()
        assert composer.get_systems() == {2}

    def test_mixed_range_returns_both_sizes(self) -> None:
        """A binary-to-ternary search returns {2, 3}."""
        composer = BladeCompositions(
            primary_elements=["Cr", "Hf", "Ta"],
            secondary_elements=[],
            primary_min=2,
            primary_max=3,
            secondary_min=0,
            secondary_max=0,
        )
        composer.generate_compositions()
        assert composer.get_systems() == {2, 3}


class TestEdgeCases:
    """Edge cases for BladeCompositions."""

    def test_single_element_pool_min_equals_max(self) -> None:
        """Single-element primary pool with min=max=1 yields one composition."""
        composer = BladeCompositions(
            primary_elements=["Cr"],
            secondary_elements=[],
            primary_min=1,
            primary_max=1,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        assert result == [["Cr"]]

    @pytest.mark.parametrize(
        ("n_elements", "k", "expected"),
        [
            (4, 2, 6),
            (4, 3, 4),
            (5, 2, 10),
        ],
    )
    def test_combination_count_matches_formula(self, n_elements: int, k: int, expected: int) -> None:
        """C(n, k) compositions are generated for n primary elements taken k at a time."""
        elements = ["A", "B", "C", "D", "E"][:n_elements]
        composer = BladeCompositions(
            primary_elements=elements,
            secondary_elements=[],
            primary_min=k,
            primary_max=k,
            secondary_min=0,
            secondary_max=0,
        )
        result = composer.generate_compositions()
        assert len(result) == expected
