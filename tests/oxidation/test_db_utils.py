"""Tests for blade.oxidation._db_utils pure helper functions."""

from __future__ import annotations

import pytest

from blade.oxidation._db_utils import (
    _fmt_stable,
    _is_element_col,
    _row_els,
    calculate_dH,
    clean_file_name,
    get_formula_counts,
    get_formula_elements,
    is_oxide_formula,
    is_single_element_formula,
    is_true_stable,
    normalize_formula,
    oxide_allowed_for_elements,
    split_list,
)


class TestNormalizeFormula:
    """normalize_formula() reduces formulas and handles bad input."""

    def test_simple_oxide_reduced(self) -> None:
        """Fe2O3 reduces to Fe2O3 (already reduced)."""
        assert normalize_formula("Fe2O3") == "Fe2O3"

    def test_reducible_formula(self) -> None:
        """Fe4O6 reduces to Fe2O3."""
        assert normalize_formula("Fe4O6") == "Fe2O3"

    def test_single_element(self) -> None:
        """Single element returns unchanged symbol."""
        assert normalize_formula("Fe") == "Fe"

    def test_invalid_formula_returns_stripped(self) -> None:
        """Unrecognised string returns stripped input."""
        result = normalize_formula("  ???  ")
        assert result == "???"


class TestCalculateDH:
    """calculate_dH() computes formation enthalpy per atom."""

    def test_pure_element_dH_is_zero(self) -> None:
        """DH for a pure element equals energy minus its own reference → 0."""
        refs = {"Fe": -8.0}
        dH = calculate_dH("Fe", -8.0, refs)
        assert dH == pytest.approx(0.0, abs=1e-10)

    def test_binary_oxide_dH_sign(self) -> None:
        """DH is negative when compound is more stable than elemental refs."""
        refs = {"Fe": -8.0, "O": -4.0}
        # Fe2O3: 5 atoms, ref_total=(2*-8+3*-4)/5=-5.6
        # dH = -7.0 - (-5.6) = -1.4
        dH = calculate_dH("Fe2O3", -7.0, refs)
        assert dH == pytest.approx(-1.4, abs=1e-9)

    def test_dH_zero_when_energy_equals_refs(self) -> None:
        """DH is zero when energy per atom matches weighted reference."""
        refs = {"Cr": -10.0, "Hf": -12.0}
        # CrHf: 2 atoms, ref per atom = (-10-12)/2 = -11
        dH = calculate_dH("CrHf", -11.0, refs)
        assert dH == pytest.approx(0.0, abs=1e-10)


class TestGetFormulaElements:
    """get_formula_elements() returns element sets."""

    def test_binary_oxide(self) -> None:
        """Fe2O3 elements are {Fe, O}."""
        assert get_formula_elements("Fe2O3") == {"Fe", "O"}

    def test_ternary(self) -> None:
        """CrHfTa has three elements."""
        assert get_formula_elements("CrHfTa") == {"Cr", "Hf", "Ta"}

    def test_single_element(self) -> None:
        """Pure element returns singleton set."""
        assert get_formula_elements("Al") == {"Al"}

    def test_invalid_returns_empty_set(self) -> None:
        """Garbage input returns empty set."""
        assert get_formula_elements("!!!") == set()


class TestIsOxideFormula:
    """is_oxide_formula() detects oxide compounds."""

    def test_fe2o3_is_oxide(self) -> None:
        """Fe2O3 is an oxide."""
        assert is_oxide_formula("Fe2O3") is True

    def test_pure_o2_is_not_oxide(self) -> None:
        """O2 is not an oxide (only one element)."""
        assert is_oxide_formula("O2") is False

    def test_metal_without_oxygen_is_not_oxide(self) -> None:
        """CrHf has no oxygen."""
        assert is_oxide_formula("CrHf") is False

    def test_tio2_is_oxide(self) -> None:
        """TiO2 is an oxide."""
        assert is_oxide_formula("TiO2") is True


class TestOxideAllowedForElements:
    """oxide_allowed_for_elements() checks cation membership."""

    def test_allowed_oxide(self) -> None:
        """Fe2O3 is allowed when Fe is in the allowed set."""
        assert oxide_allowed_for_elements("Fe2O3", {"Fe"}) is True

    def test_disallowed_cation(self) -> None:
        """Cr2O3 is not allowed when Cr is absent from allowed elements."""
        assert oxide_allowed_for_elements("Cr2O3", {"Fe"}) is False

    def test_non_oxide_returns_false(self) -> None:
        """CrHf has no oxygen, always returns False."""
        assert oxide_allowed_for_elements("CrHf", {"Cr", "Hf"}) is False

    def test_multi_cation_all_allowed(self) -> None:
        """FeAlO3 is allowed when both Fe and Al are in allowed set."""
        assert oxide_allowed_for_elements("FeAlO3", {"Fe", "Al"}) is True


class TestSplitList:
    """split_list() parses comma-separated cell values."""

    def test_basic_split(self) -> None:
        """Comma-separated string splits into trimmed tokens."""
        assert split_list("a, b, c") == ["a", "b", "c"]

    def test_nan_returns_empty(self) -> None:
        """Float NaN returns empty list."""
        assert split_list(float("nan")) == []

    def test_empty_string_returns_empty(self) -> None:
        """Empty or whitespace-only string returns empty list."""
        assert split_list("") == []
        assert split_list("   ") == []

    def test_single_token(self) -> None:
        """Single value without comma returns single-element list."""
        assert split_list("alpha") == ["alpha"]


class TestIsTrueStable:
    """is_true_stable() converts truthy representations."""

    @pytest.mark.parametrize("val", [True, "true", "TRUE", "1", "yes", "Yes"])
    def test_truthy_values(self, val) -> None:
        """Recognized truthy strings and bool True return True."""
        assert is_true_stable(val) is True

    @pytest.mark.parametrize("val", [False, "false", "FALSE", "0", "no"])
    def test_falsy_values(self, val) -> None:
        """Boolean False and non-truthy strings return False."""
        assert is_true_stable(val) is False

    def test_nan_returns_false(self) -> None:
        """NaN returns False."""
        assert is_true_stable(float("nan")) is False


class TestIsElementCol:
    """_is_element_col() recognises valid element symbols."""

    @pytest.mark.parametrize("symbol", ["Fe", "O", "Cr", "Hf", "Ta", "B"])
    def test_valid_element_symbols(self, symbol: str) -> None:
        """Real element symbols return True."""
        assert _is_element_col(symbol) is True

    @pytest.mark.parametrize("col", ["energy", "formula", "dH", "source", "Xx"])
    def test_non_element_columns(self, col: str) -> None:
        """Non-element column names return False."""
        assert _is_element_col(col) is False


class TestCleanFileName:
    """clean_file_name() sanitises names for the filesystem."""

    def test_spaces_become_underscores(self) -> None:
        """Spaces are replaced with underscores."""
        assert clean_file_name("my formula") == "my_formula"

    def test_reserved_chars_replaced(self) -> None:
        """Windows-reserved characters become underscores."""
        assert "<" not in clean_file_name("a<b>c")
        assert ">" not in clean_file_name("a<b>c")

    def test_blank_returns_unknown_parent(self) -> None:
        """Blank string returns 'unknown_parent'."""
        assert clean_file_name("") == "unknown_parent"

    def test_nan_string_returns_unknown_parent(self) -> None:
        """'nan' and 'NaN' return 'unknown_parent'."""
        assert clean_file_name("nan") == "unknown_parent"
        assert clean_file_name("NaN") == "unknown_parent"


class TestGetFormulaCounts:
    """get_formula_counts() returns element → count mapping."""

    def test_fe2o3_counts(self) -> None:
        """Fe2O3 yields {Fe: 2.0, O: 3.0}."""
        counts = get_formula_counts("Fe2O3")
        assert counts["Fe"] == pytest.approx(2.0)
        assert counts["O"] == pytest.approx(3.0)

    def test_single_element_count_is_one(self) -> None:
        """Pure Al yields {Al: 1.0}."""
        counts = get_formula_counts("Al")
        assert counts["Al"] == pytest.approx(1.0)

    def test_invalid_returns_empty_dict(self) -> None:
        """Garbage returns empty dict."""
        assert get_formula_counts("!!!") == {}


class TestIsSingleElementFormula:
    """is_single_element_formula() checks purity of formula."""

    def test_pure_iron_matches_fe(self) -> None:
        """'Fe' matches element 'Fe'."""
        assert is_single_element_formula("Fe", "Fe") is True

    def test_oxide_does_not_match_single_element(self) -> None:
        """Fe2O3 is not a single-element formula for Fe."""
        assert is_single_element_formula("Fe2O3", "Fe") is False

    def test_wrong_element(self) -> None:
        """'Fe' does not match element 'Cr'."""
        assert is_single_element_formula("Fe", "Cr") is False


class TestFmtStable:
    """_fmt_stable() formats stability flags as TRUE/FALSE/empty."""

    @pytest.mark.parametrize("val", [True, "true", "1", "yes"])
    def test_truthy_becomes_TRUE(self, val) -> None:
        """Truthy values format as 'TRUE'."""
        assert _fmt_stable(val) == "TRUE"

    @pytest.mark.parametrize("val", [False, "false", "0", "no"])
    def test_falsy_becomes_FALSE(self, val) -> None:
        """Falsy values format as 'FALSE'."""
        assert _fmt_stable(val) == "FALSE"

    def test_nan_becomes_empty_string(self) -> None:
        """NaN formats as empty string."""
        assert _fmt_stable(float("nan")) == ""


class TestRowEls:
    """_row_els() returns a frozenset of formula elements."""

    def test_binary_oxide(self) -> None:
        """Fe2O3 yields frozenset({Fe, O})."""
        assert _row_els("Fe2O3") == frozenset({"Fe", "O"})

    def test_ternary(self) -> None:
        """CrHfTa yields frozenset of three elements."""
        assert _row_els("CrHfTa") == frozenset({"Cr", "Hf", "Ta"})

    def test_invalid_returns_empty_frozenset(self) -> None:
        """Garbage returns empty frozenset."""
        assert _row_els("!!!") == frozenset()
