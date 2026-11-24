"""Unit tests for blade.analysis._sublattice — pure-math helpers only, no real data."""

from __future__ import annotations

import numpy as np

from blade.analysis._stability_tables import _cef_enthalpy, _ideal_mixing_2subl
from blade.analysis._sublattice import (
    cross_endmember_fracs,
    detect_sublattice_metals,
    parse_subl_fracs,
    phase_stoich_from_mults,
    product_simplex_grid,
)

# ---------------------------------------------------------------------------
# parse_subl_fracs
# ---------------------------------------------------------------------------


def test_parse_subl_fracs_max():
    result = parse_subl_fracs("sqs_lev=1_a_Ti=0.5,a_Nb=0.5,b_Al=1.0,c_C=1.0")
    assert result == {"a": {"Ti": 0.5, "Nb": 0.5}, "b": {"Al": 1.0}, "c": {"C": 1.0}}


def test_parse_subl_fracs_empty():
    assert parse_subl_fracs("") == {}


def test_parse_subl_fracs_single():
    result = parse_subl_fracs("sqs_lev=0_a_Ti=1.0,b_Al=1.0")
    assert result == {"a": {"Ti": 1.0}, "b": {"Al": 1.0}}


# ---------------------------------------------------------------------------
# product_simplex_grid
# ---------------------------------------------------------------------------


def test_product_simplex_grid_shape_and_sum():
    subl_metals = [["Ti", "Nb"], ["Al"]]
    subl_mults = [2, 1]
    y_overall, (y_a, y_b) = product_simplex_grid(subl_metals, subl_mults, y_step=0.5)

    # Binary Ti-Nb at step 0.5: [1,0], [0.5,0.5], [0,1] → 3 points
    # Pure Al: [1] → 1 point
    assert y_a.shape == (3, 2)
    assert y_b.shape == (1, 1)
    assert y_overall.shape == (3, 3)  # 3 total metals (Ti, Nb, Al)

    # Each row of y_overall sums to 1
    np.testing.assert_allclose(y_overall.sum(axis=1), np.ones(3), atol=1e-12)


def test_product_simplex_grid_values():
    subl_metals = [["Ti", "Nb"], ["Al"]]
    subl_mults = [2, 1]
    y_overall, _ = product_simplex_grid(subl_metals, subl_mults, y_step=0.5)

    # Each row sums to 1 (already tested above); check each row is achievable
    # The 3 rows correspond to the 3 simplex points on the Ti-Nb binary; each paired with
    # pure Al.  Columns are [Ti, Nb, Al].  simplex_grid_nd ordering is ascending y[0],
    # so rows are [0,1,…], [0.5,0.5,…], [1,0,…] for the y_a_grid.
    # Row 0: Ti=0, Nb=2/3, Al=1/3
    # Row 1: Ti=1/3, Nb=1/3, Al=1/3
    # Row 2: Ti=2/3, Nb=0, Al=1/3
    expected = np.array([[0, 2 / 3, 1 / 3], [1 / 3, 1 / 3, 1 / 3], [2 / 3, 0, 1 / 3]])
    np.testing.assert_allclose(y_overall, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# phase_stoich_from_mults
# ---------------------------------------------------------------------------


def test_phase_stoich_max():
    assert abs(phase_stoich_from_mults([2, 1], 1) - 1 / 3) < 1e-12


def test_phase_stoich_diboride():
    assert abs(phase_stoich_from_mults([1], 2) - 2.0) < 1e-12


# ---------------------------------------------------------------------------
# _cef_enthalpy
# ---------------------------------------------------------------------------


def test_cef_enthalpy_2x1():
    # h_ref[0,0] = 3.0 (Ti-Al endpoint), h_ref[1,0] = 5.0 (Nb-Al endpoint)
    h_ref = np.array([[3.0], [5.0]])  # shape (2, 1)
    y_a = np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])  # 3 points
    y_b = np.array([[1.0]])  # 1 point (pure Al)

    h0 = _cef_enthalpy(y_a, y_b, h_ref)
    # H[k] = y_a_Ti*h[0,0] + y_a_Nb*h[1,0]
    # = [1*3+0*5, 0.5*3+0.5*5, 0*3+1*5] = [3, 4, 5]
    np.testing.assert_allclose(h0, [3.0, 4.0, 5.0], atol=1e-12)


def test_cef_enthalpy_2x2():
    h_ref = np.array([[2.0, 4.0], [6.0, 8.0]])  # shape (2, 2)
    y_a = np.array([[1.0, 0.0], [0.0, 1.0]])  # pure Ti, pure Nb
    y_b = np.array([[1.0, 0.0], [0.0, 1.0]])  # pure Al, pure Si

    h0 = _cef_enthalpy(y_a, y_b, h_ref)
    # 2×2 product: Ti-Al=2, Ti-Si=4, Nb-Al=6, Nb-Si=8
    np.testing.assert_allclose(h0, [2.0, 4.0, 6.0, 8.0], atol=1e-12)


# ---------------------------------------------------------------------------
# _ideal_mixing_2subl
# ---------------------------------------------------------------------------


def test_ideal_mixing_2subl_pure_endpoints_zero():
    # Pure endpoints have no mixing entropy
    y_a = np.array([[1.0, 0.0], [0.0, 1.0]])
    y_b = np.array([[1.0]])
    s_mix = _ideal_mixing_2subl(y_a, y_b, [2, 1])
    # Both rows: S_a = 0 (pure), S_b = 0 (pure Al) → 0
    np.testing.assert_allclose(s_mix, [0.0, 0.0], atol=1e-12)


def test_ideal_mixing_2subl_binary_shape():
    y_a = np.array([[0.5, 0.5]])
    y_b = np.array([[1.0]])
    s_mix = _ideal_mixing_2subl(y_a, y_b, [2, 1])
    # S_a = 0.5*ln(0.5) + 0.5*ln(0.5) = ln(0.5)
    # S_b = 0 (pure)
    # s_mix = (2*ln(0.5) + 1*0) / 3
    expected = 2 * np.log(0.5) / 3
    np.testing.assert_allclose(s_mix[0], expected, atol=1e-12)


# ---------------------------------------------------------------------------
# detect_sublattice_metals
# ---------------------------------------------------------------------------


def test_detect_sublattice_metals(tmp_path):
    # Create mock sqs_lev= dirs
    (tmp_path / "sqs_lev=0_a_Ti=1.0,b_Al=1.0,c_C=1.0").mkdir()
    (tmp_path / "sqs_lev=1_a_Ti=0.5,a_Nb=0.5,b_Al=1.0,c_C=1.0").mkdir()
    (tmp_path / "not_a_sqs_dir").mkdir()

    subl, is_multi = detect_sublattice_metals(tmp_path, phase_element="C")

    assert is_multi is True
    assert set(subl.keys()) == {"a", "b"}
    assert set(subl["a"]) == {"Ti", "Nb"}
    assert subl["b"] == ["Al"]


def test_detect_sublattice_metals_single_subl(tmp_path):
    (tmp_path / "sqs_lev=0_a_Ti=1.0,b_B=1.0").mkdir()
    (tmp_path / "sqs_lev=1_a_Ti=0.5,a_Nb=0.5,b_B=1.0").mkdir()

    subl, is_multi = detect_sublattice_metals(tmp_path, phase_element="B")

    # B is excluded; only 'a' sublattice has variable metals
    assert is_multi is False
    assert set(subl["a"]) == {"Ti", "Nb"}
    assert "b" not in subl


# ---------------------------------------------------------------------------
# cross_endmember_fracs
# ---------------------------------------------------------------------------


def test_cross_endmember_fracs():
    endpoints = cross_endmember_fracs([["Ti", "Nb"], ["Al", "Si"]], [2, 1])
    assert len(endpoints) == 4  # 2×2 cross product
    # Ti-Al: Ti=2/3, Al=1/3
    ti_al = next(e for e in endpoints if abs(e.get("Ti", 0) - 2 / 3) < 1e-9)
    assert abs(ti_al["Al"] - 1 / 3) < 1e-9


# ---------------------------------------------------------------------------
# BUG 12 — formula-prefixed sqs_lev= dir names
# ---------------------------------------------------------------------------


def test_detect_sublattice_metals_formula_prefix(tmp_path):
    # Database stage copies dirs as "<formula>_sqs_lev=..." — must still parse.
    (tmp_path / "Ti2AlC_sqs_lev=0_a_Ti=1.0,b_Al=1.0,c_C=1.0").mkdir()
    (tmp_path / "TiNbAlC_sqs_lev=1_a_Ti=0.5,a_Nb=0.5,b_Al=1.0,c_C=1.0").mkdir()

    subl, is_multi = detect_sublattice_metals(tmp_path, phase_element="C")

    assert is_multi is True
    assert set(subl["a"]) == {"Ti", "Nb"}
    assert subl["b"] == ["Al"]


def test_detect_sublattice_metals_no_prefix_unchanged(tmp_path):
    # Original format (no formula prefix) must still work after the fix.
    (tmp_path / "sqs_lev=0_a_Ti=1.0,b_Al=1.0,c_C=1.0").mkdir()
    (tmp_path / "sqs_lev=1_a_Ti=0.5,a_Nb=0.5,b_Al=1.0,c_C=1.0").mkdir()

    subl, is_multi = detect_sublattice_metals(tmp_path, phase_element="C")

    assert is_multi is True
    assert set(subl["a"]) == {"Ti", "Nb"}
    assert subl["b"] == ["Al"]


# ---------------------------------------------------------------------------
# BUG 13 — _cef_enthalpy shape coverage (1+2 and asymmetric)
# ---------------------------------------------------------------------------


def test_cef_enthalpy_1x2():
    # 1 M-metal, 2 A-metals
    h_ref = np.array([[2.0, 4.0]])  # shape (1, 2): Ti-Al=2, Ti-Si=4
    y_a = np.array([[1.0]])  # 1 pt (pure Ti)
    y_b = np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])  # 3 pts

    h0 = _cef_enthalpy(y_a, y_b, h_ref)
    # H[k] = y_b_Al*h[0,0] + y_b_Si*h[0,1]
    np.testing.assert_allclose(h0, [2.0, 3.0, 4.0], atol=1e-12)
    assert h0.shape == (3,)


def test_cef_enthalpy_2x1_shape():
    # Existing 2+1 coverage (matches test_cef_enthalpy_2x1 above but explicit shape check)
    h_ref = np.array([[3.0], [5.0]])  # (2,1)
    y_a = np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])
    y_b = np.array([[1.0]])
    h0 = _cef_enthalpy(y_a, y_b, h_ref)
    assert h0.shape == (3,)
    np.testing.assert_allclose(h0, [3.0, 4.0, 5.0], atol=1e-12)
