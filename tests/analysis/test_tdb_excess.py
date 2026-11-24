"""Tests for _load_tdb_excess and use_tdb_excess path in _build_phase_data_2subl."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from blade.analysis._stability_tables import _build_phase_data_2subl, _load_tdb_excess

# ---------------------------------------------------------------------------
# Minimal TDB fixtures
# ---------------------------------------------------------------------------

_TDB_2PLUS1 = """\
PHASE MAX1_3 % 3 2 1 1  !
CONSTITUENT MAX1_3 :NB,TI:AL:C:!
PARAMETER G(MAX1_3,TI:AL:C;0) 298.15 -100000.0 ; 10000.0 N !
PARAMETER G(MAX1_3,NB:AL:C;0) 298.15 -200000.0 ; 10000.0 N !
PARAMETER L(MAX1_3,NB,TI:AL:C;0) 298.15 -96485.0 ; 10000.0 N !
PARAMETER L(MAX1_3,NB,TI:AL:C;1) 298.15 +96485.0 ; 10000.0 N !
"""

_TDB_1PLUS2 = """\
PHASE MAX1_3 % 3 2 1 1  !
CONSTITUENT MAX1_3 :TI:AL,SI:C:!
PARAMETER G(MAX1_3,TI:AL:C;0) 298.15 -100000.0 ; 10000.0 N !
PARAMETER G(MAX1_3,TI:SI:C;0) 298.15 -150000.0 ; 10000.0 N !
PARAMETER L(MAX1_3,TI:AL,SI:C;0) 298.15 +96485.0 ; 10000.0 N !
"""

_TDB_2PLUS2 = """\
PHASE MAX1_4 % 3 2 1 1  !
CONSTITUENT MAX1_4 :NB,TI:AL,SI:C:!
PARAMETER G(MAX1_4,TI:AL:C;0) 298.15 -100000.0 ; 10000.0 N !
PARAMETER G(MAX1_4,NB:AL:C;0) 298.15 -200000.0 ; 10000.0 N !
PARAMETER G(MAX1_4,TI:SI:C;0) 298.15 -150000.0 ; 10000.0 N !
PARAMETER G(MAX1_4,NB:SI:C;0) 298.15 -250000.0 ; 10000.0 N !
PARAMETER L(MAX1_4,NB,TI:AL:C;0) 298.15 -96485.0 ; 10000.0 N !
PARAMETER L(MAX1_4,TI:AL,SI:C;0) 298.15 +96485.0 ; 10000.0 N !
PARAMETER L(MAX1_4,NB:AL,SI:C;0) 298.15 -96485.0 ; 10000.0 N !
PARAMETER L(MAX1_4,NB,TI:AL,SI:C;0) 298.15 +96485.0 ; 10000.0 N !
"""

_TDB_NO_L = """\
PHASE MAX1_2 % 3 2 1 1  !
CONSTITUENT MAX1_2 :NB:AL:C:!
PARAMETER G(MAX1_2,NB:AL:C;0) 298.15 -300000.0 ; 10000.0 N !
"""


@pytest.fixture()
def tdb_2p1(tmp_path):
    p = tmp_path / "MAX1_3.tdb"
    p.write_text(_TDB_2PLUS1)
    return p


@pytest.fixture()
def tdb_1p2(tmp_path):
    p = tmp_path / "MAX1_3.tdb"
    p.write_text(_TDB_1PLUS2)
    return p


@pytest.fixture()
def tdb_2p2(tmp_path):
    p = tmp_path / "MAX1_4.tdb"
    p.write_text(_TDB_2PLUS2)
    return p


@pytest.fixture()
def tdb_no_l(tmp_path):
    p = tmp_path / "MAX1_2.tdb"
    p.write_text(_TDB_NO_L)
    return p


# ---------------------------------------------------------------------------
# _load_tdb_excess unit tests
# ---------------------------------------------------------------------------

def test_no_l_params_returns_zeros(tdb_no_l):
    subl_metals = [["Nb"], ["Al"]]
    subl_mults = [2, 1]
    y_a = np.ones((1, 1))
    y_b = np.ones((1, 1))
    exc = _load_tdb_excess(tdb_no_l, subl_metals, subl_mults, y_a, y_b)
    assert exc.shape == (1,)
    np.testing.assert_allclose(exc, 0.0)


def test_missing_file_returns_zeros():
    subl_metals = [["Nb", "Ti"], ["Al"]]
    subl_mults = [2, 1]
    y_a = np.array([[0.5, 0.5]])
    y_b = np.ones((1, 1))
    exc = _load_tdb_excess(Path("/nonexistent/path.tdb"), subl_metals, subl_mults, y_a, y_b)
    np.testing.assert_allclose(exc, 0.0)


def test_msite_pair_l0_value(tdb_2p1):
    """L(NB,TI:AL;0) = -96485 J/mol → -1/3 eV/metal at equimolar NB=TI=0.5."""
    subl_metals = [["Nb", "Ti"], ["Al"]]
    subl_mults = [2, 1]
    # 1 a-point: Nb=0.5, Ti=0.5 → y_a_NB * y_a_TI * (y_a_NB - y_a_TI)^0 = 0.25
    # 1 b-point: Al=1.0 → y_b_Al = 1.0
    # L0 = -96485 J/mol → -1.0 eV/fu → /3 metal = -0.3333 eV/metal
    # G_exc = 0.25 * 1.0 * (-0.3333) = -0.08333 eV/metal
    y_a = np.array([[0.5, 0.5]])
    y_b = np.ones((1, 1))
    exc = _load_tdb_excess(tdb_2p1, subl_metals, subl_mults, y_a, y_b)
    # L0 contribution only at order 0: 0.5 * 0.5 * 1 * 1.0 * (-1/3) = -1/12
    # L1: 0.5 * 0.5 * (0.5 - 0.5) * 1.0 * (1/3) = 0
    expected = 0.5 * 0.5 * 1.0 * (-96485.0 / 96485.0 / 3)
    np.testing.assert_allclose(exc[0], expected, rtol=1e-6)


def test_msite_pair_l1_antisymmetry(tdb_2p1):
    """L1 term changes sign when Nb and Ti site fractions are swapped."""
    subl_metals = [["Nb", "Ti"], ["Al"]]
    subl_mults = [2, 1]
    y_b = np.ones((1, 1))
    y_a_nb_rich = np.array([[0.75, 0.25]])
    y_a_ti_rich = np.array([[0.25, 0.75]])
    exc_nb = _load_tdb_excess(tdb_2p1, subl_metals, subl_mults, y_a_nb_rich, y_b)
    exc_ti = _load_tdb_excess(tdb_2p1, subl_metals, subl_mults, y_a_ti_rich, y_b)
    # L0 is symmetric in y_i*y_j; L1 ~ y_i*y_j*(y_i-y_j) flips sign
    # So total excess won't be equal; the L1 contribution should flip
    # L1(Nb=0.75,Ti=0.25) → L1 term = 0.75*0.25*(0.75-0.25) = +0.09375 * (+1/3)
    # L1(Nb=0.25,Ti=0.75) → L1 term = 0.25*0.75*(0.25-0.75) = -0.09375 * (+1/3)
    l1_nb = 0.75 * 0.25 * (0.75 - 0.25) * (96485.0 / 96485.0 / 3)
    l1_ti = 0.25 * 0.75 * (0.25 - 0.75) * (96485.0 / 96485.0 / 3)
    assert l1_nb > 0 and l1_ti < 0
    # Difference in total excess should be 2 * l1_contribution
    np.testing.assert_allclose(exc_nb[0] - exc_ti[0], l1_nb - l1_ti, rtol=1e-6)


def test_asite_pair_l0_value(tdb_1p2):
    """L(TI:AL,SI;0) = +96485 J/mol → +1/3 eV/metal at equimolar AL=SI=0.5."""
    subl_metals = [["Ti"], ["Al", "Si"]]
    subl_mults = [2, 1]
    y_a = np.ones((1, 1))
    y_b = np.array([[0.5, 0.5]])
    exc = _load_tdb_excess(tdb_1p2, subl_metals, subl_mults, y_a, y_b)
    # y_a_Ti=1, y_b_Al=0.5, y_b_Si=0.5, (y_b_Al - y_b_Si)^0 = 1
    # G_exc = 1.0 * 0.5 * 0.5 * 1.0 * (1/3) = 1/12
    expected = 1.0 * 0.5 * 0.5 * (96485.0 / 96485.0 / 3)
    np.testing.assert_allclose(exc[0], expected, rtol=1e-6)


def test_2plus2_all_term_types(tdb_2p2):
    """2+2 system: M-site pair, A-site pair, and reciprocal all accumulate."""
    subl_metals = [["Nb", "Ti"], ["Al", "Si"]]
    subl_mults = [2, 1]
    y_a = np.array([[0.5, 0.5]])
    y_b = np.array([[0.5, 0.5]])
    exc = _load_tdb_excess(tdb_2p2, subl_metals, subl_mults, y_a, y_b)
    assert exc.shape == (1,)
    # Must be non-zero (all three term types contribute)
    assert abs(exc[0]) > 1e-10


def test_output_shape_product_grid(tdb_2p2):
    """Output shape = n_a_pts * n_b_pts for arbitrary grid sizes."""
    subl_metals = [["Nb", "Ti"], ["Al", "Si"]]
    subl_mults = [2, 1]
    y_a = np.array([[1.0, 0.0], [0.75, 0.25], [0.5, 0.5], [0.25, 0.75], [0.0, 1.0]])
    y_b = np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])
    exc = _load_tdb_excess(tdb_2p2, subl_metals, subl_mults, y_a, y_b)
    assert exc.shape == (5 * 3,)


def test_pure_endpoints_zero_m_site_excess(tdb_2p1):
    """At pure-metal endpoints y_i=1, y_j=0: M-site pair excess = 0."""
    subl_metals = [["Nb", "Ti"], ["Al"]]
    subl_mults = [2, 1]
    y_a_nb = np.array([[1.0, 0.0]])
    y_a_ti = np.array([[0.0, 1.0]])
    y_b = np.ones((1, 1))
    exc_nb = _load_tdb_excess(tdb_2p1, subl_metals, subl_mults, y_a_nb, y_b)
    exc_ti = _load_tdb_excess(tdb_2p1, subl_metals, subl_mults, y_a_ti, y_b)
    np.testing.assert_allclose(exc_nb, 0.0, atol=1e-12)
    np.testing.assert_allclose(exc_ti, 0.0, atol=1e-12)


def test_reciprocal_term_separable(tdb_2p2):
    """Reciprocal term is zero at pure endpoints of either sublattice."""
    subl_metals = [["Nb", "Ti"], ["Al", "Si"]]
    subl_mults = [2, 1]
    # Only reciprocal term: set up so M-site and A-site pairs also zero
    # Use Nb=1 (no M-site pair) and Al=Si=0.5 (A-site pair present for Nb only)
    # At Nb=1.0, the M-site pair contribution y_NB * y_TI = 0 → zero
    # Reciprocal: y_NB * y_TI * y_AL * y_SI = 0 → zero
    y_a_pure = np.array([[1.0, 0.0]])  # Nb only
    y_b_mid = np.array([[0.5, 0.5]])
    exc = _load_tdb_excess(tdb_2p2, subl_metals, subl_mults, y_a_pure, y_b_mid)
    # Reciprocal and M-site pair both zero; only A-site pair for Nb survives
    # L(NB:AL,SI;0) = -96485 → -1/3 eV/fu → G = 1.0 * 0.5*0.5 * (-1/3) = -1/12
    expected_asite = 1.0 * 0.5 * 0.5 * (-96485.0 / 96485.0 / 3)
    np.testing.assert_allclose(exc[0], expected_asite, rtol=1e-6)


# ---------------------------------------------------------------------------
# _build_phase_data_2subl integration: use_tdb_excess changes phase_H0
# ---------------------------------------------------------------------------

def _make_sqs_dirs(tmp_path, subl_metals, subl_mults, h_ref_values):
    """Create minimal SQS structure dirs for a 2-sublattice system."""
    import textwrap

    blade_dir = tmp_path / "blade"
    blade_dir.mkdir(parents=True, exist_ok=True)

    metals_a, metals_b = subl_metals
    # Build all cross-product endpoints
    for i, ma in enumerate(metals_a):
        for j, mb in enumerate(metals_b):
            fracs_a = ",".join(f"a_{m}={1.0 if k == i else 0.0}" for k, m in enumerate(metals_a))
            fracs_b = ",".join(f"b_{m}={1.0 if k == j else 0.0}" for k, m in enumerate(metals_b))
            dname = f"sqs_lev=0_{fracs_a},{fracs_b},c_C=1.0"
            d = blade_dir / dname
            d.mkdir()
            (d / "energy").write_text(f"{h_ref_values[i][j] * 3}\n")  # e_total; n_fu = mult_a+mult_b
            # CONTCAR: 2 M + 1 A + 1 C atoms
            n_ma = 2 if subl_mults[0] == 2 else 1
            n_mb = 1
            poscar = textwrap.dedent(f"""\
                dummy
                1.0
                3.0 0.0 0.0
                0.0 3.0 0.0
                0.0 0.0 3.0
                {ma} {mb} C
                {n_ma} {n_mb} 1
                Direct
            """)
            (d / "CONTCAR").write_text(poscar)
    return blade_dir


def test_use_tdb_excess_false_unchanged(tmp_path):
    """use_tdb_excess=False → phase_H0 equals ideal CEF (no excess)."""
    subl_mults = [2, 1]
    # 1 b-metal means no A-site mixing
    subl_metals_1b = [["Nb", "Ti"], ["Al"]]
    h_ref_1b = [[-2.0], [-1.8]]
    blade_dir = _make_sqs_dirs(tmp_path / "s1", subl_metals_1b, subl_mults, h_ref_1b)

    tdb_path = tmp_path / "s1" / "dummy.tdb"
    tdb_path.write_text(_TDB_2PLUS1)

    pd_no_exc = _build_phase_data_2subl(
        blade_dir, subl_metals_1b, subl_mults, "C", 1 / 3, 3, 0.5, None, True,
        use_tdb_excess=False, tdb_path=tdb_path,
    )
    pd_with_exc = _build_phase_data_2subl(
        blade_dir, subl_metals_1b, subl_mults, "C", 1 / 3, 3, 0.5, None, True,
        use_tdb_excess=True, tdb_path=tdb_path,
    )

    # With excess active: phase_H0 must differ (L params are non-zero at mixed pts)
    # At pure endpoints they're equal; at mixed pts they differ.
    h0_no = pd_no_exc["phase_H0"]
    h0_ex = pd_with_exc["phase_H0"]
    assert h0_no.shape == h0_ex.shape
    # g_ref_grid should be identical (always ideal CEF surface)
    np.testing.assert_allclose(pd_no_exc["g_ref_grid"], pd_with_exc["g_ref_grid"])
    # phase_H0 must differ somewhere for mixed compositions
    assert not np.allclose(h0_no, h0_ex), "use_tdb_excess=True should change phase_H0"


def test_use_tdb_excess_no_tdb_path_unchanged(tmp_path):
    """use_tdb_excess=True with tdb_path=None → phase_H0 unchanged (graceful skip)."""
    subl_metals = [["Nb", "Ti"], ["Al"]]
    subl_mults = [2, 1]
    h_ref = [[-2.0], [-1.8]]
    blade_dir = _make_sqs_dirs(tmp_path / "s2", subl_metals, subl_mults, h_ref)

    pd_no_path = _build_phase_data_2subl(
        blade_dir, subl_metals, subl_mults, "C", 1 / 3, 3, 0.5, None, True,
        use_tdb_excess=True, tdb_path=None,
    )
    pd_no_exc = _build_phase_data_2subl(
        blade_dir, subl_metals, subl_mults, "C", 1 / 3, 3, 0.5, None, True,
        use_tdb_excess=False, tdb_path=None,
    )
    np.testing.assert_allclose(pd_no_path["phase_H0"], pd_no_exc["phase_H0"])


def test_real_tdb_alti_nb(tmp_path):
    """Smoke test against real AlNbTi MAX1_3 TDB from test data."""
    tdb_real = Path(
        "/Users/chasekatz/Desktop/School/Research/BLADE/Files/Comps_MAX_test"
        "/AlNbTi/MAX1_3/MAX1_3.tdb"
    )
    if not tdb_real.exists():
        pytest.skip("test data not present")

    subl_metals = [["Nb", "Ti"], ["Al"]]
    subl_mults = [2, 1]
    y_a = np.array([[0.5, 0.5]])
    y_b = np.ones((1, 1))

    exc = _load_tdb_excess(tdb_real, subl_metals, subl_mults, y_a, y_b)
    assert exc.shape == (1,)
    # L0(NB,TI:AL;0)=-16029.1 J/mol, L1=-2701.6 J/mol
    # At Nb=Ti=0.5: L0_term = 0.25*(−16029.1/96485/3) = -0.01389 eV/metal-atom
    # L1_term = 0.25*(0.5-0.5)*... = 0
    expected = 0.25 * (-16029.1 / 96485.0 / 3)
    np.testing.assert_allclose(exc[0], expected, rtol=1e-4)
