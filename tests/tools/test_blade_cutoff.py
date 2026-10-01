"""Tests for BladeCutoff neighbor-shell distance computation."""

from __future__ import annotations

import numpy as np
import pytest

from blade.tools.blade_cutoff import BladeCutoff


@pytest.fixture
def cutoff() -> BladeCutoff:
    """A default BladeCutoff instance."""
    return BladeCutoff()


@pytest.fixture
def cubic_lattice(cutoff: BladeCutoff) -> np.ndarray:
    """3×3 lattice matrix for a cubic unit cell with a = 3.0 Å."""
    return cutoff.lattice_from_params(3.0, 3.0, 3.0, 90.0, 90.0, 90.0)


@pytest.fixture
def hedb1_lattice(cutoff: BladeCutoff) -> np.ndarray:
    """3×3 lattice matrix for the hexagonal HEDB1 prototype."""
    return cutoff.lattice_from_params(3.0624, 3.0624, 3.2392, 90.0, 90.0, 120.0)


class TestLatticeFromParamsCubic:
    """lattice_from_params() for cubic symmetry (a=b=c, all angles 90°)."""

    def test_output_shape_is_3x3(self, cubic_lattice: np.ndarray) -> None:
        """Returns a (3, 3) array for any valid input."""
        assert cubic_lattice.shape == (3, 3)

    def test_cubic_diagonal_is_a(self, cutoff: BladeCutoff) -> None:
        """Diagonal elements equal a for a cubic cell."""
        a = 4.05
        lat = cutoff.lattice_from_params(a, a, a, 90.0, 90.0, 90.0)
        assert lat[0, 0] == pytest.approx(a, abs=1e-10)
        assert lat[1, 1] == pytest.approx(a, abs=1e-10)
        assert lat[2, 2] == pytest.approx(a, abs=1e-10)

    def test_cubic_off_diagonal_near_zero(self, cutoff: BladeCutoff) -> None:
        """Off-diagonal elements are negligible for a cubic cell."""
        a = 3.5
        lat = cutoff.lattice_from_params(a, a, a, 90.0, 90.0, 90.0)
        assert lat[0, 1] == pytest.approx(0.0, abs=1e-10)
        assert lat[0, 2] == pytest.approx(0.0, abs=1e-10)
        assert lat[1, 2] == pytest.approx(0.0, abs=1e-10)

    def test_first_row_is_a_zero_zero(self, cubic_lattice: np.ndarray) -> None:
        """First lattice vector lies along x: [a, 0, 0]."""
        assert cubic_lattice[0, 1] == pytest.approx(0.0, abs=1e-10)
        assert cubic_lattice[0, 2] == pytest.approx(0.0, abs=1e-10)


class TestLatticeFromParamsHexagonal:
    """lattice_from_params() for hexagonal symmetry (gamma=120°)."""

    def test_hexagonal_a_component(self, cutoff: BladeCutoff) -> None:
        """First vector magnitude equals a for the hexagonal cell."""
        a = 3.0624
        lat = cutoff.lattice_from_params(a, a, 3.2392, 90.0, 90.0, 120.0)
        assert lat[0, 0] == pytest.approx(a, abs=1e-10)

    def test_hexagonal_b_y_positive(self, hedb1_lattice: np.ndarray) -> None:
        """The b-vector y-component is positive for gamma=120°."""
        assert hedb1_lattice[1, 1] > 0.0

    def test_hexagonal_c_z_positive(self, hedb1_lattice: np.ndarray) -> None:
        """The c-vector z-component is positive."""
        assert hedb1_lattice[2, 2] > 0.0


class TestReadCoords:
    """read_coords() parses ATAT fractional coordinate strings."""

    def test_two_atom_basis_shape(self, cutoff: BladeCutoff) -> None:
        """A two-line coordinate string produces shape (2, 3)."""
        coords = "0.0 0.0 0.0 a\n0.333333 0.666667 0.5 B"
        frac = cutoff.read_coords(coords)
        assert frac.shape == (2, 3)

    def test_three_atom_basis_shape(self, cutoff: BladeCutoff) -> None:
        """A three-line coordinate string produces shape (3, 3)."""
        coords = "0.0 0.0 0.0 a\n0.333333 0.666667 0.5 B\n0.666667 0.333333 0.5 B"
        frac = cutoff.read_coords(coords)
        assert frac.shape == (3, 3)

    def test_extra_tokens_ignored(self, cutoff: BladeCutoff) -> None:
        """Labels after the three fractional coordinates are silently ignored."""
        coords = "0.1 0.2 0.3 sublattice_label extra_token"
        frac = cutoff.read_coords(coords)
        assert frac[0, 0] == pytest.approx(0.1, abs=1e-10)
        assert frac[0, 1] == pytest.approx(0.2, abs=1e-10)
        assert frac[0, 2] == pytest.approx(0.3, abs=1e-10)

    def test_first_atom_at_origin(self, cutoff: BladeCutoff) -> None:
        """Origin atom fractional coordinates are parsed as (0, 0, 0)."""
        coords = "0.000000 0.000000 0.000000 a\n0.333333 0.666667 0.500000 B"
        frac = cutoff.read_coords(coords)
        assert frac[0] == pytest.approx([0.0, 0.0, 0.0], abs=1e-6)


class TestGetShells:
    """get_shells() returns sorted unique neighbor-shell distances."""

    def test_first_shell_is_positive(self, cutoff: BladeCutoff, cubic_lattice: np.ndarray) -> None:
        """First-neighbor shell distance is positive."""
        frac = cutoff.read_coords("0.0 0.0 0.0")
        shells = cutoff.get_shells(cubic_lattice, frac, rep=(3, 3, 3))
        assert shells[0] > 0.0

    def test_cubic_first_shell_equals_a(self, cutoff: BladeCutoff) -> None:
        """First shell equals lattice constant a for a simple cubic cell."""
        a = 3.0
        lat = cutoff.lattice_from_params(a, a, a, 90.0, 90.0, 90.0)
        frac = cutoff.read_coords("0.0 0.0 0.0")
        shells = cutoff.get_shells(lat, frac, rep=(3, 3, 3))
        assert shells[0] == pytest.approx(a, abs=1e-4)

    def test_shells_are_sorted_ascending(self, cutoff: BladeCutoff, cubic_lattice: np.ndarray) -> None:
        """Returned shell distances are in non-decreasing order."""
        frac = cutoff.read_coords("0.0 0.0 0.0")
        shells = cutoff.get_shells(cubic_lattice, frac, rep=(3, 3, 3))
        assert all(shells[i] <= shells[i + 1] for i in range(len(shells) - 1))

    def test_hexagonal_first_shell_positive(self, cutoff: BladeCutoff, hedb1_lattice: np.ndarray) -> None:
        """First neighbor shell distance is positive for the HEDB1 hexagonal lattice."""
        frac = cutoff.read_coords("0.000000 0.000000 0.000000 a\n0.333333 0.666667 0.500000 B\n0.666667 0.333333 0.500000 B")
        shells = cutoff.get_shells(hedb1_lattice, frac, rep=(4, 3, 2))
        assert shells[0] > 0.0

    def test_multiple_shells_returned(self, cutoff: BladeCutoff, cubic_lattice: np.ndarray) -> None:
        """At least two distinct neighbor shells are identified for a cubic cell."""
        frac = cutoff.read_coords("0.0 0.0 0.0")
        shells = cutoff.get_shells(cubic_lattice, frac, rep=(3, 3, 3))
        assert len(shells) >= 2
