"""Neighbor-shell cutoff distance computation for ATAT mcsqs.

This module provides :class:`BladeCutoff`, which converts lattice parameters
into a 3×3 lattice matrix, reads fractional coordinates from an ATAT
coordinate string, and computes neighbor-shell distances via a supercell
construction. The resulting shell distances are used to set the ``-2``,
``-3``, and ``-4`` cutoff arguments for ``corrdump`` and ``mcsqs``.

Example::

    from blade.tools.blade_cutoff import BladeCutoff

    cutoff = BladeCutoff()
    lattice = cutoff.lattice_from_params(3.3, 3.3, 5.2, 90, 90, 120)
    coords = "0.0 0.0 0.0\\n0.333333 0.666667 0.5"
    frac = cutoff.read_coords(coords)
    shells = cutoff.get_shells(lattice, frac, rep=(4, 3, 2))
    print(f"1NN = {shells[0]:.4f}, 2NN = {shells[1]:.4f}")
"""

from __future__ import annotations

import re

import numpy as np

__author__ = "Chase Katz"


def _build_supercell(frac: np.ndarray, rep: tuple[int, int, int]) -> np.ndarray:
    """Tile fractional coordinates into a supercell.

    Args:
        frac (np.ndarray): Fractional coordinates, shape ``(N, 3)``.
        rep (tuple[int, int, int]): Repetitions along each lattice vector.

    Returns:
        np.ndarray: Supercell fractional coordinates, shape ``(N * nx*ny*nz, 3)``.
    """
    nx, ny, nz = rep
    shifts = np.array(
        [[i, j, k] for i in range(nx) for j in range(ny) for k in range(nz)],
        dtype=float,
    )
    return (frac[None, :, :] + shifts[:, None, :]).reshape(-1, 3)


def _min_image(df: np.ndarray, rep: tuple[int, int, int]) -> np.ndarray:
    """Apply minimum-image convention in fractional coordinates.

    Args:
        df (np.ndarray): Fractional displacement vectors.
        rep (tuple[int, int, int]): Supercell repetitions along each axis.

    Returns:
        np.ndarray: Wrapped displacement vectors.
    """
    rep_arr = np.array(rep, dtype=float)
    return df - rep_arr * np.round(df / rep_arr)


class BladeCutoff:
    """Compute neighbor-shell distances for ATAT cutoff parameters.

    Used internally by :class:`~blade.tools.blade_sqsgen.BladeSQS` to derive
    the ``-2``, ``-3``, and ``-4`` cutoffs passed to ``corrdump`` and
    ``mcsqs``.
    """

    def __init__(self) -> None:
        """Initialize BladeCutoff."""

    def lattice_from_params(
        self,
        a: float,
        b: float,
        c: float,
        alpha: float,
        beta: float,
        gamma: float,
    ) -> np.ndarray:
        """Convert lattice parameters to a 3×3 lattice matrix in Ångströms.

        Uses the standard crystallographic convention where the **a** vector
        lies along x, **b** is in the xy-plane, and **c** has a general
        orientation.

        Args:
            a (float): Lattice parameter *a* in Ångströms.
            b (float): Lattice parameter *b* in Ångströms.
            c (float): Lattice parameter *c* in Ångströms.
            alpha (float): Angle between **b** and **c** in degrees.
            beta (float): Angle between **a** and **c** in degrees.
            gamma (float): Angle between **a** and **b** in degrees.

        Returns:
            np.ndarray: 3×3 lattice matrix where rows are lattice vectors,
            shape ``(3, 3)``.
        """
        alpha_r = np.radians(alpha)
        beta_r = np.radians(beta)
        gamma_r = np.radians(gamma)

        ax = a
        bx = b * np.cos(gamma_r)
        by = b * np.sin(gamma_r)
        cx = c * np.cos(beta_r)
        cy = c * (np.cos(alpha_r) - np.cos(beta_r) * np.cos(gamma_r)) / np.sin(gamma_r)
        cz = np.sqrt(c**2 - cx**2 - cy**2)

        return np.array([
            [ax, 0.0, 0.0],
            [bx, by, 0.0],
            [cx, cy, cz],
        ])

    def read_coords(self, coord_string: str) -> np.ndarray:
        """Parse an ATAT fractional coordinate string into an array.

        Each line of ``coord_string`` must begin with three whitespace-separated
        floats representing fractional coordinates. Additional tokens (e.g.,
        sublattice labels) on each line are ignored.

        Args:
            coord_string (str): Multiline string of fractional coordinates,
                one atom per line (e.g., the ``coords`` field of a
                ``phases_dict``).

        Returns:
            np.ndarray: Fractional coordinate array, shape ``(N, 3)``.
        """
        frac = []
        for line in coord_string.strip().splitlines():
            vals = re.split(r"\s+", line.strip())
            frac.append([float(vals[0]), float(vals[1]), float(vals[2])])
        return np.array(frac)

    def get_shells(
        self,
        lattice: np.ndarray,
        frac: np.ndarray,
        rep: tuple[int, int, int],
        tol: float = 1e-4,
    ) -> np.ndarray:
        """Compute normalized neighbor-shell distances for a supercell.

        Builds a supercell from ``frac`` and ``rep``, computes all pairwise
        distances under the minimum-image convention, and groups them into
        shells. Shell distances are normalized by the first-neighbor distance
        so they can be used as dimensionless cutoff ratios.

        Args:
            lattice (np.ndarray): 3×3 lattice matrix (rows = lattice vectors).
            frac (np.ndarray): Fractional coordinates, shape ``(N, 3)``.
            rep (tuple[int, int, int]): Supercell repetitions ``(nx, ny, nz)``.
            tol (float, optional): Distance tolerance for shell grouping in
                Ångströms. Defaults to ``1e-4``.

        Returns:
            np.ndarray: Sorted array of normalized shell distances, starting
            from 1.0 (first-neighbor shell).
        """
        frac_sc = _build_supercell(frac, rep)
        n = len(frac_sc)
        dists: list[float] = []

        for i in range(n):
            df = frac_sc[i + 1:] - frac_sc[i]
            df = _min_image(df, rep)
            cart = df @ lattice
            r = np.linalg.norm(cart, axis=1)
            dists.extend(float(val) for val in r if val > 1e-6)

        dists_sorted = np.sort(dists)

        shells: list[float] = []
        for r in dists_sorted:
            if not shells or abs(r - shells[-1]) > tol:
                shells.append(r)

        shells_arr = np.array(shells)
        return shells_arr / shells_arr[0]
