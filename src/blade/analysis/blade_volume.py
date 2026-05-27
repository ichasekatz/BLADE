"""Volume and lattice parameter extraction from POSCAR files.

This module provides :class:`BLADEVolume`, which scans a composition
directory for ``POSCAR`` files, parses lattice parameters and atom counts,
and returns a :class:`pandas.DataFrame` suitable for downstream analysis.
SQS level and fractional composition metadata are extracted automatically
from the directory-name convention used by BLADE (``sqs_lev=N_a_X=...``).

Example::

    from blade.analysis.blade_volume import BLADEVolume
    from pathlib import Path

    vol = BLADEVolume()
    df = vol.scan_poscars(Path("CrHfTa"))
    print(df[["phase_folder", "sqs_level", "volume_per_atom_A3"]].head())
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass

__author__ = "Chase Katz"


class BLADEVolume:
    """Extract volume and lattice parameters from BLADE POSCAR trees.

    Traverses a per-composition directory that follows the BLADE directory
    convention (``<phase>/<sqs_lev=N_...>/POSCAR``) and collects lattice
    geometry, atom counts, and SQS metadata into a single DataFrame.

    Attributes:
        data (pd.DataFrame | None): Populated after calling
            :meth:`scan_poscars`.
    """

    def __init__(self) -> None:
        """Initialize BLADEVolume."""
        self.data: pd.DataFrame | None = None

    def parse_sqs_meta(self, poscar_path: Path) -> tuple[int | None, dict[str, float]]:
        """Extract SQS level and fractional composition from a POSCAR path.

        Walks up the parent directories looking for a folder name matching
        the pattern ``sqs_lev=<N>`` optionally followed by
        ``_a_<element>=<fraction>`` tokens.

        Args:
            poscar_path (Path): Path to a ``POSCAR`` file inside a BLADE
                directory tree.

        Returns:
            tuple[int | None, dict[str, float]]: A two-element tuple:

            - ``sqs_level``: Integer SQS level, or ``None`` if not found.
            - ``a_fracs``: Dict mapping element symbol to its fractional
              composition on the *a* sublattice (e.g., ``{"Cr": 0.5,
              "Hf": 0.25, "Ta": 0.25}``).
        """
        sqs_level: int | None = None
        a_fracs: dict[str, float] = {}

        for parent in poscar_path.parents:
            if "sqs_lev=" in parent.name:
                match = re.search(r"sqs_lev=(\d+)", parent.name)
                if match:
                    sqs_level = int(match.group(1))
                for el, val in re.findall(r"a_([A-Za-z]+)=([0-9]*\.?[0-9]+)", parent.name):
                    a_fracs[el] = float(val)
                break

        return sqs_level, a_fracs

    def poscar_lattice_and_counts(
        self, poscar_path: Path
    ) -> tuple[np.ndarray, int, dict[str, int]]:
        """Parse a POSCAR file and return its lattice matrix and atom counts.

        Supports both the old VASP 4 format (counts on line 6, no element
        line) and the VASP 5 format (element symbols on line 6, counts on
        line 7).

        Args:
            poscar_path (Path): Path to the ``POSCAR`` (or ``CONTCAR``) file.

        Returns:
            tuple[np.ndarray, int, dict[str, int]]: A three-element tuple:

            - ``lattice``: 3×3 lattice matrix in Ångströms, shape ``(3, 3)``.
            - ``natoms``: Total number of atoms.
            - ``counts_map``: Dict mapping element symbol to atom count.
              Empty if the file uses VASP 4 format (no element symbols).
        """
        with open(poscar_path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        scale = float(lines[1])
        lattice = (
            np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)], dtype=float)
            * scale
        )

        i = 5
        toks = lines[i].split()
        if self._all_int(toks):
            elems: list[str] = []
            counts = [int(x) for x in toks]
            i += 1
        else:
            elems = toks
            i += 1
            counts = [int(x) for x in lines[i].split()]
            i += 1

        if i < len(lines) and lines[i].lower().startswith("s"):
            i += 1

        natoms = int(sum(counts))
        counts_map = {e: int(c) for e, c in zip(elems, counts)} if elems else {}
        return lattice, natoms, counts_map

    def cellpar_from_lattice(
        self, lattice: np.ndarray
    ) -> tuple[float, float, float, float, float, float]:
        """Compute lattice parameters from a 3×3 lattice matrix.

        Args:
            lattice (np.ndarray): 3×3 lattice matrix (rows = lattice vectors).

        Returns:
            tuple[float, float, float, float, float, float]:
            ``(a, b, c, alpha, beta, gamma)`` where lengths are in Ångströms
            and angles are in degrees.
        """
        a_vec, b_vec, c_vec = lattice

        def _angle(u: np.ndarray, v: np.ndarray) -> float:
            cos_val = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            cos_val = max(-1.0, min(1.0, cos_val))
            return float(math.degrees(math.acos(cos_val)))

        return (
            float(np.linalg.norm(a_vec)),
            float(np.linalg.norm(b_vec)),
            float(np.linalg.norm(c_vec)),
            _angle(b_vec, c_vec),
            _angle(a_vec, c_vec),
            _angle(a_vec, b_vec),
        )

    def scan_poscars(self, comp_dir: Path) -> pd.DataFrame:
        """Scan a composition directory and collect POSCAR geometry data.

        Recursively searches ``comp_dir`` for ``POSCAR`` files, parses each
        one, and returns a DataFrame with one row per POSCAR.

        Args:
            comp_dir (Path): Root directory for a single composition
                (e.g., ``Path("CrHfTa")``).

        Returns:
            pd.DataFrame: DataFrame with columns:

            - ``composition_folder`` (str)
            - ``phase_folder`` (str)
            - ``sqs_level`` (int | None)
            - ``sqs_a_fracs_json`` (str): JSON-encoded fractional compositions.
            - ``poscar_path`` (str)
            - ``volume_A3`` (float)
            - ``natoms`` (int)
            - ``volume_per_atom_A3`` (float | None)
            - ``a_A``, ``b_A``, ``c_A`` (float): Lattice lengths in Ångströms.
            - ``alpha_deg``, ``beta_deg``, ``gamma_deg`` (float): Angles in degrees.
            - ``poscar_counts_json`` (str): JSON-encoded element counts.
        """
        rows: list[dict] = []
        comp_name = comp_dir.name
        print(f"Checking for POSCARs in: {comp_dir}")

        for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
            phase_name = phase_dir.name
            for poscar_path in sorted(phase_dir.rglob("POSCAR")):
                print(f"  Found POSCAR: {poscar_path}")
                try:
                    lattice, natoms, counts_map = self.poscar_lattice_and_counts(poscar_path)
                except Exception as e:
                    print(f"  Read failed: {poscar_path} -> {e}")
                    continue

                sqs_level, a_fracs = self.parse_sqs_meta(poscar_path)
                vol = float(abs(np.linalg.det(lattice)))
                vpa = vol / natoms if natoms else None
                a, b, c, alpha, beta, gamma = self.cellpar_from_lattice(lattice)

                rows.append({
                    "composition_folder": comp_name,
                    "phase_folder": phase_name,
                    "sqs_level": sqs_level,
                    "sqs_a_fracs_json": json.dumps(a_fracs, sort_keys=True),
                    "poscar_path": str(poscar_path),
                    "volume_A3": vol,
                    "natoms": natoms,
                    "volume_per_atom_A3": vpa,
                    "a_A": a,
                    "b_A": b,
                    "c_A": c,
                    "alpha_deg": alpha,
                    "beta_deg": beta,
                    "gamma_deg": gamma,
                    "poscar_counts_json": json.dumps(counts_map, sort_keys=True),
                })

        self.data = pd.DataFrame(rows)
        return self.data

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _all_int(tokens: list[str]) -> bool:
        """Return True if every token in the list is a valid integer string."""
        try:
            [int(t) for t in tokens]
            return True
        except ValueError:
            return False
