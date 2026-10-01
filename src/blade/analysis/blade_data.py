"""Structural and energetic data extraction from POSCAR and energy files.

Provides :class:`BLADEData`, which scans a composition directory for
``POSCAR`` and ``energy`` files and returns a :class:`pandas.DataFrame`.
``BLADEVolume`` is a deprecated alias for backward compatibility.
"""

from __future__ import annotations

import json
import math
import re
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from pathlib import Path

    pass

__author__ = "Chase Katz"

# ---------------------------------------------------------------------------
# POSCAR format constants (line indices are 0-based after stripping blanks)
# ---------------------------------------------------------------------------
_poscar_scale_line: int = 1
"""Index of the universal scale-factor line in a POSCAR file."""

_poscar_lattice_start: int = 2
"""Index of the first lattice-vector line in a POSCAR file."""

_poscar_lattice_end: int = 5
"""One past the index of the last lattice-vector line (exclusive slice end)."""

_poscar_species_line: int = 5
"""Index of the species/count line; may be element symbols or integer counts."""


class BLADEData:
    """Extract structural and energetic data from BLADE POSCAR/energy trees.

    Recursively scans a composition directory for ``POSCAR`` files,
    parses lattice parameters and atom counts, reads adjacent ``energy``
    files, and stores results in a :class:`pandas.DataFrame` accessible
    via :attr:`data`.
    """

    def __init__(self) -> None:
        """Initialize BLADEData with an empty data store.

        Attributes:
            data: Populated by :meth:`scan_poscars`; ``None`` until then.
        """
        self.data: pd.DataFrame | None = None

    def parse_sqs_meta(self, poscar_path: Path) -> tuple[int | None, dict[str, float]]:
        """Extract SQS level and fractional composition from a POSCAR path.

        Walks up parent directories looking for a folder matching
        ``sqs_lev=<N>`` optionally followed by ``_a_<element>=<fraction>``
        tokens.

        Args:
            poscar_path: Path to the ``POSCAR`` file whose parent directories
                are searched.

        Returns:
            A 2-tuple ``(sqs_level, a_fracs)`` where ``sqs_level`` is the
            integer extracted from the ``sqs_lev=<N>`` folder name (or
            ``None`` if not found), and ``a_fracs`` maps element symbols to
            their sublattice fractions parsed from ``a_<el>=<val>`` tokens.
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

    def poscar_lattice_and_counts(self, poscar_path: Path) -> tuple[np.ndarray, int, dict[str, int]]:
        """Parse a POSCAR and return its lattice matrix and atom counts.

        Handles both VASP 4 (counts on line 6) and VASP 5 (elements on
        line 6, counts on line 7) formats.  ``counts_map`` is empty for
        VASP 4 files because element labels are absent.

        Args:
            poscar_path: Path to the ``POSCAR`` file to parse.

        Returns:
            A 3-tuple ``(lattice, natoms, counts_map)`` where ``lattice`` is
            a ``(3, 3)`` float array of row vectors scaled by the universal
            scale factor, ``natoms`` is the total atom count, and
            ``counts_map`` maps element symbol to atom count (empty for
            VASP 4 format).
        """
        with poscar_path.open() as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        scale = float(lines[_poscar_scale_line])
        lattice = (
            np.array(
                [[float(x) for x in lines[i].split()] for i in range(_poscar_lattice_start, _poscar_lattice_end)],
                dtype=float,
            )
            * scale
        )

        i = _poscar_species_line
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
        counts_map = {e: int(c) for e, c in zip(elems, counts, strict=False)} if elems else {}
        return lattice, natoms, counts_map

    def read_energy(self, poscar_path: Path) -> float | None:
        """Read total energy in eV from an ``energy`` file next to a POSCAR.

        Args:
            poscar_path: Path to the ``POSCAR`` file; the ``energy`` file is
                expected in the same directory.

        Returns:
            The energy as a float if the file exists and is parseable;
            ``None`` otherwise.
        """
        energy_path = poscar_path.parent / "energy"
        if not energy_path.exists():
            return None
        try:
            text = energy_path.read_text().strip()
            return float(text.split()[0])
        except Exception:
            return None

    def cellpar_from_lattice(self, lattice: np.ndarray) -> tuple[float, float, float, float, float, float]:
        """Return ``(a, b, c, alpha, beta, gamma)`` from a 3×3 lattice matrix.

        Args:
            lattice: A ``(3, 3)`` array whose rows are the lattice vectors
                **a**, **b**, **c** in Ångströms.

        Returns:
            A 6-tuple ``(a, b, c, alpha, beta, gamma)`` where the lengths
            are in Ångströms and the angles are in degrees.
        """
        a_vec, b_vec, c_vec = lattice

        def _angle(u: np.ndarray, v: np.ndarray) -> float:
            """Compute the angle between vectors *u* and *v* in degrees.

            Args:
                u: First vector.
                v: Second vector.

            Returns:
                Angle in degrees, clamped to ``[0, 180]`` to guard against
                floating-point values outside ``[-1, 1]`` for ``acos``.
            """
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
        """Scan a composition directory and collect structural and energy data.

        Recursively searches ``comp_dir`` for ``POSCAR`` files and reads the
        corresponding ``energy`` file when present. Returns a DataFrame with
        one row per POSCAR and columns:

        - ``composition_folder``, ``phase_folder`` (str)
        - ``sqs_level`` (int | None)
        - ``sqs_a_fracs_json`` (str): JSON-encoded sublattice fractions.
        - ``poscar_path`` (str)
        - ``volume_A3`` (float), ``natoms`` (int),
          ``volume_per_atom_A3`` (float | None)
        - ``a_A``, ``b_A``, ``c_A`` (float): lengths in Å.
        - ``alpha_deg``, ``beta_deg``, ``gamma_deg`` (float)
        - ``poscar_counts_json`` (str): JSON-encoded element counts.
        - ``energy_eV`` (float | None),
          ``energy_per_atom_eV`` (float | None)

        Args:
            comp_dir: Root directory for one composition, containing one
                sub-directory per phase.

        Returns:
            DataFrame with one row per successfully parsed POSCAR.  Also
            stored as :attr:`data`.
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

                energy = self.read_energy(poscar_path)
                epa = energy / natoms if (energy is not None and natoms) else None

                rows.append(
                    {
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
                        "energy_eV": energy,
                        "energy_per_atom_eV": epa,
                    }
                )

        self.data = pd.DataFrame(rows)
        return self.data

    @staticmethod
    def _all_int(tokens: list[str]) -> bool:
        """Return ``True`` if every token in *tokens* can be cast to ``int``.

        Used to distinguish VASP 4 POSCAR format (counts-only species line)
        from VASP 5 (element-symbol line followed by counts).

        Args:
            tokens: List of whitespace-split string tokens to test.

        Returns:
            ``True`` when all tokens are valid integers; ``False`` otherwise.
        """
        try:
            [int(t) for t in tokens]
            return True
        except ValueError:
            return False
