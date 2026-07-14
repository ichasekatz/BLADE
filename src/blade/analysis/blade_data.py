"""Structural and energetic data extraction from POSCAR and energy files.

Provides :class:`BLADEData`, which scans a composition directory for
``POSCAR`` and ``energy`` files and returns a :class:`pandas.DataFrame`.
``BLADEVolume`` is a deprecated alias for backward compatibility.
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


class BLADEData:
    """Extract structural and energetic data from BLADE POSCAR/energy trees."""

    def __init__(self) -> None:
        self.data: pd.DataFrame | None = None

    def parse_sqs_meta(self, poscar_path: Path) -> tuple[int | None, dict[str, float]]:
        """Extract SQS level and fractional composition from a POSCAR path.

        Walks up parent directories looking for a folder matching
        ``sqs_lev=<N>`` optionally followed by ``_a_<element>=<fraction>`` tokens.
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

        Handles both VASP 4 (counts on line 6) and VASP 5 (elements on line 6,
        counts on line 7) formats. ``counts_map`` is empty for VASP 4 files.
        """
        with open(poscar_path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        scale = float(lines[1])
        lattice = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)], dtype=float) * scale

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

    def read_energy(self, poscar_path: Path) -> float | None:
        """Read total energy in eV from an ``energy`` file next to a POSCAR."""
        energy_path = poscar_path.parent / "energy"
        if not energy_path.exists():
            return None
        try:
            text = energy_path.read_text().strip()
            return float(text.split()[0])
        except Exception:
            return None

    def cellpar_from_lattice(self, lattice: np.ndarray) -> tuple[float, float, float, float, float, float]:
        """Return ``(a, b, c, alpha, beta, gamma)`` from a 3×3 lattice matrix."""
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
        """Scan a composition directory and collect structural and energy data.

        Recursively searches ``comp_dir`` for ``POSCAR`` files and reads the
        corresponding ``energy`` file when present. Returns a DataFrame with
        one row per POSCAR and columns:

        - ``composition_folder``, ``phase_folder`` (str)
        - ``sqs_level`` (int | None)
        - ``sqs_a_fracs_json`` (str): JSON-encoded sublattice fractions.
        - ``poscar_path`` (str)
        - ``volume_A3`` (float), ``natoms`` (int), ``volume_per_atom_A3`` (float | None)
        - ``a_A``, ``b_A``, ``c_A`` (float): lengths in Å.
        - ``alpha_deg``, ``beta_deg``, ``gamma_deg`` (float)
        - ``poscar_counts_json`` (str): JSON-encoded element counts.
        - ``energy_eV`` (float | None), ``energy_per_atom_eV`` (float | None)
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
        try:
            [int(t) for t in tokens]
            return True
        except ValueError:
            return False
