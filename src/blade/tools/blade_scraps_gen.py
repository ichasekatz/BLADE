"""SCRAPS-based SQS generation — drop-in replacement for BladeSQS.

``ScrapsSQSGen`` keeps ATAT's ``sqs2tdb -mk`` to scaffold the
``sqsdb_lev=*`` directory tree (so BladeTDBGen / MaterialsFramework work
unchanged) and replaces ``mcsqs`` with SCRAPS in each composition directory.
``SCRAPS.vasp`` is converted to ATAT ``bestsqs.out`` format after each run.

Placeholder element trick:
    Each (sublattice_letter, rank) slot is assigned a distinct dummy element
    from ``_PLACEHOLDER_POOL`` so the output can be unambiguously re-labeled
    as ATAT labels (``a_A``, ``a_B``, ``b_A``, …) without position-matching.

MPI warning:
    macOS + Anaconda: ``which mpirun`` often returns Hydra (MPICH), not the
    Homebrew OpenMPI that built ``scraps``.  Always resolve via
    ``_find_mpirun()`` from ``scraps_ase.py`` — never use
    ``shutil.which("mpirun")`` alone.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
from ase import Atoms
from ase.cell import Cell

_PLACEHOLDER_POOL = [
    "Ni", "Co", "Cu", "Ag", "Au", "Pt", "Pd",
    "Rh", "Ir", "Os", "Ru", "Re",
]


# ------------------------------------------------------------------
# Primitive-cell helpers
# ------------------------------------------------------------------

def _prim_cell_matrix(phases_dict: dict) -> np.ndarray:
    return Cell.fromcellpar([
        phases_dict["a"],     phases_dict["b"],    phases_dict["c"],
        phases_dict["alpha"], phases_dict["beta"],  phases_dict["gamma"],
    ]).array


def _parse_coords(coords_str: str) -> tuple[list[list[float]], list[str]]:
    positions, labels = [], []
    for line in coords_str.strip().splitlines():
        parts = line.strip().split()
        if len(parts) >= 4:
            positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
            labels.append(parts[3])
    return positions, labels


def _build_prim_atoms(phases_dict: dict, placeholder: str = "X") -> tuple[Atoms, list[str]]:
    cell = _prim_cell_matrix(phases_dict)
    pos, labels = _parse_coords(phases_dict["coords"])
    # Use a distinct placeholder element per sublattice letter so the POSCAR
    # has separate species blocks for each sublattice.  If all variable sites
    # used the same element, ASE would group them together and SCRAPS BASIS
    # indices would point to the wrong atoms.
    letter_to_placeholder: dict[str, str] = {}
    pool_offset = 0
    for lbl in labels:
        if lbl.islower() and len(lbl) == 1 and lbl not in letter_to_placeholder:
            letter_to_placeholder[lbl] = _PLACEHOLDER_POOL[pool_offset % len(_PLACEHOLDER_POOL)]
            pool_offset += 1
    symbols = [
        letter_to_placeholder[lbl] if (lbl.islower() and len(lbl) == 1) else lbl
        for lbl in labels
    ]
    return Atoms(symbols=symbols, scaled_positions=pos, cell=cell, pbc=True), labels


def _parse_sqsdb_comp(dirname: str) -> dict[str, list[float]]:
    return {
        m.group(1): [float(x) for x in m.group(2).split(",")]
        for m in re.finditer(r"_([a-z])=([\d.,]+)", dirname)
    }


def _var_sublattices(site_labels: list[str]) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for i, lbl in enumerate(site_labels):
        if lbl.islower() and len(lbl) == 1:
            result.setdefault(lbl, []).append(i)
    return result


def _fixed_sites(site_labels: list[str]) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for i, lbl in enumerate(site_labels):
        if not (lbl.islower() and len(lbl) == 1):
            result.setdefault(lbl, []).append(i)
    return result


# ------------------------------------------------------------------
# SCRAPS INPUT config builder
# ------------------------------------------------------------------

def _build_scraps_config(
    site_labels: list[str],
    sqsdb_comp: dict[str, list[float]],
    cell_dim: tuple[int, int, int],
    fix_multibasis_sublattice: bool = True,
) -> dict:
    """Build SCRAPS INPUT parameters for one sqsdb_lev composition.

    Returns a dict with keys:
        species, element_counts, sublattice_blocks, spectator_pairs,
        placeholder_to_atat
    """
    var_subs = _var_sublattices(site_labels)
    fix_sites = _fixed_sites(site_labels)
    cell_vol = cell_dim[0] * cell_dim[1] * cell_dim[2]

    species: list[str] = []
    element_counts: list[int] = []
    sublattice_blocks: list[tuple] = []
    spectator_pairs: list[tuple] = []
    placeholder_to_atat: dict[str, str] = {}

    pool_idx = 0
    rank_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for letter in sorted(var_subs.keys()):
        basis_idxs = var_subs[letter]
        n_sites = len(basis_idxs) * cell_vol
        fracs = sqsdb_comp.get(letter, [1.0])
        pure = len(fracs) == 1

        if pure:
            ph = _PLACEHOLDER_POOL[pool_idx % len(_PLACEHOLDER_POOL)]
            pool_idx += 1
            sp_idx = len(species)
            species.append(ph)
            element_counts.append(int(round(n_sites)))
            placeholder_to_atat[ph] = f"{letter}_A"
            for bi in basis_idxs:
                spectator_pairs.append((bi, sp_idx))
        else:
            if len(basis_idxs) == 1 or not fix_multibasis_sublattice:
                # Single-basis sublattice OR fix disabled → standard SUBLATTICE block.
                sub_sp_idxs: list[int] = []
                for rank, frac in enumerate(fracs):
                    ph = _PLACEHOLDER_POOL[pool_idx % len(_PLACEHOLDER_POOL)]
                    pool_idx += 1
                    sp_idx = len(species)
                    species.append(ph)
                    n = int(round(frac * n_sites))
                    element_counts.append(n)
                    placeholder_to_atat[ph] = f"{letter}_{rank_labels[rank]}"
                    sub_sp_idxs.append(sp_idx)
                sub_id = len(sublattice_blocks)
                sublattice_blocks.append((sub_id, list(basis_idxs), sub_sp_idxs))
            else:
                # Multi-basis variable sublattice with fix enabled: one spectator per
                # basis atom cycling through fraction-ranked species.  Avoids SCRAPS's
                # uniformity check failure when basis atoms are 0.5c apart (e.g. 2d
                # Wyckoff in MAX phases).  Elements are present but not WC-optimised.
                for rank, bi in enumerate(basis_idxs):
                    species_rank = rank % len(fracs)
                    ph = _PLACEHOLDER_POOL[pool_idx % len(_PLACEHOLDER_POOL)]
                    pool_idx += 1
                    sp_idx = len(species)
                    species.append(ph)
                    element_counts.append(cell_vol)
                    placeholder_to_atat[ph] = f"{letter}_{rank_labels[species_rank]}"
                    spectator_pairs.append((bi, sp_idx))

    for elem, b_idxs in sorted(fix_sites.items()):
        n_sites = len(b_idxs) * cell_vol
        sp_idx = len(species)
        species.append(elem)
        element_counts.append(int(round(n_sites)))
        placeholder_to_atat[elem] = elem
        for bi in b_idxs:
            spectator_pairs.append((bi, sp_idx))

    return {
        "species": species,
        "element_counts": element_counts,
        "sublattice_blocks": sublattice_blocks,
        "spectator_pairs": spectator_pairs,
        "placeholder_to_atat": placeholder_to_atat,
    }


# ------------------------------------------------------------------
# SCRAPS.vasp → ATAT bestsqs.out
# ------------------------------------------------------------------

def _write_bestsqs(
    scraps_vasp: Path,
    prim_cell: np.ndarray,
    species: list[str],
    placeholder_to_atat: dict[str, str],
    output_path: Path,
    read_scraps_output,
) -> None:
    """Convert SCRAPS.vasp to ATAT bestsqs.out format."""
    atoms = read_scraps_output(scraps_vasp, species=species)
    sup_cell = atoms.cell.array
    cart_pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()

    inv_prim = np.linalg.inv(prim_cell)
    sup_prim = sup_cell @ inv_prim
    pos_prim = cart_pos @ inv_prim

    lines: list[str] = []
    for v in prim_cell:
        lines.append(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")
    for v in sup_prim:
        lines.append(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")
    for pos, sym in zip(pos_prim, syms):
        label = placeholder_to_atat.get(sym, sym)
        lines.append(f"{pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f} {label}")

    output_path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------
# ScrapsSQSGen
# ------------------------------------------------------------------

class ScrapsSQSGen:
    """SQS generation using SCRAPS instead of ATAT corrdump + mcsqs.

    Keeps ATAT's ``sqs2tdb -mk`` to scaffold the ``sqsdb_lev=*`` directory
    structure so BladeTDBGen works unchanged, then replaces mcsqs with
    SCRAPS in each composition directory.

    Args:
        phases_dict: Phase prototype dict (same format as BladeSQS).
        sqsgen_levels: List of level dicts (same as BladeSQS).
        level: Composition mesh level.
        len_comp: Number of components.
        skip_existing_sqs: Skip if ``bestsqs.out`` already exists.
        sublattice_map: Per-sublattice element assignment (same as BladeSQS).
        sqsgen_in: Verbatim ``sqsgen.in`` override (bypasses auto-generation).
        scraps_bin: Path to compiled ``scraps`` binary.
        scraps_tools: Path to ``scraps-perpair/tools/`` (for ``scraps_ase`` imports).
        mpirun: Path to mpirun. Auto-resolved via ``_find_mpirun()`` if None.
        ranks: Number of MPI ranks.
        auto_budget: SCRAPS search budget (1=basic, 2=thorough, 3=exhaustive).
        max_shellnum: Maximum shell number for Warren-Cowley targets.
    """

    def __init__(
        self,
        phases_dict: dict,
        sqsgen_levels: list[dict],
        level: int,
        len_comp: int,
        skip_existing_sqs: bool = False,
        sublattice_map: Optional[dict] = None,
        sqsgen_in: Optional[str] = None,
        fixed_compositions: Optional[dict] = None,
        scraps_bin: Optional[Path] = None,
        scraps_tools: Optional[Path] = None,
        mpirun: Optional[str] = None,
        ranks: int = 4,
        auto_budget: int = 2,
        max_shellnum: int = 2,
        fix_multibasis_sublattice: bool = True,
    ) -> None:
        self.phases_dict              = phases_dict
        self.sqsgen_levels            = sqsgen_levels
        self.level                    = level
        self.len_comp                 = len_comp
        self.skip_existing            = skip_existing_sqs
        self.sublattice_map           = sublattice_map
        self.sqsgen_in                = sqsgen_in
        self.fixed_compositions       = fixed_compositions or {}
        self.scraps_bin               = Path(scraps_bin) if scraps_bin else None
        self.scraps_tools             = Path(scraps_tools) if scraps_tools else None
        self.ranks                    = ranks
        self.auto_budget              = auto_budget
        self.max_shellnum             = max_shellnum
        self.fix_multibasis_sublattice = fix_multibasis_sublattice

        _, self._site_labels = _build_prim_atoms(phases_dict)
        self._prim_cell      = _prim_cell_matrix(phases_dict)

        # Lazy-import scraps_ase after sys.path is extended
        self._make_scraps_inputs = None
        self._read_scraps_output = None
        self._find_mpirun        = None
        self._mpirun             = mpirun  # resolved on first use

    def _load_scraps_ase(self) -> None:
        if self._make_scraps_inputs is not None:
            return
        if self.scraps_tools is not None and str(self.scraps_tools) not in sys.path:
            sys.path.insert(0, str(self.scraps_tools))
        from scraps_ase import make_scraps_inputs, read_scraps_output, _find_mpirun
        self._make_scraps_inputs = make_scraps_inputs
        self._read_scraps_output = read_scraps_output
        self._find_mpirun        = _find_mpirun
        if self._mpirun is None:
            self._mpirun = _find_mpirun()

    # ── Public interface (mirrors BladeSQS) ─────────────────────────

    def sqs_struct(self) -> tuple[str, str]:
        """Build ``sqsgen.in`` + ``rndstr.skel`` text via BladeSQS helpers."""
        from blade.tools.blade_sqsgen import BladeSQS
        _sqs = BladeSQS(
            phases_dict=self.phases_dict,
            sqsgen_levels=self.sqsgen_levels,
            level=self.level,
            len_comp=self.len_comp,
            skip_existing_sqs=self.skip_existing,
            sublattice_map=self.sublattice_map,
            sqsgen_in=self.sqsgen_in,
            fixed_compositions=self.fixed_compositions,
        )
        return _sqs.sqs_struct()

    def sqs_gen(
        self,
        phase: dict,
        paths: list,
        iter: float,
        params: dict,
    ) -> None:
        """Generate SQS structures via SCRAPS for all sqsdb_lev compositions.

        1. Writes ``rndstr.skel`` + ``sqsgen.in`` (same as BladeSQS).
        2. Runs ``sqs2tdb -mk`` to scaffold ``sqsdb_lev=*`` directories.
        3. For each mixed-composition ``sqsdb_lev`` dir, runs SCRAPS and
           converts ``SCRAPS.vasp`` → ``bestsqs.out``.

        Args:
            phase  : phase entry dict (must contain ``"lattice"``).
            paths  : ``[path0, path1, path2]`` path bundle.
            iter   : unused (SCRAPS uses ``auto_budget`` / early-stop).
            params : must contain ``"super_cell_size"`` ``tuple[int,int,int]``.
        """
        self._load_scraps_ase()

        if self.scraps_bin is None or not self.scraps_bin.exists():
            raise FileNotFoundError(
                f"SCRAPS binary not found at {self.scraps_bin}. "
                "Pass scraps_bin= or run build.sh in scraps-perpair/."
            )
        if self._mpirun is None:
            raise RuntimeError("No mpirun found. Set MPIRUN env var or install OpenMPI.")

        lattice  = phase["lattice"]
        cell_dim = params["super_cell_size"]
        dir_name = Path(paths[1]) / "Files" / "SQS" / f"{lattice}_{self.len_comp}"

        if self.skip_existing and dir_name.exists():
            print(f"Skipping SQS scaffold for {lattice}_{self.len_comp}: exists")
        else:
            if dir_name.exists():
                shutil.rmtree(dir_name)
            dir_name.mkdir(parents=True, exist_ok=True)
            sqsgen, rndstr = self.sqs_struct()
            (dir_name / "rndstr.skel").write_text(rndstr)
            (dir_name / "sqsgen.in").write_text(sqsgen)

            result = subprocess.run(
                ["sqs2tdb", "-mk"],
                cwd=dir_name, capture_output=True, text=True, check=False,
            )
            print(result.stdout)
            if result.stderr:
                print("sqs2tdb stderr:", result.stderr)

        prim_atoms, _ = _build_prim_atoms(self.phases_dict)

        for sqsdir in sorted(dir_name.glob("sqsdb_lev=*/")):
            self._run_scraps_in_dir(sqsdir, prim_atoms, cell_dim)

        self._write_objective_summary(dir_name)

    # ── Private helpers ─────────────────────────────────────────────

    def _write_endmember_bestsqs(
        self,
        sqsdir: Path,
        cell_dim: tuple[int, int, int],
    ) -> None:
        """Write bestsqs.out in prim-vector-unit format for a pure endmember.

        Pure-composition dirs don't need SCRAPS (atom order is trivial), but
        ``_poscar_from_bestsqs`` requires a ``bestsqs.out`` in our custom
        prim-vector-unit format.  Any pre-existing file from a mcsqs run uses
        Cartesian supercell vectors and would produce a wrong POSCAR, so we
        always regenerate here.

        Supercell is tiled as ``cell_dim[0] × cell_dim[1] × cell_dim[2]``
        unit cells.  All variable-sublattice sites receive the ``_A`` label
        (single occupant); fixed sites keep their element symbol.

        Args:
            sqsdir   : ``sqsdb_lev=*`` directory.
            cell_dim : Supercell dimensions ``(n0, n1, n2)``.
        """
        prim_cell = self._prim_cell
        pos_list, site_labels = _parse_coords(self.phases_dict["coords"])
        n0, n1, n2 = int(cell_dim[0]), int(cell_dim[1]), int(cell_dim[2])

        # Diagonal supercell matrix in prim-vector units
        sup_prim = np.diag([float(n0), float(n1), float(n2)])

        lines: list[str] = []
        for v in prim_cell:
            lines.append(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")
        for v in sup_prim:
            lines.append(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")

        for i in range(n0):
            for j in range(n1):
                for k in range(n2):
                    for (pos, lbl) in zip(pos_list, site_labels):
                        px = pos[0] + i
                        py = pos[1] + j
                        pz = pos[2] + k
                        # Variable sublattice → single occupant a_A, b_A, etc.
                        if lbl.islower() and len(lbl) == 1:
                            atat_label = f"{lbl}_A"
                        else:
                            atat_label = lbl  # fixed element (e.g., B)
                        lines.append(f"{px:.10f} {py:.10f} {pz:.10f} {atat_label}")

        (sqsdir / "bestsqs.out").write_text("\n".join(lines) + "\n")
        print(f"  Endmember bestsqs.out written: {sqsdir.name}")

    def _run_scraps_in_dir(
        self,
        sqsdir: Path,
        prim_atoms: Atoms,
        cell_dim: tuple[int, int, int],
    ) -> None:
        sqsdb_comp = _parse_sqsdb_comp(sqsdir.name)

        all_fracs = [f for fracs in sqsdb_comp.values() for f in fracs]
        non_zero  = [f for f in all_fracs if f > 0.0]
        if non_zero and all(f == 1.0 for f in non_zero):
            # Always regenerate endmember bestsqs.out in our prim-vector-unit
            # format regardless of any pre-existing mcsqs file (wrong format).
            self._write_endmember_bestsqs(sqsdir, cell_dim)
            return

        if self.skip_existing and (sqsdir / "bestsqs.out").exists():
            print(f"Skipping SCRAPS for {sqsdir.name}: bestsqs.out exists")
            return

        if not (sqsdir / "rndstr.in").exists():
            print(f"No rndstr.in in {sqsdir}, skipping")
            return

        cfg = _build_scraps_config(
            site_labels=self._site_labels,
            sqsdb_comp=sqsdb_comp,
            cell_dim=cell_dim,
            fix_multibasis_sublattice=self.fix_multibasis_sublattice,
        )

        n_basis  = len(self._site_labels)
        cell_vol = cell_dim[0] * cell_dim[1] * cell_dim[2]
        expected = n_basis * cell_vol
        total    = sum(cfg["element_counts"])
        if total != expected:
            print(
                f"WARNING: element_counts sum {total} != expected {expected} "
                f"in {sqsdir.name} — rounding mismatch, skipping"
            )
            return

        work_dir = sqsdir / "_scraps_run"
        work_dir.mkdir(exist_ok=True)

        self._make_scraps_inputs(
            work_dir=work_dir,
            atoms=prim_atoms,
            element_counts=cfg["element_counts"],
            cell_dim=cell_dim,
            max_shellnum=self.max_shellnum,
            auto_budget=self.auto_budget,
            sublattices=cfg["sublattice_blocks"] or None,
            spectator_basis_species=cfg["spectator_pairs"] or None,
        )

        scraps_dest = work_dir / "scraps"
        shutil.copy(self.scraps_bin, scraps_dest)
        os.chmod(scraps_dest, 0o755)

        log_path = work_dir / "scraps.log"
        cmd = [self._mpirun, "--oversubscribe", "-np", str(self.ranks), "./scraps"]
        print(f"SCRAPS {sqsdir.name}: {' '.join(cmd)}")
        t0 = time.time()
        with log_path.open("w") as fh:
            rc = subprocess.run(
                cmd, cwd=work_dir, stdout=fh, stderr=subprocess.STDOUT
            ).returncode
        elapsed = time.time() - t0

        scraps_vasp = work_dir / "SCRAPS.vasp"
        if rc != 0 or not scraps_vasp.exists():
            print(f"  SCRAPS failed (rc={rc}) for {sqsdir.name} ({elapsed:.1f}s)")
            return

        fitness = None
        for line in log_path.read_text().splitlines():
            if line.startswith("FINAL fitness"):
                try:
                    fitness = float(line.split("=")[-1].strip())
                except Exception:
                    pass
                break
        print(f"  fitness={fitness}  wall={elapsed:.1f}s")

        _write_bestsqs(
            scraps_vasp=scraps_vasp,
            prim_cell=self._prim_cell,
            species=cfg["species"],
            placeholder_to_atat=cfg["placeholder_to_atat"],
            output_path=sqsdir / "bestsqs.out",
            read_scraps_output=self._read_scraps_output,
        )
        shutil.copy(scraps_vasp, sqsdir / "SCRAPS.vasp")
        shutil.copy(log_path, sqsdir / "scraps.log")
        print(f"  bestsqs.out written → {sqsdir / 'bestsqs.out'}")

    def _write_objective_summary(self, parent_dir: Path) -> None:
        output_file = parent_dir / "objective_functions.txt"
        lines = ["folder\tscraps_fitness"]
        for sqsdir in sorted(parent_dir.glob("sqsdb_lev=*/")):
            log = sqsdir / "scraps.log"
            fitness = None
            if log.exists():
                for line in log.read_text().splitlines():
                    if line.startswith("FINAL fitness"):
                        try:
                            fitness = float(line.split("=")[-1].strip())
                        except Exception:
                            pass
                        break
            if fitness is not None:
                lines.append(f"{sqsdir.name}\t{fitness}")
        output_file.write_text("\n".join(lines))
