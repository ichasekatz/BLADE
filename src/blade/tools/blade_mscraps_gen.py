"""m_scraps-based SQS generation — drop-in replacement for BladeSQS.

``MScrapsSQSGen`` mirrors ``ScrapsSQSGen`` exactly, replacing the ``scraps_ase``
script import with the ``m_scraps`` Python package.  The additional m_scraps
keywords — SRO targets, higher-order figure correlations, and MC kernel
selection — are exposed as constructor parameters.

Uses the same binary (``SCRAPs/scraps``) and same ``bestsqs.out`` format as
``ScrapsSQSGen``, so it is interchangeable with the SCRAPS integration.
"""

from __future__ import annotations

import contextlib
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.cell import Cell

_placeholder_pool = [
    "Ni", "Co", "Cu", "Ag", "Au", "Pt", "Pd", "Rh", "Ir", "Os", "Ru", "Re",
]


# ------------------------------------------------------------------
# Primitive-cell helpers (shared with blade_scraps_gen)
# ------------------------------------------------------------------


def _prim_cell_matrix(phases_dict: dict) -> np.ndarray:
    """Return 3×3 cell matrix from phases_dict lattice parameters."""
    return Cell.fromcellpar(
        [
            phases_dict["a"],
            phases_dict["b"],
            phases_dict["c"],
            phases_dict["alpha"],
            phases_dict["beta"],
            phases_dict["gamma"],
        ]
    ).array


def _parse_coords(coords_str: str) -> tuple[list[list[float]], list[str]]:
    """Parse fractional coordinate block into positions and site labels."""
    positions, labels = [], []
    for line in coords_str.strip().splitlines():
        parts = line.strip().split()
        if len(parts) >= 4:
            positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
            labels.append(parts[3])
    return positions, labels


def _build_prim_atoms(phases_dict: dict) -> tuple[Atoms, list[str]]:
    """Build ASE Atoms for the primitive cell with placeholder species.

    Args:
        phases_dict: Phase prototype dict with lattice params and ``coords``.

    Returns:
        (atoms, site_labels) — atoms has one placeholder element per variable
        sublattice letter; site_labels is the raw label list from coords.
    """
    cell = _prim_cell_matrix(phases_dict)
    pos, labels = _parse_coords(phases_dict["coords"])
    letter_to_placeholder: dict[str, str] = {}
    pool_offset = 0
    for lbl in labels:
        if lbl.islower() and len(lbl) == 1 and lbl not in letter_to_placeholder:
            letter_to_placeholder[lbl] = _placeholder_pool[pool_offset % len(_placeholder_pool)]
            pool_offset += 1
    symbols = [
        letter_to_placeholder[lbl] if (lbl.islower() and len(lbl) == 1) else lbl
        for lbl in labels
    ]
    return Atoms(symbols=symbols, scaled_positions=pos, cell=cell, pbc=True), labels


def _parse_sqsdb_comp(dirname: str) -> dict[str, list[float]]:
    """Parse composition fractions from a ``sqsdb_lev=*`` directory name."""
    return {
        m.group(1): [float(x) for x in m.group(2).split(",")]
        for m in re.finditer(r"_([a-z])=([\d.,]+)", dirname)
    }


def _var_sublattices(site_labels: list[str]) -> dict[str, list[int]]:
    """Map sublattice letter → site indices for variable sublattices."""
    result: dict[str, list[int]] = {}
    for i, lbl in enumerate(site_labels):
        if lbl.islower() and len(lbl) == 1:
            result.setdefault(lbl, []).append(i)
    return result


def _fixed_sites(site_labels: list[str]) -> dict[str, list[int]]:
    """Map element symbol → site indices for fixed (non-variable) sites."""
    result: dict[str, list[int]] = {}
    for i, lbl in enumerate(site_labels):
        if not (lbl.islower() and len(lbl) == 1):
            result.setdefault(lbl, []).append(i)
    return result


def _resolve_mpirun(hint: str | None = None) -> str | None:
    """Resolve the mpirun executable, preferring an explicit hint or MPIRUN env var.

    Args:
        hint: Explicit path provided by the caller; returned as-is if given.

    Returns:
        Resolved mpirun path string, or None if not found.
    """
    if hint:
        return hint
    env_val = os.environ.get("MPIRUN")
    if env_val:
        return env_val
    return shutil.which("mpirun")


# ------------------------------------------------------------------
# SCRAPS INPUT config builder (identical logic to blade_scraps_gen)
# ------------------------------------------------------------------


def _build_scraps_config(
    site_labels: list[str],
    sqsdb_comp: dict[str, list[float]],
    cell_dim: tuple[int, int, int],
    fix_multibasis_sublattice: bool = True,
) -> dict:
    """Build SCRAPS INPUT parameters for one sqsdb_lev composition.

    Args:
        site_labels: Site label list from the primitive cell coords block.
        sqsdb_comp: Parsed composition dict from the sqsdb directory name.
        cell_dim: Supercell tile dimensions ``(n0, n1, n2)``.
        fix_multibasis_sublattice: Use per-atom spectators for multi-basis
            variable sublattices to avoid SCRAPS uniformity check failures.

    Returns:
        Dict with keys: species, element_counts, sublattice_blocks,
        spectator_pairs, placeholder_to_atat.
    """
    var_subs = _var_sublattices(site_labels)
    fix_sites_map = _fixed_sites(site_labels)
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
            ph = _placeholder_pool[pool_idx % len(_placeholder_pool)]
            pool_idx += 1
            sp_idx = len(species)
            species.append(ph)
            element_counts.append(int(round(n_sites)))
            placeholder_to_atat[ph] = f"{letter}_A"
            for bi in basis_idxs:
                spectator_pairs.append((bi, sp_idx))
        elif len(basis_idxs) == 1 or not fix_multibasis_sublattice:
            sub_sp_idxs: list[int] = []
            for rank, frac in enumerate(fracs):
                ph = _placeholder_pool[pool_idx % len(_placeholder_pool)]
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
            # Multi-basis variable sublattice: per-atom spectators to avoid
            # SCRAPS uniformity check failure on 2d Wyckoff positions.
            for rank, bi in enumerate(basis_idxs):
                species_rank = rank % len(fracs)
                ph = _placeholder_pool[pool_idx % len(_placeholder_pool)]
                pool_idx += 1
                sp_idx = len(species)
                species.append(ph)
                element_counts.append(cell_vol)
                placeholder_to_atat[ph] = f"{letter}_{rank_labels[species_rank]}"
                spectator_pairs.append((bi, sp_idx))

    for elem, b_idxs in sorted(fix_sites_map.items()):
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
    """Convert SCRAPS.vasp to ATAT bestsqs.out format.

    Args:
        scraps_vasp: Path to the ``SCRAPS.vasp`` output file.
        prim_cell: 3×3 primitive cell matrix in Å.
        species: Species list matching the SCRAPS INPUT ordering.
        placeholder_to_atat: Map from placeholder element → ATAT label.
        output_path: Destination ``bestsqs.out`` path.
        read_scraps_output: Callable from m_scraps or scraps_ase.
    """
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
    for pos, sym in zip(pos_prim, syms, strict=False):
        label = placeholder_to_atat.get(sym, sym)
        lines.append(f"{pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f} {label}")

    output_path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------
# MScrapsSQSGen
# ------------------------------------------------------------------


class MScrapsSQSGen:
    """SQS generation using the m_scraps Python package instead of ATAT mcsqs.

    Mirrors ``ScrapsSQSGen`` but imports from the ``m_scraps`` package (a
    proper Python package rather than a tools-directory script), and exposes
    m_scraps-specific parameters: SRO targets, higher-order figure correlations,
    and MC kernel selection.

    The SCRAPS binary (``SCRAPs/scraps``) is the same binary used by SCRAPS.
    ``m_scraps`` is a Python wrapper that adds the ``SUBLATTICE``, ``SRO``,
    and ``FIGURES`` INPUT keywords.

    Args:
        phases_dict: Phase prototype dict (same format as BladeSQS).
        sqsgen_levels: List of level dicts (same as BladeSQS).
        level: Composition mesh level.
        len_comp: Number of components.
        skip_existing_sqs: Skip if ``bestsqs.out`` already exists.
        sublattice_map: Per-sublattice element assignment (same as BladeSQS).
        sqsgen_in: Verbatim ``sqsgen.in`` override.
        fixed_compositions: Fixed-composition sublattice overrides.
        scraps_bin: Path to the compiled ``scraps`` binary.
        mpirun: mpirun executable. Resolved from ``MPIRUN`` env var or PATH if None.
        ranks: Number of MPI ranks.
        auto_budget: SCRAPS search budget (1=basic, 2=thorough, 3=exhaustive).
        max_shellnum: Maximum shell number for Warren-Cowley pair targets.
        fix_multibasis_sublattice: Use per-atom spectators for multi-basis variable
            sublattices; avoids SCRAPS uniformity check failure on 2d Wyckoff sites.
        sro_targets: Per-pair SRO targets as list of ``(i, j, shell, alpha)`` tuples
            where ``i``, ``j`` are species indices (0-based) and ``alpha`` is the
            Warren-Cowley parameter target.
        figure_cutoffs: Cutoff distances (Å) for higher-order figure correlations.
        figure_species: Species indices included in figure correlation targets.
        figure_objective: Figure objective function — ``"minimize"`` or ``"match"``.
        figure_wd: Distance weight for figure objective.
        figure_wn: Neighbour-count weight for figure objective.
        figure_wr: Rank weight for figure objective.
        figure_tol: Convergence tolerance for figure objective.
        mc_kernel: MC kernel — ``"sa"`` (simulated annealing) or ``"pt"``
            (parallel tempering). Passed to SCRAPS via ``MC_KERNEL`` INPUT keyword.
    """

    def __init__(
        self,
        phases_dict: dict,
        sqsgen_levels: list[dict],
        level: int,
        len_comp: int,
        skip_existing_sqs: bool = False,
        sublattice_map: dict | None = None,
        sqsgen_in: str | None = None,
        fixed_compositions: dict | None = None,
        scraps_bin: Path | None = None,
        mpirun: str | None = None,
        ranks: int = 4,
        auto_budget: int = 2,
        max_shellnum: int = 2,
        fix_multibasis_sublattice: bool = True,
        sro_targets: list[tuple[int, int, int, float]] | None = None,
        figure_cutoffs: list[float] | None = None,
        figure_species: list[int] | None = None,
        figure_objective: str = "minimize",
        figure_wd: float = 1.0,
        figure_wn: float = 0.75,
        figure_wr: float = 1.0,
        figure_tol: float = 1e-4,
        mc_kernel: str = "sa",
    ) -> None:
        self.phases_dict = phases_dict
        self.sqsgen_levels = sqsgen_levels
        self.level = level
        self.len_comp = len_comp
        self.skip_existing = skip_existing_sqs
        self.sublattice_map = sublattice_map
        self.sqsgen_in = sqsgen_in
        self.fixed_compositions = fixed_compositions or {}
        self.scraps_bin = Path(scraps_bin) if scraps_bin else None
        self.ranks = ranks
        self.auto_budget = auto_budget
        self.max_shellnum = max_shellnum
        self.fix_multibasis_sublattice = fix_multibasis_sublattice
        self.sro_targets = sro_targets
        self.figure_cutoffs = figure_cutoffs
        self.figure_species = figure_species
        self.figure_objective = figure_objective
        self.figure_wd = figure_wd
        self.figure_wn = figure_wn
        self.figure_wr = figure_wr
        self.figure_tol = figure_tol
        self.mc_kernel = mc_kernel

        _, self._site_labels = _build_prim_atoms(phases_dict)
        self._prim_cell = _prim_cell_matrix(phases_dict)

        self._make_scraps_inputs = None
        self._read_scraps_output = None
        self._mpirun: str | None = _resolve_mpirun(mpirun)

    def _load_m_scraps(self) -> None:
        """Lazy-import m_scraps package functions on first use."""
        if self._make_scraps_inputs is not None:
            return
        from m_scraps import make_scraps_inputs, read_scraps_output

        self._make_scraps_inputs = make_scraps_inputs
        self._read_scraps_output = read_scraps_output

    # ── Public interface (mirrors BladeSQS / ScrapsSQSGen) ──────────

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
        params: dict,
    ) -> None:
        """Generate SQS structures via m_scraps for all sqsdb_lev compositions.

        1. Writes ``rndstr.skel`` + ``sqsgen.in`` (same as BladeSQS).
        2. Runs ``sqs2tdb -mk`` to scaffold ``sqsdb_lev=*`` directories.
        3. For each mixed-composition ``sqsdb_lev`` dir, runs m_scraps and
           converts ``SCRAPS.vasp`` → ``bestsqs.out``.

        Args:
            phase  : phase entry dict (must contain ``"lattice"``).
            paths  : ``[path0, path1, path2]`` path bundle.
            params : must contain ``"super_cell_size"`` ``tuple[int,int,int]``.
        """
        self._load_m_scraps()

        if self.scraps_bin is None or not self.scraps_bin.exists():
            raise FileNotFoundError(
                f"SCRAPS binary not found at {self.scraps_bin}. "
                "Build it from the m_scraps repo: cd SCRAPs && make."
            )
        if self._mpirun is None:
            raise RuntimeError(
                "No mpirun found. Set the MPIRUN environment variable or install OpenMPI."
            )

        lattice = phase["lattice"]
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
                cwd=dir_name,
                capture_output=True,
                text=True,
                check=False,
            )
            print(result.stdout)
            if result.stderr:
                print("sqs2tdb stderr:", result.stderr)

        prim_atoms, _ = _build_prim_atoms(self.phases_dict)

        for sqsdir in sorted(dir_name.glob("sqsdb_lev=*/")):
            self._run_mscraps_in_dir(sqsdir, prim_atoms, cell_dim)

        self._write_objective_summary(dir_name)

    # ── Private helpers ─────────────────────────────────────────────

    def _write_endmember_bestsqs(
        self,
        sqsdir: Path,
        cell_dim: tuple[int, int, int],
    ) -> None:
        """Write bestsqs.out in prim-vector-unit format for a pure endmember.

        Pure-composition dirs don't need SCRAPS (atom order is trivial), but
        ``_poscar_from_bestsqs`` requires ``bestsqs.out`` in prim-vector format.

        Args:
            sqsdir   : ``sqsdb_lev=*`` directory.
            cell_dim : Supercell tile dimensions ``(n0, n1, n2)``.
        """
        prim_cell = self._prim_cell
        pos_list, site_labels = _parse_coords(self.phases_dict["coords"])
        n0, n1, n2 = int(cell_dim[0]), int(cell_dim[1]), int(cell_dim[2])
        sup_prim = np.diag([float(n0), float(n1), float(n2)])

        lines: list[str] = []
        for v in prim_cell:
            lines.append(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")
        for v in sup_prim:
            lines.append(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}")

        for i in range(n0):
            for j in range(n1):
                for k in range(n2):
                    for pos, lbl in zip(pos_list, site_labels, strict=False):
                        px, py, pz = pos[0] + i, pos[1] + j, pos[2] + k
                        atat_label = f"{lbl}_A" if (lbl.islower() and len(lbl) == 1) else lbl
                        lines.append(f"{px:.10f} {py:.10f} {pz:.10f} {atat_label}")

        (sqsdir / "bestsqs.out").write_text("\n".join(lines) + "\n")
        (sqsdir / "wait").touch()
        print(f"  Endmember bestsqs.out written: {sqsdir.name}")

    def _run_mscraps_in_dir(
        self,
        sqsdir: Path,
        prim_atoms: Atoms,
        cell_dim: tuple[int, int, int],
    ) -> None:
        """Run m_scraps in one sqsdb_lev directory and convert output to bestsqs.out.

        Intended as a single-directory worker called from ``sqs_gen``.

        Args:
            sqsdir     : ``sqsdb_lev=*`` directory.
            prim_atoms : ASE Atoms for the primitive cell (with placeholder species).
            cell_dim   : Supercell tile dimensions ``(n0, n1, n2)``.
        """
        sqsdb_comp = _parse_sqsdb_comp(sqsdir.name)

        all_fracs = [f for fracs in sqsdb_comp.values() for f in fracs]
        non_zero = [f for f in all_fracs if f > 0.0]
        if non_zero and all(f == 1.0 for f in non_zero):
            self._write_endmember_bestsqs(sqsdir, cell_dim)
            return

        if self.skip_existing and (sqsdir / "bestsqs.out").exists():
            print(f"Skipping m_scraps for {sqsdir.name}: bestsqs.out exists")
            (sqsdir / "wait").touch()
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

        n_basis = len(self._site_labels)
        cell_vol = cell_dim[0] * cell_dim[1] * cell_dim[2]
        expected = n_basis * cell_vol
        total = sum(cfg["element_counts"])
        if total != expected:
            print(
                f"WARNING: element_counts sum {total} != expected {expected} "
                f"in {sqsdir.name} — rounding mismatch, skipping"
            )
            return

        # MC kernel is passed as an extra INPUT line; not a direct make_scraps_inputs param.
        extra_lines: str = f"MC_KERNEL {self.mc_kernel}\n" if self.mc_kernel else ""

        work_dir = sqsdir / "_mscraps_run"
        work_dir.mkdir(exist_ok=True)

        figure_cutoffs = self.figure_cutoffs or None
        figure_kwargs: dict = {}
        if figure_cutoffs:
            figure_kwargs = {
                "figure_cutoffs": figure_cutoffs,
                "figure_species": self.figure_species or None,
                "figure_objective": self.figure_objective,
                "figure_wd": self.figure_wd,
                "figure_wn": self.figure_wn,
                "figure_wr": self.figure_wr,
                "figure_tol": self.figure_tol,
            }
        self._make_scraps_inputs(
            work_dir=work_dir,
            atoms=prim_atoms,
            element_counts=cfg["element_counts"],
            cell_dim=cell_dim,
            max_shellnum=self.max_shellnum,
            auto_budget=self.auto_budget,
            sublattices=cfg["sublattice_blocks"] or None,
            spectator_basis_species=cfg["spectator_pairs"] or None,
            sro_targets=self.sro_targets or None,
            extra_lines=extra_lines or None,
            **figure_kwargs,
        )

        scraps_dest = work_dir / "scraps"
        shutil.copy(self.scraps_bin, scraps_dest)
        scraps_dest.chmod(0o755)

        log_path = work_dir / "scraps.log"
        cmd = [self._mpirun, "--oversubscribe", "-np", str(self.ranks), "./scraps"]
        print(f"m_scraps {sqsdir.name}: {' '.join(cmd)}")
        t0 = time.time()
        with log_path.open("w") as fh:
            rc = subprocess.run(
                cmd, cwd=work_dir, stdout=fh, stderr=subprocess.STDOUT, check=False
            ).returncode
        elapsed = time.time() - t0

        scraps_vasp = work_dir / "SCRAPS.vasp"
        if rc != 0 or not scraps_vasp.exists():
            print(f"  m_scraps failed (rc={rc}) for {sqsdir.name} ({elapsed:.1f}s)")
            return

        fitness = None
        for line in log_path.read_text().splitlines():
            if line.startswith("FINAL fitness"):
                with contextlib.suppress(Exception):
                    fitness = float(line.split("=")[-1].strip())
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
        (sqsdir / "wait").touch()
        print(f"  bestsqs.out written → {sqsdir / 'bestsqs.out'}")

    def _write_objective_summary(self, parent_dir: Path) -> None:
        """Write a tab-separated fitness summary for all completed sqsdb_lev dirs.

        Args:
            parent_dir: Parent SQS directory containing ``sqsdb_lev=*`` subdirs.
        """
        output_file = parent_dir / "objective_functions.txt"
        lines = ["folder\tscraps_fitness"]
        for sqsdir in sorted(parent_dir.glob("sqsdb_lev=*/")):
            log = sqsdir / "scraps.log"
            fitness = None
            if log.exists():
                for line in log.read_text().splitlines():
                    if line.startswith("FINAL fitness"):
                        with contextlib.suppress(Exception):
                            fitness = float(line.split("=")[-1].strip())
                        break
            if fitness is not None:
                lines.append(f"{sqsdir.name}\t{fitness}")
        output_file.write_text("\n".join(lines))
