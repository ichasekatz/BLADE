"""SQS input generation and ATAT mcsqs execution.

This module provides :class:`BladeSQS`, which writes ATAT input files
(``rndstr.skel``, ``sqsgen.in``), runs ``sqs2tdb -mk`` to populate
``sqsdb_lev=*`` sub-directories, and then executes ``corrdump`` followed
by parallel ``mcsqs`` instances for each composition directory.

Cutoff distances for ``corrdump``/``mcsqs`` are derived automatically from
the lattice parameters using :class:`~blade.tools.blade_cutoff.BladeCutoff`.

Example::

    from blade.tools.blade_sqsgen import BladeSQS

    sqs = BladeSQS(
        phases_dict=phases["PHASE1"],
        sqsgen_levels=sqsgen_levels,
        level=5,
        len_comp=3,
        skip_existing_sqs=True,
    )
    sqs.sqs_gen(phase=phase_list[0], paths=paths, params=mcsqs_params)
"""

from __future__ import annotations

import math
import shutil
import subprocess
from pathlib import Path

from blade.tools._mcsqs import (
    _run_mcsqs_in_dir,
    _write_objective_summary,
    monitor_bestcorr,
    monitor_bestcorr_parallel,
    read_objective,
)
from blade.tools._sqsgen_text import _build_sqsgen_text
from blade.tools.blade_cutoff import BladeCutoff

__author__ = "Chase Katz"

__all__ = [
    "BladeSQS",
    "read_objective",
    "monitor_bestcorr",
    "monitor_bestcorr_parallel",
]

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

# How many neighbor shells to print during SQS setup (informational only).
_n_shells_to_print: int = 10


class BladeSQS:
    """Generate SQS inputs and run ATAT mcsqs for a phase prototype.

    Handles the full SQS generation sub-workflow:

    1. Build ``rndstr.skel`` and ``sqsgen.in`` from the phase prototype.
    2. Run ``sqs2tdb -mk`` to create ``sqsdb_lev=*`` sub-directories.
    3. Run ``corrdump`` then parallel ``mcsqs`` in each sub-directory,
       with an automatic timeout via a ``stopsqs`` sentinel file.
    4. Summarize objective-function values across all sub-directories.

    Attributes:
        phases_dict (dict): Phase prototype with keys ``"a"``, ``"b"``,
            ``"c"``, ``"alpha"``, ``"beta"``, ``"gamma"``, ``"vectors"``,
            and ``"coords"``.
        sqsgen_levels (list[dict]): Composition level seeds for ``sqsgen.in``.
            Each entry must have ``"level"``, ``"letter"``, and
            ``"compositions"`` keys.  Each seed composition is automatically
            expanded to all canonical (sorted-descending) compositions that
            share its denominator — e.g., ``[0.75, 0.25]`` also generates
            ``[0.5, 0.5]`` for a binary system.  Pure endmembers are not
            expanded.
        level (int): Highest sqsgen level to include (inclusive).
        len_comp (int): Number of elements in the target system.
        skip_existing_sqs (bool): Skip directories that already contain a
            ``bestcorr.out`` file.
        sqsgen_in (str | None): Optional verbatim ``sqsgen.in`` content.
            When set, ``_build_sqsgen_text`` is bypassed entirely.
        sublattice_map (dict[str, list[str]] | None): Per-sublattice active
            species lists.  When provided, compositions for constrained
            sublattices are zero-padded to ``len_comp`` so that ternary (or
            higher-order) sqsdb entries correctly fix inactive species at 0.
    """

    def __init__(
        self,
        phases_dict: dict,
        sqsgen_levels: list[dict],
        level: int,
        len_comp: int,
        skip_existing_sqs: bool = False,
        sublattice_map: dict[str, list[str]] | None = None,
        sqsgen_in: str | None = None,
        fixed_compositions: dict[str, list[float]] | None = None,
    ) -> None:
        """Initialize BladeSQS.

        Args:
            phases_dict (dict): Phase prototype dictionary. Required keys:
                ``"a"``, ``"b"``, ``"c"``, ``"alpha"``, ``"beta"``,
                ``"gamma"``, ``"vectors"`` (lattice vector string),
                ``"coords"`` (fractional coordinate string with sublattice
                labels).
            sqsgen_levels (list[dict]): Ordered list of composition level
                definitions for ``sqsgen.in``. Each entry must contain:

                - ``"level"`` (int): Level index.
                - ``"letter"`` (list[str]): Sublattice letter(s).
                - ``"compositions"`` (list[list[float]]): Fractional
                  compositions for that sublattice.

            level (int): Maximum sqsgen level to include (e.g., ``5`` will
                include levels ``[0, 1, 2, 3, 4, 5]``).
            len_comp (int): Number of elements in the chemical system.
                Controls which composition branches of ``sqsgen.in`` are
                written.
            skip_existing_sqs (bool, optional): If ``True``, skip SQS
                generation for directories that already exist, and skip
                ``mcsqs`` runs for directories that already have a
                ``bestcorr.out`` file. Defaults to ``False``.
            sqsgen_in (str | None, optional): Verbatim ``sqsgen.in`` content.
                When provided, bypasses ``_build_sqsgen_text`` entirely and
                writes this string directly to ``sqsgen.in``.  Useful for
                multi-sublattice phases that require hand-crafted composition
                lines.  Defaults to ``None``.
            sublattice_map (dict[str, list[str]] | None, optional): Maps
                sublattice letters to their active element lists.  When set,
                ``_build_sqsgen_text`` pads each constrained sublattice's
                compositions with trailing zeros up to ``len_comp``, producing
                sqsdb entries like ``a=0.5,0.5,0.0`` for a binary-on-sublattice
                phase in a ternary system.  Defaults to ``None``.
        """
        self.phases_dict = phases_dict
        self.a = phases_dict["a"]
        self.b = phases_dict["b"]
        self.c = phases_dict["c"]
        self.alpha = phases_dict["alpha"]
        self.beta = phases_dict["beta"]
        self.gamma = phases_dict["gamma"]
        self.vectors = phases_dict["vectors"]
        self.unit_cell = phases_dict["coords"]
        self.sqsgen_levels = sqsgen_levels
        self.level = level
        self.len_comp = len_comp
        self.skip_existing_sqs = skip_existing_sqs
        self.sublattice_map = sublattice_map
        self.sqsgen_in = sqsgen_in
        # Merge any sublattices marked "Constant" in sublattice_map into fixed_compositions
        # so they are excluded from cross-product permutation in sqsgen.in.
        merged_fixed = dict(fixed_compositions or {})
        if sublattice_map:
            for phase_map in sublattice_map.values():
                constant_letters = phase_map.get("Constant", [])
                for letter in constant_letters:
                    if letter not in merged_fixed:
                        merged_fixed[letter] = []  # placeholder; caller must also set fixed_compositions
        self.fixed_compositions = merged_fixed

        self.sqsgen_text: str = ""
        self.rndstr: str = ""

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def sqs_struct(self) -> tuple[str, str]:
        """Build the ``sqsgen.in`` text and the ``rndstr.skel`` text.

        Selects which composition levels to include based on :attr:`level`
        and :attr:`len_comp`, then formats the ATAT input files.

        Returns:
            tuple[str, str]: A two-element tuple ``(sqsgen_text, rndstr_text)``
            where ``sqsgen_text`` is the content for ``sqsgen.in`` and
            ``rndstr_text`` is the content for ``rndstr.skel``.
        """
        rndstr_header = f"{self.a} {self.b} {self.c} {self.alpha} {self.beta} {self.gamma}\n{self.vectors.strip()}"
        print(rndstr_header)

        sqsgen = (
            self.sqsgen_in
            if self.sqsgen_in is not None
            else _build_sqsgen_text(
                self.sqsgen_levels,
                self.level,
                self.len_comp,
                self.unit_cell,
                self.fixed_compositions,
            )
        )
        rndstr = rndstr_header.strip() + "\n" + self.unit_cell.strip()
        print(rndstr)

        self.sqsgen_text = sqsgen
        self.rndstr = rndstr
        return sqsgen, rndstr

    def sqs_gen(
        self,
        phase: dict,
        paths: list[str | Path],
        params: dict,
    ) -> None:
        """Generate SQS folders and run ``corrdump`` + ``mcsqs``.

        For each ``sqsdb_lev=*`` sub-directory created by ``sqs2tdb -mk``:

        1. Computes cutoff distances from lattice parameters.
        2. Runs ``corrdump`` to generate cluster correlations.
        3. Spawns ``params["parallel_runs"]`` parallel ``mcsqs`` processes.
        4. Stops them after ``params["time"]`` seconds by writing a
           ``stopsqs`` sentinel file.
        5. Writes ``objective_functions.txt`` summarizing all runs.

        Args:
            phase (dict): Single phase entry from the phase list. Must contain
                a ``"lattice"`` key.
            paths (list[str | Path]): Three-element path bundle
                (see :class:`BladeTDBGen` for the convention).
            params (dict): ``mcsqs`` run parameters. Required keys:

                - ``"super_cell_size"`` (tuple[int, int, int])
                - ``"parallel_runs"`` (int)
                - ``"time"`` (float): positive run duration in seconds
                - ``"2"``, ``"3"``, ``"4"`` (int): neighbor-shell indices
                - ``"wr"``, ``"wn"``, ``"wd"`` (float): mcsqs weights
        """
        if float(params.get("time", 0)) <= 0:
            raise ValueError('params["time"] must be positive')

        dir_name = Path(paths[1]) / "Files" / "SQS" / f"{phase['lattice']}_{self.len_comp}"

        if not self.skip_existing_sqs and dir_name.exists():
            print(f"Removing existing SQS directory: {dir_name}")
            shutil.rmtree(dir_name)

        if self.skip_existing_sqs and dir_name.exists():
            print(f"Skipping SQS generation for {phase['lattice']}_{self.len_comp}: found existing folder at {dir_name}")
        else:
            dir_name.mkdir(parents=True, exist_ok=True)
            sqsgen, rndstr = self.sqs_struct()

            (dir_name / "rndstr.skel").write_text(rndstr)
            print(f"File created at: {dir_name / 'rndstr.skel'}")

            (dir_name / "sqsgen.in").write_text(sqsgen)
            print(f"File created at: {dir_name / 'sqsgen.in'}")

            result = subprocess.run(
                ["sqs2tdb", "-mk"],
                cwd=dir_name,
                capture_output=True,
                text=True,
                check=False,
            )
            print(result.stdout)
            if result.stderr:
                print("Error:", result.stderr)

        parent_dir = Path(paths[1]) / "Files" / "SQS" / f"{phase['lattice']}_{self.len_comp}"

        cutoff = BladeCutoff()
        lattice = cutoff.lattice_from_params(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        frac = cutoff.read_coords(self.unit_cell)
        shells = cutoff.get_shells(lattice, frac, params["super_cell_size"])

        print("Neighbor shells:")
        for i, s in enumerate(shells[:_n_shells_to_print], 1):
            print(f"  {i}NN = {s:.6f} Å")
        print(f"Bond length (1NN): {shells[0]:.6f} Å")

        cutoff_mode = params.get("cutoff_mode", "nn")

        def _resolve_cutoff(n: float) -> float:
            if cutoff_mode == "angstrom":
                return float(n)
            lo = int(n) - 1
            frac = n - int(n)
            if frac == 0:
                return float(shells[lo])
            return float(shells[lo] + frac * (shells[lo + 1] - shells[lo]))

        cutoff_dict: dict[str, float] = {
            "-2": _resolve_cutoff(params["2"]),
        }
        if params.get("3", 0):
            cutoff_dict["-3"] = _resolve_cutoff(params["3"])
        if params.get("4", 0):
            cutoff_dict["-4"] = _resolve_cutoff(params["4"])
        print(
            f"Derived cutoffs: {cutoff_dict['-2']:.5f} (pairs)"
            + (f", {cutoff_dict['-3']:.5f} (triplets)" if "-3" in cutoff_dict else "")
            + (f", {cutoff_dict['-4']:.5f} (quadruplets)" if "-4" in cutoff_dict else "")
        )

        n_atoms = len(self.unit_cell.strip().splitlines()) * math.prod(params["super_cell_size"])

        for sqsdir in parent_dir.glob("sqsdb_lev=*/"):
            _run_mcsqs_in_dir(sqsdir, n_atoms, cutoff_dict, params, self.skip_existing_sqs)

        _write_objective_summary(parent_dir)

    def rename_files(
        self,
        specific_phase: dict,
        paths: list[str | Path],
        sqsgen_levels2: list[dict],
    ) -> None:
        """Rename ``sqsdb_lev=*`` folders to include fixed-sublattice labels.

        Also appends the corresponding sublattice composition to each line
        of ``sqsgen.in`` so that future ``sqs2tdb`` calls include the fixed
        species.

        Args:
            specific_phase (dict): Phase entry (must contain ``"lattice"`` key).
            paths (list[str | Path]): Three-element path bundle.
            sqsgen_levels2 (list[dict]): Fixed-sublattice definitions. Each
                entry must have ``"letter"`` and ``"compositions"`` keys.
        """
        folder = Path(paths[1]) / "Files" / "SQS" / f"{specific_phase['lattice']}_{self.len_comp}"

        for level_def in sqsgen_levels2:
            letter = level_def["letter"]
            compositions = level_def["compositions"]
            suffix = f"_{letter}={compositions}"

            for sqsdir in folder.glob("sqsdb_lev=*/"):
                if sqsdir.is_dir() and suffix not in sqsdir.name:
                    new_path = sqsdir.parent / f"{sqsdir.name}{suffix}"
                    sqsdir.rename(new_path)
                    print(f"Renamed {sqsdir} -> {new_path}")

            sqsgen_path = folder / "sqsgen.in"
            if sqsgen_path.exists():
                lines = sqsgen_path.read_text().splitlines()
                new_lines = [
                    line + f"\t\t{letter}={compositions}" if line.strip() and suffix not in line else line for line in lines
                ]
                sqsgen_path.write_text("\n".join(new_lines) + "\n")
                print(f"Updated {sqsgen_path}")
