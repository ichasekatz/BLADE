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
        phases_dict=phases["HEDB1"],
        sqsgen_levels=sqsgen_levels,
        level=5,
        len_comp=3,
        skip_existing_sqs=True,
    )
    sqs.sqs_gen(phase=phase_list[0], paths=paths, iter=1_000_000, params=mcsqs_params)
"""

from __future__ import annotations

import math
import re
import shutil
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import gcd
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import TYPE_CHECKING

from blade.tools.blade_cutoff import BladeCutoff

if TYPE_CHECKING:
    pass

__author__ = "Chase Katz"


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
                        merged_fixed[letter] = []   # placeholder; caller must also set fixed_compositions
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
        rndstr_header = (
            f"{self.a} {self.b} {self.c} {self.alpha} {self.beta} {self.gamma}\n"
            f"{self.vectors.strip()}"
        )
        print(rndstr_header)

        sqsgen = self.sqsgen_in if self.sqsgen_in is not None else self._build_sqsgen_text()
        rndstr = rndstr_header.strip() + "\n" + self.unit_cell.strip()
        print(rndstr)

        self.sqsgen_text = sqsgen
        self.rndstr = rndstr
        return sqsgen, rndstr

    def read_objective(self, bestcorr_path: str | Path) -> float | None:
        """Parse the objective-function value from a ``bestcorr.out`` file.

        Args:
            bestcorr_path (str | Path): Path to a ``bestcorr.out`` file.

        Returns:
            float | None: The objective-function value, or ``None`` if the
            file does not exist or the value cannot be parsed.
        """
        bestcorr_path = Path(bestcorr_path)
        if not bestcorr_path.exists():
            return None
        text = bestcorr_path.read_text(errors="ignore")
        match = re.search(
            r"Objective_function\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
            text,
        )
        return float(match.group(1)) if match else None

    def monitor_bestcorr(
        self,
        sqsdir: str | Path,
        stop_event: threading.Event,
        interval: int = 5,
    ) -> None:
        """Monitor a single ``bestcorr.out`` file and log objective changes.

        Runs in a background thread. Polls ``bestcorr.out`` every
        ``interval`` seconds and writes new objective values to
        ``objective_history.txt`` in the same directory.

        Args:
            sqsdir (str | Path): Directory containing ``bestcorr.out``.
            stop_event (threading.Event): Setting this event stops the monitor.
            interval (int, optional): Poll interval in seconds.
                Defaults to ``5``.
        """
        bestcorr_path = Path(sqsdir) / "bestcorr.out"
        history_path = Path(sqsdir) / "objective_history.txt"
        start_time = time.time()
        last_objective = None

        with open(history_path, "w") as f:
            f.write("time_seconds\tobjective\n")
            while not stop_event.is_set():
                objective = self.read_objective(bestcorr_path)
                if objective is not None and objective != last_objective:
                    elapsed = time.time() - start_time
                    print(f"{Path(sqsdir).name}: time={elapsed:.1f}s objective={objective}")
                    f.write(f"{elapsed:.2f}\t{objective}\n")
                    f.flush()
                    last_objective = objective
                time.sleep(interval)

    def monitor_bestcorr_parallel(
        self,
        sqsdir: str | Path,
        stop_event: threading.Event,
        interval: int = 5,
    ) -> None:
        """Monitor all ``bestcorr*.out`` files in a parallel-run directory.

        Runs in a background thread. Polls every ``bestcorr*.out`` file found
        in ``sqsdir`` and logs new objective values to
        ``objective_history.txt``.

        Args:
            sqsdir (str | Path): Directory containing ``bestcorr*.out`` files
                from parallel ``mcsqs -ip=N`` runs.
            stop_event (threading.Event): Setting this event stops the monitor.
            interval (int, optional): Poll interval in seconds.
                Defaults to ``5``.
        """
        sqsdir = Path(sqsdir)
        history_path = sqsdir / "objective_history.txt"
        start_time = time.time()
        last_objectives: dict[str, float] = {}

        with open(history_path, "w") as f:
            f.write("time_seconds\tid\tobjective\n")
            while not stop_event.is_set():
                for bestcorr_path in sorted(sqsdir.glob("bestcorr*.out")):
                    objective = self.read_objective(bestcorr_path)
                    if objective is None:
                        continue
                    stem = bestcorr_path.stem.replace("bestcorr", "")
                    file_id = stem if stem else "main"
                    key = bestcorr_path.name
                    if last_objectives.get(key) != objective:
                        elapsed = time.time() - start_time
                        print(
                            f"{sqsdir.name}: id={file_id} "
                            f"time={elapsed:.1f}s objective={objective}"
                        )
                        f.write(f"{elapsed:.2f}\t{file_id}\t{objective}\t")
                        f.flush()
                        last_objectives[key] = objective
                time.sleep(interval)

    def sqs_gen(
        self,
        phase: dict,
        paths: list[str | Path],
        iter: float,
        params: dict,
    ) -> None:
        """Generate SQS folders and run ``corrdump`` + ``mcsqs``.

        For each ``sqsdb_lev=*`` sub-directory created by ``sqs2tdb -mk``:

        1. Computes cutoff distances from lattice parameters.
        2. Runs ``corrdump`` to generate cluster correlations.
        3. Spawns ``params["parallel_runs"]`` parallel ``mcsqs`` processes.
        4. Stops them after ``params["time"]`` seconds (if ``use_time=True``)
           or after ``iter`` steps, by writing a ``stopsqs`` sentinel file.
        5. Writes ``objective_functions.txt`` summarizing all runs.

        Args:
            phase (dict): Single phase entry from the phase list. Must contain
                a ``"lattice"`` key.
            paths (list[str | Path]): Three-element path bundle
                (see :class:`BladeTDBGen` for the convention).
            iter (float): Maximum ``mcsqs`` iteration count (used when
                ``params["use_time"]`` is ``False``).
            params (dict): ``mcsqs`` run parameters. Required keys:

                - ``"super_cell_size"`` (tuple[int, int, int])
                - ``"parallel_runs"`` (int)
                - ``"use_time"`` (bool)
                - ``"time"`` (float): seconds (when ``use_time=True``)
                - ``"2"``, ``"3"``, ``"4"`` (int): neighbor-shell indices
                - ``"wr"``, ``"wn"``, ``"wd"`` (float): mcsqs weights
        """
        dir_name = Path(paths[1]) / "Files" / "SQS" / f"{phase['lattice']}_{self.len_comp}"

        if not self.skip_existing_sqs and dir_name.exists():
            print(f"Removing existing SQS directory: {dir_name}")
            shutil.rmtree(dir_name)

        if self.skip_existing_sqs and dir_name.exists():
            print(
                f"Skipping SQS generation for {phase['lattice']}_{self.len_comp}: "
                f"found existing folder at {dir_name}"
            )
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
        lattice = cutoff.lattice_from_params(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
        frac = cutoff.read_coords(self.unit_cell)
        shells = cutoff.get_shells(lattice, frac, params["super_cell_size"])

        print("Neighbor shells:")
        for i, s in enumerate(shells[:10], 1):
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

        n_atoms = (
            len(self.unit_cell.strip().splitlines())
            * math.prod(params["super_cell_size"])
        )

        for sqsdir in parent_dir.glob("sqsdb_lev=*/"):
            self._run_mcsqs_in_dir(sqsdir, n_atoms, cutoff_dict, iter, params)

        self._write_objective_summary(parent_dir)

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
                    line + f"\t\t{letter}={compositions}"
                    if line.strip() and suffix not in line
                    else line
                    for line in lines
                ]
                sqsgen_path.write_text("\n".join(new_lines) + "\n")
                print(f"Updated {sqsgen_path}")

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _build_sqsgen_text(self) -> str:
        def _trim_comp(vals: list[float]) -> list[float]:
            vals = list(vals)
            vals = vals + [0.0] * max(0, self.len_comp - len(vals))
            while len(vals) > 1 and vals[-1] == 0.0:
                vals.pop()
            return vals

        def _fmt_comp(vals: list[float]) -> str:
            return ",".join(map(str, vals))

        def _is_pure_one(vals: list[float]) -> bool:
            return vals == [1.0]

        def _lcm(a: int, b: int) -> int:
            return a * b // gcd(a, b)

        def _composition_denominator(vals: list[float]) -> int:
            denom = 1
            for f in vals:
                frac = Fraction(float(f)).limit_denominator(100)
                denom = _lcm(denom, frac.denominator)
            return denom

        var_letters: set[str] = set()

        for line in self.unit_cell.strip().splitlines():
            parts = line.strip().split()
            if len(parts) >= 4:
                label = parts[3]
                if label.islower() and len(label) == 1:
                    var_letters.add(label)

        sqsgen = ""

        if self.level >= 1 and self.len_comp == 1:
            indices = [0]
        elif self.level >= 3 and self.len_comp == 2:
            indices = [0, 1, 2, 5]
        else:
            indices = list(range(self.level + 1))

        best_level_for_comp: dict[tuple[float, ...], int] = {}

        for i in indices:
            level_def = self.sqsgen_levels[i]
            level_num = level_def["level"]

            for comp in level_def["compositions"]:
                non_zero = [float(f) for f in comp if float(f) > 0.0]

                if non_zero == [1.0]:
                    vals = _trim_comp(list(comp))
                    key = tuple(vals)

                    if key not in best_level_for_comp or level_num < best_level_for_comp[key]:
                        best_level_for_comp[key] = level_num

                    continue

                vals = _trim_comp([float(f) for f in comp])
                key = tuple(vals)

                if key not in best_level_for_comp or level_num < best_level_for_comp[key]:
                    best_level_for_comp[key] = level_num

        comp_entries: list[tuple[int, list[float]]] = [
            (level_num, list(comp))
            for comp, level_num in best_level_for_comp.items()
        ]

        comp_entries = sorted(comp_entries, key=lambda x: (x[0], x[1]))

        all_letters: list[str] = []

        for i in indices:
            level_def = self.sqsgen_levels[i]
            for letter in level_def["letter"]:
                # Skip letters with a fixed composition — they are appended separately
                if letter in self.fixed_compositions:
                    continue
                if letter in var_letters and letter not in all_letters:
                    all_letters.append(letter)

        if not all_letters:
            # All variable sublattices are fixed — emit a single level=0 line
            # so sqs2tdb -mk can create the sqsdb directory.
            if self.fixed_compositions:
                fixed_suffix = "".join(
                    f"\t\t{letter}={_fmt_comp(comp)}"
                    for letter, comp in self.fixed_compositions.items()
                )
                return f"level=0{fixed_suffix}\n"
            return ""

        if len(all_letters) == 1:
            letter = all_letters[0]

            for level_num, vals in comp_entries:
                if level_num > self.level:
                    continue

                line = f"level={level_num}\t\t{letter}={_fmt_comp(vals)}"
                if self.fixed_compositions:
                    line += "".join(
                        f"\t\t{fl}={_fmt_comp(fc)}"
                        for fl, fc in self.fixed_compositions.items()
                    )
                sqsgen += line + "\n"

            return sqsgen

        endmember_line = "level=0"
        for letter in all_letters:
            endmember_line += f"\t\t{letter}=1.0"
        sqsgen += endmember_line + "\n"

        non_endmember_entries = [
            (level_num, vals)
            for level_num, vals in comp_entries
            if level_num > 0 and not _is_pure_one(vals)
        ]

        written: set[str] = set()

        for level_num, vals in non_endmember_entries:
            if level_num > self.level:
                continue

            for active_letter in all_letters:
                line = f"level={level_num}"

                for letter in all_letters:
                    if letter == active_letter:
                        line += f"\t\t{letter}={_fmt_comp(vals)}"
                    else:
                        line += f"\t\t{letter}=1.0"

                if line not in written:
                    sqsgen += line + "\n"
                    written.add(line)

        for k in range(2, len(all_letters) + 1):
            for active_letters in combinations(all_letters, k):
                for combo in product(non_endmember_entries, repeat=k):
                    combo_levels = [entry[0] for entry in combo]
                    combo_vals = [entry[1] for entry in combo]

                    combined_level = max(combo_levels)

                    if combined_level > self.level:
                        continue

                    line = f"level={combined_level}"

                    for letter in all_letters:
                        if letter in active_letters:
                            idx = active_letters.index(letter)
                            line += f"\t\t{letter}={_fmt_comp(combo_vals[idx])}"
                        else:
                            line += f"\t\t{letter}=1.0"

                    if line not in written:
                        sqsgen += line + "\n"
                        written.add(line)

        # Append fixed-composition sublattices to every generated line
        if self.fixed_compositions:
            fixed_suffix = "".join(
                f"\t\t{letter}={_fmt_comp(comp)}"
                for letter, comp in self.fixed_compositions.items()
            )
            sqsgen = "\n".join(
                line + fixed_suffix if line.strip() else line
                for line in sqsgen.splitlines()
            ) + "\n"

        return sqsgen

    def _run_mcsqs_in_dir(
        self,
        sqsdir: Path,
        n_atoms: int,
        cutoff_dict: dict[str, float],
        iter: float,
        params: dict,
    ) -> None:
        """Run ``corrdump`` then parallel ``mcsqs`` in a single sqsdb directory.

        Args:
            sqsdir (Path): ``sqsdb_lev=*`` sub-directory.
            n_atoms (int): Supercell size argument for ``mcsqs -n``.
            cutoff_dict (dict[str, float]): Cutoff distances keyed by
                ``"-2"``, ``"-3"``, ``"-4"``.
            iter (float): Maximum mcsqs iteration count
                (used when ``use_time=False``).
            params (dict): mcsqs run parameters (see :meth:`sqs_gen`).
        """
        if self.skip_existing_sqs and (sqsdir / "bestcorr.out").exists():
            print(f"Skipping mcsqs for {sqsdir}: bestcorr.out already exists.")
            return

        folder_name = sqsdir.name
        sublattice_fracs = re.findall(r'_([a-z])=([\d.,]+)', folder_name)
        all_fracs = [
            float(x)
            for _, comp_str in sublattice_fracs
            for x in comp_str.split(",")
        ]
        non_zero = [f for f in all_fracs if f > 0.0]
        if non_zero and all(f == 1.0 for f in non_zero):
            print(f"Skipping pure-species directory: {sqsdir}")
            return

        print(f"Running corrdump in {folder_name}")
        try:
            subprocess.run(
                [
                    "corrdump",
                    "-l=rndstr.in",
                    "-ro",
                    "-noe",
                    "-nop",
                    "-clus",
                    f"-2={cutoff_dict['-2']}",
                    *([f"-3={cutoff_dict['-3']}"] if "-3" in cutoff_dict else []),
                    *([f"-4={cutoff_dict['-4']}"] if "-4" in cutoff_dict else []),
                ],
                cwd=sqsdir,
                check=True,
            )
        except subprocess.CalledProcessError:
            print(f"corrdump failed in {sqsdir}, skipping.")
            return

        print(f"Running mcsqs with {n_atoms} atoms in {folder_name}")
        n_parallel = params["parallel_runs"]
        stopsqs_path = sqsdir / "stopsqs"
        if stopsqs_path.exists():
            stopsqs_path.unlink()

        stop_monitor = threading.Event()
        monitor_thread = threading.Thread(
            target=self.monitor_bestcorr_parallel,
            args=(sqsdir, stop_monitor, 5),
            daemon=True,
        )
        monitor_thread.start()

        processes = [
            self._spawn_mcsqs(sqsdir, ip, n_atoms, cutoff_dict, iter, params)
            for ip in range(1, n_parallel + 1)
        ]

        try:
            if params["use_time"]:
                self._wait_with_timeout(sqsdir, processes, stopsqs_path, params)
            for ip, p in processes:
                ret = p.wait()
                status = "OK" if ret == 0 else f"exit code {ret}"
                print(f"mcsqs -ip={ip} finished ({status}) in {sqsdir}")

            bestcorr_files = list(sqsdir.glob("bestcorr*.out"))
            if bestcorr_files:
                print(f"Running mcsqs -best in {sqsdir}")
                subprocess.run(["mcsqs", "-best"], cwd=sqsdir, check=True)
            else:
                print(f"No bestcorr*.out files in {sqsdir}; skipping mcsqs -best")

        except subprocess.CalledProcessError:
            print(f"mcsqs -best failed in {sqsdir}, continuing.")
        finally:
            stop_monitor.set()
            monitor_thread.join()
            if stopsqs_path.exists():
                stopsqs_path.unlink()

    def _spawn_mcsqs(
        self,
        sqsdir: Path,
        ip: int,
        n_atoms: int,
        cutoff_dict: dict[str, float],
        iter: float,
        params: dict,
    ) -> tuple[int, subprocess.Popen]:
        """Start a single ``mcsqs -ip=N`` process.

        Args:
            sqsdir (Path): Working directory for the process.
            ip (int): Parallel instance index.
            n_atoms (int): Supercell size.
            cutoff_dict (dict[str, float]): Cutoff distances.
            iter (float): Max iterations (when ``use_time=False``).
            params (dict): mcsqs parameters.

        Returns:
            tuple[int, subprocess.Popen]: ``(ip, process)`` pair.
        """
        cmd = [
            "mcsqs",
            f"-n={n_atoms}",
            "-l=rndstr.in",
            f"-2={cutoff_dict['-2']:.5f}",
            *([f"-3={cutoff_dict['-3']:.5f}"] if "-3" in cutoff_dict else []),
            *([f"-4={cutoff_dict['-4']:.5f}"] if "-4" in cutoff_dict else []),
            f"-wr={params['wr']}",
            f"-wn={params['wn']}",
            f"-wd={params['wd']}",
            f"-ip={ip}",
        ]
        if not params["use_time"]:
            cmd.append(f"-ms={iter}")

        print(f"Starting mcsqs -ip={ip} in {sqsdir}")
        p = subprocess.Popen(cmd, cwd=sqsdir)
        return ip, p

    def _wait_with_timeout(
        self,
        sqsdir: Path,
        processes: list[tuple[int, subprocess.Popen]],
        stopsqs_path: Path,
        params: dict,
    ) -> None:
        """Wait for time-limited mcsqs runs and stop them via stopsqs.

        Args:
            sqsdir (Path): Directory where ``stopsqs`` is written.
            processes (list[tuple[int, subprocess.Popen]]): Running processes.
            stopsqs_path (Path): Path to write the ``stopsqs`` sentinel file.
            params (dict): Must contain ``"time"``.
        """
        start_time = time.time()
        while time.time() - start_time < params["time"]:
            if stopsqs_path.exists():
                print(f"Detected existing stopsqs at {stopsqs_path}; stopping early")
                break
            if all(p.poll() is not None for _, p in processes):
                print("All mcsqs processes finished before time limit")
                break
            time.sleep(1)

        if not stopsqs_path.exists():
            stopsqs_path.touch()
            print(f"Wrote stopsqs at {stopsqs_path} after {params['time']}s")

        for ip, p in processes:
            if p.poll() is None:
                p.kill()
                print(f"Killed mcsqs -ip={ip} in {sqsdir}")

    def _write_objective_summary(self, parent_dir: Path) -> None:
        """Write ``objective_functions.txt`` summarizing all sqsdb runs.

        Args:
            parent_dir (Path): Directory containing ``sqsdb_lev=*`` folders.
        """
        output_file = parent_dir / "objective_functions.txt"
        lines = ["folder\tbestcorr_path\tobjective"]

        for sqsdir in parent_dir.glob("sqsdb_lev=*/"):
            best_path = None
            objective = None
            for candidate in sqsdir.rglob("bestcorr.out"):
                objective = self.read_objective(candidate)
                best_path = candidate
                if objective is None:
                    print(f"Could not parse objective in {candidate}")

            if best_path is None:
                print(f"No bestcorr.out found in {sqsdir}")
            elif objective is not None:
                folder = best_path.parent.name
                print(f"{folder}: {objective}")
                lines.append(f"{folder}\t{best_path}\t{objective}")

        output_file.write_text("\n".join(lines))
