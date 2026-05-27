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
    sqs.sqs_gen(phase=phase_list[0], paths=paths, iter=1_000_000, params=mcsqs_params)
"""

from __future__ import annotations

import math
import re
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
        sqsgen_levels (list[dict]): Composition levels for ``sqsgen.in``.
            Each entry must have ``"level"``, ``"letter"``, and
            ``"compositions"`` keys.
        level (int): Highest sqsgen level to include (inclusive).
        len_comp (int): Number of elements in the target system.
        skip_existing_sqs (bool): Skip directories that already contain a
            ``bestcorr.out`` file.
    """

    def __init__(
        self,
        phases_dict: dict,
        sqsgen_levels: list[dict],
        level: int,
        len_comp: int,
        skip_existing_sqs: bool = False,
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

        sqsgen = self._build_sqsgen_text()
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
        dir_name = Path(paths[1]) / "SQS" / f"{phase['lattice']}_{self.len_comp}"

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

        parent_dir = Path(paths[1]) / "SQS" / f"{phase['lattice']}_{self.len_comp}"

        cutoff = BladeCutoff()
        lattice = cutoff.lattice_from_params(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
        frac = cutoff.read_coords(self.unit_cell)
        shells = cutoff.get_shells(lattice, frac, params["super_cell_size"])

        print("Neighbor shells:")
        for i, s in enumerate(shells[:10], 1):
            print(f"  {i}NN = {s:.6f}")

        cutoff_dict = {
            "-2": shells[params["2"] - 1],
            "-3": shells[params["3"] - 1],
            "-4": shells[params["4"] - 1],
        }
        print(
            f"Derived cutoffs: {cutoff_dict['-2']:.5f} (pairs), "
            f"{cutoff_dict['-3']:.5f} (triplets), "
            f"{cutoff_dict['-4']:.5f} (quadruplets)"
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
        folder = Path(paths[1]) / "SQS" / f"{specific_phase['lattice']}_{self.len_comp}"

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
        """Format the ``sqsgen.in`` content from the level definitions.

        Returns:
            str: Multi-line string suitable for writing to ``sqsgen.in``.
        """
        sqsgen = ""
        if self.level >= 1 and self.len_comp == 1:
            indices = [0]
        elif self.level >= 3 and self.len_comp == 2:
            indices = [0, 1, 2, 5]
        else:
            indices = list(range(self.level + 1))

        for i in indices:
            level_def = self.sqsgen_levels[i]
            level_num = level_def["level"]
            letters = level_def["letter"]
            for comp in level_def["compositions"]:
                comp_text = ",".join(map(str, comp))
                line = f"level={level_num}"
                for letter in letters:
                    line += f"\t\t{letter}={comp_text}"
                sqsgen += line + "\n"

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
        comp_str = folder_name.split("=")[-1]
        fractions = [float(x) for x in comp_str.split(",")]

        if fractions == [1.0, 0.0]:
            print(f"Skipping pure-species directory: {sqsdir}")
            return

        print(f"Running corrdump in {fractions}")
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
                    f"-3={cutoff_dict['-3']}",
                    f"-4={cutoff_dict['-4']}",
                ],
                cwd=sqsdir,
                check=True,
            )
        except subprocess.CalledProcessError:
            print(f"corrdump failed in {sqsdir}, skipping.")
            return

        print(f"Running mcsqs with {n_atoms} atoms in {fractions}")
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
            f"-3={cutoff_dict['-3']:.5f}",
            f"-4={cutoff_dict['-4']:.5f}",
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
