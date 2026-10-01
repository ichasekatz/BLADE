"""Low-level mcsqs execution helpers for BLADE SQS generation.

Standalone functions extracted from :class:`~blade.tools.blade_sqsgen.BladeSQS`
that handle objective monitoring, corrdump/mcsqs subprocess management, and
objective-function summarisation.  Imported by ``blade_sqsgen`` — do not call
directly unless you know what you are doing.
"""

from __future__ import annotations

import re
import subprocess
import threading
import time
from pathlib import Path

__author__ = "Chase Katz"

# Poll cadence (seconds) for the mcsqs timeout loop in _wait_with_timeout.
_mcsqs_poll_interval: int = 1


def read_objective(bestcorr_path: str | Path) -> float | None:
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
    if re.search(r"Objective_function\s*=\s*Perfect_match", text, re.IGNORECASE):
        return float("-inf")
    match = re.search(
        r"Objective_function\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
        text,
    )
    return float(match.group(1)) if match else None


def monitor_bestcorr(
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

    with history_path.open("w") as f:
        f.write("time_seconds\tobjective\n")
        while not stop_event.is_set():
            objective = read_objective(bestcorr_path)
            if objective is not None and objective != last_objective:
                elapsed = time.time() - start_time
                print(f"{Path(sqsdir).name}: time={elapsed:.1f}s objective={objective}")
                f.write(f"{elapsed:.2f}\t{objective}\n")
                f.flush()
                last_objective = objective
            time.sleep(interval)


def monitor_bestcorr_parallel(
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

    with history_path.open("w") as f:
        f.write("time_seconds\tid\tobjective\n")
        while not stop_event.is_set():
            for bestcorr_path in sorted(sqsdir.glob("bestcorr*.out")):
                objective = read_objective(bestcorr_path)
                if objective is None:
                    continue
                stem = bestcorr_path.stem.replace("bestcorr", "")
                file_id = stem if stem else "main"
                key = bestcorr_path.name
                if last_objectives.get(key) != objective:
                    elapsed = time.time() - start_time
                    print(f"{sqsdir.name}: id={file_id} time={elapsed:.1f}s objective={objective}")
                    f.write(f"{elapsed:.2f}\t{file_id}\t{objective}\t")
                    f.flush()
                    last_objectives[key] = objective
            time.sleep(interval)


def _run_mcsqs_in_dir(
    sqsdir: Path,
    n_atoms: int,
    cutoff_dict: dict[str, float],
    params: dict,
    skip_existing_sqs: bool = False,
) -> None:
    """Run ``corrdump`` then parallel ``mcsqs`` in a single sqsdb directory.

    Args:
        sqsdir (Path): ``sqsdb_lev=*`` sub-directory.
        n_atoms (int): Supercell size argument for ``mcsqs -n``.
        cutoff_dict (dict[str, float]): Cutoff distances keyed by
            ``"-2"``, ``"-3"``, ``"-4"``.
        params (dict): mcsqs run parameters (see
            :meth:`~blade.tools.blade_sqsgen.BladeSQS.sqs_gen`).
        skip_existing_sqs (bool, optional): If ``True``, skip directories
            that already contain a ``bestcorr.out`` file.  Defaults to
            ``False``.
    """
    if skip_existing_sqs and (sqsdir / "bestcorr.out").exists():
        print(f"Skipping mcsqs for {sqsdir}: bestcorr.out already exists.")
        return

    folder_name = sqsdir.name
    sublattice_fracs = re.findall(r"_([a-z])=([\d.,]+)", folder_name)
    all_fracs = [float(x) for _, comp_str in sublattice_fracs for x in comp_str.split(",")]
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
        target=monitor_bestcorr_parallel,
        args=(sqsdir, stop_monitor, 5),
        daemon=True,
    )
    monitor_thread.start()

    processes = [_spawn_mcsqs(sqsdir, ip, n_atoms, cutoff_dict, params) for ip in range(1, n_parallel + 1)]

    try:
        _wait_with_timeout(sqsdir, processes, stopsqs_path, params)
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
    sqsdir: Path,
    ip: int,
    n_atoms: int,
    cutoff_dict: dict[str, float],
    params: dict,
) -> tuple[int, subprocess.Popen]:
    """Start a single ``mcsqs -ip=N`` process.

    Args:
        sqsdir (Path): Working directory for the process.
        ip (int): Parallel instance index.
        n_atoms (int): Supercell size.
        cutoff_dict (dict[str, float]): Cutoff distances.
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
    print(f"Starting mcsqs -ip={ip} in {sqsdir}")
    p = subprocess.Popen(cmd, cwd=sqsdir)
    return ip, p


def _wait_with_timeout(
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
        time.sleep(_mcsqs_poll_interval)

    if not stopsqs_path.exists():
        stopsqs_path.touch()
        print(f"Wrote stopsqs at {stopsqs_path} after {params['time']}s")

    for ip, p in processes:
        if p.poll() is None:
            p.kill()
            print(f"Killed mcsqs -ip={ip} in {sqsdir}")


def _write_objective_summary(parent_dir: Path) -> None:
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
            objective = read_objective(candidate)
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
