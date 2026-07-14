"""TDB fitting driver for high-throughput CALPHAD database generation.

This module provides :class:`BladeTDBGen`, which orchestrates
MaterialsFramework's :class:`~materialsframework.tools.sqs2tdb.Sqs2tdb`
workflow across many chemical compositions.  The class stores configuration
in :meth:`__init__` and executes the fitting loop only when :meth:`fit` is
called — no side effects on construction.

Example::

    from blade.tools.blade_tdb_gen import BladeTDBGen

    gen = BladeTDBGen(
        phases=phase_list,
        liquid=False,
        paths=[path0, path1, path2],
        composition_list=composition_list,
        level=5,
        skip_existing=True,
        terms_in={"PHASE1": "1,0\n2,0\n2,1\n2,2\n3,0\n"},
    )
    gen.fit()
"""

from __future__ import annotations

import os
import shutil
import threading
import time
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

__author__ = "Chase Katz"


class BladeTDBGen:
    """Run Sqs2tdb fitting across many chemical compositions.

    For each composition this class:

    1. Copies SQS directories from the BLADE staging tree into the
       MaterialsFramework ``sqsdb`` tree.
    2. Creates a per-composition output directory.
    3. Invokes :class:`~materialsframework.tools.sqs2tdb.Sqs2tdb` with a
       :class:`~materialsframework.calculators.GraceCalculator` backend.

    Construction is side-effect-free — the fitting loop only runs when
    :meth:`fit` is called explicitly.

    Attributes:
        phases (list[dict]): Phase prototype definitions (must include
            ``"lattice"`` key).
        liquid (bool): Whether to include the LIQUID phase in each fit.
        path0 (Path): Root working directory (restored after each composition).
        path1 (Path): BLADE staging directory containing generated SQS trees.
        path2 (Path): MaterialsFramework ``sqsdb`` base directory.
        composition_list (list[list[str]]): Chemical systems to process.
        level (int): SQS level passed to
            :meth:`~materialsframework.tools.sqs2tdb.Sqs2tdb.fit`.
        sqsgen_levels2 (list[dict]): Fixed-sublattice species definitions
            (e.g., a fixed sublattice in multi-sublattice systems).
            Auto-derived from phase coords if not provided.
        skip_existing (bool): If ``True``, skip compositions that already
            have a ``.tdb`` file in their output directory.
        terms_in (dict[str, str] | None): Per-phase ``terms.in`` content,
            keyed by lattice base name (e.g., ``"PHASE1"``).  Phases absent
            from the dict use the ATAT-generated template unchanged.
        sublattice_map (dict[str, dict[str, list[str]]] | None): Per-phase
            sublattice species assignments.  Outer key is lattice base name;
            inner key is sublattice letter; value is list of element symbols.
            Phases absent from the dict use the ATAT default.
    """

    def __init__(
        self,
        phases: list[dict],
        liquid: bool,
        paths: Sequence[str | Path],
        composition_list: list[list[str]],
        level: int,
        sqsgen_levels2: list[dict] | None = None,
        skip_existing: bool = False,
        phases_dict: dict[str, dict] | None = None,
        output_dir: str | Path | None = None,
        terms_in: dict[str, str] | None = None,
        mult_in: dict[str, str] | None = None,
        sublattice_map: dict[str, dict[str, list[str]]] | None = None,
        tdb_params: dict | None = None,
        fixed_compositions: dict[str, list[float]] | None = None,
    ) -> None:
        """Initialize BladeTDBGen.

        Args:
            phases (list[dict]): Phase prototype definitions. Each dict must
                contain at least a ``"lattice"`` key
                (e.g., ``[{"lattice": "PHASE1", ...}]``).
            liquid (bool): If ``True``, append ``"LIQUID"`` to every phase
                list passed to Sqs2tdb.
            paths (Sequence[str | Path]): Three-element path bundle:

                - ``paths[0]``: Root directory; the working directory is
                  restored here after each composition.
                - ``paths[1]``: BLADE staging directory (contains
                  ``SQS/<lattice>_<n>/`` trees).
                - ``paths[2]``: MaterialsFramework ``sqsdb`` base directory.

            composition_list (list[list[str]]): Chemical systems to process,
                each given as a list of element symbols
                (e.g., ``[["Cr", "Hf", "Ta"]]``).
            level (int): SQS level controlling which sqsdb sub-directories
                are included in the fit.
            sqsgen_levels2 (list[dict] | None, optional): Fixed-sublattice
                species definitions.  Each dict must have at least
                ``"element"``, ``"letter"``, ``"compositions"``, and
                ``"count"`` keys.  If ``None``, auto-derived from the
                ``"coords"`` strings of each phase in *phases*.
            skip_existing (bool, optional): Skip a composition if its output
                directory already contains a ``.tdb`` file.
                Defaults to ``False``.
            output_dir (str | Path | None, optional): Root directory under
                which per-composition output folders are created.  Defaults
                to ``paths[1]`` when ``None``.
            terms_in (dict[str, str] | None, optional): Per-phase
                ``terms.in`` content written between the two ``sqs2tdb -fit``
                calls.  Keys are lattice base names (e.g., ``"PHASE1"``);
                values are the full file content as a string.  Phases absent
                from the dict use the ATAT-generated template unchanged.
                Defaults to ``None``.
            mult_in (dict[str, str] | None, optional): Per-phase ``mult.in``
                content written between the two ``sqs2tdb -fit`` calls,
                after ``terms.in``.  Keys are lattice base names; values are
                the full file content as a string.  Phases absent from the
                dict keep the ATAT-generated ``mult.in`` unchanged.
                Defaults to ``None``.
            sublattice_map (dict[str, dict[str, list[str]]] | None, optional):
                Per-phase sublattice species assignments written to
                ``species.in`` between the two ``sqs2tdb -cp`` calls.
                Outer key is lattice base name; inner key is sublattice letter;
                value is list of element symbols
                (e.g., ``{"PHASE1": {"a": ["Cr", "Hf"], "b": ["Al"]}}``).
                Phases absent from the dict use the ATAT default (all species
                on all variable sublattices).  Defaults to ``None``.
            tdb_params (dict | None, optional): Calculator and fit parameters
                forwarded to :class:`~materialsframework.tools.sqs2tdb.Sqs2tdb`.
                Supported keys:

                - ``"fmax"`` (float): Force convergence criterion in eV/Å.
                  Default ``0.005``.
                - ``"verbose"`` (bool): Print relaxation output. Default ``True``.
                - ``"calculator"`` (str): Device string passed to
                  ``GraceCalculator`` (e.g. ``"cuda"`` or ``"cpu"``).
                  Default ``"cuda"``.
                - ``"t_min"`` (float): Lower temperature bound for CALPHAD
                  fit in K. Default ``298.15``.
                - ``"t_max"`` (float): Upper temperature bound for CALPHAD
                  fit in K. Default ``10000.0``.
                - ``"sro"`` (bool): Include short-range order correction.
                  Default ``False``.
                - ``"bv"`` (float): Energy bump value. Default ``1e-3``.
                - ``"phonon"`` (bool): Include phonon end-member contributions.
                  Default ``False``.
                - ``"open_calphad"`` (bool): Write OpenCalphad-compliant
                  ``.tdb`` output. Default ``False``.
                - ``"terms"`` (str | None): Additional interaction terms
                  string. Default ``None``.

                Unset keys fall back to their listed defaults. Defaults to
                ``None`` (all defaults active).
        """
        self.phases = phases
        self.liquid = liquid
        self.path0 = Path(paths[0])
        self.path1 = Path(paths[1])
        self.path2 = Path(paths[2])
        self.composition_list = composition_list
        self.level = level
        if sqsgen_levels2 is None:
            if phases_dict is not None:
                self.sqsgen_levels2 = {
                    p["lattice"]: self._auto_sqsgen_levels2([phases_dict[p["lattice"]]])
                    for p in phases
                    if p.get("lattice") in phases_dict
                }
            else:
                self.sqsgen_levels2 = self._auto_sqsgen_levels2(phases)
        else:
            self.sqsgen_levels2 = sqsgen_levels2
        self.skip_existing = skip_existing
        self.output_dir = Path(output_dir) if output_dir is not None else self.path1
        self.terms_in = terms_in
        self.mult_in = mult_in
        self.sublattice_map = sublattice_map
        self.tdb_params = tdb_params
        self.fixed_compositions = fixed_compositions or {}
        self.composition_elements: dict[tuple, list[str]] = {}  # (level, frac_str) → elements

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def fit(self) -> None:
        """Run the SQS-to-TDB fitting loop over all compositions.

        For each unique system size, copies the corresponding SQS directories
        from the BLADE staging tree (``path1/Files/SQS/<lattice>_<n>/``) into
        the MaterialsFramework sqsdb tree (``path2/<lattice>_<n>/``).  It then
        iterates over every composition in :attr:`composition_list`, creates a
        per-composition output folder, and calls
        :meth:`~materialsframework.tools.sqs2tdb.Sqs2tdb.fit`.

        The current working directory is temporarily changed to the
        per-composition folder and restored to :attr:`path0` afterward, as
        required by the MaterialsFramework fitting workflow.
        """
        for comp in self.composition_list:
            length = len(comp)
            comp_name = "".join(comp)
            comp_dir = self.output_dir / comp_name

            if self.skip_existing:
                comp_dir.mkdir(parents=True, exist_ok=True)
                continue
            elif comp_dir.exists():
                shutil.rmtree(comp_dir)

            for phase in self.phases:
                lattice = phase["lattice"]
                src = self.path1 / "Files" / "SQS" / f"{lattice}_{length}"
                dst = self.path2 / f"{lattice}_{length}"
                if src.exists():
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
                else:
                    print(f"Warning: SQS source not found, skipping copy: {src}")

            comp_dir.mkdir(parents=True, exist_ok=True)
            phase_list = self._build_phase_list(length)
            os.chdir(comp_dir)
            _has_constant = bool(self.composition_elements) or (
                self.sublattice_map and any(
                    bool(phase_map.get("Constant"))  # only True if Constant list is non-empty
                    for phase_map in self.sublattice_map.values()
                )
            )
            try:
                self._run_sqsfit(comp, phase_list, _has_constant)
            finally:
                os.chdir(self.path0)
                tdb_exists = any(comp_dir.glob("*.tdb"))
                if self.fixed_compositions and _has_constant and tdb_exists:
                    self._filter_permutation_dirs(comp_dir, phase_list)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _watch_and_delete(self, comp_dir: Path, phase_list: list[str],
                          stop_event: threading.Event) -> None:
        """Background thread: delete non-matching sqs_lev dirs as soon as they appear."""
        import re as _re

        desired: dict[str, float] = {}
        for phase_map in (self.sublattice_map or {}).values():
            for letter, elements in phase_map.items():
                if letter == "Constant" or not isinstance(elements, list):
                    continue
                if letter not in self.fixed_compositions:
                    continue
                for el, frac in zip(elements, self.fixed_compositions[letter]):
                    desired[f"{letter}_{el}"] = round(float(frac), 8)

        if not desired:
            return

        print(f"[watcher] desired={desired}")

        def _matches(name: str) -> bool:
            out = {}
            for m in _re.finditer(r'([a-z]_[A-Z][a-z]?)=([\d.]+)', name):
                out[m.group(1)] = round(float(m.group(2)), 8)
            result = all(abs(out.get(k, -1) - v) < 1e-6 for k, v in desired.items())
            if not result:
                misses = {k: (out.get(k, "MISSING"), v) for k, v in desired.items()
                          if abs(out.get(k, -1) - v) >= 1e-6}
                print(f"[watcher] MISMATCH {name[:60]}... misses={misses}")
            return result

        seen: set[str] = set()
        while not stop_event.is_set():
            for phase in phase_list:
                phase_dir = comp_dir / phase
                if not phase_dir.exists():
                    continue
                for d in list(phase_dir.iterdir()):
                    if not d.is_dir() or d.name in seen:
                        continue
                    # Wait for element-specific rename (c_C= suffix signals completion)
                    if "c_C=" not in d.name:
                        continue
                    seen.add(d.name)
                    if not _matches(d.name):
                        shutil.rmtree(d, ignore_errors=True)
                        print(f"[watcher] Deleted {d.name}")
                    else:
                        print(f"[watcher] Keeping {d.name}")
            time.sleep(0.01)

    def _filter_permutation_dirs(self, comp_dir: Path, phase_list: list[str]) -> None:
        """Delete sqs_lev dirs whose element-fraction assignments don't match
        the desired ordering from sublattice_map + fixed_compositions.

        Parses each dir name to extract element→fraction pairs and compares
        against the desired mapping. Keeps only the exact match.
        """
        import re as _re

        # Build desired element→fraction dict from sublattice_map + fixed_compositions
        desired: dict[str, float] = {}
        for phase_map in (self.sublattice_map or {}).values():
            for letter, elements in phase_map.items():
                if letter == "Constant":
                    continue
                if not isinstance(elements, list) or letter not in self.fixed_compositions:
                    continue
                for el, frac in zip(elements, self.fixed_compositions[letter]):
                    desired[f"{letter}_{el}"] = round(float(frac), 8)

        if not desired:
            return

        print(f"[filter] desired mapping: {desired}")

        def _parse_dir_fracs(name: str) -> dict[str, float]:
            out = {}
            for m in _re.finditer(r'([a-z]_[A-Z][a-z]?)=([\d.]+)', name):
                out[m.group(1)] = round(float(m.group(2)), 8)
            return out

        def _matches(name: str) -> bool:
            parsed = _parse_dir_fracs(name)
            return all(
                abs(parsed.get(k, -1) - v) < 1e-6
                for k, v in desired.items()
            )

        for phase in phase_list:
            phase_dir = comp_dir / phase
            print(f"[filter] scanning {phase_dir} (exists={phase_dir.exists()})")
            if not phase_dir.exists():
                continue
            kept = deleted = 0
            for d in list(phase_dir.iterdir()):
                if not d.is_dir():
                    continue
                m = _matches(d.name)
                print(f"[filter]   {'KEEP' if m else 'DEL '} {d.name}")
                if m:
                    kept += 1
                else:
                    shutil.rmtree(d)
                    deleted += 1
            print(f"[filter] done: kept={kept} deleted={deleted}")

    @staticmethod
    def _auto_sqsgen_levels2(phases: list[dict]) -> list[dict]:
        """Derive fixed-sublattice definitions from phase coords.

        Parses each phase's ``coords`` string and identifies sites whose
        label is an element symbol (starts with an uppercase letter) rather
        than a lowercase sublattice letter.  Returns one entry per unique
        fixed element, with the next available lowercase letter that is not
        already used as a variable-sublattice label.

        Args:
            phases (list[dict]): Phase prototype dicts, each containing a
                ``"coords"`` key with fractional-coordinate lines of the form
                ``"x y z LABEL"``.

        Returns:
            list[dict]: Entries for fixed sublattices, each with keys
            ``"element"``, ``"letter"``, ``"compositions"``, and ``"count"``.
            Returns ``[]`` if no fixed sites are found.
        """
        used_var_letters: set[str] = set()
        fixed_counts: dict[str, int] = {}

        for phase in phases:
            coords = phase.get("coords", "")
            for line in coords.strip().splitlines():
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                label = parts[3]
                if label.islower() and len(label) == 1:
                    used_var_letters.add(label)
                elif label and label[0].isupper():
                    fixed_counts[label] = fixed_counts.get(label, 0) + 1

        if not fixed_counts:
            return []

        result: list[dict] = []
        alphabet = iter("abcdefghijklmnopqrstuvwxyz")
        for element in sorted(fixed_counts):
            while True:
                letter = next(alphabet)
                if letter not in used_var_letters:
                    break
            result.append({
                "element": element,
                "letter": letter,
                "compositions": "1.0",
                "count": str(fixed_counts[element]),
            })
        return result

    def _build_phase_list(self, length: int) -> list[str]:
        """Construct the phase identifier list for a given system size.

        Appends ``_<length>`` to each base phase name and optionally adds
        ``"LIQUID"`` when :attr:`liquid` is ``True``.

        Args:
            length (int): Number of elements in the target chemical system.

        Returns:
            list[str]: Phase identifiers ready for
            :meth:`~materialsframework.tools.sqs2tdb.Sqs2tdb.fit`
            (e.g., ``["FCC_3", "BCC_3", "LIQUID"]`` or ``["PHASE1_2", "LIQUID"]``).
        """
        phase_list = [f"{phase['lattice']}_{length}" for phase in self.phases]
        if self.liquid:
            phase_list.append("LIQUID")
        return phase_list

    def _refit_tdb(self, comp: list[str], phase_list: list[str]) -> None:
        """Refit only: update terms.in if changed, then rerun sqs2tdb -fit and -tdb.

        Used when skip_existing=True to regenerate the TDB without re-relaxing structures.
        """
        import subprocess as _sp

        for lattice in phase_list:
            lat_dir = Path(lattice)
            if not lat_dir.exists():
                continue
            # Write updated terms.in if provided
            lattice_base = lattice.rsplit("_", 1)[0]
            if self.terms_in and lattice_base in self.terms_in:
                (lat_dir / "terms.in").write_text(self.terms_in[lattice_base])
                print(f"Updated terms.in for {lattice}")
            # Rerun sqs2tdb -fit
            r = _sp.run(["sqs2tdb", "-fit"], cwd=lat_dir, capture_output=True, text=True)
            if r.stdout.strip():
                print(r.stdout.strip()[:200])
        # Rerun sqs2tdb -tdb at comp level
        r = _sp.run(["sqs2tdb", "-tdb"], capture_output=True, text=True)
        if r.stdout.strip():
            print(r.stdout.strip()[:200])

    def _run_sqsfit(self, comp: list[str], phase_list: list[str],
                    use_filter: bool = False) -> None:
        """Invoke Sqs2tdb for a single composition.

        Constructs a :class:`~materialsframework.calculators.GraceCalculator`
        and a :class:`~materialsframework.tools.sqs2tdb.Sqs2tdb` instance,
        then calls :meth:`~materialsframework.tools.sqs2tdb.Sqs2tdb.fit`.

        Args:
            comp (list[str]): Element symbols for the target chemical system.
            phase_list (list[str]): Phase identifiers to include in the fit.
        """
        from materialsframework.calculators.registry import get_calculator
        from materialsframework.tools.sqs2tdb import Sqs2tdb

        p = self.tdb_params or {}
        calc_name = p.get("calculator", "grace")
        calc_kwargs = p.get("calculator_kwargs", {"steps": 1000, "device": "cuda"})

        # backward compat: if calculator was set to a device string, assume grace
        if calc_name in ("cpu", "cuda", "mps"):
            calc_kwargs = dict(calc_kwargs)
            calc_kwargs.setdefault("device", calc_name)
            calc_name = "grace"

        calc = get_calculator(calc_name, **calc_kwargs)
        sqs = Sqs2tdb(
            fmax=p.get("fmax", 0.005),
            verbose=p.get("verbose", True),
            track_trajectory=p.get("track_trajectory", True),
            calculator=calc,
        )
        # Strip BLADE-internal keys ("Constant") before passing to MaterialsFramework
        clean_sublattice_map = None
        if self.sublattice_map:
            clean_sublattice_map = {
                phase: {k: v for k, v in phase_map.items() if k != "Constant"}
                for phase, phase_map in self.sublattice_map.items()
            }

        if use_filter:
            import re as _re

            if self.composition_elements:
                # Per-composition element spec: (level, frac_str) → allowed elements list.
                # dir_filter keeps a dir only if its elements match the spec for that level/frac.
                # Dirs whose (level, frac) has no spec are kept (no restriction).
                comp_el = self.composition_elements

                def _dir_filter(subdir: Path) -> bool:
                    name = subdir.name
                    level_m = _re.match(r'sqs_lev=(\d+)_', name)
                    if not level_m:
                        return True
                    level_n = int(level_m.group(1))
                    by_letter: dict[str, list[tuple[str, float]]] = {}
                    for m in _re.finditer(r'([a-z])_([A-Z][a-z]?)=([\d.]+)', name):
                        letter, el, frac = m.group(1), m.group(2), float(m.group(3))
                        by_letter.setdefault(letter, []).append((el, frac))
                    a_pairs = sorted(by_letter.get('a', []), key=lambda x: -x[1])
                    if not a_pairs:
                        return True
                    a_fracs = [str(round(f, 6)) for _, f in a_pairs]
                    while len(a_fracs) > 1 and a_fracs[-1] == "0.0":
                        a_fracs.pop()
                    frac_key = ",".join(a_fracs)
                    key = (level_n, frac_key)
                    if key in comp_el:
                        actual = sorted(el for el, _ in a_pairs)
                        expected = sorted(comp_el[key])
                        if len(actual) < len(expected):
                            # More elements specified than fracs — allow any subset (permutations)
                            return set(actual).issubset(set(expected))
                        return actual == expected
                    return True

                sqs.dir_filter = _dir_filter

            elif self.fixed_compositions:
                desired: dict[str, float] = {}
                for phase_map in (self.sublattice_map or {}).values():
                    for letter, elements in phase_map.items():
                        if letter == "Constant" or not isinstance(elements, list):
                            continue
                        if letter not in self.fixed_compositions:
                            continue
                        for el, frac in zip(elements, self.fixed_compositions[letter]):
                            desired[f"{letter}_{el}"] = round(float(frac), 8)

                def _dir_filter(subdir: Path) -> bool:  # type: ignore[no-redef]
                    name = subdir.name
                    parsed = {}
                    for m in _re.finditer(r'([a-z]_[A-Z][a-z]?)=([\d.]+)', name):
                        parsed[m.group(1)] = round(float(m.group(2)), 8)
                    return all(abs(parsed.get(k, -1) - v) < 1e-6 for k, v in desired.items())

                sqs.dir_filter = _dir_filter

        smap_els = [
            el for phase_map in (self.sublattice_map or {}).values()
            for letter, elems in phase_map.items()
            if letter != "Constant" and isinstance(elems, list)
            for el in elems
            if el not in comp
        ]
        full_species = list(comp) + list(dict.fromkeys(smap_els))

        sqs.fit(
            paths=self.path2,
            sqsgen_levels2=self.sqsgen_levels2,
            species=full_species,
            lattices=phase_list,
            level=self.level,
            t_min=p.get("t_min", 298.15),
            t_max=p.get("t_max", 10000.0),
            sro=p.get("sro", False),
            bv=p.get("bv", 5e-3),
            phonon=p.get("phonon", False),
            open_calphad=p.get("open_calphad", False),
            terms=p.get("terms", None),
            terms_in=self.terms_in,
            mult_in=self.mult_in,
            sublattice_map=clean_sublattice_map,
            skip_existing=self.skip_existing,
        )
