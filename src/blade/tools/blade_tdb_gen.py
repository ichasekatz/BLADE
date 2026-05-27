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
        sqsgen_levels2=sqsgen_levels2,
        skip_existing=True,
    )
    gen.fit()
"""

from __future__ import annotations

import os
import shutil
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
        skip_existing (bool): If ``True``, skip compositions that already
            have a ``.tdb`` file in their output directory.
    """

    def __init__(
        self,
        phases: list[dict],
        liquid: bool,
        paths: Sequence[str | Path],
        composition_list: list[list[str]],
        level: int,
        sqsgen_levels2: list[dict],
        skip_existing: bool = False,
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
            sqsgen_levels2 (list[dict]): Fixed-sublattice species definitions.
                Each dict must have at least ``"element"``, ``"letter"``,
                ``"compositions"``, and ``"count"`` keys.
            skip_existing (bool, optional): Skip a composition if its output
                directory already contains a ``.tdb`` file.
                Defaults to ``False``.
        """
        self.phases = phases
        self.liquid = liquid
        self.path0 = Path(paths[0])
        self.path1 = Path(paths[1])
        self.path2 = Path(paths[2])
        self.composition_list = composition_list
        self.level = level
        self.sqsgen_levels2 = sqsgen_levels2
        self.skip_existing = skip_existing

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def fit(self) -> None:
        """Run the SQS-to-TDB fitting loop over all compositions.

        For each unique system size, this method first copies the
        corresponding SQS directories from the BLADE staging tree
        (``path1/SQS/<lattice>_<n>/``) into the MaterialsFramework sqsdb
        tree (``path2/<lattice>_<n>/``).  It then iterates over every
        composition in :attr:`composition_list`, creates a per-composition
        output folder under ``path1``, and calls
        :meth:`~materialsframework.tools.sqs2tdb.Sqs2tdb.fit`.

        The current working directory is temporarily changed to the
        per-composition folder and restored to :attr:`path0` afterward, as
        required by the MaterialsFramework fitting workflow.
        """
        unique_lengths = sorted({len(comp) for comp in self.composition_list})

        for length in unique_lengths:
            for phase in self.phases:
                lattice = phase["lattice"]
                src = self.path1 / "SQS" / f"{lattice}_{length}"
                dst = self.path2 / f"{lattice}_{length}"
                if src.exists():
                    shutil.copytree(src, dst, dirs_exist_ok=True)
                    print(f"Copied {src} -> {dst}")
                else:
                    print(f"Skipping missing folder: {src}")

        for comp in self.composition_list:
            length = len(comp)
            comp_name = "".join(comp)
            comp_dir = self.path1 / comp_name

            if self.skip_existing:
                existing_tdbs = list(comp_dir.glob("*.tdb"))
                if existing_tdbs:
                    print(
                        f"Skipping {comp_name}: found existing TDB file(s): "
                        f"{[p.name for p in existing_tdbs]}"
                    )
                    continue

            comp_dir.mkdir(parents=True, exist_ok=True)
            phase_list = self._build_phase_list(length)
            os.chdir(comp_dir)
            self._run_sqsfit(comp, phase_list)
            os.chdir(self.path0)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

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

    def _run_sqsfit(self, comp: list[str], phase_list: list[str]) -> None:
        """Invoke Sqs2tdb for a single composition.

        Constructs a :class:`~materialsframework.calculators.GraceCalculator`
        and a :class:`~materialsframework.tools.sqs2tdb.Sqs2tdb` instance,
        then calls :meth:`~materialsframework.tools.sqs2tdb.Sqs2tdb.fit`.

        Args:
            comp (list[str]): Element symbols for the target chemical system.
            phase_list (list[str]): Phase identifiers to include in the fit.
        """
        from materialsframework.calculators import GraceCalculator
        from materialsframework.tools.sqs2tdb import Sqs2tdb

        calc = GraceCalculator(steps=1000, device="cuda")
        sqs = Sqs2tdb(fmax=0.005, verbose=True, calculator=calc)
        sqs.fit(
            paths=self.path2,
            sqsgen_levels2=self.sqsgen_levels2,
            species=comp,
            lattices=phase_list,
            level=self.level,
        )
