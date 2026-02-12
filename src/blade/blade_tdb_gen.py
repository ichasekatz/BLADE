"""
This module defines the `BladeTDBGen` class, which includes tools for running SQS-to-TDB fitting over
a list of chemical systems.

The `BladeTDBGen` class automates the execution of MaterialsFramework's `Sqs2tdb` fitting workflow across
multiple chemical compositions. For each composition, the class creates a dedicated output directory,
constructs the appropriate phase list, and runs the SQS-to-TDB fitting pipeline at a specified SQS level.
This functionality is intended for high-throughput thermodynamic database generation workflows in
multicomponent materials design.
"""

import os
from pathlib import Path

from materialsframework.calculators import ORBCalculator as Calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb


class BladeTDBGen:
    """
    Run Sqs2tdb fitting over many compositions.

    For each composition, this class:
      - creates an output folder under path2/<composition>
      - constructs a phase list based on the system size (e.g., ["HEDB1_3", "HEDB2_3"])
      - optionally includes the LIQUID phase
      - runs Sqs2tdb.fit() using a GraceCalculator

    Notes:
        This workflow may write files into the current working directory. The working directory is
        reset after each composition to avoid persistent side effects.
    """
    def __init__(self, phases, liquid, paths, composition_list, level):
        """
        Initializes the BladeTDBGen object.

        Args:
            phases (list[str]): Base phase names to include in the fitting workflow (e.g., ["HEDB1", "HEDB2"]).
                Each phase name is expanded with a system-size suffix during execution.
            liquid (bool): If True, include the LIQUID phase in the fitting.
            paths (list[str | Path]): Path bundle used by the workflow. Expected indices:
                  - paths[0]: root or original working directory
                  - paths[2]: base output directory for per-composition folders
            composition_list (list[list[str]]): List of chemical systems to process, where each system
                is represented by a list of element symbols. Example: [["Cr", "Hf"], ["Cr", "Hf", "Ta"]]
            level (int): SQS level passed through to Sqs2tdb.fit().
        """
        self.phases = phases
        self.liquid = liquid
        self.path0 = paths[0]
        self.path2 = paths[2]
        self.composition_list = composition_list
        self.level = level

        def sqsfit_func(comp, phases, level):
            """
            Run Sqs2tdb fitting for a single composition.

            This helper function constructs a Sqs2tdb object with a GraceCalculator backend and executes the
            fitting routine for the specified composition and phase list.

            Args:
                comp (list[str]): Element symbols defining the target chemical system.
                phases (list[str]): Phase or lattice identifiers to include in the fit.
                level (int): SQS level controlling the database depth and fitting configuration.
            """
            calc = Calculator(steps=1000, device="cpu")

            sqs = Sqs2tdb(fmax=0.001, verbose=True, calculator=calc)

            sqs.fit(
                species=comp,
                lattices=phases,
                level=level,
            )

        def directory(phases, liquid, length):
            """
            Construct a list of phase names for a given system size.

            Each base phase name is expanded by appending an underscore and the number of elements in the
            chemical system. Optionally, the LIQUID phase is appended.

            Args:
                phases (list[str]): Base phase names (e.g., ["HEDB1", "HEDB2"]).
                liquid (bool): Whether to include the LIQUID phase.
                length (int): Number of elements in the chemical system.

            Returns:
                list[str]: Expanded list of phase identifiers suitable for Sqs2tdb fitting.
            """
            phase_list = []
            for phase in phases:
                phase_list += [f"{phase}_{length}"]
            if liquid:
                phase_list.append("LIQUID")
            return phase_list

        """
        Loop over all compositions and run SQS-to-TDB fitting.

        For each composition, an output directory is created, the corresponding phase list is
        generated based on system size, and the Sqs2tdb fitting workflow is executed.
        """
        for comp in composition_list:
            length = len(comp)
            file = f"{self.path2}{''.join(comp)}"
            Path(file).mkdir(parents=True, exist_ok=True)
            phase_list = directory(self.phases, self.liquid, length)
            os.chdir(file)
            sqsfit_func(comp, phase_list, self.level)
            os.chdir(self.path0)
