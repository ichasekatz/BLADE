import os
from pathlib import Path
from materialsframework.calculators import GraceCalculator as Calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb
from blade.blade_compositions import BladeCompositions
from blade.blade_sqs import BladeSQS

class BladeTDBGen:
    def __init__(self, phases, liquid, paths, composition_list, level):
        self.phases = phases
        self.liquid = liquid
        self.path0 = paths[0]
        self.path2 = paths[2]
        self.composition_list = composition_list
        self.level = level

        def sqsfit_func(comp, phases, level):
            # Initialize the wrapper
            sqs = Sqs2tdb(fmax=0.001, verbose=True, calculator=Calculator(device="cpu"))

            # Fit SQS
            sqs.fit(
                species=comp,
                lattices=phases,
                level=level,
            )

        def directory(phases, liquid, length):
            phase_list = []
            for phase in phases:
                phase_list += [f"{phase}_{length}"]
            if liquid:
                phase_list.append("LIQUID")
            return phase_list

        for comp in composition_list:
            length = len(comp)
            # Create and move to output folder
            file = f"{self.path2}{''.join(comp)}"
            Path(file).mkdir(parents=True, exist_ok=True)

            phase_list = directory(self.phases, self.liquid, length)
            os.chdir(file)
            sqsfit_func(comp, phase_list, self.level)

            # Change back to main folder
            os.chdir(self.path0)
