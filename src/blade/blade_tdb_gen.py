import os
from pathlib import Path
from materialsframework.calculators import GraceCalculator as Calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb
from blade.blade_compositions import BladeCompositions
from blade.blade_sqs import BladeSQS

class BladeTDBGen:
    def __init__(self, phases, liquid, paths, composition_settings, sqs_gen_settings,):
        self.phases = phases
        self.liquid = liquid
        self.path0 = paths[0]
        self.path1 = paths[1]
        self.path2 = paths[2]

        self.transition_metals = composition_settings[0]
        self.rare_earths = composition_settings[1]
        self.system_size = composition_settings[2]
        self.tm_element_range = composition_settings[3]
        self.re_element_range = composition_settings[4]
        self.allow_lower_order = composition_settings[5]

        self.a = sqs_gen_settings[0]
        self.b = sqs_gen_settings[1]
        self.c = sqs_gen_settings[2]
        self.alpha = sqs_gen_settings[3]
        self.beta = sqs_gen_settings[4]
        self.gamma = sqs_gen_settings[5]
        self.rndstr = sqs_gen_settings[6]
        self.sqsgen_levels = sqs_gen_settings[7]
        self.level = sqs_gen_settings[8]
        self.time = sqs_gen_settings[9]
        
        compositions = BladeCompositions(
            self.transition_metals,
            self.rare_earths,
            self.system_size,
            tm_min=self.tm_element_range[0],
            tm_max=self.tm_element_range[1],
            re_min=self.re_element_range[0],
            re_max=self.re_element_range[1],
            allow_lower_order=self.allow_lower_order,
        )

        composition_list = compositions.generate_compositions()

        unique_len_comps = compositions.get_systems()

        print("Compositions: ", composition_list)
        print("Total # compositions: ", len(composition_list))
        print("Unique length compositions: ", unique_len_comps)

        os.chdir(self.path1)

        sqs_gen = BladeSQS(self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.rndstr, self.sqsgen_levels, self.level)

        for phase in self.phases:
            sqs_gen.sqs_gen(unique_len_comps, phase, self.path1, self.time)

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
