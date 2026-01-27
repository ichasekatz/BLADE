import itertools
import os


class BladePhaseDiagram:
    def __init__(self, liquid, path2, phase):

        self.liquid = liquid
        self.path2 = path2
        self.phase = phase

    def directory(self, comp):
        length = len(comp)
        # Create and move to output folder
        file = f'{self.path2}{"".join(comp)}'
        os.makedirs(file, exist_ok=True)

        # Specify inputs
        if self.liquid:
            phases = [(self.phase+"_"+str(length)), "LIQUID"]
        else:
            phases = [(self.phase+"_"+str(length))]

        return file, phases

    def file_names(self, comp):
        # Editing names
        elements = [el.upper() for el in comp]
        file_names = ["_".join(p) for p in itertools.permutations(elements)]
        print(file_names)

        return file_names, elements
