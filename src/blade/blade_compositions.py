import itertools


class BladeCompositions:
    def __init__(self, transition_metals, rare_earths, system_size, tm_min, tm_max, re_min, re_max, allow_lower_order):
        self.transition_metals = transition_metals
        self.rare_earths = rare_earths
        self.system_size = system_size
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.re_min = re_min
        self.re_max = re_max
        self.allow_lower_order = allow_lower_order

    def generate_compositions(self):
        # Generate compositions

        # Generate all possible transition metal compositions
        tm_elements = list(range(self.tm_min, self.tm_max + 1))
        re_elements = list(range(self.re_min, self.re_max + 1))

        tm_combos = []
        re_combos = []
        combined_comps = []
        compositions = []

        # Generate all possible combinations of transition metals
        for i in tm_elements:
            # If only transition metals are 0, add all rare earth combinations
            if i == 0:
                for j in re_elements:
                    # If both are 0, only add Boron
                    if j == 0:
                        compositions += [[""]]
                    else:
                        combined_comps += [list(c) for c in itertools.combinations(self.rare_earths, j)]
            if i != 0:
                tm_combos += [list(c) for c in itertools.combinations(self.transition_metals, i)]

        # Generate all possible combinations of rare earths
        for j in re_elements:
            # If rare earths are 0, add all transition metal combinations
            if j == 0:
                combined_comps += tm_combos
            else:
                re_combos += [list(c) for c in itertools.combinations(self.rare_earths, j)]

        # Combine transition metal and rare earth combinations
        for tm_comp in tm_combos:
            for re_comp in re_combos:
                combined_comps.append(list(tm_comp) + re_comp)

        # Remove compositions that exceed system size
        combined_comps = [c for c in combined_comps if len(c) <= self.system_size]

        # Add Boron to each combination, remove lower order compositions if needed, and sort
        for i in combined_comps:
            compositions += [sorted(i)]
        if not self.allow_lower_order:
            compositions = [c for c in compositions if len(c) == self.system_size]
        compositions = sorted(compositions)
        print(compositions)
        print(f"Total compositions: {len(compositions)}")

        # Ensure no duplicate compositions
        unique_comps = {tuple(sorted(c)) for c in compositions}
        print("Unique:", len(unique_comps))

        self.compositions = compositions
        return compositions

    def get_systems(self):
        compositions = self.compositions
        # Count unique system sizes
        len_comps = []
        for i in compositions:
            len_comps += [len(i)]
        print(len_comps)
        unique_len_comps = set(len_comps) if len(len_comps) >= 2 else {len_comps[0]}
        print("Unique Systems:", unique_len_comps)
        return unique_len_comps
