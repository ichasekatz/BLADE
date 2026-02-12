"""

"""
import numpy as np
import math
import itertools


class BladeCutoff:
    def __init__(self, phases_dict, tol=1e-3, supercell=2):

        self.a = phases_dict["a"]
        self.b = phases_dict["b"]
        self.c = phases_dict["c"]
        self.alpha = phases_dict["alpha"]
        self.beta = phases_dict["beta"]
        self.gamma = phases_dict["gamma"]

        self.unit_cell = phases_dict["coords"]

        self.tol = tol
        self.supercell = supercell

        self.frac, self.species = self._parse_coords()

        self.A = self._lattice_vectors()

    def _parse_coords(self):
        frac = []
        species = []

        for line in self.unit_cell.strip().splitlines():
            toks = line.split()
            if len(toks) < 4:
                continue

            frac.append([float(toks[0]), float(toks[1]), float(toks[2])])
            species.append(toks[3])

        return np.array(frac, dtype=float), species

    def _lattice_vectors(self):

        alpha = math.radians(self.alpha)
        beta = math.radians(self.beta)
        gamma = math.radians(self.gamma)

        a1 = np.array([self.a, 0.0, 0.0])
        a2 = np.array([self.b * math.cos(gamma), self.b * math.sin(gamma), 0.0])

        cx = self.c * math.cos(beta)
        cy = self.c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
        cz = math.sqrt(max(self.c**2 - cx**2 - cy**2, 0.0))

        a3 = np.array([cx, cy, cz])

        return np.vstack([a1, a2, a3]).T

    def compute_shells(self):

        cart = self.frac @ self.A.T
        dists = []

        R = self.supercell
        translations = list(itertools.product(range(-R, R+1), repeat=3))

        for i in range(len(cart)):
            ri = cart[i]
            for j in range(len(self.frac)):
                fj = self.frac[j]
                for t in translations:
                    if i == j and t == (0,0,0):
                        continue

                    rj = (fj + np.array(t)) @ self.A.T
                    d = np.linalg.norm(rj - ri)

                    if d > 1e-8:
                        dists.append(d)

        dists = np.array(dists)
        dists.sort()

        shells = [dists[0]]

        for d in dists[1:]:
            if abs(d - shells[-1]) > self.tol:
                shells.append(d)

        return shells

    def get_cutoffs(self, pair_shells=3):

        shells = self.compute_shells()

        if len(shells) < pair_shells + 2:
            raise RuntimeError(
                "Not enough shells detected — increase supercell to 3."
            )

        cut2 = 0.5 * (shells[pair_shells-1] + shells[pair_shells])
        cut3 = 0.5 * (shells[pair_shells] + shells[pair_shells+1])

        return cut2, cut3, shells
    



