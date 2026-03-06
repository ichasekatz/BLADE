"""

"""
import math
import os
from pathlib import Path

import numpy as np
import re

def build_supercell(frac, rep):
    nx,ny,nz = rep
    shifts = np.array([[i,j,k] for i in range(nx)
                            for j in range(ny)
                            for k in range(nz)],dtype=float)
    return (frac[None,:,:] + shifts[:,None,:]).reshape(-1,3)


def min_image(df,rep):
    rep=np.array(rep,dtype=float)
    return df - rep*np.round(df/rep)

class BladeCutoff:
    def __init__(self):
        self.self = self

    def lattice_from_params(self, a,b,c,alpha,beta,gamma):
        """Convert a,b,c,alpha,beta,gamma -> 3x3 lattice (Å)."""
        alpha=np.radians(alpha)
        beta=np.radians(beta)
        gamma=np.radians(gamma)

        ax = a
        ay = 0
        az = 0

        bx = b*np.cos(gamma)
        by = b*np.sin(gamma)
        bz = 0

        cx = c*np.cos(beta)
        cy = c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cz = np.sqrt(c**2 - cx**2 - cy**2)

        return np.array([
            [ax,ay,az],
            [bx,by,bz],
            [cx,cy,cz]
        ])


    def read_coords(self, coord_string):
        lines = [l.strip() for l in coord_string.strip().split("\n")]
        frac=[]
        for l in lines:
            vals = re.split(r"\s+", l)
            frac.append([float(vals[0]), float(vals[1]), float(vals[2])])
        return np.array(frac)


    def get_shells(self, lattice,frac,rep,tol=1e-4):
        frac_sc = build_supercell(frac,rep)
        n=len(frac_sc)
        dists=[]

        for i in range(n):
            df = frac_sc[i+1:] - frac_sc[i]
            df = min_image(df,rep)
            cart = df @ lattice
            r = np.linalg.norm(cart,axis=1)
            for val in r:
                if val>1e-6:
                    dists.append(val)

        dists=np.sort(dists)

        shells=[]
        for r in dists:
            if len(shells)==0 or abs(r-shells[-1])>tol:
                shells.append(r)

        shells = np.array(shells)
        shells = shells / shells[0]

        return shells
    
    # def derive_cutoffs0(self, num_elements, shell_distances):
    #     shells = np.array(shell_distances, dtype=float)
    #     shells = shells / shells[0]

    #     if num_elements == 2:
    #         pair_shells = 3
    #     elif num_elements == 3:
    #         pair_shells = 2
    #     else:
    #         pair_shells = 2

    #     pair_shells = min(pair_shells, len(shells)-1)

    #     r2 = (shells[pair_shells-1] + shells[pair_shells]) / 2
    #     r3 = (shells[0] + shells[1]) / 2
    #     r4 = shells[0] * 1.01

    #     return {
    #         "-2": round(float(r2), 1),
    #         "-3": round(float(r3), 1),
    #         "-4": round(float(r4), 1),
    #     }

    def derive_cutoffs(self, num_elements, shells):
            shell_distances = np.asarray(shells)
            idx = np.where(np.abs(shell_distances - 2) < 1e-3)[0]
            idx = idx[0]
            print(idx)
            if idx == 3:
                if num_elements == 2:
                    cut2 = shell_distances[idx]
                    cut3 = math.ceil(((shell_distances[idx-1])+(shell_distances[idx-2]))/2*10)/10
                    cut4 = math.floor(shell_distances[idx-3]*10)/10
                else:
                    cut2 = math.floor(shell_distances[idx-1]*10)/10
                    cut3 = math.floor(shell_distances[idx-3]*10)/10
                    cut4 = math.floor(shell_distances[idx-3]*10)/10
            elif idx == 4:
                if num_elements == 2:
                    cut2 = shell_distances[idx]
                    cut3 = math.ceil(((shell_distances[idx-2])+(shell_distances[idx-3]))/2*10)/10
                    cut4 = math.ceil(((shell_distances[idx-2])+(shell_distances[idx-3]))/2*10)/10
                else:
                    cut2 = math.ceil(((shell_distances[idx-1])+(shell_distances[idx-2]))/2*10)/10
                    cut3 = math.ceil(((shell_distances[idx-2])+(shell_distances[idx-3]))/2*10)/10
                    cut4 = math.ceil(((shell_distances[idx-3])+(shell_distances[idx-4]))/2*10)/10
            elif idx == 5:
                if num_elements == 2:
                    cut2 = shell_distances[idx]
                    cut3 = math.ceil(((shell_distances[idx-2])+(shell_distances[idx-3]))/2*10)/10
                    cut4 = math.ceil(((shell_distances[idx-3])+(shell_distances[idx-4]))/2*10)/10
                else:
                    cut2 = math.ceil(((shell_distances[idx-1])+(shell_distances[idx-2]))/2*10)/10
                    cut3 = math.ceil(((shell_distances[idx-3])+(shell_distances[idx-4]))/2*10)/10
                    cut4 = math.ceil(((shell_distances[idx-4])+(shell_distances[idx-5]))/2*10)/10
            return {
                "-2": cut2,
                "-3": cut3,
                "-4": cut4
            }



    # def derive_cutoffs1(self, num_elements, n_sites, nn_list, max_pair_shells=6):
    #     shells = np.array(nn_list, dtype=float)
    #     shells = np.sort(shells)

    #     if len(shells) < 2:
    #         raise ValueError("Need at least 1NN and 2NN values.")

    #     shells = shells / shells[0]

    #     q = int(num_elements)
    #     if q < 2:
    #         raise ValueError("num_elements must be >= 2.")

    #     n_sites = int(n_sites)
    #     if n_sites < 2:
    #         raise ValueError("n_sites must be >= 2.")

    #     def midpoint_cutoff(k):
    #         k = min(k, len(shells) - 1)
    #         return 0.5 * (shells[k - 1] + shells[k])

    #     r3 = midpoint_cutoff(1)          # between 1NN and 2NN
    #     r4 = 1.01                        # strict 1NN

    #     # estimate "difficulty": independent correlations grow like (q-1)^m
    #     # pairs: (q-1)^2, triplets: (q-1)^3, quadruplets: (q-1)^4
    #     # we keep -3/-4 local, so only pairs drive the shell depth decision.
    #     pair_state_dim = (q - 1) ** 2

    #     # heuristic capacity: you can match O(n_sites) constraints; pairs add ~pair_state_dim per shell
    #     # so number of pair shells you can afford is roughly proportional to n_sites / pair_state_dim
    #     # clamp to [1, max_pair_shells]
    #     k2 = int(np.floor(n_sites / max(1, 8 * pair_state_dim))) + 1
    #     k2 = max(1, min(k2, max_pair_shells))

    #     # hard caps by component count (keeps SQS solvable in high-q systems)
    #     if q == 2:
    #         k2 = max(k2, 3)              # binaries usually can afford 3 shells
    #     elif q == 3:
    #         k2 = min(k2, 2)              # ternaries: keep pairs local
    #     else:
    #         k2 = min(k2, 2)              # 4+ elements: almost always 1–2 shells

    #     # if we don't have enough shells provided, reduce k2
    #     k2 = min(k2, len(shells) - 1)

    #     r2 = midpoint_cutoff(k2)

    #     return {
    #         "-2": round(float(r2), 2),
    #         "-3": round(float(r3), 2),
    #         "-4": round(float(r4), 2),
    #         "pair_shells_used": k2,
    #     }
    
    # def derive_cutoffs2(self, pair, triplet, sqsdir):
    #     PAIR_SHELLS = pair
    #     TRIPLET_SHELLS = triplet
    #     EPS = 1e-6

    #     sqsdir_path = Path(sqsdir)

    #     # find clusters file
    #     candidates = [sqsdir_path / "clusters.out"]
    #     clusters_path = None
    #     for c in candidates:
    #         if c.exists():
    #             clusters_path = c
    #             break
    #     if clusters_path is None:
    #         clus_like = sorted(sqsdir_path.glob("*clus*"), key=lambda p: p.stat().st_mtime, reverse=True)
    #         if clus_like:
    #             clusters_path = clus_like[0]
    #         else:
    #             raise FileNotFoundError(f"Could not find clusters file in {sqsdir}")

    #     # read non-empty lines
    #     lines = [ln.strip() for ln in clusters_path.read_text(errors="ignore").splitlines() if ln.strip()]

    #     pair_diams = []
    #     trip_diams = []

    #     i = 0
    #     n = len(lines)
    #     while i < n:
    #         # block header
    #         try:
    #             multiplicity = int(float(lines[i]))  # sometimes "6" sometimes "6.0000"
    #             diameter = float(lines[i + 1])
    #             npoints = int(float(lines[i + 2]))
    #         except (ValueError, IndexError) as e:
    #             raise ValueError(f"Failed parsing cluster header near line {i+1} in {clusters_path.name}: {e}")

    #         # collect by npoints
    #         if npoints == 2:
    #             pair_diams.append(diameter)
    #         elif npoints == 3:
    #             trip_diams.append(diameter)

    #         # skip the point lines
    #         i += 3 + npoints

    #     pair_shells = sorted(set(pair_diams))
    #     trip_shells = sorted(set(trip_diams))

    #     if len(pair_shells) < PAIR_SHELLS:
    #         raise ValueError(f"Not enough pair shells: found {len(pair_shells)}, need {PAIR_SHELLS}")
    #     if len(trip_shells) < TRIPLET_SHELLS:
    #         raise ValueError(f"Not enough triplet shells: found {len(trip_shells)}, need {TRIPLET_SHELLS}")

    #     cut2 = pair_shells[PAIR_SHELLS - 1] + EPS
    #     cut3 = trip_shells[TRIPLET_SHELLS - 1] + EPS

    #     print(
    #         f"Derived cutoffs from {clusters_path.name}: "
    #         f"-2={cut2:.6f}, -3={cut3:.6f} "
    #         f"(pairshell={PAIR_SHELLS}, tripshell={TRIPLET_SHELLS})"
    #     )

    #     return cut2, cut3



    def derive_cutoffs2(self, sqsdir):
        path = os.path.join(sqsdir, "clusters.out")
        data = []

        with open(path) as f:
            lines = [l.strip() for l in f if l.strip()]

        i = 0
        n = len(lines)

        while i < n:
            try:
                int(lines[i].split()[0])
            except:
                i += 1
                continue

            if i+2 >= n:
                break

            try:
                r = float(lines[i+1].split()[0])
                order = int(lines[i+2].split()[0])
            except:
                i += 1
                continue

            if order == 2:
                kind = "pair"
            elif order == 3:
                kind = "triplet"
            elif order == 4:
                kind = "quadruplet"
            else:
                kind = f"order_{order}"

            data.append((kind, r))
            i += 3 + order

        return data



