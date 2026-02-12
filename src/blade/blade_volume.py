"""

"""

import numpy as np
import re
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd


class BLADEVolume:
    def __init__(self):
        self.data = None
    def parse_sqs_meta(self, poscar_path: Path):
        sqs_level = None
        a_fracs = {}
        for p in poscar_path.parents:
            name = p.name
            if "sqs_lev=" in name:
                m = re.search(r"sqs_lev=(\d+)", name)
                if m:
                    sqs_level = int(m.group(1))
                for el, val in re.findall(r"a_([A-Za-z]+)=([0-9]*\.?[0-9]+)", name):
                    a_fracs[el] = float(val)
                break
        return sqs_level, a_fracs


    def _all_int(self, tokens):
        try:
            [int(t) for t in tokens]
            return True
        except Exception:
            return False


    def poscar_lattice_and_counts(self, poscar_path: Path):
        with open(poscar_path, "r") as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        scale = float(lines[1])
        lattice = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)], dtype=float) * scale

        i = 5
        toks = lines[i].split()
        if self._all_int(toks):
            elems = []
            counts = [int(x) for x in toks]
            i += 1
        else:
            elems = toks
            i += 1
            counts = [int(x) for x in lines[i].split()]
            i += 1

        if i < len(lines) and lines[i].lower().startswith("s"):
            i += 1

        natoms = int(sum(counts))
        counts_map = {}
        if elems:
            for e, c in zip(elems, counts):
                counts_map[e] = int(c)

        return lattice, natoms, counts_map


    def cellpar_from_lattice(self, lattice):
        a_vec, b_vec, c_vec = lattice
        a = float(np.linalg.norm(a_vec))
        b = float(np.linalg.norm(b_vec))
        c = float(np.linalg.norm(c_vec))

        def ang(u, v):
            cuv = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            cuv = max(-1.0, min(1.0, cuv))
            return float(math.degrees(math.acos(cuv)))

        alpha = ang(b_vec, c_vec)
        beta = ang(a_vec, c_vec)
        gamma = ang(a_vec, b_vec)
        return a, b, c, alpha, beta, gamma


    def scan_poscars(self, comp_dir: Path):
        rows = []
        comp_name = comp_dir.name
        print(f"Checking for POSCARs in: {comp_dir}")

        for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
            phase_name = phase_dir.name

            for poscar_path in sorted(phase_dir.rglob("POSCAR")):
                print(f"  Found POSCAR: {poscar_path}")

                try:
                    lattice, natoms, counts_map = self.poscar_lattice_and_counts(poscar_path)
                except Exception as e:
                    print(f"  Read failed: {poscar_path} -> {e}")
                    continue

                sqs_level, a_fracs = self.parse_sqs_meta(poscar_path)
                vol = float(abs(np.linalg.det(lattice)))
                vpa = vol / natoms if natoms else None
                a, b, c, alpha, beta, gamma = self.cellpar_from_lattice(lattice)

                rows.append(
                    {
                        "composition_folder": comp_name,
                        "phase_folder": phase_name,
                        "sqs_level": sqs_level,
                        "sqs_a_fracs_json": json.dumps(a_fracs, sort_keys=True),
                        "poscar_path": str(poscar_path),
                        "volume_A3": vol,
                        "natoms": natoms,
                        "volume_per_atom_A3": vpa,
                        "a_A": a,
                        "b_A": b,
                        "c_A": c,
                        "alpha_deg": alpha,
                        "beta_deg": beta,
                        "gamma_deg": gamma,
                        "poscar_counts_json": json.dumps(counts_map, sort_keys=True),
                    }
                )

        return pd.DataFrame(rows)
