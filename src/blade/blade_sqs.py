"""
This module defines the `BladeSQS` class for generating ATAT SQS inputs and running mcsqs workflows.

The `BladeSQS` class writes ATAT input files (e.g., rndstr.skel and sqsgen.in), runs sqs2tdb setup,
and then executes corrdump and mcsqs for each generated composition directory. It also provides logic
for choosing a supercell size that matches requested fractional compositions on the 'a' sublattice.
"""

import math
import subprocess
import threading
from fractions import Fraction
from pathlib import Path


class BladeSQS:
    """
    Generate SQS inputs and run ATAT mcsqs for a given phase prototype.

    The phase prototype is specified via a phases_dict containing lattice parameters and a unit-cell
    coordinate string. The class can:
      - build ATAT rndstr inputs (rndstr.skel)
      - write sqsgen.in from preset levels
      - run sqs2tdb -mk to create sqsdb_lev=* folders
      - compute a suitable n_atoms for mcsqs based on target fractions
      - run corrdump + mcsqs with an automatic timeout via stopsqs
    """
    def __init__(self, phases_dict, sqsgen_levels, level):
        """
        Initializes the BladeSQS object.

        Args:
            phases_dict (dict[str, object]): Dictionary describing a phase prototype.
                Must include: "a", "b", "c", "alpha", "beta", "gamma", and "coords".
            sqsgen_levels (list[str]): Lines for sqsgen.in, where higher indices correspond to higher optimization levels.
            level (int): Maximum sqsgen level to include (inclusive). For example, level=2 will include levels [0, 1, 2].
        """
        self.a = phases_dict["a"]
        self.b = phases_dict["b"]
        self.c = phases_dict["c"]
        self.alpha = phases_dict["alpha"]
        self.beta = phases_dict["beta"]
        self.gamma = phases_dict["gamma"]
        self.unit_cell = phases_dict["coords"]
        self.sqsgen_levels = sqsgen_levels
        self.level = level

    def sqs_struct(self):
        """
        Build the sqsgen.in text and the rndstr.skel text.

        Returns:
            tuple[str, str]:
                (sqsgen_text, rndstr_text)
        """
        rndstr1 = f"""
        {self.a} {self.b} {self.c} {self.alpha} {self.beta} {self.gamma}
        1 0 0
        0 1 0
        0 0 1
        """
        sqsgen = ""
        for i in range(self.level + 1):
            sqsgen += self.sqsgen_levels[i] + "\n"

        rndstr = rndstr1.strip() + "\n" + self.unit_cell.strip()
        print(rndstr)
        self.sqsgen_text = sqsgen
        self.rndstr = rndstr

        return sqsgen, rndstr

    def supercell_size(self, fractions, min_a_sites=16, max_den=96):
        """
        Determine a supercell size compatible with target fractions.

        This computes:
          - the number of 'a' sites needed to represent fractions exactly (up to max_den)
          - a prototype replication factor k such that the total sites follow the (a_sites + b_sites) ratio
          - integer species counts on the 'a' sublattice

        Args:
            fractions (list[float]): Species fractions on the 'a' sublattice.
            min_a_sites (int, optional): Minimum number of 'a' sites to enforce for quality. Defaults to 16.
            max_den (int, optional): Maximum denominator used when rationalizing fractions. Defaults to 96.

        Returns:
            tuple[int, list[int]]:
                (n_total, counts) where:
                  - n_total is the total number of atoms (a + b sublattices)
                  - counts are integer counts per species on the 'a' sublattice

        Raises:
            ValueError:
                If the unit cell contains no 'a' sites.
        """
        fracs = [Fraction(float(f)).limit_denominator(max_den) for f in fractions]
        denom_lcm = math.lcm(*[f.denominator for f in fracs])

        a_sites = 0
        b_sites = 0
        sites = 0
        for line in self.unit_cell.splitlines():
            parts = line.split()
            if len(parts) != 4:
                continue
            sites += 1
            if parts[-1] == "a":
                a_sites += 1
            else:
                b_sites += 1

        if a_sites == 0:
            raise ValueError("No 'a' sites found in unit_cell")

        n_a0 = math.lcm(a_sites, denom_lcm)

        if n_a0 < min_a_sites:
            mult = math.ceil(min_a_sites / n_a0)
            n_a = n_a0 * mult
        else:
            n_a = n_a0

        k = n_a // a_sites

        n_total = k * (a_sites + b_sites)

        counts = [int(n_a * f) for f in fracs]
        diff = n_a - sum(counts)
        counts[0] += diff

        return n_total, counts

    def sqs_gen(self, unique_len_comps, phase, path1, iter):
        """
        Generate SQS folders and run corrdump + mcsqs for each composition.

        For each `length` in unique_len_comps, this will:
          1) create <path1>/<phase>_<length>/
          2) write rndstr.skel and sqsgen.in
          3) run `sqs2tdb -mk` to create sqsdb_lev=* subfolders
          4) for each sqsdb_lev=* folder, run corrdump then mcsqs -n=<n_atoms>
             with a timeout implemented by creating a `stopsqs` file

        Args:
            unique_len_comps (iterable[int]): Set or list of system sizes (e.g., {2, 3}).
            phase (str): Phase name used to construct output folder names.
            path1 (str | Path): Root directory where phase folders are created.
            time (float): Seconds to wait before touching `stopsqs` to stop mcsqs.
        """
        for length in unique_len_comps:
            dir_name = Path(path1) / (phase + "_" + str(length))
            dir_name.mkdir(parents=True, exist_ok=True)
            file_path = dir_name / "rndstr.skel"
            sqsgen, rndstr = self.sqs_struct()
            with file_path.open("w") as f:
                f.write(rndstr)
            print(f"File created at: {file_path}")

            dir_name = Path(path1) / (phase + "_" + str(length))
            dir_name.mkdir(parents=True, exist_ok=True)
            file_path = dir_name / "sqsgen.in"
            with file_path.open("w") as f:
                f.write(sqsgen)
            print(f"File created at: {file_path}")

            cmd = ["sqs2tdb", "-mk"]
            result = subprocess.run(cmd, cwd=dir_name, capture_output=True, text=True, check=False)
            print(result.stdout)
            if result.stderr:
                print("Error:", result.stderr)

            parent_dir = Path(path1) / (phase + "_" + str(length))
            for sqsdir in parent_dir.glob("sqsdb_lev=*/"):
                folder_name = sqsdir.name
                folder_lev = folder_name.split("=")[1]
                folder_lev = folder_lev.split("_")[0]
                comp_str = folder_name.split("=")[-1]  # e.g. "0.75,0.125,0.125"
                fractions = [float(x) for x in comp_str.split(",")]
                if len(fractions) == 1:  # Only one species, no mixing
                    print(fractions)
                    print(f"Skipping pure species directory: {sqsdir}")
                    continue
                n_atoms, counts = self.supercell_size(fractions)
                print(f"Running corrdump in {fractions}")
                try:
                    subprocess.run(
                        [
                            "corrdump",
                            "-l=rndstr.in",
                            "-ro",
                            "-noe",
                            "-nop",
                            "-clus",
                            "-2=1,2,3",
                            "-3=1",
                        ],
                        cwd=sqsdir,
                        check=True,
                    )

                    print(f"Running mcsqs with {n_atoms} atoms in {fractions}")
                    subprocess.run(
                        ["mcsqs", f"-n={n_atoms}", f"-ms={iter}"],
                        cwd=sqsdir,
                        check=True,
                    )

                except subprocess.CalledProcessError:
                    print(f"corrdump or mcsqs failed in {sqsdir}, continuing.")
                    continue