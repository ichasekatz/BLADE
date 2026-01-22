import os
import subprocess
import math
from pathlib import Path
import threading
from fractions import Fraction

def supercell_size(self, fractions):
    fracs = [Fraction(str(f)).limit_denominator() for f in fractions]
    denom_lcm = math.lcm(*[f.denominator for f in fracs])

    # Count only actual atomic site lines: "x y z label"
    a_sites = 0
    b_sites = 0
    sites = 0
    for line in self.rndstr.splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        label = parts[-1]
        sites += 1
        if label == "B":
            b_sites += 1
        else:
            a_sites += 1

    base_sites = a_sites
    count_B = b_sites

    # Choose a number of 'a' sites that can represent the fractions exactly
    n_a = math.lcm(base_sites, denom_lcm)

    if n_a >= 100:
        n_a = base_sites * len(fracs)

    n_total = n_a + count_B
    n_total = ((n_total + (sites-1)) // sites) * sites

    n_a = n_total - count_B

    counts = [int(n_a * f) for f in fracs]
    diff = n_a - sum(counts)
    counts[0] += diff

    return n_total, counts

def sqs_struct(self):
    rndstr1 = f"""
    {self.a} {self.b} {self.c} {self.alpha} {self.beta} {self.gamma}
    {self.a} 0 0
    {self.a/2} {self.a*math.sqrt(3)/2} 0
    0 0 {self.c}
    """

    sqsgen = ""
    for i in range((self.level+1)):
        sqsgen += self.sqsgen_levels[i] + "\n"

    rndstr = rndstr1.strip() + "\n" + self.unit_cell.strip()
    self.sqsgen_text = sqsgen
    self.rndstr = rndstr

    return sqsgen, rndstr

class BladeSQS:
    def __init__(self, a, b , c, alpha, beta, gamma, unit_cell, sqsgen_levels, level):

        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.unit_cell = unit_cell
        self.sqsgen_levels = sqsgen_levels
        self.level = level

    def sqs_struct(self):
        rndstr1 = f"""
        {self.a} {self.b} {self.c} {self.alpha} {self.beta} {self.gamma}
        {self.a} 0 0
        {self.a/2} {self.a*math.sqrt(3)/2} 0
        0 0 {self.c}
        """

        sqsgen = ""
        for i in range((self.level+1)):
            sqsgen += self.sqsgen_levels[i] + "\n"

        rndstr = rndstr1.strip() + "\n" + self.unit_cell.strip()

        self.sqs_gen_text = sqsgen
        self.rndstr = rndstr

        return sqsgen, rndstr

    def supercell_size(self, fractions):
        fracs = [Fraction(str(f)).limit_denominator() for f in fractions]
        denom_lcm = math.lcm(*[f.denominator for f in fracs])

        # Count only actual atomic site lines: "x y z label"
        a_sites = 0
        b_sites = 0
        sites = 0
        for line in self.rndstr.splitlines():
            parts = line.split()
            if len(parts) != 4:
                continue
            label = parts[-1]
            sites += 1
            if label == "B":
                b_sites += 1
            else:
                a_sites += 1

        base_sites = a_sites
        count_B = b_sites

        # Choose a number of 'a' sites that can represent the fractions exactly
        n_a = math.lcm(base_sites, denom_lcm)

        if n_a >= 100:
            n_a = base_sites * len(fracs)

        n_total = n_a + count_B
        n_total = ((n_total + (sites-1)) // sites) * sites

        n_a = n_total - count_B

        counts = [int(n_a * f) for f in fracs]
        diff = n_a - sum(counts)
        counts[0] += diff

        return n_total, counts


    def sqs_gen(self, unique_len_comps, phase, path1, time):
        # Generate sqs
        for length in unique_len_comps:
            dir_name = os.path.join(path1, (phase+"_"+str(length)))
            os.makedirs(dir_name, exist_ok=True)
            file_path = os.path.join(dir_name, "rndstr.skel")
            sqsgen, rndstr = self.sqs_struct()
            with open(file_path, "w") as f:
                f.write(rndstr)
            print(f"File created at: {file_path}")

            dir_name = os.path.join(path1, (phase+"_"+str(length)))
            os.makedirs(dir_name, exist_ok=True)
            file_path = os.path.join(dir_name, "sqsgen.in")
            with open(file_path, "w") as f:
                f.write(sqsgen)
            print(f"File created at: {file_path}")

            cmd = ['sqs2tdb', '-mk']
            # Run the command in the target directory
            result = subprocess.run(cmd, cwd=dir_name, capture_output=True, text=True)
            # Output results (optional)
            print(result.stdout)
            if result.stderr:
                print("Error:", result.stderr)

            parent_dir = Path(os.path.join(path1, (phase+"_"+str(length))))
            for sqsdir in parent_dir.glob('sqsdb_lev=*/'):
                folder_name = sqsdir.name
                folder_lev = folder_name.split('=')[1]
                folder_lev = folder_lev.split('_')[0]
                comp_str = folder_name.split('=')[-1]   # e.g. "0.75,0.125,0.125"
                fractions = [float(x) for x in comp_str.split(',')]
                # Skip pure-species (end-member) folders
                if len(fractions) == 1:    # Only one species, no mixing
                    print(fractions)
                    print(f"Skipping pure species directory: {sqsdir}")
                    continue
                N_atoms, counts = self.supercell_size(fractions)

                print(f"Running corrdump in {fractions}")
                try:
                    subprocess.run(
                        [
                            "corrdump",
                            "-l=rndstr.in", "-ro", "-noe", "-nop", "-clus",
                            "-2=1,2,3", "-3=1"
                        ],
                        cwd=sqsdir, check=True
                    )

                    print(f"Running mcsqs with {N_atoms} atoms in {fractions}")
                    def make_stopsqs(sqsdir, time, fractions):
                        if not completed.is_set():
                            Path(sqsdir, "stopsqs").touch()
                            print(f"touch stopsqs in {fractions} after {time} seconds")
                    # Start the timer BEFORE you launch mcsqs
                    completed = threading.Event()
                    timer = threading.Timer(time, make_stopsqs(sqsdir, time, fractions))
                    timer.start()

                    try:
                        subprocess.run(
                            ["mcsqs", f"-n={N_atoms}"],
                            cwd=sqsdir,
                            check=True
                        )
                        completed.set()  # Mark process as completed, so timer won't fire
                    except subprocess.CalledProcessError:
                        print(f"corrdump or mcsqs failed in {sqsdir}, continuing.")
                    finally:
                        timer.cancel()
                except subprocess.CalledProcessError:
                    print(f"corrdump or mcsqs failed in {sqsdir}, continuing.")
                    continue
