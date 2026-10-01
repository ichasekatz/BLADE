"""SQS generation with BladeSQS.

Special quasi-random structures (SQS) are supercell models that reproduce the
pair- and multi-body correlation functions of a random alloy up to a chosen
cutoff. BLADE uses ATAT (mcsqs / corrdump / sqs2tdb) to generate one SQS per
composition level defined in sqsgen_levels.

When to use this script
-----------------------
Run this before TDB fitting (03_tdb_generation.py).  You need it whenever:
  - a new phase prototype is added to phases_dict, or
  - the element set or desired composition grid changes, or
  - the existing bestcorr.out files were produced with different mcsqs settings.

Prerequisites
-------------
  - ATAT binaries (mcsqs, corrdump, sqs2tdb) on PATH.
  - The BLADE and sqsdb directories must exist (path1, path2 below).
"""

from pathlib import Path

from blade.tools.blade_sqsgen import BladeSQS

# ---------------------------------------------------------------------------
# Paths — only sqsdb_dir needs to be set explicitly; blade_root is auto-detected
# ---------------------------------------------------------------------------
blade_root = Path(__file__).parent.parent
sqsdb_dir = Path("/path/to/PhaseForge/atat/data/sqsdb")  # adjust this
paths = [blade_root.parent, blade_root, sqsdb_dir]

# ---------------------------------------------------------------------------
# AlB2-type diboride phase prototype (HEDB1 standard example)
# Two sublattices: metal site (a) and boron site (B fixed).
# ---------------------------------------------------------------------------
phases: dict[str, dict] = {
    "HEDB1": {
        "a": 3.14,  # Å — average over active metals; refine as needed
        "b": 3.14,
        "c": 3.53,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"  # metal site — mixed
            "0.333333 0.666667 0.500000 B\n"  # boron — fixed
            "0.666667 0.333333 0.500000 B\n"
        ),
    },
}

# Phase list: one entry per (phase prototype, sqsdb lattice label, supercell).
# supercell_size controls total atom count: unit_cell_sites × product(dims).
# (2, 2, 2) → 24 atoms for a 3-site AlB2 cell; increase for better SQS quality.
phase_list = [
    {
        "generator_name": "HEDB1",
        "lattice": "HEDB1",
        "supercell_size": (2, 2, 2),
    },
]

# ---------------------------------------------------------------------------
# Composition levels
# Level 0: endmember (pure-A on the metal site) — required baseline.
# Level 1: equimolar binary alloy.
# ---------------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]], "letter": ["a"]},
    {"level": 1, "compositions": [[0.5, 0.5]], "letter": ["a"]},
]

# Maximum level index to write into sqsgen.in (inclusive).
level = 1

# ---------------------------------------------------------------------------
# mcsqs run parameters
# ---------------------------------------------------------------------------
mcsqs_params = {
    "time": 30,  # seconds per sqsdb directory
    "cutoff_mode": "nn",  # "nn" = nearest-neighbour shell index
    "2": 4,  # pair cutoff: 4th-NN shell
    "3": 3,  # triplet cutoff
    "4": 2,  # quadruplet cutoff
    "wr": 20,
    "wn": 0.75,
    "wd": 1,
    "parallel_runs": 8,
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    for phase in phase_list:
        # Binary system: len_comp=2 (one metal + one endmember placeholder).
        for len_comp in [2]:
            sqs = BladeSQS(
                phases_dict=phases[phase["generator_name"]],
                sqsgen_levels=sqsgen_levels,
                level=level,
                len_comp=len_comp,
                skip_existing_sqs=False,
            )
            params = mcsqs_params | {"super_cell_size": phase["supercell_size"]}
            sqs.sqs_gen(phase=phase, paths=paths, params=params)


if __name__ == "__main__":
    main()
