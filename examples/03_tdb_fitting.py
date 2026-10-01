"""Demonstrate BladeTDBGen: drive MLIP-backed CALPHAD TDB fitting.

BladeTDBGen orchestrates MaterialsFramework's Sqs2tdb workflow across
many chemical compositions.  The constructor is side-effect-free; all
computation runs when .fit() is called.

tdb_params controls the GraceCalculator (or other MLIP) device, the
relaxation convergence criterion (fmax), the temperature range for
CALPHAD fitting, and CALPHAD model flags (sro, phonon, bv, etc.).

Adjust the path variables before running.
"""

from __future__ import annotations

from pathlib import Path

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_tdb_gen import BladeTDBGen

# ── Paths ─────────────────────────────────────────────────────────────────────
blade_root = Path(__file__).parent.parent
sqsdb_dir = Path("/path/to/PhaseForge/atat/data/sqsdb")  # adjust this

paths = [blade_root.parent, blade_root, sqsdb_dir]
comps_dir = blade_root / "Files" / "Comps"

# ── Phase and composition definitions ────────────────────────────────────────
phases: dict[str, dict] = {
    "HEDB1": {
        "a": 3.0624,
        "b": 3.0624,
        "c": 3.2392,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": ("0.000000 0.000000 0.000000 a\n0.333333 0.666667 0.500000 B\n0.666667 0.333333 0.500000 B\n"),
    },
}

phase_list = [
    {"generator_name": "HEDB", "lattice": "HEDB1", "supercell_size": (4, 3, 2)},
]

# ── MLIP and fitting parameters ───────────────────────────────────────────────
# calculator: "orb" | "grace" | "mace" | "uma" | "chgnet" (see MaterialsFramework)
tdb_params: dict = {
    "fmax": 1e-4,
    "verbose": True,
    "calculator": "orb",
    "calculator_kwargs": {"steps": 1000, "device": "cpu"},
    "t_min": 298.15,
    "t_max": 10000.0,
    "sro": False,
    "bv": 1e-3,
    "phonon": False,
    "open_calphad": False,
    "track_trajectory": False,
    "terms": None,
}

# Per-phase CALPHAD interaction parameter model (ATAT terms_in format).
terms_in: dict[str, str] = {
    "HEDB1": "1,0:1,0\n2,2:1,0\n",
    "HEDB1_2": "1,0:1,0\n2,2:1,0\n",
}


def run_tdb(elements: list[str]) -> None:
    """Fit TDB databases for all compositions of a binary system.

    Args:
        elements: Primary elements to combine, e.g. ["Cr", "Hf"].

    Returns:
        None. TDB files are written to Files/Comps/<system>/<phase>/.
    """
    composer = BladeCompositions(
        primary_elements=elements,
        secondary_elements=[],
        primary_min=2,
        primary_max=2,
        secondary_min=0,
        secondary_max=0,
    )
    composition_list = composer.generate_compositions()

    gen = BladeTDBGen(
        phases=phase_list,
        phases_dict=phases,
        liquid=False,
        paths=paths,
        composition_list=composition_list,
        level=1,
        skip_existing=False,
        refit_existing=False,
        output_dir=comps_dir,
        terms_in=terms_in,
        tdb_params=tdb_params,
    )
    gen.fit()

    # Output TDB files are at:
    for comp in composition_list:
        system = "".join(comp)
        tdb_path = comps_dir / system
        if tdb_path.exists():
            tdbs = list(tdb_path.rglob("*.tdb"))
            print(f"{system}: {len(tdbs)} TDB file(s) — {[t.name for t in tdbs]}")


if __name__ == "__main__":
    run_tdb(["Cr", "Hf"])
