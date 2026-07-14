"""Refit TDB from existing relaxed structures with updated terms.in / mult.in.

Skips sqs2tdb -cp and MLIP relaxation entirely — assumes sqs_lev=* dirs and
energy files already exist in comp_dir. Writes new terms.in / mult.in, runs
sqs2tdb -fit twice, then sqs2tdb -tdb.

Usage:
    pixi run python examples/extra/refit_tdb.py
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from pycalphad import Database

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"

# Root output directory containing per-composition folders (e.g. Files/Comps/NbVW/)
comps_dir = path1 / "Files" / "Comps"

# ------------------------------------------------------------------
# Plot settings — set update_plots=True to regenerate plots after refit
# ------------------------------------------------------------------
update_plots = True
plot_temps = [300, 1000, 2000, 3000, 4000]
phase_T_range = (300, 4500, 50)
fixed_elements = {""}  # fixed sublattice species (excluded from plot axes), {""} if none
metal_fraction = 1  # fraction of atoms that are metal (1/3 for MB2)

# ------------------------------------------------------------------
# Fit parameters (must match original run)
# ------------------------------------------------------------------
t_min = False  # 298.15
t_max = False  # 10000.0
bv = 5e-3
sro = True

# ------------------------------------------------------------------
# New terms.in / mult.in content
# Keys are lattice base names (e.g. "HEDB1", "BCCmcsqs").
# Set to None to use sqs2tdb's default from the first -fit pass.
# ------------------------------------------------------------------
terms_in: dict[str, str] | None = {
    "HEDB1_1": "1,0:1,0\n2,2:1,0\n",  # unary   (e.g. CrB)
    "HEDB1_2": "1,0:1,0\n2,2:1,0\n",  # binary   (e.g. CrHfB)
    "HEDB1_3": "1,0:1,0\n2,2:1,0\n3,0:1,0\n",  # ternary  (e.g. CrHfNbB)
    "BCC1_3": "1,0\n2,2\n3,0\n",  # ternary  (e.g. CrHfNb)
}

mult_in: dict[str, str] | None = None
# mult_in = {"HEDB1": "a=1\tb=2"}

# ------------------------------------------------------------------
# Which compositions to refit (None = all found in comps_dir)
# ------------------------------------------------------------------
compositions: list[str] | None = None
# compositions = ["NbVW", "CrHf"]


# ------------------------------------------------------------------
# Refit loop
# ------------------------------------------------------------------
def run(cmd: list[str], cwd: Path) -> None:
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout.strip())
    if result.returncode not in (0, 1):  # sqs2tdb exits 1 on first -fit pass (expected)
        print(f"  stderr: {result.stderr.strip()}")


comp_dirs = sorted(
    [d for d in comps_dir.iterdir() if d.is_dir()] if compositions is None else [comps_dir / c for c in compositions]
)

for comp_dir in comp_dirs:
    if not comp_dir.is_dir():
        print(f"Skipping missing: {comp_dir}")
        continue

    print(f"\n{'='*50}")
    print(f"Refitting: {comp_dir.name}")

    # Each lattice is a subdirectory containing sqs_lev=* dirs
    lattice_dirs = [d for d in comp_dir.iterdir() if d.is_dir() and not d.name.startswith(".")]

    if not lattice_dirs:
        print("  No lattice dirs found, skipping")
        continue

    original_dir = Path.cwd()
    os.chdir(comp_dir)

    try:
        for lat_dir in sorted(lattice_dirs):
            lattice_base = lat_dir.name.rsplit("_", 1)[0]

            # Check sqs_lev dirs with energy files exist
            sqs_dirs = list(lat_dir.glob("sqs_lev=*/energy"))
            if not sqs_dirs:
                print(f"  {lat_dir.name}: no energy files found, skipping")
                continue

            print(f"  Lattice: {lat_dir.name} ({len(sqs_dirs)} energy files)")

            # Build -fit args
            fit_args = ["sqs2tdb", "-fit"]
            if t_min:
                fit_args.append(f"-Tl={t_min}")
            if t_max:
                fit_args.append(f"-Tu={t_max}")
            if bv:
                fit_args.append(f"-bv={bv}")
            if sro:
                fit_args.append("-sro")

            # First -fit pass — generates default terms.in
            print("    sqs2tdb -fit (pass 1)")
            run(fit_args, cwd=lat_dir)

            # Write custom terms.in — full name (e.g. "HEDB1_2") takes priority over base ("HEDB1")
            if terms_in:
                content = terms_in.get(lat_dir.name) or terms_in.get(lattice_base)
                if content:
                    (lat_dir / "terms.in").write_text(content)
                    print(f"    Wrote terms.in [{lat_dir.name}]: {content.strip()!r}")

            # Write custom mult.in — same priority logic
            if mult_in:
                m_content = mult_in.get(lat_dir.name) or mult_in.get(lattice_base)
                if m_content:
                    (lat_dir / "mult.in").write_text(m_content)
                    print(f"    Wrote mult.in [{lat_dir.name}]")

            # Second -fit pass — uses updated terms.in / mult.in
            print("    sqs2tdb -fit (pass 2)")
            run(fit_args, cwd=lat_dir)

        # Generate TDB
        print("  sqs2tdb -tdb")
        run(["sqs2tdb", "-tdb"], cwd=comp_dir)

        tdbs = list(comp_dir.glob("*.tdb"))
        if tdbs:
            print(f"  TDB written: {[t.name for t in tdbs]}")
        else:
            print("  Warning: no TDB file found after -tdb")

        # Regenerate plots
        if update_plots and tdbs:
            from blade.analysis.blade_visual import BladeVisualizer

            viz = BladeVisualizer()

            for tdb_path in tdbs:
                try:
                    tdb = Database(str(tdb_path))
                except Exception as e:
                    print(f"  Failed to load TDB for plotting: {e}")
                    continue

                # Derive metals and fixed_species from TDB elements
                all_els = [el for el in tdb.elements if el != "VA"]
                metals = sorted(el.title() for el in all_els if el.upper() not in {f.upper() for f in fixed_elements})
                fixed_sp = {el.title(): (1.0 - metal_fraction) / max(len(fixed_elements), 1) for el in fixed_elements}
                phases = list(tdb.phases.keys())
                comp_name = "".join(metals)

                if len(metals) < 2:
                    continue

                print(f"  Plotting: {comp_name}  metals={metals}  phases={phases}")

                viz.plot_gibbs_energy(
                    tdb=tdb,
                    metals=metals,
                    phase=phases[0],
                    fixed_species=fixed_sp,
                    temperatures=plot_temps,
                    output_path=comp_dir / f"{comp_name}_Gibbs_Energy.png",
                )
                viz.plot_gibbs_mixing(
                    tdb=tdb,
                    metals=metals,
                    phase=phases[0],
                    fixed_species=fixed_sp,
                    temperatures=plot_temps,
                    output_path=comp_dir / f"{comp_name}_Gibbs_Mixing.png",
                )
                if len(metals) == 2:
                    viz.plot_binary_phase_diagram(
                        tdb=tdb,
                        metals=metals,
                        phases=phases,
                        fixed_species=fixed_sp,
                        temperature_range=phase_T_range,
                        output_path=comp_dir / f"{comp_name}_Phase_Diagram.png",
                    )
                print(f"  Plots updated for {comp_name}")

    finally:
        os.chdir(original_dir)

print("\nDone.")
