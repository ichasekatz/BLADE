"""Example: structure and phase-diagram visualization with BladeVisualizer.

This script shows how to use BladeVisualizer for five common tasks:

1. Combining CONTCAR structures from multiple SQS levels into one image.
2. Plotting Gibbs energy vs metal-site composition at multiple temperatures.
3. Plotting Gibbs energy of mixing (ΔG_mix) vs composition.
4. Plotting a pseudo-binary CALPHAD phase diagram (T vs X).
5. Stitching phase-diagram PNGs side by side (after generating them above).
6. Building a relaxation movie from trajectory files (requires ffmpeg for MP4).

Prerequisites:
    - matplotlib, ase, Pillow, imageio, and pycalphad must be installed.
    - CONTCAR / PNG / TDB files must exist at the paths used below.
    - ffmpeg is optional (needed for MP4 output only).

    python examples/04_visualization.py
"""

from pathlib import Path

from pycalphad import Database

from blade.analysis.blade_visual import BladeVisualizer

viz = BladeVisualizer()

# ------------------------------------------------------------------
# Paths — adjust to your environment
# ------------------------------------------------------------------
path1 = Path("/Users/chasekatz/Desktop/School/Research/BLADE/Files/Comps/")

# ------------------------------------------------------------------
# Example 1: combine CONTCAR structures for one composition + phase
# ------------------------------------------------------------------
print("=" * 60)
print("Example 1: Combined CONTCAR structures")
print("=" * 60)

comp_name = "CrHf"
phase_name = "ALLOY1_1_2"
comp_dir = path1 / comp_name
phase_dir = comp_dir / phase_name

contcars = sorted(phase_dir.glob("sqs_lev=*/CONTCAR"))
if contcars:
    out = comp_dir / f"Combined_CONTCARs_{comp_name}_{phase_name}.png"
    viz.contcar(contcars, save=out)
    print(f"Saved {len(contcars)} structures -> {out}")
else:
    print(f"No CONTCARs found in {phase_dir} — skipping.")

# ------------------------------------------------------------------
# Example 2: Gibbs energy vs composition (single-phase, multi-T)
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 2: Gibbs energy vs composition")
print("=" * 60)

comp_name = "CrHf"
phase_name = "ALLOY2_1_2"

# fixed_species: elements held at constant total mole fraction
# e.g. {"Zr": 1/2} for ALLOY3 in 03_tdb_generation.py, where Zr is fixed at 50% and CrHf varies from 0-50%
fixed_species = {"Zr": 1 / 2}
# fixed_species = {}

tdb_files = sorted((path1 / comp_name).glob("*.tdb"))
if tdb_files:
    tdb = Database(str(tdb_files[0]))
    viz.plot_gibbs_energy(
        tdb=tdb,
        metals=["Cr", "Hf"],
        phase=phase_name,
        fixed_species=fixed_species,
        temperatures=[300, 1000, 2000, 3000, 4000],
        output_path=path1 / comp_name / f"{comp_name}_Gibbs_Energy.png",
    )
else:
    print(f"No TDB files found in {path1 / comp_name} — skipping.")

# ------------------------------------------------------------------
# Example 3: Gibbs energy of mixing (ΔG_mix vs composition)
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 3: Gibbs energy of mixing")
print("=" * 60)

comp_name = "CrHf"
phase_name = "ALLOY1_1_2"

if tdb_files:
    viz.plot_gibbs_mixing(
        tdb=tdb,
        metals=["Cr", "Hf"],
        phase=phase_name,
        fixed_species=fixed_species,
        temperatures=[300, 1000, 2000, 3000, 4000],
        output_path=path1 / comp_name / f"{comp_name}_Gibbs_Mixing.png",
    )
else:
    print(f"No TDB files found in {path1 / comp_name} — skipping.")

# ------------------------------------------------------------------
# Example 4: pseudo-binary CALPHAD phase diagram (T vs X)
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 4: Pseudo-binary phase diagram")
print("=" * 60)

comp_name = "CrHf"
phase_name = "ALLOY1_1_2"

if tdb_files:
    viz.plot_binary_phase_diagram(
        tdb=tdb,
        metals=["Cr", "Hf"],
        phases=[phase_name],
        fixed_species=fixed_species,
        temperature_range=(300, 4500, 50),
        output_path=path1 / comp_name / f"{comp_name}_Phase_Diagram.png",
    )
else:
    print(f"No TDB files found in {path1 / comp_name} — skipping.")

# ------------------------------------------------------------------
# Example 5: combine phase diagram PNGs (generated above)
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 5: Combined phase diagrams")
print("=" * 60)

pngs = sorted(path1.glob("**/*_Phase_Diagram.png"))
if pngs:
    out_pd = path1 / "Combined_Phase_Diagrams.png"
    viz.phase_diagram(pngs, save=out_pd)
    print(f"Combined {len(pngs)} phase diagrams -> {out_pd}")
else:
    print("No *_Phase_Diagram.png files found — skipping.")

# ------------------------------------------------------------------
# Example 6: relaxation movie (GIF + optional MP4)
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 6: Relaxation movie")
print("=" * 60)

composition_list = [["Cr", "Hf"]]

viz.make_combined_relaxation_movie(
    composition_list=composition_list,
    path1=path1,
    traj_name="relaxation_live.xyz",
    fps=10,
)
