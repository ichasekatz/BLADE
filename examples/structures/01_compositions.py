"""Example: generating element combinations with BladeCompositions.

This script shows how to enumerate chemical systems for a high-throughput
search. BladeCompositions accepts separate pools of primary and secondary
elements, plus constraints on how many of each to include per system.

Run this script without any dependencies on ATAT, MaterialsFramework, or
a GPU — it only uses the standard library and itertools.

    python examples/01_compositions.py
"""

from blade.tools.blade_compositions import BladeCompositions

# ------------------------------------------------------------------
# Example 1: pure ternary transition-metal systems
# ------------------------------------------------------------------
print("=" * 60)
print("Example 1: Ternary TM systems (no secondary elements)")
print("=" * 60)

composer = BladeCompositions(
    primary_elements=["Cr", "Hf", "Ta", "Ti", "Mo"],
    secondary_elements=[],
    primary_min=3,
    primary_max=3,
    secondary_min=0,
    secondary_max=0,
)
comps = composer.generate_compositions()
systems = composer.get_systems()

print(f"Total compositions: {len(comps)}")
print(f"System sizes:       {systems}")
print(f"First 5:            {comps[:5]}")

# ------------------------------------------------------------------
# Example 2: binary and ternary systems
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 2: Binary and ternary TM systems")
print("=" * 60)

composer2 = BladeCompositions(
    primary_elements=["Zr", "Hf", "Ta", "Cr"],
    secondary_elements=[],
    primary_min=2,
    primary_max=3,
    secondary_min=0,
    secondary_max=0,
)
comps2 = composer2.generate_compositions()
systems2 = composer2.get_systems()

print(f"Total compositions: {len(comps2)}")
print(f"System sizes:       {systems2}")
print(f"All compositions:   {comps2}")

# ------------------------------------------------------------------
# Example 3: mixed TM + rare-earth ternary systems
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 3: Two TMs + one rare earth (ternary)")
print("=" * 60)

composer3 = BladeCompositions(
    primary_elements=["Zr", "Hf", "Ta"],
    secondary_elements=["Y", "La"],
    primary_min=2,
    primary_max=2,
    secondary_min=1,
    secondary_max=1,
)
comps3 = composer3.generate_compositions()

print(f"Total compositions: {len(comps3)}")
print(f"All compositions:   {comps3}")
