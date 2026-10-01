"""Demonstrate BladeCompositions: enumerate N-component element systems.

BladeCompositions generates all valid combinations of primary (and optional
secondary) elements subject to min/max count constraints. The resulting list
feeds directly into BladeSQS and BladeTDBGen.
"""

from __future__ import annotations

from blade.tools.blade_compositions import BladeCompositions

# ── Element pool ──────────────────────────────────────────────────────────────
primary_elements = ["Cr", "Hf", "Mo", "Ta", "Zr"]

if __name__ == "__main__":
    # Binary: exactly 2 primary elements, no secondary
    binary = BladeCompositions(
        primary_elements=primary_elements,
        secondary_elements=[],
        primary_min=2,
        primary_max=2,
        secondary_min=0,
        secondary_max=0,
    )
    binary_comps = binary.generate_compositions()
    print(f"Binary  — systems: {binary.get_systems()}, total: {len(binary_comps)}")
    print(f"  {binary_comps}\n")

    # Ternary: exactly 3 primary elements
    ternary = BladeCompositions(
        primary_elements=primary_elements,
        secondary_elements=[],
        primary_min=3,
        primary_max=3,
        secondary_min=0,
        secondary_max=0,
    )
    ternary_comps = ternary.generate_compositions()
    print(f"Ternary — systems: {ternary.get_systems()}, total: {len(ternary_comps)}")
    print(f"  {ternary_comps[:5]}{'...' if len(ternary_comps) > 5 else ''}\n")

    # Mixed binary + ternary, filter to binaries only
    mixed = BladeCompositions(
        primary_elements=["Cr", "Hf", "Mo", "Ta"],
        secondary_elements=[],
        primary_min=2,
        primary_max=3,
        secondary_min=0,
        secondary_max=0,
    )
    all_comps = mixed.generate_compositions()
    binaries = [c for c in all_comps if len(c) == 2]
    print(f"Mixed   — all: {len(all_comps)}, binaries only: {len(binaries)}")
    print(f"  {binaries}")
