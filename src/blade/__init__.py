"""BLADE — high-throughput materials discovery and thermodynamic database generation.

Supports any crystal structure with fixed and variable sublattice sites.

Classes are lazily imported on first access to avoid loading heavy dependencies
(torch, materialsframework, pymatgen) at import time.

Example::

    from blade import BladeCompositions, BladeSQS, BladeTDBGen
    from blade import BladeVisualizer, BLADEVolume
"""

from __future__ import annotations

import importlib

__version__ = "0.2.0"
__author__ = "Chase Katz"

_CLASS_MAP: dict[str, tuple[str, str]] = {
    "BladeCompositions": ("blade.tools.blade_compositions", "BladeCompositions"),
    "BladeSQS": ("blade.tools.blade_sqsgen", "BladeSQS"),
    "BladeTDBGen": ("blade.tools.blade_tdb_gen", "BladeTDBGen"),
    "BladeCutoff": ("blade.tools.blade_cutoff", "BladeCutoff"),
    "BladeVisualizer": ("blade.analysis.blade_visual", "BladeVisualizer"),
    "BLADEVolume": ("blade.analysis.blade_volume", "BLADEVolume"),
}

__all__ = list(_CLASS_MAP)


def __getattr__(name: str) -> type:
    """Lazily import and return a class by name.

    Args:
        name (str): The class name to look up.

    Returns:
        type: The requested class.

    Raises:
        AttributeError: If ``name`` is not a known BLADE class.
    """
    if name in _CLASS_MAP:
        module_path, class_name = _CLASS_MAP[name]
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
