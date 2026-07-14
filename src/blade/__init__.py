"""BLADE — Batch Lattice Analysis and Discovery Engine.

High-level lazy imports for the most commonly used classes.

Example::

    from blade.tools import BladeCompositions, BladeSQS, BladeTDBGen
    from blade.analysis import BladeVisualizer, BLADEData
    from blade.oxidation import OxideCompositions, OxideDatabase
"""

from __future__ import annotations

import importlib

__version__ = "1.6"
__author__ = "Chase Katz"

_MAP: dict[str, tuple[str, str]] = {
    # tools
    "BladeCompositions": ("blade.tools.blade_compositions", "BladeCompositions"),
    "BladeSQS":          ("blade.tools.blade_sqsgen",       "BladeSQS"),
    "BladeTDBGen":       ("blade.tools.blade_tdb_gen",      "BladeTDBGen"),
    "BladeCutoff":       ("blade.tools.blade_cutoff",       "BladeCutoff"),
    "ScrapsSQSGen":      ("blade.tools.blade_scraps_gen",   "ScrapsSQSGen"),
    # analysis
    "BladeVisualizer":   ("blade.analysis.blade_visual",    "BladeVisualizer"),
    "BLADEData":         ("blade.analysis.blade_data",      "BLADEData"),
    "BLADEVolume":       ("blade.analysis.blade_data",      "BLADEVolume"),
    # oxidation
    "OxideCompositions": ("blade.oxidation.compositions",   "OxideCompositions"),
    "OxideDatabase":     ("blade.oxidation.database",       "OxideDatabase"),
}

__all__ = list(_MAP)


def __getattr__(name: str) -> type:
    if name in _MAP:
        module_path, class_name = _MAP[name]
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
