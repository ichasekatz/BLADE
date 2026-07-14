"""BLADE analysis subpackage — structure visualization and data extraction.

Classes are lazily imported on first access.

Example::

    from blade.analysis import BladeVisualizer, BLADEData
    from blade.analysis import BLADEVolume  # deprecated alias for BLADEData
"""

from __future__ import annotations

import importlib

__author__ = "Chase Katz"

_ANALYSIS_MAP: dict[str, tuple[str, str]] = {
    "BladeVisualizer": ("blade.analysis.blade_visual", "BladeVisualizer"),
    "BLADEData": ("blade.analysis.blade_data", "BLADEData"),
    "BLADEVolume": ("blade.analysis.blade_data", "BLADEVolume"),  # deprecated alias
}

__all__ = list(_ANALYSIS_MAP)


def __getattr__(name: str) -> type:
    if name in _ANALYSIS_MAP:
        module_path, class_name = _ANALYSIS_MAP[name]
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
