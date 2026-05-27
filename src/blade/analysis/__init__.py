"""BLADE analysis subpackage — structure visualization and volume extraction.

Classes are lazily imported on first access.

Example::

    from blade.analysis import BladeVisualizer, BLADEVolume
"""

from __future__ import annotations

import importlib

__author__ = "Chase Katz"

_ANALYSIS_MAP: dict[str, tuple[str, str]] = {
    "BladeVisualizer": ("blade.analysis.blade_visual", "BladeVisualizer"),
    "BLADEVolume": ("blade.analysis.blade_volume", "BLADEVolume"),
}

__all__ = list(_ANALYSIS_MAP)


def __getattr__(name: str) -> type:
    """Lazily import and return an analysis class by name.

    Args:
        name (str): The class name to look up.

    Returns:
        type: The requested class.

    Raises:
        AttributeError: If ``name`` is not found in the analysis map.
    """
    if name in _ANALYSIS_MAP:
        module_path, class_name = _ANALYSIS_MAP[name]
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
