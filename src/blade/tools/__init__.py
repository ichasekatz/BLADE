"""BLADE tools subpackage — composition generation, SQS input creation, and TDB fitting.

Classes are lazily imported on first access.

Example::

    from blade.tools import BladeCompositions, BladeSQS, BladeTDBGen
"""

from __future__ import annotations

import importlib

__author__ = "Chase Katz"

_TOOLS_MAP: dict[str, tuple[str, str]] = {
    "BladeCompositions": ("blade.tools.blade_compositions", "BladeCompositions"),
    "BladeSQS": ("blade.tools.blade_sqsgen", "BladeSQS"),
    "BladeTDBGen": ("blade.tools.blade_tdb_gen", "BladeTDBGen"),
    "BladeCutoff": ("blade.tools.blade_cutoff", "BladeCutoff"),
}

__all__ = list(_TOOLS_MAP)


def __getattr__(name: str) -> type:
    """Lazily import and return a tool class by name.

    Args:
        name (str): The class name to look up.

    Returns:
        type: The requested class.

    Raises:
        AttributeError: If ``name`` is not found in the tools map.
    """
    if name in _TOOLS_MAP:
        module_path, class_name = _TOOLS_MAP[name]
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
