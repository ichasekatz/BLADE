"""BLADE oxidation subpackage — oxide composition enumeration and database construction.

Classes are lazily imported on first access.

Example::

    from blade.oxidation import OxideCompositions, OxideDatabase
"""

from __future__ import annotations

import importlib

__author__ = "Chase Katz"

_OXIDATION_MAP: dict[str, tuple[str, str]] = {
    "OxideCompositions": ("blade.oxidation.compositions", "OxideCompositions"),
    "OxideDatabase":     ("blade.oxidation.database",     "OxideDatabase"),
}

__all__ = list(_OXIDATION_MAP)


def __getattr__(name: str) -> type:
    if name in _OXIDATION_MAP:
        module_path, class_name = _OXIDATION_MAP[name]
        module = importlib.import_module(module_path)
        return getattr(module, class_name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
