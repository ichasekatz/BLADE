"""Standalone utility helpers for the OxideDatabase pipeline.

All functions are pure (no ``self`` dependency) and can be imported
independently of :class:`blade.oxidation.database.OxideDatabase`.
"""

from __future__ import annotations

import re

import pandas as pd
from pymatgen.core import Composition
from pymatgen.core import Element as _PMGElement

__all__ = [
    "normalize_formula",
    "calculate_dH",
    "get_formula_elements",
    "is_oxide_formula",
    "oxide_allowed_for_elements",
    "split_list",
    "is_true_stable",
    "_is_element_col",
    "clean_file_name",
    "get_formula_counts",
    "is_single_element_formula",
    "_fmt_stable",
    "_row_els",
]


def normalize_formula(formula: str) -> str:
    """Return the reduced formula string, falling back to the raw string on error.

    Args:
        formula: Any formula string recognised by pymatgen.

    Returns:
        Reduced formula (e.g. ``"Fe2O3"``), or the stripped input on failure.
    """
    try:
        return Composition(str(formula)).reduced_formula
    except Exception:
        return str(formula).strip()


def calculate_dH(formula: str, energy_per_atom: float, refs: dict) -> float:
    """Compute formation enthalpy per atom relative to elemental references.

    Args:
        formula: Chemical formula of the compound.
        energy_per_atom: Total DFT/MLIP energy divided by number of atoms (eV/atom).
        refs: Mapping of element symbol → reference energy per atom (eV/atom).

    Returns:
        Formation enthalpy in eV/atom.
    """
    comp = Composition(formula)
    total_atoms = comp.num_atoms
    ref_total = sum(amount * refs[el] for el, amount in comp.get_el_amt_dict().items() if el in refs)
    return energy_per_atom - (ref_total / total_atoms)


def get_formula_elements(formula: str) -> set[str]:
    """Return the set of element symbols present in *formula*.

    Args:
        formula: Chemical formula string.

    Returns:
        Set of element symbol strings; empty set on parse failure.
    """
    try:
        return set(Composition(str(formula)).get_el_amt_dict().keys())
    except Exception:
        return set()


def is_oxide_formula(formula: str) -> bool:
    """Return ``True`` if *formula* contains oxygen and at least one other element.

    Args:
        formula: Chemical formula string.

    Returns:
        ``True`` for multi-element oxide formulas.
    """
    elems = get_formula_elements(formula)
    return "O" in elems and len(elems) >= 2


def oxide_allowed_for_elements(formula: str, allowed_elements: set[str]) -> bool:
    """Return ``True`` if *formula* is an oxide whose non-oxygen elements are all in *allowed_elements*.

    Args:
        formula: Chemical formula string.
        allowed_elements: Set of permissible non-oxygen element symbols.

    Returns:
        ``True`` when the formula is an oxide and its cation elements are a
        subset of *allowed_elements*.
    """
    elems = get_formula_elements(formula)
    if "O" not in elems:
        return False
    return (elems - {"O"}).issubset(allowed_elements)


def split_list(text) -> list[str]:
    """Split a comma-separated string into a list of stripped tokens.

    Args:
        text: Raw cell value (may be ``NaN`` or empty).

    Returns:
        List of non-empty stripped substrings.
    """
    if pd.isna(text) or str(text).strip() == "":
        return []
    return [x.strip() for x in str(text).split(",") if x.strip()]


def is_true_stable(value) -> bool:
    """Interpret a variety of truthy representations as ``True``.

    Args:
        value: Raw stability value (bool, string, ``NaN``, …).

    Returns:
        ``True`` only when *value* unambiguously represents stability.
    """
    if pd.isna(value):
        return False
    if value is True:
        return True
    return str(value).strip().lower() in {"true", "1", "yes"}


def _is_element_col(col: str) -> bool:
    """Return ``True`` if *col* is a valid pymatgen element symbol.

    Args:
        col: Column name to test.

    Returns:
        ``True`` when pymatgen recognises *col* as an element symbol.
    """
    try:
        _PMGElement(col)
        return True
    except Exception:
        return False


def clean_file_name(name: str) -> str:
    """Sanitise *name* for use as a file-system basename.

    Replaces Windows-reserved characters and spaces; falls back to
    ``"unknown_parent"`` for blank or ``"nan"`` inputs.

    Args:
        name: Raw name string.

    Returns:
        File-system-safe string with spaces replaced by underscores.
    """
    name = str(name)
    if name.strip() == "" or name.lower() == "nan":
        name = "unknown_parent"
    name = re.sub(r'[<>:"/\\|?*]', "_", name)
    return name.replace(" ", "_")


def get_formula_counts(formula: str) -> dict[str, float]:
    """Return element → atom-count mapping for *formula*.

    Args:
        formula: Chemical formula string.

    Returns:
        Dict mapping element symbol to its count; empty dict on failure.
    """
    try:
        return Composition(str(formula)).get_el_amt_dict()
    except Exception:
        return {}


def is_single_element_formula(formula: str, element: str) -> bool:
    """Return ``True`` if *formula* contains only *element*.

    Args:
        formula: Chemical formula string.
        element: Element symbol to test for.

    Returns:
        ``True`` when the formula reduces to a single-element system matching *element*.
    """
    try:
        return set(Composition(str(formula)).get_el_amt_dict().keys()) == {element}
    except Exception:
        return False


def _fmt_stable(val) -> str:
    """Format a stability value as ``"TRUE"``, ``"FALSE"``, or ``""``.

    Args:
        val: Raw stability value.

    Returns:
        ``"TRUE"`` / ``"FALSE"`` string, or empty string for missing values.
    """
    if pd.isna(val):
        return ""
    return "TRUE" if str(val).strip().lower() in {"true", "1", "yes"} else "FALSE"


def _row_els(formula: str) -> frozenset:
    """Return a frozenset of element symbols in *formula*.

    Args:
        formula: Chemical formula string.

    Returns:
        Frozenset of element symbols; empty frozenset on failure.
    """
    try:
        return frozenset(Composition(str(formula)).get_el_amt_dict().keys())
    except Exception:
        return frozenset()
