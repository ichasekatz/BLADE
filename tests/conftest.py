"""Shared fixtures for the BLADE test suite."""

from __future__ import annotations

import pytest


@pytest.fixture
def hedb1_phases_dict() -> dict:
    """HEDB1 phase parameters for the hexagonal AlB2-prototype structure."""
    return {
        "a": 3.0624,
        "b": 3.0624,
        "c": 3.2392,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": ("0.000000 0.000000 0.000000 a\n0.333333 0.666667 0.500000 B\n0.666667 0.333333 0.500000 B"),
    }


@pytest.fixture
def hedb1_phase_entry() -> dict:
    """Single HEDB1 phase entry dict."""
    return {"generator_name": "HEDB", "lattice": "HEDB1", "supercell_size": (2, 2, 2)}


@pytest.fixture
def sqsgen_levels() -> list[dict]:
    """Two-level SQS generation spec for a binary system."""
    return [
        {"level": 0, "compositions": [[1.0, 0.0]], "letter": ["a"]},
        {"level": 1, "compositions": [[0.5, 0.5]], "letter": ["a"]},
    ]
