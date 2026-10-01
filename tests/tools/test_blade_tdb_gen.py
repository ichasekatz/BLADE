"""Tests for BladeTDBGen."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

pytest.importorskip("materialsframework", reason="materialsframework not installed")

from blade.tools.blade_tdb_gen import BladeTDBGen  # noqa: E402


@pytest.fixture
def minimal_gen(tmp_path: Path, hedb1_phases_dict: dict, hedb1_phase_entry: dict) -> BladeTDBGen:
    """A BladeTDBGen constructed with minimal required arguments."""
    sqsdb = tmp_path / "sqsdb"
    sqsdb.mkdir()
    return BladeTDBGen(
        phases=[hedb1_phase_entry],
        liquid=False,
        paths=[tmp_path, tmp_path / "blade", sqsdb],
        composition_list=[["Cr", "Hf"]],
        level=1,
        phases_dict={"HEDB1": hedb1_phases_dict},
        skip_existing=False,
        refit_existing=False,
        output_dir=tmp_path / "Comps",
        terms_in=None,
        tdb_params={"calculator": "orb", "fmax": 1e-2, "verbose": False},
    )


class TestBladeTDBGenConstructor:
    """BladeTDBGen construction is side-effect-free."""

    def test_stores_composition_list(self, minimal_gen: BladeTDBGen) -> None:
        """composition_list is stored verbatim on the instance."""
        assert minimal_gen.composition_list == [["Cr", "Hf"]]

    def test_stores_level(self, minimal_gen: BladeTDBGen) -> None:
        """Level is stored as provided."""
        assert minimal_gen.level == 1

    def test_stores_skip_existing_false(self, minimal_gen: BladeTDBGen) -> None:
        """skip_existing=False is stored correctly."""
        assert minimal_gen.skip_existing is False

    def test_stores_liquid_false(self, minimal_gen: BladeTDBGen) -> None:
        """liquid=False is stored correctly."""
        assert minimal_gen.liquid is False

    def test_no_files_created_on_construction(self, tmp_path: Path, minimal_gen: BladeTDBGen) -> None:
        """Constructor creates no files or directories."""
        comps_dir = tmp_path / "Comps"
        assert not comps_dir.exists()

    def test_skip_existing_true_stored(self, tmp_path: Path, hedb1_phases_dict: dict, hedb1_phase_entry: dict) -> None:
        """skip_existing=True is stored correctly."""
        sqsdb = tmp_path / "sqsdb2"
        sqsdb.mkdir()
        gen = BladeTDBGen(
            phases=[hedb1_phase_entry],
            liquid=False,
            paths=[tmp_path, tmp_path / "blade", sqsdb],
            composition_list=[["Cr", "Hf"]],
            level=1,
            phases_dict={"HEDB1": hedb1_phases_dict},
            skip_existing=True,
        )
        assert gen.skip_existing is True


@pytest.mark.integration
def test_fit_produces_tdb(tmp_path: Path, hedb1_phases_dict: dict, hedb1_phase_entry: dict) -> None:
    """fit() writes a .tdb file for each composition (requires MLIP + ATAT)."""
    if not shutil.which("sqs2tdb"):
        pytest.skip("sqs2tdb not on PATH")

    sqsdb = tmp_path / "sqsdb"
    sqsdb.mkdir()
    comps_dir = tmp_path / "Comps"
    gen = BladeTDBGen(
        phases=[hedb1_phase_entry],
        liquid=False,
        paths=[tmp_path, tmp_path / "blade", sqsdb],
        composition_list=[["Cr", "Hf"]],
        level=1,
        phases_dict={"HEDB1": hedb1_phases_dict},
        skip_existing=False,
        refit_existing=False,
        output_dir=comps_dir,
        terms_in=None,
        tdb_params={"calculator": "orb", "fmax": 1e-2, "verbose": False},
    )
    gen.fit()
    assert any(comps_dir.rglob("*.tdb"))
