"""Tests for BladeSQS."""

from __future__ import annotations

import shutil

import pytest

from blade.tools.blade_sqsgen import BladeSQS


class TestBladeSQSConstructor:
    """BladeSQS constructor stores parameters without side effects."""

    def test_stores_phases_dict(self, hedb1_phases_dict: dict, sqsgen_levels: list[dict]) -> None:
        """phases_dict is stored verbatim on the instance."""
        sqs = BladeSQS(
            phases_dict=hedb1_phases_dict,
            sqsgen_levels=sqsgen_levels,
            level=1,
            len_comp=2,
            skip_existing_sqs=False,
        )
        assert sqs.phases_dict is hedb1_phases_dict

    def test_stores_sqsgen_levels(self, hedb1_phases_dict: dict, sqsgen_levels: list[dict]) -> None:
        """sqsgen_levels is stored verbatim on the instance."""
        sqs = BladeSQS(
            phases_dict=hedb1_phases_dict,
            sqsgen_levels=sqsgen_levels,
            level=1,
            len_comp=2,
            skip_existing_sqs=False,
        )
        assert sqs.sqsgen_levels is sqsgen_levels

    def test_stores_level(self, hedb1_phases_dict: dict, sqsgen_levels: list[dict]) -> None:
        """Level is stored as provided."""
        sqs = BladeSQS(
            phases_dict=hedb1_phases_dict,
            sqsgen_levels=sqsgen_levels,
            level=3,
            len_comp=2,
            skip_existing_sqs=False,
        )
        assert sqs.level == 3

    def test_stores_len_comp(self, hedb1_phases_dict: dict, sqsgen_levels: list[dict]) -> None:
        """len_comp is stored as provided."""
        sqs = BladeSQS(
            phases_dict=hedb1_phases_dict,
            sqsgen_levels=sqsgen_levels,
            level=1,
            len_comp=3,
            skip_existing_sqs=False,
        )
        assert sqs.len_comp == 3

    def test_skip_existing_sqs_true(self, hedb1_phases_dict: dict, sqsgen_levels: list[dict]) -> None:
        """skip_existing_sqs=True is accepted and stored."""
        sqs = BladeSQS(
            phases_dict=hedb1_phases_dict,
            sqsgen_levels=sqsgen_levels,
            level=1,
            len_comp=2,
            skip_existing_sqs=True,
        )
        assert sqs.skip_existing_sqs is True

    def test_skip_existing_sqs_false(self, hedb1_phases_dict: dict, sqsgen_levels: list[dict]) -> None:
        """skip_existing_sqs=False is accepted and stored."""
        sqs = BladeSQS(
            phases_dict=hedb1_phases_dict,
            sqsgen_levels=sqsgen_levels,
            level=1,
            len_comp=2,
            skip_existing_sqs=False,
        )
        assert sqs.skip_existing_sqs is False

    def test_default_skip_existing_sqs_is_false(self, hedb1_phases_dict: dict, sqsgen_levels: list[dict]) -> None:
        """skip_existing_sqs defaults to False when not provided."""
        sqs = BladeSQS(
            phases_dict=hedb1_phases_dict,
            sqsgen_levels=sqsgen_levels,
            level=1,
            len_comp=2,
        )
        assert sqs.skip_existing_sqs is False


@pytest.mark.integration
def test_sqs_gen_creates_directories(
    tmp_path,
    hedb1_phases_dict: dict,
    hedb1_phase_entry: dict,
    sqsgen_levels: list[dict],
) -> None:
    """sqs_gen() creates sqsdb_lev=* directories under the blade staging area (requires ATAT on PATH)."""
    if not shutil.which("mcsqs"):
        pytest.skip("mcsqs not on PATH")

    sqsdb = tmp_path / "sqsdb"
    sqsdb.mkdir()
    paths = [tmp_path, tmp_path / "blade", sqsdb]

    sqs = BladeSQS(
        phases_dict=hedb1_phases_dict,
        sqsgen_levels=sqsgen_levels,
        level=1,
        len_comp=2,
        skip_existing_sqs=False,
    )
    params = {
        "time": 5,
        "cutoff_mode": "nn",
        "2": 3,
        "3": 2,
        "4": 2,
        "parallel_runs": 1,
        "super_cell_size": (2, 2, 2),
        "wr": 1.0,
        "wn": 1.0,
        "wd": 1.0,
    }
    sqs.sqs_gen(phase=hedb1_phase_entry, paths=paths, params=params)

    blade_sqs_dir = tmp_path / "blade" / "Files" / "SQS"
    lev_dirs = list(blade_sqs_dir.glob("**/sqsdb_lev=*"))
    assert len(lev_dirs) > 0
