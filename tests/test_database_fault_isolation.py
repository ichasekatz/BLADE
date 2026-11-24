"""Verify database stage fault isolation: a single failed MP query or ORB
relaxation does not abort the full pipeline run.
"""

from __future__ import annotations

import importlib.util
import logging
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Load database_framework from examples/ without installing it as a package.
_fw_path = Path(__file__).parents[1] / "examples" / "oxidation" / "database_framework.py"
_spec = importlib.util.spec_from_file_location("blade_database_framework", _fw_path)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["blade_database_framework"] = _mod
_spec.loader.exec_module(_mod)
OxidationPipeline = _mod.OxidationPipeline
PipelineConfig = _mod.PipelineConfig


def _make_config(tmp_path: Path, **overrides) -> PipelineConfig:
    defaults = {
        "files_dir": tmp_path,
        "run_compositions": False,
        "run_download": False,
        "run_energies": False,
        "run_database": False,
    }
    defaults.update(overrides)
    return PipelineConfig(**defaults)


def test_mp_query_failure_continues(tmp_path, caplog):
    """A failed MP query for one chemsys logs a warning and skips that system."""
    cfg = _make_config(tmp_path, run_download=True)
    pipeline = OxidationPipeline(cfg)

    # Provide a minimal composition_list.xlsx with two systems.
    import pandas as pd

    comp_df = pd.DataFrame(
        [
            {"elements": "Ti, Al, C, O", "composition": "TiAlCO", "metal_composition": "TiAl", "both_composition": "TiAlC"},
            {"elements": "Nb, Al, C, O", "composition": "NbAlCO", "metal_composition": "NbAl", "both_composition": "NbAlC"},
        ]
    )
    comp_df.to_excel(tmp_path / "composition_list.xlsx", index=False)
    (tmp_path / cfg.poscar_folder).mkdir(parents=True, exist_ok=True)

    call_count = 0

    def _fake_search(**kwargs):
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            raise ConnectionError("simulated MP API failure")
        return []

    fake_mpr = MagicMock()
    fake_mpr.__enter__ = lambda s: s
    fake_mpr.__exit__ = MagicMock(return_value=False)
    fake_mpr.materials.summary.search.side_effect = _fake_search

    fake_mp_api_client = MagicMock()
    fake_mp_api_client.MPRester = lambda key: fake_mpr
    with (
        caplog.at_level(logging.WARNING, logger="blade_database_framework"),
        patch.dict(os.environ, {"MP_API_KEY": "fake-key"}),
        patch.dict(sys.modules, {
            "mp_api": MagicMock(),
            "mp_api.client": fake_mp_api_client,
            "pymatgen": MagicMock(),
            "pymatgen.core": MagicMock(),
            "pymatgen.io": MagicMock(),
            "pymatgen.io.vasp": MagicMock(),
        }),
    ):
        pipeline.download_poscars()

    assert call_count == 2, "both systems queried"
    assert any("MP query failed" in r.message for r in caplog.records), "failure logged as WARNING"


def test_orb_failure_continues(tmp_path, caplog):
    """A failed ORB relaxation for one POSCAR logs a warning and continues."""
    cfg = _make_config(tmp_path, run_energies=True, mlip="orb")
    pipeline = OxidationPipeline(cfg)

    mlip_dir = tmp_path / cfg.mlip_folder
    good = mlip_dir / "good" / "POSCAR"
    bad = mlip_dir / "bad" / "POSCAR"
    good.parent.mkdir(parents=True)
    bad.parent.mkdir(parents=True)
    good.write_text("dummy")
    bad.write_text("dummy")

    relaxed: list[str] = []

    def _fake_calculate(structure_dir, relax=True):
        name = structure_dir.name
        if name == "bad":
            raise RuntimeError("simulated ORB crash")
        relaxed.append(name)
        (structure_dir / "energy").write_text("-10.0\n")

    fake_solver = MagicMock()
    fake_solver._calculate.side_effect = _fake_calculate
    fake_calc = MagicMock()

    fake_mf_registry = MagicMock()
    fake_mf_registry.get_calculator = MagicMock(return_value=fake_calc)
    fake_mf_sqs2tdb = MagicMock()
    fake_mf_sqs2tdb.Sqs2tdb = MagicMock(return_value=fake_solver)
    with (
        caplog.at_level(logging.WARNING, logger="blade_database_framework"),
        patch.dict(sys.modules, {
            "materialsframework": MagicMock(),
            "materialsframework.calculators": MagicMock(),
            "materialsframework.calculators.registry": fake_mf_registry,
            "materialsframework.tools": MagicMock(),
            "materialsframework.tools.sqs2tdb": fake_mf_sqs2tdb,
        }),
        patch.object(pipeline, "_prepare_mlip_workspace", return_value=None),
        patch.object(pipeline, "_write_energy_per_atom", return_value=None),
    ):
        pipeline.calculate_energies()

    assert "good" in relaxed, "successful POSCAR was processed"
    assert any("ORB relaxation failed" in r.message for r in caplog.records), "failure logged as WARNING"


def test_stage_failure_raises_at_end(tmp_path):
    """A crashed stage logs and raises RuntimeError after all stages attempt."""
    cfg = _make_config(tmp_path, run_compositions=True, run_download=True)
    pipeline = OxidationPipeline(cfg)

    ran: list[str] = []

    def _fail():
        raise ValueError("stage 1 boom")

    def _ok():
        ran.append("download")

    pipeline.generate_compositions = _fail
    pipeline.download_poscars = _ok

    with pytest.raises(RuntimeError, match="stage.*failed"):
        pipeline.run()

    assert "download" in ran, "second stage ran despite first failing"
