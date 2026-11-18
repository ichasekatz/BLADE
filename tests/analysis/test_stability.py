"""Tests for blade.analysis stability: RK fit ordering, competitor loading, oxide toggle."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from blade.analysis._stability_tables import _load_xlsx_fixed_rows, build_phase_data
from blade.analysis.blade_stability import (
    BladeStability,
    StabilityConfigError,
    StabilityThresholds,
    _run_stability_grid,
    _stability_worker,
)

# Synthetic Cr-Cu diboride family, energies in eV per formula unit (1 metal + 2 B).
E_CR, E_CU, E_B = -9.5, -4.0, -6.7  # elemental references, eV/atom
H_CRB2 = E_CR + 2 * E_B - 3 * 0.35  # ΔH_f = -0.35 eV/atom (stable diboride)
H_CUB2 = E_CU + 2 * E_B + 3 * 0.30  # ΔH_f = +0.30 eV/atom (unstable diboride)
L0, L1 = -0.20, 0.05  # RK parameters for the (Cr, Cu) pair, z = y_Cr - y_Cu
E_CRBO3_PER_ATOM = -8.0  # strongly bound fictitious borate
MU_O2_REF = -4.96  # eV/O atom, ORB O2 reference


def _h_ss(y_cr: float) -> float:
    """Exact synthetic SS enthalpy per formula unit at metal-site fraction y_Cr."""
    y_cu = 1.0 - y_cr
    return y_cr * H_CRB2 + y_cu * H_CUB2 + y_cr * y_cu * (L0 + L1 * (y_cr - y_cu))


@pytest.fixture
def system_dir(tmp_path: Path) -> Path:
    """System dir with SQS structures and a competitor xlsx."""
    sdir = tmp_path / "system_structures" / "CrCuB"
    blade = sdir / "blade"
    for n_cr in range(5):  # 4 metal sites, 8 B sites
        y_cr = n_cr / 4
        name = f"sqs_lev=0_a_Cr={y_cr},a_Cu={1 - y_cr},b_B=1.0"
        d = blade / name
        d.mkdir(parents=True)
        (d / "energy").write_text(f"{4 * _h_ss(y_cr):.10f}\n")
        (d / "CONTCAR").write_text(f"x\n1.0\n1 0 0\n0 1 0\n0 0 1\nCr Cu B\n{n_cr} {4 - n_cr} 8\n")
    comp = sdir / "parent_excel" / "MP_ORB"
    comp.mkdir(parents=True)
    pd.DataFrame(
        [
            {"formula": "Cr", "energy": E_CR, "source": "MP+ORB"},
            {"formula": "Cu", "energy": E_CU + 0.5, "source": "MP+ORB"},  # higher polymorph
            {"formula": "Cu", "energy": E_CU, "source": "MP+ORB"},
            {"formula": "B", "energy": E_B, "source": "MP+ORB"},
            {"formula": "O", "energy": MU_O2_REF, "source": "MP+ORB"},
            {"formula": "CrBO3", "energy": E_CRBO3_PER_ATOM, "source": "MP+ORB"},
            {"formula": "AlB2", "energy": -7.0, "source": "MP+ORB"},  # foreign element
            {"formula": "CrB2", "energy": H_CRB2 / 3, "source": "BLADE"},
        ]
    ).to_excel(comp / "competitors.xlsx", index=False)
    return sdir


class TestRkFitOrdering:
    """The RK fit must not depend on the metal order detect_metals returns."""

    @pytest.mark.parametrize("metals", [["Cr", "Cu"], ["Cu", "Cr"]])
    def test_fit_reproduces_sqs_energies(self, system_dir: Path, metals: list[str]) -> None:
        """For either metal order, the fitted H(y) equals the exact synthetic H(y)."""
        pdata = build_phase_data(system_dir / "blade", metals, "B", 2.0, rk_order=1, y_step=0.25)
        y_cr = pdata["y_grid"][:, metals.index("Cr")]
        expected = np.array([_h_ss(y) for y in y_cr])
        assert pdata["phase_H0"] == pytest.approx(expected, abs=1e-9)


class TestCompetitorLoading:
    """_load_xlsx_fixed_rows keeps the right rows for each oxide setting."""

    def test_oxides_excluded(self, system_dir: Path) -> None:
        """With oxides excluded: elements kept, O2/oxides/foreign/BLADE rows dropped."""
        rows = _load_xlsx_fixed_rows(system_dir / "parent_excel/MP_ORB", ["Cr", "Cu"], "B", exclude_oxides=True)
        assert sorted(r["phase_id"] for r in rows) == ["MP_B", "MP_Cr", "MP_Cu"]
        assert all(r["O"] == 0 for r in rows)

    def test_oxides_included(self, system_dir: Path) -> None:
        """With oxides included: CrBO3 is loaded with its O count; O2 itself never is."""
        rows = _load_xlsx_fixed_rows(system_dir / "parent_excel/MP_ORB", ["Cr", "Cu"], "B", exclude_oxides=False)
        by_id = {r["phase_id"]: r for r in rows}
        assert sorted(by_id) == ["MP_B", "MP_Cr", "MP_CrBO3", "MP_Cu"]
        assert by_id["MP_CrBO3"]["O"] == 3.0
        assert by_id["MP_CrBO3"]["phase_el"] == 1.0

    def test_lowest_energy_polymorph_kept(self, system_dir: Path) -> None:
        """Duplicate formulas resolve to the lowest energy, not the first row."""
        rows = _load_xlsx_fixed_rows(system_dir / "parent_excel/MP_ORB", ["Cr", "Cu"], "B")
        assert {r["phase_id"]: r["energy_eV_per_atom"] for r in rows}["MP_Cu"] == pytest.approx(E_CU)

    def test_oxide_energy_counts_oxygen_atoms(self, system_dir: Path) -> None:
        """Fixed-phase formula energy of CrBO3 includes all 5 atoms."""
        pdata = build_phase_data(
            system_dir / "blade", ["Cr", "Cu"], "B", 2.0, rk_order=1, y_step=0.25,
            competitors_dir=system_dir / "parent_excel/MP_ORB", exclude_oxides=False,
        )  # fmt: skip
        idx = list(pdata["phase_ids"]).index("MP_CrBO3")
        assert pdata["fixed_energy_formula"][idx] == pytest.approx(5 * E_CRBO3_PER_ATOM)
        assert pdata["fixed_O"][idx] == 3.0


def _grid(system_dir: Path, *, competitors: bool, exclude_oxides: bool, mu: tuple[float, ...] = ()) -> pd.DataFrame:
    pdata = build_phase_data(
        system_dir / "blade", ["Cr", "Cu"], "B", 2.0, rk_order=1, y_step=0.25,
        competitors_dir=system_dir / "parent_excel/MP_ORB" if competitors else None,
        exclude_oxides=exclude_oxides,
    )  # fmt: skip
    return _run_stability_grid(pdata, np.array([300.0]), "CrCuB", mu)


def _dg_at(df: pd.DataFrame, x_cr: float, mu: float | None = None) -> float:
    sel = np.isclose(df["x_Cr"], x_cr)
    if mu is not None:
        sel &= np.isclose(df["mu_O"], mu)
    return float(df.loc[sel, "dG_hull"].iloc[0])


class TestStabilityPhysics:
    """dG_hull responds correctly to elemental and oxide competitors."""

    def test_unstable_diboride_needs_elemental_competitors(self, system_dir: Path) -> None:
        """CuB2 (ΔH_f > 0) is only flagged unstable once Cu + 2B can compete."""
        without = _grid(system_dir, competitors=False, exclude_oxides=True)
        with_el = _grid(system_dir, competitors=True, exclude_oxides=True)
        assert _dg_at(without, 0.0) == pytest.approx(0.0, abs=1e-9)
        assert _dg_at(with_el, 0.0) == pytest.approx(3 * 0.30, abs=1e-6)
        assert _dg_at(with_el, 1.0) == pytest.approx(0.0, abs=1e-6)  # CrB2 stays stable

    def test_excluded_oxides_ignore_borate(self, system_dir: Path) -> None:
        """With oxides excluded CrB2 is stable and rows carry mu_O = NaN."""
        df = _grid(system_dir, competitors=True, exclude_oxides=True)
        assert _dg_at(df, 1.0) == pytest.approx(0.0, abs=1e-6)
        assert df["mu_O"].isna().all()
        assert not df["oxides_included"].any()

    def test_included_oxides_depend_on_mu_o(self, system_dir: Path) -> None:
        """At O2-rich μ_O the borate destabilises CrB2; at very low μ_O it cannot."""
        df = _grid(system_dir, competitors=True, exclude_oxides=False, mu=(-12.0, MU_O2_REF))
        assert df["oxides_included"].all()
        # Cr + 2B → CrBO3 + B at mu_O: Δ = (5 e_ox − 3 μ + e_B) − H_CrB2
        expected = H_CRB2 - (5 * E_CRBO3_PER_ATOM - 3 * MU_O2_REF + E_B)
        assert _dg_at(df, 1.0, MU_O2_REF) == pytest.approx(expected, abs=1e-6)
        assert _dg_at(df, 1.0, -12.0) == pytest.approx(0.0, abs=1e-6)


class TestConfigAndCache:
    """Oxide config validation and settings-aware result caching."""

    def test_including_oxides_requires_mu_o(self, tmp_path: Path) -> None:
        """exclude_oxides=False without mu_o_values is a config error."""
        with pytest.raises(StabilityConfigError, match="mu_o_values"):
            BladeStability(tmp_path, tmp_path / "out", exclude_oxides=False)

    def test_cache_invalidated_when_oxide_setting_changes(self, system_dir: Path, tmp_path: Path) -> None:
        """skip_existing reuses a parquet only if it was made with the same oxide setting."""
        out = tmp_path / "out"
        out.mkdir()

        def run(exclude: bool, mu: tuple[float, ...]) -> pd.DataFrame:
            args = (system_dir, system_dir / "blade", "parent_excel/MP_ORB", exclude, mu, "B", 2.0,
                    np.array([300.0]), 0.25, 1, out, True, StabilityThresholds())  # fmt: skip
            df, err = _stability_worker(args)
            assert err is None
            return df

        excluded = run(True, ())
        included = run(False, (MU_O2_REF,))
        assert not excluded["oxides_included"].any()
        assert included["oxides_included"].all()
        assert _dg_at(included, 1.0) > _dg_at(excluded, 1.0) + 1.0


def test_end_to_end_run(system_dir: Path, tmp_path: Path) -> None:
    """BladeStability.run (spawned worker pool) writes a combined parquet."""
    out = BladeStability(
        system_dir.parent, tmp_path / "Stability", competitors_subdir="parent_excel/MP_ORB",
        t_min=300, t_max=400, t_step=100, composition_step=0.25, rk_order=1, n_workers=1,
    ).run()  # fmt: skip
    df = pd.read_parquet(out)
    assert set(df["system"]) == {"CrCuB"}
    assert len(df) == 2 * 5
    assert (df["n_competitors"] == 3).all()


class TestParetoColumns:
    """Per-row single-phase flags and formation energies feeding the Pareto analysis."""

    @pytest.fixture
    def grid(self, system_dir: Path) -> pd.DataFrame:
        """Fine binary grid with elemental competitors at two temperatures."""
        pdata = build_phase_data(
            system_dir / "blade", ["Cr", "Cu"], "B", 2.0, rk_order=1, y_step=0.05,
            competitors_dir=system_dir / "parent_excel/MP_ORB",
        )  # fmt: skip
        return _run_stability_grid(pdata, np.array([300.0, 2000.0]), "CrCuB")

    def test_formation_energies_bracket_dg_hull(self, grid: pd.DataFrame) -> None:
        """G_ss_formation − G_formation equals dG_hull everywhere."""
        diff = grid["G_ss_formation"] - grid["G_formation"]
        assert diff.to_numpy() == pytest.approx(grid["dG_hull"].to_numpy(), abs=1e-9)

    def test_flag_hierarchy(self, grid: pd.DataFrame) -> None:
        """Unique ⇒ strict ⇒ retained, and strict rows are on the hull within tolerance."""
        assert not (grid["is_single_phase_strict"] & ~grid["is_parent_structure_retained"]).any()
        assert not (grid["is_unique_homogeneous"] & ~grid["is_single_phase_strict"]).any()
        assert (grid.loc[grid["is_single_phase_strict"], "dG_hull"] <= 1e-3).all()

    def test_endpoints_never_flagged(self, grid: pd.DataFrame) -> None:
        """Pure endpoints are not mixed alloys, so no flag is set there."""
        ends = grid[np.isclose(grid["x_Cr"], 0) | np.isclose(grid["x_Cr"], 1)]
        assert not ends[["is_parent_structure_retained", "is_single_phase_strict", "is_unique_homogeneous"]].any().any()

    def test_cu_rich_decomposes(self, grid: pd.DataFrame) -> None:
        """Near CuB2 (ΔH_f > 0) metal leaves the AlB₂ phase: ss_fraction < 1, not retained."""
        row = grid[np.isclose(grid["x_Cr"], 0.05) & (grid["T"] == 300)].iloc[0]
        assert row["ss_fraction"] < 0.95
        assert row["n_active_fixed"] >= 1
        assert not row["is_parent_structure_retained"]

    def test_cu_expelled_even_when_cr_rich(self, grid: pd.DataFrame) -> None:
        """With CuB2 0.9 eV/f.u. above Cu + 2B, Cr0.9Cu0.1 still expels Cu at 2000 K.

        0.1 × 0.9 decomposition penalty outweighs mixing enthalpy (−0.014) plus
        entropy (−0.056 eV), so ss_fraction = 0.9 (all Cr stays in the AlB₂ phase).
        """
        row = grid[np.isclose(grid["x_Cr"], 0.9) & (grid["T"] == 2000)].iloc[0]
        assert row["ss_fraction"] == pytest.approx(0.9)
        assert not row["is_parent_structure_retained"]

    def test_cr_rich_is_single_phase_without_elements(self, system_dir: Path) -> None:
        """Against SQS points only, Cr-rich alloys are one homogeneous AlB₂ phase at high T."""
        pdata = build_phase_data(system_dir / "blade", ["Cr", "Cu"], "B", 2.0, rk_order=1, y_step=0.05)
        grid = _run_stability_grid(pdata, np.array([2000.0]), "CrCuB")
        row = grid[np.isclose(grid["x_Cr"], 0.9)].iloc[0]
        assert row["ss_fraction"] == pytest.approx(1.0)
        assert row["is_single_phase_strict"]
        assert row["is_unique_homogeneous"]


def test_run_pareto_toggle(system_dir: Path, tmp_path: Path) -> None:
    """run_pareto=True writes the raw/summary database and the Pareto outputs."""
    out = tmp_path / "Stability"
    BladeStability(
        system_dir.parent, out, competitors_subdir=None,  # SQS-only hull → some strict-stable alloys
        t_min=1000, t_max=2000, t_step=100, composition_step=0.05, rk_order=1, n_workers=1,
        run_pareto=True, pareto_options={"top_n": 5, "verbose_ranking": False},
    ).run()  # fmt: skip
    raw = pd.read_parquet(out / "stability_raw.parquet")
    summary = pd.read_parquet(out / "stability_summary.parquet")
    assert {"dG_above_hull", "ss_fraction", "is_single_phase_strict", "G_ss_formation"} <= set(raw.columns)
    assert len(summary) == raw["composition"].nunique()
    metrics = pd.read_parquet(out / "pareto_analysis" / "pareto_metrics.parquet")
    assert (metrics["pareto_rank"] == 1).any()
    assert (out / "pareto_analysis" / "Diboride_Synthesizability_Ranking.xlsx").exists()
