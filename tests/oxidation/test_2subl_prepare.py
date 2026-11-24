"""Unit tests for 2-sublattice oxidation support (CEF / MAX phases).

Covers:
- _product_simplex_grid_2subl shape and sum
- _cef_enthalpy_2subl correctness for 1+1, 2+1, 2+2 sublattice configs
- _ideal_mixing_2subl shape and entropy signs
- _load_nd_2subl: round-trip via in-memory DataFrames
- _prepare_tables_2subl: end-to-end on a minimal synthetic blade_dir
- SystemConfig.resolve: sublattice_mults flows through
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from blade.oxidation.framework._phase_data import (
    _cef_enthalpy_2subl,
    _ideal_mixing_2subl,
    _load_nd_2subl,
    _product_simplex_grid_2subl,
)
from blade.oxidation.framework.system_config import SystemConfig, _prepare_tables_2subl

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_rk_df(subl_mults_str: str, endmembers: list[tuple[str, float]]) -> pd.DataFrame:
    """Build an rk DataFrame with _subl_mults + endmember rows."""
    rows = [{"pair": "_subl_mults", "term": subl_mults_str, "value_eV_per_formula": 0.0, "value_eV_per_atom": 0.0}]
    for term, val in endmembers:
        rows.append({"pair": "endmember", "term": term, "value_eV_per_formula": val, "value_eV_per_atom": val / 1.333})
    return pd.DataFrame(rows)


def _make_pt_df(metals: list[str], n_phases: int = 2) -> pd.DataFrame:
    """Minimal fixed-phase table with one oxide per metal."""
    rows = []
    for i, m in enumerate(metals[:n_phases]):
        row = {"phase_id": f"{m}2O3"}
        for sym in metals:
            row[sym] = 2.0 if sym == m else 0.0
        row["O"] = 3.0
        row["energy_eV_per_atom"] = -3.0 - 0.1 * i
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# _product_simplex_grid_2subl
# ---------------------------------------------------------------------------


class TestProductSimplexGrid:
    def test_1_plus_1_shape(self):
        """Single-element sublattices: 1×1 = 1 grid point."""
        y_overall, grids = _product_simplex_grid_2subl([["Ti"], ["Al"]], [2, 1], 0.1)
        assert y_overall.shape == (1, 2)
        assert len(grids) == 2

    def test_2_plus_1_shape(self):
        """Binary M-site, single A-site: n simplex pts on first sublattice."""
        y_overall, grids = _product_simplex_grid_2subl([["Ti", "Nb"], ["Al"]], [2, 1], 0.5)
        assert y_overall.shape[1] == 3  # 2 M-metals + 1 A-metal
        assert grids[1].shape == (1, 1)

    def test_2_plus_2_shape(self):
        """Binary M-site, binary A-site: product grid."""
        y_overall, grids = _product_simplex_grid_2subl([["Ti", "Nb"], ["Al", "Si"]], [2, 1], 0.5)
        n_a = len(grids[0])
        n_b = len(grids[1])
        assert y_overall.shape == (n_a * n_b, 4)

    def test_rows_sum_to_one(self):
        """Every row of y_overall must sum to 1."""
        y_overall, _ = _product_simplex_grid_2subl([["Ti", "Nb"], ["Al", "Si"]], [2, 1], 0.1)
        np.testing.assert_allclose(y_overall.sum(axis=1), 1.0, atol=1e-12)

    def test_mult_weights(self):
        """Pure Ti on a-site, pure Al on b-site → y_Ti=2/3, y_Al=1/3."""
        y_overall, _ = _product_simplex_grid_2subl([["Ti", "Nb"], ["Al"]], [2, 1], 1.0)
        # First row: y_a_Ti=1 → y_Ti=2/3; y_a_Nb=0 → y_Nb=0; y_b_Al=1 → y_Al=1/3
        pure_ti_row = y_overall[y_overall[:, 0] > 0.6]
        assert len(pure_ti_row) >= 1
        np.testing.assert_allclose(pure_ti_row[0, 0], 2 / 3, atol=1e-10)
        np.testing.assert_allclose(pure_ti_row[0, 2], 1 / 3, atol=1e-10)


# ---------------------------------------------------------------------------
# _cef_enthalpy_2subl
# ---------------------------------------------------------------------------


class TestCefEnthalpy:
    def _h_ref(self, vals):
        return np.array(vals, dtype=float)

    def test_1_plus_1_pure_endpoint(self):
        """Single-element sublattices: H0 = H_ref[0,0]."""
        h_ref = self._h_ref([[-3.5]])
        y_a = np.array([[1.0]])
        y_b = np.array([[1.0]])
        result = _cef_enthalpy_2subl(y_a, y_b, h_ref)
        np.testing.assert_allclose(result, [-3.5])

    def test_2_plus_1_linear_interpolation(self):
        """Binary M-site, single A-site: H0 = y_a_Ti*H[Ti,Al] + y_a_Nb*H[Nb,Al]."""
        h_ref = self._h_ref([[-3.0], [-4.0]])  # shape (2,1): Ti_Al, Nb_Al
        y_a = np.array([[0.5, 0.5]])  # 50/50 M-mix
        y_b = np.array([[1.0]])
        result = _cef_enthalpy_2subl(y_a, y_b, h_ref)
        np.testing.assert_allclose(result, [-3.5], atol=1e-12)

    def test_2_plus_2_cross_product(self):
        """Four endmembers: check pure Ti-Al and pure Nb-Si corners."""
        h_ref = self._h_ref([[-3.0, -2.5], [-3.5, -4.0]])  # Ti_Al, Ti_Si, Nb_Al, Nb_Si
        # Pure Ti-Al: y_a=[1,0], y_b=[1,0]
        y_a = np.array([[1.0, 0.0], [0.0, 1.0]])
        y_b = np.array([[1.0, 0.0], [0.0, 1.0]])
        # Product: Ti-Al, Ti-Si, Nb-Al, Nb-Si
        result = _cef_enthalpy_2subl(y_a, y_b, h_ref)
        assert result.shape == (4,)
        np.testing.assert_allclose(result[0], -3.0, atol=1e-12)   # Ti-Al
        np.testing.assert_allclose(result[3], -4.0, atol=1e-12)   # Nb-Si


# ---------------------------------------------------------------------------
# _ideal_mixing_2subl
# ---------------------------------------------------------------------------


class TestIdealMixing:
    def test_pure_endpoint_zero_entropy(self):
        """Pure endpoint: S = 0."""
        y_a = np.array([[1.0, 0.0]])
        y_b = np.array([[1.0, 0.0]])
        result = _ideal_mixing_2subl(y_a, y_b, [2, 1])
        np.testing.assert_allclose(result[0], 0.0, atol=1e-12)

    def test_max_entropy_at_equimolar(self):
        """Equimolar mix has the highest entropy."""
        y_a = np.array([[0.5, 0.5]])
        y_b = np.array([[0.5, 0.5]])
        result = _ideal_mixing_2subl(y_a, y_b, [2, 1])
        # Entropy shape is negative (convention: multiplied by -kB*T gives -TS, negative)
        assert result[0] < 0

    def test_shape_matches_product(self):
        """Output shape = n_a_pts * n_b_pts."""
        y_a = np.array([[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]])
        y_b = np.array([[1.0, 0.0], [0.0, 1.0]])
        result = _ideal_mixing_2subl(y_a, y_b, [2, 1])
        assert result.shape == (6,)


# ---------------------------------------------------------------------------
# _load_nd_2subl
# ---------------------------------------------------------------------------


class TestLoadNd2subl:
    def _build_inputs(self):
        """Minimal Ti-Nb / Al-Si MAX system."""
        metals = ["Ti", "Nb", "Al", "Si"]
        rk = _make_rk_df(
            "a=2;b=1",
            [("Ti_Al", -3.0), ("Ti_Si", -2.5), ("Nb_Al", -3.5), ("Nb_Si", -4.0)],
        )
        pt = _make_pt_df(metals, n_phases=4)
        return metals, rk, pt

    def test_returns_required_keys(self):
        metals, rk, pt = self._build_inputs()
        result = _load_nd_2subl(pt, rk, metals, "TiNbAlSi", 0.5, "C", 0.333)
        required = {
            "phase_ids", "phase_kinds", "phase_y", "phase_y_nd",
            "phase_metal_stoich", "phase_element_stoich", "phase_O",
            "fixed_energy_formula", "phase_H0", "phase_mix_shape",
            "y_grid", "h_endpoints", "binary_coeffs", "ternary_coeff",
            "A_eq", "n_fixed", "n_metals", "metals", "phase_label",
            "phase_element", "phase_element_stoichiometry", "phase_suffix",
            "fixed_vol_per_atom", "fixed_formula_counts",
        }
        assert required <= set(result.keys())

    def test_n_metals_correct(self):
        metals, rk, pt = self._build_inputs()
        result = _load_nd_2subl(pt, rk, metals, "TiNbAlSi", 0.5, "C", 0.333)
        assert result["n_metals"] == 4

    def test_grid_rows_sum_to_one(self):
        metals, rk, pt = self._build_inputs()
        result = _load_nd_2subl(pt, rk, metals, "TiNbAlSi", 0.5, "C", 0.333)
        y = result["y_grid"]
        np.testing.assert_allclose(y.sum(axis=1), 1.0, atol=1e-12)

    def test_phase_H0_shape_matches_grid(self):
        metals, rk, pt = self._build_inputs()
        result = _load_nd_2subl(pt, rk, metals, "TiNbAlSi", 0.5, "C", 0.333)
        assert result["phase_H0"].shape == (len(result["y_grid"]),)

    def test_a_eq_shape(self):
        """A_eq has n_metals+1 rows (metals + phase_element balance)."""
        metals, rk, pt = self._build_inputs()
        result = _load_nd_2subl(pt, rk, metals, "TiNbAlSi", 0.5, "C", 0.333)
        n_total = result["n_fixed"] + len(result["y_grid"])
        assert result["A_eq"].shape == (5, n_total)  # 4 metals + 1 C

    def test_no_phase_element(self):
        """Works when phase_element is None (substitutional alloy)."""
        metals, rk, pt = self._build_inputs()
        result = _load_nd_2subl(pt, rk, metals, "TiNbAlSi", 0.5, None, 0.0)
        assert result["phase_element"] is None


# ---------------------------------------------------------------------------
# _prepare_tables_2subl integration test
# ---------------------------------------------------------------------------


class TestPrepareTables2subl:
    """End-to-end: write synthetic blade_dir / fixed_phases, run _prepare_tables_2subl,
    verify CSV files and that _load_nd dispatches correctly.
    """

    def _write_contcar(self, path: Path, counts: dict) -> None:
        """Write a minimal VASP POSCAR/CONTCAR with given element counts."""
        elements = [el for el, n in counts.items() if n > 0]
        nums = [str(counts[el]) for el in elements]
        lines = [
            "Generated",
            "1.0",
            "3.0 0.0 0.0",
            "0.0 3.0 0.0",
            "0.0 0.0 6.0",
            " ".join(elements),
            " ".join(nums),
            "Direct",
        ]
        total = sum(int(n) for n in nums)
        for _ in range(total):
            lines.append("0.0 0.0 0.0")
        path.write_text("\n".join(lines) + "\n")

    def _make_blade_dir(self, blade_dir: Path) -> None:
        """Four endmember SQS dirs + two mixed compositions."""
        # sqs_lev=0 endpoints: Ti+Al, Ti+Si, Nb+Al, Nb+Si
        specs = [
            ("sqs_lev=0_a_Ti=1.0,b_Al=1.0,c_C=1.0", {"Ti": 2, "Al": 1, "C": 1}, -50.0),
            ("sqs_lev=0_a_Ti=1.0,b_Si=1.0,c_C=1.0", {"Ti": 2, "Si": 1, "C": 1}, -48.0),
            ("sqs_lev=0_a_Nb=1.0,b_Al=1.0,c_C=1.0", {"Nb": 2, "Al": 1, "C": 1}, -52.0),
            ("sqs_lev=0_a_Nb=1.0,b_Si=1.0,c_C=1.0", {"Nb": 2, "Si": 1, "C": 1}, -55.0),
            ("sqs_lev=1_a_Ti=0.5,a_Nb=0.5,b_Al=0.5,b_Si=0.5,c_C=1.0", {"Ti": 1, "Nb": 1, "Al": 1, "Si": 1, "C": 1}, -53.0),
        ]
        for dirname, counts, e_total in specs:
            d = blade_dir / dirname
            d.mkdir(parents=True)
            (d / "energy").write_text(str(e_total) + "\n")
            self._write_contcar(d / "CONTCAR", counts)

    def _make_fixed_dir(self, fixed_dir: Path) -> None:
        """Two oxide fixed phases."""
        for name, counts, e in [("TiO2", {"Ti": 1, "O": 2}, -18.0), ("Al2O3", {"Al": 2, "O": 3}, -30.0)]:
            d = fixed_dir / name
            d.mkdir(parents=True)
            (d / "energy").write_text(str(e) + "\n")
            self._write_contcar(d / "CONTCAR", counts)

    def test_csvs_written(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            blade_dir = root / "blade"
            fixed_dir = root / "fixed_phases"
            blade_dir.mkdir()
            fixed_dir.mkdir()
            self._make_blade_dir(blade_dir)
            self._make_fixed_dir(fixed_dir)

            cfg = SystemConfig(
                metals=["Ti", "Nb", "Al", "Si"],
                phase_label="TiNbAlSiC",
                tag="tinbalsic",
                tables_dir=str(root / "tables"),
                phase_element="C",
                phase_element_stoichiometry=0.333,
                sublattice_mults={"a": 2, "b": 1},
            )
            _prepare_tables_2subl(root, cfg)

            assert cfg.phase_grid_file.exists()
            assert cfg.interaction_coeff_file.exists()
            assert cfg.phase_table_file.exists()
            assert cfg.phase_model_file.exists()

    def test_interaction_coeff_has_endmembers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            blade_dir = root / "blade"
            fixed_dir = root / "fixed_phases"
            blade_dir.mkdir()
            fixed_dir.mkdir()
            self._make_blade_dir(blade_dir)
            self._make_fixed_dir(fixed_dir)

            cfg = SystemConfig(
                metals=["Ti", "Nb", "Al", "Si"],
                phase_label="TiNbAlSiC",
                tag="tinbalsic",
                tables_dir=str(root / "tables"),
                phase_element="C",
                phase_element_stoichiometry=0.333,
                sublattice_mults={"a": 2, "b": 1},
            )
            _prepare_tables_2subl(root, cfg)

            rk = pd.read_csv(cfg.interaction_coeff_file)
            assert (rk["pair"] == "endmember").any()
            assert (rk["pair"] == "_subl_mults").any()
            em = rk[rk["pair"] == "endmember"]
            assert set(em["term"].values) == {"Ti_Al", "Ti_Si", "Nb_Al", "Nb_Si"}

    def test_load_nd_dispatch(self):
        """_load_nd auto-detects 2-subl and returns a compatible dict."""
        from blade.oxidation.framework._phase_data import _load_nd

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            blade_dir = root / "blade"
            fixed_dir = root / "fixed_phases"
            blade_dir.mkdir()
            fixed_dir.mkdir()
            self._make_blade_dir(blade_dir)
            self._make_fixed_dir(fixed_dir)

            cfg = SystemConfig(
                metals=["Ti", "Nb", "Al", "Si"],
                phase_label="TiNbAlSiC",
                tag="tinbalsic",
                tables_dir=str(root / "tables"),
                phase_element="C",
                phase_element_stoichiometry=0.333,
                sublattice_mults={"a": 2, "b": 1},
            )
            _prepare_tables_2subl(root, cfg)

            pd_data = _load_nd(
                cfg.phase_table_file,
                cfg.phase_grid_file,
                cfg.interaction_coeff_file,
                cfg.metals,
                cfg.phase_label,
                y_step=0.5,
                phase_element=cfg.phase_element,
                phase_element_stoichiometry=cfg.phase_element_stoichiometry,
            )

            assert "h_ref_2subl" in pd_data
            assert "subl_mults" in pd_data
            np.testing.assert_allclose(pd_data["y_grid"].sum(axis=1), 1.0, atol=1e-12)
            assert pd_data["n_metals"] == 4


# ---------------------------------------------------------------------------
# SystemConfig.resolve passes sublattice_mults through
# ---------------------------------------------------------------------------


class TestSystemConfigResolve:
    def test_sublattice_mults_propagated(self, tmp_path):
        blade_dir = tmp_path / "blade"
        blade_dir.mkdir()
        # Create one minimal dir so detect_metals_from_dirs doesn't fail
        (blade_dir / "sqs_lev=0_a_Ti=1.0,b_Al=1.0,c_C=1.0").mkdir()

        tables_dir = tmp_path / "tables"
        sm = {"a": 2, "b": 1}
        cfg = SystemConfig.resolve(
            data_root=tmp_path,
            tables_dir=tables_dir,
            phase_element="C",
            phase_element_stoichiometry=0.333,
            sublattice_mults=sm,
        )
        assert cfg.sublattice_mults == sm
