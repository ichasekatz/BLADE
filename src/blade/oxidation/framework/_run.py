#!/usr/bin/env python3
"""Phase oxidation analysis — single entry point.

Edit the CONFIG block below, then run:
    python run.py

Works with any binary or ternary mixed-phase system.
Set DATA_ROOT to a folder with the configured structure subdirs.
The script auto-detects M1 and M2 from the directory names and prepares tables
on the first run. Subsequent runs reuse cached tables.

If tables already exist (e.g. the committed movb_*.csv / hfzrb_*.csv files),
set DATA_ROOT = None and specify M1, M2 directly, or let the script infer them
from existing CSV column names.
"""

import sys
from pathlib import Path

import numpy as np

# ============================================================
# CONFIG — edit this block
# ============================================================

# Path to raw calculation folder (must contain the configured subdirs below).
# Set to None if tables already exist and you only want to run analysis.
DATA_ROOT = None
BLADE_SUBDIR = "blade"
FIXED_PHASES_SUBDIR = "fixed_phases"

# Where processed CSV tables live / will be written
TABLES_DIR = Path(__file__).parent.parent / "tables"

# Where figures are saved
FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Element overrides — leave None to auto-detect from DATA_ROOT or table columns
# Binary:  M1="Mo", M2="V",  M3=None
# Ternary: M1="Cr", M2="Mo", M3="V"  (or pass M1=["Cr","Mo","V"])
M1 = None  # e.g. "Mo"  or list ["Cr","Mo","V"]
M2 = None  # e.g. "V"   (ignored when M1 is a list)
M3 = None  # e.g. "Cr"  (set for ternary)
PHASE_ELEMENT = None
PHASE_ELEMENT_STOICHIOMETRY = 0.0

# Redlich-Kister polynomial order (used only when preparing new tables)
RK_ORDER = 3

# ------ Which analyses to run ------
RUN_OXIDATION_SCAN = True
RUN_MUO_X_MAP = True
RUN_MUO_T_MAP = True
RUN_FIXED_PHASE_MAP = True  # control: fixed-phase composition vs flexible
RUN_TERNARY_3D_MAP = True  # 3D (composition × μO) map, ternary systems only

# If True, skip LP solving when analysis CSVs already exist — regenerate figures only.
SKIP_IF_ANALYSIS_EXISTS = True
ALLOW_CALCULATIONS = True

# ------ Per-plot on/off (within each analysis) ------
# Oxidation scan
PLOT_SCAN_NONPHASE = True  # phase fraction vs μO for fixed phases
PLOT_SCAN_PHASE_Y = True  # flexible phase composition vs μO
# muO-x / muO-T assemblage maps
PLOT_MUO_X_REGION = True  # assemblage region map (muO-x)
PLOT_MUO_T_REGION = True  # assemblage region map (muO-T)
PLOT_MUO_T_SCALARS = True  # scalar maps (phase fraction, absorbed O, etc.)
# Fixed-phase map
PLOT_FIXED_REGION = True  # assemblage region map
PLOT_FIXED_ONSET = True  # oxidation onset comparison
# Ternary 3D
PLOT_3D_EQUILIBRIUM = True  # equilibrium phase composition in (y_M1, μO, y_M2)
PLOT_3D_BOUNDARY = True  # isosurface oxidation boundary
PLOT_3D_ONSET_SURF = True  # 3D onset surface over simplex
PLOT_3D_TERNARY_DIAG = True  # 2D ternary diagram coloured by oxidation resistance

# ------ Ternary 3D map grid ------
TERNARY_3D_COMP_STEP = 0.05  # ternary simplex step (0.05 → 231 pts)
TERNARY_3D_MU_O = None  # None = reuse MAP_X_MU_O; or set e.g. np.arange(-10,-4,0.1)

# ------ Phase pseudo-grid ------
Y_STEP = 0.01  # spacing in y (M1 fraction on metal sublattice)
ACTIVE_THRESHOLD = 1e-9  # phase fractions below this treated as zero
PLOT_THRESHOLD = 1e-4  # threshold for including phase in plots
INCLUDE_0P01_TO_0P05_COMPONENTS = True

# ------ Oxidation scan (1-D vs oxygen axis) ------
SCAN_X = 0.50  # initial x_M1
SCAN_T = [500]  # K
SCAN_MU_O = np.arange(-8.8, -4.0 + 0.001, 0.03)  # eV/O
SCAN_OXYGEN_MODE = "muO"  # "muO" | "log10pO2" | "pO2"
SCAN_WRITE_ALL_PHASES = False
SCAN_MAKE_PLOTS = True
SCAN_PLOT_ONLY_ACTIVE = True

# ------ muO-x phase map ------
MAP_X_T = 700  # K — used for onset/fixed maps
MAP_X_T_VALUES = [700]  # sweep temps for muO-x + ternary-3D
ANIM_FPS = 2  # GIF/MP4 frame rate
MAP_X_X_VALUES = np.arange(0.0, 1.0 + 0.005, 0.01)  # initial x_M1 grid
MAP_X_MU_O = np.arange(-8.8, -8.0 + 0.001, 0.005)  # eV/O grid
MAP_X_MAKE_PLOTS = True
MAP_X_WRITE_WIDE = False

# ------ muO-T phase map ------
MAP_T_X = 0.50  # initial x_M1 (fixed)
MAP_T_MU_O = np.arange(-8.8, -4.0 + 0.01, 0.05)
MAP_T_T_VALUES = np.arange(500, 1501, 100)  # K grid
MAP_T_MAKE_PLOTS = True
MAP_T_WRITE_WIDE = False

# ------ Fixed-composition control map ------
# Grid reuses MAP_X_T, MAP_X_X_VALUES, MAP_X_MU_O
FIXED_OXIDATION_THRESHOLD = 1e-8  # absorbed O threshold for onset detection

# ------ Ternary composition split ------
# When scanning x_M1 in a ternary (M1-M2-M3) system, the remaining (1-x_M1)
# is split between M2 and M3 according to these ratios (must sum to 1).
# Default 0.5/0.5 = equimolar M2/M3 pseudobinary join.
TERNARY_M2_FRAC = 0.5  # fraction of (1-x_M1) assigned to M2

# ------ Plot style ------
REGION_LABEL_MODE = "id"  # "id" (integer) or "full" (assemblage string)
REGION_LABEL_FONTSIZE = 14
BOUNDARY_LW = 0.8

# ============================================================
# END CONFIG
# ============================================================

sys.path.insert(0, str(Path(__file__).parent))
from system_config import SystemConfig, prepare_tables  # noqa: E402
from thermodynamics import (  # noqa: E402
    KB,
    assign_region_ids,
    build_assemblage_labels,
    format_exact_phase_fraction_line,
    format_phase_detail_line,
    grid_edges,
    ideal_mixing_nd,
    load_phase_data,
    log10_po2_from_mu_o,
    mu_o_from_log10_po2,
    muggianu_energy_nd,
    plot_region_map,
    solve_grand_lp,
    solve_grand_lp_batch,
)

# ------------------------------------------------------------------
# System setup
# ------------------------------------------------------------------


def setup_system() -> SystemConfig:
    if isinstance(M1, list):
        metals = M1
    elif M1 and M2 and M3:
        metals = sorted([M1, M2, M3])
    elif M1 and M2:
        metals = [M1, M2]
    else:
        metals = None
    cfg = SystemConfig.resolve(
        DATA_ROOT,
        TABLES_DIR,
        metals,
        blade_subdir=BLADE_SUBDIR,
        fixed_phases_subdir=FIXED_PHASES_SUBDIR,
        phase_element=PHASE_ELEMENT,
        phase_element_stoichiometry=PHASE_ELEMENT_STOICHIOMETRY,
    )
    if not cfg.tables_ready():
        if DATA_ROOT is None or not Path(DATA_ROOT).is_dir():
            raise FileNotFoundError(
                f"Tables not found in {TABLES_DIR} and DATA_ROOT is not a valid directory.\n"
                f"Set DATA_ROOT to a folder with {BLADE_SUBDIR!r} + {FIXED_PHASES_SUBDIR!r} subdirs, "
                f"or set M1, M2 (and M3 for ternary) and point TABLES_DIR at existing CSV files."
            )
        missing = ", ".join(p.name for p in cfg.missing_table_files())
        if missing:
            print(f"[setup] Missing tables: {missing}")
        print(f"[setup] Tables not found — preparing from {DATA_ROOT} ...")
        cfg.blade_subdir = BLADE_SUBDIR
        cfg.fixed_phases_subdir = FIXED_PHASES_SUBDIR
        prepare_tables(Path(DATA_ROOT), cfg, RK_ORDER)
    return cfg


def _csv_complete(path: "Path", expected_rows: int, tol: float = 0.99) -> bool:
    """Return True if CSV exists and has at least tol × expected_rows data rows."""
    try:
        import pandas as _pd

        if not path.exists():
            return False
        df = _pd.read_csv(path)
        return len(df) >= int(expected_rows * tol)
    except Exception:
        return False


def _component_presence_threshold() -> float:
    return 0.01 if INCLUDE_0P01_TO_0P05_COMPONENTS else 0.05


def _component_cache_matches(path: "Path") -> bool:
    try:
        import pandas as _pd

        values = _pd.read_csv(path, usecols=["component_presence_threshold"])
        return not values.empty and np.allclose(
            values["component_presence_threshold"].to_numpy(dtype=float),
            _component_presence_threshold(),
        )
    except Exception:
        return False


def load_phases(cfg: SystemConfig) -> dict:
    return load_phase_data(
        cfg.phase_table_file,
        cfg.phase_grid_file,
        cfg.interaction_coeff_file,
        cfg.metals,
        phase_label=cfg.phase_label,
        y_step=Y_STEP,
        phase_element=cfg.phase_element,
        phase_element_stoichiometry=cfg.phase_element_stoichiometry,
    )


def _make_conservation_rhs(x_M1: float, n_metals: int) -> np.ndarray:
    """Build elemental conservation RHS vector for given x_M1.

    The final entry is the configured flexible-phase element stoichiometry.
    """
    if n_metals == 2:
        return np.array([x_M1, 1.0 - x_M1, PHASE_ELEMENT_STOICHIOMETRY])
    else:
        x_rest = 1.0 - x_M1
        return np.array([x_M1, TERNARY_M2_FRAC * x_rest, (1.0 - TERNARY_M2_FRAC) * x_rest, PHASE_ELEMENT_STOICHIOMETRY])


# ------------------------------------------------------------------
# Oxygen axis helpers
# ------------------------------------------------------------------


def _axis_to_mu_o(ov, T, mode):
    if mode == "muO":
        mu_o = float(ov)
        log10p = log10_po2_from_mu_o(T, mu_o)
        return mu_o, log10p, 10.0**log10p
    elif mode == "log10pO2":
        log10p = float(ov)
        mu_o = mu_o_from_log10_po2(T, log10p)
        return mu_o, log10p, 10.0**log10p
    elif mode == "pO2":
        po2 = float(ov)
        log10p = np.log10(po2)
        return mu_o_from_log10_po2(T, log10p), log10p, po2
    raise ValueError(f"Unknown SCAN_OXYGEN_MODE: {mode!r}")


# ------------------------------------------------------------------
# LP solve loop — shared by scan and map analyses
# ------------------------------------------------------------------


def _solve_row(A_eq, b_eq, grand, phase_ids, phase_kinds, phase_y, phase_O, phase_mask, phase_label, metals=None):
    """Solve LP for one state; return processed result dict."""
    n = len(phase_ids)
    amounts, omega, ok = solve_grand_lp(A_eq, b_eq, grand)
    if not ok:
        return dict(
            ok=False,
            amounts=np.full(n, np.nan),
            omega=np.nan,
            fracs=np.full(n, np.nan),
            absorbed_O=np.nan,
            phase_amt=np.nan,
            avg_y=np.nan,
            y_min=np.nan,
            y_max=np.nan,
            n_active_phase=0,
            exact_label="no feasible assemblage",
            family_label="no feasible assemblage",
        )
    total_amt = np.nansum(amounts)
    fracs = amounts / total_amt if total_amt > 0 else amounts.copy()
    fracs[np.abs(fracs) < ACTIVE_THRESHOLD] = 0.0
    absorbed_O = np.nansum(amounts * phase_O)
    bamt = np.nansum(amounts[phase_mask])
    avg_y = np.nansum(amounts[phase_mask] * phase_y[phase_mask]) / bamt if bamt > ACTIVE_THRESHOLD else np.nan
    ab = phase_mask & (fracs > ACTIVE_THRESHOLD)
    y_min = float(np.nanmin(phase_y[ab])) if np.any(ab) else np.nan
    y_max = float(np.nanmax(phase_y[ab])) if np.any(ab) else np.nan
    n_ab = int(np.sum(ab))
    metal_prefix = "".join(metals or [])
    suffix = phase_label[len(metal_prefix) :] if phase_label.startswith(metal_prefix) else ""
    el, fl = build_assemblage_labels(
        phase_ids,
        phase_kinds,
        fracs,
        ACTIVE_THRESHOLD,
        phase_label,
        metals=metals or [],
        phase_suffix=suffix,
        family_threshold=(0.01 if INCLUDE_0P01_TO_0P05_COMPONENTS else 0.05),
        family_values=amounts,
        family_threshold_inclusive=INCLUDE_0P01_TO_0P05_COMPONENTS,
    )
    return dict(
        ok=True,
        amounts=amounts,
        omega=omega,
        fracs=fracs,
        absorbed_O=absorbed_O,
        phase_amt=bamt,
        avg_y=avg_y,
        y_min=y_min,
        y_max=y_max,
        n_active_phase=n_ab,
        exact_label=el,
        family_label=fl,
    )


def _process_row(amounts, omega, ok, phase_ids, phase_kinds, phase_y, phase_O, phase_mask, phase_label, metals=None):
    """Post-LP processing — same logic as _solve_row but accepts pre-solved (amounts,omega,ok)."""
    n = len(phase_ids)
    if not ok:
        return dict(
            ok=False,
            amounts=np.full(n, np.nan),
            omega=np.nan,
            fracs=np.full(n, np.nan),
            absorbed_O=np.nan,
            phase_amt=np.nan,
            avg_y=np.nan,
            y_min=np.nan,
            y_max=np.nan,
            n_active_phase=0,
            exact_label="no feasible assemblage",
            family_label="no feasible assemblage",
        )
    total_amt = np.nansum(amounts)
    fracs = amounts / total_amt if total_amt > 0 else amounts.copy()
    fracs[np.abs(fracs) < ACTIVE_THRESHOLD] = 0.0
    absorbed_O = np.nansum(amounts * phase_O)
    bamt = np.nansum(amounts[phase_mask])
    avg_y = np.nansum(amounts[phase_mask] * phase_y[phase_mask]) / bamt if bamt > ACTIVE_THRESHOLD else np.nan
    ab = phase_mask & (fracs > ACTIVE_THRESHOLD)
    y_min = float(np.nanmin(phase_y[ab])) if np.any(ab) else np.nan
    y_max = float(np.nanmax(phase_y[ab])) if np.any(ab) else np.nan
    n_ab = int(np.sum(ab))
    metal_prefix = "".join(metals or [])
    suffix = phase_label[len(metal_prefix) :] if phase_label.startswith(metal_prefix) else ""
    el, fl = build_assemblage_labels(
        phase_ids,
        phase_kinds,
        fracs,
        ACTIVE_THRESHOLD,
        phase_label,
        metals=metals or [],
        phase_suffix=suffix,
        family_threshold=(0.01 if INCLUDE_0P01_TO_0P05_COMPONENTS else 0.05),
        family_values=amounts,
        family_threshold_inclusive=INCLUDE_0P01_TO_0P05_COMPONENTS,
    )
    return dict(
        ok=True,
        amounts=amounts,
        omega=omega,
        fracs=fracs,
        absorbed_O=absorbed_O,
        phase_amt=bamt,
        avg_y=avg_y,
        y_min=y_min,
        y_max=y_max,
        n_active_phase=n_ab,
        exact_label=el,
        family_label=fl,
    )


def _build_region_details(
    region_ids,
    exact_labels,
    wide_f,
    phase_ids,
    phase_kinds,
    metals,
    phase_suffix="",
    presence_values=None,
    presence_threshold=0.0,
    presence_threshold_inclusive=True,
):
    from collections import Counter, defaultdict

    from thermodynamics import (
        _phase_comp_range_label,
        _phase_component_signature,
        _short_label,
    )

    details = {}
    min_frac = 1e-12
    presence = np.asarray(wide_f, dtype=float) if presence_values is None else np.asarray(presence_values, dtype=float)
    fixed_cols = defaultdict(list)
    phase_cols = defaultdict(list)
    for idx, (pid, kind) in enumerate(zip(phase_ids, phase_kinds)):
        if kind == "phase":
            phase_cols[_phase_component_signature(pid, metals)].append(idx)
        else:
            fixed_cols[_short_label(pid)].append(idx)
    for rid in np.unique(region_ids):
        mask = region_ids == rid
        if not np.any(mask):
            continue
        exact_vals = [str(v) for v in np.asarray(exact_labels, dtype=object)[mask].tolist()]
        exact_label = Counter(exact_vals).most_common(1)[0][0] if exact_vals else ""
        ranges = []
        for family, cols in fixed_cols.items():
            present_vals = np.nansum(np.asarray(presence[mask][:, cols], dtype=float), axis=1)
            present_max = np.nanmax(present_vals, initial=0.0)
            if (presence_threshold_inclusive and present_max < presence_threshold - 1e-12) or (
                not presence_threshold_inclusive and present_max <= presence_threshold + 1e-12
            ):
                continue
            vals = np.nansum(np.asarray(wide_f[mask][:, cols], dtype=float), axis=1)
            vals = vals[~np.isnan(vals)]
            vals = vals[vals > min_frac]
            if len(vals) == 0:
                continue
            ranges.append((family, float(np.min(vals)), float(np.max(vals))))
        for _, cols in phase_cols.items():
            present_vals = np.nansum(np.asarray(presence[mask][:, cols], dtype=float), axis=1)
            present_max = np.nanmax(present_vals, initial=0.0)
            if (presence_threshold_inclusive and present_max < presence_threshold - 1e-12) or (
                not presence_threshold_inclusive and present_max <= presence_threshold + 1e-12
            ):
                continue
            cell_values = np.asarray(wide_f[mask][:, cols], dtype=float)
            vals = np.nansum(cell_values, axis=1)
            vals = vals[~np.isnan(vals)]
            vals = vals[vals > min_frac]
            if len(vals) == 0:
                continue
            active_local = np.any(cell_values > min_frac, axis=0)
            active_cols = [col for col, active in zip(cols, active_local) if active]
            name = _phase_comp_range_label(
                np.asarray(phase_ids, dtype=object)[active_cols],
                metals,
                phase_suffix=phase_suffix,
            )
            ranges.append((name, float(np.min(vals)), float(np.max(vals))))
        details[int(rid)] = {"exact_label": exact_label, "phase_ranges": ranges}
    return details


# ------------------------------------------------------------------
# 1. Oxidation scan
# ------------------------------------------------------------------


def run_oxidation_scan(cfg: SystemConfig, pd_data: dict) -> None:
    import pandas as pd

    M1, tag = cfg.metals[0], cfg.tag

    if SCAN_OXYGEN_MODE == "muO":
        oxygen_axis = SCAN_MU_O
        oxygen_label = r"$\mu_O^*$ (eV per O atom)"
    elif SCAN_OXYGEN_MODE == "log10pO2":
        oxygen_axis = np.arange(-220, 20.1, 2.0)
        oxygen_label = r"$\log_{10}(p_{\mathrm{O}_2}/\mathrm{bar})$"
    else:
        oxygen_axis = np.logspace(-220, 20, 121)
        oxygen_label = r"$p_{\mathrm{O}_2}$ (bar)"

    phase_ids = pd_data["phase_ids"]
    phase_kinds = pd_data["phase_kinds"]
    phase_y = pd_data["phase_y"]
    phase_metal = pd_data["phase_metal_stoich"]  # (n_phases, n_metals)
    phase_element_stoich = pd_data["phase_element_stoich"]
    phase_O = pd_data["phase_O"]
    fixed_ef = pd_data["fixed_energy_formula"]
    phase_H0 = pd_data["phase_H0"]
    phase_mix = pd_data["phase_mix_shape"]
    A_eq = pd_data["A_eq"]
    n_fixed = pd_data["n_fixed"]
    n_metals = pd_data["n_metals"]
    n_phases = len(phase_ids)
    phase_mask = phase_kinds == "phase"
    b_eq = _make_conservation_rhs(SCAN_X, n_metals)

    total = len(SCAN_T) * len(oxygen_axis)
    scan_csv = TABLES_DIR / f"{tag}_scan_summary.csv"
    if SKIP_IF_ANALYSIS_EXISTS and _csv_complete(scan_csv, total):
        long_csv = TABLES_DIR / f"{tag}_scan_phase_fractions_long.csv"
        long_df = None
        if long_csv.exists():
            try:
                candidate = pd.read_csv(long_csv)
                required = {"T_K", "kind", "phase_id", "oxygen_axis", "phase_fraction_formula"}
                if required.issubset(candidate.columns):
                    long_df = candidate
            except (OSError, ValueError):
                pass
        if not SCAN_MAKE_PLOTS or long_df is not None:
            print(f"  {total} grid points — data cached, regenerating figures")
            sum_df = pd.read_csv(scan_csv)
            if SCAN_MAKE_PLOTS:
                _plot_scan(sum_df, long_df, cfg, oxygen_label)
            return
        if not ALLOW_CALCULATIONS:
            raise RuntimeError(f"plot-only mode requires the oxidation-scan phase cache: {long_csv}")

    if not ALLOW_CALCULATIONS:
        raise RuntimeError(f"plot-only mode requires the oxidation-scan cache: {scan_csv}")

    print(f"  {total} grid points ...")
    summary_rows, long_rows = [], []
    done = 0

    for T in SCAN_T:
        phase_G = phase_H0 + KB * T * phase_mix
        for ov in oxygen_axis:
            mu_o, log10p, po2 = _axis_to_mu_o(ov, T, SCAN_OXYGEN_MODE)
            grand = np.concatenate([fixed_ef - phase_O[:n_fixed] * mu_o, phase_G])
            r = _solve_row(
                A_eq,
                b_eq,
                grand,
                phase_ids,
                phase_kinds,
                phase_y,
                phase_O,
                phase_mask,
                cfg.phase_label,
                metals=cfg.metals,
            )
            amounts, fracs = r["amounts"], r["fracs"]

            summary_rows.append(
                {
                    f"x_{M1}": SCAN_X,
                    "T_K": T,
                    "oxygen_axis": float(ov),
                    "muO_eV_per_O": mu_o,
                    "log10pO2_bar": log10p,
                    "pO2_bar": po2,
                    "Omega_eV_per_initial_formula": r["omega"],
                    "absorbed_O_atoms": r["absorbed_O"],
                    "phase_amount_formula_units": r["phase_amt"],
                    f"average_phase_y_{M1}": r["avg_y"],
                    f"phase_y_min_{M1}": r["y_min"],
                    f"phase_y_max_{M1}": r["y_max"],
                    "n_active_phase_pseudophases": r["n_active_phase"],
                    "active_phases": " + ".join(
                        phase_ids[fracs > ACTIVE_THRESHOLD].tolist() if not np.all(np.isnan(fracs)) else []
                    ),
                }
            )
            for ip in range(n_phases):
                frac = float(fracs[ip]) if not np.isnan(fracs[ip]) else 0.0
                if SCAN_WRITE_ALL_PHASES or frac > ACTIVE_THRESHOLD:
                    long_rows.append(
                        {
                            f"x_{M1}": SCAN_X,
                            "T_K": T,
                            "oxygen_axis": float(ov),
                            "muO_eV_per_O": mu_o,
                            "log10pO2_bar": log10p,
                            "pO2_bar": po2,
                            "phase_id": phase_ids[ip],
                            "kind": phase_kinds[ip],
                            f"y_{M1}": phase_y[ip],
                            "amount_formula_units": float(amounts[ip]) if not np.isnan(amounts[ip]) else np.nan,
                            "phase_fraction_formula": frac,
                            **{f"{m}_stoich": phase_metal[ip, k] for k, m in enumerate(cfg.metals)},
                            "phase_element_stoich": phase_element_stoich[ip],
                            "O_stoich": phase_O[ip],
                            "grand_eV_per_formula": float(grand[ip]),
                        }
                    )
            done += 1
            if done % max(1, total // 10) == 0 or done == total:
                print(f"    {done}/{total}")

    sum_df = pd.DataFrame(summary_rows)
    long_df = pd.DataFrame(long_rows)
    prefix = TABLES_DIR / f"{tag}_scan"
    sum_df.to_csv(str(prefix) + "_summary.csv", index=False)
    long_df.to_csv(str(prefix) + "_phase_fractions_long.csv", index=False)
    print(f"  Wrote {prefix}_summary.csv + _phase_fractions_long.csv")
    if SCAN_MAKE_PLOTS:
        _plot_scan(sum_df, long_df, cfg, oxygen_label)


def _plot_scan(sum_df, long_df, cfg: SystemConfig, oxygen_label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from thermodynamics import _phase_comp_label

    M1 = cfg.metals[0]
    sys_str = "–".join(cfg.metals + ([cfg.phase_element] if cfg.phase_element else []))
    phase_prefix = "".join(cfg.metals)
    phase_suffix = cfg.phase_label[len(phase_prefix) :] if cfg.phase_label.startswith(phase_prefix) else ""
    scan_root = FIGURES_DIR / ("".join(cfg.metals) + (cfg.phase_element or "")) / "oxidation_scan"
    scan_root.mkdir(parents=True, exist_ok=True)

    for T in SCAN_T:
        smask = np.abs(sum_df["T_K"] - T) < 1e-12
        sub = sum_df[smask].sort_values("oxygen_axis")
        oa = sub["oxygen_axis"].values
        scan_dir = scan_root / f"T{int(T)}"
        scan_dir.mkdir(parents=True, exist_ok=True)

        # Fixed-phase fractions
        lmask = (np.abs(long_df["T_K"] - T) < 1e-12) & (long_df["kind"] == "fixed")
        if PLOT_SCAN_NONPHASE and lmask.any():
            fig, ax = plt.subplots(figsize=(10, 6))
            for pid, grp in long_df[lmask].groupby("phase_id"):
                grp = grp.sort_values("oxygen_axis")
                fracs = grp["phase_fraction_formula"].values
                if SCAN_PLOT_ONLY_ACTIVE and fracs.max() < PLOT_THRESHOLD:
                    continue
                ax.plot(grp["oxygen_axis"].values, fracs, lw=1.8, label=pid.split("_")[-1])
            ax.set_xlabel(oxygen_label)
            ax.set_ylabel("fixed-phase fraction")
            ax.set_title(f"{sys_str} fixed phases, $x_{{{M1}}}$={SCAN_X:.4g}, T={T} K")
            ax.set_ylim(-0.02, 1.02)
            ax.grid(True)
            ax.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=9)
            fig.tight_layout()
            fig.savefig(scan_dir / "fixed_phase_fraction.png", dpi=150)
            plt.close(fig)

        # Flexible-phase composition
        avg_y = sub[f"average_phase_y_{M1}"].values.copy()
        y_min = sub[f"phase_y_min_{M1}"].values.copy()
        y_max = sub[f"phase_y_max_{M1}"].values.copy()
        ba = sub["phase_amount_formula_units"].values
        avg_y[ba <= PLOT_THRESHOLD] = np.nan
        y_min[ba <= PLOT_THRESHOLD] = np.nan
        y_max[ba <= PLOT_THRESHOLD] = np.nan

        phase_tmask = (np.abs(long_df["T_K"] - T) < 1e-12) & (long_df["kind"] == "phase")
        bmask = phase_tmask & (long_df["phase_fraction_formula"] > PLOT_THRESHOLD)
        if not PLOT_SCAN_PHASE_Y:
            continue
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(oa, avg_y, "k-", lw=2.2, label="avg")
        ax.plot(oa, y_min, "k--", lw=1.4, label="min")
        ax.plot(oa, y_max, "k:", lw=1.8, label="max")
        ax.set_xlabel(oxygen_label)
        ax.set_ylabel(f"$x_{{{M1}}}$ in phase")
        ax.set_title(f"{sys_str} phase composition, x={SCAN_X:.4g}, T={T} K")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        fig.savefig(scan_dir / f"average_phase_y_{M1}.png", dpi=150)
        plt.close(fig)

        fig, ax_frac = plt.subplots(figsize=(10, 6))
        if bmask.any():
            sb = long_df[phase_tmask].copy()
            sb["display_phase"] = [
                _phase_comp_label(
                    str(pid),
                    cfg.metals,
                    phase_suffix=phase_suffix,
                )
                for pid in sb["phase_id"].values
            ]
            grouped = sb.groupby(["display_phase", "oxygen_axis"], as_index=False)["phase_fraction_formula"].sum()
            for label, grp in grouped.groupby("display_phase", sort=True):
                grp = grp.sort_values("oxygen_axis")
                if grp["phase_fraction_formula"].max() <= PLOT_THRESHOLD:
                    continue
                ax_frac.plot(
                    grp["oxygen_axis"].values,
                    grp["phase_fraction_formula"].values,
                    lw=1.8,
                    label=label,
                )
        ax_frac.set_xlabel(oxygen_label)
        ax_frac.set_ylabel("flexible-phase fraction")
        ax_frac.set_title(f"{sys_str} phase fractions, x={SCAN_X:.4g}, T={T} K")
        ax_frac.set_ylim(-0.02, 1.02)
        ax_frac.grid(True)
        if bmask.any():
            ax_frac.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=9)
        fig.tight_layout()
        fig.savefig(scan_dir / "phase_fraction.png", dpi=150)
        plt.close(fig)

    print(f"  Figures → {FIGURES_DIR}")


# ------------------------------------------------------------------
# 2. muO-x phase map
# ------------------------------------------------------------------


def _plot_muO_x_cache(
    cfg,
    pd_data,
    T,
    x_values,
    mu_o_vals,
    exact_l,
    fam_l,
    wide_f,
    wide_amounts,
) -> None:
    if not MAP_X_MAKE_PLOTS:
        return

    phase_ids = pd_data["phase_ids"]
    phase_kinds = pd_data["phase_kinds"]
    region_ids, rlabels = assign_region_ids(fam_l.tolist())
    region_grid = region_ids.reshape(len(x_values), len(mu_o_vals))
    sys_str = "–".join(cfg.metals + ([cfg.phase_element] if cfg.phase_element else []))
    base = FIGURES_DIR / ("".join(cfg.metals) + (cfg.phase_element or "")) / "muO_x_phase_map"
    out_dir = base / f"T{T}" if len(MAP_X_T_VALUES) > 1 else base
    m3 = cfg.metals[2] if len(cfg.metals) > 2 else None

    region_details = _build_region_details(
        region_ids,
        exact_l,
        wide_f,
        phase_ids,
        phase_kinds,
        cfg.metals,
        pd_data.get("phase_suffix", ""),
        presence_values=wide_amounts,
        presence_threshold=(0.01 if INCLUDE_0P01_TO_0P05_COMPONENTS else 0.05),
        presence_threshold_inclusive=INCLUDE_0P01_TO_0P05_COMPONENTS,
    )
    cell_texts = np.array(
        [
            format_exact_phase_fraction_line(
                row,
                phase_ids,
                phase_kinds,
                cfg.metals,
                phase_suffix=pd_data.get("phase_suffix", ""),
            )
            for row in wide_f
        ],
        dtype=object,
    ).reshape(len(x_values), len(mu_o_vals))
    png = plot_region_map(
        mu_o_vals,
        x_values,
        region_grid,
        rlabels,
        T,
        cfg.metals[0],
        cfg.metals[1],
        out_dir,
        REGION_LABEL_MODE,
        REGION_LABEL_FONTSIZE,
        BOUNDARY_LW,
        M3=m3,
        region_details=region_details,
        cell_text_grid=cell_texts,
        system_label=sys_str,
    )
    print(f"  Figure → {png}")


def run_muO_x_map(cfg: SystemConfig, pd_data: dict) -> None:
    import pandas as pd

    M1, tag = cfg.metals[0], cfg.tag
    x_values = MAP_X_X_VALUES
    mu_o_vals = MAP_X_MU_O
    T = MAP_X_T

    phase_ids = pd_data["phase_ids"]
    phase_kinds = pd_data["phase_kinds"]
    phase_y = pd_data["phase_y"]
    phase_O = pd_data["phase_O"]
    fixed_ef = pd_data["fixed_energy_formula"]
    phase_H0 = pd_data["phase_H0"]
    phase_mix = pd_data["phase_mix_shape"]
    A_eq = pd_data["A_eq"]
    n_fixed = pd_data["n_fixed"]
    n_phases = len(phase_ids)
    phase_mask = phase_kinds == "phase"
    phase_G = phase_H0 + KB * T * phase_mix

    n_x, n_mu = len(x_values), len(mu_o_vals)
    n_states = n_x * n_mu
    cache_npz = TABLES_DIR / f"{tag}_muO_x_phase_map_T{int(T)}_arrays.npz"
    if SKIP_IF_ANALYSIS_EXISTS and cache_npz.exists():
        try:
            with np.load(cache_npz, allow_pickle=True) as cached:
                cache_matches = (
                    np.array_equal(cached["x_values"], np.asarray(x_values, dtype=float))
                    and np.array_equal(cached["mu_o_values"], np.asarray(mu_o_vals, dtype=float))
                    and float(cached["temperature"]) == float(T)
                    and np.array_equal(cached["phase_ids"].astype(str), phase_ids.astype(str))
                )
                if cache_matches:
                    exact_l = cached["exact_labels"].astype(object)
                    fam_l = cached["family_labels"].astype(object)
                    wide_f = cached["phase_fractions"]
                    wide_amounts = cached["phase_amounts"]
                else:
                    exact_l = None
            if exact_l is not None:
                print(f"  {n_states} grid points — exact data cached, skipping LP")
                _plot_muO_x_cache(
                    cfg,
                    pd_data,
                    T,
                    x_values,
                    mu_o_vals,
                    exact_l,
                    fam_l,
                    wide_f,
                    wide_amounts,
                )
                return
        except (OSError, KeyError, ValueError):
            pass

    if not ALLOW_CALCULATIONS:
        raise RuntimeError(f"plot-only mode requires the exact muO-x cache: {cache_npz}")

    print(f"  {n_states} grid points ({n_x} x * {n_mu} muO) at T={T} K ...")

    omega_g = np.full(n_states, np.nan)
    abs_O_g = np.full(n_states, np.nan)
    phase_amt_g = np.full(n_states, np.nan)
    phase_frac_g = np.full(n_states, np.nan)
    avgy_g = np.full(n_states, np.nan)
    ymin_g = np.full(n_states, np.nan)
    ymax_g = np.full(n_states, np.nan)
    n_active_phase_g = np.zeros(n_states, dtype=int)
    exact_l = np.empty(n_states, dtype=object)
    fam_l = np.empty(n_states, dtype=object)
    wide_f = np.zeros((n_states, n_phases))
    wide_amounts = np.zeros((n_states, n_phases))

    for ix, x in enumerate(x_values):
        b_eq = _make_conservation_rhs(x, pd_data["n_metals"])
        # Build every objective at once, then use deterministic cold LP solves.
        fixed_grand = fixed_ef[np.newaxis, :] - mu_o_vals[:, np.newaxis] * phase_O[:n_fixed][np.newaxis, :]
        grands = np.hstack(
            [
                fixed_grand,
                np.broadcast_to(phase_G, (len(mu_o_vals), len(phase_G))),
            ]
        )
        results = solve_grand_lp_batch(A_eq, b_eq, grands)
        for imu, (amounts, omega, ok) in enumerate(results):
            row = ix * n_mu + imu
            r = _process_row(
                amounts,
                omega,
                ok,
                phase_ids,
                phase_kinds,
                phase_y,
                phase_O,
                phase_mask,
                cfg.phase_label,
                metals=cfg.metals,
            )
            wide_f[row] = r["fracs"] if r["ok"] else 0.0
            wide_amounts[row] = r["amounts"] if r["ok"] else 0.0
            omega_g[row] = r["omega"]
            abs_O_g[row] = r["absorbed_O"]
            phase_amt_g[row] = r["phase_amt"]
            phase_frac_g[row] = np.nansum(r["fracs"][phase_mask]) if r["ok"] else np.nan
            avgy_g[row] = r["avg_y"]
            ymin_g[row] = r["y_min"]
            ymax_g[row] = r["y_max"]
            n_active_phase_g[row] = r["n_active_phase"]
            exact_l[row] = r["exact_label"]
            fam_l[row] = r["family_label"]
        if (ix + 1) % max(1, n_x // 10) == 0 or ix == n_x - 1:
            print(f"    {ix + 1}/{n_x} x done")

    state_x = np.repeat(x_values, n_mu)
    state_mu = np.tile(mu_o_vals, n_x)

    summary_df = pd.DataFrame(
        {
            f"x_{M1}": state_x,
            "T_K": T,
            "muO_eV_per_O": state_mu,
            "Omega_eV_per_initial_formula": omega_g,
            "absorbed_O_atoms": abs_O_g,
            "phase_amount_formula_units": phase_amt_g,
            "phase_fraction": phase_frac_g,
            f"average_phase_y_{M1}": avgy_g,
            f"phase_y_min_{M1}": ymin_g,
            f"phase_y_max_{M1}": ymax_g,
            "n_active_phase_pseudophases": n_active_phase_g,
            "assemblage_exact": exact_l,
            "assemblage_family": fam_l,
            "component_presence_threshold": _component_presence_threshold(),
        }
    )
    prefix = TABLES_DIR / f"{tag}_muO_x_phase_map"
    summary_df.to_csv(str(prefix) + "_summary.csv", index=False)
    np.savez_compressed(
        cache_npz,
        x_values=np.asarray(x_values, dtype=float),
        mu_o_values=np.asarray(mu_o_vals, dtype=float),
        temperature=float(T),
        phase_ids=np.asarray(phase_ids, dtype=str),
        exact_labels=np.asarray(exact_l, dtype=str),
        family_labels=np.asarray(fam_l, dtype=str),
        phase_fractions=wide_f,
        phase_amounts=wide_amounts,
    )
    print(f"  Wrote {prefix}_summary.csv")

    if MAP_X_WRITE_WIDE:
        wide_df = pd.DataFrame(wide_f, columns=phase_ids.tolist())
        wide_df.insert(0, f"x_{M1}", state_x)
        wide_df.insert(1, "T_K", T)
        wide_df.insert(2, "muO_eV_per_O", state_mu)
        wide_df.to_csv(str(prefix) + "_phase_fractions_wide.csv", index=False)

    _plot_muO_x_cache(
        cfg,
        pd_data,
        T,
        x_values,
        mu_o_vals,
        exact_l,
        fam_l,
        wide_f,
        wide_amounts,
    )


# ------------------------------------------------------------------
# 3. muO-T phase map
# ------------------------------------------------------------------


def run_muO_T_map(cfg: SystemConfig, pd_data: dict) -> None:
    import matplotlib
    import pandas as pd

    matplotlib.use("Agg")
    M1, M2, tag = cfg.metals[0], cfg.metals[1], cfg.tag
    sys_str = "–".join(cfg.metals + ([cfg.phase_element] if cfg.phase_element else []))
    mu_o_vals = MAP_T_MU_O
    T_vals = MAP_T_T_VALUES
    x = MAP_T_X

    phase_ids = pd_data["phase_ids"]
    phase_kinds = pd_data["phase_kinds"]
    phase_y = pd_data["phase_y"]
    phase_O = pd_data["phase_O"]
    fixed_ef = pd_data["fixed_energy_formula"]
    phase_H0 = pd_data["phase_H0"]
    phase_mix = pd_data["phase_mix_shape"]
    A_eq = pd_data["A_eq"]
    n_fixed = pd_data["n_fixed"]
    phase_mask = phase_kinds == "phase"
    b_eq = _make_conservation_rhs(x, pd_data["n_metals"])

    n_T, n_mu = len(T_vals), len(mu_o_vals)
    n_states = n_T * n_mu

    # Skip LP if CSV already complete — regenerate figures only
    prefix = TABLES_DIR / f"{tag}_muO_T_phase_map"
    summary_csv = Path(str(prefix) + "_summary.csv")
    if SKIP_IF_ANALYSIS_EXISTS and _csv_complete(summary_csv, n_states) and _component_cache_matches(summary_csv):
        print(f"  {n_states} grid points — data cached, regenerating figures")
        summary_df = pd.read_csv(summary_csv)
        fam_l = np.array(summary_df["assemblage_family"].tolist(), dtype=object)
        region_ids, rlabels = assign_region_ids(fam_l.tolist())
        region_grid = region_ids.reshape(n_T, n_mu)
        exact_l = np.array(summary_df["assemblage_exact"].tolist(), dtype=object)
        wide_f = np.zeros((n_states, len(phase_ids)))
        region_details = _build_region_details(
            region_ids, exact_l, wide_f, phase_ids, phase_kinds, cfg.metals, pd_data.get("phase_suffix", "")
        )
        cell_texts = np.empty((n_T, n_mu), dtype=object)
        cell_texts[:] = ""
        if MAP_T_MAKE_PLOTS:
            out_dir = FIGURES_DIR / ("".join(cfg.metals) + (cfg.phase_element or "")) / "muO_T_phase_map"
            out_dir.mkdir(parents=True, exist_ok=True)
            _plot_muO_T(
                mu_o_vals,
                T_vals,
                region_grid,
                rlabels,
                x,
                M1,
                M2,
                out_dir,
                sys_str,
                region_details=region_details,
                cell_text_grid=cell_texts,
            )
        return

    if not ALLOW_CALCULATIONS:
        raise RuntimeError(f"plot-only mode requires the muO-T cache: {summary_csv}")

    print(f"  {n_states} grid points ({n_T} T * {n_mu} muO) at x_{M1}={x} ...")

    omega_g = np.full(n_states, np.nan)
    abs_O_g = np.full(n_states, np.nan)
    phase_amt_g = np.full(n_states, np.nan)
    phase_frac_g = np.full(n_states, np.nan)
    avgy_g = np.full(n_states, np.nan)
    n_active_phase_g = np.zeros(n_states, dtype=int)
    exact_l = np.empty(n_states, dtype=object)
    fam_l = np.empty(n_states, dtype=object)
    wide_f = np.zeros((n_states, len(phase_ids)))
    wide_amounts = np.zeros((n_states, len(phase_ids)))

    for iT, T in enumerate(T_vals):
        phase_G = phase_H0 + KB * T * phase_mix
        for imu, mu_o in enumerate(mu_o_vals):
            row = iT * n_mu + imu
            grand = np.concatenate([fixed_ef - phase_O[:n_fixed] * mu_o, phase_G])
            r = _solve_row(
                A_eq,
                b_eq,
                grand,
                phase_ids,
                phase_kinds,
                phase_y,
                phase_O,
                phase_mask,
                cfg.phase_label,
                metals=cfg.metals,
            )
            omega_g[row] = r["omega"]
            abs_O_g[row] = r["absorbed_O"]
            phase_amt_g[row] = r["phase_amt"]
            phase_frac_g[row] = np.nansum(r["fracs"][phase_mask]) if r["ok"] else np.nan
            avgy_g[row] = r["avg_y"]
            n_active_phase_g[row] = r["n_active_phase"]
            exact_l[row] = r["exact_label"]
            fam_l[row] = r["family_label"]
            wide_f[row] = r["fracs"] if r["ok"] else 0.0
            wide_amounts[row] = r["amounts"] if r["ok"] else 0.0
        if (iT + 1) % max(1, n_T // 5) == 0 or iT == n_T - 1:
            print(f"    {iT + 1}/{n_T} T done")

    region_ids, rlabels = assign_region_ids(fam_l.tolist())
    region_grid = region_ids.reshape(n_T, n_mu)
    state_T = np.repeat(T_vals, n_mu)
    state_mu = np.tile(mu_o_vals, n_T)
    region_details = _build_region_details(
        region_ids,
        exact_l,
        wide_f,
        phase_ids,
        phase_kinds,
        cfg.metals,
        pd_data.get("phase_suffix", ""),
        presence_values=wide_amounts,
        presence_threshold=(0.01 if INCLUDE_0P01_TO_0P05_COMPONENTS else 0.05),
        presence_threshold_inclusive=INCLUDE_0P01_TO_0P05_COMPONENTS,
    )
    cell_texts = np.array(
        [
            format_exact_phase_fraction_line(
                row, phase_ids, phase_kinds, cfg.metals, phase_suffix=pd_data.get("phase_suffix", "")
            )
            for row in wide_f
        ],
        dtype=object,
    ).reshape(n_T, n_mu)

    summary_df = pd.DataFrame(
        {
            f"x_{M1}": x,
            "T_K": state_T,
            "muO_eV_per_O": state_mu,
            "Omega_eV_per_initial_formula": omega_g,
            "absorbed_O_atoms": abs_O_g,
            "phase_amount_formula_units": phase_amt_g,
            "phase_fraction": phase_frac_g,
            f"average_phase_y_{M1}": avgy_g,
            "n_active_phase_pseudophases": n_active_phase_g,
            "assemblage_exact": exact_l,
            "assemblage_family": fam_l,
            "component_presence_threshold": _component_presence_threshold(),
        }
    )
    prefix = TABLES_DIR / f"{tag}_muO_T_phase_map"
    summary_df.to_csv(str(prefix) + "_summary.csv", index=False)
    print(f"  Wrote {prefix}_summary.csv")

    if MAP_T_MAKE_PLOTS:
        out_dir = FIGURES_DIR / ("".join(cfg.metals) + (cfg.phase_element or "")) / "muO_T_phase_map"
        out_dir.mkdir(parents=True, exist_ok=True)
        _plot_muO_T(
            mu_o_vals,
            T_vals,
            region_grid,
            rlabels,
            x,
            M1,
            M2,
            out_dir,
            sys_str,
            region_details=region_details,
            cell_text_grid=cell_texts,
        )
        if PLOT_MUO_T_SCALARS:
            _plot_scalar_maps_T(mu_o_vals, T_vals, summary_df, M1, M2, out_dir, sys_str)


def _plot_muO_T(
    mu_o_vals, T_vals, region_grid, rlabels, x, M1, M2, out_dir, sys_str=None, region_details=None, cell_text_grid=None
):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from thermodynamics import (
        add_region_annotation,
        region_annotation_text,
        separate_region_annotations,
        write_region_map_html,
    )

    n_reg = len(rlabels)
    cmap = plt.get_cmap("tab20", max(n_reg, 1))
    colors = {r + 1: cmap(r) for r in range(n_reg)}
    mu_edges = grid_edges(mu_o_vals)
    T_edges = grid_edges(T_vals)

    if REGION_LABEL_MODE == "id":
        fig, (ax, ax_leg) = plt.subplots(
            1,
            2,
            figsize=(16, 8),
            gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05},
        )
        ax_leg.axis("off")
    else:
        fig, ax = plt.subplots(figsize=(13, 8))
        ax_leg = None

    for iT in range(len(T_vals)):
        for imu in range(len(mu_o_vals)):
            rid = region_grid[iT, imu]
            rect = plt.Rectangle(
                (mu_edges[imu], T_edges[iT]),
                mu_edges[imu + 1] - mu_edges[imu],
                T_edges[iT + 1] - T_edges[iT],
                color=colors.get(rid, (0.85, 0.85, 0.85, 1.0)),
                linewidth=0,
            )
            ax.add_patch(rect)
    for iT in range(len(T_vals)):
        for imu in range(len(mu_o_vals) - 1):
            if region_grid[iT, imu] != region_grid[iT, imu + 1]:
                ax.plot([mu_edges[imu + 1]] * 2, [T_edges[iT], T_edges[iT + 1]], "k-", lw=BOUNDARY_LW)
    for iT in range(len(T_vals) - 1):
        for imu in range(len(mu_o_vals)):
            if region_grid[iT, imu] != region_grid[iT + 1, imu]:
                ax.plot([mu_edges[imu], mu_edges[imu + 1]], [T_edges[iT + 1]] * 2, "k-", lw=BOUNDARY_LW)
    annotation_artists = []
    for rid, lbl in enumerate(rlabels, start=1):
        mask = region_grid == rid
        if not np.any(mask):
            continue
        rows, cols = np.where(mask)
        r = max(0, min(len(T_vals) - 1, int(np.round(np.median(rows)))))
        c = max(0, min(len(mu_o_vals) - 1, int(np.round(np.median(cols)))))
        annotation = region_annotation_text(rid, lbl, region_details, REGION_LABEL_MODE)
        annotation_artists.append(
            add_region_annotation(ax, mu_o_vals[c], T_vals[r], annotation, REGION_LABEL_FONTSIZE, colors[rid])
        )
    ax.set_xlim(mu_o_vals[0], mu_o_vals[-1])
    ax.set_ylim(T_vals[0], T_vals[-1])
    ax.set_xlabel(r"$\mu_O$ (eV per O atom)", fontsize=13)
    ax.set_ylabel("T (K)", fontsize=13)
    title_str = sys_str if sys_str else f"{M1}–{M2}"
    ax.set_title(f"{title_str} assemblage map, $x_{{{M1}}}$={x}", fontsize=14)
    separate_region_annotations(ax, annotation_artists)

    def _format_detail(detail):
        if not detail:
            return ""
        phase_ranges = detail.get("phase_ranges") or []
        if phase_ranges:
            return format_phase_detail_line(phase_ranges)
        exact = detail.get("exact_label")
        return exact or ""

    handles = []
    for rid, lbl in enumerate(rlabels, start=1):
        if not np.any(region_grid == rid):
            continue
        detail = _format_detail(region_details.get(rid) if region_details else None)
        label = f"{rid}: {detail if detail else lbl}"
        handles.append(mpatches.Patch(color=colors[rid], label=label))
    leg_fs = max(5, min(8, int(110 / max(len(handles), 1))))
    if ax_leg is not None:
        ax_leg.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            fontsize=leg_fs,
            framealpha=0.9,
            handlelength=1.2,
            borderaxespad=0,
            title="Region key",
            title_fontsize=leg_fs + 1,
        )

    fig.savefig(out_dir / "assemblage_region_map.png", dpi=150, bbox_inches="tight")
    try:
        # Exact hover text is prepared by callers that have pointwise fractions.
        _cell_texts = cell_text_grid
        if _cell_texts is None:
            _cell_texts = np.full(region_grid.shape, "", dtype=object)
        write_region_map_html(
            mu_o_vals,
            T_vals,
            region_grid,
            rlabels,
            out_dir / "assemblage_region_map.html",
            title=f"{title_str} assemblage map, x_{M1}={x}",
            x_label="mu_O (eV per O atom)",
            y_label="T (K)",
            region_details=region_details,
            cell_text_grid=_cell_texts,
            region_label_mode=REGION_LABEL_MODE,
            region_label_fontsize=REGION_LABEL_FONTSIZE,
        )
    except Exception as e:
        print(f"[html skipped] {e}")
    plt.close(fig)
    print(f"  Figure → {out_dir / 'assemblage_region_map.png'}")


def _plot_scalar_maps_T(mu_o_vals, T_vals, df, M1, M2, out_dir, sys_str=None):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if sys_str is None:
        sys_str = f"{M1}–{M2}"
    scalars = {
        f"average_phase_y_{M1}": f"avg phase $x_{{{M1}}}$",
        "phase_fraction": "total phase fraction",
    }
    for col, label in scalars.items():
        if col not in df.columns:
            continue
        vals = df[col].values.reshape(len(T_vals), len(mu_o_vals))
        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.pcolormesh(mu_o_vals, T_vals, vals, shading="nearest", cmap="viridis")
        fig.colorbar(im, ax=ax, label=label)
        ax.set_xlabel(r"$\mu_O$ (eV/O)", fontsize=12)
        ax.set_ylabel("T (K)", fontsize=12)
        ax.set_title(f"{sys_str}: {label}", fontsize=13)
        fig.tight_layout()
        fig.savefig(out_dir / f"{col}.png", dpi=150)
        plt.close(fig)


# ------------------------------------------------------------------
# 4. Fixed-phase control map
# ------------------------------------------------------------------


def run_fixed_phase_map(cfg: SystemConfig, pd_data: dict) -> None:
    import matplotlib
    import matplotlib.pyplot as plt
    import pandas as pd

    matplotlib.use("Agg")
    M1, M2, tag = cfg.metals[0], cfg.metals[1], cfg.tag
    sys_str = "–".join(cfg.metals + ([cfg.phase_element] if cfg.phase_element else []))
    x_values = MAP_X_X_VALUES
    mu_o_vals = MAP_X_MU_O
    T = MAP_X_T

    prefix = TABLES_DIR / f"{tag}_fixed_phase_muO_x_phase_map"
    onset_csv = TABLES_DIR / f"{tag}_fixed_phase_muO_x_phase_map_oxidation_onset.csv"
    n_x_expected = len(MAP_X_X_VALUES)

    fixed_cache_matches = False
    if SKIP_IF_ANALYSIS_EXISTS and _csv_complete(onset_csv, n_x_expected):
        try:
            cached_onset = pd.read_csv(onset_csv)
            expected_mu_step = float(mu_o_vals[1] - mu_o_vals[0]) if len(mu_o_vals) > 1 else np.nan
            fixed_cache_matches = (
                np.allclose(cached_onset[f"x_{M1}"], x_values, atol=1e-10, rtol=0.0)
                and np.allclose(cached_onset["T_K"], T, atol=1e-10, rtol=0.0)
                and np.allclose(cached_onset["muO_min"], mu_o_vals[0], atol=1e-10, rtol=0.0)
                and np.allclose(cached_onset["muO_max"], mu_o_vals[-1], atol=1e-10, rtol=0.0)
                and np.allclose(cached_onset["muO_step"], expected_mu_step, atol=1e-10, rtol=0.0, equal_nan=True)
            )
        except (KeyError, OSError, TypeError, ValueError):
            fixed_cache_matches = False
    if fixed_cache_matches:
        print("  [fixed-phase] onset CSV exists — regenerating figures only")
        _run_fixed_phase_plots(cfg, M1, M2, tag, sys_str, T, prefix, plt, pd)
        return

    if not ALLOW_CALCULATIONS:
        raise RuntimeError(f"plot-only mode requires the fixed-phase cache: {onset_csv}")

    # Load N-metal phase data for fixed-phase control calculation
    # Use load_phase_data to get all necessary arrays including Muggianu coefficients
    ctrl_data = load_phases(cfg)
    metals_local = ctrl_data["metals"]
    n_metals_local = ctrl_data["n_metals"]
    fixed_ids = ctrl_data["phase_ids"][: ctrl_data["n_fixed"]]
    fixed_metal_s = ctrl_data["phase_metal_stoich"][: ctrl_data["n_fixed"]]
    fixed_phase_element = ctrl_data["phase_element_stoich"][: ctrl_data["n_fixed"]]
    fixed_O = ctrl_data["phase_O"][: ctrl_data["n_fixed"]]
    fixed_ef = ctrl_data["fixed_energy_formula"]
    n_fixed = ctrl_data["n_fixed"]
    h_ep = ctrl_data["h_endpoints"]
    bc_coeffs = ctrl_data["binary_coeffs"]
    tc_coeff = ctrl_data["ternary_coeff"]

    n_x, n_mu = len(x_values), len(mu_o_vals)
    n_states = n_x * n_mu
    print(f"  {n_states} grid points ({n_x} x * {n_mu} muO) at T={T} K ...")

    omega_g = np.full(n_states, np.nan)
    phase_amt_g = np.full(n_states, np.nan)
    phase_frac_g = np.full(n_states, np.nan)
    exact_l = np.empty(n_states, dtype=object)
    fam_l = np.empty(n_states, dtype=object)
    onset = np.full(n_x, np.nan)

    for ix, x in enumerate(x_values):
        # Build fixed-composition phase at x_M1=x (remaining split by TERNARY_M2_FRAC)
        b_eq_full = _make_conservation_rhs(x, n_metals_local)
        # Phase composition vector: same split as b_eq_full
        y_phase = b_eq_full[:n_metals_local] / b_eq_full[:n_metals_local].sum()
        y_phase_mat = y_phase[np.newaxis, :]

        phase_H = float(muggianu_energy_nd(y_phase_mat, h_ep, bc_coeffs, tc_coeff)[0])
        phase_G = phase_H + KB * T * float(ideal_mixing_nd(y_phase_mat)[0])

        if n_metals_local == 2:
            bid = f"{cfg.phase_label}_y={y_phase[0]:.4f}"
        else:
            bid = (
                cfg.phase_label
                + "_"
                + "_".join(f"{metals_local[k]}={y_phase[k]:.4f}" for k in range(n_metals_local - 1))
            )

        pids = np.concatenate([fixed_ids, [bid]])
        pkinds = np.concatenate([np.full(n_fixed, "fixed"), ["phase"]])
        pmetal = np.vstack([fixed_metal_s, y_phase[np.newaxis, :]])
        phase_element_values = np.concatenate([fixed_phase_element, [cfg.phase_element_stoichiometry]])
        pO = np.concatenate([fixed_O, [0.0]])
        A_eq = np.vstack([pmetal.T, phase_element_values[np.newaxis, :]])
        b_eq = b_eq_full
        bm = pkinds == "phase"
        onset_found = False

        for imu, mu_o in enumerate(mu_o_vals):
            row = ix * n_mu + imu
            grand = np.concatenate([fixed_ef - pO[:n_fixed] * mu_o, [phase_G]])
            r = _solve_row(
                A_eq,
                b_eq,
                grand,
                pids,
                pkinds,
                np.concatenate([np.full(n_fixed, np.nan), [x]]),
                pO,
                bm,
                cfg.phase_label,
            )
            omega_g[row] = r["omega"]
            phase_amt_g[row] = r["phase_amt"]
            phase_frac_g[row] = np.nansum(r["fracs"][bm]) if r["ok"] else np.nan
            exact_l[row] = r["exact_label"]
            fam_l[row] = r["family_label"]
            if r["ok"] and not onset_found:
                non_b_O = np.sum(r["amounts"][~bm] * pO[~bm])
                if non_b_O > FIXED_OXIDATION_THRESHOLD:
                    onset[ix] = mu_o
                    onset_found = True
        if (ix + 1) % max(1, n_x // 10) == 0 or ix == n_x - 1:
            print(f"    {ix+1}/{n_x} x done")

    region_ids, rlabels = assign_region_ids(fam_l.tolist())
    state_x = np.repeat(x_values, n_mu)
    state_mu = np.tile(mu_o_vals, n_x)

    summary_df = pd.DataFrame(
        {
            f"x_{M1}": state_x,
            "T_K": T,
            "muO_eV_per_O": state_mu,
            "Omega_eV_per_initial_formula": omega_g,
            "phase_amount_formula_units": phase_amt_g,
            "phase_fraction": phase_frac_g,
            "assemblage_exact": exact_l,
            "assemblage_family": fam_l,
        }
    )
    mu_step = float(mu_o_vals[1] - mu_o_vals[0]) if len(mu_o_vals) > 1 else np.nan
    onset_df = pd.DataFrame(
        {
            f"x_{M1}": x_values,
            "T_K": T,
            "muO_min": float(mu_o_vals[0]),
            "muO_max": float(mu_o_vals[-1]),
            "muO_step": mu_step,
            "onset_muO_fixed_phase_eV": onset,
        }
    )
    prefix = TABLES_DIR / f"{tag}_fixed_phase_muO_x_phase_map"
    summary_df.to_csv(str(prefix) + "_summary.csv", index=False)
    onset_df.to_csv(str(prefix) + "_oxidation_onset.csv", index=False)
    print(f"  Wrote {prefix}_summary.csv + _oxidation_onset.csv")

    # Compare onset to flexible-phase model if that summary exists
    flex_onset = None
    flex_file = TABLES_DIR / f"{tag}_muO_x_phase_map_summary.csv"
    if flex_file.exists():
        flex_df = pd.read_csv(flex_file)
        flex_onset = []
        for x in x_values:
            sub = flex_df[np.abs(flex_df[f"x_{M1}"] - x) < 1e-12].sort_values("muO_eV_per_O")
            found = np.nan
            for _, row in sub.iterrows():
                if row.get("absorbed_O_atoms", 0) > FIXED_OXIDATION_THRESHOLD:
                    found = row["muO_eV_per_O"]
                    break
            flex_onset.append(found)
        pd.DataFrame(
            {
                f"x_{M1}": x_values,
                "onset_muO_fixed_phase_eV": onset,
                "onset_muO_flexible_phase_eV": flex_onset,
                "onset_shift_eV": np.array(flex_onset, dtype=float) - onset,
            }
        ).to_csv(str(prefix) + "_onset_vs_flexible.csv", index=False)

    _run_fixed_phase_plots(cfg, M1, M2, tag, sys_str, T, prefix, plt, pd)


def _run_fixed_phase_plots(cfg, M1, M2, tag, sys_str, T, prefix, plt, pd):
    """Regenerate fixed-phase figures from existing CSVs."""
    onset_csv = prefix.parent / f"{tag}_fixed_phase_muO_x_phase_map_oxidation_onset.csv"
    flex_csv = prefix.parent / f"{tag}_fixed_phase_muO_x_phase_map_onset_vs_flexible.csv"

    if not PLOT_FIXED_ONSET or not onset_csv.exists():
        return

    out_dir = FIGURES_DIR / ("".join(cfg.metals) + (cfg.phase_element or "")) / "fixed_phase_muO_x_phase_map"
    out_dir.mkdir(parents=True, exist_ok=True)

    df_on = pd.read_csv(onset_csv)
    x_col = f"x_{M1}"
    y_fixed = "onset_muO_fixed_phase_eV"
    x_vals = df_on[x_col].values
    onset = df_on[y_fixed].values

    flex_onset = None
    if flex_csv.exists():
        df_fl = pd.read_csv(flex_csv)
        if "onset_muO_flexible_phase_eV" in df_fl.columns:
            flex_onset = df_fl["onset_muO_flexible_phase_eV"].values

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x_vals, onset, "o-", ms=4, lw=1.8, label="fixed-composition phase")
    if flex_onset is not None:
        ax.plot(x_vals, np.array(flex_onset, dtype=float), "s--", ms=4, lw=1.8, label="flexible-composition phase")
    ax.set_xlabel(f"initial $x_{{{M1}}}$", fontsize=12)
    ax.set_ylabel(r"oxidation onset $\mu_O$ (eV/O)", fontsize=12)
    ax.set_title(f"Oxidation onset comparison, T={T} K", fontsize=13)
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(out_dir / "oxidation_onset_comparison.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------
# 5. Ternary 3D map  (ternary systems only)
# ------------------------------------------------------------------


def _simplex_prism_wireframe(mu_min, mu_max, color="rgba(40,40,40,0.7)", width=3):
    """Triangular prism wireframe for ternary simplex × μO.

    Simplex vertices in (X, Z) = (1,0), (0,1), (0,0). μO on Y axis.
    Returns a list of plotly Scatter3d traces.
    """
    import plotly.graph_objects as go

    verts = [(1.0, 0.0), (0.0, 1.0), (0.0, 0.0)]
    traces = []
    kw = dict(mode="lines", line=dict(color=color, width=width), showlegend=False)
    # Top and bottom triangles
    for mu in [mu_min, mu_max]:
        xs = [v[0] for v in verts] + [verts[0][0]]
        zs = [v[1] for v in verts] + [verts[0][1]]
        traces.append(go.Scatter3d(x=xs, y=[mu] * 4, z=zs, **kw))
    # Vertical edges
    for x0, z0 in verts:
        traces.append(go.Scatter3d(x=[x0, x0], y=[mu_min, mu_max], z=[z0, z0], **kw))
    return traces


def _ternary_composition_grid(comp_step):
    n_steps = int(round(1.0 / comp_step))
    comp_pts = []
    for i in range(n_steps + 1):
        x1 = i / n_steps
        for j in range(n_steps - i + 1):
            x2 = j / n_steps
            x3 = 1.0 - x1 - x2
            if x3 < -1e-10:
                break
            comp_pts.append((x1, x2, max(0.0, x3)))
    return comp_pts


def _calculate_ternary_3d(cfg, pd_data, T, mu_o_vals, comp_step):
    phase_kinds = pd_data["phase_kinds"]
    phase_O = pd_data["phase_O"]
    fixed_ef = pd_data["fixed_energy_formula"]
    phase_G = pd_data["phase_H0"] + KB * T * pd_data["phase_mix_shape"]
    A_eq = pd_data["A_eq"]
    n_fixed = pd_data["n_fixed"]
    phase_mask = phase_kinds == "phase"
    phase_y_nd = pd_data["phase_y_nd"]
    comp_pts = _ternary_composition_grid(comp_step)

    n_comp = len(comp_pts)
    n_mu = len(mu_o_vals)
    n_pts = n_comp * n_mu
    print(f"  [3D ternary] {n_pts} grid pts ({n_comp} comp × {n_mu} μO) at T={T} K ...")
    s_x1 = np.empty(n_pts)
    s_x2 = np.empty(n_pts)
    s_mu = np.empty(n_pts)
    s_bfrac = np.full(n_pts, np.nan)
    s_by1 = np.full(n_pts, np.nan)
    s_by2 = np.full(n_pts, np.nan)
    onset_mu = np.full(n_comp, np.nan)

    row = 0
    for ic, (x1, x2, x3) in enumerate(comp_pts):
        b_eq = np.array([x1, x2, x3, cfg.phase_element_stoichiometry])
        onset_found = False
        for mu_o in mu_o_vals:
            grand = np.concatenate([fixed_ef - phase_O[:n_fixed] * mu_o, phase_G])
            amounts, _, ok = solve_grand_lp(A_eq, b_eq, grand)
            s_x1[row], s_x2[row], s_mu[row] = x1, x2, mu_o
            if ok:
                total = np.nansum(amounts)
                fracs = amounts / total if total > 0 else amounts.copy()
                fracs[np.abs(fracs) < ACTIVE_THRESHOLD] = 0.0
                phase_amount = float(np.nansum(amounts[phase_mask]))
                s_bfrac[row] = float(np.nansum(fracs[phase_mask]))
                if phase_amount > ACTIVE_THRESHOLD:
                    s_by1[row] = float(np.nansum(amounts[phase_mask] * phase_y_nd[phase_mask, 0])) / phase_amount
                    s_by2[row] = float(np.nansum(amounts[phase_mask] * phase_y_nd[phase_mask, 1])) / phase_amount
                if not onset_found and np.nansum(amounts * phase_O) > FIXED_OXIDATION_THRESHOLD:
                    onset_mu[ic] = mu_o
                    onset_found = True
            row += 1
        if (ic + 1) % max(1, n_comp // 10) == 0 or ic == n_comp - 1:
            print(f"    {ic + 1}/{n_comp} compositions done")

    return comp_pts, s_x1, s_x2, s_mu, s_bfrac, s_by1, s_by2, onset_mu


def run_ternary_3d_map(cfg: SystemConfig, pd_data: dict) -> None:
    """3D plot: composition simplex × μO, coloured by dominant phase fraction.

    X = x_M1 (initial), Y = μO (eV/O), Z = x_M2 (initial).
    x_M3 = 1 - x_M1 - x_M2 enforced on the simplex.

    Generates:
      - Scatter: each LP solution coloured by total phase fraction.
      - Surface: oxidation onset μO as a function of composition.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    T = MAP_X_T  # define first — used in out_dir path
    M1, M2, M3 = cfg.metals[0], cfg.metals[1], cfg.metals[2]
    tag = cfg.tag
    sys_str = "–".join(cfg.metals + ([cfg.phase_element] if cfg.phase_element else []))
    _tern_base = FIGURES_DIR / ("".join(cfg.metals) + (cfg.phase_element or "")) / "ternary_3d_map"
    out_dir = _tern_base / f"T{T}" if len(MAP_X_T_VALUES) > 1 else _tern_base
    mu_o_vals = TERNARY_3D_MU_O if TERNARY_3D_MU_O is not None else MAP_X_MU_O
    comp_step = TERNARY_3D_COMP_STEP
    scatter_csv = TABLES_DIR / f"{tag}_ternary_3d_scatter_T{int(T)}.csv"
    onset_csv_path = TABLES_DIR / f"{tag}_ternary_3d_onset_T{int(T)}.csv"
    n_steps_exp = int(round(1.0 / comp_step))
    n_comp_expected = sum(range(n_steps_exp + 2))
    n_scatter_expected = n_comp_expected * len(mu_o_vals)
    cache_ready = False
    df_sc = None
    df_on = None
    if (
        SKIP_IF_ANALYSIS_EXISTS
        and _csv_complete(scatter_csv, n_scatter_expected)
        and _csv_complete(onset_csv_path, n_comp_expected)
    ):
        try:
            df_sc = pd.read_csv(scatter_csv)
            df_on = pd.read_csv(onset_csv_path)
            expected_comp = np.asarray(_ternary_composition_grid(comp_step), dtype=float)
            expected_x1 = np.repeat(expected_comp[:, 0], len(mu_o_vals))
            expected_x2 = np.repeat(expected_comp[:, 1], len(mu_o_vals))
            expected_mu = np.tile(np.asarray(mu_o_vals, dtype=float), len(expected_comp))
            cache_ready = (
                np.allclose(df_sc[f"x_{M1}_initial"], expected_x1, atol=1e-10, rtol=0.0)
                and np.allclose(df_sc[f"x_{M2}_initial"], expected_x2, atol=1e-10, rtol=0.0)
                and np.allclose(df_sc["muO_eV_per_O"], expected_mu, atol=1e-10, rtol=0.0)
                and np.allclose(df_on[f"x_{M1}"], expected_comp[:, 0], atol=1e-10, rtol=0.0)
                and np.allclose(df_on[f"x_{M2}"], expected_comp[:, 1], atol=1e-10, rtol=0.0)
                and np.allclose(df_on[f"x_{M3}"], expected_comp[:, 2], atol=1e-10, rtol=0.0)
            )
        except (KeyError, TypeError, ValueError):
            cache_ready = False
    if cache_ready:
        print("  [ternary 3D] exact per-temperature data cached, skipping LP")
        s_x1 = df_sc[f"x_{M1}_initial"].to_numpy(dtype=float)
        s_x2 = df_sc[f"x_{M2}_initial"].to_numpy(dtype=float)
        s_mu = df_sc["muO_eV_per_O"].to_numpy(dtype=float)
        s_bfrac = df_sc["phase_fraction"].to_numpy(dtype=float)
        s_by1 = df_sc[f"phase_y_{M1}"].to_numpy(dtype=float)
        s_by2 = df_sc[f"phase_y_{M2}"].to_numpy(dtype=float)
        comp_pts = list(zip(df_on[f"x_{M1}"], df_on[f"x_{M2}"], df_on[f"x_{M3}"]))
        onset_mu = df_on["onset_muO_eV"].to_numpy(dtype=float)
    else:
        if not ALLOW_CALCULATIONS:
            raise RuntimeError("plot-only mode requires both ternary caches: " f"{scatter_csv} and {onset_csv_path}")
        comp_pts, s_x1, s_x2, s_mu, s_bfrac, s_by1, s_by2, onset_mu = _calculate_ternary_3d(
            cfg, pd_data, T, mu_o_vals, comp_step
        )
        pd.DataFrame(
            {
                f"x_{M1}_initial": s_x1,
                f"x_{M2}_initial": s_x2,
                "muO_eV_per_O": s_mu,
                "phase_fraction": s_bfrac,
                f"phase_y_{M1}": s_by1,
                f"phase_y_{M2}": s_by2,
            }
        ).to_csv(scatter_csv, index=False)
        pd.DataFrame(
            {
                f"x_{M1}": [p[0] for p in comp_pts],
                f"x_{M2}": [p[1] for p in comp_pts],
                f"x_{M3}": [p[2] for p in comp_pts],
                "onset_muO_eV": onset_mu,
            }
        ).to_csv(onset_csv_path, index=False)

    if not any((PLOT_3D_EQUILIBRIUM, PLOT_3D_BOUNDARY, PLOT_3D_ONSET_SURF, PLOT_3D_TERNARY_DIAG)):
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    n_comp = len(comp_pts)

    import plotly.graph_objects as go

    valid = ~np.isnan(s_bfrac)
    beq_ok = valid & ~np.isnan(s_by1)  # points where phase equilibrium composition is defined
    on_valid = ~np.isnan(onset_mu)

    # ---- Plot 1: equilibrium phase composition in (y_M1, μO, y_M2) space ----
    # Each point: where does the equilibrium phase sit on the composition triangle at each μO?
    # CURVED because thermodynamics shifts the equilibrium composition continuously.
    mu_min, mu_max = float(mu_o_vals[0]), float(mu_o_vals[-1])

    pf1 = go.Figure(
        [
            go.Scatter3d(
                x=s_by1[beq_ok],
                y=s_mu[beq_ok],
                z=s_by2[beq_ok],
                mode="markers",
                marker=dict(
                    size=3,
                    color=s_bfrac[beq_ok],
                    colorscale="Viridis",
                    opacity=0.6,
                    cmin=0,
                    cmax=1,
                    colorbar=dict(title="phase<br>fraction"),
                ),
                hovertemplate=(
                    f"phase {M1}=%{{x:.3f}}<br>μO=%{{y:.3f}} eV<br>"
                    f"phase {M2}=%{{z:.3f}}<br>fraction=%{{marker.color:.3f}}<extra></extra>"
                ),
            )
        ]
        + _simplex_prism_wireframe(mu_min, mu_max)
    )
    pf1.update_layout(
        title=(
            f"{sys_str} — equilibrium phase composition vs μO, T={T} K<br>"
            f"<sub>x/y/z = phase {M1}/{M2}/{M3} fractions; colour = total phase fraction</sub>"
        ),
        scene=dict(
            xaxis=dict(title=f"phase x_{M1}", range=[0, 1]),
            yaxis=dict(title="μO (eV/O)"),
            zaxis=dict(title=f"phase x_{M2}", range=[0, 1]),
        ),
        width=950,
        height=800,
    )
    if PLOT_3D_EQUILIBRIUM:
        pf1.write_html(str(out_dir / "equilibrium_phase_composition_3d.html"))

    if PLOT_3D_EQUILIBRIUM:
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection="3d")
        sc = ax.scatter(
            s_by1[beq_ok],
            s_mu[beq_ok],
            s_by2[beq_ok],
            c=s_bfrac[beq_ok],
            cmap="viridis",
            s=5,
            alpha=0.5,
            linewidths=0,
        )
        fig.colorbar(sc, ax=ax, label="phase fraction", shrink=0.6, pad=0.1)
        ax.set_xlabel(f"phase {M1} fraction", fontsize=10, labelpad=7)
        ax.set_ylabel("μO (eV/O)", fontsize=10, labelpad=7)
        ax.set_zlabel(f"phase {M2} fraction", fontsize=10, labelpad=7)
        ax.set_title(f"{sys_str} — equilibrium phase composition\nT={T} K", fontsize=12)
        fig.tight_layout()
        fig.savefig(out_dir / "equilibrium_phase_composition_3d.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # ---- Plot 2: phase fraction isosurface in (x_M1_initial, x_M2_initial, μO) space ----
    # Curved surfaces at phase_fraction = 0.1, 0.5, 0.9 show the oxidation boundary.
    pf2 = go.Figure(_simplex_prism_wireframe(mu_min, mu_max))
    for threshold, color, name in [
        (0.9, "green", "phase fraction = 0.9 (mostly intact)"),
        (0.5, "orange", "phase fraction = 0.5 (half oxidized)"),
        (0.1, "red", "phase fraction = 0.1 (mostly oxidized)"),
    ]:
        mask = valid & (np.abs(s_bfrac - threshold) < 0.08)
        if mask.sum() < 5:
            continue
        pf2.add_trace(
            go.Scatter3d(
                x=s_x1[mask],
                y=s_mu[mask],
                z=s_x2[mask],
                mode="markers",
                marker=dict(size=3, color=color, opacity=0.5),
                name=name,
                hovertemplate=(
                    f"init x_{M1}=%{{x:.2f}}<br>μO=%{{y:.3f}} eV<br>" f"init x_{M2}=%{{z:.2f}}<extra></extra>"
                ),
            )
        )
    pf2.update_layout(
        title=(
            f"{sys_str} — oxidation boundary surfaces, T={T} K<br>"
            f"<sub>x/z = initial {M1}/{M2} fraction; y = μO; "
            f"surfaces at phase fraction = 0.1, 0.5, 0.9</sub>"
        ),
        scene=dict(
            xaxis=dict(title=f"initial x_{M1}", range=[0, 1]),
            yaxis=dict(title="μO (eV/O)"),
            zaxis=dict(title=f"initial x_{M2}", range=[0, 1]),
        ),
        width=950,
        height=800,
        legend=dict(x=0.01, y=0.99),
    )
    if PLOT_3D_BOUNDARY:
        pf2.write_html(str(out_dir / "oxidation_boundary_3d.html"))

    # ---- Plot 3: onset surface (curved surface over ternary triangle) ----
    if PLOT_3D_ONSET_SURF and on_valid.sum() >= 3:
        x1_on = np.array([comp_pts[i][0] for i in range(n_comp)])[on_valid]
        x2_on = np.array([comp_pts[i][1] for i in range(n_comp)])[on_valid]
        mu_on = onset_mu[on_valid]

        fig2 = plt.figure(figsize=(12, 9))
        ax2 = fig2.add_subplot(111, projection="3d")
        surf = ax2.plot_trisurf(x1_on, x2_on, mu_on, cmap="plasma", linewidth=0.3, edgecolor="grey", alpha=0.85)
        fig2.colorbar(surf, ax=ax2, label="onset μO (eV/O)", shrink=0.6, pad=0.1)
        ax2.set_xlabel(f"initial x_{M1}", fontsize=10, labelpad=7)
        ax2.set_ylabel(f"initial x_{M2}", fontsize=10, labelpad=7)
        ax2.set_zlabel("onset μO (eV/O)", fontsize=10, labelpad=7)
        ax2.set_title(f"{sys_str} — oxidation onset surface\nT={T} K", fontsize=12)
        fig2.tight_layout()
        fig2.savefig(out_dir / "onset_surface_3d.png", dpi=150, bbox_inches="tight")
        plt.close(fig2)

        # Base triangle outline on the simplex floor

        base_kw = dict(mode="lines", line=dict(color="rgba(40,40,40,0.6)", width=3), showlegend=False)
        mu_floor = float(mu_on.min()) - 0.05
        base_xs = [1, 0, 0, 1]
        base_ys = [0, 1, 0, 0]
        ps = go.Figure(
            [
                go.Mesh3d(
                    x=x1_on,
                    y=x2_on,
                    z=mu_on,
                    intensity=mu_on,
                    colorscale="Plasma",
                    colorbar=dict(title="onset μO<br>(eV/O)"),
                    opacity=0.9,
                    hovertemplate=(
                        f"init x_{M1}=%{{x:.3f}}<br>init x_{M2}=%{{y:.3f}}<br>" f"onset=%{{z:.3f}} eV<extra></extra>"
                    ),
                ),
                go.Scatter3d(x=base_xs, y=base_ys, z=[mu_floor] * 4, **base_kw),
            ]
        )
        ps.update_layout(
            title=f"{sys_str} — oxidation onset surface, T={T} K",
            scene=dict(
                xaxis_title=f"initial x_{M1}",
                yaxis_title=f"initial x_{M2}",
                zaxis_title="onset μO (eV/O)",
            ),
            width=950,
            height=800,
        )
        ps.write_html(str(out_dir / "onset_surface_3d.html"))

    # ---- Ternary onset diagram ----
    if PLOT_3D_TERNARY_DIAG:
        _add_ternary_diagram(
            out_dir, comp_pts, onset_mu, mu_o_vals, comp_step, M1, M2, M3, sys_str, T, tag, plt, go, np
        )

    print(f"  Figures → {out_dir}")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------


def _compile_system_animations(cfg: SystemConfig) -> None:
    """Build GIF + MP4 for each analysis type across all temperature frames."""
    import shutil
    import subprocess

    sys_name = "".join(cfg.metals) + (cfg.phase_element or "")
    sys_dir = FIGURES_DIR / sys_name
    fps = ANIM_FPS

    anim_specs = [
        ("muO_x_phase_map", "assemblage_region_map.png", "muO_x_animation"),
        ("ternary_3d_map", "ternary_onset_diagram.png", "ternary_onset_animation"),
        ("ternary_3d_map", "equilibrium_phase_composition_3d.png", "equilibrium_phase_composition_3d_animation"),
        ("ternary_3d_map", "onset_surface_3d.png", "onset_surface_3d_animation"),
    ]

    for subfolder, filename, out_stem in anim_specs:
        base = sys_dir / subfolder
        # Collect frames sorted by temperature
        frames = []
        for T in sorted(MAP_X_T_VALUES):
            p = base / f"T{T}" / filename
            if p.exists():
                frames.append(p)
        if len(frames) < 2:
            continue

        out_gif = base / f"{out_stem}.gif"
        out_mp4 = base / f"{out_stem}.mp4"

        # Build GIF with imageio or PIL — resize all frames to same size
        try:
            import imageio.v2 as imageio
            from PIL import Image as _PIL

            raw = [imageio.imread(str(f)) for f in frames]
            h0, w0 = raw[0].shape[:2]
            imgs = [
                np.array(_PIL.fromarray(img).resize((w0, h0), _PIL.LANCZOS)) if img.shape[:2] != (h0, w0) else img
                for img in raw
            ]
            imageio.mimsave(str(out_gif), imgs, fps=fps, loop=0)
        except ImportError:
            from PIL import Image

            pil_imgs = [Image.open(str(f)).convert("RGBA") for f in frames]
            w0, h0 = pil_imgs[0].size
            pil_imgs = [im.resize((w0, h0), Image.LANCZOS) for im in pil_imgs]
            pil_imgs[0].save(str(out_gif), save_all=True, append_images=pil_imgs[1:], duration=int(1000 / fps), loop=0)

        # Convert PNG frames directly to MP4 via ffmpeg if available. This is
        # much more legible than transcoding from the GIF.
        if shutil.which("ffmpeg"):
            import tempfile

            with tempfile.TemporaryDirectory() as td:
                tmp = Path(td)
                try:
                    from PIL import Image as _PIL

                    for i, f in enumerate(frames):
                        im = _PIL.open(str(f)).convert("RGB")
                        if im.size != (w0, h0):
                            im = im.resize((w0, h0), _PIL.LANCZOS)
                        im.save(tmp / f"frame_{i:04d}.png")
                    subprocess.run(
                        [
                            "ffmpeg",
                            "-y",
                            "-framerate",
                            str(fps),
                            "-i",
                            str(tmp / "frame_%04d.png"),
                            "-vf",
                            "scale=trunc(iw/2)*2:trunc(ih/2)*2",
                            "-c:v",
                            "libx264",
                            "-preset",
                            "slow",
                            "-crf",
                            "18",
                            "-pix_fmt",
                            "yuv420p",
                            str(out_mp4),
                        ],
                        capture_output=True,
                    )
                except Exception:
                    pass
        print(
            f"  {sys_name} {subfolder}: {len(frames)} frames → {out_gif.name}"
            + (f" + {out_mp4.name}" if out_mp4.exists() else " (no ffmpeg → GIF only)")
        )


def main():
    print("=== Phase oxidation analysis ===")
    cfg = setup_system()
    sys_str = "-".join(cfg.metals + ([cfg.phase_element] if cfg.phase_element else []))
    print(f"System : {sys_str}  (tag={cfg.tag})")
    print(f"Tables : {cfg.tables}")

    pd_data = load_phases(cfg)
    n_phase = len(pd_data["y_grid"])
    print(f"Phases : {pd_data['n_fixed']} fixed + {n_phase} phase pseudo-states\n")

    if RUN_OXIDATION_SCAN:
        print("--- Oxidation scan ---")
        run_oxidation_scan(cfg, pd_data)

    # Temperature-independent analyses run once at MAP_X_T
    if RUN_MUO_T_MAP:
        print(f"\n--- muO-T map (x_{cfg.metals[0]}={MAP_T_X}) ---")
        run_muO_T_map(cfg, pd_data)

    if RUN_FIXED_PHASE_MAP:
        print(f"\n--- Fixed-composition control map (T={MAP_X_T} K) ---")
        run_fixed_phase_map(cfg, pd_data)

    # Temperature sweep: muO-x map and ternary-3D run for each T in MAP_X_T_VALUES
    import sys as _sys

    _mod = _sys.modules[__name__]
    _orig_T = MAP_X_T
    for _T in MAP_X_T_VALUES:
        _mod.MAP_X_T = _T
        if RUN_MUO_X_MAP:
            print(f"\n--- muO-x map (T={_T} K) ---")
            run_muO_x_map(cfg, pd_data)
        if RUN_TERNARY_3D_MAP and cfg.n_metals == 3:
            print(f"\n--- Ternary 3D composition–μO map (T={_T} K) ---")
            run_ternary_3d_map(cfg, pd_data)
    _mod.MAP_X_T = _orig_T

    # Compile per-system animations after all temperatures
    if len(MAP_X_T_VALUES) > 1:
        _compile_system_animations(cfg)

    print("\nDone.")


if __name__ == "__main__":
    main()


# TERNARY DIAGRAM PATCH — appended
def _add_ternary_diagram(out_dir, comp_pts, onset_mu, mu_o_vals, comp_step, M1, M2, M3, sys_str, T, tag, plt, go, np):
    mu_baseline = float(mu_o_vals[0])
    x1_all = np.array([p[0] for p in comp_pts])
    x2_all = np.array([p[1] for p in comp_pts])
    x3_all = np.array([p[2] for p in comp_pts])
    resist = onset_mu - mu_baseline
    max_r = float(np.nanmax(resist)) if np.any(~np.isnan(resist)) else 6.0
    resist_filled = np.where(np.isnan(resist), max_r, resist)

    ptern = go.Figure(
        go.Scatterternary(
            a=x1_all,
            b=x2_all,
            c=x3_all,
            mode="markers",
            marker=dict(
                size=12,
                color=resist_filled,
                colorscale="RdYlGn",
                cmin=0,
                cmax=max(max_r, 0.01),
                colorbar=dict(title=f"onset − ({mu_baseline:.0f}) eV"),
                symbol="circle",
                line=dict(color="grey", width=0.5),
            ),
            hovertemplate=(
                f"x_{M1}=%{{a:.3f}}<br>x_{M2}=%{{b:.3f}}<br>x_{M3}=%{{c:.3f}}<br>"
                f"resistance=%{{marker.color:.3f}} eV<extra></extra>"
            ),
        )
    )
    ptern.update_layout(
        title=f"{sys_str} — oxidation resistance, T={T} K (green=resistant, red=early oxidation)",
        ternary=dict(aaxis=dict(title=f"x_{M1}"), baxis=dict(title=f"x_{M2}"), caxis=dict(title=f"x_{M3}")),
        width=800,
        height=700,
    )
    ptern.write_html(str(out_dir / "ternary_onset_diagram.html"))

    sqrt3_2 = np.sqrt(3) / 2
    x_cart = x2_all + 0.5 * x3_all  # M1 bottom-left, M2 bottom-right, M3 top
    y_cart = sqrt3_2 * x3_all
    fig_t, ax_t = plt.subplots(figsize=(9, 8))
    ax_t.set_aspect("equal")
    ax_t.axis("off")
    tri = np.array([[0.5, sqrt3_2], [0.0, 0.0], [1.0, 0.0], [0.5, sqrt3_2]])
    ax_t.plot(tri[:, 0], tri[:, 1], "k-", lw=1.5)
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MplPath
    from matplotlib.tri import Triangulation

    tri = Triangulation(x_cart, y_cart)
    hb = ax_t.tripcolor(tri, resist_filled, shading="flat", cmap="RdYlGn", edgecolors="black", linewidth=0.4)
    # Clip to triangle boundary
    _tri_path = MplPath([(0.0, 0.0), (1.0, 0.0), (0.5, sqrt3_2), (0.0, 0.0)])
    hb.set_clip_path(PathPatch(_tri_path, transform=ax_t.transData))
    fig_t.colorbar(
        hb, ax=ax_t, fraction=0.03, pad=0.02, label=f"onset − ({mu_baseline:.0f}) eV  (oxidation resistance)"
    )
    off = 0.04
    ax_t.set_ylim(-0.08, sqrt3_2 + 0.15)  # headroom so top vertex label clears title
    ax_t.text(-off, -off, M1, ha="right", va="top", fontsize=13, fontweight="bold")
    ax_t.text(0.5, sqrt3_2 + off, M3, ha="center", va="bottom", fontsize=13, fontweight="bold")
    ax_t.text(1 + off, -off, M2, ha="left", va="top", fontsize=13, fontweight="bold")
    # Title below triangle — avoids overlap with top vertex label
    ax_t.text(
        0.5,
        -0.06,
        f"{sys_str}  |  T={T} K  |  green=resistant  red=early oxidation",
        ha="center",
        va="top",
        transform=ax_t.transAxes,
        fontsize=10,
    )
    fig_t.tight_layout()
    fig_t.savefig(out_dir / "ternary_onset_diagram.png", dpi=180, bbox_inches="tight")
    plt.close(fig_t)
