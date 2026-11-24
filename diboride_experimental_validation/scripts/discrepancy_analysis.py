"""Why do the calculations favour mixing while Ness et al. (2026) observe lamellar separation?

Reads (never writes) the BLADE stability results, the SQS structures, and the processed /
digitized experimental and DFT data; writes into simulation_comparison/discrepancy/ and figures/.

Sections (matching discrepancy_analysis.md):
  1  x = 0.5 equilibrium at every temperature: ΔG_hull, phase fractions, AlB2 phases, flags, Hessian
  2  audit of the Fig. 2 white regions
  3  SQS / reference-state / normalization / RK-fit audit
  4  common tangent, binodal and spinodal check (AlB2-only model)
  5  eutectic vs solid-state separation: required solid mixing enthalpies
  6  ZrB2-NbB2 / ZrB2-TaB2 vs Otani (and Vorotilo) at matched conditions
  7  DFT vs ORB vs Ness thermodynamic parameters for ZrB2-TaB2
  8  six-system validation table

Usage (BLADE repository root):
    .venv/bin/python diboride_experimental_validation/scripts/discrepancy_analysis.py
"""

from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
VAL, REPO = HERE.parent, HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

from blade.analysis._stability_tables import _scan_sqs_dirs, build_phase_data, detect_metals  # noqa: E402
from blade.oxidation.framework._lp import solve_grand_lp  # noqa: E402
from compare_simulation import (  # noqa: E402
    COMPETITORS, KB, ORDER, RK_ORDER, SIM_DIR, STRUCT_DIR, SYSTEMS, X_EXP, intervals, lower_hull_gap, model_gss,
)

OUT = VAL / "simulation_comparison" / "discrepancy"
FIG = VAL / "figures"
R = 8.314462618e-3          # kJ/(mol K)
EV = 96.48533212            # kJ/mol per eV
ATOMS_PER_FU = 3            # MB2


def per_atom(e_fu):
    return e_fu / ATOMS_PER_FU


# ── 1. x = 0.5 equilibrium at every temperature ───────────────────────────────
def section1() -> pd.DataFrame:
    rows = []
    for sid in ORDER:
        tag, m4, m5 = SYSTEMS[sid]
        sim = pd.read_parquet(SIM_DIR / f"{tag}_stability.parquet")
        xcol = f"x_{m5}"
        sdir = STRUCT_DIR / tag
        metals = detect_metals(sdir / "blade", "B")
        pdat = build_phase_data(sdir / "blade", metals, "B", 2.0, rk_order=RK_ORDER, y_step=0.05,
                                competitors_dir=[sdir / c for c in COMPETITORS if (sdir / c).is_dir()])
        nf = pdat["n_fixed"]
        metal_content = pdat["A_eq"][: len(metals), :].sum(axis=0)
        y = np.array([X_EXP if m == m5 else 1 - X_EXP for m in metals])
        for t in np.sort(sim["T"].unique()):
            st = sim[sim["T"] == t].set_index(sim[sim["T"] == t][xcol].round(6))
            r = st.loc[round(X_EXP, 6)]
            d2 = (st.loc[0.55, "G_ss"] - 2 * r["G_ss"] + st.loc[0.45, "G_ss"]) / 0.05**2
            grand = np.concatenate([pdat["fixed_energy_formula"], pdat["phase_H0"] + KB * t * pdat["phase_mix_shape"]])
            amt, obj, ok = solve_grand_lp(pdat["A_eq"], np.concatenate([y, [2.0]]), grand)
            phases = []
            for k in np.flatnonzero(amt > 1e-8):
                frac = amt[k] * metal_content[k]  # fraction of metal atoms (sums to 1)
                if k >= nf:
                    phases.append(("AlB2-SS", float(pdat["y_grid"][k - nf][metals.index(m5)]), frac))
                elif pdat["fixed_is_parent"][k]:
                    phases.append(("AlB2-SQS:" + str(pdat["phase_ids"][k]), float(pdat["fixed_y"][k][metals.index(m5)]), frac))
                else:
                    phases.append((str(pdat["phase_ids"][k]), np.nan, frac))
            alb2 = [p for p in phases if p[0].startswith("AlB2")]
            rows.append({
                "system_id": sid, "T_K": t, "x_groupV": X_EXP,
                "dG_hull_eV_per_fu": r["dG_hull"], "dG_hull_meV_per_atom": per_atom(r["dG_hull"]) * 1000,
                "n_phases": len(phases), "n_active_AlB2": len(alb2),
                "AlB2_compositions_x": ";".join(f"{p[1]:.3f}" for p in alb2),
                "phase_list": "; ".join(f"{p[0]} (x={p[1]:.3f})" if not np.isnan(p[1]) else p[0] for p in phases),
                "phase_metal_fractions": ";".join(f"{p[2]:.4f}" for p in phases),
                "ss_fraction_stored": r["ss_fraction"],
                "is_single_phase_strict": bool(r["is_single_phase_strict"]),
                "is_unique_homogeneous": bool(r["is_unique_homogeneous"]),
                "hessian_eigenvalue_eV_per_fu": d2,
                "hessian_note": "binary: 1x1 Hessian = d2G_ss/dx2 (central difference, h = 0.05, stored G_ss)",
                "lp_matches_stored": bool(ok and abs(obj - r["G_lp"]) < 1e-6),
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s1_x050_equilibrium_by_temperature.csv", index=False)
    return df


# ── 2. Fig. 2 white regions ───────────────────────────────────────────────────
def section2() -> pd.DataFrame:
    rows = []
    for sid in ORDER:
        tag, _, m5 = SYSTEMS[sid]
        d = pd.read_parquet(SIM_DIR / f"{tag}_stability.parquet")
        g = d["dG_hull"]
        rows.append({
            "system_id": sid, "n_points": len(g), "n_nan": int(g.isna().sum()),
            "n_negative_roundoff": int((g < 0).sum()), "min_negative_eV": float(g.min()),
            "n_exact_zero": int((g == 0).sum()), "n_0_to_1e-9": int(((g > 0) & (g < 1e-9)).sum()),
            "n_above_1e-9": int((g >= 1e-9).sum()), "max_meV_per_fu": float(g.max() * 1000),
            "T_of_points_above_1e-9": ";".join(f"{t:.0f}" for t in sorted(d.loc[g >= 1e-9, "T"].unique())),
            "negative_points_all_strict": bool(d.loc[g < 0, "is_single_phase_strict"].all()),
            "explanation": "white = ΔG_hull slightly negative (LP round-off ~1e-15 eV) below the lowest contour level 0 "
                           "(contourf extend='max' leaves values below vmin unfilled); colour = 0 to 1e-9. Both mean on the hull.",
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s2_fig2_white_region_audit.csv", index=False)
    return df


# ── 3. SQS / reference / normalization / RK audit ─────────────────────────────
def cell_params(contcar: Path) -> dict:
    lines = contcar.read_text().splitlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(v) for v in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    la, lb, lc = np.linalg.norm(lat, axis=1)
    ang = lambda u, v: np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)))
    return {"|a1|": la, "|a2|": lb, "|a3|": lc, "alpha": ang(lat[1], lat[2]), "beta": ang(lat[0], lat[2]),
            "gamma": ang(lat[0], lat[1]), "volume_A3": abs(np.linalg.det(lat))}


def section3() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, fits = [], []
    for sid in ORDER:
        tag, m4, m5 = SYSTEMS[sid]
        bdir = STRUCT_DIR / tag / "blade"
        metals = detect_metals(bdir, "B")
        for d in sorted(bdir.iterdir()):
            cp, ep = d / "CONTCAR", d / "energy"
            if not (cp.exists() and ep.exists()):
                continue
            lines = cp.read_text().splitlines()
            counts: dict[str, int] = {}
            for s, c in zip(lines[5].split(), map(int, lines[6].split()), strict=True):
                counts[s] = counts.get(s, 0) + c
            n_m = sum(counts.get(m, 0) for m in metals)
            e = float(ep.read_text().split()[0])
            m_name = {m: float(v) for m, v in re.findall(r"a_([A-Z][a-z]?)=([0-9.]+)", d.name)}
            cell = cell_params(cp)
            rows.append({
                "system_id": sid, "structure": d.name, "x_from_dirname": m_name.get(m5, 0.0),
                "x_from_CONTCAR": counts.get(m5, 0) / n_m, "n_metal": n_m, "n_B": counts.get("B", 0),
                "B_per_metal": counts.get("B", 0) / n_m, "other_species": ";".join(s for s in counts if s not in metals + ["B"]),
                "energy_total_eV": e, "energy_eV_per_fu": e / n_m, "energy_eV_per_atom": e / sum(counts.values()),
                "volume_A3_per_fu": cell["volume_A3"] / n_m,
                "cell_alpha": cell["alpha"], "cell_beta": cell["beta"], "cell_gamma": cell["gamma"],
            })
    sqs = pd.DataFrame(rows)
    # mixing enthalpy vs same-supercell endpoints; Vegard deviation of volume
    out = []
    for sid, g in sqs.groupby("system_id", sort=False):
        e0 = g.loc[g.x_from_CONTCAR == 0, "energy_eV_per_fu"].iloc[0]
        e1 = g.loc[g.x_from_CONTCAR == 1, "energy_eV_per_fu"].iloc[0]
        v0 = g.loc[g.x_from_CONTCAR == 0, "volume_A3_per_fu"].iloc[0]
        v1 = g.loc[g.x_from_CONTCAR == 1, "volume_A3_per_fu"].iloc[0]
        g = g.copy()
        x = g.x_from_CONTCAR
        g["dH_mix_meV_per_fu"] = (g.energy_eV_per_fu - ((1 - x) * e0 + x * e1)) * 1000
        g["dH_mix_meV_per_atom"] = g.dH_mix_meV_per_fu / ATOMS_PER_FU
        g["volume_vegard_dev_pct"] = (g.volume_A3_per_fu / ((1 - x) * v0 + x * v1) - 1) * 100
        out.append(g)
        # RK fit audit: refit by least squares with the same basis and report residuals
        mix = g[(x > 0) & (x < 1)]
        xs = mix.x_from_CONTCAR.to_numpy()
        h = mix.dH_mix_meV_per_fu.to_numpy()
        for order in (0, 1, 2, 3):
            A = np.column_stack([xs * (1 - xs) * (1 - 2 * xs) ** k for k in range(order + 1)])  # z = x_IV − x_V
            coef, *_ = np.linalg.lstsq(A, h, rcond=None)
            res = h - A @ coef
            fits.append({"system_id": sid, "rk_order": order, "n_points": len(xs),
                         **{f"L{k}_meV_per_fu": coef[k] for k in range(order + 1)},
                         "rms_residual_meV_per_fu": float(np.sqrt(np.mean(res**2))),
                         "max_abs_residual_meV_per_fu": float(np.abs(res).max()),
                         "basis": "ΔH = x(1−x)·Σ L_k (x_IV − x_V)^k, x = Group-V fraction"})
    sqs = pd.concat(out, ignore_index=True)
    sqs.to_csv(OUT / "s3_sqs_audit.csv", index=False)
    fitdf = pd.DataFrame(fits)
    fitdf.to_csv(OUT / "s3_rk_fit_audit.csv", index=False)
    return sqs, fitdf


# ── 4. common tangent / binodal / spinodal ────────────────────────────────────
def section4() -> pd.DataFrame:
    rows = []
    temps = np.arange(200.0, 4000.0 + 1e-9, 25.0)
    for sid in ORDER:
        tag, _, m5 = SYSTEMS[sid]
        xv, h0, mix = model_gss(tag, m5)
        inner = (xv > 0.005) & (xv < 0.995)
        hmix = h0 - (h0[0] + (h0[-1] - h0[0]) * xv)
        d2h = np.gradient(np.gradient(hmix, xv), xv)
        for t in temps:
            g = h0 + KB * t * mix
            gap = lower_hull_gap(xv, g)
            d2 = np.gradient(np.gradient(g, xv), xv)
            ints = intervals(gap > 1e-6, xv)
            rows.append({"system_id": sid, "T_K": t,
                         "max_height_above_AlB2_hull_meV_per_fu": float(gap.max() * 1000),
                         "min_d2G_dx2_eV_per_fu": float(d2[inner].min()),
                         "x_at_min_d2G": float(xv[inner][np.argmin(d2[inner])]),
                         "binodal_intervals": ";".join(f"{a:.3f}-{b:.3f}" for a, b in ints),
                         "nontrivial_common_tangent": bool(len(ints) > 0),
                         "min_enthalpic_d2H_dx2_eV_per_fu": float(d2h[inner].min())})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s4_common_tangent_check.csv", index=False)
    return df


# ── 5. eutectic vs solid-state: required solid mixing enthalpy ────────────────
def required_regular_omega(x_sol: float, t: float) -> float:
    """Symmetric regular solution: binodal at x_sol at temperature t → Ω (kJ/mol)."""
    return R * t * np.log((1 - x_sol) / x_sol) / (1 - 2 * x_sol)


def section5() -> pd.DataFrame:
    rows = []
    for t in (2000.0, 2500.0, 3000.0, 3284.0):
        om = 2 * R * t  # critical temperature Tc = Ω / 2R
        rows.append({"case": "solid miscibility gap reaching T (Tc = Ω/2R)", "T_K": t, "x_constraint": "",
                     "required_Omega_kJ_per_mol_fu": om, "required_dHmix_0.5_meV_per_fu": om / 4 / EV * 1000,
                     "required_dHmix_0.5_meV_per_atom": om / 4 / EV * 1000 / ATOMS_PER_FU,
                     "assumption": "symmetric regular solution, ideal configurational entropy on metal sublattice, no excess entropy"})
    for label, x_sol, temps in (("Otani FZ ZrB2-TaB2 limit 0.08-0.10 (growth T unknown; near solidus)", (0.08, 0.10), (2500.0, 3000.0, 3284.0)),
                                ("Otani FZ ZrB2-NbB2 limit 0.10-0.125 (near solidus)", (0.10, 0.125), (2500.0, 3000.0, 3301.0)),
                                ("Vorotilo ZrB2-TaB2 max 0.25-0.27 at 2173 K (hot pressing; equilibrium not shown)", (0.25, 0.27), (2173.15,))):
        for t in temps:
            for xs in x_sol:
                om = required_regular_omega(xs, t)
                rows.append({"case": label, "T_K": t, "x_constraint": xs, "required_Omega_kJ_per_mol_fu": om,
                             "required_dHmix_0.5_meV_per_fu": om / 4 / EV * 1000,
                             "required_dHmix_0.5_meV_per_atom": om / 4 / EV * 1000 / ATOMS_PER_FU,
                             "assumption": "symmetric regular solution inverted from a binodal composition"})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s5_required_mixing_enthalpy.csv", index=False)
    return df


# ── 6. Otani / Vorotilo matched comparison ────────────────────────────────────
def section6(bnd: pd.DataFrame) -> pd.DataFrame:
    sol = pd.read_csv(VAL / "processed_data" / "solubility_limits.csv")
    for c in ("x_solute_lower_bound", "x_solute_upper_bound", "temperature_K"):
        sol[c] = pd.to_numeric(sol[c], errors="coerce")
    rows = []
    for rid, sid in (("L01", "ZrNb"), ("L02", "ZrTa"), ("L07", "ZrNb"), ("L08", "ZrTa"), ("L09", "ZrTa")):
        r = sol[sol.record_id == rid].iloc[0]
        t = r.temperature_K
        if np.isnan(t):
            tcomp = [2500.0, 3000.0, 3300.0]
            note = "experimental temperature not stated; model evaluated at 2500/3000/3300 K (AlB2-only model; >2500 K extrapolated)"
        else:
            tcomp = [t]
            note = "matched temperature"
        for tc in tcomp:
            b = bnd[(bnd.system_id == sid) & (np.isclose(bnd.T_K, round(tc / 25) * 25))]
            gaps = b.binodal_intervals.iloc[0] if len(b) else ""
            rows.append({"record_id": rid, "system_id": sid, "experiment": r.bound_type, "method": r.method,
                         "x_exp_lower": r.x_solute_lower_bound, "x_exp_upper": r.x_solute_upper_bound,
                         "T_exp_K": t, "T_model_K": tc,
                         "model_ZrB2_side_solubility": "complete (no AlB2 gap)" if not gaps else gaps,
                         "agreement": "no" if not gaps else "check", "note": note})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s6_otani_vorotilo_matched.csv", index=False)
    return df


# ── 7. ZrB2-TaB2: DFT vs ORB vs Ness ──────────────────────────────────────────
def ness_candidates(x: np.ndarray, o1: float = 234.2, o2: float = -18.3) -> dict[str, np.ndarray]:
    """Plausible subregular conventions for Ω1, Ω2 (kJ/mol), x = x_Ta. Ness 2026 does not state the form."""
    return {
        "RK: x(1-x)[Ω1 + Ω2(x_Zr - x_Ta)]": x * (1 - x) * (o1 + o2 * (1 - 2 * x)),
        "RK: x(1-x)[Ω1 + Ω2(x_Ta - x_Zr)]": x * (1 - x) * (o1 + o2 * (2 * x - 1)),
        "Margules: x(1-x)[Ω1 x_Ta + Ω2 x_Zr]": x * (1 - x) * (o1 * x + o2 * (1 - x)),
        "Margules: x(1-x)[Ω1 x_Zr + Ω2 x_Ta]": x * (1 - x) * (o1 * (1 - x) + o2 * x),
    }


def binodal(x: np.ndarray, h_kj: np.ndarray, t: float) -> list[tuple[float, float]]:
    g = h_kj + R * t * (np.where(x > 0, x * np.log(np.where(x > 0, x, 1)), 0) + np.where(x < 1, (1 - x) * np.log(np.where(x < 1, 1 - x, 1)), 0))
    return intervals(lower_hull_gap(x, g) > 1e-6, x)


def section7(sqs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    dft = pd.read_csv(VAL / "raw_data" / "digitized_figures" / "user_dft_ZrTaB2_mixing_enthalpy.csv")
    dft = dft[(dft.x_nominal > 0) & (dft.x_nominal < 1)]
    orb = sqs[(sqs.system_id == "ZrTa") & (sqs.x_from_CONTCAR > 0) & (sqs.x_from_CONTCAR < 1)]
    rows = []
    for _, d in dft.iterrows():
        o = orb[np.isclose(orb.x_from_CONTCAR, d.x_nominal)]
        x = d.x_nominal
        tang = x * (1 - x) * (-9.306 + 4.016 * (2 * x - 1)) / EV * 1000  # kJ/mol per MB2 -> meV/MB2
        rows.append({"x_Ta": x,
                     "DFT_user_meV_per_atom": d.dH_mix_meV_per_atom, "DFT_user_meV_per_fu": d.dH_mix_meV_per_atom * 3,
                     "DFT_user_uncertainty_meV_per_atom": d.y_uncertainty_meV_per_atom,
                     "ORB_SQS_meV_per_atom": o.dH_mix_meV_per_atom.iloc[0] if len(o) else np.nan,
                     "ORB_SQS_meV_per_fu": o.dH_mix_meV_per_fu.iloc[0] if len(o) else np.nan,
                     "ORB_minus_DFT_meV_per_atom": (o.dH_mix_meV_per_atom.iloc[0] - d.dH_mix_meV_per_atom) if len(o) else np.nan,
                     "ORB_over_DFT": (o.dH_mix_meV_per_atom.iloc[0] / d.dH_mix_meV_per_atom) if len(o) else np.nan,
                     "Tang2025_LDA_RK_meV_per_atom": tang / 3, "Tang2025_LDA_RK_meV_per_fu": tang})
    cmp = pd.DataFrame(rows)
    # RK fits to user DFT and ORB on the same basis
    for name, col in (("DFT_user", "DFT_user_meV_per_fu"), ("ORB_SQS", "ORB_SQS_meV_per_fu")):
        xs, h = cmp.x_Ta.to_numpy(), cmp[col].to_numpy()
        A = np.column_stack([xs * (1 - xs) * (2 * xs - 1) ** k for k in range(3)])
        c, *_ = np.linalg.lstsq(A, h, rcond=None)
        cmp.attrs[name] = c
    pd.DataFrame([{"dataset": k, "L0_meV_per_fu": v[0], "L1_meV_per_fu": v[1], "L2_meV_per_fu": v[2],
                   "L0_kJ_per_mol_fu": v[0] / 1000 * EV, "L1_kJ_per_mol_fu": v[1] / 1000 * EV,
                   "basis": "ΔH = x(1−x)Σ L_k (x_Ta − x_Zr)^k (Tang convention)"} for k, v in cmp.attrs.items()]
                 ).to_csv(OUT / "s7_rk_fits_ZrTaB2.csv", index=False)
    cmp.attrs = {}
    cmp.to_csv(OUT / "s7_dft_vs_orb_ZrTaB2.csv", index=False)

    # Ness Ω under candidate conventions and bases
    x = np.linspace(0, 1, 2001)
    rows = []
    for form, h in ness_candidates(x).items():
        for basis, scale in (("per mol MB2 (= per mol metal)", 1.0), ("per mol atoms", 3.0)):
            hfu = h * scale  # kJ per mol MB2
            i5 = 1000
            b3000 = binodal(x, hfu, 3000.0)
            b2173 = binodal(x, hfu, 2173.15)
            rows.append({"form": form, "Omega_basis": basis,
                         "dHmix_0.5_kJ_per_mol_fu": hfu[i5], "dHmix_0.5_meV_per_fu": hfu[i5] / EV * 1000,
                         "dHmix_0.5_meV_per_atom": hfu[i5] / EV * 1000 / 3,
                         "positive_at_x0.5": bool(hfu[i5] > 0),
                         "min_dHmix_kJ_per_mol_fu": hfu.min(), "x_at_min": x[np.argmin(hfu)],
                         "binodal_3000K": ";".join(f"{a:.4f}-{b:.4f}" for a, b in b3000),
                         "binodal_2173K": ";".join(f"{a:.4f}-{b:.4f}" for a, b in b2173)})
    ness = pd.DataFrame(rows)
    ness.to_csv(OUT / "s7_ness_omega_conventions.csv", index=False)
    return cmp, ness


# ── 8. six-system table ───────────────────────────────────────────────────────
def section8(s1, s4, sqs) -> pd.DataFrame:
    exp = pd.read_csv(VAL / "processed_data" / "experimental_observations.csv")
    th = pd.read_csv(VAL / "processed_data" / "thermochemical_data.csv")
    wen = pd.read_csv(VAL / "raw_data" / "historical_experiments" / "wen2020_extraction.csv")
    rows = []
    for sid in ORDER:
        a = s1[s1.system_id == sid]
        c = s4[s4.system_id == sid]
        e = exp[exp.system_id == sid].iloc[0]
        q = sqs[(sqs.system_id == sid) & np.isclose(sqs.x_from_CONTCAR, 0.5)]
        tang = th[(th.compound_or_system == sid) & th.quantity.str.startswith("DFT ΔHmix(0.5) = L0/4")]
        w = wen[(wen.system == sid) & wen.quantity.str.startswith("ΔHmix")]
        rows.append({
            "system_id": sid,
            "exp_observation": f"lamellae λ={e.lamellar_length_scale} µm after laser melt",
            "exp_phase_separated": True,
            "pred_class_x0.5_200-2500K": "; ".join(f"{k}:{v}" for k, v in a.apply(
                lambda r: "homogeneous" if r.is_unique_homogeneous else ("2+ AlB2" if r.n_active_AlB2 >= 2 else "other phases"), axis=1).value_counts().items()),
            "max_dG_hull_x0.5_meV_per_fu": a.dG_hull_eV_per_fu.max() * 1000,
            "min_hessian_x0.5_eV_per_fu": a.hessian_eigenvalue_eV_per_fu.min(),
            "AlB2_gap_any_T_200-4000K": bool(c.nontrivial_common_tangent.any()),
            "AlB2_gap_T_range_K": (f"{c[c.nontrivial_common_tangent].T_K.min():.0f}-{c[c.nontrivial_common_tangent].T_K.max():.0f}"
                                   if c.nontrivial_common_tangent.any() else "none"),
            "dHmix_0.5_ORB_meV_per_atom": q.dH_mix_meV_per_atom.iloc[0],
            "dHmix_0.5_Tang_LDA_meV_per_atom": tang.value.iloc[0] / EV * 1000 / 3 if len(tang) else np.nan,
            "dHmix_0.5_Wen_PBE_meV_per_atom": float(w.value.iloc[0]) / EV * 1000 if len(w) else np.nan,
            "Wen_molten_salt_1423K_result": w.notes.iloc[0] if len(w) else "",
            "verdict": "MISMATCH (model homogeneous; experiment separated)",
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s8_six_system_validation_table.csv", index=False)
    return df


# ── Figures ───────────────────────────────────────────────────────────────────
def figures(sqs, cmp, ness, req) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def save(fig, name):
        for ext in ("png", "pdf"):
            fig.savefig(FIG / f"{name}.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig)

    # fig9: ZrB2-TaB2 ΔHmix per atom: user DFT, ORB, Tang, Wen, Vorotilo
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x = np.linspace(0, 1, 401)
    ax.plot(cmp.x_Ta, cmp.ORB_SQS_meV_per_atom, "s-", color="tab:blue", label="ORB-relaxed SQS (this work)")
    ax.errorbar(cmp.x_Ta, cmp.DFT_user_meV_per_atom, yerr=cmp.DFT_user_uncertainty_meV_per_atom, fmt="o-", color="tab:orange",
                label="DFT (user; digitized from mixing_enthalpy.png)")
    ax.plot(x, x * (1 - x) * (-9.306 + 4.016 * (2 * x - 1)) / EV * 1000 / 3, "--", color="tab:green", label="Tang 2025 LDA (RK, Table 2)")
    ax.plot([0.5], [-0.873 / EV * 1000], "D", color="tab:purple", label="Wen 2020 PBE (x=0.5; per-atom basis inferred)")
    ax.plot([0.5], [-8.0], "^", color="tab:brown", label="Vorotilo 2022 DFT (x=0.5)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("x in Ta$_x$Zr$_{1-x}$B$_2$")
    ax.set_ylabel("ΔH$_{mix}$ (meV/atom)   [×3 = meV/MB$_2$]")
    ax.set_title("ZrB$_2$–TaB$_2$ mixing enthalpy: ORB vs DFT (all relative to AlB$_2$ end members)", fontsize=9)
    ax.legend(fontsize=7)
    save(fig, "fig9_ZrTaB2_dft_vs_orb")

    # fig10: Ness Ω candidates vs calculations (separate scale)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for form, h in ness_candidates(x).items():
        a1.plot(x, h / EV * 1000 / 3, label=f"{form} (Ω per mol MB$_2$)")
    a1.plot(cmp.x_Ta, cmp.DFT_user_meV_per_atom, "o", color="k", ms=3, label="DFT (user)")
    a1.plot(cmp.x_Ta, cmp.ORB_SQS_meV_per_atom, "s", color="tab:blue", ms=3, label="ORB SQS")
    a1.axhline(0, color="k", lw=0.8)
    a1.set_xlabel("x$_{Ta}$")
    a1.set_ylabel("ΔH$_{mix}$ (meV/atom)")
    a1.set_title("Ness Ω₁=234.2, Ω₂=−18.3 kJ/mol under candidate conventions\n(per-mol-atoms basis would be ×3)", fontsize=8)
    a1.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
    r = req[req.case.str.startswith("solid miscibility")]
    a2.plot(r.T_K, r["required_dHmix_0.5_meV_per_atom"], "k-o", label="gap reaching T (regular, Tc=Ω/2R)")
    for lab, mk in (("Otani FZ ZrB2-TaB2", "v"), ("Otani FZ ZrB2-NbB2", "^"), ("Vorotilo", "D")):
        d = req[req.case.str.startswith(lab)]
        a2.plot(d.T_K, d["required_dHmix_0.5_meV_per_atom"], mk, label=f"{lab} solubility → required (regular)")
    a2.axhspan(-21, -2, color="tab:blue", alpha=0.15, label="calculated range at x=0.5 (ORB, DFT)")
    a2.axhline(0, color="k", lw=0.8)
    a2.set_xlabel("T (K)")
    a2.set_ylabel("ΔH$_{mix}$(0.5) (meV/atom)")
    a2.set_title("Mixing enthalpy required by observed immiscibility (regular-solution estimate)", fontsize=8)
    a2.legend(fontsize=6)
    fig.tight_layout()
    save(fig, "fig10_ness_omega_and_required_enthalpy")

    # fig11: SQS audit — volume Vegard deviation and ΔHmix per atom, all systems
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    for sid, g in sqs.groupby("system_id", sort=False):
        g = g.sort_values("x_from_CONTCAR")
        lab = f"{SYSTEMS[sid][1]}B$_2$–{SYSTEMS[sid][2]}B$_2$"
        a1.plot(g.x_from_CONTCAR, g.dH_mix_meV_per_atom, "o-", ms=3, label=lab)
        a2.plot(g.x_from_CONTCAR, g.volume_vegard_dev_pct, "o-", ms=3, label=lab)
    a1.axhline(0, color="k", lw=0.8)
    a2.axhline(0, color="k", lw=0.8)
    a1.set_xlabel("x (Group-V)")
    a1.set_ylabel("ΔH$_{mix}$ SQS (meV/atom)")
    a2.set_xlabel("x (Group-V)")
    a2.set_ylabel("volume deviation from Vegard (%)")
    a1.legend(fontsize=7)
    a1.set_title("SQS mixing enthalpy (same-supercell end-member references)", fontsize=8)
    a2.set_title("Relaxed SQS volume vs Vegard's law", fontsize=8)
    fig.tight_layout()
    save(fig, "fig11_sqs_audit")


def main() -> None:
    warnings.filterwarnings("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    s1 = section1()
    s2 = section2()
    sqs, fits = section3()
    s4 = section4()
    req = section5()
    s6 = section6(s4)
    cmp, ness = section7(sqs)
    s8 = section8(s1, s4, sqs)
    figures(sqs, cmp, ness, req)
    pd.set_option("display.width", 220)
    print(s2[["system_id", "n_nan", "n_negative_roundoff", "n_above_1e-9", "max_meV_per_fu", "T_of_points_above_1e-9"]].to_string(index=False))
    print(sqs.groupby("system_id")[["B_per_metal", "cell_alpha", "cell_beta", "cell_gamma", "volume_vegard_dev_pct"]].agg(["min", "max"]).round(3).to_string())
    print(fits[fits.rk_order == 3][["system_id", "n_points", "rms_residual_meV_per_fu", "max_abs_residual_meV_per_fu"]].round(3).to_string(index=False))
    print(s4.groupby("system_id").agg(gap=("nontrivial_common_tangent", "any"), min_d2=("min_d2G_dx2_eV_per_fu", "min"),
                                      min_d2H=("min_enthalpic_d2H_dx2_eV_per_fu", "min")).round(4).to_string())
    print(cmp.round(2).to_string(index=False))
    print(ness[["form", "Omega_basis", "dHmix_0.5_meV_per_atom", "positive_at_x0.5", "x_at_min", "binodal_3000K", "binodal_2173K"]].round(2).to_string(index=False))
    print(req.round(1).to_string(index=False))
    print(s6.to_string(index=False))
    print(s8.drop(columns=["verdict"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
