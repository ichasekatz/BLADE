"""Build the machine-readable deliverables of final_results/ from the analysis outputs (no new calculations).

CALCULATION_MANIFEST.csv     every MLIP relaxation (Stage 1, 2, 4) per model: structure, role, system, atoms, status,
                             convergence, steps, energy, result path; plus the BLADE rerun record
MODEL_COMPARISON.csv         per model: controlling energetic quantities and benchmark errors vs published DFT / MP
SIX_SYSTEM_VALIDATION.csv    per system: stoichiometric, experimental, vacancy, equilibrium and robustness summary
SIX_SYSTEM_VALIDATION_numeric.csv  normalised indicators used for fig14
PHASE_EQUILIBRIUM_RESULTS.csv  critical boron losses, phases, compositions, fractions (closed system) + open-system
                             boron-chemical-potential thresholds
UNCERTAINTY_ANALYSIS.csv     flip thresholds, model spreads, reference deviations, numerical checks, missing physics
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
VAL = FR.parent
SD = FR / "supplementary_data"
MB = VAL / "validation_workflow" / "mlip_benchmark"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
SYSTEMS = ["TiNb", "TiTa", "ZrNb", "ZrTa", "HfNb", "HfTa"]


def manifest() -> pd.DataFrame:
    rows = []
    for stage, f in (("Stage1 pilot", MB / "pilot_manifest.csv"), ("Stage2", MB / "stage2_manifest.csv"), ("Stage4", SD / "stage4_manifest.csv")):
        man = pd.read_csv(f)
        for r in man.itertuples():
            for m in MODELS:
                p = MB / "results" / m / r.structure_id / "result.json"
                if p.exists():
                    j = json.loads(p.read_text())
                    rows.append({"stage": stage, "model": m, "structure_id": r.structure_id, "role": r.role, "system": r.system, "n_atoms": j["n_atoms"],
                                 "status": j["status"], "converged": j["converged"], "n_steps": j["n_steps"], "fmax_final_eV_A": j["fmax_final_eV_A"],
                                 "energy_eV": j["energy_eV"], "protocol": "FIRE+FrechetCellFilter fmax 0.001 eV/A <=1000 steps",
                                 "result": str(p.relative_to(VAL))})
                elif m == "orb" and pd.notna(getattr(r, "orb_energy_eV", np.nan)):
                    rows.append({"stage": stage, "model": m, "structure_id": r.structure_id, "role": r.role, "system": r.system, "n_atoms": r.n_atoms,
                                 "status": "ok (ORB screen / BLADE)", "converged": True, "energy_eV": r.orb_energy_eV,
                                 "protocol": "BLADE/materialsframework ORB relax, same settings", "result": "validation_workflow/results/orb or BLADE Files"})
    df = pd.DataFrame(rows)
    df.to_csv(FR / "CALCULATION_MANIFEST.csv", index=False)
    return df


def model_comparison() -> pd.DataFrame:
    p = pd.read_csv(SD / "x01_params.csv")
    em = pd.read_csv(SD / "x01_endmember_vacancy.csv").drop_duplicates(["end_member", "model"])
    vv = pd.read_csv(SD / "x01_vacancy_interactions.csv")
    z = pd.read_csv(SD / "z01_ce_fits.csv") if (SD / "z01_ce_fits.csv").exists() else pd.DataFrame()
    m02 = pd.read_csv(MB / "results" / "m02_end_member_benchmarks.csv")
    hull = pd.read_csv(MB / "results" / "m06_hull_consistency.csv")
    rows = []
    for m in MODELS:
        r = {"model": m}
        for e in ("TiB2", "ZrB2", "HfB2", "NbB2", "TaB2"):
            q = em[(em.model == m) & (em.end_member == e)]
            r[f"E_v_{e}_eV"] = float(q.E_v_Brich_eV.iloc[0]) if len(q) else np.nan
            q = vv[(vv.model == m) & (vv.end_member == e) & (vv.arrangement == "k4_cluster")]
            r[f"cluster_vs_spread_{e}_eV_per_vac"] = float(q.interaction_eV_per_vacancy.iloc[0]) if len(q) else np.nan
        r["TaB2_E_v_error_vs_DFT_eV (Ektarawong -0.389)"] = r["E_v_TaB2_eV"] + 0.389
        for sid in SYSTEMS:
            q = p[(p.model == m) & (p.system == sid) & (p.status == "ok")]
            r[f"dHmix_x0.5_{sid}_meV_atom"] = float(q["dH_mix_x0.5_meV_atom"].iloc[0]) if len(q) else np.nan
            r[f"slope_{sid}_eV_per_nn"] = float(q.slope_eV_per_nn.iloc[0]) if len(q) else np.nan
        mm = m02[m02.model == m].set_index("end_member")
        for e in ("ZrB2", "NbB2", "TaB2"):
            if e in mm.index:
                r[f"dHf_{e}_kJ_mol_fu"] = mm.loc[e, "dHf_kJ_mol_fu"]
                r[f"c_err_vs_Otani_{e}_pct"] = mm.loc[e, "c_err_pct"]
        h = hull[(hull.model == m) & hull.phase.str.startswith("ref_mp") & ~hull.phase.str.endswith("_B")]
        r["max_abs_boride_dHf_dev_vs_MP_eV_atom"] = float((h.dHf_eV_atom - h.MP_dHf_eV_atom).abs().max()) if len(h) else np.nan
        for e in ("TaB2", "NbB2"):
            q = z[(z.model == m) & (z.end_member == e)] if len(z) else z
            if len(q):
                r[f"ordered_{e}-x_min_dE_eV_fu"] = q.min_dE_eV_fu.iloc[0]
                r[f"ordered_{e}-x_below_MB2+B"] = q["ordered_phase_below_MB2+B"].iloc[0]
                r[f"CE_{e}_LOO_meV_per_vac"] = 1000 * q.loo_cv_eV_per_vacancy.iloc[0]
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(FR / "MODEL_COMPARISON.csv", index=False)
    return df


def phase_results() -> pd.DataFrame:
    k = pd.read_csv(SD / "x01_critical.csv")
    o = pd.read_csv(SD / "x01_open_scan.csv")
    first = []
    for key, g in o.groupby(["system", "model", "scenario", "T_K"]):
        g = g.sort_values("dmu_B_eV")
        non = g[g["class"] != "homogeneous"]
        first.append(dict(zip(["system", "model", "scenario", "T_K"], key), open_dmu_B_threshold_eV=float(non.dmu_B_eV.iloc[0]) if len(non) else np.nan,
                          open_first_assemblage=non["class"].iloc[0] if len(non) else "homogeneous to 1 eV",
                          open_phases_at_threshold=non.phases.iloc[0] if len(non) else ""))
    of = pd.DataFrame(first)
    k2 = k.merge(of, on=["system", "model", "scenario", "T_K"], how="left")
    k2.loc[k2.x != 0.5, ["open_dmu_B_threshold_eV", "open_first_assemblage", "open_phases_at_threshold"]] = np.nan
    k2.to_csv(FR / "PHASE_EQUILIBRIUM_RESULTS.csv", index=False)
    return k2


def stoich_text(sid: str) -> str:
    """Assemblage at exactly B/M = 2 (delta = 0), x = 0.5, own energetics, all temperatures, from the closed scan."""
    c = pd.read_csv(SD / "x01_closed_scan.csv", usecols=["system", "model", "scenario", "T_K", "x", "delta", "class", "second_phase"])
    g = c[(c.system == sid) & (c.scenario == "own") & (c.x == 0.5) & (c.delta == 0.0)]
    if g.empty:
        return "n/a"
    bad = g[g["class"] != "homogeneous"]
    n = g.model.nunique()
    if bad.empty:
        return f"single homogeneous AlB2 phase in all {n} models"
    return (f"homogeneous in {n - bad.model.nunique()}/{n} models; exceptions: "
            + "; ".join(f"{m} ({', '.join(sorted(set(b['class'])))})" for m, b in bad.groupby("model")))


def six_system(k2: pd.DataFrame) -> pd.DataFrame:
    p = pd.read_csv(SD / "x01_params.csv")
    p = p[p.status == "ok"]
    s8 = pd.read_csv(VAL / "simulation_comparison" / "discrepancy" / "s8_six_system_validation_table.csv").set_index("system_id")
    p01 = pd.read_csv(VAL / "validation_workflow" / "phase_fraction_analysis" / "results" / "p01_ness_lamellae.csv")
    fl = pd.read_csv(SD / "x02_flip_thresholds.csv") if (SD / "x02_flip_thresholds.csv").exists() else pd.DataFrame()
    name = {"TiNb": "TiB2-NbB2", "TiTa": "TiB2-TaB2", "ZrNb": "ZrB2-NbB2", "ZrTa": "ZrB2-TaB2", "HfNb": "HfB2-NbB2", "HfTa": "HfB2-TaB2"}
    lam = {"TiNb": 0.242, "TiTa": 0.346, "ZrNb": 0.382, "ZrTa": 0.414, "HfNb": 0.404, "HfTa": 0.631}
    wen = {"TiNb": "multiple phases", "TiTa": "single phase", "ZrNb": "multiple phases", "ZrTa": "multiple phases", "HfNb": "multiple phases", "HfTa": "multiple phases"}
    other = {"TiNb": "none beyond Wen", "TiTa": "none beyond Wen", "ZrNb": "Otani: XRD broadening grows on slow cooling; FZ crystals <= ~10 mol% Nb; Levashov 1900 C: no SS at x=0.5; Eckert: continuous SS (low resolution)",
             "ZrTa": "Otani: broadening grows on slow cooling; FZ <= ~8 mol% Ta; Vorotilo 1900 C 10 min: <= 25-27 at% Ta; Ness calorimetric Omega (unpublished)",
             "HfNb": "Eckert: continuous SS (low resolution)", "HfTa": "none beyond Wen"}
    rows, num = [], []
    for sid in SYSTEMS:
        q = p[p.system == sid]
        e = p01[p01.system == name[sid]]
        g = k2[(k2.system == sid) & (k2.x == 0.5) & (k2.T_K == 2000.0)]
        own, anc = g[g.scenario == "own"], g[g.scenario == "DFT-anch"]
        cls = own.class_at_5.value_counts()
        fown = fl[(fl.system == sid) & (fl.scenario == "own")] if len(fl) else fl
        fanc = fl[(fl.system == sid) & (fl.scenario == "DFT-anch")] if len(fl) else fl
        rows.append({
            "system": name[sid],
            "stoich_dHmix_x0.5_meV_atom_DFT (Tang / Wen)": f"{s8.loc[sid, 'dHmix_0.5_Tang_LDA_meV_per_atom']:.1f} / {s8.loc[sid, 'dHmix_0.5_Wen_PBE_meV_per_atom']:.1f}",
            "stoich_dHmix_x0.5_meV_atom_MLIP_range": f"{q['dH_mix_x0.5_meV_atom'].min():.1f} to {q['dH_mix_x0.5_meV_atom'].max():.1f}" if len(q) else "n/a",
            "stoich_prediction (x=0.5, B/M=2, own, 1500-3000 K)": stoich_text(sid),
            "exp_segregation_Ness": "yes: lamellae (SE contrast; EDS unresolved); no diffraction, no boron analysis",
            "exp_lamellar_minority_area_fraction_ROI": "; ".join(f"{min(v, 1 - v):.2f}" for v in e.A_light_otsu) if len(e) else "n/a",
            "exp_wavelength_um (Ness published)": lam[sid], "exp_Wen_1423K": wen[sid], "exp_other": other[sid],
            "evidence_solid_state_gap": "moderate (Otani cooling-rate broadening)" if sid in ("ZrNb", "ZrTa") else "none direct",
            "vacancy_E_v_MB2_eV_range (B-rich)": f"{q.dEv_MB2_eV.min():.2f} to {q.dEv_MB2_eV.max():.2f}" if len(q) else "n/a",
            "vacancy_E_v_hostB2_eV_range": f"{q.dEv_hostB2_eV.min():.2f} to {q.dEv_hostB2_eV.max():.2f}" if len(q) else "n/a",
            "coupling_slope_eV_per_nn_range": f"{q.slope_eV_per_nn.min():.3f} to {q.slope_eV_per_nn.max():.3f}" if len(q) else "n/a",
            "delta_onset_pct_own (2000 K)": f"{100*own.delta_onset.min():.1f}-{100*own.delta_onset.max():.1f}" if len(own) else "n/a",
            "delta_5_pct_own": f"{100*own.delta_5.min():.1f}-{100*own.delta_5.max():.1f}" if len(own) else "n/a",
            "delta_30_pct_own": f"{100*own.delta_30.min():.1f}-{100*own.delta_30.max():.1f}" if len(own) and own.delta_30.notna().any() else ">20 (some/all models)",
            "delta_30_pct_DFT-anch": f"{100*anc.delta_30.min():.1f}-{100*anc.delta_30.max():.1f}" if len(anc) and anc.delta_30.notna().any() else ("n/a (no DFT anchor)" if not sid.endswith("Ta") else ">20"),
            "product_at_delta5_own (models)": "; ".join(f"{k}: {v}" for k, v in cls.items()),
            "product_at_delta5_DFT-anch": "; ".join(f"{k}: {v}" for k, v in anc.class_at_5.value_counts().items()) if len(anc) else "n/a",
            "second_AlB2_B_per_M_range": f"{g.second_BM_at_5.min():.2f}-{g.second_BM_at_5.max():.2f}" if g.second_BM_at_5.notna().any() else "n/a",
            "flip_d_m_eV_own_range": f"{fown.d_m_flip.min():+.2f} to {fown.d_m_flip.max():+.2f}" if len(fown) and fown.d_m_flip.notna().any() else "n/a",
            "flip_d_m_eV_DFT-anch_range": f"{fanc.d_m_flip.min():+.2f} to {fanc.d_m_flip.max():+.2f}" if len(fanc) and fanc.d_m_flip.notna().any() else ("no flip within +-1.5 eV" if len(fanc) else "n/a"),
            "DFT_anchor_available": "TaB2 V_B (Ektarawong 2023)" if sid.endswith("Ta") else "none for NbB2",
        })
        agree = cls.max() / cls.sum() if len(cls) else np.nan
        num.append({"system": sid,
                    "stoich. miscibility (|DFT dH|/40)": min(abs(s8.loc[sid, "dHmix_0.5_Tang_LDA_meV_per_atom"]) / 40, 1),
                    "exp. segregation observed": 1.0,
                    "vacancy coupling (|slope|/0.2)": float(abs(q.slope_eV_per_nn).mean() / 0.2) if len(q) else np.nan,
                    "B-loss sensitivity (1-δ5/0.1)": float(max(0, 1 - own.delta_5.median() / 0.10)) if len(own) else np.nan,
                    "model agreement on product": float(agree),
                    "DFT anchor available": 1.0 if sid.endswith("Ta") else 0.0})
    df = pd.DataFrame(rows)
    df.to_csv(FR / "SIX_SYSTEM_VALIDATION.csv", index=False)
    pd.DataFrame(num).to_csv(SD / "SIX_SYSTEM_VALIDATION_numeric.csv", index=False)
    return df


def uncertainty() -> pd.DataFrame:
    rows = []
    p = pd.read_csv(SD / "x01_params.csv")
    p = p[p.status == "ok"]
    for sid in SYSTEMS:
        q = p[p.system == sid]
        for col, lab in (("dEv_MB2_eV", "Group-V V_B energy (eV)"), ("dEv_hostB2_eV", "host V_B energy (eV)"), ("slope_eV_per_nn", "coupling slope (eV/nn)"),
                         ("W_MB2_eV_fu", "W Group-V (eV/fu)"), ("W_hostB2_eV_fu", "W host (eV/fu)"), ("dH_mix_x0.5_meV_atom", "dHmix(0.5) (meV/atom)")):
            if len(q):
                rows.append({"category": "model disagreement (5 MLIPs)", "system": sid, "quantity": lab, "min": q[col].min(), "max": q[col].max(),
                             "spread": q[col].max() - q[col].min(), "note": "spread across models; not a calibrated confidence interval"})
    em = pd.read_csv(SD / "x01_endmember_vacancy.csv").drop_duplicates(["end_member", "model"])
    t = em[em.end_member == "TaB2"]
    rows.append({"category": "reference deviation", "system": "TaB2", "quantity": "V_B energy - DFT (-0.389 eV)", "min": (t.E_v_Brich_eV + 0.389).min(),
                 "max": (t.E_v_Brich_eV + 0.389).max(), "spread": np.nan, "note": "Ektarawong et al. 2023 (GGA, bcc Ta + alpha-B references)"})
    hull = pd.read_csv(MB / "results" / "m06_hull_consistency.csv")
    h = hull[hull.phase.str.startswith("ref_mp") & ~hull.phase.str.endswith("_B")].copy()
    h["dev"] = h.dHf_eV_atom - h.MP_dHf_eV_atom
    rows.append({"category": "reference deviation", "system": "Nb-B, Ta-B borides", "quantity": "boride formation energy - MP GGA (eV/atom)",
                 "min": h.dev.min(), "max": h.dev.max(), "spread": np.nan, "note": "same structures; MP GGA is itself a DFT approximation"})
    fl = SD / "x02_flip_thresholds.csv"
    if fl.exists():
        f = pd.read_csv(fl)
        for r in f.itertuples():
            for knob, unit in (("d_m", "eV"), ("d_h", "eV"), ("slope", "x"), ("w", "x"), ("comp", "eV/atom")):
                v = getattr(r, f"{knob}_flip")
                rows.append({"category": "flip threshold (x=0.5, 2000 K, B/M 1.95)", "system": r.system, "quantity": f"{knob} ({unit}) [{r.model}, {r.scenario}]",
                             "min": v, "max": v, "spread": np.nan, "note": f"reference: {r.reference_class}; flips to {getattr(r, f'{knob}_flip_to')}"})
    rows += [{"category": "numerical", "system": "all", "quantity": "site-to-site spread at equal n_solute (eV)", "min": 0.06, "max": 0.19, "spread": np.nan,
              "note": "replicate V_B sites, Stage 1"},
             {"category": "numerical", "system": "all", "quantity": "relaxation convergence", "min": np.nan, "max": np.nan, "spread": np.nan,
              "note": "see CALCULATION_MANIFEST (converged flag, fmax)"},
             {"category": "missing physics", "system": "all", "quantity": "liquid phase", "min": np.nan, "max": np.nan, "spread": np.nan,
              "note": "lamellae form on solidification; model is solid-only (T <= 3000 K)"},
             {"category": "missing physics", "system": "all", "quantity": "defect vibrational / electronic entropy", "min": np.nan, "max": np.nan, "spread": np.nan,
              "note": "not computed; Tang QHA vibrational mixing correction applied as a scenario only"},
             {"category": "missing physics", "system": "all", "quantity": "vacancy-vacancy interaction beyond y = 0.083", "min": np.nan, "max": np.nan,
              "spread": np.nan, "note": "mean-field W fitted to k <= 8 cells; AlB2 vacancy fraction capped at 0.125"}]
    df = pd.DataFrame(rows)
    df.to_csv(FR / "UNCERTAINTY_ANALYSIS.csv", index=False)
    return df


def main() -> None:
    print(len(manifest()), "calculation records")
    model_comparison()
    k2 = phase_results()
    six_system(k2)
    uncertainty()
    print("deliverables written")


if __name__ == "__main__":
    main()
