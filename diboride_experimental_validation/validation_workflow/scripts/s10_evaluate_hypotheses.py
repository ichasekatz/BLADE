"""Step 10: apply the pass/fail criteria (docs/WORKFLOW.md section 7) to whatever results exist.

Every verdict carries its evidence level:
  ORB-screen   from ORB (BLADE settings) only - can REJECT or flag, but SUPPORTED needs DFT confirmation
  literature   from extracted literature values (Tang 2025 QHA)
  DFT          from the prepared DFT calculations (results/dft/*.csv once they exist)
Statuses: SUPPORTED / REJECTED / INCONCLUSIVE / PENDING (inputs missing) / PENDING_EXPERIMENT.

Output: results/s10_hypothesis_status.csv (+ printed table)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

WF = Path(__file__).resolve().parents[1]
RES = WF / "results"
SOL = ("ZrNb", "ZrTa")
THR = {"H1_rec_strong": 0.30, "H1_rec_weak": 0.10, "H2_margin_meV": 2.0, "H4_bm_max": 1.95, "T_min": 2000.0}


def read(name: str) -> pd.DataFrame | None:
    p = RES / name
    return pd.read_csv(p) if p.exists() else None


def h1(sid: str, vac, rec, summ) -> dict:
    if vac is None or rec is None or summ is None:
        return {"status": "PENDING", "evidence": "vacancy relaxations / framework not finished"}
    r = float(rec.loc[rec.system == sid, "dE_rec_eV"].iloc[0])
    v = vac[(vac.system == sid) & (vac.defect == "V_B") & np.isclose(vac.x_solute, 0.5) & (vac.condition == "B-rich") & (vac.T_K_for_mu == 0)]
    slope = float(np.polyfit(v.local_n_solute, v.E_f_eV, 1)[0]) if len(v) >= 3 else np.nan
    s = summ[(summ.system == sid) & (summ.T_K >= THR["T_min"]) & (summ.B_per_metal >= 1.90) & (summ.B_per_metal < 2.0)]
    gap = {ps: bool((g.AlB2_dx > 0.1).any()) for ps, g in s.groupby("parameter_set")}
    s_all = summ[(summ.system == sid) & (summ.T_K >= 1500) & (summ.B_per_metal >= 1.80)]
    any_gap_all = bool((s_all.AlB2_dx > 0.1).any())
    edge = bool(s.get("y_at_grid_edge", pd.Series(dtype=bool)).any())
    stoich = summ[(summ.system == sid) & (summ.T_K >= THR["T_min"]) & np.isclose(summ.B_per_metal, 2.0) & np.isclose(summ.x_overall, 0.5)]
    gap_at_stoich = {ps: bool((g.AlB2_dx > 0.1).any()) for ps, g in stoich.groupby("parameter_set")}
    ev = (f"dE_rec = {r:+.3f} eV (V_B in MB2 minus ZrB2, equal mu_B); dE_v/dn_solute at x=0.5 = {slope:+.3f} eV per neighbour "
          f"(6x = {6 * slope:+.3f} eV); framework two-AlB2 (dx>0.1) at T>=2000 K, 1.90<=B/M<2: {gap}; "
          f"also at exact B/M=2, x=0.5 (vacancies + free B): {gap_at_stoich}; vacancy fraction at grid edge: {edge}")
    em = vac[(vac.system == sid) & np.isclose(vac.x_solute, 1.0) & (vac.defect == "V_B") & (vac.condition == "B-rich") & (vac.T_K_for_mu == 0)]
    ef_em = float(em.E_f_eV.min()) if len(em) else np.nan
    ev += f"; E_f(V_B) in pure MB2 at the B-RICH limit = {ef_em:+.3f} eV"
    if ef_em < 0:
        ev += (" (<0: stoichiometric MB2 unstable to B loss next to elemental B; for TaB2 this AGREES with published DFT, "
               "-0.389 eV, Ektarawong et al. 2023, which also finds TaB2-x (x~0.17-0.25) on the hull)")
    if not any_gap_all and abs(r) < THR["H1_rec_weak"] and abs(6 * slope) < THR["H1_rec_weak"]:
        st = "REJECTED"
    elif not any_gap_all:
        st = "REJECTED"  # coupling exists but does not open an AlB2 gap anywhere in the scanned window
    elif edge:
        st = "INCONCLUSIVE (vacancy fraction at grid edge: model extrapolation, needs multi-vacancy DFT)"
    elif all(gap.values()) and r <= -THR["H1_rec_strong"]:
        st = "SUPPORTED (ORB screen; needs DFT P1)"
    else:
        st = "INCONCLUSIVE"
    return {"status": st, "evidence": ev, "dE_rec_eV": r, "slope_eV_per_nn": slope}


def h4(sid: str, summ) -> dict:
    if summ is None:
        return {"status": "PENDING", "evidence": "framework not run"}
    s = summ[(summ.system == sid) & np.isclose(summ.x_overall, 0.5) & (summ.T_K >= THR["T_min"])]
    hit = s[(s.B_per_metal >= THR["H4_bm_max"]) & (s.B_per_metal < 2.0) & s.M_partition_to_other]
    single = s[(s.B_per_metal >= 1.90) & (s.B_per_metal < 2.0) & (s.n_AlB2 == 1) & (s.other_phases.fillna("") == "")]
    by = {ps: bool(len(hit[hit.parameter_set == ps])) for ps in s.parameter_set.unique()}
    big = s[s.metal_fraction_in_other >= 0.3] if "metal_fraction_in_other" in s else s.iloc[0:0]
    bm_needed = {ps: (float(g.B_per_metal.max()) if len(g) else None) for ps, g in big.groupby("parameter_set")}
    frac195 = s[np.isclose(s.B_per_metal, 1.95)].groupby("parameter_set").metal_fraction_in_other.max().round(3).to_dict() \
        if "metal_fraction_in_other" in s else {}
    ev = (f"x=0.5, T>=2000 K: Group-V-rich non-AlB2 phase at 1.95<=B/M<2 in parameter sets {by}; "
          f"cases where AlB2 alone absorbs the B deficit (B/M>=1.90): {len(single)}; metal fraction in the Group-V-rich phase at "
          f"B/M=1.95: {frac195}; highest B/M giving a lamella-scale (>=30% of metal) second phase: {bm_needed}; phases: "
          + "; ".join(sorted(set(hit.other_phases.dropna())))[:300])
    if by and all(by.values()):
        st = "SUPPORTED (ORB screen; needs DFT P1 competitor + vacancy energies)"
    elif len(single) and not any(by.values()):
        st = ("REJECTED (ORB screen; conditional: small B deficits are absorbed by vacancy-rich AlB2 because ORB vacancy energies are low - "
              "reopens if DFT P1 gives higher V_B energies)")
    else:
        st = "INCONCLUSIVE"
    return {"status": st, "evidence": ev}


def h2(sid: str, conf) -> dict:
    if conf is None:
        return {"status": "PENDING", "evidence": "configuration relaxations not finished"}
    c = conf[conf.system == sid]
    # genuinely segregated (clustered) arrangements only; single-layer superstructures are ORDERED, not clustered
    seg_keys = ("AABB", "2+2", "column", "c-superlattice 2", "c-superlattice 4")
    seg = c[c.config.str.contains("|".join(k.replace("+", r"\+") for k in seg_keys))]
    rnd = c[c.config.str.startswith("quasi-random")]
    if seg.empty or rnd.empty:
        return {"status": "PENDING", "evidence": "segregated or random cells missing"}
    sig = float(np.nanmax([c.sigma_sqs_meV_atom.max(), 0.0]))
    margin = max(2 * sig, THR["H2_margin_meV"])
    best = seg.loc[seg.minus_random_meV_atom.idxmin()]
    lser = c[c.cell.str.startswith("2x2x8")].sort_values("config")
    ordered = c[c.config.str.contains("period 1|MP ordered|superlattice 1")]
    ev = (f"lowest segregated config: {best.config} at x={best.x_solute} -> {best.minus_random_meV_atom:+.2f} meV/atom vs same-cell random "
          f"(margin {margin:.1f}); superlattice dH (meV/atom): " + ", ".join(f"{r.config.split()[-2]}:{r.dH_mix_meV_atom:.1f}" for r in lser.itertuples())
          + f"; ordered vs random: " + ", ".join(f"{r.minus_random_meV_atom:+.1f}" for r in ordered.itertuples()))
    if best.minus_random_meV_atom < -margin:
        st = "SUPPORTED (ORB screen; needs DFT P3)"
    elif (seg.minus_random_meV_atom > margin).all():
        st = "REJECTED (ORB screen; DFT P3 superlattices confirm)"
    else:
        st = "INCONCLUSIVE"
    return {"status": st, "evidence": ev, "best_segregated_minus_random_meV": float(best.minus_random_meV_atom)}


def h3(sid: str, req) -> list[dict]:
    rows = []
    if req is None:
        return [{"hypothesis": "H3-vib", "status": "PENDING"}, {"hypothesis": "H3-el", "status": "PENDING"}]
    r = req[(req.system == sid) & (req.basis.str.startswith("DFT Tang static (LDA)")) & ~req.basis.str.contains("QHA")
            & (req.criterion == "x0.5_two_phase") & (req.T_K == 2500.0)]
    need = float(r["dS_ex_required_x0.5_kB_per_atom"].iloc[0])
    tang_vib = {"ZrNb": +0.14, "ZrTa": +0.43}[sid]  # k_B/atom at 2500 K, s9 of the discrepancy analysis
    rows.append({"hypothesis": "H3-vib", "level": "literature (Tang 2025 QHA)",
                 "status": "REJECTED (literature; independent ORB/DFT phonon check P4 pending)" if tang_vib >= 0 else "INCONCLUSIVE",
                 "evidence": f"required dS_ex(0.5) at 2500 K = {need:+.3f} k_B/atom; Tang QHA dS_vib,ex = {tang_vib:+.2f} k_B/atom (opposite sign)"})
    dg_need = float(r["dg_EF_required_states_per_eV_atom"].iloc[0])
    dos = read("dft_dos_excess.csv")
    if dos is None:
        rows.append({"hypothesis": "H3-el", "level": "-", "status": "PENDING",
                     "evidence": f"needs dg(E_F) excess at x=0.5 <= {dg_need:.2f} states/eV/atom (plus the vibrational offset); DOS inputs prepared (P3)"})
    else:
        dg = float(dos.loc[dos.system == sid, "dg_EF_excess"].iloc[0])
        st = "SUPPORTED" if dg <= dg_need else ("REJECTED" if dg > 0.5 * dg_need else "INCONCLUSIVE")
        rows.append({"hypothesis": "H3-el", "level": "DFT", "status": st, "evidence": f"dg = {dg:.3f} vs required {dg_need:.3f}"})
    return rows


def main() -> None:
    vac, rec, summ = read("s04_vacancy_formation_energies.csv"), read("s04_reciprocal_energy.csv"), read("s08_summary.csv")
    conf, req = read("s05_configuration_energies.csv"), read("s06_free_energy_requirements.csv")
    out = []
    for sid in SOL:
        out.append({"system": sid, "hypothesis": "H1 B-vacancy coupling opens an AlB2 gap", "level": "ORB-screen", **h1(sid, vac, rec, summ)})
        out.append({"system": sid, "hypothesis": "H4 B loss -> Group-V-rich M3B4/M5B6/MB + Zr-rich AlB2", "level": "ORB-screen", **h4(sid, summ)})
        out.append({"system": sid, "hypothesis": "H2 SQS misses coherent clustering", "level": "ORB-screen", **h2(sid, conf)})
        for r in h3(sid, req):
            out.append({"system": sid, **r})
        out.append({"system": sid, "hypothesis": "H5 kinetic segregation in a miscible solid", "level": "experiment",
                    "status": "PENDING_EXPERIMENT",
                    "evidence": "decisive test is experimental: anneal a lamellar sample below the solidus (homogenises -> H5; coarsens with "
                                "fixed tie-line compositions -> equilibrium gap); TEM-EDS/APT lamella compositions"})
    df = pd.DataFrame(out)
    df.to_csv(RES / "s10_hypothesis_status.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 200)
    print(df[["system", "hypothesis", "level", "status", "evidence"]].to_string(index=False))


if __name__ == "__main__":
    main()
