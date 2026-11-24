"""Step 4: chemical-potential limits and vacancy formation energies (ORB screening level).

Conventions (documented in docs/WORKFLOW.md section 3):
* Energies are total ORB energies (eV) with BLADE's relaxation settings. Element references: mp-160 alpha-B,
  mp-131 Zr, mp-75 Nb, mp-569794 Ta (the references BLADE uses), relaxed in s03.
* Host (Zr1-x Mx)B2 is treated as a pseudo-compound A B2 with mu_A = (1-x) mu_Zr + x mu_M and
  mu_A + 2 mu_B = e_host (energy per formula unit of the pristine host).
* The metal exchange potential dmu = mu_M - mu_Zr is set by the host phase itself:
  dmu(x, T) = d e_AlB2/dx (from a Redlich-Kister fit of BLADE's SQS energies) + k_B T ln(x / (1-x)).
  T = 0 K (enthalpy only) and T = 2500 K (ideal configurational term added) are reported.
* mu_B is the only free variable. Its range is bounded by
    B-rich:  mu_B <= e(B) and every B-rich competitor (ZrB12 ...),
    B-poor:  every metal-rich competitor P (MB, M3B4, M5B6, M2B3, M3B2, M2B, metals):
             nZr mu_Zr + nM mu_M + nB mu_B <= E_P  -> lower bound on mu_B because nB < 2 (nZr + nM).
  Phases with B/M = 2 (ZrB2, MB2, ordered ZrMB4) do not bound mu_B; they test the stoichiometric stability of
  the host (reported separately as host_above_B2_hull).
* E_f(V_B) = E(host - B) - E(host) + mu_B ;  E_f(V_X) = E(host - X) - E(host) + mu_X  (X = Zr, Nb, Ta).
* The chemical-potential-independent reciprocal quantity
  dE_rec = [E(MB2 - B) - E(MB2)] - [E(ZrB2 - B) - E(ZrB2)] compares B-vacancy energies at equal mu_B.

Outputs: results/s04_mu_B_limits.csv, results/s04_vacancy_formation_energies.csv, results/s04_orb_hull_binaries.csv
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

WF = Path(__file__).resolve().parents[1]
RES = WF / "results"
KB = 8.617333262e-5
SOL = {"ZrNb": "Nb", "ZrTa": "Ta"}
REF_ID = {"B": "comp_mp-160_B", "Zr": "comp_mp-131_Zr", "Nb": "comp_mp-75_Nb", "Ta": "comp_mp-569794_Ta"}


def composition(formula: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}


def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    e = pd.read_csv(RES / "s03_orb_energies.csv")
    man = pd.read_csv(WF / "inputs" / "structure_manifest.csv")
    audit = pd.read_csv(RES / "s00_existing_structure_audit.csv")
    return e, man, audit


def host_energies(audit: pd.DataFrame) -> pd.DataFrame:
    """Pristine BLADE SQS energies per formula unit MB2 (reused, not recomputed)."""
    h = audit[(audit.origin == "SQS (BLADE/ATAT)") & audit.path.str.contains("system_structures/(?:NbZrB|TaZrB)/blade", regex=True)]
    rows = []
    for r in h.itertuples():
        sid = "ZrNb" if "NbZrB" in r.path else "ZrTa"
        rows.append({"system": sid, "x_solute": round(r.x_solute, 4), "E_cell": r.energy_eV_cell, "n_atoms": r.n_atoms,
                     "e_fu": r.energy_eV_cell / (r.n_atoms / 3), "path": r.path})
    return pd.DataFrame(rows).drop_duplicates(["system", "x_solute"]).sort_values(["system", "x_solute"])


def rk_dedx(hosts: pd.DataFrame, sid: str, x: float) -> tuple[float, float]:
    """d e_fu/dx at x from a 3-term RK fit of the mixing energy; also returns e_fu(x) from the fit."""
    h = hosts[hosts.system == sid]
    xs, es = h.x_solute.values, h.e_fu.values
    e0, e1 = es[xs == 0][0], es[xs == 1][0]
    mix = es - (1 - xs) * e0 - xs * e1
    a = np.vstack([xs * (1 - xs) * (1 - 2 * xs) ** k for k in range(3)]).T
    lk = np.linalg.lstsq(a[(xs > 0) & (xs < 1)], mix[(xs > 0) & (xs < 1)], rcond=None)[0]
    f = lambda t: (1 - t) * e0 + t * e1 + sum(l * t * (1 - t) * (1 - 2 * t) ** k for k, l in enumerate(lk))
    return (f(x + 1e-5) - f(x - 1e-5)) / 2e-5, f(x)


def competitor_table(e: pd.DataFrame) -> pd.DataFrame:
    c = e[e.structure_id.str.startswith("comp_")].copy()
    c["formula"] = c.structure_id.str.split("_").str[2]
    c["comp"] = c.formula.map(composition)
    rows = []
    for r in c.itertuples():
        # number of formula units in the relaxed cell
        n_fu = r.n_atoms / sum(r.comp.values())
        rows.append({"structure_id": r.structure_id, "formula": r.formula, "E_fu": r.energy_eV / n_fu, **{f"n_{k}": v for k, v in r.comp.items()},
                     "converged": r.converged})
    t = pd.DataFrame(rows).fillna(0)
    for el in ("Zr", "Nb", "Ta", "B"):
        if f"n_{el}" not in t:
            t[f"n_{el}"] = 0
    return t


def binary_hull(comp: pd.DataFrame, refs: dict[str, float]) -> pd.DataFrame:
    """Formation energies and hull distances of the M-B binaries (ORB), for documentation."""
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.core import Composition

    out = []
    for m in ("Zr", "Nb", "Ta"):
        sub = comp[(comp[[f"n_{o}" for o in ("Zr", "Nb", "Ta") if o != m]].sum(axis=1) == 0)]
        entries = [PDEntry(Composition({m: r[f"n_{m}"], "B": r.n_B}), r.E_fu, name=r.structure_id) for _, r in sub.iterrows()
                   if r[f"n_{m}"] + r.n_B > 0]
        pd_ = PhaseDiagram(entries)
        for en in entries:
            out.append({"system": f"{m}-B", "structure_id": en.name, "formula": en.composition.reduced_formula,
                        "x_B": en.composition.get_atomic_fraction("B"),
                        "dHf_eV_atom": pd_.get_form_energy_per_atom(en), "e_above_hull_eV_atom": pd_.get_e_above_hull(en)})
    return pd.DataFrame(out).sort_values(["system", "x_B"])


def mu_b_limits(comp: pd.DataFrame, refs: dict[str, float], sid: str, x: float, e_fu: float, dmu: float) -> dict:
    m = SOL[sid]
    lower, upper, lower_by, upper_by = -np.inf, refs["B"], "", "B (mp-160)"
    for r in comp.itertuples():
        nzr, nm, nb = getattr(r, "n_Zr"), getattr(r, f"n_{m}"), r.n_B
        other = sum(getattr(r, f"n_{o}") for o in ("Zr", "Nb", "Ta") if o not in ("Zr", m))
        if other or (nzr + nm == 0):
            continue
        if (x == 0.0 and nm > 0) or (x == 1.0 and nzr > 0):  # end member: only phases of its own metal
            continue
        nmet = nzr + nm
        coef = nb - 2 * nmet
        # nZr(muA - x dmu) + nM(muA + (1-x) dmu) + nB muB <= E ; muA = e_fu - 2 muB
        rhs = r.E_fu - nmet * e_fu - dmu * (nm * (1 - x) - nzr * x)
        if abs(coef) < 1e-9:
            continue
        bound = rhs / coef
        if coef < 0 and bound > lower:  # metal-rich phase -> lower bound
            lower, lower_by = bound, r.formula
        elif coef > 0 and bound < upper:  # B-rich phase -> upper bound
            upper, upper_by = bound, r.formula
    return {"mu_B_rich": upper, "B_rich_set_by": upper_by, "mu_B_poor": lower, "B_poor_set_by": lower_by,
            "mu_B_window_eV": upper - lower}


def main() -> None:
    e, man, audit = load()
    hosts = host_energies(audit)
    comp = competitor_table(e)
    refs = {el: float(e.loc[e.structure_id == rid, "energy_eV_atom"].iloc[0]) for el, rid in REF_ID.items()
            if (e.structure_id == rid).any()}
    binary_hull(comp, refs).to_csv(RES / "s04_orb_hull_binaries.csv", index=False)

    lim_rows, vac_rows = [], []
    for sid in SOL:
        for x in (0.0, 0.25, 0.5, 0.75, 1.0):
            hrow = hosts[(hosts.system == sid) & np.isclose(hosts.x_solute, x)].iloc[0]
            for t in (0.0, 2500.0):
                dmu_h, e_fit = rk_dedx(hosts, sid, min(max(x, 1e-4), 1 - 1e-4))
                dmu = dmu_h + (KB * t * np.log(x / (1 - x)) if 0 < x < 1 and t > 0 else 0.0)
                lim = mu_b_limits(comp, refs, sid, x, hrow.e_fu, dmu)
                lim_rows.append({"system": sid, "x_solute": x, "T_K": t, "dmu_M_minus_Zr_eV": dmu, "e_host_fu": hrow.e_fu, **lim})
                mu = {"B": None}
                for cond, mub in (("B-rich", lim["mu_B_rich"]), ("B-poor", lim["mu_B_poor"])):
                    mu_a = hrow.e_fu - 2 * mub
                    mu_zr, mu_m = mu_a - x * dmu, mu_a + (1 - x) * dmu
                    sub = man[(man.family == "V") & (man.system.isin([sid, "ZrB2"] if x == 0 else [sid])) & np.isclose(man.x_solute, x)]
                    for r in sub.itertuples():
                        er = e[e.structure_id == r.structure_id]
                        if er.empty:
                            continue
                        species = r.defect.split("_")[1]
                        mu_x = {"B": mub, "Zr": mu_zr}.get(species, mu_m)
                        ef = float(er.energy_eV.iloc[0]) - hrow.E_cell + mu_x
                        vac_rows.append({"structure_id": r.structure_id, "system": sid, "x_solute": x, "T_K_for_mu": t,
                                         "condition": cond, "defect": r.defect, "local_n_solute": r.local_n_solute,
                                         "local_env_def": r.local_env_def, "mu_used_eV": mu_x, "E_f_eV": ef,
                                         "relax_converged": bool(er.converged.iloc[0])})
    lim_df = pd.DataFrame(lim_rows)
    lim_df.to_csv(RES / "s04_mu_B_limits.csv", index=False)
    vac = pd.DataFrame(vac_rows)
    if not vac.empty:
        vac.to_csv(RES / "s04_vacancy_formation_energies.csv", index=False)
        # chemical-potential-independent reciprocal energy (end members)
        def dv(sid_: str, x: float) -> float:
            sub = vac[(vac.system == sid_) & np.isclose(vac.x_solute, x) & (vac.defect == "V_B") & (vac.condition == "B-rich") & (vac.T_K_for_mu == 0)]
            return float(sub.E_f_eV.iloc[0] - sub.mu_used_eV.iloc[0]) if len(sub) else np.nan
        rec = [{"system": s, "dE_rec_eV": dv(s, 1.0) - dv(s, 0.0),
                "meaning": "E(V_B in MB2) - E(V_B in ZrB2) at equal mu_B; < 0 -> vacancies prefer Group-V-rich regions"} for s in SOL]
        pd.DataFrame(rec).to_csv(RES / "s04_reciprocal_energy.csv", index=False)
        print(pd.DataFrame(rec).to_string(index=False))
        piv = vac[vac.T_K_for_mu == 0].pivot_table(index=["system", "x_solute", "defect", "local_n_solute"], columns="condition", values="E_f_eV")
        print(piv.round(3).to_string())
    print(lim_df.round(4).to_string(index=False))
    vv = vacancy_interaction(e, man, audit)
    if vv is not None:
        vv.to_csv(RES / "s04_vacancy_interaction.csv", index=False)
        print(vv.round(4).to_string(index=False))


def vacancy_interaction(e: pd.DataFrame, man: pd.DataFrame, audit: pd.DataFrame) -> pd.DataFrame | None:
    """W (eV per formula unit) in dG_excess = W y'^2 from k-vacancy cells: excess = E(k) - E0 - k e1, y' = k / 96."""
    y = man[man.family == "Y"]
    if y.empty or not y.structure_id.isin(e.structure_id).all():
        return None
    hosts = host_energies(audit)
    single = {"ZrB2": ("ZrNb", 0.0, "V_ZrB2_x0.00_VB_n0a"), "NbB2": ("ZrNb", 1.0, "V_ZrNb_x1.00_VB_n6a"),
              "TaB2": ("ZrTa", 1.0, "V_ZrTa_x1.00_VB_n6a")}
    rows = []
    for r in y.itertuples():
        sid, x, one = single[r.system]
        e0 = float(hosts[(hosts.system == sid) & np.isclose(hosts.x_solute, x)].E_cell.iloc[0])
        e1 = float(e.loc[e.structure_id == one, "energy_eV"].iloc[0]) - e0
        k = int(r.defect.split()[0])
        ek = float(e.loc[e.structure_id == r.structure_id, "energy_eV"].iloc[0]) - e0
        excess = ek - k * e1
        yv = k / 96
        rows.append({"end_member": r.system, "structure_id": r.structure_id, "k_vacancies": k, "arrangement": r.defect.split()[-1],
                     "y_Va": yv, "excess_eV_cell": excess, "excess_per_vacancy_eV": excess / k,
                     "W_eV_per_fu": excess / 48 / yv**2})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
