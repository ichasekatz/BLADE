"""Step 8: evaluate the thermodynamic framework (thermo_framework.py) with ORB-screening parameters.

Parameter sets (all from files in this workflow; BLADE outputs are read, never written):
  ORB          AlB2 end members + RK mixing from BLADE's ORB SQS energies; vacancy energies from s03/s04;
               competitors = every ORB-relaxed MP phase in Zr-M-B (s03).
  ORB-vac+DFT  same, but the stoichiometric mixing RK replaced by DFT (Tang static L0*, L1*), i.e. ORB's
               2-2.5x over-binding removed (vacancy and competitor energies still ORB).
Scans overall x = 0.05..0.95, B/M = 2.00, 1.98, 1.95, 1.90, 1.80, 1.60 and T = 300, 1500, 2000, 2500, 3000 K.

Optional liquid (--liquid): pseudo-binary B/M = 2 only, dS_fus scanned over 1R, 2R, 3R per atom (no measured
values), W_L = 0; solidus T(x) compared with Eckert incipient-melting data (typed constraints, s07).

Outputs: results/s08_equilibria_<set>.csv, results/s08_summary.csv (+ results/s08_solidus_scan.csv with --liquid)
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from thermo_framework import KB, AlB2Params, LiquidParams, System  # noqa: E402

WF = Path(__file__).resolve().parents[1]
RES = WF / "results"
SOL = {"ZrNb": "Nb", "ZrTa": "Ta"}
EV_PER_KJMOL = 1 / 96.485
TANG = {"ZrNb": [-2.515, 2.418], "ZrTa": [-9.306, 4.016]}  # kJ/mol MB2
TM_RUDY_C = {"ZrB2": 3245.0, "NbB2": 3036.0, "TaB2": 3037.0}  # Eckert Figs. 42/44/46 labels (Rudy & Windisch)


def composition(formula: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}


def build(sid: str, mixing: str) -> System:
    m = SOL[sid]
    e = pd.read_csv(RES / "s03_orb_energies.csv").set_index("structure_id")
    audit = pd.read_csv(RES / "s00_existing_structure_audit.csv")
    h = audit[(audit.system == sid) & (audit.origin == "SQS (BLADE/ATAT)") & audit.path.str.contains(
        f"system_structures/{'NbZrB' if sid == 'ZrNb' else 'TaZrB'}/blade")].drop_duplicates("x_solute").sort_values("x_solute")
    xs, efu = h.x_solute.values, h.energy_eV_atom.values * 3
    e_b = float(e.loc["comp_mp-160_B", "energy_eV_atom"])
    if mixing == "ORB":
        mix = efu - (1 - xs) * efu[0] - xs * efu[-1]
        inner = (xs > 0) & (xs < 1)
        a = np.vstack([xs * (1 - xs) * (2 * xs - 1) ** k for k in range(3)]).T[inner]
        rk = list(np.linalg.lstsq(a, mix[inner], rcond=None)[0])
    else:
        rk = [v * EV_PER_KJMOL for v in TANG[sid]]
    vac = pd.read_csv(RES / "s04_vacancy_formation_energies.csv")
    v = vac[(vac.system == sid) & (vac.defect == "V_B") & (vac.condition == "B-rich") & (vac.T_K_for_mu == 0)]
    dev = lambda x: float(v[np.isclose(v.x_solute, x)].E_f_eV.mean())
    mid = v[np.isclose(v.x_solute, 0.5)]
    fit = np.polyfit(mid.local_n_solute, mid.E_f_eV, 1) if len(mid) >= 3 else None
    ev_env = (lambda n, f=fit: float(np.polyval(f, n))) if fit is not None else None
    w = {"ZrB2": 0.0, f"{m}B2": 0.0}
    vvp = RES / "s04_vacancy_interaction.csv"
    if vvp.exists():
        vv = pd.read_csv(vvp)
        vv = vv[vv.arrangement == "spread"]
        for em in w:  # least-squares W through the origin over k = 2, 4, 8
            g = vv[vv.end_member == em]
            if len(g):
                w[em] = float((g.excess_eV_cell / 48 * g.y_Va**2).sum() / (g.y_Va**4).sum())
    alb2 = AlB2Params(g_zrb=efu[0], g_mb=efu[-1], dEv_zr=dev(0.0), dEv_m=dev(1.0), e_b=e_b, rk=rk, ev_env=ev_env,
                      w_vv_zr=w["ZrB2"], w_vv_m=w[f"{m}B2"])
    rows = []
    for sidx, r in e.iterrows():
        if not sidx.startswith("comp_"):
            continue
        formula = sidx.split("_")[2]
        c = composition(formula)
        if set(c) - {"Zr", m, "B"}:
            continue
        n_fu = r.n_atoms / sum(c.values())
        if formula in ("ZrB2", f"{m}B2"):
            continue  # represented by the AlB2 phase (BLADE SQS end members)
        rows.append({"name": f"{formula} ({sidx.split('_')[1]})", "n_Zr": c.get("Zr", 0), "n_M": c.get(m, 0), "n_B": c.get("B", 0),
                     "G": r.energy_eV / n_fu})
    return System(metal=m, alb2=alb2, compounds=pd.DataFrame(rows))


def scan(system: System, label: str, sid: str) -> pd.DataFrame:
    out = []
    for t in (300.0, 1500.0, 2000.0, 2500.0, 3000.0):
        for bm in (2.0, 1.98, 1.95, 1.90, 1.80, 1.60):
            for x in np.round(np.arange(0.05, 0.951, 0.05), 3):
                ph = system.equilibrium(x, bm, t)
                ph = ph.assign(system=sid, parameter_set=label, T_K=t, B_per_metal=bm, x_overall=x)
                out.append(ph)
    return pd.concat(out, ignore_index=True)


def summarise(eq: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, g in eq.groupby(["system", "parameter_set", "T_K", "B_per_metal", "x_overall"]):
        al = g[g.phase == "AlB2"]
        other = g[~g.phase.isin(["AlB2"])]
        rows.append(dict(zip(["system", "parameter_set", "T_K", "B_per_metal", "x_overall"], key),
                         n_AlB2=len(al), AlB2_x=";".join(f"{v:.3f}" for v in al.x_M), AlB2_yVa=";".join(f"{v:.4f}" for v in al.y_Va),
                         AlB2_dx=float(al.x_M.max() - al.x_M.min()) if len(al) > 1 else 0.0,
                         other_phases=";".join(f"{p}[x={xm:.2f}]" for p, xm in zip(other.phase, other.x_M)),
                         M_partition_to_other=bool(len(other) and (other.x_M > g.x_overall.iloc[0] + 0.05).any()),
                         metal_fraction_in_other=float(other.metal_fraction.sum()) if len(other) else 0.0,
                         y_at_grid_edge=bool(len(al) and (al.y_Va > 0.30).any())))
    return pd.DataFrame(rows)


def solidus_scan(sid: str) -> pd.DataFrame:
    rows = []
    for mult in (1.0, 2.0, 3.0):
        ds = mult * 3 * KB  # per formula unit (3 atoms), k_B per atom * mult
        base = build(sid, "ORB")
        base.liquid = LiquidParams(tm_zrb2=TM_RUDY_C["ZrB2"] + 273.15, tm_mb2=TM_RUDY_C[f"{SOL[sid]}B2"] + 273.15,
                                   ds_fus_zrb2=ds, ds_fus_mb2=ds)
        base.y_grid = np.array([0.0])
        for x in np.round(np.arange(0.1, 0.91, 0.1), 2):
            lo, hi = 2500.0, 3700.0
            for _ in range(25):
                mid = 0.5 * (lo + hi)
                ph = base.equilibrium(x, 2.0, mid)
                lo, hi = (lo, mid) if (ph.phase == "LIQUID").any() else (mid, hi)
            rows.append({"system": sid, "dS_fus_kB_per_atom": mult, "x": x, "solidus_K": hi})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--liquid", action="store_true")
    args = ap.parse_args()
    allq = []
    for sid in SOL:
        for label in ("ORB", "ORB-vac+DFT"):
            eq = scan(build(sid, "ORB" if label == "ORB" else "DFT"), label, sid)
            eq.to_csv(RES / f"s08_equilibria_{sid}_{label}.csv", index=False)
            allq.append(eq)
    summ = summarise(pd.concat(allq))
    summ.to_csv(RES / "s08_summary.csv", index=False)
    pd.set_option("display.width", 250)
    flag = summ[(summ.AlB2_dx > 0.1) | summ.M_partition_to_other]
    print(f"{len(summ)} equilibria; {int((summ.AlB2_dx > 0.1).sum())} with two AlB2 phases (dx > 0.1); "
          f"{int(summ.M_partition_to_other.sum())} with Group-V partitioning into a non-AlB2 phase")
    print(summ[summ.x_overall.isin([0.25, 0.5, 0.75]) & summ.T_K.isin([2000.0, 3000.0])].to_string(index=False, max_colwidth=70))
    if args.liquid:
        sol = pd.concat([solidus_scan(s) for s in SOL])
        sol.to_csv(RES / "s08_solidus_scan.csv", index=False)
        print(sol.to_string(index=False))


if __name__ == "__main__":
    main()
