"""Item 1: boron loss required for phase separation, per model, from the completed Stage 2/3 energetics (no new
calculations).

Free energy: m06 `build_system` (each model's own references, RK mixing, dilute vacancy energies, environment slope,
vacancy-vacancy W, line compounds incl. metal-rich borides, ordered ZrMB4 and vacancy-ordered cells). Closed system:
fixed overall metal composition x and boron content B/M = 2(1 - delta); LP Gibbs minimisation with full Zr, M and B
mass balance.
Scenarios: "own" (unmodified MLIP) for both systems; "(y<=0.125)" variants cap the AlB2 vacancy fraction at the DFT
stability limit of TaB2-x so the vacancy model is never extrapolated (deficit beyond it must form a real phase); "DFT-anchored" for ZrTa only (TaB2 dilute V_B set to the
published -0.389 eV, Ektarawong et al. 2023; everything else the model's own).

Per (model, system, scenario, T, x) the scan over delta gives:
  delta_onset   smallest delta at which the equilibrium is no longer a single AlB2 phase (homogeneous diboride
                destabilised); 0 if already multiphase at B/M = 2 (pure-B by-product phases are ignored)
  delta_5       second-phase metal fraction >= 5 %  (lab-XRD / SEM-detectable, assumption)
  delta_30      second-phase metal fraction >= 30 % (substantial, lamella-scale)
  second phase  identity at delta_5 (AlB2 / metal-rich boride / vacancy-ordered cell / ordered ZrMB4) and its
                composition and B/M
"Second phase" = all metal not in the majority AlB2 phase. Validity of the vacancy model (largest vacancy fraction
y in any AlB2 phase, per B site): <= 0.083 calibrated (k = 8 cells), <= 0.125 within the DFT stability range of
TaB2-x (x <= 0.25), larger = extrapolated (flagged, not interpreted).
T range: 300-3000 K; the solid-only model has no liquid and the measured solidus of the 50 % alloys is ~3200-3500 K
(Eckert), so 3000 K is the upper limit used.
Outputs: results/p02_delta_scan_raw.csv, results/p02_critical_deficiency.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
MB = HERE.parent / "mlip_benchmark"
sys.path.insert(0, str(MB / "scripts"))
sys.path.insert(0, str(HERE.parent / "scripts"))
from m06_stage3_thermo import MODELS, build_system, energies  # noqa: E402
from s11_vacancy_sensitivity import perturbed  # noqa: E402

DFT_TAB2_VB = -0.389
DELTAS = [0.0, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.125, 0.15, 0.175, 0.20]
TEMPS = [300.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
XS = [0.10, 0.25, 0.50, 0.75, 0.90]
Y_CAL, Y_DFT = 8 / 96, 0.125


def summarise(ph: pd.DataFrame, comp: pd.DataFrame) -> dict:
    kind = dict(zip(comp.name, comp.kind))
    ph = ph[~ph.phase.astype(str).str.endswith(("_B", "B (mp-160)", "B (mp-161)", "B (mp-1193675)", "B (mp-22046)"))].copy()
    ph = ph[ph.x_M.notna() | (ph.phase == "AlB2")]
    tot = ph.metal_fraction.sum()
    if tot <= 0:
        return {"n_phases": 0}
    al = ph[ph.phase == "AlB2"].sort_values("metal_fraction", ascending=False)
    main = al.iloc[0] if len(al) else None
    others = ph.drop(index=main.name) if main is not None else ph
    second = others.sort_values("metal_fraction", ascending=False).iloc[0] if len(others) else None
    if second is None:
        ident = "none"
    elif second.phase == "AlB2":
        ident = "AlB2 (second solid solution)"
    else:
        ident = {"ordered_vacancy": "vacancy-ordered AlB2 cell", "ordered_metal": "ordered ZrMB4"}.get(kind.get(second.phase), "metal-rich boride")
    ymax = float(al.y_Va.max()) if len(al) else 0.0
    return {"n_phases": len(ph), "second_metal_frac": float(others.metal_fraction.sum() / tot), "second_identity": ident,
            "second_phase_name": "" if second is None else str(second.phase),
            "main_x": float(main.x_M) if main is not None else np.nan, "main_B_per_M": float(2 * (1 - main.y_Va)) if main is not None else np.nan,
            "second_x": np.nan if second is None else float(second.x_M),
            "second_B_per_M": np.nan if second is None or second.phase != "AlB2" else float(2 * (1 - second.y_Va)),
            "y_max": ymax, "validity": "calibrated" if ymax <= Y_CAL + 1e-9 else ("DFT range" if ymax <= Y_DFT + 1e-9 else "extrapolated")}


def main() -> None:
    e, E = energies()
    natoms = dict(zip(e.structure_id, e.n_atoms))
    em = pd.read_csv(MB / "results" / "m02_end_member_benchmarks.csv").set_index(["model", "end_member"])
    rows = []
    for sid in ("ZrNb", "ZrTa"):
        for m in MODELS:
            base = build_system(m, sid, E, natoms)
            scen = {"own": base}
            if sid == "ZrTa":
                s = perturbed(base, d_m=float(DFT_TAB2_VB - em.loc[(m, "TaB2"), "E_f_VB_Brich_eV"]))
                s.compounds = base.compounds
                scen["DFT-anchored"] = s
            # capped variants: AlB2 vacancy fraction limited to the DFT stability range of TaB2-x (y <= 0.125); any further
            # deficit must go into a real competing phase -> no extrapolation of the vacancy model
            for nm in list(scen):
                c = type(scen[nm])(metal=scen[nm].metal, alb2=scen[nm].alb2, compounds=scen[nm].compounds,
                                   x_grid=scen[nm].x_grid, y_grid=scen[nm].y_grid[scen[nm].y_grid <= Y_DFT + 1e-12])
                scen[nm + " (y<=0.125)"] = c
            for sname, s in scen.items():
                for t in TEMPS:
                    for x in XS:
                        for d in DELTAS:
                            ph = s.equilibrium(x, 2 * (1 - d), t)
                            rows.append({"system": sid, "model": m, "scenario": sname, "T_K": t, "x": x, "delta": d,
                                         **summarise(ph, base.compounds)})
                print(f"{sid} {m} done", flush=True)
    raw = pd.DataFrame(rows)
    raw.to_csv(HERE / "results" / "p02_delta_scan_raw.csv", index=False)
    crit = []
    for key, g in raw.groupby(["system", "model", "scenario", "T_K", "x"]):
        g = g.sort_values("delta")
        multi = g[g.n_phases > 1]

        def first(cond):
            h = g[cond]
            return float(h.delta.iloc[0]) if len(h) else np.nan

        d5, d30 = first(g.second_metal_frac >= 0.05), first(g.second_metal_frac >= 0.30)
        at5 = g[g.delta == d5].iloc[0] if d5 == d5 else None
        crit.append(dict(zip(["system", "model", "scenario", "T_K", "x"], key),
                         delta_onset=float(multi.delta.iloc[0]) if len(multi) else np.nan, delta_5=d5, delta_30=d30,
                         second_at_5=None if at5 is None else at5.second_identity, second_name_at_5=None if at5 is None else at5.second_phase_name,
                         main_x_at_5=None if at5 is None else at5.main_x, second_x_at_5=None if at5 is None else at5.second_x,
                         second_BM_at_5=None if at5 is None else at5.second_B_per_M, validity_at_5=None if at5 is None else at5.validity,
                         validity_at_30=(g[g.delta == d30].iloc[0].validity if d30 == d30 else None)))
    c = pd.DataFrame(crit)
    c.to_csv(HERE / "results" / "p02_critical_deficiency.csv", index=False)
    pd.set_option("display.width", 250)
    s = c[c.x == 0.5]
    print(s.pivot_table(index=["system", "scenario", "model"], columns="T_K", values="delta_5").to_string())


if __name__ == "__main__":
    main()
