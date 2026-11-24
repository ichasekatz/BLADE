"""Item 5: predicted diffraction signatures of the candidate assemblages, from each model's relaxed structures.

Lattice model (per model, structurally supported range only):
  a(x, y), c(x, y) = host(x) + y * [(1 - x) * d(ZrB2) + x * d(MB2)]
  host(x): relaxed SQS hosts at x = 0, 0.25, 0.5, 0.75, 1 (linear interpolation in x)
  d(.)   : (a, c)(k = 8 spread vacancies) - (a, c)(host) divided by y = 8/96, i.e. calibrated only for y <= 0.083
Candidate assemblages at x = 0.5, 2000 K (compositions from the capped p02 scan, y <= 0.125):
  homogeneous stoichiometric (Zr0.5M0.5)B2
  two AlB2 solid solutions: main (x, B/M) + second (x, B/M) as predicted by each model at the delta where it appears
  AlB2 + metal-rich boride: reflections of the model's relaxed Nb3B4, Nb5B6, NbB, Nb2B3, Ta3B4, Ta5B6, TaB cells
Wavelengths: Cu Ka1 1.5406 A (Eckert: Debye-Scherrer, CuKa; Otani: not stated -> assumed Cu Ka1) and Co Ka1 1.78897 A
(Vorotilo). Ness reported no diffraction of the melts.
Outputs: results/p05_alb2_peak_positions.csv, results/p05_boride_reflections.csv, results/p05_lattice_model.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
MB = HERE.parent / "mlip_benchmark"
sys.path.insert(0, str(MB / "scripts"))
from m07_diffraction_signatures import structure  # noqa: E402

MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
LAMBDA = {"CuKa1": 1.5406, "CoKa1": 1.78897}
HKL = [(0, 0, 1), (1, 0, 0), (1, 0, 1), (0, 0, 2), (1, 1, 0), (1, 0, 2), (1, 1, 1), (2, 0, 0), (2, 0, 1)]


def cellac(model: str, sid: str) -> tuple[float, float]:
    st = structure(model, sid)
    return float(np.mean(st.lattice.abc[:2]) / 4), float(st.lattice.abc[2] / 3)


def two_theta(a: float, c: float, hkl, lam: float) -> float:
    h, k, l = hkl
    inv_d2 = 4 / 3 * (h * h + h * k + k * k) / a**2 + l * l / c**2
    return float(2 * np.degrees(np.arcsin(lam * np.sqrt(inv_d2) / 2)))


def lattice_model(model: str, sid: str):
    m = sid[2:]
    xs = [0.0, 0.25, 0.5, 0.75, 1.0]
    hosts = ["host_ZrB2", f"host_{sid}_x0.25", f"host_{sid}_x0.50", f"host_{sid}_x0.75", f"host_{m}B2"]
    ac = np.array([cellac(model, h) for h in hosts])
    dz = (np.array(cellac(model, "Y_ZrB2_k8_spread")) - ac[0]) / (8 / 96)
    dm = (np.array(cellac(model, f"Y_{m}B2_k8_spread")) - ac[-1]) / (8 / 96)

    def f(x: float, y: float) -> tuple[float, float]:
        a = np.interp(x, xs, ac[:, 0]) + y * ((1 - x) * dz[0] + x * dm[0])
        c = np.interp(x, xs, ac[:, 1]) + y * ((1 - x) * dz[1] + x * dm[1])
        return float(a), float(c)

    return f, ac, dz, dm


def main() -> None:
    crit = pd.read_csv(HERE / "results" / "p02_critical_deficiency.csv")
    raw = pd.read_csv(HERE / "results" / "p02_delta_scan_raw.csv")
    peaks, lat = [], []
    for sid in ("ZrNb", "ZrTa"):
        for m in MODELS:
            f, ac, dz, dm = lattice_model(m, sid)
            lat.append({"system": sid, "model": m, "a_x0": ac[0, 0], "c_x0": ac[0, 1], "a_x1": ac[-1, 0], "c_x1": ac[-1, 1],
                        "da_dy_ZrB2": dz[0], "dc_dy_ZrB2": dz[1], "da_dy_MB2": dm[0], "dc_dy_MB2": dm[1],
                        "c_change_pct_MB2_at_y0.083": 100 * dm[1] * (8 / 96) / ac[-1, 1]})
            scen_list = ["own (y<=0.125)"] + (["DFT-anchored (y<=0.125)"] if sid == "ZrTa" else [])
            for scen in scen_list:
                g = raw[(raw.system == sid) & (raw.model == m) & (raw.scenario == scen) & (raw.T_K == 2000.0) & (raw.x == 0.5)]
                two = g[(g.second_identity == "AlB2 (second solid solution)") & (g.second_metal_frac >= 0.05)].sort_values("delta")
                cases = {"homogeneous (Zr0.5M0.5)B2": [(0.5, 0.0, 1.0)]}
                if len(two):
                    r = two.iloc[0]
                    cases[f"two AlB2 at delta={r.delta:.3f}"] = [(r.main_x, 1 - r.main_B_per_M / 2, 1 - r.second_metal_frac),
                                                               (r.second_x, 1 - r.second_B_per_M / 2, r.second_metal_frac)]
                for case, phases in cases.items():
                    for lam_name, lam in LAMBDA.items():
                        pos = [[two_theta(*f(x, y), hkl, lam) for hkl in HKL] for x, y, _ in phases]
                        for j, hkl in enumerate(HKL):
                            row = {"system": sid, "model": m, "scenario": scen, "case": case, "wavelength": lam_name,
                                   "hkl": "".join(map(str, hkl)), "two_theta_phase1": pos[0][j]}
                            if len(phases) > 1:
                                row.update({"two_theta_phase2": pos[1][j], "splitting_deg": abs(pos[1][j] - pos[0][j]),
                                            "phase1_x": phases[0][0], "phase1_BM": 2 * (1 - phases[0][1]), "phase2_x": phases[1][0],
                                            "phase2_BM": 2 * (1 - phases[1][1]), "phase2_metal_frac": phases[1][2]})
                            peaks.append(row)
    pd.DataFrame(lat).to_csv(HERE / "results" / "p05_lattice_model.csv", index=False)
    pk = pd.DataFrame(peaks)
    pk.to_csv(HERE / "results" / "p05_alb2_peak_positions.csv", index=False)
    # metal-rich boride reflections (Cu Ka1) not within 0.3 deg of any AlB2 reflection of the homogeneous alloy, all models
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    calc = XRDCalculator(wavelength="CuKa1")
    br = []
    for m in MODELS:
        for ref in ("ref_mp-10255_Nb3B4", "ref_mp-1102394_Nb5B6", "ref_mp-2580_NbB", "ref_mp-1080021_Nb2B3", "ref_mp-10142_Ta3B4",
                    "ref_mp-28629_Ta5B6", "ref_mp-1097_TaB"):
            pat = calc.get_pattern(structure(m, ref), two_theta_range=(20, 80))
            sp = pat.y.argsort()[::-1][:6]
            for i in sp:
                br.append({"model": m, "phase": ref.split("_", 2)[-1], "two_theta_CuKa1": round(float(pat.x[i]), 3),
                           "rel_intensity": round(float(pat.y[i]), 1)})
    pd.DataFrame(br).to_csv(HERE / "results" / "p05_boride_reflections.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(lat).round(4).to_string(index=False))
    s = pk[(pk.wavelength == "CuKa1") & pk.case.str.startswith("two")]
    print(s.pivot_table(index=["system", "model", "scenario", "case"], columns="hkl", values="splitting_deg").round(2).to_string())


if __name__ == "__main__":
    main()
