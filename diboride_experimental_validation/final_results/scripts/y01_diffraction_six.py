"""Six-system diffraction predictions (Cu Ka1 1.5406 A; Co Ka1 1.78897 A for comparison with Vorotilo).

Same lattice model as validation_workflow/phase_fraction_analysis/scripts/p05_diffraction.py, extended to the Ti and Hf
hosts with the Stage 4 relaxations:
  a(x, y), c(x, y) = host(x) + y * [(1 - x) * d(host B2) + x * d(MB2)]
  host(x): each model's relaxed end members and x = 0.25 / 0.5 / 0.75 SQS (linear interpolation in x)
  d(.)   : [(a, c)(k = 8 spread vacancies) - (a, c)(end member)] / (8/96)  (calibrated to y <= 0.083)
Phases: the equilibrium assemblage at delta_5 (x = 0.5, 2000 K) from x01_critical.csv, scenarios own and DFT-anch.
Reported per system / model / scenario: 2-theta of (001), (100), (101), (002), (110) for each AlB2 phase and the
splitting between the two AlB2 phases, or, for a metal-rich boride product, the strongest boride reflection and its
distance from the nearest AlB2 reflection.
Output: supplementary_data/y01_diffraction_six.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
SD = FR / "supplementary_data"
WF = FR.parent / "validation_workflow"
MB = WF / "mlip_benchmark"
sys.path.insert(0, str(MB / "scripts"))
from m07_diffraction_signatures import structure as m07_structure  # noqa: E402

MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
LAM = {"CuKa1": 1.5406, "CoKa1": 1.78897}
HKL = [(0, 0, 1), (1, 0, 0), (1, 0, 1), (0, 0, 2), (1, 1, 0)]


def structure(model: str, sid: str):
    from pymatgen.core import Structure

    p = MB / "results" / model / sid / "CONTCAR"
    if p.exists():
        return Structure.from_file(p)
    return m07_structure(model, sid)


def ac(model: str, sid: str) -> tuple[float, float]:
    """Per-formula-unit a and c of a 4x4x3 AlB2 supercell (144 atoms)."""
    st = structure(model, sid)
    return float(np.mean(st.lattice.abc[:2]) / 4), float(st.lattice.abc[2] / 3)


def two_theta(a: float, c: float, hkl, lam: float) -> float:
    h, k, l = hkl
    inv_d2 = 4 / 3 * (h * h + h * k + k * k) / a**2 + l * l / c**2
    return float(2 * np.degrees(np.arcsin(lam * np.sqrt(inv_d2) / 2)))


def lattice_model(model: str, sid: str):
    host, m = sid[:2], sid[2:]
    pts = {0.0: f"host_{host}B2", 0.25: f"host_{sid}_x0.25", 0.5: f"host_{sid}_x0.50", 0.75: f"host_{sid}_x0.75", 1.0: f"host_{m}B2"}
    xs = sorted(pts)
    A = {x: ac(model, pts[x]) for x in xs}
    d = {}
    for em in (f"{host}B2", f"{m}B2"):
        a0, c0 = ac(model, f"host_{em}")
        a8, c8 = ac(model, f"Y_{em}_k8_spread")
        d[em] = ((a8 - a0) / (8 / 96), (c8 - c0) / (8 / 96))

    def f(x: float, y: float) -> tuple[float, float]:
        a = np.interp(x, xs, [A[v][0] for v in xs])
        c = np.interp(x, xs, [A[v][1] for v in xs])
        da = (1 - x) * d[f"{host}B2"][0] + x * d[f"{m}B2"][0]
        dc = (1 - x) * d[f"{host}B2"][1] + x * d[f"{m}B2"][1]
        return a + y * da, c + y * dc

    return f


def boride_peak(model: str, name: str, lam: float) -> float | None:
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    try:
        st = structure(model, f"ref_{name}")
    except Exception:
        return None
    pat = XRDCalculator(wavelength=lam).get_pattern(st, two_theta_range=(20, 80))
    return float(pat.x[int(np.argmax(pat.y))])


def main() -> None:
    k = pd.read_csv(SD / "x01_critical.csv")
    k = k[(k.x == 0.5) & (k.T_K == 2000.0) & (k.scenario.isin(["own", "DFT-anch"]))]
    rows = []
    for r in k.itertuples():
        try:
            f = lattice_model(r.model, r.system)
        except Exception as e:  # missing structure: report, do not impute
            rows.append({"system": r.system, "model": r.model, "scenario": r.scenario, "status": f"missing structure: {e}"})
            continue
        lam = LAM["CuKa1"]
        a_h, c_h = f(0.5, 0.0)
        base = {"system": r.system, "model": r.model, "scenario": r.scenario, "status": "ok", "class_at_delta5": r.class_at_5,
                "a_homog_A": a_h, "c_homog_A": c_h}
        for hkl in HKL:
            base[f"2th_homog_{''.join(map(str, hkl))}"] = two_theta(a_h, c_h, hkl, lam)
        if r.class_at_5 == "two AlB2" and r.second_BM_at_5 == r.second_BM_at_5:
            y2 = 1 - r.second_BM_at_5 / 2
            a1, c1 = f(r.main_x_at_5, 0.0)
            a2, c2 = f(r.second_x_at_5, min(y2, 0.125))
            base.update({"phase1_x": r.main_x_at_5, "phase2_x": r.second_x_at_5, "phase2_BM": r.second_BM_at_5,
                         "c_change_phase2_vs_1_pct": 100 * (c2 / c1 - 1), "a_change_phase2_vs_1_pct": 100 * (a2 / a1 - 1)})
            for hkl in HKL:
                base[f"split_{''.join(map(str, hkl))}_deg"] = two_theta(a2, c2, hkl, lam) - two_theta(a1, c1, hkl, lam)
        elif r.class_at_5 == "metal-rich boride" and isinstance(r.second_phase_at_5, str) and r.second_phase_at_5:
            pk = boride_peak(r.model, r.second_phase_at_5, lam)
            alb = [two_theta(a_h, c_h, hkl, lam) for hkl in HKL]
            base.update({"boride": r.second_phase_at_5, "boride_strongest_2th": pk,
                         "boride_min_distance_to_AlB2_deg": None if pk is None else float(min(abs(pk - v) for v in alb))})
        rows.append(base)
    out = pd.DataFrame(rows)
    out.to_csv(SD / "y01_diffraction_six.csv", index=False)
    print(out[[c for c in ("system", "model", "scenario", "class_at_delta5", "split_001_deg", "split_002_deg", "split_100_deg",
                           "split_101_deg", "boride", "boride_strongest_2th", "boride_min_distance_to_AlB2_deg") if c in out]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
