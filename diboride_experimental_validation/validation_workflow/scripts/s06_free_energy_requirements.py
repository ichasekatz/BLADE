"""Step 6: how much extra excess free energy would reverse the predicted complete miscibility? (item 4)

Model: G(x,T) = sum_k L_k x(1-x)(1-2x)^k [Tang convention, (x_M - x_Zr) = 2x-1 with sign handled below]
                + w_add(T) x(1-x) + RT[x ln x + (1-x) ln(1-x)],   per mol MB2 (metal-site entropy only).
For each basis and T, bisection finds the smallest w_add (kJ/mol MB2, symmetric) that
  (a) opens ANY miscibility gap,
  (b) puts x = 0.5 inside the binodal (two-phase at the Ness composition),
  (c) gives a binodal at least as wide as Otani's broadening range (ZrNb 0.2-0.8, ZrTa 0.15-0.9 used as an
      upper-bound target; Otani does not show this is the equilibrium gap).
Requirements are converted to an excess entropy (if w_add = -T dS_L: dS_ex(0.5) = -w_add/(4T) per mol MB2,
then /3R -> k_B per atom) and to an electronic DOS change dg(E_F) = dS / ((pi^2/3) k_B^2 T).
Bases: ORB (RK fit to BLADE SQS energies), DFT Tang static (L0*, L1*), DFT user (ZrTa only), and
Tang static + Tang QHA vibrational L0(T) (net requirement after the computed vibrational term).

Output: results/s06_free_energy_requirements.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

WF = Path(__file__).resolve().parents[1]
VAL = WF.parent
RES = WF / "results"
R = 8.314462618e-3  # kJ/(mol K)
KB_EV = 8.617333262e-5
X = np.linspace(1e-4, 1 - 1e-4, 2001)
# Tang convention: dH = x_M x_Zr sum L_k (x_M - x_Zr)^k ; with x = x_M: (x_M - x_Zr) = 2x - 1
BASES = {
    ("ZrNb", "DFT Tang static (LDA)"): [-2.515, 2.418],
    ("ZrTa", "DFT Tang static (LDA)"): [-9.306, 4.016],
}
QHA = {"ZrNb": (-8.8268e3, -6.8045e-2, -1.6066e-6, 8.5529e-3), "ZrTa": (-2.1210e4, -2.4561e-1, -6.3207e-6, 3.1318e-2)}
TANG_STATIC_L0 = {"ZrNb": -2.515, "ZrTa": -9.306}
OTANI_RANGE = {"ZrNb": (0.2, 0.8), "ZrTa": (0.15, 0.9)}


def orb_and_user_bases() -> None:
    # user DFT and ORB RK fits for ZrTa already exist (simulation_comparison/discrepancy/s7_rk_fits_ZrTaB2.csv)
    f = pd.read_csv(VAL / "simulation_comparison" / "discrepancy" / "s7_rk_fits_ZrTaB2.csv")
    for r in f.itertuples():
        name = "DFT user (PBE, digitized)" if r.dataset == "DFT_user" else "ORB (BLADE SQS)"
        BASES[("ZrTa", name)] = [r.L0_kJ_per_mol_fu, r.L1_kJ_per_mol_fu, r.L2_meV_per_fu * 96.485 / 1000]
    # ORB ZrNb: fit BLADE SQS energies (audit) in the same convention
    a = pd.read_csv(RES / "s00_existing_structure_audit.csv")
    h = a[(a.system == "ZrNb") & (a.origin == "SQS (BLADE/ATAT)") & a.path.str.contains("system_structures/NbZrB/blade")]
    h = h.drop_duplicates("x_solute").sort_values("x_solute")
    e_fu = h.energy_eV_atom.values * 3
    xs = h.x_solute.values
    mix = (e_fu - (1 - xs) * e_fu[0] - xs * e_fu[-1]) * 96.485  # kJ/mol MB2
    inner = (xs > 0) & (xs < 1)
    A = np.vstack([xs * (1 - xs) * (2 * xs - 1) ** k for k in range(3)]).T[inner]
    BASES[("ZrNb", "ORB (BLADE SQS)")] = list(np.linalg.lstsq(A, mix[inner], rcond=None)[0])


def g_curve(lk: list[float], w_add: float, t: float) -> np.ndarray:
    hx = sum(l * X * (1 - X) * (2 * X - 1) ** k for k, l in enumerate(lk))
    return hx + w_add * X * (1 - X) + R * t * (X * np.log(X) + (1 - X) * np.log(1 - X))


def binodal(g: np.ndarray) -> tuple[float, float] | None:
    """Widest interval where the lower convex hull departs from G (two-phase region)."""
    hull: list[int] = []  # Andrew monotone chain, lower hull (X is sorted)
    for i in range(len(X)):
        while len(hull) >= 2:
            i1, i2 = hull[-2], hull[-1]
            if (X[i2] - X[i1]) * (g[i] - g[i1]) - (g[i2] - g[i1]) * (X[i] - X[i1]) <= 0:
                hull.pop()
            else:
                break
        hull.append(i)
    gaps = [(X[i], X[j]) for i, j in zip(hull[:-1], hull[1:]) if j - i > 3]
    if not gaps:
        return None
    return max(gaps, key=lambda p: p[1] - p[0])


def required(lk: list[float], t: float, test) -> float:
    lo, hi = -200.0, 400.0
    if not test(binodal(g_curve(lk, hi, t))):
        return np.nan
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if test(binodal(g_curve(lk, mid, t))) else (mid, hi)
    return hi


def main() -> None:
    orb_and_user_bases()
    rows = []
    for (sid, basis), lk in BASES.items():
        for t in (1500.0, 2000.0, 2500.0, 3000.0, 3300.0):
            variants = [(basis, list(lk), 0.0)]
            if basis.startswith("DFT Tang"):
                a, b, c, d = QHA[sid]
                l0_t = (a + b * t**2 + c * t**3 + d * t**2 * np.log(t)) / 1000
                variants.append((basis + " + Tang QHA vib", [l0_t] + list(lk[1:]), l0_t - TANG_STATIC_L0[sid]))
            for name, l, vib_shift in variants:
                tests = {"any_gap": lambda b: b is not None,
                         "x0.5_two_phase": lambda b: b is not None and b[0] < 0.5 < b[1],
                         "otani_range": lambda b, r=OTANI_RANGE[sid]: b is not None and b[0] <= r[0] and b[1] >= r[1]}
                for crit, fn in tests.items():
                    w = required(l, t, fn)
                    ds = -w / (4 * t) * 1000  # J/(mol MB2 K) at x = 0.5
                    ds_kb_atom = ds / (3 * R * 1000)
                    rows.append({"system": sid, "basis": name, "T_K": t, "criterion": crit,
                                 "w_add_required_kJ_per_mol_MB2": w, "dH_equiv_at_x0.5_meV_atom": w / 4 / 3 / 96.485 * 1000,
                                 "dS_ex_required_x0.5_kB_per_atom": ds_kb_atom,
                                 "dg_EF_required_states_per_eV_atom": ds_kb_atom / ((np.pi**2 / 3) * KB_EV * t),
                                 "vibrational_shift_already_included_kJ": vib_shift})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "s06_free_energy_requirements.csv", index=False)
    pd.set_option("display.width", 220)
    print(df[df.T_K.isin([2000.0, 3000.0])].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
