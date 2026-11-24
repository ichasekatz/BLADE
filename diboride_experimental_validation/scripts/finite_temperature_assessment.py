"""Finite-temperature contributions to Group IV–V diboride mixing: what is known and what is required.

Uses only extracted literature values (Tang 2025 Table 2 static L0*, Table 4 QHA L0(T)) and the user's
DFT (ZrB2-TaB2 RK fit, s7_rk_fits_ZrTaB2.csv). Writes simulation_comparison/discrepancy/
s9_finite_temperature_assessment.csv and figures/fig12_finite_temperature_L0.(png|pdf).

Conventions: L0 in kJ per mol MB2 (Tang), ΔH(x=0.5) = L0/4. A symmetric regular solution has a miscibility gap
at T when L0_eff(T) > 2RT. The excess-entropy parameter S_L = −dL0/dT; ΔS_ex(0.5) = S_L/4 per mol MB2.
Electronic term: S_el = (π²/3) k_B² T g(E_F); the excess Δg(E_F) needed is reported, not assumed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

VAL = Path(__file__).resolve().parent.parent
OUT = VAL / "simulation_comparison" / "discrepancy"
R = 8.314462618  # J/(mol K)
KB_EV = 8.617333262e-5

# Tang 2025 Table 4: L0(T) = A + B T^2 + C T^3 + D T^2 ln T  (J/mol MB2)  — tang2025_dft_extraction.csv T075–T098
QHA = {"HfNb": (-1.6953e4, -1.2778e-1, -3.2267e-6, 1.6228e-2), "HfTa": (-2.6841e4, -2.4549e-1, -6.3012e-6, 3.1303e-2),
       "TiNb": (-2.8557e4, -1.0937e-1, -2.5720e-6, 1.3810e-2), "ZrNb": (-8.8268e3, -6.8045e-2, -1.6066e-6, 8.5529e-3),
       "TiTa": (-5.3691e4, -2.2006e-1, -5.5094e-6, 2.8034e-2), "ZrTa": (-2.1210e4, -2.4561e-1, -6.3207e-6, 3.1318e-2)}
STATIC = {"HfNb": -8.813, "HfTa": -14.12, "TiNb": -20.12, "ZrNb": -2.515, "TiTa": -42.21, "ZrTa": -9.306}  # kJ/mol, Table 2


def l0_qha(sid: str, t: np.ndarray) -> np.ndarray:
    a, b, c, d = QHA[sid]
    return (a + b * t**2 + c * t**3 + d * t**2 * np.log(t)) / 1000.0


def main() -> None:
    rows = []
    for sid in QHA:
        for t in (1000.0, 2000.0, 2500.0, 3000.0):
            l0 = float(l0_qha(sid, np.array([t]))[0])
            s_l = -float((l0_qha(sid, np.array([t + 1.0])) - l0_qha(sid, np.array([t - 1.0])))[0]) / 2.0 * 1000  # J/mol K
            # required: L0_eff(T) >= 2RT with static enthalpy L0* → S_L,required = −(2RT − L0*)/T
            s_req = -(2 * R * t - STATIC[sid] * 1000) / t
            rows.append({
                "system_id": sid, "T_K": t, "L0_static_kJ": STATIC[sid], "L0_QHA_kJ": l0, "gap_threshold_2RT_kJ": 2 * R * t / 1000,
                "S_L_QHA_J_per_molK": s_l, "dS_ex_0.5_QHA_kB_per_atom": s_l / 4 / R / 3,
                "S_L_required_J_per_molK": s_req, "dS_ex_0.5_required_kB_per_atom": s_req / 4 / R / 3,
                "electronic_dg_EF_required_states_per_eV_per_atom": (s_req / 4 / R / 3) / ((np.pi**2 / 3) * KB_EV * t),
                "QHA_predicts_gap": bool(l0 > 2 * R * t / 1000),
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s9_finite_temperature_assessment.csv", index=False)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = np.linspace(300, 3300, 200)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for sid in QHA:
        ax.plot(t, l0_qha(sid, t), label=f"{sid} (Tang QHA)")
    ax.plot(t, 2 * R * t / 1000, "k--", lw=1.5, label="gap threshold 2RT (regular solution)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("L$_0$(T) (kJ / mol MB$_2$)")
    ax.set_title("Tang 2025 quasi-harmonic L$_0$(T): vibrations make mixing MORE favourable\n(electronic free energy not included)", fontsize=9)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(VAL / "figures" / f"fig12_finite_temperature_L0.{ext}", dpi=200, bbox_inches="tight")
    pd.set_option("display.width", 250)
    print(df[df.T_K.isin([2500.0, 3000.0])].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
