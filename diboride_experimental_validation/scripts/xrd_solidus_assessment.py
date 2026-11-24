"""Quantitative use of the digitized Otani 2009 Fig. 1 and Eckert 1966 Figs. 41-46.

1. Otani (10-11) peak broadening -> equivalent metal-composition spread dx. The excess FWHM over the end-member
   width is divided by |d(2theta)/dx| from Otani's own lattice constants. ZrB2-TiB2 (known gap) and ZrB2-HfB2
   (complete solubility) are the controls.
2. Eckert lattice parameters vs Otani at matched x: are the two data sets consistent?
3. XRD resolving power: separation of (10-11), (002), (110) peaks for two phases x0 +/- dx/2 (Eckert/Otani slopes).
4. Eckert melting: incipient melting/collapse vs the linear line between literature end members, and the offset
   Eckert's own method shows on pure NbB2 (Fig. 46). Is a eutectic-like minimum resolvable?

Assumptions (stated in output): Cu K-alpha1 radiation (lambda = 1.5406 A; Otani does not state it), size/strain
broadening equal to the end members, and a uniform composition distribution of width dx (FWHM ~ |d2theta/dx| dx).
Writes simulation_comparison/discrepancy/s11_*.csv, s12_*.csv, s13_*.csv and figures/fig13_xrd_solidus.(png|pdf).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

VAL = Path(__file__).resolve().parent.parent
DIG = VAL / "raw_data" / "digitized_figures"
OUT = VAL / "simulation_comparison" / "discrepancy"
LAMBDA_CU_KA1_A = 1.5406


def two_theta_deg(a_a: np.ndarray, c_a: np.ndarray, hkl: tuple[int, int, int]) -> np.ndarray:
    """Bragg angle 2theta (deg) of a hexagonal reflection for lattice constants a, c in Angstrom."""
    h, k, l = hkl
    inv_d2 = 4.0 / 3.0 * (h * h + h * k + k * k) / a_a**2 + l * l / c_a**2
    return 2.0 * np.degrees(np.arcsin(LAMBDA_CU_KA1_A * np.sqrt(inv_d2) / 2.0))


def otani_series(d: pd.DataFrame, system: str, series: str) -> pd.DataFrame:
    g = d[(d.system == system) & (d.series == series)][["x_solute", "y_value"]].dropna()
    return g.groupby("x_solute", as_index=False).y_value.mean().sort_values("x_solute")


def otani_composition_spread(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for system in ("ZrB2-NbB2", "ZrB2-TaB2", "ZrB2-HfB2", "ZrB2-TiB2"):
        a = otani_series(d, system, "a lattice constant")
        c = otani_series(d, system, "c lattice constant")
        # end-member lattice constants and linear (Vegard) slopes from Otani's own end points
        a0, a1 = a.y_value.iloc[0] * 10, a.y_value.iloc[-1] * 10
        c0, c1 = c.y_value.iloc[0] * 10, c.y_value.iloc[-1] * 10
        for label in ("(10-11) FWHM arc-melted", "(10-11) FWHM water-cooled"):
            w = otani_series(d, system, label)
            if w.empty:
                continue
            w_end = 0.5 * (w.y_value.iloc[0] + w.y_value.iloc[-1])  # instrumental + size/strain baseline
            for x, fwhm in zip(w.x_solute, w.y_value):
                if x in (0.0, 1.0):
                    continue
                aa = a0 + (a1 - a0) * np.array([x - 0.005, x + 0.005])
                cc = c0 + (c1 - c0) * np.array([x - 0.005, x + 0.005])
                slope = abs(float(np.diff(two_theta_deg(aa, cc, (1, 0, 1)))[0])) / 0.01  # deg per unit x
                excess_lin = max(fwhm - w_end, 0.0)  # Lorentzian-like (linear) subtraction
                excess_quad = float(np.sqrt(max(fwhm**2 - w_end**2, 0.0)))  # Gaussian-like (quadrature)
                rows.append({
                    "system": system, "cooling": "slow (Cu block)" if "arc-melted" in label else "water-cooled",
                    "x_solute": x, "fwhm_deg": fwhm, "end_member_fwhm_deg": round(w_end, 3),
                    "d2theta_dx_deg_101": round(slope, 3),
                    "dx_equiv_linear": round(excess_lin / slope, 3), "dx_equiv_quadrature": round(excess_quad / slope, 3),
                    "fwhm_uncertainty_deg": 0.01,
                })
    return pd.DataFrame(rows)


def eckert_vs_otani(d_ot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fig, system in ((41, "ZrB2-NbB2"), (43, "ZrB2-TaB2")):
        e = pd.read_csv(DIG / f"eckert1966_fig{fig}.csv")
        for par in ("a", "c"):
            pts = e[e.series == par]
            ot = otani_series(d_ot, system, f"{par} lattice constant")
            for x, y in zip(pts.x_solute, pts.y_value):
                y_ot = float(np.interp(x, ot.x_solute, ot.y_value)) * 10
                rows.append({"system": system, "parameter": par, "x_solute": round(x, 3), "eckert_A": round(y, 4),
                             "otani_interp_A": round(y_ot, 4), "eckert_minus_otani_A": round(y - y_ot, 4)})
    return pd.DataFrame(rows)


def peak_separation(d_ot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for system in ("ZrB2-NbB2", "ZrB2-TaB2"):
        a = otani_series(d_ot, system, "a lattice constant")
        c = otani_series(d_ot, system, "c lattice constant")
        a0, a1, c0, c1 = a.y_value.iloc[0] * 10, a.y_value.iloc[-1] * 10, c.y_value.iloc[0] * 10, c.y_value.iloc[-1] * 10
        for dx in (0.05, 0.1, 0.2, 0.4, 0.6):
            xs = np.array([0.5 - dx / 2, 0.5 + dx / 2])
            aa, cc = a0 + (a1 - a0) * xs, c0 + (c1 - c0) * xs
            row = {"system": system, "x0": 0.5, "dx_between_phases": dx}
            for hkl in ((1, 0, 1), (0, 0, 2), (1, 1, 0), (1, 1, 2)):
                row[f"sep_2theta_{''.join(map(str, hkl))}_deg"] = round(abs(float(np.diff(two_theta_deg(aa, cc, hkl))[0])), 3)
            rows.append(row)
    return pd.DataFrame(rows)


def eckert_melting() -> pd.DataFrame:
    rows = []
    for fig, system in ((42, "ZrB2-NbB2"), (44, "ZrB2-TaB2"), (46, "HfB2-NbB2")):
        e = pd.read_csv(DIG / f"eckert1966_fig{fig}.csv")
        ends = e[e.series == "Rudy and Windisch"].set_index("x_solute").y_value
        t0, t1 = float(ends.loc[0.0]), float(ends.loc[1.0])
        for series in ("incipient melting", "collapse"):
            for x, t in zip(e[e.series == series].x_solute, e[e.series == series].y_value):
                lin = t0 + (t1 - t0) * x
                rows.append({"system": system, "series": series, "x_solute": round(x, 3), "T_C": t,
                             "linear_endmember_T_C": round(lin, 0), "T_minus_linear_C": round(t - lin, 0),
                             "T_minus_lower_endmember_C": round(t - min(t0, t1), 0),
                             "T_uncertainty_C": 10})
    return pd.DataFrame(rows)


def main() -> None:
    d_ot = pd.read_csv(DIG / "otani2009_fig1.csv")
    s11 = otani_composition_spread(d_ot)
    s11.to_csv(OUT / "s11_otani_xrd_composition_spread.csv", index=False)
    s12 = pd.concat([eckert_vs_otani(d_ot), peak_separation(d_ot)], keys=["eckert_vs_otani", "peak_separation"],
                    names=["table"]).reset_index(level=0)
    s12.to_csv(OUT / "s12_lattice_and_peak_separation.csv", index=False)
    s13 = eckert_melting()
    s13.to_csv(OUT / "s13_eckert_melting_vs_endmembers.csv", index=False)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    style = {"ZrB2-NbB2": "C0", "ZrB2-TaB2": "C3", "ZrB2-HfB2": "C2", "ZrB2-TiB2": "k"}
    for (system, cooling), g in s11.groupby(["system", "cooling"]):
        ls = "-" if cooling.startswith("slow") else "--"
        ax1.errorbar(g.x_solute, g.dx_equiv_linear, yerr=[np.zeros(len(g)), g.dx_equiv_quadrature - g.dx_equiv_linear], fmt="o" + ls,
                     color=style[system], ms=4, capsize=2, label=f"{system} ({cooling})")
    ax1.set_xlabel("x in Zr$_{1-x}$M$_x$B$_2$")
    ax1.set_ylabel("equivalent composition spread Δx")
    ax1.set_title("Otani (10-11) excess broadening as composition spread\nmarker: linear baseline subtraction; bar top: quadrature subtraction",
                  fontsize=9)
    ax1.legend(fontsize=7)
    for system, g in s13.groupby("system"):
        for series, m in (("incipient melting", "o"), ("collapse", "^")):
            h = g[g.series == series]
            ax2.plot(h.x_solute, h.T_minus_linear_C, m, color={"ZrB2-NbB2": "C0", "ZrB2-TaB2": "C3", "HfB2-NbB2": "C1"}[system],
                     mfc="none" if series == "collapse" else None, label=f"{system} {series}")
    ax2.axhline(0, color="k", lw=0.8)
    ax2.axhspan(-109, -57, color="grey", alpha=0.25, label="Eckert's own offset on pure NbB$_2$ (Fig. 46)")
    ax2.set_xlabel("x (Group V metal fraction)")
    ax2.set_ylabel("T − linear(end members) (°C)")
    ax2.set_title("Eckert 1966 melting points relative to literature end members", fontsize=9)
    ax2.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(VAL / "figures" / f"fig13_xrd_solidus.{ext}", dpi=200, bbox_inches="tight")

    pd.set_option("display.width", 250)
    print(s11.to_string(index=False))
    print(s12.to_string(index=False))
    print(s13.groupby(["system", "series"]).T_minus_linear_C.agg(["min", "max", "mean"]).to_string())


if __name__ == "__main__":
    main()
