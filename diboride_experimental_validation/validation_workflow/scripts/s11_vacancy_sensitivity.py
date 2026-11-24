"""Step 11: sensitivity of the boron-vacancy mechanism (H1) to the ORB vacancy energetics (no new calculations).

Uses the framework (thermo_framework.py) with parameters built exactly as in s08, then perturbs:
  d_m    shift of the B-vacancy energy in the Group-V diboride (eV)      [ORB TaB2/NbB2 values are suspect]
  d_zr   shift of the B-vacancy energy in ZrB2 (eV)
  The local-environment energies E_v(n) at x = 0.5 shift consistently: E_v(n) += d_zr + (d_m - d_zr) n / 6.
  slope  multiplier on the local-environment dependence (0 = no local preference beyond linear CEF)
  W      multiplier on the vacancy-vacancy interaction W (both end members)
  comp   shift of all metal-rich boride energies (eV/atom), testing H4 competition

Outcome classes at fixed overall composition (closed system: x_overall, B/M):
  HOMOGENEOUS          one AlB2 phase (vacancy concentration simply adjusts), maybe with free B
  AlB2_GAP             two AlB2 phases with |dx| > 0.05 (metal partitioning coupled to B content)  <- H1
  AlB2+BORIDE          one AlB2 + a non-AlB2 metal-rich boride (ordered B-deficient structure)     <- H4
  AlB2_GAP+BORIDE      both
Also:
  * critical d_m: smallest upward shift of the Group-V vacancy energy that removes AlB2_GAP;
  * open-system equilibria vs mu_B (grand potential per metal atom, lower convex hull in x);
  * spinodal map: smallest eigenvalue of the Hessian of G(x, y) per metal atom (local instability vs nucleation).

Outputs: results/s11_sensitivity_*.csv, figures/fig_s11_*.png
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s08_framework_equilibria as s08  # noqa: E402
from thermo_framework import KB, System  # noqa: E402

WF = Path(__file__).resolve().parents[1]
RES = WF / "results"
FIG = WF / "figures"
SOL = {"ZrNb": "Nb", "ZrTa": "Ta"}
METAL_RICH = ("B (", "Zr (", "Nb (", "Ta (")  # names not counted as borides


def perturbed(base: System, d_m: float = 0.0, d_zr: float = 0.0, slope: float = 1.0, w: float = 1.0, comp: float = 0.0) -> System:
    a = base.alb2
    env0 = a.ev_env
    lin0 = lambda n: a.dEv_zr + (a.dEv_m - a.dEv_zr) * n / 6.0

    def env(n: int, env0=env0, lin0=lin0) -> float:
        e = env0(n) if env0 else lin0(n)
        # keep the mean level, scale only the deviation of the n-dependence from the end-member line
        e = lin0(n) + slope * (e - lin0(n)) if slope != 1.0 else e
        return e + d_zr + (d_m - d_zr) * n / 6.0

    new = replace(a, dEv_m=a.dEv_m + d_m, dEv_zr=a.dEv_zr + d_zr, ev_env=env, w_vv_zr=a.w_vv_zr * w, w_vv_m=a.w_vv_m * w)
    c = base.compounds.copy()
    if comp:
        boride = (c.n_B > 0) & (c.n_Zr + c.n_M > 0)
        c.loc[boride, "G"] += comp * (c.loc[boride, "n_Zr"] + c.loc[boride, "n_M"] + c.loc[boride, "n_B"])
    return System(metal=base.metal, alb2=new, compounds=c, x_grid=np.linspace(0, 1, 101), y_grid=np.r_[0.0, np.geomspace(1e-5, 0.35, 70)])


def classify(ph: pd.DataFrame) -> tuple[str, float, float]:
    al = ph[ph.phase == "AlB2"]
    bor = ph[(ph.phase != "AlB2") & ~ph.phase.str.startswith(METAL_RICH) & (ph.x_M.notna())]
    dx = float(al.x_M.max() - al.x_M.min()) if len(al) > 1 else 0.0
    gap = dx > 0.05
    tag = ("AlB2_GAP+BORIDE" if gap and len(bor) else "AlB2_GAP" if gap else "AlB2+BORIDE" if len(bor) else "HOMOGENEOUS")
    ymax = float(al.y_Va.max()) if len(al) else np.nan
    return tag, dx, ymax


def critical_shift(base: System, t: float, bm: float, x: float = 0.5) -> float:
    """Smallest d_m in [0, 2] eV that removes the two-AlB2 field (bisection on the class)."""
    has_gap = lambda d: classify(perturbed(base, d_m=d).equilibrium(x, bm, t))[0].startswith("AlB2_GAP")
    if not has_gap(0.0):
        return 0.0
    if has_gap(2.0):
        return np.inf
    lo, hi = 0.0, 2.0
    for _ in range(12):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if has_gap(mid) else (lo, mid)
    return hi


def open_system(base: System, mu_b: float, t: float) -> pd.DataFrame:
    """Grand potential per metal atom, Omega(x) = min_y G(x,y) - mu_B * 2(1-y); lower hull over x incl. compounds."""
    a = base.alb2
    xs = np.linspace(0, 1, 201)
    ys = np.r_[0.0, np.geomspace(1e-6, 0.35, 400)]
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    om = a.gibbs(xx, yy, t) - mu_b * 2 * (1 - yy)
    k = om.argmin(axis=1)
    pts = [(x, om[i, k[i]], "AlB2", ys[k[i]]) for i, x in enumerate(xs)]
    for r in base.compounds.itertuples():
        nmet = r.n_Zr + r.n_M
        if nmet == 0:
            continue
        pts.append((r.n_M / nmet, (r.G - mu_b * r.n_B) / nmet, r.name, np.nan))
    pts.sort(key=lambda p: (p[0], p[1]))
    hull: list = []
    for p in pts:
        while len(hull) >= 2 and (hull[-1][0] - hull[-2][0]) * (p[1] - hull[-2][1]) - (hull[-1][1] - hull[-2][1]) * (p[0] - hull[-2][0]) <= 0:
            hull.pop()
        if hull and abs(hull[-1][0] - p[0]) < 1e-12:
            if p[1] < hull[-1][1]:
                hull[-1] = p
            continue
        hull.append(p)
    # phases in equilibrium at x = 0.5: the hull segment spanning 0.5
    for i in range(len(hull) - 1):
        if hull[i][0] <= 0.5 <= hull[i + 1][0]:
            seg = [hull[i], hull[i + 1]] if hull[i + 1][0] - hull[i][0] > 0.011 else [hull[i]]
            return pd.DataFrame([{"phase": s[2], "x_M": s[0], "y_Va": s[3]} for s in seg])
    return pd.DataFrame()


def hessian_min_eig(base: System, t: float, x: float, y: float, h: float = 2e-3) -> float:
    g = lambda xx, yy: float(base.alb2.gibbs(np.array([xx]), np.array([yy]), t)[0])
    gxx = (g(x + h, y) - 2 * g(x, y) + g(x - h, y)) / h**2
    gyy = (g(x, y + h) - 2 * g(x, y) + g(x, y - h)) / h**2
    gxy = (g(x + h, y + h) - g(x + h, y - h) - g(x - h, y + h) + g(x - h, y - h)) / (4 * h * h)
    return float(np.linalg.eigvalsh(np.array([[gxx, gxy], [gxy, gyy]]))[0])


def main() -> None:
    rows_map, rows_crit, rows_oat, rows_open, rows_sp = [], [], [], [], []
    for sid in SOL:
        for pset in ("ORB", "DFT"):
            base = s08.build(sid, pset)
            label = "ORB" if pset == "ORB" else "ORB-vac+DFT"
            # 1. 2-D map d_m x d_zr at T = 2000 K, B/M = 1.95, x = 0.5
            for d_m in np.round(np.arange(-0.2, 1.61, 0.1), 2):
                for d_zr in (-0.5, -0.25, 0.0, 0.25, 0.5):
                    tag, dx, ym = classify(perturbed(base, d_m=d_m, d_zr=d_zr).equilibrium(0.5, 1.95, 2000.0))
                    rows_map.append({"system": sid, "parameter_set": label, "d_m_eV": d_m, "d_zr_eV": d_zr, "class": tag, "AlB2_dx": dx, "y_max": ym})
            # 2. critical d_m vs T and B/M
            for t in (1500.0, 2000.0, 2500.0, 3000.0):
                for bm in (2.0, 1.98, 1.95, 1.90):
                    rows_crit.append({"system": sid, "parameter_set": label, "T_K": t, "B_per_metal": bm,
                                      "critical_d_m_eV": critical_shift(base, t, bm)})
            # 3. one-at-a-time sensitivities at T = 2000 K, B/M = 1.95
            for knob, vals in (("slope", (0.0, 0.5, 1.0, 1.5)), ("w", (0.5, 1.0, 2.0, 4.0)), ("comp", (-0.05, -0.02, 0.0, 0.02, 0.05)),
                               ("d_zr", (-0.5, 0.0, 0.5)), ("d_m", (0.0, 0.3, 0.6, 0.9))):
                for v in vals:
                    tag, dx, ym = classify(perturbed(base, **{knob: v}).equilibrium(0.5, 1.95, 2000.0))
                    rows_oat.append({"system": sid, "parameter_set": label, "knob": knob, "value": v, "class": tag, "AlB2_dx": dx, "y_max": ym})
            # 4. open system vs mu_B (between B-rich and B-poor limits of s04 at x = 0.5)
            lim = pd.read_csv(RES / "s04_mu_B_limits.csv")
            l5 = lim[(lim.system == sid) & np.isclose(lim.x_solute, 0.5) & (lim.T_K == 0)].iloc[0]
            for d_m in (0.0, 0.5):
                b = perturbed(base, d_m=d_m)
                for t in (2000.0, 2500.0):
                    for mu in np.linspace(l5.mu_B_poor - 0.1, l5.mu_B_rich, 12):
                        ph = open_system(b, mu, t)
                        rows_open.append({"system": sid, "parameter_set": label, "d_m_eV": d_m, "T_K": t, "mu_B": mu,
                                          "mu_B_minus_Brich": mu - l5.mu_B_rich,
                                          "phases": ";".join(f"{p}[x={xm:.2f},y={yv:.3f}]" for p, xm, yv in zip(ph.phase, ph.x_M, ph.y_Va))})
            # 5. spinodal map (Hessian) at T = 2000 K
            for x in np.round(np.arange(0.05, 0.96, 0.05), 2):
                for y in (0.0005, 0.005, 0.02, 0.05, 0.1):
                    rows_sp.append({"system": sid, "parameter_set": label, "T_K": 2000.0, "x": x, "y_Va": y,
                                    "min_eig_eV": hessian_min_eig(base, 2000.0, x, y)})
    pd.DataFrame(rows_map).to_csv(RES / "s11_sensitivity_map_dm_dzr.csv", index=False)
    crit = pd.DataFrame(rows_crit)
    crit.to_csv(RES / "s11_sensitivity_critical_shift.csv", index=False)
    pd.DataFrame(rows_oat).to_csv(RES / "s11_sensitivity_one_at_a_time.csv", index=False)
    pd.DataFrame(rows_open).to_csv(RES / "s11_open_system_muB.csv", index=False)
    pd.DataFrame(rows_sp).to_csv(RES / "s11_spinodal_hessian.csv", index=False)
    pd.set_option("display.width", 220)
    print(crit.pivot_table(index=["system", "parameter_set", "B_per_metal"], columns="T_K", values="critical_d_m_eV").round(2).to_string())
    print(pd.DataFrame(rows_oat).to_string(index=False))
    plot()


def plot() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m = pd.read_csv(RES / "s11_sensitivity_map_dm_dzr.csv")
    cls = {"HOMOGENEOUS": 0, "AlB2+BORIDE": 1, "AlB2_GAP": 2, "AlB2_GAP+BORIDE": 3}
    cmap = matplotlib.colors.ListedColormap(["#d9d9d9", "#9ecae1", "#fb6a4a", "#a50f15"])
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), sharex=True, sharey=True)
    for ax, ((sid, ps), g) in zip(axes.ravel(), m.groupby(["system", "parameter_set"])):
        piv = g.assign(c=g["class"].map(cls)).pivot(index="d_zr_eV", columns="d_m_eV", values="c")
        ax.imshow(piv.values, origin="lower", aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5,
                  extent=[piv.columns.min() - 0.05, piv.columns.max() + 0.05, piv.index.min() - 0.125, piv.index.max() + 0.125])
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.axhline(0, color="k", lw=0.8, ls="--")
        ax.set_title(f"{sid}, {ps}: x=0.5, B/M=1.95, 2000 K", fontsize=9)
    for ax in axes[1]:
        ax.set_xlabel("shift of Group-V B-vacancy energy d_m (eV)")
    for ax in axes[:, 0]:
        ax.set_ylabel("shift of ZrB$_2$ B-vacancy energy d_Zr (eV)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=cmap(i)) for i in range(4)]
    fig.legend(handles, list(cls), loc="lower center", ncol=4, fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "fig_s11_sensitivity_map.png", dpi=200)

    c = pd.read_csv(RES / "s11_sensitivity_critical_shift.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, (sid, g) in zip(axes, c.groupby("system")):
        for (ps, bm), h in g.groupby(["parameter_set", "B_per_metal"]):
            ax.plot(h.T_K, h.critical_d_m_eV.replace(np.inf, 2.0), "o-" if ps == "ORB" else "s--", label=f"{ps}, B/M={bm}")
        ax.set_title(f"{sid}: upward shift of the Group-V V_B energy that removes the AlB$_2$ gap", fontsize=8)
        ax.set_xlabel("T (K)")
    axes[0].set_ylabel("critical d_m (eV)  (2.0 = not removed within 2 eV)")
    axes[1].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "fig_s11_critical_shift.png", dpi=200)


if __name__ == "__main__":
    main()
