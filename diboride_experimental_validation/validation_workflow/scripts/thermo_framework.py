"""Stand-alone thermodynamic framework for (Zr,M)-B (M = Nb or Ta) that can include boron non-stoichiometry,
competing solid phases and (optionally) a liquid, without touching BLADE outputs.

Phases
------
AlB2  compound-energy formalism (Zr,M)1(B,Va)2, per formula unit (1 metal site + 2 B sites):
      G = sum_ij y_i y'_j G_ij
          + y_Zr y_M [ y'_B L_B(x) + y'_Va L_Va(x) ] + y_Zr y_M y'_B y'_Va L_rec(T)
          + kT [ sum_i y_i ln y_i + 2 sum_j y'_j ln y'_j ] + G_extra(x, T)
      G_iB  = energy of iB2 (eV per formula unit)
      G_iVa = G_iB - 2 E_B + 2 dE_v(i), dE_v(i) = E(host - B) - E(host) + E_B (= E_f(V_B) at mu_B = E_B)
      L_B(x)= sum_k L_k (2x - 1)^k  (Tang convention, x = y_M); L_Va = L_B (assumption, flagged)
      W_vv(x) y'^2 = [(1-x) W_Zr + x W_M] y'^2 is a mean-field vacancy-vacancy term fitted to multi-vacancy cells
      (it leaves the dilute single-vacancy energy unchanged).
      L_rec(T) is fitted so the dilute vacancy energy at x = 0.5 matches the local-environment average
      E_v,eff(T) = -kT ln sum_n P(n | x) exp(-E_v(n)/kT)  (binomial P for a random alloy).
Line compounds: every competitor relaxed in s03 (metals, B polymorphs, MB, M3B4, M5B6, M2B3, M3B2, M2B, ZrB12,
      ordered ZrMB4 ...) with G = E (0 K) + optional G_extra(T).
Liquid (optional, pseudo-binary ZrB2-MB2 only, B/M = 2):
      G_L(x) = sum_i x_i [G_iB + dS_fus,i (T_m,i - T)] + kT[x ln x + (1-x) ln(1-x)] + x(1-x) W_L
      T_m from Rudy & Windisch; dS_fus and W_L are INPUTS (no measured values available) to be scanned.

Equilibrium: Gibbs-energy minimisation as a linear programme over a discretised AlB2 (x, y_Va) grid + compounds
(+ liquid grid), element balances for Zr, M, B; HiGHS. Distinct AlB2 phases are identified by clustering the
active grid points (|dx| > 0.02 or |dy| > 0.01 starts a new phase).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import pandas as pd
from scipy.optimize import linprog

KB = 8.617333262e-5  # eV/K


@dataclass
class AlB2Params:
    g_zrb: float  # eV per formula unit ZrB2
    g_mb: float  # eV per formula unit MB2
    dEv_zr: float  # eV, B-vacancy energy in ZrB2 relative to elemental B
    dEv_m: float
    e_b: float  # eV/atom elemental B reference
    rk: list[float]  # eV per formula unit, Tang convention in (2x - 1)
    ev_env: Callable[[int], float] | None = None  # E_v(n_solute among 6 nn) relative to E_B, eV
    l_va_equals_l_b: bool = True
    w_vv_zr: float = 0.0  # eV/fu, B-vacancy - B-vacancy interaction W y'^2 in ZrB2 (fitted from Y series)
    w_vv_m: float = 0.0
    g_extra: Callable[[np.ndarray, float], np.ndarray] | None = None  # (x, T) -> eV/fu

    def g_iva(self) -> tuple[float, float]:
        return self.g_zrb - 2 * self.e_b + 2 * self.dEv_zr, self.g_mb - 2 * self.e_b + 2 * self.dEv_m

    def l_rec(self, t: float) -> float:
        """Reciprocal parameter so that the CEF dilute vacancy energy at x = 0.5 equals E_v,eff(0.5, T)."""
        if self.ev_env is None:
            return 0.0
        from math import comb

        n = np.arange(7)
        p = np.array([comb(6, k) * 0.5**6 for k in n])
        ev = np.array([self.ev_env(int(k)) for k in n])
        tt = max(t, 1.0)
        ev_eff = -KB * tt * np.log(np.sum(p * np.exp(-(ev - ev.min()) / (KB * tt)))) + ev.min()
        linear = 0.5 * (self.dEv_zr + self.dEv_m)
        # dG/dy at y->0, per B removed: lin + x(1-x)(L_rec)/2 at x = 0.5 -> (1/4) L_rec / 2
        return 8.0 * (ev_eff - linear)

    def gibbs(self, x: np.ndarray, y: np.ndarray, t: float) -> np.ndarray:
        gzv, gmv = self.g_iva()
        lb = sum(l * (2 * x - 1) ** k for k, l in enumerate(self.rk))
        lv = lb if self.l_va_equals_l_b else 0.0
        g = ((1 - x) * (1 - y) * self.g_zrb + x * (1 - y) * self.g_mb + (1 - x) * y * gzv + x * y * gmv
             + x * (1 - x) * ((1 - y) * lb + y * lv) + x * (1 - x) * y * (1 - y) * self.l_rec(t))
        with np.errstate(divide="ignore", invalid="ignore"):
            s_m = np.where((x > 0) & (x < 1), x * np.log(x) + (1 - x) * np.log(1 - x), 0.0)
            s_b = np.where((y > 0) & (y < 1), y * np.log(y) + (1 - y) * np.log(1 - y), 0.0)
        g = g + KB * t * (s_m + 2 * s_b) + y**2 * ((1 - x) * self.w_vv_zr + x * self.w_vv_m)
        if self.g_extra is not None:
            g = g + self.g_extra(x, t)
        return g


@dataclass
class LiquidParams:
    tm_zrb2: float  # K
    tm_mb2: float
    ds_fus_zrb2: float  # eV/K per formula unit
    ds_fus_mb2: float
    w_l: float = 0.0  # eV/fu


@dataclass
class System:
    metal: str  # "Nb" or "Ta"
    alb2: AlB2Params
    compounds: pd.DataFrame  # columns: name, n_Zr, n_M, n_B, G (eV per formula as given)
    liquid: LiquidParams | None = None
    x_grid: np.ndarray = field(default_factory=lambda: np.linspace(0, 1, 101))
    y_grid: np.ndarray = field(default_factory=lambda: np.r_[0.0, np.geomspace(1e-4, 0.35, 48)])

    def columns(self, t: float) -> tuple[np.ndarray, np.ndarray, list[dict]]:
        """Return composition matrix A (3 x N: Zr, M, B), Gibbs vector g, and labels."""
        cols, g, lab = [], [], []
        xx, yy = np.meshgrid(self.x_grid, self.y_grid, indexing="ij")
        xx, yy = xx.ravel(), yy.ravel()
        gg = self.alb2.gibbs(xx, yy, t)
        for x, y, gv in zip(xx, yy, gg):
            cols.append([1 - x, x, 2 * (1 - y)])
            g.append(gv)
            lab.append({"phase": "AlB2", "x": x, "y_Va": y, "n_metal": 1.0})
        for r in self.compounds.itertuples():
            cols.append([r.n_Zr, r.n_M, r.n_B])
            g.append(r.G)
            lab.append({"phase": r.name, "x": r.n_M / (r.n_Zr + r.n_M) if r.n_Zr + r.n_M else np.nan, "y_Va": np.nan, "n_metal": float(r.n_Zr + r.n_M)})
        if self.liquid is not None:
            lq = self.liquid
            for x in self.x_grid:
                gl = ((1 - x) * (self.alb2.g_zrb + lq.ds_fus_zrb2 * (lq.tm_zrb2 - t)) + x * (self.alb2.g_mb + lq.ds_fus_mb2 * (lq.tm_mb2 - t))
                      + x * (1 - x) * lq.w_l)
                if 0 < x < 1:
                    gl += KB * t * (x * np.log(x) + (1 - x) * np.log(1 - x))
                cols.append([1 - x, x, 2.0])
                g.append(gl)
                lab.append({"phase": "LIQUID", "x": x, "y_Va": np.nan, "n_metal": 1.0})
        return np.array(cols).T, np.array(g), lab

    def equilibrium(self, x_overall: float, b_per_metal: float, t: float, method: str = "highs") -> pd.DataFrame:
        a, g, lab = self.columns(t)
        b = np.array([1 - x_overall, x_overall, b_per_metal])
        res = linprog(g, A_eq=a, b_eq=b, bounds=(0, None), method=method)
        if not res.success:
            return pd.DataFrame([{"status": res.message}])
        act = [(i, v) for i, v in enumerate(res.x) if v > 1e-7]
        rows = [{**lab[i], "amount_fu": v} for i, v in act]
        df = pd.DataFrame(rows)
        # merge adjacent AlB2 grid points into phases
        out = []
        for ph, grp in df.groupby("phase", sort=False):
            if ph not in ("AlB2", "LIQUID"):
                out.append({"phase": ph, "x_M": grp.x.iloc[0], "y_Va": np.nan, "amount_fu": grp.amount_fu.sum(),
                            "metal_fraction": float((grp.amount_fu * grp.n_metal).sum())})
                continue
            grp = grp.sort_values("x")
            clusters, cur = [], [grp.iloc[0]]
            for _, r in list(grp.iterrows())[1:]:
                prev = cur[-1]
                if abs(r.x - prev.x) > 0.02 or (ph == "AlB2" and abs(r.y_Va - prev.y_Va) > 0.01):
                    clusters.append(cur)
                    cur = [r]
                else:
                    cur.append(r)
            clusters.append(cur)
            for c in clusters:
                c = pd.DataFrame(c)
                w = c.amount_fu / c.amount_fu.sum()
                out.append({"phase": ph, "x_M": float((w * c.x).sum()), "y_Va": float((w * c.y_Va).sum()) if ph == "AlB2" else np.nan,
                            "amount_fu": float(c.amount_fu.sum()), "metal_fraction": float(c.amount_fu.sum())})
        o = pd.DataFrame(out)
        o["G_total_eV"] = float(res.fun)
        return o
