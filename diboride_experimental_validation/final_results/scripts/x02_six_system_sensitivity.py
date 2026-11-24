"""Six-system sensitivity: how much must one energetic term change to reverse the predicted assemblage?

Reference: x = 0.5, T = 2000 K, B/M = 1.95 (delta = 2.5 %), AlB2 vacancy fraction <= 0.125, each model's own energetics
(and the DFT-anchored TaB2 V_B for the Ta systems). Classes: homogeneous | two AlB2 | metal-rich boride |
vacancy-ordered MB2-x | ordered metal. One term varied at a time:
  d_m    Group-V B-vacancy energy shift (eV)              -1.5 ... +1.5
  d_h    host (Ti/Zr/Hf)B2 B-vacancy energy shift (eV)    -1.5 ... +1.5
  slope  multiplier on the local-environment coupling      0 ... 2
  w      multiplier on the vacancy-vacancy interaction     0.1 ... 10
  comp   shift of all metal-rich borides (eV/atom; ordered cells not shifted) -0.10 ... +0.10
The smallest change that flips the class is reported (None = no flip in range).
Output: supplementary_data/x02_flip_thresholds.csv
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FR / "scripts"))
from x01_six_system_thermo import MODELS, build, energies, solve, summarise  # noqa: E402

X0, T0, BM0 = 0.5, 2000.0, 1.95


def perturb(sysm, d_m=0.0, d_h=0.0, slope=1.0, w=1.0, comp=0.0):
    a = sysm.alb2
    env0 = a.ev_env
    lin = lambda n: a.dEv_zr + (a.dEv_m - a.dEv_zr) * n / 6.0

    def env(n, env0=env0, lin=lin):
        e = env0(n)
        e = lin(n) + slope * (e - lin(n))
        return e + d_h + (d_m - d_h) * n / 6.0

    new = replace(a, dEv_m=a.dEv_m + d_m, dEv_zr=a.dEv_zr + d_h, ev_env=env, w_vv_zr=a.w_vv_zr * w, w_vv_m=a.w_vv_m * w)
    c = sysm.compounds.copy()
    if comp:
        mrb = (c.kind == "reference") & (c.n_B > 0) & (c.n_Zr + c.n_M > 0)
        c.loc[mrb, "G"] += comp * (c.loc[mrb, "n_Zr"] + c.loc[mrb, "n_M"] + c.loc[mrb, "n_B"])
    s = type(sysm)(metal=sysm.metal, alb2=new, compounds=c, x_grid=sysm.x_grid, y_grid=sysm.y_grid)
    return s


def cls(s) -> str:
    return summarise(solve(s, X0, BM0, T0), s.compounds)["class"]


def main() -> None:
    E, nat = energies()
    rows = []
    grids = {"d_m": np.round(np.arange(-1.5, 1.51, 0.05), 3), "d_h": np.round(np.arange(-1.5, 1.51, 0.1), 3),
             "slope": np.round(np.arange(0.0, 2.01, 0.1), 2), "w": np.geomspace(0.1, 10, 25), "comp": np.round(np.arange(-0.10, 0.101, 0.005), 3)}
    for sid in ["TiNb", "TiTa", "ZrNb", "ZrTa", "HfNb", "HfTa"]:
        for model in MODELS:
            for scen in ["own"] + (["DFT-anch"] if sid.endswith("Ta") else []):
                base = build(model, sid, E, nat, scen)
                if base is None:
                    continue
                ref = cls(base)
                out = {"system": sid, "model": model, "scenario": scen, "reference_class": ref}
                for knob, grid in grids.items():
                    flips = []
                    for v in grid:
                        c = cls(perturb(base, **{knob: float(v)}))
                        if c != ref:
                            dist = abs(np.log(v)) if knob == "w" else (abs(v - 1.0) if knob == "slope" else abs(v))
                            flips.append((dist, float(v), c))
                    b = min(flips) if flips else None
                    out[f"{knob}_flip"] = None if b is None else b[1]
                    out[f"{knob}_flip_to"] = None if b is None else b[2]
                rows.append(out)
                print(out, flush=True)
    pd.DataFrame(rows).to_csv(FR / "supplementary_data" / "x02_flip_thresholds.csv", index=False)


if __name__ == "__main__":
    main()
