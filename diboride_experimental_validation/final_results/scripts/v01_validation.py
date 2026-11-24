"""Validation checks for the final results (no new relaxations).

1. Phase-equilibrium mass balance: for every system/model/scenario at x = 0.25/0.5/0.75, delta = 0-20 %, T = 2000 K,
   rebuild the overall Zr(host)/M/B content from the reported phases and compare with the imposed composition.
2. Grid convergence: repeat a subset with a 2x finer x grid and 2x denser vacancy grid; compare the second-phase metal
   fraction and the assemblage class.
3. Relaxation convergence: fraction of MLIP relaxations that reached fmax <= 0.001 eV/A, from CALCULATION_MANIFEST.csv.
Output: supplementary_data/v01_mass_balance.csv, v01_grid_convergence.csv, v01_relaxation_convergence.csv, and a
printed summary.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
SD = FR / "supplementary_data"
sys.path.insert(0, str(FR / "scripts"))
from x01_six_system_thermo import MODELS, Y_CAP, build, energies, solve, summarise  # noqa: E402

SYSTEMS = ["TiNb", "TiTa", "ZrNb", "ZrTa", "HfNb", "HfTa"]


def composition(ph: pd.DataFrame, comp: pd.DataFrame) -> tuple[float, float, float]:
    """Return (host metal, Group-V metal, B) per formula-unit basis of the LP (sum of metal = 1)."""
    cmap = comp.set_index("name")
    h = m = b = 0.0
    for r in ph.itertuples():
        if r.phase == "AlB2":
            h += r.amount_fu * (1 - r.x_M)
            m += r.amount_fu * r.x_M
            b += r.amount_fu * 2 * (1 - r.y_Va)
        else:
            c = cmap.loc[r.phase]
            h += r.amount_fu * c.n_Zr
            m += r.amount_fu * c.n_M
            b += r.amount_fu * c.n_B
    return h, m, b


def main() -> None:
    E, nat = energies()
    mb, gc = [], []
    for sid in SYSTEMS:
        for model in MODELS:
            for scen in ["own"] + (["DFT-anch"] if sid.endswith("Ta") else []):
                s = build(model, sid, E, nat, scen)
                if s is None:
                    continue
                for x in (0.25, 0.5, 0.75):
                    for d in (0.0, 0.01, 0.025, 0.05, 0.10, 0.20):
                        ph = solve(s, x, 2 * (1 - d), 2000.0)
                        if "phase" not in ph.columns:
                            mb.append({"system": sid, "model": model, "scenario": scen, "x": x, "delta": d, "status": "LP failed"})
                            continue
                        h, m, b = composition(ph, s.compounds)
                        mb.append({"system": sid, "model": model, "scenario": scen, "x": x, "delta": d, "status": "ok",
                                   "metal_sum_err": abs(h + m - 1), "x_err": abs(m / (h + m) - x), "BM_err": abs(b / (h + m) - 2 * (1 - d))})
                # grid convergence at x = 0.5, 2000 K
                fine = type(s)(metal=s.metal, alb2=s.alb2, compounds=s.compounds, x_grid=np.linspace(0, 1, 201),
                               y_grid=np.r_[0.0, np.geomspace(1e-5, Y_CAP, 90)])
                for d in (0.01, 0.025, 0.05, 0.10):
                    a = summarise(solve(s, 0.5, 2 * (1 - d), 2000.0), s.compounds)
                    f = summarise(solve(fine, 0.5, 2 * (1 - d), 2000.0), fine.compounds)
                    gc.append({"system": sid, "model": model, "scenario": scen, "delta": d, "class_default": a["class"], "class_fine": f["class"],
                               "frac_default": a["second_metal_frac"], "frac_fine": f["second_metal_frac"],
                               "abs_frac_diff": abs(a["second_metal_frac"] - f["second_metal_frac"])})
            print(sid, model, "checked", flush=True)
    mb, gc = pd.DataFrame(mb), pd.DataFrame(gc)
    mb.to_csv(SD / "v01_mass_balance.csv", index=False)
    gc.to_csv(SD / "v01_grid_convergence.csv", index=False)
    ok = mb[mb.status == "ok"]
    print(f"mass balance: {len(ok)} states, LP failures {int((mb.status != 'ok').sum())}, max |dx| {ok.x_err.max():.2e}, "
          f"max |dB/M| {ok.BM_err.max():.2e}, max |metal sum - 1| {ok.metal_sum_err.max():.2e}")
    print(f"grid convergence: {len(gc)} states, class changed in {int((gc.class_default != gc.class_fine).sum())}, "
          f"median |d frac| {gc.abs_frac_diff.median():.3f}, max {gc.abs_frac_diff.max():.3f}")
    man = FR / "CALCULATION_MANIFEST.csv"
    if man.exists():
        c = pd.read_csv(man)
        r = c[c.fmax_final_eV_A.notna()].groupby("model").agg(n=("structure_id", "size"), converged=("converged", "mean"),
                                                                max_fmax=("fmax_final_eV_A", "max"), max_steps=("n_steps", "max"))
        r.to_csv(SD / "v01_relaxation_convergence.csv")
        print(r.to_string())


if __name__ == "__main__":
    main()
