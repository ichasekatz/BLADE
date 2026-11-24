"""Stage 1 projection (approximate; NOT the Stage 3 refit): what do the other MLIPs' vacancy energies imply?

For each model, the B-vacancy offsets relative to ORB are inserted into the ORB-built framework (s11 `perturbed`):
  d_zr  = E_v,model(ZrB2) - E_v,ORB(ZrB2)
  d_m   = E_v,model(MB2)  - E_v,ORB(MB2)
  slope = (dE_v/dn at x = 0.5)_model / (dE_v/dn)_ORB
Everything else (stoichiometric mixing, vacancy-vacancy W, competitor energies) stays ORB. W and the competitors
are Stage 2 quantities. Outcome classes as in s11. Output: mlip_benchmark/results/m03_projection.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

MB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MB.parent / "scripts"))
import s08_framework_equilibria as s08  # noqa: E402
from s11_vacancy_sensitivity import classify, perturbed  # noqa: E402


def main() -> None:
    em = pd.read_csv(MB / "results" / "m02_end_member_benchmarks.csv").set_index(["model", "end_member"])
    tr = pd.read_csv(MB / "results" / "m02_environment_trends.csv").set_index(["model", "system"])
    rows = []
    for sid, host in (("ZrNb", "NbB2"), ("ZrTa", "TaB2")):
        for pset in ("ORB", "DFT"):
            base = s08.build(sid, pset)
            for m in ("orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"):
                d_zr = em.loc[(m, "ZrB2"), "E_f_VB_Brich_eV"] - em.loc[("orb", "ZrB2"), "E_f_VB_Brich_eV"]
                d_m = em.loc[(m, host), "E_f_VB_Brich_eV"] - em.loc[("orb", host), "E_f_VB_Brich_eV"]
                slope = tr.loc[(m, sid), "slope_eV_per_solute_nn"] / tr.loc[("orb", sid), "slope_eV_per_solute_nn"]
                sysm = perturbed(base, d_m=d_m, d_zr=d_zr, slope=slope)
                for t in (1500.0, 2000.0, 2500.0):
                    for bm in (2.0, 1.98, 1.95, 1.90, 1.80):
                        ph = sysm.equilibrium(0.5, bm, t)
                        tag, dx, ym = classify(ph)
                        other = ph[(ph.phase != "AlB2") & ph.x_M.notna()]
                        rows.append({"system": sid, "mixing_basis": "ORB" if pset == "ORB" else "DFT (Tang)", "model_vacancy_energies": m,
                                     "d_zr_eV": d_zr, "d_m_eV": d_m, "slope_ratio": slope, "T_K": t, "B_per_metal": bm, "class": tag,
                                     "AlB2_dx": dx, "AlB2_x": ";".join(f"{v:.2f}" for v in ph[ph.phase == "AlB2"].x_M),
                                     "other_phases": ";".join(other.phase),
                                     "metal_fraction_other": float(other.metal_fraction.sum()) if len(other) else 0.0})
    df = pd.DataFrame(rows)
    df.to_csv(MB / "results" / "m03_projection.csv", index=False)
    pd.set_option("display.width", 220)
    piv = df[df.T_K == 2000.0].pivot_table(index=["system", "mixing_basis", "model_vacancy_energies"], columns="B_per_metal",
                                          values="class", aggfunc="first")
    print(piv.to_string())


if __name__ == "__main__":
    main()
