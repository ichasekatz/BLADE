"""Item 3: how much must individual energetic terms change to reverse each model's predicted phase assemblage?

Reference condition: x = 0.5, T = 2000 K, B/M = 1.95 (delta = 2.5 %), AlB2 vacancy fraction capped at y <= 0.125
(the DFT stability range of TaB2-x; no extrapolation). Assemblage classes: homogeneous / two AlB2 / metal-rich boride.
For each model ("own" energetics; ZrTa also "DFT-anchored" TaB2 V_B = -0.389 eV) the smallest change that flips the
class is found by scanning:
  d_m   shift of the Group-V B-vacancy energy (eV per vacancy)       -1.5 ... +1.5
  W     multiplier on the vacancy-vacancy interaction (both end members) 0.1 ... 10
  comp  shift of every metal-rich boride energy (eV per atom; ordered ZrMB4 / vacancy cells NOT shifted) -0.10 ... +0.10
Measured uncertainties (independent references, NOT calibrated confidence intervals):
  vacancy energy: model - DFT for TaB2 (Ektarawong 2023); model spread for NbB2 (no DFT value)
  W: model spread (no reference)
  boride formation energies: model - Materials Project (GGA) per atom, from m06_hull_consistency.csv
Outputs: results/p03_flip_thresholds.csv, results/p03_term_uncertainties.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
MB = HERE.parent / "mlip_benchmark"
sys.path.insert(0, str(MB / "scripts"))
sys.path.insert(0, str(HERE.parent / "scripts"))
sys.path.insert(0, str(HERE / "scripts"))
from m06_stage3_thermo import MODELS, build_system, energies  # noqa: E402
from p02_boron_loss_scan import DFT_TAB2_VB, Y_DFT, summarise  # noqa: E402
from s11_vacancy_sensitivity import perturbed  # noqa: E402

X0, T0, BM0 = 0.5, 2000.0, 1.95


def capped(s, base, comp_shift: float = 0.0):
    """Cap the AlB2 vacancy fraction and (optionally) shift ONLY the metal-rich borides (kind == reference with metal
    and boron), never the ordered ZrMB4 or vacancy-ordered AlB2 cells."""
    c = base.compounds.copy()
    if comp_shift:
        mrb = (c.kind == "reference") & (c.n_B > 0) & (c.n_Zr + c.n_M > 0)
        c.loc[mrb, "G"] += comp_shift * (c.loc[mrb, "n_Zr"] + c.loc[mrb, "n_M"] + c.loc[mrb, "n_B"])
    s.compounds = c
    s.y_grid = s.y_grid[s.y_grid <= Y_DFT + 1e-12]
    return s


def cls(s, comp) -> str:
    r = summarise(s.equilibrium(X0, BM0, T0), comp)
    if r.get("n_phases", 0) <= 1:
        return "homogeneous"
    return {"AlB2 (second solid solution)": "two AlB2"}.get(r["second_identity"], r["second_identity"])


def main() -> None:
    e, E = energies()
    natoms = dict(zip(e.structure_id, e.n_atoms))
    em = pd.read_csv(MB / "results" / "m02_end_member_benchmarks.csv").set_index(["model", "end_member"])
    rows = []
    for sid in ("ZrNb", "ZrTa"):
        for m in MODELS:
            base = build_system(m, sid, E, natoms)
            scen = {"own": 0.0}
            if sid == "ZrTa":
                scen["DFT-anchored"] = float(DFT_TAB2_VB - em.loc[(m, "TaB2"), "E_f_VB_Brich_eV"])
            for sname, d0 in scen.items():
                ref = cls(capped(perturbed(base, d_m=d0), base), base.compounds)
                out = {"system": sid, "model": m, "scenario": sname, "reference_class": ref}
                for knob, grid in (("d_m", np.round(np.arange(-1.5, 1.51, 0.05), 3)), ("w", np.geomspace(0.1, 10, 31)),
                                   ("comp", np.round(np.arange(-0.10, 0.101, 0.005), 3))):
                    flips = []
                    for v in grid:
                        kw = {"d_m": d0 + (v if knob == "d_m" else 0.0)}
                        if knob == "w":
                            kw["w"] = float(v)
                        c = cls(capped(perturbed(base, **kw), base, comp_shift=float(v) if knob == "comp" else 0.0), base.compounds)
                        if c != ref:
                            flips.append((abs(np.log(v)) if knob == "w" else abs(v), v, c))
                    best = min(flips) if flips else None
                    out[f"{knob}_flip_value"] = None if best is None else best[1]
                    out[f"{knob}_flip_to"] = None if best is None else best[2]
                rows.append(out)
                print(out, flush=True)
    pd.DataFrame(rows).to_csv(HERE / "results" / "p03_flip_thresholds.csv", index=False)

    # measured term uncertainties
    par = pd.read_csv(MB / "results" / "m06_model_parameters.csv")
    hull = pd.read_csv(MB / "results" / "m06_hull_consistency.csv")
    unc = []
    for sid, host in (("ZrNb", "NbB2"), ("ZrTa", "TaB2")):
        p = par[par.system == sid]
        v = p.dEv_MB2_eV
        unc.append({"system": sid, "term": f"V_B energy in {host} (eV)", "model_min": v.min(), "model_max": v.max(), "model_spread": v.max() - v.min(),
                    "independent_reference": "-0.389 (DFT, Ektarawong 2023)" if host == "TaB2" else "none (qualitative: ordered B vacancies stable, Dahlqvist 2022)",
                    "max_abs_dev_from_reference": float((v - DFT_TAB2_VB).abs().max()) if host == "TaB2" else np.nan})
        w = p.W_MB2_eV_fu
        unc.append({"system": sid, "term": f"vacancy-vacancy W in {host} (eV/fu)", "model_min": w.min(), "model_max": w.max(), "model_spread": w.max() - w.min(),
                    "independent_reference": "none", "max_abs_dev_from_reference": np.nan})
        s = p.slope_eV_per_nn
        unc.append({"system": sid, "term": "V_B environment slope (eV per Group-V neighbour)", "model_min": s.min(), "model_max": s.max(),
                    "model_spread": s.max() - s.min(), "independent_reference": "none (sign supported by Dahlqvist 2022 group trend)",
                    "max_abs_dev_from_reference": np.nan})
    b = hull[hull.phase.str.startswith("ref_mp") & hull.phase.str.contains("B")].copy()
    b = b[~b.phase.str.endswith("_B")]
    b["dev"] = b.dHf_eV_atom - b.MP_dHf_eV_atom
    for mm in ("Nb", "Ta"):
        g = b[b.system == f"{mm}-B"]
        unc.append({"system": f"{mm}-B borides", "term": "metal-rich boride formation energy, model - MP GGA (eV/atom)",
                    "model_min": g.dev.min(), "model_max": g.dev.max(), "model_spread": g.dev.max() - g.dev.min(),
                    "independent_reference": "Materials Project GGA (same structures)", "max_abs_dev_from_reference": float(g.dev.abs().max())})
    u = pd.DataFrame(unc)
    u.to_csv(HERE / "results" / "p03_term_uncertainties.csv", index=False)
    pd.set_option("display.width", 250)
    print(u.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
