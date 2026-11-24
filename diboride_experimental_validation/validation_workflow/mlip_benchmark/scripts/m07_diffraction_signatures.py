"""Which measurements distinguish the candidate decomposition products? (item 9)

From one model's relaxed structures (default MatterSim; ORB also written), compute Cu Ka1 powder patterns
(pymatgen XRDCalculator) of:
  A  homogeneous (Zr0.5M0.5)B2                       relaxed x = 0.5 SQS host
  B  two AlB2 solid solutions (ORB-type prediction)  relaxed x = 0.25 and x = 0.75 hosts as stand-ins for the Zr-rich
                                                     and Group-V-rich phases, plus the Group-V-rich B-deficient MB2-d
                                                     (k = 8 spread cell)
  C  AlB2 + metal-rich boride (four-model prediction) x = 0.25 host + M3B4, M5B6, MB (relaxed reference cells)
Reports the strongest reflections of each phase, and the reflections of the non-AlB2 borides lying more than
0.3 deg from every AlB2 reflection with >= 5 % relative intensity: the diffraction fingerprint of scenario C.
Also reports a and c of each relaxed phase (vacancy-induced change of c), for comparison with Otani/Eckert.

Outputs: mlip_benchmark/results/m07_xrd_peaks_<model>.csv, m07_distinguishing_reflections_<model>.csv,
         m07_lattice_<model>.csv, figures/fig_m07_xrd_<model>.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

MB = Path(__file__).resolve().parents[1]
OUT = MB / "results"


def structure(model: str, sid: str):
    from pymatgen.core import Structure

    if model == "orb":
        man = pd.concat([pd.read_csv(MB / "pilot_manifest.csv"), pd.read_csv(MB / "stage2_manifest.csv")]).set_index("structure_id")
        wf = MB.parent
        cand = [wf / "results" / "orb" / sid.replace("ref_", "comp_") / "CONTCAR", wf / "results" / "orb" / sid / "CONTCAR",
                Path(man.loc[sid, "start_geometry"])]
        p = next(c for c in cand if c.exists())
    else:
        p = MB / "results" / model / sid / "CONTCAR"
    return Structure.from_file(p)


def main(model: str = "mattersim") -> None:
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    calc = XRDCalculator(wavelength="CuKa1")
    rows, lat = [], []
    for sid_, m in (("ZrNb", "Nb"), ("ZrTa", "Ta")):
        phases = {
            "A: (Zr0.5M0.5)B2 homogeneous": f"host_{sid_}_x0.50",
            "B/C: Zr-rich (Zr0.75M0.25)B2": f"host_{sid_}_x0.25",
            "B: M-rich (Zr0.25M0.75)B2": f"host_{sid_}_x0.75",
            f"B: {m}B2-d (B/M=1.83)": f"Y_{m}B2_k8_spread",
            f"C: {m}3B4": {"Nb": "ref_mp-10255_Nb3B4", "Ta": "ref_mp-10142_Ta3B4"}[m],
            f"C: {m}5B6": {"Nb": "ref_mp-1102394_Nb5B6", "Ta": "ref_mp-28629_Ta5B6"}[m],
            f"C: {m}B": {"Nb": "ref_mp-2580_NbB", "Ta": "ref_mp-1097_TaB"}[m],
        }
        if m == "Nb":
            phases["C: Nb2B3"] = "ref_mp-1080021_Nb2B3"
        pats = {}
        for label, sid in phases.items():
            try:
                st = structure(model, sid)
            except (FileNotFoundError, StopIteration, KeyError):
                continue
            if "B2" in label and "B4" not in label:
                a = float(np.mean(st.lattice.abc[:2])) / 4
                c = st.lattice.abc[2] / 3
                lat.append({"system": sid_, "phase": label, "a_A": a, "c_A": c})
            pat = calc.get_pattern(st, two_theta_range=(20, 80))
            pats[label] = pat
            order = np.argsort(pat.y)[::-1][:8]
            for i in order:
                rows.append({"system": sid_, "phase": label, "two_theta": round(float(pat.x[i]), 3), "rel_intensity": round(float(pat.y[i]), 1),
                             "hkl": ",".join("".join(str(v) for v in h["hkl"]) for h in pat.hkls[i])})
        alb2_peaks = np.concatenate([pats[k].x[pats[k].y >= 5] for k in pats if "B2" in k and "B4" not in k])
        dist = []
        for k, pat in pats.items():
            if not k.startswith("C: ") or "B2" in k:
                continue
            for x, y, h in zip(pat.x, pat.y, pat.hkls):
                if y >= 10 and np.min(np.abs(alb2_peaks - x)) > 0.3:
                    dist.append({"system": sid_, "phase": k, "two_theta": round(float(x), 3), "rel_intensity": round(float(y), 1),
                                 "hkl": ",".join("".join(str(v) for v in hh["hkl"]) for hh in h),
                                 "min_separation_from_AlB2_deg": round(float(np.min(np.abs(alb2_peaks - x))), 2)})
        pd.DataFrame(dist).to_csv(OUT / f"m07_distinguishing_reflections_{model}_{sid_}.csv", index=False)
        plot(pats, sid_, model)
    pd.DataFrame(rows).to_csv(OUT / f"m07_xrd_peaks_{model}.csv", index=False)
    pd.DataFrame(lat).to_csv(OUT / f"m07_lattice_{model}.csv", index=False)
    print(pd.DataFrame(lat).round(4).to_string(index=False))
    for sid_ in ("ZrNb", "ZrTa"):
        d = pd.read_csv(OUT / f"m07_distinguishing_reflections_{model}_{sid_}.csv")
        print(sid_, d.sort_values("rel_intensity", ascending=False).head(8).to_string(index=False))


def plot(pats: dict, sid_: str, model: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (k, p) in enumerate(pats.items()):
        ax.vlines(p.x, i * 110, i * 110 + p.y, lw=1.2, color=f"C{i}")
        ax.text(80.5, i * 110 + 20, k, fontsize=7, va="bottom")
    ax.set_xlim(20, 80)
    ax.set_yticks([])
    ax.set_xlabel("2θ (Cu Kα1, deg)")
    ax.set_title(f"{sid_}: simulated powder patterns of candidate products ({model}-relaxed structures)", fontsize=9)
    fig.tight_layout()
    (MB / "figures").mkdir(exist_ok=True)
    fig.savefig(MB / "figures" / f"fig_m07_xrd_{model}_{sid_}.png", dpi=180, bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "mattersim")
