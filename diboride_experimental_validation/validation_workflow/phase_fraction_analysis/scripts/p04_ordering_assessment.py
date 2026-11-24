"""Item 4: is an order-disorder transition plausible? Uses only existing configurations (m06_ordering_energies.csv).

Available configurations per model (x = 0.5 unless stated):
  metal:      192-atom quasi-random, alternating c-layers (= MP ordered ZrMB4, P6/mmm), AABB slab, BLADE 144-atom SQS
  vacancy:    k = 4 / 8 B vacancies spread, k = 4 clustered, k = 8 confined to one B plane (ZrB2, NbB2, TaB2)
  coupled:    one V_B in Zr|Zr, Zr|M, M|M layers of the AABB slab
NOT available: the DFT cluster-expansion ground states of TaB2-x (Ektarawong 2023, x = 0.167, 0.2, 0.25; structures
in that paper's SI, not accessed) and any ordered NbB2-x structure (Dahlqvist 2022, Nb16B28 etc.).

What can be stated without a configurational free-energy model:
  dE_ord   = E(ordered) - E(random), per metal atom (eV): the ordering energy of the one ordered candidate.
  T_cross  = |dE_ord| / (k_B ln 2): the temperature at which the fully ordered state and the ideal random solution
             have equal free energy (ordered entropy 0, random entropy k_B ln 2 per metal). NOT a transition
             temperature. Short-range order lowers the disordered free energy and partial disorder lowers the ordered
             one, so the true T_c can differ substantially; T_cross is only an energy-scale indicator.
Output: results/p04_ordering_assessment.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
MB = HERE.parent / "mlip_benchmark"
KB = 8.617333262e-5


def main() -> None:
    o = pd.read_csv(MB / "results" / "m06_ordering_energies.csv")
    rows = []
    met = o[o.type == "metal arrangement"].pivot_table(index=["model", "system"], columns="config", values="dH_mix_meV_atom")
    for (m, sid), r in met.iterrows():
        rand = r["random (192-atom quasi-random)"]
        for cfg, lab in (("metal-ordered alternating c-layers", "metal ordering (ZrMB4 c-layers)"), ("segregated AABB slab", "metal segregation (AABB)")):
            de_atom = (r[cfg] - rand) / 1000  # eV per atom (3 atoms per formula unit)
            de_metal = 3 * de_atom
            rows.append({"model": m, "system": sid, "kind": lab, "dE_vs_random_meV_atom": 1000 * de_atom, "dE_per_metal_eV": de_metal,
                         "favoured": bool(de_metal < -0.002 * 3),
                         "T_cross_K (energy-scale indicator, not Tc)": abs(de_metal) / (KB * np.log(2)) if de_metal < 0 else np.nan})
    vac = o[o.type == "vacancy arrangement"]
    for r in vac.itertuples():
        rows.append({"model": r.model, "system": r.system, "kind": f"vacancy arrangement {r.config}",
                     "dE_vs_random_meV_atom": np.nan, "dE_per_metal_eV": np.nan,
                     "interaction_eV_per_vacancy": r.interaction_eV_per_vacancy,
                     "favoured": bool(r.config in ("k4_cluster", "k8_plane") and r.interaction_eV_per_vacancy < -0.02)})
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results" / "p04_ordering_assessment.csv", index=False)
    pd.set_option("display.width", 220)
    print(df[df.kind.str.startswith("metal")].round(3).to_string(index=False))
    v = df[df.kind.str.startswith("vacancy")]
    print(v.pivot_table(index=["system", "kind"], columns="model", values="interaction_eV_per_vacancy").round(3).to_string())


if __name__ == "__main__":
    main()
