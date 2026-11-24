"""Step 5: random vs layered / clustered / ordered arrangements at identical composition (H2).

dH_mix(config) = [E(config) - (1-x) E(ZrB2) - x E(MB2)] / n_atoms, with end members relaxed in the SAME
supercell and with the same settings (C family: 4x4x4, 192 atoms; L family: 2x2x8, 96 atoms). Also listed:
BLADE's 144-atom SQS (4x4x3) and the MP ordered ZrMB4 cell (6 atoms), each with its own end-member reference.
SQS sampling noise sigma_SQS = |dH(qsqs1) - dH(qsqs2)| / sqrt(2), pooled per system.

Outputs: results/s05_configuration_energies.csv, figures/fig_s05_configurations.(png|pdf)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

WF = Path(__file__).resolve().parents[1]
RES = WF / "results"
SOL = {"ZrNb": "Nb", "ZrTa": "Ta"}


def main() -> None:
    e = pd.read_csv(RES / "s03_orb_energies.csv").set_index("structure_id")
    man = pd.read_csv(WF / "inputs" / "structure_manifest.csv").set_index("structure_id")
    audit = pd.read_csv(RES / "s00_existing_structure_audit.csv")
    rows = []
    for sid, m in SOL.items():
        for fam, end0, end1 in (("C", f"C_{sid}_x0.00_end", f"C_{sid}_x1.00_end"), ("L", f"L_{sid}_end0", f"L_{sid}_end1")):
            if end0 not in e.index or end1 not in e.index:
                continue
            e0, e1 = e.loc[end0, "energy_eV_atom"], e.loc[end1, "energy_eV_atom"]
            sub = man[(man.family == fam) & (man.system == sid) & (man.x_solute > 0) & (man.x_solute < 1)]
            for s, r in sub.iterrows():
                if s not in e.index:
                    continue
                dh = e.loc[s, "energy_eV_atom"] - (1 - r.x_solute) * e0 - r.x_solute * e1
                rows.append({"system": sid, "x_solute": r.x_solute, "cell": r.cell, "structure_id": s, "config": r.config,
                             "dH_mix_meV_atom": 1000 * dh, "converged": bool(e.loc[s, "converged"])})
        # BLADE 144-atom SQS (reused) and MP ordered ZrMB4
        h = audit[(audit.system == sid) & (audit.origin == "SQS (BLADE/ATAT)") & (audit.n_atoms == 144)]
        h = h[h.path.str.contains(f"system_structures/{'NbZrB' if sid == 'ZrNb' else 'TaZrB'}/blade")]
        e0b = h[h.x_solute == 0].energy_eV_atom.iloc[0]
        e1b = h[h.x_solute == 1].energy_eV_atom.iloc[0]
        for r in h[(h.x_solute > 0) & (h.x_solute < 1)].itertuples():
            rows.append({"system": sid, "x_solute": r.x_solute, "cell": "BLADE SQS 4x4x3 (144)", "structure_id": r.path.split("/")[-1],
                         "config": "SQS (BLADE, reused)", "dH_mix_meV_atom": 1000 * (r.energy_eV_atom - (1 - r.x_solute) * e0b - r.x_solute * e1b),
                         "converged": True})
        mp = {"ZrNb": "comp_mp-1215211_ZrNbB4", "ZrTa": "comp_mp-1215209_ZrTaB4"}[sid]
        if mp in e.index:
            rows.append({"system": sid, "x_solute": 0.5, "cell": "MP ordered (6)", "structure_id": mp, "config": "MP ordered ZrMB4 (c-layer period 1)",
                         "dH_mix_meV_atom": 1000 * (e.loc[mp, "energy_eV_atom"] - 0.5 * e0b - 0.5 * e1b), "converged": bool(e.loc[mp, "converged"])})
    df = pd.DataFrame(rows)
    if df.empty:
        print("no configuration results yet")
        return
    # SQS noise and comparison with the same-cell quasi-random mean
    out = []
    for (sid, x), g in df.groupby(["system", "x_solute"]):
        q = g[g.config.str.startswith("quasi-random")].dH_mix_meV_atom
        ref = q.mean() if len(q) else g[g.config.str.startswith("SQS")].dH_mix_meV_atom.mean()
        sig = abs(q.iloc[0] - q.iloc[1]) / np.sqrt(2) if len(q) == 2 else np.nan
        g = g.assign(ref_random_meV_atom=ref, sigma_sqs_meV_atom=sig, minus_random_meV_atom=g.dH_mix_meV_atom - ref)
        out.append(g)
    df = pd.concat(out).sort_values(["system", "x_solute", "dH_mix_meV_atom"])
    df.to_csv(RES / "s05_configuration_energies.csv", index=False)
    pd.set_option("display.width", 220)
    print(df.round(2).to_string(index=False))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=False)
    marks = {"quasi": ("o", "C0"), "SQS": ("s", "k"), "c-layered period 1": ("D", "C3"), "MP ordered": ("*", "C3"),
             "c-layered AABB": ("v", "C2"), "a-rows alternating": ("<", "C1"), "a-rows 2+2": (">", "C1"),
             "single c-layer": ("D", "C3"), "single a-row": ("<", "C1"), "2x2 column": ("P", "C4"), "c-superlattice": ("^", "C5")}
    for ax, sid in zip(axes, SOL):
        g = df[df.system == sid]
        for r in g.itertuples():
            mk, col = next((v for k, v in marks.items() if r.config.startswith(k)), ("x", "grey"))
            ax.plot(r.x_solute + (0.01 if r.cell.startswith("2x2x8") else 0), r.dH_mix_meV_atom, mk, color=col, ms=7, mfc="none" if "2x2x8" in r.cell else col)
        for k, (mk, col) in marks.items():
            ax.plot([], [], mk, color=col, label=k)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"{sid[:2]}B$_2$–{sid[2:]}B$_2$ (ORB, same-cell end-member reference)", fontsize=9)
        ax.set_xlabel("x (solute)")
        ax.set_ylabel("ΔH$_{mix}$ (meV/atom)")
    axes[1].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    (WF / "figures").mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(WF / "figures" / f"fig_s05_configurations.{ext}", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
