"""Mixing enthalpy of the relaxed SQS structures behind the BLADE predictions (six pseudobinaries).

Reads the raw SQS energies (<system>/blade/*/energy + CONTCAR) — the inputs to the Redlich-Kister
fit — and reports ΔH_mix(x) = E(x) − [(1−x)·E(MIV B2) + x·E(MV B2)] in eV per MB2 formula unit,
plus the fitted RK value at x = 0.5. Positive values favour phase separation; negative values
favour mixing/ordering. Writes simulation_comparison/sqs_mixing_enthalpy.csv and
figures/fig6_sqs_mixing_enthalpy.(png|pdf).

Usage (from the BLADE repository root):
    .venv/bin/python diboride_experimental_validation/scripts/sqs_mixing_enthalpy.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
VAL, REPO = HERE.parent, HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))

from blade.analysis._stability_tables import _scan_sqs_dirs, build_phase_data, detect_metals  # noqa: E402
from compare_simulation import ORDER, RK_ORDER, STRUCT_DIR, SYSTEMS  # noqa: E402


def main() -> None:
    rows = []
    for sid in ORDER:
        tag, m4, m5 = SYSTEMS[sid]
        bdir = STRUCT_DIR / tag / "blade"
        metals = detect_metals(bdir, "B")
        sqs = _scan_sqs_dirs(bdir, metals, metals + ["B", "O"], [re.compile(rf"a_{m}=([0-9.]+)") for m in metals],
                             [f"y_{m}_metal_site" for m in metals], len(metals))
        e = {round(r[f"y_{m5}_metal_site"], 4): (r["energy_eV_per_formula"], r["phase_id"]) for r in sqs}
        e_iv, e_v = e[0.0][0], e[1.0][0]
        for x, (en, pid) in sorted(e.items()):
            rows.append({"system_id": sid, "x_groupV": x, "dH_mix_eV_per_fu": en - ((1 - x) * e_iv + x * e_v),
                         "kind": "SQS (relaxed, MLIP energy)", "structure_dir": pid})
        pdat = build_phase_data(bdir, metals, "B", 2.0, rk_order=RK_ORDER, y_step=0.01)
        xv = pdat["y_grid"][:, metals.index(m5)]
        i = int(np.argmin(np.abs(xv - 0.5)))
        rows.append({"system_id": sid, "x_groupV": 0.5, "dH_mix_eV_per_fu": pdat["phase_H0"][i] - 0.5 * (e_iv + e_v),
                     "kind": "RK fit", "structure_dir": ""})
    df = pd.DataFrame(rows)
    out = VAL / "simulation_comparison" / "sqs_mixing_enthalpy.csv"
    df.to_csv(out, index=False)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    for sid in ORDER:
        d = df[(df.system_id == sid) & (df.kind.str.startswith("SQS"))].sort_values("x_groupV")
        ax.plot(d.x_groupV, d.dH_mix_eV_per_fu * 1000, "o-", ms=3, label=f"{SYSTEMS[sid][1]}B$_2$–{SYSTEMS[sid][2]}B$_2$")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("x (Group-V MB$_2$ on metal sublattice)")
    ax.set_ylabel("ΔH$_{mix}$ (meV / MB$_2$)")
    ax.set_title("Mixing enthalpy of relaxed SQS (inputs to BLADE model)\n>0 favours separation; all six are <0 (mixing favoured)", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(VAL / "figures" / f"fig6_sqs_mixing_enthalpy.{ext}", dpi=200, bbox_inches="tight")
    print("wrote", out)
    literature_comparison(df)


EV_TO_KJ = 96.485  # kJ/mol per eV/f.u.


def element_energy(sym: str) -> float:
    """ORB energy per atom of the pure element, from any system's ORB/ folder (eV/atom)."""
    for d in sorted(STRUCT_DIR.glob(f"*/ORB/{sym}_mp-*_{sym}")):
        lines = (d / "CONTCAR").read_text().splitlines()
        n = sum(map(int, lines[6].split()))
        return float((d / "energy").read_text().split()[0]) / n
    raise FileNotFoundError(f"no ORB reference structure for {sym}")


def literature_comparison(df: pd.DataFrame) -> None:
    """ΔH_mix(0.5) vs DFT/calorimetry and end-member ΔfH vs calorimetry (fig7, fig8)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    th = pd.read_csv(VAL / "processed_data" / "thermochemical_data.csv")
    rows = []
    for sid in ORDER:
        ours = df[(df.system_id == sid) & (df.kind.str.startswith("SQS")) & (df.x_groupV == 0.5)].dH_mix_eV_per_fu.iloc[0]
        rows.append({"system_id": sid, "source": "this work: ORB-relaxed SQS", "dH_mix_0.5_eV_per_fu": ours, "kind": "MLIP"})
        t = th[(th.compound_or_system == sid) & th.quantity.str.startswith("DFT ΔHmix(0.5) = L0/4")]
        if len(t):
            rows.append({"system_id": sid, "source": "Tang 2025 DFT (LDA, SQS)", "dH_mix_0.5_eV_per_fu": t.value.iloc[0] / EV_TO_KJ, "kind": "DFT"})
        v = th[(th.compound_or_system == sid) & (th.source_id == "S16")]
        if len(v):
            rows.append({"system_id": sid, "source": "Vorotilo 2022 DFT", "dH_mix_0.5_eV_per_fu": v.value.iloc[0] * 3, "kind": "DFT"})
        c = th[(th.compound_or_system == sid) & (th.record_id == "T08")]
        if len(c):
            rows.append({"system_id": sid, "source": "Ness calorimetry Ω1/4 (per mol MB2 assumed)", "dH_mix_0.5_eV_per_fu": c.value.iloc[0] / EV_TO_KJ, "kind": "experiment"})
    lit = pd.DataFrame(rows)
    lit.to_csv(VAL / "simulation_comparison" / "mixing_enthalpy_vs_literature.csv", index=False)

    # End-member formation enthalpies: ORB vs oxide melt solution calorimetry.
    form = []
    eb = element_energy("B")
    for sid in ORDER:
        tag, m4, m5 = SYSTEMS[sid]
        bdir = STRUCT_DIR / tag / "blade"
        metals = detect_metals(bdir, "B")
        sqs = _scan_sqs_dirs(bdir, metals, metals + ["B", "O"], [re.compile(rf"a_{m}=([0-9.]+)") for m in metals],
                             [f"y_{m}_metal_site" for m in metals], len(metals))
        for m in (m4, m5):
            if any(f["compound"] == f"{m}B2" for f in form):
                continue
            e = next(r["energy_eV_per_formula"] for r in sqs if abs(r[f"y_{m}_metal_site"] - 1) < 1e-9)
            dfh = (e - element_energy(m) - 2 * eb) * EV_TO_KJ
            exp = th[(th.compound_or_system == f"{m}B2") & th.quantity.str.startswith("standard enthalpy of formation")]
            form.append({"compound": f"{m}B2", "dfH_ORB_kJ_per_mol": dfh,
                         "dfH_OMSC_kJ_per_mol": exp.value.iloc[0] if len(exp) else np.nan,
                         "OMSC_uncertainty_2SE": exp.uncertainty.iloc[0] if len(exp) else np.nan})
    fdf = pd.DataFrame(form)
    fdf["difference_ORB_minus_OMSC"] = fdf.dfH_ORB_kJ_per_mol - fdf.dfH_OMSC_kJ_per_mol
    fdf.to_csv(VAL / "simulation_comparison" / "formation_enthalpy_vs_calorimetry.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.5, 4))
    srcs = list(dict.fromkeys(lit.source))
    w = 0.8 / len(srcs)
    for k, src in enumerate(srcs):
        d = lit[lit.source == src].set_index("system_id").reindex(ORDER)
        ax.bar(np.arange(len(ORDER)) + (k - (len(srcs) - 1) / 2) * w, d["dH_mix_0.5_eV_per_fu"] * 1000, w, label=src)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(ORDER)), [f"{SYSTEMS[s][1]}B$_2$–{SYSTEMS[s][2]}B$_2$" for s in ORDER], fontsize=8)
    ax.set_ylabel("ΔH$_{mix}$(x=0.5) (meV / MB$_2$)")
    ax.set_title("Mixing enthalpy at x = 0.5: calculations vs calorimetry (only ZrB$_2$–TaB$_2$ measured)", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(VAL / "figures" / f"fig7_mixing_enthalpy_vs_literature.{ext}", dpi=200, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.errorbar(fdf.dfH_OMSC_kJ_per_mol, fdf.dfH_ORB_kJ_per_mol, xerr=fdf.OMSC_uncertainty_2SE, fmt="o")
    for _, r in fdf.iterrows():
        ax.annotate(r.compound, (r.dfH_OMSC_kJ_per_mol, r.dfH_ORB_kJ_per_mol), xytext=(4, 4), textcoords="offset points", fontsize=8)
    lim = [min(fdf.dfH_OMSC_kJ_per_mol.min(), fdf.dfH_ORB_kJ_per_mol.min()) - 20, max(fdf.dfH_OMSC_kJ_per_mol.max(), fdf.dfH_ORB_kJ_per_mol.max()) + 20]
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.set_xlabel("ΔfH, calorimetry (Ness 2025, kJ/mol)")
    ax.set_ylabel("ΔfH, ORB (this work, kJ/mol)")
    ax.set_title("End-member formation enthalpies", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(VAL / "figures" / f"fig8_formation_enthalpy_vs_calorimetry.{ext}", dpi=200, bbox_inches="tight")
    print(lit.to_string(index=False))
    print(fdf.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
