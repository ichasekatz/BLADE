"""Stage 2 analysis + Stage 3 thermodynamics with each model's internally consistent energetics.

For every model in {orb, mace-mpa, mace-omat, mattersim, sevennet} and each system (ZrNb, ZrTa):
  * energies   Stage 1 (pilot) + Stage 2 results of that model only (ORB: the stored ORB screen)
  * mixing     RK (L0, L1) fit to the model's own x = 0.25/0.5/0.75 SQS hosts vs its own end members
  * vacancies  dE_v (ZrB2, MB2), environment slope at x = 0.5, vacancy-vacancy W from k = 2/4/8 spread cells
  * phases     AlB2 CEF (Zr,M)(B,Va)2 with ideal entropy on both sublattices + W y^2 term (thermo_framework);
               line compounds: alpha-B, Zr, bcc Nb, bcc Ta, Nb2B3, Nb3B4, Nb5B6, NbB, Ta3B4, Ta5B6, TaB (model's own),
               ordered ZrMB4 (= alternating c-layer 192-atom cell), and the vacancy-ORDERED cells
               (k = 4 clustered, k = 8 single-plane) as candidate ordered B-deficient AlB2 phases.
Checks written: per-model M-B hull vs MP, 0 K B-rich stability of B-deficient cells, metal / vacancy / coupled
ordering energies. Equilibria: x = 0.05..0.95, T = 1000..3000 K, B/M = 2.00..1.80, closed system (mass balance).
Classes: HOMOGENEOUS | AlB2_GAP | AlB2+BORIDE | ORDERED_VAC | ORDERED_ZrMB4 (combinations joined by '+').

Outputs: mlip_benchmark/results/m06_*.csv, figures/fig_m06_*.png
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MB = Path(__file__).resolve().parents[1]
WF = MB.parent
sys.path.insert(0, str(WF / "scripts"))
from thermo_framework import AlB2Params, System  # noqa: E402

OUT = MB / "results"
FIG = MB / "figures"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
SYS = {"ZrNb": "Nb", "ZrTa": "Ta"}
EM_V1 = {"ZrB2": "V_ZrB2_x0.00_VB_n0a", "NbB2": "V_ZrNb_x1.00_VB_n6a", "TaB2": "V_ZrTa_x1.00_VB_n6a"}
REFS = {"B": "ref_mp-160_B", "Zr": "ref_mp-131_Zr", "Nb": "ref_mp-75_Nb", "Ta": "ref_mp-50_Ta"}
BORIDES = {"Nb": ["ref_mp-1080021_Nb2B3", "ref_mp-10255_Nb3B4", "ref_mp-1102394_Nb5B6", "ref_mp-2580_NbB"],
           "Ta": ["ref_mp-10142_Ta3B4", "ref_mp-28629_Ta5B6", "ref_mp-1097_TaB"]}
TANG = {"ZrNb": [-2.515, 2.418], "ZrTa": [-9.306, 4.016]}  # kJ/mol MB2 (DFT mixing basis for sensitivity)


def formula_counts(sid: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", sid.split("_")[-1])}


def energies() -> tuple[pd.DataFrame, dict]:
    man = pd.concat([pd.read_csv(MB / "pilot_manifest.csv"), pd.read_csv(MB / "stage2_manifest.csv")], ignore_index=True)
    rows = []
    for r in man.itertuples():
        rows.append({"model": "orb", "structure_id": r.structure_id, "energy_eV": r.orb_energy_eV, "n_atoms": r.n_atoms, "converged": True})
        for m in MODELS[1:]:
            p = MB / "results" / m / r.structure_id / "result.json"
            if p.exists():
                j = json.loads(p.read_text())
                rows.append({"model": m, "structure_id": r.structure_id, "energy_eV": j["energy_eV"], "n_atoms": j["n_atoms"],
                             "converged": j["converged"]})
    e = pd.DataFrame(rows)
    return e, {(r.model, r.structure_id): r.energy_eV for r in e.itertuples()}


def model_params(m: str, sid: str, E: dict, natoms: dict, mixing: str = "own") -> dict | None:
    mm = SYS[sid]
    need = ["host_ZrB2", f"host_{mm}B2", f"host_{sid}_x0.25", f"host_{sid}_x0.50", f"host_{sid}_x0.75", REFS["B"]]
    if any((m, k) not in E for k in need):
        return None
    e_b = E[(m, REFS["B"])] / natoms[REFS["B"]]
    efu = {0.0: E[(m, "host_ZrB2")] / 48, 1.0: E[(m, f"host_{mm}B2")] / 48}
    for x in (0.25, 0.5, 0.75):
        efu[x] = E[(m, f"host_{sid}_x{x:.2f}")] / 48
    xs = np.array([0.25, 0.5, 0.75])
    mix = np.array([efu[x] - (1 - x) * efu[0.0] - x * efu[1.0] for x in xs])
    a = np.vstack([xs * (1 - xs) * (2 * xs - 1) ** k for k in range(2)]).T
    rk = list(np.linalg.lstsq(a, mix, rcond=None)[0]) if mixing == "own" else [v / 96.485 for v in TANG[sid]]
    ev = {em: E[(m, EM_V1[em])] - E[(m, f"host_{em}")] + e_b for em in ("ZrB2", f"{mm}B2") if (m, EM_V1[em]) in E}
    env = [(n, E[(m, k)] - E[(m, f"host_{sid}_x0.50")] + e_b) for n in range(7) for k in (f"V_{sid}_x0.50_VB_n{n}a",) if (m, k) in E]
    slope_fit = np.polyfit([p[0] for p in env], [p[1] for p in env], 1)
    w = {}
    for em in ("ZrB2", f"{mm}B2"):
        e1 = E[(m, EM_V1[em])] - E[(m, f"host_{em}")]
        pts = []
        for k in (2, 4, 8):
            key = (m, f"Y_{em}_k{k}_spread")
            if key in E:
                exc = E[key] - E[(m, f"host_{em}")] - k * e1
                pts.append((k / 96, exc / 48))
        w[em] = float(sum(y2 * ex for y2, ex in ((p[0] ** 2, p[1]) for p in pts)) / sum(p[0] ** 4 for p in pts)) if pts else 0.0
    return {"e_b": e_b, "efu": efu, "rk": rk, "dEv_zr": ev["ZrB2"], "dEv_m": ev[f"{mm}B2"], "slope_fit": slope_fit,
            "w_zr": w["ZrB2"], "w_m": w[f"{mm}B2"], "mix_meV_atom": dict(zip(xs, 1000 * mix / 3))}


def build_system(m: str, sid: str, E: dict, natoms: dict, mixing: str = "own", include_ordered: bool = True) -> System | None:
    p = model_params(m, sid, E, natoms, mixing)
    if p is None:
        return None
    mm = SYS[sid]
    f = p["slope_fit"]
    alb2 = AlB2Params(g_zrb=p["efu"][0.0], g_mb=p["efu"][1.0], dEv_zr=p["dEv_zr"], dEv_m=p["dEv_m"], e_b=p["e_b"], rk=p["rk"],
                      ev_env=lambda n, f=f: float(np.polyval(f, n)), w_vv_zr=p["w_zr"], w_vv_m=p["w_m"])
    rows = []
    for ref in [REFS["B"], REFS["Zr"], REFS[mm]] + BORIDES[mm]:
        if (m, ref) not in E:
            continue
        c = formula_counts(ref)
        nfu = natoms[ref] / sum(c.values())
        rows.append({"name": ref.split("_", 1)[1], "n_Zr": c.get("Zr", 0), "n_M": c.get(mm, 0), "n_B": c.get("B", 0), "G": E[(m, ref)] / nfu,
                     "kind": "reference"})
    if include_ordered:
        key = (m, f"C_{sid}_x0.50_clay1")
        if key in E:  # 192 atoms = 32 ZrMB4 formula units
            rows.append({"name": f"ordered Zr{mm}B4 (c-layer)", "n_Zr": 1, "n_M": 1, "n_B": 4, "G": E[key] / 32, "kind": "ordered_metal"})
        for em, nzr, nm in (("ZrB2", 48, 0), (f"{mm}B2", 0, 48)):
            for tag, k in (("k4_cluster", 4), ("k8_plane", 8)):
                key = (m, f"Y_{em}_{tag}")
                if key in E:
                    rows.append({"name": f"{em}-{k}V ordered ({tag})", "n_Zr": nzr, "n_M": nm, "n_B": 96 - k, "G": E[key], "kind": "ordered_vacancy"})
    return System(metal=mm, alb2=alb2, compounds=pd.DataFrame(rows), x_grid=np.linspace(0, 1, 101),
                  y_grid=np.r_[0.0, np.geomspace(1e-5, 0.35, 60)])


def classify(ph: pd.DataFrame, comp: pd.DataFrame) -> str:
    kind = dict(zip(comp.name, comp.kind))
    al = ph[ph.phase == "AlB2"]
    tags = []
    if len(al) > 1 and al.x_M.max() - al.x_M.min() > 0.05:
        tags.append("AlB2_GAP")
    for p in ph.phase:
        k = kind.get(p)
        if k == "ordered_vacancy":
            tags.append("ORDERED_VAC")
        elif k == "ordered_metal":
            tags.append("ORDERED_ZrMB4")
        elif k == "reference" and not p.endswith(("_B", "_Zr", "_Nb", "_Ta")):
            tags.append("AlB2+BORIDE")
    return "+".join(sorted(set(tags))) if tags else ("HOMOGENEOUS" if len(al) else "NO_AlB2")


def checks(e: pd.DataFrame, E: dict, natoms: dict) -> None:
    """Per-model hull consistency, 0 K B-rich stability of B-deficient cells, ordering energies."""
    mp = pd.read_csv(WF / "inputs" / "competitors" / "manifest.csv").set_index("material_id")
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.core import Composition

    hull_rows, vac0_rows, ord_rows, par_rows = [], [], [], []
    for m in MODELS:
        for mm in ("Nb", "Ta"):
            ents = []
            for ref in [REFS["B"], REFS[mm]] + BORIDES[mm]:
                if (m, ref) in E:
                    c = formula_counts(ref)
                    nfu = natoms[ref] / sum(c.values())
                    ents.append(PDEntry(Composition({k: v * nfu for k, v in c.items()}), E[(m, ref)], name=ref))
            if (m, f"host_{mm}B2") in E:
                ents.append(PDEntry(Composition({mm: 48, "B": 96}), E[(m, f"host_{mm}B2")], name=f"host_{mm}B2"))
            extra = []
            for tag, k in (("k2_spread", 2), ("k4_spread", 4), ("k8_spread", 8), ("k4_cluster", 4), ("k8_plane", 8)):
                key = (m, f"Y_{mm}B2_{tag}")
                if key in E:
                    extra.append(PDEntry(Composition({mm: 48, "B": 96 - k}), E[key], name=f"Y_{mm}B2_{tag}"))
            if (m, EM_V1[f"{mm}B2"]) in E:
                extra.append(PDEntry(Composition({mm: 48, "B": 95}), E[(m, EM_V1[f"{mm}B2"])], name=EM_V1[f"{mm}B2"]))
            if len(ents) < 3:
                continue
            pdg = PhaseDiagram(ents)
            for en in ents:
                mid = en.name.split("_")[1] if en.name.startswith("ref_") else None
                hull_rows.append({"model": m, "system": f"{mm}-B", "phase": en.name, "dHf_eV_atom": pdg.get_form_energy_per_atom(en),
                                  "e_above_hull_eV_atom": pdg.get_e_above_hull(en),
                                  "MP_dHf_eV_atom": float(mp.loc[mid, "mp_formation_energy_eV_atom"]) if mid in mp.index else np.nan,
                                  "MP_e_hull_eV_atom": float(mp.loc[mid, "mp_e_above_hull_eV_atom"]) if mid in mp.index else np.nan})
            mu = E[(m, REFS["B"])] / natoms[REFS["B"]]
            for en in extra:
                k = 96 - int(en.composition["B"])
                d0 = en.energy + k * mu - E[(m, f"host_{mm}B2")]
                vac0_rows.append({"model": m, "end_member": f"{mm}B2", "cell": en.name, "k_vacancies": k,
                                  "dE_vs_stoich_plus_alphaB_eV": d0, "per_vacancy_eV": d0 / k,
                                  "e_above_own_hull_meV_atom": 1000 * pdg.get_e_above_hull(en, allow_negative=True),
                                  "verdict": ("B-deficient cell BELOW the MB2 + alpha-B tie line at 0 K (stoichiometric MB2 not the 0 K "
                                              "ground state)") if d0 < 0 else "stoichiometric MB2 stable vs this cell at 0 K"})
        for sid, mm in SYS.items():
            p = model_params(m, sid, E, natoms)
            if p is None:
                continue
            par_rows.append({"model": m, "system": sid, "L0_eV_fu": p["rk"][0], "L1_eV_fu": p["rk"][1],
                             "dH_mix_x0.5_meV_atom": p["mix_meV_atom"][0.5], "dEv_ZrB2_eV": p["dEv_zr"], "dEv_MB2_eV": p["dEv_m"],
                             "slope_eV_per_nn": p["slope_fit"][0], "W_ZrB2_eV_fu": p["w_zr"], "W_MB2_eV_fu": p["w_m"]})
            # ordering / segregation energies vs the model's own per-atom end members (exact for MLIPs)
            e0, e1 = E[(m, "host_ZrB2")] / 144, E[(m, f"host_{mm}B2")] / 144
            for tag, lab in (("qsqs1", "random (192-atom quasi-random)"), ("clay1", "metal-ordered alternating c-layers"), ("clay2", "segregated AABB slab")):
                key = (m, f"C_{sid}_x0.50_{tag}")
                if key in E:
                    ord_rows.append({"model": m, "system": sid, "type": "metal arrangement", "config": lab,
                                     "dH_mix_meV_atom": 1000 * (E[key] / 192 - 0.5 * e0 - 0.5 * e1)})
            ord_rows.append({"model": m, "system": sid, "type": "metal arrangement", "config": "BLADE 144-atom SQS",
                             "dH_mix_meV_atom": 1000 * (E[(m, f"host_{sid}_x0.50")] / 144 - 0.5 * e0 - 0.5 * e1)})
            for em in ("ZrB2", f"{mm}B2"):
                if em == "ZrB2" and sid == "ZrTa":
                    continue
                e1v = E[(m, EM_V1[em])] - E[(m, f"host_{em}")]
                for tag, k in (("k4_spread", 4), ("k4_cluster", 4), ("k8_spread", 8), ("k8_plane", 8)):
                    key = (m, f"Y_{em}_{tag}")
                    if key in E:
                        ord_rows.append({"model": m, "system": em, "type": "vacancy arrangement", "config": tag,
                                         "interaction_eV_per_vacancy": (E[key] - E[(m, f"host_{em}")] - k * e1v) / k})
            hc = (m, f"C_{sid}_x0.50_clay2")
            if hc in E:
                mu = E[(m, REFS["B"])] / natoms[REFS["B"]]
                for tag, n in (("ZrZr", 0), ("ZrM", 3), ("MM", 6)):
                    key = (m, f"S_{sid}_VB_{tag}")
                    if key in E:
                        ord_rows.append({"model": m, "system": sid, "type": "coupled metal-vacancy", "config": f"V_B in {tag} layer of AABB slab",
                                         "E_v_Brich_eV": E[key] - E[hc] + mu, "n_solute": n})
    pd.DataFrame(hull_rows).to_csv(OUT / "m06_hull_consistency.csv", index=False)
    pd.DataFrame(vac0_rows).to_csv(OUT / "m06_zeroK_Bdeficient_stability.csv", index=False)
    pd.DataFrame(ord_rows).to_csv(OUT / "m06_ordering_energies.csv", index=False)
    pd.DataFrame(par_rows).to_csv(OUT / "m06_model_parameters.csv", index=False)


def phase_maps(E: dict, natoms: dict) -> pd.DataFrame:
    rows = []
    for m in MODELS:
        for sid in SYS:
            for mixing in ("own", "DFT"):
                sysm = build_system(m, sid, E, natoms, mixing)
                if sysm is None:
                    continue
                for t in np.arange(1000.0, 3001.0, 250.0):
                    for bm in (2.0, 1.98, 1.95, 1.90, 1.80):
                        for x in np.round(np.arange(0.05, 0.951, 0.05), 2):
                            ph = sysm.equilibrium(x, bm, t)
                            al = ph[ph.phase == "AlB2"]
                            other = ph[ph.phase != "AlB2"]
                            rows.append({"model": m, "system": sid, "mixing": mixing, "T_K": t, "B_per_metal": bm, "x_overall": x,
                                         "class": classify(ph, sysm.compounds),
                                         "AlB2_x": ";".join(f"{v:.3f}" for v in al.x_M),
                                         "AlB2_B_per_metal": ";".join(f"{2 * (1 - v):.3f}" for v in al.y_Va),
                                         "AlB2_metal_fraction": ";".join(f"{v:.3f}" for v in al.metal_fraction),
                                         "other_phases": ";".join(other.phase),
                                         "other_metal_fraction": ";".join(f"{v:.3f}" for v in other.metal_fraction)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "m06_phase_maps.csv", index=False)
    return df


def sensitivity(E: dict, natoms: dict) -> pd.DataFrame:
    """Shift of the Group-V vacancy energy (d_m) needed to flip each model's assemblage at x = 0.5, 2000 K."""
    from s11_vacancy_sensitivity import perturbed

    rows = []
    for m in MODELS:
        for sid in SYS:
            base = build_system(m, sid, E, natoms)
            if base is None:
                continue
            for bm in (1.98, 1.95, 1.90):
                for d in np.round(np.arange(-1.4, 1.41, 0.1), 2):
                    s = perturbed(base, d_m=float(d))
                    ph = s.equilibrium(0.5, bm, 2000.0)
                    rows.append({"model": m, "system": sid, "B_per_metal": bm, "d_m_eV": d, "class": classify(ph, base.compounds)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "m06_sensitivity_dm.csv", index=False)
    return df


def plot_maps(df: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(exist_ok=True)
    classes = sorted(df["class"].unique())
    cmap = plt.get_cmap("tab10")
    cidx = {c: i for i, c in enumerate(classes)}
    for bm in (1.95, 1.90):
        fig, axes = plt.subplots(2, 5, figsize=(17, 6.5), sharex=True, sharey=True)
        for i, sid in enumerate(SYS):
            for j, m in enumerate(MODELS):
                ax = axes[i, j]
                g = df[(df.model == m) & (df.system == sid) & (df.mixing == "own") & (df.B_per_metal == bm)]
                if g.empty:
                    ax.set_axis_off()
                    continue
                ax.scatter(g.x_overall, g.T_K, c=[cmap(cidx[c] % 10) for c in g["class"]], s=28, marker="s")
                ax.set_title(f"{sid} | {m}", fontsize=8)
                if i == 1:
                    ax.set_xlabel("x (Nb or Ta)")
                if j == 0:
                    ax.set_ylabel("T (K)")
        handles = [plt.Line2D([], [], marker="s", ls="", color=cmap(cidx[c] % 10), label=c) for c in classes]
        fig.legend(handles=handles, loc="lower center", ncol=min(len(classes), 6), fontsize=7)
        fig.suptitle(f"Equilibrium phase assemblage, B/M = {bm} (each model's own energetics)", fontsize=10)
        fig.tight_layout(rect=(0, 0.06, 1, 0.96))
        fig.savefig(FIG / f"fig_m06_phase_maps_BM{bm:.2f}.png", dpi=180)


def main() -> None:
    e, E = energies()
    natoms = dict(zip(e.structure_id, e.n_atoms))
    checks(e, E, natoms)
    maps = phase_maps(E, natoms)
    plot_maps(maps)
    sens = sensitivity(E, natoms)
    pd.set_option("display.width", 240)
    print(pd.read_csv(OUT / "m06_model_parameters.csv").round(4).to_string(index=False))
    print(pd.read_csv(OUT / "m06_zeroK_Bdeficient_stability.csv").round(3).to_string(index=False))
    print(pd.read_csv(OUT / "m06_ordering_energies.csv").round(3).to_string(index=False))
    s = maps[(maps.mixing == "own") & (maps.x_overall == 0.5)]
    print(s.pivot_table(index=["system", "model", "B_per_metal"], columns="T_K", values="class", aggfunc="first").to_string())
    print(sens.pivot_table(index=["system", "model", "B_per_metal"], columns="d_m_eV", values="class", aggfunc="first").iloc[:, ::4].to_string())


if __name__ == "__main__":
    main()
