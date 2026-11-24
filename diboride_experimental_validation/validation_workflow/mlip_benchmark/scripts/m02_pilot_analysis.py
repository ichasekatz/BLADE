"""Stage 1 pilot analysis: ORB vs MACE-MPA vs MACE-OMAT vs MatterSim vs SevenNet on identical structures.

Everything is built WITHIN each model (no cross-model total-energy comparison):
  mu_B(B-rich)   = E_model(alpha-B, mp-160) / 12
  mu_B(B-poor)   end members: max over the model's own metal-rich phases P (n_M, n_B) and the metal:
                 mu_B >= (E_P - n_M e_host) / (n_B - 2 n_M),  e_host = E_host / 48 (energy per MB2)
  E_f(V_B)       = E(host - B) - E(host) + mu_B          (eV per vacancy)
  E_v            = E(host - B) - E(host) + mu_B(B-rich)  (used for environment trends / reciprocal energy)
  dE_rec         = E_v(MB2) - E_v(ZrB2)                  (independent of mu_B)
  dHf(MB2)       = [E_host - 48 E_M - 96 E_B] / 144 per atom, elements = model's own Zr/bcc-Nb/bcc-Ta/alpha-B
  e_hull(MB2)    model's own M-B convex hull (MB2 host + metal-rich borides + metal + alpha-B)
Benchmarks: experimental dHf (Ness 2025 OMSC, impurity-corrected; processed_data/thermochemical_data.csv),
MP (GGA) formation energies of the same MP structures (inputs/competitors/manifest.csv), Otani 2009 end-member
lattice constants (digitized, raw_data/digitized_figures/otani2009_fig1.csv).

Outputs: mlip_benchmark/results/m02_*.csv, mlip_benchmark/figures/fig_m02_*.png
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

MB = Path(__file__).resolve().parents[1]
WF = MB.parent
VAL = WF.parent
OUT = MB / "results"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
METAL = {"ZrB2": "Zr", "NbB2": "Nb", "TaB2": "Ta"}
MREF = {"Zr": "ref_mp-131_Zr", "Nb": "ref_mp-75_Nb", "Ta": "ref_mp-50_Ta"}
MRICH = {"NbB2": ["ref_mp-1080021_Nb2B3", "ref_mp-10255_Nb3B4", "ref_mp-1102394_Nb5B6", "ref_mp-2580_NbB", "ref_mp-75_Nb"],
         "TaB2": ["ref_mp-10142_Ta3B4", "ref_mp-28629_Ta5B6", "ref_mp-1097_TaB", "ref_mp-50_Ta"], "ZrB2": ["ref_mp-131_Zr"]}
KJMOL_PER_EV = 96.485


def comp_of(sid: str, n_atoms: int) -> tuple[int, int]:
    """(n_metal, n_B) of a reference cell from its formula in the id."""
    import re

    formula = sid.split("_")[-1]
    c = {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}
    nfu = n_atoms / sum(c.values())
    nm = sum(v for k, v in c.items() if k != "B") * nfu
    return int(round(nm)), int(round(c.get("B", 0) * nfu))


def load() -> pd.DataFrame:
    man = pd.read_csv(MB / "pilot_manifest.csv")
    rows = []
    for r in man.itertuples():
        rows.append({"model": "orb", "structure_id": r.structure_id, "energy_eV": r.orb_energy_eV, "converged": True, "status": "ok"})
        for m in MODELS[1:]:
            p = OUT / m / r.structure_id / "result.json"
            if p.exists():
                j = json.loads(p.read_text())
                rows.append({"model": m, "structure_id": r.structure_id, "energy_eV": j["energy_eV"], "converged": j["converged"],
                             "status": j["status"], "n_steps": j["n_steps"], "max_disp_A": j["max_atomic_displacement_A"],
                             "vol_change_pct": 100 * (j["volume_A3"] / j["volume_start_A3"] - 1), "cell_lengths": j["cell_lengths"],
                             "cell_angles": j["cell_angles"]})
    e = pd.DataFrame(rows).merge(man[["structure_id", "role", "system", "x_solute", "n_solute", "n_atoms"]], on="structure_id")
    return e


def orb_cell(sid: str) -> list[float] | None:
    """ORB relaxed cell lengths for hosts (BLADE CONTCAR) for the lattice benchmark."""
    man = pd.read_csv(MB / "pilot_manifest.csv").set_index("structure_id")
    from ase.io import read

    at = read(man.loc[sid, "start_geometry"], format="vasp")
    return [float(v) for v in at.cell.lengths()]


def main() -> None:
    e = load()
    E = {(r.model, r.structure_id): r.energy_eV for r in e.itertuples()}
    natoms = dict(zip(e.structure_id, e.n_atoms))
    vac_rows, em_rows = [], []
    thermo = pd.read_csv(VAL / "processed_data" / "thermochemical_data.csv")
    exp_hf = {r.compound_or_system: (r.value, r.uncertainty) for r in thermo[thermo.quantity.str.contains("standard enthalpy of formation", na=False)
                                                                              & thermo.compound_or_system.isin(METAL)].itertuples()}
    mpman = pd.read_csv(WF / "inputs" / "competitors" / "manifest.csv").set_index("formula")
    otani = pd.read_csv(VAL / "raw_data" / "digitized_figures" / "otani2009_fig1.csv")
    ot = {}
    for host, sysname, x in (("ZrB2", "ZrB2-NbB2", 0.0), ("NbB2", "ZrB2-NbB2", 1.0), ("TaB2", "ZrB2-TaB2", 1.0)):
        g = otani[(otani.system == sysname) & np.isclose(otani.x_solute, x)]
        ot[host] = (float(g[g.series == "a lattice constant"].y_value.iloc[0]) * 10, float(g[g.series == "c lattice constant"].y_value.iloc[0]) * 10)
    for m in MODELS:
        if (m, "ref_mp-160_B") not in E or (m, "host_ZrB2") not in E:
            continue
        mu_rich = E[(m, "ref_mp-160_B")] / natoms["ref_mp-160_B"]
        ev = {}
        for host in METAL:
            hid = f"host_{host}"
            if (m, hid) not in E:
                continue
            e_fu = E[(m, hid)] / 48
            # B-poor limit from the model's own metal-rich phases
            lows = []
            for ref in MRICH[host]:
                if (m, ref) not in E:
                    continue
                nm, nb = comp_of(ref, natoms[ref])
                lows.append(((E[(m, ref)] - nm * e_fu) / (nb - 2 * nm), ref.split("_")[-1]))
            mu_poor, poor_by = max(lows) if lows else (np.nan, "")
            # own hull: is MB2 stable against the metal-rich phases + B ?  (e_hull of MB2 on the model's M-B hull)
            from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
            from pymatgen.core import Composition

            el = METAL[host]
            ents = [PDEntry(Composition({el: 48, "B": 96}), E[(m, hid)], name=hid), PDEntry(Composition({"B": 12}), E[(m, "ref_mp-160_B")], name="B")]
            for ref in MRICH[host]:
                if (m, ref) in E:
                    nm, nb = comp_of(ref, natoms[ref])
                    ents.append(PDEntry(Composition({el: nm, "B": nb} if nb else {el: nm}), E[(m, ref)], name=ref))
            pdg = PhaseDiagram(ents)
            ehull = pdg.get_e_above_hull(ents[0])
            dhf = (E[(m, hid)] - 48 * E[(m, MREF[el])] / natoms[MREF[el]] - 96 * mu_rich) / 144
            vid = {"ZrB2": "V_ZrB2_x0.00_VB_n0a", "NbB2": "V_ZrNb_x1.00_VB_n6a", "TaB2": "V_ZrTa_x1.00_VB_n6a"}[host]
            ev_h = E.get((m, vid), np.nan) - E[(m, hid)] + mu_rich
            ev[host] = ev_h
            cell = orb_cell(hid) if m == "orb" else e[(e.model == m) & (e.structure_id == hid)].cell_lengths.iloc[0]
            a, c = cell[0] / 4, cell[2] / 3
            em_rows.append({"model": m, "end_member": host, "dHf_eV_atom": dhf, "dHf_kJ_mol_fu": dhf * 3 * KJMOL_PER_EV,
                            "exp_dHf_kJ_mol (Ness2025 OMSC)": exp_hf.get(host, (np.nan,))[0], "exp_unc_2SE": exp_hf.get(host, (np.nan, np.nan))[1],
                            "MP_GGA_dHf_eV_atom": float(mpman.loc[host, "mp_formation_energy_eV_atom"]) if host in mpman.index else np.nan,
                            "e_hull_MB2_own_hull_eV_atom": ehull, "mu_B_window_eV": mu_rich - mu_poor, "B_poor_set_by": poor_by,
                            "E_f_VB_Brich_eV": ev_h, "E_f_VB_Bpoor_eV": E.get((m, vid), np.nan) - E[(m, hid)] + mu_poor,
                            "a_A": a, "c_A": c, "a_exp_Otani_A": ot[host][0], "c_exp_Otani_A": ot[host][1],
                            "a_err_pct": 100 * (a / ot[host][0] - 1), "c_err_pct": 100 * (c / ot[host][1] - 1)})
        for sysx, solute_host in (("ZrNb_x0.50", "NbB2"), ("ZrTa_x0.50", "TaB2")):
            hid = f"host_{sysx}"
            if (m, hid) not in E:
                continue
            sub = e[(e.model == m) & (e.role == "V_B") & (e.system == sysx)]
            for r in sub.itertuples():
                vac_rows.append({"model": m, "system": sysx[:4], "structure_id": r.structure_id, "n_solute": r.n_solute,
                                 "E_v_Brich_eV": r.energy_eV - E[(m, hid)] + mu_rich, "converged": r.converged})
        for host, v in ev.items():
            vac_rows.append({"model": m, "system": host, "structure_id": f"end member {host}", "n_solute": 6 if host != "ZrB2" else 0,
                             "E_v_Brich_eV": v, "converged": True})
    em = pd.DataFrame(em_rows)
    vac = pd.DataFrame(vac_rows)
    em.to_csv(OUT / "m02_end_member_benchmarks.csv", index=False)
    vac.to_csv(OUT / "m02_vacancy_energies.csv", index=False)

    # environment trends and cross-model agreement
    trend = []
    for (m, sysn), g in vac[vac.system.isin(["ZrNb", "ZrTa"])].groupby(["model", "system"]):
        slope, icpt = np.polyfit(g.n_solute, g.E_v_Brich_eV, 1)
        rep = g[g.n_solute == 3].E_v_Brich_eV
        host = "NbB2" if sysn == "ZrNb" else "TaB2"
        evm = vac[(vac.model == m) & (vac.system == host)].E_v_Brich_eV
        evz = vac[(vac.model == m) & (vac.system == "ZrB2")].E_v_Brich_eV
        trend.append({"model": m, "system": sysn, "slope_eV_per_solute_nn": slope, "range_n0_to_n6_eV": 6 * slope,
                      "replicate_n3_spread_eV": float(rep.max() - rep.min()) if len(rep) > 1 else np.nan,
                      "dE_rec_eV": float(evm.iloc[0] - evz.iloc[0]) if len(evm) and len(evz) else np.nan,
                      "n_points": len(g)})
    tr = pd.DataFrame(trend)
    # rank agreement with ORB (Spearman on the n-ordered environment energies)
    from scipy.stats import spearmanr

    rk = []
    for sysn in ("ZrNb", "ZrTa"):
        base = vac[(vac.model == "orb") & (vac.system == sysn)].set_index("structure_id").E_v_Brich_eV
        for m in MODELS[1:]:
            other = vac[(vac.model == m) & (vac.system == sysn)].set_index("structure_id").E_v_Brich_eV
            common = base.index.intersection(other.index)
            if len(common) >= 4:
                rho = spearmanr(base[common], other[common]).statistic
                rk.append({"model": m, "system": sysn, "spearman_vs_orb": rho,
                           "mean_abs_diff_vs_orb_eV": float((other[common] - base[common]).abs().mean()),
                           "mean_diff_vs_orb_eV": float((other[common] - base[common]).mean())})
    tr = tr.merge(pd.DataFrame(rk), on=["model", "system"], how="left") if rk else tr
    tr.to_csv(OUT / "m02_environment_trends.csv", index=False)

    # relaxation / physical-consistency checks
    chk = e[e.model != "orb"][["model", "structure_id", "role", "status", "converged", "n_steps", "max_disp_A", "vol_change_pct", "cell_angles"]].copy()
    chk["hexagonal_angles_ok"] = chk.cell_angles.apply(lambda a: isinstance(a, list) and abs(a[0] - 90) < 1 and abs(a[1] - 90) < 1 and abs(a[2] - 120) < 1)
    chk["hexagonal_angles_ok"] = chk.hexagonal_angles_ok.astype(object)
    chk.loc[chk.role == "reference", "hexagonal_angles_ok"] = None
    chk.drop(columns="cell_angles").to_csv(OUT / "m02_relaxation_checks.csv", index=False)

    pd.set_option("display.width", 250)
    print(em.round(3).to_string(index=False))
    print(tr.round(3).to_string(index=False))
    bad = chk[(chk.status != "ok") | (~chk.converged) | (chk.hexagonal_angles_ok == False)]  # noqa: E712
    print(f"relaxation issues: {len(bad)}")
    if len(bad):
        print(bad.to_string(index=False))
    plot(vac, em)


def plot(vac: pd.DataFrame, em: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (MB / "figures").mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, sysn in zip(axes[:2], ("ZrNb", "ZrTa")):
        for i, m in enumerate(MODELS):
            g = vac[(vac.model == m) & (vac.system == sysn)].sort_values("n_solute")
            if g.empty:
                continue
            ax.plot(g.n_solute + (i - 2) * 0.06, g.E_v_Brich_eV, "o-", label=m, ms=5)
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xlabel("Nb/Ta among the 6 metal neighbours of the B site")
        ax.set_ylabel("E_f(V_B) at B-rich limit (eV)")
        ax.set_title(f"x = 0.5 {sysn[:2]}B$_2$–{sysn[2:]}B$_2$ SQS (each model's own references)", fontsize=9)
    axes[0].legend(fontsize=7)
    ax = axes[2]
    hosts = ["ZrB2", "NbB2", "TaB2"]
    w = 0.16
    for i, m in enumerate(MODELS):
        g = em[em.model == m].set_index("end_member").reindex(hosts)
        ax.bar(np.arange(3) + (i - 2) * w, g.E_f_VB_Brich_eV, w, label=m)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(range(3), hosts)
    ax.set_ylabel("E_f(V_B), B-rich limit (eV)")
    ax.set_title("Pure diborides: < 0 means MB$_2$ would lose B next to elemental B", fontsize=9)
    fig.tight_layout()
    fig.savefig(MB / "figures" / "fig_m02_pilot_vacancy_energies.png", dpi=200)


if __name__ == "__main__":
    main()
