"""Six-system nonstoichiometric thermodynamics (closed and open systems), all five MLIPs.

Free energy per formula unit of the AlB2 phase (Zr|Ti|Hf, M)(B, Va)2 (thermo_framework.AlB2Params):
  - end members and RK (L0, L1) stoichiometric mixing from each model's own x = 0.25/0.5/0.75 hosts ("own"), or
    Tang 2025 LDA L0*, L1* ("DFT mixing" sensitivity); ideal metal-sublattice entropy
  - dilute B-vacancy energies of both end members, local-environment slope from the x = 0.5 SQS vacancy series,
    mean-field vacancy-vacancy interaction W y^2 per end member (k = 2/4/8 cells), ideal vacancy entropy
  - optional vibrational correction to the mixing: Tang QHA L0(T) - L0(0) (Table 4) added as x(1-x) dL0_vib(T)
    (scenario "+vib"); electronic free energy: no data -> not included (stated)
Line compounds (each model's own energies): alpha-B, both metals, all stable metal-rich borides of both metals,
ordered ZrMB4-type c-layer cell, vacancy-ordered cells (k4 cluster, k8 plane, and every enumerated ordered MB2-x cell
of TaB2/NbB2 from Stage 4).
Validity: AlB2 vacancy fraction y <= 0.125 (DFT stability range of TaB2-x); beyond that the deficit must form a real
phase. Temperatures 1500-3000 K (no liquid; solidus ~3200-3500 K).
Scenarios: own | DFT-anch (Ta systems: TaB2 dilute V_B = -0.389 eV) | own+vib | DFT-mix (Tang mixing).
Closed system: x in {0.25, 0.5, 0.75}, B/M = 2(1 - delta), delta 0-20 %.
Open system: x = 0.5, boron chemical potential mu_B from the alpha-B limit downward (reservoir), grand potential
per metal atom, lower convex hull in x with all phases.
Outputs (final_results/supplementary_data): x01_closed_scan.csv, x01_critical.csv, x01_open_scan.csv, x01_params.csv
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
VAL = FR.parent
WF = VAL / "validation_workflow"
MB = WF / "mlip_benchmark"
sys.path.insert(0, str(WF / "scripts"))
sys.path.insert(0, str(WF / "phase_fraction_analysis" / "scripts"))
from thermo_framework import AlB2Params, System  # noqa: E402

SD = FR / "supplementary_data"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
DFT_TAB2_VB = -0.389
Y_CAP = 0.125
TANG = {"TiNb": [-20.12, 4.249], "TiTa": [-42.21, 5.988], "ZrNb": [-2.515, 2.418], "ZrTa": [-9.306, 4.016],
        "HfNb": [-8.813, 5.557], "HfTa": [-14.12, 7.290]}  # kJ/mol MB2, our convention (2x - 1) with x = Group-V fraction
QHA = {"TiNb": (-2.8557e4, -1.0937e-1, -2.5720e-6, 1.3810e-2), "TiTa": (-5.3691e4, -2.2006e-1, -5.5094e-6, 2.8034e-2),
       "ZrNb": (-8.8268e3, -6.8045e-2, -1.6066e-6, 8.5529e-3), "ZrTa": (-2.1210e4, -2.4561e-1, -6.3207e-6, 3.1318e-2),
       "HfNb": (-1.6953e4, -1.2778e-1, -3.2267e-6, 1.6228e-2), "HfTa": (-2.6841e4, -2.4549e-1, -6.3012e-6, 3.1303e-2)}
REF_METAL = {"Zr": "ref_mp-131_Zr", "Ti": "ref_mp-72_Ti", "Hf": "ref_mp-103_Hf", "Nb": "ref_mp-75_Nb", "Ta": "ref_mp-50_Ta"}
BORIDES = {"Nb": ["ref_mp-1080021_Nb2B3", "ref_mp-10255_Nb3B4", "ref_mp-1102394_Nb5B6", "ref_mp-2580_NbB"],
           "Ta": ["ref_mp-10142_Ta3B4", "ref_mp-28629_Ta5B6", "ref_mp-1097_TaB"],
           "Ti": ["ref_mp-7857_TiB", "ref_mp-1025170_Ti3B4"], "Zr": [], "Hf": []}


def spec(sid: str) -> dict:
    host, m = sid[:2], sid[2:]
    zr = host == "Zr"
    return {"host": host, "m": m,
            "em_host": f"host_{host}B2", "em_m": f"host_{m}B2",
            "vh": "V_ZrB2_x0.00_VB_n0a" if zr else f"V_{host}B2_VB",
            "vm": {"Nb": "V_ZrNb_x1.00_VB_n6a", "Ta": "V_ZrTa_x1.00_VB_n6a"}[m],
            "x": {x: f"host_{sid}_x{x:.2f}" for x in (0.25, 0.5, 0.75)},
            "env": [f"V_{sid}_x0.50_VB_n{n}a" for n in range(7)],
            "clay1": f"C_{sid}_x0.50_clay1"}


def energies() -> tuple[dict, dict]:
    mans = [pd.read_csv(MB / "pilot_manifest.csv"), pd.read_csv(MB / "stage2_manifest.csv"), pd.read_csv(SD / "stage4_manifest.csv")]
    man = pd.concat(mans, ignore_index=True).drop_duplicates("structure_id")
    E, nat = {}, {}
    for r in man.itertuples():
        nat[r.structure_id] = int(r.n_atoms)
        if r.role == "host" and pd.notna(getattr(r, "orb_energy_eV", np.nan)):
            E[("orb", r.structure_id)] = float(r.orb_energy_eV)  # BLADE's own ORB host energy (as Stage 1-3)
        elif pd.notna(getattr(r, "orb_energy_eV", np.nan)):
            E[("orb", r.structure_id)] = float(r.orb_energy_eV)
        for m in MODELS:
            p = MB / "results" / m / r.structure_id / "result.json"
            if p.exists() and (m, r.structure_id) not in E:
                j = json.loads(p.read_text())
                if j.get("status") == "ok":
                    E[(m, r.structure_id)] = float(j["energy_eV"])
    return E, nat


def counts(sid: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", sid.split("_")[-1])}


def params(model: str, sid: str, E: dict, nat: dict, mixing: str = "own") -> dict | None:
    s = spec(sid)
    need = [s["em_host"], s["em_m"], *s["x"].values(), "ref_mp-160_B", s["vh"], s["vm"]]
    if any((model, k) not in E for k in need):
        return None
    e_b = E[(model, "ref_mp-160_B")] / nat["ref_mp-160_B"]
    efu0, efu1 = E[(model, s["em_host"])] / 48, E[(model, s["em_m"])] / 48
    xs = np.array([0.25, 0.5, 0.75])
    mix = np.array([E[(model, s["x"][x])] / 48 - (1 - x) * efu0 - x * efu1 for x in xs])
    a = np.vstack([xs * (1 - xs) * (2 * xs - 1) ** k for k in range(2)]).T
    rk = list(np.linalg.lstsq(a, mix, rcond=None)[0]) if mixing == "own" else [v / 96.485 for v in TANG[sid]]
    dev = {k: E[(model, s[v])] - E[(model, s[em])] + e_b for k, v, em in (("h", "vh", "em_host"), ("m", "vm", "em_m"))}
    env = [(n, E[(model, k)] - E[(model, s["x"][0.5])] + e_b) for n, k in enumerate(s["env"]) if (model, k) in E]
    slope = np.polyfit([p[0] for p in env], [p[1] for p in env], 1) if len(env) >= 3 else None
    w = {}
    for tag, em, v1 in (("h", f"{s['host']}B2", "vh"), ("m", f"{s['m']}B2", "vm")):
        e1 = E[(model, s[v1])] - E[(model, s["em_host" if tag == "h" else "em_m"])]
        pts = []
        for k in (2, 4, 8):
            key = (model, f"Y_{em}_k{k}_spread")
            if key in E:
                pts.append((k / 96, (E[key] - E[(model, s["em_host" if tag == "h" else "em_m"])] - k * e1) / 48))
        w[tag] = float(sum(y**2 * ex for y, ex in pts) / sum(y**4 for y, _ in pts)) if pts else np.nan
    return {"e_b": e_b, "efu0": efu0, "efu1": efu1, "rk": rk, "dEv_h": dev["h"], "dEv_m": dev["m"], "slope": slope, "w_h": w["h"], "w_m": w["m"],
            "mix_x0.5_meV_atom": 1000 * mix[1] / 3, "n_env": len(env)}


def build(model: str, sid: str, E: dict, nat: dict, scenario: str = "own") -> System | None:
    s = spec(sid)
    p = params(model, sid, E, nat, "DFT" if scenario == "DFT-mix" else "own")
    if p is None or p["slope"] is None or np.isnan(p["w_h"]) or np.isnan(p["w_m"]):
        return None
    dEv_m = DFT_TAB2_VB if (scenario.startswith("DFT-anch") and s["m"] == "Ta") else p["dEv_m"]
    shift = dEv_m - p["dEv_m"]
    f = p["slope"]
    g_extra = None
    if scenario.endswith("+vib"):
        a_, b_, c_, d_ = QHA[sid]
        l0 = lambda t: (a_ + b_ * t**2 + c_ * t**3 + d_ * t**2 * np.log(max(t, 1.0))) / 1000 / 96.485
        g_extra = lambda x, t, l0=l0: x * (1 - x) * (l0(t) - l0(1.0))
    alb2 = AlB2Params(g_zrb=p["efu0"], g_mb=p["efu1"], dEv_zr=p["dEv_h"], dEv_m=dEv_m, e_b=p["e_b"], rk=p["rk"],
                      ev_env=lambda n, f=f, sh=shift: float(np.polyval(f, n)) + sh * n / 6, w_vv_zr=p["w_h"], w_vv_m=p["w_m"], g_extra=g_extra)
    rows = []
    for ref in ["ref_mp-160_B", REF_METAL[s["host"]], REF_METAL[s["m"]]] + BORIDES[s["host"]] + BORIDES[s["m"]]:
        if (model, ref) not in E:
            continue
        c = counts(ref)
        nfu = nat[ref] / sum(c.values())
        rows.append({"name": ref.split("_", 1)[1], "n_Zr": c.get(s["host"], 0), "n_M": c.get(s["m"], 0), "n_B": c.get("B", 0),
                     "G": E[(model, ref)] / nfu, "kind": "reference"})
    if (model, s["clay1"]) in E:
        rows.append({"name": f"ordered {s['host']}{s['m']}B4", "n_Zr": 1, "n_M": 1, "n_B": 4, "G": E[(model, s["clay1"])] / 32, "kind": "ordered_metal"})
    for em, nh, nm in ((f"{s['host']}B2", 48, 0), (f"{s['m']}B2", 0, 48)):
        for tag, k in (("k4_cluster", 4), ("k8_plane", 8)):
            if (model, f"Y_{em}_{tag}") in E:
                rows.append({"name": f"{em}-{k}V ({tag})", "n_Zr": nh, "n_M": nm, "n_B": 96 - k, "G": E[(model, f"Y_{em}_{tag}")], "kind": "ordered_vacancy"})
    # enumerated ordered MB2-x cells (Stage 4), only for the Group-V end member of this system
    for (mm, sid_), en in E.items():
        if mm != model or not sid_.startswith(f"O_{s['m']}B2_"):
            continue
        n = int(re.search(r"_N(\d+)k(\d+)_", sid_).group(1))
        k = int(re.search(r"_N(\d+)k(\d+)_", sid_).group(2))
        corr = shift * k if (scenario.startswith("DFT-anch") and s["m"] == "Ta") else 0.0  # same DFT anchoring per vacancy
        rows.append({"name": sid_, "n_Zr": 0, "n_M": n, "n_B": 2 * n - k, "G": en + corr, "kind": "ordered_vacancy"})
    ys = np.r_[0.0, np.geomspace(1e-5, Y_CAP, 45)]
    return System(metal=s["m"], alb2=alb2, compounds=pd.DataFrame(rows), x_grid=np.linspace(0, 1, 101), y_grid=ys)


DX_MERGE, DBM_MERGE = 0.05, 0.05  # AlB2 grid points closer than this in x and B/M are one phase (LP grid discretisation)


def merge_alb2(ph: pd.DataFrame) -> pd.DataFrame:
    """Merge AlB2 'phases' that are neighbouring grid points of one solution (|dx| <= 0.05 and |dB/M| <= 0.05)."""
    al = ph[ph.phase == "AlB2"].sort_values("metal_fraction", ascending=False)
    if len(al) < 2:
        return ph
    groups: list[list] = []
    for i, r in al.iterrows():
        for g in groups:
            x0, y0 = g[0][1].x_M, g[0][1].y_Va
            if abs(r.x_M - x0) <= DX_MERGE and abs(2 * (r.y_Va - y0)) <= DBM_MERGE:
                g.append((i, r))
                break
        else:
            groups.append([(i, r)])
    rows = []
    for g in groups:
        w = np.array([r.metal_fraction for _, r in g])
        base = g[0][1].copy()
        if w.sum() > 0:
            base["x_M"] = float(np.dot(w, [r.x_M for _, r in g]) / w.sum())
            base["y_Va"] = float(np.dot(w, [r.y_Va for _, r in g]) / w.sum())
        base["metal_fraction"] = float(w.sum())
        rows.append(base)
    return pd.concat([ph[ph.phase != "AlB2"], pd.DataFrame(rows)], ignore_index=True)


LP_LOG: list[dict] = []  # solver retries / failures, written to supplementary_data/x01_lp_log.csv


def solve(sysm: System, x: float, b_per_metal: float, t: float, tag: str = "") -> pd.DataFrame:
    """System.equilibrium with HiGHS (default), retried with dual simplex and interior point if the solve fails.

    Mass balance and the objective are identical for every method; only the numerical path differs. States that
    fail with all three methods are returned as a one-row frame with 'status' and are reported, never imputed.
    """
    msgs = []
    for method in ("highs", "highs-ds", "highs-ipm"):
        ph = sysm.equilibrium(x, b_per_metal, t, method=method)
        if "phase" in ph.columns:
            if msgs:
                LP_LOG.append({"tag": tag, "x": x, "B_per_M": b_per_metal, "T_K": t, "outcome": f"solved with {method}", "messages": " | ".join(msgs)})
            return ph
        msgs.append(f"{method}: {ph.status.iloc[0]}")
    LP_LOG.append({"tag": tag, "x": x, "B_per_M": b_per_metal, "T_K": t, "outcome": "failed", "messages": " | ".join(msgs)})
    return pd.DataFrame([{"status": " | ".join(msgs)}])


def summarise(ph: pd.DataFrame, comp: pd.DataFrame) -> dict:
    if "phase" not in ph.columns:  # LP failed with every method
        return {"class": "LP failed", "second_metal_frac": np.nan, "second_phase": "", "main_x": np.nan, "main_B_per_M": np.nan,
                "second_x": np.nan, "second_B_per_M": np.nan, "y_max": np.nan}
    kind = dict(zip(comp.name, comp.kind))
    ph = ph[~ph.phase.astype(str).str.contains(r"_B$|mp-160_B|^B ")].copy()
    ph = ph[ph.x_M.notna() | (ph.phase == "AlB2")]
    ph = merge_alb2(ph)
    tot = ph.metal_fraction.sum()
    al = ph[ph.phase == "AlB2"].sort_values("metal_fraction", ascending=False)
    if tot <= 0 or al.empty:
        main = ph.sort_values("metal_fraction", ascending=False).iloc[0] if len(ph) else None
    else:
        main = al.iloc[0]
    others = ph.drop(index=main.name) if main is not None else ph
    second = others.sort_values("metal_fraction", ascending=False).iloc[0] if len(others) else None
    if second is None:
        ident = "homogeneous"
    elif second.phase == "AlB2":
        ident = "two AlB2"
    else:
        ident = {"ordered_vacancy": "vacancy-ordered MB2-x", "ordered_metal": "ordered metal (ZrMB4-type)"}.get(kind.get(second.phase), "metal-rich boride")
    return {"class": ident, "second_metal_frac": float(others.metal_fraction.sum() / tot) if tot > 0 else np.nan,
            "second_phase": "" if second is None else str(second.phase), "main_x": float(main.x_M) if main is not None else np.nan,
            "main_B_per_M": float(2 * (1 - main.y_Va)) if main is not None and main.phase == "AlB2" else np.nan,
            "second_x": np.nan if second is None else float(second.x_M),
            "second_B_per_M": np.nan if second is None or second.phase != "AlB2" else float(2 * (1 - second.y_Va)),
            "y_max": float(al.y_Va.max()) if len(al) else np.nan}


def open_scan(sysm: System, t: float, dmu_list) -> list[dict]:
    """Grand potential per metal atom at fixed mu_B = mu_B(alpha-B) - dmu; lower hull over x; phases at x = 0.5."""
    a = sysm.alb2
    xs = np.linspace(0, 1, 201)
    ys = sysm.y_grid
    xx, yy = np.meshgrid(xs, ys, indexing="ij")
    out = []
    for dmu in dmu_list:
        mu = a.e_b - dmu
        om = a.gibbs(xx, yy, t) - mu * 2 * (1 - yy)
        k = om.argmin(axis=1)
        pts = [(x, om[i, k[i]], "AlB2", ys[k[i]]) for i, x in enumerate(xs)]
        for r in sysm.compounds.itertuples():
            nm = r.n_Zr + r.n_M
            if nm:
                pts.append((r.n_M / nm, (r.G - mu * r.n_B) / nm, r.name, np.nan))
        pts.sort(key=lambda p: (p[0], p[1]))
        hull: list = []
        for p in pts:
            if hull and abs(hull[-1][0] - p[0]) < 1e-12:
                if p[1] < hull[-1][1]:
                    hull[-1] = p
                    while len(hull) >= 3 and (hull[-2][0] - hull[-3][0]) * (hull[-1][1] - hull[-3][1]) - (hull[-2][1] - hull[-3][1]) * (hull[-1][0] - hull[-3][0]) <= 0:
                        hull.pop(-2)
                continue
            while len(hull) >= 2 and (hull[-1][0] - hull[-2][0]) * (p[1] - hull[-2][1]) - (hull[-1][1] - hull[-2][1]) * (p[0] - hull[-2][0]) <= 0:
                hull.pop()
            hull.append(p)
        seg = None
        for i in range(len(hull) - 1):
            if hull[i][0] <= 0.5 <= hull[i + 1][0]:
                seg = [hull[i]] if abs(hull[i][0] - 0.5) < 1e-9 else ([hull[i + 1]] if abs(hull[i + 1][0] - 0.5) < 1e-9 else [hull[i], hull[i + 1]])
                break
        if seg is None:
            continue
        al_single = len(seg) == 1 and seg[0][2] == "AlB2"
        two_al = len(seg) == 2 and seg[0][2] == "AlB2" and seg[1][2] == "AlB2" and abs(seg[1][0] - seg[0][0]) > 0.05
        cls = "homogeneous" if al_single else ("two AlB2" if two_al else ("phase change" if len(seg) == 1 else "AlB2 + other"))
        out.append({"dmu_B_eV": dmu, "class": cls, "phases": ";".join(f"{p[2]}[x={p[0]:.2f},y={p[3]:.3f}]" for p in seg),
                    "y_at_x0.5": float(seg[0][3]) if al_single else np.nan})
    return out


def export_tables(E: dict, nat: dict) -> None:
    """Per-structure tables for figures: mixing enthalpies, vacancy vs environment, vacancy interactions, ordering."""
    mix, env, vv, ordr, emv = [], [], [], [], []
    for sid in ["TiNb", "TiTa", "ZrNb", "ZrTa", "HfNb", "HfTa"]:
        s = spec(sid)
        for model in MODELS:
            if any((model, k) not in E for k in (s["em_host"], s["em_m"], "ref_mp-160_B")):
                continue
            e0, e1 = E[(model, s["em_host"])] / 144, E[(model, s["em_m"])] / 144
            eB = E[(model, "ref_mp-160_B")] / nat["ref_mp-160_B"]
            for x, k in s["x"].items():
                if (model, k) in E:
                    mix.append({"system": sid, "model": model, "x": x, "dH_mix_meV_atom": 1000 * (E[(model, k)] / 144 - (1 - x) * e0 - x * e1)})
            for n, k in enumerate(s["env"]):
                if (model, k) in E and (model, s["x"][0.5]) in E:
                    env.append({"system": sid, "model": model, "n_solute": n, "E_v_Brich_eV": E[(model, k)] - E[(model, s["x"][0.5])] + eB})
            for tag in ("qsqs1", "clay1", "clay2"):
                k = f"C_{sid}_x0.50_{tag}"
                if (model, k) in E:
                    ordr.append({"system": sid, "model": model, "config": tag, "dH_mix_meV_atom": 1000 * (E[(model, k)] / 192 - 0.5 * e0 - 0.5 * e1)})
            for lab, n in (("HH", 0), ("HM", 3), ("MM", 6), ("ZrZr", 0), ("ZrM", 3)):
                k = f"S_{sid}_VB_{lab}"
                hc = (model, f"C_{sid}_x0.50_clay2")
                if (model, k) in E and hc in E:
                    ordr.append({"system": sid, "model": model, "config": f"V_B in {lab} layer", "n_solute": n, "E_v_Brich_eV": E[(model, k)] - E[hc] + eB})
        for em, vkey, hkey in ((f"{s['host']}B2", s["vh"], s["em_host"]), (f"{s['m']}B2", s["vm"], s["em_m"])):
            for model in MODELS:
                if (model, vkey) not in E or (model, hkey) not in E:
                    continue
                eB = E[(model, "ref_mp-160_B")] / nat["ref_mp-160_B"]
                e1 = E[(model, vkey)] - E[(model, hkey)]
                emv.append({"end_member": em, "model": model, "E_v_Brich_eV": e1 + eB})
                for tag in ("k2_spread", "k4_spread", "k8_spread", "k4_cluster", "k8_plane"):
                    key = (model, f"Y_{em}_{tag}")
                    if key in E:
                        k = int(tag[1])
                        vv.append({"end_member": em, "model": model, "arrangement": tag, "k": k,
                                   "interaction_eV_per_vacancy": (E[key] - E[(model, hkey)] - k * e1) / k})
    pd.DataFrame(mix).drop_duplicates().to_csv(SD / "x01_mix_points.csv", index=False)
    pd.DataFrame(env).drop_duplicates().to_csv(SD / "x01_vacancy_env.csv", index=False)
    pd.DataFrame(vv).drop_duplicates().to_csv(SD / "x01_vacancy_interactions.csv", index=False)
    pd.DataFrame(ordr).drop_duplicates().to_csv(SD / "x01_ordering.csv", index=False)
    pd.DataFrame(emv).drop_duplicates().to_csv(SD / "x01_endmember_vacancy.csv", index=False)


def main() -> None:
    E, nat = energies()
    export_tables(E, nat)
    if len(sys.argv) > 1 and sys.argv[1] == "tables":
        return
    systems = ["TiNb", "TiTa", "ZrNb", "ZrTa", "HfNb", "HfTa"]
    prow, closed, crit, openr = [], [], [], []
    deltas = [0.0, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.125, 0.15, 0.175, 0.20]
    for sid in systems:
        for model in MODELS:
            p = params(model, sid, E, nat)
            if p is None:
                prow.append({"system": sid, "model": model, "status": "missing inputs"})
                continue
            prow.append({"system": sid, "model": model, "status": "ok", "L0_eV_fu": p["rk"][0], "L1_eV_fu": p["rk"][1], "dH_mix_x0.5_meV_atom": p["mix_x0.5_meV_atom"],
                         "dEv_hostB2_eV": p["dEv_h"], "dEv_MB2_eV": p["dEv_m"], "slope_eV_per_nn": None if p["slope"] is None else p["slope"][0],
                         "n_env_points": p["n_env"], "W_hostB2_eV_fu": p["w_h"], "W_MB2_eV_fu": p["w_m"]})
            scens = ["own", "own+vib", "DFT-mix"] + (["DFT-anch", "DFT-anch+vib"] if sid.endswith("Ta") else [])
            for scen in scens:
                sysm = build(model, sid, E, nat, scen)
                if sysm is None:
                    continue
                for t in (1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0):
                    for x in (0.25, 0.5, 0.75):
                        recs = []
                        for d in deltas:
                            r = summarise(solve(sysm, x, 2 * (1 - d), t, f"{sid}/{model}/{scen}"), sysm.compounds)
                            recs.append({"system": sid, "model": model, "scenario": scen, "T_K": t, "x": x, "delta": d, **r})
                        closed += recs
                        g = pd.DataFrame(recs)
                        first = lambda c: float(g[c].delta.iloc[0]) if c.any() else np.nan
                        multi = g["class"] != "homogeneous"
                        d5 = first(g.second_metal_frac >= 0.05)
                        at5 = g[g.delta == d5].iloc[0] if d5 == d5 else None
                        crit.append({"system": sid, "model": model, "scenario": scen, "T_K": t, "x": x,
                                     "delta_onset": float(g[multi].delta.iloc[0]) if multi.any() else np.nan,
                                     "delta_5": d5, "delta_30": first(g.second_metal_frac >= 0.30),
                                     "class_at_5": None if at5 is None else at5["class"], "second_phase_at_5": None if at5 is None else at5.second_phase,
                                     "main_x_at_5": None if at5 is None else at5.main_x, "second_x_at_5": None if at5 is None else at5.second_x,
                                     "second_BM_at_5": None if at5 is None else at5.second_B_per_M,
                                     "frac_at_delta_0.05": float(g[g.delta == 0.05].second_metal_frac.iloc[0])})
                if scen in ("own", "DFT-anch"):
                    for t in (2000.0, 3000.0):
                        for r in open_scan(sysm, t, np.round(np.arange(0.0, 1.01, 0.05), 3)):
                            openr.append({"system": sid, "model": model, "scenario": scen, "T_K": t, **r})
            print(sid, model, "done", flush=True)
    pd.DataFrame(prow).to_csv(SD / "x01_params.csv", index=False)
    pd.DataFrame(closed).to_csv(SD / "x01_closed_scan.csv", index=False)
    pd.DataFrame(crit).to_csv(SD / "x01_critical.csv", index=False)
    pd.DataFrame(openr).to_csv(SD / "x01_open_scan.csv", index=False)
    pd.DataFrame(LP_LOG, columns=["tag", "x", "B_per_M", "T_K", "outcome", "messages"]).to_csv(SD / "x01_lp_log.csv", index=False)
    n_fail = sum(r["outcome"] == "failed" for r in LP_LOG)
    print(f"LP: {len(LP_LOG) - n_fail} states needed a fallback solver, {n_fail} failed with all methods", flush=True)


if __name__ == "__main__":
    main()
