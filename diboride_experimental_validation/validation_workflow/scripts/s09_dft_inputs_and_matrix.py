"""Step 9: prepare (never launch) DFT and phonon inputs for the confirmation subset, and write the prioritised
calculation matrix (item 7).

DFT inputs: vasp_inputs/<tier>/<calc_id>/{POSCAR, INCAR, POTCAR.spec, README}. POTCARs are NOT generated
(licensed); POTCAR.spec lists the PAW symbols (MP defaults: Zr_sv, Nb_pv, Ta_pv, B). k-points via KSPACING so
every cell has the same reciprocal-space density. Starting geometries are the ORB-relaxed CONTCARs where they
exist (else the unrelaxed POSCAR / BLADE CONTCAR).

INCAR presets (PBE; same for every structure in a comparison):
  relax_fixed_cell  ISIF = 2 for point defects at the relaxed host lattice (dilute-defect convention)
  relax_full        ISIF = 3 for configurations / competitors (ENCUT 520 eV, 1.3x the B default, for Pulay stress)
  static_dos        ISMEAR = -5 tetrahedron, NEDOS 3001, LORBIT 11, denser KSPACING
Phonons: phonopy_disp.conf only (finite displacement 0.01 A), not run.

Cost estimates are rough: core-hours ~ 400 * (N_atoms/144)^2.5 for a 144-atom PBE relaxation (about 60 ionic
steps), static/DOS 0.15x, phonons (#displacements x static). They are for planning only.

Outputs: vasp_inputs/..., results/s09_calculation_matrix.csv
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

WF = Path(__file__).resolve().parents[1]
BLADE = WF.parents[1]
OUT = WF / "vasp_inputs"
POTCAR = {"Zr": "Zr_sv", "Nb": "Nb_pv", "Ta": "Ta_pv", "B": "B"}
COMMON = {"PREC": "Accurate", "ENCUT": 520, "EDIFF": "1E-6", "ISPIN": 1, "LREAL": "Auto", "LWAVE": ".FALSE.", "LCHARG": ".FALSE.",
          "KGAMMA": ".TRUE.", "NELM": 120, "ALGO": "Normal"}
PRESETS = {
    "relax_fixed_cell": {**COMMON, "IBRION": 2, "ISIF": 2, "NSW": 200, "EDIFFG": -0.01, "ISMEAR": 1, "SIGMA": 0.1, "KSPACING": 0.25},
    "relax_full": {**COMMON, "IBRION": 2, "ISIF": 3, "NSW": 200, "EDIFFG": -0.01, "ISMEAR": 1, "SIGMA": 0.1, "KSPACING": 0.25},
    "static_dos": {**COMMON, "IBRION": -1, "NSW": 0, "ISMEAR": -5, "NEDOS": 3001, "LORBIT": 11, "KSPACING": 0.15, "LCHARG": ".TRUE."},
}


def cost(n_atoms: int, kind: str) -> float:
    base = 400 * (n_atoms / 144) ** 2.5
    return max(1.0, round({"relax_fixed_cell": base, "relax_full": base * 1.3, "static_dos": base * 0.15}.get(kind, base), 0))


def geometry(struct_id: str) -> Path | None:
    for p in (WF / "results" / "orb" / struct_id / "CONTCAR", WF / "inputs" / "structures" / struct_id / "POSCAR"):
        if p.exists():
            return p
    return None


def blade_host(sid: str, x: float) -> Path:
    m = {"ZrNb": "Nb", "ZrTa": "Ta"}[sid]
    names = {0.0: "ZrB2_sqs_lev=0_a_Zr=1.0,b_B=1.0", 0.25: f"Zr3{m}B8_sqs_lev=2_a_{m}=0.25,a_Zr=0.75,b_B=1.0",
             0.5: f"Zr{m}B4_sqs_lev=1_a_{m}=0.5,a_Zr=0.5,b_B=1.0", 0.75: f"Zr{m}3B8_sqs_lev=2_a_{m}=0.75,a_Zr=0.25,b_B=1.0",
             1.0: f"{m}B2_sqs_lev=0_a_{m}=1.0,b_B=1.0"}
    return BLADE / "Files" / "system_structures" / ("NbZrB" if sid == "ZrNb" else "TaZrB") / "blade" / names[x] / "CONTCAR"


def matrix() -> list[dict]:
    """Rows: calc_id, priority, hypothesis, system, x, structure, source, n_atoms, calc_type, preset, cost, status, decides."""
    rows = []
    def add(prio, hyp, sid, x, sid_struct, src, n, ctype, preset, decides, status="prepared (not launched)"):
        rows.append({"priority": prio, "hypothesis": hyp, "system": sid, "x_solute": x, "structure": sid_struct, "geometry_source": str(src) if src else "",
                     "n_atoms": n, "calc_type": ctype, "incar_preset": preset, "est_core_hours": cost(n, preset) if preset in PRESETS else None,
                     "status": status, "decides": decides})
    for sid in ("ZrNb", "ZrTa"):
        m = sid[2:]
        # P1: H1 vacancy energetics (decisive quantity: dE_rec and local slope)
        for x, host in ((0.0, "ZrB2"), (1.0, f"{m}B2")):
            if x == 0 and sid == "ZrTa":
                continue  # ZrB2 shared with ZrNb
            vid = "V_ZrB2_x0.00_VB_n0a" if x == 0 else f"V_{sid}_x1.00_VB_n6a"
            add(1, "H1", sid if x else "ZrB2", x, f"{host} host 4x4x3", blade_host(sid, x), 144, "DFT relax (host)", "relax_full", "reference for E_f(V_B)")
            add(1, "H1", sid if x else "ZrB2", x, vid, geometry(vid), 143, "DFT relax V_B", "relax_fixed_cell", "dE_rec = E_v(MB2) - E_v(ZrB2)")
        refs = {"ZrNb": ["mp-160_B", "mp-131_Zr", "mp-75_Nb"], "ZrTa": ["mp-50_Ta"]}[sid]
        for ref in refs:
            p = WF / "inputs" / "competitors" / ref / "POSCAR"
            n = sum(int(v) for v in p.read_text().splitlines()[6].split()) if p.exists() else 0
            add(1, "H1/H4", sid, None, ref, p, n, "DFT relax element reference", "relax_full", "chemical-potential reference (bcc Ta, not beta-Ta)")
        for yid in ([f"Y_ZrB2_k4_spread"] if sid == "ZrNb" else []) + [f"Y_{m}B2_k4_spread", f"Y_{m}B2_k8_spread"]:
            add(1, "H1/H4", "ZrB2" if "ZrB2" in yid else sid, 0.0 if "ZrB2" in yid else 1.0, yid, geometry(yid), 144 - int(yid.split("_k")[1][0]),
                "DFT relax multi-V_B", "relax_fixed_cell", "vacancy-vacancy interaction W -> B-deficient homogeneity range of MB2")
        add(1, "H1", sid, 0.5, f"{sid} SQS host x=0.5", blade_host(sid, 0.5), 144, "DFT relax (host)", "relax_full", "reference for E_f(V_B) vs n")
        for n in (0, 3, 6):
            vid = f"V_{sid}_x0.50_VB_n{n}a"
            add(1, "H1", sid, 0.5, vid, geometry(vid), 143, "DFT relax V_B", "relax_fixed_cell", "slope dE_v/dn_solute")
        for comp in ({"ZrNb": ["mp-1080021_Nb2B3", "mp-10255_Nb3B4", "mp-1102394_Nb5B6", "mp-2580_NbB"],
                      "ZrTa": ["mp-10142_Ta3B4", "mp-28629_Ta5B6", "mp-1097_TaB"]}[sid]):
            p = WF / "inputs" / "competitors" / comp / "POSCAR"
            n = sum(int(v) for v in p.read_text().splitlines()[6].split()) if p.exists() else 0
            add(1, "H1/H4", sid, 1.0, comp, p, n, "DFT relax competitor", "relax_full", "B-poor mu_B limit; H4 phase energy")
        # P2: H4 framework re-run with DFT energies (no new structures beyond P1) -> analysis only
        add(2, "H4", sid, 0.5, "framework s08 with DFT parameters", None, 0, "analysis", "-", "Group-V-rich M3B4/M5B6 + Zr-rich AlB2 at B/M 1.95-1.98?",
            status="analysis step after P1")
        # P3: H2 configurations (2x2x8 superlattices, cheapest consistent cell)
        for s in (f"L_{sid}_end0", f"L_{sid}_end1", f"L_{sid}_p1", f"L_{sid}_p2", f"L_{sid}_p4"):
            add(3, "H2", sid, 0.5 if "_p" in s else (0.0 if s.endswith("0") else 1.0), s, geometry(s), 96, "DFT relax configuration", "relax_full",
                "coherent segregation energy vs layer thickness")
        for s in (f"C_{sid}_x0.50_qsqs1", f"C_{sid}_x0.50_alay2"):
            add(4, "H2", sid, 0.5, s, geometry(s), 192, "DFT relax configuration", "relax_full", "in-plane segregation vs SQS (only if ORB H2 inconclusive)")
        # P3: H3-el DOS
        for x in (0.0, 0.25, 0.5, 0.75, 1.0):
            if x == 0.0 and sid == "ZrTa":
                continue
            add(3, "H3-el", sid if x else "ZrB2", x, f"BLADE SQS x={x}", blade_host(sid, x), 144, "DFT static DOS", "static_dos", "dg(E_F) excess at x=0.5")
        # P4: H3-vib phonons
        add(4, "H3-vib", sid, 0.5, f"ordered ZrMB4 + end members (primitive)", None, 6, "phonons (ORB first, DFT check)", "-", "sign of dS_vib,ex",
            status="phonopy conf prepared; ORB phonons recommended first")
        add(5, "H3-vib", sid, 0.5, "SQS x=0.5 harmonic phonons", blade_host(sid, 0.5), 144, "phonons DFT (864 displ. w/o symmetry)", "-",
            "independent check of Tang QHA sign", status="not recommended unless ORB phonons disagree with Tang")
        # P5: liquid
        add(5, "H6-liquid", sid, 0.5, "AIMD liquid (Zr,M)B2 at ~3300 K", None, 200, "DFT-MD", "-", "liquid mixing W_L for eutectic temperature",
            status="only if H1/H4 open a solid gap")
    df = pd.DataFrame(rows)
    df.insert(0, "calc_id", [f"M{i:03d}" for i in range(1, len(df) + 1)])
    return df.to_dict("records")


def write_vasp(rows: list[dict]) -> None:
    for r in rows:
        if r["incar_preset"] not in PRESETS or not r["geometry_source"]:
            continue
        d = OUT / f"P{r['priority']}_{r['hypothesis'].replace('/', '-')}" / f"{r['calc_id']}_{r['structure'].replace(' ', '_').replace('=', '')}"
        d.mkdir(parents=True, exist_ok=True)
        from ase.io import read, write

        atoms = read(r["geometry_source"], format="vasp")
        write(d / "POSCAR", atoms, format="vasp", direct=True, sort=True)  # species grouped -> one POTCAR block each
        species = (d / "POSCAR").read_text().splitlines()[5].split()
        (d / "POTCAR.spec").write_text(" ".join(POTCAR[s] for s in species) + "\n")
        (d / "INCAR").write_text("".join(f"{k} = {v}\n" for k, v in PRESETS[r["incar_preset"]].items()))
        (d / "README").write_text(f"{r['calc_id']} | {r['hypothesis']} | {r['calc_type']} | decides: {r['decides']}\n"
                                  f"geometry: {r['geometry_source']}\nNOT LAUNCHED. Build POTCAR from POTCAR.spec (PBE PAW, 54 set).\n")
    # phonopy template for end members / ordered cells
    ph = OUT / "P4_H3-vib"
    ph.mkdir(parents=True, exist_ok=True)
    (ph / "phonopy_disp.conf").write_text("# phonopy -d phonopy_disp.conf -c POSCAR  (primitive AlB2 / ZrMB4 cells)\n"
                                         "DIM = 4 4 3\nDISPLACEMENT_DISTANCE = 0.01\nPM = AUTO\n")
    (ph / "README").write_text("Harmonic phonons for ZrB2, NbB2, TaB2 (3-atom) and ordered ZrNbB4/ZrTaB4 (6-atom, DIM 4 4 2).\n"
                               "Run with ORB first (cheap); DFT only for 2 cells if ORB contradicts Tang 2025's sign of the excess vibrational entropy.\n")


def main() -> None:
    rows = matrix()
    write_vasp(rows)
    df = pd.DataFrame(rows)
    df.to_csv(WF / "results" / "s09_calculation_matrix.csv", index=False)
    pd.set_option("display.width", 250)
    print(df[["calc_id", "priority", "hypothesis", "system", "x_solute", "structure", "n_atoms", "calc_type", "est_core_hours"]].to_string(index=False))
    print("total estimated core-hours by priority:")
    print(df.groupby("priority").est_core_hours.sum().to_string())


if __name__ == "__main__":
    main()
