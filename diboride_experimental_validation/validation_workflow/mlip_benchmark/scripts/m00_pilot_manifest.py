"""Stage 1 pilot manifest: the exact starting geometries used for the ORB screen, reused for every MLIP.

Large cells (BLADE 4x4x3 AlB2 supercells, 144 / 143 atoms):
  hosts      ZrB2, NbB2, TaB2 and the x = 0.5 SQS of ZrB2-NbB2 and ZrB2-TaB2 (BLADE ORB-relaxed CONTCARs, the same
             geometry the vacancy cells were cut from)
  vacancies  one V_B in each pure diboride; V_B at every n_solute = 0..6 (+ replicate n = 3) in both x = 0.5 SQS
Small reference cells (MP POSCARs, as for ORB): alpha-B (mp-160), Zr (mp-131), bcc Nb (mp-75), bcc Ta (mp-50),
  and the metal-rich borides that set the B-poor limit (Nb2B3, Nb3B4, Nb5B6, NbB, Ta3B4, Ta5B6, TaB).
Output: mlip_benchmark/pilot_manifest.csv (structure_id, role, system, x_solute, n_solute, n_atoms, start_geometry, orb_result)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

WF = Path(__file__).resolve().parents[2]
BLADE = WF.parents[1]
MB = WF / "mlip_benchmark"
HOST = {"ZrB2": ("NbZrB", "ZrB2_sqs_lev=0_a_Zr=1.0,b_B=1.0"), "NbB2": ("NbZrB", "NbB2_sqs_lev=0_a_Nb=1.0,b_B=1.0"),
        "TaB2": ("TaZrB", "TaB2_sqs_lev=0_a_Ta=1.0,b_B=1.0"), "ZrNb_x0.50": ("NbZrB", "ZrNbB4_sqs_lev=1_a_Nb=0.5,a_Zr=0.5,b_B=1.0"),
        "ZrTa_x0.50": ("TaZrB", "ZrTaB4_sqs_lev=1_a_Ta=0.5,a_Zr=0.5,b_B=1.0")}
REFS = ["mp-160_B", "mp-131_Zr", "mp-75_Nb", "mp-50_Ta", "mp-1080021_Nb2B3", "mp-10255_Nb3B4", "mp-1102394_Nb5B6", "mp-2580_NbB",
        "mp-10142_Ta3B4", "mp-28629_Ta5B6", "mp-1097_TaB"]


def main() -> None:
    rows = []
    for hid, (sysdir, name) in HOST.items():
        p = BLADE / "Files" / "system_structures" / sysdir / "blade" / name
        rows.append({"structure_id": f"host_{hid}", "role": "host", "system": hid, "x_solute": 0.5 if "x0.50" in hid else None,
                     "n_solute": None, "start_geometry": str(p / "CONTCAR"), "orb_energy_eV": float((p / "energy").read_text().split()[0])})
    man = pd.read_csv(WF / "inputs" / "structure_manifest.csv")
    orb = pd.read_csv(WF / "results" / "s03_orb_energies.csv").set_index("structure_id")
    pick = man[(man.family == "V") & (man.defect == "V_B") & (man.x_solute.isin([0.0, 0.5, 1.0]))]
    for r in pick.itertuples():
        host = "ZrB2" if r.system == "ZrB2" else ({"ZrNb": "NbB2", "ZrTa": "TaB2"}[r.system] if r.x_solute == 1.0 else f"{r.system}_x0.50")
        rows.append({"structure_id": r.structure_id, "role": "V_B", "system": host, "x_solute": r.x_solute, "n_solute": r.local_n_solute,
                     "start_geometry": str(WF / r.poscar), "orb_energy_eV": float(orb.loc[r.structure_id, "energy_eV"])})
    for ref in REFS:
        rows.append({"structure_id": f"ref_{ref}", "role": "reference", "system": ref.split("_")[1], "x_solute": None, "n_solute": None,
                     "start_geometry": str(WF / "inputs" / "competitors" / ref / "POSCAR"),
                     "orb_energy_eV": float(orb.loc[f"comp_{ref}", "energy_eV"])})
    df = pd.DataFrame(rows)
    from ase.io import read

    df["n_atoms"] = [len(read(p, format="vasp")) for p in df.start_geometry]
    df.to_csv(MB / "pilot_manifest.csv", index=False)
    print(df[["structure_id", "role", "system", "n_solute", "n_atoms"]].to_string(index=False))
    print(f"{len(df)} structures, {int((df.n_atoms > 100).sum())} large cells")


if __name__ == "__main__":
    main()
