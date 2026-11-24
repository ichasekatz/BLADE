"""Stage 2 (reduced) manifest: structures added to the Stage 1 pilot, same starting geometries as the ORB screen.

  Y   multi-vacancy cells in ZrB2, NbB2, TaB2 (k = 2, 4, 8 spread; k = 4 clustered; k = 8 confined to one B plane)
      -> each model's own vacancy-vacancy interaction W and vacancy ordering / clustering energies
  H   stoichiometric SQS hosts at x = 0.25 and 0.75 (BLADE ORB-relaxed CONTCARs)
      -> each model's own stoichiometric mixing energy (RK fit with the Stage 1 x = 0, 0.5, 1 hosts)
  V   B vacancies at x = 0.25 (n_solute = 0 and 4) and x = 0.75 (n = 2 and 6): lowest/highest environment present
      -> coupling across composition
  C   x = 0.5 metal configurations in the 4x4x4 cell: quasi-random, alternating c-layers (= ordered ZrMB4), AABB slab
      -> metal ordering vs random vs segregation (per-atom end-member references from the Stage 1 hosts: exact for
         perfect periodic crystals with an MLIP, which has no k-point sampling)
  S   one B vacancy in the AABB cell in a Zr|Zr, Zr|M or M|M boron layer -> coupled metal-vacancy arrangements
Output: mlip_benchmark/stage2_manifest.csv (same columns as pilot_manifest.csv)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

MB = Path(__file__).resolve().parents[1]
WF = MB.parent
BLADE = WF.parents[1]


def main() -> None:
    man = pd.read_csv(WF / "inputs" / "structure_manifest.csv").set_index("structure_id")
    orb = pd.read_csv(WF / "results" / "s03_orb_energies.csv").set_index("structure_id")
    rows = []

    def add(sid: str, role: str, system: str, x, n, geom: Path, e_orb: float) -> None:
        rows.append({"structure_id": sid, "role": role, "system": system, "x_solute": x, "n_solute": n, "start_geometry": str(geom),
                     "orb_energy_eV": e_orb})

    for em in ("ZrB2", "NbB2", "TaB2"):
        for tag in ("k2_spread", "k4_spread", "k8_spread", "k4_cluster", "k8_plane"):
            sid = f"Y_{em}_{tag}"
            add(sid, "multi_V_B", em, None, None, WF / man.loc[sid, "poscar"], orb.loc[sid, "energy_eV"])
    for sysn, sysdir, m in (("ZrNb", "NbZrB", "Nb"), ("ZrTa", "TaZrB", "Ta")):
        for x, name in ((0.25, f"Zr3{m}B8_sqs_lev=2_a_{m}=0.25,a_Zr=0.75,b_B=1.0"), (0.75, f"Zr{m}3B8_sqs_lev=2_a_{m}=0.75,a_Zr=0.25,b_B=1.0")):
            p = BLADE / "Files" / "system_structures" / sysdir / "blade" / name
            add(f"host_{sysn}_x{x:.2f}", "host", f"{sysn}_x{x:.2f}", x, None, p / "CONTCAR", float((p / "energy").read_text().split()[0]))
        for sid, x, n in ((f"V_{sysn}_x0.25_VB_n0a", 0.25, 0), (f"V_{sysn}_x0.25_VB_n4a", 0.25, 4),
                          (f"V_{sysn}_x0.75_VB_n2a", 0.75, 2), (f"V_{sysn}_x0.75_VB_n6a", 0.75, 6)):
            add(sid, "V_B", f"{sysn}_x{x:.2f}", x, n, WF / man.loc[sid, "poscar"], orb.loc[sid, "energy_eV"])
        for tag in ("qsqs1", "clay1", "clay2"):
            sid = f"C_{sysn}_x0.50_{tag}"
            add(sid, "config", sysn, 0.5, None, WF / man.loc[sid, "poscar"], orb.loc[sid, "energy_eV"])
        for tag, n in (("ZrZr", 0), ("ZrM", 3), ("MM", 6)):
            sid = f"S_{sysn}_VB_{tag}"
            add(sid, "coupled_V_B", sysn, 0.5, n, WF / man.loc[sid, "poscar"], orb.loc[sid, "energy_eV"])
    df = pd.DataFrame(rows)
    from ase.io import read

    df["n_atoms"] = [len(read(p, format="vasp")) for p in df.start_geometry]
    df.to_csv(MB / "stage2_manifest.csv", index=False)
    print(df.groupby("role").agg(n=("structure_id", "size"), atoms=("n_atoms", "sum")).to_string())
    print(f"{len(df)} structures per model")


if __name__ == "__main__":
    main()
