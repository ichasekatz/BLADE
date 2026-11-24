"""Stage 4 structures: complete the six Ness systems + ordered B-vacancy search. POSCARs and a manifest only.

Same generators and conventions as Stage 1-3 (validation_workflow/scripts/s02_build_structures.py): starting
geometries are BLADE's ORB-relaxed 144-atom SQS (4x4x3 AlB2) cells; vacancies are cut from them.

A  TiB2-NbB2, TiB2-TaB2, HfB2-NbB2, HfB2-TaB2 (systems without vacancy data):
   hosts x = 0.25, 0.5, 0.75; V_B at x = 0.5 for every local environment present (n_solute = 0..6, replicate at 3);
   V_B at x = 0.25 / 0.75 for the lowest and highest environment; 4x4x4 metal configurations at x = 0.5
   (quasi-random, alternating c-layers, AABB slab); V_B in the host|host, host|M, M|M layers of the slab.
B  TiB2 and HfB2 end members: host, single V_B, k = 2/4/8 spread, k = 4 clustered, k = 8 single-plane vacancies.
C  Competitors: Ti (mp-72), Hf (mp-103) [BLADE's references], TiB, Ti3B4 (stable MP borides, relaxed earlier in
   bfix_rerun/MP_binary_borides). Hf has no stable boride other than HfB2 in MP.
D  Ordered B-vacancy enumeration in TaB2 and NbB2 (DFT: TaB2-x on the hull for x = 0.167-0.25, Ektarawong 2023;
   ordered NbB2-y stable up to ~10 % B sites, Dahlqvist 2022). All Hermite-normal-form supercells of the AlB2
   primitive cell with N = 2..6 formula units, k = 1 vacancy (and N = 6, 8 with k = 2), giving B/M = 1.5-1.917;
   duplicates removed with pymatgen StructureMatcher. Lattice from the ORB-relaxed BLADE end member.
Output: final_results/supplementary_data/stage4_structures/<id>/POSCAR, stage4_manifest.csv
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read, write

FR = Path(__file__).resolve().parents[1]
BLADE = FR.parents[1]
WF = FR.parent / "validation_workflow"
sys.path.insert(0, str(WF / "scripts"))
from s02_build_structures import alb2_supercell, b_site_env, decorate, m_site_env, quasi_random  # noqa: E402

OUT = FR / "supplementary_data" / "stage4_structures"
SYSTEMS = {"TiNb": ("NbTiB", "Ti", "Nb"), "TiTa": ("TaTiB", "Ti", "Ta"), "HfNb": ("HfNbB", "Hf", "Nb"), "HfTa": ("HfTaB", "Hf", "Ta")}
rows: list[dict] = []


def sqs_dir(sysdir: str, host: str, m: str, x: float) -> Path:
    base = BLADE / "Files" / "system_structures" / sysdir / "blade"
    if x == 0.0:
        return next(base.glob(f"{host}B2_sqs_lev=0_*"))
    if x == 1.0:
        return next(base.glob(f"{m}B2_sqs_lev=0_*"))
    lev = 1 if x == 0.5 else 2
    return next(d for d in base.glob(f"*_sqs_lev={lev}_*") if f"a_{m}={x}" in d.name)


def save(at: Atoms, sid: str, **meta) -> None:
    d = OUT / sid
    d.mkdir(parents=True, exist_ok=True)
    order = {"B": 2}
    at = at[np.argsort([order.get(s, 0 if s in ("Ti", "Hf", "Zr") else 1) for s in at.symbols], kind="stable")]
    write(d / "POSCAR", at, format="vasp", direct=True, sort=False)
    rows.append({"structure_id": sid, "n_atoms": len(at), "start_geometry": str(d / "POSCAR"), **meta})


def build_systems() -> None:
    for sid, (sysdir, host, m) in SYSTEMS.items():
        for x in (0.25, 0.5, 0.75):
            p = sqs_dir(sysdir, host, m, x)
            rows.append({"structure_id": f"host_{sid}_x{x:.2f}", "role": "host", "system": f"{sid}_x{x:.2f}", "x_solute": x,
                         "start_geometry": str(p / "CONTCAR"), "orb_energy_eV": float((p / "energy").read_text().split()[0]), "n_atoms": 144})
        for x in (0.25, 0.5, 0.75):
            hostat = read(sqs_dir(sysdir, host, m, x) / "CONTCAR", format="vasp")
            env = b_site_env(hostat, m)
            picks = []
            if x == 0.5:
                for n, g in env.groupby("n_solute_nn6"):
                    picks.append((n, int(g.atom_index.iloc[0]), "a"))
                    if n == 3 and len(g) > 1:
                        picks.append((n, int(g.atom_index.iloc[len(g) // 2]), "b"))
            else:
                for n in (env.n_solute_nn6.min(), env.n_solute_nn6.max()):
                    picks.append((int(n), int(env[env.n_solute_nn6 == n].atom_index.iloc[0]), "a"))
            for n, idx, tag in picks:
                v = hostat.copy()
                del v[idx]
                save(v, f"V_{sid}_x{x:.2f}_VB_n{n}{tag}", role="V_B", system=f"{sid}_x{x:.2f}", x_solute=x, n_solute=n)
        # metal configurations in 4x4x4
        h0 = read(sqs_dir(sysdir, host, m, 0.0) / "CONTCAR", format="vasp")
        h1 = read(sqs_dir(sysdir, host, m, 1.0) / "CONTCAR", format="vasp")
        a = (np.linalg.norm(h0.cell[0]) + np.linalg.norm(h1.cell[0])) / 8
        c = (np.linalg.norm(h0.cell[2]) + np.linalg.norm(h1.cell[2])) / 6
        base, ijk = alb2_supercell(a, c, (4, 4, 4))
        _, spin, _ = quasi_random(base, 0.5, seed=150)
        cfgs = {"qsqs1": spin < 0, "clay1": ijk[:, 2] % 2 == 1, "clay2": ijk[:, 2] // 2 == 1}
        for tag, mask in cfgs.items():
            at = decorate(base, m, mask)
            at.set_chemical_symbols([host if s == "Zr" else s for s in at.get_chemical_symbols()])
            save(at, f"C_{sid}_x0.50_{tag}", role="config", system=sid, x_solute=0.5)
            if tag == "clay2":
                b_idx = [i for i, s in enumerate(at.symbols) if s == "B"]
                zb = at.get_scaled_positions()[b_idx, 2]
                for lab, k, n in (("HH", 0, 0), ("HM", 1, 3), ("MM", 2, 6)):
                    cand = [b for b, z in zip(b_idx, zb) if abs(z - (k + 0.5) / 4) < 0.02]
                    v = at.copy()
                    del v[cand[0]]
                    save(v, f"S_{sid}_VB_{lab}", role="coupled_V_B", system=sid, x_solute=0.5, n_solute=n)


def build_end_members() -> None:
    for em, sysdir in (("TiB2", "NbTiB"), ("HfB2", "HfNbB")):
        p = next((BLADE / "Files" / "system_structures" / sysdir / "blade").glob(f"{em}_sqs_lev=0_*"))
        rows.append({"structure_id": f"host_{em}", "role": "host", "system": em, "x_solute": 0.0, "start_geometry": str(p / "CONTCAR"),
                     "orb_energy_eV": float((p / "energy").read_text().split()[0]), "n_atoms": 144})
        hostat = read(p / "CONTCAR", format="vasp")
        b_idx = np.array([i for i, s in enumerate(hostat.symbols) if s == "B"])
        d = hostat.get_all_distances(mic=True)[np.ix_(b_idx, b_idx)]
        v = hostat.copy()
        del v[int(b_idx[0])]
        save(v, f"V_{em}_VB", role="V_B", system=em, x_solute=0.0, n_solute=0)
        for k in (2, 4, 8):
            ch = [0]
            while len(ch) < k:
                ch.append(int(np.argmax(d[:, ch].min(axis=1))))
            v = hostat.copy()
            del v[[int(b_idx[c]) for c in ch]]
            save(v, f"Y_{em}_k{k}_spread", role="multi_V_B", system=em, k_vacancies=k)
        near = list(np.argsort(d[0])[:4])
        v = hostat.copy()
        del v[[int(b_idx[c]) for c in near]]
        save(v, f"Y_{em}_k4_cluster", role="multi_V_B", system=em, k_vacancies=4)
        zb = hostat.get_scaled_positions()[b_idx, 2]
        layer = np.where(np.abs(zb - zb[0]) < 0.02)[0]
        dl = d[np.ix_(layer, layer)]
        pl = [0]
        while len(pl) < 8:
            pl.append(int(np.argmax(dl[:, pl].min(axis=1))))
        v = hostat.copy()
        del v[[int(b_idx[layer[c]]) for c in pl]]
        save(v, f"Y_{em}_k8_plane", role="multi_V_B", system=em, k_vacancies=8)


def build_competitors() -> None:
    comp = [("Ti", BLADE / "Files" / "MaterialsProject_Comps_ORB" / "Ti" / "mp-72_Ti"),
            ("Hf", BLADE / "Files" / "MaterialsProject_Comps_ORB" / "Hf" / "mp-103_Hf")]
    for f in sorted((BLADE / "bfix_rerun" / "MP_binary_borides" / "B-Ti").iterdir()):
        comp.append((f.name.split("_")[-1], f))
    for name, d in comp:
        src = d / "POSCAR" if (d / "POSCAR").exists() else d / "CONTCAR"
        sid = f"ref_{d.name}"
        rows.append({"structure_id": sid, "role": "reference", "system": name, "start_geometry": str(src),
                     "orb_energy_eV": float((d / "energy").read_text().split()[0]) if (d / "energy").exists() else np.nan,
                     "n_atoms": len(read(src, format="vasp"))})


def hnf(n: int):
    for a in range(1, n + 1):
        if n % a:
            continue
        for c in range(1, n // a + 1):
            if (n // a) % c:
                continue
            f = n // a // c
            for b in range(c):
                for d in range(f):
                    for e in range(f):
                        yield np.array([[a, 0, 0], [b, c, 0], [d, e, f]])


def build_ordered_vacancies() -> None:
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.io.ase import AseAtomsAdaptor

    sm = StructureMatcher(ltol=0.1, stol=0.1, angle_tol=3, primitive_cell=False, scale=False)
    for em, sysdir in (("TaB2", "TaTiB"), ("NbB2", "NbTiB")):
        h = read(next((BLADE / "Files" / "system_structures" / sysdir / "blade").glob(f"{em}_sqs_lev=0_*")) / "CONTCAR", format="vasp")
        a, c = float(np.mean(h.cell.lengths()[:2]) / 4), float(h.cell.lengths()[2] / 3)
        prim = Atoms(em[:2] + "B2", scaled_positions=[(0, 0, 0), (1 / 3, 2 / 3, 0.5), (2 / 3, 1 / 3, 0.5)],
                     cell=[[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]], pbc=True)
        from ase.build import make_supercell

        for n, k in ((2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (6, 2), (8, 2)):
            uniq = []
            for H in hnf(n):
                sc = make_supercell(prim, H)
                if min(sc.cell.lengths()) < 2.5:
                    continue
                b_idx = [i for i, s in enumerate(sc.symbols) if s == "B"]
                for combo in itertools.combinations(b_idx, k):
                    if k == 2 and combo[0] != b_idx[0]:
                        continue  # translation: fix the first vacancy
                    v = sc.copy()
                    del v[list(combo)]
                    st = AseAtomsAdaptor.get_structure(v)
                    if any(sm.fit(st, u) for u in uniq):
                        continue
                    uniq.append(st)
                    if k == 1:
                        break  # one vacancy: every site of this supercell is equivalent by translation/symmetry of prim
            for i, st in enumerate(uniq):
                at = AseAtomsAdaptor.get_atoms(st)
                y = k / (2 * n)
                save(at, f"O_{em}_N{n}k{k}_{i:02d}", role="ordered_V_B", system=em, k_vacancies=k, n_fu=n, y_Va=y,
                     B_per_M=2 * (1 - y))
            print(em, n, k, len(uniq), flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_end_members()
    build_systems()
    build_competitors()
    build_ordered_vacancies()
    df = pd.DataFrame(rows)
    df.to_csv(FR / "supplementary_data" / "stage4_manifest.csv", index=False)
    print(df.groupby("role").agg(n=("structure_id", "size"), atoms=("n_atoms", "sum")).to_string())


if __name__ == "__main__":
    main()
