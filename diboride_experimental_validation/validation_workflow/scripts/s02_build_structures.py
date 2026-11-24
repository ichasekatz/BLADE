"""Step 2: build every structure in the calculation matrix (POSCARs only; nothing is relaxed here).

Families (all AlB2-derived, metal sublattice ordering only; B sublattice full unless a vacancy is stated):

V   Vacancies in the EXISTING relaxed BLADE 144-atom SQS (4x4x3 AlB2 cells, 48 M + 96 B) at x = 0, 0.25, 0.5,
    0.75, 1. V_B sites are chosen to span the number of solute (Nb/Ta) atoms among the 6 metal neighbours of the
    B site (n_sol = 0..6). V_Zr / V_M sites are chosen at the lowest and highest solute count among the 8
    nearest metal neighbours. The pristine host energies are reused from BLADE (not recomputed).
C   Configurations in a common 4x4x4 cell (64 M + 128 B = 192 atoms) at x = 0.25, 0.5, 0.75: two
    quasi-random (SQS-like, pair-correlation optimised) realisations, c-layered (period 1 and 2 at x = 0.5;
    single solute layer at 0.25/0.75), a-layered rows (period 1 and 2), and a compact column cluster
    (x = 0.25/0.75). End members are built in the same cell.
L   c-superlattice series in a 2x2x8 cell (32 M + 64 B = 96 atoms) at x = 0.5 with 1, 2 and 4 consecutive
    layers of each metal (period 2, 4, 8 layers), plus end members in the same cell.
S   Vacancy-segregation test: the 4x4x4 c-layered AABB cell at x = 0.5 with one B vacancy in the B layer
    between two Zr layers, at a Zr|M interface, or between two M layers.

Writes inputs/structures/<id>/POSCAR and inputs/structure_manifest.csv. Deterministic (fixed seeds).
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read, write

BLADE = Path(__file__).resolve().parents[3]
WF = Path(__file__).resolve().parents[1]
OUTS = WF / "inputs" / "structures"
SOLUTE = {"ZrNb": "Nb", "ZrTa": "Ta"}
SYSDIR = {"ZrNb": "NbZrB", "ZrTa": "TaZrB"}
SQS_NAME = {  # BLADE folder names (x = solute fraction)
    0.0: "ZrB2_sqs_lev=0_a_Zr=1.0,b_B=1.0", 0.25: "Zr3{M}B8_sqs_lev=2_a_{M}=0.25,a_Zr=0.75,b_B=1.0",
    0.5: "Zr{M}B4_sqs_lev=1_a_{M}=0.5,a_Zr=0.5,b_B=1.0", 0.75: "Zr{M}3B8_sqs_lev=2_a_{M}=0.75,a_Zr=0.25,b_B=1.0",
    1.0: "{M}B2_sqs_lev=0_a_{M}=1.0,b_B=1.0"}
# ORB-relaxed end-member lattice constants (A) from the BLADE SQS cells (a = |a1|/4, c = |a3|/3), used as
# starting geometry only (every structure is fully relaxed afterwards).
rows: list[dict] = []


def blade_sqs(sid: str, x: float) -> Atoms:
    m = SOLUTE[sid]
    return read(BLADE / "Files" / "system_structures" / SYSDIR[sid] / "blade" / SQS_NAME[x].format(M=m) / "CONTCAR")


def lattice_constants(sid: str, x: float) -> tuple[float, float]:
    at = blade_sqs(sid, x)
    return float(np.linalg.norm(at.cell[0]) / 4), float(np.linalg.norm(at.cell[2]) / 3)


def alb2_supercell(a: float, c: float, n: tuple[int, int, int]) -> tuple[Atoms, np.ndarray]:
    """AlB2 supercell with all metal sites 'Zr'. Returns atoms and integer (i, j, k) cell index of each metal."""
    prim = Atoms("ZrB2", scaled_positions=[(0, 0, 0), (1 / 3, 2 / 3, 0.5), (2 / 3, 1 / 3, 0.5)],
                 cell=[[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]], pbc=True)
    sc = prim.repeat(n)
    idx = np.array([np.round(p * n).astype(int) % n for p, s in zip(sc.get_scaled_positions(), sc.symbols) if s == "Zr"])
    return sc, idx


def metal_indices(at: Atoms) -> np.ndarray:
    return np.array([i for i, s in enumerate(at.symbols) if s != "B"])


def pair_shells(at: Atoms, metals: np.ndarray, rmax: float = 6.5) -> list[np.ndarray]:
    """Boolean adjacency matrices of the metal sublattice, one per distinct neighbour distance (rounded 0.05 A)."""
    d = at.get_all_distances(mic=True)[np.ix_(metals, metals)]
    shells = sorted({round(v / 0.05) * 0.05 for v in d.ravel() if 0.1 < v < rmax})
    return [np.isclose(np.round(d / 0.05) * 0.05, s) for s in shells]


def quasi_random(at: Atoms, x: float, seed: int, n_steps: int = 40000) -> Atoms:
    """SQS-like decoration: Metropolis swaps minimising the deviation of shell pair correlations from random."""
    rng = np.random.default_rng(seed)
    metals = metal_indices(at)
    shells = pair_shells(at, metals)
    n_sol = int(round(x * len(metals)))
    spin = np.full(len(metals), 1.0)
    spin[rng.choice(len(metals), n_sol, replace=False)] = -1.0
    target = (1 - 2 * x) ** 2
    weights = np.array([1.0 / (k + 1) for k in range(len(shells))])

    def cost(s: np.ndarray) -> float:
        return float(sum(w * abs((s @ sh @ s) / sh.sum() - target) for w, sh in zip(weights, shells)))

    cur, temp = cost(spin), 0.05
    for step in range(n_steps):
        i, j = rng.choice(np.where(spin > 0)[0]), rng.choice(np.where(spin < 0)[0])
        spin[i], spin[j] = spin[j], spin[i]
        new = cost(spin)
        if new <= cur or rng.random() < np.exp((cur - new) / temp):
            cur = new
        else:
            spin[i], spin[j] = spin[j], spin[i]
        temp = max(1e-4, 0.05 * (1 - step / n_steps))
    out = at.copy()
    return out, spin, cur


def decorate(at: Atoms, solute: str, mask: np.ndarray) -> Atoms:
    out = at.copy()
    metals = metal_indices(out)
    syms = out.get_chemical_symbols()
    for k, i in enumerate(metals):
        syms[i] = solute if mask[k] else "Zr"
    out.set_chemical_symbols(syms)
    return out


def save(at: Atoms, sid_: str, **meta) -> None:
    folder = OUTS / sid_
    folder.mkdir(parents=True, exist_ok=True)
    at = at[np.argsort([{"Zr": 0, "Nb": 1, "Ta": 1, "B": 2}[s] for s in at.symbols], kind="stable")]
    write(folder / "POSCAR", at, format="vasp", direct=True, sort=False)
    sym = at.get_chemical_symbols()
    n_m = sum(s != "B" for s in sym)
    rows.append({"structure_id": sid_, "n_atoms": len(at), "formula": at.get_chemical_formula(mode="metal"),
                 "B_per_metal": sym.count("B") / n_m, "poscar": str((folder / "POSCAR").relative_to(WF)), **meta})


def b_site_env(at: Atoms, solute: str) -> pd.DataFrame:
    metals = metal_indices(at)
    b_idx = np.array([i for i, s in enumerate(at.symbols) if s == "B"])
    d = at.get_all_distances(mic=True)[np.ix_(b_idx, metals)]
    nn = np.argsort(d, axis=1)[:, :6]
    n_sol = np.array([sum(at.symbols[metals[j]] == solute for j in row) for row in nn])
    return pd.DataFrame({"atom_index": b_idx, "n_solute_nn6": n_sol})


def m_site_env(at: Atoms, solute: str) -> pd.DataFrame:
    metals = metal_indices(at)
    d = at.get_all_distances(mic=True)[np.ix_(metals, metals)]
    nn = np.argsort(d, axis=1)[:, 1:9]  # 6 in-plane + 2 out-of-plane
    n_sol = np.array([sum(at.symbols[metals[j]] == solute for j in row) for row in nn])
    return pd.DataFrame({"atom_index": metals, "species": [at.symbols[i] for i in metals], "n_solute_nn8": n_sol})


def build_vacancies() -> None:
    done_zrb2 = False
    for sid, m in SOLUTE.items():
        for x in (0.0, 0.25, 0.5, 0.75, 1.0):
            if x == 0.0 and done_zrb2:
                continue
            host = blade_sqs(sid, x)
            sys_label = "ZrB2" if x == 0.0 else sid
            parent = str((BLADE / "Files" / "system_structures" / SYSDIR[sid] / "blade" / SQS_NAME[x].format(M=m)).relative_to(BLADE.parent))
            env_b = b_site_env(host, m)
            picks = []
            for n, g in env_b.groupby("n_solute_nn6"):
                picks.append((n, int(g.atom_index.iloc[0]), "a"))
                if x == 0.5 and n == 3 and len(g) > 1:  # replicate site: sampling noise at fixed n_sol
                    picks.append((n, int(g.atom_index.iloc[len(g) // 2]), "b"))
            for n, idx, tag in picks:
                v = host.copy()
                del v[idx]
                save(v, f"V_{sys_label}_x{x:.2f}_VB_n{n}{tag}", family="V", system=sys_label, x_solute=x,
                     cell="BLADE SQS 4x4x3 (144)", config="SQS (BLADE)", defect="V_B", local_n_solute=n,
                     local_env_def="solute among 6 metal nn of B site", parent=parent,
                     hypothesis="H1 vacancy coupling; H1b B-deficient decomposition")
            env_m = m_site_env(host, m)
            for spec in sorted(env_m.species.unique()):
                g = env_m[env_m.species == spec].sort_values(["n_solute_nn8", "atom_index"])
                for tag, row in (("min", g.iloc[0]), ("max", g.iloc[-1])):
                    if tag == "max" and row.n_solute_nn8 == g.iloc[0].n_solute_nn8:
                        continue
                    v = host.copy()
                    del v[int(row.atom_index)]
                    save(v, f"V_{sys_label}_x{x:.2f}_V{spec}_{tag}n{int(row.n_solute_nn8)}", family="V", system=sys_label,
                         x_solute=x, cell="BLADE SQS 4x4x3 (144)", config="SQS (BLADE)", defect=f"V_{spec}",
                         local_n_solute=int(row.n_solute_nn8), local_env_def="solute among 8 metal nn of metal site",
                         parent=parent, hypothesis="H1 vacancy coupling (metal sublattice)")
            if x == 0.0:
                done_zrb2 = True


def build_configurations() -> None:
    n = (4, 4, 4)
    for sid, m in SOLUTE.items():
        a0, c0 = lattice_constants(sid, 0.0)
        a1, c1 = lattice_constants(sid, 1.0)
        for x in (0.0, 0.25, 0.5, 0.75, 1.0):
            a, c = a0 + (a1 - a0) * x, c0 + (c1 - c0) * x
            base, ijk = alb2_supercell(a, c, n)
            meta = dict(family="C", system=sid, x_solute=x, cell="4x4x4 (192)", defect="none",
                        hypothesis="H2 SQS vs segregated/ordered arrangements")
            if x in (0.0, 1.0):
                save(decorate(base, m, np.full(len(ijk), x == 1.0)), f"C_{sid}_x{x:.2f}_end", config="end member", **meta)
                continue
            for seed in (1, 2):
                _, spin, cst = quasi_random(base, x, seed=seed * 100 + int(100 * x))
                save(decorate(base, m, spin < 0), f"C_{sid}_x{x:.2f}_qsqs{seed}", config=f"quasi-random #{seed}",
                     sqs_cost=round(cst, 4), **meta)
            i, k = ijk[:, 0], ijk[:, 2]
            j = ijk[:, 1]
            if x == 0.5:
                configs = {"clay1": k % 2 == 1, "clay2": k // 2 == 1, "alay1": i % 2 == 1, "alay2": i // 2 == 1}
            else:
                one = (k == 0) if x == 0.25 else (k != 0)
                row = (i == 0) if x == 0.25 else (i != 0)
                col = ((i < 2) & (j < 2)) if x == 0.25 else ~((i < 2) & (j < 2))
                configs = {"clay": one, "alay": row, "column": col}
            label = {"clay1": "c-layered period 1 (=MP ZrMB4 order)", "clay2": "c-layered AABB (slab, 2 layers)",
                     "alay1": "a-rows alternating", "alay2": "a-rows 2+2 (slab normal to a)",
                     "clay": "single c-layer ordered", "alay": "single a-row plane", "column": "2x2 column cluster"}
            for key, mask in configs.items():
                assert abs(mask.mean() - x) < 1e-9, (key, mask.mean(), x)
                save(decorate(base, m, mask), f"C_{sid}_x{x:.2f}_{key}", config=label[key], **meta)


def build_superlattices() -> None:
    n = (2, 2, 8)
    for sid, m in SOLUTE.items():
        a0, c0 = lattice_constants(sid, 0.0)
        a1, c1 = lattice_constants(sid, 1.0)
        meta = dict(family="L", system=sid, cell="2x2x8 (96)", defect="none",
                    hypothesis="H2 coherent segregation energy vs layer thickness")
        for x in (0.0, 1.0):
            base, ijk = alb2_supercell(a0 + (a1 - a0) * x, c0 + (c1 - c0) * x, n)
            save(decorate(base, m, np.full(len(ijk), x == 1.0)), f"L_{sid}_end{int(x)}", x_solute=x, config="end member", **meta)
        base, ijk = alb2_supercell((a0 + a1) / 2, (c0 + c1) / 2, n)
        for p in (1, 2, 4):
            save(decorate(base, m, (ijk[:, 2] // p) % 2 == 1), f"L_{sid}_p{p}", x_solute=0.5,
                 config=f"c-superlattice {p}+{p} layers", **meta)


def build_segregation_test() -> None:
    for sid, m in SOLUTE.items():
        a0, c0 = lattice_constants(sid, 0.0)
        a1, c1 = lattice_constants(sid, 1.0)
        base, ijk = alb2_supercell((a0 + a1) / 2, (c0 + c1) / 2, (4, 4, 4))
        at = decorate(base, m, ijk[:, 2] // 2 == 1)  # metal layers k=0,1 Zr; k=2,3 solute
        # B layer at z = (k + 0.5)/4 lies between metal layers k and k+1 (mod 4)
        b_idx = [i for i, s in enumerate(at.symbols) if s == "B"]
        zb = at.get_scaled_positions()[b_idx, 2]
        for tag, k in (("ZrZr", 0), ("ZrM", 1), ("MM", 2)):
            cand = [b for b, z in zip(b_idx, zb) if abs(z - (k + 0.5) / 4) < 0.02]
            v = at.copy()
            del v[cand[0]]
            save(v, f"S_{sid}_VB_{tag}", family="S", system=sid, x_solute=0.5, cell="4x4x4 (192)",
                 config="c-layered AABB", defect="V_B", local_n_solute={"ZrZr": 0, "ZrM": 3, "MM": 6}[tag],
                 local_env_def="solute among 6 metal nn of B site", parent="C_" + sid + "_x0.50_clay2",
                 hypothesis="H1 vacancy segregation to M-rich regions")


def build_vacancy_series() -> None:
    """Y family: k = 2, 4, 8 B vacancies in the BLADE end-member cells (spread: greedy max-min distance) plus
    k = 4 clustered (nearest B neighbours). Gives the B-sublattice vacancy-vacancy interaction L(B,Va) per metal,
    needed so the CEF does not extrapolate dilute energies linearly to large vacancy fractions."""
    for sid, x, label in (("ZrNb", 0.0, "ZrB2"), ("ZrNb", 1.0, "NbB2"), ("ZrTa", 1.0, "TaB2")):
        host = blade_sqs(sid, x)
        b_idx = np.array([i for i, s in enumerate(host.symbols) if s == "B"])
        d = host.get_all_distances(mic=True)[np.ix_(b_idx, b_idx)]
        for k in (2, 4, 8):
            chosen = [0]
            while len(chosen) < k:
                chosen.append(int(np.argmax(d[:, chosen].min(axis=1))))
            v = host.copy()
            del v[[int(b_idx[c]) for c in chosen]]
            save(v, f"Y_{label}_k{k}_spread", family="Y", system=label, x_solute=x, cell="BLADE SQS 4x4x3 (144)", config="end member",
                 defect=f"{k} V_B spread", parent=label, hypothesis="H1/H4 vacancy-vacancy interaction (CEF L(B,Va))")
        # k = 8 vacancies ordered in ONE B plane (spread within that plane): vacancy-plane ordering, a precursor of the
        # M3B4 / M5B6 intergrowth motifs; compared with k8_spread (3-D random-like) at identical composition.
        zb = host.get_scaled_positions()[b_idx, 2]
        layer = np.where(np.abs(zb - zb[0]) < 0.02)[0]
        dl = d[np.ix_(layer, layer)]
        pl = [0]
        while len(pl) < 8:
            pl.append(int(np.argmax(dl[:, pl].min(axis=1))))
        v = host.copy()
        del v[[int(b_idx[layer[c]]) for c in pl]]
        save(v, f"Y_{label}_k8_plane", family="Y", system=label, x_solute=x, cell="BLADE SQS 4x4x3 (144)", config="end member",
             defect="8 V_B plane", parent=label, hypothesis="vacancy ordering (B-plane) vs dispersed vacancies")
        near = list(np.argsort(d[0])[:4])
        v = host.copy()
        del v[[int(b_idx[c]) for c in near]]
        save(v, f"Y_{label}_k4_cluster", family="Y", system=label, x_solute=x, cell="BLADE SQS 4x4x3 (144)", config="end member",
             defect="4 V_B clustered", parent=label, hypothesis="H1/H4 vacancy clustering / ordering")


def main() -> None:
    OUTS.mkdir(parents=True, exist_ok=True)
    build_vacancy_series()
    build_vacancies()
    build_configurations()
    build_superlattices()
    build_segregation_test()
    df = pd.DataFrame(rows)
    df.to_csv(WF / "inputs" / "structure_manifest.csv", index=False)
    print(df.groupby(["family", "system"]).size().to_string())
    print(f"total {len(df)} structures; atoms: {df.n_atoms.sum()}")


if __name__ == "__main__":
    main()
