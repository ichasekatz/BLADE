"""Ordered B-vacancy ground states and a finite-temperature check by cluster expansion + Monte Carlo.

Data: Stage 4 enumerated ordered MB2-x cells (TaB2, NbB2; N = 2-8 formula units, k = 1-2 vacancies, 223 symmetry-
distinct structures each), the 144-atom single / multi-vacancy cells, and the pristine end member, all relaxed by each
model with the common protocol.
1. 0 K formation energy of MB2-y relative to MB2 + alpha-B (each model's own references):
      dE(y) = [E_cell - n_fu e_MB2 + k E_B] / n_fu      (eV per formula unit; y = k / (2 n_fu))
   Lower convex hull of dE(y) between y = 0 and the vacancy-richest cell (B/M >= 1.5). Points below 0 at y > 0 mean
   a B-deficient ordered phase is more stable than MB2 + B at 0 K (as the DFT cluster expansion of Ektarawong 2023
   reports for TaB2-x, x = 0.167-0.25). Published DFT benchmark is qualitative here (their structures not available).
2. Lattice-gas cluster expansion on the B sublattice (vacancy = 1): E_cell = n_fu J0 + k J1 + sum_s J_s N_s, N_s =
   number of vacancy-vacancy pairs in shell s (periodic images included), shells from the ideal (unrelaxed) lattice,
   up to 6.0 A. Ordinary least squares; leave-one-out CV error reported. Only accepted if LOO-CV < 10 meV per
   vacancy-containing cell per vacancy.
3. Canonical Metropolis Monte Carlo (vacancy-B swaps) at fixed y = 1/12 and 1/8 on an 8x8x6 primitive supercell
   (768 B sites), T = 300-3000 K: energy and heat capacity C = var(E)/(k_B T^2). A heat-capacity maximum indicates an
   order-disorder crossover OF THE FITTED MODEL; it is reported as such, not as a measured or DFT transition
   temperature (relaxation, CE truncation and the absence of vibrational entropy limit it).
Outputs: supplementary_data/z01_ordered_hull.csv, z01_ce_fits.csv, z01_mc.csv
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

FR = Path(__file__).resolve().parents[1]
SD = FR / "supplementary_data"
MB = FR.parent / "validation_workflow" / "mlip_benchmark"
MODELS = ["orb", "mace-mpa", "mace-omat", "mattersim", "sevennet"]
KB = 8.617333262e-5
RCUT = 6.0


def energy(model: str, sid: str) -> float | None:
    p = MB / "results" / model / sid / "result.json"
    if p.exists():
        j = json.loads(p.read_text())
        return float(j["energy_eV"]) if j.get("status") == "ok" else None
    return None


def host_energy(model: str, em: str) -> float | None:
    if model == "orb":
        man = pd.read_csv(MB / "pilot_manifest.csv").set_index("structure_id")
        return float(man.loc[f"host_{em}", "orb_energy_eV"])
    return energy(model, f"host_{em}")


def vacancy_sites(poscar: Path, em: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (cell, cartesian vacancy positions) from the unrelaxed POSCAR: ideal B sites minus occupied ones."""
    from ase.io import read

    at = read(poscar, format="vasp")
    cell = at.cell.array
    metals = at.positions[np.array([s != "B" for s in at.symbols])]
    bpos = at.positions[np.array([s == "B" for s in at.symbols])]
    # ideal B sites: each metal at R has B at R + (1/3, 2/3, 1/2) and R + (2/3, 1/3, 1/2) of the primitive cell
    # orientation-independent primitive vectors from the metal sublattice (cells may be Niggli-reduced / rotated)
    d = metals[None, :, :] - metals[:, None, :]
    frac = np.linalg.solve(cell.T, d.reshape(-1, 3).T).T
    frac -= np.round(frac)
    dc = (frac @ cell).reshape(len(metals), len(metals), 3)[0]
    ln = np.linalg.norm(dc, axis=1)
    vecs = dc[(ln > 0.5) & (ln < 4.0)]
    ext = []  # add short lattice vectors of the cell itself (single-layer / single-column cells have no metal image)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                v = np.array([i, j, k]) @ cell
                if 0.5 < np.linalg.norm(v) < 4.0:
                    ext.append(v)
    vecs = np.vstack([vecs] + ([np.array(ext)] if ext else []))
    L = np.linalg.norm(vecs, axis=1)

    def partners(i):
        cosv = vecs @ vecs[i] / (L * L[i])
        return np.where((np.abs(L - L[i]) < 0.05) & (np.abs(cosv - 0.5) < 0.03))[0]

    inpl = [i for i in range(len(vecs)) if len(partners(i)) >= 1]
    i1 = min(inpl, key=lambda i: L[i])
    a1, a2 = vecs[i1], vecs[partners(i1)[0]]
    nrm = np.cross(a1, a2)
    nrm /= np.linalg.norm(nrm)
    proj = np.abs(vecs @ nrm)
    cz = float(np.min(proj[proj > 0.5]))
    a3 = cz * nrm
    t1, t2 = (a1 + a2) / 3 + a3 / 2, 2 * (a1 + a2) / 3 + a3 / 2
    ideal = np.vstack([metals + t1, metals + t2])
    vac = []
    for p in ideal:
        f = np.linalg.solve(cell.T, (bpos - p).T).T
        f -= np.round(f)
        if np.min(np.linalg.norm(f @ cell, axis=1)) > 0.6:
            fp = np.linalg.solve(cell.T, p)
            fp -= np.floor(fp)
            if not any(np.linalg.norm(((fp - q) - np.round(fp - q)) @ cell) < 0.3 for q in [np.linalg.solve(cell.T, v) for v in vac]):
                vac.append(fp @ cell)
    return cell, np.array(vac)


def pair_shells(cell: np.ndarray, vac: np.ndarray, shells: np.ndarray) -> np.ndarray:
    counts = np.zeros(len(shells))
    if len(vac) == 0:
        return counts
    rng = range(-3, 4)
    imgs = np.array([[i, j, k] for i in rng for j in rng for k in rng]) @ cell
    for i in range(len(vac)):
        for j in range(len(vac)):
            d = np.linalg.norm(vac[j] - vac[i] + imgs, axis=1)
            d = d[(d > 0.1) & (d < RCUT)]
            for dd in d:
                s = np.argmin(np.abs(shells - dd))
                if abs(shells[s] - dd) < 0.12:
                    counts[s] += 0.5  # each pair counted twice over (i, j)
    return counts


def main() -> None:
    man = pd.read_csv(SD / "stage4_manifest.csv")
    ords = man[man.role == "ordered_V_B"]
    hull_rows, fit_rows, mc_rows = [], [], []
    for em in ("TaB2", "NbB2"):
        geo = {}
        for r in ords[ords.system == em].itertuples():
            geo[r.structure_id] = vacancy_sites(Path(r.start_geometry), em)
        alld = np.concatenate([np.linalg.norm(v[1][:, None] - v[1][None], axis=2).ravel() for v in geo.values() if len(v[1]) > 1] + [np.array([])])
        # shells: distinct B-B distances of the ideal lattice up to RCUT, from the largest cell
        big = max(geo, key=lambda s: len(geo[s][1]) + 0 * len(s))
        cellb, _ = geo[big]
        shells = np.array(sorted({round(x, 1) for x in np.linspace(1.7, RCUT, 1)}))
        # build shells from an ideal reference B lattice of a 4x4x3 MB2 cell
        from ase.io import read

        ref = read(next(Path("/home/ichasekatz/BLADE/BLADE/Files/system_structures").glob(f"*/blade/{em}_sqs_lev=0_*")) / "CONTCAR", format="vasp")
        bp = ref.positions[np.array([s == "B" for s in ref.symbols])]
        dd = ref.get_all_distances(mic=True)[np.ix_(np.array([s == "B" for s in ref.symbols]), np.array([s == "B" for s in ref.symbols]))]
        shells = np.array(sorted({round(float(x), 1) for x in dd.ravel() if 0.5 < x < RCUT}))
        for model in MODELS:
            e0 = host_energy(model, em)
            eb = energy(model, "ref_mp-160_B") if model != "orb" else float(pd.read_csv(MB / "pilot_manifest.csv").set_index("structure_id").loc["ref_mp-160_B", "orb_energy_eV"])
            if e0 is None or eb is None:
                continue
            efu, eB = e0 / 48, eb / 12
            X, Y, recs = [], [], []
            for r in ords[ords.system == em].itertuples():
                e = energy(model, r.structure_id)
                if e is None:
                    continue
                n, k = int(r.n_fu), int(r.k_vacancies)
                de = (e - n * efu + k * eB) / n
                recs.append({"model": model, "end_member": em, "structure_id": r.structure_id, "n_fu": n, "k": k, "y": k / (2 * n),
                             "B_per_M": 2 - k / n, "dE_eV_fu": de, "dE_per_vacancy_eV": de * n / k})
                cell, vac = geo[r.structure_id]
                X.append(np.r_[n, k, pair_shells(cell, vac, shells)])
                Y.append(e - n * efu)  # energy relative to the stoichiometric reference; J0 absorbs residual
            if not recs:
                continue
            h = pd.DataFrame(recs)
            # lower hull in (y, dE)
            pts = sorted([(0.0, 0.0, "MB2")] + [(r.y, r.dE_eV_fu, r.structure_id) for r in h.itertuples()])
            hull: list = []
            for p in pts:
                while len(hull) >= 2 and (hull[-1][0] - hull[-2][0]) * (p[1] - hull[-2][1]) - (hull[-1][1] - hull[-2][1]) * (p[0] - hull[-2][0]) <= 0:
                    hull.pop()
                hull.append(p)
            on = {p[2] for p in hull}
            h["on_0K_hull"] = h.structure_id.isin(on) & (h.dE_eV_fu < 0)  # stable only if below the MB2 + B tie line
            best = h.loc[h.groupby("y").dE_eV_fu.idxmin()]
            hull_rows.append(h)
            # CE fit (n_fu coefficient absorbs per-cell reference; k = point term; pair shells)
            X, Y = np.array(X, float), np.array(Y, float)
            keep = X[:, 2:].sum(axis=0) > 0
            Xf = np.c_[X[:, :2], X[:, 2:][:, keep]]
            coef, *_ = np.linalg.lstsq(Xf, Y, rcond=None)
            loo = []
            for i in range(len(Y)):
                msk = np.arange(len(Y)) != i
                c_, *_ = np.linalg.lstsq(Xf[msk], Y[msk], rcond=None)
                loo.append((Xf[i] @ c_ - Y[i]) / max(X[i, 1], 1))
            cv = float(np.sqrt(np.mean(np.square(loo))))
            fit_rows.append({"model": model, "end_member": em, "n_structures": len(Y), "n_pair_shells": int(keep.sum()),
                             "J_point_eV": coef[1], "loo_cv_eV_per_vacancy": cv, "accepted": cv < 0.010,
                             "shells_A": ";".join(f"{s:.1f}" for s in shells[keep]), "J_pairs_eV": ";".join(f"{v:.4f}" for v in coef[2:]),
                             "min_dE_eV_fu": float(h.dE_eV_fu.min()), "y_of_min": float(h.loc[h.dE_eV_fu.idxmin(), "y"]),
                             "ordered_phase_below_MB2+B": bool((h.dE_eV_fu < 0).any()), "n_cells_below_MB2+B": int((h.dE_eV_fu < 0).sum()),
                             "dE_per_vac_spread_eV": float(h.dE_per_vacancy_eV.max() - h.dE_per_vacancy_eV.min())})
            print(em, model, f"{len(Y)} cells, CV {1000*cv:.1f} meV/vac, min dE {h.dE_eV_fu.min():+.3f} eV/fu at y={fit_rows[-1]['y_of_min']:.3f}", flush=True)
            if cv < 0.010:
                mc_rows += monte_carlo(model, em, coef[2:], shells[keep], ref)
    if hull_rows:
        pd.concat(hull_rows).to_csv(SD / "z01_ordered_hull.csv", index=False)
    pd.DataFrame(fit_rows).to_csv(SD / "z01_ce_fits.csv", index=False)
    pd.DataFrame(mc_rows).to_csv(SD / "z01_mc.csv", index=False)


def monte_carlo(model, em, J, shells, ref, seed=7) -> list[dict]:
    """Canonical lattice-gas MC on the B sublattice of an 8x8x6 primitive supercell built from the relaxed end member."""
    from ase import Atoms

    a = float(np.mean(ref.cell.lengths()[:2]) / 4)
    c = float(ref.cell.lengths()[2] / 3)
    prim = Atoms("XB2", scaled_positions=[(0, 0, 0), (1 / 3, 2 / 3, 0.5), (2 / 3, 1 / 3, 0.5)],
                 cell=[[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]], pbc=True).repeat((8, 8, 6))
    b = np.array([s == "B" for s in prim.get_chemical_symbols()])
    pos = prim.positions[b]
    d = prim.get_all_distances(mic=True)[np.ix_(b, b)]
    nb = len(pos)
    nbr = []  # per site: list of (neighbour, J)
    for i in range(nb):
        lst = []
        for s, j in zip(shells, J):
            idx = np.where(np.abs(d[i] - s) < 0.12)[0]
            lst += [(int(q), float(j)) for q in idx if q != i]
        nbr.append(lst)
    rng = np.random.default_rng(seed)
    out = []
    for y in (1 / 12, 1 / 8):
        nv = int(round(y * nb))
        occ = np.zeros(nb, int)
        occ[rng.choice(nb, nv, replace=False)] = 1

        def local(i):
            return sum(jj * occ[q] for q, jj in nbr[i])

        for t in (3000.0, 2500.0, 2000.0, 1500.0, 1000.0, 600.0, 300.0):  # annealing schedule, continue from previous
            beta = 1 / (KB * t)
            es = []
            e = 0.5 * sum(occ[i] * local(i) for i in range(nb))
            for step in range(60 * nb):
                v = rng.choice(np.where(occ == 1)[0])
                bsite = rng.choice(np.where(occ == 0)[0])
                de = -local(v) + local(bsite) - sum(jj for q, jj in nbr[bsite] if q == v)
                if de <= 0 or rng.random() < np.exp(-beta * de):
                    occ[v], occ[bsite] = 0, 1
                    e += de
                if step > 20 * nb and step % nb == 0:
                    es.append(e)
            es = np.array(es)
            out.append({"model": model, "end_member": em, "y": y, "T_K": t, "E_per_vacancy_eV": float(es.mean() / nv),
                        "C_per_vacancy_kB": float(es.var() / (KB * t) ** 2 / nv)})
    return out


if __name__ == "__main__":
    main()
