"""Step 0: audit existing computational data for reusable Zr-Nb-B / Zr-Ta-B structures.

Scans every BLADE structure tree (read-only) for relaxed cells whose species are a subset of {Zr, Nb, B} or
{Zr, Ta, B}. Each cell is classified by metal composition, B/metal ratio (vacancy / off-stoichiometry), origin
(SQS, MP-ordered, competitor), and a geometry fingerprint. The fingerprint separates distinct SQS realizations
from copies made by ``reuse_suborder_relaxations``. Nothing is written outside validation_workflow/results/.

Output: results/s00_existing_structure_audit.csv, results/s00_audit_summary.md
"""

from __future__ import annotations

import hashlib
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

BLADE = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "results"
ROOTS = ["Files/system_structures", "Files/system_structures1", "Files/system_structures_4metal",
         "Files/system_structures_5metal", "Files/Comps", "Files/Comps1", "Files/MaterialsProject_Comps_ORB",
         "Files/SQS", "../SQS"]
SYSTEMS = {"ZrNb": {"Zr", "Nb", "B"}, "ZrTa": {"Zr", "Ta", "B"}}


def read_poscar(path: Path) -> tuple[np.ndarray, list[str], list[int], np.ndarray] | None:
    """Minimal VASP5 POSCAR/CONTCAR reader (direct or Cartesian). Returns lattice, species, counts, frac coords."""
    try:
        lines = path.read_text().splitlines()
        scale = float(lines[1].split()[0])
        lat = np.array([[float(v) for v in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
        species = lines[5].split()
        counts = [int(v) for v in lines[6].split()]
        mode_line = 7
        if lines[mode_line].strip()[0] in "sS":
            mode_line += 1
        cart = lines[mode_line].strip()[0] in "cCkK"
        n = sum(counts)
        xyz = np.array([[float(v) for v in lines[mode_line + 1 + i].split()[:3]] for i in range(n)])
        frac = np.linalg.solve(lat.T, (xyz * scale).T).T if cart else xyz
        return lat, species, counts, frac % 1.0
    except (ValueError, IndexError, OSError):
        return None


def fingerprint(lat: np.ndarray, species: list[str], counts: list[int], frac: np.ndarray) -> str:
    """Order-invariant fingerprint: sorted metal-metal nearest-neighbour pair counts + rounded cell metric."""
    labels = np.repeat(species, counts)
    metals = [i for i, s in enumerate(labels) if s != "B"]
    cart = frac @ lat
    pairs: Counter[str] = Counter()
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) @ lat
    for ii, i in enumerate(metals):
        for j in metals[ii + 1:]:
            d = np.min(np.linalg.norm(cart[j] - cart[i] + shifts, axis=1))
            if d < 3.6:  # in-plane (~a) and out-of-plane (~c) metal neighbours
                tag = "ab" if d < 3.3 else "c"
                pairs["-".join(sorted((labels[i], labels[j]))) + tag] += 1
    metric = np.round(np.sort(np.linalg.norm(lat, axis=1)), 1)
    key = f"{sorted(pairs.items())}|{metric.tolist()}|{sorted(zip(species, counts))}"
    return hashlib.sha1(key.encode()).hexdigest()[:12]


def classify(path: Path, species: list[str], counts: list[int]) -> str:
    p = str(path)
    if "sqs_lev" in p:
        return "SQS (BLADE/ATAT)"
    if re.search(r"mp-\d+", p):
        return "MP structure (ORB-relaxed)"
    return "other"


def main() -> None:
    rows = []
    for root in ROOTS:
        base = (BLADE / root).resolve()
        if not base.exists():
            continue
        for contcar in base.rglob("CONTCAR"):
            parsed = read_poscar(contcar)
            if parsed is None:
                continue
            lat, species, counts, frac = parsed
            spset = set(species)
            for sid, allowed in SYSTEMS.items():
                if not spset <= allowed or spset == {"B"}:
                    continue
                comp = dict(zip(species, counts))
                n_metal = sum(v for k, v in comp.items() if k != "B")
                energy_file = contcar.parent / "energy"
                energy = float(energy_file.read_text().split()[0]) if energy_file.exists() else np.nan
                solute = "Nb" if sid == "ZrNb" else "Ta"
                rows.append({
                    "system": sid, "path": str(contcar.parent.relative_to(BLADE.parent)),
                    "origin": classify(contcar, species, counts), "n_atoms": sum(counts),
                    "n_Zr": comp.get("Zr", 0), f"n_solute": comp.get(solute, 0), "n_B": comp.get("B", 0),
                    "x_solute": comp.get(solute, 0) / n_metal if n_metal else np.nan,
                    "B_per_metal": comp.get("B", 0) / n_metal if n_metal else np.nan,
                    "has_B_vacancy_or_offstoich": bool(n_metal and abs(comp.get("B", 0) / n_metal - 2) > 1e-6),
                    "energy_eV_cell": energy, "energy_eV_atom": energy / sum(counts) if energy == energy else np.nan,
                    "fingerprint": fingerprint(lat, species, counts, frac),
                })
    df = pd.DataFrame(rows).drop_duplicates(subset=["system", "path"])
    df["n_copies_same_fingerprint"] = df.groupby("fingerprint").path.transform("count")
    df = df.sort_values(["system", "x_solute", "origin", "fingerprint", "path"])
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "s00_existing_structure_audit.csv", index=False)

    uniq = df.drop_duplicates("fingerprint")
    lines = ["# Step 0 audit: reusable Zr–Nb–B / Zr–Ta–B structures", "",
             f"Scanned roots: {', '.join(ROOTS)}. {len(df)} cells found, {len(uniq)} distinct geometries "
             "(fingerprint = metal-pair neighbour counts + cell metric).", ""]
    for sid in SYSTEMS:
        u = uniq[uniq.system == sid]
        lines.append(f"## {sid}")
        lines.append("")
        lines.append("| origin | x_solute | n_atoms | B/M | distinct realizations | copies on disk | E (eV/atom) |")
        lines.append("|---|---|---|---|---|---|---|")
        for (orig, x, n, bm), g in u.groupby(["origin", "x_solute", "n_atoms", "B_per_metal"]):
            copies = int(df[df.fingerprint.isin(g.fingerprint)].shape[0])
            e = ", ".join(f"{v:.5f}" for v in g.energy_eV_atom.dropna().unique()[:4])
            lines.append(f"| {orig} | {x:.3f} | {n} | {bm:.3f} | {len(g)} | {copies} | {e} |")
        lines.append("")
    vac = df[df.has_B_vacancy_or_offstoich & (df.n_Zr + df.n_solute > 0) & (df.n_B > 0)]
    lines.append(f"Boron-deficient / off-stoichiometric metal-boride cells found: **{len(vac)}**.")
    lines.append(vac[["system", "path", "B_per_metal"]].to_string(index=False) if len(vac) else "")
    (OUT / "s00_audit_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
