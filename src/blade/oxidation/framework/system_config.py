"""System configuration and table preparation for N-metal phase models.

Handles binary and ternary mixed-phase systems uniformly.
Metal detection comes from the configured structure subdirectory, and
fixed-phase data comes from a configurable companion subdirectory.
"""

import csv
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List

import numpy as np


def _stoichiometric_suffix(element: str | None, stoichiometry: float) -> str:
    if not element or stoichiometry <= 0:
        return ""
    coefficient = "" if abs(float(stoichiometry) - 1.0) < 1e-12 else f"{stoichiometry:g}"
    return f"{element}{coefficient}"


@dataclass
class SystemConfig:
    metals: List[str]
    phase_label: str
    tag: str
    tables_dir: str  # absolute path
    phase_element: str | None = None
    phase_element_stoichiometry: float = 0.0
    blade_subdir: str = "blade"
    fixed_phases_subdir: str = "fixed_phases"

    @property
    def n_metals(self) -> int:
        return len(self.metals)

    @property
    def is_binary(self) -> bool:
        return self.n_metals == 2

    @property
    def is_ternary(self) -> bool:
        return self.n_metals == 3

    # ---- File paths -------------------------------------------------------
    @property
    def tables(self) -> Path:
        return Path(self.tables_dir)

    @property
    def phase_table_file(self) -> Path:
        return self.tables / f"{self.tag}_phase_table.csv"

    @property
    def phase_grid_file(self) -> Path:
        return self.tables / f"{self.tag}_phase_energy_points.csv"

    @property
    def interaction_coeff_file(self) -> Path:
        return self.tables / f"{self.tag}_phase_interaction_coefficients.csv"

    @property
    def phase_fit_file(self) -> Path:
        return self.tables / f"{self.tag}_phase_fit_points.csv"

    @property
    def phase_grid_0k_file(self) -> Path:
        return self.tables / f"{self.tag}_phase_grid_0K.csv"

    @property
    def phase_model_file(self) -> Path:
        return self.tables / f"{self.tag}_phase_model.txt"

    @property
    def config_file(self) -> Path:
        return self.tables / f"{self.tag}_system_config.json"

    def tables_ready(self) -> bool:
        return self.phase_table_file.exists() and self.phase_grid_file.exists() and self.interaction_coeff_file.exists()

    def missing_table_files(self) -> list[Path]:
        files = [self.phase_table_file, self.phase_grid_file, self.interaction_coeff_file]
        return [p for p in files if not p.exists()]

    def save(self) -> None:
        self.tables.mkdir(parents=True, exist_ok=True)
        with open(self.config_file, "w") as f:
            json.dump(asdict(self), f, indent=2)

    # ---- Construction -----------------------------------------------------
    @classmethod
    def resolve(
        cls,
        data_root,
        tables_dir,
        metals=None,
        blade_subdir: str = "blade",
        fixed_phases_subdir: str = "fixed_phases",
        phase_element: str | None = None,
        phase_element_stoichiometry: float = 0.0,
    ) -> "SystemConfig":
        """Find or build SystemConfig. Priority:
        1. Existing *_system_config.json matching metals (if given)
        2. Auto-detect from data_root/blade_subdir/
        3. Auto-detect from existing *_phase_table.csv columns
        """
        tables_dir = Path(tables_dir)
        tables_dir.mkdir(parents=True, exist_ok=True)
        blade_subdir = str(blade_subdir)
        fixed_phases_subdir = str(fixed_phases_subdir)
        phase_element = str(phase_element).strip() if phase_element else None
        phase_element_stoichiometry = float(phase_element_stoichiometry) if phase_element else 0.0
        if phase_element and phase_element_stoichiometry <= 0:
            raise ValueError("phase_element_stoichiometry must be positive when phase_element is set")

        # Normalise metals argument
        if metals is not None:
            metals = sorted(metals) if not isinstance(metals, list) else metals

        # 1. Try existing config JSON
        for cfg_file in sorted(tables_dir.glob("*_system_config.json")):
            with open(cfg_file) as f:
                data = json.load(f)
            cfg_metals = data.get("metals", [])
            if metals and sorted(metals) != sorted(cfg_metals):
                continue
            data["tables_dir"] = str(tables_dir)
            data.setdefault("blade_subdir", blade_subdir)
            data.setdefault("fixed_phases_subdir", fixed_phases_subdir)
            # Run-file chemistry is authoritative, including for old caches.
            data["phase_element"] = phase_element
            data["phase_element_stoichiometry"] = phase_element_stoichiometry
            suffix = _stoichiometric_suffix(phase_element, phase_element_stoichiometry)
            data["tag"] = "".join(m.lower() for m in cfg_metals) + (phase_element or "").lower()
            data["phase_label"] = "".join(cfg_metals) + suffix
            return cls(**data)

        # 2. Auto-detect from DATA_ROOT
        if data_root is not None:
            blade = Path(data_root) / blade_subdir
            if blade.is_dir():
                detected = cls._detect_metals_from_dirs(blade, {phase_element} if phase_element else set())
                if len(detected) >= 2:
                    metals = metals or detected

        # 3. Auto-detect from phase table columns
        if not metals:
            for tbl in sorted(tables_dir.glob("*_phase_table.csv")):
                detected = cls._detect_metals_from_table(tbl, {phase_element} if phase_element else set())
                if len(detected) >= 2:
                    metals = detected
                    break

        if not metals or len(metals) < 2:
            raise ValueError(
                f"Cannot determine metals. Set DATA_ROOT to a folder with {blade_subdir!r} subdirs, "
                "or set metals explicitly."
            )

        suffix = _stoichiometric_suffix(phase_element, phase_element_stoichiometry)
        tag = "".join(m.lower() for m in metals) + (phase_element or "").lower()
        phase_label = "".join(metals) + suffix
        return cls(
            metals=metals,
            phase_label=phase_label,
            tag=tag,
            tables_dir=str(tables_dir),
            phase_element=phase_element,
            phase_element_stoichiometry=phase_element_stoichiometry,
            blade_subdir=blade_subdir,
            fixed_phases_subdir=fixed_phases_subdir,
        )

    # ---- Detection helpers ------------------------------------------------
    @staticmethod
    def _detect_metals_from_dirs(blade_dir: Path, excluded_elements=None) -> list:
        """Parse a_X=N tokens from blade/ subdir names → sorted metal list."""
        pattern = re.compile(r"a_([A-Za-z]{1,3})=[0-9.]+")
        metals: set = set()
        excluded = {str(e).upper() for e in (excluded_elements or set())}
        for d in blade_dir.iterdir():
            if not d.is_dir():
                continue
            for m in pattern.finditer(d.name):
                elem = m.group(1)
                if elem.upper() not in excluded:
                    metals.add(elem)
        return sorted(metals)

    @staticmethod
    def _detect_metals_from_table(phase_table_file: Path, excluded_elements=None) -> list:
        """Read column names from *_phase_table.csv → metal columns."""
        non_metal = {"phase_id", "O", "energy_eV_per_atom"}
        non_metal.update(excluded_elements or set())
        with open(phase_table_file, newline="") as f:
            header = next(csv.reader(f))
        return sorted(c for c in header if c not in non_metal)


# ---------------------------------------------------------------------------
# Table preparation — unified for N metals
# ---------------------------------------------------------------------------


def prepare_tables(data_root: Path, cfg: SystemConfig, rk_order: int = 3) -> None:
    """Read raw structure folders, fit the interaction model, and write CSVs."""
    data_root = Path(data_root)
    blade_dir = data_root / cfg.blade_subdir
    fixed_dir = data_root / cfg.fixed_phases_subdir
    metals = cfg.metals
    n_metals = cfg.n_metals
    phase_element = cfg.phase_element
    phase_stoich = cfg.phase_element_stoichiometry
    elements = metals + ([phase_element] if phase_element else []) + ["O"]

    cfg.tables.mkdir(parents=True, exist_ok=True)

    # ---- Parse SQS structures -------------------------------------------
    # Each structure: {y_metals: array(n_metals), energy, total_atoms, formula_units}
    all_rows = []
    metal_patterns = [re.compile(rf"a_{m}=([0-9.]+)") for m in metals]

    for d in blade_dir.iterdir():
        if not d.is_dir():
            continue
        ep = d / "energy"
        cp = d / "CONTCAR"
        if not ep.exists() or not cp.exists():
            continue

        # Parse metal counts from directory name
        metal_counts = []
        for pat in metal_patterns:
            m = pat.search(d.name)
            metal_counts.append(float(m.group(1)) if m else 0.0)
        total_metal = sum(metal_counts)
        if total_metal <= 0:
            continue

        y = [c / total_metal for c in metal_counts]  # fractional occupancies

        counts = _read_contcar_counts(cp, elements)
        n_atoms = sum(counts.values())
        n_fu = sum(counts[m] for m in metals)  # formula units
        if n_fu <= 0 or n_atoms <= 0:
            continue

        e_total = _read_scalar(ep)
        row = {"phase_id": d.name}
        for i, metal in enumerate(metals):
            row[f"y_{metal}_metal_site"] = y[i]
        row["energy_eV_per_formula"] = e_total / n_fu
        row["energy_eV_per_atom"] = e_total / n_atoms
        row["source_total_energy_eV"] = e_total
        row["source_total_atoms"] = int(n_atoms)
        row["source_formula_units"] = int(n_fu)
        all_rows.append(row)

    if not all_rows:
        raise RuntimeError(
            f"No valid structures in {blade_dir}. " f"Expected directories containing a_{metals[0]}=3,a_{metals[1]}=1"
        )

    y_cols = [f"y_{m}_metal_site" for m in metals]
    all_rows.sort(key=lambda r: [r[c] for c in y_cols])

    phase_grid_keys = (
        ["phase_id"]
        + y_cols
        + [
            "energy_eV_per_formula",
            "energy_eV_per_atom",
            "source_total_energy_eV",
            "source_total_atoms",
            "source_formula_units",
        ]
    )
    _write_csv(cfg.phase_grid_file, phase_grid_keys, all_rows)
    print(f"  Wrote {cfg.phase_grid_file}  ({len(all_rows)} structures)")

    # ---- Endpoint energies ----------------------------------------------
    h_ep = np.zeros(n_metals)
    for i, m in enumerate(metals):
        eye = [0.0] * n_metals
        eye[i] = 1.0
        eps = [
            r["energy_eV_per_formula"]
            for r in all_rows
            if all(abs(r[y_cols[k]] - eye[k]) < 1e-10 for k in range(n_metals))
        ]
        if not eps:
            raise ValueError(f"Missing pure endpoint y_{m}=1 in {blade_dir} data")
        h_ep[i] = eps[0]

    # ---- Global simultaneous fit of all Muggianu parameters ------------
    # Uses ALL SQS data (binary sub-system + ternary) in a single lstsq.
    # Each structure contributes to EVERY pair's RK parameters, not just
    # the binary sub-system it nominally belongs to.
    #
    # Model (linear in all L parameters):
    #   H_excess = Σ_{i<j} y_i·y_j · Σ_k L_k^{ij}·(y_i-y_j)^k
    #            + y_0·y_1·y_2 · L^{012}   (ternary only)
    #
    # Endpoints are excluded from the fit (H_excess = 0 by construction).
    binary_coeffs: dict = {}
    ternary_coeff = 0.0
    rk_rows = []

    pairs = [(i, j) for i in range(n_metals) for j in range(i + 1, n_metals)]
    n_pair_params = len(pairs) * (rk_order + 1)
    has_ternary = n_metals == 3
    n_params = n_pair_params + (1 if has_ternary else 0)

    basis_rows, targets = [], []
    for r in all_rows:
        y = np.array([r[y_cols[k]] for k in range(n_metals)])
        # Skip pure endpoints — they pin H_ep, not interaction params
        if np.any(y >= 1.0 - 1e-10):
            continue
        h_lin = float(np.dot(y, h_ep))
        h_exc = r["energy_eV_per_formula"] - h_lin

        basis = []
        for i, j in pairs:
            yi, yj = y[i], y[j]
            z = yi - yj
            z_pow = 1.0
            for _ in range(rk_order + 1):
                basis.append(yi * yj * z_pow)
                z_pow *= z
        if has_ternary:
            basis.append(y[0] * y[1] * y[2])

        basis_rows.append(basis)
        targets.append(h_exc)

    if basis_rows:
        A = np.array(basis_rows)
        b = np.array(targets)
        coeffs_all, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    else:
        coeffs_all = np.zeros(n_params)

    # Unpack fitted coefficients
    idx = 0
    for i, j in pairs:
        c_ij = coeffs_all[idx : idx + rk_order + 1]
        binary_coeffs[(i, j)] = c_ij
        pair_str = f"{metals[i]}-{metals[j]}"
        for k, c in enumerate(c_ij):
            rk_rows.append(
                {
                    "pair": pair_str,
                    "term": f"L{k}",
                    "value_eV_per_formula": float(c),
                    "value_eV_per_atom": float(c) / (1.0 + phase_stoich),
                }
            )
        idx += rk_order + 1

    if has_ternary:
        ternary_coeff = float(coeffs_all[idx])
        rk_rows.append(
            {
                "pair": "-".join(metals),
                "term": "L0",
                "value_eV_per_formula": ternary_coeff,
                "value_eV_per_atom": ternary_coeff / (1.0 + phase_stoich),
            }
        )
        print(
            f'  Ternary L^{{{",".join(metals)}}} = {ternary_coeff:.6f} eV  '
            f"(global fit over {len(basis_rows)} structures)"
        )

    _write_csv(cfg.interaction_coeff_file, ["pair", "term", "value_eV_per_formula", "value_eV_per_atom"], rk_rows)
    print(f"  Wrote {cfg.interaction_coeff_file}")

    # ---- RK fit points -------------------------------------------------
    fit_rows = []
    for r in all_rows:
        y_vec = np.array([r[y_cols[k]] for k in range(n_metals)])
        h_fit = muggianu_energy_nd_scalar(y_vec, h_ep, binary_coeffs, ternary_coeff)
        fit_rows.append(
            {
                "phase_id": r["phase_id"],
                **{y_cols[k]: r[y_cols[k]] for k in range(n_metals)},
                "energy_eV_per_formula": r["energy_eV_per_formula"],
                "rk_fit_eV_per_formula": h_fit,
                "residual_eV_per_formula": r["energy_eV_per_formula"] - h_fit,
            }
        )
    _write_csv(
        cfg.phase_fit_file,
        ["phase_id"] + y_cols + ["energy_eV_per_formula", "rk_fit_eV_per_formula", "residual_eV_per_formula"],
        fit_rows,
    )
    print(f"  Wrote {cfg.phase_fit_file}")

    # ---- Fixed phases ---------------------------------------------------
    phase_rows = []
    if not fixed_dir.is_dir():
        raise FileNotFoundError(
            f"Fixed-phase directory not found: {fixed_dir}. "
            f"Set fixed_phases_subdir in run_simple.py to the correct folder name."
        )
    for d in fixed_dir.iterdir():
        if not d.is_dir():
            continue
        ep = d / "energy"
        cp = d / "CONTCAR"
        if not ep.exists() or not cp.exists():
            continue
        counts = _read_contcar_counts(cp, elements)
        n_atoms = sum(counts.values())
        if n_atoms <= 0:
            continue
        reduced = _reduce_counts(counts)
        e_total = _read_scalar(ep)

        # Exclude line compounds incompatible with the configured flexible phase.
        total_metal = sum(reduced[m] for m in metals)
        phase_element_count = reduced[phase_element] if phase_element else 0.0
        o_count = reduced["O"]
        if phase_element_count > 0 and o_count == 0 and total_metal > 0:
            element_per_metal = phase_element_count / total_metal
            if abs(element_per_metal - phase_stoich) > 0.3:
                continue

        row = {"phase_id": d.name}
        for m in metals:
            row[m] = reduced[m]
        if phase_element:
            row[phase_element] = reduced[phase_element]
        row["O"] = reduced["O"]
        row["energy_eV_per_atom"] = e_total / n_atoms
        phase_rows.append(row)

    if not phase_rows:
        raise RuntimeError(f"No valid fixed phases in {fixed_dir}")
    phase_rows.sort(key=lambda r: r["phase_id"])
    phase_keys = ["phase_id"] + metals + ([phase_element] if phase_element else []) + ["O", "energy_eV_per_atom"]
    _write_csv(cfg.phase_table_file, phase_keys, phase_rows)
    print(f"  Wrote {cfg.phase_table_file}  ({len(phase_rows)} phases)")

    # ---- Human-readable model summary ----------------------------------
    lines = [
        f'Phase: {"_".join(f"{m}_y{i}" for i,m in enumerate(metals))} target',
        "",
        "RK/Muggianu model:",
        "  H(y) = sum_i y_i * H_i + sum_{i<j} y_i*y_j * sum_k L_k^{ij}*(y_i-y_j)^k",
    ]
    if n_metals == 3:
        lines.append("       + y0*y1*y2 * L^{012}")
    lines += ["", "Endpoints (eV/formula):"] + [f"  H_{m} = {h_ep[i]}" for i, m in enumerate(metals)]
    lines += ["", "G(y,T) = H(y) + k_B T * sum_i y_i*ln(y_i)", "", "k_B = 8.617333262145e-5 eV/K"]
    cfg.phase_model_file.write_text("\n".join(lines) + "\n")
    print(f"  Wrote {cfg.phase_model_file}")

    cfg.save()
    print(f"  Wrote {cfg.config_file}")


# ---------------------------------------------------------------------------
# Muggianu scalar helper (for fitting)
# ---------------------------------------------------------------------------


def muggianu_energy_nd_scalar(y_vec, h_ep, binary_coeffs, ternary_coeff):
    """Scalar version for single composition point."""
    H = float(np.dot(y_vec, h_ep))
    for (i, j), L in binary_coeffs.items():
        yi, yj = y_vec[i], y_vec[j]
        z, zpow, poly = yi - yj, 1.0, 0.0
        for c in L:
            poly += c * zpow
            zpow *= z
        H += yi * yj * poly
    if len(y_vec) >= 3 and ternary_coeff != 0.0:
        H += ternary_coeff * y_vec[0] * y_vec[1] * y_vec[2]
    return H


# ---------------------------------------------------------------------------
# Binary RK fitting
# ---------------------------------------------------------------------------


def _fit_rk_binary(points, h_end0, h_end1, order):
    """Fit RK coefficients for a binary pair.

    points: list of (y0, energy_eV_per_formula) where y0 = fraction of metal 0.
    h_end0, h_end1: endpoint energies at y0=0 and y0=1 respectively.
    """
    rows = []
    for y0, e in points:
        if y0 <= 1e-10 or y0 >= 1.0 - 1e-10:
            continue
        z = 2.0 * y0 - 1.0
        pre = y0 * (1.0 - y0)
        basis = [pre * z**k for k in range(order + 1)]
        excess = e - ((1.0 - y0) * h_end0 + y0 * h_end1)
        rows.append({"basis": basis, "excess": excess})

    if not rows or len(rows) < order + 1:
        # Fall back to fewer terms
        actual_order = max(0, len(rows) - 1)
        if actual_order < order:
            return _fit_rk_binary(points, h_end0, h_end1, actual_order)

    n = order + 1
    A = [[sum(r["basis"][i] * r["basis"][j] for r in rows) for j in range(n)] for i in range(n)]
    b = [sum(r["basis"][i] * r["excess"] for r in rows) for i in range(n)]
    return _solve_linear(A, b)


def _solve_linear(A, b):
    n = len(b)
    A = [row[:] for row in A]
    b = b[:]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(A[r][col]))
        A[col], A[pivot] = A[pivot], A[col]
        b[col], b[pivot] = b[pivot], b[col]
        if abs(A[col][col]) < 1e-14:
            raise ValueError("Singular matrix in RK fit")
        for row in range(col + 1, n):
            f = A[row][col] / A[col][col]
            A[row] = [A[row][j] - f * A[col][j] for j in range(n)]
            b[row] -= f * b[col]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]
    return x


# ---------------------------------------------------------------------------
# File I/O helpers
# ---------------------------------------------------------------------------


def _read_scalar(path: Path) -> float:
    return float(path.read_text().strip().split()[0])


def _read_contcar_counts(path: Path, elements: list) -> dict:
    lines = path.read_text().splitlines()
    if len(lines) < 7:
        raise ValueError(f"CONTCAR too short: {path}")
    syms = lines[5].split()
    counts = list(map(int, lines[6].split()))
    result = {e: 0 for e in elements}
    for s, c in zip(syms, counts):
        if s in result:
            result[s] += c
    return result


def _gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _reduce_counts(counts: dict) -> dict:
    vals = [v for v in counts.values() if v > 0]
    if not vals:
        return counts
    g = vals[0]
    for v in vals[1:]:
        g = _gcd(g, v)
    if g <= 1:
        return counts
    return {k: v // g for k, v in counts.items()}


def _write_csv(path: Path, fieldnames: list, rows: list) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
