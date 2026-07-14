"""Core thermodynamic model for N-metal mixed-phase oxidation.

Supports binary and ternary mixed-phase systems.
Grand-potential LP minimization over an N-metal phase simplex
grid plus fixed line-compound phases.

Phase free energy (Muggianu model):
    G(y, T) = H_Muggianu(y) + k_B T * sum_i y_i * ln(y_i)
"""

import re

import numpy as np
import pandas as pd
from scipy.optimize import linprog

KB = 8.617333262145e-5
R = 8.31446261815324
EV_KJ_PER_MOL = 96.48533212331002

_re_phase_y = re.compile(r"y=([0-9.]+)$")
_re_phase_ynd = re.compile(r"([A-Z][a-z]?)=([0-9.]+)")


# ---------------------------------------------------------------------------
# Simplex grid
# ---------------------------------------------------------------------------


def simplex_grid_nd(n_metals: int, step: float = 0.01) -> np.ndarray:
    """Composition grid on the (n_metals-1)-simplex. Returns (n_pts, n_metals)."""
    if n_metals < 2:
        raise ValueError(f"simplex_grid_nd needs at least 2 metals; got {n_metals}")
    if n_metals == 2:
        y0 = np.arange(0.0, 1.0 + step / 2.0, step)
        if abs(y0[-1] - 1.0) > 1e-12:
            y0 = np.append(y0, 1.0)
        return np.column_stack([y0, 1.0 - y0])
    elif n_metals == 3:
        pts = []
        n = int(round(1.0 / step))
        for i in range(n + 1):
            y0 = i / n
            for j in range(n - i + 1):
                y1 = j / n
                y2 = 1.0 - y0 - y1
                if y2 < -1e-10:
                    break
                pts.append([y0, y1, max(0.0, y2)])
        return np.array(pts)

    n = int(round(1.0 / step))
    if n <= 0:
        raise ValueError(f"Invalid simplex step: {step}")

    pts = []

    def _fill(prefix, remaining_units, remaining_dims):
        if remaining_dims == 1:
            pts.append(prefix + [remaining_units / n])
            return
        for units in range(remaining_units + 1):
            _fill(prefix + [units / n], remaining_units - units, remaining_dims - 1)

    _fill([], n, n_metals)
    return np.array(pts)


# ---------------------------------------------------------------------------
# N-metal thermodynamic functions
# ---------------------------------------------------------------------------


def muggianu_energy_nd(
    y_mat: np.ndarray, h_endpoints: np.ndarray, binary_coeffs: dict, ternary_coeff: float = 0.0
) -> np.ndarray:
    """Muggianu mixing enthalpy for an N-metal phase family (eV/formula unit).

    y_mat:         (n_pts, n_metals), rows sum to 1
    h_endpoints:   (n_metals,) endpoint energies
    binary_coeffs: {(i,j): RK_array} for each pair i<j
    ternary_coeff: single L^{012} ternary interaction
    """
    H = y_mat @ h_endpoints
    for (i, j), L in binary_coeffs.items():
        yi, yj = y_mat[:, i], y_mat[:, j]
        z, z_pow, poly = yi - yj, np.ones(len(y_mat)), np.zeros(len(y_mat))
        for c in L:
            poly += c * z_pow
            z_pow *= z
        H += yi * yj * poly
    if y_mat.shape[1] >= 3 and ternary_coeff != 0.0:
        H += ternary_coeff * y_mat[:, 0] * y_mat[:, 1] * y_mat[:, 2]
    return H


def ideal_mixing_nd(y_mat: np.ndarray) -> np.ndarray:
    """N-component ideal mixing entropy shape: Σ_i y_i*ln(y_i). Returns (n_pts,)."""
    s = np.zeros(len(y_mat))
    for j in range(y_mat.shape[1]):
        y = y_mat[:, j]
        m = y > 0.0
        s[m] += y[m] * np.log(y[m])
    return s


# ---------------------------------------------------------------------------
# LP solver
# ---------------------------------------------------------------------------


def solve_grand_lp_batch(A_eq, b_eq, grand_stack):
    """Solve an objective stack without per-batch thread-pool overhead.

    Each row still uses the exact same cold LP solve as :func:`solve_grand_lp`,
    preserving deterministic phase amounts for degenerate optima.
    """
    return [solve_grand_lp(A_eq, b_eq, grand) for grand in grand_stack]


def solve_grand_lp(A_eq, b_eq, grand):
    """Minimize grand @ n s.t. A_eq @ n = b_eq, n >= 0. Returns (amounts, obj, ok)."""
    bounds = [(0.0, None)] * len(grand)
    opts = {"disp": False}
    # revised simplex is ~1.5× faster for small dense LPs; fall back to HiGHS if removed
    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = linprog(grand, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="revised simplex", options=opts)
    except Exception:
        result = linprog(grand, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs", options=opts)
    if result.status == 0:
        x = np.maximum(result.x, 0.0)
        x[x < 1e-10] = 0.0
        return x, float(result.fun), True
    return np.full(len(grand), np.nan), np.nan, False


# ---------------------------------------------------------------------------
# Oxygen thermochemistry
# ---------------------------------------------------------------------------


def _oxygen_shomate_coeffs(T):
    if 100 <= T < 700:
        return 31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0
    elif 700 <= T < 2000:
        return 30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0
    elif 2000 <= T <= 6000:
        return 20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185, 0.0
    else:
        raise ValueError(f"Shomate O2 valid 100-6000 K; got T={T}")


def _oxygen_delta_g0(T):
    A, B, C, D, E, F, G, Hc = _oxygen_shomate_coeffs(T)
    t = T / 1000.0
    Hr = A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F - Hc
    S = A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G
    T0 = 298.15
    A0, B0, C0, D0, E0, _, G0, _ = _oxygen_shomate_coeffs(T0)
    t0 = T0 / 1000.0
    S0 = A0 * np.log(t0) + B0 * t0 + C0 * t0**2 / 2 + D0 * t0**3 / 3 - E0 / (2 * t0**2) + G0
    return Hr - T * S / 1000.0 + T0 * S0 / 1000.0


def mu_o_from_log10_po2(T, log10_po2, mu_o_offset=0.0):
    return (_oxygen_delta_g0(T) + R * T * np.log(10.0) * log10_po2 / 1000.0) / (2.0 * EV_KJ_PER_MOL) + mu_o_offset


def log10_po2_from_mu_o(T, mu_o, mu_o_offset=0.0):
    return (2.0 * (mu_o - mu_o_offset) * EV_KJ_PER_MOL - _oxygen_delta_g0(T)) / (R * T * np.log(10.0) / 1000.0)


# ---------------------------------------------------------------------------
# Phase data loading — unified N-metal
# ---------------------------------------------------------------------------


def load_phase_data(
    phase_table_file,
    phase_points_file,
    rk_coeff_file,
    metals_or_M1,
    M2=None,
    y_column=None,
    phase_label=None,
    y_step=0.01,
    phase_element=None,
    phase_element_stoichiometry=0.0,
):
    """Load phase data for binary or N-metal system.

    The flexible phase's non-metal element and stoichiometry are explicit inputs.
    """
    if isinstance(metals_or_M1, (list, tuple)):
        metals = list(metals_or_M1)
        if phase_label is None:
            phase_label = "".join(metals)
    else:
        metals = [metals_or_M1, M2]
        if phase_label is None:
            phase_label = f"{metals_or_M1}{M2}"
    return _load_nd(
        phase_table_file,
        phase_points_file,
        rk_coeff_file,
        metals,
        phase_label,
        y_step,
        phase_element,
        phase_element_stoichiometry,
    )


def _load_nd(
    phase_table_file,
    phase_points_file,
    rk_coeff_file,
    metals,
    phase_label,
    y_step,
    phase_element,
    phase_element_stoichiometry,
):
    n = len(metals)
    pt = pd.read_csv(phase_table_file)
    bt = pd.read_csv(phase_points_file)
    rk = pd.read_csv(rk_coeff_file)

    # Fixed phases
    fids = pt["phase_id"].astype(str).values
    fmetal = np.column_stack([pt[m].astype(float).values for m in metals])
    phase_element_stoichiometry = float(phase_element_stoichiometry) if phase_element else 0.0
    fphase_element = (
        pt[phase_element].astype(float).values if phase_element and phase_element in pt.columns else np.zeros(len(pt))
    )
    fO = pt["O"].astype(float).values
    fatoms = fmetal.sum(axis=1) + fphase_element + fO
    fef = pt["energy_eV_per_atom"].astype(float).values * fatoms
    n_fixed = len(fids)

    # Simplex grid
    y_grid = simplex_grid_nd(n, y_step)  # (n_pts, n)
    n_pts = len(y_grid)

    # Endpoint energies from the phase grid CSV
    y_cols = [f"y_{m}_metal_site" for m in metals]
    if all(c in bt.columns for c in y_cols):
        by = bt[y_cols].astype(float).values
    else:
        # Legacy binary CSV: only y_M1_metal_site present; derive y_M2 = 1 - y_M1
        old = [c for c in bt.columns if c.startswith("y_") and "metal" in c]
        y0 = bt[old[0]].astype(float).values
        by = np.column_stack([y0, 1.0 - y0])
    be = bt["energy_eV_per_formula"].astype(float).values

    h_ep = np.zeros(n)
    for i, m in enumerate(metals):
        eye = np.zeros(n)
        eye[i] = 1.0
        mask = np.all(np.abs(by - eye) < 1e-10, axis=1)
        if not np.any(mask):
            raise ValueError(f"Missing pure endpoint y_{m}=1 in {phase_points_file}")
        h_ep[i] = float(be[mask][0])

    # Interaction coefficients
    bc: dict = {}
    tc = 0.0
    if "pair" in rk.columns:
        for pstr, grp in rk.groupby("pair"):
            nd = pstr.count("-")
            if nd == 1:
                a, b = pstr.split("-")
                if a in metals and b in metals:
                    i, j = metals.index(a), metals.index(b)
                    if i > j:
                        i, j = j, i
                    bc[(i, j)] = grp.sort_values("term")["value_eV_per_formula"].values
            elif nd >= 2:
                tc = float(grp["value_eV_per_formula"].iloc[0])
    else:
        bc[(0, 1)] = rk["value_eV_per_formula"].values

    phase_H0 = muggianu_energy_nd(y_grid, h_ep, bc, tc)
    phase_mix = ideal_mixing_nd(y_grid)

    # Phase IDs for the phase grid
    if n == 2:
        phase_ids = np.array([f"{phase_label}_y={y[0]:.4f}" for y in y_grid])
    else:
        phase_ids = np.array(
            [phase_label + "_" + "_".join(f"{metals[k]}={y_grid[r,k]:.4f}" for k in range(n - 1)) for r in range(n_pts)]
        )

    phase_ids = np.concatenate([fids, phase_ids])
    phase_kinds = np.concatenate([np.full(n_fixed, "fixed"), np.full(n_pts, "phase")])
    pmetal = np.vstack([fmetal, y_grid])
    pphase_element = np.concatenate(
        [
            fphase_element,
            phase_element_stoichiometry * np.ones(n_pts),
        ]
    )
    pO = np.concatenate([fO, np.zeros(n_pts)])
    A_eq = np.vstack([pmetal.T, pphase_element[np.newaxis, :]])
    phase_y = np.concatenate([np.full(n_fixed, np.nan), y_grid[:, 0]])
    nan_f = np.full((n_fixed, n), np.nan)
    phase_y_nd = np.vstack([nan_f, y_grid])

    suffix = phase_label[len("".join(metals)) :] if phase_label.startswith("".join(metals)) else ""

    return {
        "phase_ids": phase_ids,
        "phase_kinds": phase_kinds,
        "phase_y": phase_y,
        "phase_y_nd": phase_y_nd,
        "phase_metal_stoich": pmetal,
        "phase_element_stoich": pphase_element,
        "phase_O": pO,
        "fixed_energy_formula": fef,
        "phase_H0": phase_H0,
        "phase_mix_shape": phase_mix,
        "y_grid": y_grid,
        "h_endpoints": h_ep,
        "binary_coeffs": bc,
        "ternary_coeff": tc,
        "A_eq": A_eq,
        "n_fixed": n_fixed,
        "n_metals": n,
        "metals": metals,
        "phase_label": phase_label,
        "phase_element": phase_element,
        "phase_element_stoichiometry": phase_element_stoichiometry,
        "phase_suffix": suffix,
    }


# ---------------------------------------------------------------------------
# Assemblage labeling
# ---------------------------------------------------------------------------


def _short_label(pid):
    return pid.split("_")[-1]


def _phase_composition(pid, metals):
    """Extract the metal-sublattice composition encoded in a pseudo-phase ID."""
    n = len(metals)
    if n == 2:
        m = _re_phase_y.search(pid)
        if not m:
            return None
        y0 = float(m.group(1))
        return np.array([y0, max(0.0, 1.0 - y0)], dtype=float)
    else:
        matches = {s: float(v) for s, v in _re_phase_ynd.findall(pid) if s in metals}
        fracs = [matches.get(m, 0.0) for m in metals[:-1]]
        fracs.append(max(0.0, 1.0 - sum(fracs)))
        return np.asarray(fracs, dtype=float)


def _phase_comp_label_rounded(pid, metals, decimals=1, phase_suffix=""):
    """Label a flexible phase by the elements present above the cutoff."""
    sfx = phase_suffix
    raw = _phase_composition(pid, metals)
    if raw is None:
        return pid
    present = [metal for metal, frac in zip(metals, raw) if frac > 0.01 + 1e-12]
    if len(present) == 1:
        return f"{present[0]}{sfx}"
    return "(" + "".join(present) + ")" + sfx


def _phase_comp_label(pid, metals, phase_suffix=""):
    """Format phase composition label.

    Components at or below 0.01 are omitted. Above that cutoff, composition
    values do not change the phase name; only the set of present elements does.
    """
    sfx = phase_suffix
    fracs = _phase_composition(pid, metals)
    if fracs is None:
        return pid
    present = [metal for metal, frac in zip(metals, fracs) if frac > 0.01 + 1e-12]
    if len(present) == 1:
        return f"{present[0]}{sfx}"
    return "(" + "".join(present) + ")" + sfx


def _phase_component_signature(pid, metals):
    """Displayed formula for elements above the region-presence cutoff."""
    fracs = _phase_composition(pid, metals)
    if fracs is None:
        return str(pid).split("_")[0]
    present = [m for m, f in zip(metals, fracs) if f > 0.01 + 1e-12]
    base = str(pid).split("_")[0]
    metal_prefix = "".join(metals)
    suffix = base[len(metal_prefix) :] if base.startswith(metal_prefix) else ""
    if len(present) == 1:
        return f"{present[0]}{suffix}"
    return f"({''.join(present)}){suffix}"


def _phase_comp_range_label(pids, metals, phase_suffix=""):
    """Label a region's flexible phase by its component signature."""
    compositions = [comp for comp in (_phase_composition(str(pid), metals) for pid in pids) if comp is not None]
    if not compositions:
        return str(pids[0]).split("_")[0] if len(pids) else ""
    comp = np.vstack(compositions)
    present = [metal for i, metal in enumerate(metals) if float(np.max(comp[:, i])) > 0.01 + 1e-12]
    if len(present) == 1:
        return f"{present[0]}{phase_suffix}"
    return f"({''.join(present)}){phase_suffix}"


def format_phase_fraction_comp(name, lo, hi, low_cutoff=0.01, low_mid_cutoff=0.01, range_cutoff=0.02):
    """Format a phase fraction with its composition label.

    Rules:
    - an entire phase range < 0.01 → omitted
    - narrow/constant values → round-half-up to two decimals
    - a range endpoint <= 0.05 → endpoint ``0``
    - an entirely >= 0.95 phase → phase name without a ``1`` prefix
    - lo <= 0.05 and hi >= 0.95 → ``0-1{name}``
    - hi - lo <= 0.02        → rounded average
    - otherwise              → rounded range
    """
    lo = min(1.0, max(0.0, float(lo)))
    hi = min(1.0, max(0.0, float(hi)))

    def key_number(value):
        return np.floor(float(value) * 100.0 + 0.5 + 1e-12) / 100.0

    def key_value(value):
        return f"{key_number(value):.2f}"

    if hi < low_cutoff - 1e-12:
        return ""
    if lo <= 0.05 + 1e-12 and key_number(hi) >= 0.95 - 1e-12:
        return f"0-1{name}"
    if key_number(lo) >= 0.95 - 1e-12:
        return name
    crosses_presence_cutoff = lo < low_cutoff - 1e-12 <= hi
    if (hi - lo) <= range_cutoff + 1e-12 and not crosses_presence_cutoff:
        avg = 0.5 * (lo + hi)
        if key_number(avg) >= 0.95 - 1e-12:
            return name
        return f"{key_value(avg)}{name}"
    lo_text = "0" if lo <= 0.05 + 1e-12 else key_value(lo)
    hi_text = "1" if key_number(hi) >= 0.95 - 1e-12 else key_value(hi)
    return f"{lo_text}-{hi_text}{name}"


def _format_exact_fraction_value(frac, significant_digits=2):
    """Format exact hover fractions with fixed-point significant digits."""
    frac = float(frac)
    if frac == 0.0:
        return "0"
    magnitude = int(np.floor(np.log10(abs(frac))))
    for digits in range(significant_digits, significant_digits + 8):
        decimals = max(0, digits - magnitude - 1)
        txt = f"{frac:.{decimals}f}".rstrip("0").rstrip(".")
        # Do not make an included >0.01 value look like the excluded cutoff.
        if frac > 0.01 and float(txt or 0.0) <= 0.01:
            continue
        return txt if txt else "0"
    return f"{frac:.10f}".rstrip("0").rstrip(".")


def _clean_phase_label(name: str) -> str:
    """Remove minor-element annotations (<=0.05) from a mixed-phase label.

    The formula suffix is preserved while minor components are removed.
    """
    import re

    m = re.match(r"\(([^)]+)\)(.*)$", name)
    if not m:
        return name
    inner, suffix = m.groups()
    # Parse tokens: optional fraction then metal symbol
    tok_re = re.compile(r"(<0\.\d+|0\.0\d*)?([A-Z][a-z]?)")
    kept = []
    dominant = None
    for frac, metal in tok_re.findall(inner):
        if frac == "":
            # dominant metal (no fraction shown) — always keep
            dominant = metal
        elif frac.startswith("<") or float(frac) <= 0.05:
            continue  # drop
        else:
            kept.append(f"{frac}{metal}")
    if dominant is None:
        return name
    if not kept:
        return f"{dominant}{suffix}"
    return "(" + dominant + "".join(kept) + ")" + suffix


def format_phase_detail_line(phase_ranges):
    """Format region phase ranges into one compact hover/legend line.

    Composition labels already apply the shared minor-component rules. Entries
    with the same displayed composition are merged into one phase range.
    """
    if not phase_ranges:
        return ""
    # Merge pseudo-phases that format to the same displayed composition.
    from collections import defaultdict

    merged: dict[str, list] = defaultdict(lambda: [1.0, 0.0])
    for name, lo, hi in phase_ranges:
        clean = str(name)
        merged[clean][0] = min(merged[clean][0], float(lo))
        merged[clean][1] = max(merged[clean][1], float(hi))
    ordered = sorted(
        merged.items(),
        key=lambda kv: 0.5 * (kv[1][0] + kv[1][1]),
        reverse=True,
    )
    parts = []
    for clean_name, (lo, hi) in ordered:
        formatted = format_phase_fraction_comp(clean_name, lo, hi)
        if formatted:
            parts.append(formatted)
    return " + ".join(parts)


def format_exact_phase_fraction_line(fracs, phase_ids, phase_kinds, metals, min_frac=1e-12, phase_suffix=""):
    """Format exact phase fractions for a single state.

    Returns a compact exact-fraction string for one equilibrium state.
    """
    from collections import defaultdict

    totals = defaultdict(float)
    for pid, kind, frac in zip(phase_ids, phase_kinds, fracs):
        if np.isnan(frac) or float(frac) <= min_frac:
            continue
        name = _phase_comp_label(pid, metals, phase_suffix=phase_suffix) if kind == "phase" else _short_label(pid)
        totals[name] += float(frac)
    if not totals:
        return ""
    ordered = sorted(totals.items(), key=lambda x: x[1], reverse=True)
    return " + ".join(f"{_format_exact_fraction_value(f)}{name}" for name, f in ordered)


def format_phase_fraction_summary_line(summary):
    """Format a saved `phase_fraction_summary` string into exact hover text."""
    if not isinstance(summary, str) or not summary.strip():
        return ""
    pairs = []
    for tok in summary.split("|"):
        tok = tok.strip()
        if "=" not in tok:
            continue
        name, val = tok.rsplit("=", 1)
        try:
            pairs.append((name.strip(), float(val.strip())))
        except ValueError:
            continue
    if not pairs:
        return ""
    pairs.sort(key=lambda x: x[1], reverse=True)
    return " + ".join(f"{_format_exact_fraction_value(frac)}{name}" for name, frac in pairs)


def build_assemblage_labels(
    phase_ids,
    phase_kinds,
    fractions,
    threshold,
    phase_label,
    metals=None,
    phase_suffix="",
    family_threshold=None,
    family_values=None,
    family_threshold_inclusive=True,
):
    """Build (exact_label, family_label) for an assemblage.

    exact_label      — full detail, every phase above threshold.
    family_label     — used for region colouring; phases below family_threshold
                       are ignored so trace amounts don't create spurious regions.
    family_threshold — defaults to max(threshold, 0.01) if None.
    family_values    — raw phase amounts used for region presence. Defaults
                       to normalized fractions for backward compatibility.
    """
    if metals is None:
        metals = []
    if family_threshold is None:
        family_threshold = max(threshold, 0.01)

    # Exact label: everything above the LP threshold
    active_exact = fractions > threshold
    aids_e, akinds_e = phase_ids[active_exact], phase_kinds[active_exact]
    if len(aids_e) == 0:
        return "no feasible assemblage", "no feasible assemblage"

    exact = []
    for pid, kind in zip(aids_e, akinds_e):
        exact.append(
            _phase_comp_label(pid, metals, phase_suffix=phase_suffix) if kind == "phase" else _short_label(pid)
        )
    exact_label = " + ".join(exact)

    # Family label: ignore phases below family_threshold to avoid trace-phase splits
    presence_values = fractions if family_values is None else np.asarray(family_values)
    if family_threshold_inclusive:
        active_family = presence_values >= family_threshold - 1e-12
    else:
        active_family = presence_values > family_threshold + 1e-12

    aids_f, akinds_f = phase_ids[active_family], phase_kinds[active_family]
    ffixed, plbls = [], []
    for pid, kind in zip(aids_f, akinds_f):
        if kind == "phase":
            plbls.append(_phase_component_signature(pid, metals))
        else:
            ffixed.append(_short_label(pid))

    fphase = sorted(set(plbls))

    family_parts = sorted(set(ffixed + fphase))
    family_label = " + ".join(family_parts) if family_parts else exact_label
    return exact_label, family_label


def assign_region_ids(family_labels):
    labels, idx = [], {}
    ids = np.zeros(len(family_labels), dtype=int)
    for i, lbl in enumerate(family_labels):
        if lbl not in idx:
            idx[lbl] = len(labels) + 1
            labels.append(lbl)
        ids[i] = idx[lbl]
    return ids, labels


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def grid_edges(values):
    v = np.asarray(values, dtype=float)
    if len(v) == 1:
        return np.array([v[0] - 0.5, v[0] + 0.5])
    mids = 0.5 * (v[:-1] + v[1:])
    return np.concatenate([[v[0] - 0.5 * (v[1] - v[0])], mids, [v[-1] + 0.5 * (v[-1] - v[-2])]])


def region_annotation_text(rid, label, region_details=None, mode="id"):
    """Return either a region number or its formatted phase summary."""
    normalized = str(mode).strip().lower()
    if normalized in {"id", "ids", "number", "numbers"}:
        return str(rid)
    if normalized not in {"phase", "phases", "full"}:
        raise ValueError("region_label_mode must be 'id' or 'phases', " f"not {mode!r}")
    detail = region_details.get(rid) if region_details else None
    text = ""
    if detail:
        phase_ranges = detail.get("phase_ranges") or []
        if phase_ranges:
            text = format_phase_detail_line(phase_ranges)
        else:
            text = detail.get("exact_label") or ""
    text = text or str(label)
    # In-graph summaries are vertical lists: one complete component per row.
    return "\n".join(text.split(" + "))


def add_region_annotation(ax, x, y, text, fontsize, region_color="white"):
    """Add a bordered region label for later collision adjustment."""
    from matplotlib.colors import to_rgb

    darker_color = tuple(0.82 * channel for channel in to_rgb(region_color))
    return ax.text(
        x,
        y,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight="bold",
        zorder=5,
        bbox=dict(fc=darker_color, ec="black", lw=0.8, alpha=1.0, pad=2.4),
    )


def separate_region_annotations(ax, texts):
    """Greedily move labels in display space to reduce visible overlap."""
    if len(texts) < 2:
        return
    canvas = ax.figure.canvas
    canvas.draw()
    renderer = canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer=renderer)
    placed = []
    offsets = [(0, 0)]
    for radius in range(8, 129, 8):
        offsets.extend(
            [
                (0, radius),
                (0, -radius),
                (radius, radius),
                (-radius, radius),
                (radius, -radius),
                (-radius, -radius),
                (radius, 0),
                (-radius, 0),
            ]
        )
    for text in texts:
        original = ax.transData.transform(text.get_position())
        best = None
        for dx, dy in offsets:
            candidate = original + np.array([dx, dy], dtype=float)
            text.set_position(ax.transData.inverted().transform(candidate))
            box = text.get_window_extent(renderer=renderer).expanded(1.10, 1.16)
            outside = (
                max(0.0, axes_box.x0 - box.x0)
                + max(0.0, box.x1 - axes_box.x1)
                + max(0.0, axes_box.y0 - box.y0)
                + max(0.0, box.y1 - axes_box.y1)
            )
            overlap = sum(
                max(0.0, min(box.x1, other.x1) - max(box.x0, other.x0))
                * max(0.0, min(box.y1, other.y1) - max(box.y0, other.y0))
                for other in placed
            )
            score = overlap + outside * 1_000.0 + 0.01 * (abs(dx) + abs(dy))
            if best is None or score < best[0]:
                best = (score, candidate, box)
            if overlap == 0.0 and outside == 0.0:
                break
        _, candidate, box = best
        text.set_position(ax.transData.inverted().transform(candidate))
        placed.append(box)


def write_region_map_html(
    x_values,
    y_values,
    region_id_grid,
    region_labels,
    html_path,
    title,
    x_label,
    y_label,
    region_details=None,
    cell_text_grid=None,
    region_label_mode="id",
    region_label_fontsize=11,
):
    """Write a hoverable HTML version of a region map.

    The HTML shows the same region partition as the PNG and exposes the
    per-region phase-fraction text on hover.
    """
    import hashlib

    import matplotlib.pyplot as plt
    import plotly.graph_objects as go
    from matplotlib.colors import to_hex, to_rgb

    x_values = np.asarray(x_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    region_id_grid = np.asarray(region_id_grid, dtype=int)
    if region_id_grid.shape != (len(y_values), len(x_values)):
        raise ValueError(
            "region_id_grid shape must be (len(y_values), len(x_values)) "
            f"but got {region_id_grid.shape} for {len(y_values)}x{len(x_values)}"
        )

    def _label_color(label):
        cmap = plt.get_cmap("tab20")
        h = int(hashlib.md5(str(label).encode()).hexdigest(), 16)
        return to_hex(cmap(h % 20 / 20))

    def _format_detail(detail):
        if not detail:
            return ""
        phase_ranges = detail.get("phase_ranges") or []
        if phase_ranges:
            return format_phase_detail_line(phase_ranges)
        exact = detail.get("exact_label")
        return exact or ""

    n_reg = len(region_labels)
    colors = [_label_color(lbl) for lbl in region_labels]
    if n_reg == 0:
        colors = ["#d9d9d9"]
        n_reg = 1

    colorscale = []
    for i, color in enumerate(colors, start=1):
        lo = (i - 1) / n_reg
        hi = i / n_reg
        colorscale.append([lo, color])
        colorscale.append([hi, color])

    hover_text = np.full(region_id_grid.shape, "", dtype=object)
    if cell_text_grid is not None:
        hover_text[:] = np.asarray(cell_text_grid, dtype=object)
    else:
        for rid, lbl in enumerate(region_labels, start=1):
            mask = region_id_grid == rid
            detail = _format_detail(region_details.get(rid) if region_details else None)
            text = f"Region {rid}: {lbl}"
            if detail:
                text = f"{text}<br>{detail}"
            hover_text[mask] = text
    fig = go.Figure(
        data=go.Heatmap(
            x=x_values,
            y=y_values,
            z=region_id_grid,
            text=hover_text,
            hovertemplate="%{text}<br>" + x_label + "=%{x:.3f}<br>" + y_label + "=%{y:.3f}<extra></extra>",
            colorscale=colorscale,
            zmin=0.5,
            zmax=n_reg + 0.5,
            showscale=False,
            hoverongaps=False,
        )
    )

    annotations = []
    for rid, lbl in enumerate(region_labels, start=1):
        mask = region_id_grid == rid
        if not np.any(mask):
            continue
        rows, cols = np.where(mask)
        r = max(0, min(len(y_values) - 1, int(np.round(np.median(rows)))))
        c = max(0, min(len(x_values) - 1, int(np.round(np.median(cols)))))
        annotations.append(
            dict(
                x=float(x_values[c]),
                y=float(y_values[r]),
                text=region_annotation_text(rid, lbl, region_details, region_label_mode).replace("\n", "<br>"),
                showarrow=False,
                font=dict(size=region_label_fontsize, color="black"),
                bgcolor=(
                    "rgba("
                    + ",".join(str(round(255 * 0.82 * channel)) for channel in to_rgb(colors[rid - 1]))
                    + ",1.0)"
                ),
                bordercolor="black",
                borderwidth=1,
                borderpad=3,
            )
        )

    fig.update_layout(
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        annotations=annotations,
        margin=dict(l=70, r=30, t=60, b=60),
        hovermode="closest",
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)
    fig.write_html(str(html_path), include_plotlyjs=True, full_html=True)
    return html_path


def plot_region_map(
    mu_o_values,
    x_values,
    region_id_grid,
    region_labels,
    T_value,
    M1,
    M2,
    plot_output_dir,
    region_label_mode="id",
    region_label_fontsize=14,
    boundary_lw=0.8,
    M3=None,
    region_details=None,
    cell_text_grid=None,
    system_label=None,
):
    """muO–x assemblage map with legend panel. Saves PNG and HTML."""
    import matplotlib

    matplotlib.use("Agg")
    import hashlib
    from pathlib import Path

    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    n_reg = len(region_labels)
    # Deterministic color: hash of label string → same assemblage always same color
    # across temperatures, systems, and runs.
    cmap = plt.get_cmap("tab20")  # 20 fixed slots

    def _label_color(label):
        h = int(hashlib.md5(label.encode()).hexdigest(), 16)
        return cmap(h % 20 / 20)

    colors = {r + 1: _label_color(lbl) for r, lbl in enumerate(region_labels)}

    show_region_key = str(region_label_mode).strip().lower() in {"id", "ids", "number", "numbers"}
    if show_region_key:
        fig, (ax, ax_leg) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05})
        ax_leg.axis("off")
    else:
        fig, ax = plt.subplots(figsize=(13, 8))
        ax_leg = None

    # Fast rendering: build RGBA array and draw with pcolormesh (single draw call)
    color_arr = np.array(
        [colors.get(int(rid), (0.9, 0.9, 0.9, 1.0)) for rid in region_id_grid.ravel()], dtype=float
    ).reshape(len(x_values), len(mu_o_values), 4)
    xe = grid_edges(x_values)
    me = grid_edges(mu_o_values)
    ax.pcolormesh(me, xe, color_arr, shading="flat")

    # Draw each unlike-neighbour interface once. Contouring numeric region IDs
    # stacks lines when adjacent IDs differ by more than one.
    if n_reg > 1:
        from matplotlib.collections import LineCollection

        segments = []
        vertical = region_id_grid[:, :-1] != region_id_grid[:, 1:]
        for ix, imu in np.argwhere(vertical):
            segments.append(
                [
                    (me[imu + 1], xe[ix]),
                    (me[imu + 1], xe[ix + 1]),
                ]
            )
        horizontal = region_id_grid[:-1, :] != region_id_grid[1:, :]
        for ix, imu in np.argwhere(horizontal):
            segments.append(
                [
                    (me[imu], xe[ix + 1]),
                    (me[imu + 1], xe[ix + 1]),
                ]
            )
        if segments:
            ax.add_collection(
                LineCollection(
                    segments,
                    colors="black",
                    linewidths=boundary_lw,
                    capstyle="butt",
                    joinstyle="miter",
                    zorder=3,
                )
            )

    annotation_artists = []
    for rid, lbl in enumerate(region_labels, start=1):
        mask = region_id_grid == rid
        if not np.any(mask):
            continue
        rows, cols = np.where(mask)
        r = max(0, min(len(x_values) - 1, int(np.round(np.median(rows)))))
        c = max(0, min(len(mu_o_values) - 1, int(np.round(np.median(cols)))))
        annotation = region_annotation_text(rid, lbl, region_details, region_label_mode)
        annotation_artists.append(
            add_region_annotation(ax, mu_o_values[c], x_values[r], annotation, region_label_fontsize, colors[rid])
        )

    ax.set_xlim(float(mu_o_values[0]), float(mu_o_values[-1]))
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel(r"$\mu_O$ (eV per O atom)", fontsize=13)
    ylabel = f"initial $x_{{{M1}}}$" + (f"  ($x_{{{M2}}}$=$x_{{{M3}}}$)" if M3 else "")
    ax.set_ylabel(ylabel, fontsize=13)
    sys_str = system_label or ("–".join([M1, M2] + ([M3] if M3 else [])))
    ax.set_title(f"{sys_str} assemblage boundaries, T={T_value} K", fontsize=14)
    ax.tick_params(labelsize=11)
    separate_region_annotations(ax, annotation_artists)

    def _format_detail(detail):
        if not detail:
            return ""
        phase_ranges = detail.get("phase_ranges") or []
        if phase_ranges:
            return format_phase_detail_line(phase_ranges)
        exact = detail.get("exact_label")
        return exact or ""

    def _legend_label(rid, lbl):
        if region_details is None or rid not in region_details:
            return lbl
        detail = _format_detail(region_details[rid])
        return detail if detail else lbl

    handles = [
        mpatches.Patch(color=colors[rid], label=f"{rid}: {_legend_label(rid, lbl)}")
        for rid, lbl in enumerate(region_labels, start=1)
        if np.any(region_id_grid == rid)
    ]
    fs = max(5, min(8, int(110 / max(len(handles), 1))))
    if ax_leg is not None:
        ax_leg.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.0, 1.0),
            fontsize=fs,
            framealpha=0.9,
            handlelength=1.2,
            borderaxespad=0,
            title="Region key",
            title_fontsize=fs + 1,
        )

    out = Path(plot_output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "assemblage_region_map.png", dpi=150, bbox_inches="tight")
    try:
        write_region_map_html(
            mu_o_values,
            x_values,
            region_id_grid,
            region_labels,
            out / "assemblage_region_map.html",
            title=f"{sys_str} assemblage boundaries, T={T_value} K",
            x_label="mu_O (eV per O atom)",
            y_label=(f"initial x_{M1}" + (f" (x_{M2}=x_{M3})" if M3 else "")),
            region_details=region_details,
            cell_text_grid=cell_text_grid,
            region_label_mode=region_label_mode,
            region_label_fontsize=region_label_fontsize,
        )
    except Exception as e:
        print(f"[html skipped] {e}")
    plt.close(fig)
    return str(out / "assemblage_region_map.png")
