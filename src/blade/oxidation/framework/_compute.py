"""Computation helpers for SingleCompositionAnalyzer.

Contains elemental-conservation utilities and all LP-scan routines
(run_scan, run_muO_T_map, run_muO_x_map) — extracted from single.py
to keep that module slim.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path

# ---------------------------------------------------------------- b_eq helpers


def _conservation_rhs(self) -> np.ndarray:
    """Build the elemental conservation vector from explicit settings.

    Args:
        self: SingleCompositionAnalyzer instance.

    Returns:
        1-D array concatenating metal fractions and the phase-element stoichiometry.
    """
    return np.concatenate([self._comp, [self.phase_element_stoichiometry]])


def _family_threshold(self) -> float:
    """Return the minimum phase-fraction for family-membership decisions.

    Args:
        self: SingleCompositionAnalyzer instance.

    Returns:
        0.01 when include_0p01_to_0p05_components is True, else 0.05.
    """
    return 0.01 if self.include_0p01_to_0p05_components else 0.05


def _amounts_from_fractions(self, fracs: np.ndarray, b_eq: np.ndarray) -> np.ndarray:
    """Recover the LP amount scale from normalized cached fractions.

    Args:
        self: SingleCompositionAnalyzer instance.
        fracs: Normalized phase-fraction vector (sums to 1).
        b_eq: Elemental conservation right-hand-side vector.

    Returns:
        Rescaled amount vector consistent with the conservation constraints.
    """
    conserved = np.asarray(self.pd_data["A_eq"] @ fracs, dtype=float)
    b_eq = np.asarray(b_eq, dtype=float)
    valid = (np.abs(b_eq) > 1e-12) & (np.abs(conserved) > 1e-12)
    if not np.any(valid):
        return np.asarray(fracs, dtype=float)
    scale = float(np.median(b_eq[valid] / conserved[valid]))
    return np.asarray(fracs, dtype=float) * scale


# ---------------------------------------------------------------- analyses


def run_scan(
    self,
    T_values: list[float] | np.ndarray,
    mu_O_values: np.ndarray,
    skip_if_exists: bool = True,
    active_threshold: float = 1e-9,
    plot_threshold: float = 1e-4,
) -> Path:
    """Run a 1-D oxidation scan: phase fractions vs μO at each temperature.

    Args:
        self: SingleCompositionAnalyzer instance.
        T_values: Sequence of temperatures (K) at which to run the scan.
        mu_O_values: 1-D array of oxygen chemical potentials (eV/O).
        skip_if_exists: Re-use on-disk caches when available.
        active_threshold: Phase fractions below this are zeroed out.
        plot_threshold: Groups with max fraction below this are not plotted.

    Returns:
        Path to the oxidation_scan output directory.
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import pandas as pd

    from .thermodynamics import (
        KB,
        _phase_comp_label,
        _short_label,
        solve_grand_lp,
    )

    mpl.use("Agg")

    self._ensure_path()
    pd_d = self.pd_data
    b_eq = self._conservation_rhs()
    out = self.out_dir / "oxidation_scan"
    cache_root = self.out_dir / "cache" / "oxidation_scan"
    cache_root.mkdir(parents=True, exist_ok=True)
    mu_vals = np.asarray(mu_O_values)
    T_list = list(T_values)
    n_phases = len(pd_d["phase_ids"])

    for T in T_list:
        temp_out = out / f"T{int(T)}"
        temp_out.mkdir(parents=True, exist_ok=True)
        cache = cache_root / f"T{int(T)}.npz"
        legacy_csv = self.out_dir / f"scan_{int(T)}K.csv"
        wide_f = None

        if skip_if_exists and cache.exists():
            try:
                cached = np.load(cache)
                candidate = np.asarray(cached["fractions"], dtype=float)
                grid_matches = "mu_values" not in cached.files or np.array_equal(cached["mu_values"], mu_vals)
                if candidate.shape == (len(mu_vals), n_phases) and grid_matches:
                    wide_f = candidate
                    if "mu_values" not in cached.files:
                        np.savez_compressed(
                            cache,
                            fractions=wide_f,
                            mu_values=mu_vals,
                        )
            except Exception:
                wide_f = None

        # Convert an existing wide CSV cache once without recalculating LPs.
        if wide_f is None and skip_if_exists and legacy_csv.exists():
            try:
                legacy = pd.read_csv(legacy_csv)
                if len(legacy) == len(mu_vals):
                    wide_f = np.column_stack(
                        [
                            (legacy[str(pid)].to_numpy(dtype=float) if str(pid) in legacy.columns else np.zeros(len(mu_vals)))
                            for pid in pd_d["phase_ids"]
                        ]
                    )
                    np.savez_compressed(
                        cache,
                        fractions=wide_f,
                        mu_values=mu_vals,
                    )
            except Exception:
                wide_f = None

        if wide_f is None:
            if not self._run_calculations:
                raise RuntimeError(f"plot-only mode requires the single-scan cache: {cache}")
            phase_G = pd_d["phase_H0"] + KB * float(T) * pd_d["phase_mix_shape"]
            wide_f = np.zeros((len(mu_vals), n_phases), dtype=float)
            for imu, mu_o in enumerate(mu_vals):
                grand = np.concatenate(
                    [
                        pd_d["fixed_energy_formula"] - pd_d["phase_O"][: pd_d["n_fixed"]] * mu_o,
                        phase_G,
                    ]
                )
                amounts, _, ok = solve_grand_lp(pd_d["A_eq"], b_eq, grand)
                if ok:
                    total = float(np.nansum(amounts))
                    fracs = amounts / total if total > 0 else amounts.copy()
                    fracs[np.abs(fracs) < active_threshold] = 0.0
                    wide_f[imu] = fracs
            np.savez_compressed(cache, fractions=wide_f, mu_values=mu_vals)

        if not self._run_plots:
            continue

        grouped_cols = {}
        for idx, (pid, kind) in enumerate(zip(pd_d["phase_ids"], pd_d["phase_kinds"], strict=False)):
            label = (
                _phase_comp_label(
                    str(pid),
                    self.metals,
                    phase_suffix=pd_d.get("phase_suffix", ""),
                )
                if kind == "phase"
                else _short_label(str(pid))
            )
            grouped_cols.setdefault(label, []).append(idx)

        active_groups = []
        for label, cols in grouped_cols.items():
            curve = np.sum(wide_f[:, cols], axis=1)
            if np.nanmax(curve) > plot_threshold:
                active_groups.append((label, curve))
        if active_groups:
            fig, ax = plt.subplots(figsize=(10, 6))
            for label, curve in active_groups:
                ax.plot(mu_vals, curve, lw=1.8, label=label)
            ax.set_xlabel(r"$\mu_O$ (eV per O atom)", fontsize=13)
            ax.set_ylabel("Phase fraction", fontsize=13)
            ax.set_title(
                f"{''.join(self.metals)}{self.phase_element or ''} scan | {self._x_label} | T={int(T)} K",
                fontsize=13,
            )
            ax.set_xlim(float(mu_vals[0]), float(mu_vals[-1]))
            ax.set_ylim(-0.02, 1.02)
            ax.grid(True, alpha=0.35)
            ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=9, borderaxespad=0)
            png = temp_out / "phase_fractions.png"
            fig.savefig(png, dpi=170, bbox_inches="tight")
            plt.close(fig)
            print(f"  scan {int(T)}K → {png}")
    return out


def run_muO_T_map(
    self,
    T_values: np.ndarray,
    mu_O_values: np.ndarray,
    skip_if_exists: bool = True,
    active_threshold: float = 1e-9,
    boundary_lw: float = 0.8,
    region_label_fontsize: int | None = None,
) -> Path:
    """Compute and plot the μO vs T assemblage map at this fixed composition.

    Args:
        self: SingleCompositionAnalyzer instance.
        T_values: 1-D array of temperatures (K).
        mu_O_values: 1-D array of oxygen chemical potentials (eV/O).
        skip_if_exists: Re-use on-disk caches when available.
        active_threshold: Phase fractions below this are zeroed out.
        boundary_lw: Line width for region-boundary lines.
        region_label_fontsize: Override for in-map label font size.

    Returns:
        Path to the muO_T_phase_map output directory.
    """
    from framework._run import _build_region_details

    from .thermodynamics import (
        KB,
        assign_region_ids,
        build_assemblage_labels,
        format_exact_phase_fraction_line,
        solve_grand_lp,
    )

    self._ensure_path()
    pd_d = self.pd_data
    b_eq = self._conservation_rhs()
    out = self.out_dir
    map_out = out / "muO_T_phase_map"
    cache_dir = out / "cache" / "muO_T_phase_map"
    cache_dir.mkdir(parents=True, exist_ok=True)
    T_vals = np.asarray(T_values)
    mu_vals = np.asarray(mu_O_values)
    M1 = self.metals[0]
    M2 = self.metals[1] if len(self.metals) > 1 else ""
    sys_str = "–".join(self.metals + ([self.phase_element] if self.phase_element else []))
    n_T, n_mu = len(T_vals), len(mu_vals)
    n_states = n_T * n_mu
    if region_label_fontsize is None:
        region_label_fontsize = self.region_label_fontsize

    fractions_cache = cache_dir / "phase_fractions.npz"
    legacy_cache = out / "muO_T_map_phase_fractions.npz"
    cache_source = fractions_cache if fractions_cache.exists() else legacy_cache
    if skip_if_exists and cache_source.exists():
        try:
            cached = np.load(cache_source)
            wide_f = cached["fractions"]
            if wide_f.shape != (n_states, len(pd_d["phase_ids"])):
                raise ValueError("cached phase-fraction matrix has the wrong shape")
            if "mu_values" in cached.files and not np.array_equal(cached["mu_values"], mu_vals):
                raise ValueError("cached muO grid does not match")
            if "temperature_values" in cached.files and not np.array_equal(cached["temperature_values"], T_vals):
                raise ValueError("cached temperature grid does not match")
            wide_amounts = np.vstack([_amounts_from_fractions(self, fracs, b_eq) for fracs in wide_f])
            rebuilt = [
                build_assemblage_labels(
                    pd_d["phase_ids"],
                    pd_d["phase_kinds"],
                    fracs,
                    active_threshold,
                    pd_d.get("phase_label", ""),
                    metals=self.metals,
                    phase_suffix=pd_d.get("phase_suffix", ""),
                    family_threshold=self._family_threshold,
                    family_values=wide_amounts[i],
                    family_threshold_inclusive=self.include_0p01_to_0p05_components,
                )
                for i, fracs in enumerate(wide_f)
            ]
            exact_l = np.array([item[0] for item in rebuilt], dtype=object)
            fam_l = np.array([item[1] for item in rebuilt], dtype=object)
            stale = np.mean(fam_l == "no feasible assemblage") > 0.9
            if not stale:
                if cache_source != fractions_cache or "mu_values" not in cached.files or "temperature_values" not in cached.files:
                    np.savez_compressed(
                        fractions_cache,
                        fractions=wide_f,
                        mu_values=mu_vals,
                        temperature_values=T_vals,
                    )
                if not self._run_plots:
                    return map_out
                region_ids, rlabels = assign_region_ids(fam_l.tolist())
                region_grid = region_ids.reshape(n_T, n_mu)
                region_details = _build_region_details(
                    region_ids,
                    exact_l,
                    wide_f,
                    pd_d["phase_ids"],
                    pd_d["phase_kinds"],
                    self.metals,
                    pd_d.get("phase_suffix", ""),
                    presence_values=wide_amounts,
                    presence_threshold=self._family_threshold,
                    presence_threshold_inclusive=self.include_0p01_to_0p05_components,
                )
                cell_texts = np.array(
                    [
                        format_exact_phase_fraction_line(
                            wide_f[i],
                            pd_d["phase_ids"],
                            pd_d["phase_kinds"],
                            self.metals,
                            phase_suffix=pd_d.get("phase_suffix", ""),
                        )
                        for i in range(n_states)
                    ],
                    dtype=object,
                ).reshape(n_T, n_mu)
                self._plot_muO_T(
                    mu_vals,
                    T_vals,
                    region_grid,
                    rlabels,
                    region_details,
                    cell_texts,
                    map_out,
                    sys_str,
                    M1,
                    M2,
                    boundary_lw,
                    region_label_fontsize,
                    pd_d["phase_ids"],
                    pd_d["phase_kinds"],
                )
                return map_out
        except Exception:
            pass

    if not self._run_calculations:
        raise RuntimeError(f"plot-only mode requires the single muO-T cache: {fractions_cache}")

    print(f"  muO-T map: {n_states} LP solves …")
    fam_l = np.empty(n_states, dtype=object)
    exact_l = np.empty(n_states, dtype=object)
    wide_f = np.zeros((n_states, len(pd_d["phase_ids"])))
    wide_amounts = np.zeros_like(wide_f)

    for iT, T in enumerate(T_vals):
        phase_G = pd_d["phase_H0"] + KB * float(T) * pd_d["phase_mix_shape"]
        for imu, mu_o in enumerate(mu_vals):
            row = iT * n_mu + imu
            grand = np.concatenate(
                [
                    pd_d["fixed_energy_formula"] - pd_d["phase_O"][: pd_d["n_fixed"]] * mu_o,
                    phase_G,
                ]
            )
            amounts, _, ok = solve_grand_lp(pd_d["A_eq"], b_eq, grand)
            if ok:
                total = float(np.nansum(amounts))
                fracs = amounts / total if total > 0 else amounts.copy()
                fracs[np.abs(fracs) < active_threshold] = 0.0
                el, fl = build_assemblage_labels(
                    pd_d["phase_ids"],
                    pd_d["phase_kinds"],
                    fracs,
                    active_threshold,
                    pd_d.get("phase_label", ""),
                    metals=self.metals,
                    phase_suffix=pd_d.get("phase_suffix", ""),
                    family_threshold=self._family_threshold,
                    family_values=amounts,
                    family_threshold_inclusive=self.include_0p01_to_0p05_components,
                )
                wide_f[row] = fracs
                wide_amounts[row] = amounts
            else:
                fl = el = "no feasible assemblage"
            fam_l[row] = fl
            exact_l[row] = el
        if (iT + 1) % max(1, n_T // 5) == 0 or iT == n_T - 1:
            print(f"    {iT + 1}/{n_T} T done")

    region_ids, rlabels = assign_region_ids(fam_l.tolist())
    region_grid = region_ids.reshape(n_T, n_mu)
    region_details = _build_region_details(
        region_ids,
        exact_l,
        wide_f,
        pd_d["phase_ids"],
        pd_d["phase_kinds"],
        self.metals,
        pd_d.get("phase_suffix", ""),
        presence_values=wide_amounts,
        presence_threshold=self._family_threshold,
        presence_threshold_inclusive=self.include_0p01_to_0p05_components,
    )
    cell_texts = np.array(
        [
            format_exact_phase_fraction_line(
                wide_f[i],
                pd_d["phase_ids"],
                pd_d["phase_kinds"],
                self.metals,
                phase_suffix=pd_d.get("phase_suffix", ""),
            )
            for i in range(n_states)
        ],
        dtype=object,
    ).reshape(n_T, n_mu)

    np.savez_compressed(
        fractions_cache,
        fractions=wide_f,
        mu_values=mu_vals,
        temperature_values=T_vals,
    )
    if not self._run_plots:
        return map_out
    self._plot_muO_T(
        mu_vals,
        T_vals,
        region_grid,
        rlabels,
        region_details,
        cell_texts,
        map_out,
        sys_str,
        M1,
        M2,
        boundary_lw,
        region_label_fontsize,
        pd_d["phase_ids"],
        pd_d["phase_kinds"],
    )
    return map_out


def run_muO_x_map(
    self,
    T: float,
    x_values: np.ndarray,
    mu_O_values: np.ndarray,
    skip_if_exists: bool = True,
) -> Path:
    """Compute and plot the μO vs composition map at fixed T.

    Args:
        self: SingleCompositionAnalyzer instance.
        T: Temperature in Kelvin.
        x_values: 1-D array of first-metal mole fractions.
        mu_O_values: 1-D array of oxygen chemical potentials (eV/O).
        skip_if_exists: Re-use on-disk caches when available.

    Returns:
        Path to the muO_x_phase_map/T<int(T)> output directory.
    """
    from framework._run import _build_region_details

    from .thermodynamics import (
        KB,
        assign_region_ids,
        build_assemblage_labels,
        format_exact_phase_fraction_line,
        plot_region_map,
        solve_grand_lp,
    )

    self._ensure_path()
    pd_d = self.pd_data
    x_vals = np.asarray(x_values)
    mu_vals = np.asarray(mu_O_values)
    M1 = self.metals[0]
    M2 = self.metals[1] if len(self.metals) > 1 else ""
    M3 = self.metals[2] if len(self.metals) > 2 else None
    out = self.out_dir
    map_out = out / "muO_x_phase_map" / f"T{int(T)}"
    cache_dir = out / "cache" / "muO_x_phase_map" / f"T{int(T)}"
    cache_dir.mkdir(parents=True, exist_ok=True)
    fractions_cache = cache_dir / "phase_fractions.npz"
    legacy_cache = out / f"muO_x_map_T{int(T)}_phase_fractions.npz"
    cache_source = fractions_cache if fractions_cache.exists() else legacy_cache
    n_expected = len(x_vals) * len(mu_vals)

    if skip_if_exists and cache_source.exists():
        try:
            cached = np.load(cache_source)
            wide_f = cached["fractions"]
            if wide_f.shape != (n_expected, len(pd_d["phase_ids"])):
                raise ValueError("cached phase-fraction matrix has the wrong shape")
            if "mu_values" in cached.files and not np.array_equal(cached["mu_values"], mu_vals):
                raise ValueError("cached muO grid does not match")
            if "x_values" in cached.files and not np.array_equal(cached["x_values"], x_vals):
                raise ValueError("cached composition grid does not match")
            if "temperature" in cached.files and not np.isclose(float(cached["temperature"]), float(T)):
                raise ValueError("cached temperature does not match")
            if (
                cache_source != fractions_cache
                or "mu_values" not in cached.files
                or "x_values" not in cached.files
                or "temperature" not in cached.files
            ):
                np.savez_compressed(
                    fractions_cache,
                    fractions=wide_f,
                    mu_values=mu_vals,
                    x_values=x_vals,
                    temperature=float(T),
                )
            if not self._run_plots:
                return map_out
            rebuilt = []
            wide_amounts = np.zeros_like(wide_f)
            for row_index, fracs in enumerate(wide_f):
                ix = row_index // len(mu_vals)
                x = float(x_vals[ix])
                n_m = len(self.metals)
                comp = np.array([x, 1.0 - x]) if n_m == 2 else np.array([x] + [(1.0 - x) / (n_m - 1)] * (n_m - 1))
                b_eq = np.concatenate([comp, [self.phase_element_stoichiometry]])
                wide_amounts[row_index] = _amounts_from_fractions(self, fracs, b_eq)
                rebuilt.append(
                    build_assemblage_labels(
                        pd_d["phase_ids"],
                        pd_d["phase_kinds"],
                        fracs,
                        1e-9,
                        pd_d.get("phase_label", ""),
                        metals=self.metals,
                        phase_suffix=pd_d.get("phase_suffix", ""),
                        family_threshold=self._family_threshold,
                        family_values=wide_amounts[row_index],
                        family_threshold_inclusive=self.include_0p01_to_0p05_components,
                    )
                )
            exact_l = np.array([item[0] for item in rebuilt], dtype=object)
            fam_l = np.array([item[1] for item in rebuilt], dtype=object)
            region_ids, rlabels = assign_region_ids(fam_l.tolist())
            region_grid = region_ids.reshape(len(x_vals), len(mu_vals))
            region_details = _build_region_details(
                region_ids,
                exact_l,
                wide_f,
                pd_d["phase_ids"],
                pd_d["phase_kinds"],
                self.metals,
                pd_d.get("phase_suffix", ""),
                presence_values=wide_amounts,
                presence_threshold=self._family_threshold,
                presence_threshold_inclusive=self.include_0p01_to_0p05_components,
            )
            cell_texts = np.array(
                [
                    format_exact_phase_fraction_line(
                        row,
                        pd_d["phase_ids"],
                        pd_d["phase_kinds"],
                        self.metals,
                        phase_suffix=pd_d.get("phase_suffix", ""),
                    )
                    for row in wide_f
                ],
                dtype=object,
            ).reshape(len(x_vals), len(mu_vals))
            plot_region_map(
                mu_vals,
                x_vals,
                region_grid,
                rlabels,
                T,
                M1,
                M2,
                map_out,
                M3=M3,
                region_label_mode=self.region_label_mode,
                region_label_fontsize=self.region_label_fontsize,
                region_details=region_details,
                cell_text_grid=cell_texts,
                system_label=self.sys_cfg.phase_label,
            )
            return map_out
        except Exception:
            pass

    if not self._run_calculations:
        raise RuntimeError(f"plot-only mode requires the single muO-x cache: {fractions_cache}")

    print(f"  muO-x map T={int(T)}K: {n_expected} LP solves …")
    phase_G = pd_d["phase_H0"] + KB * float(T) * pd_d["phase_mix_shape"]
    wide_f = np.zeros((n_expected, len(pd_d["phase_ids"])))
    exact_l = np.empty(n_expected, dtype=object)
    fam_l = np.empty(n_expected, dtype=object)
    wide_amounts = np.zeros_like(wide_f)
    row_index = 0
    for x in x_vals:
        n_m = len(self.metals)
        comp = (
            np.array([float(x), 1.0 - float(x)])
            if n_m == 2
            else np.array([float(x)] + [(1.0 - float(x)) / (n_m - 1)] * (n_m - 1))
        )
        b_eq = np.concatenate([comp, [self.phase_element_stoichiometry]])
        for mu_o in mu_vals:
            grand = np.concatenate(
                [
                    pd_d["fixed_energy_formula"] - pd_d["phase_O"][: pd_d["n_fixed"]] * mu_o,
                    phase_G,
                ]
            )
            amounts, _, ok = solve_grand_lp(pd_d["A_eq"], b_eq, grand)
            if ok:
                total = float(np.nansum(amounts))
                fracs = amounts / total if total > 0 else amounts.copy()
                fracs[np.abs(fracs) < 1e-9] = 0.0
                el, fl = build_assemblage_labels(
                    pd_d["phase_ids"],
                    pd_d["phase_kinds"],
                    fracs,
                    1e-9,
                    pd_d.get("phase_label", ""),
                    metals=self.metals,
                    phase_suffix=pd_d.get("phase_suffix", ""),
                    family_threshold=self._family_threshold,
                    family_values=amounts,
                    family_threshold_inclusive=self.include_0p01_to_0p05_components,
                )
            else:
                el = fl = "no feasible assemblage"
            exact_l[row_index] = el
            fam_l[row_index] = fl
            if ok:
                wide_f[row_index] = fracs
                wide_amounts[row_index] = amounts
            row_index += 1
    np.savez_compressed(
        fractions_cache,
        fractions=wide_f,
        mu_values=mu_vals,
        x_values=x_vals,
        temperature=float(T),
    )
    if not self._run_plots:
        return map_out
    region_ids, rlabels = assign_region_ids(fam_l.tolist())
    region_grid = region_ids.reshape(len(x_vals), len(mu_vals))
    region_details = _build_region_details(
        region_ids,
        exact_l,
        wide_f,
        pd_d["phase_ids"],
        pd_d["phase_kinds"],
        self.metals,
        pd_d.get("phase_suffix", ""),
        presence_values=wide_amounts,
        presence_threshold=self._family_threshold,
        presence_threshold_inclusive=self.include_0p01_to_0p05_components,
    )
    cell_texts = np.array(
        [
            format_exact_phase_fraction_line(
                row, pd_d["phase_ids"], pd_d["phase_kinds"], self.metals, phase_suffix=pd_d.get("phase_suffix", "")
            )
            for row in wide_f
        ],
        dtype=object,
    ).reshape(len(x_vals), len(mu_vals))
    plot_region_map(
        mu_vals,
        x_vals,
        region_grid,
        rlabels,
        T,
        M1,
        M2,
        map_out,
        M3=M3,
        region_label_mode=self.region_label_mode,
        region_label_fontsize=self.region_label_fontsize,
        region_details=region_details,
        cell_text_grid=cell_texts,
        system_label=self.sys_cfg.phase_label,
    )
    return map_out
