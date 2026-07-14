"""SingleCompositionAnalyzer — all graphs for one specific composition.

Usage in run_simple.py or a notebook::

    from framework import SingleCompositionAnalyzer
    import numpy as np

    ana = SingleCompositionAnalyzer(
        system_dir = Path("/path/to/CrHfB"),
        tables_dir = Path("tables"),
        figures_dir = Path("figures"),
        metals     = ["Cr", "Hf"],          # auto-detected if None
        composition = 0.3,                  # binary: x_M1 scalar
        #composition = [0.3, 0.5, 0.2],    # ternary: [x_M1, x_M2, x_M3]
    )
    ana.run(
        T_values   = np.arange(200, 2001, 10),
        mu_O_values = np.arange(-10, -4, 0.1),
        scan_T     = [700, 1000, 1300],
    )
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from .utils import prepare_system_tables_dir, system_key


class SingleCompositionAnalyzer:
    """Run all standard analyses for one specific composition.

    Parameters
    ----------
    system_dir   : raw GRACE calculation folder (blade/ + grace/ subdirs)
    tables_dir   : root containing one CSV folder per modeled system
    figures_dir  : where figures are saved
    metals       : element list e.g. ['Cr','Hf']. Auto-detected if None.
    composition  : scalar x_M1 for binary; list [x_M1,x_M2,...] for n-ary.
                   Fractions must sum to 1.
    rk_order     : RK polynomial order for table preparation
    y_step       : flexible-phase composition grid step
    mixed_phase_subdir  : subfolder name for SQS/GRACE outputs
    fixed_phases_subdir : subfolder name for reference phase calculations
    """

    def __init__(
        self,
        system_dir: str | Path | None = None,
        tables_dir: str | Path = "tables",
        figures_dir: str | Path = "figures",
        metals: list[str] | None = None,
        composition: float | list[float] = 0.5,
        rk_order: int = 3,
        y_step: float = 0.01,
        mixed_phase_subdir: str = "blade",
        fixed_phases_subdir: str = "ORB",
        phase_element: str | None = None,
        phase_element_stoichiometry: float = 0.0,
        include_0p01_to_0p05_components: bool = True,
        region_label_mode: str = "id",
        region_label_fontsize: int = 12,
    ):
        self.system_dir = Path(system_dir) if system_dir else None
        self.tables_root = Path(tables_dir)
        self.figures_dir = Path(figures_dir)
        self.rk_order = rk_order
        self.y_step = y_step
        self.mixed_phase_subdir = mixed_phase_subdir
        self.fixed_phases_subdir = fixed_phases_subdir
        self.phase_element = str(phase_element).strip() if phase_element else None
        self.phase_element_stoichiometry = float(phase_element_stoichiometry) if self.phase_element else 0.0
        self.include_0p01_to_0p05_components = bool(include_0p01_to_0p05_components)
        mode = str(region_label_mode).strip().lower()
        if mode in {"id", "ids", "number", "numbers"}:
            self.region_label_mode = "id"
        elif mode in {"phase", "phases", "full"}:
            self.region_label_mode = "phases"
        else:
            raise ValueError("region_label_mode must be 'id' or 'phases'")
        self.region_label_fontsize = int(region_label_fontsize)
        if self.region_label_fontsize <= 0:
            raise ValueError("region_label_fontsize must be positive")
        if self.phase_element and self.phase_element_stoichiometry <= 0:
            raise ValueError("phase_element_stoichiometry must be positive when phase_element is set")

        self._ensure_path()

        # Resolve metals
        if metals is None:
            metals = self._detect_metals()
        self.metals = metals
        self.n_metals = len(metals)
        self.tables_dir = prepare_system_tables_dir(self.tables_root, self.metals, self.phase_element)

        # Normalise composition
        if np.isscalar(composition):
            x = float(composition)
            if self.n_metals == 2:
                self._comp = np.array([x, 1.0 - x])
            else:
                rest = (1.0 - x) / (self.n_metals - 1)
                self._comp = np.array([x] + [rest] * (self.n_metals - 1))
        else:
            self._comp = np.asarray(composition, dtype=float)
            if abs(self._comp.sum() - 1.0) > 1e-6:
                raise ValueError("Composition fractions must sum to 1.0")

        self._x_label = "  ".join(f"x_{m}={self._comp[i]:.3g}" for i, m in enumerate(self.metals))
        self._x_str = "_".join(f"{m}{self._comp[i]:.2f}".replace(".", "p") for i, m in enumerate(self.metals))

        self._sys_cfg = None
        self._pd_data = None
        self._run_calculations = True
        self._run_plots = True

    # ---------------------------------------------------------------- setup

    def _ensure_path(self):
        root = Path(__file__).parent
        for p in [str(root), str(root / "python")]:
            if p not in sys.path:
                sys.path.insert(0, p)

    def _detect_metals(self) -> list[str]:
        self._ensure_path()
        from .system_config import SystemConfig

        if self.system_dir and (self.system_dir / self.mixed_phase_subdir).is_dir():
            return SystemConfig._detect_metals_from_dirs(
                self.system_dir / self.mixed_phase_subdir, {self.phase_element} if self.phase_element else set()
            )
        raise ValueError("system_dir with mixed_phase_subdir required for auto-detection")

    def prepare(self) -> "SingleCompositionAnalyzer":
        """Prepare tables (skips if already exist). Returns self."""
        self._ensure_path()
        from .system_config import SystemConfig, prepare_tables

        self._sys_cfg = SystemConfig.resolve(
            self.system_dir,
            self.tables_dir,
            self.metals,
            blade_subdir=self.mixed_phase_subdir,
            fixed_phases_subdir=self.fixed_phases_subdir,
            phase_element=self.phase_element,
            phase_element_stoichiometry=self.phase_element_stoichiometry,
        )
        if not self._sys_cfg.tables_ready():
            print(f"  [{self._sys_cfg.tag}] preparing tables …")
            prepare_tables(self.system_dir, self._sys_cfg, self.rk_order)
        else:
            print(f"  [{self._sys_cfg.tag}] tables exist — skipping")
        return self

    def load(self) -> dict:
        """Load equilibrium phase data from the prepared tables."""
        self._ensure_path()
        from .thermodynamics import load_phase_data

        if self._sys_cfg is None:
            self.prepare()
        pd_d = load_phase_data(
            self._sys_cfg.phase_table_file,
            self._sys_cfg.phase_grid_file,
            self._sys_cfg.interaction_coeff_file,
            self._sys_cfg.metals,
            phase_label=self._sys_cfg.phase_label,
            y_step=self.y_step,
            phase_element=self._sys_cfg.phase_element,
            phase_element_stoichiometry=self._sys_cfg.phase_element_stoichiometry,
        )
        self._pd_data = pd_d
        return self._pd_data

    @property
    def pd_data(self) -> dict:
        if self._pd_data is None:
            self.load()
        return self._pd_data

    @property
    def sys_cfg(self):
        if self._sys_cfg is None:
            self.prepare()
        return self._sys_cfg

    @property
    def out_dir(self) -> Path:
        sys_name = system_key(self.metals, self.phase_element)
        d = self.figures_dir / sys_name / "single_composition" / self._x_str
        d.mkdir(parents=True, exist_ok=True)
        return d

    # ---------------------------------------------------------------- b_eq

    def _conservation_rhs(self) -> np.ndarray:
        """Build the elemental conservation vector from explicit settings."""
        return np.concatenate([self._comp, [self.phase_element_stoichiometry]])

    @property
    def _family_threshold(self) -> float:
        return 0.01 if self.include_0p01_to_0p05_components else 0.05

    def _amounts_from_fractions(self, fracs: np.ndarray, b_eq: np.ndarray) -> np.ndarray:
        """Recover the LP amount scale from normalized cached fractions."""
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
        """1D oxidation scan: phase fractions vs μO at each T. Returns out_dir."""
        import matplotlib
        import matplotlib.pyplot as plt
        import pandas as pd
        from thermodynamics import (
            KB,
            _phase_comp_label,
            _short_label,
            solve_grand_lp,
        )

        matplotlib.use("Agg")

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
                                (
                                    legacy[str(pid)].to_numpy(dtype=float)
                                    if str(pid) in legacy.columns
                                    else np.zeros(len(mu_vals))
                                )
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
            for idx, (pid, kind) in enumerate(zip(pd_d["phase_ids"], pd_d["phase_kinds"])):
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
        """μO vs T assemblage map at this fixed composition."""
        from framework._run import _build_region_details
        from thermodynamics import (
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
                wide_amounts = np.vstack([self._amounts_from_fractions(fracs, b_eq) for fracs in wide_f])
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
                    if (
                        cache_source != fractions_cache
                        or "mu_values" not in cached.files
                        or "temperature_values" not in cached.files
                    ):
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
                print(f"    {iT+1}/{n_T} T done")

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

    def _plot_muO_T(
        self,
        mu_vals,
        T_vals,
        region_grid,
        rlabels,
        region_details,
        cell_texts,
        out_dir,
        sys_str,
        M1,
        M2,
        boundary_lw,
        region_label_fontsize,
        phase_ids,
        phase_kinds,
    ):
        import hashlib as _hl

        import matplotlib
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt
        from thermodynamics import (
            add_region_annotation,
            format_phase_detail_line,
            grid_edges,
            region_annotation_text,
            separate_region_annotations,
            write_region_map_html,
        )

        matplotlib.use("Agg")
        out_dir.mkdir(parents=True, exist_ok=True)
        _cmap = plt.get_cmap("tab20")
        # Deterministic hash color — same assemblage always same color
        colors = {r + 1: _cmap(int(_hl.md5(lbl.encode()).hexdigest(), 16) % 20 / 20) for r, lbl in enumerate(rlabels)}
        mu_edges = grid_edges(mu_vals)
        T_edges = grid_edges(T_vals)

        if self.region_label_mode == "id":
            fig, (ax, ax_leg) = plt.subplots(
                1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05}
            )
            ax_leg.axis("off")
        else:
            fig, ax = plt.subplots(figsize=(13, 8))
            ax_leg = None

        for iT in range(len(T_vals)):
            for imu in range(len(mu_vals)):
                rid = region_grid[iT, imu]
                ax.add_patch(
                    plt.Rectangle(
                        (mu_edges[imu], T_edges[iT]),
                        mu_edges[imu + 1] - mu_edges[imu],
                        T_edges[iT + 1] - T_edges[iT],
                        color=colors.get(rid, (0.85, 0.85, 0.85, 1.0)),
                        linewidth=0,
                    )
                )
        for iT in range(len(T_vals)):
            for imu in range(len(mu_vals) - 1):
                if region_grid[iT, imu] != region_grid[iT, imu + 1]:
                    ax.plot([mu_edges[imu + 1]] * 2, [T_edges[iT], T_edges[iT + 1]], "k-", lw=boundary_lw)
        for iT in range(len(T_vals) - 1):
            for imu in range(len(mu_vals)):
                if region_grid[iT, imu] != region_grid[iT + 1, imu]:
                    ax.plot([mu_edges[imu], mu_edges[imu + 1]], [T_edges[iT + 1]] * 2, "k-", lw=boundary_lw)
        annotation_artists = []
        for rid, lbl in enumerate(rlabels, start=1):
            mask = region_grid == rid
            if not np.any(mask):
                continue
            rows, cols = np.where(mask)
            r = max(0, min(len(T_vals) - 1, int(np.round(np.median(rows)))))
            c = max(0, min(len(mu_vals) - 1, int(np.round(np.median(cols)))))
            annotation = region_annotation_text(rid, lbl, region_details, self.region_label_mode)
            region_color = colors[rid]

            lightness = 0.5   # 0.0 = unchanged, larger = lighter
            annotation_color = (
                region_color[0] + (1.0 - region_color[0]) * lightness,
                region_color[1] + (1.0 - region_color[1]) * lightness,
                region_color[2] + (1.0 - region_color[2]) * lightness,
                region_color[3],
            )

            annotation_artists.append(
                add_region_annotation(
                    ax,
                    mu_vals[c],
                    T_vals[r],
                    annotation,
                    9,
                    annotation_color,
                )
            )
        ax.set_xlim(mu_vals[0], mu_vals[-1])
        ax.set_ylim(T_vals[0], T_vals[-1])
        ax.tick_params(axis="x", labelsize=16)
        ax.tick_params(axis="y", labelsize=16)
        ax.set_xlabel(r"$\mu_O$ (eV per O atom)", fontsize=20)
        ax.set_ylabel("T (K)", fontsize=20)
        ax.set_title(f"{sys_str} | {self._x_label}", fontsize=24)
        separate_region_annotations(ax, annotation_artists)

        def _fmt(detail):
            if not detail:
                return ""
            pr = detail.get("phase_ranges") or []
            return format_phase_detail_line(pr) if pr else (detail.get("exact_label") or "")

        handles = []
        for rid, lbl in enumerate(rlabels, start=1):
            if not np.any(region_grid == rid):
                continue
            detail = _fmt(region_details.get(rid) if region_details else None)
            handles.append(mpatches.Patch(color=colors[rid], label=f"{rid}: {detail if detail else lbl}"))
        leg_fs = max(5, min(8, int(110 / max(len(handles), 1))))
        if ax_leg is not None:
            ax_leg.legend(
                handles=handles,
                loc="upper left",
                bbox_to_anchor=(0.0, 1.0),
                fontsize=leg_fs,
                framealpha=0.9,
                handlelength=1.2,
                borderaxespad=0,
                title="Region key",
                title_fontsize=leg_fs + 5,
            )

        png = out_dir / "muO_T_map.png"
        fig.savefig(png, dpi=150, bbox_inches="tight")
        try:
            write_region_map_html(
                mu_vals,
                T_vals,
                region_grid,
                rlabels,
                out_dir / "muO_T_map.html",
                title=f"{sys_str} | {self._x_label}",
                x_label="mu_O (eV/O)",
                y_label="T (K)",
                region_details=region_details,
                cell_text_grid=cell_texts,
                region_label_mode=self.region_label_mode,
                region_label_fontsize=region_label_fontsize,
            )
        except Exception as e:
            print(f"  [html skipped] {e}")
        plt.close(fig)
        print(f"  Figure → {png}")

    def run_muO_x_map(
        self,
        T: float,
        x_values: np.ndarray,
        mu_O_values: np.ndarray,
        skip_if_exists: bool = True,
    ) -> Path:
        """μO vs x map at fixed T showing where this composition lies."""
        from framework._run import _build_region_details
        from thermodynamics import (
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
                    wide_amounts[row_index] = self._amounts_from_fractions(fracs, b_eq)
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

    # ---------------------------------------------------------------- run all

    def run(
        self,
        T_values: np.ndarray | None = None,
        mu_O_values: np.ndarray | None = None,
        scan_T: list[float] | None = None,
        x_values: np.ndarray | None = None,
        skip_if_exists: bool = True,
        run_calculations: bool = True,
        run_plots: bool = True,
    ) -> Path:
        """Run all analyses for this composition. Returns output directory.

        When both stage flags are true, all numerical caches are completed
        before a second, cache-only plotting pass starts.
        """
        import traceback

        T_vals = np.asarray(T_values if T_values is not None else np.arange(200, 2001, 10))
        mu_vals = np.asarray(mu_O_values if mu_O_values is not None else np.arange(-10, -4, 0.1))
        s_T = scan_T if scan_T is not None else [int(T_vals[len(T_vals) // 2])]
        x_vals = np.asarray(x_values if x_values is not None else np.arange(0.0, 1.01, 0.01))

        self.prepare()
        self.load()

        def execute_pass(allow_calculations: bool, allow_plots: bool, use_cache: bool) -> None:
            self._run_calculations = allow_calculations
            self._run_plots = allow_plots
            failures = []
            for label, fn in [
                ("scan", lambda: self.run_scan(s_T, mu_vals, use_cache)),
                ("muO-T map", lambda: self.run_muO_T_map(T_vals, mu_vals, use_cache)),
                ("muO-x map", lambda: [self.run_muO_x_map(T, x_vals, mu_vals, use_cache) for T in s_T]),
            ]:
                try:
                    fn()
                except Exception as e:
                    print(f"  [{label}] FAILED: {e}")
                    traceback.print_exc()
                    failures.append(f"{label}: {e}")
            if failures:
                raise RuntimeError("; ".join(failures))

        if run_calculations and run_plots:
            print("\n=== Single-composition calculation pass ===")
            execute_pass(True, False, skip_if_exists)
            print("\n=== Single-composition plotting pass (cached data only) ===")
            execute_pass(False, True, True)
        elif run_calculations or run_plots:
            execute_pass(run_calculations, run_plots, True if run_plots else skip_if_exists)
        else:
            print("Single-composition oxidation calculations and plots are both disabled.")

        if run_plots:
            print(f"\n  All figures → {self.out_dir}")
        return self.out_dir
