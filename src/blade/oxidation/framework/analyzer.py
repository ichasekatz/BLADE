"""SystemAnalyzer — runs all analyses for one system."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from .config import Config
from .utils import (
    csv_has_rows,
    fmt_frac,
    make_animation,
    normalize_region_labels,
    phase_short,
    system_key,
)


class SystemAnalyzer:
    """Run all configured analyses for one system.

    Usage::

        cfg = Config(...)
        analyzer = SystemAnalyzer(system_dir, config=cfg)
        analyzer.run()
    """

    def __init__(self, system_dir: Path, config: Config):
        self.system_dir = Path(system_dir)
        self.config = config

        # Lazy-loaded after tables are prepared
        self._sys_cfg = None
        self._pd_data = None

    # ---------------------------------------------------------------- setup

    def _ensure_imports(self):
        root = Path(__file__).parent
        for p in [str(root), str(root / "python")]:
            if p not in sys.path:
                sys.path.insert(0, p)

    def prepare(self, metals: list[str]) -> bool:
        """Load SystemConfig and prepare tables. Returns False on failure."""
        self._ensure_imports()
        from .system_config import SystemConfig, prepare_tables

        try:
            self._sys_cfg = SystemConfig.resolve(
                self.system_dir,
                self.config.tables_dir,
                metals,
                blade_subdir=self.config.blade_subdir,
                fixed_phases_subdir=self.config.fixed_phases_subdir,
                phase_element=self.config.phase_element,
                phase_element_stoichiometry=self.config.phase_element_stoichiometry,
            )
            cfg = self.config
            if not self._sys_cfg.tables_ready() or not cfg.skip_if_tables_exist:
                prepare_tables(self.system_dir, self._sys_cfg, cfg.rk_order)
            return True
        except Exception as e:
            print(f"  table prep failed: {e}")
            return False

    def load_phases(self):
        """Load pd_data from tables. Must call prepare() first."""
        self._ensure_imports()
        from ._run import load_phases as _load_phases

        return _load_phases(self._sys_cfg)

    def run_standard(self) -> None:
        """Run standard analyses via MapsRunner (no direct run.py calls)."""
        from .maps import MapsRunner

        any_standard = any(
            [
                self.config.run_oxidation_scan,
                self.config.run_muO_x_map,
                self.config.run_muO_T_map,
                self.config.run_fixed_phase_map,
                self.config.run_ternary_3d_map,
            ]
        )
        if any_standard:
            runner = MapsRunner(self.system_dir, self._sys_cfg.metals, self.config)
            runner.run(self._sys_cfg, self._pd_data)

    def run_composition_slices(self) -> None:
        if not self.config.run_composition_slice_maps:
            return
        if self._pd_data is None:
            self._pd_data = self.load_phases()
        self._run_composition_slice_maps(self._sys_cfg, self._pd_data)
        # muO-T slices: ~62k LP per spec — only run when explicitly enabled
        if getattr(self.config, "run_composition_slice_muT_maps", False):
            self._run_composition_slice_muT_maps(self._sys_cfg, self._pd_data)

    def run_onset_auc(self) -> None:
        if not self.config.run_onset_auc_diagrams:
            return
        if self._pd_data is None:
            self._pd_data = self.load_phases()
        self._run_onset_auc_diagrams(self._sys_cfg, self._pd_data)

    def run(self) -> None:
        """Run all enabled analyses in sequence."""
        if self._pd_data is None:
            self._pd_data = self.load_phases()
        import traceback

        failures = []
        for label, fn in (
            ("standard analyses", self.run_standard),
            ("composition slices", self.run_composition_slices),
            ("onset AUC", self.run_onset_auc),
        ):
            try:
                fn()
            except Exception as e:
                print(f"  [{label}] FAILED: {e}")
                traceback.print_exc()
                failures.append(f"{label}: {e}")
        if failures:
            raise RuntimeError("; ".join(failures))

    @staticmethod
    def _cache_columns_match(csv_path, expected_columns) -> bool:
        """Return true only when cached coordinates exactly match the requested grid."""
        import pandas as pd

        try:
            columns = list(expected_columns)
            df = pd.read_csv(csv_path, usecols=columns)
            if len(df) != len(next(iter(expected_columns.values()))):
                return False
            return all(
                np.allclose(
                    df[column].to_numpy(dtype=float),
                    np.asarray(expected, dtype=float),
                    atol=1e-10,
                    rtol=0.0,
                    equal_nan=True,
                )
                for column, expected in expected_columns.items()
            )
        except (OSError, KeyError, TypeError, ValueError):
            return False

    def _slice_cache_matches(self, csv_path, T, spec, metals, x_values, mu_values) -> bool:
        if not (
            self.config.skip_if_analysis_exists
            and csv_has_rows(csv_path, len(x_values) * len(mu_values))
            and self._component_cache_matches(csv_path)
        ):
            return False
        compositions = np.asarray(
            [self._composition_from_slice(float(x), len(metals), spec["axis"], spec["remainder"]) for x in x_values],
            dtype=float,
        )
        expected = {
            "T_K": np.full(len(x_values) * len(mu_values), float(T)),
            "muO_eV_per_O": np.tile(mu_values, len(x_values)),
            f"x_{metals[spec['axis']]}_axis": np.repeat(x_values, len(mu_values)),
        }
        expected.update({f"x_{metal}": np.repeat(compositions[:, i], len(mu_values)) for i, metal in enumerate(metals)})
        return self._cache_columns_match(csv_path, expected)

    def _slice_muT_cache_matches(self, csv_path, comp, axis_m, x, metals, T_values, mu_values) -> bool:
        if not (
            self.config.skip_if_analysis_exists
            and csv_has_rows(csv_path, len(T_values) * len(mu_values))
            and self._component_cache_matches(csv_path)
        ):
            return False
        n_states = len(T_values) * len(mu_values)
        expected = {
            "T_K": np.repeat(T_values, len(mu_values)),
            "muO_eV_per_O": np.tile(mu_values, len(T_values)),
            f"x_{axis_m}_axis": np.full(n_states, float(x)),
        }
        expected.update({f"x_{metal}": np.full(n_states, float(comp[i])) for i, metal in enumerate(metals)})
        if not self._cache_columns_match(csv_path, expected):
            return False
        try:
            import pandas as pd

            columns = set(pd.read_csv(csv_path, nrows=0).columns)
            return {
                "phase_fraction",
                "parent_phase_fraction",
                "oxide_phase_fraction",
                "assemblage_exact",
                "phase_fraction_summary",
            }.issubset(columns)
        except (OSError, ValueError):
            return False

    def _onset_cache_matches(self, csv_path, comp_grid, comp_step, mu_values, T, metals) -> bool:
        if not (self.config.skip_if_analysis_exists and csv_has_rows(csv_path, len(comp_grid))):
            return False
        compositions = np.asarray(comp_grid, dtype=float)
        n_rows = len(compositions)
        mu_step = float(mu_values[1] - mu_values[0]) if len(mu_values) > 1 else np.nan
        expected = {
            "T_K": np.full(n_rows, float(T)),
            "comp_step": np.full(n_rows, float(comp_step)),
            "muO_min": np.full(n_rows, float(mu_values[0])),
            "muO_max": np.full(n_rows, float(mu_values[-1])),
            "muO_step": np.full(n_rows, mu_step),
        }
        expected.update({f"x_{metal}": compositions[:, i] for i, metal in enumerate(metals)})
        if not self._cache_columns_match(csv_path, expected):
            return False
        try:
            import pandas as pd

            columns = set(pd.read_csv(csv_path, nrows=0).columns)
            return {
                "feasible_muO_max_eV",
                "first_infeasible_muO_eV",
                "feasible_muO_points",
                "parent_phase_fraction_auc_eV",
                "oxide_phase_fraction_auc_eV",
            }.issubset(columns)
        except (OSError, ValueError):
            return False

    # ---------------------------------------------------------------- internals

    def _composition_slice_specs(self, metals):
        from itertools import product

        n = len(metals)
        metal_lookup = {metal.lower(): i for i, metal in enumerate(metals)}
        priority = [str(element).strip() for element in self.config.slice_axis_priority if str(element).strip()]
        priority_axes = [metal_lookup[element.lower()] for element in priority if element.lower() in metal_lookup]
        if priority and not priority_axes:
            raise ValueError(f"none of the slice_axis_priority elements {priority} " f"are present in {metals}")
        selector = priority_axes[0] if priority else self.config.slice_axis
        if selector is None:
            axes = list(range(n))
        elif isinstance(selector, str):
            axes = [i for i, metal in enumerate(metals) if metal.lower() == selector.strip().lower()]
            if not axes:
                raise ValueError(f"slice_axis {selector!r} is not one of {metals}")
        else:
            axis_index = int(selector)
            if axis_index < 0 or axis_index >= n:
                raise ValueError(f"slice_axis index {axis_index} is outside 0..{n - 1}")
            axes = [axis_index]
        if n == 2:
            return [
                {
                    "axis": axis,
                    "remainder": [1.0],
                    "name": f"axis_{metals[axis]}",
                    "label": f"x_{metals[axis]}",
                }
                for axis in axes
            ]
        specs = []
        ratio_values = sorted(
            {
                round(float(value), 12)
                for value in self.config.slice_remainder_ratios
                if -1e-12 <= float(value) <= 1.0 + 1e-12
            }
        )
        if not ratio_values:
            raise ValueError("slice_remainder_ratios must contain values from 0 to 1")
        for axis in axes:
            others = [i for i in range(n) if i != axis]
            if n == 3:
                for ratio in ratio_values:
                    rem = [float(ratio), 1.0 - float(ratio)]
                    name = (
                        f"axis_{metals[axis]}_"
                        f"{metals[others[0]]}{fmt_frac(rem[0])}_"
                        f"{metals[others[1]]}{fmt_frac(rem[1])}"
                    )
                    label = (
                        f"x_{metals[axis]}, remainder "
                        f"{rem[0]:.2f}{metals[others[0]]}/"
                        f"{rem[1]:.2f}{metals[others[1]]}"
                    )
                    specs.append({"axis": axis, "remainder": rem, "name": name, "label": label})
            else:
                remainder_grid = [
                    list(values) for values in product(ratio_values, repeat=n - 1) if abs(sum(values) - 1.0) < 1e-9
                ]
                if not remainder_grid:
                    remainder_grid = [[1.0 / (n - 1)] * (n - 1)]
                for rem in remainder_grid:
                    suffix = "_".join(f"{metals[index]}{fmt_frac(value)}" for index, value in zip(others, rem))
                    split = "/".join(f"{value:.2f}{metals[index]}" for index, value in zip(others, rem))
                    specs.append(
                        {
                            "axis": axis,
                            "remainder": rem,
                            "name": f"axis_{metals[axis]}_{suffix}",
                            "label": f"x_{metals[axis]}, remainder {split}",
                        }
                    )
        return specs

    @staticmethod
    def _composition_from_slice(x_axis: float, n_metals: int, axis: int, remainder_weights):
        comp = np.zeros(n_metals, dtype=float)
        comp[axis] = x_axis
        rest = max(0.0, 1.0 - x_axis)
        others = [i for i in range(n_metals) if i != axis]
        weights = np.asarray(remainder_weights, dtype=float)
        weights = weights / weights.sum() if weights.sum() > 0 else np.ones(len(others)) / len(others)
        for i, w in zip(others, weights):
            comp[i] = rest * w
        return comp

    # A metal is absent (→ new region) if its fraction stays below this.
    _COMPONENT_THRESHOLD = 0.01

    def _component_cache_matches(self, csv_path) -> bool:
        try:
            import pandas as pd

            values = pd.read_csv(csv_path, usecols=["component_presence_threshold"])[
                "component_presence_threshold"
            ].to_numpy(dtype=float)
            expected = 0.01 if self.config.include_0p01_to_0p05_components else 0.05
            return len(values) > 0 and np.allclose(values, expected)
        except Exception:
            return False

    def _coarse_family_label(self, phase_ids, phase_kinds, values, phase_y_nd, metals):
        """Region key from compounds present above the configured raw amount."""
        from thermodynamics import _phase_component_signature

        threshold = 0.01 if self.config.include_0p01_to_0p05_components else 0.05
        active = (
            values >= threshold - 1e-12 if self.config.include_0p01_to_0p05_components else values > threshold + 1e-12
        )
        fixed = sorted(
            {phase_short(pid) for pid, kind in zip(phase_ids[active], phase_kinds[active]) if kind != "phase"}
        )
        phase_signatures = sorted(
            {
                _phase_component_signature(pid, metals)
                for pid, kind in zip(phase_ids[active], phase_kinds[active])
                if kind == "phase"
            }
        )
        parts = fixed + phase_signatures
        return " + ".join(parts) if parts else "no feasible assemblage"

    @staticmethod
    def _merge_same_component_adjacent_regions(labels: list, nrows: int, ncols: int) -> list:
        """Flood-fill adjacent cells whose phase composition ranges touch or overlap.

        Two cells connect when:
          - same fixed phases
          - same metals present in the flexible phase
          - every metal's composition range overlaps or touches
            (e.g. 0.00-0.05 connects 0.05-0.10 at 0.05)

        BFS expands via direct neighbours so transitivity holds: A→B→C merge
        even if A and C don't directly touch.
        """
        import re
        from collections import Counter, deque

        _PHASE_RE = re.compile(r"\(([^)]+)\)(?:[A-Z][A-Za-z0-9.]*)?")
        _COMPONENT_RE = re.compile(r"(?:(\d*\.?\d+)(?:-(\d*\.?\d+))?)?([A-Z][a-z]?)")
        _TOL = 1e-9

        def _parse_match(match):
            tokens = _COMPONENT_RE.findall(match.group(1))
            elements = [element for _, _, element in tokens]
            # Parentheses containing oxygen belong to fixed formulas such as
            # Cr2(MoO4)3, not to a variable-composition phase signature.
            if len(elements) < 2 or "O" in elements:
                return None
            out = {}
            for lo_text, hi_text, element in tokens:
                if lo_text:
                    lo = float(lo_text)
                    hi = float(hi_text) if hi_text else lo
                else:
                    lo, hi = 0.0, 1.0
                out[element] = (lo, hi)
            return out or None

        def _parse_phase(lbl):
            for match in _PHASE_RE.finditer(lbl):
                parsed = _parse_match(match)
                if parsed is not None:
                    return parsed
            return None

        def _fixed(lbl):
            without_phase = _PHASE_RE.sub(
                lambda match: "" if _parse_match(match) is not None else match.group(0),
                lbl,
            )
            parts = sorted(part.strip() for part in without_phase.split("+") if part.strip())
            return "|".join(parts)

        def _connects(la, lb):
            if _fixed(la) != _fixed(lb):
                return False
            ra, rb = _parse_phase(la), _parse_phase(lb)
            if ra is None and rb is None:
                return la == lb
            if ra is None or rb is None:
                return False
            if set(ra) != set(rb):
                return False
            for metal in ra:
                lo_a, hi_a = ra[metal]
                lo_b, hi_b = rb[metal]
                if hi_a < lo_b - _TOL or hi_b < lo_a - _TOL:
                    return False
            return True

        grid = np.array(labels, dtype=object).reshape(nrows, ncols)
        visited = np.zeros((nrows, ncols), dtype=bool)
        result = grid.copy()

        for r0 in range(nrows):
            for c0 in range(ncols):
                if visited[r0, c0]:
                    continue
                visited[r0, c0] = True
                component = [(r0, c0)]
                queue = deque([(r0, c0)])
                while queue:
                    r, c = queue.popleft()
                    for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                        r2, c2 = r + dr, c + dc
                        if 0 <= r2 < nrows and 0 <= c2 < ncols and not visited[r2, c2]:
                            if _connects(grid[r, c], grid[r2, c2]):
                                visited[r2, c2] = True
                                component.append((r2, c2))
                                queue.append((r2, c2))

                canonical = Counter(grid[r, c] for r, c in component).most_common(1)[0][0]
                for r, c in component:
                    result[r, c] = canonical

        return result.ravel().tolist()

    def _solve_metrics(self, pd_data, comp, T, mu_o, phase_G):
        self._ensure_imports()
        from thermodynamics import KB, build_assemblage_labels, solve_grand_lp  # noqa

        phase_O = pd_data["phase_O"]
        n_fixed = pd_data["n_fixed"]
        fixed_ef = pd_data["fixed_energy_formula"]
        grand = np.concatenate([fixed_ef - phase_O[:n_fixed] * mu_o, phase_G])
        b_eq = np.concatenate([comp, [pd_data.get("phase_element_stoichiometry", 0.0)]])
        amounts, omega, ok = solve_grand_lp(pd_data["A_eq"], b_eq, grand)
        n = len(pd_data["phase_ids"])
        threshold = self.config.active_threshold
        if not ok:
            return {
                "ok": False,
                "amounts": np.full(n, np.nan),
                "fracs": np.full(n, np.nan),
                "omega": np.nan,
                "absorbed_O": np.nan,
                "phase_fraction": np.nan,
                "parent_phase_fraction": np.nan,
                "oxide_phase_fraction": np.nan,
                "exact_label": "no feasible assemblage",
                "family_label": "no feasible assemblage",
            }
        total = float(np.nansum(amounts))
        raw_fracs = amounts / total if total > 0 else amounts.copy()
        fracs = raw_fracs.copy()
        fracs[np.abs(fracs) < threshold] = 0.0
        phase_mask = pd_data["phase_kinds"] == "phase"
        phase_y_nd = np.asarray(pd_data["phase_y_nd"], dtype=float)
        # Preserve phase identity by element set, not by its initial ratios.
        parent_mask = phase_mask & np.all(phase_y_nd > 0.01, axis=1)
        oxide_mask = np.asarray(phase_O, dtype=float) > 1e-12
        exact_label, family_label = build_assemblage_labels(
            pd_data["phase_ids"],
            pd_data["phase_kinds"],
            fracs,
            threshold,
            pd_data.get("phase_label", ""),
            metals=pd_data.get("metals", []),
            phase_suffix=pd_data.get("phase_suffix", ""),
            family_threshold=(0.01 if self.config.include_0p01_to_0p05_components else 0.05),
            family_values=amounts,
            family_threshold_inclusive=self.config.include_0p01_to_0p05_components,
        )
        return {
            "ok": True,
            "amounts": amounts,
            "fracs": fracs,
            "omega": float(omega),
            "absorbed_O": float(np.nansum(amounts * phase_O)),
            "phase_fraction": float(np.nansum(fracs[phase_mask])),
            "parent_phase_fraction": float(np.nansum(raw_fracs[parent_mask])),
            "oxide_phase_fraction": float(np.nansum(raw_fracs[oxide_mask])),
            "exact_label": exact_label,
            "family_label": family_label,
        }

    def _find_external_diagram(self, metals, T=None) -> Path | None:
        base = self.config.ternary_diagrams_root / system_key(metals)
        if not base.exists():
            return None
        # Prefer clean GIF
        clean_gif = base / f"{system_key(metals)}_phase_evolution_clean.gif"
        if clean_gif.exists():
            return clean_gif
        if T is not None:
            for p in [
                base / "per_temperature" / f"{int(round(T))}K.png",
                base / "per_temperature" / f"{int(round(T))}K_gibbs.png",
            ]:
                if p.exists():
                    return p
        for p in [
            base / f"{system_key(metals)}_phase_map.png",
            base / f"{system_key(metals)}_phase_diagram.png",
            base / f"{system_key(metals)}_ternary.png",
        ]:
            if p.exists():
                return p
        # Case-insensitive fallback: prefer files containing "phase" and "diagram"
        pngs = sorted(base.glob("*.png"))
        for p in pngs:
            n = p.name.lower()
            if "phase" in n and "diagram" in n:
                return p
        for p in pngs:
            if "phase" in p.name.lower():
                return p
        return pngs[0] if pngs else None

    @staticmethod
    def _read_diagram_image(diag: Path, T=None, t_start: float = 200.0, gif_t_step: float = 10.0):
        if diag.suffix.lower() == ".gif":
            from PIL import Image as _PILImg

            im = _PILImg.open(diag)
            frame_idx = (
                max(0, int(round((float(T) - t_start) / gif_t_step)))
                if T is not None
                else getattr(im, "n_frames", 1) - 1
            )
            im.seek(min(frame_idx, getattr(im, "n_frames", 1) - 1))
            return np.array(im.convert("RGB"))
        import matplotlib.pyplot as _plt

        return _plt.imread(diag)

    def _plot_region_map_png(
        self,
        mu_values,
        x_values,
        region_grid,
        region_labels,
        sys_cfg,
        T,
        spec,
        out_dir,
        region_details=None,
        cell_text_grid=None,
        y_label=None,
        y_axis_label=None,
    ) -> Path:
        import matplotlib

        matplotlib.use("Agg")
        import hashlib as _hl

        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt

        self._ensure_imports()
        from thermodynamics import (
            add_region_annotation,
            format_phase_detail_line,
            grid_edges,
            region_annotation_text,
            separate_region_annotations,
            write_region_map_html,
        )

        cfg = self.config
        out_dir.mkdir(parents=True, exist_ok=True)
        cmap = plt.get_cmap("tab20")

        def _lc(lbl):
            return cmap(int(_hl.md5(lbl.encode()).hexdigest(), 16) % 20 / 20)

        colors = {r + 1: _lc(lbl) for r, lbl in enumerate(region_labels)}
        xe = grid_edges(x_values)
        me = grid_edges(mu_values)

        if cfg.region_label_mode == "id":
            fig, (ax, ax_leg) = plt.subplots(
                1, 2, figsize=(17, 8), gridspec_kw={"width_ratios": [3.2, 1.2], "wspace": 0.05}
            )
            ax_leg.axis("off")
        else:
            fig, ax = plt.subplots(figsize=(13, 8))
            ax_leg = None
        for ix in range(len(x_values)):
            for imu in range(len(mu_values)):
                rid = int(region_grid[ix, imu])
                ax.add_patch(
                    plt.Rectangle(
                        (me[imu], xe[ix]),
                        me[imu + 1] - me[imu],
                        xe[ix + 1] - xe[ix],
                        color=colors.get(rid, (0.9, 0.9, 0.9, 1.0)),
                        linewidth=0,
                    )
                )
        for ix in range(len(x_values)):
            for imu in range(len(mu_values) - 1):
                if region_grid[ix, imu] != region_grid[ix, imu + 1]:
                    ax.plot([me[imu + 1]] * 2, [xe[ix], xe[ix + 1]], "k-", lw=cfg.boundary_lw)
        for ix in range(len(x_values) - 1):
            for imu in range(len(mu_values)):
                if region_grid[ix, imu] != region_grid[ix + 1, imu]:
                    ax.plot([me[imu], me[imu + 1]], [xe[ix + 1]] * 2, "k-", lw=cfg.boundary_lw)
        annotation_artists = []
        for rid, _ in enumerate(region_labels, start=1):
            mask = region_grid == rid
            if not np.any(mask):
                continue
            rows, cols = np.where(mask)
            r = max(0, min(len(x_values) - 1, int(np.round(np.median(rows)))))
            c = max(0, min(len(mu_values) - 1, int(np.round(np.median(cols)))))
            annotation = region_annotation_text(
                rid,
                region_labels[rid - 1],
                region_details,
                cfg.region_label_mode,
            )
            annotation_artists.append(
                add_region_annotation(
                    ax,
                    mu_values[c],
                    x_values[r],
                    annotation,
                    cfg.region_label_fontsize,
                    colors[rid],
                )
            )
        axis_metal = sys_cfg.metals[spec["axis"]] if sys_cfg else ""
        resolved_y_label = y_label if y_label else f"initial x_{axis_metal}"
        ax.set_xlim(float(mu_values[0]), float(mu_values[-1]))
        ax.set_ylim(float(x_values[0]), float(x_values[-1]))
        ax.set_xlabel(r"$\mu_O$ (eV per O atom)", fontsize=13)
        ax.set_ylabel(resolved_y_label, fontsize=13)
        title_T = f" | T={T} K" if T is not None else ""
        ax.set_title(f"{sys_cfg.phase_label} | {spec['label']}{title_T}", fontsize=13)
        separate_region_annotations(ax, annotation_artists)

        def _format_detail(detail):
            if not detail:
                return ""
            phase_ranges = detail.get("phase_ranges") or []
            if phase_ranges:
                return format_phase_detail_line(phase_ranges)
            exact = detail.get("exact_label")
            return exact or ""

        handles = []
        for rid, lbl in enumerate(region_labels, start=1):
            if not np.any(region_grid == rid):
                continue
            detail = _format_detail(region_details.get(rid) if region_details else None)
            label = f"{rid}: {detail if detail else lbl}"
            handles.append(mpatches.Patch(color=colors[rid], label=label))
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
        out = out_dir / "assemblage_region_map.png"
        fig.savefig(out, dpi=170, bbox_inches="tight")
        try:
            write_region_map_html(
                mu_values,
                x_values,
                region_grid,
                region_labels,
                out_dir / "assemblage_region_map.html",
                title=f"{sys_cfg.phase_label} | {spec['label']}{title_T}",
                x_label="mu_O (eV per O atom)",
                y_label=resolved_y_label,
                region_details=region_details,
                cell_text_grid=cell_text_grid,
                region_label_mode=cfg.region_label_mode,
                region_label_fontsize=cfg.region_label_fontsize,
            )
        except Exception as e:
            print(f"[html skipped] {e}")
        plt.close(fig)
        return out

    def _plot_side_by_side(self, src_map, metals, T, spec, out_dir, composition=None) -> Path | None:
        if len(metals) < 2 or not self.config.make_side_by_side_miscibility:
            return None
        diag = self._find_external_diagram(metals, T)
        if not src_map.exists():
            return None
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(18, 7), gridspec_kw={"wspace": 0.04})
        axes[0].imshow(plt.imread(src_map))
        axes[0].axis("off")
        if diag is not None:
            is_gif = Path(diag).suffix.lower() == ".gif"
            if len(metals) == 3:
                # Ternary: show diagram with composition scan-line overlay
                sqrt3_2 = np.sqrt(3.0) / 2.0
                axes[1].imshow(self._read_diagram_image(diag, T), extent=(0, 1, 0, sqrt3_2), origin="upper")
                axes[1].set_xlim(0, 1)
                axes[1].set_ylim(0, sqrt3_2)
                axes[1].axis("off")
                x0 = float(self.config.slice_x_values[0])
                x1 = float(self.config.slice_x_values[-1])
                c0 = self._composition_from_slice(x0, 3, spec["axis"], spec["remainder"])
                c1 = self._composition_from_slice(x1, 3, spec["axis"], spec["remainder"])
                xa, ya = c0[1] + 0.5 * c0[2], sqrt3_2 * c0[2]
                xb, yb = c1[1] + 0.5 * c1[2], sqrt3_2 * c1[2]
                axes[1].plot([xa, xb], [ya, yb], color="black", lw=5, alpha=0.85)
                axes[1].plot([xa, xb], [ya, yb], color="white", lw=2.5, alpha=0.95)
                off = 0.05
                axes[1].text(0.50, sqrt3_2 + off, metals[2], ha="center", va="bottom", fontsize=12, fontweight="bold")
                axes[1].text(-off, -off * 0.4, metals[0], ha="right", va="top", fontsize=12, fontweight="bold")
                axes[1].text(1.0 + off, -off * 0.4, metals[1], ha="left", va="top", fontsize=12, fontweight="bold")
            else:
                # Binary (or higher-order): show diagram with T line + element labels.
                # BLADEvisualizer: x-axis = x_{M2} (0=pure M1 left, 1=pure M2 right)
                # y-axis: Temperature (K), bottom=T_min, top=T_max
                # data x range ≈ transAxes [0.12, 0.92], y range ≈ [0.12, 0.88]
                axes[1].imshow(self._read_diagram_image(diag, T), aspect="auto")
                axes[1].axis("off")
                if len(metals) == 2:
                    M1, M2 = metals[0], metals[1]
                    # Element labels at x-axis ends (x_M2: left=0=M1, right=1=M2)
                    x_left, x_right = 0.12, 0.92  # data x bounds in transAxes
                    y_xlabel = 0.08  # just below data area
                    for x_pos, ha, label in [
                        (x_left, "center", M1),
                        (x_right, "center", M2),
                    ]:
                        axes[1].text(
                            x_pos,
                            y_xlabel,
                            label,
                            transform=axes[1].transAxes,
                            ha=ha,
                            va="top",
                            fontsize=13,
                            fontweight="bold",
                            bbox=dict(fc="white", ec="none", alpha=0.7, pad=1),
                        )
                    if not is_gif and T is not None:
                        t0 = float(getattr(self.config, "phase_diagram_t_start", 300))
                        t1 = float(getattr(self.config, "phase_diagram_t_end", 4500))
                        y_data_lo, y_data_hi = 0.12, 0.88
                        if t1 > t0:
                            raw = min(1.0, max(0.0, (float(T) - t0) / (t1 - t0)))
                            y_line = y_data_lo + raw * (y_data_hi - y_data_lo)
                            axes[1].plot(
                                [x_left, x_right],
                                [y_line, y_line],
                                transform=axes[1].transAxes,
                                color="black",
                                lw=2.8,
                                solid_capstyle="round",
                            )
                            axes[1].plot(
                                [x_left, x_right],
                                [y_line, y_line],
                                transform=axes[1].transAxes,
                                color="white",
                                lw=1.2,
                                solid_capstyle="round",
                            )
                    if composition is not None:
                        comp = np.asarray(composition, dtype=float)
                        if comp.shape == (2,):
                            x_second = min(1.0, max(0.0, float(comp[1])))
                            x_line = x_left + x_second * (x_right - x_left)
                            y_data_lo, y_data_hi = 0.12, 0.88
                            axes[1].plot(
                                [x_line, x_line],
                                [y_data_lo, y_data_hi],
                                transform=axes[1].transAxes,
                                color="black",
                                lw=3.2,
                                solid_capstyle="round",
                            )
                            axes[1].plot(
                                [x_line, x_line],
                                [y_data_lo, y_data_hi],
                                transform=axes[1].transAxes,
                                color="white",
                                lw=1.4,
                                solid_capstyle="round",
                            )
            # Temperature stamp on all diagram types
            axes[1].text(
                0.50,
                0.02,
                (f"T = {int(T)} K  |  " if T is not None else "") + spec["label"],
                transform=axes[1].transAxes,
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight="bold",
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=2),
            )
        else:
            axes[1].axis("off")
            axes[1].text(
                0.5,
                0.5,
                "No external miscibility diagram found",
                transform=axes[1].transAxes,
                ha="center",
                va="center",
                fontsize=11,
            )
        fig.suptitle(f"{''.join(metals)}{self.config.phase_element or ''} | {spec['label']}", fontsize=13)
        out = out_dir / "assemblage_with_miscibility.png"
        fig.savefig(out, dpi=170, bbox_inches="tight")
        plt.close(fig)
        return out

    def _run_composition_slice_maps(self, sys_cfg, pd_data) -> None:
        import pandas as pd

        self._ensure_imports()

        from thermodynamics import KB, assign_region_ids, format_phase_fraction_summary_line

        cfg = self.config
        metals = sys_cfg.metals
        sys_name = system_key(metals, cfg.phase_element)
        specs = self._composition_slice_specs(metals)
        n_metals = len(metals)
        x_values = np.asarray(cfg.slice_x_values, dtype=float)
        mu_values = np.asarray(cfg.slice_mu_O, dtype=float)

        print(f"--- Composition slice maps ({len(specs)} slices) ---")
        slice_states = {
            spec["name"]: {
                "root": (cfg.figures_dir / sys_name / "composition_slice_maps" / spec["name"]),
                "map_frames": [],
                "panel_frames": [],
            }
            for spec in specs
        }
        reused_caches = 0
        calculated_caches = 0

        # Temperature-first scheduling gives every slice a chart before moving
        # to the next temperature. This also makes interrupted batch runs resume
        # across all compositions instead of repeatedly starting at slice one.
        for T in cfg.slice_T_values:
            for spec in specs:
                state = slice_states[spec["name"]]
                slice_root = state["root"]
                plots_root = slice_root / "muO_x_plots"
                plots_dir = plots_root / f"T{int(T)}"
                plots_dir.mkdir(parents=True, exist_ok=True)
                map_png = plots_dir / "assemblage_region_map.png"
                map_html = plots_dir / "assemblage_region_map.html"
                panel_png = plots_dir / "assemblage_with_miscibility.png"
                plot_settings_path = plots_dir / "plot_settings.json"
                expected_plot_settings = {
                    "annotation_style_version": 6,
                    "region_label_mode": cfg.region_label_mode,
                    "region_label_fontsize": cfg.region_label_fontsize,
                    "boundary_lw": cfg.boundary_lw,
                }
                try:
                    import json

                    plot_settings_match = (
                        plot_settings_path.exists()
                        and json.loads(plot_settings_path.read_text()) == expected_plot_settings
                    )
                except (OSError, ValueError, TypeError):
                    plot_settings_match = False

                # Reuse and reorganize figures made by the previous flat layout.
                legacy_files = {
                    plots_root / f"assemblage_region_map_{int(T)}K.png": map_png,
                    plots_root / f"assemblage_region_map_{int(T)}K.html": map_html,
                    plots_root / f"assemblage_with_miscibility_{int(T)}K.png": panel_png,
                }
                for old_path, new_path in legacy_files.items():
                    if old_path.exists() and not new_path.exists():
                        old_path.replace(new_path)
                csv_path = cfg.tables_dir / f"{sys_cfg.tag}_slice_{spec['name']}_T{int(T)}.csv"
                cache_ready = self._slice_cache_matches(csv_path, T, spec, metals, x_values, mu_values)
                if cache_ready and not cfg.run_plots:
                    reused_caches += 1
                    continue
                if not cache_ready and not cfg.run_calculations:
                    raise RuntimeError(f"plot-only mode requires the composition-slice cache: {csv_path}")
                if (
                    cache_ready
                    and plot_settings_match
                    and map_png.exists()
                    and map_html.exists()
                    and panel_png.exists()
                ):
                    reused_caches += 1
                    state["map_frames"].append((float(T), map_png))
                    state["panel_frames"].append((float(T), panel_png))
                    continue

                if cache_ready:
                    reused_caches += 1
                    df = pd.read_csv(csv_path)
                else:
                    calculated_caches += 1
                    phase_G = pd_data["phase_H0"] + KB * float(T) * pd_data["phase_mix_shape"]
                    rows = []
                    for x in x_values:
                        comp = self._composition_from_slice(float(x), n_metals, spec["axis"], spec["remainder"])
                        for mu_o in mu_values:
                            r = self._solve_metrics(pd_data, comp, float(T), float(mu_o), phase_G)
                            family = (
                                self._coarse_family_label(
                                    pd_data["phase_ids"],
                                    pd_data["phase_kinds"],
                                    r["amounts"],
                                    pd_data["phase_y_nd"],
                                    metals,
                                )
                                if r["ok"]
                                else "no feasible assemblage"
                            )
                            rows.append(
                                {
                                    "T_K": float(T),
                                    "muO_eV_per_O": float(mu_o),
                                    f"x_{metals[spec['axis']]}_axis": float(x),
                                    **{f"x_{m}": float(comp[i]) for i, m in enumerate(metals)},
                                    "Omega_eV_per_initial_formula": r["omega"],
                                    "absorbed_O_atoms": r["absorbed_O"],
                                    "phase_fraction": r["phase_fraction"],
                                    "parent_phase_fraction": r["parent_phase_fraction"],
                                    "oxide_phase_fraction": r["oxide_phase_fraction"],
                                    "assemblage_exact": r["exact_label"],
                                    "phase_fraction_summary": self._phase_fraction_summary(
                                        r["fracs"],
                                        pd_data["phase_ids"],
                                        pd_data["phase_kinds"],
                                        metals,
                                        pd_data.get("phase_suffix", ""),
                                    ),
                                    "assemblage_family_coarse": family,
                                    "component_presence_threshold": (
                                        0.01 if cfg.include_0p01_to_0p05_components else 0.05
                                    ),
                                }
                            )
                    df = pd.DataFrame(rows)
                    df.to_csv(csv_path, index=False)

                if not cfg.run_plots:
                    continue

                labels = normalize_region_labels(df["assemblage_family_coarse"].astype(str).tolist())
                labels = self._merge_same_component_adjacent_regions(labels, len(x_values), len(mu_values))
                region_ids, region_labels = assign_region_ids(labels)
                region_grid = region_ids.reshape(len(x_values), len(mu_values))
                region_details = self._build_region_details_from_csv(
                    df,
                    region_ids,
                    pd_data["phase_ids"],
                    pd_data["phase_kinds"],
                    metals,
                    pd_data.get("phase_suffix", ""),
                )
                cell_texts = np.array(
                    [format_phase_fraction_summary_line(s) for s in df["phase_fraction_summary"].astype(str).tolist()],
                    dtype=object,
                ).reshape(len(x_values), len(mu_values))
                if not plot_settings_match or not (map_png.exists() and map_html.exists()):
                    generated_map = self._plot_region_map_png(
                        mu_values,
                        x_values,
                        region_grid,
                        region_labels,
                        sys_cfg,
                        int(T),
                        spec,
                        plots_dir,
                        region_details=region_details,
                        cell_text_grid=cell_texts,
                    )
                    if generated_map and generated_map.exists() and generated_map != map_png:
                        generated_map.replace(map_png)
                    html_src = plots_dir / "assemblage_region_map.html"
                    if html_src.exists() and html_src != map_html:
                        html_src.replace(map_html)
                if map_png.exists() and (not plot_settings_match or not panel_png.exists()):
                    generated_panel = self._plot_side_by_side(map_png, metals, int(T), spec, plots_dir)
                    if generated_panel and generated_panel.exists() and generated_panel != panel_png:
                        generated_panel.replace(panel_png)
                if map_png.exists() and map_html.exists():
                    plot_settings_path.write_text(json.dumps(expected_plot_settings, indent=2, sort_keys=True) + "\n")
                if map_png.exists():
                    state["map_frames"].append((float(T), map_png))
                if panel_png.exists():
                    state["panel_frames"].append((float(T), panel_png))

        for spec in specs:
            state = slice_states[spec["name"]]
            slice_root = state["root"]
            map_frames = [path for _, path in sorted(state["map_frames"])]
            panel_frames = [path for _, path in sorted(state["panel_frames"])]
            if cfg.run_animations:
                make_animation(
                    map_frames,
                    slice_root / "muO_x.gif",
                    slice_root / "muO_x.mp4",
                    cfg.anim_fps,
                    cfg.mp4_crf,
                    cfg.mp4_preset,
                )
                if panel_frames:
                    make_animation(
                        panel_frames,
                        slice_root / "muO_x_miscibility.gif",
                        slice_root / "muO_x_miscibility.mp4",
                        cfg.anim_fps,
                        cfg.mp4_crf,
                        cfg.mp4_preset,
                    )
        print(
            f"  Composition slice data: reused {reused_caches} cache(s); "
            f"calculated {calculated_caches} missing cache(s)"
        )

    @staticmethod
    def _phase_fraction_summary(fracs, phase_ids, phase_kinds, metals, phase_suffix=""):
        from collections import defaultdict

        from thermodynamics import _phase_comp_label, _short_label

        totals = defaultdict(float)
        min_frac = 1e-12
        any_valid = False
        for pid, kind, frac in zip(phase_ids, phase_kinds, fracs):
            if np.isnan(frac) or float(frac) <= min_frac:
                continue
            any_valid = True
            name = _phase_comp_label(pid, metals, phase_suffix=phase_suffix) if kind == "phase" else _short_label(pid)
            totals[name] += float(frac)

        pairs = sorted(totals.items(), key=lambda x: x[1], reverse=True)
        if not pairs and len(fracs) > 0 and any_valid:
            idx = int(np.nanargmax(np.nan_to_num(fracs, nan=-1.0)))
            if idx >= 0:
                pid = phase_ids[idx]
                kind = phase_kinds[idx]
                name = (
                    _phase_comp_label(pid, metals, phase_suffix=phase_suffix) if kind == "phase" else _short_label(pid)
                )
                pairs = [(name, float(fracs[idx]))]
        return " | ".join(f"{name}={frac:.6f}" for name, frac in pairs)

    @staticmethod
    def _parse_phase_fraction_summary(summary: str) -> list[tuple[str, float]]:
        if not isinstance(summary, str) or not summary.strip():
            return []
        out = []
        for tok in summary.split("|"):
            tok = tok.strip()
            if "=" not in tok:
                continue
            name, val = tok.rsplit("=", 1)
            try:
                out.append((name.strip(), float(val.strip())))
            except ValueError:
                continue
        return out

    def _build_region_details_from_csv(self, df, region_ids, phase_ids, phase_kinds, metals, phase_suffix=""):
        from collections import Counter, defaultdict

        from thermodynamics import _phase_comp_label, _short_label

        details = {}
        min_frac = 1e-12
        exact_labels = (
            df["assemblage_exact"].astype(str).values
            if "assemblage_exact" in df.columns
            else np.array([""] * len(df), dtype=object)
        )
        frac_summaries = (
            df["phase_fraction_summary"].astype(str).values
            if "phase_fraction_summary" in df.columns
            else np.array([""] * len(df), dtype=object)
        )
        family_aliases = defaultdict(list)
        for pid, kind in zip(phase_ids, phase_kinds):
            family = _phase_comp_label(pid, metals, phase_suffix=phase_suffix) if kind == "phase" else _short_label(pid)
            family_aliases[family].append(family)
            family_aliases[family].append(str(pid).split("_")[-1])
        for rid in np.unique(region_ids):
            mask = region_ids == rid
            if not np.any(mask):
                continue
            exact_vals = [v for v in exact_labels[mask].tolist() if v]
            exact_label = Counter(exact_vals).most_common(1)[0][0] if exact_vals else ""
            phase_maps = [dict(self._parse_phase_fraction_summary(s)) for s in frac_summaries[mask].tolist()]
            ranges = []
            seen = set()
            for pid, kind in zip(phase_ids, phase_kinds):
                name = (
                    _phase_comp_label(pid, metals, phase_suffix=phase_suffix) if kind == "phase" else _short_label(pid)
                )
                if name in seen:
                    continue
                seen.add(name)
                aliases = list(dict.fromkeys(family_aliases.get(name, [name])))  # deduplicated
                vals = []
                for m in phase_maps:
                    total = 0.0
                    found = False
                    for alias in aliases:
                        if alias in m:
                            total += float(m[alias])
                            found = True
                    if found:
                        vals.append(min(total, 1.0))  # clamp to [0,1]
                vals = [v for v in vals if v > min_frac]
                display_threshold = 0.01 if self.config.include_0p01_to_0p05_components else 0.05
                if vals and max(vals) >= display_threshold - 1e-12:
                    ranges.append((name, float(min(vals)), float(max(vals))))
            details[int(rid)] = {
                "exact_label": exact_label,
                "phase_ranges": ranges,
            }
        return details

    def _run_composition_slice_muT_maps(self, sys_cfg, pd_data) -> None:
        """Generate μO–T assemblage maps at fixed compositions along each slice.

        Perpendicular to the μO–x slices: here composition is fixed at each
        sample point along the slice axis, and T varies.  Saves PNGs + HTMLs
        into composition_slice_maps/{spec}/{x_label}/muO_T/ and places
        miscibility-gap diagrams beside them.
        """
        import pandas as pd

        self._ensure_imports()
        from thermodynamics import KB, assign_region_ids, format_phase_fraction_summary_line

        cfg = self.config
        metals = sys_cfg.metals
        sys_name = system_key(metals, cfg.phase_element)
        specs = self._composition_slice_specs(metals)
        n_metals = len(metals)
        x_values = np.asarray(cfg.slice_muT_x_values, dtype=float)
        mu_values = np.asarray(cfg.slice_mu_O, dtype=float)
        T_values = np.asarray(cfg.slice_T_values, dtype=float)

        print(
            f"--- Composition slice muO–T maps ({len(specs)} slices, "
            f"{len(x_values)} x × {len(T_values)} T × {len(mu_values)} muO) ---"
        )
        for spec in specs:
            slice_root = cfg.figures_dir / sys_name / "composition_slice_maps" / spec["name"]
            plots_dir = slice_root / "muO_T_plots"
            plots_dir.mkdir(parents=True, exist_ok=True)
            map_frames, panel_frames = [], []

            for x in x_values:
                x_round = round(float(x), 6)  # avoid 0.30000000001
                x_disp = f"{x_round:.2f}"  # "0.30"
                x_str = x_disp.replace(".", "p")  # "0p30"
                comp = self._composition_from_slice(float(x), n_metals, spec["axis"], spec["remainder"])
                axis_m = metals[spec["axis"]]
                csv_path = cfg.tables_dir / f"{sys_cfg.tag}_slice_muT_{spec['name']}_{x_str}.csv"
                if self._slice_muT_cache_matches(
                    csv_path,
                    comp,
                    axis_m,
                    x,
                    metals,
                    T_values,
                    mu_values,
                ):
                    df = pd.read_csv(csv_path)
                else:
                    if not cfg.run_calculations:
                        raise RuntimeError(f"plot-only mode requires the composition-slice muO-T cache: {csv_path}")
                    rows = []
                    for T in T_values:
                        phase_G = pd_data["phase_H0"] + KB * float(T) * pd_data["phase_mix_shape"]
                        for mu_o in mu_values:
                            r = self._solve_metrics(pd_data, comp, float(T), float(mu_o), phase_G)
                            family = (
                                self._coarse_family_label(
                                    pd_data["phase_ids"],
                                    pd_data["phase_kinds"],
                                    r["amounts"],
                                    pd_data["phase_y_nd"],
                                    metals,
                                )
                                if r["ok"]
                                else "no feasible assemblage"
                            )
                            rows.append(
                                {
                                    "T_K": float(T),
                                    "muO_eV_per_O": float(mu_o),
                                    f"x_{axis_m}_axis": float(x),
                                    **{f"x_{m}": float(comp[i]) for i, m in enumerate(metals)},
                                    "phase_fraction": r["phase_fraction"],
                                    "parent_phase_fraction": r["parent_phase_fraction"],
                                    "oxide_phase_fraction": r["oxide_phase_fraction"],
                                    "assemblage_exact": r["exact_label"],
                                    "phase_fraction_summary": self._phase_fraction_summary(
                                        r["fracs"],
                                        pd_data["phase_ids"],
                                        pd_data["phase_kinds"],
                                        metals,
                                        pd_data.get("phase_suffix", ""),
                                    ),
                                    "assemblage_family_coarse": family,
                                    "component_presence_threshold": (
                                        0.01 if cfg.include_0p01_to_0p05_components else 0.05
                                    ),
                                }
                            )
                    df = pd.DataFrame(rows)
                    df.to_csv(csv_path, index=False)

                if not cfg.run_plots:
                    continue

                labels = normalize_region_labels(df["assemblage_family_coarse"].astype(str).tolist())
                labels = self._merge_same_component_adjacent_regions(labels, len(T_values), len(mu_values))
                region_ids, region_labels = assign_region_ids(labels)
                region_grid = region_ids.reshape(len(T_values), len(mu_values))
                region_details = self._build_region_details_from_csv(
                    df,
                    region_ids,
                    pd_data["phase_ids"],
                    pd_data["phase_kinds"],
                    metals,
                    pd_data.get("phase_suffix", ""),
                )
                cell_texts = np.array(
                    [format_phase_fraction_summary_line(s) for s in df["phase_fraction_summary"].astype(str)],
                    dtype=object,
                ).reshape(len(T_values), len(mu_values))

                # Modify spec label to include current composition
                spec_with_x = dict(spec)
                spec_with_x["label"] = f"x_{axis_m}={x_disp}"

                map_png = self._plot_region_map_png(
                    mu_values,
                    T_values,
                    region_grid,
                    region_labels,
                    sys_cfg,
                    None,
                    spec_with_x,
                    plots_dir,
                    region_details=region_details,
                    cell_text_grid=cell_texts,
                    y_label="T (K)",
                )
                # Rename with composition suffix
                if map_png and map_png.exists():
                    dst = plots_dir / f"assemblage_region_map_x{x_str}.png"
                    map_png.rename(dst)
                    map_png = dst
                    html_src = plots_dir / "assemblage_region_map.html"
                    if html_src.exists():
                        html_src.rename(plots_dir / f"assemblage_region_map_x{x_str}.html")
                panel_png = self._plot_side_by_side(
                    map_png,
                    metals,
                    None,
                    spec_with_x,
                    plots_dir,
                    composition=comp,
                )
                if panel_png and panel_png.exists():
                    dst = plots_dir / f"assemblage_with_miscibility_x{x_str}.png"
                    panel_png.rename(dst)
                    panel_png = dst
                if map_png:
                    map_frames.append(map_png)
                if panel_png:
                    panel_frames.append(panel_png)

            # GIF/MP4 sit outside muO_T_plots/, labelled muO_T.*
            if cfg.run_animations:
                make_animation(
                    map_frames,
                    slice_root / "muO_T.gif",
                    slice_root / "muO_T.mp4",
                    cfg.anim_fps,
                    cfg.mp4_crf,
                    cfg.mp4_preset,
                )
                if panel_frames:
                    make_animation(
                        panel_frames,
                        slice_root / "muO_T_miscibility.gif",
                        slice_root / "muO_T_miscibility.mp4",
                        cfg.anim_fps,
                        cfg.mp4_crf,
                        cfg.mp4_preset,
                    )

    def _composition_grid_for_onset(self, n_metals: int):
        self._ensure_imports()
        from thermodynamics import simplex_grid_nd

        cfg = self.config
        step = (
            cfg.onset_comp_step_binary
            if n_metals == 2
            else cfg.onset_comp_step_ternary if n_metals == 3 else cfg.onset_comp_step_higher
        )
        grid = simplex_grid_nd(n_metals, step)
        while len(grid) > cfg.onset_max_comp_points and step < 0.5:
            step *= 2
            grid = simplex_grid_nd(n_metals, step)
        return grid, step

    @staticmethod
    def _state_key(comp, mu_o):
        """Stable coordinate key for an already calculated equilibrium state."""
        return tuple(round(float(x), 10) for x in comp) + (round(float(mu_o), 10),)

    def _load_reusable_scalar_states(self, sys_cfg, T):
        """Read onset inputs from existing full-grid calculation tables."""
        import pandas as pd

        metals = sys_cfg.metals
        reusable = {}
        paths = sorted(self.config.tables_dir.glob(f"{sys_cfg.tag}_slice_*_T{int(T)}.csv"))
        if len(metals) == 2:
            paths.append(self.config.tables_dir / f"{sys_cfg.tag}_muO_x_phase_map_T{int(T)}_cache.csv")

        required = {
            "muO_eV_per_O",
            "phase_fraction",
            "parent_phase_fraction",
            "oxide_phase_fraction",
            "absorbed_O_atoms",
        }
        for path in paths:
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path)
            except Exception:
                continue
            if not required.issubset(df.columns):
                continue

            metal_cols = [f"x_{m}" for m in metals]
            has_full_composition = all(c in df.columns for c in metal_cols)
            if not has_full_composition and len(metals) != 2:
                continue

            for values in df.to_dict("records"):
                if has_full_composition:
                    comp = [values[c] for c in metal_cols]
                else:
                    x0 = float(values[f"x_{metals[0]}"])
                    comp = [x0, 1.0 - x0]
                reusable[self._state_key(comp, values["muO_eV_per_O"])] = (
                    float(values["phase_fraction"]),
                    float(values["parent_phase_fraction"]),
                    float(values["oxide_phase_fraction"]),
                    float(values["absorbed_O_atoms"]),
                )
        return reusable

    def _run_onset_auc_diagrams(self, sys_cfg, pd_data) -> None:
        import pandas as pd

        self._ensure_imports()
        from thermodynamics import KB

        cfg = self.config
        metals = sys_cfg.metals
        n_metals = len(metals)
        comp_grid, comp_step = self._composition_grid_for_onset(n_metals)
        mu_values = np.asarray(cfg.onset_auc_mu_O, dtype=float)
        sys_name = system_key(metals, cfg.phase_element)
        out_root = cfg.figures_dir / sys_name / "onset_auc"
        out_root.mkdir(parents=True, exist_ok=True)

        print(f"--- Onset diagrams ({len(comp_grid)} comps, step={comp_step:g}) ---")
        onset_frames = []
        all_summary = []
        for T in cfg.onset_auc_T_values:
            csv_path = cfg.tables_dir / f"{sys_cfg.tag}_onset_auc_T{int(T)}.csv"
            if self._onset_cache_matches(csv_path, comp_grid, comp_step, mu_values, T, metals):
                df = pd.read_csv(csv_path)
            else:
                if not cfg.run_calculations:
                    raise RuntimeError(f"plot-only mode requires a current onset cache: {csv_path}")
                phase_G = pd_data["phase_H0"] + KB * float(T) * pd_data["phase_mix_shape"]
                reusable = self._load_reusable_scalar_states(sys_cfg, T)
                reused_states = 0
                solved_states = 0
                rows = []
                for comp in comp_grid:
                    phase_fraction_curve, parent_fraction_curve, oxide_fraction_curve, sampled_mu = [], [], [], []
                    onset_mu = np.nan
                    first_infeasible_mu = np.nan
                    for mu_o in mu_values:
                        cached = reusable.get(self._state_key(comp, mu_o))
                        if cached is not None:
                            phase_fraction, parent_fraction, oxide_fraction, absorbed_o = cached
                            ok = not np.isnan(phase_fraction)
                            reused_states += 1
                        else:
                            r = self._solve_metrics(pd_data, comp, float(T), float(mu_o), phase_G)
                            phase_fraction = r["phase_fraction"] if r["ok"] else np.nan
                            parent_fraction = r["parent_phase_fraction"] if r["ok"] else np.nan
                            oxide_fraction = r["oxide_phase_fraction"] if r["ok"] else np.nan
                            absorbed_o = r["absorbed_O"] if r["ok"] else np.nan
                            ok = r["ok"]
                            solved_states += 1
                        if not ok:
                            first_infeasible_mu = float(mu_o)
                            break
                        phase_fraction_curve.append(phase_fraction)
                        parent_fraction_curve.append(parent_fraction)
                        oxide_fraction_curve.append(oxide_fraction)
                        sampled_mu.append(float(mu_o))
                        if np.isnan(onset_mu) and absorbed_o > cfg.onset_threshold:
                            onset_mu = float(mu_o)
                    curve = np.asarray(phase_fraction_curve, dtype=float)
                    parent_curve = np.asarray(parent_fraction_curve, dtype=float)
                    oxide_curve = np.asarray(oxide_fraction_curve, dtype=float)
                    sampled_mu_array = np.asarray(sampled_mu, dtype=float)
                    if len(curve) == 0:
                        auc = np.nan
                        parent_auc = np.nan
                        oxide_auc = np.nan
                        feasible_mu_max = np.nan
                    else:
                        auc = 0.0 if len(curve) == 1 else float(np.trapezoid(curve, sampled_mu_array))
                        parent_auc = (
                            0.0 if len(parent_curve) == 1 else float(np.trapezoid(parent_curve, sampled_mu_array))
                        )
                        oxide_auc = 0.0 if len(oxide_curve) == 1 else float(np.trapezoid(oxide_curve, sampled_mu_array))
                        feasible_mu_max = float(sampled_mu_array[-1])
                    rows.append(
                        {
                            "T_K": float(T),
                            **{f"x_{m}": float(comp[i]) for i, m in enumerate(metals)},
                            "onset_muO_eV": onset_mu,
                            "phase_fraction_auc_eV": auc,
                            "parent_phase_fraction_auc_eV": parent_auc,
                            "oxide_phase_fraction_auc_eV": oxide_auc,
                            "phase_element": cfg.phase_element or "",
                            "feasible_muO_max_eV": feasible_mu_max,
                            "first_infeasible_muO_eV": first_infeasible_mu,
                            "feasible_muO_points": len(sampled_mu),
                            "comp_step": float(comp_step),
                            "muO_min": float(mu_values[0]),
                            "muO_max": float(mu_values[-1]),
                            "muO_step": float(mu_values[1] - mu_values[0]) if len(mu_values) > 1 else np.nan,
                        }
                    )
                df = pd.DataFrame(rows)
                df.to_csv(csv_path, index=False)
                print(f"    T={int(T)} K: reused {reused_states} states, " f"solved {solved_states} new states")
            all_summary.append(df)

            if not cfg.run_plots:
                continue

            t_dir = out_root / f"T{int(T)}"
            t_dir.mkdir(parents=True, exist_ok=True)
            if n_metals == 2:
                onset_png = self._plot_onset_auc_binary(df, sys_cfg, int(T), t_dir)
            elif n_metals == 3:
                onset_png = self._plot_onset_auc_ternary(df, sys_cfg, int(T), t_dir)
            else:
                onset_png = None
                (t_dir / "README.txt").write_text(
                    f"{''.join(metals)}{cfg.phase_element or ''} has {n_metals} metals.\n"
                    "n-ary onset written to CSV; no static simplex plot above ternary.\n"
                )
            if onset_png:
                onset_frames.append(onset_png)

        if all_summary:
            import pandas as pd

            pd.concat(all_summary, ignore_index=True).to_csv(
                cfg.tables_dir / f"{sys_cfg.tag}_onset_auc_all_temperatures.csv",
                index=False,
            )
        if cfg.run_animations:
            make_animation(
                onset_frames,
                out_root / "onset_diagram.gif",
                out_root / "onset_diagram.mp4",
                cfg.anim_fps,
                cfg.mp4_crf,
                cfg.mp4_preset,
            )

    @staticmethod
    def _plot_onset_auc_binary(df, sys_cfg, T, out_dir):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        m0 = sys_cfg.metals[0]
        x = df[f"x_{m0}"].values
        order = np.argsort(x)
        x = x[order]
        onset = df["onset_muO_eV"].values[order]
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(x, onset, "o-", lw=1.8, ms=4)
        ax.set_xlabel(f"x_{m0}")
        ax.set_ylabel(r"onset $\mu_O$ (eV/O)")
        ax.set_title(f"{sys_cfg.phase_label} onset | T={T} K")
        ax.grid(True, alpha=0.35)
        onset_png = out_dir / "onset_diagram.png"
        fig.savefig(onset_png, dpi=180, bbox_inches="tight")
        plt.close(fig)
        return onset_png

    @staticmethod
    def _ternary_frame(ax, metals):
        sqrt3_2 = np.sqrt(3.0) / 2.0
        tri = np.array([[0.5, sqrt3_2], [0.0, 0.0], [1.0, 0.0], [0.5, sqrt3_2]])
        ax.plot(tri[:, 0], tri[:, 1], "k-", lw=1.5)
        off = 0.06
        ax.text(-off, -off * 0.7, metals[0], ha="right", va="top", fontsize=12, fontweight="bold")
        ax.text(0.5, sqrt3_2 + off, metals[2], ha="center", va="bottom", fontsize=12, fontweight="bold")
        ax.text(1 + off, -off * 0.7, metals[1], ha="left", va="top", fontsize=12, fontweight="bold")
        ax.set_xlim(-0.08, 1.08)
        ax.set_ylim(-0.08, sqrt3_2 + 0.14)
        ax.set_aspect("equal")
        ax.axis("off")

    @staticmethod
    def _plot_ternary_scalar(df, sys_cfg, T, out_dir, value_col, title, color_label, filename, cmap="viridis"):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as tri_mod

        metals = sys_cfg.metals
        x2 = df[f"x_{metals[1]}"].values
        x3 = df[f"x_{metals[2]}"].values
        vals = df[value_col].astype(float).values
        xcart = x2 + 0.5 * x3
        ycart = (np.sqrt(3.0) / 2.0) * x3
        triang = tri_mod.Triangulation(xcart, ycart)
        tri_vals = vals[triang.triangles]
        triang.set_mask(np.any(np.isnan(tri_vals), axis=1))
        fig, ax = plt.subplots(figsize=(8.5, 7.8))
        SystemAnalyzer._ternary_frame(ax, metals)
        valid = ~np.isnan(vals)
        if valid.sum() >= 3:
            v_min, v_max = float(np.nanmin(vals)), float(np.nanmax(vals))
            if v_max > v_min:
                pc = ax.tricontourf(triang, vals, levels=30, cmap=cmap, vmin=v_min, vmax=v_max)
            else:
                pc = ax.tripcolor(triang, vals, shading="flat", cmap=cmap)
        else:
            pc = ax.tripcolor(triang, vals, shading="flat", cmap=cmap)
        fig.colorbar(pc, ax=ax, fraction=0.04, pad=0.03, label=color_label)
        ax.set_title(f"{sys_cfg.phase_label} {title} | T={T} K", fontsize=12)
        out = out_dir / filename
        fig.savefig(out, dpi=190, bbox_inches="tight")
        plt.close(fig)
        return out

    @classmethod
    def _plot_onset_auc_ternary(cls, df, sys_cfg, T, out_dir):
        return cls._plot_ternary_scalar(
            df,
            sys_cfg,
            T,
            out_dir,
            "onset_muO_eV",
            "oxidation onset",
            r"onset $\mu_O$ (eV/O)",
            "onset_diagram.png",
            cmap="RdYlGn",
        )
