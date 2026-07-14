"""BatchRunner — iterates over all system folders in system_structures_root."""

from __future__ import annotations

import sys
import traceback
from dataclasses import replace
from pathlib import Path

from .analyzer import SystemAnalyzer
from .config import Config
from .utils import (
    make_animation,
    prepare_system_tables_dir,
    system_key,
    system_tag,
    tile_images,
)


class BatchRunner:
    """Run analyses across all system folders.

    Usage::

        from ._batch_internal import BatchRunner; from .config import Config
        cfg = Config(system_structures_root="...", ...)
        runner = BatchRunner(cfg)
        tags_info = runner.run()   # returns list for Compiler
    """

    def __init__(self, config: Config):
        self.config = config

    def _ensure_imports(self):
        framework_dir = str(Path(__file__).parent)
        if framework_dir in sys.path:
            sys.path.remove(framework_dir)
        sys.path.insert(0, framework_dir)

    def _detect_metals(self, system_dir: Path) -> list[str] | None:
        self._ensure_imports()
        from .system_config import SystemConfig

        try:
            return SystemConfig._detect_metals_from_dirs(
                system_dir / self.config.blade_subdir,
                {self.config.phase_element} if self.config.phase_element else set(),
            )
        except Exception as e:
            print(f"  metal detection failed: {e}")
            return None

    def _check_l0_filter(self, tag: str, config: Config | None = None) -> bool:
        """Return True if system should be included (L0 ≤ 0 or filter disabled)."""
        import pandas as pd

        cfg = config or self.config
        if not cfg.filter_no_miscibility_gap:
            return True
        rk_file = cfg.tables_dir / f"{tag}_phase_interaction_coefficients.csv"
        try:
            rk_df = pd.read_csv(rk_file)
            l0_rows = rk_df[rk_df["term"] == "L0"]
            if "pair" in rk_df.columns:
                l0_rows = l0_rows[l0_rows["pair"].str.count("-") == 1]
            l0_vals = l0_rows["value_eV_per_formula"].astype(float).values
            max_l0 = float(l0_vals.max()) if len(l0_vals) > 0 else 0.0
            if max_l0 > 0:
                pairs = (
                    l0_rows.loc[l0_rows["value_eV_per_formula"] > 0, "pair"].tolist()
                    if "pair" in rk_df.columns
                    else ["L0"]
                )
                print(f"  L0={max_l0:+.4f} eV > 0 ({', '.join(pairs)}) → miscibility gap → filtered")
                return False
            print(f"  L0={max_l0:+.4f} eV ≤ 0 → included")
        except Exception as e:
            print(f"  L0 check failed ({e}) → included by default")
        return True

    def _compile_slice_video_grid(self, metals: list[str]) -> None:
        if not self.config.compile_slice_video_grids:
            return
        cfg = self.config
        sys_name = system_key(metals, cfg.phase_element)
        base = cfg.figures_dir / sys_name / "composition_slice_maps"
        if not base.exists():
            return
        grid_dir = base / "_grid"
        frames = []
        for T in cfg.slice_T_values:
            imgs = sorted(base.glob(f"*/T{int(T)}/assemblage_with_miscibility.png"))
            if not imgs:
                imgs = sorted(base.glob(f"*/T{int(T)}/assemblage_region_map.png"))
            out = tile_images(imgs, grid_dir / f"T{int(T)}_grid.png", cols=3)
            if out:
                frames.append(out)
        make_animation(
            frames,
            grid_dir / "slice_maps_grid.gif",
            grid_dir / "slice_maps_grid.mp4",
            cfg.anim_fps,
            cfg.mp4_crf,
            cfg.mp4_preset,
        )

    def _compile_global_slice_video_grid(self, tags_info: list) -> None:
        if not self.config.compile_slice_video_grids:
            return
        cfg = self.config
        global_dir = cfg.figures_dir / "compiled_slice_map_video_grid"
        frames = []
        systems = [sn for _, _, _, sn in tags_info][: cfg.max_global_video_grid_items]
        for T in cfg.slice_T_values:
            imgs = [cfg.figures_dir / sn / "composition_slice_maps" / "_grid" / f"T{int(T)}_grid.png" for sn in systems]
            out = tile_images(imgs, global_dir / f"T{int(T)}_grid.png", cols=3)
            if out:
                frames.append(out)
        make_animation(
            frames,
            global_dir / "compiled_slice_maps_grid.gif",
            global_dir / "compiled_slice_maps_grid.mp4",
            cfg.anim_fps,
            cfg.mp4_crf,
            cfg.mp4_preset,
        )

    def run(self) -> list:
        """Run all systems. Returns tags_info list for Compiler."""
        cfg = self.config
        cfg.tables_dir.mkdir(parents=True, exist_ok=True)
        cfg.figures_dir.mkdir(parents=True, exist_ok=True)
        self._ensure_imports()
        import pandas as pd

        from .system_config import SystemConfig, prepare_tables

        root = cfg.system_structures_root
        all_dirs = sorted(root.iterdir()) if root.exists() else []
        if cfg.systems:
            requested = set(cfg.systems)
            available = {directory.name for directory in all_dirs if directory.is_dir()}
            missing = sorted(requested - available)
            if missing:
                print(f"Requested system folders not found: {', '.join(missing)}")
            all_dirs = [directory for directory in all_dirs if directory.name in requested]
        blade_subdir = cfg.blade_subdir
        all_systems = [d for d in all_dirs if d.is_dir() and (d / blade_subdir).is_dir()]
        skipped_unary = [d.name for d in all_dirs if d.is_dir() and not (d / blade_subdir).is_dir()]

        print(f"Found {len(all_systems)} systems, {len(skipped_unary)} single-phase folders (skipped).")
        if skipped_unary:
            print(f"Skipped single-phase folders: {', '.join(skipped_unary)}\n")

        done, failed, skipped, filtered = [], [], [], []

        for i, system_dir in enumerate(all_systems):
            name = system_dir.name
            if name in cfg.skip_systems:
                print(f"[{name}] skipped (in skip_systems)")
                skipped.append(name)
                continue

            metals = self._detect_metals(system_dir)
            if metals is None or len(metals) < 2:
                print(f"[{name}] expected 2+ metals, got {metals} — skipping")
                failed.append((name, f"detected {metals}"))
                continue
            if cfg.elements and not set(metals).issubset(set(cfg.elements)):
                print(f"[{name}] skipped (contains metals outside selected elements {cfg.elements})")
                skipped.append(name)
                continue

            m1, m2 = metals[0], metals[1]
            m3 = metals[2] if len(metals) >= 3 else None
            tag = system_tag(metals, cfg.phase_element)
            table_dir = prepare_system_tables_dir(cfg.tables_dir, metals, cfg.phase_element)
            system_cfg = replace(cfg, tables_dir=table_dir)

            # Ensure tables exist for L0 check
            cfg_tmp = SystemConfig.resolve(
                system_dir,
                table_dir,
                metals,
                blade_subdir=cfg.blade_subdir,
                fixed_phases_subdir=cfg.fixed_phases_subdir,
                phase_element=cfg.phase_element,
                phase_element_stoichiometry=cfg.phase_element_stoichiometry,
            )
            if not cfg_tmp.tables_ready():
                if not cfg.run_calculations:
                    missing = ", ".join(str(path.name) for path in cfg_tmp.missing_table_files())
                    message = f"plot-only mode requires prepared equilibrium tables: {missing}"
                    print(f"[{name}] {message} — skipping")
                    failed.append((name, message))
                    continue
                try:
                    missing = ", ".join(str(p.name) for p in cfg_tmp.missing_table_files())
                    if missing:
                        print(f"[{name}] missing tables: {missing} — rebuilding")
                    prepare_tables(system_dir, cfg_tmp, cfg.rk_order)
                except Exception as e:
                    print(f"[{name}] table prep failed: {e} — skipping")
                    failed.append((name, str(e)))
                    continue
            rk_file = cfg_tmp.interaction_coeff_file

            if not self._check_l0_filter(tag, system_cfg):
                l0_val = 0.0
                try:
                    rk_df = pd.read_csv(rk_file)
                    l0_rows = rk_df[rk_df["term"] == "L0"]
                    if "pair" in rk_df.columns:
                        l0_rows = l0_rows[l0_rows["pair"].str.count("-") == 1]
                    l0_val = float(l0_rows["value_eV_per_formula"].max())
                except Exception:
                    pass
                filtered.append((name, l0_val))
                continue

            if cfg.skip_if_tables_exist and cfg_tmp.tables_ready():
                print(f"[{name}]  tables cached → analysis only")

            metals_str = "/".join(filter(None, [m1, m2, m3]))
            print(f"\n{'='*60}")
            print(f"  {name}  ({i+1}/{len(all_systems)})  [{metals_str}]")
            print(f"{'='*60}")

            try:
                analyzer = SystemAnalyzer(system_dir, config=system_cfg)
                ok = analyzer.prepare(metals)
                if not ok:
                    failed.append((name, "prepare failed"))
                    continue
                analyzer.run()
                self._compile_slice_video_grid(metals)
                done.append((name, m1, m2, m3))
            except Exception as e:
                print(f"\n  ERROR in {name}: {e}")
                traceback.print_exc()
                failed.append((name, str(e)))

        # Summary
        print(f"\n{'='*60}")
        print("BATCH COMPLETE")
        print(f"  Analysed: {len(done)}/{len(all_systems)}")
        print(f"  Filtered (L0>0): {len(filtered)}")
        print(f"  Failed:   {len(failed)}")
        print(f"  Skipped:  {len(skipped)}")
        if filtered:
            print("\nFiltered (miscibility gap):")
            for name, l0 in sorted(filtered, key=lambda x: -x[1]):
                print(f"  {name}  L0={l0:+.4f} eV")
        if failed:
            print("\nFailed:")
            for name, err in failed:
                print(f"  {name}: {err}")
        print(f"\nTables  → {cfg.tables_dir}")
        print(f"Figures → {cfg.figures_dir}")

        if failed:
            names = ", ".join(name for name, _ in failed)
            print(f"\nWarning: oxidation batch failed for {len(failed)} system(s): {names}")

        tags_info = [
            (
                system_tag(list(filter(None, [m1, m2, m3])), cfg.phase_element),
                m1,
                m2,
                system_key(list(filter(None, [m1, m2, m3])), cfg.phase_element),
            )
            for _, m1, m2, m3 in done
        ]
        self._compile_global_slice_video_grid(tags_info)
        return tags_info
