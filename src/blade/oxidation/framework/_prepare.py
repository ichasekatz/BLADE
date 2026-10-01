"""Preparation helpers for SingleCompositionAnalyzer.

Contains path/environment setup, metal auto-detection, table preparation,
and phase-data loading — extracted from single.py to keep that module slim.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from blade.oxidation.framework.single import SingleCompositionAnalyzer


def _ensure_path(self) -> None:
    """Insert the framework root and its python/ sub-directory onto sys.path.

    Args:
        self: SingleCompositionAnalyzer instance.
    """
    root = Path(__file__).parent
    for p in [str(root), str(root / "python")]:
        if p not in sys.path:
            sys.path.insert(0, p)


def _detect_metals(self) -> list[str]:
    """Auto-detect metal elements from the mixed-phase subdirectory structure.

    Args:
        self: SingleCompositionAnalyzer instance.

    Returns:
        Ordered list of metal element symbols found in the directory tree.

    Raises:
        ValueError: If system_dir or the mixed_phase_subdir is not available.
    """
    _ensure_path(self)
    from .system_config import SystemConfig

    if self.system_dir and (self.system_dir / self.mixed_phase_subdir).is_dir():
        return SystemConfig._detect_metals_from_dirs(
            self.system_dir / self.mixed_phase_subdir, {self.phase_element} if self.phase_element else set()
        )
    raise ValueError("system_dir with mixed_phase_subdir required for auto-detection")


def prepare(self) -> SingleCompositionAnalyzer:
    """Prepare equilibrium tables; skips if they already exist.

    Args:
        self: SingleCompositionAnalyzer instance.

    Returns:
        self, for method chaining.
    """
    _ensure_path(self)
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
    """Load equilibrium phase data from the prepared tables.

    Args:
        self: SingleCompositionAnalyzer instance.

    Returns:
        Phase-data dictionary as returned by thermodynamics.load_phase_data.
    """
    _ensure_path(self)
    from .thermodynamics import load_phase_data

    if self._sys_cfg is None:
        prepare(self)
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
