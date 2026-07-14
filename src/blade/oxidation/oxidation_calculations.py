"""BLADE interface to the generalized grand-potential oxidation backend."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class OxidationConfig:
    """Paths and chemistry shared by oxidation calculations and plots."""

    files_dir: Path
    # None uses the framework bundled beside this module in BLADE.
    framework_root: Path | None = None
    phase_element: str | None = None
    phase_element_stoichiometry: float = 0.0
    mixed_phase_subdir: str = "blade"
    fixed_phases_subdir: str = "ORB"
    region_label_mode: str = "phases"
    region_label_fontsize: int = 7
    slice_axis: int | str | None = None
    slice_axis_priority: list[str] | None = None

    def __post_init__(self) -> None:
        self.files_dir = Path(self.files_dir).expanduser().resolve()
        self.framework_root = (
            Path(self.framework_root).expanduser().resolve() if self.framework_root is not None else None
        )
        self.phase_element = str(self.phase_element).strip() if self.phase_element else None
        self.phase_element_stoichiometry = float(self.phase_element_stoichiometry) if self.phase_element else 0.0
        if self.phase_element and self.phase_element_stoichiometry <= 0:
            raise ValueError("phase_element_stoichiometry must be positive when phase_element is set")
        mode = str(self.region_label_mode).strip().lower()
        if mode not in {"id", "phases"}:
            raise ValueError("region_label_mode must be 'id' or 'phases'")
        self.region_label_mode = mode
        self.region_label_fontsize = int(self.region_label_fontsize)

    @property
    def structures_dir(self) -> Path:
        return self.files_dir / "system_structures"

    @property
    def phase_diagrams_dir(self) -> Path:
        return self.files_dir / "Phase_Diagrams"

    @property
    def tables_dir(self) -> Path:
        return self.files_dir / "oxidation" / "tables"

    @property
    def outputs_dir(self) -> Path:
        return self.files_dir / "oxidation" / "figures"

    def load_backend(self) -> None:
        backend_root = self.framework_root or Path(__file__).resolve().parent
        if not (backend_root / "framework" / "__init__.py").is_file():
            raise FileNotFoundError(f"No bundled oxidation framework found under {backend_root}")
        root = str(backend_root)
        if root not in sys.path:
            sys.path.insert(0, root)
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.outputs_dir.mkdir(parents=True, exist_ok=True)


class OxidationCalculator:
    """Configure and run cached batch or fixed-composition equilibrium solves."""

    def __init__(self, config: OxidationConfig):
        self.config = config
        self.config.load_backend()

    def batch_config(self, **overrides: Any):
        """Build the bundled oxidation BatchConfig with BLADE path defaults."""
        from framework import BatchConfig

        settings: dict[str, Any] = {
            "phase_element": self.config.phase_element,
            "phase_element_stoichiometry": self.config.phase_element_stoichiometry,
            "system_structures_root": self.config.structures_dir,
            "tables_dir": self.config.tables_dir,
            "figures_dir": self.config.outputs_dir,
            "phase_diagrams_root": self.config.phase_diagrams_dir,
            "mixed_phase_subdir": self.config.mixed_phase_subdir,
            "fixed_phases_subdir": self.config.fixed_phases_subdir,
            "skip_if_tables_exist": True,
            "skip_if_analysis_exists": True,
            "region_label_mode": self.config.region_label_mode,
            "region_label_fontsize": self.config.region_label_fontsize,
            "slice_axis": self.config.slice_axis,
            "slice_axis_priority": self.config.slice_axis_priority or [],
        }
        settings.update(overrides)
        return BatchConfig(**settings)

    def run_batch(self, **overrides: Any) -> list:
        """Run enabled analyses, reusing every compatible table and grid cache."""
        from framework import BatchRunner

        return BatchRunner(self.batch_config(**overrides)).run()

    def single(
        self,
        system: str,
        composition: float | list[float],
        metals: list[str] | None = None,
        **overrides: Any,
    ):
        """Create a configured SingleCompositionAnalyzer for one BLADE system."""
        from framework import SingleCompositionAnalyzer

        settings: dict[str, Any] = {
            "system_dir": self.config.structures_dir / system,
            "tables_dir": self.config.tables_dir,
            "figures_dir": self.config.outputs_dir,
            "metals": metals,
            "composition": composition,
            "mixed_phase_subdir": self.config.mixed_phase_subdir,
            "fixed_phases_subdir": self.config.fixed_phases_subdir,
            "phase_element": self.config.phase_element,
            "phase_element_stoichiometry": self.config.phase_element_stoichiometry,
            "region_label_mode": self.config.region_label_mode,
            "region_label_fontsize": self.config.region_label_fontsize,
        }
        settings.update(overrides)
        return SingleCompositionAnalyzer(**settings)
