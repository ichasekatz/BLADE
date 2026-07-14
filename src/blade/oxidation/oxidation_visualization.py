"""Cached visualization entry points for BLADE oxidation analyses."""

from __future__ import annotations

from typing import Any

from .oxidation_calculations import OxidationCalculator, OxidationConfig


class OxidationVisualizer:
    """Regenerate oxidation plots from compatible BLADE/Files caches."""

    def __init__(self, config: OxidationConfig):
        self.config = config
        self.calculator = OxidationCalculator(config)

    def render_cached_batch(self, **overrides: Any) -> list:
        """Render enabled maps while requiring existing prepared system tables."""
        if not any(self.config.tables_dir.glob("*/")):
            raise FileNotFoundError(f"No prepared oxidation tables exist in {self.config.tables_dir}")
        settings: dict[str, Any] = {
            "skip_if_tables_exist": True,
            "skip_if_analysis_exists": True,
            "region_label_mode": self.config.region_label_mode,
            "region_label_fontsize": self.config.region_label_fontsize,
        }
        settings.update(overrides)
        return self.calculator.run_batch(**settings)

    def set_region_labels(self, mode: str, fontsize: int | None = None) -> None:
        """Switch future region maps between numeric IDs and formatted phases."""
        normalized = str(mode).strip().lower()
        if normalized not in {"id", "phases"}:
            raise ValueError("mode must be 'id' or 'phases'")
        self.config.region_label_mode = normalized
        if fontsize is not None:
            if int(fontsize) <= 0:
                raise ValueError("fontsize must be positive")
            self.config.region_label_fontsize = int(fontsize)
