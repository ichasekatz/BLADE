"""BLADE interface to the generalized grand-potential oxidation backend."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

_default_region_label_fontsize: int = 7


@dataclass
class OxidationConfig:
    """Paths and chemistry shared by oxidation calculations and plots.

    Attributes:
        files_dir: Root output directory for all BLADE oxidation outputs.
        framework_root: Path to the bundled oxidation framework package.
            ``None`` uses the framework that lives beside this module.
        phase_element: Element symbol whose stoichiometry is fixed across the
            composition sweep (e.g. ``"B"`` for borides).
        phase_element_stoichiometry: Atoms-per-formula-unit for
            *phase_element*; must be positive when *phase_element* is set.
        mixed_phase_subdir: Sub-directory name for BLADE/MLIP-relaxed phases.
        fixed_phases_subdir: Sub-directory name for fixed reference phases.
        region_label_mode: How to label stability regions — ``"id"`` for
            numeric IDs, ``"phases"`` for phase formulae.
        region_label_fontsize: Font size (pt) for region labels in plots.
        slice_axis: Composition axis to fix when projecting high-dimensional
            diagrams to 2-D.
        slice_axis_priority: Ordered list of axis candidates when
            *slice_axis* is ``None``.
    """

    files_dir: Path
    # None uses the framework bundled beside this module in BLADE.
    framework_root: Path | None = None
    phase_element: str | None = None
    phase_element_stoichiometry: float = 0.0
    mixed_phase_subdir: str = "blade"
    fixed_phases_subdir: str = "ORB"
    region_label_mode: str = "phases"
    region_label_fontsize: int = _default_region_label_fontsize
    slice_axis: int | str | None = None
    slice_axis_priority: list[str] | None = None

    def __post_init__(self) -> None:
        """Coerce, validate, and normalise all config fields after construction.

        Raises:
            ValueError: If *phase_element* is set but
                *phase_element_stoichiometry* is not positive.
            ValueError: If *region_label_mode* is not ``"id"`` or
                ``"phases"``.
        """
        self.files_dir = Path(self.files_dir).expanduser().resolve()
        self.framework_root = Path(self.framework_root).expanduser().resolve() if self.framework_root is not None else None
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
        """Absolute path to the ``system_structures`` sub-directory."""
        return self.files_dir / "system_structures"

    @property
    def phase_diagrams_dir(self) -> Path:
        """Absolute path to the ``Phase_Diagrams`` sub-directory."""
        return self.files_dir / "Phase_Diagrams"

    @property
    def tables_dir(self) -> Path:
        """Absolute path to the ``oxidation/tables`` sub-directory."""
        return self.files_dir / "oxidation" / "tables"

    @property
    def outputs_dir(self) -> Path:
        """Absolute path to the ``oxidation/figures`` sub-directory."""
        return self.files_dir / "oxidation" / "figures"

    def load_backend(self) -> None:
        """Inject the bundled framework into ``sys.path`` and create output dirs.

        Inserts *framework_root* (or the directory containing this module) at
        the front of :data:`sys.path` so that ``import framework`` resolves to
        the bundled oxidation backend.  Also creates :attr:`tables_dir` and
        :attr:`outputs_dir`.

        Raises:
            FileNotFoundError: If no valid ``framework/__init__.py`` is found
                under the resolved backend root.
        """
        backend_root = self.framework_root or Path(__file__).resolve().parent
        if not (backend_root / "framework" / "__init__.py").is_file():
            raise FileNotFoundError(f"No bundled oxidation framework found under {backend_root}")
        root = str(backend_root)
        if root not in sys.path:
            sys.path.insert(0, root)
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.outputs_dir.mkdir(parents=True, exist_ok=True)


class OxidationCalculator:
    """Configure and run cached batch or fixed-composition equilibrium solves.

    Wraps the bundled ``framework.BatchRunner`` and
    ``framework.SingleCompositionAnalyzer`` with BLADE path conventions so
    callers only need to supply an :class:`OxidationConfig`.

    Args:
        config: Fully initialised configuration object.  ``load_backend()``
            is called automatically during construction.
    """

    def __init__(self, config: OxidationConfig) -> None:
        """Store config and load the framework backend.

        Args:
            config: Oxidation configuration; :meth:`OxidationConfig.load_backend`
                is invoked immediately.
        """
        self.config = config
        self.config.load_backend()

    def batch_config(self, **overrides: Any):
        """Build a ``BatchConfig`` populated with BLADE path defaults.

        Keyword arguments in *overrides* are merged on top of the defaults
        derived from :attr:`config`, allowing per-call customisation without
        mutating the shared config.

        Args:
            **overrides: Any ``BatchConfig`` field to override.

        Returns:
            A fully configured ``BatchConfig`` instance ready to pass to
            :class:`framework.BatchRunner`.
        """
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
        """Run all enabled analyses, reusing compatible table and grid caches.

        Constructs a :class:`framework.BatchRunner` from :meth:`batch_config`
        and delegates to its ``run()`` method.

        Args:
            **overrides: Forwarded verbatim to :meth:`batch_config`.

        Returns:
            List of result objects returned by the framework runner.
        """
        from framework import BatchRunner

        return BatchRunner(self.batch_config(**overrides)).run()

    def single(
        self,
        system: str,
        composition: float | list[float],
        metals: list[str] | None = None,
        **overrides: Any,
    ):
        """Create a ``SingleCompositionAnalyzer`` for one BLADE system.

        Args:
            system: Name of the system sub-directory under
                :attr:`~OxidationConfig.structures_dir`
                (e.g. ``"HfCrB2-O"``).
            composition: Scalar or list of metal fractions passed directly
                to ``SingleCompositionAnalyzer``.
            metals: Ordered list of metal element symbols.  ``None`` lets
                the analyzer infer them from the structure files.
            **overrides: Extra keyword arguments forwarded to
                ``SingleCompositionAnalyzer``.

        Returns:
            A configured ``SingleCompositionAnalyzer`` instance ready to
            call ``.run()`` on.
        """
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
