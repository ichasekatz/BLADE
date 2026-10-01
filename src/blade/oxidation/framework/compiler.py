"""Compiler — post-processing: combined figures, onset plots, miscibility gaps."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .combined_figures import compile_combined_figures
from .miscibility import plot_miscibility_gaps
from .onset_plots import compile_onset_lines, compile_ternary_onset_grid

if TYPE_CHECKING:
    from .config import Config


class Compiler:
    """Compile post-processing figures after BatchRunner finishes.

    Usage::

        from .compiler import Compiler; from .config import Config
        cfg = Config(...)
        comp = Compiler(cfg)
        comp.compile(tags_info)   # tags_info from BatchRunner.run()
    """

    def __init__(self, config: Config):
        self.config = config

    def compile(self, tags_info: list) -> None:
        """Run all compilation steps.

        Args:
            tags_info: List of (tag, m1, m2, sys_name) tuples from BatchRunner.
        """
        self.compile_onset_lines(tags_info)
        self.compile_ternary_onset_grid()
        self.plot_miscibility_gaps()
        self.compile_combined_figures()
        print("Post-processing done.")

    def compile_onset_lines(self, tags_info: list) -> None:
        """Copy assemblage maps + build combined onset line plots.

        Args:
            tags_info: List of (tag, m1, m2, sys_name) tuples from BatchRunner.
        """
        compile_onset_lines(self.config, tags_info)

    def compile_ternary_onset_grid(self) -> None:
        """Build a grid of ternary onset plots across all systems."""
        compile_ternary_onset_grid(self.config)

    def plot_miscibility_gaps(self) -> None:
        """Compute and plot miscibility gaps for all ternary TDB systems."""
        plot_miscibility_gaps(self.config)

    def compile_combined_figures(self) -> None:
        """Build per-system three-panel combined figures and an overview grid."""
        compile_combined_figures(self.config)
