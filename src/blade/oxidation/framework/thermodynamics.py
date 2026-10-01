"""Core thermodynamic model for N-metal mixed-phase oxidation.

Supports binary and ternary mixed-phase systems.
Grand-potential LP minimization over an N-metal phase simplex
grid plus fixed line-compound phases.

Phase free energy (Muggianu model):
    G(y, T) = H_Muggianu(y) + k_B T * sum_i y_i * ln(y_i)

This module re-exports all public names from focused sibling modules:
  _grid        — simplex_grid_nd, grid_edges
  _mixing      — muggianu_energy_nd, ideal_mixing_nd
  _lp          — solve_grand_lp_batch, solve_grand_lp
  _oxygen      — mu_o_from_log10_po2, log10_po2_from_mu_o
  _phase_data  — load_phase_data
  _labels      — assemblage labeling and region-map plotting
"""

from __future__ import annotations

from ._grid import grid_edges, simplex_grid_nd
from ._labels import (
    _clean_phase_label,
    _format_exact_fraction_value,
    _phase_comp_label,
    _phase_comp_label_rounded,
    _phase_comp_range_label,
    _phase_component_signature,
    _phase_composition,
    _short_label,
    add_region_annotation,
    assign_region_ids,
    build_assemblage_labels,
    format_exact_phase_fraction_line,
    format_phase_detail_line,
    format_phase_fraction_comp,
    format_phase_fraction_summary_line,
    plot_region_map,
    region_annotation_text,
    separate_region_annotations,
    write_region_map_html,
)
from ._lp import solve_grand_lp, solve_grand_lp_batch
from ._mixing import ideal_mixing_nd, muggianu_energy_nd
from ._oxygen import (
    EV_KJ_PER_MOL,
    R,
    _oxygen_delta_g0,
    _oxygen_shomate_coeffs,
    log10_po2_from_mu_o,
    mu_o_from_log10_po2,
)
from ._phase_data import _load_nd, load_phase_data

# Module-level constant not consumed by any sibling module
KB = 8.617333262145e-5  # eV / K (Boltzmann constant)

__all__ = [
    # constants
    "KB",
    "R",
    "EV_KJ_PER_MOL",
    # _grid
    "simplex_grid_nd",
    "grid_edges",
    # _mixing
    "muggianu_energy_nd",
    "ideal_mixing_nd",
    # _lp
    "solve_grand_lp_batch",
    "solve_grand_lp",
    # _oxygen
    "_oxygen_shomate_coeffs",
    "_oxygen_delta_g0",
    "mu_o_from_log10_po2",
    "log10_po2_from_mu_o",
    # _phase_data
    "load_phase_data",
    "_load_nd",
    # _labels
    "_short_label",
    "_phase_composition",
    "_phase_comp_label_rounded",
    "_phase_comp_label",
    "_phase_component_signature",
    "_phase_comp_range_label",
    "format_phase_fraction_comp",
    "_format_exact_fraction_value",
    "_clean_phase_label",
    "format_phase_detail_line",
    "format_exact_phase_fraction_line",
    "format_phase_fraction_summary_line",
    "build_assemblage_labels",
    "assign_region_ids",
    "region_annotation_text",
    "add_region_annotation",
    "separate_region_annotations",
    "write_region_map_html",
    "plot_region_map",
]
