#!/usr/bin/env python3
"""Run the complete BLADE -> database -> oxidation workflow from TOML input."""

from __future__ import annotations

import argparse
import ast
import importlib.util
import itertools
import os
import re
import shutil
import sys
import tomllib
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

# Optional example-prototype lattice presets used by automatic averaging.
# Unlisted prototypes/elements use periodic-table references, and TOML
# phase.element_lattice_parameters entries may override or extend all presets.
_LATTICE_PARAMETER_LIBRARY: dict[str, dict[str, dict[str, float]]] = {
    "HEDB": {
        "Cr": {"a": 2.969, "c": 3.066},
        "Hf": {"a": 3.141, "c": 3.470},
        "Mo": {"a": 3.053, "c": 3.169},
        "Nb": {"a": 3.086, "c": 3.269},
        "Ta": {"a": 3.098, "c": 3.227},
        "Ti": {"a": 3.028, "c": 3.228},
        "V": {"a": 2.998, "c": 3.057},
        "W": {"a": 3.020, "c": 3.137},
        "Zr": {"a": 3.169, "c": 3.530},
    },
    "MAX": {
        "Cr": {"a": 2.86, "c": 12.80},
        "Hf": {"a": 3.28, "c": 14.36},
        "Mo": {"a": 2.96, "c": 13.20},
        "Nb": {"a": 3.10, "c": 13.80},
        "Ta": {"a": 3.08, "c": 14.10},
        "Ti": {"a": 3.04, "c": 13.60},
        "V": {"a": 2.91, "c": 13.20},
        "W": {"a": 2.91, "c": 13.80},
        "Zr": {"a": 3.32, "c": 14.57},
    },
}


@lru_cache(maxsize=1)
def _periodic_lattice_parameters() -> dict[str, dict[str, float]]:
    """Build reference a/b/c values for all 118 element symbols.

    ASE reference-state cells are used when available. Elements without a
    bulk reference cell use an isotropic diameter from ASE's covalent radius;
    this includes estimated values for short-lived synthetic elements.
    """
    from ase.data import chemical_symbols, covalent_radii, reference_states

    table = {}
    for atomic_number, symbol in enumerate(chemical_symbols[1:119], start=1):
        state = reference_states[atomic_number]
        if state and state.get("a") is not None:
            a = float(state["a"])
            b = a * float(state.get("b/a", 1.0))
            c = a * float(state.get("c/a", 1.0))
        else:
            a = b = c = 2.0 * float(covalent_radii[atomic_number])
        table[symbol] = {"a": a, "b": b, "c": c}
    return table


def _range(values: list[float], *, integer: bool = False) -> np.ndarray:
    if len(values) != 3:
        raise ValueError(f"range must be [start, stop, step], got {values!r}")
    start, stop, step = values
    if step <= 0:
        raise ValueError(f"range step must be positive, got {step}")
    dtype = int if integer else float
    return np.arange(start, stop + step * 0.001, step, dtype=dtype)


def _module_available(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ModuleNotFoundError, ValueError):
        return False


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _assignment_names(statement: ast.stmt) -> set[str]:
    targets = []
    if isinstance(statement, ast.Assign):
        targets = statement.targets
    elif isinstance(statement, ast.AnnAssign):
        targets = [statement.target]
    else:
        return set()
    names = set()
    for target in targets:
        if isinstance(target, ast.Name):
            names.add(target.id)
    return names


def _run_configured_script(path: Path, overrides: dict[str, Any]) -> None:
    """Execute a legacy example after replacing its top-level settings."""
    tree = ast.parse(path.read_text(), filename=str(path))
    tree.body = [statement for statement in tree.body if not (_assignment_names(statement) & overrides.keys())]
    ast.fix_missing_locations(tree)
    namespace = {
        "__name__": "__main__",
        "__file__": str(path),
        "__package__": None,
        **overrides,
    }
    exec(compile(tree, str(path), "exec"), namespace)


def _copy_nested_mapping(value: Any, name: str) -> dict:
    """Copy a TOML mapping while rejecting ambiguous non-table values."""
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a TOML table")
    return {str(key): dict(item) if isinstance(item, dict) else item for key, item in value.items()}


def _database_api_key(database: dict[str, Any]) -> str:
    """Resolve an MP key from an explicit value or configured environment.

    Older input files stored the literal key in api_key_env. Preserve that
    behavior when the named variable is absent, while preferring the explicit
    api_key field or a real environment variable.
    """
    explicit_key = str(database.get("api_key", "")).strip()
    if explicit_key:
        return explicit_key

    env_name = str(database.get("api_key_env", "MP_API_KEY")).strip()
    env_value = os.environ.get(env_name, "").strip()
    if env_value:
        return env_value
    if env_name and env_name != "MP_API_KEY":
        return env_name
    return ""


def _system_key(value: str | list[str]) -> str:
    """Return an order-independent key such as CrHf for a chemical system."""
    if isinstance(value, str):
        elements = re.findall(r"[A-Z][a-z]?", value)
        if "".join(elements) != value:
            raise ValueError(f"invalid system name {value!r}")
    else:
        elements = [str(element) for element in value]
    if not elements or len(set(elements)) != len(elements):
        raise ValueError(f"system must contain unique element symbols: {value!r}")
    return "".join(sorted(elements))


def _system_elements(value: str) -> list[str]:
    elements = re.findall(r"[A-Z][a-z]?", value)
    if "".join(elements) != value:
        raise ValueError(f"invalid system name {value!r}")
    return elements


class FullFrameworkPipeline:
    STAGE_ORDER = ("tdb", "phase_visualization", "database", "oxidation")

    def __init__(self, input_path: Path):
        self.input_path = input_path.resolve()
        with self.input_path.open("rb") as handle:
            self.settings = tomllib.load(handle)
        paths = self.settings["paths"]
        self.blade_root = Path(paths["blade_root"]).expanduser().resolve()
        self.files_dir = Path(paths.get("files_dir", self.blade_root / "Files")).expanduser().resolve()
        self.sqsdb_dir = Path(paths["sqsdb_dir"]).expanduser().resolve()
        self.scraps_repo = (
            Path(paths.get("scraps_repo", self.blade_root.parent / "SCRAPS" / "scraps-perpair")).expanduser().resolve()
        )
        source_root = str(self.blade_root / "src")
        if source_root not in sys.path:
            sys.path.insert(0, source_root)
        self.examples_dir = self.blade_root / "examples"
        self.oxidation_framework_dir = self.blade_root / "src" / "blade" / "oxidation" / "framework"

    def enabled(self, stage: str) -> bool:
        stages = self.settings.get("stages", {})
        aliases = {
            "tdb": ("sqs_generation", "tdb_fitting"),
            "phase_plots": ("phase_diagrams",),
            "phase_grid": ("phase_diagrams",),
            "phase_visualization": ("phase_diagrams",),
            "database": ("oxide_database",),
            "oxidation": ("oxidation_graphs",),
        }
        if stage in aliases and any(name in stages for name in aliases[stage]):
            return any(bool(stages[name]) for name in aliases[stage])
        if stage in stages:
            return bool(stages[stage])
        # Backward compatibility for pre-1.6 stage names.
        legacy = {
            "sqs_generation": "tdb",
            "tdb_fitting": "tdb",
            "phase_diagrams": "phase_plots",
            "phase_visualization": "phase_plots",
            "oxide_database": "database",
            "oxidation_graphs": "oxidation",
        }
        return bool(stages.get(legacy.get(stage, stage), False))

    def _sqs_method(self) -> str:
        method = str(self.settings.get("tdb", {}).get("sqs_method", "mcsqs")).lower()
        if method not in {"mcsqs", "scraps"}:
            raise ValueError("tdb.sqs_method must be 'mcsqs' or 'scraps'")
        return method

    def _tdb_driver(self) -> Path:
        mode = str(self.settings.get("tdb", {}).get("prototype_driver", "standard")).lower()
        drivers = {"standard": "hedb", "multibasis": "max"}
        if mode not in drivers:
            raise ValueError("tdb.prototype_driver must be 'standard' or 'multibasis'")
        # Legacy module names identify bundled implementations, not chemistry.
        family = drivers[mode]
        suffix = "_scraps" if self._sqs_method() == "scraps" else ""
        return self.examples_dir / "structures" / f"tdb_gen_{family}{suffix}.py"

    def _phase_plot_systems(self) -> set[str]:
        """Return explicit plot systems or reproduce the TDB composition set."""
        configured = self.settings.get("phase_plots", {}).get("systems")
        if configured is not None:
            if not isinstance(configured, list) or not configured:
                raise ValueError("phase_plots.systems must be a non-empty list")
            return {_system_key(str(system)) for system in configured}

        cfg = self.settings["tdb"]
        primary = list(cfg.get("primary_elements", []))
        secondary = list(cfg.get("secondary_elements", []))
        systems = set()
        for primary_count in range(cfg["primary_min"], cfg["primary_max"] + 1):
            for primary_group in itertools.combinations(primary, primary_count):
                for secondary_count in range(cfg.get("secondary_min", 0), cfg.get("secondary_max", 0) + 1):
                    for secondary_group in itertools.combinations(secondary, secondary_count):
                        elements = [*primary_group, *secondary_group]
                        if elements:
                            systems.add(_system_key(elements))
        return systems

    def _input_bundle(self, source: dict[str, Any], *, name: str) -> dict[str, Any]:
        """Normalize model-input and sublattice-site settings from TOML."""
        bundle = {
            key: _copy_nested_mapping(source.get(key), f"{name}.{key}")
            for key in (
                "terms_in",
                "mult_in",
                "sublattice_map",
                "sqsgen_in",
                "fixed_compositions",
            )
        }
        sites = _copy_nested_mapping(source.get("sites"), f"{name}.sites")
        for phase_name, phase_sites in sites.items():
            if not isinstance(phase_sites, dict):
                raise ValueError(f"{name}.sites.{phase_name} must be a TOML table")
            phase_map = bundle["sublattice_map"].setdefault(phase_name, {})
            if not isinstance(phase_map, dict):
                raise ValueError(f"{name}.sublattice_map.{phase_name} must be a table")
            constant = list(phase_map.get("Constant", []))
            for letter, site in phase_sites.items():
                if not isinstance(site, dict):
                    raise ValueError(f"{name}.sites.{phase_name}.{letter} must be a TOML table")
                elements = [str(element) for element in site.get("elements", [])]
                if not elements:
                    raise ValueError(f"{name}.sites.{phase_name}.{letter}.elements cannot be empty")
                phase_map[str(letter)] = elements
                if site.get("constant", False):
                    composition = [float(value) for value in site.get("composition", [])]
                    if not composition:
                        composition = [1.0 / len(elements)] * len(elements)
                    bundle["fixed_compositions"][str(letter)] = composition
                    if letter not in constant:
                        constant.append(str(letter))
            if constant:
                phase_map["Constant"] = constant

        # A fixed composition always means fixed element order. Add the
        # internal Constant marker wherever that sublattice is defined.
        for phase_map in bundle["sublattice_map"].values():
            if not isinstance(phase_map, dict):
                continue
            constant = list(phase_map.get("Constant", []))
            for letter in bundle["fixed_compositions"]:
                if letter in phase_map and letter not in constant:
                    constant.append(letter)
            if constant:
                phase_map["Constant"] = constant
        return bundle

    def _system_overrides(self) -> dict[str, dict[str, Any]]:
        overrides = {}
        for raw_name, source in self.settings.get("tdb_system_overrides", {}).items():
            if not isinstance(source, dict):
                raise ValueError(f"tdb_system_overrides.{raw_name} must be a table")
            key = _system_key(raw_name)
            if key in overrides:
                raise ValueError(f"duplicate system override for {key}")
            bundle = self._input_bundle(source, name=f"tdb_system_overrides.{raw_name}")
            if bundle["sqsgen_in"]:
                raise ValueError(
                    f"tdb_system_overrides.{raw_name}.sqsgen_in is ambiguous because "
                    "SQS templates are shared by system size; use sites or "
                    "fixed_compositions for system-specific fractions"
                )
            overrides[key] = bundle
        return overrides

    def _phase_lattice(self, cfg: dict, phase: dict, *, emit: bool = True) -> dict[str, float]:
        base = {
            "a": float(phase["a"]),
            "b": float(phase.get("b", phase["a"])),
            "c": float(phase["c"]),
        }
        if not phase.get("use_average_lattice", False):
            return base

        configured_table = phase.get("element_lattice_parameters", {})
        if not isinstance(configured_table, dict):
            raise ValueError("phase.element_lattice_parameters must be a TOML table")
        family = str(phase.get("lattice_parameter_family", phase.get("generator_name", ""))).strip().upper()
        parameter_table = {element: dict(values) for element, values in _periodic_lattice_parameters().items()}
        for element, values in _LATTICE_PARAMETER_LIBRARY.get(family, {}).items():
            family_values = dict(values)
            family_values.setdefault("b", family_values.get("a"))
            parameter_table.setdefault(element, {}).update(family_values)
        for element, values in configured_table.items():
            if not isinstance(values, dict):
                raise ValueError(f"phase.element_lattice_parameters.{element} must be a TOML table")
            parameter_table.setdefault(str(element), {}).update(values)
        average_spec = phase.get("average_elements")
        axis_sources = {"a": "a", "b": "b", "c": "c"}
        if average_spec is None:
            elements = []
        elif isinstance(average_spec, dict):
            elements = [str(element) for element in average_spec.get("elements", [])]
            for axis in axis_sources:
                source = str(average_spec.get(axis, axis)).strip().lower()
                if source not in axis_sources:
                    raise ValueError(f"phase.average_elements.{axis} must be 'a', 'b', or 'c'")
                axis_sources[axis] = source
        elif isinstance(average_spec, list):
            # Backward-compatible shorthand; dictionary form is preferred.
            elements = [str(element) for element in average_spec]
        else:
            raise ValueError("phase.average_elements must be a TOML inline table")
        if not elements:
            elements = list(
                dict.fromkeys(list(cfg.get("primary_elements", [])) + list(cfg.get("secondary_elements", [])))
            )
        if not elements:
            raise ValueError("phase.use_average_lattice requires at least one element")

        averaged = {}
        for axis in ("a", "b", "c"):
            source_axis = axis_sources[axis]
            values = []
            for element in elements:
                item = parameter_table.get(element)
                if item is None:
                    raise ValueError(f"Unknown element symbol in lattice average: {element}")
                values.append(float(item[source_axis]))
            averaged[axis] = float(np.mean(values))
        if emit:
            mapping = ", ".join(f"{axis}:{axis_sources[axis]}" for axis in ("a", "b", "c"))
            print(
                f"{family or 'custom'} lattice average from "
                + ", ".join(elements)
                + ": "
                + ", ".join(f"{axis}={averaged[axis]:.6g}" for axis in ("a", "b", "c"))
                + f"  ({mapping})"
            )
        return averaged

    @staticmethod
    def _validate_fixed_sites(bundle: dict[str, Any], name: str, errors: list[str]) -> None:
        phase_maps = bundle["sublattice_map"]
        for letter, composition in bundle["fixed_compositions"].items():
            try:
                values = [float(value) for value in composition]
            except (TypeError, ValueError):
                errors.append(f"{name}.fixed_compositions.{letter} must be numeric")
                continue
            if not np.isclose(sum(values), 1.0, atol=1e-8):
                errors.append(f"{name}.fixed_compositions.{letter} sums to " f"{sum(values):g}, not 1")
            for phase_name, phase_map in phase_maps.items():
                if not isinstance(phase_map, dict) or letter not in phase_map:
                    continue
                elements = phase_map[letter]
                if len(elements) != len(values):
                    errors.append(
                        f"{name}: {phase_name}.{letter} has {len(elements)} elements "
                        f"but its fixed composition has {len(values)} fractions"
                    )

    def validate(self, *, check_dependencies: bool = True) -> list[str]:
        errors = []
        if not self.blade_root.is_dir():
            errors.append(f"paths: BLADE root does not exist: {self.blade_root}")
        if self.enabled("tdb") and not self.sqsdb_dir.is_dir():
            errors.append(f"tdb: sqsdb directory does not exist: {self.sqsdb_dir}")
        required_scripts = {
            "tdb": self._tdb_driver(),
            "phase_visualization": self.examples_dir / "extra" / "phase_visualization.py",
            "database": self.examples_dir / "oxidation" / "database_framework.py",
        }
        for stage, script in required_scripts.items():
            if self.enabled(stage) and not script.is_file():
                errors.append(f"{stage}: missing script {script}")
        if self.enabled("oxidation") and not (self.oxidation_framework_dir / "__init__.py").is_file():
            errors.append(f"oxidation: missing bundled framework under {self.oxidation_framework_dir}")

        tdb = self.settings.get("tdb", {})
        levels = {entry.get("level") for entry in self.settings.get("tdb_sqs_levels", [])}
        if self.enabled("tdb") and tdb.get("level") not in levels:
            errors.append(f"tdb: selected level {tdb.get('level')} is absent from tdb_sqs_levels")
        if self.enabled("tdb"):
            try:
                inputs = self._input_bundle(self.settings.get("tdb_inputs", {}), name="tdb_inputs")
                system_overrides = self._system_overrides()
                self._phase_lattice(tdb, self.settings["phase"], emit=False)
                allowed = set(tdb.get("primary_elements", [])) | set(tdb.get("secondary_elements", []))
                for system_name, override in system_overrides.items():
                    unknown = set(_system_elements(system_name)) - allowed
                    if unknown:
                        errors.append(
                            f"tdb_system_overrides.{system_name}: elements are not in "
                            f"tdb element pools: {sorted(unknown)}"
                        )
                    self._validate_fixed_sites(override, f"tdb_system_overrides.{system_name}", errors)
                self._validate_fixed_sites(inputs, "tdb_inputs", errors)
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"tdb inputs: {exc}")
        phase_plots = self.settings.get("phase_plots", {})
        phase_grid = self.settings.get("phase_grid", {})
        if self.enabled("phase_plots"):
            try:
                if not self._phase_plot_systems():
                    errors.append("phase_plots: no systems selected")
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"phase_plots: {exc}")
        if self.enabled("phase_plots") and self.enabled("phase_grid"):
            plot_range = (
                phase_plots.get("t_start"),
                phase_plots.get("t_end"),
                phase_plots.get("gif_t_step"),
            )
            grid_range = (
                phase_grid.get("t_start"),
                phase_grid.get("t_end"),
                phase_grid.get("t_step"),
            )
            if plot_range != grid_range:
                errors.append(
                    f"phase_grid: temperature range {grid_range} does not match " f"phase GIF range {plot_range}"
                )
            if phase_plots.get("output_folder", "Phase_Diagrams") != phase_grid.get(
                "phase_diagrams_folder", "Phase_Diagrams"
            ):
                errors.append("phase_grid: input folder does not match phase_plots output folder")
        oxidation = self.settings.get("oxidation", {})
        phase_element = oxidation.get("phase_element", "").strip()
        stoichiometry = float(oxidation.get("phase_element_stoichiometry", 0.0))
        if self.enabled("phase_plots") and phase_element and stoichiometry > 0:
            expected_metal_fraction = 1.0 / (1.0 + stoichiometry)
            actual_metal_fraction = float(phase_plots.get("metal_fraction", 1.0))
            if not np.isclose(actual_metal_fraction, expected_metal_fraction):
                errors.append(
                    "phase_plots: metal_fraction does not match oxidation "
                    f"phase stoichiometry; expected {expected_metal_fraction:g}"
                )
        database = self.settings.get("database", {})
        fixed_elements = {str(item) for item in database.get("fixed_elements", [])}
        if (
            self.enabled("database")
            and self.enabled("oxidation")
            and phase_element
            and phase_element not in fixed_elements
        ):
            errors.append(
                f"oxidation: phase_element {phase_element!r} is absent from "
                f"database.fixed_elements {sorted(fixed_elements)}"
            )
        single = self.settings.get("oxidation_single", {})
        composition = single.get("composition", [])
        if composition and not np.isclose(sum(composition), 1.0, atol=1e-8):
            errors.append(f"oxidation_single: composition sums to {sum(composition):g}, not 1")
        oxidation_mode = str(oxidation.get("mode", "batch")).strip().lower()
        range_settings = {
            "batch": (
                ("oxidation_batch", "temperature_range"),
                ("oxidation_batch", "mu_o_range"),
                ("oxidation_batch", "scan_mu_o_range"),
                ("oxidation_batch", "map_x_range"),
            ),
            "single": (
                ("oxidation_single", "temperature_range"),
                ("oxidation_single", "mu_o_range"),
                ("oxidation_single", "x_range"),
            ),
        }
        if oxidation_mode not in range_settings:
            errors.append("oxidation.mode must be 'batch' or 'single'")
        for section, key in range_settings.get(oxidation_mode, ()):
            try:
                _range(self.settings[section][key])
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"{section}.{key}: {exc}")
        if not check_dependencies:
            return errors

        dependencies = {
            "tdb": ["blade", "pycalphad"],
            "phase_plots": ["blade", "pycalphad", "scipy", "matplotlib", "imageio"],
            "phase_grid": ["PIL", "matplotlib"],
            "database": ["blade", "pandas", "openpyxl"],
            "oxidation": ["blade", "scipy", "pandas", "matplotlib", "plotly"],
        }
        for stage, modules in dependencies.items():
            if not self.enabled(stage):
                continue
            for module in modules:
                if not _module_available(module):
                    errors.append(f"{stage}: Python dependency {module!r} is unavailable")
        if self.enabled("sqs_generation"):
            if self._sqs_method() == "mcsqs" and shutil.which("mcsqs") is None:
                errors.append("tdb: mcsqs is not on PATH")
            elif self._sqs_method() == "scraps":
                scraps = self.settings.get("tdb_scraps", {})
                scraps_bin = Path(scraps.get("bin", self.scraps_repo / "SCRAPs" / "scraps"))
                scraps_tools = Path(scraps.get("tools", self.scraps_repo / "tools"))
                if not scraps_bin.is_file():
                    errors.append(f"tdb: SCRAPS executable does not exist: {scraps_bin}")
                if not scraps_tools.is_dir():
                    errors.append(f"tdb: SCRAPS tools directory does not exist: {scraps_tools}")
                if not _module_available("ase"):
                    errors.append("tdb: SCRAPS requires Python dependency 'ase'")
        if self.enabled("phase_grid") and shutil.which("ffmpeg") is None:
            errors.append("phase_grid: ffmpeg is not on PATH")
        if self.enabled("database") and database.get("run_download", False):
            if not _module_available("mp_api"):
                errors.append("database: mp-api is unavailable")
            if not _database_api_key(database):
                errors.append(
                    "database: no Materials Project key; set database.api_key "
                    "or the environment variable named by database.api_key_env"
                )
        if self.enabled("database") and database.get("run_energies", False):
            if not _module_available("materialsframework"):
                errors.append("database: materialsframework is unavailable")
        return errors

    def describe(self) -> None:
        print(f"Input:       {self.input_path}")
        print(f"BLADE:       {self.blade_root}")
        print(f"Files:       {self.files_dir}")
        print(f"Oxidation:   {self.oxidation_framework_dir}")
        for stage in (
            "sqs_generation",
            "tdb_fitting",
            "phase_diagrams",
            "oxide_database",
            "oxidation_graphs",
        ):
            print(f"  {stage:12s} {'enabled' if self.enabled(stage) else 'disabled'}")

    def run_tdb(self) -> None:
        cfg = self.settings["tdb"]
        phase = self.settings["phase"]
        fit = dict(self.settings["tdb_fit"])
        fit["calculator"] = cfg["mlip"]
        fit["calculator_kwargs"] = dict(self.settings.get("tdb_mlip_kwargs", {}))
        sqs = dict(self.settings["tdb_sqs"])
        sqs["2"] = sqs.pop("pair_cutoff")
        sqs["3"] = sqs.pop("triplet_cutoff")
        sqs["4"] = sqs.pop("quadruplet_cutoff")
        phase_key = phase["key"]
        lattice = self._phase_lattice(cfg, phase)
        phases = {
            phase_key: {
                **lattice,
                "alpha": phase.get("alpha", 90),
                "beta": phase.get("beta", 90),
                "gamma": phase.get("gamma", 120),
                "vectors": phase["vectors"],
                "coords": phase["coords"],
            }
        }
        phase_list = [
            {
                "generator_name": phase["generator_name"],
                "lattice": phase_key,
                "supercell_size": tuple(phase["supercell_size"]),
            }
        ]
        inputs = self._input_bundle(self.settings.get("tdb_inputs", {}), name="tdb_inputs")
        terms = phase.get("terms", "").strip()
        if terms and phase_key not in inputs["terms_in"]:
            inputs["terms_in"][phase_key] = terms
        system_overrides = self._system_overrides()

        sqsgen_levels = [
            {
                **entry,
                "letter": list(entry["letter"]),
                "compositions": [list(values) for values in entry["compositions"]],
            }
            for entry in self.settings["tdb_sqs_levels"]
        ]
        selected_level = next(entry for entry in sqsgen_levels if entry["level"] == cfg["level"])
        for override in system_overrides.values():
            for letter, composition in override["fixed_compositions"].items():
                values = [float(value) for value in composition]
                if values not in selected_level["compositions"]:
                    selected_level["compositions"].append(values)
                if letter not in selected_level["letter"]:
                    selected_level["letter"].append(letter)
        scraps = self.settings.get("tdb_scraps", {})
        overrides = {
            "path0": self.blade_root.parent,
            "path1": self.blade_root,
            "path2": self.sqsdb_dir,
            "paths": [self.blade_root.parent, self.blade_root, self.sqsdb_dir],
            "level": cfg["level"],
            "run_sqs": self.enabled("sqs_generation"),
            "skip_existing_sqs": cfg["skip_existing_sqs"],
            "run_tdb": self.enabled("tdb_fitting"),
            "skip_existing_tdb": cfg["skip_existing_tdb"],
            "refit_existing_tdb": cfg.get("refit_existing_tdb", False),
            "skip_existing_plots": cfg.get("skip_existing_plots", False),
            "generate_gibbs_energy": self.enabled("tdb_fitting") and cfg.get("generate_gibbs_energy", True),
            "generate_gibbs_mixing": self.enabled("tdb_fitting") and cfg.get("generate_gibbs_mixing", True),
            "generate_phase_diagram": self.enabled("tdb_fitting") and cfg.get("generate_phase_diagram", True),
            "generate_combined_phase_diagram": self.enabled("tdb_fitting")
            and cfg.get("generate_combined_phase_diagram", True),
            "generate_contcar_plots": self.enabled("tdb_fitting") and cfg.get("generate_contcar_plots", True),
            "mlip": cfg["mlip"],
            "mlip_kwargs": dict(self.settings.get("tdb_mlip_kwargs", {})),
            "tdb_params": fit,
            "terms_in": inputs["terms_in"] or None,
            "mult_in": inputs["mult_in"] or None,
            "sublattice_map": inputs["sublattice_map"] or None,
            "sqsgen_in": inputs["sqsgen_in"] or None,
            "fixed_compositions": inputs["fixed_compositions"] or None,
            "system_overrides": system_overrides or None,
            "run_movie": self.enabled("tdb_fitting") and cfg.get("run_movie", False),
            "primary_elements": list(cfg["primary_elements"]),
            "secondary_elements": list(cfg.get("secondary_elements", [])),
            "primary_min": cfg["primary_min"],
            "primary_max": cfg["primary_max"],
            "secondary_min": cfg.get("secondary_min", 0),
            "secondary_max": cfg.get("secondary_max", 0),
            "phases": phases,
            "phase_list": phase_list,
            "liquid": cfg.get("liquid", False),
            "sqsgen_levels": sqsgen_levels,
            "mcsqs_params": sqs,
            "SCRAPS_REPO": self.scraps_repo,
            "SCRAPS_BIN": Path(scraps.get("bin", self.scraps_repo / "SCRAPs" / "scraps")),
            "SCRAPS_TOOLS": Path(scraps.get("tools", self.scraps_repo / "tools")),
            "scraps_ranks": int(scraps.get("ranks", 10)),
            "auto_budget": int(scraps.get("auto_budget", 3)),
            "max_shellnum": int(scraps.get("max_shellnum", 4)),
            "fix_multibasis_sublattice": bool(
                scraps.get(
                    "fix_multibasis_sublattice",
                    str(phase.get("generator_name", "")).upper() == "MAX",
                )
            ),
            "comps_dir": self.files_dir / cfg.get("comps_folder", "Comps"),
        }
        _run_configured_script(self._tdb_driver(), overrides)

    def run_phase_visualization(self) -> None:
        plot_cfg = self.settings["phase_plots"]
        grid_cfg = self.settings["phase_grid"]
        overrides = {
            "run_phase_plots": self.enabled("phase_plots"),
            "run_phase_grid": self.enabled("phase_grid"),
            "path0": self.blade_root.parent,
            "path1": self.blade_root,
            "comps_dir": self.files_dir / plot_cfg.get("comps_folder", "Comps"),
            "output_dir": self.files_dir / plot_cfg.get("output_folder", "Phase_Diagrams"),
            "gif_dir": self.files_dir / grid_cfg.get("phase_diagrams_folder", "Phase_Diagrams"),
            "t_start": plot_cfg["t_start"],
            "t_end": plot_cfg["t_end"],
            "t_step": plot_cfg["t_step"],
            "skip": plot_cfg["skip_missing_calculations"],
            "make_gif": plot_cfg["make_gif"],
            "gif_t_step": plot_cfg["gif_t_step"],
            "gif_fps": plot_cfg["gif_fps"],
            "fixed_elements": set(plot_cfg.get("fixed_elements", [])),
            "metal_fraction": plot_cfg["metal_fraction"],
            "n_grid": plot_cfg["n_grid"],
            "two_phase_threshold": plot_cfg["two_phase_threshold"],
            "pressure": plot_cfg["pressure"],
            "panel_w": grid_cfg["panel_width"],
            "panel_h": grid_cfg["panel_height"],
            "n_cols": grid_cfg.get("columns"),
            "fps": grid_cfg["fps"],
            "bg_color": tuple(grid_cfg["background_rgb"]),
            "header_h": grid_cfg["header_height"],
            "cbar_w": grid_cfg["colorbar_width"],
            "t_start_gif": grid_cfg["t_start"],
            "selected_systems": self._phase_plot_systems(),
            "t_end_gif": grid_cfg["t_end"],
            "grid_t_step": grid_cfg["t_step"],
        }
        _run_configured_script(self.examples_dir / "extra" / "phase_visualization.py", overrides)

    def run_database(self) -> None:
        cfg = dict(self.settings["database"])
        api_key = _database_api_key(cfg)
        cfg.pop("api_key", None)
        cfg.pop("api_key_env", None)
        if api_key:
            os.environ["MP_API_KEY"] = api_key
        cfg["files_dir"] = self.files_dir
        cfg["fixed_elements"] = frozenset(cfg.get("fixed_elements", []))
        cfg["mlip_kwargs"] = dict(self.settings.get("database_mlip_kwargs", {}))
        cfg["fallback_refs"] = dict(self.settings.get("database_fallback_refs", {}))
        module = _load_module(
            self.examples_dir / "oxidation" / "database_framework.py",
            "blade_database_framework",
        )
        module.OxidationPipeline(module.PipelineConfig(**cfg)).run()

    def run_oxidation(self) -> None:
        from blade.oxidation import OxidationCalculator, OxidationConfig

        cfg = self.settings["oxidation"]
        phase_element = cfg.get("phase_element", "").strip() or None
        oxidation_cfg = OxidationConfig(
            files_dir=self.files_dir,
            phase_element=phase_element,
            phase_element_stoichiometry=cfg.get("phase_element_stoichiometry", 0.0),
            mixed_phase_subdir=cfg.get("mixed_phase_subdir", "blade"),
            fixed_phases_subdir=cfg.get("fixed_phases_subdir", "ORB"),
            region_label_mode=cfg.get("region_label_mode", "phases"),
            region_label_fontsize=cfg.get("region_label_fontsize", 7),
            slice_axis_priority=list(cfg.get("slice_axis_priority", [])),
        )
        calculator = OxidationCalculator(oxidation_cfg)
        mode = cfg.get("mode", "batch").strip().lower()
        if mode == "batch":
            batch = dict(self.settings["oxidation_batch"])
            batch["temperature_values"] = _range(batch.pop("temperature_range"), integer=True)
            batch["mu_O_values"] = _range(batch.pop("mu_o_range"))
            batch["scan_mu_O"] = _range(batch.pop("scan_mu_o_range"))
            batch["map_x_x_values"] = _range(batch.pop("map_x_range"))
            calculator.run_batch(**batch)
        elif mode == "single":
            single = dict(self.settings["oxidation_single"])
            system = single.pop("system")
            metals = list(single.pop("metals"))
            composition = list(single.pop("composition"))
            temperatures = _range(single.pop("temperature_range"), integer=True)
            mu_o = _range(single.pop("mu_o_range"))
            scan_t = list(single.pop("scan_temperatures"))
            x_values = _range(single.pop("x_range"))
            skip = bool(single.pop("skip_if_exists", True))
            run_calculations = bool(single.pop("run_calculations", True))
            run_plots = bool(single.pop("run_plots", True))
            analyzer = calculator.single(system=system, metals=metals, composition=composition, **single)
            analyzer.run(
                T_values=temperatures,
                mu_O_values=mu_o,
                scan_T=scan_t,
                x_values=x_values,
                skip_if_exists=skip,
                run_calculations=run_calculations,
                run_plots=run_plots,
            )
        else:
            raise ValueError("oxidation.mode must be 'batch' or 'single'")

    def run(self, *, dry_run: bool = False) -> None:
        self.describe()
        methods = {
            "tdb": self.run_tdb,
            "phase_visualization": self.run_phase_visualization,
            "database": self.run_database,
            "oxidation": self.run_oxidation,
        }
        for stage in self.STAGE_ORDER:
            if not self.enabled(stage):
                continue
            print(f"\n=== {stage} ===")
            if dry_run:
                print("dry-run: not executed")
            else:
                methods[stage]()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=Path(__file__).with_name("full_framework.toml"))
    parser.add_argument("--check", action="store_true", help="validate configuration and dependencies")
    parser.add_argument("--dry-run", action="store_true", help="show enabled stages without running them")
    oxidation_mode = parser.add_mutually_exclusive_group()
    oxidation_mode.add_argument("--oxidation-plot-only", action="store_true", help=argparse.SUPPRESS)
    oxidation_mode.add_argument("--oxidation-onset-only", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    pipeline = FullFrameworkPipeline(args.input)
    if args.oxidation_plot_only:
        pipeline.settings["stages"] = {"oxidation": True}
        batch = pipeline.settings.setdefault("oxidation_batch", {})
        batch["run_calculations"] = False
        batch["run_plots"] = True
    elif args.oxidation_onset_only:
        pipeline.settings["stages"] = {"oxidation": True}
        batch = pipeline.settings.setdefault("oxidation_batch", {})
        batch.update(
            {
                "run_calculations": True,
                "run_plots": True,
                "run_scan": False,
                "run_muO_x_map": False,
                "run_onset": True,
                "run_3d_onset": False,
                "run_animations": False,
                "run_composition_slices": False,
                "run_composition_slice_muT": False,
                "calculation_first": False,
            }
        )
    errors = pipeline.validate(check_dependencies=args.check)
    if errors:
        print("Configuration errors:")
        for error in errors:
            print(f"  - {error}")
        return 2
    if args.check:
        pipeline.describe()
        print("Configuration and dependencies are valid.")
        return 0
    try:
        pipeline.run(dry_run=args.dry_run)
    except RuntimeError as error:
        if not args.oxidation_plot_only:
            raise
        print(f"\nPartial oxidation plotting completed with unavailable caches: {error}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
