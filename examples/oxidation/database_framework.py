"""Run oxidation examples 01-04 as one resumable BLADE workflow.

Set the Materials Project key before running the download stage:

    export MP_API_KEY="your-key"
    python examples/oxidation/00_full_pipeline.py

The original 01-04 examples remain useful for running individual stages.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from blade.oxidation.compositions import OxideCompositions
from blade.oxidation.database import OxideDatabase


@dataclass
class PipelineConfig:
    files_dir: Path = Path("/Users/chasekatz/Desktop/School/Research/BLADE/Files")
    primary_elements: list[str] = field(default_factory=lambda: ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"])
    secondary_elements: list[str] = field(default_factory=list)
    fixed_elements: frozenset[str] = frozenset({"B"})
    oxygen_element: str = "O"
    primary_min: int = 0
    primary_max: int = 3
    secondary_min: int = 0
    secondary_max: int = 0

    stable_only: bool = True
    skip_blade_formulas: bool = True
    skip_existing_poscars: bool = True
    skip_existing_energies: bool = True

    mlip: str = "orb"
    mlip_kwargs: dict = field(default_factory=lambda: {"steps": 1000, "device": "cpu"})
    mlip_label: str = "MP+ORB"
    mlip_folder: str = "MaterialsProject_Comps_ORB"
    poscar_folder: str = "MaterialsProject_Comps_POSCARs"
    blade_comp_folders: list[str] = field(default_factory=lambda: ["Comps"])

    run_compositions: bool = True
    run_download: bool = True
    run_energies: bool = True
    run_database: bool = True

    fallback_refs: dict[str, float] = field(
        default_factory=lambda: {
            "C": -9.200,
            "Cr": -9.632,
            "Hf": -9.956,
            "Mo": -10.850,
            "Nb": -10.094,
            "O": -4.949,
            "Ta": -11.853,
            "Ti": -7.897,
            "V": -9.080,
            "W": -12.960,
            "Zr": -8.547,
            "B": -6.680,
        }
    )

    def __post_init__(self) -> None:
        self.files_dir = Path(self.files_dir).expanduser().resolve()
        self.files_dir.mkdir(parents=True, exist_ok=True)


class OxidationPipeline:
    """Orchestrate composition, MP, MLIP, and database stages in order."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    @property
    def poscar_dir(self) -> Path:
        return self.config.files_dir / self.config.poscar_folder

    @property
    def mlip_dir(self) -> Path:
        return self.config.files_dir / self.config.mlip_folder

    def generate_compositions(self) -> None:
        cfg = self.config
        OxideCompositions(
            files_dir=cfg.files_dir,
            primary_elements=cfg.primary_elements,
            secondary_elements=cfg.secondary_elements,
            oxygen_element=cfg.oxygen_element,
            primary_min=cfg.primary_min,
            primary_max=cfg.primary_max,
            secondary_min=cfg.secondary_min,
            secondary_max=cfg.secondary_max,
            include_no_oxygen=True,
            include_fixed_oxygen=True,
            include_added_oxygen=True,
            include_fixed=True,
            fixed_elements=cfg.fixed_elements,
        ).run()

    def _blade_formulas(self) -> set[str]:
        if not self.config.skip_blade_formulas:
            return set()
        blade_xlsx = self.config.files_dir / "blade_generated_data.xlsx"
        if not blade_xlsx.exists():
            return set()
        from pymatgen.core import Composition

        formulas = set()
        frame = pd.read_excel(blade_xlsx)
        for formula in frame.get("formula", pd.Series(dtype=str)).dropna():
            try:
                formulas.add(Composition(str(formula)).reduced_formula)
            except Exception:
                continue
        return formulas

    def download_poscars(self) -> None:
        from mp_api.client import MPRester
        from pymatgen.core import Composition
        from pymatgen.io.vasp import Poscar

        cfg = self.config
        api_key = os.environ.get("MP_API_KEY", "").strip()
        if not api_key:
            raise RuntimeError("Set MP_API_KEY before running the Materials Project stage")

        composition_xlsx = cfg.files_dir / "composition_list.xlsx"
        if not composition_xlsx.exists():
            raise FileNotFoundError(f"Missing composition table: {composition_xlsx}")
        self.poscar_dir.mkdir(parents=True, exist_ok=True)
        blade_formulas = self._blade_formulas()

        frame = pd.read_excel(composition_xlsx)
        systems = []
        for _, row in frame.iterrows():
            elements = [item.strip() for item in str(row["elements"]).split(",") if item.strip()]
            if elements:
                systems.append(
                    {
                        "composition": row["composition"],
                        "metal_composition": row.get("metal_composition", ""),
                        "both_composition": row.get("both_composition", ""),
                        "chemsys": "-".join(sorted(elements)),
                    }
                )
        systems_df = pd.DataFrame(systems).drop_duplicates(subset=["chemsys"])
        fields = [
            "material_id",
            "formula_pretty",
            "energy_per_atom",
            "formation_energy_per_atom",
            "energy_above_hull",
            "is_stable",
            "volume",
            "density",
            "nsites",
            "symmetry",
            "structure",
        ]
        summary_rows = []
        with MPRester(api_key) as mpr:
            for _, system in systems_df.iterrows():
                chemsys = system["chemsys"]
                print(f"Searching Materials Project: {chemsys}")
                try:
                    docs = mpr.materials.summary.search(
                        chemsys=chemsys,
                        is_stable=True if cfg.stable_only else None,
                        fields=fields,
                    )
                except Exception as exc:
                    print(f"  Failed {chemsys}: {exc}")
                    continue
                for doc in docs:
                    reduced = Composition(doc.formula_pretty).reduced_formula
                    if reduced in blade_formulas:
                        continue
                    material_dir = self.poscar_dir / chemsys / f"{doc.material_id}_{doc.formula_pretty}"
                    poscar_path = material_dir / "POSCAR"
                    material_dir.mkdir(parents=True, exist_ok=True)
                    if not (cfg.skip_existing_poscars and poscar_path.exists()):
                        Poscar(doc.structure).write_file(poscar_path)
                    symmetry = doc.symmetry
                    summary_rows.append(
                        {
                            "input_composition": system["composition"],
                            "input_metal_composition": system["metal_composition"],
                            "input_both_composition": system["both_composition"],
                            "chemsys": chemsys,
                            "material_id": str(doc.material_id),
                            "formula": doc.formula_pretty,
                            "energy_per_atom": doc.energy_per_atom,
                            "formation_energy_per_atom": doc.formation_energy_per_atom,
                            "energy_above_hull": doc.energy_above_hull,
                            "is_stable": doc.is_stable,
                            "volume": doc.volume,
                            "volume_per_atom": doc.volume / doc.nsites if doc.nsites else None,
                            "density": doc.density,
                            "nsites": doc.nsites,
                            "crystal_system": getattr(symmetry, "crystal_system", None),
                            "spacegroup_number": getattr(symmetry, "number", None),
                            "spacegroup_symbol": getattr(symmetry, "symbol", None),
                            "poscar_path": str(poscar_path),
                        }
                    )
        pd.DataFrame(summary_rows).to_excel(self.poscar_dir / "materials_project_summary.xlsx", index=False)

    def _prepare_mlip_workspace(self) -> None:
        """Mirror downloaded POSCARs into the potential-specific workspace."""
        for source in self.poscar_dir.rglob("POSCAR"):
            relative = source.relative_to(self.poscar_dir)
            destination = self.mlip_dir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                shutil.copy2(source, destination)
        summary = self.poscar_dir / "materials_project_summary.xlsx"
        if summary.exists():
            shutil.copy2(summary, self.mlip_dir / summary.name)

    @staticmethod
    def _write_energy_per_atom(structure_dir: Path) -> None:
        energy_path = structure_dir / "energy"
        poscar_path = structure_dir / "POSCAR"
        if not energy_path.exists() or not poscar_path.exists():
            return
        energy_text = energy_path.read_text().strip()
        if not energy_text:
            return
        lines = poscar_path.read_text().splitlines()
        tokens = lines[5].split()
        counts = (
            [int(token) for token in tokens]
            if all(token.isdigit() for token in tokens)
            else [int(token) for token in lines[6].split()]
        )
        energy_per_atom = float(energy_text.split()[0]) / sum(counts)
        (structure_dir / "energy_per_atom").write_text(f"{energy_per_atom:.16f}\n")

    def calculate_energies(self) -> None:
        from materialsframework.calculators.registry import get_calculator
        from materialsframework.tools.sqs2tdb import Sqs2tdb

        cfg = self.config
        self._prepare_mlip_workspace()
        calculator = get_calculator(cfg.mlip, **cfg.mlip_kwargs)
        solver = Sqs2tdb(fmax=1e-3, verbose=True, track_trajectory=False, calculator=calculator)
        poscars = sorted(self.mlip_dir.rglob("POSCAR"))
        print(f"Found {len(poscars)} POSCARs for {cfg.mlip} relaxation")
        for poscar in poscars:
            structure_dir = poscar.parent
            if cfg.skip_existing_energies and (structure_dir / "energy").exists():
                self._write_energy_per_atom(structure_dir)
                continue
            print(f"  Relaxing: {structure_dir.relative_to(self.mlip_dir)}")
            try:
                solver._calculate(structure_dir, relax=True)
                self._write_energy_per_atom(structure_dir)
            except Exception as exc:
                print(f"  Failed {structure_dir.name}: {exc}")

    def build_database(self) -> None:
        cfg = self.config
        database = OxideDatabase(
            files_dir=cfg.files_dir,
            mlip_sources=[(cfg.mlip_label, cfg.mlip_folder)],
            fallback_refs=cfg.fallback_refs,
            poscar_folder=cfg.poscar_folder,
            fixed_elements=cfg.fixed_elements,
        )
        blade_dirs = [cfg.files_dir / folder for folder in cfg.blade_comp_folders]
        database.scan_blade_data(
            blade_comp_dir=blade_dirs[0],
            mlip_ref_label=cfg.mlip_label,
            mlip_ref_folder=self.mlip_dir,
            oxygen_element=cfg.oxygen_element,
        )
        database.run()
        database.collect_structures(
            blade_comp_dirs=blade_dirs,
            mlip_ref_folder=self.mlip_dir,
            oxygen_element=cfg.oxygen_element,
        )

    def run(self) -> None:
        stages = [
            ("01 compositions", self.config.run_compositions, self.generate_compositions),
            ("02 Materials Project POSCARs", self.config.run_download, self.download_poscars),
            ("03 MLIP energies", self.config.run_energies, self.calculate_energies),
            ("04 oxidation database", self.config.run_database, self.build_database),
        ]
        for name, enabled, function in stages:
            if not enabled:
                print(f"Skipping {name}")
                continue
            print(f"\n--- {name} ---")
            function()


if __name__ == "__main__":
    OxidationPipeline(PipelineConfig()).run()
