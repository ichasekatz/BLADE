"""OxideDatabase — merge BLADE + Materials Project data and build per-parent Excel files.

Sections executed by OxideDatabase.run():
    1  Scan single-element CONTCAR/energy files to build per-source element refs
    2  Merge BLADE + MP data → all_energy_data.xlsx
    3  Oxide matching + append oxide/MP rows → blade_oxide_data.xlsx
    4  Per-system energy/dH files → system_energy_dH/<system>.xlsx
    5  Per-parent Excel files → parent_formula_excel_files/<formula>.xlsx
    6  Pure-element energy table → pure_element_energies.xlsx
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition, Element as _PMGElement, Structure


# ==============================================================================
# Module-level helpers (no self needed)
# ==============================================================================

def normalize_formula(formula: str) -> str:
    try:
        return Composition(str(formula)).reduced_formula
    except Exception:
        return str(formula).strip()


def calculate_dH(formula: str, energy_per_atom: float, refs: dict) -> float:
    comp = Composition(formula)
    total_atoms = comp.num_atoms
    ref_total = sum(amount * refs[el] for el, amount in comp.get_el_amt_dict().items() if el in refs)
    return energy_per_atom - (ref_total / total_atoms)


def get_formula_elements(formula: str) -> set[str]:
    try:
        return set(Composition(str(formula)).get_el_amt_dict().keys())
    except Exception:
        return set()


def is_oxide_formula(formula: str) -> bool:
    elems = get_formula_elements(formula)
    return "O" in elems and len(elems) >= 2


def oxide_allowed_for_elements(formula: str, allowed_elements: set[str]) -> bool:
    elems = get_formula_elements(formula)
    if "O" not in elems:
        return False
    return (elems - {"O"}).issubset(allowed_elements)


def split_list(text) -> list[str]:
    if pd.isna(text) or str(text).strip() == "":
        return []
    return [x.strip() for x in str(text).split(",") if x.strip()]


def is_true_stable(value) -> bool:
    if pd.isna(value):
        return False
    if value is True:
        return True
    return str(value).strip().lower() in {"true", "1", "yes"}


def _is_element_col(col: str) -> bool:
    try:
        _PMGElement(col)
        return True
    except Exception:
        return False


def clean_file_name(name: str) -> str:
    name = str(name)
    if name.strip() == "" or name.lower() == "nan":
        name = "unknown_parent"
    name = re.sub(r'[<>:"/\\|?*]', "_", name)
    return name.replace(" ", "_")


def get_formula_counts(formula: str) -> dict[str, float]:
    try:
        return Composition(str(formula)).get_el_amt_dict()
    except Exception:
        return {}


def is_single_element_formula(formula: str, element: str) -> bool:
    try:
        return set(Composition(str(formula)).get_el_amt_dict().keys()) == {element}
    except Exception:
        return False


def _fmt_stable(val) -> str:
    if pd.isna(val):
        return ""
    return "TRUE" if str(val).strip().lower() in {"true", "1", "yes"} else "FALSE"


def _row_els(formula: str) -> frozenset:
    try:
        return frozenset(Composition(str(formula)).get_el_amt_dict().keys())
    except Exception:
        return frozenset()


# ==============================================================================
# OxideDatabase class
# ==============================================================================

class OxideDatabase:
    """Orchestrate all sections of the oxide database build pipeline."""

    def __init__(
        self,
        files_dir: Path,
        mlip_sources: list[tuple[str, str]],
        fallback_refs: dict[str, float],
        include_mlip: bool = True,
        include_mp: bool = True,
        poscar_folder: str = "MaterialsProject_Comps_POSCARs",
        fixed_elements: frozenset[str] = frozenset({"B"}),
    ):
        self.files_dir = files_dir
        self.mlip_sources = mlip_sources
        self.fallback_refs = fallback_refs
        self.include_mlip = include_mlip
        self.include_mp = include_mp
        self.poscar_folder = poscar_folder
        self.fixed_elements = fixed_elements

        # Derived paths
        mp_folder_base = mlip_sources[0][1] if mlip_sources else "MaterialsProject_Comps"
        self.mp_comps_dir = files_dir / mp_folder_base

        # Build {label: Path} for each MLIP source that actually exists on disk
        self.mlip_dirs: dict[str, Path] = {
            label: files_dir / folder
            for label, folder in mlip_sources
            if (files_dir / folder).exists()
        }

        # Which MLIP source to use for BLADE/fixed-element refs
        self._mlip_ref_label = "MP+GRACE"

        # Shared state populated during run()
        self.mp_summary_df: pd.DataFrame = pd.DataFrame()
        self.blade_df: pd.DataFrame = pd.DataFrame()
        self.combined_df: pd.DataFrame = pd.DataFrame()
        self.element_refs: dict[str, float] = {}
        self._all_element_refs: dict[str, dict[str, float]] = {}
        self.mlip_lookups: dict[str, dict] = {}
        self._mlip_lookups_all: dict[str, dict] = {}

        # Populated in _build_parent_excel
        self._s6_poscar_index: dict[str, Path] = {}
        self._s6_contcar_index: dict[str, Path] = {}

        print(f"MP base folder : {mp_folder_base}")
        for label, path in self.mlip_dirs.items():
            print(f"  {label:<15}: {path.name}")

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def scan_blade_data(
        self,
        blade_comp_dir: Path,
        mlip_ref_label: str = "MP+GRACE",
        mlip_ref_folder: Path | None = None,
        secondary_elements: list[str] | None = None,
        oxygen_element: str = "O",
    ) -> None:
        """Scan BLADE CONTCARs, compute dH, write blade_generated_data.xlsx.

        Call this before run() to populate blade_generated_data.xlsx from a
        BLADE Comps directory. Moved from OxideCompositions so all database
        steps are co-located with the database script.
        """
        from pymatgen.core import Composition as _Comp, Structure as _St

        non_primary = set(secondary_elements or []) | self.fixed_elements | {oxygen_element}
        element_refs = self._scan_blade_element_refs(mlip_ref_label, mlip_ref_folder)
        print(f"\nElement refs ({mlip_ref_label if element_refs else 'MP DFT fallback'}): {element_refs}")

        rows = []
        all_elements: set[str] = set()

        for contcar_path in Path(blade_comp_dir).rglob("CONTCAR"):
            sqs_folder = contcar_path.parent
            energy_path = sqs_folder / "energy"
            try:
                structure = _St.from_file(contcar_path)
            except Exception as e:
                print(f"Skipping {contcar_path}: {e}")
                continue

            comp_dict = structure.composition.get_el_amt_dict()
            all_elements.update(comp_dict.keys())
            n_atoms = len(structure)
            formula = structure.composition.reduced_formula

            _metals_in = sorted(el for el in comp_dict if el not in non_primary)
            _fixed_in  = sorted(el for el in comp_dict if el in non_primary)
            if len(_metals_in) == 1:
                parent_system = "".join(_metals_in + _fixed_in)
            else:
                _comp_dir_name = sqs_folder.parent.parent.name
                try:
                    _all_metals = sorted(_Comp(_comp_dir_name).get_el_amt_dict().keys())
                except Exception:
                    _all_metals = _metals_in
                parent_system = "".join(_all_metals + _fixed_in)

            row: dict = {
                "formula": formula, "parent_system": parent_system,
                "energy": None, "dH": None, "entry_id": sqs_folder.name,
                "file_path": str(contcar_path.resolve()),
                "is_stable": None, "source": "BLADE", "volume": structure.volume,
            }
            for el, amount in comp_dict.items():
                row[el] = amount

            if energy_path.exists():
                try:
                    total_energy = float(energy_path.read_text().strip().split()[0])
                    epa = total_energy / n_atoms
                    row["energy"] = epa
                    try:
                        from blade.oxidation.compositions import calculate_dH
                        row["dH"] = calculate_dH(formula, epa, element_refs)
                    except Exception:
                        row["dH"] = None
                except Exception:
                    pass
            rows.append(row)

        blade_df = pd.DataFrame(rows)
        element_cols = sorted(all_elements)
        for el in element_cols:
            if el not in blade_df.columns:
                blade_df[el] = 0
        blade_df[element_cols] = blade_df[element_cols].fillna(0)
        ordered_cols = (
            ["formula", "parent_system"] + element_cols
            + ["energy", "dH", "entry_id", "file_path", "is_stable", "source", "volume"]
        )
        blade_df = blade_df[ordered_cols]
        blade_df.to_excel(self.files_dir / "blade_generated_data.xlsx", index=False)
        print(blade_df)

    def _scan_blade_element_refs(
        self,
        mlip_ref_label: str,
        mlip_ref_folder: Path | None,
    ) -> dict[str, float]:
        """Scan MLIP folder for pure-element reference energies."""
        element_refs: dict[str, float] = {}
        if mlip_ref_folder and Path(mlip_ref_folder).exists():
            from pymatgen.core import Structure as _St
            for _contcar in Path(mlip_ref_folder).rglob("CONTCAR"):
                _ef = _contcar.parent / "energy"
                if not _ef.exists():
                    continue
                try:
                    _s = _St.from_file(_contcar)
                    _els = set(_s.composition.get_el_amt_dict().keys())
                    if len(_els) != 1:
                        continue
                    _el = next(iter(_els))
                    _epa = float(_ef.read_text().strip().split()[0]) / len(_s)
                    if _el not in element_refs or _epa < element_refs[_el]:
                        element_refs[_el] = _epa
                except Exception:
                    continue
            if element_refs:
                print(f"Element refs from {mlip_ref_label} ({len(element_refs)} elements):")
                for _el, _v in sorted(element_refs.items()):
                    print(f"  {_el}: {_v:.6f} eV/atom")
        return element_refs

    def collect_structures(
        self,
        blade_comp_dirs: Path | list[Path],
        mlip_ref_folder: Path | None = None,
        output_dir: Path | None = None,
        oxygen_element: str = "O",
    ) -> None:
        """Collect CONTCAR + energy files per system into output_dir.

        For each system in system_energy_dH/:
          - Copies MLIP-relaxed MP structures (lowest energy per formula)
          - Copies BLADE SQS structures matching the system
          - Copies TDB files for the system
          - Copies parent_formula Excel files

        Call after run().
        """
        import shutil
        from pymatgen.core import Composition as _Comp, Structure as _St

        system_dH_dir    = self.files_dir / "system_energy_dH"
        parent_excel_dir = self.files_dir / "parent_formula_excel_files"
        out = Path(output_dir) if output_dir else self.files_dir / "system_structures"
        if isinstance(blade_comp_dirs, (str, Path)):
            blade_comp_dirs = [Path(blade_comp_dirs)]
        blade_comp_dirs = [Path(d) for d in blade_comp_dirs]
        out.mkdir(parents=True, exist_ok=True)

        def _normalize(f: str) -> str:
            try:
                return _Comp(str(f)).reduced_formula
            except Exception:
                return str(f).strip()

        def _els(f: str) -> frozenset:
            try:
                return frozenset(_Comp(str(f)).get_el_amt_dict().keys())
            except Exception:
                return frozenset()

        def _read_energy(p: Path) -> float | None:
            try:
                return float(p.read_text().strip().split()[0])
            except Exception:
                return None

        def _copy(contcar: Path, energy: Path | None, dest: Path) -> None:
            dest.mkdir(parents=True, exist_ok=True)
            shutil.copy2(contcar, dest / "CONTCAR")
            if energy and energy.exists():
                shutil.copy2(energy, dest / "energy")

        # Build stable IDs from MP summary
        stable_ids: set[str] = set()
        grace_base = mlip_ref_folder or (self.files_dir / "MaterialsProject_Comps_GRACE")
        _summary = grace_base / "materials_project_summary.xlsx"
        if _summary.exists():
            _df = pd.read_excel(_summary)
            _df["energy_above_hull"] = pd.to_numeric(_df["energy_above_hull"], errors="coerce")
            _df["is_stable"] = _df["is_stable"].apply(lambda v: str(v).strip().lower() in {"true","1","yes"})
            stable_ids = set(_df[_df["is_stable"] | (_df["energy_above_hull"]==0)]["material_id"].astype(str).str.strip())
        print(f"Stable MP entries: {len(stable_ids)}")

        # MLIP structure index
        mlip_index: dict[str, dict] = {}
        for contcar in sorted(grace_base.rglob("CONTCAR")):
            entry_id = contcar.parent.name
            mp_id = entry_id.split("_")[0]
            if stable_ids and mp_id not in stable_ids:
                continue
            try:
                s = _St.from_file(str(contcar))
                f = s.composition.reduced_formula
                epa = (_read_energy(contcar.parent/"energy") or 0) / len(s)
            except Exception:
                continue
            if f not in mlip_index or epa < mlip_index[f]["epa"]:
                mlip_index[f] = {"entry_id": entry_id, "contcar": contcar,
                                  "energy": contcar.parent/"energy" if (contcar.parent/"energy").exists() else None,
                                  "epa": epa}

        # POSCAR fallback
        poscar_base = self.files_dir / "MaterialsProject_Comps"
        poscar_index: dict[str, dict] = {}
        if poscar_base.exists():
            for poscar in sorted(poscar_base.rglob("POSCAR")):
                entry_id = poscar.parent.name
                mp_id = entry_id.split("_")[0]
                if stable_ids and mp_id not in stable_ids:
                    continue
                try:
                    f = _St.from_file(str(poscar)).composition.reduced_formula
                except Exception:
                    continue
                if f not in poscar_index:
                    poscar_index[f] = {"entry_id": entry_id, "poscar": poscar}

        # BLADE structure index — scan all provided comp dirs
        blade_index: dict[str, list[dict]] = {}
        for blade_comp_dir in blade_comp_dirs:
            for contcar in sorted(blade_comp_dir.rglob("CONTCAR")):
                try:
                    f = _St.from_file(str(contcar)).composition.reduced_formula
                except Exception:
                    continue
                blade_index.setdefault(f, []).append({
                    "entry_id": contcar.parent.name,
                    "contcar": contcar,
                    "energy": contcar.parent/"energy" if (contcar.parent/"energy").exists() else None,
                    "comp_dir": contcar.parent.parent.parent.name,
                })

        print(f"MLIP index: {len(mlip_index)} formulas | BLADE index: {len(blade_index)} formulas")

        for xlsx in sorted(p for p in system_dH_dir.glob("*.xlsx") if not p.name.startswith("~$")):
            system = xlsx.stem
            print(f"\nSystem: {system}")
            sys_els = _els(system)
            sys_out = out / system

            df = pd.read_excel(xlsx)
            if "formula" not in df.columns:
                continue

            # MLIP-relaxed structures → grace/ folder; MP POSCAR fallback → mp/ folder
            grace_label = Path(grace_base).name.replace("MaterialsProject_Comps_", "") or "grace"
            n_mlip = n_mp = n_miss = 0
            for raw in df["formula"].dropna().unique():
                f = _normalize(str(raw))
                hit = mlip_index.get(f)
                if hit:
                    _copy(hit["contcar"], hit["energy"], sys_out/grace_label/f"{f}_{hit['entry_id']}")
                    n_mlip += 1
                else:
                    fb = poscar_index.get(f)
                    if fb:
                        _copy(fb["poscar"], None, sys_out/"mp"/f"{f}_{fb['entry_id']}")
                        n_mp += 1
                    else:
                        n_miss += 1

            # BLADE structures
            n_blade = 0
            seen: set = set()
            for f, entries in blade_index.items():
                f_els = _els(f)
                if not f_els or not (f_els - {oxygen_element}).issubset(sys_els):
                    continue
                for hit in entries:
                    if hit["contcar"] in seen:
                        continue
                    comp_els = _els(hit["comp_dir"])
                    if comp_els and not comp_els.issubset(sys_els):
                        continue
                    seen.add(hit["contcar"])
                    _copy(hit["contcar"], hit["energy"], sys_out/"blade"/f"{f}_{hit['entry_id']}")
                    n_blade += 1

            # TDB files — search all comp dirs
            n_tdb = 0
            for blade_comp_dir in blade_comp_dirs:
                for comp_dir in blade_comp_dir.iterdir():
                    if not comp_dir.is_dir() or not _els(comp_dir.name).issubset(sys_els):
                        continue
                    for tdb in comp_dir.glob("*.tdb"):
                        (sys_out/"tdb").mkdir(parents=True, exist_ok=True)
                        shutil.copy2(tdb, sys_out/"tdb"/tdb.name)
                        n_tdb += 1

            # System energy Excel
            if xlsx.exists():
                shutil.copy2(xlsx, sys_out/f"{system}_system_energy.xlsx")

            # Parent formula Excels
            n_parent = 0
            if parent_excel_dir.exists():
                for src in parent_excel_dir.iterdir():
                    if not src.is_dir():
                        continue
                    for pxlsx in src.glob("*.xlsx"):
                        if pxlsx.name.startswith("~$"):
                            continue
                        if not (_els(pxlsx.stem) - {oxygen_element}).issubset(sys_els):
                            continue
                        dest = sys_out/"parent_excel"/src.name
                        dest.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(pxlsx, dest/pxlsx.name)
                        n_parent += 1

            print(f"  {grace_label}: {n_mlip} | mp: {n_mp} | blade: {n_blade} | tdb: {n_tdb} | parent_excel: {n_parent} | missing: {n_miss}")

        print(f"\nDone. Output: {out.resolve()}")

    def run(self) -> None:
        self.files_dir.mkdir(parents=True, exist_ok=True)
        self._scan_refs()
        self._build_all_energy_data()
        self._build_oxide_data()
        self._build_system_energy_dh()
        self._build_parent_excel()
        self._build_pure_element_energies()

    # ------------------------------------------------------------------
    # Section 1: Scan element refs from MLIP CONTCAR/energy files
    # ------------------------------------------------------------------

    def _scan_refs(self) -> None:
        # Build stable mp-id set from summary so element refs use only stable structures
        stable_ids: set[str] = set()
        mp_summary_path = self.mp_comps_dir / "materials_project_summary.xlsx"
        if mp_summary_path.exists():
            try:
                _smdf = pd.read_excel(mp_summary_path)
                stable_ids = set(
                    _smdf[
                        _smdf["is_stable"].apply(lambda v: str(v).strip().lower() in {"true", "1", "yes"})
                        | (pd.to_numeric(_smdf.get("energy_above_hull", pd.Series(dtype=float)), errors="coerce") == 0)
                    ]["material_id"].astype(str).str.strip()
                )
            except Exception:
                pass

        self._all_element_refs = {}
        for lbl, mdir in self.mlip_dirs.items():
            scanned = self._scan_element_refs_folder(mdir, stable_ids=stable_ids or None)
            if scanned:
                # Sanity-check: replace values >5x outside the fallback range with fallback
                for el, val in list(scanned.items()):
                    fb = self.fallback_refs.get(el)
                    if fb is not None and abs(val - fb) > abs(fb) * 5:
                        print(
                            f"  WARNING: {el} scanned ref {val:.3f} looks bad "
                            f"(fallback={fb:.3f}), using fallback"
                        )
                        scanned[el] = fb
                    elif fb is None and (abs(val) > 1000 or val > 0):
                        print(f"  WARNING: {el} scanned ref {val:.3f} removed (no fallback)")
                        del scanned[el]
                self._all_element_refs[lbl] = scanned
                print(f"Element refs scanned from {lbl}: {len(scanned)} elements")

        # Primary refs for BLADE dH (falls back to hardcoded if not found)
        self.element_refs = self._all_element_refs.get(self._mlip_ref_label, self.fallback_refs)

    @staticmethod
    def _scan_element_refs_folder(folder: Path,
                                   stable_ids: set[str] | None = None) -> dict[str, float]:
        best: dict[str, float] = {}
        if not folder.exists():
            return best
        for contcar in folder.rglob("CONTCAR"):
            if stable_ids is not None:
                mp_id = contcar.parent.name.split("_")[0]
                if mp_id not in stable_ids:
                    continue
            ef = contcar.parent / "energy"
            if not ef.exists():
                continue
            try:
                s = Structure.from_file(contcar)
                els = set(s.composition.get_el_amt_dict().keys())
                if len(els) != 1:
                    continue
                el = next(iter(els))
                epa = float(ef.read_text().strip().split()[0]) / len(s)
                if el not in best or epa < best[el]:
                    best[el] = epa
            except Exception:
                continue
        return best

    # ------------------------------------------------------------------
    # Section 2: Merge BLADE + MP data → all_energy_data.xlsx
    # ------------------------------------------------------------------

    def _build_all_energy_data(self) -> None:
        blade_xlsx = self.files_dir / "blade_generated_data.xlsx"
        mp_summary_xlsx = self.mp_comps_dir / "materials_project_summary.xlsx"
        all_energy_data_xlsx = self.files_dir / "all_energy_data.xlsx"

        blade_df = pd.read_excel(blade_xlsx)

        # Load mp_summary_df whenever MP DFT rows or MLIP rows are needed
        need_mp_summary = self.include_mp or self.include_mlip
        if need_mp_summary and mp_summary_xlsx.exists():
            self.mp_summary_df = pd.read_excel(mp_summary_xlsx)

        # Build CONTCAR index for MLIP sources (needed even when include_mp=False)
        mlip_contcar_idx: dict[str, dict[str, Path]] = {
            label: self._build_contcar_index(path)
            for label, path in self.mlip_dirs.items()
        } if self.include_mlip else {}

        # Pre-build POSCAR index for MP DFT file paths
        _mp_poscar_dir = self.files_dir / self.poscar_folder
        _idx_base_poscar: dict[str, Path] = {}
        if self.include_mp:
            _poscar_source = str(_mp_poscar_dir) if _mp_poscar_dir.exists() else "MLIP folders"
            _poscar_scan_dirs = [_mp_poscar_dir] if _mp_poscar_dir.exists() else list(self.mlip_dirs.values())
            for scan_dir in _poscar_scan_dirs:
                for p in scan_dir.rglob("POSCAR"):
                    key = p.parent.name.split("_")[0]
                    if key not in _idx_base_poscar:
                        _idx_base_poscar[key] = p
            print(f"  Base POSCAR index: {len(_idx_base_poscar)} entries (from {_poscar_source})")

        rows = []
        all_elements: set[str] = set()

        if need_mp_summary and not self.mp_summary_df.empty:
            for _, mp_row in self.mp_summary_df.iterrows():
                material_id = str(mp_row.get("material_id", "")).strip()
                formula_dft = str(mp_row.get("formula", ""))

                # MP DFT row — only when include_mp=True
                if self.include_mp:
                    poscar_path = str(_idx_base_poscar[material_id].resolve()) if material_id in _idx_base_poscar else ""
                    vol = ""
                    if poscar_path:
                        try:
                            s = Structure.from_file(poscar_path)
                            vol = s.volume
                            all_elements.update(s.composition.get_el_amt_dict().keys())
                        except Exception:
                            pass
                    dft_row: dict = {
                        "formula": formula_dft,
                        "energy": mp_row.get("energy_per_atom"),
                        "dH": mp_row.get("formation_energy_per_atom"),
                        "entry_id": material_id,
                        "file_path": poscar_path,
                        "is_stable": _fmt_stable(mp_row.get("is_stable")),
                        "source": "Materials Project",
                        "volume": vol,
                    }
                    try:
                        dft_row.update(Composition(formula_dft).get_el_amt_dict())
                    except Exception:
                        pass
                    rows.append(dft_row)

                # MLIP rows — only when include_mlip=True
                if self.include_mlip:
                    for label, idx in mlip_contcar_idx.items():
                        if material_id in idx:
                            r = self._mlip_row_from_contcar(
                                idx[material_id], material_id, mp_row, label, all_elements
                            )
                            if r:
                                rows.append(r)

        mp_df = pd.DataFrame(rows)

        fixed_cols = ["formula", "energy", "dH", "entry_id", "file_path", "is_stable", "source", "volume"]
        all_element_cols = sorted(
            set(c for c in blade_df.columns if c not in fixed_cols + ["phase_label", "experiment_dH", "parent_system"])
            | all_elements
        )

        # Rename phase_label → formula in blade_df
        if "phase_label" in blade_df.columns:
            blade_df = blade_df.rename(columns={"phase_label": "formula"})
        if "formula" in blade_df.columns and blade_df.columns.tolist().count("formula") > 1:
            blade_df = blade_df.loc[:, ~blade_df.columns.duplicated(keep="first")]

        blade_df["is_stable"] = blade_df["is_stable"].apply(_fmt_stable) if "is_stable" in blade_df.columns else ""

        for el in all_element_cols:
            if el not in blade_df.columns:
                blade_df[el] = 0
            if el not in mp_df.columns:
                mp_df[el] = 0

        blade_df[all_element_cols] = blade_df[all_element_cols].fillna(0)
        mp_df[all_element_cols] = mp_df[all_element_cols].fillna(0)

        # Deduplicate BLADE rows by formula — keep lowest-energy SQS per formula
        blade_df["energy"] = pd.to_numeric(blade_df["energy"], errors="coerce")
        blade_df = (
            blade_df
            .sort_values("energy", na_position="last")
            .drop_duplicates(subset=["formula"], keep="first")
            .reset_index(drop=True)
        )

        ordered_cols = (
            ["formula"]
            + all_element_cols
            + ["energy", "dH", "entry_id", "file_path", "is_stable", "source", "volume"]
        )
        ordered_cols = [c for c in ordered_cols if c in blade_df.columns or c in mp_df.columns]

        combined_df = pd.concat(
            [blade_df[[c for c in ordered_cols if c in blade_df.columns]],
             mp_df[[c for c in ordered_cols if c in mp_df.columns]]],
            ignore_index=True,
        )
        combined_df = (
            combined_df
            .sort_values(
                ["source", "formula"],
                key=lambda s: s.str.lower() if s.dtype == object else s,
            )
            .reset_index(drop=True)
        )
        combined_df.to_excel(all_energy_data_xlsx, index=False)
        print(f"Saved combined file: {all_energy_data_xlsx.resolve()}")
        print(f"BLADE rows: {len(blade_df)}  |  MP rows: {len(mp_df)}  |  Total: {len(combined_df)}")

        # Store for later sections
        self.blade_df = blade_df
        self.combined_df = combined_df

    def _build_contcar_index(self, folder: Path | None) -> dict[str, Path]:
        if folder is None or not folder.exists():
            return {}
        print(f"Indexing CONTCARs in {folder} ...")
        idx: dict[str, Path] = {}
        for c in folder.rglob("CONTCAR"):
            key = c.parent.name.split("_")[0]
            if key not in idx:
                idx[key] = c
        print(f"  Found {len(idx)} CONTCARs")
        return idx

    def _mlip_row_from_contcar(
        self,
        contcar: Path,
        material_id: str,
        mp_row: pd.Series,
        source_label: str,
        comp_elements: set,
    ) -> dict | None:
        energy_file = contcar.parent / "energy"
        src_refs = self._all_element_refs.get(source_label, self.element_refs)
        try:
            s = Structure.from_file(contcar)
            comp_dict = s.composition.get_el_amt_dict()
            comp_elements.update(comp_dict.keys())
            n = len(s)
            formula = s.composition.reduced_formula
            if energy_file.exists():
                epa = float(energy_file.read_text().strip().split()[0]) / n
                dh = calculate_dH(formula, epa, src_refs)
            else:
                epa = mp_row.get("energy_per_atom")
                dh = mp_row.get("formation_energy_per_atom")
            row = {
                "formula": formula,
                "energy": epa,
                "dH": dh,
                "entry_id": material_id,
                "file_path": str(contcar.resolve()),
                "is_stable": _fmt_stable(mp_row.get("is_stable")),
                "source": source_label,
                "volume": s.volume,
            }
            row.update(Composition(formula).get_el_amt_dict())
            return row
        except Exception as e:
            print(f"  Error reading {contcar}: {e}")
            return None

    # ------------------------------------------------------------------
    # Section 3: Oxide matching + append oxide/MP rows → blade_oxide_data.xlsx
    # ------------------------------------------------------------------

    def _build_oxide_data(self) -> None:
        # Controls which rows get appended under each BLADE parent
        include_mp_dft = True
        oxide_column = "stable_mp_oxides"  # or "possible_mp_oxides"

        blade_oxide_data_xlsx = self.files_dir / "blade_oxide_data.xlsx"

        blade_df = pd.read_excel(self.files_dir / "blade_generated_data.xlsx")
        blade_df = blade_df.drop_duplicates(subset=["formula"], keep="first").reset_index(drop=True)
        mp_summary = pd.read_excel(self.mp_comps_dir / "materials_project_summary.xlsx")

        _meta_cols = {
            "formula", "phase_label", "energy", "dH", "entry_id", "file_path",
            "is_stable", "source", "volume", "parent_formula", "parent_entry_id",
        }
        blade_element_cols = sorted(
            [c for c in blade_df.columns if c not in _meta_cols and _is_element_col(c)],
            key=lambda e: _PMGElement(e).Z,
        )
        metal_elements = [el for el in blade_element_cols if el not in self.fixed_elements | {"O"}]
        drop_cols = {"phase_label", "formula"}

        def get_blade_allowed_elements(row: pd.Series) -> set[str]:
            allowed: set[str] = set()
            for el in metal_elements:
                if el in row.index and pd.notna(row[el]) and float(row[el]) > 0:
                    allowed.add(el)
            for _fe in self.fixed_elements:
                if _fe in row.index and pd.notna(row[_fe]) and float(row[_fe]) > 0:
                    allowed.add(_fe)
            return allowed

        def make_row(d: dict, role: str, parent_formula: str, parent_entry_id, source_label: str) -> dict:
            out = {k: v for k, v in d.items() if k not in drop_cols}
            out["row_role"] = role
            out["parent_formula"] = parent_formula
            out["parent_entry_id"] = parent_entry_id
            out["source"] = source_label
            return out

        # Index MP summary by normalized formula
        mp_summary["formula_norm"] = mp_summary["formula"].astype(str).apply(normalize_formula)
        mp_oxides = mp_summary[mp_summary["formula"].apply(is_oxide_formula)].copy()
        mp_oxides["energy_above_hull"] = pd.to_numeric(
            mp_oxides.get("energy_above_hull", None), errors="coerce"
        )
        mp_oxides["formation_energy_per_atom"] = pd.to_numeric(
            mp_oxides.get("formation_energy_per_atom", None), errors="coerce"
        )

        # Best-per-formula for MP DFT
        sort_col = "energy_above_hull" if "energy_above_hull" in mp_oxides.columns else "formation_energy_per_atom"
        mp_best = (
            mp_oxides.sort_values(sort_col)
            .drop_duplicates("formula_norm", keep="first")
            .set_index("formula_norm")
        )

        # Stable material_id set for MLIP lookup filtering
        stable_mp_ids: set[str] = set(
            self.mp_summary_df[
                self.mp_summary_df["is_stable"].apply(
                    lambda v: str(v).strip().lower() in {"true", "1", "yes"}
                )
                | (
                    pd.to_numeric(
                        self.mp_summary_df.get("energy_above_hull", pd.Series(dtype=float)),
                        errors="coerce",
                    ) == 0
                )
            ]["material_id"].astype(str).str.strip()
        )

        # Stable-only MLIP lookups — used for dH refs and system_energy_dH
        self.mlip_lookups = {
            label: self._build_mlip_lookup(path, self._all_element_refs.get(label, self.element_refs), stable_mp_ids)
            for label, path in self.mlip_dirs.items()
            if self.include_mlip
        }
        # Unfiltered MLIP lookups — used by _apply_mlip_values so all structures get correct source/energy
        self._mlip_lookups_all = {
            label: self._build_mlip_lookup(path, self._all_element_refs.get(label, self.element_refs), stable_ids=None)
            for label, path in self.mlip_dirs.items()
            if self.include_mlip
        }

        # All element cols for output ordering (re-derive from combined)
        all_element_cols = sorted(
            el for el in self.combined_df.columns
            if _is_element_col(el)
        )

        expanded_rows: list[dict] = []
        missing_rows: list[dict] = []

        for _, blade_row in blade_df.iterrows():
            if str(blade_row.get("source", "")).upper() != "BLADE":
                continue

            parent_formula = str(blade_row.get("formula", blade_row.get("phase_label", "")))
            parent_entry_id = blade_row.get("entry_id", "")
            allowed = get_blade_allowed_elements(blade_row)

            # Compute oxide summary metadata for this parent
            possible = mp_oxides[mp_oxides["formula"].apply(lambda f: oxide_allowed_for_elements(f, allowed))]
            stable = possible[possible["is_stable"] == True]
            possible_formulas = sorted(possible["formula"].dropna().unique())
            stable_formulas = sorted(stable["formula"].dropna().unique())

            # Lowest eHull/dH/energy per source
            ehull_cols: dict[str, tuple] = {}
            dH_cols: dict[str, tuple] = {}
            energy_cols: dict[str, tuple] = {}
            if not possible.empty:
                if possible["energy_above_hull"].notna().any():
                    b = possible.sort_values("energy_above_hull").iloc[0]
                    ehull_cols["mp_dft"] = (b["formula"], b["energy_above_hull"])
                if possible["formation_energy_per_atom"].notna().any():
                    b = possible.sort_values("formation_energy_per_atom").iloc[0]
                    dH_cols["mp_dft"] = (b["formula"], b["formation_energy_per_atom"])
                if possible["energy_per_atom"].notna().any():
                    b = possible.sort_values("energy_per_atom").iloc[0]
                    energy_cols["mp_dft"] = (b["formula"], b["energy_per_atom"])
            for lbl, lkp in self.mlip_lookups.items():
                lkp_key = lbl.replace("+", "_").replace(" ", "_").lower()
                mlip_oxides_dH = [
                    (fml, lkp[normalize_formula(fml)]["dH"])
                    for fml in possible_formulas
                    if normalize_formula(fml) in lkp and lkp[normalize_formula(fml)].get("dH") is not None
                ]
                mlip_oxides_E = [
                    (fml, lkp[normalize_formula(fml)]["energy"])
                    for fml in possible_formulas
                    if normalize_formula(fml) in lkp and lkp[normalize_formula(fml)].get("energy") is not None
                ]
                if mlip_oxides_dH:
                    dH_cols[lkp_key] = min(mlip_oxides_dH, key=lambda x: x[1])
                if mlip_oxides_E:
                    energy_cols[lkp_key] = min(mlip_oxides_E, key=lambda x: x[1])

            # Non-oxide sub-compositions: BLADE rows + pure element references
            non_oxide_seen: set[str] = set()
            non_oxide_rows: list[dict] = []

            for el in sorted(allowed):
                if el in self.element_refs and el not in non_oxide_seen and el != parent_formula:
                    non_oxide_rows.append({"formula": el, "dH": None, "is_stable": True})
                    non_oxide_seen.add(el)

            for _, br in blade_df.iterrows():
                bf = str(br.get("formula", br.get("phase_label", "")))
                if bf == parent_formula or bf in non_oxide_seen:
                    continue
                ba = get_blade_allowed_elements(br)
                if ba and ba.issubset(allowed):
                    dh_val = pd.to_numeric(br.get("dH"), errors="coerce")
                    stable_val = str(br.get("is_stable", "")).strip().upper() == "TRUE"
                    non_oxide_rows.append({
                        "formula": bf,
                        "dH": dh_val if not pd.isna(dh_val) else None,
                        "is_stable": stable_val,
                    })
                    non_oxide_seen.add(bf)

            # Add stable MP non-oxide compounds whose elements ⊆ allowed
            for _, _mp_r in mp_summary.iterrows():
                _mf = str(_mp_r.get("formula", "")).strip()
                if not _mf or _mf in non_oxide_seen:
                    continue
                if is_oxide_formula(_mf):
                    continue
                try:
                    _mels = set(Composition(_mf).get_el_amt_dict().keys())
                except Exception:
                    continue
                if not _mels.issubset(allowed):
                    continue
                if not is_true_stable(_mp_r.get("is_stable")):
                    continue
                non_oxide_rows.append({
                    "formula": _mf,
                    "dH": pd.to_numeric(_mp_r.get("formation_energy_per_atom"), errors="coerce"),
                    "is_stable": True,
                })
                non_oxide_seen.add(_mf)

            non_oxide_formulas = [r["formula"] for r in non_oxide_rows]
            dh_ranked = [r for r in non_oxide_rows if r["dH"] is not None]
            lowest_dH_non_oxide, lowest_dH_non_oxide_val = None, None
            if dh_ranked:
                best_nr = min(dh_ranked, key=lambda r: r["dH"])
                lowest_dH_non_oxide, lowest_dH_non_oxide_val = best_nr["formula"], best_nr["dH"]
            # Use dH as proxy for eHull (BLADE does not have hull distance)
            lowest_ehull_non_oxide, lowest_ehull_non_oxide_val = lowest_dH_non_oxide, lowest_dH_non_oxide_val

            # Non-oxide + non-boride: pure metals / metal alloys only (no B, no O)
            metal_only = {el for el in allowed if el in metal_elements}
            nnb_rows: list[dict] = []
            nnb_seen: set[str] = set()
            for el in sorted(metal_only):
                if el in self.element_refs and el != parent_formula:
                    nnb_rows.append({"formula": el, "dH": None, "is_stable": True})
                    nnb_seen.add(el)
            for _, br in blade_df.iterrows():
                bf = str(br.get("formula", br.get("phase_label", "")))
                if bf == parent_formula or bf in nnb_seen:
                    continue
                ba = get_blade_allowed_elements(br)
                ba_metals = {e for e in ba if e in metal_elements}
                if ba_metals and ba_metals.issubset(metal_only) and not (ba & self.fixed_elements) and "O" not in ba:
                    dh_v = pd.to_numeric(br.get("dH"), errors="coerce")
                    st_v = str(br.get("is_stable", "")).strip().upper() == "TRUE"
                    nnb_rows.append({
                        "formula": bf,
                        "dH": dh_v if not pd.isna(dh_v) else None,
                        "is_stable": st_v,
                    })
                    nnb_seen.add(bf)

            # Also add stable MP metal-only compounds (no B, no O) from mp_summary
            for _, _mp_r in mp_summary.iterrows():
                _mf = str(_mp_r.get("formula", "")).strip()
                if not _mf or _mf in nnb_seen:
                    continue
                try:
                    _mels = set(Composition(_mf).get_el_amt_dict().keys())
                except Exception:
                    continue
                if (_mels & self.fixed_elements) or "O" in _mels:
                    continue
                if not _mels.issubset(metal_only) or len(_mels) < 2:
                    continue
                if not is_true_stable(_mp_r.get("is_stable")):
                    continue
                nnb_rows.append({"formula": _mf, "dH": pd.to_numeric(_mp_r.get("formation_energy_per_atom"), errors="coerce"), "is_stable": True})
                nnb_seen.add(_mf)

            nnb_formulas = [r["formula"] for r in nnb_rows]
            nnb_stable_formulas = [r["formula"] for r in nnb_rows if r["is_stable"]]
            nnb_dh_ranked = [r for r in nnb_rows if r["dH"] is not None]
            lowest_dH_nnb, lowest_dH_nnb_val = None, None
            lowest_ehull_nnb, lowest_ehull_nnb_val = None, None
            if nnb_dh_ranked:
                nnb_best = min(nnb_dh_ranked, key=lambda r: r["dH"])
                lowest_dH_nnb, lowest_dH_nnb_val = nnb_best["formula"], nnb_best["dH"]
                lowest_ehull_nnb, lowest_ehull_nnb_val = lowest_dH_nnb, lowest_dH_nnb_val

            # BLADE parent row
            br_out = {k: v for k, v in blade_row.to_dict().items() if k not in drop_cols}
            br_out["parent_formula"] = parent_formula
            br_out["possible_non_oxides"] = ", ".join(sorted(non_oxide_formulas))
            br_out["n_non_oxides"] = len(non_oxide_formulas)
            br_out["possible_non_oxide_non_fixed"] = ", ".join(sorted(nnb_formulas))
            br_out["n_non_oxide_non_fixed"] = len(nnb_formulas)
            br_out["possible_mp_oxides"] = ", ".join(possible_formulas)
            br_out["stable_mp_oxides"] = ", ".join(stable_formulas)
            br_out["n_possible_oxides"] = len(possible_formulas)
            br_out["n_stable_oxides"] = len(stable_formulas)
            expanded_rows.append(br_out)

            for oxide_formula in split_list(blade_row.get(oxide_column, "")):
                if not oxide_allowed_for_elements(oxide_formula, allowed):
                    continue
                fkey = normalize_formula(oxide_formula)

                if include_mp_dft and fkey in mp_best.index:
                    expanded_rows.append(make_row(
                        mp_best.loc[fkey].to_dict(), "mp_dft_oxide",
                        parent_formula, parent_entry_id, "Materials Project",
                    ))

                for lbl, lkp in self.mlip_lookups.items():
                    if fkey in lkp:
                        role = f"mp_{lbl.replace('+', '').replace(' ', '_').lower()}_oxide"
                        expanded_rows.append(make_row(
                            lkp[fkey], role,
                            parent_formula, parent_entry_id, lbl,
                        ))

                found_any = (fkey in mp_best.index) or any(fkey in lkp for lkp in self.mlip_lookups.values())
                if not found_any:
                    missing_rows.append({"parent_formula": parent_formula, "oxide_formula": oxide_formula})

        expanded_df = pd.DataFrame(expanded_rows)
        missing_df = pd.DataFrame(missing_rows)

        blade_oxide_drop = drop_cols | {"is_stable", "parent_entry_id", "row_role"}
        for col in blade_oxide_drop:
            if col in expanded_df.columns:
                expanded_df = expanded_df.drop(columns=[col])
        for col in ("formula_norm",):
            if col in expanded_df.columns:
                expanded_df = expanded_df.drop(columns=[col])

        # Rename parent_formula → formula
        if "parent_formula" in expanded_df.columns:
            if "formula" in expanded_df.columns:
                expanded_df = expanded_df.drop(columns=["formula"])
            expanded_df = expanded_df.rename(columns={"parent_formula": "formula"})

        # Column order: formula → elements → energy → dH → entry_id → file_path → is_stable → source → volume → rest
        el_cols_out = [c for c in all_element_cols if c in expanded_df.columns]
        front = (
            ["formula"]
            + el_cols_out
            + ["energy", "dH", "entry_id", "file_path", "is_stable", "source", "volume"]
        )
        front = [c for c in front if c in expanded_df.columns]
        rest = [c for c in expanded_df.columns if c not in front]
        expanded_df = expanded_df[front + rest]

        with pd.ExcelWriter(blade_oxide_data_xlsx, engine="openpyxl") as writer:
            expanded_df.to_excel(writer, sheet_name="expanded_data", index=False)
            if not missing_df.empty:
                missing_df.to_excel(writer, sheet_name="missing_oxides", index=False)

        print(f"Saved: {blade_oxide_data_xlsx.resolve()}")
        print(f"Expanded rows: {len(expanded_df)}  |  Missing: {len(missing_df)}")

    @staticmethod
    def _build_mlip_lookup(
        folder: Path | None,
        refs: dict[str, float],
        stable_ids: set[str] | None = None,
    ) -> dict[str, dict]:
        if folder is None or not folder.exists():
            return {}
        lookup: dict[str, dict] = {}
        for contcar in folder.rglob("CONTCAR"):
            mp_id = contcar.parent.name.split("_")[0]
            if stable_ids is not None and mp_id not in stable_ids:
                continue
            ef = contcar.parent / "energy"
            if not ef.exists():
                continue
            try:
                s = Structure.from_file(contcar)
                n = len(s)
                fn = s.composition.reduced_formula
                epa = float(ef.read_text().strip().split()[0]) / n
                dh = calculate_dH(fn, epa, refs)
                fkey = normalize_formula(fn)
                if fkey not in lookup or epa < lookup[fkey]["energy"]:
                    lookup[fkey] = {
                        "energy": epa,
                        "dH": dh,
                        "file_path": str(contcar.resolve()),
                        "formula_norm": fkey,
                    }
            except Exception:
                continue
        return lookup

    # ------------------------------------------------------------------
    # Section 4: Per-system energy/dH files
    # ------------------------------------------------------------------

    def _build_system_energy_dh(self) -> None:
        blade_raw = pd.read_excel(self.files_dir / "blade_generated_data.xlsx")
        blade_raw["formula"] = blade_raw["formula"].astype(str).str.strip()
        blade_raw["parent_system"] = (
            blade_raw["parent_system"].astype(str).str.strip()
            if "parent_system" in blade_raw.columns
            else blade_raw["formula"]
        )

        # Build system → element set mapping
        blade_systems: dict[str, frozenset] = {}
        for _, r in blade_raw.iterrows():
            sys = str(r["parent_system"]).strip()
            if sys and sys.lower() != "nan":
                try:
                    els = frozenset(Composition(sys).get_el_amt_dict().keys())
                except Exception:
                    try:
                        els = frozenset(Composition(str(r["formula"])).get_el_amt_dict().keys())
                    except Exception:
                        continue
                if sys not in blade_systems:
                    blade_systems[sys] = els

        # MP base: stable structures + all pure elements
        mp_base = self.mp_summary_df.copy()
        mp_base["energy_above_hull"] = pd.to_numeric(
            mp_base.get("energy_above_hull", pd.Series(dtype=float)), errors="coerce"
        )
        mp_base["energy_per_atom"] = pd.to_numeric(
            mp_base.get("energy_per_atom", pd.Series(dtype=float)), errors="coerce"
        )
        mp_base["formation_energy_per_atom"] = pd.to_numeric(
            mp_base.get("formation_energy_per_atom", pd.Series(dtype=float)), errors="coerce"
        )
        mp_base["_is_stable"] = (
            mp_base["is_stable"].apply(lambda v: str(v).strip().lower() in {"true", "1", "yes"})
            if "is_stable" in mp_base.columns
            else False
        )
        mp_base["_is_pure_el"] = mp_base["formula"].apply(
            lambda f: len(get_formula_elements(str(f))) == 1
        )
        mp_base = (
            mp_base[mp_base["_is_stable"] | (mp_base["energy_above_hull"] == 0) | mp_base["_is_pure_el"]]
            .sort_values("energy_above_hull", na_position="last")
            .drop_duplicates(subset=["formula"], keep="first")
            .reset_index(drop=True)
        )
        cs_mp = mp_base[["formula"]].copy()
        cs_mp["energy_mp"] = mp_base["energy_per_atom"].values
        cs_mp["dH_mp"] = mp_base["formation_energy_per_atom"].values
        # Pure elements: normalize formula to element symbol (e.g. O2 → O), set fallback energy, dH=0
        _pure_mask = cs_mp["formula"].apply(lambda f: len(get_formula_elements(str(f))) == 1)
        cs_mp.loc[_pure_mask, "formula"] = cs_mp.loc[_pure_mask, "formula"].apply(
            lambda f: next(iter(get_formula_elements(str(f))), f)
        )
        cs_mp.loc[_pure_mask, "energy_mp"] = cs_mp.loc[_pure_mask, "formula"].apply(
            lambda f: self.fallback_refs.get(f)
        )
        cs_mp.loc[_pure_mask, "dH_mp"] = 0.0

        # BLADE lookup: formula → lowest-energy SQS entry
        blade_lkp: dict[str, dict] = {}
        for _, r in blade_raw.iterrows():
            bf = normalize_formula(str(r.get("formula", "")))
            be = pd.to_numeric(r.get("energy"), errors="coerce")
            bd = pd.to_numeric(r.get("dH"), errors="coerce")
            if not bf or bf == "nan":
                continue
            if bf not in blade_lkp or (
                pd.notna(be) and (blade_lkp[bf]["energy"] is None or be < blade_lkp[bf]["energy"])
            ):
                blade_lkp[bf] = {"energy": be, "dH": bd}

        # Append BLADE-only formulas not in MP rows
        cs_mp_norms = set(cs_mp["formula"].apply(normalize_formula))
        blade_only = [
            {"formula": bf, "energy_mp": None, "dH_mp": None}
            for bf in blade_lkp
            if bf not in cs_mp_norms
        ]
        if blade_only:
            cs_mp = pd.concat([cs_mp, pd.DataFrame(blade_only)], ignore_index=True)

        cs_mp["energy_blade"] = cs_mp["formula"].apply(
            lambda f: blade_lkp.get(normalize_formula(f), {}).get("energy")
        )
        cs_mp["dH_blade"] = cs_mp["formula"].apply(
            lambda f: 0.0 if len(_row_els(f)) == 1
            else blade_lkp.get(normalize_formula(f), {}).get("dH")
        )

        # Add MLIP columns
        for lbl, _ in self.mlip_sources:
            lkp = self.mlip_lookups.get(lbl, {})
            col_key = lbl.replace("+", "").replace(" ", "_").lower()
            cs_mp[f"energy_{col_key}"] = cs_mp["formula"].apply(
                lambda f: lkp.get(normalize_formula(f), {}).get("energy")
            )
            cs_mp[f"dH_{col_key}"] = cs_mp["formula"].apply(
                lambda f: 0.0 if len(_row_els(f)) == 1
                else lkp.get(normalize_formula(f), {}).get("dH")
            )

        # Keep only stable MP structures (BLADE rows always kept)
        stable_mp_norms = set(
            self.mp_summary_df[
                self.mp_summary_df["is_stable"].apply(
                    lambda v: str(v).strip().lower() in {"true", "1", "yes"}
                )
                | (
                    pd.to_numeric(
                        self.mp_summary_df.get("energy_above_hull", pd.Series(dtype=float)),
                        errors="coerce",
                    ) == 0
                )
            ]["formula"].apply(normalize_formula)
        )
        blade_norms = set(blade_lkp.keys())
        cs_mp = cs_mp[
            cs_mp["formula"].apply(normalize_formula).isin(stable_mp_norms | blade_norms)
        ].copy().reset_index(drop=True)

        cs_mp["_els"] = cs_mp["formula"].apply(_row_els)

        sys_dir = self.files_dir / "system_energy_dH"
        sys_dir.mkdir(parents=True, exist_ok=True)

        def _sys_row_order(f: str) -> int:
            els = get_formula_elements(str(f))
            if len(els) == 1:
                return 0
            if (els & self.fixed_elements) and "O" not in els:
                return 1
            if "O" in els:
                return 2
            return 3

        for sys, sys_els in sorted(blade_systems.items()):
            mask = cs_mp["_els"].apply(
                lambda e: bool(e) and (e == {"O"} or (bool(e - {"O"}) and (e - {"O"}).issubset(sys_els)))
            )
            sys_df = cs_mp[mask].drop(columns=["_els"]).copy()
            sys_df = sys_df.drop_duplicates(subset="formula", keep="first")
            if "is_stable" in sys_df.columns:
                sys_df = sys_df[
                    sys_df["formula"].apply(lambda f: len(get_formula_elements(str(f))) == 1)
                    | sys_df["is_stable"].apply(
                        lambda v: True
                        if (pd.isna(v) or str(v).strip() in ("", "nan", "None"))
                        else str(v).strip().lower() in {"true", "1", "yes"}
                    )
                ].copy()

            sys_df["_ord"] = sys_df["formula"].apply(_sys_row_order)
            sys_df = (
                sys_df.sort_values(["_ord", "formula"])
                .drop(columns=["_ord"])
                .reset_index(drop=True)
            )
            if sys_df.empty:
                continue
            sys_df.to_excel(sys_dir / f"{str(sys).replace('/', '_')}.xlsx", index=False)

        print(f"Saved {len(blade_systems)} system files → {sys_dir}")

    # ------------------------------------------------------------------
    # Section 5: Per-parent Excel files
    # ------------------------------------------------------------------

    def _build_parent_excel(self) -> None:
        out_root = self.mp_comps_dir
        out_root.mkdir(parents=True, exist_ok=True)

        # Build POSCAR index — prefer dedicated POSCARs folder (original MP structures)
        _poscar_scan_root = self.files_dir / self.poscar_folder
        if not _poscar_scan_root.exists():
            _poscar_scan_root = out_root  # fallback to GRACE folder
        print(f"Indexing POSCARs in {_poscar_scan_root} for section 6 ...")
        self._s6_poscar_index = {}
        for p in _poscar_scan_root.rglob("POSCAR"):
            key = p.parent.name.split("_")[0]
            if key not in self._s6_poscar_index:
                self._s6_poscar_index[key] = p
        print(f"  Found {len(self._s6_poscar_index)} POSCARs")

        # Build CONTCAR index from all MLIP dirs
        self._s6_contcar_index = {}
        for lbl, mdir in self.mlip_dirs.items():
            for c in mdir.rglob("CONTCAR"):
                key = c.parent.name.split("_")[0]
                if key not in self._s6_contcar_index:
                    self._s6_contcar_index[key] = c
        print(f"  Found {len(self._s6_contcar_index)} CONTCARs")

        blade_df = pd.read_excel(self.files_dir / "blade_generated_data.xlsx")
        mp_df = pd.read_excel(self.mp_comps_dir / "materials_project_summary.xlsx")

        blade_df.columns = blade_df.columns.str.strip()
        mp_df.columns = mp_df.columns.str.strip()

        formula_col = "parent_formula" if "parent_formula" in blade_df.columns else "formula"
        blade_df = blade_df.drop_duplicates(subset=[formula_col], keep="first").copy()

        _meta_cols_s6 = {
            "formula", "phase_label", "energy", "dH", "entry_id", "file_path",
            "is_stable", "source", "volume", "parent_formula", "parent_entry_id",
        }
        all_element_order = sorted(
            {c for c in blade_df.columns if c not in _meta_cols_s6 and _is_element_col(c)} | {"O"},
            key=lambda e: _PMGElement(e).Z,
        )
        metal_elements = [el for el in all_element_order if el not in self.fixed_elements | {"O"}]

        base_cols_after_elements = ["energy", "dH", "entry_id", "file_path", "is_stable", "source", "volume"]

        output_dir = self.files_dir / "parent_formula_excel_files"
        source_dirs = {"MP": output_dir / "MP"}
        for lbl in self.mlip_dirs:
            source_dirs[lbl] = output_dir / lbl.replace("+", "_")
        for d in source_dirs.values():
            d.mkdir(parents=True, exist_ok=True)

        mp_df["formula"] = mp_df["formula"].astype(str)
        mp_oxides_df = mp_df[mp_df["formula"].apply(is_oxide_formula)].copy()
        if "is_stable" in mp_oxides_df.columns:
            mp_oxides_df = mp_oxides_df[mp_oxides_df["is_stable"].apply(is_true_stable)].copy()

        def get_blade_allowed_elements(row: pd.Series) -> set[str]:
            allowed: set[str] = set()
            for el in metal_elements:
                if el in row.index and pd.notna(row[el]) and float(row[el]) > 0:
                    allowed.add(el)
            for _fe in self.fixed_elements:
                if _fe in row.index and pd.notna(row[_fe]) and float(row[_fe]) > 0:
                    allowed.add(_fe)
            return allowed

        def get_value(row: pd.Series, possible_cols: list[str], default=""):
            for col in possible_cols:
                if col in row.index and pd.notna(row[col]):
                    return row[col]
            return default

        def get_mp_structure_path(row: pd.Series, poscar_only: bool = False) -> str:
            entry_id = str(get_value(row, ["entry_id", "material_id"], "")).strip()
            if not poscar_only and entry_id and entry_id.lower() != "nan":
                hit = self._s6_contcar_index.get(entry_id)
                if hit:
                    return str(hit.resolve())
            if entry_id and entry_id.lower() != "nan":
                hit = self._s6_poscar_index.get(entry_id)
                if hit:
                    return str(hit.resolve())
            raw_path = get_value(row, ["file_path", "poscar_path"], "")
            if raw_path:
                p = Path(str(raw_path))
                if p.exists():
                    return str(p.resolve())
            return ""

        def standardize_blade_row(row: pd.Series) -> dict:
            formula = str(get_value(row, ["formula", "phase_label"], "")).strip()
            new: dict = {"formula": formula}
            for el in all_element_order:
                new[el] = row[el] if el in row.index and pd.notna(row[el]) else get_formula_counts(formula).get(el, 0)
            new["energy"] = get_value(row, ["energy", "energy_per_atom"], "")
            new["dH"] = get_value(row, ["dH", "formation_energy_per_atom"], "")
            new["entry_id"] = get_value(row, ["entry_id", "material_id"], "")
            new["file_path"] = get_value(row, ["file_path"], "")
            new["is_stable"] = True
            new["source"] = get_value(row, ["source"], "BLADE")
            new["volume"] = get_value(row, ["volume"], "")
            return new

        def standardize_mp_row(row: pd.Series) -> dict:
            formula = str(get_value(row, ["formula", "phase_label"], "")).strip()
            counts = get_formula_counts(formula)
            new: dict = {"formula": formula}
            for el in all_element_order:
                new[el] = row[el] if el in row.index and pd.notna(row[el]) else counts.get(el, 0)
            new["energy"] = get_value(row, ["energy", "energy_per_atom"], "")
            new["dH"] = get_value(row, ["dH", "formation_energy_per_atom"], "")
            new["entry_id"] = get_value(row, ["entry_id", "material_id"], "")
            new["file_path"] = get_mp_structure_path(row)
            new["is_stable"] = is_true_stable(get_value(row, ["is_stable"], True))
            new["source"] = get_value(row, ["source"], "Materials Project")
            new["volume"] = get_value(row, ["volume"], "")
            return new

        def standardize_pure_element_row(row: pd.Series, element: str) -> dict:
            new = standardize_mp_row(row)
            new["formula"] = element
            new["energy"] = self.fallback_refs.get(element, new["energy"])
            new["dH"] = 0
            new["is_stable"] = True
            for el in all_element_order:
                new[el] = 1 if el == element else 0
            if new["source"] == "":
                new["source"] = "Materials Project"
            return new

        def make_blank_pure_element_row(element: str) -> dict:
            new: dict = {"formula": element}
            for el in all_element_order:
                new[el] = 1 if el == element else 0
            new.update({
                "energy": self.fallback_refs.get(element, ""),
                "dH": 0,
                "entry_id": "",
                "file_path": "",
                "is_stable": True,
                "source": "Materials Project",
                "volume": "",
            })
            return new

        def is_metal_only_formula(formula: str, allowed_els: set[str]) -> bool:
            try:
                els = set(Composition(str(formula)).get_el_amt_dict().keys())
                return len(els) > 1 and not (els & (self.fixed_elements | {"O"})) and els.issubset(allowed_els)
            except Exception:
                return False

        def _is_pure_element(formula: str) -> bool:
            try:
                return len(Composition(formula).get_el_amt_dict()) == 1
            except Exception:
                return False

        def _apply_mlip_values(df: pd.DataFrame, src_label: str,
                               stable_only: bool = True) -> pd.DataFrame:
            """Replace energy/dH/file_path/source for non-BLADE rows using MLIP lookup.
            Pure element rows always keep dH=0 (they are the reference state)."""
            lkp = (self.mlip_lookups if stable_only else self._mlip_lookups_all).get(src_label, {})
            if not lkp:
                return df
            df = df.copy()
            for idx, row in df.iterrows():
                if str(row.get("source", "")) == "BLADE":
                    continue
                fkey = normalize_formula(str(row.get("formula", "")))
                if fkey in lkp:
                    df.at[idx, "energy"] = lkp[fkey].get("energy")
                    df.at[idx, "file_path"] = lkp[fkey].get("file_path", row.get("file_path", ""))
                    df.at[idx, "source"] = src_label
                    if not _is_pure_element(str(row.get("formula", ""))):
                        df.at[idx, "dH"] = lkp[fkey].get("dH")
            return df

        def _fix_mp_paths(df: pd.DataFrame) -> pd.DataFrame:
            """For MP DFT subfolder: replace CONTCAR paths with POSCAR paths."""
            df = df.copy()
            for idx, row in df.iterrows():
                eid = str(row.get("entry_id", "")).strip()
                if eid and eid.lower() != "nan":
                    hit = self._s6_poscar_index.get(eid)
                    if hit:
                        df.at[idx, "file_path"] = str(hit.resolve())
            return df

        def _filter_stable(df: pd.DataFrame) -> pd.DataFrame:
            """Keep only stable MP rows + all BLADE rows."""
            return df[
                (df["is_stable"].astype(str).str.upper() == "TRUE")
                | (df["source"].astype(str).str.strip() == "BLADE")
            ].copy()

        def _row_order(formula: str) -> int:
            els = get_formula_elements(formula)
            if len(els) == 1:
                return 0  # pure element
            if (els & self.fixed_elements) and "O" not in els:
                return 1  # fixed-element compound (boride/carbide/etc.)
            if "O" in els:
                return 2  # oxide
            return 3  # metal compound

        # Build pure element row lookup from mp_df
        pure_element_rows: dict[str, pd.Series] = {}
        for el in all_element_order:
            el_rows = mp_df[mp_df["formula"].apply(lambda f: is_single_element_formula(f, el))].copy()
            if el_rows.empty:
                continue
            if "is_stable" in el_rows.columns:
                stable = el_rows[el_rows["is_stable"].apply(is_true_stable)]
                if not stable.empty:
                    el_rows = stable
            if "energy_above_hull" in el_rows.columns:
                el_rows["energy_above_hull"] = pd.to_numeric(el_rows["energy_above_hull"], errors="coerce")
                el_rows = el_rows.sort_values("energy_above_hull", na_position="last")
            pure_element_rows[el] = el_rows.iloc[0]

        # Print element ref comparison table
        print(f"\n{'Element':<8} {'Hardcoded DFT (eV/atom)':>24} {'MP energy_per_atom (eV/atom)':>30} {'Diff':>10}")
        print("-" * 76)
        for el in sorted(self.element_refs):
            hard = self.fallback_refs.get(el)
            if hard is None:
                continue
            mp_row = pure_element_rows.get(el)
            mp_e = (
                float(mp_row["energy_per_atom"])
                if mp_row is not None and pd.notna(mp_row.get("energy_per_atom"))
                else None
            )
            diff = (mp_e - hard) if mp_e is not None else None
            mp_str = f"{mp_e:.6f}" if mp_e is not None else "N/A"
            diff_str = f"{diff:+.6f}" if diff is not None else "N/A"
            print(f"{el:<8} {hard:>24.6f} {mp_str:>30} {diff_str:>10}")

        missing_pure_elements: set[str] = set()
        missing_mp_contcars: list[dict] = []

        for _, blade_row in blade_df.iterrows():
            allowed = get_blade_allowed_elements(blade_row)
            if not allowed:
                continue

            required_pure = set(allowed) | {"O"}
            parent_formula = str(get_value(blade_row, ["parent_formula", "formula", "phase_label"], "parent"))

            rows: list[dict] = []

            for el in all_element_order:
                if el not in required_pure:
                    continue
                if el in pure_element_rows:
                    rows.append(standardize_pure_element_row(pure_element_rows[el], el))
                else:
                    rows.append(make_blank_pure_element_row(el))
                    missing_pure_elements.add(el)

            rows.append(standardize_blade_row(blade_row))

            # Include all BLADE sub-compositions whose metal elements are a strict subset
            parent_formula_val = str(get_value(blade_row, ["formula", "phase_label"], ""))
            for _, sub_row in blade_df.iterrows():
                sub_formula = str(get_value(sub_row, ["formula", "phase_label"], ""))
                if sub_formula == parent_formula_val:
                    continue
                sub_allowed = get_blade_allowed_elements(sub_row)
                if sub_allowed and sub_allowed < allowed:
                    rows.append(standardize_blade_row(sub_row))

            # Stable MP oxides
            for _, oxide_row in mp_oxides_df[
                mp_oxides_df["formula"].apply(lambda f: oxide_allowed_for_elements(f, allowed))
            ].iterrows():
                rows.append(standardize_mp_row(oxide_row))

            # Non-oxide non-boride MP compounds (metal alloys) — stable only
            metal_allowed = {el for el in allowed if el in metal_elements}
            alloy_mask = mp_df["formula"].apply(lambda f: is_metal_only_formula(f, metal_allowed))
            alloy_df = mp_df[alloy_mask]
            if "is_stable" in alloy_df.columns:
                alloy_df = alloy_df[alloy_df["is_stable"].apply(is_true_stable)]
            for _, alloy_row in alloy_df.iterrows():
                rows.append(standardize_mp_row(alloy_row))

            out = pd.DataFrame(rows)
            for el in all_element_order:
                out[el] = pd.to_numeric(out[el], errors="coerce").fillna(0)

            element_cols_to_keep = [
                el for el in all_element_order
                if el in out.columns and out[el].sum() != 0
            ]
            final_cols = ["formula", *element_cols_to_keep, *base_cols_after_elements]
            for col in final_cols:
                if col not in out.columns:
                    out[col] = ""
            out = out[final_cols].copy()

            # Dedup by formula+source: sort energy → volume → entry_id number, keep best
            def _eid_num(eid):
                try:
                    return int(str(eid).split("-")[-1])
                except Exception:
                    return 999_999_999

            out["_sort_stable"] = out["is_stable"].apply(
                lambda v: 0 if is_true_stable(v) else 1)
            out["_sort_vol"] = pd.to_numeric(out.get("volume", pd.Series(dtype=float)), errors="coerce")
            out["_sort_eid"] = out["entry_id"].apply(_eid_num)
            out["_sort_e"] = pd.to_numeric(out["energy"], errors="coerce")
            out = (
                out.sort_values(
                    ["_sort_stable", "_sort_e", "_sort_vol", "_sort_eid"],
                    na_position="last"
                )
                .drop_duplicates(subset=["formula", "source"], keep="first")
                .drop(columns=["_sort_stable", "_sort_vol", "_sort_eid", "_sort_e"])
                .reset_index(drop=True)
            )
            out["is_stable"] = out["is_stable"].apply(
                lambda x: ""
                if (pd.isna(x) or str(x).strip() in ("", "nan", "None"))
                else "TRUE" if is_true_stable(x) else "FALSE"
            )

            out["_order"] = out["formula"].apply(_row_order)
            out = out.sort_values(["_order", "formula"]).drop(columns=["_order"]).reset_index(drop=True)

            mp_missing = out[
                (out["source"].astype(str).str.strip() == "Materials Project")
                & (out["file_path"].astype(str).str.strip() == "")
            ]
            for _, mr in mp_missing.iterrows():
                missing_mp_contcars.append({
                    "parent_formula": parent_formula,
                    "phase_label": mr.get("phase_label", ""),
                    "formula": mr.get("formula", ""),
                    "entry_id": mr.get("entry_id", ""),
                })

            fname = f"{clean_file_name(parent_formula)}.xlsx"
            for src_label, src_dir in source_dirs.items():
                src_out = _filter_stable(out.copy())
                if src_label == "MP":
                    src_out = _fix_mp_paths(src_out)
                else:
                    src_out = _apply_mlip_values(src_out, src_label)
                if not src_out.empty:
                    src_out.to_excel(src_dir / fname, index=False)

        print(f"Done. Wrote files to: {output_dir.resolve()}")
        print(f"Unique BLADE formulas written: {len(blade_df)}")

        if missing_pure_elements:
            print(f"Warning: Missing pure element refs (blank rows added): {sorted(missing_pure_elements)}")

        if missing_mp_contcars:
            missing_path = output_dir / "missing_mp_contcar_paths.xlsx"
            pd.DataFrame(missing_mp_contcars).to_excel(missing_path, index=False)
            print(f"Warning: Some MP rows have blank file_path. Report: {missing_path.resolve()}")
        else:
            print("All Materials Project rows have CONTCAR paths.")

        # Store for _build_pure_element_energies
        self._pure_element_rows = pure_element_rows
        self._all_element_order = all_element_order
        self._mp_df_s6 = mp_df

    # ------------------------------------------------------------------
    # Section 6: Pure-element energy table
    # ------------------------------------------------------------------

    def _build_pure_element_energies(self) -> None:
        all_element_order = self._all_element_order
        pure_element_rows = self._pure_element_rows
        mp_df = self._mp_df_s6

        all_ref_elements = sorted(set(list(all_element_order) + ["O"]))
        el_ref_rows = []
        for el in all_ref_elements:
            mp_er = pure_element_rows.get(el)
            if mp_er is None:
                mp_all = mp_df[mp_df["formula"].apply(lambda f: is_single_element_formula(f, el))].copy()
                if not mp_all.empty:
                    if "energy_above_hull" in mp_all.columns:
                        mp_all["energy_above_hull"] = pd.to_numeric(mp_all["energy_above_hull"], errors="coerce")
                        mp_all = mp_all.sort_values("energy_above_hull", na_position="last")
                    mp_er = mp_all.iloc[0]
            mid = mp_er["material_id"] if mp_er is not None else ""
            # Use hardcoded fallback refs for MP DFT energy (consistent reference)
            mp_e = self.fallback_refs.get(el)
            row = {"element": el, "material_id_mp_dft": mid, "energy_mp_dft": mp_e}
            for lbl, refs in self._all_element_refs.items():
                col = lbl.replace("+", "_").replace(" ", "_").lower()
                row[f"energy_{col}"] = refs.get(el)
            el_ref_rows.append(row)

        el_ref_df = pd.DataFrame(el_ref_rows)
        el_ref_xlsx = self.files_dir / "pure_element_energies.xlsx"
        el_ref_df.to_excel(el_ref_xlsx, index=False)
        print(f"Saved pure element energies: {el_ref_xlsx}")
