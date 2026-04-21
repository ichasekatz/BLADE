"""

"""

from pathlib import Path
import shutil
import os
import subprocess
import sys
from itertools import permutations
from math import prod

import numpy as np
from materialsframework.calculators import GraceCalculator as Calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar


class BladeTDBGen:
    """

    """
    def __init__(self, phases, liquid, paths, composition_list, level, skip_existing=False):
        """

        """
        self.phases = phases
        self.liquid = liquid
        self.path0 = Path(paths[0])
        self.path2 = Path(paths[1])
        self.composition_list = composition_list
        self.level = level
        self.skip_existing = skip_existing

    @staticmethod
    def normalize_fractions(fractions, n_elements):
        fractions = list(fractions)
        if len(fractions) > n_elements:
            return None
        while len(fractions) < n_elements:
            fractions.append(0.0)
        return fractions

    @staticmethod
    def expand_all_cases(levels, elements):
        cases = []
        seen = set()

        for entry in levels:
            level = entry["level"]

            for raw_fractions in entry["compositions"]:
                fractions = BladeTDBGen.normalize_fractions(raw_fractions, len(elements))

                if fractions is None:
                    print(
                        f"Skipping level {level}: too many fractions "
                        f"{raw_fractions} for {len(elements)} elements"
                    )
                    continue

                # add all unique permutations of the composition
                for perm in set(permutations(fractions)):
                    perm = list(perm)
                    key = (level, tuple(perm))

                    if key in seen:
                        continue

                    seen.add(key)
                    cases.append(
                        {
                            "level": level,
                            "fractions": perm,
                        }
                    )

        # optional: stable ordering for readability
        cases.sort(key=lambda x: (x["level"], x["fractions"]))
        return cases

    @staticmethod
    def composition_string(elements, fractions):
        return "".join(f"{el}{frac:g}" for el, frac in zip(elements, fractions))

    @staticmethod
    def filename_string(elements, fractions):
        return "".join(f"{el}{str(frac).replace('.', 'p')}" for el, frac in zip(elements, fractions))

    @staticmethod
    def poscar_filename(elements, fractions, phase_name):
        return f"POSCAR_{BladeTDBGen.filename_string(elements, fractions)}_{phase_name.upper()}"

    @staticmethod
    def folder_name(elements, fractions, level):
        parts = [f"sqs_lev={level}"]
        comp_parts = [f"{el}={frac:g}" for el, frac in zip(elements, fractions)]
        parts.append("a_" + ",a_".join(comp_parts))
        return "_".join(parts)

    @staticmethod
    def is_pure_endmember(fractions, tol=1e-12):
        nonzero = sum(1 for x in fractions if abs(x) > tol)
        return nonzero <= 1

    @staticmethod
    def get_phase_supercell_size(phase, default_supercell_size=(2, 2, 2)):
        return tuple(phase.get("supercell_size", default_supercell_size))

    @staticmethod
    def infer_phase_multiplicities_from_generated_poscars(
        phases,
        cases,
        elements,
        workdir,
        default_supercell_size=(2, 2, 2),
    ):
        phase_multiplicities = {}

        for phase in phases:
            found = False
            phase_supercell_size = BladeTDBGen.get_phase_supercell_size(phase, default_supercell_size)
            n_super = prod(phase_supercell_size)

            for case in cases:
                poscar_path = workdir / BladeTDBGen.poscar_filename(
                    elements,
                    case["fractions"],
                    phase["generator_name"],
                )

                if not poscar_path.exists():
                    continue

                structure = Poscar.from_file(poscar_path).structure
                n_sites = len(structure)

                if n_sites % n_super != 0:
                    raise ValueError(
                        f"Cannot infer mult for {phase['lattice']} from {poscar_path}: "
                        f"{n_sites} sites is not divisible by supercell product {n_super} "
                        f"for supercell_size={phase_supercell_size}"
                    )

                mult = n_sites // n_super
                phase_multiplicities[phase["lattice"]] = mult
                print(
                    f"Using generated POSCAR for {phase['lattice']}: "
                    f"{poscar_path} ({n_sites} sites -> mult={mult}, "
                    f"supercell_size={phase_supercell_size})"
                )
                found = True
                break

            if not found:
                raise FileNotFoundError(
                    f"Could not infer mult for {phase['lattice']} from generated POSCAR files"
                )

        return phase_multiplicities

    @staticmethod
    def write_species_in(base_dir, elements, phases):
        content = f"a={','.join(elements)}\n"
        for phase in phases:
            lattice_dir = base_dir / phase["lattice"]
            lattice_dir.mkdir(parents=True, exist_ok=True)
            (lattice_dir / "species.in").write_text(content)

    @staticmethod
    def write_mult_in(base_dir, phase_multiplicities):
        for lattice, mult in phase_multiplicities.items():
            lattice_dir = base_dir / lattice
            lattice_dir.mkdir(parents=True, exist_ok=True)
            (lattice_dir / "mult.in").write_text(f"a={mult}\n")

    @staticmethod
    def write_terms_in(base_dir, phases, terms=None):
        for phase in phases:
            lattice = phase["lattice"]
            lattice_dir = base_dir / lattice
            lattice_dir.mkdir(parents=True, exist_ok=True)

            if terms is not None:
                content = terms
            else:
                if lattice in ["BCC_A2", "FCC_A1", "HCP_A3"]:
                    content = "1,0\n2,1\n"
                else:
                    content = "1,0:1,0\n2,0:1,0\n"

            (lattice_dir / "terms.in").write_text(content)

    @staticmethod
    def parse_atat_str(path):
        lines = [line.rstrip() for line in Path(path).read_text().splitlines() if line.strip()]
        if len(lines) < 7:
            raise ValueError(f"{path} does not look like a valid ATAT str.out file")

        coord_sys = np.array([[float(x) for x in lines[i].split()] for i in range(3)], dtype=float)
        supercell = np.array([[float(x) for x in lines[i].split()] for i in range(3, 6)], dtype=float)

        species = []
        coords = []
        for line in lines[6:]:
            parts = line.split()
            if len(parts) < 4:
                continue
            coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
            species.append(parts[3])

        coords = np.array(coords, dtype=float)
        return coord_sys, supercell, coords, species

    @staticmethod
    def write_atat_str(path, coord_sys, supercell, coords, species):
        with Path(path).open("w") as f:
            for row in coord_sys:
                f.write(" ".join(f"{x:.6f}" for x in row) + "\n")
            for row in supercell:
                f.write(" ".join(f"{x:.6f}" for x in row) + "\n")
            for xyz, sp in zip(coords, species):
                f.write(" ".join(f"{x:.6f}" for x in xyz) + f" {sp}\n")

    @staticmethod
    def wrap_frac_near_reference(frac, ref_frac):
        return frac - np.round(frac - ref_frac)

    @staticmethod
    def write_relaxed_str_out_from_template(template_str, contcar, output_str):
        coord_sys, supercell, template_coords, template_species = BladeTDBGen.parse_atat_str(template_str)
        relaxed = Structure.from_file(contcar)

        relaxed_species = [site.species_string for site in relaxed]
        if len(relaxed_species) != len(template_species):
            raise ValueError("Template str.out and CONTCAR have different numbers of atoms")

        if relaxed_species != template_species:
            raise ValueError("Species order in CONTCAR does not match template str.out")

        supercell_inv = np.linalg.inv(supercell)
        template_frac = template_coords @ supercell_inv

        relaxed_frac = np.array(relaxed.frac_coords, dtype=float)
        relaxed_frac_unwrapped = np.array(
            [BladeTDBGen.wrap_frac_near_reference(f, r) for f, r in zip(relaxed_frac, template_frac)],
            dtype=float,
        )

        new_coords = relaxed_frac_unwrapped @ supercell
        BladeTDBGen.write_atat_str(output_str, coord_sys, supercell, new_coords, template_species)

    @staticmethod
    def calculate_all_structures(s2t, cases, elements, phases, workdir):
        for case in cases:
            level = case["level"]
            fractions = case["fractions"]
            comp_name = BladeTDBGen.folder_name(elements, fractions, level)

            for phase in phases:
                folder = workdir / phase["lattice"] / comp_name

                if not (folder / "POSCAR").exists():
                    print(f"Skipping missing POSCAR in {folder}")
                    continue

                try:
                    print(f"\nCalculating {folder}")
                    s2t._calculate(folder)

                    template_str = folder / "str_template.out"
                    contcar = folder / "CONTCAR"
                    str_out = folder / "str.out"

                    if template_str.exists() and contcar.exists():
                        BladeTDBGen.write_relaxed_str_out_from_template(template_str, contcar, str_out)
                        print(f"Wrote str.out in {folder}")
                    else:
                        print(f"Missing template or CONTCAR in {folder}")

                    print("Files now in folder:")
                    print(sorted(p.name for p in folder.iterdir()))

                except Exception as exc:
                    print(f"Failed calculation for {folder}: {exc}")

    @staticmethod
    def check_required_files(cases, elements, phases, workdir):
        required = ["POSCAR", "energy", "CONTCAR", "str.out", "force.out", "stress.out"]
        missing_any = False

        for case in cases:
            level = case["level"]
            fractions = case["fractions"]
            comp_name = BladeTDBGen.folder_name(elements, fractions, level)

            for phase in phases:
                folder = workdir / phase["lattice"] / comp_name
                missing = [name for name in required if not (folder / name).exists()]

                if missing:
                    missing_any = True
                    print(f"Missing files in {folder}: {', '.join(missing)}")

        return not missing_any

    @staticmethod
    def run_fit_model(elements, phases, sqsgen_levels, params):
        fit_cmd = (
            "from materialsframework.tools.sqs2tdb import Sqs2tdb; "
            "from materialsframework.calculators import GraceCalculator as Calculator; "
            f"s=Sqs2tdb(fmax={params['fmax']}, verbose={params['verbose']}, calculator=Calculator(device='{params['calculator']}')); "
            f"s.species={elements!r}; "
            f"s.lattices={[phase['lattice'] for phase in phases]!r}; "
            f"s.level={max(entry['level'] for entry in sqsgen_levels)!r}; "
            f"s.t_min={params['t_min']}; "
            f"s.t_max={params['t_max']}; "
            f"s.sro={params['sro']}; "
            f"s.bv={params['bv']}; "
            f"s.phonon={params['phonon']}; "
            f"s.open_calphad={params['open_calphad']}; "
            f"s.terms={params['terms']}; "
            "s._fit_model()"
        )

        result = subprocess.run(
            [sys.executable, "-c", fit_cmd],
            capture_output=True,
            text=True,
        )

        print("\n_fit_model return code:", result.returncode)
        print("FIT STDOUT:")
        print(result.stdout)
        print("FIT STDERR:")
        print(result.stderr)

        return result.returncode == 0

    @staticmethod
    def _sqs_source_case_name(fractions, level):
        labels = [chr(ord("a") + i) for i in range(len(fractions))]
        comp_parts = [f"{labels[i]}={fractions[i]:g}" for i in range(len(fractions))]
        return f"sqs_lev={level}_" + ",".join(comp_parts)

    @staticmethod
    def _replace_variable_labels_in_text(text, elements):
        labels = [chr(ord("a") + i) for i in range(len(elements))]
        tmp_tokens = {}

        for i, lbl in enumerate(labels):
            tmp = f"__TMPLBL{i}__"
            tmp_tokens[lbl] = tmp
            text = text.replace(f" {lbl}\n", f" {tmp}\n")
            text = text.replace(f" {lbl} ", f" {tmp} ")
            text = text.replace(f"{lbl}=", f"{tmp}=")

        for i, lbl in enumerate(labels):
            el = elements[i]
            tmp = tmp_tokens[lbl]
            text = text.replace(f" {tmp}\n", f" {el}\n")
            text = text.replace(f" {tmp} ", f" {el} ")
            text = text.replace(f"{tmp}=", f"{el}=")

        return text

    @staticmethod
    def _rewrite_folder_labels_to_elements(folder, elements):
        label_to_element = {chr(ord("a") + i): elements[i] for i in range(len(elements))}

        for fname in ["POSCAR", "str_template.out", "str.out", "species.in"]:
            path = folder / fname
            if not path.exists():
                continue

            text = path.read_text()
            text = BladeTDBGen._replace_variable_labels_in_text(text, elements)

            if fname == "POSCAR":
                lines = text.splitlines()

                if len(lines) >= 8:
                    atom_lines = lines[8:]

                    grouped = {}
                    species_order = []

                    for line in atom_lines:
                        parts = line.split()
                        if len(parts) < 4:
                            continue

                        coords = parts[:-1]
                        sp = parts[-1]
                        sp = label_to_element.get(sp, sp)

                        if sp not in grouped:
                            grouped[sp] = []
                            species_order.append(sp)

                        grouped[sp].append(" ".join(coords + [sp]))

                    new_atom_lines = []
                    for sp in species_order:
                        new_atom_lines.extend(grouped[sp])

                    lines[5] = " ".join(species_order)
                    lines[6] = " ".join(str(len(grouped[sp])) for sp in species_order)
                    lines = lines[:8] + new_atom_lines

                    text = "\n".join(lines) + "\n"

            path.write_text(text)

    def copy_sqs_folders_into_workdir(self, cases, elements, phases, workdir):
        sqs_root = self.path2 / "SQS"

        for phase in phases:
            lattice = phase["lattice"]
            n_elements = len(elements)

            src_phase_dir = sqs_root / f"{lattice}_{n_elements}" / lattice
            dst_phase_dir = workdir / lattice
            dst_phase_dir.mkdir(parents=True, exist_ok=True)

            if not src_phase_dir.exists():
                print(f"Missing SQS source phase directory: {src_phase_dir}")
                continue

            for case in cases:
                level = case["level"]
                fractions = case["fractions"]

                src_case_name = BladeTDBGen._sqs_source_case_name(fractions, level)
                src = src_phase_dir / src_case_name

                dst_case_name = BladeTDBGen.folder_name(elements, fractions, level)
                dst = dst_phase_dir / dst_case_name

                if not src.exists():
                    print(f"Missing SQS source folder: {src}")
                    continue

                if not dst.exists():
                    shutil.copytree(src, dst)

                BladeTDBGen._rewrite_folder_labels_to_elements(dst, elements)

                # also copy one top-level POSCAR for multiplicity inference
                src_poscar = src / "POSCAR"
                if src_poscar.exists():
                    top_poscar = workdir / BladeTDBGen.poscar_filename(
                        elements,
                        fractions,
                        phase["generator_name"],
                    )
                    shutil.copy2(src_poscar, top_poscar)

                    text = top_poscar.read_text()
                    text = BladeTDBGen._replace_variable_labels_in_text(text, elements)

                    lines = text.splitlines()
                    if len(lines) >= 6:
                        counts = lines[6].split()
                        lines[5] = " ".join(elements[:len(counts)])
                        text = "\n".join(lines) + "\n"

                    top_poscar.write_text(text)

    def fit_tdb(self, cases, elements, phases, workdir, sqsgen_levels, params, default_supercell_size=(2, 2, 2)):
        self.copy_sqs_folders_into_workdir(cases, elements, phases, workdir)

        BladeTDBGen.write_species_in(workdir, elements, phases)

        phase_multiplicities = BladeTDBGen.infer_phase_multiplicities_from_generated_poscars(
            phases=phases,
            cases=cases,
            elements=elements,
            workdir=workdir,
            default_supercell_size=default_supercell_size,
        )
        BladeTDBGen.write_mult_in(workdir, phase_multiplicities)
        BladeTDBGen.write_terms_in(workdir, phases, terms=None)

        s2t = Sqs2tdb(
            fmax=params['fmax'],
            verbose=params['verbose'],
            calculator=Calculator(device=params['calculator']),
        )

        BladeTDBGen.calculate_all_structures(s2t, cases, elements, phases, workdir)

        fit_cases = [case for case in cases if not BladeTDBGen.is_pure_endmember(case["fractions"])]

        print("\nChecking required files...")
        all_good = BladeTDBGen.check_required_files(fit_cases, elements, phases, workdir)

        if not all_good:
            print("\nAborting fit because some folders are missing required files.")
            return

        print("\nPhase multiplicities:")
        for lattice, mult in phase_multiplicities.items():
            print(f"{lattice}: a={mult}")

        print("\nFolders included in fit:")
        for case in fit_cases:
            level = case["level"]
            fractions = case["fractions"]
            comp_name = BladeTDBGen.folder_name(elements, fractions, level)
            for phase in phases:
                folder = workdir / phase["lattice"] / comp_name
                print(folder)

        print("\nFitting CALPHAD model...")

        s2t.species = elements
        s2t.lattices = [phase["lattice"] for phase in phases]
        s2t.level = max(entry["level"] for entry in sqsgen_levels)
        s2t.t_min = params['t_min']
        s2t.t_max = params['t_max']
        s2t.sro = params['sro']
        s2t.bv = params['bv']
        s2t.phonon = params['phonon']
        s2t.open_calphad = params['open_calphad']
        s2t.terms = params['terms']

        original_dir = Path.cwd()
        workdir = workdir.resolve()
        tdb_name = "_".join(sorted(el.upper() for el in elements)) + ".tdb"
        tdb_path = workdir / tdb_name
        moved = []

        try:
            os.chdir(workdir)

            for case in cases:
                if BladeTDBGen.is_pure_endmember(case["fractions"]):
                    level = case["level"]
                    comp_name = BladeTDBGen.folder_name(elements, case["fractions"], level)

                    for phase in phases:
                        src = Path(phase["lattice"]) / comp_name
                        if src.exists():
                            dst = Path(phase["lattice"]) / f"__skip__{comp_name}"
                            src.rename(dst)
                            moved.append((dst, src))

            fit_success = BladeTDBGen.run_fit_model(elements, phases, sqsgen_levels, params)

            if not fit_success:
                print("\nSkipping -tdb because _fit_model() failed.")
                return

            s2t._run_command("sqs2tdb", ["-tdb"])

            print(f"Expected TDB path: {tdb_path}")
            print(f"TDB exists: {tdb_path.exists()}")
            print("Files in workdir after -tdb:")
            print(sorted(p.name for p in Path('.').iterdir()))

            if not tdb_path.exists():
                print("\nSkipping success message because TDB file was not created.")
                return

        except Exception as exc:
            print(f"Failed to fit TDB: {exc}")

        finally:
            for dst, src in moved:
                if dst.exists():
                    dst.rename(src)
            os.chdir(original_dir)

        if tdb_path.exists():
            print(f"\nGenerated TDB: {tdb_path}")
        else:
            print(f"\nTDB was not created: {tdb_path}")

    def run_workflow(
        self,
        sqsgen_levels,
        elements,
        phase_dicts,
        workdir,
        params,
        default_supercell_size=(2, 2, 2),
    ):
        workdir.mkdir(parents=True, exist_ok=True)

        cases = BladeTDBGen.expand_all_cases(sqsgen_levels, elements)

        print("\nExpanded compositions:")
        for case in cases:
            print(
                f"level={case['level']} "
                f"fractions={case['fractions']} "
                f"name={BladeTDBGen.composition_string(elements, case['fractions'])}"
            )

        self.fit_tdb(
            cases=cases,
            elements=elements,
            phases=phase_dicts,
            workdir=workdir,
            sqsgen_levels=sqsgen_levels,
            params=params,
            default_supercell_size=default_supercell_size,
        )

    def get_composition_dir(self, elements):
        return self.path2 / "".join(elements)

    @staticmethod
    def get_tdb_name(elements):
        return "_".join(sorted(el.upper() for el in elements)) + ".tdb"

    def get_tdb_path(self, elements):
        return self.get_composition_dir(elements) / BladeTDBGen.get_tdb_name(elements)

    def should_skip_composition(self, elements):
        tdb_path = self.get_tdb_path(elements)
        return self.skip_existing and tdb_path.exists()

    def run_single_composition(
        self,
        sqsgen_levels,
        elements,
        phase_dicts,
        params,
        default_supercell_size=(2, 2, 2),
    ):
        workdir = self.get_composition_dir(elements)
        workdir.mkdir(parents=True, exist_ok=True)

        tdb_path = self.get_tdb_path(elements)

        if self.should_skip_composition(elements):
            print(f"Skipping {elements} because TDB already exists: {tdb_path}")
            return

        print(f"\nRunning composition: {elements}")
        print(f"Workdir: {workdir}")
        print(f"Expected TDB: {tdb_path}")

        self.run_workflow(
            sqsgen_levels=sqsgen_levels,
            elements=elements,
            phase_dicts=phase_dicts,
            workdir=workdir,
            params=params,
            default_supercell_size=default_supercell_size,
        )

    def run_all_compositions(
        self,
        sqsgen_levels,
        phase_dicts,
        params = {
            "fmax": 1e-4,
            "verbose": True,
            "calculator": Calculator(device="cuda"),
            "t_min": 298.15,
            "t_max": 10000.0,
            "sro": False,
            "bv": 1e-3,
            "phonon": False,
            "open_calphad": False,
            "terms": None,
        },
        default_supercell_size=(2, 2, 2),
    ):
        original_dir = Path.cwd()

        try:
            for comp in self.composition_list:
                elements = list(comp)

                try:
                    self.run_single_composition(
                        sqsgen_levels=sqsgen_levels,
                        elements=elements,
                        phase_dicts=phase_dicts,
                        params=params,
                        default_supercell_size=default_supercell_size,
                    )
                except Exception as exc:
                    print(f"Failed for composition {elements}: {exc}")

        finally:
            os.chdir(original_dir)