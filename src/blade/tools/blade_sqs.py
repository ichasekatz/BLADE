"""

"""

from pathlib import Path
import shutil
from itertools import permutations

import numpy as np
from pymatgen.core import Structure, Composition, Lattice
from pymatgen.io.vasp import Poscar
from sqsgenerator import optimize, parse_config, to_pymatgen
from sqsgenerator.core import LogLevel

class BladeSQS:
    """
    Generate SQS inputs and run ATAT mcsqs for a given phase prototype.

    The phase prototype is specified via a phases_dict containing lattice parameters and a unit-cell
    coordinate string. The class can:
      - build ATAT rndstr inputs (rndstr.skel)
      - write sqsgen.in from preset levels
      - run sqs2tdb -mk to create sqsdb_lev=* folders
      - compute a suitable n_atoms for mcsqs based on target fractions
      - run corrdump + mcsqs with an automatic timeout via stopsqs
    """
    def __init__(self, phases_dict, sqsgen_levels, level):
        self.phases_dict = phases_dict
        self.a = phases_dict["a"]
        self.b = phases_dict["b"]
        self.c = phases_dict["c"]
        self.alpha = phases_dict["alpha"]
        self.beta = phases_dict["beta"]
        self.gamma = phases_dict["gamma"]
        self.vectors = phases_dict["vectors"]
        self.unit_cell = phases_dict["coords"]
        self.sqsgen_levels = sqsgen_levels
        self.level = level

    @staticmethod
    def placeholder_elements(n):
        pool = ["Xe", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru"]
        if n > len(pool):
            raise ValueError(f"Need {n} placeholder elements, but only {len(pool)} are defined.")
        return pool[:n]

    @staticmethod
    def placeholder_letter_map(elements):
        return {el: chr(ord("a") + i) for i, el in enumerate(elements)}

    @staticmethod
    def inverse_placeholder_letter_map(elements):
        return {chr(ord("a") + i): el for i, el in enumerate(elements)}

    @staticmethod
    def lowercase_labels(n):
        if n > 26:
            raise ValueError("More than 26 variable elements are not supported with single-letter lowercase labels.")
        return [chr(ord("a") + i) for i in range(n)]

    @staticmethod
    def relabel_poscar_species(poscar_path, elements):
        """
        Convert placeholder elements back to lowercase labels
        in the created POSCAR files.

        Example:
            Xe -> a
            Kr -> b
            B  -> B   (unchanged because not in placeholder list)
        """
        repl = {el: chr(ord("a") + i) for i, el in enumerate(elements)}

        with open(poscar_path, "r") as f:
            lines = f.readlines()

        if len(lines) < 8:
            raise ValueError(f"POSCAR too short: {poscar_path}")

        # Update header species line, but only for placeholder elements.
        species_line = lines[5].split()
        new_species_line = [repl.get(sp, sp) for sp in species_line]
        lines[5] = " ".join(new_species_line) + "\n"

        # Optional: update comment line in the same spirit.
        comment_parts = lines[0].split()
        updated_comment_parts = []
        for token in comment_parts:
            replaced = token
            for el, lbl in repl.items():
                if token.startswith(el):
                    suffix = token[len(el):]
                    replaced = f"{lbl}{suffix}"
                    break
            updated_comment_parts.append(replaced)
        lines[0] = " ".join(updated_comment_parts) + "\n"

        # Update site labels only if the species is one of the placeholder elements.
        for i in range(8, len(lines)):
            parts = lines[i].split()
            if len(parts) >= 4:
                site = parts[-1]
                if site in repl:
                    parts[-1] = repl[site]
                    lines[i] = " ".join(parts) + "\n"

        with open(poscar_path, "w") as f:
            f.writelines(lines)

    @staticmethod
    def relabel_str_template(template_path, elements):
        repl = {el: chr(ord("a") + i) for i, el in enumerate(elements)}

        with open(template_path, "r") as f:
            lines = f.readlines()

        # First 6 lines are lattice/supercell information
        # Atomic labels start on line 7
        for i in range(6, len(lines)):
            parts = lines[i].split()

            if len(parts) >= 4:
                site = parts[-1]

                # Replace Xe -> a, Kr -> b, etc.
                # Keep B or any other fixed species unchanged
                parts[-1] = repl.get(site, site)

                lines[i] = " ".join(parts) + "\n"

        with open(template_path, "w") as f:
            f.writelines(lines)

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

        for entry in levels:
            level = entry["level"]

            for raw_fractions in entry["compositions"]:
                fractions = BladeSQS.normalize_fractions(raw_fractions, len(elements))

                if fractions is None:
                    print(
                        f"Skipping level {level}: too many fractions "
                        f"{raw_fractions} for {len(elements)} elements"
                    )
                    continue

                cases.append(
                    {
                        "level": level,
                        "fractions": fractions,
                    }
                )

        return cases

    @staticmethod
    def composition_string(elements, fractions):
        return "".join(f"{el}{frac:g}" for el, frac in zip(elements, fractions))

    @staticmethod
    def filename_string(elements, fractions):
        return "".join(
            f"{el}{str(frac).replace('.', 'p')}"
            for el, frac in zip(elements, fractions)
        )

    @staticmethod
    def poscar_filename(elements, fractions, phase_name):
        repl = BladeSQS.placeholder_letter_map(elements)
        parts = []
        for el, frac in zip(elements, fractions):
            parts.append(f"{repl[el]}{str(frac).replace('.', 'p')}")
        return f"POSCAR_{''.join(parts)}_{phase_name.upper()}"

    @staticmethod
    def folder_name(elements, fractions, level):
        repl = BladeSQS.placeholder_letter_map(elements)
        parts = [f"sqs_lev={level}"]
        comp_parts = [f"{repl[el]}={frac:g}" for el, frac in zip(elements, fractions)]
        parts.append(",".join(comp_parts))
        return "_".join(parts)

    @staticmethod
    def get_phase_supercell_size(phase, default_supercell_size=(2, 2, 2)):
        return tuple(phase.get("supercell_size", default_supercell_size))

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
    def build_initial_template_str_from_poscar(poscar_path, output_path):
        structure = Structure.from_file(poscar_path)

        lattice = np.array(structure.lattice.matrix, dtype=float)
        coord_sys = lattice.copy()
        supercell = np.eye(3, dtype=float)

        coords = np.array(structure.frac_coords, dtype=float)
        species = [site.species_string for site in structure]

        BladeSQS.write_atat_str(output_path, coord_sys, supercell, coords, species)

    @staticmethod
    def apply_lowercase_permutation_to_file(file_path, mapping):
        """
        Replace only lowercase variable labels in a file.
        Uppercase fixed species like B remain unchanged.

        mapping example:
            {"a": "b", "b": "a", "c": "c"}
        """
        file_path = Path(file_path)

        with open(file_path, "r") as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            parts = line.split()

            if not parts:
                new_lines.append(line)
                continue

            replaced_any = False
            new_parts = []
            for token in parts:
                new_token = mapping.get(token, token)
                if new_token != token:
                    replaced_any = True
                new_parts.append(new_token)

            if replaced_any:
                new_lines.append(" ".join(new_parts) + "\n")
                continue

            if len(parts) >= 4:
                site = parts[-1]
                if site in mapping:
                    parts[-1] = mapping[site]
                    new_lines.append(" ".join(parts) + "\n")
                    continue

            new_lines.append(line)

        with open(file_path, "w") as f:
            f.writelines(new_lines)

    @staticmethod
    def remap_folder_name_from_mapping(folder_name, mapping):
        """
        Example:
            folder_name = "sqs_lev=1_a=0.75,b=0.25,c=0"
            mapping = {"a": "b", "b": "a", "c": "c"}

        returns:
            "sqs_lev=1_a=0.25,b=0.75,c=0"
        """
        marker = "sqs_lev="
        if not folder_name.startswith(marker):
            return folder_name

        parts = folder_name.split("_", 2)
        if len(parts) < 3:
            return folder_name

        prefix = f"{parts[0]}_{parts[1]}"
        comp_part = parts[2]

        assignments = comp_part.split(",")
        old_values = {}
        label_order = []

        for item in assignments:
            label, value = item.split("=", 1)
            old_values[label] = value
            label_order.append(label)

        inverse_mapping = {v: k for k, v in mapping.items()}

        new_assignments = []
        for new_label in label_order:
            source_label = inverse_mapping[new_label]
            new_value = old_values[source_label]
            new_assignments.append(f"{new_label}={new_value}")

        return f"{prefix}_{','.join(new_assignments)}"

    @staticmethod
    def duplicate_structure_directories_with_permutations(workdir, n_elements):
        """
        Duplicate each sqs_lev=* directory for all lowercase-label permutations.

        Folder names are rewritten to reflect the actual remapped composition.
        Only lowercase labels are remapped in file contents.
        Uppercase fixed species like B remain unchanged.
        """
        workdir = Path(workdir)
        labels = BladeSQS.lowercase_labels(n_elements)
        all_perms = list(permutations(labels))

        for phase_dir in workdir.iterdir():
            if not phase_dir.is_dir():
                continue

            src_dirs = [p for p in phase_dir.iterdir() if p.is_dir() and p.name.startswith("sqs_lev=")]

            for src_dir in src_dirs:
                for perm in all_perms:
                    mapping = dict(zip(labels, perm))

                    # Skip identity permutation
                    if all(mapping[k] == k for k in labels):
                        continue

                    dst_name = BladeSQS.remap_folder_name_from_mapping(src_dir.name, mapping)
                    dst_dir = src_dir.parent / dst_name

                    if dst_dir == src_dir:
                        continue

                    # Avoid overwriting if another permutation already created this exact composition folder
                    if dst_dir.exists():
                        continue

                    shutil.copytree(src_dir, dst_dir)

                    for candidate in ["POSCAR", "str_template.out"]:
                        candidate_path = dst_dir / candidate
                        if candidate_path.exists():
                            BladeSQS.apply_lowercase_permutation_to_file(candidate_path, mapping)

                    print(f"Created remapped directory: {dst_dir}")

    def generate_custom_sqs(
        self,
        composition,
        supercell_size,
        shell_weights,
        iterations,
    ):
        comp = Composition(composition)

        # Build the physical lattice from the phase parameters
        base_lattice = Lattice.from_parameters(
            a=self.a,
            b=self.b,
            c=self.c,
            alpha=self.alpha,
            beta=self.beta,
            gamma=self.gamma,
        )

        # Apply the prototype vector matrix on top of the physical lattice
        vec_matrix = np.array([
            [float(x) for x in line.split()]
            for line in self.vectors.strip().splitlines()
        ], dtype=float)

        lattice = vec_matrix @ np.array(base_lattice.matrix, dtype=float)

        variable_coords = []
        fixed_coords = []
        fixed_species = []

        for line in self.unit_cell.strip().splitlines():
            parts = line.split()
            if len(parts) < 4:
                raise ValueError("Each coordinate line must include x y z and a site label.")

            xyz = [float(parts[0]), float(parts[1]), float(parts[2])]
            site = parts[3]

            if site.islower():
                variable_coords.append(xyz)
            else:
                fixed_coords.append(xyz)
                fixed_species.append(site)

        if not variable_coords:
            raise ValueError("No lowercase variable sites found.")

        n_variable_per_cell = len(variable_coords)
        n_variable_sites = int(n_variable_per_cell * np.prod(supercell_size))

        frac_dict = comp.fractional_composition.as_dict()
        composition_dict = {
            el: int(round(frac * n_variable_sites))
            for el, frac in frac_dict.items()
        }

        diff = n_variable_sites - sum(composition_dict.values())
        if diff != 0:
            first_key = next(iter(composition_dict))
            composition_dict[first_key] += diff

        configuration = {
            "structure": {
                "lattice": lattice,
                "coords": variable_coords,
                "species": ["Xe"] * len(variable_coords),
                "supercell": supercell_size,
            },
            "iterations": iterations,
            "shell_weights": shell_weights,
            "composition": composition_dict,
            "iteration_mode": "random",
        }

        parsed = parse_config(configuration)
        results = optimize(parsed, level=LogLevel.warn)

        variable_structure = to_pymatgen(results.best().structure())
        objective = results.best().objective

        # Add back fixed sites across the full supercell
        full_lattice = variable_structure.lattice
        full_species = [site.species_string for site in variable_structure]
        full_coords = [site.frac_coords for site in variable_structure]

        nx, ny, nz = supercell_size
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    shift = np.array([i, j, k], dtype=float)
                    for sp, fc in zip(fixed_species, fixed_coords):
                        new_fc = (np.array(fc) + shift) / np.array([nx, ny, nz], dtype=float)
                        full_species.append(sp)
                        full_coords.append(new_fc)

        structure = Structure(
            lattice=full_lattice,
            species=full_species,
            coords=full_coords,
            coords_are_cartesian=False,
        ).get_sorted_structure()

        return structure, objective

    def generate_all_poscars(
        self,
        cases,
        elements,
        phases,
        workdir,
        default_supercell_size=(2, 2, 2),
        shell_weights=None,
        iterations=2000,
    ):
        if shell_weights is None:
            shell_weights = {
                1: 1.0,
                2: 0.5,
            }

        for case in cases:
            level = case["level"]
            fractions = case["fractions"]
            composition = BladeSQS.composition_string(elements, fractions)

            crystal_structure = phases["generator_name"]
            phase_supercell_size = BladeSQS.get_phase_supercell_size(phases, default_supercell_size)

            print(
                f"\nGenerating level {level}: "
                f"{composition} {crystal_structure.upper()}"
            )

            try:
                structure, objective = self.generate_custom_sqs(
                    composition=composition,
                    supercell_size=phase_supercell_size,
                    shell_weights=shell_weights,
                    iterations=iterations,
                )

                filename = workdir / BladeSQS.poscar_filename(
                    elements,
                    fractions,
                    crystal_structure,
                )

                filename.parent.mkdir(parents=True, exist_ok=True)

                print(f"Composition: {composition}")
                print(f"Supercell size: {phase_supercell_size}")
                print(f"Objective: {objective}")
                print(f"Number of atoms: {len(structure)}")

                Poscar(structure).write_file(filename)
                print(f"Saved to {filename}")

            except Exception as exc:
                print(
                    f"Failed for {composition} "
                    f"{crystal_structure.upper()} with supercell_size={phase_supercell_size}: {exc}"
                )

    @staticmethod
    def prepare_structure_directories(cases, elements, phases, workdir):
        for case in cases:
            level = case["level"]
            fractions = case["fractions"]
            comp_name = BladeSQS.folder_name(elements, fractions, level)

            src = workdir / BladeSQS.poscar_filename(
                elements,
                fractions,
                phases["generator_name"],
            )
            dst_dir = workdir / phases["lattice"] / comp_name
            dst = dst_dir / "POSCAR"
            template_str = dst_dir / "str_template.out"

            if not src.exists():
                print(f"Skipping missing file: {src}")
                continue

            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(src, dst)
            BladeSQS.build_initial_template_str_from_poscar(dst, template_str)
            BladeSQS.relabel_str_template(template_str, elements)

            print(f"Level {level}: copied {src} -> {dst}")
            print(f"Level {level}: wrote template {template_str}")

    @staticmethod
    def rewrite_all_poscars_to_letters(workdir, elements):
        workdir = Path(workdir)
        for poscar_path in workdir.rglob("POSCAR*"):
            if poscar_path.is_file():
                BladeSQS.relabel_poscar_species(poscar_path, elements)

    def sqs_gen(
        self,
        len_comp,
        specific_phase,
        work_path,
        iterations=2000,
        shell_weights=None,
        default_supercell_size=(2, 2, 2),
    ):
        workdir = Path(work_path) / "SQS" / f"{specific_phase['lattice']}_{len_comp}"
        workdir.mkdir(parents=True, exist_ok=True)

        elements = self.placeholder_elements(len_comp)
        cases = self.expand_all_cases(self.sqsgen_levels[: (self.level + 1)], elements)

        print("\nExpanded compositions:")
        for case in cases:
            print(
                f"level={case['level']} "
                f"fractions={case['fractions']} "
                f"name={self.composition_string(elements, case['fractions'])}"
            )

        self.generate_all_poscars(
            cases=cases,
            elements=elements,
            phases=specific_phase,
            workdir=workdir,
            default_supercell_size=default_supercell_size,
            shell_weights=shell_weights,
            iterations=iterations,
        )

        self.prepare_structure_directories(
            cases=cases,
            elements=elements,
            phases=specific_phase,
            workdir=workdir,
        )

        self.rewrite_all_poscars_to_letters(workdir, elements)

        self.duplicate_structure_directories_with_permutations(
            workdir=workdir,
            n_elements=len_comp,
        )