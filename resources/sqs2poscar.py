import numpy as np
from pathlib import Path

INPUT_SQS_PATH = r"/Users/chasekatz/Desktop/School/Research/PhaseForge/PhaseForge/atat/data/sqsdb/FCC_A1/sqsdb_lev=3_a=0.33333,0.33333,0.33333/bestsqs.out"
OUTPUT_POSCAR_PATH = r"/Users/chasekatz/Desktop/School/Research/BLADE/BLADE/resources/POSCARs/POSCARFCC"
POSCAR_COMMENT = "SQS structure (ATAT -> POSCAR)"
SCALE_FACTOR = 1.0
REPEAT_XYZ = (1, 1, 1)

def inverse_matrix_3x3(m: np.ndarray) -> np.ndarray:
    det = np.linalg.det(m)
    if abs(det) < 1e-14:
        raise ValueError("Singular lattice matrix (det ~ 0). Cannot invert.")
    return np.linalg.inv(m)

def unique_preserve_order(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

def read_bestsqs_cpp_compatible(path: str):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SQS file not found:\n{path}")

    lines = [l.strip() for l in path.read_text().splitlines() if l.strip()]
    if len(lines) < 7:
        raise ValueError("File does not look like a valid ATAT mcsqs bestsqs.out (too few lines).")

    def parse3(i):
        toks = lines[i].split()
        if len(toks) < 3:
            raise ValueError(f"Cannot parse 3 floats on line {i+1}: {lines[i]}")
        return [float(toks[0]), float(toks[1]), float(toks[2])]

    vec1 = np.array([parse3(0), parse3(1), parse3(2)], dtype=float)
    vec2 = np.array([parse3(3), parse3(4), parse3(5)], dtype=float)

    coords_basis = []
    species = []

    for i in range(6, len(lines)):
        toks = lines[i].split()
        if len(toks) < 4:
            raise ValueError(f"Bad atom line on line {i+1}: {lines[i]}")
        x, y, z = float(toks[0]), float(toks[1]), float(toks[2])
        elem = toks[3]
        coords_basis.append([x, y, z])
        species.append(elem)

    if not coords_basis:
        raise ValueError("No atomic positions found in bestsqs.out")

    return vec1, vec2, np.array(coords_basis, dtype=float), species

def convert_bestsqs_to_poscar(vec1, vec2, atom_basis, species):
    latvec = vec2 @ vec1
    atom_cart = atom_basis @ vec1

    elem_order = unique_preserve_order(species)
    counts = [species.count(el) for el in elem_order]

    sorted_cart = []
    sorted_species = []

    for el in elem_order:
        for i, sp in enumerate(species):
            if sp == el:
                sorted_cart.append(atom_cart[i])
                sorted_species.append(sp)

    sorted_cart = np.array(sorted_cart, dtype=float)

    if len(sorted_cart) != len(species):
        raise ValueError("Sorting error: atom count mismatch.")

    latvec_inv = inverse_matrix_3x3(latvec)
    direct = sorted_cart @ latvec_inv
    direct = direct - np.floor(direct)

    return latvec, elem_order, counts, direct, sorted_species

def repeat_supercell(latvec, species, direct_coords, reps=(1, 1, 1)):
    rx, ry, rz = reps

    if rx < 1 or ry < 1 or rz < 1:
        raise ValueError("Repeat values must be positive integers.")

    new_latvec = latvec.copy()
    new_latvec[0] *= rx
    new_latvec[1] *= ry
    new_latvec[2] *= rz

    scale = np.array([rx, ry, rz], dtype=float)
    new_species = []
    new_direct = []

    for ix in range(rx):
        for iy in range(ry):
            for iz in range(rz):
                shift = np.array([ix, iy, iz], dtype=float)
                for sp, coord in zip(species, direct_coords):
                    new_coord = (coord + shift) / scale
                    new_direct.append(new_coord)
                    new_species.append(sp)

    return new_latvec, new_species, np.array(new_direct, dtype=float)

def write_poscar(path: str, comment: str, scale: float, latvec: np.ndarray,
                 elem_order, counts, direct_coords: np.ndarray):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        f.write(f"{comment}\n")
        f.write(f"{scale:.16g}\n")

        for i in range(3):
            f.write(f"{latvec[i,0]:16.10f} {latvec[i,1]:16.10f} {latvec[i,2]:16.10f}\n")

        f.write(" ".join(elem_order) + "\n")
        f.write(" ".join(str(c) for c in counts) + "\n")
        f.write("Direct\n")

        for c in direct_coords:
            f.write(f"{c[0]:16.10f} {c[1]:16.10f} {c[2]:16.10f}\n")

print("Reading ATAT bestsqs.out (C++-compatible parser)...")
vec1, vec2, atom_basis, species = read_bestsqs_cpp_compatible(INPUT_SQS_PATH)
print(f"Atoms read: {len(species)}")

latvec, elem_order, counts, direct, sorted_species = convert_bestsqs_to_poscar(
    vec1,
    vec2,
    atom_basis,
    species
)

expanded_latvec, expanded_species, expanded_direct = repeat_supercell(
    latvec,
    sorted_species,
    direct,
    REPEAT_XYZ
)

expanded_elem_order = unique_preserve_order(expanded_species)
expanded_counts = [expanded_species.count(el) for el in expanded_elem_order]

print("Writing POSCAR...")
write_poscar(
    OUTPUT_POSCAR_PATH,
    POSCAR_COMMENT,
    SCALE_FACTOR,
    expanded_latvec,
    expanded_elem_order,
    expanded_counts,
    expanded_direct
)

print("\nSUCCESS")
print("Input :", INPUT_SQS_PATH)
print("Output:", OUTPUT_POSCAR_PATH)
print("Elements:", " ".join(expanded_elem_order))
print("Counts  :", " ".join(map(str, expanded_counts)))
print("Repeat  :", REPEAT_XYZ)