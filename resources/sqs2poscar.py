import numpy as np
from collections import OrderedDict
from pathlib import Path

# ============================================================
# USER SETTINGS  (EDIT ONLY THESE TWO LINES)
# ============================================================

# sqsdb_lev=0_a=1
# sqsdb_lev=1_a=0.5,0.5
# sqsdb_lev=2_a=0.75,0.25
# sqsdb_lev=3_a=0.33333,0.33333,0.33333
# sqsdb_lev=4_a=0.5,0.25,0.25
# sqsdb_lev=5_a=0.625,0.375
# sqsdb_lev=5_a=0.875,0.125

INPUT_SQS_PATH = r"/Users/chasekatz/Desktop/School/Research/PhaseForge/PhaseForge/atat/data/sqsdb/FCC_A1/sqsdb_lev=3_a=0.33333,0.33333,0.33333/bestsqs.out"
OUTPUT_POSCAR_PATH = r"/Users/chasekatz/Desktop/School/Research/BLADE/BLADE/resources/POSCARs/POSCARFCC"
POSCAR_COMMENT = "SQS structure (ATAT -> POSCAR)"
SCALE_FACTOR = 1.0  # POSCAR universal scale

# ============================================================
# Helpers
# ============================================================

def inverse_matrix_3x3(m: np.ndarray) -> np.ndarray:
    """Exact analog of the C++ inverse_matrix() with determinant check."""
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

# ============================================================
# Read ATAT bestsqs.out using the same logic as the C++ code
# Format expected (ATAT 3.0 mcsqs bestsqs.out):
#   line 1-3: basis coordinate system (Vec1)
#   line 4-6: lattice vectors in basis coordinates (Vec2)
#   line 7- : atomic positions in basis coordinates + element label
# ============================================================

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

    # Vec1: basis vectors (3 lines)
    vec1 = np.array([parse3(0), parse3(1), parse3(2)], dtype=float)  # 3x3

    # Vec2: lattice vectors in basis coordinates (3 lines)
    vec2 = np.array([parse3(3), parse3(4), parse3(5)], dtype=float)  # 3x3

    coords_basis = []
    species = []
    for i in range(6, len(lines)):
        toks = lines[i].split()
        if len(toks) < 4:
            raise ValueError(f"Bad атом line (need x y z elem) on line {i+1}: {lines[i]}")
        x, y, z = float(toks[0]), float(toks[1]), float(toks[2])
        elem = toks[3]
        coords_basis.append([x, y, z])
        species.append(elem)

    if not coords_basis:
        raise ValueError("No atomic positions found in bestsqs.out")

    return vec1, vec2, np.array(coords_basis, dtype=float), species

# ============================================================
# Convert using the same math as the C++ code
# - Lattice in Cartesian: LatVec = Vec2 * Vec1  (with the same index convention)
# - Atom Cartesian: AtomCart = AtomBasis * Vec1
# - Sort atoms by first-appearance element order
# - Convert to Direct using inverse(LatVec)
# ============================================================

def convert_bestsqs_to_poscar(vec1, vec2, atom_basis, species):
    # C++ computes:
    # LatVec[i][j] = vec2[i][0]*vec1[0][j] + vec2[i][1]*vec1[1][j] + vec2[i][2]*vec1[2][j]
    # That is plain matrix multiplication: LatVec = vec2 @ vec1
    latvec = vec2 @ vec1  # 3x3

    # Atom Cartesian:
    # atom_cart = atom_basis @ vec1
    atom_cart = atom_basis @ vec1

    # Element ordering by first appearance (same as C++)
    elem_order = unique_preserve_order(species)

    # Counts
    counts = [species.count(el) for el in elem_order]

    # Sort atoms by elem_order (stable within each element like C++)
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

    # Direct coords:
    # C++ uses custom inverse and then:
    # direct_x = x*Inv[0][0] + y*Inv[1][0] + z*Inv[2][0]
    # This corresponds to: direct = cart @ inv(latvec)
    latvec_inv = inverse_matrix_3x3(latvec)
    direct = sorted_cart @ latvec_inv

    # Wrap to [0,1) for VASP friendliness (C++ did not, but it's safe & avoids negatives)
    direct = direct - np.floor(direct)

    return latvec, elem_order, counts, direct

# ============================================================
# Write POSCAR
# ============================================================

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

# ============================================================
# Main
# ============================================================

print("Reading ATAT bestsqs.out (C++-compatible parser)...")
vec1, vec2, atom_basis, species = read_bestsqs_cpp_compatible(INPUT_SQS_PATH)
print(f"Atoms read: {len(species)}")

latvec, elem_order, counts, direct = convert_bestsqs_to_poscar(vec1, vec2, atom_basis, species)

print("Writing POSCAR...")
write_poscar(
    OUTPUT_POSCAR_PATH,
    POSCAR_COMMENT,
    SCALE_FACTOR,
    latvec,
    elem_order,
    counts,
    direct
)

print("\nSUCCESS")
print("Input :", INPUT_SQS_PATH)
print("Output:", OUTPUT_POSCAR_PATH)
print("Elements:", " ".join(elem_order))
print("Counts  :", " ".join(map(str, counts)))


# INPUT_SQS_PATH = r"/Users/chasekatz/Desktop/School/Research/PhaseForge/PhaseForge/atat/data/sqsdb/FCC1_2/sqsdb_lev=3_a=0.33333,0.33333,0.33333/bestsqs.out"
# OUTPUT_POSCAR_PATH = r"/Users/chasekatz/Desktop/School/Research/BLADE/BLADE/resources/POSCARs/POSCAR2"

# POSCAR_COMMENT = "SQS structure (ATAT -> POSCAR)"

# # ============================================================
# # Read ATAT SQS file
# # ============================================================

# def read_sqs(sqs_path):
#     sqs_path = Path(sqs_path)

#     if not sqs_path.exists():
#         raise FileNotFoundError(f"SQS file not found:\n{sqs_path}")

#     with open(sqs_path) as f:
#         lines = [l.strip() for l in f if l.strip()]

#     if len(lines) < 4:
#         raise ValueError("File does not look like a valid ATAT SQS structure")

#     # First 3 lines = lattice vectors (Cartesian Å)
#     lattice = np.array([
#         list(map(float, lines[0].split())),
#         list(map(float, lines[1].split())),
#         list(map(float, lines[2].split()))
#     ], dtype=float)

#     coords = []
#     species = []

#     for line in lines[3:]:
#         toks = line.split()
#         if len(toks) < 4:
#             continue
#         coords.append([float(toks[0]), float(toks[1]), float(toks[2])])
#         species.append(toks[3])

#     if len(coords) == 0:
#         raise ValueError("No atomic coordinates found in SQS file")

#     return lattice, np.array(coords), species


# # ============================================================
# # Write POSCAR (no species reordering)
# # ============================================================

# def write_poscar(poscar_path, lattice, coords, species, comment):

#     poscar_path = Path(poscar_path)
#     poscar_path.parent.mkdir(parents=True, exist_ok=True)

#     # Preserve first-appearance order of species
#     species_map = OrderedDict()
#     for i, sp in enumerate(species):
#         if sp not in species_map:
#             species_map[sp] = []
#         species_map[sp].append(coords[i])

#     unique_species = list(species_map.keys())
#     counts = [len(species_map[sp]) for sp in unique_species]

#     with open(poscar_path, "w") as f:

#         # header
#         f.write(comment + "\n")
#         f.write("1.0\n")

#         # lattice
#         for vec in lattice:
#             f.write(f"{vec[0]:16.10f} {vec[1]:16.10f} {vec[2]:16.10f}\n")

#         # species and counts
#         f.write(" ".join(unique_species) + "\n")
#         f.write(" ".join(map(str, counts)) + "\n")

#         # fractional coords
#         f.write("Direct\n")

#         for sp in unique_species:
#             for c in species_map[sp]:
#                 f.write(f"{c[0]:16.10f} {c[1]:16.10f} {c[2]:16.10f}\n")


# # ============================================================
# # Main
# # ============================================================


# print("Reading SQS structure...")
# lattice, coords, species = read_sqs(INPUT_SQS_PATH)

# print(f"Atoms read: {len(species)}")
# print("Writing POSCAR...")

# write_poscar(
#     OUTPUT_POSCAR_PATH,
#     lattice,
#     coords,
#     species,
#     POSCAR_COMMENT
# )

# print("\nSUCCESS")
# print("Input :", INPUT_SQS_PATH)
# print("Output:", OUTPUT_POSCAR_PATH)

