import matplotlib.pyplot as plt
from ase.io import read
from ase.visualize.plot import plot_atoms
from PIL import Image
import numpy as np
from ase import Atoms


def read_poscar_inline_symbols(poscar_path):
    with open(poscar_path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    scale = float(lines[1])

    # cell vectors (3 lines)
    cell = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)], dtype=float)
    cell *= scale

    # After cell: either a counts line, or a symbols line then counts line.
    # Your example: counts directly: "48 7 49"
    i = 5
    toks = lines[i].split()

    def all_int(ts):
        try:
            [int(x) for x in ts]
            return True
        except Exception:
            return False

    if all_int(toks):
        counts = [int(x) for x in toks]
        i += 1
    else:
        # symbols line then counts line
        i += 1
        counts = [int(x) for x in lines[i].split()]
        i += 1

    # Optional "Selective dynamics"
    if lines[i].lower().startswith("s"):
        i += 1

    # Coord mode
    mode = lines[i].lower()
    direct = mode.startswith("d")
    i += 1

    n_atoms = sum(counts)

    symbols = []
    positions = []
    for ln in lines[i : i + n_atoms]:
        parts = ln.split()
        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
        sym = parts[-1]  # last token is element, e.g. B / Cr / Hf
        positions.append([x, y, z])
        symbols.append(sym)

    positions = np.array(positions, dtype=float)

    # Convert Direct -> Cartesian if needed
    if direct:
        positions = positions @ cell

    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)

class BLADEVisualizer:
    """
    Minimal POSCAR visualizer.
    """

    def poscar(self, poscars, save=None):
        """
        Visualize a POSCAR/CONTCAR.

        Parameters
        ----------
        poscar : str or Path
            Path to POSCAR or CONTCAR
        save : str or Path, optional
            If given, save image to this path
        """
        fig, axes = plt.subplots(1, len(poscars), figsize=(4 * len(poscars), 4))

        if len(poscars) == 1:
            axes = [axes]

        for ax, p in zip(axes, poscars):
            atoms = read_poscar_inline_symbols(p)
            plot_atoms(atoms, ax, rotation="65x,45y,0z", show_unit_cell=True)
            ax.set_axis_off()

        if save is not None:
            fig.savefig(save, dpi=200, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()

    def phase_diagram(self, images, save):
            """
            Combine PNG images side-by-side and save.

            Parameters
            ----------
            images : list[str | Path]
                Paths to image files
            save : str | Path
                Output image path
            """
            imgs = [Image.open(p) for p in images]

            widths, heights = zip(*(img.size for img in imgs))
            total_width = sum(widths)
            max_height = max(heights)

            combined = Image.new("RGB", (total_width, max_height), "white")

            x_offset = 0
            for img in imgs:
                combined.paste(img, (x_offset, 0))
                x_offset += img.width

            combined.save(save)

