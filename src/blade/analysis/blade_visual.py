"""
This module defines the `BladeVisualizer` class, which is used for visualizing atomic structures and phase diagrams.

The `BladeVisualizer` class provides helper functions for reading POSCAR files with inline element symbols
and for visualizing atomic structures and phase diagrams using ASE and Matplotlib. It is primarily intended
for lightweight inspection and comparison of SQS and related structures generated within the BLADE workflow.
"""

import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.visualize.plot import plot_atoms
from PIL import Image


class BLADEVisualizer:
    """
    Visualization utilities for BLADE structures and phase diagrams.

    This class provides convenience methods for:
      - rendering one or more POSCAR structures side-by-side
      - combining multiple phase diagram images into a single figure

    The visualizations are intended for quick inspection and reporting rather than publication-quality rendering.
    """
    def poscar(self, poscars, save=None):
        """
        Visualize one or more POSCAR structures.

        Each POSCAR is read using inline element symbols and rendered as a 2D projection using ASE's plotting
        utilities. Multiple structures are displayed side-by-side for easy comparison.

        Args:
            poscars (list[str | Path]): List of paths to POSCAR files to visualize.
            save (str | Path, optional): If provided, save the resulting figure to this path. If None, the figure
                is displayed interactively.
        """
        fig, axes = plt.subplots(1, len(poscars), figsize=(4 * len(poscars), 4))

        if len(poscars) == 1:
            axes = [axes]

        for ax, p in zip(axes, poscars):
            atoms = self.read_poscar_inline_symbols(p)
            plot_atoms(atoms, ax, rotation="65x,45y,0z", show_unit_cell=True)
            ax.set_axis_off()

        if save is not None:
            fig.savefig(save, dpi=200, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()

    def phase_diagram(self, images, save):
        """
        Combine multiple phase diagram images into a single image.

        This method concatenates a list of image files horizontally into one combined image and saves
        the result. All images are aligned at the top and padded with a white background as needed.

        Args:
            images (list[str | Path]): Paths to image files to combine.
            save (str | Path): Output path for the combined image.
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

    def read_poscar_inline_symbols(self, poscar_path):
        """
        Read a POSCAR file that includes inline element symbols per atom.

        This function parses a VASP POSCAR-style file where each atomic coordinate line ends with an
        explicit element symbol. It supports both Direct and Cartesian coordinate formats and returns an
        ASE Atoms object with periodic boundary conditions enabled.

        Args:
            poscar_path (str | Path): Path to the POSCAR file to be read.

        Returns:
            ase.Atoms: ASE Atoms object containing atomic symbols, Cartesian positions, lattice vectors,
                and periodic boundary conditions.
        """
        with open(poscar_path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        scale = float(lines[1])

        cell = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)], dtype=float)
        cell *= scale

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
            i += 1
            counts = [int(x) for x in lines[i].split()]
            i += 1

        if lines[i].lower().startswith("s"):
            i += 1

        mode = lines[i].lower()
        direct = mode.startswith("d")
        i += 1

        n_atoms = sum(counts)

        symbols = []
        positions = []
        for ln in lines[i : i + n_atoms]:
            parts = ln.split()
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            sym = parts[-1]
            positions.append([x, y, z])
            symbols.append(sym)

        positions = np.array(positions, dtype=float)

        if direct:
            positions = positions @ cell

        return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
