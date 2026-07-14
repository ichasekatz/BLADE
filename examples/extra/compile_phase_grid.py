"""Compile all phase diagram GIFs into synchronized grid videos (MP4).

Each GIF is resized to panel_w x panel_h and tiled into an n_cols x n_rows grid.
All GIFs should have the same number of frames.

Usage:
    pixi run python examples/extra/compile_phase_grid.py
"""

from __future__ import annotations

import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0    = Path("/Users/chasekatz/Desktop/School/Research")
path1    = path0 / "BLADE"
gif_dir  = path1 / "Files" / "Phase_Diagrams"

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------
panel_w  = 150
panel_h  = 130
n_cols   = None
fps      = 20

bg_color = (255, 255, 255)
header_h = 40
cbar_w   = 60

# Temperature range — must match plot_phase.py gif settings
t_start_gif = 200
t_end_gif   = 2000
gif_t_step  = 10

temperatures = np.arange(t_start_gif, t_end_gif + gif_t_step, gif_t_step)
print(f"Temperature frames expected: {len(temperatures)}")

cmap = mcm.plasma
norm = mcolors.Normalize(vmin=t_start_gif, vmax=t_end_gif)


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def load_font(size: int):
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def parse_elements(gif_path: Path) -> list[str]:
    """Extract element list from filename.

    Example:
        CrHfNb_phase_evolution_clean.gif -> ['Cr', 'Hf', 'Nb']
    """
    stem = (
        gif_path.stem
        .replace("_phase_evolution_clean", "")
        .replace("_phase_evolution", "")
    )
    return re.findall(r"[A-Z][a-z]?", stem)


def overlay_element_triangle(
    img: Image.Image,
    elements: list[str],
    panel_w: int,
    panel_h: int,
) -> Image.Image:
    """Draw element labels in the top-right of panel."""
    if len(elements) != 3:
        return img

    draw = ImageDraw.Draw(img)
    font = load_font(11)

    cx = panel_w - 28
    cy = 14

    positions = {
        elements[2]: (cx, cy),           # top
        elements[0]: (cx - 16, cy + 16), # bottom-left
        elements[1]: (cx + 16, cy + 16), # bottom-right
    }

    for el, (x, y) in positions.items():
        draw.text((x + 1, y + 1), el, fill=(0, 0, 0), font=font, anchor="mm")
        draw.text((x, y), el, fill=(255, 255, 255), font=font, anchor="mm")

    return img


def get_gif_frame(reader: Image.Image, frame_index: int) -> Image.Image:
    """Get GIF frame safely using PIL.

    If requested frame is beyond available frames, use the final frame.
    """
    max_frame = getattr(reader, "n_frames", 1) - 1
    reader.seek(min(frame_index, max_frame))
    return reader.convert("RGB")


# ------------------------------------------------------------------
# Find all GIFs and group by element count
# ------------------------------------------------------------------
all_gifs = sorted(gif_dir.rglob("*_phase_evolution_clean.gif"))

print(f"Found {len(all_gifs)} clean GIFs")

if not all_gifs:
    raise SystemExit("No clean GIFs found — run plot_phase.py first.")

groups: dict[int, list[Path]] = defaultdict(list)

for p in all_gifs:
    els = parse_elements(p)
    groups[len(els)].append(p)

size_names = {
    2: "binary",
    3: "ternary",
    4: "quaternary",
    5: "quinary",
    6: "senary",
    7: "septenary",
    8: "octonary",
    9: "nonary",
    10: "denary",
}

# ------------------------------------------------------------------
# Main writer
# ------------------------------------------------------------------
def write_grid(gif_paths: list[Path], out_mp4: Path, label: str) -> None:
    n = len(gif_paths)

    nc = n_cols if n_cols is not None else math.ceil(math.sqrt(n))
    nr = math.ceil(n / nc)

    gw = nc * panel_w + cbar_w
    gh = nr * panel_h + header_h
    cbar_h = nr * panel_h

    print(f"\n{label}: {n} GIFs  grid {nc}×{nr}  {gw}×{gh}px")

    # Open GIFs with PIL
    readers: list[Image.Image] = []
    n_frames_list: list[int] = []

    for p in gif_paths:
        im = Image.open(p)
        nf = getattr(im, "n_frames", 1)
        readers.append(im)
        n_frames_list.append(nf)

    unique_counts = sorted(set(n_frames_list))
    expected_frames = len(temperatures)

    print(f"  frame counts found: {unique_counts}")
    print(f"  expected frames from temperature range: {expected_frames}")

    if len(unique_counts) != 1:
        print("  WARNING: not all GIFs have the same frame count.")
        print("  Shorter GIFs will hold their last frame.")

    # Use expected temperature frames unless every GIF is shorter
    n_frames = min(expected_frames, max(n_frames_list))

    print(f"  using {n_frames} frames")

    # Colorbar strip
    cbar_arr = np.zeros((cbar_h, cbar_w, 3), dtype=np.uint8)

    for row_px in range(cbar_h):
        frac = 1.0 - row_px / cbar_h
        rgba = cmap(frac)
        cbar_arr[row_px, :] = [
            int(rgba[0] * 255),
            int(rgba[1] * 255),
            int(rgba[2] * 255),
        ]

    ffmpeg_cmd = [
        "ffmpeg",
        "-y",
        "-f", "rawvideo",
        "-vcodec", "rawvideo",
        "-s", f"{gw}x{gh}",
        "-pix_fmt", "rgb24",
        "-r", str(fps),
        "-i", "-",
        "-vcodec", "libx264",
        "-crf", "23",
        "-pix_fmt", "yuv420p",
        str(out_mp4),
    ]

    proc = subprocess.Popen(
        ffmpeg_cmd,
        stdin=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )

    if proc.stdin is None:
        raise RuntimeError("Could not open ffmpeg stdin.")

    try:
        for f in range(n_frames):
            T = temperatures[f] if f < len(temperatures) else t_end_gif

            if f % 50 == 0:
                print(f"  Frame {f}/{n_frames}  T={int(T)} K")

            canvas = np.full((gh, gw, 3), bg_color, dtype=np.uint8)

            # Header
            rgba = cmap(norm(T))
            hdr_color = [
                int(rgba[0] * 255),
                int(rgba[1] * 255),
                int(rgba[2] * 255),
            ]

            canvas[:header_h, :] = hdr_color

            hdr_img = Image.fromarray(canvas[:header_h].copy())
            draw = ImageDraw.Draw(hdr_img)

            draw.text(
                (gw // 2, header_h // 2),
                f"{label}  |  T = {int(T)} K",
                fill=(255, 255, 255),
                anchor="mm",
                font=load_font(24),
            )

            canvas[:header_h] = np.array(hdr_img)

            # Panels
            for idx, (reader, gif_path) in enumerate(zip(readers, gif_paths)):
                row = idx // nc
                col = idx % nc

                els = parse_elements(gif_path)

                try:
                    img = get_gif_frame(reader, f)
                    img = img.resize((panel_w, panel_h), Image.LANCZOS)
                    img = overlay_element_triangle(img, els, panel_w, panel_h)
                    panel = np.array(img)

                except Exception:
                    panel = np.full(
                        (panel_h, panel_w, 3),
                        bg_color,
                        dtype=np.uint8,
                    )

                y0 = header_h + row * panel_h
                x0 = col * panel_w

                canvas[y0:y0 + panel_h, x0:x0 + panel_w] = panel

            # Grid lines
            for row in range(nr):
                for col in range(nc):
                    y0 = header_h + row * panel_h
                    x0 = col * panel_w

                    canvas[y0, x0:x0 + panel_w] = [80, 80, 80]
                    canvas[y0 + panel_h - 1, x0:x0 + panel_w] = [80, 80, 80]
                    canvas[y0:y0 + panel_h, x0] = [80, 80, 80]
                    canvas[y0:y0 + panel_h, x0 + panel_w - 1] = [80, 80, 80]

            # Colorbar
            canvas[header_h:header_h + cbar_h, nc * panel_w:] = cbar_arr

            marker_y = header_h + int((1.0 - norm(T)) * cbar_h)
            marker_y = max(header_h, min(header_h + cbar_h - 1, marker_y))

            canvas[marker_y - 1:marker_y + 2, nc * panel_w:] = [255, 255, 255]
            canvas[marker_y, nc * panel_w:] = [0, 0, 0]

            proc.stdin.write(canvas.tobytes())

    finally:
        proc.stdin.close()
        proc.wait()

        for r in readers:
            r.close()

    print(f"  Saved: {out_mp4.name}  ({out_mp4.stat().st_size / 1e6:.1f} MB)")


# ------------------------------------------------------------------
# Write one grid video per system size
# ------------------------------------------------------------------
for n_el in sorted(groups.keys()):
    gifs = sorted(groups[n_el])
    size_label = size_names.get(n_el, f"{n_el}-component")
    mp4_out = gif_dir / f"grid_{size_label}.mp4"

    write_grid(gifs, mp4_out, size_label)

print("\nAll grid videos saved.")