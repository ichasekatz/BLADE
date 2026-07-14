"""Shared utilities: animation, image tiling, CSV validation."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def system_tables_dir(tables_root: str | Path, metals, phase_element=None) -> Path:
    """Return the table folder dedicated to one modeled composition system."""
    root = Path(tables_root)
    name = system_key(metals, phase_element)
    return root if root.name == name else root / name


def prepare_system_tables_dir(tables_root: str | Path, metals, phase_element=None) -> Path:
    """Create a system table folder and migrate matching legacy flat files."""
    root = Path(tables_root)
    destination = system_tables_dir(root, metals, phase_element)
    destination.mkdir(parents=True, exist_ok=True)
    if destination == root:
        return destination

    tag = system_tag(metals, phase_element)
    for source in root.glob(f"{tag}_*"):
        if not source.is_file():
            continue
        target = destination / source.name
        if target.exists():
            continue
        source.replace(target)
        print(f"  Moved legacy table: {source.name} -> {destination.name}/")
    return destination


def make_animation(
    frame_paths, out_gif: Path, out_mp4: Path, fps: int = 2, mp4_crf: int = 18, mp4_preset: str = "slow"
) -> None:
    frame_paths = [Path(p) for p in frame_paths if Path(p).exists()]
    if len(frame_paths) < 2:
        return
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    try:
        import imageio.v2 as imageio
        from PIL import Image
    except ImportError:
        print("  animation skipped: imageio/Pillow missing")
        return

    raw = [Image.open(p).convert("RGB") for p in frame_paths]
    w0, h0 = raw[0].size
    frames = [im.resize((w0, h0), Image.LANCZOS) if im.size != (w0, h0) else im for im in raw]
    imageio.mimsave(str(out_gif), [np.array(im) for im in frames], fps=fps, loop=0)

    try:
        import shutil
        import subprocess
        import tempfile

        if not shutil.which("ffmpeg"):
            return
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            for i, im in enumerate(frames):
                im.save(tmp / f"frame_{i:04d}.png")
            subprocess.run(
                [
                    "ffmpeg",
                    "-y",
                    "-framerate",
                    str(fps),
                    "-i",
                    str(tmp / "frame_%04d.png"),
                    "-vf",
                    "scale=trunc(iw/2)*2:trunc(ih/2)*2",
                    "-c:v",
                    "libx264",
                    "-preset",
                    mp4_preset,
                    "-crf",
                    str(mp4_crf),
                    "-pix_fmt",
                    "yuv420p",
                    str(out_mp4),
                ],
                capture_output=True,
                check=False,
            )
    except Exception as e:
        print(f"  mp4 skipped for {out_mp4.name}: {e}")


def tile_images(image_paths, out_path: Path, cols: int = 4) -> Path | None:
    from PIL import Image, ImageOps

    image_paths = [Path(p) for p in image_paths if Path(p).exists()]
    if not image_paths:
        return None
    import math

    imgs = [Image.open(p).convert("RGB") for p in image_paths]
    thumb_w = min(720, max(360, imgs[0].size[0] // 2))
    thumb_h = int(thumb_w * imgs[0].size[1] / imgs[0].size[0])
    thumbs = [ImageOps.contain(im, (thumb_w, thumb_h), Image.LANCZOS) for im in imgs]
    rows = int(math.ceil(len(thumbs) / cols))
    canvas = Image.new("RGB", (cols * thumb_w, rows * thumb_h), "white")
    for idx, im in enumerate(thumbs):
        canvas.paste(im, ((idx % cols) * thumb_w, (idx // cols) * thumb_h))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(out_path)
    return out_path


def csv_has_rows(path: Path, expected_rows: int) -> bool:
    if not path.exists():
        return False
    try:
        import pandas as pd

        return len(pd.read_csv(path, usecols=[0])) >= expected_rows
    except Exception:
        return False


def fmt_frac(value: float) -> str:
    return f"{value:.2f}".rstrip("0").rstrip(".").replace(".", "p")


def _fmt_val(v: float) -> str:
    """Format a fraction value: strip trailing zeros; edge values → '0' or '1'."""
    if v <= 0.0:
        return "0"
    if v >= 1.0:
        return "1"
    return f"{v:.2f}".rstrip("0").rstrip(".")


def fmt_range(value_min: float, value_max: float, step: float) -> str:
    import math

    lo = max(0.0, math.floor((value_min + 1e-12) / step) * step)
    hi = min(1.0, math.ceil((value_max - 1e-12) / step) * step)
    # Snap near-zero lo to 0, near-one hi to 1 (0.01 edge tolerance)
    if lo <= 0.01:
        lo = 0.0
    if hi >= 0.99:
        hi = 1.0
    if abs(hi - lo) < step * 0.5:
        return _fmt_val(lo)
    return f"{_fmt_val(lo)}-{_fmt_val(hi)}"


def normalize_region_labels(labels: list[str]) -> list[str]:
    """Clean up labels read from cached CSVs.

    Collapses degenerate ranges where both values are identical:
      '1.00-1.00Hf' → '1.00Hf'
      '0.05-0.05Cr' → '0.05Cr'
    Also collapses near-pure mixed-phase compositions while preserving any
    configured formula suffix.
    """
    import re

    _DUPE_RANGE = re.compile(r"(\d+\.\d+)-\1(?=[A-Z])")
    _PHASE_SEG = re.compile(r"\(([^)]+)\)([A-Z][A-Za-z0-9.]*)?")
    _TOK = re.compile(r"(?:(\d+\.\d+(?:-\d+\.\d+)?)\s*)?([A-Z][a-z]?)")

    def _collapse_phase(segment: str, suffix: str = "") -> str:
        toks = []
        for raw in segment.split():
            m = _TOK.fullmatch(raw)
            if not m:
                continue
            frac_s, metal = m.groups()
            if frac_s is None:
                return f"({segment}){suffix}"
            if "-" in frac_s:
                lo_s, hi_s = frac_s.split("-", 1)
                lo_v = float(lo_s) if lo_s not in ("0", "") else 0.0
                hi_v = float(hi_s) if hi_s != "1" else 1.0
            else:
                lo_v = hi_v = float(frac_s) if frac_s not in ("0", "1") else (0.0 if frac_s == "0" else 1.0)
            toks.append((metal, lo_v, hi_v))
        if not toks:
            return f"({segment}){suffix}"
        # Pure end-member: any metal's hi ≥ 0.99
        for metal, lo_v, hi_v in toks:
            if hi_v >= 0.99:
                return f"{metal}{suffix}"
        # Wide range spanning nearly 0→1: collapse to pure end-member of dominant metal
        for metal, lo_v, hi_v in toks:
            if lo_v <= 0.01 and hi_v >= 0.99:
                return f"{metal}{suffix}"
        return f"({segment}){suffix}"

    out = []
    for lbl in labels:
        lbl = _DUPE_RANGE.sub(r"\1", lbl)
        lbl = _PHASE_SEG.sub(lambda m: _collapse_phase(m.group(1), m.group(2) or ""), lbl)
        out.append(lbl)
    return out


def stoichiometric_suffix(element: str | None, stoichiometry: float) -> str:
    """Return a formula suffix such as ``X2`` from explicit phase settings."""
    if not element or stoichiometry <= 0:
        return ""
    value = float(stoichiometry)
    coefficient = "" if abs(value - 1.0) < 1e-12 else f"{value:g}"
    return f"{element}{coefficient}"


def system_key(metals, phase_element: str | None = None) -> str:
    """Stable system folder name; stoichiometry is intentionally omitted."""
    return "".join(m.title() for m in metals) + (phase_element or "")


def phase_formula(metals, phase_element: str | None, stoichiometry: float) -> str:
    return "".join(m.title() for m in metals) + stoichiometric_suffix(phase_element, stoichiometry)


def system_tag(metals, phase_element: str | None = None) -> str:
    return "".join(m.lower() for m in metals) + (phase_element or "").lower()


def phase_short(pid: str) -> str:
    return str(pid).split("_")[-1]
