"""Quantitative image analysis of Ness et al. (2026) Figure 5 (six 50 mol% Group IV-Group V diboride melt pools).

Source image: ness_images/main/fig-004.png, extracted losslessly with `pdfimages -png -f 3 -l 8` from the paper PDF
(embedded JPEG, 300 ppi in the PDF; 2075 x 932 px). The original SEM files are not available; this is a published,
JPEG-compressed, rescaled figure, so all results carry that limitation.

Imaging (from the paper): secondary-electron (SE) micrographs of polished, 1 nm Au-coated melt pools. SE contrast is
mixed topographic + compositional; the authors identify light lamellae as (Group V)B2-rich (HfB2-NbB2: HfB2-rich)
from brightness and low-resolution EDS. EDS "lacked the resolution to map the lamellae".

Method per panel:
  1. Scale: the black scale bar inside the white label box (bottom right) is detected as the longest dark horizontal
     run; um/px = bar length (um, from the label) / run length (px).
  2. Regions of interest (ROI): rectangles chosen by eye to contain only the lamellar (eutectic) microstructure,
     excluding large non-lamellar grains, cracks, voids and labels. Coordinates below (figure pixels) and overlays
     in figures/ make the choice auditable.
  3. Phase segmentation: Gaussian smoothing (sigma 1 px), Otsu threshold; the light-phase area fraction A_L is
     reported for the Otsu threshold and for thresholds shifted by +-10 % of the grey-level range (sensitivity).
  4. Characteristic wavelength: radially averaged 2-D power spectrum of the mean-subtracted, Hann-windowed ROI;
     peak in the RDF -> lambda = 1/k (um). Compared with Ness's published values as a method check.
Outputs: results/p01_ness_lamellae.csv, figures/fig_p01_roi_overlays.png
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parents[1]
IMG = HERE / "ness_images" / "main" / "fig-004.png"
NESS_LAMBDA = {"TiB2-NbB2": 0.242, "TiB2-TaB2": 0.346, "ZrB2-NbB2": 0.382, "ZrB2-TaB2": 0.414, "HfB2-NbB2": 0.404,
               "HfB2-TaB2": 0.631}  # um, published Acta Mater. Table 2 (the SSRN preprint lists TiB2-TaB2 as 0.520)
# panel boxes (x0, y0, x1, y1) in figure pixels and the scale-bar length printed on each panel
PANELS = {
    "TiB2-NbB2": ((0, 0, 657, 442), 1.0), "ZrB2-NbB2": ((671, 0, 1329, 442), 5.0), "HfB2-NbB2": ((1342, 0, 2075, 442), 5.0),
    "TiB2-TaB2": ((0, 450, 657, 932), 5.0), "ZrB2-TaB2": ((671, 450, 1329, 932), 4.0), "HfB2-TaB2": ((1342, 450, 2075, 932), 5.0),
}
# lamellar ROIs (figure pixels), chosen by eye; see overlay figure
ROIS = {  # v2: lamella-only areas read from the v1 overlay (large grains, blobs, cracks and labels excluded)
    "TiB2-NbB2": [(170, 70, 430, 330)],
    "ZrB2-NbB2": [(936, 102, 1031, 176), (725, 255, 785, 340)],
    "HfB2-NbB2": [(1630, 108, 1700, 210), (1445, 245, 1555, 325)],
    "TiB2-TaB2": [(165, 610, 290, 690)],
    "ZrB2-TaB2": [(760, 660, 868, 741)],
    "HfB2-TaB2": [(1560, 480, 1740, 640)],
}


def scale_um_per_px(gray: np.ndarray, box, bar_um: float) -> tuple[float, int]:
    x0, y0, x1, y1 = box
    sub = gray[y0 + int(0.75 * (y1 - y0)):y1, x0 + int(0.55 * (x1 - x0)):x1]
    best = 0
    for row in sub:
        dark = row < 60
        run = cur = 0
        for d in dark:
            cur = cur + 1 if d else 0
            run = max(run, cur)
        best = max(best, run)
    return bar_um / best, best


def otsu(v: np.ndarray) -> float:
    hist, edges = np.histogram(v, bins=256, range=(0, 256))
    p = hist / hist.sum()
    w = np.cumsum(p)
    mu = np.cumsum(p * np.arange(256))
    sb = (mu[-1] * w - mu) ** 2 / (w * (1 - w) + 1e-12)
    return float(np.argmax(sb))


def wavelength(roi: np.ndarray, um_px: float) -> float:
    a = roi - roi.mean()
    a = a * np.outer(np.hanning(a.shape[0]), np.hanning(a.shape[1]))
    ps = np.abs(np.fft.fftshift(np.fft.fft2(a))) ** 2
    ny, nx = ps.shape
    yy, xx = np.indices(ps.shape)
    ky = (yy - ny // 2) / ny
    kx = (xx - nx // 2) / nx
    k = np.hypot(kx, ky)  # cycles per px
    bins = np.linspace(0, 0.5, 80)
    idx = np.digitize(k.ravel(), bins)
    rdf = np.array([ps.ravel()[idx == i].mean() if np.any(idx == i) else 0 for i in range(1, len(bins))])
    kc = 0.5 * (bins[1:] + bins[:-1])
    valid = (kc > 3.0 / min(ny, nx)) & (kc < 0.25) & (rdf > 0)  # beyond ROI-size and JPEG/pixel-noise frequencies
    lk, lr = np.log(kc[valid]), np.log(rdf[valid])
    bg = np.polyval(np.polyfit(lk, lr, 1), lk)  # aperiodic power-law background (as in Ness's peak/background separation)
    kpk = kc[valid][np.argmax(lr - bg)]
    return um_px / kpk


def main() -> None:
    from scipy.ndimage import gaussian_filter

    gray = np.asarray(Image.open(IMG).convert("L"), dtype=float)
    rows = []
    for sysn, (box, bar_um) in PANELS.items():
        um_px, bar_px = scale_um_per_px(gray, box, bar_um)
        for i, (x0, y0, x1, y1) in enumerate(ROIS[sysn]):
            roi = gaussian_filter(gray[y0:y1, x0:x1], 1.0)
            t = otsu(roi)
            sep = roi[roi > t].mean() - roi[roi <= t].mean()  # separation of the two class means
            fr = {f"A_light_{tag}": float((roi > t + d * sep).mean()) for tag, d in (("otsu", 0.0), ("minus10", -0.10), ("plus10", 0.10))}
            fr["class_separation_grey"] = float(sep)
            # bimodality: fraction of pixels within +-10 % of the threshold (high = weak contrast, unreliable fraction)
            fr["frac_pixels_near_threshold"] = float((np.abs(roi - t) < 0.10 * sep).mean())
            rows.append({"system": sysn, "roi": i, "x0": x0, "y0": y0, "x1": x1, "y1": y1, "scale_bar_um": bar_um, "scale_bar_px": bar_px,
                         "um_per_px": um_px, "roi_area_um2": (x1 - x0) * (y1 - y0) * um_px**2, "otsu_threshold": t, **fr,
                         "lambda_fft_um": wavelength(gray[y0:y1, x0:x1], um_px), "lambda_ness_um": NESS_LAMBDA[sysn],
                         "lambda_over_px": NESS_LAMBDA[sysn] / um_px})
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results" / "p01_ness_lamellae.csv", index=False)
    pd.set_option("display.width", 220)
    print(df[["system", "roi", "um_per_px", "roi_area_um2", "A_light_otsu", "A_light_minus10", "A_light_plus10", "lambda_fft_um",
              "lambda_ness_um", "lambda_over_px"]].round(3).to_string(index=False))
    overlay(gray, df)


def overlay(gray, df) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from scipy.ndimage import gaussian_filter

    fig, axes = plt.subplots(2, 1, figsize=(16, 15))
    axes[0].imshow(gray, cmap="gray")
    for r in df.itertuples():
        axes[0].add_patch(Rectangle((r.x0, r.y0), r.x1 - r.x0, r.y1 - r.y0, fill=False, ec="red", lw=1.5))
        axes[0].text(r.x0, r.y0 - 4, f"{r.system}#{r.roi}", color="red", fontsize=7)
    axes[0].set_title("Ness 2026 Fig. 5 with lamellar ROIs (red)")
    axes[0].axis("off")
    seg = np.zeros_like(gray)
    for r in df.itertuples():
        roi = gaussian_filter(gray[r.y0:r.y1, r.x0:r.x1], 1.0)
        seg[r.y0:r.y1, r.x0:r.x1] = np.where(roi > r.otsu_threshold, 2, 1)
    axes[1].imshow(gray, cmap="gray")
    axes[1].imshow(np.ma.masked_equal(seg, 0), cmap="coolwarm", alpha=0.55, vmin=1, vmax=2)
    axes[1].set_title("Otsu segmentation inside ROIs (red = light phase, blue = dark phase)")
    axes[1].axis("off")
    (HERE / "figures").mkdir(exist_ok=True)
    fig.savefig(HERE / "figures" / "fig_p01_roi_overlays.png", dpi=130, bbox_inches="tight")


if __name__ == "__main__":
    main()
