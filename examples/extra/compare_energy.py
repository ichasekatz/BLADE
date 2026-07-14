"""Compare energy per atom between two folders of relaxed structures.

Scans both folders for POSCAR/CONTCAR + energy files, matches by reduced formula,
and plots parity + residual histogram. Works with any folder structure —
Materials Project downloads, BLADE Comps dirs, or custom folders.

Optionally compares folder_1 against DFT values from a materials_project_summary.xlsx.

Output:
    <output_dir>/energy_comparison.png
    <output_dir>/energy_comparison.xlsx

Usage:
    pixi run python examples/extra/compare_energy.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymatgen.core import Structure

# ------------------------------------------------------------------
# Paths — adjust to your environment
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
files_dir = path1 / "Files"

folder_1 = files_dir / "MaterialsProject_Comps_GRACE"
folder_2 = files_dir / "MaterialsProject_Comps_ORB"
label_1  = folder_1.name.replace("MaterialsProject_Comps_", "") or "Folder1"
label_2  = folder_2.name.replace("MaterialsProject_Comps_", "") or "Folder2"

output_dir = files_dir
output_dir.mkdir(parents=True, exist_ok=True)

# Optional: path to MP summary xlsx for DFT comparison (set None to skip)
dft_xlsx = folder_1 / "materials_project_summary.xlsx"

# ------------------------------------------------------------------
# Scanner: works for any folder with POSCAR/CONTCAR + energy files
# ------------------------------------------------------------------
def scan_energies(folder: Path, label: str) -> pd.DataFrame:
    rows: list[dict] = []
    for poscar_path in sorted(folder.rglob("POSCAR")):
        subdir = poscar_path.parent
        energy_file = subdir / "energy"
        if not energy_file.exists():
            continue

        struct_path = subdir / "CONTCAR"
        if not struct_path.exists():
            struct_path = poscar_path
        try:
            s = Structure.from_file(struct_path)
        except Exception as e:
            print(f"  Skipping {poscar_path}: {e}")
            continue

        try:
            total = float(energy_file.read_text().strip().split()[0])
            epa = total / len(s)
        except Exception:
            continue

        rows.append({
            "formula": s.composition.reduced_formula,
            f"epa_{label}": epa,
            "n_atoms": len(s),
            "n_elements": len(s.composition.get_el_amt_dict()),
            "folder": folder.name,
            "subdir": str(subdir),
        })

    df = pd.DataFrame(rows)
    print(f"{label}: {len(df)} structures with energy files in {folder.name}")
    return df


def plot_parity_and_residual(
    x: pd.Series, y: pd.Series,
    xlabel: str, ylabel: str,
    title_parity: str, title_residual: str,
    color_by: pd.Series | None = None,
    out_png: Path | None = None,
    out_xlsx: Path | None = None,
) -> None:
    diff = y - x
    mae  = diff.abs().mean()
    rmse = (diff ** 2).mean() ** 0.5
    mean_d = diff.mean()
    print(f"MAE:  {mae:.4f} eV/atom")
    print(f"RMSE: {rmse:.4f} eV/atom")
    print(f"Mean: {mean_d:+.4f} eV/atom")

    cmap = plt.get_cmap("tab10")
    if color_by is not None:
        vals = sorted(color_by.unique())
        colors = {v: cmap(i % 10) for i, v in enumerate(vals)}
    else:
        colors = None

    _pct_lo, _pct_hi = 2, 98
    all_vals = pd.concat([x, y]).dropna()
    _lo = np.percentile(all_vals, _pct_lo) - 0.5
    _hi = np.percentile(all_vals, _pct_hi) + 0.5
    in_view = x.between(_lo, _hi) & y.between(_lo, _hi)
    n_out = (~in_view).sum()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    if colors and color_by is not None:
        for v in vals:
            mask = (color_by == v) & in_view
            ax.scatter(x[mask], y[mask], color=colors[v], alpha=0.7, s=30,
                       label=f"{v}-element", zorder=3)
            out_mask = (color_by == v) & ~in_view
            if out_mask.any():
                ax.scatter(x[out_mask].clip(_lo, _hi), y[out_mask].clip(_lo, _hi),
                           color=colors[v], alpha=0.4, s=60, marker="^",
                           edgecolors="k", linewidths=0.5, zorder=4)
    else:
        ax.scatter(x[in_view], y[in_view], alpha=0.6, s=30, zorder=3)
    lims = [_lo, _hi]
    ax.plot(lims, lims, "k--", lw=1, zorder=2)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(title_parity, fontsize=14)
    if colors:
        ax.legend(fontsize=10, framealpha=0.8)
    note = f"\n(▲ {n_out} outliers clipped)" if n_out else ""
    ax.text(0.05, 0.95,
            f"MAE = {mae:.4f} eV/atom\nRMSE = {rmse:.4f} eV/atom\nn = {len(x)}{note}",
            transform=ax.transAxes, va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    ax = axes[1]
    _d_lo = np.percentile(diff.dropna(), _pct_lo)
    _d_hi = np.percentile(diff.dropna(), _pct_hi)
    ax.hist(diff.clip(_d_lo, _d_hi), bins=40, color="steelblue", edgecolor="white", alpha=0.85)
    ax.axvline(0,      color="k",   lw=1.5, linestyle="--", label="zero")
    ax.axvline(mean_d, color="red", lw=1.5, linestyle="--", label=f"mean = {mean_d:+.4f}")
    ax.set_xlabel(f"{ylabel.split('(')[0].strip()} − {xlabel.split('(')[0].strip()} (eV/atom)", fontsize=13)
    ax.set_ylabel("Count", fontsize=13)
    ax.set_title(title_residual, fontsize=14)
    ax.legend(fontsize=10)

    plt.tight_layout()
    if out_png:
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        print(f"Saved: {out_png}")
    plt.close()

    if out_xlsx is not None:
        result = pd.DataFrame({"x_epa": x, "y_epa": y, "diff_epa": diff})
        result.to_excel(out_xlsx, index=False)
        print(f"Saved: {out_xlsx}")


# ==============================================================================
# Section 1: folder_1 vs DFT (optional)
# ==============================================================================
df1 = scan_energies(folder_1, label_1)

if dft_xlsx and dft_xlsx.exists():
    print(f"\n--- {label_1} vs DFT ---")
    mp_df = pd.read_excel(dft_xlsx)
    mp_df["formula"] = mp_df["formula"].astype(str).str.strip()
    mp_df["epa_dft"] = pd.to_numeric(mp_df.get("energy_per_atom", mp_df.get("energy_per_atom_dft")), errors="coerce")

    df1_dedup_dft = df1.sort_values(f"epa_{label_1}").drop_duplicates("formula")
    mp_dedup = mp_df[["formula", "epa_dft"]].dropna().drop_duplicates("formula")
    merged = df1_dedup_dft.merge(mp_dedup, on="formula", how="inner")
    merged = merged.dropna(subset=["epa_dft", f"epa_{label_1}"])
    print(f"Matched {len(merged)} unique formulas")

    plot_parity_and_residual(
        x=merged["epa_dft"],
        y=merged[f"epa_{label_1}"],
        xlabel="MP DFT energy (eV/atom)",
        ylabel=f"{label_1} energy (eV/atom)",
        title_parity=f"Energy: {label_1} vs DFT",
        title_residual=f"Residuals: {label_1} − DFT",
        color_by=merged["n_elements"],
        out_png=output_dir / f"{label_1}_vs_DFT.png",
        out_xlsx=output_dir / f"{label_1}_vs_DFT.xlsx",
    )
else:
    print(f"\nSection 1 skipped — no DFT xlsx at {dft_xlsx}")

# ==============================================================================
# Section 2: folder_1 vs folder_2
# ==============================================================================
if folder_2.exists():
    print(f"\n--- {label_1} vs {label_2} ---")
    df2 = scan_energies(folder_2, label_2)

    # Deduplicate by formula (keep lowest energy per formula) to avoid many-to-many join
    df1_dedup = df1.sort_values(f"epa_{label_1}").drop_duplicates("formula")
    df2_dedup = df2.sort_values(f"epa_{label_2}").drop_duplicates("formula")

    merged2 = df1_dedup.merge(
        df2_dedup[["formula", f"epa_{label_2}", "n_elements"]],
        on="formula", how="inner", suffixes=("", "_2"),
    )
    merged2 = merged2.dropna(subset=[f"epa_{label_1}", f"epa_{label_2}"])
    print(f"Matched {len(merged2)} unique formulas")

    n_el_col = "n_elements" if "n_elements" in merged2.columns else "n_elements_2"

    plot_parity_and_residual(
        x=merged2[f"epa_{label_1}"],
        y=merged2[f"epa_{label_2}"],
        xlabel=f"{label_1} energy (eV/atom)",
        ylabel=f"{label_2} energy (eV/atom)",
        title_parity=f"Energy: {label_2} vs {label_1}",
        title_residual=f"Residuals: {label_2} − {label_1}",
        color_by=merged2[n_el_col],
        out_png=output_dir / f"{label_1}_vs_{label_2}.png",
        out_xlsx=output_dir / f"{label_1}_vs_{label_2}.xlsx",
    )
else:
    print(f"\nSection 2 skipped — {folder_2} not found.")
