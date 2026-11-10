"""Plot the inlier residual histograms from the outlier experiment.

Usage: python tools/plot_residuals.py results/seed42/residual_hist.csv docs/residuals.png
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main(csv_path: str, out_path: str) -> None:
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    width = float(rows[1]["bin_start_px"]) - float(rows[0]["bin_start_px"])
    bins = [float(r["bin_start_px"]) for r in rows[:-1]]  # last row = everything beyond the range
    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=150)
    styles = {"squared": "#d95f02", "huber": "#7570b3", "cauchy": "#1b9e77"}
    for name, color in styles.items():
        share = [float(r[name]) * 100 for r in rows[:-1]]
        tail = float(rows[-1][name]) * 100
        ax.step(bins, share, where="post", color=color, lw=1.8,
                label=f"{name} loss  ({tail:.1f}% of clean matches > {bins[-1] + width:.0f} px)")
    ax.set_xlabel("reprojection error of the clean matches after BA (px)")
    ax.set_ylabel("share of clean matches per 0.1 px (%)")
    ax.set_title("5% bad matches: a squared loss lets them pull the solution; robust losses ignore them")
    ax.title.set_fontsize(10)
    ax.set_xlim(0, bins[-1] + width)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path)
    print("wrote", out_path)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
