#!/usr/bin/env python3
"""The confrontation assay, composed from the plate photographs.

Each plate carries its own controls: the two strains growing alone at the top and
the same two strains confronted below, so the comparison the scoring rests on is
visible within a single image. Obverse and reverse are both shown because several
phenotypes - diffusible pigment and vegetative-mycelium colour - are only readable
from the underside.

The dish is located by a Hough circle fit rather than a hand-drawn box, so the
crops are concentric and the scale bar follows from the fitted radius.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from PIL import Image, ImageOps
from skimage import feature, transform

ROOT = Path(__file__).resolve().parents[1]
PLATES = ROOT / "data" / "plates"
OUT = ROOT / "results" / "figures"

PLATE_MM = 90.0          # outer diameter of a standard Petri dish
MARGIN = 1.01            # crop just outside the fitted rim

# (obverse, reverse, strain A, strain B)
PAIRS = [
    ("20250111_123839.jpg", "20250111_123846.jpg", "S524", "S1350"),
    ("20250111_123904.jpg", "20250111_123910.jpg", "S713'W", "S1350"),
    ("20250111_123746.jpg", "20250111_123754.jpg", "MS3 11", "S1350"),
]

# where S1350 sits in the confronted pair, as a fraction of the cropped plate
MARKS = [(0.635, 0.515), (0.607, 0.500), (0.600, 0.490)]


def use_arial() -> str:
    for p in ("/mnt/c/Windows/Fonts/arial.ttf", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"):
        if Path(p).exists():
            fm.fontManager.addfont(p)
            bold = p.replace("arial.ttf", "arialbd.ttf").replace("Arial.ttf", "Arial_Bold.ttf")
            if Path(bold).exists():
                fm.fontManager.addfont(bold)
            return fm.FontProperties(fname=p).get_name()
    return "DejaVu Sans"


def fit_circle(im: Image.Image, s: int = 8) -> tuple[int, int, int]:
    g = im.convert("L")
    W, H = g.size
    sm = np.asarray(g.resize((W // s, H // s)), dtype=float) / 255.0
    edges = feature.canny(sm, sigma=3.0, low_threshold=0.06, high_threshold=0.18)
    radii = np.arange(int(0.30 * W / s), int(0.52 * W / s), 2)
    h = transform.hough_circle(edges, radii)
    _, cx, cy, r = transform.hough_circle_peaks(h, radii, total_num_peaks=1)
    return int(cx[0] * s), int(cy[0] * s), int(r[0] * s)


def plate(path: Path, mirror: bool = False):
    im = ImageOps.exif_transpose(Image.open(path)).convert("RGB")
    cx, cy, r = fit_circle(im)
    half = int(r * MARGIN)
    crop = im.crop((cx - half, cy - half, cx + half, cy + half))
    if mirror:
        crop = ImageOps.mirror(crop)
    return np.asarray(crop), (2 * r) / PLATE_MM, crop.size[0]


def main() -> None:
    font = use_arial()
    plt.rcParams.update({"font.family": font, "savefig.dpi": 300})

    fig, axes = plt.subplots(2, 3, figsize=(7.3, 5.15))
    fig.subplots_adjust(left=0.035, right=0.999, top=0.955, bottom=0.004,
                        wspace=0.012, hspace=0.012)

    tags = [["a", "b", "c"], ["d", "e", "f"]]
    for col, (front, back, a_lab, b_lab) in enumerate(PAIRS):
        for row, (fn, mir) in enumerate(((front, False), (back, True))):
            img, ppm, side = plate(PLATES / fn, mirror=mir)
            ax = axes[row][col]
            ax.imshow(img)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)

            ax.text(0.022, 0.978, tags[row][col], transform=ax.transAxes,
                    fontsize=9, fontweight="bold", va="top", ha="left", color="white",
                    path_effects=_halo())
            if row == 0:
                ax.set_title(f"{a_lab} × {b_lab}", fontsize=8, pad=3, color="black")
                # one arrowhead on S1350 in the confronted pair
                # short stub from directly below, so it clears the marker-pen
                # labels written to the right of the colony
                mx, my = MARKS[col]
                ax.annotate("", xy=(mx * side, (my + 0.030) * side),
                            xytext=(mx * side, (my + 0.115) * side),
                            arrowprops=dict(arrowstyle="-|>", color="white", lw=1.2,
                                            mutation_scale=8,
                                            path_effects=_halo(1.6)))
            if col == 0:
                ax.set_ylabel("obverse" if row == 0 else "reverse", fontsize=7.5,
                              color="black", labelpad=2)
            if row == 0 and col == 0:
                w = 10 * ppm
                x0, y0 = side * 0.055, side * 0.945
                ax.add_patch(Rectangle((x0, y0), w, side * 0.017, color="white",
                                       ec="black", lw=0.4, zorder=6))
                ax.text(x0 + w / 2, y0 - side * 0.012, "10 mm", ha="center", va="bottom",
                        fontsize=7, color="white", zorder=7, path_effects=_halo())

    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"figure_plates.{ext}", dpi=400, bbox_inches="tight",
                    facecolor="white", pad_inches=0.012)
    plt.close(fig)
    print(f"  wrote figure_plates ({font})")


def _halo(lw: float = 1.8):
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=lw, foreground="black")]


if __name__ == "__main__":
    main()
