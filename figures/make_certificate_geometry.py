#!/usr/bin/env python3
"""Conceptual geometry of the infinite-lattice a posteriori certificate.

All drawing coordinates are illustrative. They are not numerical certificate
data, the 109-chart atlas, an energy surface, or computed transition paths.
The inclusion condition in panel (b) is the dimension-independent Hilbert-ball
test used in the article; only one transverse direction is drawn.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch
import numpy as np

OUT = Path(__file__).resolve().parent
plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.linewidth": .7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "path",
    "savefig.facecolor": "white",
})

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#007D59"
INK = "#222222"
GRAY = "#888888"


def line(ax, x, y, **kwargs):
    ax.plot(x, y, solid_capstyle="round", **kwargs)


def label(ax, text, xy, xytext, color=INK, ha="left"):
    ax.annotate(text, xy=xy, xytext=xytext, color=color,
                ha=ha, va="center", fontsize=8.5,
                arrowprops={"arrowstyle": "-", "lw": .6,
                            "color": color, "shrinkB": 3},
                annotation_clip=False)


fig, (a, b) = plt.subplots(1, 2, figsize=(6.0, 3.55),
                          gridspec_kw={"width_ratios": [1.03, 1]})
fig.subplots_adjust(left=.025, right=.99, top=.875, bottom=.18, wspace=.14)

# (a) A projected coordinate picture. Rectangle outlines are coordinate
# domains only; their full transverse fibers are suppressed in this panel.
a.set_xlim(-1.05, 1.15)
a.set_ylim(-.94, .98)
a.set_aspect("equal")
a.axis("off")
fig.text(.028, .94, "(a) Exact cuts and coordinate coverage",
         fontsize=10, va="top", fontweight="bold")

charts = [(-.98, -.68, 1.18, 1.36),
          (-.35, -.74, 1.22, 1.48),
          (.25, -.64, .85, 1.28)]
for x, y, w, hgt in charts:
    a.add_patch(Rectangle((x, y), w, hgt, fill=False, ec="#A5A5A5",
                          lw=.75, ls=(0, (2.2, 2.2)), zorder=1))

# A_0 is a bounded coordinate rectangle; the two exact linear cuts trim it.
a.add_patch(Rectangle((-.85, -.55), 1.85, 1.1, facecolor="#F4F4F4",
                      edgecolor="#A0A0A0", lw=.7, zorder=0))
clip_x = (.8 - .55) / .65
domain = [(-.85, -.55), (clip_x, -.55), (1, -.15),
          (1, .15), (clip_x, .55), (-.85, .55)]
a.add_patch(Polygon(domain, closed=True, facecolor="#E7F2EC",
                    edgecolor=INK, lw=1.0, zorder=2))

xx = np.linspace(.09, 1.1, 100)
for sign in [1, -1]:
    yy = sign * (.8 - .65 * xx)
    line(a, xx, yy, color=BLUE, lw=1.05, ls=(0, (4, 2)), zorder=3)
    sx, sy = .65, sign * (.8 - .65 * .65)
    # Only the section of the exact cap in this reaction-coordinate plane
    # is drawn, along the exact plane ell=0.
    tt = np.linspace(-.105, .105, 2)
    tangent = np.array([1, -sign*.65])
    tangent /= np.linalg.norm(tangent)
    line(a, sx + tt*tangent[0], sy + tt*tangent[1],
         color=BLUE, lw=3.1, zorder=5)
    a.plot(sx, sy, marker="D", markersize=4.6, color=ORANGE,
           markeredgecolor="white", markeredgewidth=.4, zorder=7)
    a.text(sx-.025, sy+sign*.115, r"$s_+$" if sign==1 else r"$s_-$",
           color=ORANGE, ha="center", va="center", fontsize=10)

# Coordinate face a_1=.2 of the left chart lies strictly inside a neighbor.
line(a, [.2, .2], [-.27, .27], color=GREEN, lw=2,
     ls=(0, (1.7, 1.7)), zorder=5)
a.plot(-.46, 0, "o", color=INK, markersize=4, zorder=6)
a.text(-.48, -.13, r"$m_0$", fontsize=10, ha="center")
a.text(-.35, .34, r"$\mathcal{D}$", fontsize=12, color=INK)

# This is a schematic projection of a matching upper path, not an orbit.
path = np.array([[-.46, 0], [-.18, .025], [.13, .10], [.40, .25],
                 [.65, .3775], [.92, .70]])
line(a, path[:, 0], path[:, 1], color=ORANGE, lw=1.4, zorder=6)
a.add_patch(FancyArrowPatch(path[-2], path[-1], arrowstyle="-|>",
                            mutation_scale=8, color=ORANGE, lw=1.4,
                            zorder=6))
a.text(.73, .81, r"toward $\mathcal{K}$", fontsize=8.3, color=ORANGE,
       ha="center")
label(a, "covered coordinate face", (.2, -.16), (-.99, -.84), color=GREEN)
label(a, "retained old wall", (-.85, .18), (-1.00, .84))
a.text(1.03, .40, r"$\ell_+=0$", color=BLUE, fontsize=8.5,
       rotation=-33, ha="center")
a.text(1.03, -.40, r"$\ell_-=0$", color=BLUE, fontsize=8.5,
       rotation=33, ha="center")
a.text(.91, -.78, r"$(a_1,a_2)$", fontsize=9, ha="center")

# (b) A coordinate/transverse slice. The orange segment is the full allowed
# low-sublevel ball in a one-direction depiction. The strict norm inequality
# below certifies inclusion in ALL transverse directions in Y.
b.set_xlim(-.98, 1.22)
b.set_ylim(-1.18, 1.35)
b.axis("off")
fig.text(.55, .94, "(b) Entire low-sublevel fibers",
         fontsize=10, va="top", fontweight="bold")

xi, xj = -.78, -.05
right_i, right_j = .25, 1.12
yi, yj = 0, .16
Ri, Rj, rlow = .85, .90, .38
b.add_patch(Rectangle((xi, yi-Ri), right_i-xi, 2*Ri,
                      facecolor="#F4F4F4", edgecolor=GRAY,
                      lw=.8, ls=(0, (2.2, 2.2)), zorder=1))
b.add_patch(Rectangle((xj, yj-Rj), right_j-xj, 2*Rj,
                      facecolor="#EDF5FA", edgecolor=BLUE,
                      lw=1.0, zorder=2))
b.text(-.61, .63, r"$U_i$", color=GRAY, fontsize=11)
b.text(.75, .82, r"$U_j$", color=BLUE, fontsize=11)

astar = right_i
line(b, [astar, astar], [yi-Ri, yi+Ri], color=GRAY, lw=1,
     ls=(0, (3, 2)), zorder=3)
line(b, [astar, astar], [yi-rlow, yi+rlow], color=ORANGE,
     lw=5, zorder=4)
line(b, [astar-.07, astar+.07], [yi-rlow, yi-rlow],
     color=ORANGE, lw=1.3, zorder=4)
line(b, [astar-.07, astar+.07], [yi+rlow, yi+rlow],
     color=ORANGE, lw=1.3, zorder=4)
b.plot(astar, yi, marker="o", markersize=3, color=ORANGE, zorder=5)
b.plot(astar, yj, marker="+", markersize=6, markeredgewidth=1.2,
       color=BLUE, zorder=5)
label(b, r"$a_*\in\partial A_i\cap A_j$", (astar, .85),
      (-.84, 1.13), ha="left")
label(b, r"$\|w\|\leq\bar{r}_i$", (astar, -.38),
      (-.86, -.54), color=ORANGE)
b.text(.39, -.04, r"$y_i$", color=ORANGE, fontsize=9)
b.text(.40, .23, r"$y_j$", color=BLUE, fontsize=9)

# Radius of the neighboring ball and the transverse slice indication.
b.add_patch(FancyArrowPatch((.91, yj), (.91, yj+Rj), arrowstyle="<->",
                            mutation_scale=7, color=BLUE, lw=.8, zorder=5))
b.text(.98, yj+Rj/2, r"$R_j$", color=BLUE, fontsize=9,
       va="center", ha="left")
line(b, [.55, 1.10], [-1.02, -1.02], color=INK, lw=.7)
b.add_patch(FancyArrowPatch((1.03, -1.02), (1.16, -1.02),
                            arrowstyle="-|>", mutation_scale=6,
                            color=INK, lw=.7))
b.text(.54, -1.16, "one coordinate direction", fontsize=8.2, ha="center")
b.text(-.88, -.98, "one direction in $Y$", fontsize=8.2,
       rotation=90, va="bottom", ha="center")

fig.text(.50, .105,
         r"$\bar{r}_i=\min(R_i,r_i^+),\qquad"
         r"\|y_i-y_j\|+\bar{r}_i<R_j$",
         ha="center", va="center", fontsize=10)
fig.text(.50, .041,
         "Conceptual schematic: no calibrated coordinates or certified path data.",
         ha="center", va="center", fontsize=8.4, color="#555555")

for suffix in ["pdf", "svg", "png"]:
    fig.savefig(OUT / f"certificate_geometry.{suffix}", dpi=300,
                metadata={"Creator": "Matplotlib; make_certificate_geometry.py"})
plt.close(fig)
