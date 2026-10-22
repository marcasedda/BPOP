"""Draw the conceptual B-POP workflow described in ../BPOP_WORKFLOW.md.

Run with Python + Matplotlib; the only figure artifact is BPOP_WORKFLOW.pdf.
Arrows show logical dependencies, rather than the order of every source call.
"""

import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCRATCH = ROOT / "tmp" / "bpop-workflow-figure"
SCRATCH.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(SCRATCH / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon


INK = "#253448"
BLUE = "#1f77b4"
ORANGE = "#b9681d"
GREEN = "#26764b"
GRAY = "#657180"
fig, ax = plt.subplots(figsize=(18, 14))
fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
ax.set(xlim=(0, 18), ylim=(0, 14))
ax.axis("off")
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
checks = []


def text(x, y, value, size=11, **kwargs):
    kwargs.setdefault("color", INK)
    artist = ax.text(x, y, value, fontsize=size, **kwargs)
    return artist


def box(x, y, w, h, heading, body, color=GRAY, fill="#f5f7fa", size=11):
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.015,rounding_size=0.09",
        linewidth=1.3, edgecolor=color, facecolor=fill, zorder=2,
    )
    ax.add_patch(patch)
    heading_artist = text(x + 0.18, y + h - 0.19, heading, 12.5,
                          va="top", weight="bold", zorder=3)
    body_artist = text(x + 0.18, y + h - 0.53, body, size,
                       va="top", linespacing=1.5, zorder=3)
    checks.append((patch, [heading_artist, body_artist]))


def decision(cx, cy, w, h, label, color):
    patch = Polygon([(cx, cy + h / 2), (cx + w / 2, cy),
                     (cx, cy - h / 2), (cx - w / 2, cy)],
                    facecolor="white", edgecolor=color, linewidth=1.4, zorder=2)
    ax.add_patch(patch)
    artist = text(cx, cy, label, 11, ha="center", va="center", zorder=3)
    checks.append((patch, [artist]))


def arrow(points, color=GRAY, dashed=False):
    for a, b in zip(points[:-2], points[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=1.35,
                linestyle="--" if dashed else "-", zorder=1)
    ax.add_patch(FancyArrowPatch(
        points[-2], points[-1], arrowstyle="-|>", mutation_scale=13,
        linewidth=1.35, color=color, linestyle="--" if dashed else "-",
        shrinkA=0, shrinkB=2, zorder=1,
    ))


# Shared setup: candidate fractions describe simulation coverage.
box(0.45, 12.1, 10.6, 1.2,
    "1  Configure and compile the C++ population simulator",
    "Precomputed SEVN / MOBSE stellar catalogues + cluster tables + cosmic-age grid\n"
    "Set input_params.h: paths, candidate counts, cluster physics, spins; build BPOP.x")
box(0.45, 10.65, 10.6, 1.05,
    "2  Allocate candidates and sample their formation conditions",
    r"Isolated / young / globular / nuclear channels; draw $z_{\rm for}$ and $Z$, then group by metallicity")
arrow([(5.75, 12.1), (5.75, 11.7)])
arrow([(2.55, 10.65), (2.55, 10.2)], BLUE)
arrow([(8.125, 10.65), (8.125, 10.2)], ORANGE)

# Isolated lane.
box(0.45, 9.1, 4.2, 1.1, "3a  Isolated binary evolution",
    "Draw a complete upstream binary record:\ncomponent masses, orbit and delay; assign spins",
    BLUE, "#eef5fb", size=10.7)
box(0.45, 7.35, 4.2, 1.3, "4a  Use the imported merger delay",
    r"$t_{\rm merge}=t_{\rm for}+\Delta t_{\rm input}$" + "\n"
    "Keep correlations from stellar binary evolution.", BLUE, "#eef5fb")
arrow([(2.55, 9.1), (2.55, 8.65)], BLUE)
decision(2.55, 6.5, 3.6, 1.05, "Merger before\nthe present age?", BLUE)
arrow([(2.55, 7.35), (2.55, 7.025)], BLUE)
box(0.45, 5.0, 4.2, 0.9, "5a  Isolated merger properties",
    "Masses, spins, remnant and effective spin", BLUE, "#eef5fb", size=10.5)
arrow([(2.55, 5.975), (2.55, 5.9)], BLUE)
text(2.82, 5.94, "Yes", 9.5, color=BLUE, va="center")
arrow([(2.55, 5.0), (2.55, 3.0)], BLUE)

# Dynamical lane: two distinct retention tests and evolving cluster state.
box(5.65, 9.1, 4.95, 1.1, "3b  Young / globular / nuclear clusters",
    "Draw cluster mass and radius; estimate escape speed.\n"
    "Sample BHs retained after natal kicks; assign spins.", ORANGE, "#fff5eb", size=10.7)
box(5.65, 7.35, 4.95, 1.3, "4b  Assemble and evolve a dynamical binary",
    "Birth / sinking → binary formation → pairing / hardening\n"
    "→ GW inspiral; check ejection versus in-cluster merger.\n"
    "Update the cluster as time advances.", ORANGE, "#fff5eb", size=10.6)
arrow([(8.125, 9.1), (8.125, 8.65)], ORANGE)
decision(8.125, 6.5, 3.6, 1.05, "Merger before\nthe present age?", ORANGE)
arrow([(8.125, 7.35), (8.125, 7.025)], ORANGE)
box(5.65, 4.75, 4.95, 1.1, "5b  Calculate the merger remnant",
    "Numerical-relativity fits → mass, spin and GW recoil.\n"
    "Record the dynamical chain step.", ORANGE, "#fff5eb", size=10.7)
arrow([(8.125, 5.975), (8.125, 5.85)], ORANGE)
text(8.39, 5.92, "Yes", 9.5, color=ORANGE, va="center")
decision(8.125, 3.8, 4.95, 1.05,
         "Retained after recoil, and\nanother merger possible?", ORANGE)
arrow([(8.125, 4.75), (8.125, 4.325)], ORANGE)
arrow([(10.6, 3.8), (11.65, 3.8), (11.65, 8.0), (10.6, 8.0)], ORANGE)
text(11.37, 5.85, "Yes: reuse remnant + draw companion", 10,
     color=ORANGE, rotation=90, ha="center", va="center")
arrow([(8.125, 3.275), (8.125, 3.0)], ORANGE)
text(8.29, 3.15, "No: end chain", 9.5, color=ORANGE)

# Both failed-clock branches feed the non-merger catalogue.
arrow([(4.35, 6.5), (5.1, 6.5), (5.1, 3.0)], GRAY)
arrow([(6.325, 6.5), (5.1, 6.5)], GRAY)
text(4.6, 6.65, "No", 9.5, ha="center")
text(5.9, 6.65, "No", 9.5, ha="center")
text(4.87, 4.5, "Non-mergers", 9.5, rotation=90, ha="center", va="center")

# Interface groups all products without implying that every branch writes each.
box(0.45, 1.35, 10.6, 1.65, "6  Archive the C++ catalogues and run configuration",
    "", fill="#f5f7fa")
text(0.68, 2.36, "Catalogue.txt", 11.2, weight="bold")
text(0.68, 2.08, "Isolated mergers + selected terminal chain outcomes", 10.5)
text(6.05, 2.36, "Catalogue_multiple_dyn.txt", 11.2, weight="bold")
text(6.05, 2.08, "Intermediate dynamical chain records", 10.5)
text(0.68, 1.76, "Larger_than_tH.txt", 11.2, weight="bold")
text(0.68, 1.48, "Failed / unfinished candidates; efficiency denominators", 10.5)
text(6.05, 1.76, "retention_Z<Z>.dat", 11.2, weight="bold")
text(6.05, 1.48, "Natal-kick cumulative retention tables", 10.5)

# Right lane: catalogue weighting remains a separate modelling stage.
arrow([(11.05, 2.2), (12.2, 2.2), (12.2, 12.7), (12.75, 12.7)], GREEN)
text(12.04, 10.1, "Catalogue inputs", 9.5, color=GREEN,
     rotation=90, ha="center", va="center")
box(12.75, 12.1, 4.8, 1.2, "7  Set Python cosmological models",
    "Channel formation histories, metallicity probabilities,\n"
    "normalizations, external isolated efficiency; Planck18.", GREEN, "#eff8f2", size=10.6)
box(12.75, 9.7, 4.8, 1.85, "8  Estimate merger efficiency by metallicity",
    r"$\eta_c(Z)$: mergers per $M_\odot$ of formed stars" + "\n"
    "Isolated: external efficiency table.\n"
    "Clusters: merger / non-merger yields, BH and binary\n"
    "fractions, and natal-kick retention.", GREEN, "#eff8f2", size=10.8)
arrow([(15.15, 12.1), (15.15, 11.55)], GREEN)
box(12.75, 6.95, 4.8, 2.2, "9  Weight formation histories and delays",
    "Sum metallicity bins and integrate formation epochs:\n"
    r"$\mathcal{R}_c(z)=10^9\sum_k\int dz_f\,\psi_c\,\eta_{c,k}\,P_{c,k}\,q_{c,k}\,|dt/dz_f|$" + "\n"
    r"Implemented delay-bin weight: $q=dN/(N\,\Delta t_{\rm merger})$." + "\n"
    r"Output: channel rates $\mathcal{R}_c(z)$ [Gpc$^{-3}$ yr$^{-1}$]." + "\n"
    "See the normalization note below.", GREEN, "#eff8f2", size=10.6)
arrow([(15.15, 9.7), (15.15, 9.15)], GREEN)
box(12.75, 4.75, 4.8, 1.65, "10  Convert rates to observer counts",
    r"$\Delta N_c\simeq\mathcal{R}_c(z)\,\Delta V_{\rm com}\,T_{\rm obs}/(1+z)$" + "\n"
    "Comoving volume × observer interval;\n"
    "the redshift factor converts source to observer time.", GREEN, "#eff8f2", size=10.8)
arrow([(15.15, 6.95), (15.15, 6.4)], GREEN)
box(12.75, 2.65, 4.8, 1.55, "11  Draw mock cosmic event catalogues",
    "Match channel, metallicity and merger-redshift bins;\n"
    "sample catalogue rows with replacement.\n"
    "Write channel files and Catalogue_global.txt.", GREEN, "#eff8f2", size=10.8)
arrow([(15.15, 4.75), (15.15, 4.2)], GREEN)
text(12.85, 2.37, "Sparse candidate bins can reduce written counts.\n"
     "Detector selection requires an additional model.", 10, va="top", linespacing=1.5)

# Source-review distinctions central to interpreting this checkout.
ax.plot([0.45, 17.55], [1.05, 1.05], color="#d5dce3", lw=1)
text(0.45, 0.83,
     "Reading the diagram: arrows show logical dependencies; candidate fractions set sampling coverage, not cosmic channel fractions.  BH = black hole; GW = gravitational wave.",
     10.2, va="top")
text(0.45, 0.56,
     "Rate normalization: the implemented delay-bin denominator includes non-mergers, as does the dynamical efficiency; their combined normalization needs scientific review.",
     10.2, va="top")
text(0.45, 0.29,
     r"Pop III trial: $Z=10^{-11}$, generic formation histories and zero isolated efficiency; calibrated Pop III rates remain unfinished.  Source: BPOP_WORKFLOW.md, §§1–11.",
     10.2, va="top")

# Check actual rendered geometry rather than accepting clipped text.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for patch, artists in checks:
    outer = patch.get_window_extent(renderer)
    for artist in artists:
        inner = artist.get_window_extent(renderer)
        if inner.width and inner.height and (
            inner.x0 < outer.x0 or inner.x1 > outer.x1
            or inner.y0 < outer.y0 or inner.y1 > outer.y1
        ):
            raise RuntimeError(f"Text exceeds its node: {artist.get_text()}")
for artist in ax.texts:
    extent = artist.get_window_extent(renderer)
    if not fig.bbox.contains(extent.x0, extent.y0) or not fig.bbox.contains(extent.x1, extent.y1):
        raise RuntimeError(f"Text exceeds the page: {artist.get_text()}")
visible_text = [artist for artist in ax.texts if artist.get_text()]
for i, first in enumerate(visible_text):
    for second in visible_text[i + 1:]:
        if first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)):
            raise RuntimeError(f"Text overlaps: {first.get_text()} / {second.get_text()}")

output = ROOT / "BPOP_WORKFLOW.pdf"
fig.savefig(output, bbox_inches="tight", pad_inches=0.08,
            metadata={"Title": "B-POP workflow", "Subject": "Conceptual flow from BPOP_WORKFLOW.md",
                      "Creator": "figures/plot_bpop_workflow.py", "CreationDate": None,
                      "ModDate": None})
plt.close(fig)
print(f"Wrote {output}; checked {len(ax.texts)} text elements and {len(checks)} nodes.")
