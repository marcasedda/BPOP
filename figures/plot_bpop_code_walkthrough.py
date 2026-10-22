"""Illustrate the code paths described in ../BPOP_CODE_WALKTHROUGH.md.

Panel (a): the dynamical merger clock and hierarchical loop in
src/Main_Bpop_dev.cpp. Panel (b): which C++ output columns the production
Python rate script reads, with the review's issue IDs (section 8) marked.

Run with Python + Matplotlib; writes ../BPOP_CODE_WALKTHROUGH_FIGURE.pdf.
Pass --png to also write a raster preview under ../tmp/bpop-code-figure/.
The figure is a schematic: stage widths are not to scale.
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCRATCH = ROOT / "tmp" / "bpop-code-figure"
SCRATCH.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(SCRATCH / "mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

# Reference categorical slots 1-3 (validated all-pairs) + reserved "critical".
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8984"
CPP = "#2a78d6"        # C++ simulator stages
PY = "#eb6834"         # Python rate stages
FILE = "#1baf7a"       # files passed between stages
CRIT = "#d03b3b"       # issue markers (always icon + ID, never colour alone)
TINT = {CPP: "#e6effa", PY: "#fcebe3", FILE: "#e3f5ee", MUTED: "#f3f3f1"}

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
fig, ax = plt.subplots(figsize=(16, 11.4))
fig.patch.set_facecolor(SURFACE)
fig.subplots_adjust(left=0.005, right=0.995, bottom=0.005, top=0.995)
ax.set(xlim=(0, 16), ylim=(0, 11.4))
ax.axis("off")
checks = []  # (container patch, [text artists]) for the overflow check


def text(x, y, value, size=10.5, **kw):
    kw.setdefault("color", INK)
    kw.setdefault("zorder", 4)
    return ax.text(x, y, value, fontsize=size, **kw)


def box(x, y, w, h, heading, body="", color=MUTED, size=10, center=False):
    patch = FancyBboxPatch((x, y), w, h,
                           boxstyle="round,pad=0.01,rounding_size=0.08",
                           lw=1.4, edgecolor=color, facecolor=TINT[color],
                           zorder=2)
    ax.add_patch(patch)
    arts = []
    if center:
        arts.append(text(x + w / 2, y + h / 2, heading, size, ha="center",
                         va="center", weight="bold", linespacing=1.35))
    else:
        arts.append(text(x + 0.14, y + h - 0.13, heading, size + 0.5,
                         va="top", weight="bold"))
        if body:
            arts.append(text(x + 0.14, y + h - 0.45, body, size, va="top",
                             color=INK2, linespacing=1.45))
    checks.append((patch, arts))
    return patch


def arrow(p, q, color=MUTED, dashed=False, rad=0.0, label=None, lpos=0.5,
          loff=(0, 0.14)):
    ax.add_patch(FancyArrowPatch(
        p, q, arrowstyle="-|>", mutation_scale=13, lw=1.4, color=color,
        linestyle="--" if dashed else "-", shrinkA=1, shrinkB=1,
        connectionstyle=f"arc3,rad={rad}", zorder=1))
    if label:
        x = p[0] + (q[0] - p[0]) * lpos + loff[0]
        y = p[1] + (q[1] - p[1]) * lpos + loff[1]
        text(x, y, label, 9.5, ha="center", va="center", color=INK2,
             style="italic")


def issue(x, y, tag, note, ha="left"):
    """Red '!' badge with the issue ID, plus a short ink-coloured note."""
    badge = text(x, y, f"! {tag}", 9.5, ha=ha, va="center", color="white",
                 weight="bold", zorder=6,
                 bbox=dict(boxstyle="round,pad=0.22", fc=CRIT, ec=CRIT))
    if note:
        fig.canvas.draw()
        bb = badge.get_window_extent().transformed(ax.transData.inverted())
        nx = bb.x1 + 0.1 if ha == "left" else bb.x0 - 0.1
        text(nx, y, note, 9.5, ha=ha, va="center", color=INK, zorder=6)


def panel_title(y, value):
    text(0.25, y, value, 14, weight="bold", va="center")


# ---------------------------------------------------------------------------
# (a) Dynamical merger clock and hierarchical chain (Main_Bpop_dev.cpp)
# ---------------------------------------------------------------------------
panel_title(11.05, "(a)  One dynamical source in Main_Bpop_dev.cpp: "
            "the merger clock and the hierarchical loop")
text(0.25, 10.68, "Schematic, stage widths not to scale.  time starts at "
     "t_form; each stage adds to it, and the cluster (M, R, σ, ρ, v_esc) is "
     "re-evolved after every stage (●).", 10, color=INK2, va="center")

stages = [
    ("sinking / BH birth", "max(t_DF, t_BH)", 2.05),
    ("binary forms", "min(t_capture, t_3bb)", 2.05),
    ("pairing", "t_12", 1.25),
    ("hardening", "5(m1+m2)/m_per · t_bbh", 2.35),
    ("GW inspiral", "t_GW(a_crit, e)", 1.85),
]
x0, y0, h = 0.55, 8.85, 0.85
x = x0
bounds = [x0]
for name, formula, w in stages:
    box(x, y0, w, h, f"{name}\n{formula}", color=CPP, size=9.8, center=True)
    x += w + 0.1
    bounds.append(x - 0.05)
x_end = x - 0.1

# time axis with re-evolution markers
arrow((x0, 8.6), (x_end + 0.25, 8.6), color=INK2)
for bx in bounds:
    ax.plot(bx, 8.6, "o", ms=7, color=CPP, mec=SURFACE, mew=2, zorder=5)
text(x0, 8.33, "t_form", 10, ha="center", va="center", color=INK2)
text(x_end, 8.33, "t_merge", 10, ha="center", va="center", color=INK2)

# core-collapse snap inside the binary-formation stage
xcc = bounds[1] + 1.15
ax.plot([xcc, xcc], [8.5, 9.95], color=INK2, lw=1.2, ls=(0, (3, 2)), zorder=3)
text(xcc - 0.08, 10.17, "t_form + t_cc: a stage that crosses core collapse\n"
     "restarts there and is recomputed", 9.5, ha="center", va="center",
     color=INK2)

# C1 sits over the hardening / inspiral stages (a_crit uses semi_ej)
issue(bounds[3] + 0.3, 10.1, "C1",
      "first merger: semi_ej → a_crit\nuses uninitialised mu_red")

# Decision: merges before today?
cx, cy = 12.55, y0 + h / 2
diamond = Polygon([(cx, cy + 0.62), (cx + 1.0, cy), (cx, cy - 0.62),
                   (cx - 1.0, cy)], fc="white", ec=CPP, lw=1.4, zorder=2)
ax.add_patch(diamond)
checks.append((diamond, [text(cx, cy, "time <\nt_today ?", 10,
                              ha="center", va="center", linespacing=1.3)]))
arrow((x_end, cy), (cx - 1.0, cy), color=CPP)

box(14.0, 9.7, 1.85, 0.95, "Larger_than_tH",
    "non-merger row", color=FILE, size=9.5)
arrow((cx + 1.0, cy), (14.0, 10.05), color=MUTED, label="no", lpos=0.45,
      loff=(0.0, 0.2))

box(13.75, 7.8, 2.1, 0.95, "SREM2",
    "remnant M, χ, recoil", color=CPP, size=9.5)
arrow((cx, cy - 0.62), (13.75, 8.27), color=CPP, rad=0.25, label="yes",
      lpos=0.4, loff=(-0.2, -0.05))

# Hierarchical loop
box(0.55, 6.3, 10.95, 1.45, "Hierarchical loop  (:2646)",
    "remnant → new primary;  draw a secondary with kick ≤ current v_esc;\n"
    "t_DF(r_∞) + t_12 + hardening + t_GW;  SREM2;  nrecy += 1;  "
    "one row → Catalogue_multiple_dyn.txt", color=CPP, size=10)
arrow((13.75, 7.95), (11.5, 7.5), color=CPP, label="remnant retained",
      lpos=0.5, loff=(0.15, 0.2))
arrow((6.0, 7.75), (6.0, 8.52), color=CPP, dashed=True,
      label="repeat the clock", lpos=0.5, loff=(0.95, 0.0))

# Stop conditions
text(0.6, 5.85, "The chain stops when:", 10.5, va="center", weight="bold")
chips = ["highlander (next merger after today)", "GWrecoiled",
         "ejected", "inevap / dissolved", "M > msmbhmax", "nrecy > nbhs"]
cxp = 2.6
for c in chips:
    t = text(cxp, 5.85, c, 9.5, va="center", color=INK,
             bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=MUTED,
                       lw=1.1))
    fig.canvas.draw()
    bb = t.get_window_extent().transformed(ax.transData.inverted())
    cxp = bb.x1 + 0.42

box(11.85, 6.3, 4.0, 1.0, "Catalogue.txt",
    "one terminal summary row per source\n(absolute time in col. 14)",
    color=FILE, size=9.5)
arrow((11.5, 6.6), (11.85, 6.6), color=MUTED)
text(11.68, 6.42, "on stop", 9, ha="center", va="top", color=INK2,
     style="italic")
issue(11.9, 5.85, "C4", "highlander row mixes next pair\nwith the last merger time")
issue(0.6, 5.35, "C6", "a GW-recoiled final step is written to Catalogue.txt "
      "but not to Catalogue_multiple_dyn.txt")

ax.plot([0.25, 15.75], [4.95, 4.95], color="#e1e0d9", lw=1.2, zorder=0)

# ---------------------------------------------------------------------------
# (b) C++ -> Python handoff (BPOPRate_v2025.py, lib_BPOPRate.py)
# ---------------------------------------------------------------------------
panel_title(4.65, "(b)  What the production Python stage reads, "
            "and where to check it")

files = [
    ("Catalogue.txt", "24, 25 z_merge, z_form · 1 Z · 3 env\n"
                      "14 t_merge (absolute) · 15–17 cluster M, R", None, 3.3, 0.9),
    ("Larger_than_tH.txt", "13 Z · 14 env · 22 nrecy · 23 end time · 24 t_form\n"
                           "rate stage bins non-mergers by col. 23, not 24",
     "P1", 2.2, 0.9),
    ("retention_Z*.dat", "v,  CDF_single,  CDF_mixed\n"
                         "weighted by fmix = IDmod/4, not C++ mixing",
     "P2", 1.1, 0.9),
    ("eta_isolated_binaries_SEVN", "external η_iso(Z) table", None, 0.2, 0.7),
]
for name, body, tag, y, hh in files:
    box(0.55, y, 4.3, hh, name, body, color=FILE, size=9.3)
    if tag:
        issue(4.72, y + hh - 0.2, tag, "", ha="right")

stages_py = [
    ("efficiency  η(Z, env)",
     "fbh/m★ · N_mer/(N_mer+N_non) · fbin · ⟨f_ret⟩_M", 3.0),
    ("rate convolution  R(z)",
     "∫ ψ · η · (dN/N) · P(Z|z_f) · |dt/dz| dz_f / Δt", 1.75),
    ("bsmplr  mock events",
     "N = R ΔV T_obs/(1+z);  draw rows with replacement", 0.5),
]
for name, body, y in stages_py:
    box(6.55, y, 4.8, 1.0, name, body, color=PY, size=9.3)

outs = [("efficiency_Z_SET*.txt", 3.15), ("Merger_rates_*_<env>.txt", 1.9),
        ("Catalogue_<env>.txt\nCatalogue_global.txt", 0.55)]
for name, y in outs:
    box(12.4, y, 3.45, 0.75, name, color=FILE, size=9.3, center=True)

# file -> stage arrows
arrow((4.85, 3.85), (6.55, 3.65), color=MUTED)            # Catalogue -> eta
arrow((4.85, 3.6), (6.55, 2.4), color=MUTED)              # Catalogue -> rate
arrow((4.85, 3.4), (6.55, 1.0), color=MUTED)              # Catalogue -> smplr
arrow((4.85, 2.75), (6.55, 3.4), color=MUTED)             # LtH -> eta
arrow((4.85, 2.5), (6.55, 2.15), color=MUTED)             # LtH -> rate
arrow((4.85, 1.55), (6.55, 3.2), color=MUTED)             # retention -> eta
arrow((4.85, 0.55), (6.55, 3.05), color=MUTED)            # eta_iso -> eta
# stage -> stage and stage -> output
arrow((8.95, 3.0), (8.95, 2.75), color=PY)
arrow((8.95, 1.75), (8.95, 1.5), color=PY)
for (sy, oy) in [(3.5, 3.52), (2.25, 2.27), (1.0, 0.92)]:
    arrow((11.35, sy), (12.4, oy), color=PY)

issue(9.5, 0.22, "P3", "int() truncates expected counts N ≥ 1")

# Key: swatch + label, laid out left to right, right-aligned to the page.
items = [(CPP, "C++ stage"), (PY, "Python stage"), (FILE, "file")]
kx = 9.6
for col, lab in items:
    ax.add_patch(FancyBboxPatch((kx, 4.57), 0.3, 0.16,
                                boxstyle="round,pad=0.01", fc=TINT[col],
                                ec=col, lw=1.3, zorder=3))
    t = text(kx + 0.42, 4.65, lab, 9.5, va="center", color=INK2)
    fig.canvas.draw()
    kx = t.get_window_extent().transformed(ax.transData.inverted()).x1 + 0.35
issue(kx, 4.65, "C1", "issue ID in walkthrough §8")


# ---------------------------------------------------------------------------
# Overflow check: every text artist must sit inside its container.
# ---------------------------------------------------------------------------
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
inv = ax.transData.inverted()
problems = []
for patch, arts in checks:
    pb = patch.get_window_extent(renderer)
    for a in arts:
        tb = a.get_window_extent(renderer)
        if (tb.x0 < pb.x0 - 1 or tb.x1 > pb.x1 + 1 or
                tb.y0 < pb.y0 - 1 or tb.y1 > pb.y1 + 1):
            problems.append(a.get_text().splitlines()[0])
fb = fig.bbox
for a in ax.texts:
    tb = a.get_window_extent(renderer)
    if tb.x0 < fb.x0 or tb.x1 > fb.x1 or tb.y0 < fb.y0 or tb.y1 > fb.y1:
        problems.append("off-page: " + a.get_text().splitlines()[0])
if problems:
    print("Text overflow:", *problems, sep="\n  ")

out = ROOT / "BPOP_CODE_WALKTHROUGH_FIGURE.pdf"
fig.savefig(out, facecolor=SURFACE)
print("Wrote", out)
if "--png" in sys.argv:
    png = SCRATCH / "preview.png"
    fig.savefig(png, dpi=80, facecolor=SURFACE)
    print("Preview", png)
