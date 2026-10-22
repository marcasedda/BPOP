#!/usr/bin/env python3
"""Reproduce the diagnostic PDF and individual panels from bundled aggregates."""
from pathlib import Path
import argparse
from contextlib import nullcontext
import hashlib
import json
import os

BASE = Path(__file__).resolve().parent
# All temporary state belongs in the session's project-local tmp directory.
# A copied portable bundle may override this path with MPLCONFIGDIR.
os.environ.setdefault("MPLCONFIGDIR", str(BASE.parents[1]/"tmp"/"pop3-mixing-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.text import Text
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

RUNS = ("pop3_mixing0", "pop3_mixing0.5", "pop3_mixing1.0")
LABELS = ("mixing = 0", "mixing = 0.5", "mixing = 1")
COLORS = ("C0", "C1", "C2")
ENVS = ("none", "young", "globular", "nuclear")
ENV_LABELS = ("Isolated", "Young", "Globular", "Nuclear")
AXLABELS = {
    "m1": r"Primary mass [$M_\odot$]", "m2": r"Secondary mass [$M_\odot$]",
    "remnant": r"Remnant mass [$M_\odot$]", "q": r"Mass ratio $\min(m_1,m_2)/\max(m_1,m_2)$",
    "chi_eff": r"Effective spin $\chi_{\rm eff}$", "delay": "Elapsed delay [Gyr]",
    "generation": "Recorded previous-merger count", "zmerger": "Merger redshift",
    "zformation": "Formation redshift", "eccentricity": "Binary eccentricity",
    "spin1": "Primary dimensionless spin", "remnant_spin": "Remnant dimensionless spin",
    "recoil": r"GW recoil [km s$^{-1}$]", "recoil_escape": r"GW recoil / escape speed",
    "formation_time": "Formation cosmic age [Gyr]", "merger_time": "Merger cosmic age [Gyr]",
    "kick": r"Natal kick [km s$^{-1}$]", "chain_length": "Recorded events per chain",
    "growth": "Maximum recorded remnant / first primary",
    "formation_delay": "BH formation time [Gyr]", "capture_delay": "Binary-single capture time [Gyr]",
    "threebody_delay": "Three-body formation time [Gyr]", "friction_delay": "Dynamical friction time [Gyr]",
    "hardening_delay": "Hardening time [Gyr]", "pairing_delay": "Pairing time [Gyr]",
    "inspiral_delay": "GW inspiral time [Gyr]",
}
LOG = {"m1", "m2", "remnant", "delay", "formation_time", "merger_time", "kick",
       "recoil", "recoil_escape", "generation", "chain_length", "growth",
       "formation_delay", "capture_delay", "threebody_delay", "friction_delay",
       "hardening_delay", "pairing_delay", "inspiral_delay"}

PANEL_GROUPS = (
    ("inputs", ("single_mixed_mass", "single_mixed_formation_time", "single_mixed_kick",
                "mixture_mass", "mixture_formation_time", "mixture_kick")),
    ("yields", ("terminal_yield", "recorded_events_per_candidate", "recorded_chains_per_candidate",
                "terminal_count", "recorded_event_count", "terminal_status")),
    ("terminal", ("m1", "m2", "remnant", "q", "chi_eff", "delay", "generation", "zmerger", "eccentricity")),
    ("chain_events", ("m1", "m2", "remnant", "q", "chi_eff", "delay", "generation", "recoil_escape", "eccentricity")),
    ("timescales", ("formation_delay", "capture_delay", "threebody_delay", "friction_delay",
                    "hardening_delay", "pairing_delay", "inspiral_delay", "formation_time", "merger_time")),
    ("terminal_by_environment", ENVS),
    ("first_events", ("young_m1", "globular_m1", "nuclear_m1", "young_q", "globular_q", "nuclear_q")),
    ("per_chain", ("length", "maximum_mass", "growth", "young_maximum_mass",
                   "globular_maximum_mass", "nuclear_maximum_mass")),
    ("terminal_extra", ("zformation", "recoil", "spin1", "remnant_spin", "formation_time", "merger_time")),
)
SHORT_SUMMARY = (
    "Each PDF isolates one comparison across mixing 0, 0.5, and 1. Input panels show the BH "
    "reservoirs; yield panels count terminal mergers; chain panels describe repeated mergers. "
    "Increasing mixing shifts the terminal median remnant mass from 314 to 105 to 85 solar masses. "
    "Recorded chain events fall from 24.06 to 13.70 to 4.24 million, mainly through shorter nuclear "
    "recycling histories. Terminal counts change much less. These are individual stochastic runs; "
    "chain records omit some recoil-ending events, and extreme masses can reach the imposed growth limit."
)


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split-only", action="store_true",
                        help="export individual panels while preserving the original combined PDF")
    args = parser.parse_args()
    output_dir = BASE/"panels"
    output_dir.mkdir(exist_ok=True)
    checksums = json.loads((BASE/"aggregate_checksums.json").read_text())
    for name, expected in checksums.items():
        assert checksum(BASE/name) == expected, (name, "checksum mismatch")
    table = pd.read_csv(BASE/"statistics.csv").set_index(["run", "scope", "environment", "metric"])
    summary = pd.read_csv(BASE/"summary.csv").set_index(["run", "environment"])
    statuses = pd.read_csv(BASE/"statuses.csv")
    data = np.load(BASE/"histograms.npz", allow_pickle=False)
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 9,
                         "xtick.labelsize": 9, "ytick.labelsize": 9,
                         "pdf.fonttype": 42, "axes.grid": True, "grid.alpha": .2})
    page_validation = []
    panel_validation = []
    page_notes = []

    def histogram(run, scope, env, metric):
        row = table.loc[(run, scope, env, metric)]
        key = row.histogram_key
        counts, edges = data[key], data[key+"_edges"]
        assert counts.sum()+row.underflow+row.overflow == row.n
        return counts, edges, int(row.n)

    def panel(ax, scope, env, metric):
        occupied = []
        for run, label, color in zip(RUNS, LABELS, COLORS):
            if (run, scope, env, metric) not in table.index:
                continue
            counts, edges, n = histogram(run, scope, env, metric)
            ax.stairs(counts/(max(n, 1)*np.diff(edges)), edges, label=label, color=color)
            idx = np.flatnonzero(counts)
            if len(idx):
                occupied.extend([edges[idx[0]], edges[idx[-1]+1]])
        ax.set_xlabel(AXLABELS[metric])
        ax.set_ylabel("Probability density")
        if metric in LOG:
            ax.set_xscale("symlog", linthresh=1) if metric == "generation" else ax.set_xscale("log")
        if occupied and min(occupied) < max(occupied):
            ax.set_xlim(min(occupied), max(occupied))
        if metric in ("m1", "m2", "remnant", "generation", "chain_length", "growth", "delay"):
            ax.set_yscale("log")
        ax.legend(frameon=False)

    def save(pdf, fig, note):
        page_notes.append(note)
        fig.tight_layout(h_pad=1.6, w_pad=1.3)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        outside = []
        omitted_ticks = set()
        for ax in fig.axes:
            for axis, limits in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
                for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                    if not min(limits) <= tick.get_loc() <= max(limits):
                        omitted_ticks.update((id(tick.label1), id(tick.label2)))
        for artist in fig.findobj(Text):
            if not artist.get_visible() or not artist.get_text() or id(artist) in omitted_ticks:
                continue
            box = artist.get_window_extent(renderer)
            # Matplotlib creates tick artists beyond the displayed axis limits;
            # ignore those unless they are actually drawn inside the axes box.
            if artist.axes is not None and artist.get_clip_on():
                continue
            if box.x0 < -2 or box.y0 < -2 or box.x1 > fig.bbox.x1+2 or box.y1 > fig.bbox.y1+2:
                outside.append(artist.get_text())
        page_validation.append(dict(page=len(page_validation)+1, text_outside_canvas=outside))
        assert not outside, outside
        if pdf is not None:
            pdf.savefig(fig, bbox_inches="tight", pad_inches=.04)
        group, names = PANEL_GROUPS[len(page_validation)-1]
        assert len(names) == len(fig.axes)
        figtexts = [(text, text.get_visible()) for text in fig.texts]
        for text, _ in figtexts:
            text.set_visible(False)
        for ax, name in zip(fig.axes, names):
            # Descriptions are written to Markdown; PDFs contain axes and legends.
            # The original dense page remains available as the combined PDF.
            for other in fig.axes:
                other.set_visible(other is ax)
            original_legend = ax.get_legend()
            if original_legend is None:
                ax.legend(handles=[Line2D([], [], color=color, label=label)
                                   for color, label in zip(COLORS, LABELS)],
                          frameon=False, fontsize=8)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            bounds = ax.get_tightbbox(renderer)
            crop = bounds.transformed(fig.dpi_scale_trans.inverted()).padded(.06)
            filename = f"{len(page_validation):02d}_{group}_{name}.pdf"
            path = output_dir/filename
            fig.savefig(path, bbox_inches=crop,
                        metadata={"Title": f"Pop III {group}: {name}", "Creator": "plot_diagnostics.py",
                                  "CreationDate": None, "ModDate": None})
            panel_validation.append(dict(file=filename, group=group, quantity=name,
                                         width_inches=float(crop.width), height_inches=float(crop.height),
                                         bytes=path.stat().st_size, sha256=checksum(path)))
            if original_legend is None:
                ax.get_legend().remove()
        for ax in fig.axes:
            ax.set_visible(True)
        for text, visible in figtexts:
            text.set_visible(visible)
        plt.close(fig)

    combined = nullcontext(None) if args.split_only else PdfPages(
        BASE/"diagnostics.pdf", metadata={"Title": "Pop III mixing comparison",
                                         "Creator": "plot_diagnostics.py",
                                         "CreationDate": None, "ModDate": None})
    with combined as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
        for col, metric in enumerate(("m1", "formation_time", "kick")):
            for kind, color in (("single", "C0"), ("mixed", "C2")):
                counts, edges, n = histogram(kind, "input_spectrum", "all", metric)
                axes[0, col].stairs(counts/(n*np.diff(edges)), edges, label=kind+" spectrum", color=color)
            for mix, label, color in zip((0, .5, 1), LABELS, COLORS):
                a, edges, na = histogram("single", "input_spectrum", "all", metric)
                b, _, nb = histogram("mixed", "input_spectrum", "all", metric)
                density = ((1-mix)*a/na+mix*b/nb)/np.diff(edges)
                axes[1, col].stairs(density, edges, label=label, color=color)
            for ax in axes[:, col]:
                ax.set_xscale("log")
                ax.set_yscale("log")
                ax.set_xlabel(AXLABELS[metric].replace("Primary", "Input BH"))
                ax.set_ylabel("Probability density")
                ax.legend(frameon=False)
                a, edges, _ = histogram("single", "input_spectrum", "all", metric)
                b, _, _ = histogram("mixed", "input_spectrum", "all", metric)
                idx = np.flatnonzero(a+b)
                ax.set_xlim(edges[idx[0]], edges[idx[-1]+1])
        save(pdf, fig, "Single/mixed panels compare the supplied reservoir components. Mixture panels show unconditioned combinations; actual sampling also requires natal-kick retention.")

        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
        x = np.arange(4)
        for i, (run, label, color) in enumerate(zip(RUNS, LABELS, COLORS)):
            rows = summary.loc[run].reindex(ENVS)
            axes[0, 0].bar(x+(i-1)*.24, rows.terminal_yield, .24, color=color, label=label)
            axes[0, 1].bar(x+(i-1)*.24, rows.recorded_chain_events/rows.candidates, .24, color=color)
            axes[0, 2].bar(x+(i-1)*.24, rows.recorded_chains/rows.candidates, .24, color=color)
            axes[1, 0].plot((0, .5, 1)[i], summary.loc[(run, "all"), "terminal_mergers"], "o", color=color)
            axes[1, 1].plot((0, .5, 1)[i], summary.loc[(run, "all"), "recorded_chain_events"], "o", color=color)
        for ax in axes[0]:
            ax.set_xticks(x, ENV_LABELS, rotation=15)
            ax.set_xlabel("Environment")
        axes[0, 0].set_ylabel("Terminal mergers / candidates")
        axes[0, 0].legend(frameon=False)
        axes[0, 1].set_ylabel("Recorded chain events / candidates")
        axes[0, 2].set_ylabel("Recorded chains / candidates")
        for ax in axes[1, :2]:
            ax.set_xlabel("Mixing parameter")
            ax.set_xticks([0, .5, 1])
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axes[1, 0].set_ylabel("Terminal merger count")
        axes[1, 1].set_ylabel("Recorded chain event count")
        dynamic = statuses[(statuses.scope == "terminal") & (statuses.environment != "none")]
        status_names = sorted(dynamic.status.unique())
        for i, (run, color) in enumerate(zip(RUNS, COLORS)):
            totals = dynamic[dynamic.run == run].groupby("status")["count"].sum().reindex(status_names, fill_value=0)
            axes[1, 2].bar(np.arange(len(status_names))+(i-1)*.24, totals/totals.sum(), .24, color=color)
        axes[1, 2].set_xticks(np.arange(len(status_names)), status_names, rotation=30, ha="right")
        axes[1, 2].set_ylabel("Fraction of dynamical terminal outcomes")
        save(pdf, fig, "Exact logged allocations set yield denominators. Recorded chain events and terminal outcomes are distinct populations.")

        for scope, metrics, note in (
            ("terminal", ("m1", "m2", "remnant", "q", "chi_eff", "delay", "generation", "zmerger", "eccentricity"),
             "Terminal catalogue; all environments. Delay = merger cosmic age minus formation cosmic age."),
            ("chain_event", ("m1", "m2", "remnant", "q", "chi_eff", "delay", "generation", "recoil_escape", "eccentricity"),
             "All recorded dynamical chain events; event-weighted. Recoil-ending events may be absent from the writer."),
            ("chain_event", ("formation_delay", "capture_delay", "threebody_delay", "friction_delay", "hardening_delay", "pairing_delay", "inspiral_delay", "formation_time", "merger_time"),
             "Recorded chain timescales converted from years to Gyr. Log axes omit nonpositive values; omissions are counted in statistics.csv."),
        ):
            fig, axes = plt.subplots(3, 3, figsize=(12, 10.5))
            for ax, metric in zip(axes.flat, metrics):
                panel(ax, scope, "all", metric)
            save(pdf, fig, note+" Densities use all finite values, including out-of-range values in their denominator.")

        fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
        for ax, env, name in zip(axes.flat, ENVS, ENV_LABELS):
            panel(ax, "terminal", env, "remnant")
            ax.set_ylabel(name+"\nProbability density")
        save(pdf, fig, "Terminal remnant masses separated by environment; common bin edges across all three runs.")

        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
        for col, (env, name) in enumerate(zip(ENVS[1:], ENV_LABELS[1:])):
            for row, metric in enumerate(("m1", "q")):
                panel(axes[row, col], "chain_first", env, metric)
                axes[row, col].set_ylabel(name+"\nProbability density")
        save(pdf, fig, "Recorded generation-zero events; excludes first mergers omitted by the recoil-before-write branch.")

        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
        for ax, metric in zip(axes[0], ("chain_length", "remnant", "growth")):
            panel(ax, "per_chain", "all", metric)
        for ax, env, name in zip(axes[1], ENVS[1:], ENV_LABELS[1:]):
            panel(ax, "per_chain", env, "remnant")
            ax.set_ylabel(name+"\nProbability density")
        save(pdf, fig, "One entry per chain with at least one recorded event. Remnant denotes its maximum recorded remnant mass.")

        fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
        for ax, metric in zip(axes.flat, ("zformation", "recoil", "spin1", "remnant_spin", "formation_time", "merger_time")):
            panel(ax, "terminal", "all", metric)
        save(pdf, fig, "Additional terminal distributions; all environments. Natal Fuller spins do not constrain recycled-primary spins.")

    validation = {"pages": page_validation,
                                                          "matplotlib": matplotlib.__version__,
                                                          "descriptions_in_markdown_only": True,
                                                          "aggregate_checksum_files_verified": len(checksums)}
    if not args.split_only:
        (BASE/"plot_validation.json").write_text(json.dumps(validation, indent=2)+"\n")
    validation["panels"] = panel_validation
    (output_dir/"validation.json").write_text(json.dumps(validation, indent=2)+"\n")
    assert len(SHORT_SUMMARY.split()) < 100
    index = ["# Individual diagnostic panels", "", SHORT_SUMMARY, "",
             "## Top three results", "",
             "Results below follow mixing **0 → 0.5 → 1**.", "",
             "1. **Lower remnant masses:** median terminal remnant mass falls "
             "**314 → 105 → 85 solar masses**. "
             "[Mass comparison](03_terminal_remnant.pdf).",
             "2. **Less nuclear recycling:** recorded nuclear chain events fall "
             "**23.94 → 13.60 → 4.16 million**, while nuclear terminal yields stay "
             "near **89–93%**. Mixing changes repeated growth much more than terminal "
             "merger production. [Event yields](02_yields_recorded_events_per_candidate.pdf).",
             "3. **Rare long chains dominate:** the longest **1%** of recorded chains "
             "supply **99.42% → 99.15% → 97.75%** of all recorded events; median chain "
             "length is **one** in every run. [Chain lengths](08_per_chain_length.pdf).", "",
             "Start with [terminal yield](02_yields_terminal_yield.pdf), "
             "[terminal remnant masses](03_terminal_remnant.pdf), "
             "[recorded event counts](02_yields_recorded_event_count.pdf), and "
             "[chain lengths](08_per_chain_length.pdf).", ""]
    for (group, _), note in zip(PANEL_GROUPS, page_notes):
        index.extend(["## "+group.replace("_", " ").capitalize(), "", note, ""])
        index.extend(f"- [{row['quantity'].replace('_', ' ')}]({row['file']})"
                     for row in panel_validation if row["group"] == group)
        index.append("")
    (output_dir/"README.md").write_text("\n".join(index))
    print(f"Wrote {len(panel_validation)} individual panels to {output_dir}")
    if not args.split_only:
        print(f"Wrote {BASE/'diagnostics.pdf'} ({len(page_validation)} pages)")


if __name__ == "__main__":
    main()
