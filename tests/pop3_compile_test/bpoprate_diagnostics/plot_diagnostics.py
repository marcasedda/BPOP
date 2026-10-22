"""Plot archived BPOPRate trial outputs without rerunning the simulation."""
from pathlib import Path
import hashlib
import json
import os

BASE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(BASE / "tmp" / "mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


def main():
    data = BASE / "data"
    manifest = json.loads((BASE / "input_manifest.json").read_text())
    for name, expected in manifest.items():
        assert hashlib.sha256((data / name).read_bytes()).hexdigest() == expected, name
    metadata = json.loads((data / "run_metadata.json").read_text())
    channels = ("nuclear", "globular", "young", "none")
    labels = ("Nuclear", "Globular", "Young", "Isolated")
    colors = ("#0072B2", "#D55E00", "#009E73", "#999999")
    rates = [np.loadtxt(data / f"Merger_rates_a_{c}_massive_no.txt") for c in channels]
    for rate in rates:
        assert rate.shape == (metadata["settings"]["grid_size"] - 1, 3)
        assert np.isfinite(rate).all() and (rate[:, 1:] >= 0).all()
        np.testing.assert_array_equal(rate[:, 0], rates[0][:, 0])
        assert (np.diff(rate[:, 0]) > 0).all()
    efficiency = np.loadtxt(data / "efficiency_Z_SET2a_massive_no.txt", dtype=str)
    assert tuple(efficiency[:, 2]) == channels
    eta = efficiency[:, 1].astype(float)
    assert np.isfinite(eta).all() and (eta >= 0).all()
    catalogue = np.loadtxt(data / "Catalogue_global.txt", dtype=str, skiprows=1, ndmin=2)
    assert catalogue.shape == (16, 22)
    zgw = catalogue[:, 1].astype(float)
    assert np.isfinite(zgw).all() and (zgw >= 0).all() and (zgw <= 15).all()
    counts = {c: int(np.count_nonzero(catalogue[:, 12] == c)) for c in channels}
    assert counts == dict(nuclear=3, globular=2, young=11, none=0)

    plt.rcParams.update({"font.size": 10, "axes.labelsize": 11,
                         "legend.fontsize": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.1),
                             gridspec_kw={"width_ratios": [1.65, 1, 1.2]})
    ax = axes[0]
    for rate, label, color in zip(rates[:3], labels[:3], colors[:3]):
        ax.step(rate[:, 0], rate[:, 1], where="mid", color=color,
                linewidth=0.85, alpha=0.85, label=label)
    total = np.sum([rate[:, 1] for rate in rates], axis=0)
    ax.step(rates[0][:, 0], total, where="mid", color="black", linewidth=1,
            alpha=0.65, label="Total")
    ax.set(xlabel=r"Merger redshift $z$",
           ylabel=r"Merger rate density [yr$^{-1}$ Gpc$^{-3}$]", xlim=(0, 15))
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1]
    ax.bar(np.arange(4), eta / 1e-4, color=colors, width=0.65)
    ax.set_xticks(np.arange(4), labels, rotation=25, ha="right")
    ax.set(xlabel="Environment", ylabel=r"Merger efficiency [$10^{-4}$ $M_\odot^{-1}$]")
    ax.set_ylim(0, float(eta.max() / 1e-4) * 1.18)
    for i, value in enumerate(eta / 1e-4):
        ax.annotate(f"{value:.2f}" if value else "0 (disabled)", (i, value),
                    xytext=(0, 5), textcoords="offset points", ha="center", fontsize=9)
    ax.grid(axis="y", alpha=0.2)
    ax.set_axisbelow(True)

    ax = axes[2]
    bins = np.linspace(0, 15, 11)
    samples = [zgw[catalogue[:, 12] == c] for c in channels[:3]]
    hist, _, _ = ax.hist(samples, bins=bins, stacked=True, color=colors[:3],
                          label=[f"{label} ({counts[c]})" for label, c in zip(labels[:3], channels[:3])],
                          edgecolor="white", linewidth=0.6)
    assert int(hist[-1].sum()) == len(catalogue)
    ax.set(xlabel=r"Sampled merger redshift $z$", ylabel="Sampled sources per bin", xlim=(0, 15))
    ax.set_ylim(0, max(5, float(hist[-1].max()) + 2.5))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    ax.set_axisbelow(True)

    # An annotation records trial assumptions; no figure or panel titles.
    footer = (r"Execution trial: generic SFR; $Z=10^{-11}$; isolated channel disabled. "
              r"$T_{\rm obs}=0.01$ yr; 16 sampled sources."
              "\nRaw rate bins shown; sparse input population; no uncertainty estimates.")
    fig.text(0.5, 0.035, footer, ha="center", va="bottom", fontsize=9, color="0.3")
    fig.tight_layout(rect=(0, 0.16, 1, 1), w_pad=2)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        for axis in (ax.xaxis, ax.yaxis):
            lower, upper = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks():
                if not lower <= tick.get_loc() <= upper:
                    tick.label1.set_visible(False)
                    tick.label2.set_visible(False)
    for text in fig.findobj(matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if box.width and box.height:
                assert fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1), text.get_text()
    output = BASE.parent / "bpoprate_diagnostics.pdf"
    fig.savefig(output, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    summary = {"sample_counts": counts, "rate_rows": len(rates[0]),
               "total_rate_maximum": float(total.max()),
               "grid_size": metadata["settings"]["grid_size"],
               "formation_grid_size": metadata["settings"]["formation_grid_size"],
               "checks": ["input hashes", "finite rates and efficiencies", "matching redshift grids",
                          "catalogue shape and channel counts", "histogram count conservation",
                          "all visible text within figure bounds"]}
    (BASE / "plot_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(output)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
