# B-POP workflow figure

`plot_bpop_workflow.py` generates the single-page vector PDF
[`../BPOP_WORKFLOW.pdf`](../BPOP_WORKFLOW.pdf). Its text is distilled from
[`../BPOP_WORKFLOW.md`](../BPOP_WORKFLOW.md), sections 1–11, reviewed on
2026-10-12. It requires Python and Matplotlib, with no simulation data or
project modules. All figure content is embedded in the script.

From the repository root:

```bash
rtk proxy /home/chun-lung.chan/.conda/envs/popIII/bin/python figures/plot_bpop_workflow.py
```

Any Python environment with Matplotlib can be used. Paths resolve from the
script location, so invocation from another working directory also works.
The output is always placed one directory above the script. Matplotlib cache
storage is under the repository's `tmp/bpop-workflow-figure/`.

The diagram shows logical dependencies, not every source-call ordering. It
distinguishes simulation coverage from cosmic weighting, natal-kick retention
from remnant recoil retention, and representative merger summaries from
intermediate chain records. The normalization and Pop III notes preserve the
limitations documented in the source report. It is an explanatory schematic,
not a new simulation or validation of physical rates.

The script checks text against node/page bounds and checks text overlaps.
Ghostscript can validate the resulting PDF without creating an image preview:

```bash
rtk proxy gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage BPOP_WORKFLOW.pdf
```
