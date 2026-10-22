# BPOPRate trial diagnostics

[bpoprate_diagnostics.pdf](../bpoprate_diagnostics.pdf) contains three panels:

- Raw merger rate density versus redshift, for each enabled environment and
  their total. No extra smoothing or uncertainty estimate is applied.
- Dynamical merger efficiencies at the supplied metallicity. Isolated
  efficiency is zero because that channel was disabled in the trial.
- Stacked redshift counts for all 16 sampled sources: nuclear 3, globular 2,
  young 11. Counts use ten common bins over redshift 0–15.

These are diagnostics of an execution trial with generic star-formation
histories, fixed Z=1e-11, and Tobs=0.01 yr, not a calibrated Pop III prediction.
The rate tables have 499 rows from the original 500-point grids. The sparse
input population produces strongly varying rates between bins.

This directory contains all inputs needed to reproduce the figure. Copy the
whole directory and run, from any working directory. The PDF is written beside
this directory, one level above the script:

```bash
python /path/to/diagnostics/plot_diagnostics.py
```

Dependencies used: Python 3.11, NumPy 1.26.4, and Matplotlib 3.11.1. The script
resolves input/output paths relative to its own location, checks the SHA-256
hashes in `input_manifest.json`, validates data shape and histogram counts, and
saves only PDF. Cache files stay under this directory's `tmp/`.
`plot_validation.json` records the numerical and text-layout checks. To place
Matplotlib cache files in another temporary directory, set `MPLCONFIGDIR`.

The PDF also passed Ghostscript validation:

```bash
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage ../bpoprate_diagnostics.pdf
```
