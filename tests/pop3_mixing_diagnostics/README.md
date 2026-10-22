# Pop III mixing comparison

[REPORT.md](REPORT.md) explains the results; [diagnostics.pdf](diagnostics.pdf)
contains nine pages of comparison figures. Inputs are the existing completed
`../pop3_mixing0`, `../pop3_mixing0.5`, and `../pop3_mixing1.0` experiments.
The analysis does not run or change B-POP.

For smaller figures, open the [individual-panel index](panels/README.md).
It provides a description under 100 words, four suggested starting plots, and
61 single-panel PDFs grouped by subject, together with the top three results.
Descriptions and population/normalization notes appear in Markdown only;
PDFs contain the plots, axis labels, and legends. The combined PDF is also
regenerated with this presentation.

Regenerate only the small PDFs, preserving the original combined PDF:

```bash
python /path/to/pop3_mixing_diagnostics/plot_diagnostics.py --split-only
```

## Reproduce figures and report

Copy this entire diagnostic directory, including its CSV, NPZ, JSON, and
`snapshots/` files. From any working directory:

```bash
python /path/to/pop3_mixing_diagnostics/plot_diagnostics.py
python /path/to/pop3_mixing_diagnostics/write_report.py
```

Dependencies: Python 3.11, NumPy, pandas, and Matplotlib. The project environment
is `/home/chun-lung.chan/.conda/envs/popIII/bin/python`; exact analysis and plot
versions are recorded in `manifest.json` and `plot_validation.json`.
No raw run catalogues or project Python modules are needed for these steps.
Only PDF figure output is generated. For a copied bundle, set `MPLCONFIGDIR`
to a writable temporary cache directory if its default location is unsuitable.
Within this checkout the default cache is project-local `tmp/pop3-mixing-mpl/`.

## Reanalyse the raw experiments

This requires the original full repository, shared redshift table, and the three
run directories with their raw catalogues and stellar input datasets:

```bash
rtk proxy /home/chun-lung.chan/.conda/envs/popIII/bin/python -B \
  tests/pop3_mixing_diagnostics/analyze.py --chunk-size 100000
```

This overwrites this directory's derived analysis products. It reads all chain
events in bounded chunks, independently hashes/counts the raw files, verifies
input equivalence, and saves per-chain summaries. The three chain files total
about 13.5 GB; rerendering the PDF does not rescan them. Analysis progress can
be redirected to a log under the repository's `tmp/` directory.

Then regenerate the report and PDF with the commands above. Validate the PDF:

```bash
rtk proxy gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage \
  tests/pop3_mixing_diagnostics/diagnostics.pdf
```

Run the bounded raw-sample and full aggregate checks with:

```bash
rtk proxy /home/chun-lung.chan/.conda/envs/popIII/bin/python -B \
  tests/pop3_mixing_diagnostics/check_analysis.py
```

## Files and definitions

- `summary.csv`: exact candidate allocation, terminal yield, non-merger record
  counts, recorded event/chain counts, chain lengths, mass thresholds, runtime.
- `settings.csv`: every parsed configuration macro, including inactive values.
- `statistics.csv`: finite/non-finite counts, mean, extrema, strict mass-threshold
  counts, zeros, negatives, and histogram underflow/overflow for each quantity.
- `histograms.npz`: common bin edges and exact bin counts for all plotted groups.
- `pop3_mixing*_terminal_quantiles.csv`: exact finite-value terminal percentiles.
- `pop3_mixing*_chains.npz`: one entry per chain with at least one recorded event;
  chain ID, length, maximum remnant mass, first primary mass, last remnant mass,
  first/last event cosmic ages in years, maximum generation counter, environment
  code, mapped terminal ID, and whether that ID occurs in the terminal catalogue.
  Environment codes are 1=isolated, 2=young, 3=globular, 4=nuclear (0=aggregate
  label, unused for chain arrays). There are no isolated entries in chain files.
- `statuses.csv`: exact status counts, with terminal and recorded-event scopes
  kept separate. `none`, `young`, `globular`, `nuclear` retain their source labels.
- `snapshots/`: input headers, source/Makefile snapshots, build paths, logs,
  metallicity allocation, and sampled natal-retention CDFs. Headers are the
  available current-directory evidence, not independently archived compiler inputs.
- `manifest.json`: raw source paths relative to the original repository, SHA-256,
  sizes, independent line counts, analysis versions, and chunk size.
- `aggregate_checksums.json`: checksums for portable plotting inputs; verified
  before rendering. `validation.json` records the full input/schema audit.
- `plot_validation.json`: PDF page count and rendered text-bound checks.
- `panels/validation.json`: individual PDF names, dimensions, sizes, checksums,
  and source-page text-bound checks. `panels/pdf_validation.json` records the
  single-page and Ghostscript checks for the exported files.
- `reproduction_validation.json`: historical copied-bundle exact-render
  comparison, plus the current caption-free PDF hash and Ghostscript checks.
- `analysis_checks.json`: sample chunk-independence and full aggregate checks.
- `check_analysis.py`: bounded sample checks for chunk-independent histogram
  accumulation and independent consistency checks of full per-chain aggregates.

Terminal yields use the exact logged allocations, not merger/(merger+nonmerger).
Main and chain time columns are interpreted from the current source writers.
The strict printed-age count `chain_events_before_present` can differ slightly
from the total because six-significant-digit output rounds near the present age.
Distribution densities divide by all finite values and bin width. Out-of-range
values are counted, not silently renormalized away; log panels omit nonpositive
values. Physical-bound flags are retained and reported rather than repaired.
Mass thresholds use strict `>` comparisons. All figure statistics use full data;
no display samples, event bootstrap, smoothing, or confidence intervals are used.

## PDF page guide

1. Input mass, formation-time, and natal-kick spectra; unconditioned mixtures.
2. Candidate-normalized yields, counts, and dynamical terminal statuses.
3. Terminal masses, mass ratio, effective spin, delay, recycling, redshift, eccentricity.
4. All recorded chain-event distributions, including recoil/escape speed.
5. Recorded event timescales and cosmic ages.
6. Terminal remnant mass by environment.
7. Recorded generation-zero primary masses and mass ratios by cluster type.
8. Per-chain lengths, maximum recorded masses, and growth.
9. Additional terminal redshift, spin, recoil, and cosmic-age distributions.

The recorded chain-event list can omit recoil-ending events. It and the terminal
catalogue overlap and must not be summed as disjoint event populations. Mixing
does not act on the isolated channel; the experiments have one unpaired
stochastic realization per setup. These diagnostics do not calibrate cosmic rates.
