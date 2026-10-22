# Project handoff

Updated: 2026-10-22. Repository: `/home/chun-lung.chan/project/popIII_study/BPOP`.

## Fork transfer (2026-10-22)

The working checkout is now `BPOP` on `juno_pop3`. Its remote remains
`git@github.com:marcasedda/BPOP.git`; nothing was pushed. The local fork at
`../BPOP_p3` is preserved unchanged. Its commit `447bbfd` (`added test for pop3`)
was fast-forwarded into this branch, followed by transfer of its working files.
The historical sections below describe the work as recorded in the source
clone; their old checkout paths and working-tree status are historical.

Source, configuration, scripts, and documentation are included in the transfer
commit. Raw simulation outputs, downloaded inputs, diagnostic data/figures,
snapshots, temporary environments, and compiler products are copied locally and
excluded from new commits. Root documentation PDFs are included. Existing
tracked files `build/BPOP.x`, `build/Main_Bpop_dev.o`, and `build/funcGR.o`
retain their copied modifications outside the commit; ignore rules do not
untrack them. Do not stage these artifacts with `git add -u` or `git add .`.

The destination's existing `graphify-out/` was preserved; the source graph and
checkout-specific `.git`, `.codex`, `.agents`, and `.aws` directories were not
copied. Transfer records are in local `tmp/fork-transfer-20261022/`.

Copied executables, generated build paths, and virtual environments may still
refer to `BPOP_p3`. Clean and rebuild the relevant target before future
simulation work from this checkout. Absolute paths in archived provenance
remain unchanged deliberately. This transfer does not rerun simulations,
rates, or plots, nor does it certify the existing scientific results.

Transfer validation: 15 incoming Python scripts passed AST parsing with the
`popIII` environment; all 31 files above 10 MiB were copied and ignored by Git.
The 53-file source/documentation selection is below 1 MiB in total. Source
clone status/remotes and all 24 existing destination graph file checksums were
preserved. New transfer edits pass whitespace checks; a full staged review
also reports 2,634 inherited whitespace findings in 13 unchanged incoming
files. Their formatting is preserved; details are in the local transfer log.

## Smaller mixing figures (2026-10-19)

Follow-up presentation change: the user requested descriptions in Markdown
instead of inside PDF plots, and asked for the top three results. Both the
combined PDF and all 61 individual panels are now regenerated without prose
footers. Axis labels and legends remain. `panels/README.md` holds the migrated
population/normalization notes under each topic and a linked top-three section:
terminal median remnant mass 314 -> 105 -> 85 Msun; nuclear chain event totals
23.94 -> 13.60 -> 4.16 million while terminal yields remain 89-93%; longest 1%
of recorded chains supply 99.42% / 99.15% / 97.75% of events. The short overall
description remains under 100 words. This supersedes the footer/original-PDF
preservation details in the initial split task below. Numerical aggregates are
unchanged. PDF validation checks that no description text remains in the PDFs;
prior copied-render reproduction evidence is retained as historical metadata.

The user found the combined diagnostic too large and requested smaller pieces
plus an explanation under 100 words. `plot_diagnostics.py --split-only` now
exports all 61 axes as individual single-page PDFs in
`tests/pop3_mixing_diagnostics/panels/`, preserving the original combined PDF.
The normal plotting command also emits these pieces. Each piece includes its
own legend and the source page's population/normalization note; only the selected
axes are drawn, with a cropped PDF page and wrapped context text.

`panels/README.md` contains the short overall explanation, four suggested
starting plots, and a topic-grouped complete index. Main README and report link
to it; `write_report.py` preserves the new index link on regeneration. No raw
analysis or simulations are needed. Key starting files are
`02_yields_terminal_yield.pdf`, `03_terminal_remnant.pdf`,
`02_yields_recorded_event_count.pdf`, and `08_per_chain_length.pdf`.
`panels/validation.json` records each file's dimensions, size, and checksum;
`panels/pdf_validation.json` records one-page/Ghostscript checks. The original
PDF checksum and all aggregate checksums remain unchanged. No additional image
format is generated; the plotting skill's PDF-only constraint is retained.

## Mixing-run diagnostics (2026-10-15)

Implemented the approved full comparison in `tests/pop3_mixing_diagnostics/`.
Read the existing `tests/pop3_mixing0`, `pop3_mixing0.5`, and `pop3_mixing1.0`
runs without changing their inputs, executables, or catalogues. Their current
headers differ only in mixing, each uses N=1,000,000, and local C++ main sources
plus all three stellar datasets match byte for byte. Header/compiler provenance
is limited to the available directory files and corroborating run logs.

Deliverables: `REPORT.md`, nine-page `diagnostics.pdf`, README, standalone
analysis/plot/report/check scripts, CSV tables, histogram NPZ, per-chain NPZ,
source/configuration/log snapshots, raw manifest/checksums, and validation JSONs.
The whole directory is about 7.5 MB; plot/report reproduction needs only this
bundle and Python packages. Raw catalogue analysis requires the original repo.

Key results (mixing 0 / 0.5 / 1):

- Terminal mergers: 396,658 / 363,491 / 373,946.
- Recorded chain events: 24,056,741 / 13,704,958 / 4,244,100; all 42,005,799
  records streamed, with independent newline-count reconciliation.
- Terminal median remnant masses: 313.928 / 104.848 / 85.218 Msun.
- Nuclear recorded events: 23,942,113 / 13,603,487 / 4,156,801, despite terminal
  nuclear yields remaining near 89–93%. Top 1% of recorded chains account for
  99.42% / 99.15% / 97.75% of recorded events.
- Mixing 0/0.5 reach the imposed msmbhmax=5e7 Msun growth stopping scale.

Catalogue semantics are explicitly audited: chain column 31 plus one maps to
terminal IDs; recoil exits before the chain writer can omit the final event.
Terminal and intermediate catalogues overlap and must not simply be added.
Times follow current writers (absolute cosmic ages, years); elapsed delay
subtracts formation age. A few printed ages round across the present-age cutoff.
Non-merger IDs have writer-specific offsets and are not complementary IDs.

All chain numeric columns are finite, chronological checks passed, and the
checked physical bounds passed. Non-finite terminal semimajor axes occur in
19 / 10 / 9 rows; many non-merger separations are non-finite. These columns are
counted explicitly and not used by the mass/spin/delay plots. No cosmic rates,
simulations, significance estimates, or independent-seed reruns were generated.

Validation: independent full row counts, exact terminal/log counts, 250,003-row
samples per run with 7,777-row alternate chunks, and independently grouped sample
chains compared to full arrays. Nine PDF pages passed text-bound and Ghostscript
checks. A copied bundle reproduced the identical PDF from another working
directory without raw data. `analysis_checks.json`, `plot_validation.json`, and
`reproduction_validation.json` record these checks. Temporary state is in root
`tmp/`; existing working changes were preserved.

Commands (from repository root; use the popIII Python environment):

```bash
rtk proxy /home/chun-lung.chan/.conda/envs/popIII/bin/python -B tests/pop3_mixing_diagnostics/plot_diagnostics.py
rtk proxy /home/chun-lung.chan/.conda/envs/popIII/bin/python -B tests/pop3_mixing_diagnostics/write_report.py
rtk proxy /home/chun-lung.chan/.conda/envs/popIII/bin/python -B tests/pop3_mixing_diagnostics/check_analysis.py
```

Full reanalysis uses `analyze.py --chunk-size 100000`, took about 310 seconds,
and overwrites derived products; rerun plot/report afterward. The plotting
script validates bundled aggregate checksums and emits PDF only.

## Pop III parameter recommendation audit (2026-10-14)

The user asked whether the original Pop III recommendations in the parameter
comments are applied. Compared the HEAD version of
`tests/pop3_compile_test/input_params_p3.h`, its current annotations, the local
driver/Makefile, and all three retained simulator configuration snapshots.
All 67 current macro definitions match HEAD and all three snapshots. No
parameter, source, executable, or run output was changed during this audit.

- Applied: `uppergap="no"`, `bhseed="no"`, `mixing=0.5` (one of the original
  suggested 0/0.5/1 choices), `spinlb="fuller"`, and `DynOvTot=0.9`.
- The original comment `Ex. for pop 3: DynOvTot 0.9 and pYC or pGC = DynOvTot`
  is not fully applied: current YC/GC/NC shares are 0.75/0.15/0.10. The code
  multiplies YC/GC shares by DynOvTot; they are conditional dynamical shares,
  not shares of the total population. The original example is ambiguous and
  does not establish a unique replacement configuration.
- The original Elbadry19 note suggested `sigma_metal_clu=0.5`; current value
  is 0.6, but the active Pop III branch bypasses that metallicity scatter.
  `SEVN="pop3"` selects Z=1e-11 despite `MonoZ="no"`; the explicit `mono_Z`
  override is inactive with the current continuous SFR and AS20 settings.
- The `MCL_popIII`/`RCL_popIII` fixed-cluster and `MBH_popIII` seed parameters
  are inactive with `TagR="AS20"` and `bhseed="no"`; their names alone do
  not enable the separate seed experiment.
- The root build uses `include/input_params.h`, which remains a different
  non-Pop-III configuration. Use the local test Makefile for the Pop III setup.

Validation used macro comparisons, current-source branch checks, archived logs,
and `make -C tests/pop3_compile_test -n`. The dry run schedules recompilation
after the annotation edits; no build or simulation was executed. The existing
graph was queried for orientation, but excludes the local experimental source.

## Workflow figure (2026-10-12)

The user requested session-start and a PDF explaining B-POP from the existing
`BPOP_WORKFLOW.md`. Created root `BPOP_WORKFLOW.pdf`, a single-page vector
schematic with isolated and dynamical C++ lanes, the retained-remnant feedback
loop, the four catalogue products, and the Python efficiency/rate/event stages.
It distinguishes simulation coverage from cosmic channel fractions and retains
the report's normalization, sparse-bin, detector-selection, and Pop III limits.
The report itself and all simulator/rate sources were preserved by this task.

Reproduction source is `figures/plot_bpop_workflow.py`; `figures/README.md`
documents its source and dependencies. It embeds all figure text and requires
only Python plus Matplotlib. It resolves output paths from its own location,
uses project-local `tmp/bpop-workflow-figure/` for cache/validation files, and
writes only the root PDF as a figure artifact. Reproduce with:

```bash
rtk proxy /home/chun-lung.chan/.conda/envs/popIII/bin/python figures/plot_bpop_workflow.py
rtk proxy gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage BPOP_WORKFLOW.pdf
```

Validation checked all 51 text elements for page bounds and overlaps and 17
nodes for text containment; generation also succeeded from an unrelated
project-local working directory. Ghostscript validated the PDF and extracted
its text. No raster previews, simulations, builds, or rate calculations were
performed for this figure task. The existing graph was queried for orientation;
current report/source checks supplied finer details. No graph rebuild occurred.
Existing working changes and the concurrently refreshed handoff were preserved.

## Architecture and important files

- Standard C++ simulator: `src/Main_Bpop_dev.cpp`; shared physics/sampling:
  `src/funcGR.cpp`, `include/funcGR.h`; configuration: `include/input_params.h`.
  Root Makefile delegates to `build/Makefile`, producing `build/BPOP.x`.
- Standard rates: `utils/BPOP_MERGER_RATE/BPOPRate_v2025.py` and
  `utils/BPOP_LIBPY/lib_BPOPRate.py`. Author-specific absolute paths remain;
  local Pop III repairs have not been propagated to production scripts.
- Experimental workflow: `tests/pop3_compile_test/` has its own C++ main,
  parameter header, Makefile, Python driver/library, and usage README.
- `BPOP_WORKFLOW.md`: current-source workflow, flowchart, equations, 19
  literature links, and catalogue/schema distinctions.
- `tests/pop3_compile_test/input_params_p3.h`: compact parameter comments
  explain units, effects, inactive settings, and Pop III overrides. All 67
  macro definitions and all non-comment content still match HEAD.
- `POP3_REFACTORING_PLAN.md`: proposed improvements, not implemented. Profile
  before changing RNG/bin-count/sampling behavior; preserve scientific models.
- `tests/pop3_compile_test/POP3_LABEL_REPORT.md`: older review; Python
  syntax/runtime findings predate the local repairs. Use current code/README.

## Working tree

HEAD and recorded `origin/main`: `447bbfd` (`added test for pop3`). Staging is
empty. No commit/push was made for the workflow report or parameter comments.
Preserve existing changes; do not reset or stage the working tree wholesale.

Pre-existing modifications: root `.gitignore`, root/build Makefiles, tracked
root build artifacts, local Python pair, test `.gitignore`, and test README.
This session added the workflow report and parameter comments; HANDOFF.md is
refreshed here. AGENTS.md, refactoring/label reports, diagnostic PDF and its
reproduction directory, and project-local `tmp/` also remain untracked.
The separate figure task generated untracked `BPOP_WORKFLOW.pdf` and `figures/`;
see its reproduction and validation details above.

Large local datasets are excluded from Git: `data_1e-11.dat`,
`spectrum1e-11.txt`, `spectrum_cleaned00000000001.txt` in the test directory.
The README directs readers to the maintainer for download links. Keep these,
generated catalogues, and environments out of source commits.

## Pop III simulator

- Local Makefile builds `BPOP.x` from the local main and shared `funcGR.cpp`;
  objects/dependencies live under `.build/`. Generated `build_paths.h` records
  checkout root/data directory. Header changes rebuild; compiler-flag changes
  or checkout relocation require a clean rebuild.
- Defaults: N=1000, CHUNKS="no", SEVN="pop3", uppergap="no", bhseed="no",
  spinlb="fuller", DynOvTot=0.9, pYC/pGC/pNC=0.75/0.15/0.1, mixing=0.5,
  fbin=0, CLfill="GG23". Pop III normally fixes Z=1e-11, bypassing metallicity
  scatter. See annotated header for explicit override conditions.
- Loaders accept the binary catalogue text header, retain complete correlated
  numeric rows, reject malformed/non-finite inputs, and warn/discard the mixed
  spectrum's incomplete final record.
- Use `make run`: it checks inputs, creates fresh `tmp/pop3-run.*/`, archives
  executable/configuration, and requires the completion marker. Output.log is
  renamed Output.txt after simulator file moves. Do not run directly in the
  dataset directory: the simulator moves working-directory text/data files.
  `clean` preserves datasets and completed runs.

## Python workflow and retained results

- Local driver explicitly imports `lib_BPOPRate_p3.py`. CLI selects catalogue
  directory, grid sizes, observer interval, seed, and rate-only execution.
  Each invocation creates `tmp/bpoprate-run.*/`, archives sources, and records
  input hashes/versions/assumptions. Planck18 recomputes redshift-grid ages.
- Repairs include modern SciPy Simpson calls, singleton interpolation,
  non-merger formation times from column 25, filtered efficiency indexing,
  zero-channel outputs, integer sample indices, preserved status strings,
  clamped negative smoothed rates, and preserved mono-metallicity.
- Required `--smoke-test` uses generic model `a` SFRs, fixed Z=1e-11 with
  probability one, zero isolated efficiency/fraction, and inherited dynamical
  normalizations. `models("pop3")` raises for missing calibration. Successful
  execution is not physical Pop III validation. Rate execution creates no plots.
- Retained simulator runs: `tmp/pop3-run.3vULll/` (382 main rows), `Rt5vbi/`
  (351 main rows), and `c44UBX/`. Prior validation checked finite main/chain
  catalogues, basic physical bounds, and natal retention CDFs. Records include
  `tmp/pop3_validation.json` and `tmp/pop3_first_run_hashes.json`.
- Full rate trial: `tmp/bpoprate-run.l2fkfo47/`, input from `3vULll`,
  500-point merger/formation grids, Tobs=0.01 yr, seed=12345. Four rate tables:
  499 rows, three columns, finite/nonnegative. Mock catalogues: 22 columns;
  nuclear=3, globular=2, young=11, isolated=0, global=16 events. Run.log,
  metadata, validation script, and `validation.json` remain there.
- Diagnostic PDF: `tests/pop3_compile_test/bpoprate_diagnostics.pdf`.
  Adjacent `bpoprate_diagnostics/` contains script, seven archived inputs,
  checksum manifest, README, and plot validation. Script resolves paths from
  its own location and writes PDF one level above. Figure labels the generic
  execution-trial assumptions. Obsolete scratch/runs were previously removed;
  use the retained paths above.

## Outstanding issues and future work

1. Calibrated Pop III SFR/channel normalizations and isolated efficiency are
   missing; clusters still use generic YC/GC/NC prescriptions. Supply these
   before interpreting physical rates.
2. Efficiency and the rate kernel's dN/N both contain non-merger factors.
   Review their combined scientific normalization.
3. The 1,000-candidate population is sparse; grid changes substantially affect
   rates. Sampling skips unsupported bins, so expected/written counts differ.
   Larger populations and convergence checks remain necessary.
4. Dissolved clusters can write infinite non-merger semimajor axes. The local
   rate workflow does not read those columns; physical treatment is unresolved.
5. Schema discrepancies: main `tlast_mer` is absolute cosmic age; subtract
   `tfor` for elapsed delay. Root chain writer emits 32 fields, schema says 33;
   Pop III adds BH count. Chain times are written in years where schema says
   Myr. Main catalogue summarizes terminal outcomes, not all intermediate events.
6. Kick-constrained sampling can loop indefinitely when no row is eligible;
   shared helpers repeatedly initialize RNG engines. See refactoring plan.

All requested documentation work is complete. These are future directions
requiring an appropriate task, not instructions to implement automatically.

## Graph, environment, and constraints

`graphify-out/` has 268 nodes, 479 edges, 11 communities, covering the original
19-file corpus. It excludes experimental sources/recent documents; verify
current source directly. Raw extraction flagged 155 dangling endpoints and
three collapsed relationships; exported endpoints resolve through external
nodes. Cost zeros mean unmeasured host-agent token usage. Preserve annotations.

Graphify 0.9.87 lives in `tmp/graphify-venv`; installed skill is 0.9.56.
`/home/chun-lung.chan/public_html/BPOP_p3_graphify` is a verified symlink to
local graph output, exposing artifact changes without automatic rebuilding.
The workflow review was saved in graph memory; no graph rebuild occurred.

Use `/home/chun-lung.chan/.conda/envs/popIII/bin/python` (Python 3.11) for
scientific execution. System `python3` is 3.6 and lacks newer Python APIs.
Conda supplies NumPy, SciPy, Astropy, pandas, Matplotlib, seaborn, and
scikit-learn; retained run metadata records exact versions.

Follow AGENTS.md and `/home/chun-lung.chan/.codex/RTK.md`: shell commands start
with `rtk`; use `rtk proxy` for unfiltered behavior. Temporary storage belongs
in project-local `tmp/`. Generate plots only for an explicit task requiring them.

## Validation at handoff

- `git diff --check` passes across tracked changes.
- Parameter header matches HEAD after comments/whitespace are removed.
- AST syntax checks pass for the local rate driver, library, and plot script.
- Workflow report local links and code/equation fences pass checks.
- All seven archived diagnostic input hashes match; PDF has a PDF signature.
- Retained runs/reports, graph counts, public symlink, and dataset ignore rules
  were checked. Runtime/plot results above are retained evidence, not reruns.
- No new build, simulation, rate calculation, or plot was needed for the
  documentation changes. Continuous-learning review found no new cross-project
  lesson requiring a skill-library update.

## Commands

From repository root:

```bash
rtk proxy make -C tests/pop3_compile_test
rtk proxy make -C tests/pop3_compile_test run
rtk proxy make -C tests/pop3_compile_test clean
rtk proxy tmp/graphify-venv/bin/graphify query "<question>"
rtk git status --short
```

Reproduce rates from the test directory (Bash), replacing the input path for
a new simulation:

```bash
rtk proxy bash -c 'catalogue_dir=(../../tmp/pop3-run.3vULll/SIM_*); /home/chun-lung.chan/.conda/envs/popIII/bin/python -u BPOPRate_v2025_p3.py --catalogue-dir "${catalogue_dir[0]}" --smoke-test --tobs 0.01 --seed 12345'
```

For a faster check, add `--grid-size 80 --formation-grid-size 80`; add
`--rates-only` to skip sampling. Changed grids change sparse-sample rates.
