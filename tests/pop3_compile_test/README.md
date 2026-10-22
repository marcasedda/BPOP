# Population III compile test

Build and run the local C++ Pop III variant using the shared physics helpers
and tables from this repository. Requires a full checkout, Make, and a C++17
compiler such as `g++`.

## Download the inputs

The following datasets are excluded from Git because of their size.
**Contact the repository maintainer to request download links**, then place
all three files in `tests/pop3_compile_test/` with these exact names:

| File | Approximate size |
| --- | ---: |
| `data_1e-11.dat` | 13 MiB |
| `spectrum1e-11.txt` | 135 MiB |
| `spectrum_cleaned00000000001.txt` | 27 MiB |

Compilation works without these files; `make run` requires them and checks
that they exist before starting.

## Build and run

From the repository root:

```bash
cd tests/pop3_compile_test
make clean && make
make run
```

`make` builds the local `BPOP.x`. Use `make run` for execution: it creates a
fresh output directory, protecting the input files from the simulator's bulk
file moves. Avoid running `./BPOP.x` directly in this input directory.

## Configuration

Edit `input_params_p3.h`, then run `make` to recompile. Defaults are 1,000
sources at metallicity `1e-11`, Fuller spins, and disabled chunking, upper
mass-gap sampling, and BH seeds. Other model settings retain their supplied
values.

Paths are generated automatically from the checkout location. After moving
the checkout or changing compiler flags, use `make clean && make`. To select
a compiler explicitly, use `make CXX=g++`.

## Outputs

Each run prints its directory: `../../tmp/pop3-run.XXXXXX/`.

- `Output.txt`: execution log.
- `SIM_*/`: merger catalogues and retention table.
- `BPOP.x` and both configuration headers: snapshots used for the run.

`make clean` removes local build artifacts while preserving datasets and runs.

## Known limitations

- The supplied mixed spectrum has an incomplete final row. The loader warns
  and discards it; other malformed or non-finite input stops execution.
- Dissolved clusters can produce infinite semimajor axes in
  `Larger_than_tH.txt`. Check results before scientific use.
- Pop III rate calibration remains unfinished: no Pop III star-formation law
  or isolated-binary efficiency table was supplied. The Python workflow supports
  an explicitly labelled execution trial described below. `make run` executes
  only C++.

## Run the Python rate trial

Use Python 3.9 or newer with NumPy, SciPy (including
`make_smoothing_spline`), Astropy, pandas, Matplotlib, seaborn, and scikit-learn.
The local `popIII` conda environment has these dependencies. From this directory:

```bash
conda activate popIII
catalogue_dir=(../../tmp/pop3-run.3vULll/SIM_*)
python BPOPRate_v2025_p3.py --catalogue-dir "${catalogue_dir[0]}" \
    --smoke-test --tobs 0.01
```

Replace the example simulator directory with the directory printed by your
`make run`. The shell example uses Bash. The input directory must contain
`Catalogue.txt`, `Larger_than_tH.txt`, and `retention_Z1e-11.dat`.

Each invocation creates a fresh `../../tmp/bpoprate-run.*/`. Inputs are read
in place. `run_metadata.json` records the input paths, SHA-256 hashes, grid
settings, seed, package versions, and trial assumptions. Copies of both Python
source files are archived there. No plots are generated.

The default grids have 500 points each over redshift 0–15. For a quicker trial,
add `--grid-size 80 --formation-grid-size 80`. `--rates-only` skips catalogue
sampling; `--seed` controls its random draws; `--tobs` is the observer time in
years and affects sampled counts, not merger rate density. The short observation
time above limits trial catalogue size; the code's default remains 10 years.

Full runs place efficiency tables, four environment rate tables, and five
sample catalogues under `SET*/SMPCAT_Massive_no/`. Rate tables have redshift,
total rate, and the single-metallicity contribution; rates have units
`yr^-1 Gpc^-3`. Sample catalogues have the 22 columns described in the root
README. Rate-only runs retain their tables directly in the new output directory.

`--smoke-test` explicitly chooses the existing generic model `a` star-formation
histories, fixes the population to `Z=1e-11` with probability one, and disables
the isolated channel because its Pop III efficiency is unavailable. Other
normalisations (`fbh=0.0026`, `mstar=1 Msun`, binary fractions, and channel
scalings) are inherited. These results are execution checks, not calibrated
Pop III predictions. Running without this flag reports the missing model.

## How BPOPRate works

The completed trial's three-panel diagnostic figure is
[bpoprate_diagnostics.pdf](bpoprate_diagnostics.pdf). Its script and archived
plot inputs are in [bpoprate_diagnostics/](bpoprate_diagnostics/README.md).

`BPOPRate_v2025_p3.py` orchestrates the calculation; `lib_BPOPRate_p3.py`
provides the star-formation models, interpolation, cosmological factors, and
sampler. The shared production scripts under `utils/` are unchanged.

1. `models`, `msca`, and `psi` define star-formation histories and their scaling
   for nuclear, globular, young clusters, and isolated binaries. Planck18 ages
   convert between cosmic time and redshift. The checked-in redshift table sets
   the interpolation grid; its ages are recalculated with Planck18.
2. The script reads merger metallicities, environments, cluster properties,
   formation redshifts, and absolute merger times from `Catalogue.txt`. It reads
   non-merger formation times from column 25 of `Larger_than_tH.txt`. Infinite
   semimajor axes in that file are not read by this workflow.
3. Dynamical efficiency per metallicity is
   `eta = (fbh/mstar) * Nmerge/(Nmerge+Nnonmerge) * fbin * retention_weight`.
   The retention weight is a cluster-mass-weighted combination of the single
   and binary natal-kick retention curves, using `fmix=0.5`. With the inherited
   `repetita="no"`, the extra previous-merger count is excluded. The isolated
   channel normally needs an external metallicity-dependent efficiency table.
4. For each merger-redshift bin, the script integrates over earlier formation
   redshifts. Its integrand is
   `SFR * eta * (dN/N) * P(Z|z_form) * |dt/dz_form| / delta_t`, multiplied by
   `1e9` to convert volume units from Mpc cubed to Gpc cubed. `dN` counts
   candidates that form in a given formation bin and merge in the target cosmic
   time bin; `N` includes mergers and non-mergers formed in that bin. The trial
   fixes `P(Z|z_form)=1` for the supplied mono-metallicity population.
5. `bsmplr` smooths the rate tables and `Nsam` computes expected observer counts
   as `rate * delta_comoving_volume * Tobs/(1+z)`. It draws matching simulated
   mergers with replacement in each metallicity, environment, and redshift bin,
   then writes the per-environment and combined catalogues. A fixed seed makes
   sampling reproducible for identical inputs, grids, and package versions.

The 1,000-source trial is sparse: grid resolution changes its rate estimates,
and sampling skips bins with no matching simulated merger. Spline smoothing
also redistributes rates across bins. Expected and written source counts can
differ. The inherited efficiency and `dN/N` both involve non-merger counts;
their combined normalisation needs scientific review before physical inference.
Production use requires a Pop III SFR, channel normalisations, an isolated
efficiency table, and larger populations with resolution/convergence checks.
