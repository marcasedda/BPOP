# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

B-POP is a semi-analytic population synthesizer for binary black hole (BBH) mergers in isolated binaries and in young, globular and nuclear star clusters (YC/GC/NC). This checkout also has an experimental Population III (Pop III) variant. See `AGENTS.md` for repository guidelines. `BPOP_WORKFLOW.md` explains the physics and data flow, with source line references. `HANDOFF.md` records the latest session state, the retained runs and the open issues.

## Mandatory conventions

- Prefix every shell command with `rtk`. Use `rtk proxy <cmd>` when you need unfiltered output.
- Put temporary files and run directories in the project-local `tmp/`.
- Run scientific Python with `/home/chun-lung.chan/.conda/envs/popIII/bin/python` (3.11, conda env `popIII`). The system `python3` is 3.6 and too old.
- The working tree has uncommitted changes from earlier sessions. Keep them: never reset the tree or stage everything at once.
- Do not commit generated catalogues, large datasets (`tests/pop3_compile_test/*.dat`, `spectrum*.txt`), environments or `tmp/`.
- Generate plots only when a task explicitly asks for them.

## Build and run

Standard simulator (from the repo root):

```bash
rtk proxy make                # builds build/BPOP.x (root Makefile delegates to build/Makefile)
rtk proxy make clean && rtk proxy make
rtk proxy make CXX=g++
```

All parameters are compile-time macros in `include/input_params.h`, so any parameter change requires a rebuild. `PREDIR`, `PATH` and `SINGPTH` hold absolute paths to the upstream author's SEVN/MOBSE stellar catalogues. Set them to local datasets before running. Copy `BPOP.x` into a dedicated run directory and run `./BPOP.x > Output.txt` there. Parallel runs use the scripts in `tools/`: copy them into the run directory, set `N` in `bpop_chunks.sh`, run `bpop_chunks.sh`, then run `bpop_catcatalogues.sh`. Each chunk writes `chunk.out` and needs `CHUNKS "yes"`.

Pop III variant (self-contained under `tests/pop3_compile_test/`):

```bash
rtk proxy make -C tests/pop3_compile_test            # build
rtk proxy make -C tests/pop3_compile_test run        # checks inputs, runs in fresh tmp/pop3-run.XXXXXX/
rtk proxy make -C tests/pop3_compile_test clean      # keeps datasets and completed runs
```

`make run` needs three large datasets that are not in Git (`data_1e-11.dat`, `spectrum1e-11.txt`, `spectrum_cleaned00000000001.txt`). The test Makefile generates `.build/build_paths.h` from the checkout location. Rebuild from clean after moving the checkout or changing compiler flags.

Pop III rate trial (run from `tests/pop3_compile_test/`):

```bash
rtk proxy bash -c 'd=(../../tmp/pop3-run.3vULll/SIM_*); /home/chun-lung.chan/.conda/envs/popIII/bin/python -u BPOPRate_v2025_p3.py --catalogue-dir "${d[0]}" --smoke-test --tobs 0.01 --seed 12345'
```

For a quick check, add `--grid-size 80 --formation-grid-size 80`. Add `--rates-only` to skip sampling. Each invocation writes to a fresh `tmp/bpoprate-run.*/` with archived sources and `run_metadata.json`.

There is no test framework or linter. To validate a change, build, run a small population and check the catalogues: values must be finite, the column counts and physical bounds must match `struct_*.md`, and stochastic outputs should be compared statistically. Report the settings and commands you used.

## Architecture

The pipeline has two stages.

1. **C++ simulator** (`src/Main_Bpop_dev.cpp` holds a single ~3400-line `main`; physics and sampling helpers are in `src/funcGR.cpp` / `include/funcGR.h`):
   - It allocates the `N` candidates deterministically among the isolated, YC, GC and NC channels (`DynOvTot`, `pYC`, `pGC`, `pNC`; NC gets the remainder).
   - It samples formation redshift (`SFRTYPE_*`) and metallicity, then processes candidates grouped by metallicity.
   - Isolated binaries draw correlated masses and delays from precomputed SEVN/MOBSE catalogues.
   - Dynamical binaries draw a host cluster from the tables in `include/ClustersRhM/`. The code then advances through sinking, pairing, hardening and GW inspiral. It computes the remnant mass, spin and recoil, and if the remnant is retained it loops into hierarchical mergers.
   - Outputs are `Catalogue.txt`, `Catalogue_multiple_dyn.txt`, `Larger_than_tH.txt` (non-mergers) and `retention_Z*.dat`.
   - At exit it runs `system("mv *.txt *.dat <SIM_dir>")`, which moves **every** `.txt` and `.dat` file in the working directory. Never run the binary in a directory that holds inputs or sources.
2. **Python rate stage** (`utils/BPOP_MERGER_RATE/BPOPRate_v2025.py` driver plus `utils/BPOP_LIBPY/lib_BPOPRate.py` library):
   - It computes metallicity-dependent merger efficiencies and convolves star-formation histories with the delay-time distributions. The results are merger rate densities per environment.
   - It converts rates to observer counts and resamples matching catalogue rows into mock event catalogues.
   - The production driver is configured by editing module-level variables. It reads catalogues from `../` and imports the library through a hard-coded `sys.path` entry that points at the upstream author's machine.

Pop III variant: `tests/pop3_compile_test/` has its own `Main_Bpop_dev_p3.cpp`, `input_params_p3.h`, Makefile and Python pair (`*_p3.py`). It links the shared `src/funcGR.cpp`, so changing shared helpers affects both drivers. The Pop III Python pair has local repairs and a CLI that the production scripts under `utils/` do not have.

## Known caveats (see HANDOFF.md for details)

- Pop III SFR, channel normalisations and isolated efficiency are uncalibrated. `--smoke-test` output checks that the code runs; it is not a physical prediction.
- Schema docs disagree with the code in places. The main catalogue's `tlast_mer` is absolute cosmic age. The chain writer emits 32 fields while the schema lists 33. Chain times are in years where the schema says Myr.
- Kick-constrained BH sampling can loop forever when no row is eligible. The RNG is re-seeded on every draw in `funcGR.cpp`. `POP3_REFACTORING_PLAN.md` lists proposed improvements that are not yet implemented.
- Dissolved clusters can write infinite semimajor axes to `Larger_than_tH.txt`.
