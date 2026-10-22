# Repository Guidelines

## Project Structure & Module Organization

B-POP synthesizes binary black hole populations and merger catalogues. The C++ simulator lives in `src/Main_Bpop_dev.cpp`; physics and sampling helpers are in `src/funcGR.cpp` and `include/funcGR.h`. Compile-time settings reside in `include/input_params.h`, with cluster tables under `include/ClustersRhM/`. `build/` contains the build Makefile and executable. `tools/` provides chunk execution and catalogue collation scripts. Python postprocessing lives in `utils/BPOP_MERGER_RATE/` and `utils/BPOP_LIBPY/`. Read `README.md`, `input.md`, and `struct_*.md` for usage and output schemas. `tests/pop3_compile_test/` contains experimental Population III variants.

## Build, Test, and Development Commands

- `make`: build `build/BPOP.x` using C++17 and optimization.
- `make clean && make`: rebuild all objects after configuration or compiler changes.
- `make CXX=g++`: select the compiler explicitly.
- `make install PREFIX=/chosen/environment`: install the executable and parameter snapshot; defaults to the active conda prefix.
- `./BPOP.x > Output.txt`: run from a prepared simulation directory containing the copied executable and required inputs.
- `sh bpop_chunks.sh`: launch chunks from the prepared run directory; configure the script's `N` first. Run `sh bpop_catcatalogues.sh` after completion.
- `python BPOPRate_v2025.py`: run postprocessing from its prepared directory with simulator catalogues in the parent directory and `lib_BPOPRate` on the Python import path.

## Coding Style & Naming Conventions

Match surrounding C++ indentation, typically two spaces; use four spaces for Python and tabs for Makefile recipes. Preserve existing scientific identifiers and parameter macros. Give new helpers descriptive names and document units and physical assumptions. No formatter or linter configuration is present; avoid unrelated reformatting.

## Testing Guidelines

No automated test framework, test target, or coverage threshold is configured. Treat the Population III variants as experimental fixtures. For simulator changes, compile and run a small configured population with available input catalogues. Check finite values, catalogue columns, and physical bounds against `struct_*.md`; compare stochastic results statistically. Record settings, input datasets, and validation commands.

## Commit & Pull Request Guidelines

History uses short descriptive subjects such as `Update input.md` and `fixed bugs for timing and mass selection`; no formal prefix convention is established. Keep commits focused. PRs should explain the changed behavior, scientific impact, configuration, and validation results; link relevant issues and update schema documentation when outputs change.

## Configuration & Agent Instructions

Set `PREDIR`, `SINGPTH`, and `PATH` for local datasets before compiling. Archive the parameter header with each run. Keep generated catalogues and temporary environments out of source commits. Agents must prefix shell commands with `rtk`, using `rtk proxy` for unfiltered behavior, and store temporary files in project-local `tmp/`.
