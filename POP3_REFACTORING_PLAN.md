# Incremental improvements to the Population III workflow

## Summary

Prioritize **random-number generation in C++** and **repeated catalogue scans in Python**. These are the clearest performance opportunities in the current source. This is a static review; speedups still need measurement.

Keep the existing scientific prescriptions, smoke-test assumptions, catalogue formats, and CLI defaults.

## Recommended changes, in priority order

1. **Seed the C++ random engine once per run.**  
   In [funcGR.cpp](src/funcGR.cpp#L605), `rnd()` constructs a `random_device` and seeds a new `mt19937` for every draw. Gaussian helpers do the same. Pass one run-owned engine through the sampling helpers, including `singBHt_mix`. Add an optional seed and record the resolved value. This removes repeated initialization and enables reproducible simulation runs. Random sequences will change, so validate distributions statistically.

2. **Replace repeated Python catalogue scans with bin counts.**  
   The [rate calculation](tests/pop3_compile_test/BPOPRate_v2025_p3.py#L766) repeatedly scans merger and non-merger arrays inside the merger-redshift × formation-redshift loops. For each existing formation grid, assign catalogue rows to bins once, then calculate denominator and matching-merger counts with indexed accumulation. Preserve the current moving formation grids, half-open bin boundaries, cosmic-time selection, and Simpson integration. This removes an inner loop that repeatedly touches the whole catalogue.

3. **Remove repeated array copying during Python sampling.**  
   In [lib_BPOPRate_p3.py](tests/pop3_compile_test/lib_BPOPRate_p3.py#L614), each `np.concatenate` copies the accumulated selection again. Collect index arrays in a list and concatenate once, preserving the existing reverse accumulation order. Replace repeated catalogue-string accumulation with `join` or streamed writes. Evaluate comoving volumes at all redshift edges together and reuse their differences within a sampling run.

4. **Make kick-constrained sampling terminate reliably.**  
   Both single and mixed BH samplers repeatedly draw until `kick <= vescape`; they cannot finish when no eligible row exists. Build a sorted index of kicks once per spectrum, find the eligible prefix by binary search, and uniformly select an original row from that prefix. Return an explicit “no eligible BH” result and fail with context where the caller requires one. Preserve each row’s correlated mass, delay, and kick values.

5. **Extract small components from the Pop III driver.**  
   Start with spectrum loading, BH sampling, and catalogue writing. Introduce a `SpectrumCatalogue` containing the aligned vectors and a named `BHSample` result instead of `single_bh[0..2]`. Replace manually allocated arrays with vectors or fixed-size arrays, retaining the current data layout where practical. Avoid a full simulator rewrite.

6. **Reduce output overhead and isolate file handling.**  
   Replace per-row `endl` with `'\n'`; flush at completion or explicit checkpoints. Open output files directly inside the run directory and track their paths, replacing shell commands that move every `.txt` and `.dat` file. Report missing inputs and failed writes with a nonzero exit status.

## Validation

- Capture baseline wall time and peak memory separately for C++ simulation, Python rates, and Python sampling.
- Compare optimized rate-bin counts exactly against the existing loops, including boundary values, empty bins, and the single-metallicity case.
- Compare rate tables numerically before output rounding; retain identical grid settings and integration conventions.
- Check seeded repeatability, catalogue columns, physical bounds, and sampling distributions.
- Test empty spectra, no eligible kicks, truncated input, and failed output creation.
- Implement each improvement separately so its benefit and scientific impact can be assessed.

## Defaults and boundaries

Start with items **1–3**, then address sampler termination and structural cleanup. Shared C++ helper changes require checking the production driver too. Defer parallel execution, binary input caches, runtime replacement of all parameter macros, and changes to Pop III calibration until profiling or a separate scientific review justifies them.
