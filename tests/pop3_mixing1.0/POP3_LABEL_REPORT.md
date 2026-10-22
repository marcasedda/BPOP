# Why the `pop3` label is needed

Date: 2026-10-12. Scope: the C++ trial and Python source files in
`tests/pop3_compile_test/`, compared with the shared simulator sources.
Line references below describe the source at the time of this review.

## Main finding

`pop3` is a model-selection value. In the C++ trial, `#define SEVN "pop3"`
selects a one-metallicity input population, bypasses the usual metallicity
sampling, and accommodates the supplied Population III catalogue formats.
Removing these branches or replacing the value with `"yes"` would change
which data and sampling rules the simulator uses.

The exact word `pop3` is a software convention: another name could work if all
producers and consumers were changed together. **The distinction between this
population and the other supported populations is necessary in the current
implementation.** A `_p3` filename alone cannot select it.

The label currently supports a runnable C++ trial. It does not establish a
complete Population III cluster or cosmological merger-rate model: several
corresponding Python branches and the dedicated cluster prescription remain
unfinished.

## 1. Different kinds of Population III labels

| Label | Location | Meaning |
| --- | --- | --- |
| `SEVN "pop3"` | [input_params_p3.h](input_params_p3.h), line 15 | C++ model selector, incorporated into the executable at compilation. Despite its name, `SEVN` distinguishes several catalogue modes, including `"no"`, `"yes"`, `"all"`, `"iorio"`, and `"pop3"`. |
| `pre == "pop3"` | [BPOPRate_v2025_p3.py](BPOPRate_v2025_p3.py), lines 88, 211, 236 | Intended Python selector for the input population, metallicity list, and efficiency file. Its current default is `pre = "SEVN"`. |
| `a == "pop3"` | [lib_BPOPRate_p3.py](lib_BPOPRate_p3.py), lines 69–72 | Intended Python selector for star-formation and metallicity models. This is separate from `pre`; the driver currently sets `a = "a"`. |
| `_p3` suffix | Source/header/library filenames | Identifies the experimental variant and keeps its files separate. The Makefile explicitly compiles the local main, and the main explicitly includes the local header. The suffix has no automatic physical meaning. |
| `MCL_popIII`, `RCL_popIII`, `MBH_popIII`, `aBH_popIII` | Header lines 110–114; main lines 249–254 | Parameters for optional cluster/black-hole seed prescriptions. These names also occur in the shared main; they do not independently enable Population III mode. |

There is no command-line `--pop3` option in this trial. Changing the C++ header
requires recompilation; the local Makefile tracks that dependency
([Makefile](Makefile), lines 28–29).

## 2. What the C++ label actually controls

All source references in this section are to
[Main_Bpop_dev_p3.cpp](Main_Bpop_dev_p3.cpp).

| Source lines | Behavior selected by `SEVN == "pop3"` | Why it matters |
| --- | --- | --- |
| 423–435 | Sets `nmetal = 1`, also used as `numZ`. | The supplied trial has one metallicity dataset rather than the ordinary multi-metallicity grid. Array sizes and catalogue loops must agree with that coverage. |
| 504–506 | Sets `met[0] = 1.E-11`. | Identifies the input metallicity and constructs the binary catalogue name `data_1e-11.dat`. |
| 581–583 | Sets `metdyn[0] = 1.E-11`. | Initializes the metallicity used to name the dynamical retention table. Lines 1633–1639 use it to write `retention_Z1e-11.dat`; omitting this initialization would leave the naming value indeterminate. |
| 864–865 | Uses `log10(met[0])` for isolated binaries and young-cluster sources under the current formation-time settings. | Avoids sampling an enriched-population metallicity distribution for the one-metallicity input population. |
| 912–913 | Applies the same fixed metallicity to globular/nuclear cluster sources under the current settings. | Keeps the other dynamical channels consistent with the supplied grid. |
| 1139–1150 | Recognizes and skips the supplied `BH_tot ...` binary-catalogue header. | The following reader expects 20 numeric values per record; the text header cannot be read as such a record. |
| 1459–1460 | Uses the full metallicity string for the mixed spectrum. | Opens the supplied `spectrum1e-11.txt`. |
| 1498–1499 | Uses the special cleaned-spectrum suffix `00000000001`. | Opens the supplied `spectrum_cleaned00000000001.txt`, whose naming convention differs from the mixed spectrum. |
| 1932–1937 | Contains only comments about a future Population III cluster mass/radius prescription. | Marks an intended extension point; it currently changes no cluster property. |

### Why a tiny nonzero metallicity appears

This implementation represents its Population III input population with
`Z = 1e-11`. It repeatedly takes `log10(Z)` and uses that value in filenames.
Using literal zero would give a non-finite logarithm and would request different
filenames. The code therefore requires the supplied nonzero value. Interpreting
it as a numerical proxy for a nearly pristine population is consistent with
the Population III comments, but the source does not document the stellar-track
provenance or validate that scientific approximation.

### Why filename handling needs special treatment

The ordinary spectrum-name construction removes the first two characters of
the streamed metallicity string (lines 1419–1429), matching values written as
`0.0002` and similar decimals. For `1e-11`, this operation produces `-11`.
Without the Pop III filename branches, the code would request
`spectrum-11.txt` and `spectrum_cleaned-11.txt`, neither of which is the supplied
file. The cleaned-spectrum suffix is currently hardcoded for this dataset; a
different Pop III input metallicity would require revisiting it.

### Why `MonoZ` alone is insufficient

`MonoZ` changes assigned source metallicities (lines 968–973). It does not
select a one-entry catalogue grid, initialize its entries, skip the binary
header, or fix the spectrum names. With an ordinary catalogue mode, later
source-to-grid assignment maps metallicities to that mode's available inputs
(for example, lines 1067–1072 map values below the minimum to `met[0]`). Setting
`mono_Z = 1e-11` therefore cannot substitute for selecting the Pop III inputs.

The fixed-metallicity branches are also conditional: `bigbang`, `single`, and
`grid` formation settings are checked earlier, and the later `MonoZ`/seed
override still applies. The present header uses `continuous`, `MonoZ "no"`,
and `TagR "AS20"`, so those overrides do not apply to the current trial.

## 3. What the label does not automatically select

The local Makefile links the main with shared
[../../src/funcGR.cpp](../../src/funcGR.cpp). The shared physics implementation
has no `pop3` selector. Most merger, spin, kick, and cluster-evolution routines
are reused; the supplied input distributions and separate configuration
choices determine how they are applied.

The current header separately sets `spinlb "fuller"`, `uppergap "no"`,
`bhseed "no"`, `mixing 0.5`, and `fbin 0.0`. Comments recommend some of these
settings for the Pop III project, but `SEVN "pop3"` does not enforce them.
The seed constants are not an automatic Population III seed population:
their cluster values are used by `TagR == "seed"` (main lines 1876–1878), and
their black-hole values by the matching seed branch (lines 2344–2348).

The empty block at lines 1932–1937 explicitly proposes reusing existing young
and globular cluster models for testing. The current configuration also assigns
a nuclear-cluster fraction (`pNC = 0.1`). Consequently, the label alone does
not imply that the adopted environments, masses, radii, or formation times
have been calibrated for Population III clusters. Formation times still use
the existing `func.sfr_red(...)` path (lines 820–836).

The local main also differs from the shared main in several timing and
hierarchical-merger calculations. Those changes are not guarded by `pop3`.
A comparison between the two executables would therefore change more than the
population label and should not be interpreted as a controlled test of the
label's scientific effect.

## 4. Why Python needs corresponding labels—and why it is unfinished

The rate calculation must select inputs and population weights consistent
with the C++ catalogues. The draft makes that distinction in two independent
places:

- `pre == "pop3"` selects `Z = [1.E-11]`, proposes `IBcorr = 0.185`, and sets
  `app = "_pre"` for the efficiency-file lookup. The resulting filename is
  `eta_isolated_binaries_pre` when `aCE` is empty. A comment says the efficiency
  file is still to be supplied (driver lines 239–240); this file is absent from
  the trial directory.
- `models("pop3")` returns `"pop3"` for all four star-formation and metallicity
  model entries, ordered nuclear, globular, young, isolated. Its normalization
  parameters remain zero, with a comment that the dynamical normalization
  needs discussion (library lines 69–81).

Neither selector is currently active in the driver defaults, and the driver
imports `lib_BPOPRate` from a hardcoded external directory rather than the local
`lib_BPOPRate_p3` (lines 22–23). Even after correcting those connections, simply
setting both labels would leave these concrete problems:

| Problem | Evidence and consequence |
| --- | --- |
| Driver cannot parse. | Line 88 is `elif(pre=="pop3":`, missing a closing parenthesis. Other top-level statements also have inconsistent indentation. |
| No Pop III star-formation function. | `psi()` supports `mf17` and `kr13`, then falls back to `psiNC` for every other name (library lines 138–150). `"pop3"` therefore receives the generic Gaussian rather than a dedicated prescription. |
| No Pop III metallicity weighting. | `metcor()` has no `"pop3"` case and returns an unassigned `logz_me` for that input (library lines 153–162). The driver calls this function before its Gaussian metallicity-bin calculation (lines 773–778). |
| Remaining generic rate normalization. | The driver still sets total star formation using `psi("mf17", ...)` (line 285). Its own comments identify Pop III isolated/young star formation as pending work (lines 280–282). |

Thus, the Python labels express necessary distinctions, but their downstream
implementations are incomplete. A one-metallicity Pop III rate model needs an
explicit weighting rule for that population and defined star-formation and
efficiency normalizations; selecting a string does not provide those inputs.

## 5. Output provenance

The C++ catalogue header records `Metal` and `EnvType`, but it has no explicit
Population III column (main line 1050). The `SIM_*` directory-name construction
does not include `SEVN` (lines 3365–3450 and following). Therefore, the source
label currently selects behavior without being preserved as an explicit model
label in those output names.

The run wrapper copies `input_params_p3.h` and the generated path header into
each run directory ([Makefile](Makefile), line 47). Preserve that snapshot with
the catalogues to establish that `SEVN "pop3"` was used. If outputs are later
combined across populations, explicit population metadata would help preserve
their origin; adding it would be a separate implementation task.

## 6. Evidence and validation limits

This review used direct source searches, a comparison of the local and shared
mains, input-file header inspection, and read-only checks of previous trial
outputs. The existing graphify graph was queried for context, but its corpus
does not include this experimental directory; the findings above come from the
current files.

Checks performed during this review:

- Verified that all three input files listed by the local Makefile exist, and
  inspected their first lines. The binary file has the expected text header;
  both spectra begin with numeric five-value records.
- Parsed both Python files with `ast.parse` without importing them. The driver
  fails at line 88; the local library parses successfully.
- Executed only the selected library function definitions, extracted using
  Python's AST, with standard-library substitutes for the required numerical
  operations. `models("pop3")` returned the eight intended tags and zero
  normalization parameters; `psi("pop3", 1, 1, 1, 1, 1)` returned `1.0` through
  its Gaussian fallback; `metcor("pop3", 10, 0.2)` raised `UnboundLocalError`.
  These are narrow logic checks, not a full Python pipeline run.
- Inspected existing runs `tmp/pop3-run.Rt5vbi/` and `tmp/pop3-run.3vULll/`.
  Their main catalogues contain 351 and 382 rows, respectively, with the sole
  metallicity `1e-11`. Both contain `retention_Z1e-11.dat` and a completion
  marker in their logs. These are prior runs, not simulations started by this
  review.

No simulator source, configuration, input data, or executable was changed.
No new build, simulation, rate calculation, or plot was needed for this report.
The evidence establishes how the labels operate and identifies unfinished
branches; it does not validate the Population III physics.

## Conclusion

Keep the `pop3` selection in the current C++ trial because it binds the
population to its actual metallicity coverage and input conventions. Keep the
corresponding distinction in future Python rate work, while implementing and
validating the missing prescriptions. File suffixes and comments communicate
intent; the executable branches and compatible scientific inputs determine
the resulting population.
