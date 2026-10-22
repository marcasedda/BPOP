# Pop III mixing diagnostic

This compares the existing one-million-candidate runs at mixing 0, 0.5, and 1. Only `mixing` differs among their 66 parsed parameter macros (the header guard is excluded). Their local C++ main sources and all three stellar input datasets match byte for byte. The result is a comparison of simulator populations, not calibrated cosmic rates.

[Individual figures and brief explanation](panels/README.md) · [Figures (9-page PDF)](diagnostics.pdf) · [Environment summaries](summary.csv) · [All settings](settings.csv) · [Validation](validation.json)

## Main results

| Mixing | Terminal mergers | Terminal yield | Recorded chain events | Chains with records | Mean recorded events/chain |
| --- | --- | --- | --- | --- | --- |
| 0 | 396,658 | 39.6658% | 24,056,741 | 98,750 | 243.61 |
| 0.5 | 363,491 | 36.3491% | 13,704,958 | 97,392 | 140.72 |
| 1 | 373,946 | 37.3946% | 4,244,100 | 93,862 | 45.22 |

- Mixing 0.5 changes terminal counts by -8.36% and recorded chain events by -43.03% relative to mixing 0.
- Mixing 1 changes terminal counts by -5.73% and recorded chain events by -82.36% relative to mixing 0.

These two counts measure different things. The terminal catalogue contains one written endpoint per successful candidate; a recorded chain can contribute many intermediate mergers. Long chains dominate the event-weighted distributions. Per-chain figures give every recorded chain equal weight.

## Input setup and mechanism

| Parameter | mixing 0 | mixing 0.5 | mixing 1 |
| --- | --- | --- | --- |
| SEVN | "pop3" | "pop3" | "pop3" |
| N | 1000000 | 1000000 | 1000000 |
| DynOvTot | 0.9 | 0.9 | 0.9 |
| pYC | 0.75 | 0.75 | 0.75 |
| pGC | 0.15 | 0.15 | 0.15 |
| pNC | 0.1 | 0.1 | 0.1 |
| uppergap | "no" | "no" | "no" |
| bhseed | "no" | "no" | "no" |
| SFRTYPE_ISO | "continuous" | "continuous" | "continuous" |
| SFRTYPE_CLU | "continuous" | "continuous" | "continuous" |
| mixing | 0.0 | 0.5 | 1.0 |
| fbin | 0.0 | 0.0 | 0.0 |
| TagR | "AS20" | "AS20" | "AS20" |
| CLfill | "GG23" | "GG23" | "GG23" |
| spinlb | "fuller" | "fuller" | "fuller" |

`mixing` is the probability of selecting the mixed/binary-origin reservoir, not the isolated-binary fraction and not `fbin`. A uniform draw selects the reservoir for a dynamical candidate, and its companion draws use that reservoir. Both samplers reject BHs whose natal kick exceeds the cluster escape speed. Consequently, the unconditioned 0.5 mixture shown on PDF page 1 need not equal the distribution of successful mergers. The spectrum arrays preserve the input correlations between mass, formation time, and kick.

| Reservoir | Complete records | Mean BH mass [Msun] | Max BH mass [Msun] | Mean kick [km/s] | Discarded EOF values |
| --- | --- | --- | --- | --- | --- |
| single | 551,680 | 166.267 | 599.474 | 6.033 | 0 |
| mixed | 3,626,021 | 46.663 | 330.796 | 4.733 | 2 |

The EOF handling matches the simulator: only an incomplete final group of five numbers is discarded. These are input-spectrum summaries, before retention selection. A smaller reservoir mass scale can change recoil, pairing times, retention, and the supply of companions; the runs do not isolate those mechanisms individually.

## Environment-specific terminal yields

Denominators come from `Creating N` in each run log: isolated 99,999; young 675,000; globular 135,000; nuclear 90,001. This preserves the simulator's floating-point truncation and residual allocation rather than replacing them with rounded nominal fractions.

| Environment | mixing 0 | mixing 0.5 | mixing 1 |
| --- | --- | --- | --- |
| none | 97,422 (97.42%) | 97,387 (97.39%) | 97,352 (97.35%) |
| young | 100,769 (14.93%) | 77,884 (11.54%) | 78,155 (11.58%) |
| globular | 115,777 (85.76%) | 108,182 (80.13%) | 115,083 (85.25%) |
| nuclear | 82,690 (91.88%) | 80,038 (88.93%) | 83,356 (92.62%) |

The isolated channel does not use the mixing reservoir. Its differences across these runs provide a stochastic reference, not a paired control: random seeds are unrecorded and candidate IDs are not paired between runs.

| Environment | Recorded events: mixing 0 | mixing 0.5 | mixing 1 |
| --- | --- | --- | --- |
| young | 12,559 | 10,921 | 11,846 |
| globular | 102,069 | 90,550 | 75,453 |
| nuclear | 23,942,113 | 13,603,487 | 4,156,801 |

Nuclear clusters dominate the recorded event total in all three runs. Their event count falls strongly with increasing mixing, while their terminal yield remains near 89–93%. The main change is therefore in how far a small set of nuclear chains continues recycling, rather than in whether most nuclear candidates produce a terminal merger.

## Masses, delays, and growth

| Mixing | Terminal median Mrem [Msun] | Median elapsed delay [Gyr] | Mrem > 100 | Mrem > 1,000 | Mrem > 10,000 |
| --- | --- | --- | --- | --- | --- |
| 0 | 313.928 | 0.0235 | 260,598 | 14,739 | 249 |
| 0.5 | 104.848 | 0.0438 | 192,153 | 7,511 | 201 |
| 1 | 85.218 | 0.0646 | 140,724 | 154 | 146 |

| Mixing | Median chain length | 99th percentile length | Maximum length | Largest recorded Mrem [Msun] | Chains reaching > 10,000 Msun |
| --- | --- | --- | --- | --- | --- |
| 0 | 1 | 9.0 | 319,717 | 5.00006e+07 | 249 |
| 0.5 | 1 | 7.0 | 320,470 | 5.00005e+07 | 201 |
| 1 | 1 | 2.0 | 153,739 | 6.7568e+06 | 146 |

Lengths count recorded rows, not the `Nrec` counter. Maximum masses and growth factors use every recorded event; no display subsampling is used. Source column 17 is retained as the recorded previous-merger counter and is not assumed to equal a chain's row count. The fixed-seed macros are inactive with the current AS20 / no-seed setup, so massive remnants here do not imply that the separate prescribed seed experiment was enabled.

| Mixing | Top 1% of recorded chains (rounded up) | Share of all recorded events | Chains exceeding 5e7 Msun |
| --- | --- | --- | --- |
| 0 | 988 | 99.42% | 7 |
| 0.5 | 974 | 99.15% | 5 |
| 1 | 939 | 97.75% | 0 |

The event-weighted comparisons are dominated by a very small number of long chains. The largest masses for mixing 0 and 0.5 reach the configured `msmbhmax=5e7 Msun` stopping scale: the check acts on primary mass after an event is written, so the last remnant may slightly exceed it. Those maxima are affected by the simulation's imposed growth limit and are not unconstrained physical mass predictions.

## Catalogue semantics and reconciliation

The column mapping is checked against the local C++ writers. Main column 15 and chain column 16 are absolute cosmic ages in years; elapsed delay subtracts formation age. Chain columns 9–15 are also written in years, despite the older schema's Myr labels. The chain file has 33 columns in this local variant.

Chain column 23 identifies the candidate chain. Its column 31 plus one maps to the terminal catalogue ID because `itot` is incremented before the final writer. A GW-recoil exit occurs before the chain writer; a terminal merger may therefore be absent from the chain file, or differ from its last recorded event. Other termination paths also prevent simple equality. The files must not be added together as disjoint event lists.

| Mixing | Chains mapped to terminal | Last event matches terminal | Last event differs | Dynamical terminal without chain | Chain without terminal |
| --- | --- | --- | --- | --- | --- |
| 0 | 98750 | 30560 | 68190 | 200486 | 0 |
| 0.5 | 97392 | 42214 | 55178 | 168712 | 0 |
| 1 | 93862 | 44725 | 49137 | 182732 | 0 |

The match compares both component masses, remnant mass, and cosmic age with relative tolerance 2e-5 (catalogues print about six significant digits). Non-merger record IDs have writer-specific offsets and can overlap terminal IDs; their raw record counts are reported, but they are not used as complementary denominators. `chain_events_before_present` in the CSV applies a strict comparison to the printed cosmic ages and the simulator's 13.78689 Gyr present age. A few records round across that boundary at six significant digits; they remain included in the recorded-event diagnostics and are not classified as physical future mergers.

## Data-quality audit

### Mixing 0

- Full raw chain row count: 24,056,741; streamed and independent byte-level newline counts agree.
- Within-chain time/generation ordering violations: 0.
- Terminal numeric non-finite counts (zero-based columns): `{'30': 19}`.
- Chain numeric non-finite counts: `{}`.
- Non-merger numeric non-finite counts: `{'20': 481428, '21': 481409}`.
- Terminal physical-bound flags: `{'nonpositive_components': 0, 'remnant_exceeds_total': 0, 'spin_outside_bounds': 0, 'eccentricity_outside_bounds': 0, 'negative_elapsed_delay': 0}`.
- Chain physical-bound flags: `{'nonpositive_components': 0, 'remnant_exceeds_total': 0, 'spin_outside_bounds': 0, 'eccentricity_outside_bounds': 0, 'negative_elapsed_delay': 0}`.
- Non-merger duplicate IDs: 0; terminal/non-merger ID overlap: 108707.

### Mixing 0.5

- Full raw chain row count: 13,704,958; streamed and independent byte-level newline counts agree.
- Within-chain time/generation ordering violations: 0.
- Terminal numeric non-finite counts (zero-based columns): `{'30': 10}`.
- Chain numeric non-finite counts: `{}`.
- Non-merger numeric non-finite counts: `{'20': 492202, '21': 492192}`.
- Terminal physical-bound flags: `{'nonpositive_components': 0, 'remnant_exceeds_total': 0, 'spin_outside_bounds': 0, 'eccentricity_outside_bounds': 0, 'negative_elapsed_delay': 0}`.
- Chain physical-bound flags: `{'nonpositive_components': 0, 'remnant_exceeds_total': 0, 'spin_outside_bounds': 0, 'eccentricity_outside_bounds': 0, 'negative_elapsed_delay': 0}`.
- Non-merger duplicate IDs: 0; terminal/non-merger ID overlap: 97726.

### Mixing 1

- Full raw chain row count: 4,244,100; streamed and independent byte-level newline counts agree.
- Within-chain time/generation ordering violations: 0.
- Terminal numeric non-finite counts (zero-based columns): `{'30': 9}`.
- Chain numeric non-finite counts: `{}`.
- Non-merger numeric non-finite counts: `{'20': 503442, '21': 503433}`.
- Terminal physical-bound flags: `{'nonpositive_components': 0, 'remnant_exceeds_total': 0, 'spin_outside_bounds': 0, 'eccentricity_outside_bounds': 0, 'negative_elapsed_delay': 0}`.
- Chain physical-bound flags: `{'nonpositive_components': 0, 'remnant_exceeds_total': 0, 'spin_outside_bounds': 0, 'eccentricity_outside_bounds': 0, 'negative_elapsed_delay': 0}`.
- Non-merger duplicate IDs: 0; terminal/non-merger ID overlap: 92607.

The non-finite terminal entries are in semimajor axis (one-based column 31): 19, 10, and 9 entries across mixing 0, 0.5, and 1. Non-merger separation columns 21–22 contain many non-finite values. These columns are not used for the displayed mass, spin, or delay distributions. All recorded chain numeric columns are finite, and the checked mass/spin/eccentricity/time bounds have no violations. Spin-bound checks here cover primary spin, remnant spin, and effective spin; secondary spin is not separately aggregated.

Non-finite entries are excluded only from the affected quantity, not silently from whole records. Physical-bound flags are reported rather than repaired. Logarithmic panels omit nonpositive values. Every quantity's finite count, non-finite count, zeros, negative values, and histogram underflow/overflow are recorded in `statistics.csv`. Plotted densities divide by the full finite count and bin width; a panel with omitted mass can therefore integrate to less than one. Terminal percentile tables are exact; no chain-event percentiles are inferred from coarse histograms. Per-chain length percentiles above are exact over the compact chain arrays.

## Interpretation limits

- One realization per mixing value cannot measure between-run scatter or establish ensemble significance. No bootstrap across dependent chain events is used.
- The common inputs and matching current sources support comparison, but the directory headers are the available configuration evidence, not an independently recorded compiler provenance. The logs corroborate population size, channel allocation, and named modes.
- Generic cluster structure/evolution and uniform-redshift formation sampling remain inherited assumptions. No cosmic-rate calibration or detector selection is applied.
- The reported runtime reflects completed logs and is not a controlled performance benchmark.

## Reproduction and provenance

See [README.md](README.md). `manifest.json` records raw-input SHA-256 values and line counts; `aggregate_checksums.json` verifies the portable plotting inputs. Analysis used Python 3.11.15, NumPy 1.26.4, pandas 3.0.5, and chunks of 100,000 rows. Source snapshots and logs are under `snapshots/`. Raw runs were read without modification.
