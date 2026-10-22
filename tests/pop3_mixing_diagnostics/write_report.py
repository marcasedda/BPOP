#!/usr/bin/env python3
"""Write the numerical report from the bundled aggregate data."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
RUNS = ("pop3_mixing0", "pop3_mixing0.5", "pop3_mixing1.0")
MIX = dict(zip(RUNS, (0, .5, 1)))


def markdown(columns, rows):
    return "\n".join(["| " + " | ".join(columns) + " |",
                       "| " + " | ".join(["---"]*len(columns)) + " |"] +
                      ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def main():
    summary = pd.read_csv(BASE/"summary.csv").set_index(["run", "environment"])
    stats = pd.read_csv(BASE/"statistics.csv").set_index(["run", "scope", "environment", "metric"])
    validation = json.loads((BASE/"validation.json").read_text())
    manifest = json.loads((BASE/"manifest.json").read_text())
    configs = pd.read_csv(BASE/"settings.csv")
    quantiles = {run: pd.read_csv(BASE/f"{run}_terminal_quantiles.csv").set_index(["environment", "metric"])
                 for run in RUNS}
    text = ["# Pop III mixing diagnostic", "",
            "This compares the existing one-million-candidate runs at mixing 0, 0.5, and 1. "
            "Only `mixing` differs among their 66 parsed parameter macros (the header guard is excluded). "
            "Their local C++ main sources and all three stellar input datasets match byte for byte. "
            "The result is a comparison of simulator populations, not calibrated cosmic rates.", "",
            "[Individual figures and brief explanation](panels/README.md) · "
            "[Figures (9-page PDF)](diagnostics.pdf) · [Environment summaries](summary.csv) · "
            "[All settings](settings.csv) · [Validation](validation.json)", "",
            "## Main results", ""]
    rows = []
    for run in RUNS:
        row = summary.loc[(run, "all")]
        rows.append([MIX[run], f"{int(row.terminal_mergers):,}", f"{row.terminal_yield:.4%}",
                     f"{int(row.recorded_chain_events):,}", f"{int(row.recorded_chains):,}",
                     f"{row.mean_recorded_chain_length:.2f}"])
    text += [markdown(["Mixing", "Terminal mergers", "Terminal yield", "Recorded chain events",
                       "Chains with records", "Mean recorded events/chain"], rows), ""]
    baseline = summary.loc[(RUNS[0], "all")]
    for run in RUNS[1:]:
        row = summary.loc[(run, "all")]
        text.append(f"- Mixing {MIX[run]} changes terminal counts by "
                    f"{100*(row.terminal_mergers/baseline.terminal_mergers-1):+.2f}% and recorded chain "
                    f"events by {100*(row.recorded_chain_events/baseline.recorded_chain_events-1):+.2f}% "
                    "relative to mixing 0.")
    text += ["", "These two counts measure different things. The terminal catalogue contains one "
             "written endpoint per successful candidate; a recorded chain can contribute many intermediate "
             "mergers. Long chains dominate the event-weighted distributions. Per-chain figures give every "
             "recorded chain equal weight.", "", "## Input setup and mechanism", ""]
    wanted = ("N", "SEVN", "DynOvTot", "pYC", "pGC", "pNC", "mixing", "fbin", "spinlb",
              "uppergap", "bhseed", "SFRTYPE_ISO", "SFRTYPE_CLU", "TagR", "CLfill")
    text += [markdown(["Parameter", "mixing 0", "mixing 0.5", "mixing 1"],
                      [[row.parameter, row[RUNS[0]], row[RUNS[1]], row[RUNS[2]]]
                       for _, row in configs.iterrows() if row.parameter in wanted]), "",
             "`mixing` is the probability of selecting the mixed/binary-origin reservoir, not the "
             "isolated-binary fraction and not `fbin`. A uniform draw selects the reservoir for a "
             "dynamical candidate, and its companion draws use that reservoir. Both samplers reject BHs "
             "whose natal kick exceeds the cluster escape speed. Consequently, the unconditioned "
             "0.5 mixture shown on PDF page 1 need not equal the distribution of successful mergers. "
             "The spectrum arrays preserve the input correlations between mass, formation time, and kick.", ""]
    input_rows = []
    for kind in ("single", "mixed"):
        mass = stats.loc[(kind, "input_spectrum", "all", "m1")]
        kick = stats.loc[(kind, "input_spectrum", "all", "kick")]
        input_rows.append([kind, f"{int(mass.n):,}", f"{mass['mean']:.3f}",
                           f"{mass.maximum:.3f}", f"{kick['mean']:.3f}",
                           validation["spectra"][kind]["incomplete_eof_values"]])
    text += [markdown(["Reservoir", "Complete records", "Mean BH mass [Msun]", "Max BH mass [Msun]",
                       "Mean kick [km/s]", "Discarded EOF values"], input_rows), "",
             "The EOF handling matches the simulator: only an incomplete final group of five numbers "
             "is discarded. These are input-spectrum summaries, before retention selection. "
             "A smaller reservoir mass scale can change recoil, pairing times, retention, and the "
             "supply of companions; the runs do not isolate those mechanisms individually.", "",
             "## Environment-specific terminal yields", "",
             "Denominators come from `Creating N` in each run log: isolated 99,999; young 675,000; "
             "globular 135,000; nuclear 90,001. This preserves the simulator's floating-point truncation "
             "and residual allocation rather than replacing them with rounded nominal fractions.", ""]
    yield_rows = []
    for env in ("none", "young", "globular", "nuclear"):
        yield_rows.append([env] + [f"{int(summary.loc[(run, env), 'terminal_mergers']):,} "
                                  f"({summary.loc[(run, env), 'terminal_yield']:.2%})" for run in RUNS])
    text += [markdown(["Environment", "mixing 0", "mixing 0.5", "mixing 1"], yield_rows), "",
             "The isolated channel does not use the mixing reservoir. Its differences across these "
             "runs provide a stochastic reference, not a paired control: random seeds are unrecorded "
             "and candidate IDs are not paired between runs.", "",
             markdown(["Environment", "Recorded events: mixing 0", "mixing 0.5", "mixing 1"],
                      [[env] + [f"{int(summary.loc[(run, env), 'recorded_chain_events']):,}" for run in RUNS]
                       for env in ("young", "globular", "nuclear")]), "",
             "Nuclear clusters dominate the recorded event total in all three runs. Their event "
             "count falls strongly with increasing mixing, while their terminal yield remains "
             "near 89–93%. The main change is therefore in how far a small set of nuclear chains "
             "continues recycling, rather than in whether most nuclear candidates produce a terminal merger.", "",
             "## Masses, delays, and growth", ""]
    mass_rows = []
    for run in RUNS:
        row = summary.loc[(run, "all")]
        mass_rows.append([MIX[run], f"{quantiles[run].loc[('all', 'remnant'), 'p50']:.3f}",
                          f"{quantiles[run].loc[('all', 'delay'), 'p50']:.4f}",
                          f"{int(row.terminal_above_100):,}", f"{int(row.terminal_above_1000):,}",
                          f"{int(row.terminal_above_10000):,}"])
    text += [markdown(["Mixing", "Terminal median Mrem [Msun]", "Median elapsed delay [Gyr]",
                       "Mrem > 100", "Mrem > 1,000", "Mrem > 10,000"], mass_rows), ""]
    growth_rows = []
    concentration_rows = []
    for run in RUNS:
        with np.load(BASE/f"{run}_chains.npz") as chain:
            growth_rows.append([MIX[run], f"{np.median(chain['recorded_events']):.0f}",
                                f"{np.quantile(chain['recorded_events'], .99):.1f}",
                                f"{chain['recorded_events'].max():,}",
                                f"{chain['maximum_remnant_mass'].max():.6g}",
                                int((chain["maximum_remnant_mass"] > 10000).sum())])
            lengths = np.sort(chain["recorded_events"])
            top_n = int(np.ceil(.01*len(lengths)))
            concentration_rows.append([MIX[run], top_n,
                                       f"{100*lengths[-top_n:].sum()/lengths.sum():.2f}%",
                                       int((chain["maximum_remnant_mass"] > 5e7).sum())])
    text += [markdown(["Mixing", "Median chain length", "99th percentile length", "Maximum length",
                       "Largest recorded Mrem [Msun]", "Chains reaching > 10,000 Msun"], growth_rows), "",
             "Lengths count recorded rows, not the `Nrec` counter. Maximum masses and growth factors "
             "use every recorded event; no display subsampling is used. Source column 17 is retained "
             "as the recorded previous-merger counter and is not assumed to equal a chain's row count. "
             "The fixed-seed macros are inactive with the current AS20 / no-seed setup, so massive "
             "remnants here do not imply that the separate prescribed seed experiment was enabled.", "",
             markdown(["Mixing", "Top 1% of recorded chains (rounded up)", "Share of all recorded events",
                       "Chains exceeding 5e7 Msun"], concentration_rows), "",
             "The event-weighted comparisons are dominated by a very small number of long chains. "
             "The largest masses for mixing 0 and 0.5 reach the configured `msmbhmax=5e7 Msun` "
             "stopping scale: the check acts on primary mass after an event is written, so the last "
             "remnant may slightly exceed it. Those maxima are affected by the simulation's imposed "
             "growth limit and are not unconstrained physical mass predictions.", "",
             "## Catalogue semantics and reconciliation", "",
             "The column mapping is checked against the local C++ writers. Main column 15 and chain "
             "column 16 are absolute cosmic ages in years; elapsed delay subtracts formation age. "
             "Chain columns 9–15 are also written in years, despite the older schema's Myr labels. "
             "The chain file has 33 columns in this local variant.", "",
             "Chain column 23 identifies the candidate chain. Its column 31 plus one maps to the "
             "terminal catalogue ID because `itot` is incremented before the final writer. A GW-recoil "
             "exit occurs before the chain writer; a terminal merger may therefore be absent from "
             "the chain file, or differ from its last recorded event. Other termination paths also "
             "prevent simple equality. The files must not be added together as disjoint event lists.", ""]
    reconcile = []
    for run in RUNS:
        a = validation[run]
        reconcile.append([MIX[run], a["chains_mapping_to_terminal"], a["mapped_terminal_matching_last_record"],
                          a["mapped_terminal_differing_from_last_record"], a["terminal_dyn_without_recorded_chain"],
                          a["chains_without_terminal"]])
    text += [markdown(["Mixing", "Chains mapped to terminal", "Last event matches terminal", "Last event differs",
                       "Dynamical terminal without chain", "Chain without terminal"], reconcile), "",
             "The match compares both component masses, remnant mass, and cosmic age with relative "
             "tolerance 2e-5 (catalogues print about six significant digits). Non-merger record IDs "
             "have writer-specific offsets and can overlap terminal IDs; their raw record counts "
             "are reported, but they are not used as complementary denominators. "
             "`chain_events_before_present` in the CSV applies a strict comparison to the printed "
             "cosmic ages and the simulator's 13.78689 Gyr present age. A few records round across "
             "that boundary at six significant digits; they remain included in the recorded-event "
             "diagnostics and are not classified as physical future mergers.", "",
             "## Data-quality audit", ""]
    for run in RUNS:
        a = validation[run]
        text += [f"### Mixing {MIX[run]}", "",
                 f"- Full raw chain row count: {a['chain_rows']:,}; streamed and independent byte-level newline counts agree.",
                 f"- Within-chain time/generation ordering violations: {a['chronology_errors']}.",
                 f"- Terminal numeric non-finite counts (zero-based columns): `{a['terminal_numeric_nonfinite']}`.",
                 f"- Chain numeric non-finite counts: `{a['chain_numeric_nonfinite']}`.",
                 f"- Non-merger numeric non-finite counts: `{a['nonmerger_numeric_nonfinite']}`.",
                 f"- Terminal physical-bound flags: `{a['terminal_physical']}`.",
                 f"- Chain physical-bound flags: `{a['chain_physical']}`.",
                 f"- Non-merger duplicate IDs: {a['nonmerger_duplicate_ids']}; terminal/non-merger ID overlap: {a['terminal_nonmerger_id_overlap']}.", ""]
    text += ["The non-finite terminal entries are in semimajor axis (one-based column 31): 19, 10, "
             "and 9 entries across mixing 0, 0.5, and 1. Non-merger separation columns 21–22 contain "
             "many non-finite values. These columns are not used for the displayed mass, spin, or "
             "delay distributions. All recorded chain numeric columns are finite, and the checked "
             "mass/spin/eccentricity/time bounds have no violations. Spin-bound checks here cover "
             "primary spin, remnant spin, and effective spin; secondary spin is not separately aggregated.", "",
             "Non-finite entries are excluded only from the affected quantity, not silently from "
             "whole records. Physical-bound flags are reported rather than repaired. Logarithmic panels "
             "omit nonpositive values. Every quantity's finite count, non-finite count, zeros, negative "
             "values, and histogram underflow/overflow are recorded in `statistics.csv`. Plotted densities "
             "divide by the full finite count and bin width; a panel with omitted mass can therefore "
             "integrate to less than one. Terminal percentile tables are exact; no chain-event "
             "percentiles are inferred from coarse histograms. Per-chain length percentiles above "
             "are exact over the compact chain arrays.", "", "## Interpretation limits", "",
             "- One realization per mixing value cannot measure between-run scatter or establish "
             "ensemble significance. No bootstrap across dependent chain events is used.",
             "- The common inputs and matching current sources support comparison, but the directory "
             "headers are the available configuration evidence, not an independently recorded compiler "
             "provenance. The logs corroborate population size, channel allocation, and named modes.",
             "- Generic cluster structure/evolution and uniform-redshift formation sampling remain "
             "inherited assumptions. No cosmic-rate calibration or detector selection is applied.",
             "- The reported runtime reflects completed logs and is not a controlled performance benchmark.", "",
             "## Reproduction and provenance", "",
             "See [README.md](README.md). `manifest.json` records raw-input SHA-256 values and line "
             "counts; `aggregate_checksums.json` verifies the portable plotting inputs. Analysis used "
             f"Python {manifest['python']}, NumPy {manifest['numpy']}, pandas {manifest['pandas']}, "
             f"and chunks of {manifest['chunk_size']:,} rows. Source snapshots and logs are under "
             "`snapshots/`. Raw runs were read without modification.", ""]
    (BASE/"REPORT.md").write_text("\n".join(text))
    print(BASE/"REPORT.md")


if __name__ == "__main__":
    main()
