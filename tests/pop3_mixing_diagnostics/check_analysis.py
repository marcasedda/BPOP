#!/usr/bin/env python3
"""Check bounded raw samples and the complete saved per-chain products."""
from collections import Counter
import importlib.util
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("mixing_analysis", BASE/"analyze.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def reduce_sample(frame, chunk_size):
    acc = analysis.Accumulator()
    counts, maximum, first, last = Counter(), {}, {}, {}
    for start in range(0, len(frame), chunk_size):
        block = frame.iloc[start:start+chunk_size]
        acc.group("sample", "chain_event", block[24].to_numpy(), analysis.event_metrics(block, True))
        for chain, group in block.groupby(22, sort=False):
            counts[int(chain)] += len(group)
            maximum[int(chain)] = max(maximum.get(int(chain), 0), group[25].max())
            first.setdefault(int(chain), group[0].iloc[0])
            last[int(chain)] = group[25].iloc[-1]
    return acc, counts, maximum, first, last


def main():
    summary = pd.read_csv(BASE/"summary.csv").set_index(["run", "environment"])
    results = {}
    for run in analysis.RUNS:
        path = next((analysis.PROJECT/"tests"/run).glob("SIM_*"))/"Catalogue_multiple_dyn.txt"
        # 250,003 rows is deliberately not a multiple of either chunk size.
        sample = next(analysis.frames(path, 33, (23, 24), 250003))
        one, *reference = reduce_sample(sample, len(sample))
        many, *split = reduce_sample(sample, 7777)
        assert reference == split
        for key in one.values:
            a, b = one.values[key], many.values[key]
            assert np.array_equal(a["hist"], b["hist"])
            for name in a:
                if name == "hist":
                    continue
                if name == "sum":
                    assert np.isclose(a[name], b[name], rtol=1e-12, atol=1e-8)
                else:
                    assert a[name] == b[name], (key, name)
        # Compare the actual full-pass reduction to independent pandas groups
        # for every chain wholly contained in the bounded raw prefix.
        groups = sample.groupby(22, sort=False)
        complete_ids = sample[22].unique()[:-1].astype(np.int64)
        with np.load(BASE/f"{run}_chains.npz") as chain:
            idx = np.searchsorted(chain["chain_id"], complete_ids)
            assert np.array_equal(chain["chain_id"][idx], complete_ids)
            assert np.array_equal(chain["recorded_events"][idx], groups.size().loc[complete_ids])
            assert np.array_equal(chain["maximum_remnant_mass"][idx], groups[25].max().loc[complete_ids])
            assert np.array_equal(chain["first_primary_mass"][idx], groups[0].first().loc[complete_ids])
            assert np.array_equal(chain["last_remnant_mass"][idx], groups[25].last().loc[complete_ids])
            total = summary.loc[(run, "all")]
            assert len(chain["chain_id"]) == total.recorded_chains
            assert chain["recorded_events"].sum() == total.recorded_chain_events
            assert chain["maximum_remnant_mass"].max() == total.maximum_recorded_remnant_mass
            for code, env in enumerate(analysis.ENVS):
                if code == 0:
                    continue
                mask = chain["environment_code"] == code
                assert mask.sum() == summary.loc[(run, env), "recorded_chains"]
                assert chain["recorded_events"][mask].sum() == summary.loc[(run, env), "recorded_chain_events"]
        results[run] = dict(raw_sample_rows=len(sample), alternate_chunk_size=7777,
                            complete_sample_chains_checked=len(complete_ids),
                            histogram_chunk_independence=True, full_chain_summary_consistency=True)
    (BASE/"analysis_checks.json").write_text(json.dumps(results, indent=2)+"\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
