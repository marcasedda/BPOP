#!/usr/bin/env python3
"""Stream the existing mixing experiments; never modify their inputs.

The column maps below follow Main_Bpop_dev_p3.cpp, not the older schema prose.
Run with the project's Python 3.11 environment. Plotting is a separate step.
"""
from pathlib import Path
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import platform
import re
import shutil
import time

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
PROJECT = BASE.parents[1]
ENVS = ("all", "none", "young", "globular", "nuclear")
RUNS = ("pop3_mixing0", "pop3_mixing0.5", "pop3_mixing1.0")
MASS_EDGES = np.geomspace(0.1, 1e9, 161)
TIME_EDGES = np.geomspace(1e-12, 1e5, 171)  # Gyr
EDGES = {
    "m1": MASS_EDGES, "m2": MASS_EDGES, "remnant": MASS_EDGES,
    "q": np.linspace(0, 1, 51), "chi_eff": np.linspace(-1, 1, 81),
    "spin1": np.linspace(0, 1, 81), "remnant_spin": np.linspace(0, 1, 81),
    "recoil": np.geomspace(1e-7, 1e5, 121),
    "recoil_escape": np.geomspace(1e-10, 1e5, 121),
    "delay": TIME_EDGES, "formation_time": TIME_EDGES,
    "merger_time": TIME_EDGES,
    "zmerger": np.linspace(0, 20, 81), "zformation": np.linspace(0, 20, 81),
    "generation": np.r_[-0.5, 0.5, 1.5, 2.5, np.geomspace(3.5, 1e6, 91)],
    "eccentricity": np.linspace(0, 1, 81),
    "formation_delay": TIME_EDGES, "capture_delay": TIME_EDGES,
    "threebody_delay": TIME_EDGES, "friction_delay": TIME_EDGES,
    "hardening_delay": TIME_EDGES, "pairing_delay": TIME_EDGES,
    "inspiral_delay": TIME_EDGES,
    "kick": np.geomspace(1e-7, 1e5, 121),
    "chain_length": np.r_[0.5, 1.5, 2.5, np.geomspace(3.5, 1e6, 91)],
    "growth": np.geomspace(1e-2, 1e8, 141),
}


def sha(path):
    digest = hashlib.sha256()
    lines = 0
    last = b""
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
            lines += block.count(b"\n")
            last = block[-1:]
    return {"sha256": digest.hexdigest(), "bytes": path.stat().st_size,
            "lines": lines + int(bool(last) and last != b"\n")}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def parameters(path):
    return {m.group(1): m.group(2).split("//")[0].strip()
            for m in re.finditer(r"^\s*#define\s+(\w+)\s+([^\n]+)",
                                 path.read_text(), re.M)}


def frames(path, width, strings, chunk_size, header=False):
    with path.open() as stream:
        first = stream.readline()
        if header:
            first = stream.readline()
        assert len(first.split()) == width, (path, "first-row width")
    types = {i: (str if i in strings else np.float64) for i in range(width)}
    yield from pd.read_csv(path, sep=r"\s+", header=None, names=range(width),
                           skiprows=int(header), dtype=types, chunksize=chunk_size,
                           on_bad_lines="error")


class Accumulator:
    def __init__(self):
        self.values = {}

    def add(self, run, scope, env, metric, values):
        key = (run, scope, env, metric)
        if key not in self.values:
            self.values[key] = {"n": 0, "nonfinite": 0, "zero": 0,
                                "negative": 0, "sum": 0., "minimum": None,
                                "maximum": None, "underflow": 0, "overflow": 0,
                                "hist": np.zeros(len(EDGES[metric])-1, np.int64),
                                "above_100": 0, "above_1000": 0,
                                "above_10000": 0}
        out = self.values[key]
        x = np.asarray(values)
        out["nonfinite"] += int((~np.isfinite(x)).sum())
        x = x[np.isfinite(x)]
        if not len(x):
            return
        out["n"] += len(x)
        out["zero"] += int((x == 0).sum())
        out["negative"] += int((x < 0).sum())
        out["sum"] += float(x.sum())
        out["minimum"] = float(x.min()) if out["minimum"] is None else min(out["minimum"], float(x.min()))
        out["maximum"] = float(x.max()) if out["maximum"] is None else max(out["maximum"], float(x.max()))
        edges = EDGES[metric]
        out["hist"] += np.histogram(x, edges)[0]
        out["underflow"] += int((x < edges[0]).sum())
        out["overflow"] += int((x > edges[-1]).sum())
        if metric in ("m1", "m2", "remnant"):
            for threshold in (100, 1000, 10000):
                out[f"above_{threshold}"] += int((x > threshold).sum())

    def group(self, run, scope, env, metrics):
        for group in ENVS:
            mask = np.ones(len(env), bool) if group == "all" else env == group
            if mask.any():
                for metric, values in metrics.items():
                    self.add(run, scope, group, metric, values[mask])

    def save(self):
        rows, arrays = [], {}
        for i, (key, value) in enumerate(sorted(self.values.items())):
            run, scope, env, metric = key
            edges = EDGES[metric]
            counts = value["hist"]
            assert int(counts.sum()) + value["underflow"] + value["overflow"] == value["n"]
            name = f"h{i}"
            arrays[name] = counts
            arrays[name + "_edges"] = edges
            row = dict(zip(("run", "scope", "environment", "metric"), key))
            row.update({k: v for k, v in value.items() if k not in ("hist", "sum")})
            row.update(histogram_key=name, mean=value["sum"]/value["n"] if value["n"] else None)
            rows.append(row)
        pd.DataFrame(rows).to_csv(BASE / "statistics.csv", index=False)
        np.savez_compressed(BASE / "histograms.npz", **arrays)
        return rows


def numeric_audit(frame, string_cols):
    return {str(c): int((~np.isfinite(frame[c].to_numpy())).sum())
            for c in frame if c not in string_cols and (~np.isfinite(frame[c].to_numpy())).any()}


def event_metrics(frame, chain=False):
    if chain:
        m1, m2, spin1, remnant, spin, chi, recoil = (frame[c].to_numpy() for c in (0, 1, 2, 25, 26, 27, 28))
        formation, merger, generation, esc = (frame[c].to_numpy() for c in (7, 15, 16, 29))
    else:
        m1, m2, spin1, remnant, spin, chi, recoil = (frame[c].to_numpy() for c in (5, 6, 7, 9, 10, 11, 12))
        formation, merger, generation, esc = (frame[c].to_numpy() for c in (13, 14, 2, 17))
    with np.errstate(divide="ignore", invalid="ignore"):
        result = dict(m1=m1, m2=m2, spin1=spin1, remnant=remnant,
                      remnant_spin=spin, chi_eff=chi, recoil=recoil,
                      q=np.minimum(m1, m2)/np.maximum(m1, m2),
                      delay=(merger-formation)/1e9, formation_time=formation/1e9,
                      merger_time=merger/1e9, generation=generation,
                      recoil_escape=np.where(esc > 0, recoil/esc, np.nan))
    if chain:
        result.update({name: frame[col].to_numpy()/1e9 for col, name in
                       enumerate(("formation_delay", "capture_delay", "threebody_delay",
                                  "friction_delay", "hardening_delay", "pairing_delay",
                                  "inspiral_delay"), 8)})
        result["eccentricity"] = frame[31].to_numpy()
    else:
        result.update(zmerger=frame[24].to_numpy(), zformation=frame[25].to_numpy(),
                      eccentricity=frame[29].to_numpy())
    return result


def physical_audit(metrics):
    return {"nonpositive_components": int(((metrics["m1"] <= 0) | (metrics["m2"] <= 0)).sum()),
            "remnant_exceeds_total": int((metrics["remnant"] >
                                           (metrics["m1"]+metrics["m2"])*1.0001).sum()),
            "spin_outside_bounds": int(((metrics["spin1"] < 0) | (metrics["spin1"] > 1) |
                                         (metrics["remnant_spin"] < 0) | (metrics["remnant_spin"] > 1) |
                                         (abs(metrics["chi_eff"]) > 1)).sum()),
            "eccentricity_outside_bounds": int(((metrics["eccentricity"] < 0) |
                                                (metrics["eccentricity"] > 1)).sum()),
            "negative_elapsed_delay": int((metrics["delay"] < 0).sum())}


def analyze_run(root, acc, manifest, chunk_size):
    run = root.name
    sims = list(root.glob("SIM_*"))
    assert len(sims) == 1, (root, "expected one simulation")
    sim = sims[0]
    snap = BASE / "snapshots" / run
    snap.mkdir(parents=True, exist_ok=True)
    for path in (root/"input_params_p3.h", root/"Main_Bpop_dev_p3.cpp", root/"Makefile",
                 root/".build/build_paths.h", sim/"Output.txt", sim/"met_dist.txt",
                 sim/"retention_Z1e-11.dat"):
        shutil.copy2(path, snap/path.name)
        manifest[str(path.relative_to(PROJECT))] = sha(path)
    log = (sim/"Output.txt").read_text()
    assert "Elapsed time = " in log
    allocation = [int(v) for v in re.search(r"Creating N = (.*)", log)[1].split()]
    n, niso, nyc, ngc, nnc = allocation
    candidates = dict(zip(ENVS, (n, niso, nyc, ngc, nnc)))
    logged_terminal = int(re.search(r"Number of simulated mergers = (\d+)", log)[1])
    runtime = float(re.search(r"Elapsed time = ([\d.]+)s", log)[1])
    for filename in ("data_1e-11.dat", "spectrum1e-11.txt", "spectrum_cleaned00000000001.txt",
                     "Catalogue.txt", "Larger_than_tH.txt", "Catalogue_multiple_dyn.txt"):
        path = sim/filename
        print(f"Hash/count {run}/{filename}", flush=True)
        manifest[str(path.relative_to(PROJECT))] = sha(path)

    path = sim/"Catalogue.txt"
    terminal = pd.concat(frames(path, 36, (3, 4, 18), chunk_size, True), ignore_index=True)
    assert len(terminal) == logged_terminal
    assert len(terminal) == manifest[str(path.relative_to(PROJECT))]["lines"]-1
    assert terminal[0].is_unique and set(terminal[3]) <= set(ENVS[1:])
    assert (terminal[1] == 1e-11).all()
    metrics = event_metrics(terminal)
    acc.group(run, "terminal", terminal[3].to_numpy(), metrics)
    quantiles = []
    for env in ENVS:
        mask = np.ones(len(terminal), bool) if env == "all" else terminal[3].to_numpy() == env
        for metric, values in metrics.items():
            x = values[mask & np.isfinite(values)]
            quantiles.append(dict(run=run, environment=env, metric=metric,
                                  **dict(zip(("p05", "p50", "p95", "p99"),
                                             np.quantile(x, (.05, .5, .95, .99)).tolist() if len(x) else [None]*4))))
    pd.DataFrame(quantiles).to_csv(BASE/f"{run}_terminal_quantiles.csv", index=False)
    terminal_counts = Counter(terminal[3])
    statuses = terminal.groupby([3, 18]).size()
    audit = {"terminal_numeric_nonfinite": numeric_audit(terminal, (3, 4, 18)),
             "terminal_physical": physical_audit(metrics)}
    terminal_ids = terminal[0].to_numpy(dtype=np.int64)
    assert np.equal(terminal[0], terminal_ids).all()
    terminal_index = terminal.set_index(0)

    nonmerge_counts = Counter()
    nonmerge_ids = []
    nonmerge_nonfinite = Counter()
    nonmerge_rows = 0
    for frame in frames(sim/"Larger_than_tH.txt", 25, (14,), chunk_size):
        assert set(frame[14]) <= set(ENVS[1:])
        assert (frame[13] == 1e-11).all()
        nonmerge_rows += len(frame)
        nonmerge_counts.update(frame[14])
        nonmerge_nonfinite.update(numeric_audit(frame, (14,)))
        nonmerge_ids.append(frame[0].to_numpy(dtype=np.int64))
    nonmerge_ids = np.concatenate(nonmerge_ids)
    assert nonmerge_rows == manifest[str((sim/"Larger_than_tH.txt").relative_to(PROJECT))]["lines"]
    audit["nonmerger_numeric_nonfinite"] = dict(nonmerge_nonfinite)
    audit["nonmerger_duplicate_ids"] = int(len(nonmerge_ids)-len(np.unique(nonmerge_ids)))
    audit["terminal_nonmerger_id_overlap"] = int(len(np.intersect1d(terminal_ids, nonmerge_ids)))

    # Arrays are bounded by the number of requested candidates, not event count.
    count = np.zeros(n, np.int64)
    maxmass = np.zeros(n)
    firstmass = np.full(n, np.nan)
    lastmass = np.zeros(n)
    lasttime = np.zeros(n)
    firsttime = np.full(n, np.nan)
    maxgen = np.zeros(n)
    lastm1 = np.zeros(n)
    lastm2 = np.zeros(n)
    eventid = np.full(n, -1, np.int64)
    envcode = np.full(n, -1, np.int8)
    chain_nonfinite, chain_physical = Counter(), Counter()
    chain_rows = 0
    chain_counts, chain_status = Counter(), Counter()
    chain_before = Counter()
    previous_id, previous_time, previous_generation = None, None, None
    chronology_errors = 0
    # Present-age limit used by the simulator, from its shared table.
    age_table = np.loadtxt(PROJECT/"include/redshift_time.txt")
    hubble = float(age_table[0, 1])*1e9
    for frame in frames(sim/"Catalogue_multiple_dyn.txt", 33, (23, 24), chunk_size):
        env = frame[24].to_numpy()
        assert set(env) <= set(ENVS[2:])
        ids = frame[22].to_numpy(dtype=np.int64)
        assert np.equal(ids, frame[22]).all() and ((ids >= 0) & (ids < n)).all()
        values = event_metrics(frame, True)
        acc.group(run, "chain_event", env, values)
        # First events separately expose the reservoir effect before recycling.
        first = frame[16].to_numpy() == 0
        acc.group(run, "chain_first", env[first], {k: v[first] for k, v in values.items()})
        chain_nonfinite.update(numeric_audit(frame, (23, 24)))
        chain_physical.update(physical_audit(values))
        chain_counts.update(env)
        chain_before.update(env[frame[15].to_numpy() < hubble])
        chain_status.update(zip(env, frame[23]))
        # Source writes one contiguous block per chain. Validate across chunks.
        changes = np.r_[True, ids[1:] != ids[:-1]]
        starts = np.flatnonzero(changes)
        ends = np.r_[starts[1:]-1, len(ids)-1]
        blockids = ids[starts]
        assert (np.diff(blockids) > 0).all()
        if previous_id is not None:
            assert ids[0] >= previous_id
            if ids[0] == previous_id:
                chronology_errors += int(frame[15].iloc[0] < previous_time)
                chronology_errors += int(frame[16].iloc[0] < previous_generation)
        same = ids[1:] == ids[:-1]
        chronology_errors += int(((np.diff(frame[15]) < 0) & same).sum())
        chronology_errors += int(((np.diff(frame[16]) < 0) & same).sum())
        previous_id, previous_time, previous_generation = ids[-1], frame[15].iloc[-1], frame[16].iloc[-1]
        new = count[blockids] == 0
        firstmass[blockids[new]] = frame[0].to_numpy()[starts[new]]
        firsttime[blockids[new]] = frame[15].to_numpy()[starts[new]]
        count[blockids] += ends-starts+1
        maxmass[blockids] = np.maximum(maxmass[blockids], np.maximum.reduceat(frame[25].to_numpy(), starts))
        maxgen[blockids] = np.maximum(maxgen[blockids], np.maximum.reduceat(frame[16].to_numpy(), starts))
        for target, col in ((lastmass, 25), (lasttime, 15), (lastm1, 0), (lastm2, 1), (eventid, 30)):
            target[blockids] = frame[col].to_numpy()[ends]
        codes = np.array([ENVS.index(e) for e in env[starts]], np.int8)
        assert ((envcode[blockids] == -1) | (envcode[blockids] == codes)).all()
        envcode[blockids] = codes
        chain_rows += len(frame)
        if chain_rows % (chunk_size*20) == 0:
            print(f"Streamed {run}: {chain_rows:,} chain rows", flush=True)
    assert count.sum() == chain_rows
    assert chain_rows == manifest[str((sim/"Catalogue_multiple_dyn.txt").relative_to(PROJECT))]["lines"]
    assert chronology_errors == 0
    active = np.flatnonzero(count)
    mapped_ids = eventid[active]+1  # writer increments itot before terminal output
    assert len(np.unique(mapped_ids)) == len(mapped_ids)
    present = np.isin(mapped_ids, terminal_ids)
    matched = terminal_index.reindex(mapped_ids[present])
    exact_last = np.ones(present.sum(), bool)
    for data, col in ((lastm1, 5), (lastm2, 6), (lastmass, 9), (lasttime, 14)):
        exact_last &= np.isclose(data[active[present]], matched[col].to_numpy(), rtol=2e-5, atol=1e-5)
    audit.update(chain_numeric_nonfinite=dict(chain_nonfinite), chain_physical=dict(chain_physical),
                 chain_rows=chain_rows, distinct_recorded_chains=len(active), chronology_errors=chronology_errors,
                 chains_mapping_to_terminal=int(present.sum()),
                 mapped_terminal_matching_last_record=int(exact_last.sum()),
                 mapped_terminal_differing_from_last_record=int((~exact_last).sum()),
                 chains_without_terminal=int((~present).sum()),
                 terminal_dyn_without_recorded_chain=int((terminal[3] != "none").sum()-present.sum()))
    chain_data = dict(chain_id=active, recorded_events=count[active],
                      maximum_remnant_mass=maxmass[active], first_primary_mass=firstmass[active],
                      last_remnant_mass=lastmass[active], first_event_time_yr=firsttime[active],
                      last_event_time_yr=lasttime[active], maximum_recorded_generation=maxgen[active],
                      environment_code=envcode[active], terminal_id=mapped_ids,
                      has_terminal=present)
    np.savez_compressed(BASE/f"{run}_chains.npz", **chain_data)
    chainenv = np.array(ENVS)[envcode[active]]
    acc.group(run, "per_chain", chainenv,
              dict(chain_length=count[active], remnant=maxmass[active],
                   growth=maxmass[active]/firstmass[active]))
    rows = []
    for env in ENVS:
        nt = len(terminal) if env == "all" else terminal_counts[env]
        nc = chain_rows if env == "all" else chain_counts[env]
        nm = nonmerge_rows if env == "all" else nonmerge_counts[env]
        mask = np.ones(len(active), bool) if env == "all" else chainenv == env
        rem = metrics["remnant"] if env == "all" else metrics["remnant"][terminal[3].to_numpy() == env]
        rows.append(dict(run=run, mixing=float(parameters(root/"input_params_p3.h")["mixing"]),
                         environment=env, candidates=candidates[env], terminal_mergers=nt,
                         terminal_yield=nt/candidates[env], nonmerger_records=nm,
                         recorded_chain_events=nc, recorded_chains=int(mask.sum()),
                         chain_events_before_present=sum(chain_before.values()) if env == "all" else chain_before[env],
                         mean_recorded_chain_length=float(count[active[mask]].mean()) if mask.any() else None,
                         maximum_recorded_chain_length=int(count[active[mask]].max()) if mask.any() else 0,
                         maximum_recorded_remnant_mass=float(maxmass[active[mask]].max()) if mask.any() else None,
                         terminal_above_100=int((rem > 100).sum()),
                         terminal_above_1000=int((rem > 1000).sum()),
                         terminal_above_10000=int((rem > 10000).sum()), runtime_seconds=runtime))
    status_rows = [dict(run=run, scope="terminal", environment=env, status=status, count=int(value))
                   for (env, status), value in statuses.items()]
    status_rows += [dict(run=run, scope="chain_event", environment=env, status=status, count=int(value))
                    for (env, status), value in chain_status.items()]
    print(f"Finished {run}: {len(terminal):,} terminal, {chain_rows:,} recorded events", flush=True)
    return rows, status_rows, audit


def analyze_spectra(roots, acc, manifest, chunk_size):
    audit = {}
    for kind, name in (("single", "spectrum_cleaned00000000001.txt"), ("mixed", "spectrum1e-11.txt")):
        hashes = [manifest[str(next(root.glob("SIM_*"))/name).replace(str(PROJECT)+"/", "")]["sha256"] for root in roots]
        assert len(set(hashes)) == 1, (kind, "input spectra differ")
        path = next(roots[0].glob("SIM_*"))/name
        # Match the simulator: whitespace-delimited groups of five, discarding
        # only an incomplete EOF record. Keep row correlations intact.
        tail = np.empty(0)
        metrics = {name: [] for name in ("m1", "formation_time", "kick")}
        with path.open() as stream:
            while lines := stream.readlines(8*1024*1024):
                values = np.fromstring(" ".join(lines), sep=" ")
                values = np.r_[tail, values]
                end = len(values)//5*5
                array = values[:end].reshape(-1, 5)
                assert np.isfinite(array).all()
                metrics["m1"].append(array[:, 1])
                metrics["formation_time"].append(array[:, 3]/1000.)  # input Myr -> Gyr
                metrics["kick"].append(array[:, 4])
                tail = values[end:]
        metrics = {k: np.concatenate(v) for k, v in metrics.items()}
        acc.group(kind, "input_spectrum", np.repeat("none", len(metrics["m1"])), metrics)
        audit[kind] = dict(records=len(metrics["m1"]), incomplete_eof_values=len(tail))
        for root in roots:
            mix = float(parameters(root/"input_params_p3.h")["mixing"])
            # Plotting uses these component weights; retention-conditioned draws
            # are distinct from this unconditioned input mixture.
            audit[kind][root.name+"_weight"] = mix if kind == "mixed" else 1-mix
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-size", type=int, default=100000)
    args = parser.parse_args()
    assert args.chunk_size > 0
    started = time.monotonic()
    roots = [PROJECT/"tests"/run for run in RUNS]
    configs = {root.name: parameters(root/"input_params_p3.h") for root in roots}
    config_rows = [{"parameter": key, **{run: config[key] for run, config in configs.items()}}
                   for key in configs[RUNS[0]]]
    pd.DataFrame(config_rows).to_csv(BASE/"settings.csv", index=False)
    differences = [row["parameter"] for row in config_rows if len({row[run] for run in RUNS}) > 1]
    assert differences == ["mixing"], differences
    assert len({sha(root/"Main_Bpop_dev_p3.cpp")["sha256"] for root in roots}) == 1
    acc, manifest, results, statuses, validation = Accumulator(), {}, [], [], {}
    manifest["include/redshift_time.txt"] = sha(PROJECT/"include/redshift_time.txt")
    for root in roots:
        rows, status_rows, audit = analyze_run(root, acc, manifest, args.chunk_size)
        results.extend(rows)
        statuses.extend(status_rows)
        validation[root.name] = audit
    validation["spectra"] = analyze_spectra(roots, acc, manifest, args.chunk_size)
    for name in ("data_1e-11.dat", "spectrum1e-11.txt", "spectrum_cleaned00000000001.txt"):
        assert len({value["sha256"] for key, value in manifest.items() if key.endswith("/"+name)}) == 1
    acc.save()
    pd.DataFrame(results).to_csv(BASE/"summary.csv", index=False)
    pd.DataFrame(statuses).to_csv(BASE/"statuses.csv", index=False)
    dump(BASE/"validation.json", validation)
    dump(BASE/"manifest.json", dict(created_utc=datetime.now(timezone.utc).isoformat(),
                                    python=platform.python_version(), numpy=np.__version__,
                                    pandas=pd.__version__, chunk_size=args.chunk_size,
                                    elapsed_analysis_seconds=time.monotonic()-started,
                                    parameter_differences=differences, sources=manifest))
    products = list(BASE.glob("*.csv"))+list(BASE.glob("*.npz"))+list((BASE/"snapshots").rglob("*"))
    dump(BASE/"aggregate_checksums.json", {str(p.relative_to(BASE)): sha(p)["sha256"]
                                          for p in products if p.is_file()})
    print("Analysis complete", flush=True)


if __name__ == "__main__":
    main()
