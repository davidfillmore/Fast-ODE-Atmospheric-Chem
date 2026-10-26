"""TS1 experiments: every solver configuration on every scenario, in parallel.

Each job runs one configuration on one scenario's first CELLS cells, in groups of
32 cells that share a step, scores the snapshots against the scenario's reference
over those cells, and returns one row: digits overall, the worst species, the
gate, and the operation counts. Rows go to results/ts1/runs.csv (appended);
a job whose row is already there is skipped.

    python3 experiments/ts1_experiments.py [--configs alg,restart,...] [--cells 128]
        [--scenarios all|disclosed|NAME,...] [--workers 8] [--out results/ts1/runs.csv]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# All eleven paper scenarios are distributed together. Keep the original
# disclosed/held-out grouping only for the --scenarios disclosed selector.
SCENARIO_ROOT = ROOT / "data/scenarios"
DISCLOSED = SCENARIO_ROOT
HELDOUT = SCENARIO_ROOT
DISCLOSED_NAMES = ("tier1-30s", "diurnal-30min", "diurnal-10min",
                   "diurnal-5min", "polluted-10min", "clean-10min")
SCENARIOS = {
    "tier1-30s": DISCLOSED, "diurnal-30min": DISCLOSED, "diurnal-10min": DISCLOSED,
    "diurnal-5min": DISCLOSED, "polluted-10min": DISCLOSED, "clean-10min": DISCLOSED,
    "summer-15min": HELDOUT, "tropics-20min": HELDOUT, "winter60-10min": HELDOUT,
    "polluted-48h-30min": HELDOUT, "clean-48h-5min": HELDOUT,
}

# name -> (solver kind, options, atol, list of rtol); atol None is rtol times the
# scoring floor, the algorithm's rule (Section 4)
FLOOR = 1e-14
RT_ALG = [3e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4]
RT_LOOP = [3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4]
RT_ABL = [1e-2, 2e-3, 5e-4]
CONFIGS = {
    "alg": ("algorithm", {}, None, RT_ALG),
    "alg-a20": ("algorithm", {}, 1e-20, RT_ALG),
    "alg-a15": ("algorithm", {}, 1e-15, RT_ALG),
    "alg-elementary": ("algorithm", {"controller": "elementary"}, None, RT_ABL),
    "alg-reset": ("algorithm", {"carry": False, "quantize": False}, None, RT_ABL),
    "alg-noquant": ("algorithm", {"quantize": False}, None, RT_ABL),
    "alg-nomemory": ("algorithm", {"carry_memory": False}, None, RT_ABL),
    "restart-p": ("restarting", {}, None, RT_LOOP),
    "restart": ("restarting", {}, 1e-15, RT_LOOP),
    "restart-a20": ("restarting", {}, 1e-20, [1e-2, 1e-3, 3e-4]),
    "be": ("backward-euler", {}, 1e-3, [1e-6]),
    "be-tight": ("backward-euler", {}, 1e-15, [1e-3, 1e-6]),
    # the absolute-tolerance sweep, on a subset of scenarios (--scenarios)
    "alg-a12": ("algorithm", {}, 1e-12, [1e-2, 2e-3]),
    "alg-a13": ("algorithm", {}, 1e-13, [1e-2, 2e-3]),
    "alg-a14": ("algorithm", {}, 1e-14, [1e-2, 2e-3]),
    "alg-a16": ("algorithm", {}, 1e-16, [1e-2, 2e-3]),
    "alg-a17": ("algorithm", {}, 1e-17, [1e-2, 2e-3]),
    "alg-a18": ("algorithm", {}, 1e-18, [1e-2, 2e-3]),
    "restart-a13": ("restarting", {}, 1e-13, [1e-2, 1e-3, 1e-4]),
}
MAIN = ["alg", "restart-p", "alg-a20", "alg-a15", "restart", "restart-a20", "be", "be-tight",
        "alg-elementary", "alg-reset", "alg-noquant", "alg-nomemory"]
ATOL_SWEEP = ["alg-a12", "alg-a13", "alg-a14", "alg-a16", "alg-a17", "alg-a18", "restart-a13"]
ATOL_SCENARIOS = ["diurnal-10min", "polluted-10min", "tropics-20min", "winter60-10min"]

FIELDS = ["config", "solver", "rtol", "atol", "scenario", "cells", "intervals", "digits", "worst_species",
          "worst_digits", "gate", "level2", "level3", "failed_calls", "first_failure", "steps", "accepted",
          "rejected", "rhs", "jac", "lu", "solves", "newton", "failures", "unconverged", "seconds",
          "worst5"]


def make_solver(kind, options):
    from fastode import Algorithm, RestartingLoop, BackwardEuler, ROS3
    if kind == "algorithm":
        return Algorithm(ROS3, **options)
    if kind == "restarting":
        return RestartingLoop(ROS3)
    return BackwardEuler()


def job(spec):
    config, kind, options, atol, rtol, scenario, ncells = spec
    atol = rtol * FLOOR if atol is None else atol
    from fastode.ts1 import TS1, Cells, Scenario, run_scenario, score
    t0 = time.time()
    ts1 = TS1()
    scn = Scenario(SCENARIOS[scenario] / scenario / "scenario.txt")
    cells = Cells(ts1, scn.cells, which=range(ncells))
    solver = make_solver(kind, options)
    snaps, counts, statuses = run_scenario(scn, solver, cells, ts1, rtol, atol, group_size=32)
    ref = scn.reference_snapshots if scn.snapshot_every else scn.reference_final
    D, per, gate = score(snaps, range(ncells), ts1.species, ref)
    worst = sorted(per.items(), key=lambda kv: kv[1])
    ws, wd = worst[0] if worst else ("", float("nan"))
    level = lambda L: bool(gate and D >= L and wd >= L - 1)
    c = counts.as_dict()
    return {
        "config": config, "solver": kind, "rtol": rtol, "atol": atol, "scenario": scenario,
        "cells": ncells, "intervals": scn.steps, "digits": round(D, 4), "worst_species": ws,
        "worst_digits": round(wd, 4), "gate": gate, "level2": level(2), "level3": level(3),
        "failed_calls": len(statuses), "first_failure": statuses[0][2] if statuses else "",
        **{k: c[k] for k in ("steps", "accepted", "rejected", "rhs", "jac", "lu", "solves", "newton",
                             "failures", "unconverged")},
        "seconds": round(time.time() - t0, 1),
        "worst5": json.dumps([(s, round(d, 3)) for s, d in worst[:5]]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="main", help="main, atol, or a list of names")
    ap.add_argument("--rtols", default="", help="restrict to these relative tolerances")
    ap.add_argument("--cells", type=int, default=64)
    ap.add_argument("--scenarios", default="all")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=str(ROOT / "results/ts1/runs.csv"))
    args = ap.parse_args()
    names = (list(SCENARIOS) if args.scenarios == "all" else
             list(DISCLOSED_NAMES) if args.scenarios == "disclosed" else
             args.scenarios.split(","))
    if not 1 <= args.cells <= 64:
        ap.error("the bundled paper datasets contain 64 cells; use --cells 1..64")
    unknown = set(names) - SCENARIOS.keys()
    if unknown:
        ap.error("unknown scenarios: " + ", ".join(sorted(unknown)))
    only = {float(r) for r in args.rtols.split(",") if r}
    out = Path(args.out)
    done = set()
    if out.exists():    # rows already in the output are not run again
        with open(out, newline="") as f:
            done = {(r["config"], float(r["rtol"]), r["scenario"], int(r["cells"])) for r in csv.DictReader(f)}
    specs = []
    configs = MAIN if args.configs == "main" else ATOL_SWEEP if args.configs == "atol" else args.configs.split(",")
    for config in configs:
        kind, options, atol, rtols = CONFIGS[config]
        for rtol in rtols:
            if only and rtol not in only:
                continue
            for s in (ATOL_SCENARIOS if config in ATOL_SWEEP else names):
                if (config, rtol, s, args.cells) not in done:
                    specs.append((config, kind, options, atol, rtol, s, args.cells))
    # longest scenarios first, so the pool drains evenly
    weight = {"clean-48h-5min": 6, "diurnal-5min": 3, "tier1-30s": 1}
    specs.sort(key=lambda sp: -weight.get(sp[5], 2))
    print(f"{len(specs)} jobs ({len(done)} rows already in {out.name})", flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    new = not out.exists()
    from multiprocessing import get_context
    with get_context("spawn").Pool(args.workers) as pool, open(out, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        for i, row in enumerate(pool.imap_unordered(job, specs), 1):
            w.writerow(row)
            f.flush()
            print(f"[{i}/{len(specs)}] {row['config']} rtol={row['rtol']:g} {row['scenario']}: "
                  f"{row['digits']} ({row['worst_species']} {row['worst_digits']}) "
                  f"steps/int {row['steps'] / row['intervals']:.2f} {row['seconds']}s", flush=True)


if __name__ == "__main__":
    main()
