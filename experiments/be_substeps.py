"""Backward Euler at equal accuracy: how many sub-steps per interval does it
need to reach the digits the algorithm reaches?

Each host interval is split into N equal sub-intervals, and the baseline's
backward Euler is called once on each, as a host would call it if it shortened
the chemistry step. The photolysis stays constant over the host interval, as in
every other run. Rows go to results/ts1/be_substeps.csv.

    python3 experiments/be_substeps.py [--workers 8]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))

CASES = [(s, n) for s in ("tier1-30s", "diurnal-10min", "diurnal-30min") for n in (1, 4, 16, 64)]
CASES += [("diurnal-30min", 256)]


class SubsteppedBackwardEuler:
    """The baseline's backward Euler called on N equal pieces of each interval."""

    name = "backward-euler"

    def __init__(self, n):
        from fastode import BackwardEuler
        self.n = n
        self.be = BackwardEuler()

    def new_memory(self):
        return None

    def solve(self, system, Y, dt, atol, rtol, counts, memory=None):
        status = "converged"
        for _ in range(self.n):
            Y, st = self.be.solve(system, Y, dt / self.n, atol, rtol, counts)
            if st != "converged":
                status = st
        return Y, status


def job(case):
    scenario, n = case
    from fastode.ts1 import TS1, Cells, Scenario, run_scenario, score
    from ts1_experiments import SCENARIOS
    t0 = time.time()
    ts1 = TS1()
    scn = Scenario(SCENARIOS[scenario] / scenario / "scenario.txt")
    cells = Cells(ts1, scn.cells, which=range(64))
    snaps, c, statuses = run_scenario(scn, SubsteppedBackwardEuler(n), cells, ts1, 1e-6, 1e-3, group_size=32)
    ref = scn.reference_snapshots if scn.snapshot_every else scn.reference_final
    D, per, gate = score(snaps, range(64), ts1.species, ref)
    worst = min(per, key=per.get)
    per_int = 1.0 / scn.steps / 2          # per interval and 32-cell group
    work = (10.187 * c.lu + 2.520 * c.solves + 4.785 * c.jac + 3.189 * c.rhs) * per_int
    return {"scenario": scenario, "substeps": n, "digits": round(D, 4), "worst_species": worst,
            "worst_digits": round(per[worst], 4), "level2": bool(gate and D >= 2 and per[worst] >= 1),
            "level3": bool(gate and D >= 3 and per[worst] >= 2),
            "newton_per_interval": round(c.newton * per_int, 2), "work_kmadd": round(work, 1),
            "failures": c.failures, "unconverged": c.unconverged, "failed_calls": len(statuses),
            "seconds": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    from multiprocessing import get_context
    rows = []
    cases = sorted(CASES, key=lambda c: -c[1])           # the longest first
    with get_context("spawn").Pool(args.workers) as pool:
        for row in pool.imap_unordered(job, cases):
            rows.append(row)
            print(row, flush=True)
    rows.sort(key=lambda r: (r["scenario"], r["substeps"]))
    out = ROOT / "results/ts1/be_substeps.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
