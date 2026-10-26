"""Where the algorithm's rejections fall within an interval, on TS1.

Runs the algorithm (atol = rtol * a) on two day-and-night scenarios, 64 cells in
two groups, and counts the rejected attempts on the first step of each interval
and on later steps. Rows go to results/ts1/rejections.csv.

    python3 experiments/ts1_rejections.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
import fastode.solvers as S  # noqa: E402
from fastode import Algorithm, ROS3  # noqa: E402
from fastode.ts1 import TS1, Cells, Scenario, run_scenario  # noqa: E402
from ts1_experiments import FLOOR, SCENARIOS  # noqa: E402

CASES = [("diurnal-10min", 5e-3), ("diurnal-10min", 1e-3), ("polluted-48h-30min", 5e-3)]


class Counting(Algorithm):
    """The algorithm, recording every error estimate of every call."""

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.calls = []

    def solve(self, *args, **kw):
        log = []
        norm = S.error_norm

        def spy(*a):
            value = norm(*a)
            log.append(value)
            return value

        S.error_norm = spy
        try:
            return super().solve(*args, **kw)
        finally:
            S.error_norm = norm
            self.calls.append(log)


def main():
    ts1 = TS1()
    rows = []
    for scenario, rtol in CASES:
        scn = Scenario(SCENARIOS[scenario] / scenario / "scenario.txt")
        cells = Cells(ts1, scn.cells, which=range(64))
        solver = Counting(ROS3)
        run_scenario(scn, solver, cells, ts1, rtol, rtol * FLOOR, group_size=32)
        first = later = 0
        for errs in solver.calls:
            i = 0
            while i < len(errs) and errs[i] >= 1.0:     # attempts before the first acceptance
                first += 1
                i += 1
            later += sum(e >= 1.0 for e in errs[i + 1:])
        rows.append({"scenario": scenario, "rtol": rtol, "calls": len(solver.calls),
                     "rejected_first_step": first, "rejected_later": later})
        print(f"{scenario} rtol={rtol:g}: {first} rejections on the first step of an interval, "
              f"{later} later, in {len(solver.calls)} calls", flush=True)
    out = ROOT / "results/ts1/rejections.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
