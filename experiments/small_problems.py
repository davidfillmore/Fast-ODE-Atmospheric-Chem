"""Experiments I on the small problems: the three solvers, called once per host
interval, against Radau references; digits and operation counts.

    python3 experiments/small_problems.py [--out results/small/runs.csv]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fastode import Algorithm, BackwardEuler, Counts, ROS3, RestartingLoop  # noqa: E402
from fastode import problems as P  # noqa: E402
from fastode.reference import digits, radau_reference  # noqa: E402

M_AIR = 2.5e19
TROPO_FLOOR = 1e-14 * 6.02214076e23 * 1e-6      # 1e-14 mol m-3 in molecules cm-3


def tropo_diurnal(system, m, dt=600.0):
    P.tropospheric_photolysis(system, P.daylight(m * dt))


# name -> (builder, dt, intervals, threshold, on_interval)
PROBLEMS = {
    "robertson": (lambda: P.robertson(), 1.0e3, 100, 1e-12, None),
    "chapman": (lambda: P.chapman(), 1.0e5, 100, 1e-12, None),
    "tropospheric": (lambda: P.tropospheric(), 600.0, 144, TROPO_FLOOR, None),
    "tropospheric-diurnal": (lambda: P.tropospheric(), 600.0, 144, TROPO_FLOOR, tropo_diurnal),
    "triples": (lambda: P.triples(N=10), 0.2, 100, 1e-12, None),
    "hires": (lambda: P.hires(), 321.8122 / 100, 100, 1e-12, None),
}

# name -> (solver, absolute tolerance as a multiple of the floor a; None is rtol times a,
# the algorithm's rule)
SOLVERS = {
    "algorithm": (lambda: Algorithm(ROS3), None),
    "algorithm-tiny": (lambda: Algorithm(ROS3), 1e-6),
    "algorithm-a/10": (lambda: Algorithm(ROS3), 0.1),
    "restarting": (lambda: RestartingLoop(ROS3), None),
    "restarting-a/10": (lambda: RestartingLoop(ROS3), 0.1),
    "backward-euler": (lambda: BackwardEuler(), 1.0),
}


def run(system, Y0, solver, dt, intervals, rtol, atol, on_interval):
    Y = Y0.copy()
    counts = Counts()
    mem = solver.new_memory()
    states, bad = [], 0
    for m in range(intervals):
        if on_interval is not None:
            on_interval(system, m)
        Y, status = solver.solve(system, Y, dt, atol, rtol, counts, mem)
        bad += status != "converged"
        states.append(Y[0].copy())
    return np.array(states), counts, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "results/small/runs.csv"))
    args = ap.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["problem", "solver", "rtol", "atol", "intervals", "digits", "worst", "worst_digits",
              "failed_calls", "steps", "accepted", "rejected", "rhs", "jac", "lu", "solves", "newton",
              "failures", "unconverged"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for name, (build, dt, intervals, a, on_interval) in PROBLEMS.items():
            system, Y0 = build()
            ref = radau_reference(system, Y0, dt, intervals, on_interval=on_interval)
            n = system.n
            for sname, (make, atol_factor) in SOLVERS.items():
                settings = ([(1e-6, a * atol_factor)] if sname == "backward-euler" else
                            [(r, a * (r if atol_factor is None else atol_factor)) for r in (1e-2, 1e-3, 1e-4)])
                for rtol, atol in settings:
                    system, Y0 = build()
                    X, counts, bad = run(system, Y0, make(), dt, intervals, rtol, np.full(n, atol), on_interval)
                    D, per = digits(X, ref, a)
                    j = min(per, key=per.get) if per else None
                    c = counts.as_dict()
                    row = {"problem": name, "solver": sname, "rtol": rtol, "atol": atol, "intervals": intervals,
                           "digits": round(D, 3), "worst": system.names[j] if j is not None else "",
                           "worst_digits": round(per[j], 3) if j is not None else "", "failed_calls": bad,
                           **{k: c[k] for k in ("steps", "accepted", "rejected", "rhs", "jac", "lu", "solves",
                                                "newton", "failures", "unconverged")}}
                    w.writerow(row)
                    f.flush()
                    print(f"{name:22s} {sname:15s} rtol={rtol:g}: {D:6.2f} digits "
                          f"(worst {row['worst']} {row['worst_digits']}), steps/int {c['steps'] / intervals:.2f}, "
                          f"lu/int {c['lu'] / intervals:.2f}", flush=True)


if __name__ == "__main__":
    main()
