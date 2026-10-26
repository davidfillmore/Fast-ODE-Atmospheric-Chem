"""TS1 with constrained radical families: 4, 9, and 17 species held in
quasi-steady state, on the 30-second scenario (noon photolysis throughout, where
the constraints have index one), against the full ODE.

    python3 experiments/ts1_dae.py [--cells 64] [--rtols 1e-3,1e-4]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from fastode import Algorithm, RestartingLoop, ROS3, RODAS3  # noqa: E402
from fastode.ts1 import TS1, Cells, Scenario, run_scenario, score  # noqa: E402
from dae_experiments import FAMILIES  # noqa: E402

SCN = ROOT / "data/scenarios/tier1-30s/scenario.txt"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=64)
    ap.add_argument("--rtols", default="1e-3,1e-4")
    ap.add_argument("--floor", type=float, default=1e-14, help="the absolute tolerance is rtol times this")
    args = ap.parse_args()
    ts1 = TS1()
    scn = Scenario(SCN)
    cells = Cells(ts1, scn.cells, which=range(args.cells))
    groups = (args.cells + 31) // 32
    out = ROOT / "results/dae/ts1_families.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    configs = [("ode", ROS3, Algorithm, ()), ("ode", RODAS3, Algorithm, ())]
    configs += [(f"dae-{k}", RODAS3, Algorithm, tuple(fam)) for k, fam in FAMILIES.items()]
    configs += [("dae-4", RODAS3, RestartingLoop, tuple(FAMILIES[4])), ("ode", RODAS3, RestartingLoop, ())]
    for rtol in (float(r) for r in args.rtols.split(",")):
        for form, tab, Solver, fam in configs:
            snaps, counts, statuses = run_scenario(scn, Solver(tab), cells, ts1, rtol, rtol * args.floor, group_size=32,
                                                   constrained=fam)
            D, per, gate = score(snaps, range(args.cells), ts1.species, scn.reference_final)
            worst = min(per, key=per.get)
            fam_digits = min(per[s] for s in fam if s in per) if fam else float("nan")
            row = {"rtol": rtol, "form": form, "tableau": tab.name, "solver": Solver.name, "digits": round(D, 3),
                   "worst": worst, "worst_digits": round(per[worst], 3),
                   "constrained_min_digits": round(fam_digits, 3),
                   "steps_per_interval": round(counts.steps / scn.steps / groups, 3),
                   "lu_per_interval": round(counts.lu / scn.steps / groups, 3),
                   "proj_passes_per_interval": round(counts.proj_passes / scn.steps / groups, 3),
                   "failed_calls": len(statuses), "first_failure": statuses[0][2] if statuses else ""}
            rows.append(row)
            print(row, flush=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
