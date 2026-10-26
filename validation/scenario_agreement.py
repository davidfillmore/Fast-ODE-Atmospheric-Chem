"""Check the prototypes of the baseline solvers against MICM on a short scenario.

A mini scenario is cut from a TS1 scenario (its first CELLS cells, its first STEPS
steps). MICM runs it with its Rosenbrock loop (micm_scenario, --method
rosenbrock3 or rosenbrock4dae) and with its backward Euler (micm_be_scenario,
built from validation/micm_be_scenario.cpp); the prototypes run it with the same
tolerances and one group of CELLS cells. Prints the largest relative state
difference over species above a floor, and both sets of operation counts.

    python3 validation/scenario_agreement.py SCENARIO_DIR [--cells 32] [--steps 36]
        [--rtol 1e-3] [--atol 1e-15] [--micm BIN] [--micm-be BIN]
"""
import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fastode import ROS3, RODAS3, RestartingLoop, BackwardEuler  # noqa: E402
from fastode.ts1 import TS1, Cells, Scenario, run_scenario, read_csv_table  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("scenario")
ap.add_argument("--cells", type=int, default=32)
ap.add_argument("--steps", type=int, default=36)
ap.add_argument("--rtol", type=float, default=1e-3)
ap.add_argument("--atol", type=float, default=1e-15)
ap.add_argument("--floor", type=float, default=1e-14)
ap.add_argument("--micm", required=True, help="path to the MICM Rosenbrock scenario driver")
ap.add_argument("--micm-be", required=True)
args = ap.parse_args()

src = Scenario(Path(args.scenario) / "scenario.txt")
work = Path(tempfile.mkdtemp(prefix="agree-"))
lines = Path(src.cells).read_text().splitlines()
(work / "cells.csv").write_text("\n".join(lines[: args.cells + 1]) + "\n")
(work / "photolysis.csv").write_text((src.path.parent / "photolysis.csv").read_text())
(work / "scenario.txt").write_text(
    f"cells = cells.csv\ndt = {src.dt:g}\nsteps = {args.steps}\nphotolysis = photolysis.csv\n")

ts1 = TS1()
cells = Cells(ts1, work / "cells.csv")
scn = Scenario(work / "scenario.txt")


def micm_run(binary, method):
    out = subprocess.run([binary, "--threads", "1", "--group-size", str(args.cells), "--method", method,
                          "--rtol", repr(args.rtol), "--atol", repr(args.atol), str(work / "scenario.txt")],
                         cwd=work, capture_output=True, text=True)
    stats = dict(re.findall(r"(\w+)=([0-9.eE+-]+)", out.stdout))
    header, rows = read_csv_table(work / "ts1_final.csv")
    names = header[2:]
    Y = np.array([[r[2 + names.index(s)] for s in ts1.species] for r in rows])
    return Y, stats, out.stderr.strip()


def compare(label, Y_micm, stats, solver):
    snaps, counts, statuses = run_scenario(scn, solver, cells, ts1, args.rtol, args.atol, group_size=args.cells)
    Y = snaps[args.steps]
    mask = np.abs(Y_micm) >= args.floor
    rel = np.where(mask, np.abs(Y - Y_micm) / np.maximum(np.abs(Y_micm), 1e-300), 0.0)
    j = np.unravel_index(np.argmax(rel), rel.shape)
    print(f"{label}: max relative difference {rel.max():.2e} over values >= {args.floor:g} "
          f"(worst {ts1.species[j[1]]}); statuses {statuses[:3]}")
    print(f"   MICM:      {stats}")
    print(f"   prototype: steps {counts.steps} accepted {counts.accepted} rejected {counts.rejected} "
          f"rhs {counts.rhs} jac {counts.jac} lu {counts.lu} solves {counts.solves} newton {counts.newton}")


for method, tab in (("rosenbrock3", ROS3), ("rosenbrock4dae", RODAS3)):
    Y, stats, err = micm_run(args.micm, method)
    compare(f"{method} restarting loop", Y, stats, RestartingLoop(tab))
Y, stats, err = micm_run(args.micm_be, "be")
compare("backward Euler", Y, stats, BackwardEuler())
