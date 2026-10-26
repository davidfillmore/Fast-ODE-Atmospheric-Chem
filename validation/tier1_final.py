"""Validate the TS1 mechanism and the restarting-loop prototype against MICM.

MICM's run of the Tier 1 benchmark (benchmarks/reference/README.md): 100 cells of
ts1_cells_100_seed1.csv in one group, one untimed warm-up step and 100 steps of
30 s, MICM's three-stage Rosenbrock with its default tolerances (relative 1e-6,
absolute 1e-3 mol m-3). The prototype must reproduce MICM's final state, which
checks the mechanism, the rate constants, and the loop at once.

    python3 validation/tier1_final.py [--cells 100]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fastode import RestartingLoop, ROS3, Counts  # noqa: E402
from fastode.chem import ChemSystem  # noqa: E402
from fastode.ts1 import TS1, Cells, read_csv_table  # noqa: E402

REF = Path(__file__).resolve().parents[1] / "data/validation"

ap = argparse.ArgumentParser()
ap.add_argument("--cells", type=int, default=100)
ap.add_argument("--steps", type=int, default=101)
args = ap.parse_args()

ts1 = TS1()
cells = Cells(ts1, REF / "ts1_cells_100_seed1.csv", which=range(args.cells))
k = ts1.rate_constants(cells.T, cells.P, cells.params)
system = ChemSystem(ts1.mech, k, dense=False)
solver = RestartingLoop(ROS3)
Y = cells.Y.copy()
counts = Counts()
atol = np.full(ts1.mech.n, 1e-3)
t0 = time.time()
for step in range(args.steps):
    Y, status = solver.solve(system, Y, 30.0, atol, 1e-6, counts)
    if status != "converged":
        print("step", step, status)
        break
print(f"{args.steps} steps in {time.time() - t0:.1f} s; counts {counts.as_dict()}")

header, micm = read_csv_table(REF / "micm_ts1_final_100_seed1.csv")
names = header[2:]
micm = micm[: args.cells, 2:]
ours = np.array([[Y[c, ts1.index[s]] for s in names] for c in range(args.cells)])
diff = np.abs(ours - micm) / np.maximum(np.abs(micm), 1e-30)
worst = np.unravel_index(np.argmax(diff), diff.shape)
print(f"largest relative difference {diff.max():.3e} (cell {worst[0]}, {names[worst[1]]})")
big = diff > 1e-10
print(f"entries above 1e-10: {int(big.sum())} of {diff.size}")
