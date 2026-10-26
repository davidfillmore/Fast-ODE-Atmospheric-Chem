"""Compare ROS3 with carried step and memory against the restarting loop."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fastode import Algorithm, ROS3, RestartingLoop
from fastode.ts1 import Cells, Scenario, TS1, run_scenario, score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", default="tier1-30s")
    parser.add_argument("--cells", type=int, default=64)
    parser.add_argument("--rtol", type=float, default=1e-3)
    args = parser.parse_args()
    if not 1 <= args.cells <= 64:
        parser.error("the bundled paper datasets contain 64 cells; use --cells 1..64")
    if args.rtol <= 0:
        parser.error("--rtol must be positive")
    root = Path(__file__).resolve().parents[1]
    path = root / "data/scenarios" / args.scenario / "scenario.txt"
    if not path.is_file():
        parser.error(f"unknown bundled scenario: {args.scenario}")
    scenario, ts1 = Scenario(path), TS1()
    cells = Cells(ts1, scenario.cells, which=range(args.cells))
    reference = (scenario.reference_snapshots if scenario.snapshot_every
                 else scenario.reference_final)
    failed = False
    for solver in (Algorithm(ROS3), RestartingLoop(ROS3)):
        # a = 1e-14 mol m-3 is the paper's scoring floor. run_scenario
        # retains a separate solver.new_memory() for each group across calls.
        snapshots, counts, statuses = run_scenario(
            scenario, solver, cells, ts1, args.rtol, args.rtol * 1e-14,
            group_size=32,
        )
        digits, per_species, gate = score(
            snapshots, range(args.cells), ts1.species, reference,
        )
        worst = min(per_species, key=per_species.get)
        groups = (args.cells + 31) // 32
        print(f"{solver.name}: {digits:.4f} digits; worst {worst} "
              f"{per_species[worst]:.4f}; {counts.steps / groups / scenario.steps:.3f} "
              f"steps/group/interval; {counts.rejected} rejected; "
              f"{len(statuses)} failed calls; gate={gate}")
        failed |= bool(statuses) or not gate
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
