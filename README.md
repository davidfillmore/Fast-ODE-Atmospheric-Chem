# A Fast ODE Solver for Atmospheric Chemistry

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23672473.svg)](https://doi.org/10.5281/zenodo.23672473)

The Python reference implementation accompanying the paper of the same name.
It integrates stiff atmospheric reaction systems over a sequence of host-model
intervals. `Algorithm(ROS3)` carries its proposed step and H211b controller
memory across intervals, rounds steps to divisors of the interval, and uses
`atol = rtol * a`, where `a` is the concentration floor used to score accuracy.

The numerical package is copied unchanged from the paper's implementation.
The restarting Rosenbrock and backward Euler baselines, constrained-system
extension with RODAS3, experiment drivers, tests, and all eleven TS1 scenarios
are included. The scenarios contain the **first 64 cells**, exactly the subset
used in the paper, with their reference final states and snapshots.

## Citation

Please cite the archived code when using this software:

Fillmore, D. (2026). *davidfillmore/Fast-ODE-Atmospheric-Chem: v2026.10.001* [Software]. Zenodo.
[https://doi.org/10.5281/zenodo.23672473](https://doi.org/10.5281/zenodo.23672473).

This DOI identifies release `v2026.10.001`, published on 2026-10-26.
Project institution: **Broomfield Digital LLC**.

The accompanying paper is available on Zenodo:
[A Fast ODE Solver for Atmospheric Chemistry](https://zenodo.org/records/23654022).
Machine-readable citation metadata is in [CITATION.cff](CITATION.cff).

## Install and run

Use Python 3.10 or newer, NumPy, and SciPy. From a clone of this repository:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[test]'
python -m pytest
python examples/ts1.py
```

The example runs the algorithm and the restarting loop on the same 64 TS1
cells, with `rtol = 1e-3` and `atol = 1e-17 mol m-3`. It reports accuracy
against the bundled reference, attempted steps per group per interval,
rejections, and failed calls. Try `--scenario diurnal-10min` for a day/night
cycle or `--cells 32` for one group.

The package installation includes the TS1 mechanism. Run the experiment
scripts and example from the repository clone to access the scenario datasets.
[requirements-tested.txt](requirements-tested.txt) records the versions used
for the import checks; to use those versions, install it before `.[test]`.

## Calling the algorithm

The state has shape `(number_of_cells, number_of_species)`. A group shares one
step size. A system supplies `rhs(Y)`, `jac(Y)`, and a factorization of
`alpha * M - J`; see [system.py](fastode/system.py) and
[chem.py](fastode/chem.py). For an unconstrained ODE, `M` is the identity.

```python
import numpy as np
from fastode import Algorithm, Counts, ROS3
from fastode.system import DenseSystem

system = DenseSystem(
    1,
    rhs=lambda Y: -Y,
    jac=lambda Y: np.full((Y.shape[0], 1, 1), -1.0),
)
solver = Algorithm(ROS3)
memory = solver.new_memory()
counts = Counts()
Y = np.array([[1.0]])
rtol, a = 1e-3, 1e-14
atol = np.full(system.n, rtol * a)
for interval in range(10):
    Y, status = solver.solve(system, Y, 0.1, atol, rtol, counts, memory)
    if status != "converged":
        raise RuntimeError(status)
```

Retain **the same `memory` for each group** across successive host intervals.
Create a separate memory for each independent group or trajectory. Omitting
memory restarts the step and controller on each call. The caller sets the
absolute tolerance; `Algorithm` does not choose it automatically.

The implementation is in [solvers.py](fastode/solvers.py), the ROS3 and RODAS3
coefficients in [tableaus.py](fastode/tableaus.py), and the TS1 rate constants,
scenario runner, and scoring in [ts1.py](fastode/ts1.py). `run_scenario` retains
memory for each group and updates photolysis at host-interval boundaries.
The implementation counts numerical operations. The paper's work comparisons
use these counts rather than Python execution times.

## Reproduce the experiments

A short comparison, with the same two groups of 32 cells used in the paper:

```sh
python experiments/ts1_experiments.py --configs alg,restart-p \
  --scenarios tier1-30s --rtols 1e-3 --workers 1
```

Run the full TS1 comparison and tolerance sweeps:

```sh
python experiments/ts1_experiments.py --configs main --workers 7
python experiments/ts1_experiments.py --configs atol --workers 2 \
  --out results/ts1/atol.csv
python experiments/ts1_rejections.py
python experiments/be_substeps.py --workers 7
python experiments/make_tables.py
```

These use all eleven scenarios and 64 cells by default. The main and absolute
tolerance sweeps total 588 runs. Outputs are written under `results/` and the
generated LaTeX tables under `tables/`. The sweep driver appends rows and skips
completed jobs: use a new `--out` filename to repeat a sweep. Smaller runs may
use `--cells 1..64`; aggregate paper tables require the full sweeps.

The additional small-system and constrained-system experiments are:

```sh
python experiments/small_problems.py
python experiments/dae_experiments.py all
python experiments/ts1_dae.py
python experiments/make_tables_dae.py
```

The original experiment records are in [reference-results/](reference-results/),
separate from fresh runs. To regenerate tables from those records immediately:

```sh
python experiments/make_tables.py \
  --runs reference-results/ts1/runs.csv \
  --atol reference-results/ts1/atol.csv \
  --be-substeps reference-results/ts1/be_substeps.csv
python experiments/make_tables_dae.py --results reference-results
```

## Checks and provenance

```sh
python -m pytest
python tools/check_tableaus.py
python validation/tier1_final.py
```

The tests cover order conditions, observed ODE and DAE convergence orders,
invariants, projection, step rounding, controller memory, rejections, and
minimum-step and attempt-budget limits. The last command compares the
restarting loop with the included 100-cell MICM reference. See
[validation/README.md](validation/README.md) for optional comparisons with a
separate MICM build.

[provenance/README.md](provenance/README.md) identifies the source commit and
the small changes needed to make the drivers portable.
[provenance/verification.md](provenance/verification.md) records the import
checks and their results.
[data/scenarios/README.md](data/scenarios/README.md) describes the datasets.
The original validation report is preserved in
[reference-results/validation.md](reference-results/validation.md).

This repository uses the existing [GPLv3 license](LICENSE). Attribution and
the original license and notice for the MICM-derived mechanism are retained
under [third-party/](third-party/README.md).
