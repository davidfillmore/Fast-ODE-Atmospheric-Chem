"""fastode: pure-Python prototypes for "A Fast ODE Solver for Atmospheric Chemistry".

Three solvers on one set of systems, each counting its operations:

- ``Algorithm``: ROS3 (RODAS3 with constraints) with the step carried across
  intervals and quantized to them, the H211b filter with its memory carried,
  and the absolute tolerance chosen by the caller (Section 4 of the paper).
- ``RestartingLoop``: the baseline Rosenbrock loop, restarted every interval at
  1e-6 of the interval under the elementary controller (Section 5.1).
- ``BackwardEuler``: the baseline backward Euler solver with Newton's
  iteration and CAM-Chem's step reductions (Section 5.2).

Work is counted, never timed.
"""
from .tableaus import ROS3, RODAS3, Tableau
from .counts import Counts
from .solvers import Algorithm, RestartingLoop, BackwardEuler, project

__all__ = ["ROS3", "RODAS3", "Tableau", "Counts", "Algorithm", "RestartingLoop",
           "BackwardEuler", "project"]
