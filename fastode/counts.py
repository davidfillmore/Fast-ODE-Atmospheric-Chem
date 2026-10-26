"""Operation counters. Every count is per call of the counted operation on a whole
group of cells: one factorization of a group is one factorization of each of its
cells, so counts per group and interval are counts per cell and interval."""
from __future__ import annotations

from dataclasses import dataclass, fields


@dataclass
class Counts:
    steps: int = 0          # attempted steps (accepted + rejected)
    accepted: int = 0
    rejected: int = 0
    rhs: int = 0            # right-hand-side (forcing) evaluations
    jac: int = 0            # Jacobian evaluations
    lu: int = 0             # factorizations
    solves: int = 0         # triangular-solve pairs
    newton: int = 0         # Newton iterations (backward Euler)
    failures: int = 0       # Newton convergence failures (backward Euler)
    unconverged: int = 0    # unconverged iterates accepted (backward Euler)
    proj_passes: int = 0    # projection passes (consistent initialization)
    proj_trials: int = 0    # projection line-search trials
    calls: int = 0          # solver calls (intervals)

    def add(self, other: "Counts") -> None:
        for f in fields(self):
            setattr(self, f.name, getattr(self, f.name) + getattr(other, f.name))

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}
