"""Systems: a group of C cells that share one step, with state Y of shape (C, n).

A system supplies the right-hand side, its Jacobian, and factorizations of
alpha*M - J, where M = diag(mass) holds 1 for a differential and 0 for an
algebraic (constrained) variable.
"""
from __future__ import annotations

import numpy as np

from .linalg import DenseLU, dense_matrix


class System:
    n: int
    mass: np.ndarray
    names: list

    def rhs(self, Y: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def jac(self, Y: np.ndarray):
        raise NotImplementedError

    def factor(self, alpha: float, J, algebraic_rows_only: bool = False):
        raise NotImplementedError

    @property
    def has_constraints(self) -> bool:
        return bool(np.any(self.mass == 0.0))


class DenseSystem(System):
    """A system given by functions rhs(Y) -> (C, n) and jac(Y) -> (C, n, n)."""

    def __init__(self, n, rhs, jac, mass=None, names=None):
        self.n = n
        self._rhs = rhs
        self._jac = jac
        self.mass = np.ones(n) if mass is None else np.asarray(mass, dtype=float)
        self.names = names or [f"y{i}" for i in range(n)]

    def rhs(self, Y):
        return self._rhs(Y)

    def jac(self, Y):
        return self._jac(Y)

    def factor(self, alpha, J, algebraic_rows_only=False):
        return DenseLU(dense_matrix(alpha, self.mass, J, algebraic_rows_only))
