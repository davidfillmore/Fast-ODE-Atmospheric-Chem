"""Rosenbrock tableaus in the transformed (KPP) form of the paper's eq. (ros-implform).

    (M/(h gamma) - J) K_r = F(Y + sum_s a_rs K_s) + (M/h) sum_s c_rs K_s,
    Y_new = Y + sum_r m_r K_r,     error = sum_r e_r K_r.

The coefficients are those of Sandu et al. (1997), as KPP prints them.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Tableau:
    name: str
    gamma: float
    a: np.ndarray        # (S, S), strictly lower triangular
    c: np.ndarray        # (S, S), strictly lower triangular
    m: np.ndarray        # (S,)
    e: np.ndarray        # (S,)
    order: int
    embedded: int

    @property
    def stages(self) -> int:
        return len(self.m)

    @property
    def k(self) -> float:
        """The order of the error estimate, p_hat + 1: the controller's exponent is 1/k."""
        return float(self.embedded + 1)

    @property
    def new_forcing(self) -> list[bool]:
        """Stage r evaluates F anew unless its predicted point is stage r-1's."""
        flags = [True]
        for r in range(1, self.stages):
            prev = np.concatenate([self.a[r - 1, : r - 1], [0.0]])
            flags.append(not (np.array_equal(self.a[r, :r], prev) and self.a[r, r - 1] == 0.0))
        return flags

    def standard_form(self):
        """A, Gamma, b, bhat of the standard form (eq. ros-implcoef inverted)."""
        S = self.stages
        Gamma = np.linalg.inv(np.eye(S) / self.gamma - self.c)
        A = self.a @ Gamma
        b = Gamma.T @ self.m
        bhat = b - Gamma.T @ self.e
        return A, Gamma, b, bhat

    def index1(self, q: int) -> float:
        """b^T B^{-1} c^q, which the index-1 condition of that order asks to equal one."""
        A, Gamma, b, _ = self.standard_form()
        B = A + Gamma
        return float(b @ np.linalg.solve(B, A.sum(axis=1) ** q))


def _lower(S, flat):
    M = np.zeros((S, S))
    k = 0
    for i in range(1, S):
        for j in range(i):
            M[i, j] = flat[k]
            k += 1
    return M


ROS3 = Tableau(
    name="ros3",
    gamma=0.43586652150845899941601945119356,
    a=_lower(3, [1.0, 1.0, 0.0]),
    c=_lower(3, [-0.10156171083877702091975600115545e01,
                 0.40759956452537699824805835358067e01,
                 0.92076794298330791242156818474003e01]),
    m=np.array([1.0, 0.61697947043828245592553615689730e01, -0.42772256543218573326238373806514]),
    e=np.array([0.5, -0.29079558716805469821718236208017e01, 0.22354069897811569627360909276199]),
    order=3, embedded=2,
)

RODAS3 = Tableau(
    name="rodas3",
    gamma=0.5,
    a=_lower(4, [0.0, 2.0, 0.0, 2.0, 0.0, 1.0]),
    c=_lower(4, [4.0, 1.0, -1.0, 1.0, -1.0, -8.0 / 3.0]),
    m=np.array([2.0, 0.0, 1.0, 1.0]),
    e=np.array([0.0, 0.0, 0.0, 1.0]),
    order=3, embedded=2,
)
