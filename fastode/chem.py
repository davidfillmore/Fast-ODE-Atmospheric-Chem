"""Mechanisms: a reaction list compiled to vectorized rates, forcing, and Jacobian.

A reaction with rate constant k and reactants A_1 ... A_p proceeds at the rate
r = k x_1 ... x_p; it removes one of each reactant per event and adds the yield
of each product (eq. quadform for at most two reactants). Third bodies, such as
M, are not state variables: their concentration is folded into k by the caller.
The Jacobian has a fixed sparsity pattern, and its values are a sparse linear map
of the partial derivatives of the rates.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .linalg import DenseLU, SparseLU, SparsePattern, dense_matrix
from .system import System


class Mechanism:
    def __init__(self, species, reactions, third_bodies=("M",)):
        """species: state species names; reactions: (reactants, products) with
        reactants a list of names (repeats allowed) and products a list of
        (name, yield)."""
        self.species = list(species)
        n = self.n = len(self.species)
        index = {s: i for i, s in enumerate(self.species)}
        R = self.R = len(reactions)
        lists = []
        self.third_body_power = np.zeros(R)
        for j, (reactants, _) in enumerate(reactions):
            lst = []
            for s in reactants:
                if s in third_bodies:
                    self.third_body_power[j] += 1
                else:
                    lst.append(index[s])
            lists.append(lst)
        P = self.P = max(1, max(len(lst) for lst in lists))
        self.react = np.full((R, P), n, dtype=np.int64)       # n points at a column of ones
        for j, lst in enumerate(lists):
            self.react[j, : len(lst)] = lst
        S = np.zeros((n, R))
        for j, (_, products) in enumerate(reactions):
            for i in lists[j]:
                S[i, j] -= 1.0
            for s, y in products:
                if s in third_bodies:
                    continue
                S[index[s], j] += y
        self.S = S
        self.ST = np.ascontiguousarray(S.T)

        entries: dict[tuple[int, int], int] = {}
        rows, cols = [], []

        def entry(i, k):
            key = (i, k)
            if key not in entries:
                entries[key] = len(rows)
                rows.append(i)
                cols.append(k)
            return entries[key]

        for i in range(n):
            entry(i, i)
        jp, ent, coef = [], [], []
        for j in range(R):
            nz = np.nonzero(S[:, j])[0]
            for p in range(P):
                kk = self.react[j, p]
                if kk == n:
                    continue
                for i in nz:
                    jp.append(j * P + p)
                    ent.append(entry(int(i), int(kk)))
                    coef.append(S[i, j])
        self.pattern = SparsePattern(n, np.array(rows), np.array(cols))
        G = sp.csr_matrix((coef, (jp, ent)), shape=(R * P, len(rows)))
        self.GT = G.T.tocsr()

    def _gather(self, Y):
        C = Y.shape[0]
        Yp = np.concatenate([Y, np.ones((C, 1))], axis=1)
        return Yp[:, self.react]                                  # (C, R, P)

    def rates(self, k, Y):
        return k * self._gather(Y).prod(axis=2)

    def rhs(self, k, Y):
        return self.rates(k, Y) @ self.ST

    def jac_values(self, k, Y):
        X = self._gather(Y)
        C = X.shape[0]
        d = np.empty_like(X)
        for p in range(self.P):
            prod = k.copy() if k.ndim == 2 else np.broadcast_to(k, X.shape[:2]).copy()
            for q in range(self.P):
                if q != p:
                    prod *= X[:, :, q]
            d[:, :, p] = prod
        return np.asarray((self.GT @ d.reshape(C, -1).T).T)    # (C, nnz)

    def dense_jacobian(self, vals):
        C = vals.shape[0]
        J = np.zeros((C, self.n, self.n))
        J[:, self.pattern.rows, self.pattern.cols] = vals
        return J


class ChemSystem(System):
    """A mechanism with rate constants k of shape (C, R); mass[i] = 0 holds
    species i in quasi-steady state (its row becomes the constraint F_i = 0)."""

    def __init__(self, mech: Mechanism, k: np.ndarray, mass=None, dense=None):
        self.mech = mech
        self.k = np.atleast_2d(k)
        self.n = mech.n
        self.names = mech.species
        self.mass = np.ones(self.n) if mass is None else np.asarray(mass, dtype=float)
        self.dense = (self.n <= 60) if dense is None else dense

    def rhs(self, Y):
        return self.mech.rhs(self.k, Y)

    def jac(self, Y):
        vals = self.mech.jac_values(self.k, Y)
        return self.mech.dense_jacobian(vals) if self.dense else vals

    def factor(self, alpha, J, algebraic_rows_only=False):
        if self.dense:
            return DenseLU(dense_matrix(alpha, self.mass, J, algebraic_rows_only))
        vals = self.mech.pattern.matrix_values(alpha, self.mass, J, algebraic_rows_only)
        return SparseLU(self.mech.pattern, vals)
