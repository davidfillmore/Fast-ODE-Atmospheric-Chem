"""Factorizations of the stage matrix alpha*M - J for a group of C cells.

Dense: one LU per cell (small systems). Sparse: the C cells as one
block-diagonal matrix, factored by SuperLU (TS1). The operation counts do not
depend on which is used.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class DenseLU:
    def __init__(self, A: np.ndarray):
        # A: (C, n, n)
        self._lus = [sla.lu_factor(A[c], check_finite=False) for c in range(A.shape[0])]

    def solve(self, B: np.ndarray) -> np.ndarray:
        return np.stack([sla.lu_solve(lu, B[c], check_finite=False) for c, lu in enumerate(self._lus)])


def dense_matrix(alpha: float, mass: np.ndarray, J: np.ndarray, algebraic_rows_only: bool = False):
    """alpha*diag(mass) - J, or, for the projection, diag(mass) - J with only the
    algebraic rows of J kept (alpha = 1)."""
    if algebraic_rows_only:
        J = J * (mass == 0.0)[None, :, None]
    A = -J.copy()
    idx = np.arange(J.shape[1])
    A[:, idx, idx] += alpha * mass
    return A


class SparsePattern:
    """A fixed sparsity pattern (rows, cols) that includes the diagonal."""

    def __init__(self, n: int, rows: np.ndarray, cols: np.ndarray):
        self.n = n
        self.rows = np.asarray(rows, dtype=np.int64)
        self.cols = np.asarray(cols, dtype=np.int64)
        self.nnz = len(self.rows)
        pos = {(int(r), int(c)): k for k, (r, c) in enumerate(zip(self.rows, self.cols))}
        self.diag = np.array([pos[(i, i)] for i in range(n)], dtype=np.int64)
        # column-major order of one block, for CSC
        self._order = np.lexsort((self.rows, self.cols))
        self._block_indices = self.rows[self._order]
        counts = np.bincount(self.cols, minlength=n)
        self._block_indptr = np.concatenate([[0], np.cumsum(counts)])
        self._csc_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    def csc_structure(self, C: int):
        if C not in self._csc_cache:
            n, nnz = self.n, self.nnz
            indices = np.concatenate([self._block_indices + c * n for c in range(C)])
            indptr = np.concatenate([self._block_indptr[:-1] + c * nnz for c in range(C)] + [[C * nnz]])
            self._csc_cache[C] = (indices, indptr)
        return self._csc_cache[C]

    def matrix_values(self, alpha: float, mass: np.ndarray, Jvals: np.ndarray, algebraic_rows_only: bool = False):
        vals = -Jvals
        if algebraic_rows_only:
            vals = vals * (mass[self.rows] == 0.0)[None, :]
        else:
            vals = vals.copy()
        vals[:, self.diag] += alpha * mass[None, :]
        return vals


class SparseLU:
    def __init__(self, pattern: SparsePattern, vals: np.ndarray):
        C = vals.shape[0]
        indices, indptr = pattern.csc_structure(C)
        data = vals[:, pattern._order].ravel()
        n = pattern.n * C
        A = sp.csc_matrix((data, indices, indptr), shape=(n, n))
        self._lu = spla.splu(A, permc_spec="MMD_AT_PLUS_A")
        self._shape = (C, pattern.n)

    def solve(self, B: np.ndarray) -> np.ndarray:
        return self._lu.solve(np.ascontiguousarray(B).ravel()).reshape(self._shape)
