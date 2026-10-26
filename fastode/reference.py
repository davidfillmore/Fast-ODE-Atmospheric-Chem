"""References and scoring for the small problems.

References come from SciPy's Radau method at a tight tolerance, interval by
interval, so that piecewise-constant coefficients change exactly at the interval
boundaries; Radau is an implicit Runge-Kutta method, independent of the
Rosenbrock methods under test.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp


def radau_reference(system, Y0, dt, intervals, rtol=1e-12, atol=1e-30, on_interval=None):
    """States at the ends of `intervals` intervals of length dt (one cell, ODE only)."""
    y = Y0[0].astype(float).copy()
    out = []
    for m in range(intervals):
        if on_interval is not None:
            on_interval(system, m)

        def f(t, v):
            return system.rhs(v[None, :])[0]

        def jac(t, v):
            J = system.jac(v[None, :])
            if J.ndim == 2:      # sparse values
                J = system.mech.dense_jacobian(J)
            return J[0]

        sol = solve_ivp(f, (0.0, dt), y, method="Radau", jac=jac, rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        out.append(y.copy())
    return np.array(out)


def digits(result, reference, threshold, columns=None):
    """Sandu digits overall and per variable over all samples (rows) whose
    reference value is at least the threshold; capped at 16."""
    R = np.asarray(reference)
    X = np.asarray(result)
    cols = range(R.shape[1]) if columns is None else columns
    per = {}
    all_sq = []
    for j in cols:
        mask = R[:, j] >= threshold
        if not np.any(mask):
            continue
        sq = ((X[mask, j] - R[mask, j]) / R[mask, j]) ** 2
        all_sq.extend(sq.tolist())
        ms = float(np.mean(sq))
        per[j] = 16.0 if ms <= 0.0 else min(16.0, -math.log10(math.sqrt(ms)))
    ms = float(np.mean(all_sq)) if all_sq else float("nan")
    overall = 16.0 if ms <= 0.0 else min(16.0, -math.log10(math.sqrt(ms)))
    return overall, per
