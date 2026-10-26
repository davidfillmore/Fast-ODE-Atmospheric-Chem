"""The three solvers of the paper, each called once per split interval.

- RestartingLoop: the baseline Rosenbrock loop, line for line as the baseline
  implementation (MICM at commit e15898ce): a restart at 1e-6 of the interval,
  the elementary controller, KPP's constants, the root-mean-square norm over the
  group, and, with constraints, the step change in place of the algebraic rows of
  the error estimate.
- Algorithm: the same stages and norm, with the step carried across intervals
  and quantized to them, the H211b filter with its memory carried (the
  elementary factor after a rejection), the last step of an interval landing on
  its end exactly (a last step shortened to fit is kept from the controller, and
  the next interval starts from the proposal it was shortened from), and the
  embedded estimate kept on the algebraic rows.
- BackwardEuler: the baseline backward Euler with full Newton, its convergence
  test on the correction, and CAM-Chem's step reductions.

The prototypes keep the Jacobian after a rejection (the baseline's in-place LU
re-evaluates it); the count of Jacobians then compares the methods, not the
storage scheme.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .counts import Counts
from .tableaus import ROS3, Tableau

ROUND_OFF = float(np.finfo(float).eps)
H_MIN_FRACTION = 1.0e-15
H_START_FRACTION = 1.0e-6
ERROR_MIN = 1.0e-10


def error_norm(Y, Ynew, E, atol, rtol):
    """The baseline's normalized error: the root mean square over all variables of
    all cells of e / (atol + rtol max(|y|, |y_new|)), floored at 1e-10."""
    ymax = np.maximum(np.abs(Y), np.abs(Ynew))
    s = E / (atol[None, :] + rtol * ymax)
    value = math.sqrt(float(np.sum(s * s)) / s.size)
    return max(value, ERROR_MIN) if not math.isnan(value) else value


def quantize(h, interval):
    """Round h down to a divisor of the interval (eq. ros-quantize)."""
    if h >= interval:
        return interval
    n = math.ceil(interval / h - 1.0e-9)
    return interval / max(n, 1)


def stages(system, tab: Tableau, Y, F0, LU, H, counts: Counts):
    """One Rosenbrock attempt in the transformed form with the mass matrix."""
    S = tab.stages
    new = tab.new_forcing
    mass = system.mass
    dae = system.has_constraints
    K = []
    forcing = None
    for r in range(S):
        if r == 0:
            forcing = F0
        elif new[r]:
            Yr = Y.copy()
            for s in range(r):
                if tab.a[r, s] != 0.0:
                    Yr += tab.a[r, s] * K[s]
            forcing = system.rhs(Yr)
            counts.rhs += 1
        rhs = forcing.copy()
        for s in range(r):
            coef = tab.c[r, s] / H
            if coef != 0.0:
                rhs += coef * (K[s] * mass[None, :] if dae else K[s])
        K.append(LU.solve(rhs))
        counts.solves += 1
    Ynew = Y.copy()
    E = np.zeros_like(Y)
    for r in range(S):
        if tab.m[r] != 0.0:
            Ynew += tab.m[r] * K[r]
        if tab.e[r] != 0.0:
            E += tab.e[r] * K[r]
    return Ynew, E


def clamp(system, Y):
    """Negative values to zero after a call, the differential variables only when
    there are constraints (the baseline's post-solve clamp)."""
    if system.has_constraints:
        diff = system.mass == 1.0
        Y[:, diff] = np.maximum(Y[:, diff], 0.0)
    else:
        np.maximum(Y, 0.0, out=Y)
    return Y


def project(system, Y, atol, rtol, counts: Counts, max_updates=10, tol=0.1,
            max_backtracks=24, backtrack_factor=0.5, sigma=1.0e-4):
    """Consistent initialization (Appendix F): damped Newton on G(x0, z) = 0 for the
    algebraic variables, accepted on the weighted Newton correction. Returns
    (Y, status); a failure restores Y."""
    alg = system.mass == 0.0
    if not np.any(alg):
        return Y, "converged"
    atol_a = atol[alg]

    def residual(Z):
        G = system.rhs(Z)
        counts.rhs += 1
        G[:, ~alg] = 0.0
        return G

    def finite(A):
        return bool(np.all(np.isfinite(A[:, alg])))

    def weighted(Z, D):
        z, d = Z[:, alg], D[:, alg]
        scale = atol_a[None, :] + rtol * np.maximum(np.abs(z), np.abs(z + d))
        with np.errstate(divide="ignore", invalid="ignore"):
            q = np.where(scale > 0.0, np.abs(d) / scale, np.where(d == 0.0, 0.0, np.inf))
        return float(q.max())

    original = None
    for update in range(max_updates + 1):
        G = residual(Y)
        if not finite(G):
            break
        counts.proj_passes += 1
        if float(np.abs(G[:, alg]).max()) == 0.0:
            return Y, "converged"
        J = system.jac(Y)
        counts.jac += 1
        LU = system.factor(1.0, J, algebraic_rows_only=True)
        counts.lu += 1
        delta = LU.solve(G)
        counts.solves += 1
        if not finite(delta):
            break
        q = weighted(Y, delta)
        if q <= tol:
            return Y, "converged"
        if update == max_updates:
            break
        if original is None:
            original = Y.copy()
        step = 1.0
        accepted = None
        for backtrack in range(max_backtracks + 1):
            candidate = Y.copy()
            candidate[:, alg] += step * delta[:, alg]
            step_used = step
            step *= backtrack_factor
            if not finite(candidate):
                continue
            Gc = residual(candidate)
            counts.proj_trials += 1
            if not finite(Gc):
                continue
            if float(np.abs(Gc[:, alg]).max()) == 0.0:
                accepted = candidate
                break
            trial = LU.solve(Gc)
            counts.solves += 1
            if not finite(trial):
                continue
            qt = weighted(candidate, trial)
            if qt <= tol or qt < (1.0 - sigma * step_used) * q:
                accepted = candidate
                break
        if accepted is None:
            break
        Y = accepted
    if original is not None:
        Y = original
    return Y, "projection failed"


@dataclass
class Memory:
    """What the algorithm carries from one interval to the next."""
    h_next: float = 0.0
    err_prev: float = 0.0
    h_prev: float = 0.0


class RestartingLoop:
    """The baseline Rosenbrock loop."""

    name = "restarting"

    def __init__(self, tableau: Tableau = ROS3, max_steps: int = 1000, step_change_rows: bool = True):
        """step_change_rows: replace the algebraic rows of the error estimate by the
        step change, as the baseline does; False keeps the embedded estimate (a
        control that separates this choice from the restart)."""
        self.tableau = tableau
        self.max_steps = max_steps
        self.step_change_rows = step_change_rows

    def new_memory(self):
        return None

    def solve(self, system, Y, dt, atol, rtol, counts: Counts, memory=None):
        tab = self.tableau
        counts.calls += 1
        h_min = H_MIN_FRACTION * dt
        h_max = dt
        H = min(max(h_min, H_START_FRACTION * dt), h_max)
        Y = Y.copy()
        if system.has_constraints:
            Y, status = project(system, Y, atol, rtol, counts)
            if status != "converged":
                return Y, status
        alg = system.mass == 0.0
        t = 0.0
        reject_last = reject_more = False
        nsteps = 0
        status = "running"
        while (t - dt + ROUND_OFF) <= 0.0 and status == "running":
            if nsteps > self.max_steps:
                status = "too many steps"
                break
            if (t + 0.1 * H == t) or H <= min(ROUND_OFF, h_min):
                status = "step too small"
                break
            H = min(H, abs(dt - t))
            F0 = system.rhs(Y)
            counts.rhs += 1
            J = system.jac(Y)
            counts.jac += 1
            accepted = False
            while not accepted:
                LU = system.factor(1.0 / (H * tab.gamma), J)
                counts.lu += 1
                Ynew, E = stages(system, tab, Y, F0, LU, H, counts)
                if system.has_constraints and self.step_change_rows:
                    E[:, alg] = Ynew[:, alg] - Y[:, alg]
                err = error_norm(Y, Ynew, E, atol, rtol)
                nsteps += 1
                counts.steps += 1
                if math.isnan(err) or math.isinf(err):
                    status = "nan or inf"
                    break
                fac = min(6.0, max(0.2, 0.9 / err ** (1.0 / tab.k)))
                Hnew = H * fac
                if err < 1.0 or H < h_min:
                    counts.accepted += 1
                    t += H
                    Y = Ynew
                    Hnew = max(h_min, min(Hnew, h_max))
                    if reject_last:
                        Hnew = min(Hnew, H)
                    reject_last = reject_more = False
                    H = Hnew
                    accepted = True
                else:
                    if reject_more:
                        Hnew = H * 0.1
                    reject_more = reject_last
                    reject_last = True
                    H = Hnew
                    counts.rejected += 1
        if status == "running":
            status = "converged"
        return clamp(system, Y), status


class Algorithm:
    """The integrator of Section 4."""

    name = "algorithm"

    def __init__(self, tableau: Tableau = ROS3, controller: str = "h211b", b: float = 4.0,
                 safety: float = 0.9, fac_min: float = 0.2, fac_max: float = 6.0,
                 carry: bool = True, quantize: bool = True, max_steps: int = 1000,
                 carry_memory: bool = True):
        """carry: start each interval from the step proposed at the end of the last
        (else from 1e-6 of the interval); carry_memory: keep the filter's memory across
        the boundary (else the elementary rule until a step of the interval is accepted);
        quantize: round every proposal down to a divisor of the interval."""
        self.tableau = tableau
        self.controller = controller
        self.b = b
        self.safety = safety
        self.fac_min = fac_min
        self.fac_max = fac_max
        self.carry = carry
        self.quantize = quantize
        self.max_steps = max_steps
        self.carry_memory = carry_memory

    def new_memory(self):
        return Memory()

    def factor(self, err, H, mem: Memory):
        k = self.tableau.k
        fac = self.safety / err ** (1.0 / k)
        if self.controller == "h211b" and mem.err_prev > 0.0 and mem.h_prev > 0.0:
            fac = (self.safety * (err * mem.err_prev) ** (-1.0 / (self.b * k))
                   * (H / mem.h_prev) ** (-1.0 / self.b))
        return min(self.fac_max, max(self.fac_min, fac))

    def solve(self, system, Y, dt, atol, rtol, counts: Counts, memory: Memory = None):
        tab = self.tableau
        mem = memory if memory is not None else Memory()
        if not (self.carry and self.carry_memory):
            mem.err_prev = mem.h_prev = 0.0
        counts.calls += 1
        h_min = H_MIN_FRACTION * dt
        h_max = dt
        if self.carry and mem.h_next > 0.0:
            H = mem.h_next
        else:
            H = H_START_FRACTION * dt
        H = min(max(h_min, H), h_max)
        if self.quantize:
            H = quantize(H, dt)
        proposed = H
        Y = Y.copy()
        if system.has_constraints:
            Y, status = project(system, Y, atol, rtol, counts)
            if status != "converged":
                return Y, status
        t = 0.0
        reject_last = reject_more = False
        nsteps = 0
        status = "running"
        while (t - dt + ROUND_OFF) <= 0.0 and status == "running":
            if nsteps >= self.max_steps:
                status = "too many steps"
                break
            if (t + 0.1 * H == t) or H < h_min:
                status = "step too small"
                break
            # The last step lands on the end of the interval exactly. n quantized steps
            # of dt/n can fall a few ulps short of dt, and a micro-step taken there would
            # hand its own tiny proposal to the next interval.
            last = dt - t <= H * (1.0 + 1.0e-8)
            shortened = last and dt - t < H * (1.0 - 1.0e-8)
            H_ctrl = H
            if last:
                H = dt - t
            F0 = system.rhs(Y)
            counts.rhs += 1
            J = system.jac(Y)
            counts.jac += 1
            accepted = False
            while not accepted:
                # A retry is a new attempt too. Neither the step floor nor the
                # attempt budget may be bypassed by the inner rejection loop.
                if nsteps >= self.max_steps:
                    status = "too many steps"
                    break
                if H < h_min or t + 0.1 * H == t:
                    status = "step too small"
                    break
                LU = system.factor(1.0 / (H * tab.gamma), J)
                counts.lu += 1
                Ynew, E = stages(system, tab, Y, F0, LU, H, counts)
                err = error_norm(Y, Ynew, E, atol, rtol)
                nsteps += 1
                counts.steps += 1
                if math.isnan(err) or math.isinf(err):
                    status = "nan or inf"
                    break
                if err < 1.0 and shortened:
                    # A last step shortened to fit the interval is kept from the controller:
                    # its estimate stays out of the memory, and the next interval starts
                    # from the proposal it was shortened from.
                    counts.accepted += 1
                    t = dt
                    Y = Ynew
                    proposed = H_ctrl
                    accepted = True
                elif err < 1.0:
                    fac = self.factor(err, H, mem)
                    counts.accepted += 1
                    t = dt if last else t + H
                    Y = Ynew
                    mem.err_prev = err
                    mem.h_prev = H
                    Hnew = H * fac
                    if reject_last:
                        Hnew = min(Hnew, H)
                    Hnew = max(h_min, min(Hnew, h_max))
                    if self.quantize:
                        Hnew = quantize(Hnew, dt)
                    reject_last = reject_more = False
                    proposed = Hnew
                    H = Hnew
                    accepted = True
                else:
                    # The filter weighs the last accepted estimate too, and right after a
                    # rejection its factor can be near one or above; the elementary factor
                    # is below the safety factor whenever err >= 1.
                    Hnew = H * max(self.fac_min, self.safety / err ** (1.0 / tab.k))
                    if reject_more:
                        Hnew = H * 0.1
                    reject_more = reject_last
                    reject_last = True
                    if self.quantize:
                        Hnew = quantize(Hnew, dt)
                    H = Hnew
                    last = dt - t <= H * (1.0 + 1.0e-8)
                    shortened = last and dt - t < H * (1.0 - 1.0e-8)
                    H_ctrl = H
                    if last:
                        H = dt - t
                    counts.rejected += 1
        mem.h_next = proposed
        if status == "running":
            status = "converged"
        return clamp(system, Y), status


class BackwardEuler:
    """The baseline backward Euler with Newton's iteration."""

    name = "backward-euler"

    def __init__(self, max_iterations: int = 11, reductions=(0.5, 0.5, 0.5, 0.5, 0.1),
                 small: float = 1.0e-40):
        self.max_iterations = max_iterations
        self.reductions = reductions
        self.small = small

    def new_memory(self):
        return None

    def converged(self, delta, Y, atol, rtol):
        if not np.all(np.isfinite(delta)):
            return False
        a = np.abs(delta)
        bound = np.maximum(np.maximum(self.small, atol[None, :]), rtol * np.abs(Y))
        return bool(np.all(a <= bound))

    def solve(self, system, Y, dt, atol, rtol, counts: Counts, memory=None):
        counts.calls += 1
        H = dt
        t = 0.0
        n_success = 0
        n_fail = 0
        Yn = Y.copy()
        Yn1 = Y.copy()
        status = "converged"
        while t < dt:
            ok = False
            iterations = 0
            while True:
                counts.newton += 1
                F = system.rhs(Yn1)
                counts.rhs += 1
                J = system.jac(Yn1)
                counts.jac += 1
                LU = system.factor(1.0 / H, J)
                counts.lu += 1
                delta = LU.solve(F - (Yn1 - Yn) / H)
                counts.solves += 1
                Yn1 = np.maximum(Yn1 + delta, 0.0)
                iterations += 1
                if iterations > 1:
                    ok = self.converged(delta, Yn1, atol, rtol)
                if ok or iterations >= self.max_iterations:
                    break
            counts.steps += 1
            if not ok:
                counts.rejected += 1
                counts.failures += 1
                n_success = 0
                if n_fail >= len(self.reductions):
                    t += H
                    counts.unconverged += 1
                    status = "accepted unconverged"
                    if not np.all(np.isfinite(Yn1)):
                        status = "nan or inf"
                    break
                Yn1 = Yn.copy()
                H *= self.reductions[n_fail]
                n_fail += 1
            else:
                counts.accepted += 1
                t += H
                Yn = Yn1.copy()
                n_success += 1
                if n_success >= 2:
                    n_success = 0
                    H *= 2.0
            H = min(H, dt - t)
        return clamp(system, Yn1), status


def fixed_steps(system, tab: Tableau, Y, t_end, n_steps, counts: Counts = None):
    """n_steps Rosenbrock steps of equal length, without error control, for
    convergence studies."""
    counts = counts or Counts()
    H = t_end / n_steps
    for _ in range(n_steps):
        F0 = system.rhs(Y)
        J = system.jac(Y)
        LU = system.factor(1.0 / (H * tab.gamma), J)
        Y, _ = stages(system, tab, Y, F0, LU, H, counts)
    return Y
