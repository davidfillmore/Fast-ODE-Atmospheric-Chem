"""Tests: tableau orders and index-1 conditions, convergence orders in use,
linear invariants, the projection's affine covariance, quantization, and the
step rules of the two Rosenbrock loops.

    python3 -m pytest tests
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fastode import Algorithm, Counts, ROS3, RODAS3, RestartingLoop, project  # noqa: E402
from fastode import problems as P  # noqa: E402
from fastode.solvers import fixed_steps, quantize  # noqa: E402
from fastode.system import DenseSystem  # noqa: E402


@pytest.mark.parametrize("tab", [ROS3, RODAS3])
def test_order_conditions(tab):
    A, Gamma, b, bhat = tab.standard_form()
    beta = A + Gamma
    ap = A.sum(axis=1)
    bp = beta.sum(axis=1)
    assert abs(b.sum() - 1) < 1e-12
    assert abs(b @ bp - 0.5) < 1e-12
    assert abs(b @ ap ** 2 - 1 / 3) < 1e-12
    assert abs(b @ (beta @ bp) - 1 / 6) < 1e-12
    assert abs(bhat.sum() - 1) < 1e-12 and abs(bhat @ bp - 0.5) < 1e-12


def test_index1_conditions():
    assert abs(RODAS3.index1(2) - 1) < 1e-12 and abs(RODAS3.index1(3) - 1) < 1e-12
    assert abs(ROS3.index1(2) - 1.0908767) < 1e-6


def test_forcing_reuse_flags():
    assert ROS3.new_forcing == [True, True, False]
    assert RODAS3.new_forcing == [True, False, True, True]


@pytest.mark.parametrize("tab", [ROS3, RODAS3])
def test_order_three_on_a_nonlinear_ode(tab):
    # y' = -y^3, y(0) = 1: y = 1/sqrt(1+2t). (RODAS3 integrates y' = -y^2 exactly.)
    sys_ = DenseSystem(1, lambda Y: -Y ** 3, lambda Y: (-3 * Y ** 2)[:, :, None])
    exact = 1.0 / math.sqrt(3.0)
    errs = [abs(fixed_steps(sys_, tab, np.array([[1.0]]), 1.0, n)[0, 0] - exact) for n in (20, 40, 80)]
    orders = [math.log2(errs[i] / errs[i + 1]) for i in range(2)]
    assert min(orders) > 2.8


def test_index1_orders():
    system, Y0 = P.manufactured_coupled()
    from scipy.integrate import solve_ivp
    s = solve_ivp(lambda t, y: [-y[0] + y[0] ** 2 + y[1], y[0] * (y[0] ** 2 + y[1]) - y[1]], (0, 1), [0.5, 0.5],
                  method="DOP853", rtol=1e-13, atol=1e-15)
    z_exact = s.y[0, -1] ** 2 + s.y[1, -1]

    def zerr(tab, n):
        return abs(fixed_steps(system, tab, Y0.copy(), 1.0, n)[0, 2] - z_exact)

    assert math.log2(zerr(RODAS3, 64) / zerr(RODAS3, 128)) > 2.7
    assert math.log2(zerr(ROS3, 64) / zerr(ROS3, 128)) < 2.2


def test_linear_invariant_without_clamp():
    system, Y0 = P.robertson()
    for tab in (ROS3, RODAS3):
        Y = fixed_steps(system, tab, Y0.copy(), 1.0e-2, 50)    # steps the initial transient can bear
        assert abs(Y.sum() - 1.0) < 1e-14


def test_projection_is_affine_covariant():
    system, Y0 = P.tropospheric(qssa=True)
    atol = np.full(system.n, 600.0)
    ref, st = project(system, Y0.copy(), atol, 1e-3, Counts())
    assert st == "converged"
    alg = system.mass == 0.0
    D = np.where(alg, 1e12, 1.0)
    scaled = DenseSystem(system.n, lambda Y: system.rhs(Y) * D, lambda Y: system.jac(Y) * D[None, :, None],
                         mass=system.mass)
    Ys, st2 = project(scaled, Y0.copy(), atol, 1e-3, Counts())
    assert st2 == "converged"
    assert np.allclose(Ys, ref, rtol=1e-8, atol=0)


def test_quantize():
    for h in (1.0, 7.3, 33.0, 600.0, 900.0):
        q = quantize(h, 600.0)
        n = round(600.0 / q)
        assert abs(n * q - 600.0) < 1e-9 and h / 2 < q <= h or h >= 600.0


def test_restart_costs_nine_steps():
    # a linear decay that one step could cover within the tolerance
    sys_ = DenseSystem(1, lambda Y: -1e-6 * Y, lambda Y: np.full((Y.shape[0], 1, 1), -1e-6))
    c = Counts()
    RestartingLoop(ROS3).solve(sys_, np.array([[1.0]]), 600.0, np.array([1e-10]), 1e-3, c)
    assert c.accepted >= 9
    c2 = Counts()
    alg = Algorithm(ROS3)
    mem = alg.new_memory()
    Y = np.array([[1.0]])
    for _ in range(5):
        Y, _ = alg.solve(sys_, Y, 600.0, np.array([1e-10]), 1e-3, c2, mem)
    c3 = Counts()
    alg.solve(sys_, Y, 600.0, np.array([1e-10]), 1e-3, c3, mem)
    assert c3.accepted == 1


def test_quantized_steps_land_on_the_interval_end():
    # seven steps of 600/7 sum to a few ulps short of 600; the call must still take
    # seven steps and carry 600/7, not a micro-step's proposal
    sys_ = DenseSystem(1, lambda Y: -1e-9 * Y, lambda Y: np.full((Y.shape[0], 1, 1), -1e-9))
    alg = Algorithm(ROS3, fac_max=1.0)
    mem = alg.new_memory()
    mem.h_next = 600.0 / 7
    Y = np.array([[1.0]])
    for _ in range(3):
        c = Counts()
        Y, st = alg.solve(sys_, Y, 600.0, np.array([1e-10]), 1e-3, c, mem)
        assert st == "converged" and c.steps == 7
        assert abs(mem.h_next - 600.0 / 7) < 1e-9


def test_a_rejection_shrinks_the_step():
    # a carried step far too long, with a filter memory that would ask it to grow
    sys_ = DenseSystem(1, lambda Y: -Y, lambda Y: np.full((Y.shape[0], 1, 1), -1.0))
    alg = Algorithm(ROS3)
    mem = alg.new_memory()
    mem.h_next, mem.err_prev, mem.h_prev = 1.0, 1e-8, 1.0
    c = Counts()
    alg.solve(sys_, np.array([[1.0]]), 1.0, np.array([1e-12]), 4e-3, c, mem)
    assert c.rejected == 1      # the filter's factor would have kept h = 1 and been rejected thrice


def test_a_shortened_last_step_is_kept_from_the_controller():
    # the constrained triples: the step changes inside an interval, the last step is
    # shortened to fit, and its small proposal must not be carried
    system, Y0 = P.triples(1, dae=True)
    alg = Algorithm(RODAS3)
    mem = alg.new_memory()
    Y, c = Y0.copy(), Counts()
    for _ in range(100):
        Y, st = alg.solve(system, Y, 0.2, np.full(3, 1e-18), 1e-4, c, mem)
        assert st == "converged"
    assert c.steps < 250        # 311 when the shortened step's proposal was carried


def test_algorithm_does_not_force_acceptance_below_minimum_step():
    # The first attempt already uses h_min. Its rejected retry must return an
    # error, even if a still smaller step could bypass the normal error test.
    system = DenseSystem(
        1, lambda Y: -1e15 * Y,
        lambda Y: np.full((Y.shape[0], 1, 1), -1e15),
    )
    alg = Algorithm(ROS3, max_steps=10)
    mem = alg.new_memory()
    mem.h_next = 1e-15
    initial, counts = np.array([[1.0]]), Counts()
    result, status = alg.solve(
        system, initial, 1.0, np.array([1e-30]), 1e-12, counts, mem,
    )
    assert status == "step too small"
    assert counts.steps == counts.rejected == 1
    assert counts.accepted == 0
    np.testing.assert_array_equal(result, initial)


@pytest.mark.parametrize("budget", [0, 1, 2])
def test_algorithm_attempt_budget_includes_rejected_retries(budget):
    system = DenseSystem(
        1, lambda Y: -Y, lambda Y: np.full((Y.shape[0], 1, 1), -1.0),
    )
    alg = Algorithm(ROS3, max_steps=budget)
    mem = alg.new_memory()
    mem.h_next = 1.0
    initial, counts = np.array([[1.0]]), Counts()
    result, status = alg.solve(
        system, initial, 1.0, np.array([1e-30]), 1e-13, counts, mem,
    )
    assert status == "too many steps"
    assert counts.steps == counts.rejected == budget
    assert counts.accepted == 0
    np.testing.assert_array_equal(result, initial)


@pytest.mark.parametrize(
    "budget, expected_status", [(2, "too many steps"), (10, "converged")],
)
def test_algorithm_attempt_budget_includes_accepted_steps(budget, expected_status):
    system = DenseSystem(
        1, lambda Y: -1e-6 * Y,
        lambda Y: np.full((Y.shape[0], 1, 1), -1e-6),
    )
    alg = Algorithm(ROS3, max_steps=budget, fac_max=1.0)
    mem = alg.new_memory()
    mem.h_next = 0.1
    counts = Counts()
    _, status = alg.solve(
        system, np.array([[1.0]]), 1.0, np.array([1e-10]), 1e-3,
        counts, mem,
    )
    assert status == expected_status
    assert counts.steps == counts.accepted == budget
    assert counts.rejected == 0
