"""Experiments II: constrained species.

    python3 experiments/dae_experiments.py [orders|pr|robertson|chapman|tropo|triples|vdp|akzo|
                                            projection|schur|conservation|rows|qde|index2|all]

Writes CSV files to results/dae/ and prints a summary of each.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fastode import Algorithm, Counts, ROS3, RODAS3, RestartingLoop, project  # noqa: E402
from fastode import problems as P  # noqa: E402
from fastode.reference import digits, radau_reference  # noqa: E402
from fastode.solvers import fixed_steps  # noqa: E402

OUT = ROOT / "results/dae"
OUT.mkdir(parents=True, exist_ok=True)


def write(name, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT / f"{name}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def slopes(hs, errs):
    return [math.log(errs[i] / errs[i + 1]) / math.log(hs[i] / hs[i + 1]) for i in range(len(hs) - 1)]


def run_intervals(system, Y0, solver, dt, intervals, rtol, atol, on_interval=None):
    Y = Y0.copy()
    counts = Counts()
    mem = solver.new_memory()
    states, bad = [], 0
    for m in range(intervals):
        if on_interval is not None:
            on_interval(system, m)
        Y, status = solver.solve(system, Y, dt, atol, rtol, counts, mem)
        bad += status != "converged"
        states.append(Y[0].copy())
    return np.array(states), counts, bad


# ---- (a) orders on the manufactured index-1 problems ---------------------------
def orders():
    rows = []
    T = 1.0
    exact = {"circle": np.array([math.sin(T), math.cos(T)])}
    from scipy.integrate import solve_ivp

    def red(t, y):
        z = y[0] ** 2 + y[1]
        return [-y[0] + z, y[0] * z - y[1]]

    s = solve_ivp(red, (0, T), [0.5, 0.5], method="DOP853", rtol=1e-13, atol=1e-15)
    y = s.y[:, -1]
    exact["coupled"] = np.array([y[0], y[1], y[0] ** 2 + y[1]])
    for pname, build in (("circle", P.manufactured_circle), ("coupled", P.manufactured_coupled)):
        system, Y0 = build()
        alg = system.mass == 0.0
        for tab in (ROS3, RODAS3):
            Ns = [8, 16, 32, 64, 128, 256]
            ed, ea = [], []
            for N in Ns:
                Y = fixed_steps(system, tab, Y0.copy(), T, N)[0]
                e = np.abs(Y - exact[pname])
                ed.append(e[~alg].max())
                ea.append(e[alg].max())
                rows.append({"problem": pname, "tableau": tab.name, "N": N, "err_differential": e[~alg].max(),
                             "err_algebraic": e[alg].max()})
            hs = [T / N for N in Ns]
            print(f"orders {pname:8s} {tab.name:6s} differential {np.round(slopes(hs, ed), 2)} "
                  f"algebraic {np.round(slopes(hs, ea), 2)}")
    write("orders", rows)


# ---- (b) stiff order on Prothero-Robinson ----------------------------------------
def prothero_robinson():
    rows = []
    T = 1.0
    for lam in (-1.0e2, -1.0e4, -1.0e6, -1.0e8):
        system, Y0 = P.prothero_robinson(lam)
        for tab in (ROS3, RODAS3):
            Ns = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
            errs = []
            for N in Ns:
                Y = fixed_steps(system, tab, Y0.copy(), T, N)[0]
                errs.append(abs(Y[1] - math.sin(T)))
                rows.append({"lambda": lam, "tableau": tab.name, "N": N, "h": T / N, "error": errs[-1]})
            print(f"Prothero-Robinson {tab.name:6s} errors {[f'{e:.1e}' for e in errs]}")
            print(f"                  slopes {np.round(slopes([T / n for n in Ns], errs), 2)}")
    write("prothero_robinson", rows)


# ---- (c) Robertson: the QSSA-DAE against the full ODE ------------------------------
def robertson():
    rows = []
    dt, intervals, a = 1.0e3, 100, 1e-12
    for k2 in (3e5, 3e6, 3e7, 3e8, 3e9):
        full, Y0 = P.robertson(k2)
        ref = radau_reference(full, Y0, dt, intervals)
        for rtol in ((1e-2, 1e-3, 1e-4, 1e-5, 1e-6) if k2 == 3e7 else (1e-4,)):
            for form, tab in (("ode", ROS3), ("ode", RODAS3), ("dae", RODAS3)):
                for sname, solver in (("algorithm", Algorithm(tab)), ("restarting", RestartingLoop(tab))):
                    system, Y0 = P.robertson(k2, qssa=(form == "dae"))
                    if form == "dae":
                        Y0[0, 1] = P.robertson_manifold(Y0, k2)[0]
                    X, c, bad = run_intervals(system, Y0, solver, dt, intervals, rtol, np.full(3, rtol * a))
                    D, per = digits(X, ref, a, columns=[0, 2])
                    rows.append({"k2": k2, "rtol": rtol, "form": form, "tableau": tab.name, "solver": sname,
                                 "digits_AC": round(D, 3), "digits_B": round(digits(X, ref, a, columns=[1])[0], 3),
                                 "steps": c.steps, "rejected": c.rejected, "lu": c.lu, "solves": c.solves,
                                 "rhs": c.rhs, "proj_passes": c.proj_passes, "failed_calls": bad})
                    print(f"robertson k2={k2:g} rtol={rtol:g} {form} {tab.name:6s} {sname:10s}: "
                          f"A,C {D:5.2f} B {rows[-1]['digits_B']:5.2f} steps {c.steps:5d} failed {bad}")
    write("robertson", rows)


# ---- (d) Chapman ------------------------------------------------------------------
def chapman():
    rows = []
    dt, intervals, a = 1.0e5, 100, 1e-12
    for scale in (1.0, 10.0, 100.0, 1000.0):
        full, Y0 = P.chapman(scale)
        ref = radau_reference(full, Y0, dt, intervals)
        for rtol in (1e-2, 1e-4, 1e-6):
            for form, tab in (("ode", RODAS3), ("dae", RODAS3)):
                system, Y0 = P.chapman(scale, qssa=(form == "dae"))
                X, c, bad = run_intervals(system, Y0, Algorithm(tab), dt, intervals, rtol, np.full(3, rtol * a))
                D, per = digits(X, ref, a)
                rows.append({"scale": scale, "rtol": rtol, "form": form, "digits": round(D, 3),
                             "digits_O3": round(per.get(2, float("nan")), 3), "digits_O": round(per.get(0, float("nan")), 3),
                             "steps": c.steps, "lu": c.lu, "failed_calls": bad})
                print(f"chapman scale={scale:g} rtol={rtol:g} {form}: {D:5.2f} (O {rows[-1]['digits_O']}) "
                      f"steps {c.steps}")
    write("chapman", rows)


# ---- (e) the tropospheric mechanism ------------------------------------------------
def tropospheric():
    rows = []
    dt, intervals = 600.0, 144
    a = 1e-14 * 6.02214076e17
    for j in (0.1, 0.3, 1.0, 3.0, 10.0):
        full, Y0 = P.tropospheric(j)
        ref = radau_reference(full, Y0, dt, intervals)
        for rtol in ((1e-2, 1e-3, 1e-4, 1e-5, 1e-6) if j == 1.0 else (1e-4,)):
            for form in ("ode", "dae"):
                for sname, mk in (("algorithm", Algorithm), ("restarting", RestartingLoop)):
                    system, Y0 = P.tropospheric(j, qssa=(form == "dae"))
                    X, c, bad = run_intervals(system, Y0, mk(RODAS3), dt, intervals, rtol,
                                              np.full(system.n, rtol * a))
                    D, per = digits(X, ref, a)
                    rad = [P.TROPO_SPECIES.index(s) for s in P.TROPO_RADICALS]
                    Dr, _ = digits(X, ref, a, columns=rad)
                    Ds, _ = digits(X, ref, a, columns=[i for i in range(system.n) if i not in rad])
                    rows.append({"j": j, "rtol": rtol, "form": form, "solver": sname, "digits": round(D, 3),
                                 "digits_radicals": round(Dr, 3), "digits_slow": round(Ds, 3), "steps": c.steps,
                                 "rejected": c.rejected, "lu": c.lu, "proj_passes": c.proj_passes,
                                 "failed_calls": bad})
                    print(f"tropo j={j:g} rtol={rtol:g} {form} {sname:10s}: {D:5.2f} (radicals {Dr:5.2f}, "
                          f"slow {Ds:5.2f}) steps {c.steps} failed {bad}")
    write("tropospheric", rows)


# ---- (f) triples ---------------------------------------------------------------------
def triples():
    rows = []
    dt, intervals, a = 0.2, 100, 1e-12
    for N in (1, 10, 100):
        full, Y0 = P.triples(N)
        ref = radau_reference(full, Y0, dt, intervals)
        for form in ("ode", "dae"):
            system, Y0 = P.triples(N, dae=(form == "dae"))
            tab = RODAS3
            X, c, bad = run_intervals(system, Y0, Algorithm(tab), dt, intervals, 1e-4, np.full(system.n, 1e-4 * a))
            D, per = digits(X[1:], ref[1:], a)       # past the initial layer of the full system
            rows.append({"N": N, "form": form, "digits": round(D, 3), "steps": c.steps, "lu": c.lu,
                         "unknowns": system.n, "failed_calls": bad})
            print(f"triples N={N} {form}: {D:5.2f} digits, steps {c.steps}, unknowns {system.n}")
    write("triples", rows)


# ---- (g) Van der Pol: the singular-perturbation limit ------------------------------
def van_der_pol():
    rows = []
    T = 0.5
    dae, Y0 = P.van_der_pol(0.0, dae=True)
    Yd, c, _ = run_intervals(dae, Y0, Algorithm(RODAS3, max_steps=10**7), T, 1, 1e-10, np.full(2, 1e-14))
    for k in range(0, 7):
        eps = 10.0 ** (-k)
        ode, Y0 = P.van_der_pol(eps)
        Yo, co, _ = run_intervals(ode, Y0, Algorithm(RODAS3, max_steps=10**7), T, 1, 1e-10, np.full(2, 1e-14))
        diff = abs(Yo[-1, 0] - Yd[-1, 0])
        rows.append({"eps": eps, "X_ode": Yo[-1, 0], "X_dae": Yd[-1, 0], "difference": diff, "ode_steps": co.steps,
                     "dae_steps": c.steps})
    for i in range(len(rows) - 1):
        rows[i + 1]["order"] = round(math.log10(rows[i]["difference"] / rows[i + 1]["difference"]), 4)
    for r in rows:
        print(f"vdp eps={r['eps']:.0e}: |X_ode - X_dae| = {r['difference']:.3e} order {r.get('order', '')} "
              f"steps ode {r['ode_steps']} dae {r['dae_steps']}")
    write("van_der_pol", rows)


# ---- (h) Akzo Nobel: the DAE against its exact reduction -----------------------------
def akzo():
    rows = []
    T, intervals = 180.0, 18
    red, Y0r = P.akzo(dae=False)
    ref = radau_reference(red, Y0r, T / intervals, intervals)
    ref6 = np.concatenate([ref, (P.AKZO["Ks"] * ref[:, 0] * ref[:, 3])[:, None]], axis=1)
    print(f"akzo reference y1(180) = {ref[-1, 0]:.12f} (test set 0.1150794920661702)")
    for rtol in (1e-2, 1e-4, 1e-6):
        for form in ("reduced", "dae"):
            system, Y0 = P.akzo(dae=(form == "dae"))
            tab = RODAS3
            X, c, bad = run_intervals(system, Y0, Algorithm(tab), T / intervals, intervals, rtol,
                                      np.full(system.n, rtol * 1e-12))
            D, _ = digits(X, ref6 if form == "dae" else ref, 1e-12)
            rows.append({"rtol": rtol, "form": form, "digits": round(D, 3), "y1_180": X[-1, 0], "steps": c.steps,
                         "lu": c.lu, "failed_calls": bad})
            print(f"akzo rtol={rtol:g} {form:8s}: {D:5.2f} digits, y1(180) {X[-1, 0]:.10f}, steps {c.steps}")
    write("akzo", rows)


# ---- (i) the projection: weighted correction against the residual rule -------------
def projection():
    """Each case scales every constraint row by s, starts the algebraic variables
    from a guess, and runs (W) the projection of Appendix F and (R) Newton's
    iteration, undamped, stopped instead when the largest scaled residual is
    below 1e-8 (at most eleven iterations). A case is a false convergence when the rule stops with the state
    outside the integrator's tolerance of the true solution (weighted distance
    above one), and a false failure when the rule fails although an iterate
    came within a tenth of the tolerance."""
    from fastode.system import DenseSystem

    def scaled(system, s):
        alg = system.mass == 0.0
        D = np.where(alg, s, 1.0)

        def rhs(Y):
            return system.rhs(Y) * D[None, :]

        def jac(Y):
            J = system.jac(Y)
            if J.ndim == 2:
                J = system.mech.dense_jacobian(J)
            return J * D[None, :, None]

        return DenseSystem(system.n, rhs, jac, mass=system.mass, names=system.names)

    cases = []
    sysR, Y0 = P.robertson(qssa=True); Y0 = np.array([[0.9, 0.0, 0.1]]); cases.append(("robertson", sysR, Y0, 1e-3, 1e-12))
    sysC, Y0 = P.chapman(qssa=True); cases.append(("chapman", sysC, Y0, 1e-3, 1e-12))
    sysT, Y0 = P.tropospheric(qssa=True); cases.append(("tropospheric", sysT, Y0, 1e-3, 6.02e3 * 0.1))
    sysQ, Y0 = P.triples(N=3, dae=True); Y0 = np.tile([0.0, 0.0, 0.4], (1, 3)); cases.append(("triples", sysQ, Y0, 1e-3, 1e-12))

    def cold():
        def rhs(Y):
            return np.stack([np.zeros(Y.shape[0]), 2 * Y[:, 0] - Y[:, 1] ** 2], axis=1)

        def jac(Y):
            J = np.zeros((Y.shape[0], 2, 2))
            J[:, 1, 0] = 2.0
            J[:, 1, 1] = -2 * Y[:, 1]
            return J
        return DenseSystem(2, rhs, jac, mass=[1.0, 0.0], names=["B", "C"])

    cases.append(("cold-start", cold(), np.array([[0.5, 1e-6]]), 1e-3, 1e-12))
    rows = []
    for name, system, Y0, rtol, atol in cases:
        n = system.n
        atol_v = np.full(n, atol)
        alg = system.mass == 0.0
        # the true consistent state, by the projection at tight tolerance from a good start
        truth, st = project(system, Y0.copy() if name != "cold-start" else np.array([[0.5, 0.9]]), atol_v * 1e-6,
                            1e-12, Counts(), max_updates=60)
        if name == "tropospheric":
            truth, st = project(system, Y0.copy() + 1e-30, atol_v * 1e-6, 1e-12, Counts(), max_updates=200)
        zstar = truth[0, alg]

        def dist(Y):
            z = Y[0, alg]
            return float(np.max(np.abs(z - zstar) / (atol + rtol * np.abs(zstar))))

        guesses = {"zero": 0.0, "1e-6": 1e-6, "half": 0.5, "double": 2.0, "ten": 10.0}
        for gname, gf in guesses.items():
            Yg = truth.copy()
            Yg[0, alg] = zstar * gf if gname != "zero" else 0.0
            if name == "cold-start":
                Yg[0, alg] = 1e-6 if gname in ("zero", "1e-6") else 0.9 * gf
            for e in range(-14, 15):
                s = 10.0 ** e
                ss = scaled(system, s)
                Yw, stw = project(ss, Yg.copy(), atol_v, rtol, Counts())
                dw = dist(Yw)
                # residual rule: the same damped iteration, stopped on ||G||_inf <= 1e-8
                Y = Yg.copy(); conv_r = False; best = dist(Y)
                for it in range(11):
                    G = ss.rhs(Y); G[:, ~alg] = 0.0
                    if np.all(np.isfinite(G)) and np.abs(G[:, alg]).max() <= 1e-8:
                        conv_r = True
                        break
                    J = ss.jac(Y)
                    LU = ss.factor(1.0, J, algebraic_rows_only=True)
                    d = LU.solve(G)
                    if not np.all(np.isfinite(d)):
                        break
                    Y = Y.copy(); Y[:, alg] += d[:, alg]
                    best = min(best, dist(Y))
                dr = dist(Y)
                rows.append({"problem": name, "guess": gname, "scale": s,
                             "weighted_status": stw, "weighted_dist": dw,
                             "weighted_false_conv": stw == "converged" and dw > 1.0,
                             "residual_converged": conv_r, "residual_dist": dr,
                             "residual_false_conv": conv_r and dr > 1.0,
                             "residual_false_fail": (not conv_r) and best <= 0.1})
    n = len(rows)
    print(f"projection sweep: {n} cases")
    print(f"  weighted rule: false convergences {sum(r['weighted_false_conv'] for r in rows)}, "
          f"failures {sum(r['weighted_status'] != 'converged' for r in rows)}")
    print(f"  residual rule: false convergences {sum(r['residual_false_conv'] for r in rows)}, "
          f"false failures {sum(r['residual_false_fail'] for r in rows)}")
    write("projection", rows)


# ---- (j) Schur fill on TS1 -------------------------------------------------------
FAMILIES = {
    4: ["O1D", "O", "OH", "HO2"],
    9: ["O1D", "O", "OH", "HO2", "CH3O2", "C2H5O2", "CH3CO3", "C3H7O2", "PO2"],
    17: ["O1D", "O", "OH", "HO2", "CH3O2", "C2H5O2", "CH3CO3", "C3H7O2", "PO2",
         "EO2", "ENEO2", "MACRO2", "MCO3", "RO2", "XO2", "ALKO2", "HOCH2OO"],
}


def schur():
    import scipy.sparse as sp
    from fastode.ts1 import TS1
    ts1 = TS1()
    pat = ts1.mech.pattern
    n = ts1.mech.n
    A = sp.csr_matrix((np.ones(pat.nnz), (pat.rows, pat.cols)), shape=(n, n))
    rows = []
    for size, fam in FAMILIES.items():
        z = np.array([ts1.index[s] for s in fam])
        x = np.array([i for i in range(n) if i not in set(z)])
        Axx, Axz, Azx, Azz = A[x][:, x], A[x][:, z], A[z][:, x], A[z][:, z]
        # structural inverse of the algebraic block: dense within its connected components
        ncomp, lab = sp.csgraph.connected_components(Azz, directed=False)
        Zinv = sp.csr_matrix((lab[:, None] == lab[None, :]).astype(float))
        S = (Axx + Axz @ Zinv @ Azx)
        S.data[:] = 1.0
        rows.append({"family": size, "nnz_J": pat.nnz, "nnz_Jxx": Axx.nnz, "nnz_S": S.nnz,
                     "n_x": len(x), "components": ncomp})
        print(f"schur family {size}: nnz(J) {pat.nnz}, nnz(J_xx) {Axx.nnz}, nnz(S) {S.nnz}, "
              f"n_x {len(x)}, algebraic components {ncomp}")
    write("schur", rows)


# ---- (k) Robertson's conservation-form DAE to t = 1e11 ---------------------------------
def conservation():
    rows = []
    times = np.concatenate([[0.0], 10.0 ** np.arange(-3, 12)])
    for form in ("ode", "dae"):
        if form == "dae":
            system, Y0 = P.robertson_conservation()
        else:
            system, Y0 = P.robertson()
        solver = Algorithm(RODAS3)
        mem = solver.new_memory()
        Y = Y0.copy()
        c = Counts()
        drift = 0.0
        for i in range(len(times) - 1):
            Y, st = solver.solve(system, Y, times[i + 1] - times[i], np.full(3, 1e-20), 1e-6, c, mem)
            drift = max(drift, abs(Y[0].sum() - 1.0))
        rows.append({"form": form, "A_1e11": Y[0, 0], "B_1e11": Y[0, 1], "mass_drift": drift, "steps": c.steps})
        print(f"robertson {form} to 1e11: A {Y[0, 0]:.4e} B {Y[0, 1]:.4e} max |A+B+C-1| {drift:.2e} steps {c.steps}")
    print("  test-set values at 1e11: A 2.0833e-08, B 8.3334e-14")
    write("conservation", rows)


ALL = {"orders": orders, "pr": prothero_robinson, "robertson": robertson, "chapman": chapman,
       "tropo": tropospheric, "triples": triples, "vdp": van_der_pol, "akzo": akzo, "projection": projection,
       "schur": schur, "conservation": conservation}



# ---- (l) the restarting loop's two error rows, on Robertson and the tropospheric DAE ----------
def rows_control():
    """The restarting loop with the baseline's step-change rows against the same
    loop with the embedded estimate on the algebraic rows, and the algorithm."""
    rows = []
    cases = [("robertson", lambda q: P.robertson(qssa=q), 1.0e3, 100, 1e-12, None),
             ("tropospheric", lambda q: P.tropospheric(qssa=q), 600.0, 144, 1e-14 * 6.02214076e17, None)]
    for name, build, dt, intervals, a, on in cases:
        full, Y0 = build(False)
        ref = radau_reference(full, Y0, dt, intervals)
        for rtol in (1e-2, 1e-3, 1e-4):
            for label, solver in (("algorithm", Algorithm(RODAS3)),
                                  ("restarting, step-change rows", RestartingLoop(RODAS3)),
                                  ("restarting, embedded rows", RestartingLoop(RODAS3, step_change_rows=False))):
                system, Y0 = build(True)
                if name == "robertson":
                    Y0[0, 1] = P.robertson_manifold(Y0)[0]
                atol = np.full(system.n, rtol * a)
                X, c, bad = run_intervals(system, Y0, solver, dt, intervals, rtol, atol)
                D, _ = digits(X, ref, a)
                rows.append({"problem": name, "rtol": rtol, "solver": label, "digits": round(D, 3),
                             "steps": c.steps, "rejected": c.rejected, "failed_calls": bad})
                print(f"{name} rtol={rtol:g} {label:30s}: {D:5.2f} digits, steps {c.steps}, failed {bad}")
    write("rows_control", rows)


ALL["rows"] = rows_control


# ---- (m) the quadratic systems of Bacsi and Kocsis, fixed steps against a tight reference --------
def qde():
    """QDE-2 and QDE-3, the worked examples of Bacsi and Kocsis (their eqs. 19 and 23),
    on windows free of blow-up, by fixed steps against DOP853 at 1e-14."""
    from scipy.integrate import solve_ivp
    from fastode.system import DenseSystem

    def f2(x):
        x1, x2 = x
        return np.array([-x1 * x1 + 4 * x1 * x2 - x1 + 2 * x2, x1 * x1 + 2 * x2 * x2 + x1])

    def j2(x):
        x1, x2 = x
        return np.array([[-2 * x1 + 4 * x2 - 1, 4 * x1 + 2], [2 * x1 + 1, 4 * x2]])

    def f3(x):
        x1, x2, x3 = x
        return np.array([x1 * x1 + 7 * x2 * x2 - 4 * x1 * x2 + 5 * x1,
                         -2 * x2 * x2 + x1 * x2 + 2 * x1 * x3 - 7 * x2 * x3 + 2 * x2,
                         -x2 * x2 - 7 * x3 * x3 - 4 * x2 * x3 - x3])

    def j3(x):
        x1, x2, x3 = x
        return np.array([[2 * x1 - 4 * x2 + 5, 14 * x2 - 4 * x1, 0.0],
                         [x2 + 2 * x3, -4 * x2 + x1 - 7 * x3 + 2, 2 * x1 - 7 * x2],
                         [0.0, -2 * x2 - 4 * x3, -14 * x3 - 4 * x2 - 1]])

    rows = []
    for name, f, j, x0, T in (("qde2", f2, j2, [0.2, -0.3], 1.0), ("qde3", f3, j3, [-0.2, 0.1, 0.1], 0.5)):
        n = len(x0)
        system = DenseSystem(n, lambda Y, f=f: np.array([f(y) for y in Y]),
                             lambda Y, j=j: np.array([j(y) for y in Y]))
        ref = solve_ivp(lambda t, x: f(x), (0, T), x0, method="DOP853", rtol=1e-14, atol=1e-16).y[:, -1]
        for tab in (ROS3, RODAS3):
            Ns = [8, 16, 32, 64, 128]
            errs = []
            for N in Ns:
                Y = fixed_steps(system, tab, np.array([x0], dtype=float), T, N)[0]
                errs.append(float(np.max(np.abs(Y - ref))))
                rows.append({"problem": name, "tableau": tab.name, "N": N, "error": errs[-1]})
            print(f"{name} {tab.name:6s} slopes {np.round(slopes([T / N for N in Ns], errs), 2)}")
    write("qde", rows)


ALL["qde"] = qde


# ---- (n) a deliberately index-2 formulation: how does the failure show itself? -------------
def index2():
    """t' = 1, y' = z, 0 = y - sin t: the constraint does not involve z, so G_z = 0.
    The projection from an inconsistent start, and both loops from a consistent one."""
    import warnings
    from fastode.system import DenseSystem

    def rhs(Y):
        t, y, z = Y[:, 0], Y[:, 1], Y[:, 2]
        return np.stack([np.ones_like(t), z, y - np.sin(t)], axis=1)

    def jac(Y):
        J = np.zeros((Y.shape[0], 3, 3))
        J[:, 1, 2] = 1.0
        J[:, 2, 0] = -np.cos(Y[:, 0])
        J[:, 2, 1] = 1.0
        return J

    system = DenseSystem(3, rhs, jac, mass=[1.0, 1.0, 0.0], names=["t", "y", "z"])
    atol, rtol, dt = np.full(3, 1e-10), 1e-6, 0.1
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = Counts()
        Y, st = project(system, np.array([[0.0, 0.5, 1.0]]), atol, rtol, c)
        rows.append({"case": "projection, inconsistent start", "status": st, "t_reached": 0.0, "steps": 0,
                     "y_minus_sin_t": float(Y[0, 1] - math.sin(Y[0, 0]))})
        for label, solver in (("algorithm", Algorithm(RODAS3)), ("restarting", RestartingLoop(RODAS3))):
            c = Counts()
            Y, st = solver.solve(system, np.array([[0.0, 0.0, 1.0]]), dt, atol, rtol, c, solver.new_memory())
            rows.append({"case": f"{label}, consistent start", "status": st, "t_reached": float(Y[0, 0]),
                         "steps": c.steps, "y_minus_sin_t": float(Y[0, 1] - math.sin(Y[0, 0]))})
    for r in rows:
        print(f"index 2, {r['case']}: {r['status']}, t {r['t_reached']:.4f} of {dt}, steps {r['steps']}")
    write("index2", rows)


ALL["index2"] = index2


if __name__ == "__main__":
    which = sys.argv[1:] or ["all"]
    for name in (list(ALL) if which == ["all"] else which):
        print(f"==== {name}", flush=True)
        ALL[name]()
