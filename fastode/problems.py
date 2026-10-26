"""The small test problems of Appendix D, as systems of one cell.

Chemistry problems are Mechanisms, so that holding a species in quasi-steady
state is a matter of its mass entry: the constraint is its own rate equation.
The others are given by their right-hand sides and Jacobians.
"""
from __future__ import annotations

import numpy as np

from .chem import ChemSystem, Mechanism
from .system import DenseSystem


# ---- Robertson -------------------------------------------------------------
def robertson(k2=3.0e7, qssa=False):
    mech = Mechanism(["A", "B", "C"], [
        (["A"], [("B", 1.0)]),
        (["B", "B"], [("B", 1.0), ("C", 1.0)]),
        (["B", "C"], [("A", 1.0), ("C", 1.0)]),
    ])
    k = np.array([[0.04, k2, 1.0e4]])
    mass = np.array([1.0, 0.0, 1.0]) if qssa else None
    return ChemSystem(mech, k, mass=mass), np.array([[1.0, 0.0, 0.0]])


def robertson_manifold(Y, k2=3.0e7):
    A, C = Y[:, 0], Y[:, 2]
    k1, k3 = 0.04, 1.0e4
    return (-k3 * C + np.sqrt((k3 * C) ** 2 + 4.0 * k1 * k2 * A)) / (2.0 * k2)


def robertson_conservation():
    """The test set's conservation form: the C row replaced by A + B + C = 1."""
    k1, k2, k3 = 0.04, 3.0e7, 1.0e4

    def rhs(Y):
        A, B, C = Y[:, 0], Y[:, 1], Y[:, 2]
        return np.stack([-k1 * A + k3 * B * C, k1 * A - k2 * B * B - k3 * B * C, A + B + C - 1.0], axis=1)

    def jac(Y):
        A, B, C = Y[:, 0], Y[:, 1], Y[:, 2]
        J = np.zeros((Y.shape[0], 3, 3))
        J[:, 0, 0] = -k1; J[:, 0, 1] = k3 * C; J[:, 0, 2] = k3 * B
        J[:, 1, 0] = k1; J[:, 1, 1] = -2 * k2 * B - k3 * C; J[:, 1, 2] = -k3 * B
        J[:, 2, :] = 1.0
        return J

    return DenseSystem(3, rhs, jac, mass=[1.0, 1.0, 0.0], names=["A", "B", "C"]), np.array([[1.0, 0.0, 0.0]])


# ---- Chapman, dimensionless --------------------------------------------------
def chapman(scale=1.0, qssa=False):
    """x1 = O, x2 = O2, x3 = O3; 'scale' multiplies k2 and k3 (the stiffness knob)."""
    mech = Mechanism(["O", "O2", "O3"], [
        (["O2"], [("O", 2.0)]),
        (["O", "O2"], [("O3", 1.0)]),
        (["O3"], [("O", 1.0), ("O2", 1.0)]),
        (["O", "O3"], [("O2", 2.0)]),
    ])
    k = np.array([[1.0e-10, 1.0e-4 * scale, 1.0e-3, 1.0e-4 * scale]])
    mass = np.array([0.0, 1.0, 1.0]) if qssa else None
    Y0 = np.array([[0.0, 2.0e4, 1.0]])
    if qssa:
        x2, x3 = Y0[0, 1], Y0[0, 2]
        Y0[0, 0] = (2e-10 * x2 + 1e-3 * x3) / (1e-4 * scale * (x2 + x3))
    return ChemSystem(mech, k, mass=mass), Y0


# ---- the tropospheric O3-NOx-HOx mechanism (Appendix D, Table) -----------------
TROPO_SPECIES = ["O1D", "O", "OH", "HO2", "O3", "NO", "NO2", "HNO3", "CO", "CO2", "CH4", "H2O",
                 "H2O2", "O2", "N2"]
TROPO_RADICALS = ["O1D", "O", "OH", "HO2"]
TROPO_J = {"p1": 8.0e-3, "p2": 3.0e-5, "p3": 4.0e-4}


def tropospheric(j=1.0, qssa=False, T=298.15, M=2.5e19):
    """16 species counting the third body M (folded into the rate constants),
    14 reactions; molecules cm-3 and s."""
    rx = [
        (["NO2"], [("NO", 1.0), ("O", 1.0)]),                       # p1
        (["O3"], [("O1D", 1.0), ("O2", 1.0)]),                      # p2
        (["O3"], [("O", 1.0), ("O2", 1.0)]),                        # p3
        (["O1D", "N2"], [("O", 1.0), ("N2", 1.0)]),                 # r1
        (["O1D", "O2"], [("O", 1.0), ("O2", 1.0)]),                 # r2
        (["O1D", "H2O"], [("OH", 2.0)]),                            # r3
        (["O", "O2", "M"], [("O3", 1.0)]),                          # r4
        (["O", "O3"], [("O2", 2.0)]),                               # r5
        (["O3", "NO"], [("NO2", 1.0), ("O2", 1.0)]),                # r6
        (["OH", "CO"], [("HO2", 1.0), ("CO2", 1.0)]),               # r7
        (["OH", "CH4"], [("HO2", 1.0)]),                            # r8
        (["HO2", "NO"], [("OH", 1.0), ("NO2", 1.0)]),               # r9
        (["OH", "NO2"], [("HNO3", 1.0)]),                           # r10
        (["HO2", "HO2"], [("H2O2", 1.0), ("O2", 1.0)]),             # r11
    ]
    mech = Mechanism(TROPO_SPECIES, rx)

    def arr(A, C=0.0, B=0.0):
        return A * np.exp(C / T) * (T / 300.0) ** B

    k = np.array([[TROPO_J["p1"] * j, TROPO_J["p2"] * j, TROPO_J["p3"] * j,
                   arr(2.15e-11, 110), arr(3.30e-11, 55), arr(1.63e-10, 60), arr(6.0e-34, 0, -2.4) * M,
                   arr(8.0e-12, -2060), arr(3.0e-12, -1500), arr(1.5e-13), arr(2.45e-12, -1775),
                   arr(3.3e-12, 270), arr(1.0e-11), arr(3.0e-13, 460)]])
    ppb = 1e-9 * M
    Y0 = np.zeros((1, len(TROPO_SPECIES)))
    init = {"O3": 30 * ppb, "NO": 1 * ppb, "NO2": 1 * ppb, "CO": 100 * ppb, "CH4": 1800 * ppb,
            "H2O": 0.01 * M, "O2": 0.21 * M, "N2": 0.78 * M}
    for s, v in init.items():
        Y0[0, TROPO_SPECIES.index(s)] = v
    mass = None
    if qssa:
        mass = np.ones(len(TROPO_SPECIES))
        for s in TROPO_RADICALS:
            mass[TROPO_SPECIES.index(s)] = 0.0
    return ChemSystem(mech, k, mass=mass), Y0


def tropospheric_photolysis(system, j):
    system.k[0, :3] = [TROPO_J["p1"] * j, TROPO_J["p2"] * j, TROPO_J["p3"] * j]


def daylight(t):
    """The diurnal factor: a half-sine between 06:00 and 18:00, t in seconds from 00:00."""
    hour = (t / 3600.0) % 24.0
    return np.sin(np.pi * (hour - 6.0) / 12.0) if 6.0 < hour < 18.0 else 0.0


# ---- fast-equilibrium triples ------------------------------------------------------
def triples(N=10, kappa=1.0e6, ks=1.0, dae=False):
    """N independent triples A <=> B (kappa, K_eq = 1), B -> C (ks). The DAE holds
    A and B on A - B = 0 and A + B + C = 1, with C differential."""
    if not dae:
        names, rx = [], []
        for i in range(N):
            names += [f"A{i}", f"B{i}", f"C{i}"]
            rx += [([f"A{i}"], [(f"B{i}", 1.0)]), ([f"B{i}"], [(f"A{i}", 1.0)]), ([f"B{i}"], [(f"C{i}", 1.0)])]
        k = np.tile([kappa, kappa, ks], (1, N))
        Y0 = np.tile([1.0, 0.0, 0.0], (1, N))
        return ChemSystem(Mechanism(names, rx), k), Y0

    def rhs(Y):
        F = np.empty_like(Y)
        A, B, C = Y[:, 0::3], Y[:, 1::3], Y[:, 2::3]
        F[:, 0::3] = A - B
        F[:, 1::3] = A + B + C - 1.0
        F[:, 2::3] = ks * B
        return F

    def jac(Y):
        n = 3 * N
        J = np.zeros((Y.shape[0], n, n))
        for i in range(N):
            a, b, c = 3 * i, 3 * i + 1, 3 * i + 2
            J[:, a, a] = 1.0; J[:, a, b] = -1.0
            J[:, b, a] = 1.0; J[:, b, b] = 1.0; J[:, b, c] = 1.0
            J[:, c, b] = ks
        return J

    mass = np.tile([0.0, 0.0, 1.0], N)
    Y0 = np.tile([0.5, 0.5, 0.0], (1, N))
    return DenseSystem(3 * N, rhs, jac, mass=mass), Y0


# ---- the rigor suite ----------------------------------------------------------------
def hires():
    def rhs(Y):
        y = Y.T
        F = np.stack([
            -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007,
            1.71 * y[0] - 8.75 * y[1],
            -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4],
            8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3],
            -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6],
            -280.0 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6],
            280.0 * y[5] * y[7] - 1.81 * y[6],
            -280.0 * y[5] * y[7] + 1.81 * y[6],
        ])
        return F.T

    def jac(Y):
        J = np.zeros((Y.shape[0], 8, 8))
        y6, y8 = Y[:, 5], Y[:, 7]
        J[:, 0, 0] = -1.71; J[:, 0, 1] = 0.43; J[:, 0, 2] = 8.32
        J[:, 1, 0] = 1.71; J[:, 1, 1] = -8.75
        J[:, 2, 2] = -10.03; J[:, 2, 3] = 0.43; J[:, 2, 4] = 0.035
        J[:, 3, 1] = 8.32; J[:, 3, 2] = 1.71; J[:, 3, 3] = -1.12
        J[:, 4, 4] = -1.745; J[:, 4, 5] = 0.43; J[:, 4, 6] = 0.43
        J[:, 5, 3] = 0.69; J[:, 5, 4] = 1.71; J[:, 5, 5] = -0.43 - 280.0 * y8; J[:, 5, 6] = 0.69; J[:, 5, 7] = -280.0 * y6
        J[:, 6, 5] = 280.0 * y8; J[:, 6, 6] = -1.81; J[:, 6, 7] = 280.0 * y6
        J[:, 7, 5] = -280.0 * y8; J[:, 7, 6] = 1.81; J[:, 7, 7] = -280.0 * y6
        return J

    Y0 = np.array([[1.0, 0, 0, 0, 0, 0, 0, 0.0057]])
    return DenseSystem(8, rhs, jac), Y0


AKZO = dict(k1=18.7, k2=0.58, k3=0.09, k4=0.42, K=34.4, klA=3.3, Ks=115.83, pCO2=0.9, H=737.0)


def akzo(dae=True):
    """Chemical Akzo Nobel (test set): index 1 with y6 algebraic, or its exact
    five-species reduction (y6 = Ks y1 y4 substituted) when dae is False."""
    p = AKZO

    def full(Y):
        y1, y2, y3, y4, y5 = (Y[:, i] for i in range(5))
        y6 = Y[:, 5] if dae else p["Ks"] * y1 * y4
        sq = np.sqrt(np.maximum(y2, 0.0))
        r1 = p["k1"] * y1 ** 4 * sq
        r2 = p["k2"] * y3 * y4
        r3 = p["k2"] / p["K"] * y1 * y5
        r4 = p["k3"] * y1 * y4 ** 2
        r5 = p["k4"] * y6 ** 2 * sq
        fin = p["klA"] * (p["pCO2"] / p["H"] - y2)
        cols = [-2 * r1 + r2 - r3 - r4, -0.5 * r1 - r4 - 0.5 * r5 + fin, r1 - r2 + r3, -r2 + r3 - 2 * r4,
                r2 - r3 + r5]
        if dae:
            cols.append(p["Ks"] * y1 * y4 - y6)
        return np.stack(cols, axis=1)

    n = 6 if dae else 5

    def jac(Y, eps=1e-8):
        # numerical Jacobian by central differences, adequate for a reference-checked test
        C = Y.shape[0]
        J = np.zeros((C, n, n))
        for i in range(n):
            d = np.zeros_like(Y)
            h = eps * np.maximum(1.0, np.abs(Y[:, i]))
            d[:, i] = h
            J[:, :, i] = (full(Y + d) - full(Y - d)) / (2 * h[:, None])
        return J

    y0 = [0.444, 0.00123, 0.0, 0.007, 0.0]
    if dae:
        y0.append(p["Ks"] * 0.444 * 0.007)
    mass = [1, 1, 1, 1, 1, 0] if dae else None
    return DenseSystem(n, full, jac, mass=mass), np.array([y0])


def van_der_pol(eps, dae=False):
    """X' = -Z, eps Z' = (1 - X^2) Z + X (Lienard form, fast variable negated);
    the DAE is the eps = 0 limit."""
    def rhs(Y):
        X, Z = Y[:, 0], Y[:, 1]
        G = (1.0 - X * X) * Z + X
        return np.stack([-Z, G if dae else G / eps], axis=1)

    def jac(Y):
        X, Z = Y[:, 0], Y[:, 1]
        J = np.zeros((Y.shape[0], 2, 2))
        J[:, 0, 1] = -1.0
        s = 1.0 if dae else 1.0 / eps
        J[:, 1, 0] = s * (-2.0 * X * Z + 1.0)
        J[:, 1, 1] = s * (1.0 - X * X)
        return J

    Y0 = np.array([[2.0, 2.0 / 3.0]])
    return DenseSystem(2, rhs, jac, mass=[1.0, 0.0] if dae else None, names=["X", "Z"]), Y0


def prothero_robinson(lam=-1.0e6):
    """Y' = lam (Y - sin t) + cos t, autonomized with t as the first variable."""
    def rhs(Y):
        t, y = Y[:, 0], Y[:, 1]
        return np.stack([np.ones_like(t), lam * (y - np.sin(t)) + np.cos(t)], axis=1)

    def jac(Y):
        t = Y[:, 0]
        J = np.zeros((Y.shape[0], 2, 2))
        J[:, 1, 0] = -lam * np.cos(t) - np.sin(t)
        J[:, 1, 1] = lam
        return J

    return DenseSystem(2, rhs, jac, names=["t", "Y"]), np.array([[0.0, 0.0]])


def manufactured_circle():
    """y' = z, 0 = y^2 + z^2 - 1 from (0, 1): y = sin t, z = cos t."""
    def rhs(Y):
        y, z = Y[:, 0], Y[:, 1]
        return np.stack([z, y * y + z * z - 1.0], axis=1)

    def jac(Y):
        J = np.zeros((Y.shape[0], 2, 2))
        J[:, 0, 1] = 1.0
        J[:, 1, 0] = 2 * Y[:, 0]
        J[:, 1, 1] = 2 * Y[:, 1]
        return J

    return DenseSystem(2, rhs, jac, mass=[1.0, 0.0], names=["y", "z"]), np.array([[0.0, 1.0]])


def manufactured_coupled():
    """y1' = -y1 + z, y2' = y1 z - y2, 0 = y1^2 + y2 - z from (0.5, 0.5, 0.75)."""
    def rhs(Y):
        y1, y2, z = Y[:, 0], Y[:, 1], Y[:, 2]
        return np.stack([-y1 + z, y1 * z - y2, y1 * y1 + y2 - z], axis=1)

    def jac(Y):
        y1, z = Y[:, 0], Y[:, 2]
        J = np.zeros((Y.shape[0], 3, 3))
        J[:, 0, 0] = -1.0; J[:, 0, 2] = 1.0
        J[:, 1, 0] = z; J[:, 1, 1] = -1.0; J[:, 1, 2] = y1
        J[:, 2, 0] = 2 * y1; J[:, 2, 1] = 1.0; J[:, 2, 2] = -1.0
        return J

    return DenseSystem(3, rhs, jac, mass=[1.0, 1.0, 0.0], names=["y1", "y2", "z"]), np.array([[0.5, 0.5, 0.75]])
