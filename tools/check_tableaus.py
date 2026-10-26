"""Verify the Rosenbrock tableaus of Appendix C against the order conditions
of Appendix A.

Each tableau is given in the transformed (KPP) form of eq. (ros-implform):
gamma (the shared diagonal), a_ij and c_ij for j < i (flat KPP index
(i-1)(i-2)/2 + j - 1 with 1-based i, j), the update weights m_i and the
error weights e_i. The script converts to the standard form
(Gamma = ((1/gamma) I - C)^{-1}, A = Atilde Gamma, b = Gamma^T m,
b - bhat = Gamma^T e), evaluates the order conditions through order four
with the diagonal-inclusive sums of Appendix A (beta = A + Gamma,
a'_r = sum_s a_rs, beta'_r = sum_s beta_rs), the stiff-accuracy condition
b_r = beta_Sr, R(inf) = 1 - b^T B^{-1} 1 for the method and its companion,
and the index-1 conditions b^T B^{-1} a'^q = 1 for q = 2, 3 (Theorem
thm:roche), which the local error of the algebraic variables needs at
orders three and four.

Pure Python plus numpy; run: python3 tools/check_tableaus.py
"""
import numpy as np

def unpack(S, flat):
    """Lower-triangular S x S matrix from the KPP flat layout."""
    M = np.zeros((S, S))
    k = 0
    for i in range(1, S):
        for j in range(i):
            M[i, j] = flat[k]
            k += 1
    return M

TABLEAUS = {
    "ros3": dict(
        S=3, gamma=0.43586652150845899941601945119356,
        a=[1.0, 1.0, 0.0],
        c=[-0.10156171083877702091975600115545e+01,
           0.40759956452537699824805835358067e+01,
           0.92076794298330791242156818474003e+01],
        m=[1.0, 0.61697947043828245592553615689730e+01, -0.42772256543218573326238373806514],
        e=[0.5, -0.29079558716805469821718236208017e+01, 0.22354069897811569627360909276199],
        order=3, embedded=2, source="Sandu et al. 1997"),
    "rodas3": dict(
        S=4, gamma=0.5,
        a=[0.0, 2.0, 0.0, 2.0, 0.0, 1.0],
        c=[4.0, 1.0, -1.0, 1.0, -1.0, -8.0 / 3.0],
        m=[2.0, 0.0, 1.0, 1.0],
        e=[0.0, 0.0, 0.0, 1.0],
        order=3, embedded=2, source="Sandu et al. 1997"),
}

def standard_form(t):
    S, g = t["S"], t["gamma"]
    At, C = unpack(S, t["a"]), unpack(S, t["c"])
    Ginv = np.eye(S) / g - C                 # Gamma^{-1} = (1/gamma) I - C
    Gamma = np.linalg.inv(Ginv)              # lower triangular, diagonal gamma
    A = At @ Gamma                           # strictly lower
    b = Gamma.T @ np.array(t["m"])
    bhat = b - Gamma.T @ np.array(t["e"])
    return A, Gamma, b, bhat

def conditions(A, Gamma, b):
    """Residuals of the order conditions, diagonal-inclusive convention."""
    beta = A + Gamma
    ap = A.sum(axis=1)                       # a'_r
    bp = beta.sum(axis=1)                    # beta'_r (includes gamma)
    res = {
        1: [b.sum() - 1.0],
        2: [b @ bp - 0.5],
        3: [b @ ap**2 - 1.0/3.0, b @ (beta @ bp) - 1.0/6.0],
        4: [b @ ap**3 - 0.25,
            b @ (ap * (A @ bp)) - 1.0/8.0,
            b @ (beta @ ap**2) - 1.0/12.0,
            b @ (beta @ (beta @ bp)) - 1.0/24.0],
    }
    return res

def observed_order(res, tol=1e-10):
    p = 0
    for k in (1, 2, 3, 4):
        if all(abs(r) < tol for r in res[k]):
            p = k
        else:
            break
    return p

def r_infinity(B, b):
    one = np.ones(len(b))
    return 1.0 - b @ np.linalg.solve(B, one)

def index1(B, A, b, q):
    """b^T B^{-1} a'^q, which the index-1 condition of that order asks to be 1."""
    return b @ np.linalg.solve(B, A.sum(axis=1)**q)

def main():
    print(f"{'tableau':8s} {'S':>2s} {'order':>5s} {'emb.':>4s} {'stiff.acc.':>10s} {'R(inf)':>9s} {'Rhat(inf)':>10s} {'idx1 q=2':>9s} {'idx1 q=3':>9s}  source")
    ok = True
    for name, t in TABLEAUS.items():
        A, Gamma, b, bhat = standard_form(t)
        B = A + Gamma
        p = observed_order(conditions(A, Gamma, b))
        phat = observed_order(conditions(A, Gamma, bhat))
        stiff = np.allclose(b, B[-1, :], atol=1e-12)
        rinf, rhat = r_infinity(B, b), r_infinity(B, bhat)
        flag = (p == t["order"] and phat == t["embedded"])
        ok &= flag
        print(f"{name:8s} {t['S']:2d} {p:5d} {phat:4d} {str(stiff):>10s} {rinf:9.2e} {rhat:10.3f} {index1(B, A, b, 2):9.6f} {index1(B, A, b, 3):9.6f}  {t['source']}{'' if flag else '   <-- MISMATCH'}")
    print("all orders as declared:", ok)

if __name__ == "__main__":
    main()
