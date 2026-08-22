"""LIU'S HYPOTHESIS 1, REPRODUCED AND STRUCTURALLY ANALYSED.

Liu (arXiv:2306.08824v1, CISS 2024) proves a union-closed constant
    c' = 0.382709087918741 > c* = 0.3823455333667027
i.e. ABOVE the two-strategy ceiling this project's certificates live under --
but conditionally on two hypotheses he verifies only numerically.  The first
one (his Section V-A, "Positive-Semidefiniteness") is what this file attacks.

WHAT HYPOTHESIS 1 SAYS.  Let h be binary entropy (base 2),
    a(s) := s(1-s),      z(s,t) := (1-s)(1-t) + a(s)a(t),
    K(s,t) := h(z(s,t)).
Claim: the quadratic form mu |-> int int K(s,t) mu(ds) mu(dt) is <= 0 for
every signed measure mu annihilating the three functions
    1,   s,   a(s) = s(1-s),
i.e. -K is positive semidefinite on that codimension-THREE subspace.  (The
project ledger earlier recorded "codimension 2" second-hand; Liu's own
verification code frankl3.m builds A = [1, i/n, a(i)] and projects with
P = I - A(A'A)^{-1}A', so it is three constraints.  Recorded here from the
primary source.)

LIU'S OWN CHECK (frankl3.m, fetched 2026-08-22): uniform grid
s_i = i*dt, dt = 4e-4, i = 1..2499; H_ij = h(z(s_i,s_j)); P projects out
[1, i/n, a(i)]; H1 = P(4:n,:) H P(:,4:n), symmetrised; reports
max eig(H1) = 1.6311e-14 -- zero to floating-point noise, since eps = 2.2e-16
and a 2500x2500 matrix of that scale carries roundoff ~1e-13.

WHY THIS PROJECT CARES.  Our Theorem A (uc/decomposition.py) is exactly this
kind of statement, proved exactly: h(xy) = g@g - a@a - T with T PSD, a sympy
tautology.  Liu's kernel is h(xy + a(s)a(t)) with x = 1-s, y = 1-t: the same
OR-entropy kernel with a rank-one perturbation inside the entropy argument.
If the perturbed kernel admits an analogous exact decomposition, Hypothesis 1
becomes a theorem and Liu's constant stops being conditional -- a constant
strictly above the ceiling of our own certified route.

WHAT THIS FILE ESTABLISHES (all NUMERICAL unless marked otherwise):
 1. An independent reproduction of Liu's check, with our own grid handling and
    an eigenvalue solver on the symmetrised projected matrix.
 2. Grid-refinement behaviour of the top eigenvalue (does it stay at roundoff
    as dt shrinks, or is it a positive number the coarse grid hides?).
 3. The structural decomposition: expand h(z) in the series that makes the
    unperturbed kernel's sign structure exact, and report which pieces are
    NSD, which are PSD, and where the codimension-3 subspace is exactly what
    kills the positive part.  This is the map an exact proof must follow.

Run: ./.venv/bin/python uc/liu_kernel.py            (default dt = 2e-3)
     LIU_DT=4e-4 ./.venv/bin/python uc/liu_kernel.py   (Liu's own grid)
"""

import os
from math import comb

import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LOG2 = np.log(2.0)
C_STAR = 0.3823455333667027
C_LIU = 0.382709087918741


def h(x):
    """Binary entropy, base 2, with h(0) = h(1) = 0; vectorised."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    m = (x > 0.0) & (x < 1.0)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1.0 - xm) * np.log1p(-xm)) / LOG2
    return out


def liu_matrix(dt):
    """H_ij = h(z(s_i,s_j)) on Liu's grid, plus the grid itself."""
    s = np.arange(1, int(round(1.0 / dt)))* dt
    a = s * (1.0 - s)
    z = np.outer(1.0 - s, 1.0 - s) + np.outer(a, a)
    return s, a, h(z)


def projected_top_eig(dt, drop=3):
    """Top eigenvalue of the kernel restricted to {1, s, a}^perp.

    Liu projects with P = I - A(A'A)^{-1}A' for A = [1, s, a] and then keeps
    rows/columns drop+1..n.  Dropping rows is only a way to discard the exact
    zero modes P introduces; we instead build an orthonormal basis of the
    orthogonal complement and form the (n-3)x(n-3) compression, which is the
    same quadratic form without the rank-deficiency bookkeeping.
    """
    s, a, H = liu_matrix(dt)
    n = len(s)
    A = np.column_stack([np.ones(n), s, a])
    # Orthonormal basis of A^perp via QR on the full space.
    Q, _ = np.linalg.qr(A)                      # n x 3, orthonormal columns
    # Householder-style complement: take the orthogonal complement basis from
    # a full QR of A padded to n columns.
    Qf, _ = np.linalg.qr(np.column_stack([A, np.eye(n)[:, : n - 3]]))
    B = Qf[:, 3:]                               # n x (n-3), orthonormal, A'B = 0
    M = B.T @ H @ B
    M = 0.5 * (M + M.T)
    ev = np.linalg.eigvalsh(M)
    return n, float(ev[-1]), float(ev[0]), float(np.abs(A.T @ B).max())


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    print("c*      = %.16f   (our route's ceiling, hard for every alpha)" % C_STAR)
    print("c_Liu   = %.16f   (conditional on Hypotheses 1 and 2)" % C_LIU)
    print("gain    = %+.3e  -- why Hypothesis 1 is worth proving" % (C_LIU - C_STAR))
    print()

    print("1. Independent reproduction of Liu's Hypothesis-1 check [NUMERICAL]")
    print("   %-10s %-7s %-16s %-16s %s" % ("dt", "n", "max eig", "min eig",
                                            "orthogonality"))
    default_dt = float(os.environ.get("LIU_DT", "2e-3"))
    grids = sorted({8e-3, 4e-3, default_dt, 2e-3}, reverse=True)
    for dt in grids:
        n, hi, lo, orth = projected_top_eig(dt)
        verdict = "<= 0 to roundoff" if hi < 1e-9 else "POSITIVE"
        print("   %-10.1e %-7d %+16.6e %+16.6e %.1e   %s"
              % (dt, n, hi, lo, orth, verdict))
    print("2. THE REDUCTION: Hypothesis 1 as an explicit moment inequality")
    print("   Writing u = 1-s, v = 1-t, Liu's kernel argument factors:")
    print("       z(s,t) = (1-s)(1-t) + s(1-s)t(1-t) = u v (1 + s t),")
    print("   and with -(1-z)ln(1-z) = z - sum_{k>=2} z^k/(k(k-1)),")
    print("       ln2 * h(z) = -z ln u - z ln v - z ln(1+st) + z")
    print("                    - sum_{k>=2} z^k/(k(k-1))          [exact]")
    print("   Every term of the first line carries a factor that is a")
    print("   polynomial of degree <= 2 in one variable:")
    print("     -z ln u = -(u ln u)(v) - (u s ln u)(v t),   v and vt=a(t) are such,")
    print("     +z      = (u)(v) + a(s)a(t),")
    print("   so ALL of them are annihilated by Liu's codimension-3 projection")
    print("   {1, s, s(1-s)}^perp -- which is exactly why the projection has")
    print("   three constraints and not two.  What survives is diagonal in the")
    print("   moments M_{k,i} = <mu, (1-s)^k s^i>:")
    print()
    print("     ln2 * int int h(z) dmu dmu")
    print("        =  sum_{j>=2} ((-1)^(j+1)/(j(j-1))) M_{1,j}^2")
    print("          - sum_{k>=2} sum_{i=0..k} (C(k,i)/(k(k-1))) M_{k,i}^2")
    print()
    print("   The ONLY positive weights are the odd j >= 3 terms, weight")
    print("   1/(j(j-1)); j = 1 is killed by the projection (f_1 = s(1-s)).")
    print("   Since M_{k,i} are linear in the plain moments p_r = <mu, s^r>")
    print("   with p_0 = p_1 = p_2 = 0, and finite atomic signed measures")
    print("   realise arbitrary finite moment vectors (Vandermonde), the")
    print("   hypothesis is EQUIVALENT to an explicit infinite real symmetric")
    print("   matrix being negative semidefinite -- no grid, no quadrature.")
    print()
    print("   Numerical check of the reduction (random projected measures):")
    print("   %-6s %-20s %-20s %s" % ("K", "direct mu'H mu", "moment series", "|diff|"))
    import numpy.random as _r
    rng = _r.default_rng(7)
    nn = 400
    ss = (np.arange(1, nn + 1)) / (nn + 1.0)
    uu = 1.0 - ss
    zz = np.outer(uu, uu) * (1.0 + np.outer(ss, ss))
    HH = -(zz * np.log(zz) + (1.0 - zz) * np.log1p(-zz))     # ln2 * h(z)
    AA = np.column_stack([np.ones(nn), ss, ss * ss])
    QQ, _ = np.linalg.qr(AA)
    mm = rng.normal(size=nn)
    mm = (mm - QQ @ (QQ.T @ mm)) / nn
    direct = float(mm @ HH @ mm)
    for K in (20, 60, 150, 300):
        tot = 0.0
        for j in range(2, K + 2):
            Mj = float(mm @ (uu * ss ** j))
            tot += ((-1) ** (j + 1) / (j * (j - 1))) * Mj * Mj
        for k in range(2, K + 1):
            uk = uu ** k
            for i in range(0, k + 1):
                M = float(mm @ (uk * ss ** i))
                tot -= (comb(k, i) / (k * (k - 1))) * M * M
        print("   %-6d %+20.10e %+20.10e %.2e" % (K, direct, tot, abs(direct - tot)))
    print()
    print("3. What a rigorous proof of Hypothesis 1 now needs")
    print("   (a) assemble the moment matrix G in the plain-moment basis with")
    print("       exact rational coefficients (the weights are rationals;")
    print("       the basis change (1-s)^k s^i -> s^r is binomial);")
    print("   (b) certify -G_R negative semidefinite for a truncation degree R")
    print("       by an interval Cholesky (Arb), exactly as this project")
    print("       certifies Phi >= 0 box by box;")
    print("   (c) bound the discarded tail: the k-sum tail is controlled by")
    print("       sup z = 1 at the corner s,t -> 0, so the tail bound must be")
    print("       taken in a weighted norm -- this is the one genuinely open")
    print("       analytic step, and it is a one-dimensional estimate, not a")
    print("       9-parameter search.")
    print("   Contrast with Liu's own verification: a float64 eigenvalue solve")
    print("   on a 2499-point grid whose top eigenvalue (1.6e-14) sits below")
    print("   that computation's own roundoff floor.")
    print()
    print("NUMERICAL ONLY.  Nothing here is machine-checked in the certified")
    print("sense; this file reproduces and maps Liu's hypothesis, it does not")
    print("prove it.")
