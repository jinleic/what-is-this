"""norm_unit_test — Main's ordered test: unit-test the L^2(mu) norm layer against exact
Gaussian moments, no paper reading, no conventions. Feeds monomials x^i y^j into the ACTUAL
norm machinery the S-sweep uses (the Q-family Gram / F-matrix on the GL grid, and the
closed-form K^2-style 2-D route for comparison) and compares against (2i-1)!!(2j-1)!!.

If the norm layer carries a monomial/Hermite bookkeeping defect (ratios 1, 1.5, 2.5 at
degrees 1,2,3 cited by Main as the manufacturing route), it shows here.
REPORT EVERY VALUE INCLUDING MATCHES.
"""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
import d3h_sweep49_asm as A   # provides XS, WS, UX, WIJ, G_pa_terms, mult_terms, wvec, rho_c, ctab

def norm_of_matrix(F):
    return float(np.real(np.sum(A.WIJ * np.abs(F) ** 2)))

def mono_matrix(i, j):
    """F(x,y) = x^i y^j on the grid (the L2(mu) integrand object the norms square)."""
    fx = A.XS ** i
    fy = A.XS ** j
    return np.outer(fx, fy)

# EXPECTED: ||x^i y^j||^2_{L2(mu)} = E[x^{2i}] E[y^{2j}] = (2i-1)!!(2j-1)!!
def dfact(n):
    if n <= 0:
        return 1
    return n * dfact(n - 2)

print("== L2(mu) norm unit test on the sweep's GL grid (XW=6, P=18, GM=24, N/axis=432) ==")
print("  (the weight is dmu = e^{-x^2/2}/sqrt(2pi) folded into WS; the quadrature domain is [-6,6])")
print()
mono_cases = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 2), (3, 0), (0, 3),
              (3, 3), (2, 1), (1, 2), (4, 0), (0, 4)]
for (i, j) in mono_cases:
    F = mono_matrix(i, j)
    got = norm_of_matrix(F)
    want = dfact(2 * i - 1) * dfact(2 * j - 1)
    rel = abs(got - want) / want if want else abs(got)
    print(f"  ||x^{i} y^{j}||^2: machine = {got:.10f}   exact = {want:.10f}   rel = {rel:.2e}")

print()
print("== same monomials through the W-WEIGHTED kernel-style assembly path (K_r with r=0) ==")
# The chain's norm route squares functions built as (weight * mult * K). At r = 0 the kernel is
# K_0(u,v) = phi(u) phi(v) — a product kernel; Phi-like objects built from MONOMIAL weights
# should reduce to Hermite-moment algebra. Test: feed F = psi3'(x) psi3'(y) x^i y^j * 2pi vt^2 K_0:
import d3h_sweep49_asm as AS
r0 = 0.0
s0 = 1.0
Km0 = np.exp(-(A.UX[:, None] ** 2 + A.UX[None, :] ** 2) / 2) / (2 * math.pi)
import d3h_leg2_indep as M
for (i, j) in [(1, 1), (2, 2), (3, 3)]:
    F = np.outer(A.XS ** i, A.XS ** j) * Km0
    got = norm_of_matrix(F)
    # Ground truth via 1-D moments: int x^{2i} exp(-u^2) dmu with u = vt psi3(x):
    # compute exact via mpmath quad for independence:
    import mpmath
    from mpmath import mp, mpf
    mp.dps = 30
    vtf = mpf('0.128957369921412')
    def mom(ii):
        f = lambda x: x ** (2 * ii) * mpmath.e ** (-(vtf * (x ** 3 - 3 * x) / mpmath.sqrt(6)) ** 2) * mpmath.e ** (-(x * x) / 2) / mpmath.sqrt(2 * mpmath.pi)
        return mpmath.quad(f, [-mpmath.inf, mpmath.inf])
    want = float(mom(i) * mom(j))
    rel = abs(got - want) / want
    print(f"  ||x^{i} y^{j} K_0||^2 (r=0 kernel path): machine = {got:.8e}   exact = {want:.8e}   rel = {rel:.2e}")

print()
print("== CROSS-CHECK the Gram/G-machinery itself on pure Gauss-Hermite ground ==")
# The route internally computes <Q Q> as Q^T W Q. Ground truth for Q rows: Q^0_b = q_b(u(x)):
# <q_b(u) q_{b'}(u)>_{mu} with u = vt psi3(x): again not closed form; instead verify the
# QUADRATURE layer: int x^k dmu vs exact moments on our grid (this is what the norm reduces to):
for k in range(0, 13, 2):
    got = float(np.sum(A.WS * A.XS ** k))
    want = dfact(k - 1)
    print(f"  int x^{k} dmu: machine = {got:.12f}   exact = {want:.12f}   rel = {abs(got-want)/want:.2e}")

print()
print("== K^2 route: ||2 pi K||^2 at r=0 (pristine product kernel; closed form E[K^2]) ==")
# ||2piK_0||^2 = 4 pi^2 (int K_0^2 dmu dmu): K_0(u,v) = phi(u) phi(v), u = vt psi3(x):
# E[phi(u)^2] has no closed form in x — BUT at the OPERATOR level with r=0, C_r = (pi/2) q_0 q_0:
# G_{0,1} etc reduce; skip — the monomial block above is the requested unit test.
print("(K^2 closed-form cross-check skipped: phi(u)^2 moments non-closed in x; monomial table is the ordered test)")
