"""Closed forms and a rigorous audit of Liu's proposed tail domination.

The script keeps analytic, sampled-numerical, and conjectural statements visibly
separate.  It does not touch certificate campaigns.

Run from ``math`` with::

    ./.venv/bin/python uc/liu_tail.py
"""

import os

# Keep the requested spectral diagnostic on one CPU core.  These must be set
# before importing NumPy/BLAS.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from fractions import Fraction
from math import comb

import mpmath as mp
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss

VERIFY_TERMS = 1024
VERIFY_TOL = mp.mpf("1e-12")
GENERALIZED_EIGEN_CUTOFF = 1e-12


def _xlog1m_mp(x: mp.mpf) -> mp.mpf:
    """Return (1-x) log(1-x), continuously extended at x=1."""
    return mp.mpf("0") if x == 1 else (1 - x) * mp.log1p(-x)


def base_sum_mp(x: mp.mpf) -> mp.mpf:
    """sum_{n>=2} x^n/(n(n-1)), in closed form."""
    return x + _xlog1m_mp(x)


def odd_sum_mp(x: mp.mpf) -> mp.mpf:
    """sum over odd n>=3 of x^n/(n(n-1)), in closed form."""
    return x + (_xlog1m_mp(x) - (1 + x) * mp.log1p(x)) / 2


def even_sum_mp(x: mp.mpf) -> mp.mpf:
    """sum over even n>=2 of x^n/(n(n-1)), in closed form."""
    return (_xlog1m_mp(x) + (1 + x) * mp.log1p(x)) / 2


def _xlog1m_np(x: np.ndarray) -> np.ndarray:
    """Vectorised continuous extension of (1-x) log(1-x)."""
    out = np.zeros_like(x, dtype=float)
    mask = x < 1.0
    out[mask] = (1.0 - x[mask]) * np.log1p(-x[mask])
    return out


def base_sum_np(x: np.ndarray) -> np.ndarray:
    return x + _xlog1m_np(x)


def odd_sum_np(x: np.ndarray) -> np.ndarray:
    return x + (_xlog1m_np(x) - (1.0 + x) * np.log1p(x)) / 2.0


def even_sum_np(x: np.ndarray) -> np.ndarray:
    return (_xlog1m_np(x) + (1.0 + x) * np.log1p(x)) / 2.0


def exact_symbolic_checks() -> None:
    """Machine-check the algebra used in the three resummations."""
    y = sp.symbols("y")
    j = sp.symbols("j", integer=True, positive=True)
    k = sp.symbols("k", integer=True, nonnegative=True)
    i = sp.symbols("i", integer=True)

    base = y + (1 - y) * sp.log(1 - y)
    base_at_minus_y = -y + (1 + y) * sp.log(1 + y)
    odd = y + ((1 - y) * sp.log(1 - y) - (1 + y) * sp.log(1 + y)) / 2
    even = ((1 - y) * sp.log(1 - y) + (1 + y) * sp.log(1 + y)) / 2

    # SymPy returns Piecewise because it also records the convergence domain.
    summed = sp.summation(y**j / (j * (j - 1)), (j, 2, sp.oo))
    summed_on_unit_disk = summed.args[0][0] if isinstance(summed, sp.Piecewise) else summed
    assert sp.simplify(summed_on_unit_disk - base) == 0
    assert sp.simplify(odd - (base - base_at_minus_y) / 2) == 0
    assert sp.simplify(even - (base + base_at_minus_y) / 2) == 0
    assert sp.simplify(odd + even - base) == 0

    # Check the parity coefficients exactly through a nontrivial finite order;
    # the all-orders statement already follows from F(y) +/- F(-y).
    odd_series = sp.series(odd, y, 0, 18).removeO().expand()
    even_series = sp.series(even, y, 0, 18).removeO().expand()
    for degree in range(2, 18):
        coefficient = sp.Rational(1, degree * (degree - 1))
        assert odd_series.coeff(y, degree) == (coefficient if degree % 2 else 0)
        assert even_series.coeff(y, degree) == (coefficient if degree % 2 == 0 else 0)

    binomial_sum = sp.summation(sp.binomial(k, i) * y**i, (i, 0, k))
    binomial_on_unit_interval = (
        binomial_sum.args[0][0]
        if isinstance(binomial_sum, sp.Piecewise)
        else binomial_sum
    )
    assert sp.simplify(binomial_on_unit_interval - (1 + y) ** k) == 0

    # Exact identities behind the chosen Cauchy--Schwarz weights.
    m = sp.symbols("m", integer=True, nonnegative=True)
    assert sp.simplify(
        1 / ((m + 1) * (m + 2)) - (1 / (m + 1) - 1 / (m + 2))
    ) == 0


def verify_closed_forms() -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """Compare each closed kernel with a long, direct truncated series."""
    mp.mp.dps = 80
    samples = (
        (mp.mpf("0.02"), mp.mpf("0.03")),
        (mp.mpf("0.07"), mp.mpf("0.11")),
        (mp.mpf("0.23"), mp.mpf("0.81")),
        (mp.mpf("0.61"), mp.mpf("0.74")),
        (mp.mpf("0.91"), mp.mpf("0.94")),
    )
    positive_errors: list[mp.mpf] = []
    even_errors: list[mp.mpf] = []
    k_family_errors: list[mp.mpf] = []

    for s, t in samples:
        u = 1 - s
        v = 1 - t
        y = s * t
        z = u * v * (1 + y)

        positive_closed = u * v * odd_sum_mp(y)
        positive_direct = u * v * mp.fsum(
            y**j / (j * (j - 1)) for j in range(3, VERIFY_TERMS + 1, 2)
        )
        positive_errors.append(abs(positive_closed - positive_direct))

        even_closed = u * v * even_sum_mp(y)
        even_direct = u * v * mp.fsum(
            y**j / (j * (j - 1)) for j in range(2, VERIFY_TERMS + 1, 2)
        )
        even_errors.append(abs(even_closed - even_direct))

        k_family_closed = base_sum_mp(z)
        k_family_direct = mp.fsum(
            z**k / (k * (k - 1)) for k in range(2, VERIFY_TERMS + 1)
        )
        k_family_errors.append(abs(k_family_closed - k_family_direct))

    errors = max(positive_errors), max(even_errors), max(k_family_errors)
    assert all(error < VERIFY_TOL for error in errors)
    return errors


def cs_weight(m: int) -> Fraction:
    """A positive, summable weight with an exact telescoping total of one."""
    return Fraction(1, (m + 1) * (m + 2))


def cs_coefficient(j: int, m: int) -> Fraction:
    """Coefficient supplied by weighted Cauchy--Schwarz for positive index j."""
    return Fraction((m + 1) * (m + 2), j * (j - 1))


def available_g2_weight(i: int) -> Fraction:
    """Actual k=2 coefficient binom(2,i)/(2*1), extended by zero."""
    return Fraction(comb(2, i), 2) if 0 <= i <= 2 else Fraction(0)


def first_positive_index(j0: int) -> int:
    """Smallest odd j>=3 contained in a tail starting at j0."""
    j = max(3, j0)
    return j if j % 2 else j + 1


def exact_cs_checks() -> None:
    """Check the explicit Cauchy--Schwarz constants with rational arithmetic."""
    for m in range(128):
        assert cs_weight(m) * cs_coefficient(3, m) == Fraction(1, 3 * 2)
    partial_weight = sum((cs_weight(m) for m in range(128)), Fraction(0))
    assert partial_weight == 1 - Fraction(1, 129)
    assert [available_g2_weight(i) for i in range(6)] == [
        Fraction(1, 2),
        Fraction(1),
        Fraction(1, 2),
        Fraction(0),
        Fraction(0),
        Fraction(0),
    ]


def projected_nystrom(order: int, j0: int = 3) -> dict[str, float | int]:
    """Sample the projected L2 operators using Gauss--Legendre quadrature."""
    nodes, weights = leggauss(order)
    s = (nodes + 1.0) / 2.0
    weights = weights / 2.0
    left = s[:, None]
    right = s[None, :]
    uv = (1.0 - left) * (1.0 - right)
    y = left * right
    z = uv * (1.0 + y)

    positive_scalar = odd_sum_np(y)
    for j in range(3, j0, 2):
        positive_scalar -= y**j / (j * (j - 1))
    positive = uv * positive_scalar
    negative = uv * even_sum_np(y) + base_sum_np(z)

    root_weight = np.sqrt(weights)
    positive_matrix = root_weight[:, None] * positive * root_weight[None, :]
    negative_matrix = root_weight[:, None] * negative * root_weight[None, :]

    # In weighted coordinates, these are exactly the three L2 constraints.
    constraints = root_weight[:, None] * np.column_stack(
        (np.ones(order), s, s * (1.0 - s))
    )
    q, _ = np.linalg.qr(constraints, mode="complete")
    complement = q[:, 3:]
    p_projected = complement.T @ positive_matrix @ complement
    d_projected = complement.T @ negative_matrix @ complement

    p_eigenvalues = np.linalg.eigvalsh(p_projected)
    d_eigenvalues, d_eigenvectors = np.linalg.eigh(d_projected)
    gap_eigenvalues = np.linalg.eigvalsh(d_projected - p_projected)

    # This is only a conditioning-aware diagnostic.  Modes below the stated
    # cutoff are deliberately omitted, so the ratio is not a proof.
    resolved = d_eigenvalues > GENERALIZED_EIGEN_CUTOFF
    whitener = d_eigenvectors[:, resolved] / np.sqrt(d_eigenvalues[resolved])[None, :]
    relative = np.linalg.eigvalsh(whitener.T @ p_projected @ whitener)

    return {
        "order": order,
        "resolved": int(np.count_nonzero(resolved)),
        "positive_norm": float(p_eigenvalues[-1]),
        "negative_min": float(d_eigenvalues[0]),
        "gap_min": float(gap_eigenvalues[0]),
        "relative_max": float(relative[-1]),
    }


def _scientific(x: mp.mpf) -> str:
    return mp.nstr(x, 8, min_fixed=0, max_fixed=0)


def main() -> None:
    exact_symbolic_checks()
    exact_cs_checks()

    print("LIU POSITIVE-TAIL ANALYSIS")
    print("All status labels apply to the complete claim on the same line/block.")
    print()
    print("1. Closed kernels")
    print("   Put u=1-s, v=1-t, y=st, and z=uv(1+y).")
    print("   PROVED [SymPy exact infinite summation and simplification]:")
    print("     F(x) = sum_{n>=2} x^n/[n(n-1)]")
    print("          = x + (1-x) ln(1-x), with (1-x)ln(1-x)=0 at x=1.")
    print("   PROVED [exact parity split (F(y)-F(-y))/2]:")
    print("     P(s,t) = uv * { y + [(1-y)ln(1-y)-(1+y)ln(1+y)]/2 }.")
    print("     This is sum_{odd j>=3} f_j(s)f_j(t)/[j(j-1)].")
    print("   PROVED [exact parity split (F(y)+F(-y))/2]:")
    print("     N_even(s,t) = uv * {[(1-y)ln(1-y)+(1+y)ln(1+y)]/2}.")
    print("     The signed even-j contribution is -N_even(s,t).")
    print("   PROVED [SymPy binomial theorem check plus the exact F resummation]:")
    print("     N_k(s,t) = sum_{k>=2} sum_{i=0}^k binom(k,i)")
    print("                  * g_{k,i}(s)g_{k,i}(t)/[k(k-1)]")
    print("                = F(z) = z + (1-z)ln(1-z).")
    print("     The signed k>=2 contribution is -N_k(s,t).")

    p_error, even_error, k_error = verify_closed_forms()
    print(
        "   NUMERICAL [80-digit mpmath, five fixed interior pairs, "
        f"{VERIFY_TERMS} terms]:"
    )
    print(f"     max |P_closed-P_truncated|           = {_scientific(p_error)}")
    print(f"     max |N_even_closed-N_even_truncated| = {_scientific(even_error)}")
    print(f"     max |N_k_closed-N_k_truncated|       = {_scientific(k_error)}")
    print(f"     Each sampled error is below {float(VERIFY_TOL):.0e}.")

    print()
    print("2. Weighted Cauchy--Schwarz audit")
    print("   PROVED [uniform geometric telescoping, then Cauchy--Schwarz]:")
    print("     f_j = sum_{m>=0} g_{2,j+m}.")
    print("     Choose w_m = 1/[(m+1)(m+2)], so sum_m w_m = 1.")
    print("     For every odd j>=3 and every finite signed measure mu,")
    print("       <f_j>^2/[j(j-1)]")
    print("         <= sum_{m>=0} c_{j,m}<g_{2,j+m}>^2,")
    print("       c_{j,m} = (m+1)(m+2)/[j(j-1)].")
    print("     The coefficient identities are checked exactly by Fraction/SymPy.")
    demo = ", ".join(str(cs_coefficient(3, m)) for m in range(5))
    print(f"     For j=3, c_(3,m), m=0..4, are {demo}.")

    budgets = ", ".join(str(available_g2_weight(i)) for i in range(6))
    print("   PROVED [exact binomial support and Fraction arithmetic]:")
    print(f"     available k=2 weights at i=0..5 are [{budgets}].")
    print("     Thus every available k=2 weight at i>=3 is exactly zero.")
    print("     The explicit weights already demand c_(3,0)=1/3 at i=3:")
    print("       available 0, shortfall 1/3.")
    print("   PROVED [weight-independent obstruction]:")
    print("     For arbitrary positive weights W=sum_m w_m<infinity,")
    print("       c_(j,0)=W/[j(j-1)w_0] > 1/[j(j-1)],")
    print("     because w_0<W.  For a tail starting at any finite J0, let J be")
    print("     its first odd index >=3.  At i=J the available weight is zero,")
    print("     while the infimum required shortfall is 1/[J(J-1)] (not attained).")
    print("     At J=3 this best-possible infimum is 1/6.")
    print("   PROVED [the preceding exact obstruction]:")
    print("     No finite J0 closes this proposed g_2-budget domination.")
    print("     Only the vacuous empty tail J0=infinity closes.")
    print("     Consequently there is no finite leftover set for a matrix certificate;")
    print("     this route leaves every positive index {3,5,7,...} unresolved.")

    print()
    print("3. Spectral alternative")
    print("   PROVED [compact-operator theorem on projected L2(0,1)]:")
    print("     D=N_even+N_k has a continuous PSD kernel, hence is compact on the")
    print("     infinite-dimensional codimension-three subspace H={1,s,a}^perp.")
    print("     Therefore inf_{phi in H, ||phi||=1}<phi,D phi>=0; D has no")
    print("     positive smallest eigenvalue.  Every finite positive tail has")
    print("     positive norm because f_J is not in span{1,s,a}.  Hence the")
    print("     sufficient test ||P_tail|| <= lambda_min(D) cannot close for any")
    print("     finite J0.")
    print("   NUMERICAL [float64 Gauss--Legendre Nystrom discretisation]:")
    print("     order   ||P||_H       min eig(D)_H   min eig(D-P)_H resolved  max P/D")
    for order in (32, 48, 64):
        result = projected_nystrom(order)
        print(
            f"     {order:5d}  {result['positive_norm']:.12e} "
            f"{result['negative_min']:+.3e}       {result['gap_min']:+.3e}"
            f"       {result['resolved']:3d}   {result['relative_max']:.10f}"
        )
    print(
        "     Here max P/D retains only D-eigenmodes above "
        f"{GENERALIZED_EIGEN_CUTOFF:.0e}; it is a diagnostic, not a proof."
    )
    print("     The near-one ratios and roundoff-scale minima expose no spectral margin.")
    print("   NUMERICAL [64-node discretisation, closed-form positive tails]:")
    print("     J0      ||P_{j>=J0}||_H    max P_tail/D   min eig(D-P_tail)")
    for j0 in (3, 9, 17, 33):
        result = projected_nystrom(64, j0)
        print(
            f"     {j0:2d}      {result['positive_norm']:.12e}   "
            f"{result['relative_max']:.10f}     {result['gap_min']:+.3e}"
        )

    print()
    print("FINAL STATUS")
    print("   PROVED [exact support obstruction]: the requested Cauchy--Schwarz")
    print("   tail domination does not close for any finite cutoff.")
    print("   PROVED [compactness obstruction]: the separate spectral-gap criterion")
    print("   also cannot close; its negative-side lower edge is exactly zero.")
    print("   NUMERICAL [the tables above]: relative finite-dimensional spectra are")
    print("   essentially tight at one, so the experiment supplies no hidden margin.")
    print("   PROVED [exact support obstruction]: unresolved positive set is all")
    print("   odd j>=3; no finite leftover set exists.")


if __name__ == "__main__":
    main()
