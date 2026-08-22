"""EXACT certificate that the 3D Ising layer transfer operator is not spectrally Gaussian.

`experiments/e37_spectral_gaussianity.py` establishes the result in high-precision floating
point.  This script removes the floating point entirely: every number below is an exact
rational, and every spectral statement is proved by Sylvester's law of inertia.

CONSTRUCTION (exact rationality).
Write t = tanh(K*/2) and choose t RATIONAL.  Then

    exp(K* A / 2) = prod_i (cosh(K*/2) + sinh(K*/2) X_i)
                  = (1 - t^2)^{-n/2} prod_i (1 + t X_i),

and  <k| prod_i (1 + t X_i) |l> = t^{hamming(k,l)},  an exact rational matrix P_t.
From K* = -(1/2) log tanh K we get exp(-2K*) = tanh K, hence

    exp(-2K) = tanh K* = 2t / (1 + t^2),      exp(2K) = (1 + t^2) / (2t) =: q,

which is rational.  The diagonal of B has entries b_k = sum over bonds of sigma_i sigma_j,
all of the same parity eps as the bond count, so exp(K b_k) = exp(K eps) * q^{(b_k-eps)/2}
with the first factor a GLOBAL scalar.  Therefore

    S := exp(K* A/2) exp(K B) exp(K* A/2) = (global scalar) * R,
    R := P_t diag(q^{(b_k - eps)/2}) P_t,

with R an exact rational symmetric positive definite matrix.  A global positive scalar
shifts every log-eigenvalue equally, so it cancels from the shifted log-spectrum that the
Gaussianity test examines.  Hence testing R is equivalent to testing V, exactly.

CERTIFICATION (Sylvester).
For a rational shift sigma, exact symmetric Gaussian elimination of R - sigma I yields a
diagonal D with the same inertia; the number of negative entries of D is EXACTLY the number
of eigenvalues of R below sigma.  No floating point, no error term.  From this:

  * bisection gives certified rational enclosures of individual eigenvalues;
  * and, crucially, ABSENCE of any eigenvalue in an interval is certified by two counts:
    if #{lambda < sigma_lo} == #{lambda < sigma_hi} then no eigenvalue lies in
    [sigma_lo, sigma_hi).

The second fact is what makes the certificate cheap: to refute Gaussianity we must show the
predicted eigenvalue is missing, and that is two inertia counts.
"""

from __future__ import annotations

import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ising.transfer_matrix import layer_bonds  # noqa: E402


# --------------------------------------------------------------- exact rational matrix
def build_R(n: int, bonds, t: Fraction):
    """Exact rational symmetric positive definite R, similar to the layer transfer operator
    up to a global positive scalar.  Returns (R, q, eps)."""
    q = (1 + t * t) / (2 * t)                      # = exp(2K)
    dim = 1 << n
    # diagonal of B
    b = []
    for k in range(dim):
        spins = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        b.append(sum(spins[i] * spins[j] for i, j in bonds))
    eps = b[0] % 2
    assert all((x - eps) % 2 == 0 for x in b), "bond-count parity assumption violated"
    d = [q ** ((x - eps) // 2) for x in b]

    # P_t[k][l] = t^hamming(k,l)
    tp = [t ** h for h in range(n + 1)]
    P = [[tp[bin(k ^ l).count("1")] for l in range(dim)] for k in range(dim)]

    # R = P diag(d) P   (P symmetric)
    R = [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        Pi = P[i]
        for j in range(i, dim):
            Pj = P[j]
            acc = Fraction(0)
            for k in range(dim):
                acc += Pi[k] * d[k] * Pj[k]
            R[i][j] = R[j][i] = acc
    return R, q, eps


# ------------------------------------------------------------------- exact inertia count
def count_at(R, sigma: Fraction):
    """Exact number of eigenvalues of R strictly below sigma, or None if sigma hits a
    degenerate pivot.

    Symmetric Gaussian elimination on R - sigma I.  By Sylvester's law of inertia the number
    of negative pivots equals #{lambda < sigma}.  If no zero pivot occurs then sigma is not
    an eigenvalue and the count is exact.  Returning None rather than silently perturbing is
    essential: the caller must know WHICH shift a count belongs to.
    """
    n = len(R)
    a = [[R[i][j] - (sigma if i == j else 0) for j in range(n)] for i in range(n)]
    neg = 0
    for k in range(n):
        p = a[k][k]
        if p == 0:
            return None
        if p < 0:
            neg += 1
        inv = 1 / p
        row = a[k]
        for i in range(k + 1, n):
            f = a[i][k] * inv
            if f:
                ai = a[i]
                for j in range(k, n):
                    ai[j] -= f * row[j]
    return neg


def count_inside(R, lo: Fraction, hi: Fraction):
    """A shift strictly inside (lo, hi) together with its exact count.

    Never reports a count at a shift other than the one returned, which is what makes the
    bisection sound.  Degenerate pivots are escaped by moving the trial shift, and the
    RETURNED shift is the one that was actually evaluated.
    """
    span = hi - lo
    for k in range(40):
        # k = 0 is the midpoint; later trials walk inward deterministically
        mid = lo + span / 2 if k == 0 else lo + span * Fraction(1, 2) + span * Fraction(
            (-1) ** k, 3 ** (k + 2))
        if not (lo < mid < hi):
            continue
        c = count_at(R, mid)
        if c is not None:
            return mid, c
    raise RuntimeError("no non-degenerate shift found inside the bracket")


def count_outward(R, sigma: Fraction, direction: int):
    """Exact count at sigma, or at a shift moved OUTWARD (direction -1 = down, +1 = up).

    Used for the final absence window.  Moving outward can only enlarge the tested interval,
    so an 'absent' verdict obtained this way remains valid for the original interval.
    Returns (shift_used, count).
    """
    c = count_at(R, sigma)
    if c is not None:
        return sigma, c
    step = abs(sigma) if sigma != 0 else Fraction(1)
    for k in range(1, 40):
        s = sigma + direction * step * Fraction(1, 10 ** (20 - k // 2))
        c = count_at(R, s)
        if c is not None:
            return s, c
    raise RuntimeError("no non-degenerate outward shift found")


def isolate(R, index: int, lo: Fraction, hi: Fraction, rel: Fraction):
    """Certified rational enclosure of the (0-based) index-th smallest eigenvalue.

    Invariant maintained: #{lambda < lo} <= index  and  #{lambda < hi} > index, with both
    counts taken at exactly lo and hi.  Every update uses the count at the shift that was
    actually evaluated, so the enclosure is valid with no perturbation slack.
    """
    while True:
        if lo > 0 and hi - lo <= rel * lo:
            return lo, hi
        mid, c = count_inside(R, lo, hi)
        if c <= index:
            lo = mid
        else:
            hi = mid


# ------------------------------------------------------------------------------- driver
def certify(name: str, n: int, bonds, t: Fraction, rel_exp: int = 8):
    """No logarithms anywhere: the Gaussian condition in MULTIPLICATIVE form.

    For a Gaussian, spec = { lambda_min * prod_{k in S} u_k } over subsets S, with u_k > 1.
    Sorting ascending forces lambda_1/lambda_0 = u_(1) and lambda_2/lambda_0 = u_(2), the
    two smallest single-particle factors -- because u_(1) u_(2) > u_(2) makes any other
    assignment impossible.  Since {u_(1), u_(2)} is itself a subset, the value

        lambda_0 * u_(1) * u_(2) = lambda_1 * lambda_2 / lambda_0

    MUST be an eigenvalue.  Every quantity here is a ratio of eigenvalues, so certified
    rational enclosures of the three smallest eigenvalues give a certified rational
    enclosure of the predicted one -- and two inertia counts decide whether it is there.
    """
    t0 = time.time()
    R, q, eps = build_R(n, bonds, t)
    dim = 1 << n
    rel = Fraction(1, 10 ** rel_exp)

    tr = sum(R[i][i] for i in range(dim))
    assert count_at(R, Fraction(0)) == 0, "R must be positive definite"
    encl = [isolate(R, i, Fraction(0), tr + 1, rel) for i in range(3)]
    (l0, h0), (l1, h1), (l2, h2) = encl

    # CERTIFIED DISTINCTNESS of the three lowest eigenvalues.
    # The forced-eigenvalue argument needs lambda_0 < lambda_1 < lambda_2, which excludes
    # zero modes (one zero mode gives lambda_0 = lambda_1, NOT all three equal) and a
    # repeated smallest factor.  This is weaker than full spectral simplicity, which is
    # neither assumed nor established.
    three_lowest_distinct = bool(h0 < l1 and h1 < l2)
    gap01 = float(l1 - h0) / float(l1) if three_lowest_distinct else 0.0
    gap12 = float(l2 - h1) / float(l2) if three_lowest_distinct else 0.0
    assert three_lowest_distinct, (
        "CERTIFICATE INVALID: the three smallest eigenvalue enclosures overlap. A zero mode "
        "(lambda_0 = lambda_1) or a repeated smallest factor (lambda_1 = lambda_2) would make "
        "the forced eigenvalue vacuous. Narrow the enclosures or reformulate.")

    # interval arithmetic on positive rationals: outward-rounded by construction
    pred_lo = l1 * l2 / h0
    pred_hi = h1 * h2 / l0
    s_lo, c_lo = count_outward(R, pred_lo, -1)
    s_hi, c_hi = count_outward(R, pred_hi, +1)
    absent = (c_lo == c_hi)

    return {
        "name": name,
        "n_sites": n,
        "n_bonds": len(bonds),
        "t_tanh_half_Kstar": str(t),
        "exp_2K": str(q),
        "bond_parity": eps,
        "dimension": dim,
        "relative_enclosure_width": f"1e-{rel_exp}",
        "lambda0": [float(l0), float(h0)],
        "lambda1": [float(l1), float(h1)],
        "lambda2": [float(l2), float(h2)],
        "three_lowest_distinct": three_lowest_distinct,
        "gap01_relative": gap01,
        "gap12_relative": gap12,
        "predicted_eigenvalue_interval": [float(pred_lo), float(pred_hi)],
        "predicted_interval_relative_width": float((pred_hi - pred_lo) / pred_lo),
        "window_actually_tested": [str(s_lo), str(s_hi)],
        "window_encloses_prediction": bool(s_lo <= pred_lo and s_hi >= pred_hi),
        "count_below_lo": c_lo,
        "count_below_hi": c_hi,
        "predicted_eigenvalue_absent": bool(absent),
        "elapsed_seconds": round(time.time() - t0, 2),
    }


def main() -> int:
    # tanh(K*/2) = 1/3  ->  exp(2K) = (1+1/9)/(2/3) = 5/3, a benign high-temperature point
    t = Fraction(1, 3)

    cases = [
        ("1D open chain n=4 (2D Ising strip, CONTROL)", 4, [(0, 1), (1, 2), (2, 3)]),
        ("1D open chain n=6 (2D Ising strip, CONTROL)", 6,
         [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]),
        ("2D open grid 2x3 (3D Ising layer)", 6, list(layer_bonds((2, 3), (False, False)))),
    ]

    rows = []
    for name, n, bonds in cases:
        print(f"certifying {name} ...", flush=True)
        r = certify(name, n, bonds, t)
        rows.append(r)
        print(f"   exp(2K) = {r['exp_2K']}   dim = 2^{n}   #bonds = {r['n_bonds']}")
        print(f"   predicted eigenvalue in [{r['predicted_eigenvalue_interval'][0]:.12g}, "
              f"{r['predicted_eigenvalue_interval'][1]:.12g}]  "
              f"(rel width {r['predicted_interval_relative_width']:.2e})")
        print(f"   inertia counts at the window ends: {r['count_below_lo']} and "
              f"{r['count_below_hi']}")
        print(f"   predicted eigenvalue ABSENT: {r['predicted_eigenvalue_absent']}   "
              f"({r['elapsed_seconds']} s)")

    ctrl = rows[0]
    ctrl2 = rows[1]
    test = rows[2]
    checks = [
        {
            "name": "control_chain_window_contains_eigenvalue",
            "passed": not (ctrl["predicted_eigenvalue_absent"]
                           or ctrl2["predicted_eigenvalue_absent"]),
            "detail": "for the 2D Ising strip the certified forced-value window contains at "
                      "least one eigenvalue: the inertia counts differ "
                      f"({ctrl['count_below_lo']} vs {ctrl['count_below_hi']}). This is a "
                      "consistency check ONLY: differing counts locate SOME eigenvalue in a "
                      "positive-width window; exact presence of the forced value itself is "
                      "not certified (that would need an equality certificate).",
        },
        {
            "name": "3D_layer_predicted_eigenvalue_is_certifiably_absent",
            "passed": test["predicted_eigenvalue_absent"],
            "detail": "for the 2x3 3D Ising layer the forced generator sum gamma_1 + gamma_2 "
                      "is NOT in the spectrum: the exact inertia count is identical at both "
                      f"ends of the window ({test['count_below_lo']}), so no eigenvalue lies "
                      "inside it. This is an exact rational certificate, with no floating "
                      "point in the spectral statement.",
        },
    ]

    out = {
        "provenance": {
            "script": "experiments/e38_gaussianity_certificate.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": "exact Fraction arithmetic throughout; Sylvester's law of inertia "
                         "for every spectral statement; no floating point anywhere in the "
                         "certificate. Counts are only ever reported at the shift actually "
                         "evaluated, and the final absence window is rounded OUTWARD.",
            "python": platform.python_version(),
        },
        "data": {"tanh_half_Kstar": str(t), "rows": rows},
        "checks": checks,
    }
    os.makedirs("results/spectral", exist_ok=True)
    json.dump(out, open("results/spectral/gaussianity_certificate.json", "w"), indent=1)

    print()
    for c in checks:
        print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['name']}")
        print(f"         {c['detail']}")
    bad = [c["name"] for c in checks if not c["passed"]]
    print()
    if bad:
        print(f"FAIL: {bad}")
        return 1
    print("PASS: exact rational certificate -- the 3D layer spectrum is missing an "
          "eigenvalue that any Gaussian operator on n modes would be forced to have.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
