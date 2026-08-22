"""Recurrence guesser + integer-relation detector for the zeta(5) campaign.

Two tools, with their evidence status built into the contract:

guess_recurrence(seq, order, degree)
    Exact search (FLINT integer nullspace) for polynomial-coefficient
    recurrences  sum_{i=0}^{order} C_i(n) * x_{n+i} = 0,  deg C_i <= degree.
    Any returned recurrence has been verified exactly at EVERY index the
    supplied sequence covers, so it is a proved identity for that range;
    validity for all n is CONDITIONAL until certified symbolically
    (creative telescoping or an algebraic proof).

find_integer_relation(values)
    mpmath PSLQ. Output is a CANDIDATE relation only -- never evidence.
    The returned vector merely makes the dot product small at the working
    precision; callers must certify independently.

Selftest (--selftest) validates the guesser against ground truth: fed only the
raw integers q_0..q_63 of Zudilin's sequence (arXiv:math/0206178v2), it must
re-derive recursion (1) exactly -- same primitive coefficient vector -- find
nothing at order 2 (degree <= 12) or order 3 with degree <= 8, and the PSLQ
wrapper must recover a constructed relation and reject pi.

Usage:  cd math && ./.venv/bin/python zeta5/src/recsearch.py --selftest
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from flint import fmpz_mat

sys.path.insert(0, str(Path(__file__).resolve().parent))


# --- exact integer polynomial helpers (coefficient lists, low degree first) --

def pmul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pshift(a: list[int], k: int) -> list[int]:
    """Coefficients of a(n+k)."""
    out = [0] * len(a)
    for j, c in enumerate(a):
        # c*(n+k)^j
        term = [c]
        for _ in range(j):
            term = pmul(term, [k, 1])
        for d, v in enumerate(term):
            out[d] += v
    return out


def peval(a: list[int], n: int) -> int:
    v = 0
    for c in reversed(a):
        v = v * n + c
    return v


def _primitive(vec: list[int]) -> list[int]:
    g = 0
    for v in vec:
        g = math.gcd(g, v)
    vec = [v // g for v in vec]
    lead = next(v for v in reversed(vec) if v != 0)
    return [-v for v in vec] if lead < 0 else vec


# --- recurrence guessing ------------------------------------------------------

def guess_recurrence(seq: list[Fraction | int], order: int, degree: int,
                     start: int = 1, extra: int = 8) -> list[list[list[int]]]:
    """Search for sum_i C_i(n) x_{n+i} = 0 over n = start .. len(seq)-1-order.

    Returns a list of candidates (possibly empty). Each candidate is
    [C_0, ..., C_order], integer coefficient lists (low degree first),
    primitively normalized with the leading coefficient of C_order positive.
    Every returned candidate is verified exactly on the FULL supplied range,
    not just the fitting window.
    """
    unknowns = (order + 1) * (degree + 1)
    rows = unknowns + extra
    last_fit = start + rows - 1
    if last_fit + order > len(seq) - 1 - 4:
        raise ValueError(f"need seq up to index {last_fit + order + 4}, "
                         f"have {len(seq) - 1}")
    seq = [Fraction(x) for x in seq]

    matrix: list[list[int]] = []
    for n in range(start, start + rows):
        row_frac = [Fraction(n**j) * seq[n + i]
                    for i in range(order + 1) for j in range(degree + 1)]
        den = math.lcm(*(f.denominator for f in row_frac))
        matrix.append([int(f * den) for f in row_frac])

    X, nullity = fmpz_mat(matrix).nullspace()
    out = []
    for k in range(nullity):
        vec = _primitive([int(X[r, k]) for r in range(unknowns)])
        cand = [vec[i * (degree + 1):(i + 1) * (degree + 1)]
                for i in range(order + 1)]
        full = all(sum(Fraction(peval(cand[i], n)) * seq[n + i]
                       for i in range(order + 1)) == 0
                   for n in range(start, len(seq) - order))
        if full:
            out.append(cand)
    return out


# --- integer-relation candidates ---------------------------------------------

def find_integer_relation(values: list[mp.mpf], maxcoeff: int = 10**6,
                          ) -> list[int] | None:
    """PSLQ candidate for sum_i c_i * values[i] = 0, or None.

    CANDIDATE ONLY. A non-None result is re-checked to be small at working
    precision and within the coefficient bound; that is still not a proof.
    """
    rel = mp.pslq(values, maxcoeff=maxcoeff, maxsteps=20000)
    if rel is None:
        return None
    if max(abs(c) for c in rel) > maxcoeff:
        return None
    resid = mp.fsum(c * v for c, v in zip(rel, values))
    if abs(resid) > mp.mpf(10) ** (-(mp.mp.dps // 2)):
        return None
    return rel


# --- selftest -----------------------------------------------------------------

def _target_recursion_shifted() -> list[list[int]]:
    """Recursion (1) rewritten as sum_{i=0}^{3} C_i(n) x_{n+i} = 0 (m = n+2).

    C_3(n) = (n+3)^6 a0(n+2), C_2(n) = a1(n+2), C_1(n) = -4(2n+3) a2(n+2),
    C_0(n) = -4(n+1)^4 (2n+3)(2n+1) a0(n+3);  from zudilin_rec's transcription.
    """
    a0 = [-2871, 20010, -48459, 41218]
    a1 = [2 * c for c in (-60291, -271701, 460056, 3790503, 1311365,
                          -19538418, -24317344, 36002654, 89030880, 48802112)]
    a2 = [-44541, 170716, 115771, -1182926, 647130,
          2947148, -3144314, -2617900, 3874492]

    def power(base: list[int], k: int) -> list[int]:
        out = [1]
        for _ in range(k):
            out = pmul(out, base)
        return out

    c3 = pmul(power([3, 1], 6), pshift(a0, 2))
    c2 = pshift(a1, 2)
    c1 = pmul([-4], pmul([3, 2], pshift(a2, 2)))
    c0 = pmul([-4], pmul(power([1, 1], 4),
                         pmul([3, 2], pmul([1, 2], pshift(a0, 3)))))
    deg = max(map(len, (c0, c1, c2, c3)))
    padded = [c + [0] * (deg - len(c)) for c in (c0, c1, c2, c3)]
    flat = _primitive([v for c in padded for v in c])
    return [flat[i * deg:(i + 1) * deg] for i in range(4)]


def selftest() -> int:
    from zudilin_rec import generate

    failures: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        tag = "PASS" if ok else "FAIL"
        print(f"[{tag}] {name}" + (f" -- {detail}" if detail else ""))
        if not ok:
            failures.append(name)

    q, _, _ = generate(63)

    # 1. blind re-derivation of recursion (1) from raw q_n data
    found = guess_recurrence(q, order=3, degree=9)
    check("guesser: unique order-3 degree-9 recurrence on q_n",
          len(found) == 1, f"candidates={len(found)}")
    if len(found) == 1:
        check("guesser: recovered recursion (1) exactly (primitive vector)",
              found[0] == _target_recursion_shifted())

    # 2. negative controls: nothing smaller on this sequence
    check("guesser: no order-2 recurrence, degree <= 12",
          guess_recurrence(q, order=2, degree=12) == [])
    check("guesser: no order-3 recurrence, degree <= 8",
          guess_recurrence(q, order=3, degree=8) == [])

    # 3. PSLQ wrapper: recover a constructed relation, reject pi
    mp.mp.dps = 60
    z3, z5 = mp.zeta(3), mp.zeta(5)
    x = 9 * z5 + 33 * z3 - 49          # r_1 from the paper, by construction
    rel = find_integer_relation([x, mp.mpf(1), z3, z5])
    check("pslq: recovers [1, 49, -33, -9]",
          rel in ([1, 49, -33, -9], [-1, -49, 33, 9]), f"rel={rel}")
    check("pslq: no small relation for pi over (1, zeta3, zeta5)",
          find_integer_relation([mp.pi, mp.mpf(1), z3, z5],
                                maxcoeff=10**4) is None)

    if failures:
        print(f"FAILURES: {failures}")
        return 1
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    print(__doc__)
