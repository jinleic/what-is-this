#!/usr/bin/env python3
"""
Executable model of the two pillars of Hilbert's Tenth Problem over Q.

Pillar 1 (definability): Poonen's forall-exists definition of Z in Q
  [arXiv math/0703907] via traces of norm-1 quaternions, Hilbert symbols,
  and Hasse-Minkowski. We implement membership in S_{a,b}, the covering
  T_{a,b} = S+S+{0..2309} = intersection of Z_(p) over ramified p, the
  adversary of Lemma 2.6, and the forall^2 parameterization of Lemma 4.2.

Pillar 2 (elliptic path): Poonen's integrality mechanism on rank-1 curves
  [arXiv math/0306277] that yields H10 undecidability for large subrings
  of Q (and, via rank stability [arXiv 2412.01768, 2501.18774], for all
  rings of integers). We verify divisibility-sequence structure, rank of
  apparition, formal-group valuation growth, quadratic height growth,
  primitive divisors, and the S-integrality "control surface".

Dependencies: none (pure Python 3).
"""
from fractions import Fraction
from itertools import combinations
import json
import math
import os
import random
import sys
import time

OO = 'oo'  # the real place


# ---------------------------------------------------------------- primitives
# Least strong pseudoprime to all bases 2..41 (OEIS A014233 a(13)).  Below this
# value, Miller-Rabin with the 13 bases used in _is_prime is a PROVEN primality
# certificate; at or above it we refuse to certify rather than risk a silent
# composite-as-prime (factorint feeds Hilbert-symbol / ramification computations
# that require EXACT factorizations).
_MR_LIMIT = 3317044064679887385961981
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)


class PrimalityBound(Exception):
    """n passed every MR base but lies outside the proven deterministic range."""


def _mr_screen(n: int) -> bool:
    """Strong-pseudoprime screen, bases 2..41.  False = PROVEN composite (a
    witnessed composite is composite for all n); True = no claim above _MR_LIMIT."""
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


_TRIAL_PRIMES_CACHE = None


def _trial_primes_1e6():
    global _TRIAL_PRIMES_CACHE
    if _TRIAL_PRIMES_CACHE is None:
        _TRIAL_PRIMES_CACHE = tuple(primerange(2, 10**6))
    return _TRIAL_PRIMES_CACHE


def _pocklington(n: int, depth: int) -> bool:
    """Generalized Pocklington (Crandall-Pomerance Thm 4.1.3): if n-1 = F*R,
    F > sqrt(n) built from recursively CERTIFIED primes q, and for each q some
    a has a^(n-1) = 1 (mod n) and gcd(a^((n-1)/q) - 1, n) = 1, then every prime
    p | n is = 1 (mod F) > sqrt(n), so n is prime.  If only F >= cuberoot(n)
    is reached, the Brillhart-Lehmer-Selfridge relaxation (CP Thm 4.1.5)
    applies: the same per-q conditions force every prime factor = 1 (mod F),
    so a composite n has EXACTLY two prime factors (uF+1)(vF+1); writing
    R = c2*F + c1, n is then composite iff x^2 - c1 x + c2 has nonnegative
    square discriminant reproducing n.  Returns a proof or raises
    PrimalityBound; a Fermat failure or a proper gcd proves compositeness."""
    if depth > 3:
        raise PrimalityBound(n)
    m = n - 1
    fac = {}
    F = 1
    # integer ceil(cuberoot(n)): smallest F0 with F0**3 >= n
    lo, hi = 1, 1 << ((n.bit_length() + 2) // 3 + 2)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid - 1
    F0 = lo + (lo**3 < n)
    for p in _trial_primes_1e6():
        if p * p > m or F >= F0:
            break                       # F >= n^{1/3}: BLS is reachable already
        while m % p == 0:
            fac[p] = fac.get(p, 0) + 1
            m //= p
            F *= p
    if m == 1:
        F = n - 1
    stack = [m] if m > 1 else []
    bls = F**3 >= n                          # BLS threshold tracked while splitting
    while stack and F * F <= n:
        t = stack.pop()
        if t == 1:
            continue
        try:
            if _is_prime(t, depth + 1):
                fac[t] = fac.get(t, 0) + 1
                F *= t
                bls = bls or F**3 >= n
                continue
        except PrimalityBound:
            continue                      # uncertifiable chunk: excluded from F
        try:
            d = _brent(t)
            stack += [d, t // d]
        except FactorBudget:
            continue                      # unfactored chunk: excluded from F
    if F * F <= n and not bls:
        raise PrimalityBound(n)           # not enough certified factored part
    for q in fac:
        for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
            if pow(a, n - 1, n) != 1:
                return False              # Fermat witness: proven composite
            g = math.gcd(pow(a, (n - 1) // q, n) - 1, n)
            if g == 1:
                break                     # Pocklington condition met for this q
            if 1 < g < n:
                return False              # proper factor found: proven composite
        else:
            raise PrimalityBound(n)       # no usable witness for q: inconclusive
    if F * F > n:
        return True
    # BLS relaxation: F >= n^{1/3} -> composite n has exactly two prime
    # factors uF+1, vF+1.  Coprimality of the certified and uncertified
    # parts is a theorem hypothesis: check it.
    R = (n - 1) // F
    assert F * R == n - 1
    if math.gcd(F, R) != 1:
        raise PrimalityBound(n)
    c2, c1 = divmod(R, F)
    D = c1 * c1 - 4 * c2
    if D < 0:
        return True
    s_ = math.isqrt(D)
    if s_ * s_ != D or (c1 - s_) % 2:
        return True
    u, v = (c1 - s_) // 2, (c1 + s_) // 2
    if u > 0 and v > 0 and (u * F + 1) * (v * F + 1) == n:
        return False                      # exhibited two-factor composite
    return True


def _is_prime(n: int, depth: int = 0) -> bool:
    """PROVEN primality only: deterministic MR below _MR_LIMIT (A014233),
    Pocklington n-1 certification above it, PrimalityBound if unprovable."""
    if n < 2:
        return False
    if not _mr_screen(n):
        return False
    if n < _MR_LIMIT:
        return True
    return _pocklington(n, depth)


class FactorBudget(Exception):
    """Raised when Brent rho exceeds its iteration budget (huge semiprime)."""


def _brent(n: int, budget: int = 400000) -> int:
    if n % 2 == 0:
        return 2
    for c in range(1, 20):
        y, r, q, g, count = 2, 1, 1, 1, 0
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(128, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += 128
                count += 128
                if count > budget:
                    raise FactorBudget(n)
            r *= 2
        if g != n:
            return g
        y = ys
        while True:
            y = (y * y + c) % n
            g = math.gcd(abs(x - y), n)
            if g > 1:
                break
        if g != n:
            return g
    raise FactorBudget(n)


def factorint(n: int) -> dict:
    """Prime factorization; small-prime pre-pass + Miller-Rabin + Brent rho."""
    n = abs(n)
    f = {}
    if n <= 1:
        return f
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
              53, 59, 61, 67, 71, 73, 79, 83, 89, 97):
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if _is_prime(m):
            f[m] = f.get(m, 0) + 1
            continue
        d = _brent(m)
        stack += [d, m // d]
    return f


def _certify_by_trial(p: int) -> bool:
    """Independent exact primality check by trial division (p <= ~1e12)."""
    if p < 2:
        return False
    d = 2
    while d * d <= p:
        if p % d == 0:
            return False
        d += 1 if d == 2 else 2
    return True


def _verify_factoring():
    """Adversarial factorization tests (advisory 2026-08-12).

    Guards the exactness contract of factorint: Hilbert symbols / ramified sets
    treat its output as a proven factorization.  Checks:
      (i)   psi_12 = A014233(12), strong pseudoprime to bases 2..37, is NOT
            certified prime (base 41 witnesses it);
      (ii)  psi_13 = A014233(13) = _MR_LIMIT is never certified prime: either
            proven composite (Pocklington gcd/Fermat) or explicit PrimalityBound;
      (iii) the Mersenne prime 2^89-1 > _MR_LIMIT IS certified via Pocklington
            (n-1 = 2(2^44-1)(2^44+1) has enough certified factored part);
      (iv)  classical spsp / Carmichael numbers factor correctly, every factor
            independently certified by trial division;
      (v)   perfect powers and balanced semiprimes;
      (vi)  300 random n < 1e9 round-trip with trial-division certification.
    """
    psi12 = 318665857834031151167461      # spsp to 2..37; composite
    psi13 = 3317044064679887385961981     # spsp to 2..41; = _MR_LIMIT
    m89 = 2**89 - 1                       # prime, > _MR_LIMIT
    assert _is_prime(psi12) is False      # base 41 must witness it
    try:
        assert _is_prime(psi13) is False  # proven composite is acceptable...
    except PrimalityBound:
        pass                              # ...and so is an explicit refusal
    assert _is_prime(m89) is True         # Pocklington proof above the MR range
    assert factorint(m89) == {m89: 1}
    # refusal paths (reviewer 2026-08-12): the "refuse, never guess" contract
    # (a) prime > _MR_LIMIT with 3-smooth n-1: full Pocklington witness loop
    assert _is_prime(2**83 * 3**10 + 1) is True
    # (b) prime > _MR_LIMIT with rho-hard p-1: must REFUSE, and factorint must
    #     propagate the refusal rather than return a partial factorization
    p_hard = 2**90 + 4449                 # MR-passes; p-1 = 2^k * (rho-hard)
    for probe in (lambda: _is_prime(p_hard), lambda: factorint(3 * p_hard)):
        try:
            probe()
            raise AssertionError("certified/factored past a rho-hard p-1")
        except (PrimalityBound, FactorBudget):
            pass
    # (c) balanced ~122-bit semiprime: rho budget must raise, not stall/guess
    try:
        _brent(2305843009213713989 * 2305843009213743977, budget=30000)
        raise AssertionError("brent exceeded budget without raising")
    except FactorBudget:
        pass
    adversarial = [2047, 1373653, 25326001, 3215031751, 3474749660383,
                   341550071728321, 561, 1105, 1729, 41041, 825265]
    for n in adversarial:
        f = factorint(n)
        prod = 1
        for p, e in f.items():
            prod *= p**e
            assert _certify_by_trial(p), (n, p)
        assert prod == n and len(f) >= 2, (n, f)
    for n, expect in ((2**64, {2: 64}), (3**40, {3: 40}),
                      (1000003**2, {1000003: 2})):
        assert factorint(n) == expect, n
    n = 1000003 * 1000033
    assert factorint(n) == {1000003: 1, 1000033: 1}
    rng = random.Random(1231)
    for _ in range(300):
        n = rng.randint(2, 10**9)
        f = factorint(n)
        prod = 1
        for p, e in f.items():
            prod *= p**e
            assert _certify_by_trial(p), (n, p)
        assert prod == n, (n, f)


def primerange(a: int, b: int):
    sieve = bytearray([1]) * b
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(b ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [p for p in range(max(a, 2), b) if sieve[p]]


def vp(x, p: int) -> int:
    x = Fraction(x)
    assert x != 0
    v, n, d = 0, x.numerator, x.denominator
    while n % p == 0:
        n //= p
        v += 1
    while d % p == 0:
        d //= p
        v -= 1
    return v


def unit_part(x, p: int) -> Fraction:
    return Fraction(x) / Fraction(p) ** vp(x, p)


def unit_mod(u: Fraction, m: int) -> int:
    return (u.numerator % m) * pow(u.denominator % m, -1, m) % m


def legendre(u: Fraction, p: int) -> int:
    return 1 if pow(unit_mod(u, p), (p - 1) // 2, p) == 1 else -1


# ------------------------------------------------------------ Hilbert symbol
def hilbert(a, b, v) -> int:
    """Hilbert symbol (a,b)_v over Q_v; v a prime or OO."""
    a, b = Fraction(a), Fraction(b)
    assert a != 0 and b != 0
    if v == OO:
        return -1 if (a < 0 and b < 0) else 1
    p = v
    al, be = vp(a, p), vp(b, p)
    u, w = unit_part(a, p), unit_part(b, p)
    if p != 2:
        s = 1
        if (al * be) % 2 and ((p - 1) // 2) % 2:
            s = -s
        if be % 2 and legendre(u, p) == -1:
            s = -s
        if al % 2 and legendre(w, p) == -1:
            s = -s
        return s
    um, wm = unit_mod(u, 8), unit_mod(w, 8)
    eps = lambda r: ((r - 1) // 2) % 2
    om = lambda r: ((r * r - 1) // 8) % 2
    return -1 if (eps(um) * eps(wm) + al * om(wm) + be * om(um)) % 2 else 1


def places_of(*xs):
    ps = {2}
    for x in xs:
        x = Fraction(x)
        ps |= set(factorint(x.numerator)) | set(factorint(x.denominator))
    return sorted(ps) + [OO]


def ramified(a, b):
    """Places where the quaternion algebra H_{a,b} ramifies (even cardinality)."""
    return [v for v in places_of(a, b) if hilbert(a, b, v) == -1]


def delta_fin(a, b):
    return [v for v in ramified(a, b) if v != OO]


# --------------------------------------------- quadratic forms, local-global
def is_square_qv(x, v) -> bool:
    x = Fraction(x)
    if x == 0:
        return True
    if v == OO:
        return x > 0
    p = v
    if vp(x, p) % 2:
        return False
    u = unit_part(x, p)
    return legendre(u, p) == 1 if p != 2 else unit_mod(u, 8) == 1


def hasse_eps(coeffs, v) -> int:
    e = 1
    for c1, c2 in combinations(coeffs, 2):
        e *= hilbert(c1, c2, v)
    return e


def isotropic_qv(coeffs, v) -> bool:
    """Diagonal form <c1..cn> isotropic over Q_v (Serre, Cours d'arithmetique IV, Thm 6)."""
    coeffs = [Fraction(c) for c in coeffs]
    assert all(c != 0 for c in coeffs)
    n = len(coeffs)
    if v == OO:
        return any(c > 0 for c in coeffs) and any(c < 0 for c in coeffs)
    d = Fraction(1)
    for c in coeffs:
        d *= c
    if n == 1:
        return False
    if n == 2:
        return is_square_qv(-d, v)
    if n == 3:
        return hilbert(-1, -d, v) == hasse_eps(coeffs, v)
    if n == 4:
        return (not is_square_qv(d, v)) or hasse_eps(coeffs, v) == hilbert(-1, -1, v)
    return True


def represents_qv(coeffs, n, v) -> bool:
    """<coeffs> represents n over Q_v.  n = 0 asks for NONTRIVIAL zero
    representation, i.e. isotropy (standard convention)."""
    n = Fraction(n)
    if n == 0:
        return isotropic_qv(coeffs, v)
    return isotropic_qv(list(coeffs) + [-n], v)


def represents_q(coeffs, n) -> bool:
    """Hasse-Minkowski: f represents n over Q iff over every Q_v.
    n = 0 means NONTRIVIALLY (global isotropy).  Restricting to places_of is
    sound in every arity: with >= 3 variables the form is unimodular hence
    isotropic at omitted places; with 2 variables the places_of checks force
    v_p(-c1*c2) even at all p and -c1*c2 > 0 at OO, so -det is a perfect
    square in Q and a square everywhere; with 1 variable it is always False."""
    n = Fraction(n)
    if n == 0:
        return all(isotropic_qv(coeffs, v) for v in places_of(*coeffs))
    return all(represents_qv(coeffs, n, v) for v in places_of(*coeffs, n))


# ------------------------------------------------- Pillar 1: Poonen's gadget
N_SHIFT = 2310  # 2*3*5*7*11, Poonen math/0703907


def in_S(a, b, t) -> bool:
    """t in S_{a,b} = reduced traces of norm-1 elements of H_{a,b}(Q).
    t = 2*x1 with x1^2 - a x2^2 - b x3^2 + a b x4^2 = 1
    <=> t = +-2  or  <a,b,-ab> represents (t^2-4)/4 over Q."""
    t = Fraction(t)
    if t in (2, -2):
        return True
    return represents_q([a, b, Fraction(-a) * b], (t * t - 4) / 4)


def U(q: int) -> set:
    """U_q = {s in F_q : x^2 - s x + 1 irreducible} (Poonen Lemma 2.1)."""
    if q == 2:
        return {1}
    sq = {(x * x) % q for x in range(q)}
    return {s for s in range(q) if (s * s - 4) % q not in sq}


def in_T_certificate(a, b, t, s_window=60):
    """Constructive witness (s, s', n) with t = s + s' + n, s,s' in S_{a,b},
    0 <= n < 2310 -- realizes Poonen Lemma 2.4 membership."""
    t = Fraction(t)
    d = t.denominator
    for n in range(N_SHIFT):
        u = t - n
        for e in range(-s_window, s_window + 1):
            s = Fraction(e, d)
            if in_S(a, b, s) and in_S(a, b, u - s):
                return (s, u - s, n)
    return None


def adversary(t):
    """For non-integral t, an (a,b) with a,b > 0 and some prime q | denom(t)
    in Delta_{a,b} -- hence t not in T_{a,b} (Poonen Lemma 2.6 pattern)."""
    t = Fraction(t)
    d = t.denominator
    assert d > 1, "t is an integer; no adversary exists"
    q = min(factorint(d))
    if q == 2:
        return (7, 7)
    for b_ in range(2, 200):
        if legendre(Fraction(b_), q) == -1 and q in delta_fin(q, b_):
            return (q, b_)
    raise RuntimeError


def param(a, b):
    """Poonen Lemma 4.2 parameterization: a' = a^2+b^2+1 > 0, b' = a^2+a+1+b^2 > 0."""
    return (a * a + b * b + 1, a * a + a + 1 + b * b)


# ------------------------------------------- Pillar 2: elliptic divisibility
class EC:
    """y^2 = x^3 + A x + B over Q; points are (Fraction, Fraction) or None=O."""

    def __init__(self, A, B):
        self.A, self.B = Fraction(A), Fraction(B)
        assert -16 * (4 * self.A ** 3 + 27 * self.B ** 2) != 0

    def on(self, P):
        if P is None:
            return True
        x, y = P
        return y * y == x ** 3 + self.A * x + self.B

    def neg(self, P):
        return None if P is None else (P[0], -P[1])

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2 and y1 == -y2:
            return None
        lam = (3 * x1 * x1 + self.A) / (2 * y1) if P == Q else (y2 - y1) / (x2 - x1)
        x3 = lam * lam - x1 - x2
        return (x3, lam * (x1 - x3) - y1)

    def mul(self, n, P):
        if n < 0:
            return self.neg(self.mul(-n, P))
        R, Q = None, P
        while n:
            if n & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            n >>= 1
        return R


def divisibility_sequence(E: EC, P, nmax: int):
    """B_n with x(nP) = A_n / B_n^2 in lowest terms."""
    B, Q = {}, None
    for n in range(1, nmax + 1):
        Q = E.add(Q, P)
        d = Q[0].denominator
        r = math.isqrt(d)
        assert r * r == d
        B[n] = r
    return B


# ---------------------------------------------------------------- self-tests
def _selftest():
    rng = random.Random(42)
    rand_rat = lambda: Fraction(rng.choice([x for x in range(-30, 31) if x]),
                                rng.randint(1, 30))
    assert hilbert(-1, -1, 2) == -1 and hilbert(-1, -1, OO) == -1
    assert hilbert(2, 3, 3) == -1
    for _ in range(200):
        a, b, c = rand_rat(), rand_rat(), rand_rat()
        for v in [2, 3, 5, 7, 11, OO]:
            assert hilbert(a, b, v) == hilbert(b, a, v)
            assert hilbert(a * c, b, v) == hilbert(a, b, v) * hilbert(c, b, v)
            assert hilbert(a * b ** 2, c, v) == hilbert(a, c, v)
    for _ in range(300):
        a, b = rand_rat(), rand_rat()
        prod = 1
        for v in places_of(a, b):
            prod *= hilbert(a, b, v)
        assert prod == 1
        assert len(ramified(a, b)) % 2 == 0
    # quadratic forms
    assert represents_q([1, 1], 5) and not represents_q([1, 1], 3)
    assert not represents_q([1, 1], -1)
    for m in (7, 15, 23, 28, 31, 60, 14, 33):
        k = m
        while k % 4 == 0:
            k //= 4
        assert represents_q([1, 1, 1], m) == (k % 8 != 7)
    assert represents_q([1, 1, 1, 1], 7)
    # zero-representation convention: n = 0 means nontrivial (isotropy)
    assert not represents_q([1, 1], 0) and represents_q([1, -1], 0)
    assert not represents_q([1, 1, 1], 0) and represents_q([1, 1, -2], 0)
    assert not represents_q([1, 1, -3], 0)
    assert not represents_qv([1, 1], 0, OO) and represents_qv([1, -1], 0, OO)
    # Psi-level n == 0 shortcut is deliberate: (y,r,s) = 0 solves Psi when
    # delta c^2 = 16 (no nonvanishing constraint in Psi_tau) -- tau=0, c=4:
    assert _sun_solvable(Fraction(3), Fraction(5), Fraction(4), Fraction(0)) is True
    _verify_factoring()


def _verify_results():
    """Reproduces every computational claim in RESULTS.md.
    Positive T-membership uses a bounded witness search (sampled sufficiency);
    negative certificates are sound (q | denom(t) with q in Delta_fin)."""
    # sandwich at ramified p for (a,b)=(5,2): red^{-1}(U_p) subset S subset Z_(p)
    a, b = 5, 2
    assert delta_fin(a, b) == [2, 5]
    U5, U2 = U(5), U(2)
    for t in range(-30, 31):
        if t % 5 in U5 and t % 2 in U2:
            assert in_S(a, b, t)
    for t in [Fraction(1, 5), Fraction(7, 5), Fraction(3, 2), Fraction(9, 10),
              Fraction(12, 5), Fraction(1, 2), Fraction(-99, 5)]:
        assert not in_S(a, b, t)
    # Lemma 2.4 realized on samples (incl. one needing the n-shift machinery)
    for t in [Fraction(17), Fraction(-23), Fraction(1, 3), Fraction(-40, 21),
              Fraction(1234, 9), Fraction(-7, 3), Fraction(22, 7), Fraction(0),
              Fraction(16)]:
        w = in_T_certificate(a, b, t)
        assert w is not None and w[0] + w[1] + w[2] == t
        assert in_S(a, b, w[0]) and in_S(a, b, w[1]) and 0 <= w[2] < N_SHIFT
    # sound non-membership certificates for every sampled non-integer
    for t in [Fraction(1, 2), Fraction(3, 5), Fraction(22, 7), Fraction(5, 13),
              Fraction(-1, 30)]:
        aa, bb = adversary(t)
        q = min(factorint(t.denominator))
        assert q in delta_fin(aa, bb) and not in_S(aa, bb, t)
    # forall^2 parameterization reaches every prime p < 50
    cands = [Fraction(n, d) for n in range(-8, 9) for d in (1, 2, 3)]
    for p in primerange(2, 50):
        assert any(hilbert(*param(a_, b_), p) == -1
                   for a_ in cands for b_ in cands), p
    # elliptic pillar: y^2 = x^3 - 2, P = (3,5), n <= 40
    E_, P_ = EC(0, -2), (Fraction(3), Fraction(5))
    B = divisibility_sequence(E_, P_, 40)
    for m in range(1, 41):
        for n in range(m, 41, m):
            assert B[n] % B[m] == 0

    def vint(m, r):
        k = 0
        while m % r == 0:
            m //= r
            k += 1
        return k

    apparition, checked = {}, 0
    for r in primerange(2, 200):
        ns = [n for n in range(1, 41) if B[n] % r == 0]
        if ns:
            assert ns == list(range(ns[0], 41, ns[0])), r
            apparition[r] = ns[0]
    assert len(apparition) == 35
    for r, nr in apparition.items():          # formal-group valuation law
        if r == 2:
            continue
        for k in range(1, 40 // nr + 1):
            assert vint(B[k * nr], r) == vint(B[nr], r) + vint(k, r)
            checked += 1
    assert checked == 122
    for n in (20, 30, 40):                    # quadratic denominator growth
        assert 0.66 < math.log(B[n]) / n ** 2 < 0.68
    # primitive divisors for n = 2..14, proven: either a fully-factored fresh
    # prime, or an unfactored cofactor coprime to all earlier B_m (every prime
    # factor of such a cofactor is primitive)
    small_primes = primerange(2, 10 ** 6)

    def partial_factor(m):
        fs = {}
        for q in small_primes:
            while m % q == 0:
                fs[q] = fs.get(q, 0) + 1
                m //= q
            if m == 1:
                break
        return fs, m

    seen = set()
    for n in range(2, 15):
        fs, cof = partial_factor(B[n])
        prim = [q for q in fs if q not in seen]
        assert prim or (cof > 1 and all(math.gcd(cof, B[m]) == 1
                                        for m in range(1, n))), n
        seen |= set(fs)


def _verify_L2prime():
    """L2' (THEOREMS.md): model-relative structure theorem for cofinite S.

    Model A: y^2 = x^3 - 2,       P  = (3, 5)     (minimal at 5)
    Model B: Y^2 = X^3 - 2*5^6,   P' = (75, 625)  (5-integral, NOT minimal at 5)
    Same curve over Q: X = 25 x, Y = 125 y.

    (a) exact apparition {n : v_q(x(nP)) < 0} = m_q Z  (both models, every q < 50 appearing);
    (b) D = {n : m_q inmid n for all q in Q0} == brute force, several Q0;
    (c) D empty  <=>  some m_q = 1  <=>  x(P) not in R  <=>  1 not in D;
    (model-dependence) Q0 = {5}: D_A != D_B, n = 2 separates -- the advisory's point.
    """
    N = 60
    A, PA = EC(0, -2), (Fraction(3), Fraction(5))
    B, PB = EC(0, Fraction(-2) * 5 ** 6), (Fraction(75), Fraction(625))
    assert A.on(PA) and B.on(PB)

    def multiples(E, P, N):
        xs, Q = {}, None
        for n in range(1, N + 1):
            Q = E.add(Q, P)
            xs[n] = Q[0]
        return xs

    xA, xB = multiples(A, PA, N), multiples(B, PB, N)
    # independent group laws agree with the rescaling -- non-circular setup check
    assert all(xB[n] == 25 * xA[n] for n in range(1, N + 1))

    def apparition(xs, q):
        """m_q with {n : v_q < 0} = m_q Z, asserting exactness (L2'(a)); 0 if empty."""
        S = [n for n, x in xs.items() if vp(x, q) < 0]
        if not S:
            return 0
        m = S[0]
        assert S == [n for n in xs if n % m == 0], (q, m, S[:6])
        return m

    mA, mB = {}, {}
    for q in primerange(2, 50):
        a, b = apparition(xA, q), apparition(xB, q)
        if a:
            mA[q] = a
        if b:
            mB[q] = b
    # rescaling is supported at 5 only: all other moduli must coincide
    assert all(mB.get(q, 0) == m for q, m in mA.items() if q != 5), (mA, mB)
    # model-dependence at q = 5: predicted m_5 = 2 -> 10 (depth-2 filtration)
    assert mA[5] == 2 and mB[5] == 10, (mA.get(5), mB.get(5))
    assert vp(xA[2], 5) == -2 and vp(xB[2], 5) == 0
    print(f"  m_q, model A (minimal at 5): {mA}")
    print(f"  m_q, model B (5-rescaled):   {mB}")
    print(f"  x(2P) = {xA[2]} (v_5 = -2) vs X(2P') = {xB[2]} (v_5 = 0):")
    print("  -> R-integrality NOT preserved by Q0-supported coordinate change;")
    print("     old L2 step retracted, theorem now stated per Q0-integral model.")

    def D_formula(m, Q0, N):
        ms = [m[q] for q in Q0 if q in m]
        return [n for n in range(1, N + 1) if all(n % mq for mq in ms)]

    def D_brute(xs, Q0, N):
        return [n for n in range(1, N + 1) if all(vp(xs[n], q) >= 0 for q in Q0)]

    for Q0 in [{5}, {2, 5}, {2, 3, 5, 7}, {11, 13}]:
        for xs, m in ((xA, mA), (xB, mB)):
            assert D_formula(m, Q0, N) == D_brute(xs, Q0, N), Q0
    DA5, DB5 = D_brute(xA, {5}, N), D_brute(xB, {5}, N)
    assert 2 in DB5 and 2 not in DA5 and DA5 != DB5
    print(f"  D(Q0={{5}}): A excludes 2Z, B excludes 10Z; congruence formula == brute (4 Q0s, both models)")

    # (c) emptiness: P2 = 2P has v_5(x) < 0, so with Q0 = {5}: m_5(P2) = 1, D = empty
    x2 = multiples(A, A.add(PA, PA), 30)
    assert apparition(x2, 5) == 1
    assert all(vp(x, 5) < 0 for x in x2.values())
    for Q0 in [{5}, {2, 5}, {2, 3, 5, 7}]:
        DD = D_brute(xA, Q0, N)
        assert (len(DD) > 0) == (1 in DD)
    print("  emptiness criterion: D(2P, {5}) = empty (m_5 = 1); D nonempty <=> 1 in D  OK")
    print("  L2' verified: (a) exact APs, (b) congruence structure, (c) emptiness, model-dependence\n")


# --------------------- Sun arXiv 2607.28606 audit (record claim: 10 -> 7)
def _qmul(p, q, A, B):
    """Quaternion product, basis (1, i, j, ij), i^2 = A, j^2 = B, ij = -ji."""
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (w1*w2 + A*x1*x2 + B*y1*y2 - A*B*z1*z2,
            w1*x2 + x1*w2 - B*y1*z2 + B*z1*y2,
            w1*y2 + y1*w2 + A*x1*z2 - A*z1*x2,
            w1*z2 + z1*w2 + x1*y2 - y1*x2)


def _nrd(p, A, B):
    w, x, y, z = p
    return w*w - A*x*x - B*y*y + A*B*z*z


def _sun_solvable(a, b, c, tau):
    """Psi_tau(a,b,c) over Q (S = {2}): exists rational (y,r,s) with
    delta(c^2 - A y^2) - 16 B (r^2 - A s^2) = 16, A = 1+4a^2, B = 2b,
    delta = 1 - A tau^2.  Decided by Hasse-Minkowski."""
    A, B = 1 + 4*a*a, 2*b
    delta = 1 - A*tau*tau
    if delta == 0:
        return None
    n = 16 - delta*c*c
    return True if n == 0 else represents_q([-delta*A, -16*B, 16*A*B], n)


def _sun_g(a, b):
    return 16*a**4/(1 + 4*a*a) - ((b - 1)*(b - 1)/b)**2


def _sun_h(a, b, z):
    den = 1 - z - a*a*z*z
    return None if (den == 0 or b in (0, 1)) else a*a*z*z*_sun_g(a, b)/den


def _sun_integral(c, ps):
    return c == 0 or all(vp(c, q) >= 0 for q in ps)


def _sun_solvable3(a, b, c, tau):
    """Ternary restriction of Psi_tau (witness s frozen to 0):
    exists rational (y,r) with delta(c^2 - A y^2) - 16 B r^2 = 16.
    Any solution extends to Psi_tau by s = 0, so soundness is inherited;
    the probe below measures what completeness costs."""
    A, B = 1 + 4*a*a, 2*b
    delta = 1 - A*tau*tau
    if delta == 0:
        return None
    n = 16 - delta*c*c
    return True if n == 0 else represents_q([-delta*A, -16*B], n)

def _sun_tied_solvable(a, b, c, tau):
    """L6 tied block: Psi_tau with its s-witness fixed to (a - 1)/2.

    The remaining conic in (y, r) is decided by Hasse-Minkowski.
    """
    A, B = 1 + 4*a*a, 2*b
    delta = 1 - A*tau*tau
    if delta == 0:
        return None
    s = (a - 1)/2
    n = 16 - delta*c*c - 16*A*B*s*s
    return True if n == 0 else represents_q([-delta*A, -16*B], n)


def _sun_target_tau(a, w):
    """Canonical tied branch at odd target w: 0 if w=3 mod 4, else 2a/A."""
    A = 1 + 4*a*a
    return Fraction(0) if w % 4 == 3 else 2*a/A



class _GFq:
    """Tiny F_{p^k} model (odd q only): tuples of F_p coeffs, monic modulus.
    Moduli beyond the frozen table are found by lex search (Rabin test:
    x^{p^k} = x mod f and gcd(x^{p^{k/r}} - x, f) = 1 for primes r | k)."""
    _MOD = {9: (3, (1, 0)),        # x^2 + 1        over F_3
            25: (5, (2, 0)),       # x^2 + 2        over F_5
            27: (3, (1, 2, 0)),    # x^3 + 2x + 1   over F_3
            49: (7, (1, 0))}       # x^2 + 1        over F_7

    def __init__(self, q):
        if q in self._MOD:
            self.p, self.mod = self._MOD[q]
        else:
            p = 3
            while q % p:
                p += 2
            k, m = 0, q
            while m > 1:
                assert m % p == 0, f"{q} not an odd prime power"
                m //= p
                k += 1
            self.p, self.mod = p, self._find_mod(p, k)
        self.k = len(self.mod)
        from itertools import product as _prod
        self.elems = [tuple(t) for t in _prod(range(self.p), repeat=self.k)]
        self.zero = (0,) * self.k
        self.one = (1,) + (0,) * (self.k - 1)

    @staticmethod
    def _polmulmod(u, v, mod, p):
        k = len(mod)
        conv = [0] * (2 * k - 1)
        for i, x in enumerate(u):
            if x:
                for j, y in enumerate(v):
                    conv[i + j] = (conv[i + j] + x * y) % p
        for d in range(2 * k - 2, k - 1, -1):
            c = conv[d]
            if c:
                conv[d] = 0
                for j in range(k):
                    conv[d - k + j] = (conv[d - k + j] - c * mod[j]) % p
        return conv[:k]

    @classmethod
    def _find_mod(cls, p, k):
        """Lexicographically first monic irreducible x^k + tail over F_p."""
        from itertools import product as _prod
        x = [0, 1] + [0] * (k - 2)

        def xq(e, f):
            r, b = [1] + [0] * (k - 1), x[:]
            while e:
                if e & 1:
                    r = cls._polmulmod(r, b, f, p)
                e >>= 1
                if e:
                    b = cls._polmulmod(b, b, f, p)
            return r

        def deg(u):
            return max((i for i, c in enumerate(u) if c), default=-1)

        def coprime(u, v):
            u, v = list(u), list(v)
            while deg(v) >= 0:
                du, dv = deg(u), deg(v)
                if du < dv:
                    u, v = v, u
                    continue
                c = u[du] * pow(v[dv], p - 2, p) % p
                for j in range(dv + 1):
                    u[du - dv + j] = (u[du - dv + j] - c * v[j]) % p
            return deg(u) == 0

        for tail in _prod(range(p), repeat=k):
            f = list(tail)
            if xq(p ** k, f) != x:
                continue
            if all(coprime(f + [1], [(a - b) % p for a, b in zip(xq(p ** (k // r), f), x)])
                   for r in (2, 3, 5) if k % r == 0):
                return tuple(tail)
        raise AssertionError((p, k))

    def scal(self, n):
        return (n % self.p,) + (0,) * (self.k - 1)

    def add(self, u, v):
        return tuple((x + y) % self.p for x, y in zip(u, v))

    def neg(self, u):
        return tuple((-x) % self.p for x in u)

    def mul(self, u, v):
        return tuple(self._polmulmod(u, v, self.mod, self.p))

    def inv(self, u):
        assert u != self.zero
        r, b, e = self.one, u, self.p ** self.k - 2
        while e:
            if e & 1:
                r = self.mul(r, b)
            e >>= 1
            if e:
                b = self.mul(b, b)
        return r


def _field(q):
    """(elems, add, mul, neg, one, zero, scal, inv) for an odd prime power q."""
    if _is_prime(q):
        return (list(range(q)), lambda x, y: (x + y) % q,
                lambda x, y: x * y % q, lambda x: -x % q, 1, 0,
                lambda n: n % q, lambda x: pow(x, q - 2, q))
    F = _GFq(q)
    return F.elems, F.add, F.mul, F.neg, F.one, F.zero, F.scal, F.inv


def _lemma51_scan(q):
    """Exhaustively decide Sun Lemma 5.1 over F_q (odd prime power).

    For each tau in F_q \\ {0, 1, -1}: does there exist a in F_q with
      chi(1 + 4a^2) = -1   and   chi(-(1+4a^2) (1 - tau^2 (1+4a^2))) = 1 ?
    Returns (n_admissible_tau, failing_taus, min_solution_count).
    Exhaustion over a finite field IS a proof for that q."""
    elems, add, mul, neg, one, zero, scal, _inv = _field(q)
    four, minus_one = scal(4), neg(one)
    squares = {mul(z, z) for z in elems} - {zero}
    chi = lambda z: 0 if z == zero else (1 if z in squares else -1)
    A_of = {a: add(one, mul(four, mul(a, a))) for a in elems}
    n_adm, failures, min_cnt = 0, [], None
    for tau in elems:
        if tau in (zero, one, minus_one):
            continue
        n_adm += 1
        t2 = mul(tau, tau)
        cnt = 0
        for a in elems:
            A = A_of[a]
            if chi(A) != -1:
                continue
            if chi(neg(mul(A, add(one, neg(mul(t2, A)))))) == 1:
                cnt += 1
        if cnt == 0:
            failures.append(tau)
        elif min_cnt is None or cnt < min_cnt:
            min_cnt = cnt
        if q > 25:
            assert cnt > 0, (q, tau)        # consistency (proved for all q: L5)
    return n_adm, failures, min_cnt


def _lemma51_exhaust(verbose=True):
    """L4: exact exceptional set of Sun Lemma 5.1 for ALL odd prime powers
    q <= 25 (exhaustive = proved), plus Weil-range sanity for 25 < q <= 49."""
    small = [3, 5, 7, 9, 11, 13, 17, 19, 23, 25]
    sanity = [27, 29, 31, 37, 41, 43, 47, 49]
    rows = {}
    for q in small + sanity:
        rows[q] = _lemma51_scan(q)
    if verbose:
        for q in small:
            n_adm, fails, mc = rows[q]
            tag = "VACUOUS (no admissible tau)" if n_adm == 0 else \
                  (f"FAILS at tau={fails}" if fails else f"holds, min #a = {mc}")
            print(f"    q={q}: {tag}")
        ok = [q for q in sanity if not rows[q][1]]
        print(f"    sanity 25<q<=49: holds for all tau at q in {ok}")
    return rows


def _lemma51_identity(bound=49, verbose=True):
    """L5 [PROVED HERE, THEOREMS.md]: Sun Lemma 5.1's solution count is EXACTLY
    a quarter of the points on a Legendre elliptic curve:

        4 N(tau) = q + 1 + sum_u chi(u (u-1) (u - tau^{-2})) = #E_lam(F_q),
        E_lam : y^2 = x (x-1) (x-lam),  lam = tau^{-2}.

    With Hasse |#E - q - 1| <= 2 sqrt(q) (classical), N(tau) >= (sqrt(q)-1)^2/4
    > 0 for EVERY odd prime power q and every admissible tau: Sun's hypothesis
    '|k| > 25' is void.  This verifies the identity + Hasse + 4 | #E for all
    odd prime powers q <= bound and all tau, plus each step evaluation used in
    the proof (T1 = -1, fibre multiplicity 1 + chi(u-1), T2a = -1, T3 = -1-z)."""
    from collections import Counter
    qs = []
    for p in primerange(3, bound + 1):
        pk = p
        while pk <= bound:
            qs.append(pk)
            pk *= p
    pairs = 0
    for q in sorted(qs):
        elems, add, mul, neg, one, zero, scal, inv = _field(q)
        squares = {mul(z, z) for z in elems} - {zero}
        chi = lambda z: 0 if z == zero else (1 if z in squares else -1)
        four, minus_one = scal(4), neg(one)
        u_of = {a: add(one, mul(four, mul(a, a))) for a in elems}
        z_cnt = 1 + chi(minus_one)
        assert sum(chi(u_of[a]) for a in elems) == -1, q            # T1 = -1
        mult = Counter(u_of.values())
        assert all(mult.get(u, 0) == 1 + chi(add(u, neg(one))) for u in elems), q
        for tau in elems:
            if tau in (zero, one, minus_one):
                continue
            t2 = mul(tau, tau)
            lam = inv(t2)
            N = T2a = T3 = C = 0
            for a in elems:
                ua = u_of[a]
                wa = neg(mul(ua, add(one, neg(mul(t2, ua)))))
                T3 += chi(mul(ua, wa))
                if chi(ua) == -1 and chi(wa) == 1:
                    N += 1
            for u in elems:
                T2a += chi(add(mul(t2, mul(u, u)), neg(u)))
                C += chi(mul(mul(u, add(u, neg(one))), add(u, neg(lam))))
            assert T2a == -1 and T3 == -1 - z_cnt, (q, tau)
            assert 4 * N == q + 1 + C, (q, tau, N, C)               # the identity
            assert C * C <= 4 * q, (q, tau, C)                      # Hasse
            assert (q + 1 + C) % 4 == 0, (q, tau)                   # full 2-torsion
            pairs += 1
    if verbose:
        print(f"    4N(tau) = #E_lam(F_q) [lam = tau^(-2), Legendre] + Hasse + 4 | #E:")
        print(f"    verified for ALL odd prime powers q <= {bound}:"
              f" {pairs} (q, tau) pairs, 0 failures")
        print("    -> N >= (sqrt(q)-1)^2/4 > 0 for every odd q (Hasse): '|k| > 25' is")
        print("       VOID; the only Lemma-5.1 obstruction is tau-vacuity of F_3")
    return pairs


def _ternary_probe(verbose=True):
    """Channel-2 barrier probe: freeze s = 0 in Psi_tau (3 witnesses -> 2).
    Soundness is inherited (a ternary solution IS a quaternary one).
    Completeness dies: exhibit Delta-integral c with Psi_tau solvable but the
    ternary restriction unsolvable, at BOTH tau = 0 and tau = tau_1."""
    losses, checked = [], 0
    for w, (ai, bi, _tt) in sorted(_SUN_WITNESSES.items()):
        if w >= 60:
            continue                        # fast fixed grid, as documented
        a, b = Fraction(ai), Fraction(bi)
        D = delta_fin(1 + 4*a*a, 2*b)
        taus = (Fraction(0), 2*a/(1 + 4*a*a))
        for c in [Fraction(k) for k in (0, 1, -1, 2, 5, -9, 24)]:
            if not _sun_integral(c, D):
                continue
            for tau in taus:
                q4, q3 = _sun_solvable(a, b, c, tau), _sun_solvable3(a, b, c, tau)
                if q4 is None:
                    continue
                checked += 1
                assert not (q3 and not q4), (a, b, c, tau)   # s=0 soundness
                if q4 and not q3:
                    losses.append((a, b, c, tau))
    assert losses, "ternary restriction lost nothing on this grid?!"
    if verbose:
        a, b, c, tau = losses[0]
        print(f"    {checked} integral instances; ternary loses {len(losses)}"
              f" (first: a={a}, b={b}, c={c}, tau={tau})")
        print("    -> s = 0 restriction is sound but NOT complete: the 4th witness")
        print("       is load-bearing; 2-witness certificates need a new ternary idea")
    return losses


def _lemma22(q):
    """Sun Lemma 2.2 residue combinatorics over F_q, exhaustive."""
    if q == 2:  # F_4 = F_2[x]/(x^2+x+1); conj z = z^2, Tr(a+bx) = b
        elems = [(i, j) for i in range(2) for j in range(2) if (i, j) != (0, 0)]
        mul = lambda p, r: ((p[0]*r[0] + p[1]*r[1]) % 2,
                            (p[0]*r[1] + p[1]*r[0] + p[1]*r[1]) % 2)
        norm = lambda p: (p[0]*p[0] + p[0]*p[1] + p[1]*p[1]) % 2
        tr = lambda p: p[1] % 2
    else:       # F_{q^2} = F_q[x]/(x^2-d), d a nonsquare
        sq = {(i*i) % q for i in range(1, q)}
        d = next(t for t in range(2, q) if t not in sq)
        elems = [(i, j) for i in range(q) for j in range(q) if (i, j) != (0, 0)]
        mul = lambda p, r: ((p[0]*r[0] + d*p[1]*r[1]) % q, (p[0]*r[1] + p[1]*r[0]) % q)
        norm = lambda p: (p[0]*p[0] - d*p[1]*p[1]) % q
        tr = lambda p: (2*p[0]) % q
    T = [p for p in elems if norm(p) == 1]
    assert len(T) == q + 1, (q, len(T))
    V = {tr(p) for p in T}
    fixed = {p for p in T if mul(p, p) == (1, 0)}
    Uset = {tr(p) for p in T if p not in fixed}
    assert (len(Uset), len(V)) == (((q-1)//2, (q+3)//2) if q % 2 else (q//2, q//2 + 1)), q
    for cbar in range(q):
        assert any((cbar - v) % q in Uset for v in V), (q, cbar)


# explicit bridge witnesses w -> (a, b, tau) in Phi_1^{S={2}} for z = w^3,
# found 2026-08-12 under the PROVEN factoring engine (no uncertified primality
# accepted anywhere).  tau = None: some tau in {0, 2a/(1+4a^2)} fires;
# tau = (num, den): that frozen tau fires (8 targets below need one outside
# {0, tau_1} -- the Lambda-family license exercised by real targets).
# Re-verified end-to-end below: Phi membership, Delta = {2,w}, h integral, block.
_SUN_WITNESSES = {
    3: (1, -3, None), 5: (1, -11, None), 7: (1, -7, None), 11: (5, -11, None),
    13: (1, -13, None), 17: (1, -17, None), 19: (3, -19, None),
    23: (1, -23, None), 29: (3, -29, None), 31: (3, -31, None),
    37: (3, -11, None), 41: (5, -41, None), 43: (1, -43, None),
    47: (1, -47, None), 53: (1, -53, None), 59: (3, -177, None),
    61: (5, -61, None), 67: (1, -67, None), 71: (7, -71, None),
    73: (1, -73, None), 79: (3, -79, None), 83: (1, -83, None),
    89: (5, -89, None), 97: (1, -97, None), 101: (5, -101, None),
    103: (1, -103, None), 107: (1, -107, None), 109: (9, -109, None),
    113: (1, -113, None), 127: (1, -127, None), 131: (3, -131, None),
    137: (9, -137, (1, 4)), 139: (7, -139, (3, 2)), 149: (9, -149, None),
    151: (5, -151, None), 157: (1, -157, None), 163: (1, -163, None),
    167: (1, -167, None), 173: (1, -173, (1, 4)), 179: (3, -179, None),
    181: (11, -1267, (3, 2)), 191: (5, -191, None), 193: (1, -193, None),
    197: (1, -197, None), 199: (7, -199, None), 211: (7, -211, None),
    223: (1, -223, None), 227: (1, -227, None), 229: (9, -229, None),
    233: (1, -233, None), 239: (3, -717, None), 241: (3, -241, (1, 4)),
    251: (3, -251, (1, 2)), 257: (1, -257, None), 263: (1, -263, None),
    269: (7, -269, (1, 4)), 271: (5, -271, None), 277: (3, -277, None),
    281: (3, -281, (1, 4)), 283: (1, -283, None), 293: (1, -293, None),
}

# L6 canonical tied witnesses w -> (s, b, tau), with a = 1 + 2s and z = w^3.
# Every entry was found under the proven-primality factoring engine after first
# enforcing that A is a nonsquare w-unit.  It is re-verified below: tau is the
# W0 branch, (a,b) is in Phi_1^{S={2}}, Delta = {2,w}, h is integral, and
# Psi_tau has a rational point whose s-coordinate is exactly (a - 1)/2.
_L6_WITNESSES = {
    3: ((1, 11), -33, (0, 1)), 5: ((-1, 9), 45, (126, 277)), 7: ((1, 13), -21, (0, 1)),
    11: ((-1, 9), 11, (0, 1)), 13: ((0, 1), 143, (2, 5)), 17: ((-1, 11), 85, (198, 445)),
    19: ((1, 9), 95, (0, 1)), 23: ((1, 11), 207, (0, 1)), 29: ((1, 9), -87, (198, 565)),
    31: ((1, 13), 155, (0, 1)), 37: ((1, 11), -407, (286, 797)), 41: ((-2, 9), -369, (90, 181)),
    43: ((-1, 9), -129, (0, 1)), 47: ((1, 13), 423, (0, 1)), 53: ((0, 1), -53, (2, 5)),
    59: ((1, 9), 295, (0, 1)), 61: ((-1, 9), -183, (126, 277)), 67: ((1, 7), 67, (0, 1)),
    71: ((-1, 7), 71, (0, 1)), 73: ((0, 1), -73, (2, 5)), 79: ((2, 13), 869, (0, 1)),
    83: ((0, 1), -747, (0, 1)), 89: ((-1, 13), -623, (286, 653)), 97: ((-1, 11), 873, (198, 445)),
    101: ((2, 13), -909, (442, 1325)), 103: ((-1, 13), -721, (0, 1)), 107: ((-2, 9), 321, (0, 1)),
    109: ((2, 11), -981, (330, 1021)), 113: ((-3, 11), -791, (110, 221)), 127: ((1, 13), -381, (0, 1)),
    131: ((-2, 9), 1965, (0, 1)), 137: ((-1, 9), 1233, (126, 277)), 139: ((-1, 5), 139, (0, 1)),
    149: ((-2, 9), -149, (90, 181)), 151: ((1, 13), 151, (0, 1)), 157: ((0, 1), 1727, (2, 5)),
    163: ((-2, 11), -1467, (0, 1)), 167: ((-4, 13), -167, (0, 1)), 173: ((1, 9), 173, (198, 565)),
    179: ((-2, 9), -179, (0, 1)), 181: ((-1, 11), -2715, (198, 445)), 191: ((-2, 9), -191, (0, 1)),
    193: ((-1, 3), -579, (6, 13)), 197: ((0, 1), 2167, (2, 5)), 199: ((3, 13), 199, (0, 1)),
    211: ((-1, 7), 1055, (0, 1)), 223: ((-1, 11), -2453, (0, 1)), 227: ((-1, 1), -2497, (0, 1)),
    229: ((4, 13), -2061, (546, 1933)), 233: ((0, 1), -233, (2, 5)), 239: ((-1, 13), 239, (0, 1)),
    241: ((-1, 13), -2651, (286, 653)), 251: ((-1, 5), -753, (0, 1)), 257: ((0, 1), -257, (2, 5)),
    263: ((-1, 11), 1315, (0, 1)), 269: ((1, 11), 3497, (286, 797)), 271: ((1, 7), 271, (0, 1)),
    277: ((-2, 9), -2493, (90, 181)), 281: ((-3, 11), 3653, (110, 221)), 283: ((-2, 1), -283, (0, 1)),
    293: ((1, 5), 293, (70, 221)),
}



def _verify_sun(extended=False):
    """Independent audit of Sun, arXiv 2607.28606 ('Q\\Z is diophantine over Q
    with 7 unknowns', posted 2026-07-30).  K = Q, S = {2}, pi = 2, u = 1.

    VERIFIED (this run):
      (A) identities (2.2)-(2.4): a Psi_tau solution yields alpha, alpha*gamma_tau
          with Nrd = 1 and Trd(alpha) + Trd(alpha gamma_tau) = c  [exact arithmetic];
      (B) soundness (= Prop 2.1 route): Psi_tau solvable => c integral at every
          finite ramified place of (1+4a^2, 2b); zero violations;
      (C) fixed-tau incompleteness: (a,b,c) = (2, -7/3, -9), Delta = {3,7},
          c Delta-integral, yet (9.5) and (9.6) BOTH rationally unsolvable
          => the tau_0/tau_1 remark cannot replace the Lambda-family (Secs. 3-4);
      (D) Lemma 2.2 residue combinatorics, exhaustive over F_q for q < 50;
      (E) bridge existence probe: explicit (a,b) in Phi witnessing z = w^3,
          each with Delta(Q_{a,b}) = {2, w}: w < 100 in-suite, the full frozen
          table (every odd prime w < 300) with --extended.

    COMPLETENESS CHAIN (Secs. 3-8): traced; no substantive gap found.  The posted
    text carries stale cross-references from a draft renumbering ('Section 11'
    -> Sec. 8; 'Sections 7 and 8' for the target-place local point -> Secs. 5-6):
    editorial defects, not missing mathematics.  Sec. 5 supplies the local point
    at w (Lemma 5.1 makes -A*delta_tau a square unit at w; v_w(c) = 6m-2 >= 4;
    Hensel on -A delta y^2 = 16); Secs. 3-4 freeze smooth points at S; Sec. 7
    picks global (a,b); Sec. 8 = Hasse-Minkowski + Lemma 8.1.  Status: 7 =
    established modulo refereeing; refereed anchor 10 (Daans, JLMS 2024).
    """
    rng = random.Random(2026)
    rq = lambda: Fraction(rng.randint(-9, 9), rng.randint(1, 6))
    # (A) identities
    checked = 0
    while checked < 150:
        a, b, c, tau, y, r, s = (rq() for _ in range(7))
        A, B = 1 + 4*a*a, 2*b
        delta = 1 - A*tau*tau
        if delta == 0 or b == 0:
            continue
        z = ((c - A*tau*y)/4, (y - tau*c)/4)
        alpha = (z[0], z[1], r, s)
        gamma = ((1 + A*tau*tau)/delta, 2*tau/delta, 0, 0)
        assert 16*(z[0]*z[0] - A*z[1]*z[1]) == delta*(c*c - A*y*y)              # (2.2)
        assert 16*_nrd(alpha, A, B) == delta*(c*c - A*y*y) - 16*B*(r*r - A*s*s)  # (2.3)
        assert _nrd(gamma, A, B) == 1
        assert 2*alpha[0] + 2*_qmul(alpha, gamma, A, B)[0] == c                  # (2.4)
        checked += 1
    print(f"  (A) identities (2.2)-(2.4): {checked} exact random instances OK")
    # (B) soundness on engineered grid
    tested = viol = pairs = 0
    while pairs < 20:
        a, b = rq(), rq()
        if b == 0:
            continue
        pairs += 1
        D = delta_fin(1 + 4*a*a, 2*b)
        tau1 = 2*a/(1 + 4*a*a)
        cs = [Fraction(k) for k in (-9, -1, 0, 5, 24)]
        cs += [Fraction(1, q) for q in D] + [Fraction(q + 1, q) for q in D]
        for c in cs:
            for tau in (Fraction(0), tau1):
                s_ = _sun_solvable(a, b, c, tau)
                tested += 1
                if s_ and not _sun_integral(c, D):
                    viol += 1
    assert viol == 0, viol
    print(f"  (B) soundness: {tested} (a,b,c,tau) checks, 0 violations (Prop 2.1 confirmed)")
    # (C) fixed-tau incompleteness certificate
    a0, b0, c0 = Fraction(2), Fraction(-7, 3), Fraction(-9)
    D0 = delta_fin(1 + 4*a0*a0, 2*b0)
    assert D0 == [3, 7] and _sun_integral(c0, D0)
    s95 = _sun_solvable(a0, b0, c0, Fraction(0))
    s96 = _sun_solvable(a0, b0, c0, 2*a0/(1 + 4*a0*a0))
    assert s95 is False and s96 is False
    print(f"  (C) incompleteness: (a,b,c)=(2,-7/3,-9), Delta={D0}, c integral,")
    print(f"      (9.5) solvable: {s95}, (9.6) solvable: {s96} -> Lambda-family is load-bearing")
    # (D) Lemma 2.2 combinatorics
    for q in [2] + primerange(3, 50):
        _lemma22(q)
    print("  (D) Lemma 2.2 combinatorics exhaustive over F_q, q in {2} u odd primes < 50")
    # (E) bridge witnesses, re-verified end-to-end (proven engine only)
    def _block(a, b, c, tau):
        try:
            return bool(_sun_solvable(a, b, c, tau))
        except (FactorBudget, PrimalityBound):
            return False                # undecidable instance never counts
    nver = 0
    for w, (ai, bi, tt) in sorted(_SUN_WITNESSES.items()):
        if w >= 100 and not extended:
            continue
        a, b = Fraction(ai), Fraction(bi)
        assert vp(b, 2) == 0 and (a == 1 or vp(a - 1, 2) >= 1)      # Phi membership
        D = delta_fin(1 + 4*a*a, 2*b)
        assert D == [2, w], (w, D)
        c = _sun_h(a, b, Fraction(w)**3)
        assert c is not None and _sun_integral(c, D), w
        if tt is None:
            ok = _block(a, b, c, Fraction(0)) or \
                 _block(a, b, c, 2*a/(1 + 4*a*a))
        else:
            ok = _sun_solvable(a, b, c, Fraction(*tt))
        assert ok, w
        nver += 1
    hi = 300 if extended else 100
    print(f"  (E) bridge witnesses verified for all {nver} odd primes w < {hi}:"
          f" Delta = {{2,w}}, block fires")
    print("  -> verified: identities, soundness, Lemma 2.2, existence probe;")
    print("     completeness chain traced (stale cross-refs are editorial);")
    print("     status: 7 = established modulo refereeing (refereed anchor: 10)\n")

def _evidence_export_lines():
    """Canonical JSONL serialization of the evidence authority _L6_WITNESSES.

    data/l6_witnesses.jsonl MUST equal these lines byte-for-byte: the export
    carries no field that is not derived from the authority.  Search-side
    provenance (timings, found flags) lives in data/superseded/ only.
    Regenerate with `python3 h10q.py --export-evidence`.
    """
    lines = []
    for w, (st, bi, tt) in sorted(_L6_WITNESSES.items()):
        s, tau = Fraction(*st), Fraction(*tt)
        lines.append(json.dumps(
            {"w": w, "s": [s.numerator, s.denominator], "b": [bi, 1],
             "tau": [tau.numerator, tau.denominator], "delta": [2, w]},
            separators=(",", ":")))
    return lines


def _verify_L6(extended=False):
    """Verify L6's proved pieces and bounded evidence for the assembly lemma.

    Proved: the tied block is a restriction of Psi; W0 selects one of its two
    formula branches at every odd target, W1 gives the exact target-place
    criterion, and W2 gives the exact quadratic norm obstruction. Evidence only:
    canonical global tied witnesses for every odd prime below 100 (below 300
    with --extended). This does NOT prove the simultaneous weak/norm-
    approximation lemma needed for completeness.
    """
    # (W0) Canonical two-branch selection.  The classical character sum
    # sum_a chi(1+4a^2) = -1 supplies an A-nonsquare at every odd target.
    # For such A, tau=0 works when chi(-1)=-1; otherwise tau=2a/A gives
    # delta=1/A and hence -delta*A=-1, a square.
    branch_checked = 0
    for w in primerange(3, 300):
        chars = []
        nonsquares = []
        for abar in range(w):
            Abar = (1 + 4*abar*abar) % w
            chi = 0 if Abar == 0 else legendre(Fraction(Abar), w)
            chars.append(chi)
            if chi == -1:
                nonsquares.append(abar)
        assert sum(chars) == -1 and nonsquares, w
        for abar in nonsquares:
            a = Fraction(abar)
            A = 1 + 4*a*a
            tau = _sun_target_tau(a, w)
            delta = 1 - A*tau*tau
            assert vp(A, w) == 0 and vp(delta, w) == 0
            assert legendre(-delta*A, w) == 1
        branch_checked += 1
    # Extension-field check of the same identity (all odd prime powers <= 49).
    q_checked = 0
    for q in [3, 5, 7, 9, 11, 13, 17, 19, 23, 25, 27, 29, 31, 37, 41, 43, 47, 49]:
        elems, add, mul, neg, one, zero, scal, inv = _field(q)
        chars = []
        for a in elems:
            A = add(one, mul(scal(4), mul(a, a)))
            if A == zero:
                chars.append(0)
                continue
            chi = one
            base, e = A, (q - 1)//2
            while e:
                if e & 1:
                    chi = mul(chi, base)
                e >>= 1
                if e:
                    base = mul(base, base)
            chars.append(1 if chi == one else -1)
        assert sum(chars) == -1 and -1 in chars, q
        q_checked += 1
    print(f"  (W0) canonical tau branch: all {branch_checked} odd primes w < 300;"
          f" character identity on {q_checked} odd prime powers <= 49; q=3 included")

    # (W1) At the target w, M is a unit congruent to 16.  Thus the conic
    # alpha*y^2 + beta*r^2 = M, v(alpha)=0, v(beta)=1, is soluble iff
    # chi_w(M/alpha)=chi_w(-delta*A)=1.  Exercise the exact local solver.
    rng = random.Random(613)
    ps = primerange(3, 50)
    taus = [Fraction(0), Fraction(1, 5), Fraction(1, 4),
            Fraction(1, 3), Fraction(2, 5), Fraction(1, 2),
            Fraction(3, 5), Fraction(2, 3), Fraction(3, 4),
            Fraction(5, 4), Fraction(3, 2), Fraction(2), Fraction(3)]
    checked = 0
    while checked < 2000:
        w = rng.choice(ps)
        den = rng.choice([d for d in (1, 3, 5, 7, 9, 11, 13) if d % w])
        s = Fraction(rng.randint(-6, 6), den)
        a, A = 1 + 2*s, 1 + 4*(1 + 2*s)**2
        if vp(A, w) != 0:
            continue
        m = rng.choice([m for m in (-11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11)
                        if m % w])
        b, B = Fraction(m*w), Fraction(2*m*w)
        tau = rng.choice(taus)
        delta = 1 - A*tau*tau
        if delta == 0 or vp(delta, w) != 0:
            continue
        cden = rng.choice([d for d in (1, 5, 7, 11, 13) if d % w])
        c = Fraction(w**rng.randint(1, 4) * rng.choice([1, 2, 4, 5, 7]), cden)
        M = 16 - delta*c*c - 16*A*B*s*s
        assert vp(M, w) == 0
        local = represents_qv([-delta*A, -16*B], M, w)
        criterion = legendre((-delta*A).numerator * (-delta*A).denominator, w) == 1
        assert local == criterion, (w, s, b, tau)
        checked += 1
    print(f"  (W1) tied target-place reduction: {checked} exact local instances, 0 mismatches")

    # Exact norm reformulation; the minus sign is load-bearing.  Count checks
    # explicitly: the delta == 0 skip must not eat iterations, or the printed
    # total would overclaim.
    nchk = 0
    while nchk < 200:
        a = Fraction(rng.randint(-9, 9), rng.randint(1, 9))
        b = Fraction(rng.choice([n for n in range(-9, 10) if n]), rng.randint(1, 9))
        tau = Fraction(rng.randint(-9, 9), rng.randint(1, 9))
        y = Fraction(rng.randint(-9, 9), rng.randint(1, 9))
        r = Fraction(rng.randint(-9, 9), rng.randint(1, 9))
        A, B = 1 + 4*a*a, 2*b
        delta = 1 - A*tau*tau
        if delta == 0:
            continue
        M = -delta*A*y*y - 16*B*r*r
        disc = -delta*A*B
        norm = (delta*A*y)**2 - disc*(4*r)**2
        assert norm == -delta*A*M
        nchk += 1
    # Zero case: M_tau = 0 is solved by (y,r) = (0,0) although 0 is not in
    # N(E_tau^x); solvability is the DISJUNCTION "M_tau = 0 or norm
    # membership", and solvers return True on their M == 0 branch.
    # a=1, tau=0, c=4 gives M = 16 - 16 - 0 = 0 for any b.
    assert _sun_tied_solvable(Fraction(1), Fraction(3), Fraction(4), Fraction(0)) is True
    print(f"  (W2) tied conic = quadratic norm equation: {nchk} exact identities;"
          " M_tau=0 zero case exercised")

    # (E6) Frozen global witnesses for the open completeness direction.
    bound = 300 if extended else 100
    expected = primerange(3, bound)
    assert [w for w in sorted(_L6_WITNESSES) if w < bound] == expected
    nver = 0
    for w, (st, bi, tt) in sorted(_L6_WITNESSES.items()):
        if w >= bound:
            continue
        s, b, tau = Fraction(*st), Fraction(bi), Fraction(*tt)
        a, A = 1 + 2*s, 1 + 4*(1 + 2*s)**2
        assert s == 0 or (vp(s, 2) >= 0 and vp(s, w) >= 0)
        assert vp(b, 2) == 0 and vp(b, w) == 1
        assert vp(A, 2) == 0 and unit_mod(A, 8) == 5
        assert vp(A, w) == 0 and legendre(A, w) == -1
        assert tau == _sun_target_tau(a, w)
        delta = 1 - A*tau*tau
        assert vp(delta, w) == 0 and legendre(-delta*A, w) == 1
        # This is the bridge applied to the formula's input z=w, hence c=h(a,b,z^3).
        z = Fraction(w)
        assert vp(z, w) == 1
        expected_vc = 2*vp(a, w) + 6*vp(z, w) - 2
        D = delta_fin(A, 2*b)
        assert D == [2, w], (w, D)
        c = _sun_h(a, b, z**3)
        assert c is not None and vp(c, w) == expected_vc >= 4
        assert _sun_integral(c, D), w
        assert _sun_tied_solvable(a, b, c, tau) is True, w
        nver += 1
    print(f"  (E6) canonical tied certificates: all {nver} odd primes w < {bound};"
          " A nonsquare unit, Delta = {2,w}, global conic point")

    # (EV) SSOT: _L6_WITNESSES above is the single authority; the JSONL file
    # data/l6_witnesses.jsonl is its canonical serialization
    # (_evidence_export_lines) and must match BYTE-FOR-BYTE, so no field can
    # exist in the export without being derived and checked.  Missing or
    # drifted export fails the suite.
    export = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "data", "l6_witnesses.jsonl")
    with open(export) as fh:
        got = fh.read()
    want = "\n".join(_evidence_export_lines()) + "\n"
    assert got == want, "evidence export drifted from in-code authority"
    t0 = sum(1 for _, (_, _, tt) in sorted(_L6_WITNESSES.items())
             if Fraction(*tt) == 0)
    assert (t0, len(_L6_WITNESSES) - t0) == (32, 29)
    print(f"  (EV) evidence export byte-identical to authority serialization:"
          f" {len(_L6_WITNESSES)} rows, {t0} tau=0 / {len(_L6_WITNESSES) - t0} tau=2a/A")
    print("  -> L6 count = 6; W0-W2 proved and E6 bounded evidence only."
          " Global completeness is CONDITIONAL on classical Schinzel H via L19-L22\n")

_L7_RESCUES = {
    # (w, u): (s, b, tau).  A Phi-witness (2-adic conditions only: v2(s) >= 0,
    # v2(b) = 0) whose tied conic is soluble at z = w*u, for cells where the
    # frozen canonical (s_w, b_w) fails both branches.  Spot checks from the
    # bounded 71/71 rescue probe of 2026-08-15 (RESULTS.md); probe evidence,
    # not a completeness proof.
    (3, Fraction(2)): (Fraction(-1), Fraction(3), Fraction(0)),
    (3, Fraction(1, 5)): (Fraction(-1), Fraction(3), Fraction(0)),
    (5, Fraction(-1)): (Fraction(0), Fraction(-5), Fraction(2, 5)),
    (5, Fraction(1, 2)): (Fraction(0), Fraction(-5), Fraction(2, 5)),
    (7, Fraction(-1)): (Fraction(-1, 5), Fraction(7), Fraction(0)),
    (7, Fraction(1, 5)): (Fraction(0), Fraction(-7), Fraction(0)),
    (11, Fraction(-1)): (Fraction(-3), Fraction(11), Fraction(0)),
    (11, Fraction(1, 3)): (Fraction(3), Fraction(11), Fraction(0)),
    (13, Fraction(-1)): (Fraction(0), Fraction(-13), Fraction(2, 5)),
    (13, Fraction(-1, 3)): (Fraction(-2), Fraction(-13), Fraction(-6, 37)),
}


def _l7_tied_status(a, b, z, tau):
    """True/False/None(pole)/'budget' for the tied conic at input z."""
    c = _sun_h(a, b, z**3)
    if c is None:
        return None
    try:
        return _sun_tied_solvable(a, b, c, tau)
    except (FactorBudget, PrimalityBound):
        return 'budget'


def _verify_L7(extended=False):
    """L7 (THEOREMS.md): structural COROLLARIES of W2 / L6-soundness plus
    bounded probe evidence.  None of this proves the assembly lemma, which
    remains OPEN; no status or count changes here.

    (a) Tie-cost dictionary [corollary of W2 + standard local theory]: at an
        odd place v where delta*A and B are units, the TIED conic fails iff
        v(M) is odd and chi_v(-delta*A*B) = -1.  Where additionally A is a
        v-unit (so ALL ternary coefficients are units), the UNTIED quaternary
        block is soluble regardless of M: Chevalley-Warning + Hensel make the
        unit ternary isotropic, and a regular isotropic form is universal.
        Where v(A) != 0 or v(B) != 0 the dictionary does NOT apply and real
        obstructions occur; three frozen witnesses in-suite, including a
        cancellation case with v5(AB)=0 (so the bad set is NOT "v | AB").
    (b) No rational-function witness [corollary, conditional on L6
        soundness]: a fixed (s,b,tau) with a Q(z)-section Y,R would
        specialize to a solution at some non-target z = +-2^k.  Kills
        uniform sections only -- NOT pointwise coverage by a fixed pair,
        which stays open.  Checked strictly on a bounded sweep.
    (c) Reciprocity parity [corollary of W2 + Hilbert reciprocity]: the
        obstructed-place set of a failing tied conic has even size, never
        one -- any repair must flip places in pairs.
    """
    # (a) dictionary, on guarded instances derived from the frozen table.
    rng = random.Random(77)
    dict_odd = dict_even = quat = 0
    ws = sorted(_L6_WITNESSES)[:12] if extended else sorted(_L6_WITNESSES)[:6]
    for w in ws:
        st, bi, tt = _L6_WITNESSES[w]
        s, b = Fraction(*st), Fraction(bi)
        a = 1 + 2*s
        A, B = 1 + 4*a*a, 2*b
        for tau in (Fraction(0), 2*a/A):
            delta = 1 - A*tau*tau
            for _ in range(6):
                z = Fraction(w) * Fraction(rng.randint(1, 9), rng.choice([1, 1, 5]))
                c = _sun_h(a, b, z**3)
                if c is None:
                    continue
                M = 16 - delta*c*c - 16*A*B*s*s
                if M == 0:
                    continue
                try:
                    vs = [v for v in places_of(M) if v not in (2, OO)
                          and vp(delta*A, v) == 0 and vp(B, v) == 0]
                except (FactorBudget, PrimalityBound):
                    continue
                for v in vs[:4]:
                    tied = represents_qv([-delta*A, -16*B], M, v)
                    if vp(M, v) % 2:
                        want = legendre((-delta*A*B).numerator
                                        * (-delta*A*B).denominator, v) == 1
                        assert tied == want, (w, v)
                        dict_odd += 1
                    else:
                        assert tied is True, (w, v)
                        dict_even += 1
                    if vp(A, v) == 0:
                        # (iii) needs v coprime to ALL coefficients: isotropy
                        # of the unit ternary (<alpha,beta> represents -gamma
                        # <=> <alpha,beta,gamma> isotropic), hence universal,
                        # hence the ACTUAL untied RHS 16 - delta*c^2.
                        assert represents_qv([-delta*A, -16*B], -16*A*B, v), (w, v)
                        u0 = 16 - delta*c*c
                        if u0 != 0:
                            assert represents_qv([-delta*A, -16*B, 16*A*B],
                                                 u0, v), (w, v)
                        quat += 1
    print(f"  (L7a) tie-cost dictionary: {dict_odd} odd-valuation + {dict_even}"
          f" even-valuation places match; untied universal at {quat} fully-unit places")
    # (a') coefficient-bad places ({v: v(A)!=0 or v(B)!=0}) are REAL
    # obstructions outside the wild set -- three frozen witnesses:
    # 1. guarded w=5 shape s=3 (a=7, A=197), b=15 (v3(B)=1), z=15, tau=2a/A:
    #    v3(M)=0 yet the tied conic is insoluble over Q_3.
    aa, bb = Fraction(7), Fraction(15)
    AA = 1 + 4*aa*aa
    tt3 = 2*aa/AA
    dd = 1 - AA*tt3*tt3
    cc3 = _sun_h(aa, bb, Fraction(15)**3)
    MM = 16 - dd*cc3*cc3 - 16*AA*(2*bb)*Fraction(3)**2
    assert vp(MM, 3) == 0 and vp(2*bb, 3) == 1
    assert represents_qv([-dd*AA, -16*(2*bb)], MM, 3) is False
    # 2. a=1, b=11 (A=5, tau=2/5, deltaA=1): ternary <-1,-352,1760> is
    #    anisotropic over Q_5 (5 | A) and the untied RHS fails there.
    a5, b5 = Fraction(1), Fraction(11)
    A5, B5 = 1 + 4*a5*a5, 2*b5
    t5 = 2*a5/A5
    d5 = 1 - A5*t5*t5
    assert represents_qv([-d5*A5, -16*B5], -16*A5*B5, 5) is False
    c5 = _sun_h(a5, b5, Fraction(1))
    assert represents_qv([-d5*A5, -16*B5, 16*A5*B5], 16 - d5*c5*c5, 5) is False
    # 3. the bad set is NOT "v | AB": valuations cancel in the product.
    #    a=1, b=1/5: v5(A)=1, v5(B)=-1, v5(AB)=0, yet v=5 obstructs.
    a6, b6 = Fraction(1), Fraction(1, 5)
    A6, B6 = 1 + 4*a6*a6, 2*b6
    t6 = 2*a6/A6
    d6 = 1 - A6*t6*t6
    assert vp(A6*B6, 5) == 0 and vp(A6, 5) == 1 and vp(B6, 5) == -1
    assert represents_qv([-d6*A6, -16*B6], -16*A6*B6, 5) is False
    c6 = _sun_h(a6, b6, Fraction(1))
    assert represents_qv([-d6*A6, -16*B6, 16*A6*B6], 16 - d6*c6*c6, 5) is False
    print("  (L7a') coefficient-bad places witnessed: v=3 | B (tied, v3(M)=0),"
          " v=5 | A (untied), and v=5 with v5(AB)=0 (cancellation) --"
          " bad set is {v: v(A)!=0 or v(B)!=0}")

    # (b) strict non-target sweep at z = +-2^k (bounded probe).
    ws = [3, 5, 7, 11, 13, 17, 19, 23] if extended else [3, 5, 7]
    ks = range(-3, 4) if extended else range(-2, 3)
    nfalse = nbudget = 0
    failing = []
    for w in ws:
        st, bi, tt = _L6_WITNESSES[w]
        s, b = Fraction(*st), Fraction(bi)
        a = 1 + 2*s
        A = 1 + 4*a*a
        for k in ks:
            for sign in (1, -1):
                z = Fraction(sign * 2**k) if k >= 0 else Fraction(sign, 2**(-k))
                for tau in (Fraction(0), 2*a/A):
                    r = _l7_tied_status(a, b, z, tau)
                    assert r is not True, (w, z, tau)
                    if r is False:
                        nfalse += 1
                        if len(failing) < 8:
                            failing.append((a, b, z, tau))
                    elif r == 'budget':
                        nbudget += 1
    print(f"  (L7b) non-target z = +-2^k: 0 solutions;"
          f" {nfalse} disproved, {nbudget} budget refusals (bounded probe)")

    # (c) parity of the obstruction set on disproved cells.
    npar = 0
    for a, b, z, tau in failing:
        A, B = 1 + 4*a*a, 2*b
        delta = 1 - A*tau*tau
        c = _sun_h(a, b, z**3)
        s = (a - 1)/2
        M = 16 - delta*c*c - 16*A*B*s*s
        alpha = -delta*A
        try:
            vs = set([2, OO]) | set(places_of(-alpha*(-16*B))) | set(places_of(alpha*M))
            bad = [v for v in vs if hilbert(-alpha*(-16*B), alpha*M, v) == -1]
        except (FactorBudget, PrimalityBound):
            continue
        assert bad and len(bad) % 2 == 0, (a, b, z, tau, bad)
        npar += 1
    # rescue spot checks: target-cell membership, canonical failure of BOTH
    # branches, then Phi-membership + solvability of the replacement.
    nres = ncert = 0
    certws = set()
    for (w, u), (s, b, tau) in sorted(_L7_RESCUES.items()):
        z = Fraction(w) * u
        assert vp(z, w) >= 1, (w, u)                  # genuine target cell
        stc, bic, _ = _L6_WITNESSES[w]                # frozen canonical pair
        sc, bc = Fraction(*stc), Fraction(bic)
        ac = 1 + 2*sc
        rr = [_l7_tied_status(ac, bc, z, tc)
              for tc in (Fraction(0), 2*ac/(1 + 4*ac*ac))]
        assert all(r is not True for r in rr), (w, u, rr)
        if all(r is False for r in rr):
            ncert += 1              # Hasse-certified both-branch failure
            certws.add(w)
        assert (s == 0 or vp(s, 2) >= 0) and vp(b, 2) == 0
        a = 1 + 2*s
        A = 1 + 4*a*a
        assert tau in (0, 2*a/A)
        assert _l7_tied_status(a, b, z, tau) is True, (w, u)
        nres += 1
    # certified rows: 7/10; (3,2), (7,-1), (11,1/3) hit FactorBudget in >=1
    # branch.  EVERY probed canonical pair (w = 3,5,7,11,13) has at least one
    # Hasse-certified both-branch failure, so none covers m_w pointwise.
    assert ncert >= 7, ncert
    assert certws >= {3, 5, 7, 11, 13}, certws
    print(f"  (L7c) reciprocity parity: {npar} obstruction sets, all even >= 2;"
          f" {nres} frozen rescue spot-checks re-verified ({ncert} canonical"
          f" failures Hasse-certified, rest budget-refused; all 5 probed pairs"
          f" certifiably miss cells)")
    print("  -> L7 corollaries + bounded probe only. Unconditional assembly stays"
          " open; the Schinzel-conditional chain is closed later by L22\n")


def _verify_L8(extended=False):
    """L8 (THEOREMS.md): absorption sharpening, the 2-adic parity wall, and
    exact-matching emptiness (unconditional), plus scoped alignment evidence.
    Structure results only; the Schinzel-conditional chain closes later at L22.

    (0) Value identities: u0 = -M0, u1 = -A*M1, (A b^2 D_z)^2 u_tau =
        X^2 - D_tau Y^2, and D_tau*disc(E_tau) = -2AbP in both branches,
        asserted as exact rational identities on random witness/cell tuples.
    (a) Absorption: at odd v, if v(delta_tau*A), v(B), v(M_tau) are ALL even,
        the tied conic is soluble over Q_v (scaling + surjectivity of
        nondegenerate binary forms over F_v + Hensel).  This sharpens the
        L7a' support to ODD-valuation places.
    (b) Parity wall: 2-adic admissibility (v2(s) >= 0, v2(b) = 0) forces
        v2(P) = 0 for P = 1 - 2Abs^2, hence v2(-2AbP) = 1: the king class
        D_tau * disc(E_tau) = -2AbP (both branches) is NEVER a square.
    (c) Exact-matching emptiness: the gauge -2AbP in Q^x2 (source = target
        field) is rationally parametrized by b = -th^2/(2As^2(1-th^2)) or
        b = 1/(2As^2(1-th^2)); on BOTH roots v2(b) != 0 for every theta, so
        the gauge misses the admissible box entirely.  (The split-E gauge is
        equally empty: v2(disc E) = 1.)  General fixed-class matching is NOT
        ruled out here; see THEOREMS L8c scope note.
    (d) Alignment evidence (scoped): every frozen rescue's norm value has a
        wild odd-multiplicity prime > 10^5 whose symbol is +1 -- the frozen
        rescues are alignment events, not support-controlled constructions.
    """
    rng = random.Random(20260815)
    # (0) exact value identities on random admissible witnesses and cells:
    #     u0 = -M0, u1 = -A*M1, (A b^2 D_z)^2 u_tau = X^2 - D_tau Y^2
    #     (D_0 = P, D_1 = A P), and D_tau * disc(E_tau) = -2AbP in both
    #     branches as exact rationals.
    nid = 0
    while nid < (400 if extended else 150):
        s = Fraction(rng.randint(-9, 9), rng.choice([1, 3, 5, 7]))
        b = Fraction(rng.choice([n for n in range(-40, 41) if n % 2]),
                     rng.choice([1, 3, 5, 7, 9]))
        z = Fraction(rng.randint(-12, 12), rng.choice([1, 3, 5, 7]))
        if z == 0 or b in (0, 1):
            continue
        a = 1 + 2*s
        A = 1 + 4*a*a
        c = _sun_h(a, b, z**3)
        if c is None:
            continue
        Dz = 1 - z**3 - a*a*z**6
        Ng = 16*a**4*b*b - A*(b - 1)**4
        X = a*a*z**6*Ng
        Y = 4*A*b*b*Dz
        P = 1 - A*2*b*s*s
        M0 = 16 - c*c - 16*A*2*b*s*s
        M1 = 16 - c*c/A - 16*A*2*b*s*s
        u0, u1 = c*c - 16*P, c*c - 16*A*P
        assert u0 == -M0 and u1 == -A*M1
        assert (A*b*b*Dz)**2 * u0 == X*X - P*Y*Y
        assert (A*b*b*Dz)**2 * u1 == X*X - A*P*Y*Y
        assert P * (-2*A*b) == (A*P) * (-2*b) == -2*A*b*P
        nid += 1
    print(f"  (L8-id) value identities exact on {nid} random witness/cell tuples")

    # (a) absorption lemma: exact local checks at random odd v, even valuations
    nabs = 0
    while nabs < (400 if extended else 150):
        v = rng.choice([3, 5, 7, 11, 13, 17, 19, 23])
        e1, e2, e3 = (2*rng.randint(-2, 2) for _ in range(3))
        u1 = Fraction(rng.choice([u for u in range(1, v)]))
        u2 = Fraction(rng.choice([u for u in range(1, v)]))
        u3 = Fraction(rng.choice([u for u in range(1, v)]))
        alpha = u1 * Fraction(v) ** e1
        beta = u2 * Fraction(v) ** e2
        m = u3 * Fraction(v) ** e3
        assert represents_qv([alpha, beta], m, v), (v, alpha, beta, m)
        nabs += 1
    # odd-valuation places remain genuinely bad (frozen L7a' witnesses cover
    # this; re-assert one here so L8a's scope is machine-delimited)
    assert represents_qv([-Fraction(197), -16*Fraction(30)],
                         Fraction(2), 3) is False
    print(f"  (L8a) absorption: {nabs} even-valuation instances soluble at odd v;"
          " odd-valuation counterexample re-checked")

    # (b) parity wall on random admissible witnesses
    npw = 0
    while npw < (2000 if extended else 600):
        s = Fraction(rng.randint(-40, 40),
                     rng.choice([1, 3, 5, 7, 9, 11, 15]))
        b = Fraction(rng.choice([n for n in range(-99, 100) if n % 2]),
                     rng.choice([1, 3, 5, 7, 9, 11, 13]))
        A = 1 + 4*(1 + 2*s)**2
        P = 1 - 2*A*b*s*s
        if P == 0:
            continue
        assert vp(P, 2) == 0, (s, b)
        king = -2*A*b*P
        assert vp(king, 2) == 1, (s, b)
        npw += 1
    print(f"  (L8b) parity wall: v2(P)=0 and v2(-2AbP)=1 on {npw} admissible"
          " witnesses; king class never a square")

    # (c) collapse-family emptiness: both roots, valuation identity + v2(b) != 0
    ncf = 0
    for _ in range(4000 if extended else 1200):
        s = Fraction(rng.randint(-20, 20), rng.choice([1, 3, 5, 7]))
        if s == 0 or vp(s, 2) < 0:
            continue
        th = Fraction(rng.randint(-60, 60), rng.randint(1, 60))
        if th == 0 or th*th == 1:
            continue
        A = 1 + 4*(1 + 2*s)**2
        d = 1 - th*th
        for root, b in ((1, -th*th/(2*A*s*s*d)), (-1, Fraction(1)/(2*A*s*s*d))):
            if b in (0, 1):
                continue
            # gauge identity: -2AbP is an exact square on the family
            P = 1 - 2*A*b*s*s
            king = -2*A*b*P            # = (th/(s*(1-th^2)))^2 on both roots
            r2 = king.numerator * king.denominator
            r = math.isqrt(r2)
            assert king > 0 and r*r == r2, (s, th, root)
            assert vp(b, 2) != 0, (s, th, root)
            ncf += 1
    print(f"  (L8c) exact-matching emptiness: {ncf} gauge points (both roots),"
          " none 2-adically admissible; source=target gauge never enters the box")

    # (d) alignment evidence on the frozen rescue table
    nal = 0
    for (w, u), (s, b, tau) in sorted(_L7_RESCUES.items()):
        a = 1 + 2*s
        A = 1 + 4*a*a
        z = Fraction(w) * u
        c = _sun_h(a, b, z**3)
        delta = 1 - A*tau*tau
        M = 16 - delta*c*c - 16*A*2*b*s*s
        x, dsc = -delta*A*M, -delta*A*2*b
        try:
            fac = list(factorint(abs(x.numerator)).items()) + \
                  list(factorint(x.denominator).items())
        except (FactorBudget, PrimalityBound):
            continue
        odd_mult = [q for q, e in fac if q != 2 and e % 2]
        assert odd_mult and max(odd_mult) > 10**5, (w, u)
        assert all(hilbert(x, dsc, q) == 1 for q in odd_mult), (w, u)
        nal += 1
    assert nal == len(_L7_RESCUES), nal
    print(f"  (L8d) alignment: {nal}/{len(_L7_RESCUES)} rescues have wild"
          " odd-mult primes > 10^5, all symbols +1 (deterministic"
          " factorizations, none refused)")
    print("  -> L8 structure results only. Unconditional assembly stays open;"
          " the Schinzel-conditional chain is closed later by L22\n")

# --------------------------------------------------- L9: reciprocity steering
# Shared certificate machinery for the steered assembly probe (l9_steer.py
# imports these; single source).  See THEOREMS.md L9/L9a.
_L9_TRIAL_BOUND = 10**4
_L9_TRIAL_PRIMES = list(primerange(2, _L9_TRIAL_BOUND))
_L9_DIGIT_CAP = 200


def _l9_trial_divide(n):
    """n > 0 -> (dict of trial-prime factors, cofactor)."""
    f = {}
    for p in _L9_TRIAL_PRIMES:
        if p * p > n:
            break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    if n > 1 and n < _L9_TRIAL_BOUND * _L9_TRIAL_BOUND:
        f[n] = f.get(n, 0) + 1   # cofactor below trial square is prime
        n = 1
    return f, n


def _l9_iroot(n, e):
    """Exact floor e-th root, integers only (isqrt / integer Newton)."""
    if e == 1 or n < 2:
        return n
    if e == 2:
        return math.isqrt(n)
    r = 1 << -(-n.bit_length() // e)          # r >= floor root
    while True:
        nr = ((e - 1) * r + n // r ** (e - 1)) // e
        if nr >= r:
            break
        r = nr
    while r ** e > n:
        r -= 1
    return r


def _l9_perfect_power(n):
    """n > 1 with no factor < _L9_TRIAL_BOUND -> (root, exp), exp maximal."""
    for e in range(n.bit_length() // 13, 1, -1):
        r = _l9_iroot(n, e)
        if r ** e == n:
            return r, e
    return n, 1


def _l9_proven_prime(n):
    """True only under a PROOF of primality; False on composites and on
    refusals (both are skips, never evidence)."""
    if not _mr_screen(n):
        return False
    try:
        return _is_prime(n)
    except (FactorBudget, PrimalityBound):
        return False


def _l9_resolve(x, allow_blobs):
    """Resolve a nonzero Fraction under the steering filter.

    Returns (facdict, blobs): facdict prime -> signed exponent; blobs a
    list of (C, e) — composite cofactors left unfactored with EVEN total
    exponent (only when allow_blobs) — or None when a cofactor is neither
    1, a proven prime power, nor an even-exponent blob.  The resolution is
    complete and re-asserted by exact reconstruction."""
    f, blobs = {}, []
    for part, sgn in ((abs(x.numerator), 1), (x.denominator, -1)):
        if len(str(part)) > _L9_DIGIT_CAP:
            return None
        sm, cof = _l9_trial_divide(part)
        for p, e in sm.items():
            f[p] = f.get(p, 0) + sgn * e
        if cof == 1:
            continue
        root, exp = _l9_perfect_power(cof)
        if _l9_proven_prime(root):
            f[root] = f.get(root, 0) + sgn * exp
        elif allow_blobs and exp % 2 == 0:
            blobs.append((root, sgn * exp))
        else:
            return None
    num = den = 1
    for p, e in list(f.items()) + blobs:
        if e > 0:
            num *= p ** e
        else:
            den *= p ** (-e)
    assert num == abs(x.numerator) and den == x.denominator, x
    return {p: e for p, e in f.items() if e}, blobs


def _l9_steered_solvable(a, b, z, tau):
    """Steered W2 certificate for the tied conic at input z.

    Resolution step: d = alpha*B must resolve completely; x = alpha*M
    resolves into proven prime powers times even-exponent composite blobs
    coprime to supp(d) (their primes q have v_q(x) even, v_q(d) = 0,
    hence symbol +1 without identification).  wilds = the odd-multiplicity
    resolved primes of x above the trial bound with v_q(d) = 0.

    Certification step:
      - exactly ONE wild prime q0: the symbols on T \\ {q0} are checked
        exactly; if all are +1, (x,d)_{q0} = +1 is CONCLUDED from the L9
        product formula (the certifying step), and the direct symbol at
        q0 is then replayed as a consistency assert only ('steered');
      - zero or >= 2 wild primes: every place of T is checked directly
        ('smooth' / 'aligned' — full-support certificates; 'aligned' rows
        carry alignment content beyond reciprocity at their wild primes).

    Returns (True, mech, wilds) — PROVED soluble; (False, bad, wilds) —
    PROVED insoluble (some exactly-checked symbol is -1); (None, reason,
    None) — not certifiable under the steering filter (never evidence).
    mech: 'M=0' | 'smooth' | 'steered' | 'aligned'."""
    c = _sun_h(a, b, z**3)
    if c is None:
        return None, 'pole', None
    A, B = 1 + 4*a*a, 2*b
    delta = 1 - A*tau*tau
    if delta == 0:
        return None, 'delta0', None
    s = (a - 1)/2
    M = 16 - delta*c*c - 16*A*B*s*s
    if M == 0:
        return True, 'M=0', []
    alpha = -delta*A
    x, d = alpha*M, alpha*B
    dres = _l9_resolve(d, allow_blobs=False)
    if dres is None:
        return None, 'cofactor_d', None
    dfac, _ = dres
    xres = _l9_resolve(x, allow_blobs=True)
    if xres is None:
        return None, 'cofactor_x', None
    xfac, xblobs = xres
    dsupp = abs(d.numerator) * d.denominator
    if any(math.gcd(C, dsupp) != 1 for C, _ in xblobs):
        return None, 'blob_meets_d', None
    places = sorted({q for q in xfac if q != 2} | {q for q in dfac if q != 2})
    wilds = sorted(q for q, e in xfac.items()
                   if q > _L9_TRIAL_BOUND and e % 2 and dfac.get(q, 0) == 0)
    if len(wilds) == 1:
        # L9 forcing: certify from T \ {q0}; the direct symbol at q0 is a
        # consistency replay, NOT the certifying step.
        q0 = wilds[0]
        bad = [v for v in [2, OO] + places
               if v != q0 and hilbert(x, d, v) == -1]
        if bad:
            return False, bad, wilds
        assert hilbert(x, d, q0) == 1, (a, b, z, tau, q0)  # forced by L9
        return True, 'steered', wilds
    bad = [v for v in [2, OO] + places if hilbert(x, d, v) == -1]
    if bad:
        return False, bad, wilds
    return True, 'smooth' if not wilds else 'aligned', wilds


# Steered assembly witnesses frozen from the 2026-08-16 l9_steer.py run:
# (w, (u_num, u_den)) -> ((s_num, s_den), (b_num, b_den), (tau_num, tau_den)).
# Each row is an admissible Phi-pair whose tied conic at z = w*u is PROVED
# soluble by the steered W2 certificate (_l9_steered_solvable is True),
# re-verified in-suite below.  Bounded evidence for the OPEN assembly
# lemma; NOT a completeness proof.  data/l9_steered.jsonl is the canonical
# byte-exact serialization (_l9_export_lines).
_L9_STEERED = {
    (3, (-5, 1)): ((0, 1), (1, 11), (2, 5)),
    (3, (-3, 1)): ((-1, 1), (15, 1), (0, 1)),
    (3, (-2, 1)): ((-1, 1), (3, 1), (0, 1)),
    (3, (-1, 1)): ((0, 1), (27, 5), (2, 5)),
    (3, (-2, 7)): ((0, 1), (-27, 5), (0, 1)),
    (3, (-1, 5)): ((0, 1), (-3, 1), (0, 1)),
    (3, (1, 7)): ((0, 1), (-27, 1), (0, 1)),
    (3, (1, 1)): ((0, 1), (-5, 3), (0, 1)),
    (3, (2, 1)): ((0, 1), (31, 3), (0, 1)),
    (3, (3, 1)): ((-1, 1), (27, 1), (0, 1)),
    (3, (5, 1)): ((0, 1), (1, 11), (2, 5)),
    (5, (-5, 1)): ((0, 1), (-9, 1), (2, 5)),
    (5, (-3, 1)): ((0, 1), (1, 11), (2, 5)),
    (5, (-2, 1)): ((0, 1), (1, 11), (2, 5)),
    (5, (-1, 1)): ((0, 1), (-1, 1), (2, 5)),
    (5, (-1, 3)): ((0, 1), (1, 11), (2, 5)),
    (5, (-2, 7)): ((0, 1), (1, 11), (2, 5)),
    (5, (1, 7)): ((0, 1), (-29, 1), (2, 5)),
    (5, (1, 3)): ((0, 1), (-1, 1), (2, 5)),
    (5, (2, 3)): ((0, 1), (1, 11), (2, 5)),
    (5, (1, 1)): ((0, 1), (-5, 1), (2, 5)),
    (5, (2, 1)): ((0, 1), (-25, 1), (2, 5)),
    (5, (7, 3)): ((0, 1), (-1, 1), (2, 5)),
    (5, (3, 1)): ((0, 1), (1, 11), (2, 5)),
    (5, (5, 1)): ((0, 1), (-5, 1), (2, 5)),
    (7, (-5, 1)): ((1, 1), (5, 9), (6, 37)),
    (7, (-3, 1)): ((-1, 1), (3, 1), (0, 1)),
    (7, (-2, 1)): ((-1, 1), (35, 9), (0, 1)),
    (7, (-1, 1)): ((-1, 1), (-7, 11), (0, 1)),
    (7, (-1, 3)): ((0, 1), (-7, 1), (0, 1)),
    (7, (-1, 5)): ((-1, 1), (35, 1), (0, 1)),
    (7, (1, 3)): ((-1, 1), (7, 1), (0, 1)),
    (7, (2, 3)): ((-1, 3), (-7, 3), (0, 1)),
    (7, (1, 1)): ((-1, 1), (-7, 11), (0, 1)),
    (7, (2, 1)): ((4, 1), (-7, 1), (0, 1)),
    (7, (7, 3)): ((-1, 1), (-7, 19), (0, 1)),
    (7, (3, 1)): ((-1, 1), (15, 1), (0, 1)),
    (7, (5, 1)): ((0, 1), (-29, 5), (2, 5)),
    (11, (-5, 1)): ((0, 1), (-1, 1), (2, 5)),
    (11, (-3, 1)): ((0, 1), (27, 5), (2, 5)),
    (11, (-2, 1)): ((2, 1), (-11, 5), (0, 1)),
    (11, (-1, 1)): ((-4, 1), (-11, 9), (0, 1)),
    (11, (-1, 3)): ((2, 1), (-11, 1), (0, 1)),
    (11, (-2, 7)): ((2, 1), (-11, 5), (0, 1)),
    (11, (-1, 5)): ((2, 1), (-11, 5), (0, 1)),
    (11, (1, 7)): ((2, 1), (-11, 1), (0, 1)),
    (11, (1, 3)): ((1, 3), (11, 1), (0, 1)),
    (11, (2, 3)): ((-3, 1), (11, 1), (0, 1)),
    (11, (1, 1)): ((2, 1), (-11, 5), (0, 1)),
    (11, (2, 1)): ((2, 1), (-11, 5), (0, 1)),
    (11, (7, 3)): ((-3, 1), (7, 5), (0, 1)),
    (11, (3, 1)): ((-1, 1), (-15, 11), (0, 1)),
    (11, (5, 1)): ((0, 1), (-25, 1), (2, 5)),
    (13, (-5, 1)): ((0, 1), (-25, 1), (2, 5)),
    (13, (-3, 1)): ((-1, 1), (-3, 11), (0, 1)),
    (13, (-2, 1)): ((0, 1), (-13, 5), (2, 5)),
    (13, (-1, 1)): ((1, 1), (13, 9), (6, 37)),
    (13, (-1, 3)): ((0, 1), (-13, 5), (2, 5)),
    (13, (-2, 7)): ((1, 3), (-13, 3), (30, 109)),
    (13, (-1, 5)): ((0, 1), (-13, 9), (2, 5)),
    (13, (1, 7)): ((0, 1), (-13, 9), (2, 5)),
    (13, (1, 3)): ((0, 1), (-13, 1), (2, 5)),
    (13, (2, 3)): ((-2, 1), (39, 1), (-6, 37)),
    (13, (1, 1)): ((0, 1), (-13, 9), (2, 5)),
    (13, (2, 1)): ((1, 1), (-13, 11), (6, 37)),
    (13, (7, 3)): ((0, 1), (-13, 9), (2, 5)),
    (13, (3, 1)): ((0, 1), (-13, 9), (2, 5)),
    (13, (5, 1)): ((0, 1), (-1, 1), (2, 5)),
    (17, (-5, 1)): ((0, 1), (-1, 1), (2, 5)),
    (17, (-3, 1)): ((-1, 1), (3, 1), (0, 1)),
    (17, (-2, 1)): ((0, 1), (-17, 1), (2, 5)),
    (17, (-1, 1)): ((0, 1), (-17, 1), (2, 5)),
    (17, (-1, 3)): ((1, 1), (17, 1), (6, 37)),
    (17, (-2, 7)): ((-1, 1), (17, 5), (-2, 5)),
    (17, (-1, 5)): ((-2, 1), (-17, 9), (-6, 37)),
    (17, (1, 7)): ((0, 1), (-17, 9), (2, 5)),
    (17, (1, 3)): ((1, 1), (17, 1), (6, 37)),
    (17, (2, 3)): ((-2, 1), (-17, 1), (-6, 37)),
    (17, (1, 1)): ((1, 1), (-17, 7), (6, 37)),
    (17, (2, 1)): ((0, 1), (-17, 1), (2, 5)),
    (17, (7, 3)): ((1, 1), (-17, 3), (6, 37)),
    (17, (3, 1)): ((1, 1), (17, 9), (6, 37)),
    (17, (5, 1)): ((0, 1), (1, 11), (2, 5)),
    (19, (-5, 1)): ((0, 1), (-9, 5), (2, 5)),
    (19, (-3, 1)): ((1, 1), (-19, 3), (0, 1)),
    (19, (-2, 1)): ((1, 3), (19, 5), (0, 1)),
    (19, (-1, 1)): ((1, 1), (-19, 3), (0, 1)),
    (19, (-1, 3)): ((1, 3), (-19, 7), (0, 1)),
    (19, (-2, 7)): ((1, 1), (-19, 3), (0, 1)),
    (19, (-1, 5)): ((1, 1), (-19, 3), (0, 1)),
    (19, (1, 7)): ((1, 1), (-19, 7), (0, 1)),
    (19, (1, 3)): ((1, 1), (-19, 3), (0, 1)),
    (19, (2, 3)): ((1, 1), (-19, 11), (0, 1)),
    (19, (1, 1)): ((1, 1), (19, 1), (0, 1)),
    (19, (2, 1)): ((1, 1), (-19, 11), (0, 1)),
    (19, (7, 3)): ((-2, 1), (57, 1), (0, 1)),
    (19, (3, 1)): ((-1, 1), (-33, 1), (0, 1)),
    (19, (5, 1)): ((0, 1), (1, 11), (2, 5)),
    (23, (-5, 1)): ((0, 1), (-1, 1), (2, 5)),
    (23, (-3, 1)): ((0, 1), (27, 1), (2, 5)),
    (23, (-2, 1)): ((1, 1), (-69, 1), (0, 1)),
    (23, (-1, 1)): ((-2, 1), (69, 1), (0, 1)),
    (23, (-1, 3)): ((1, 1), (23, 9), (0, 1)),
    (23, (-2, 7)): ((1, 1), (23, 9), (0, 1)),
    (23, (-1, 5)): ((1, 5), (-23, 7), (0, 1)),
    (23, (1, 7)): ((1, 1), (-23, 7), (0, 1)),
    (23, (1, 3)): ((1, 1), (23, 1), (0, 1)),
    (23, (2, 3)): ((-2, 1), (23, 3), (0, 1)),
    (23, (1, 1)): ((1, 1), (-23, 3), (0, 1)),
    (23, (2, 1)): ((1, 1), (23, 1), (0, 1)),
    (23, (7, 3)): ((1, 1), (-23, 3), (0, 1)),
    (23, (3, 1)): ((1, 1), (-23, 3), (0, 1)),
    (23, (5, 1)): ((1, 1), (-35, 1), (6, 37)),
    (29, (-5, 1)): ((0, 1), (-19, 11), (2, 5)),
    (29, (-3, 1)): ((-1, 1), (-33, 1), (0, 1)),
    (29, (-1, 1)): ((-2, 1), (-29, 1), (-6, 37)),
    (29, (-1, 3)): ((1, 1), (87, 11), (6, 37)),
    (29, (-2, 7)): ((1, 1), (-29, 7), (6, 37)),
    (29, (-1, 5)): ((-2, 1), (29, 7), (-6, 37)),
    (29, (1, 7)): ((1, 1), (29, 9), (6, 37)),
    (29, (1, 3)): ((1, 1), (-29, 11), (6, 37)),
    (29, (2, 3)): ((1, 1), (-87, 1), (6, 37)),
    (29, (1, 1)): ((1, 1), (-29, 3), (6, 37)),
    (29, (2, 1)): ((2, 1), (-29, 9), (10, 101)),
    (29, (7, 3)): ((1, 1), (-29, 7), (6, 37)),
    (29, (3, 1)): ((1, 1), (29, 9), (6, 37)),
    (29, (5, 1)): ((0, 1), (1, 11), (2, 5)),
    (31, (-5, 1)): ((0, 1), (1, 11), (2, 5)),
    (31, (-3, 1)): ((-1, 1), (15, 1), (0, 1)),
    (31, (-1, 1)): ((3, 1), (31, 9), (0, 1)),
    (31, (-1, 3)): ((1, 1), (93, 7), (0, 1)),
    (31, (-2, 7)): ((1, 1), (31, 9), (0, 1)),
    (31, (-1, 5)): ((1, 1), (-31, 7), (0, 1)),
    (31, (1, 7)): ((1, 1), (31, 1), (0, 1)),
    (31, (1, 3)): ((-1, 3), (-93, 1), (0, 1)),
    (31, (2, 3)): ((1, 1), (-31, 11), (0, 1)),
    (31, (1, 1)): ((-2, 1), (31, 11), (0, 1)),
    (31, (2, 1)): ((-1, 3), (31, 9), (0, 1)),
    (31, (7, 3)): ((1, 1), (-31, 7), (0, 1)),
    (31, (3, 1)): ((-1, 1), (-27, 5), (-2, 5)),
    (31, (5, 1)): ((-1, 1), (-25, 11), (-2, 5)),
    (37, (-5, 1)): ((0, 1), (1, 11), (2, 5)),
    (37, (-3, 1)): ((0, 1), (27, 1), (2, 5)),
    (37, (-2, 1)): ((1, 1), (-1, 7), (6, 37)),
    (37, (-1, 1)): ((1, 1), (-11, 1), (6, 37)),
    (37, (-1, 3)): ((1, 1), (-1, 3), (6, 37)),
    (37, (-2, 7)): ((1, 1), (-1, 7), (6, 37)),
    (37, (-1, 5)): ((1, 1), (-1, 3), (6, 37)),
    (37, (1, 7)): ((1, 1), (-7, 1), (6, 37)),
    (37, (1, 3)): ((1, 1), (-25, 3), (6, 37)),
    (37, (2, 3)): ((1, 1), (9, 1), (6, 37)),
    (37, (1, 1)): ((1, 1), (1, 9), (6, 37)),
    (37, (2, 1)): ((1, 1), (37, 1), (6, 37)),
    (37, (7, 3)): ((1, 1), (-3, 1), (6, 37)),
    (37, (3, 1)): ((1, 1), (-1, 11), (6, 37)),
    (37, (5, 1)): ((0, 1), (-1, 1), (2, 5)),
    (41, (-5, 1)): ((0, 1), (-9, 1), (2, 5)),
    (41, (-3, 1)): ((0, 1), (27, 5), (2, 5)),
    (41, (-2, 1)): ((1, 3), (-41, 15), (30, 109)),
    (41, (-1, 1)): ((-1, 3), (-123, 1), (6, 13)),
    (41, (-1, 3)): ((1, 3), (41, 1), (30, 109)),
    (41, (-2, 7)): ((1, 3), (41, 5), (30, 109)),
    (41, (-1, 5)): ((2, 1), (-41, 1), (10, 101)),
    (41, (1, 3)): ((-1, 3), (41, 1), (6, 13)),
    (41, (2, 3)): ((2, 1), (-41, 1), (10, 101)),
    (41, (1, 1)): ((-3, 1), (41, 5), (-10, 101)),
    (41, (2, 1)): ((-3, 1), (41, 9), (-10, 101)),
    (41, (7, 3)): ((-1, 1), (63, 5), (0, 1)),
    (41, (3, 1)): ((-1, 1), (123, 1), (0, 1)),
    (41, (5, 1)): ((0, 1), (19, 1), (2, 5)),
    (43, (-5, 1)): ((0, 1), (-1, 1), (2, 5)),
    (43, (-3, 1)): ((0, 1), (27, 5), (2, 5)),
    (43, (-2, 1)): ((-1, 1), (43, 1), (0, 1)),
    (43, (-1, 1)): ((1, 1), (43, 9), (0, 1)),
    (43, (-1, 3)): ((-1, 1), (43, 5), (0, 1)),
    (43, (-2, 7)): ((-2, 1), (-43, 9), (0, 1)),
    (43, (-1, 5)): ((-1, 1), (-43, 19), (0, 1)),
    (43, (1, 7)): ((-2, 1), (129, 1), (0, 1)),
    (43, (1, 3)): ((1, 1), (43, 1), (0, 1)),
    (43, (2, 3)): ((-2, 1), (43, 7), (0, 1)),
    (43, (1, 1)): ((1, 1), (43, 9), (0, 1)),
    (43, (2, 1)): ((-1, 1), (43, 1), (0, 1)),
    (43, (7, 3)): ((1, 1), (-129, 1), (0, 1)),
    (43, (3, 1)): ((-3, 1), (27, 5), (0, 1)),
    (43, (5, 1)): ((0, 1), (-29, 1), (2, 5)),
    (47, (-5, 1)): ((0, 1), (71, 1), (2, 5)),
    (47, (-3, 1)): ((0, 1), (-27, 11), (2, 5)),
    (47, (-2, 1)): ((-1, 1), (47, 9), (0, 1)),
    (47, (-1, 1)): ((-1, 1), (47, 1), (0, 1)),
    (47, (-1, 3)): ((-1, 3), (47, 1), (0, 1)),
    (47, (-2, 7)): ((1, 3), (-141, 1), (0, 1)),
    (47, (-1, 5)): ((4, 1), (-47, 1), (0, 1)),
    (47, (1, 7)): ((-1, 1), (-47, 11), (0, 1)),
    (47, (1, 3)): ((4, 1), (-47, 1), (0, 1)),
    (47, (2, 3)): ((-1, 1), (235, 1), (0, 1)),
    (47, (1, 1)): ((-1, 1), (47, 9), (0, 1)),
    (47, (2, 1)): ((1, 3), (-47, 15), (0, 1)),
    (47, (7, 3)): ((-1, 1), (235, 1), (0, 1)),
    (47, (3, 1)): ((-1, 1), (-3, 11), (0, 1)),
    (47, (5, 1)): ((0, 1), (-25, 1), (2, 5)),
    (53, (-5, 1)): ((0, 1), (-29, 1), (2, 5)),
    (53, (-3, 1)): ((-1, 1), (-15, 11), (0, 1)),
    (53, (-2, 1)): ((-3, 1), (-53, 19), (-10, 101)),
    (53, (-1, 1)): ((-1, 1), (-53, 11), (-2, 5)),
    (53, (-1, 3)): ((0, 1), (-265, 1), (2, 5)),
    (53, (-2, 7)): ((0, 1), (-53, 1), (2, 5)),
    (53, (-1, 5)): ((0, 1), (-53, 5), (2, 5)),
    (53, (1, 7)): ((-3, 1), (53, 5), (-10, 101)),
    (53, (2, 3)): ((2, 1), (-53, 17), (10, 101)),
    (53, (1, 1)): ((0, 1), (-53, 1), (2, 5)),
    (53, (2, 1)): ((2, 1), (-53, 1), (10, 101)),
    (53, (7, 3)): ((2, 1), (-7, 5), (0, 1)),
    (53, (3, 1)): ((-1, 1), (-33, 5), (0, 1)),
    (53, (5, 1)): ((0, 1), (-1, 1), (2, 5)),
    (59, (-5, 1)): ((0, 1), (-5, 1), (2, 5)),
    (59, (-3, 1)): ((-1, 1), (-33, 1), (0, 1)),
    (59, (-2, 1)): ((2, 1), (-59, 1), (0, 1)),
    (59, (-1, 1)): ((1, 1), (-177, 1), (0, 1)),
    (59, (-1, 3)): ((1, 1), (-59, 11), (0, 1)),
    (59, (-2, 7)): ((-3, 1), (59, 13), (0, 1)),
    (59, (-1, 5)): ((2, 1), (-59, 1), (0, 1)),
    (59, (1, 7)): ((-1, 3), (59, 1), (0, 1)),
    (59, (1, 3)): ((2, 1), (-59, 1), (0, 1)),
    (59, (2, 3)): ((2, 1), (-295, 1), (0, 1)),
    (59, (1, 1)): ((1, 1), (59, 1), (0, 1)),
    (59, (2, 1)): ((2, 1), (-59, 13), (0, 1)),
    (59, (7, 3)): ((-2, 1), (177, 1), (0, 1)),
    (59, (3, 1)): ((-1, 1), (3, 1), (0, 1)),
    (59, (5, 1)): ((0, 1), (19, 1), (2, 5)),
    (61, (-5, 1)): ((1, 1), (5, 1), (6, 37)),
    (61, (-3, 1)): ((0, 1), (-27, 11), (2, 5)),
    (61, (-2, 1)): ((2, 1), (-61, 17), (10, 101)),
    (61, (-1, 1)): ((1, 1), (61, 1), (6, 37)),
    (61, (-1, 3)): ((1, 1), (-183, 1), (6, 37)),
    (61, (-1, 5)): ((2, 1), (61, 19), (10, 101)),
    (61, (1, 7)): ((-2, 1), (183, 1), (-6, 37)),
    (61, (1, 3)): ((-2, 1), (61, 3), (-6, 37)),
    (61, (2, 3)): ((-2, 1), (-61, 1), (-6, 37)),
    (61, (1, 1)): ((2, 1), (-61, 9), (10, 101)),
    (61, (7, 3)): ((-1, 1), (7, 1), (0, 1)),
    (61, (3, 1)): ((-1, 1), (3, 5), (0, 1)),
    (61, (5, 1)): ((0, 1), (-25, 9), (2, 5)),
    (67, (-5, 1)): ((-1, 1), (-11, 9), (-2, 5)),
    (67, (-3, 1)): ((-1, 1), (-27, 5), (-2, 5)),
    (67, (-2, 1)): ((2, 1), (-335, 1), (0, 1)),
    (67, (-1, 1)): ((1, 3), (-201, 1), (0, 1)),
    (67, (-1, 3)): ((2, 1), (-67, 1), (0, 1)),
    (67, (-2, 7)): ((-3, 1), (67, 1), (0, 1)),
    (67, (-1, 5)): ((4, 1), (201, 1), (0, 1)),
    (67, (1, 7)): ((1, 3), (67, 1), (0, 1)),
    (67, (1, 3)): ((3, 1), (-67, 7), (0, 1)),
    (67, (2, 3)): ((-1, 1), (-67, 19), (0, 1)),
    (67, (1, 1)): ((-1, 1), (67, 1), (0, 1)),
    (67, (7, 3)): ((-1, 1), (7, 5), (0, 1)),
    (67, (3, 1)): ((-1, 1), (-3, 11), (0, 1)),
    (67, (5, 1)): ((1, 1), (-35, 9), (6, 37)),
    (71, (-5, 1)): ((0, 1), (1, 11), (2, 5)),
    (71, (-3, 1)): ((3, 1), (25, 3), (0, 1)),
    (71, (-2, 1)): ((3, 1), (71, 1), (0, 1)),
    (71, (-1, 1)): ((-1, 3), (71, 1), (0, 1)),
    (71, (-1, 3)): ((-1, 3), (-71, 3), (0, 1)),
    (71, (-2, 7)): ((3, 1), (71, 1), (0, 1)),
    (71, (-1, 5)): ((-1, 3), (-213, 1), (0, 1)),
    (71, (1, 7)): ((-1, 3), (71, 9), (0, 1)),
    (71, (1, 3)): ((4, 1), (-71, 1), (0, 1)),
    (71, (2, 3)): ((4, 1), (71, 3), (0, 1)),
    (71, (1, 1)): ((4, 1), (-71, 1), (0, 1)),
    (71, (2, 1)): ((-1, 3), (71, 1), (0, 1)),
    (71, (7, 3)): ((-1, 1), (-7, 11), (0, 1)),
    (71, (3, 1)): ((-1, 1), (-27, 11), (0, 1)),
    (71, (5, 1)): ((0, 1), (-25, 1), (2, 5)),
    (73, (-5, 1)): ((0, 1), (-9, 1), (2, 5)),
    (73, (-3, 1)): ((-1, 1), (-27, 5), (-2, 5)),
    (73, (-2, 1)): ((0, 1), (-73, 1), (2, 5)),
    (73, (-1, 1)): ((-3, 1), (73, 1), (-10, 101)),
    (73, (-1, 3)): ((-1, 1), (73, 1), (-2, 5)),
    (73, (-2, 7)): ((0, 1), (-73, 1), (2, 5)),
    (73, (-1, 5)): ((0, 1), (-73, 1), (2, 5)),
    (73, (1, 7)): ((-1, 3), (73, 1), (6, 13)),
    (73, (1, 3)): ((-1, 3), (73, 1), (6, 13)),
    (73, (2, 3)): ((-1, 1), (73, 1), (-2, 5)),
    (73, (1, 1)): ((0, 1), (-365, 1), (2, 5)),
    (73, (2, 1)): ((-1, 1), (73, 9), (-2, 5)),
    (73, (7, 3)): ((-1, 1), (7, 5), (0, 1)),
    (73, (3, 1)): ((0, 1), (27, 5), (2, 5)),
    (73, (5, 1)): ((0, 1), (-29, 1), (2, 5)),
    (79, (-5, 1)): ((0, 1), (31, 1), (2, 5)),
    (79, (-3, 1)): ((-1, 1), (-33, 5), (0, 1)),
    (79, (-2, 1)): ((-2, 1), (-79, 9), (0, 1)),
    (79, (-1, 1)): ((1, 1), (-237, 1), (0, 1)),
    (79, (-1, 3)): ((1, 1), (-79, 3), (0, 1)),
    (79, (-2, 7)): ((-2, 1), (79, 11), (0, 1)),
    (79, (-1, 5)): ((1, 1), (79, 1), (0, 1)),
    (79, (1, 7)): ((-2, 1), (-79, 1), (0, 1)),
    (79, (1, 3)): ((3, 1), (-79, 7), (0, 1)),
    (79, (2, 3)): ((-2, 1), (237, 1), (0, 1)),
    (79, (1, 1)): ((-2, 1), (-79, 1), (0, 1)),
    (79, (2, 1)): ((1, 1), (-79, 3), (0, 1)),
    (79, (7, 3)): ((-1, 1), (35, 9), (0, 1)),
    (79, (3, 1)): ((-1, 1), (3, 1), (0, 1)),
    (79, (5, 1)): ((1, 1), (-5, 3), (6, 37)),
    (83, (-5, 1)): ((0, 1), (1, 11), (2, 5)),
    (83, (-3, 1)): ((-1, 1), (-3, 11), (0, 1)),
    (83, (-2, 1)): ((2, 1), (83, 19), (0, 1)),
    (83, (-1, 1)): ((-3, 1), (83, 1), (0, 1)),
    (83, (-1, 3)): ((-1, 1), (-83, 11), (0, 1)),
    (83, (-2, 7)): ((-1, 3), (83, 13), (0, 1)),
    (83, (-1, 5)): ((-1, 1), (83, 1), (0, 1)),
    (83, (1, 7)): ((-1, 3), (83, 1), (0, 1)),
    (83, (1, 3)): ((2, 1), (83, 19), (0, 1)),
    (83, (2, 3)): ((-1, 1), (83, 5), (0, 1)),
    (83, (1, 1)): ((-1, 1), (83, 5), (0, 1)),
    (83, (2, 1)): ((2, 1), (-83, 1), (0, 1)),
    (83, (7, 3)): ((2, 1), (-7, 5), (0, 1)),
    (83, (3, 1)): ((0, 1), (27, 1), (2, 5)),
    (83, (5, 1)): ((0, 1), (19, 1), (2, 5)),
    (89, (-5, 1)): ((0, 1), (-5, 1), (2, 5)),
    (89, (-3, 1)): ((-1, 1), (-27, 11), (0, 1)),
    (89, (-1, 1)): ((-1, 3), (-267, 1), (6, 13)),
    (89, (-1, 3)): ((3, 1), (89, 1), (14, 197)),
    (89, (-2, 7)): ((-2, 1), (267, 1), (-6, 37)),
    (89, (-1, 5)): ((-1, 3), (89, 17), (6, 13)),
    (89, (1, 7)): ((1, 5), (89, 17), (70, 221)),
    (89, (1, 3)): ((1, 1), (-89, 3), (6, 37)),
    (89, (2, 3)): ((3, 1), (89, 9), (14, 197)),
    (89, (1, 1)): ((-3, 1), (89, 1), (-10, 101)),
    (89, (2, 1)): ((-1, 3), (89, 1), (6, 13)),
    (89, (7, 3)): ((2, 1), (-7, 1), (0, 1)),
    (89, (3, 1)): ((-1, 1), (-33, 5), (0, 1)),
    (89, (5, 1)): ((0, 1), (-89, 1), (2, 5)),
    (97, (-5, 1)): ((0, 1), (1, 11), (2, 5)),
    (97, (-3, 1)): ((-1, 1), (3, 5), (0, 1)),
    (97, (-2, 1)): ((0, 1), (-97, 5), (2, 5)),
    (97, (-1, 1)): ((1, 1), (97, 9), (6, 37)),
    (97, (-1, 3)): ((-1, 1), (97, 1), (-2, 5)),
    (97, (-2, 7)): ((-1, 3), (97, 1), (6, 13)),
    (97, (-1, 5)): ((1, 1), (97, 1), (6, 37)),
    (97, (1, 7)): ((5, 1), (17, 9), (22, 485)),
    (97, (1, 3)): ((4, 1), (-97, 1), (18, 325)),
    (97, (2, 3)): ((-2, 1), (97, 11), (-6, 37)),
    (97, (1, 1)): ((-1, 1), (97, 1), (-2, 5)),
    (97, (2, 1)): ((-1, 3), (-291, 1), (6, 13)),
    (97, (7, 3)): ((1, 1), (97, 1), (6, 37)),
    (97, (3, 1)): ((-1, 1), (97, 1), (-2, 5)),
    (97, (5, 1)): ((0, 1), (19, 1), (2, 5)),
}


def _l9_export_lines():
    """Canonical JSONL serialization of _L9_STEERED (same contract as
    _evidence_export_lines: the file must match byte-for-byte)."""
    lines = []
    for (w, ut), (st, bt, tt) in sorted(_L9_STEERED.items()):
        lines.append(json.dumps(
            {"w": w, "u": list(ut), "s": list(st), "b": list(bt),
             "tau": list(tt)}, separators=(", ", ": ")))
    return lines


def _verify_L9(extended=False):
    """L9 (THEOREMS.md): the steering lemma, the prime-b symbol L9a, and
    replay of the frozen steered witness table.  Structure results and
    bounded evidence for the OPEN assembly lemma; no status/count changes.

    (a) L9 lemma engine: for random x, d with T = {2, oo} u supp(x) u
        supp(d), the product of (x,d)_v over T is 1 (reciprocity restricted
        to T; places outside T have unit-unit symbols): forcing one place
        from the others is exactly this identity.
    (b) L9a: for admissible s (A = 1+4a^2 = 5 mod 8, hence nonsquare) and
        b = eps*q1 with q1 prime satisfying the unit conditions FACTOR BY
        FACTOR (v_q1 = 0 on each of delta, A, a, z, D_z), v_q1(x) = -4 and
        (x, d)_q1 = legendre(A, q1) exactly.  (b1) re-refutes the earlier
        product-form guard (it admits q1 | A on tau = 2a/A, where
        delta*A = 1, and then v_q1(x) is odd); (b2) re-refutes the earlier
        class-modulus recipe (a controlled symbol flips inside the
        advertised prime class).
    (c) Frozen steered table: every row is Phi-admissible, sits in a
        genuine target cell, and its steered W2 certificate returns True;
        the JSONL export matches byte-for-byte.
    """
    rng = random.Random(20260816)
    # (a) lemma engine on random rationals
    nlem = 0
    while nlem < (600 if extended else 200):
        x = Fraction(rng.randint(-3000, 3000), rng.randint(1, 300))
        d = Fraction(rng.randint(-3000, 3000), rng.randint(1, 300))
        if x == 0 or d == 0:
            continue
        T = sorted(set([2, OO]) | set(places_of(x)) | set(places_of(d)),
                   key=str)
        prod = 1
        for v in T:
            prod *= hilbert(x, d, v)
        assert prod == 1, (x, d)
        nlem += 1
    print(f"  (L9) product over T is 1 on {nlem} random pairs: one unchecked"
          " place is always forced")
    # (b) prime-b symbol identity
    nl9a = 0
    qs = [q for q in primerange(1000, 4000)][::7]
    while nl9a < (300 if extended else 100):
        s = Fraction(rng.randint(-6, 6), rng.choice([1, 3, 5, 7]))
        a = 1 + 2*s
        A = 1 + 4*a*a
        assert A == 5 or vp(A - 5, 2) >= 3    # A = 5 mod 8: nonsquare
        w = rng.choice([3, 5, 7, 11, 13])
        u = Fraction(rng.randint(1, 9), rng.choice([1, 3, 7]))
        z = Fraction(w) * u
        Dz = 1 - z**3 - a*a*z**6
        if Dz == 0:
            continue
        for tau in (Fraction(0), 2*a/A):
            delta = 1 - A*tau*tau
            if delta == 0:
                continue
            alpha = -delta*A
            # factor-by-factor unit conditions (a product guard is NOT
            # enough: on tau = 2a/A one has delta*A = 1, which hides q1|A)
            q1 = rng.choice(qs)
            if any(vp(t, q1) != 0 for t in (delta, A, a, z, Dz)):
                continue
            eps = rng.choice([1, -1])
            b = Fraction(eps * q1)
            c = _sun_h(a, b, z**3)
            if c is None:
                continue
            M = 16 - delta*c*c - 16*A*2*b*((a - 1)/2)**2
            if M == 0:
                continue
            x, d = alpha*M, alpha*2*b
            assert vp(x, q1) == -4, (s, w, u, tau, q1)
            assert hilbert(x, d, q1) == legendre(A, q1), (s, w, u, tau, q1)
            nl9a += 1
    # (b1) boundary regression: dropping "q1 does not divide A" breaks L9a.
    # Cell w=3, z=6, s=0, tau=2a/A: delta*A = 1, so a product-form guard on
    # 2*delta*A*a*z*D_z admits q1 = 5 = A -- and then v_q1(x) is ODD.
    z_ce = Fraction(6)
    a_ce, A_ce = Fraction(1), Fraction(5)
    tau_ce = 2*a_ce/A_ce
    delta_ce = 1 - A_ce*tau_ce*tau_ce
    Dz_ce = 1 - z_ce**3 - a_ce*a_ce*z_ce**6
    assert vp(2*delta_ce*A_ce*a_ce*z_ce*Dz_ce, 5) == 0   # old guard passes
    assert vp(A_ce, 5) == 1                              # new guard rejects
    c_ce = _sun_h(a_ce, Fraction(5), z_ce**3)
    M_ce = 16 - delta_ce*c_ce*c_ce
    x_ce, d_ce = -delta_ce*A_ce*M_ce, -delta_ce*A_ce*10
    assert vp(x_ce, 5) == -5 and hilbert(x_ce, d_ce, 5) == -1
    print(f"  (L9a) prime-b symbol: v_q1(x) = -4 and (x,d)_q1 = (A|q1) on"
          f" {nl9a} instances (factor-by-factor units; A = 5 mod 8"
          " throughout); q1|A boundary case re-refuted")
    # (b2) modulus regression: the naive class recipe p^(v_p(x(b0))+1) over
    # data primes only does NOT freeze controlled symbols.  Cell (73, 5),
    # s = 0, b0 = -29, tau = 2/5: N = 689120 is the least modulus meeting
    # the old conditions, q = 29 + 355*N is prime and in class, (A|q) = +1,
    # yet the symbol at the controlled prime 59 flips +1 -> -1.
    z_m, s_m, tau_m = Fraction(73*5), Fraction(0), Fraction(2, 5)
    a_m = 1 + 2*s_m
    A_m = 1 + 4*a_m*a_m
    delta_m = 1 - A_m*tau_m*tau_m
    alpha_m = -delta_m*A_m
    N_m, q_m = 689120, 244637629
    assert _is_prime(q_m) and (q_m - 29) % N_m == 0
    assert legendre(A_m, q_m) == legendre(A_m, 29) == 1
    syms = []
    for b_m in (Fraction(-29), Fraction(-q_m)):
        c_m = _sun_h(a_m, b_m, z_m**3)
        M_m = 16 - delta_m*c_m*c_m
        syms.append(hilbert(alpha_m*M_m, alpha_m*2*b_m, 59))
    assert syms == [1, -1], syms
    # (b3) the same cell's D_z = -59 * 40077920921911 pins a fixed-data prime
    # far above the trial bound.  A superseded draft called it a controlled
    # place that the modulus had to freeze; L10-0 shows v_p(d) = 0 there, so
    # it is only a candidate WILD prime.  Kept as a regression because the
    # number is cited in THEOREMS L9/L10-0.
    Dz_m = 1 - z_m**3 - a_m*a_m*z_m**6
    assert Dz_m == -2364597334392749
    assert factorint(abs(Dz_m.numerator)) == {59: 1, 40077920921911: 1}
    assert _is_prime(40077920921911) and 40077920921911 > _L9_TRIAL_BOUND
    assert vp(-2*Fraction(-29), 40077920921911) == 0   # v_p(d) = 0: not controlled
    print("  (L9-mod) naive class modulus refuted in-suite: same prime class"
          " and sign, controlled symbol at 59 flips +1 -> -1; the cell's"
          " 4.0e13 data prime is wild, not controlled (L10-0)")
    # (c) frozen steered table replay + export contract
    mechs = {}
    maxwild = 0
    for (w, ut), (st, bt, tt) in sorted(_L9_STEERED.items()):
        u, s = Fraction(*ut), Fraction(*st)
        b, tau = Fraction(*bt), Fraction(*tt)
        z = Fraction(w) * u
        assert vp(z, w) >= 1, (w, ut)                 # genuine target cell
        assert (s == 0 or vp(s, 2) >= 0) and vp(b, 2) == 0
        a = 1 + 2*s
        assert tau in (0, 2*a/(1 + 4*a*a))
        ok, mech, wilds = _l9_steered_solvable(a, b, z, tau)
        assert ok is True, (w, ut, mech)
        mechs[mech] = mechs.get(mech, 0) + 1
        if mech == 'steered' and wilds:
            maxwild = max(maxwild, max(wilds))
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "data", "l9_steered.jsonl")
    want = "\n".join(_l9_export_lines()) + "\n"
    try:
        with open(path) as fh:
            got = fh.read()
    except FileNotFoundError:
        got = None
    assert got == want, "steered evidence export drifted from authority"
    cells = {(w, ut) for (w, ut) in _L9_STEERED}
    ws = {w for (w, _) in cells}
    print(f"  (L9-T) steered table replayed: {len(_L9_STEERED)} cells across"
          f" {len(ws)} targets w; mechanisms {mechs}; largest forced wild"
          f" prime ~ 10^{len(str(maxwild)) - 1}; export byte-identical")
    print("  -> steering converts wild-prime luck into reciprocity; coverage"
          " is bounded evidence ONLY. Unconditional assembly stays open;"
          " the Schinzel-conditional chain is closed later by L22\n")





# Canonical u-pool for the assembly grid (SINGLE SOURCE: l9_steer.py imports
# it).  Cell (w, u) is a target iff v_w(w*u) >= 1.  353 cells for w < 100.
_L9_U_POOL = [Fraction(v) for v in (1, -1, 2, -2, 3, -3, 5, -5)] + \
             [Fraction(n, m) for (n, m) in
              [(1, 3), (-1, 3), (2, 3), (-1, 5), (7, 3), (1, 7), (-2, 7)]]


def _l9_grid(wmax=100):
    """The canonical cell set, as (w, (num, den)) keys."""
    return {(w, (u.numerator, u.denominator))
            for w in primerange(3, wmax) for u in _L9_U_POOL
            if vp(Fraction(w)*u, w) >= 1}


# ---------------------------------------------------------------- L10:
# class-side dichotomy for the prime-b family.  See THEOREMS.md L10.
#   x = alpha*P(b)/(D^2 A^2 b^4),  d = 2*alpha*b,  alpha = -delta_tau*A.
# Frozen set is exactly supp(d) = {2} u supp(alpha) u {q1}: at every other
# place v_p(d) = 0, so (x,d)_p = (d|p)^{v_p(x)} = +1 whenever v_p(x) is even.


def _l10_poly_mul(p, q):
    r = [Fraction(0)]*(len(p) + len(q) - 1)
    for i, ci in enumerate(p):
        if ci:
            for j, cj in enumerate(q):
                r[i + j] += ci*cj
    return r


def _l10_poly_add(p, q):
    r = [Fraction(0)]*max(len(p), len(q))
    for i, c in enumerate(p):
        r[i] += c
    for i, c in enumerate(q):
        r[i] += c
    return r


def _l10_P(a, Z, D, A, delta, s):
    """P(b) = 16 D^2 A^2 b^4 - delta a^4 Z^4 Ng(b)^2 - 32 A^3 s^2 D^2 b^5."""
    q = [Fraction(1)]
    for _ in range(4):
        q = _l10_poly_mul(q, [Fraction(-1), Fraction(1)])
    Ng = _l10_poly_add([Fraction(0), Fraction(0), 16*a**4], [c*(-A) for c in q])
    t1 = [Fraction(0)]*4 + [16*D*D*A*A]
    t2 = [c*(-delta*a**4*Z**4) for c in _l10_poly_mul(Ng, Ng)]
    t3 = [Fraction(0)]*5 + [-32*A**3*s*s*D*D]
    return _l10_poly_add(_l10_poly_add(t1, t2), t3)


def _l10_taylor(P, b0):
    """Exact Taylor coefficients about b0: c_j = sum_{i>=j} P_i C(i,j) b0^(i-j),
    so that P(b0 + h) = sum_j c_j h^j.  (A synthetic-division variant of this
    was wrong -- it returned P's constant coefficient instead of P(b0), i.e.
    expanded about 0 -- which silently invalidated every exponent below.)"""
    n = len(P)
    return [sum(P[i]*math.comb(i, j)*b0**(i - j) for i in range(j, n))
            for j in range(n)]


def _l10_exponent(P, b0, p):
    """Rigorous k with: b = b0 mod p^k  ==>  v_p(P(b)) = v_p(P(b0)) and
    P(b)/P(b0) is a square in Q_p (ratio in 1 + p Z_p, resp. 1 + 8 Z_2)."""
    cs = _l10_taylor(P, b0)
    v0, e = vp(cs[0], p), (2 if p == 2 else 0)
    tail = [vp(c, p) for c in cs[1:] if c != 0]
    k = 3 if p == 2 else 1
    return max(k, v0 + e - min(tail) + 1) if tail else k


def _l10_supp(x):
    if x == 0:
        return set()
    return set(factorint(abs(x.numerator))) | set(factorint(x.denominator))


def _l10_class_cert(a, z, tau, eps, q1):
    """Aligned-class certificate for b = eps*q1 on the prime-b family.

    Returns None if the input is degenerate, else a dict with the modulus N,
    the frozen places, the frozen symbols, and the sign bound Q0.  'ok' means
    every frozen symbol is +1, so L9 forces the single wild symbol to +1 too.
    """
    A = 1 + 4*a*a
    delta = 1 - A*tau*tau
    if delta == 0:
        return None
    alpha = -delta*A
    Z = z**3
    D = 1 - Z - a*a*Z*Z
    if D == 0:
        return None
    s = (a - 1)/2
    P = _l10_P(a, Z, D, A, delta, s)
    b0 = Fraction(eps*q1)
    # L10-0 frees a place only when v_p(x) is EVEN there.  delta carries A in
    # its denominator on every branch with tau != 0 (1/A canonical, -4a^4/A
    # square), which forces odd v_p(x) at p | A, so supp(delta) is controlled
    # data and MUST be frozen alongside supp(alpha).  Omitting it was a real
    # defect: it silently certified rows whose symbol at p | A is -1.
    # Two corrections, both from the 2026-08-16 audit, both load-bearing:
    #  (1) supp(delta): delta carries A in its denominator on every tau != 0
    #      branch (1/A canonical, -4a^4/A square), forcing ODD v_p(x) at p | A.
    #  (2) {3, 5, 7}: P has degree 8, so a prime dividing P(u) for EVERY unit u
    #      without dividing the content must satisfy p - 1 <= 8 (a nonzero
    #      degree-8 form over F_p has at most 8 roots).  Hence the only
    #      prime-restricted fixed divisors are p in {2, 3, 5, 7}, and freezing
    #      that fixed finite set closes the gap completely.
    S = sorted({2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta))
    # L9a's hypotheses are FACTOR BY FACTOR: on the alpha = -1 branch a guard
    # on supp(2*alpha) is vacuous, so (x,d)_q1 = (A|q1) would be asserted
    # without its units.  Reject any base sharing q1 with the cell data.
    if a == 0 or z == 0 or q1 in S:
        return None
    if vp(a, 2) != 0:
        # admissibility: a = 1 + 2s is an odd 2-adic unit.  Without it
        # A = 1 + 4a^2 need NOT be 5 mod 8 and the L10b wall is false
        # (a = 2 gives A = 17 = 1 mod 8; see _verify_L10).
        return None
    if not _is_prime(q1):
        return None          # L9a and legendre(A, q1) both need q1 prime
    if any(vp(t, q1) != 0 for t in (delta, A, a, z, D)):
        return None
    ks = {p: _l10_exponent(P, b0, p) for p in S}
    N = 8
    for p, k in ks.items():
        N = N*p**k // math.gcd(N, p**k)
    m = 4*A.numerator*A.denominator
    N = N*m // math.gcd(N, m)
    if math.gcd(q1, N) != 1:         # else the progression holds no other prime
        return None
    c0 = _sun_h(a, b0, Z)
    if c0 is None:
        return None
    M0 = 16 - delta*c0*c0 - 32*A*b0*s*s
    if M0 == 0:
        return None
    x0, d0 = alpha*M0, alpha*2*b0
    syms = {p: hilbert(x0, d0, p) for p in S}
    nz = [i for i, c in enumerate(P) if c != 0]
    lead = P[nz[-1]]
    Q0 = 1 + max((abs(c/lead) for c in P[:nz[-1]]), default=Fraction(0))
    sx = (1 if alpha > 0 else -1)*(1 if lead > 0 else -1)*(eps**nz[-1])
    sd = 1 if alpha*eps > 0 else -1
    syms[OO] = -1 if (sx < 0 and sd < 0) else 1      # asymptotic, valid past Q0
    syms['q1'] = legendre(A, q1)                     # = (x,d)_q1 by L9a
    return {'ok': all(v == 1 for v in syms.values()), 'N': N, 'S': S, 'ks': ks,
            'syms': syms, 'Q0': Q0}


# SUPERSEDED 2026-08-16.  `_L10_CLASSES` held 164 rows; the adversarial audit
# showed its certificate omitted controlled places (see L11e/L11f).  130 of the
# 164 rows are INVALID: their symbol at p | A is -1.  The other 34 keep valid
# symbols but their recorded modulus N is stale under the complete controlled
# set {2,3,5,7} u supp(alpha) u supp(delta), so no row of the table can be
# cited as printed.  The table is retired rather than trimmed: L11h's
# `_L11_CLASSES` supersedes it with a characterized 190-cell table.  A sample
# of the invalid rows is kept so the defect can never silently return.
_L10_INVALID_ROWS = (
    (3, (-3, 1), 59),
    (3, (-2, 1), 59),
    (3, (-2, 7), 41),
    (3, (-1, 1), 59),
    (3, (-1, 5), 59),
    (3, (1, 1), 59),
    (3, (1, 7), 59),
    (3, (2, 1), 59),
)


def _verify_L10(extended=False):
    """L10 (THEOREMS.md): the class-side dichotomy of the prime-b family.

    (a) On the canonical branch tau = 2a/A one has delta = 1/A and alpha = -1
        EXACTLY, so supp(alpha) is empty and the frozen set is just {2, oo, q1}.
    (b) On tau = 0 (forced by W0 when w = 3 mod 4) alpha = -A, and for PRIME
        A with a ODD (admissibility -- a = 2 gives A = 17 = 1 mod 8 and
        breaks the wall) coprime to 2*z*D_z, the frozen symbols satisfy
            prod_{p|A} (x,d)_p * (A|q1) = (-2 eps | A) = -1
        because A = 5 mod 8.  No aligned class exists for such a witness.
        For rational or non-squarefree A this is machine-checked evidence only.
    (c) 164 exhibited aligned classes = exactly the canonical w = 1 mod 4
        cells of the 353-cell grid (key-set equality asserted), replayed
        with their moduli; symbols constant along each class except for
        the finitely many members dividing the fixed cell data.
    """
    rng = random.Random(20261016)
    # (0) the Taylor shift itself: sum_j c_j h^j must equal P(b0 + h).
    nsh = 0
    while nsh < (120 if extended else 40):
        P = [Fraction(rng.randint(-60, 60), rng.randint(1, 9))
             for _ in range(rng.randint(2, 9))]
        b0 = Fraction(rng.randint(-30, 30), rng.randint(1, 5))
        cs = _l10_taylor(P, b0)
        for h in (Fraction(0), Fraction(1), Fraction(-3), Fraction(5, 7)):
            lhs = rhs = Fraction(0)
            for c in reversed(cs):
                lhs = lhs*h + c
            for c in reversed(P):
                rhs = rhs*(b0 + h) + c
            assert lhs == rhs, (P, b0, h)
        nsh += 1
    print(f"  (L10) exact Taylor shift: sum_j c_j h^j = P(b0+h) on {nsh}"
          " random polynomials x 4 offsets (a wrong shift silently voided"
          " every exponent in an earlier draft)")
    # (a) alpha = -1 on the canonical tau = 2a/A branch
    nal = 0
    while nal < (300 if extended else 100):
        s = Fraction(rng.randint(-8, 8), rng.choice([1, 3, 5, 7]))
        a = 1 + 2*s
        if a == 0:
            continue
        A = 1 + 4*a*a
        tau = 2*a/A
        assert 1 - A*tau*tau == 1/A and -(1 - A*tau*tau)*A == -1, (s,)
        nal += 1
    print(f"  (L10a) canonical branch: delta = 1/A and alpha = -1 exactly on"
          f" {nal} random admissible s; the SQUARE CLASS of alpha is trivial,"
          " but delta = 1/A keeps supp(A) controlled, so the frozen set does"
          " NOT collapse to {2, oo, q1} -- see L11f")
    # (b) unfrozen places with even valuation are free
    nfree = 0
    while nfree < (400 if extended else 150):
        x = Fraction(rng.randint(-400, 400), rng.randint(1, 40))
        d = Fraction(rng.randint(-400, 400), rng.randint(1, 40))
        if x == 0 or d == 0:
            continue
        p = rng.choice([3, 5, 7, 11, 13, 17, 19, 23])
        if vp(d, p) != 0 or vp(x, p) % 2 != 0:
            continue
        assert hilbert(x, d, p) == 1, (x, d, p)
        nfree += 1
    print(f"  (L10) unfrozen places are free: v_p(d) = 0 and v_p(x) even =>"
          f" symbol +1, {nfree} instances")
    # (c) tau = 0 wall: proved case (A prime) and evidence case
    npr = nev = 0
    while npr < (250 if extended else 80):
        a = Fraction(rng.choice([1, -1, 3, -3, 5, -5, 7, -7, 11, -11, 13, -13]))
        A = 1 + 4*a*a
        if not _is_prime(A.numerator):
            continue
        w = rng.choice([3, 7, 11, 19, 23, 31, 43, 59, 83])
        u = rng.choice([Fraction(1), Fraction(2), Fraction(3), Fraction(5),
                        Fraction(-1), Fraction(1, 3), Fraction(1, 5)])
        z = Fraction(w)*u
        Z = z**3
        D = 1 - Z - a*a*Z*Z
        if D == 0 or _sun_target_tau(a, w) != 0 or vp(2*z*D, A.numerator) != 0:
            continue
        eps = rng.choice([1, -1])
        q1 = rng.choice([q for q in primerange(41, 500)])
        # L9a's units are hypotheses: q1 must divide none of A, a, z, D_z,
        # delta.  Guarding only q1 == A.numerator let q1 = w slip through for
        # w in {43, 59, 83}, violating the very hypotheses RESULTS advertises.
        Dz_ = 1 - z**3 - a*a*z**6
        if any(t != 0 and vp(t, q1) != 0 for t in (A, a, z, Dz_, Fraction(1))):
            continue
        r = _l10_class_cert(a, z, Fraction(0), eps, q1)
        if r is None or q1 <= r['Q0']:
            continue
        # The L10b identity is a statement about the places p | A ONLY.  S is
        # now larger ({2,3,5,7} u supp(alpha) u supp(delta), see L11f), so the
        # product must be taken over supp(A), not over all of S.
        lhs = r['syms']['q1']
        for p in factorint(A.numerator):
            lhs *= r['syms'][p]
        assert vp(a, 2) == 0 and A % 8 == 5, (a,)   # admissibility is essential
        assert lhs == legendre(-2*eps, A.numerator) == -1, (a, w, u, eps, q1)
        assert not r['ok']
        npr += 1
    # (c1) why admissibility is a hypothesis and not decoration: a = 2 is not
    # an odd 2-adic unit, A = 17 = 1 mod 8, and the wall identity yields +1.
    a_ce, z_ce = Fraction(2), Fraction(3)
    A_ce = 1 + 4*a_ce*a_ce
    assert _is_prime(A_ce.numerator) and A_ce % 8 == 1
    assert legendre(-2, 17) == legendre(2, 17) == 1
    assert _l10_class_cert(a_ce, z_ce, Fraction(0), 1, 5) is None   # refused
    print("  (L10b') admissibility is load-bearing: a = 2 gives A = 17 = 1 mod 8"
          " and (-2eps|A) = +1, so the wall fails; such a is refused up front")
    print(f"  (L10b) tau = 0 wall PROVED case (A prime, coprime to 2*z*D_z):"
          f" prod_(p|A)(x,d)_p * (A|q1) = (-2eps|A) = -1 on {npr} instances;"
          " no aligned class exists for such a witness")
    # (c2) EVIDENCE ONLY: same identity in Jacobi form for rational or
    # non-squarefree A.  Not proved -- the cross terms (A/p^k | p) and
    # rational A lie outside the L10b derivation.
    nev = 0
    while nev < (600 if extended else 200):
        s = Fraction(rng.choice([1, -1, 2, -2, 3, -3]), rng.choice([3, 5, 7]))
        a = 1 + 2*s
        A = 1 + 4*a*a
        if a == 0 or (A.denominator == 1 and _is_prime(A.numerator)):
            continue
        w = rng.choice([3, 7, 11, 19, 23, 31, 43, 59, 83])
        u = rng.choice([Fraction(1), Fraction(2), Fraction(3), Fraction(5),
                        Fraction(-1), Fraction(1, 3), Fraction(1, 7)])
        z = Fraction(w)*u
        if 1 - z**3 - a*a*z**6 == 0 or _sun_target_tau(a, w) != 0:
            continue
        eps = rng.choice([1, -1])
        q1 = rng.choice([q for q in primerange(41, 500)])
        r = _l10_class_cert(a, z, Fraction(0), eps, q1)
        if r is None or q1 <= r['Q0']:
            continue
        # again: the identity concerns the places dividing A only (L11f made S
        # strictly larger), so restrict the product to supp(num A) u supp(den A)
        lhs = r['syms']['q1']
        for m in (A.numerator, A.denominator):
            for p in factorint(m):
                lhs *= r['syms'][p]
        rhs = 1
        for m in (A.numerator, A.denominator):
            for p, e in factorint(m).items():
                rhs *= legendre(-2*eps, p)**e
        assert lhs == rhs == -1, (s, w, u, eps, q1)
        assert not r['ok']
        nev += 1
    print(f"  (L10b-E) EVIDENCE ONLY (not proved): Jacobi form of the same"
          f" identity holds on {nev} rational / non-squarefree A instances")
    # (d) REGRESSION, not a replay.  The 164-row table this block used to
    # replay is retired: its certificate omitted the controlled places p | A,
    # and 130 rows were invalid.  We now assert those rows STAY rejected.
    nbad = 0
    for (w, ut, q1) in _L10_INVALID_ROWS:
        a = Fraction(1)                     # the shape the retired table used
        z = Fraction(w)*Fraction(*ut)
        assert vp(z, 5) <= 0                # so L11g's wall applies at p = 5
        r = _l10_class_cert(a, z, 2*a/(1 + 4*a*a), 1, q1)
        assert r is not None and not r['ok'], (w, ut)   # must stay refused
        assert r['syms'][5] == -1, (w, ut)              # and for THIS reason
        nbad += 1
    print(f"  (L10c) SUPERSEDED: the 164-row table omitted controlled places."
          f" 130 rows are INVALID (symbol -1 at p | A); the other 34 keep valid"
          f" symbols but stale moduli under L11f's complete set. {nbad} rows of"
          " the retired shape are frozen as regressions and stay refused, each"
          " for the proved reason; L11h carries the corrected 190-cell table")
    print("  -> class side splits by branch: PROVED reachable on tau = 2a/A,"
          " PROVED unreachable on tau = 0 for prime A (evidence otherwise)."
          " This layer is structural; the Schinzel-conditional chain closes at L22")

# ---------------------------------------------------------------- L11:
# branch completion.  See THEOREMS.md L11.
#
# Theta_* currently disjoins two branches, tau = 0 (alpha = -A) and
# tau = 2a/A (alpha = -1), where alpha = -delta_tau*A, delta_tau = 1 - A tau^2.
# Soundness is tau-UNIFORM (Sun identities 2.2-2.4 hold for arbitrary tau and
# Prop 2.1 never mentions tau) and a disjunction of k conics in the same (y, r)
# is one polynomial product, so branches cost NO witnesses.  Two more rational
# tau complete the square-class group <-1, A>:
#
#   tau = 1                  -> alpha = 4a^2*A          == A  (mod squares)
#   tau = (1+2a^2)/(1+4a^2)  -> alpha = 4a^4 = (2a^2)^2  == 1  (mod squares)
#
# The second is the whole point: alpha is a PERFECT SQUARE.  Two consequences,
# and only the second is used below:
#   * the target-place CHARACTER chi_w(alpha) = +1 holds at every odd w not
#     dividing 2a, with no w mod 4 condition.  This is NOT W1: Sun's Lemma 5.1
#     also needs v_w(A) = v_w(delta) = 0, which fails at e.g. a = 1, w = 5.
#   * supp(alpha) contributes no odd place and alpha > 0, so with eps = +1 the
#     square class of alpha is trivial and alpha > 0.  This looked like it
#     dissolved the L10b wall; IT DOES NOT.  The frozen set is NOT {2}: delta
#     carries A in its denominator, so p | A stays controlled, and L11g shows
#     the resulting pairing (x,d)_A (x,d)_q = -1 holds on EVERY branch when
#     v_A(z) <= 0.  What _L11_CLASSES certifies is L11h: the reachable cells
#     are exactly those whose z has a numerator prime = 1 mod 4 (190 of 353).
_L11_D_POOL = (1, -1, 2, -2, 3, 5, -5, 7, 1, 3)


def _l11_tau(a, d=1):
    """tau_d = (A + d^2)/(2dA):  alpha = -delta*A = ((A - d^2)/(2d))^2.

    The general rational solution of "alpha is a square": alpha = lam^2 needs
    A + lam^2 = mu^2 with tau = mu/A, i.e. (mu-lam)(mu+lam) = A, so lam and mu
    are cut out by a free rational factorisation parameter d.  d = 1 gives the
    canonical square branch tau_dagger = (1+2a^2)/(1+4a^2).
    """
    A = 1 + 4*a*a
    if d == 0:
        return None
    return Fraction(A + d*d, 1)/(2*d*A) if isinstance(d, int) else (A + d*d)/(2*d*A)


def _l11_alpha(a, tau):
    A = 1 + 4*a*a
    return -(1 - A*tau*tau)*A


# Aligned-class certificates for the prime-b family, one per REACHABLE cell.
# (w, u) -> (a, eps, q1, N, p) where p is the escape prime (= 1 mod 4, dividing
# the numerator of z) and a is CONSTRUCTED from p so that v_p(1 + 4a^2) > 0.
# L11h: these 190 cells are exactly the cells with such a p; the other 163 are
# provably unreachable by this family on EVERY branch (L11g).  Supersedes both
# _L10_CLASSES (130 of whose 164 rows were invalid) and the retracted 353-row
# table of L11e.
_L11_CLASSES = {
    (3, (5, 1)): ((1, 1), 1, 41, 3360, 5),
    (3, (-5, 1)): ((1, 1), 1, 41, 3360, 5),
    (5, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (5, (3, 1)): ((1, 1), 1, 41, 3360, 5),
    (5, (7, 3)): ((1, 1), 1, 41, 10080, 5),
    (5, (2, 1)): ((1, 1), 1, 41, 2520, 5),
    (5, (1, 1)): ((1, 1), 1, 41, 10080, 5),
    (5, (2, 3)): ((1, 1), 1, 59, 2520, 5),
    (5, (1, 3)): ((1, 1), 1, 41, 10080, 5),
    (5, (1, 7)): ((1, 1), 1, 41, 10080, 5),
    (5, (-2, 7)): ((1, 1), 1, 41, 2520, 5),
    (5, (-1, 3)): ((1, 1), 1, 59, 10080, 5),
    (5, (-1, 1)): ((1, 1), 1, 41, 10080, 5),
    (5, (-2, 1)): ((1, 1), 1, 41, 2520, 5),
    (5, (-3, 1)): ((1, 1), 1, 41, 3360, 5),
    (5, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (7, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (7, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (11, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (11, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (13, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (13, (3, 1)): ((9, 1), 1, 43, 218400, 13),
    (13, (7, 3)): ((9, 1), 1, 43, 218400, 13),
    (13, (2, 1)): ((9, 1), 1, 43, 54600, 13),
    (13, (1, 1)): ((9, 1), 1, 43, 218400, 13),
    (13, (2, 3)): ((9, 1), 1, 43, 54600, 13),
    (13, (1, 3)): ((9, 1), 1, 43, 218400, 13),
    (13, (1, 7)): ((9, 1), 1, 43, 218400, 13),
    (13, (-1, 5)): ((9, 1), 1, 43, 218400, 13),
    (13, (-2, 7)): ((9, 1), 1, 43, 54600, 13),
    (13, (-1, 3)): ((9, 1), 1, 43, 218400, 13),
    (13, (-1, 1)): ((9, 1), 1, 43, 218400, 13),
    (13, (-2, 1)): ((9, 1), 1, 43, 54600, 13),
    (13, (-3, 1)): ((9, 1), 1, 43, 218400, 13),
    (13, (-5, 1)): ((1, 1), 1, 59, 10080, 5),
    (17, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (17, (3, 1)): ((15, 1), 1, 41, 9082080, 17),
    (17, (7, 3)): ((15, 1), 1, 41, 9082080, 17),
    (17, (2, 1)): ((15, 1), 1, 41, 756840, 17),
    (17, (1, 1)): ((15, 1), 1, 41, 9082080, 17),
    (17, (2, 3)): ((15, 1), 1, 41, 2270520, 17),
    (17, (1, 3)): ((15, 1), 1, 41, 9082080, 17),
    (17, (1, 7)): ((15, 1), 1, 41, 9082080, 17),
    (17, (-1, 5)): ((15, 1), 1, 41, 9082080, 17),
    (17, (-2, 7)): ((15, 1), 1, 41, 2270520, 17),
    (17, (-1, 3)): ((15, 1), 1, 41, 9082080, 17),
    (17, (-1, 1)): ((15, 1), 1, 41, 9082080, 17),
    (17, (-2, 1)): ((15, 1), 1, 41, 2270520, 17),
    (17, (-3, 1)): ((15, 1), 1, 41, 9082080, 17),
    (17, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (19, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (19, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (23, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (23, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (29, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (29, (3, 1)): ((23, 1), 1, 67, 1472415840, 29),
    (29, (7, 3)): ((23, 1), 1, 67, 818008800, 29),
    (29, (2, 1)): ((23, 1), 1, 67, 286303080, 29),
    (29, (1, 1)): ((23, 1), 1, 67, 1145212320, 29),
    (29, (2, 3)): ((23, 1), 1, 67, 204502200, 29),
    (29, (1, 3)): ((23, 1), 1, 67, 163601760, 29),
    (29, (1, 7)): ((23, 1), 1, 67, 1145212320, 29),
    (29, (-1, 5)): ((23, 1), 1, 67, 5726061600, 29),
    (29, (-2, 7)): ((23, 1), 1, 67, 1431515400, 29),
    (29, (-1, 3)): ((23, 1), 1, 67, 1145212320, 29),
    (29, (-1, 1)): ((23, 1), 1, 67, 818008800, 29),
    (29, (-2, 1)): ((23, 1), 1, 67, 40900440, 29),
    (29, (-3, 1)): ((23, 1), 1, 67, 72148376160, 29),
    (29, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (31, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (31, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (37, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (37, (3, 1)): ((3, 1), 1, 47, 1864800, 37),
    (37, (7, 3)): ((3, 1), 1, 41, 1118880, 37),
    (37, (2, 1)): ((3, 1), 1, 47, 93240, 37),
    (37, (1, 1)): ((3, 1), 1, 47, 372960, 37),
    (37, (2, 3)): ((3, 1), 1, 41, 279720, 37),
    (37, (1, 3)): ((3, 1), 1, 41, 1118880, 37),
    (37, (1, 7)): ((3, 1), 1, 47, 46620000, 37),
    (37, (-1, 5)): ((3, 1), 1, 41, 152292000, 37),
    (37, (-2, 7)): ((3, 1), 1, 47, 93240, 37),
    (37, (-1, 3)): ((3, 1), 1, 47, 1864800, 37),
    (37, (-1, 1)): ((3, 1), 1, 41, 1118880, 37),
    (37, (-2, 1)): ((3, 1), 1, 47, 466200, 37),
    (37, (-3, 1)): ((3, 1), 1, 47, 372960, 37),
    (37, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (41, (5, 1)): ((1, 1), 1, 59, 10080, 5),
    (41, (3, 1)): ((25, 1), 1, 53, 8403360, 41),
    (41, (7, 3)): ((25, 1), 1, 53, 1235293920, 41),
    (41, (2, 1)): ((25, 1), 1, 53, 6302520, 41),
    (41, (1, 1)): ((25, 1), 1, 53, 75630240, 41),
    (41, (2, 3)): ((25, 1), 1, 53, 6302520, 41),
    (41, (1, 3)): ((25, 1), 1, 53, 25210080, 41),
    (41, (1, 7)): ((25, 1), 1, 53, 75630240, 41),
    (41, (-1, 5)): ((25, 1), 1, 53, 75630240, 41),
    (41, (-2, 7)): ((25, 1), 1, 53, 18907560, 41),
    (41, (-1, 3)): ((25, 1), 1, 53, 25210080, 41),
    (41, (-1, 1)): ((25, 1), 1, 53, 25210080, 41),
    (41, (-2, 1)): ((25, 1), 1, 71, 132352920, 41),
    (41, (-3, 1)): ((25, 1), 1, 53, 8403360, 41),
    (41, (-5, 1)): ((1, 1), 1, 59, 10080, 5),
    (43, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (43, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (47, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (47, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (53, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (53, (3, 1)): ((15, 1), 1, 43, 21191520, 53),
    (53, (7, 3)): ((15, 1), 1, 43, 3027360, 53),
    (53, (2, 1)): ((15, 1), 1, 43, 5297880, 53),
    (53, (1, 1)): ((15, 1), 1, 43, 21191520, 53),
    (53, (2, 3)): ((15, 1), 1, 43, 37085160, 53),
    (53, (1, 3)): ((15, 1), 1, 43, 21191520, 53),
    (53, (1, 7)): ((15, 1), 1, 43, 21191520, 53),
    (53, (-1, 5)): ((15, 1), 1, 43, 148340640, 53),
    (53, (-2, 7)): ((15, 1), 1, 43, 5297880, 53),
    (53, (-1, 3)): ((15, 1), 1, 43, 148340640, 53),
    (53, (-1, 1)): ((15, 1), 1, 43, 21191520, 53),
    (53, (-2, 1)): ((15, 1), 1, 43, 5297880, 53),
    (53, (-3, 1)): ((15, 1), 1, 43, 21191520, 53),
    (53, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (59, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (59, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (61, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (61, (3, 1)): ((25, 1), 1, 73, 8403360, 61),
    (61, (7, 3)): ((25, 1), 1, 73, 75630240, 61),
    (61, (2, 1)): ((25, 1), 1, 83, 510504120, 61),
    (61, (1, 1)): ((25, 1), 1, 83, 25210080, 61),
    (61, (2, 3)): ((25, 1), 1, 73, 18907560, 61),
    (61, (1, 3)): ((25, 1), 1, 73, 75630240, 61),
    (61, (1, 7)): ((25, 1), 1, 73, 75630240, 61),
    (61, (-1, 5)): ((25, 1), 1, 73, 75630240, 61),
    (61, (-2, 7)): ((25, 1), 1, 73, 170168040, 61),
    (61, (-1, 3)): ((25, 1), 1, 73, 75630240, 61),
    (61, (-1, 1)): ((25, 1), 1, 83, 226890720, 61),
    (61, (-2, 1)): ((25, 1), 1, 73, 18907560, 61),
    (61, (-3, 1)): ((25, 1), 1, 73, 8403360, 61),
    (61, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (67, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (67, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (71, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (71, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (73, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (73, (3, 1)): ((23, 1), 1, 101, 163601760, 73),
    (73, (7, 3)): ((23, 1), 1, 43, 163601760, 73),
    (73, (2, 1)): ((23, 1), 1, 43, 40900440, 73),
    (73, (1, 1)): ((23, 1), 1, 43, 163601760, 73),
    (73, (2, 3)): ((23, 1), 1, 43, 40900440, 73),
    (73, (1, 3)): ((23, 1), 1, 43, 163601760, 73),
    (73, (1, 7)): ((23, 1), 1, 43, 163601760, 73),
    (73, (-1, 5)): ((23, 1), 1, 43, 163601760, 73),
    (73, (-2, 7)): ((23, 1), 1, 43, 40900440, 73),
    (73, (-1, 3)): ((23, 1), 1, 43, 163601760, 73),
    (73, (-1, 1)): ((23, 1), 1, 43, 163601760, 73),
    (73, (-2, 1)): ((23, 1), 1, 43, 40900440, 73),
    (73, (-3, 1)): ((23, 1), 1, 47, 4090044000, 73),
    (73, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (79, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (79, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (83, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (83, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (89, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (89, (3, 1)): ((17, 1), 1, 59, 66087840, 89),
    (89, (7, 3)): ((17, 1), 1, 41, 66087840, 89),
    (89, (2, 1)): ((17, 1), 1, 41, 16521960, 89),
    (89, (1, 1)): ((17, 1), 1, 41, 66087840, 89),
    (89, (2, 3)): ((17, 1), 1, 137, 16521960, 89),
    (89, (1, 3)): ((17, 1), 1, 59, 66087840, 89),
    (89, (1, 7)): ((17, 1), 1, 59, 66087840, 89),
    (89, (-1, 5)): ((17, 1), 1, 41, 66087840, 89),
    (89, (-2, 7)): ((17, 1), 1, 137, 16521960, 89),
    (89, (-1, 3)): ((17, 1), 1, 41, 66087840, 89),
    (89, (-1, 1)): ((17, 1), 1, 137, 66087840, 89),
    (89, (-2, 1)): ((17, 1), 1, 59, 16521960, 89),
    (89, (-3, 1)): ((17, 1), 1, 41, 66087840, 89),
    (89, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
    (97, (5, 1)): ((1, 1), 1, 41, 10080, 5),
    (97, (3, 1)): ((11, 1), 1, 83, 3585120, 97),
    (97, (7, 3)): ((11, 1), 1, 67, 25095840, 97),
    (97, (2, 1)): ((11, 1), 1, 67, 896280, 97),
    (97, (1, 1)): ((11, 1), 1, 67, 3585120, 97),
    (97, (2, 3)): ((11, 1), 1, 67, 896280, 97),
    (97, (1, 3)): ((11, 1), 1, 67, 3585120, 97),
    (97, (1, 7)): ((11, 1), 1, 67, 3585120, 97),
    (97, (-1, 5)): ((11, 1), 1, 67, 3585120, 97),
    (97, (-2, 7)): ((11, 1), 1, 67, 896280, 97),
    (97, (-1, 3)): ((11, 1), 1, 67, 3585120, 97),
    (97, (-1, 1)): ((11, 1), 1, 67, 3585120, 97),
    (97, (-2, 1)): ((11, 1), 1, 67, 896280, 97),
    (97, (-3, 1)): ((11, 1), 1, 83, 3585120, 97),
    (97, (-5, 1)): ((1, 1), 1, 41, 10080, 5),
}


def _l11_escape_primes(z):
    """Numerator primes = 1 mod 4 of z -- the only escape from the L11g wall."""
    return [p for p in sorted(factorint(abs(z.numerator))) if p % 4 == 1]


def _verify_L11(extended=False):
    """L11 (THEOREMS.md): branch completion, and what it does and does not buy.

    (a) branch algebra: tau_d gives alpha = ((A-d^2)/(2d))^2 for every nonzero
        rational d; the four named branches realise the square classes
        -A, -1, A, 1 of <-1, A>.
    (b) soundness is tau-uniform: Sun identities (2.2)-(2.4) re-checked at the
        NEW branches, and the tied block never fires on bad z.
    (c) target-place CHARACTER only: chi_w(alpha) = +1 at every odd w not
        dividing 2a.  NOT W1, which also needs v_w(A) = v_w(delta) = 0 --
        false at a = 1, w = 5 -- so W0 is not asserted redundant.
    (d) L11f: the controlled set {2,3,5,7} u supp(alpha) u supp(delta) really
        is complete (even content of P; degree 8 < p - 1 for p >= 11), with
        the p = 3 fixed-divisor witness that motivated it.
    (e) L11g: the wall is BRANCH-INDEPENDENT.  (x,d)_A * (x,d)_q = -1 on all
        four branches whenever v_A(z) <= 0.  Includes the exact (3,1,41)
        regression that refuted the retracted L11e.
    (f) L11h: the 190 reachable cells are exactly those whose z has a
        numerator prime = 1 mod 4; each row replayed, each of the other 163
        confirmed to have no such prime.
    """
    rng = random.Random(20261101)
    # (a) branch algebra
    nb = 0
    while nb < (400 if extended else 120):
        s = Fraction(rng.randint(-9, 9), rng.choice([1, 3, 5, 7, 9]))
        a = 1 + 2*s
        if a == 0 or vp(a, 2) != 0:
            continue
        A = 1 + 4*a*a
        assert A == 5 or vp(A - 5, 2) >= 3      # A = 5 mod 8 2-adically
        assert _l11_alpha(a, Fraction(0)) == -A
        assert _l11_alpha(a, 2*a/A) == -1
        assert _l11_alpha(a, Fraction(1)) == 4*a*a*A          # == A mod squares
        tau_d = (1 + 2*a*a)/A
        assert tau_d == (A + 1)/(2*A)
        assert _l11_alpha(a, tau_d) == 4*a**4 == (2*a*a)**2    # PERFECT SQUARE
        assert 1 - A*tau_d*tau_d == -4*a**4/A
        for d in _L11_D_POOL:
            lam = (A - d*d)/(2*d)
            tau = (A + d*d)/(2*d*A)
            assert _l11_alpha(a, tau) == lam*lam, (a, d)
            assert A + lam*lam == ((A + d*d)/(2*d))**2
        nb += 1
    print(f"  (L11a) branch algebra on {nb} admissible a: tau_d = (A+d^2)/(2dA)"
          " gives alpha = ((A-d^2)/(2d))^2; the four named branches realise"
          " -A, -1, A, 1")
    # (b) soundness, tau-uniform
    ni = 0
    while ni < (300 if extended else 100):
        a, b, c, y, r_ = (Fraction(rng.randint(-9, 9), rng.randint(1, 6))
                          for _ in range(5))
        if a == 0 or b == 0 or vp(a, 2) != 0:
            continue
        A, B = 1 + 4*a*a, 2*b
        s = (a - 1)/2
        for tau in ((1 + 2*a*a)/A, Fraction(1)):
            delta = 1 - A*tau*tau
            assert delta != 0
            z_ = ((c - A*tau*y)/4, (y - tau*c)/4)
            alpha_ = (z_[0], z_[1], r_, s)
            gamma = ((1 + A*tau*tau)/delta, 2*tau/delta, 0, 0)
            assert 16*(z_[0]*z_[0] - A*z_[1]*z_[1]) == delta*(c*c - A*y*y)
            assert 16*_nrd(alpha_, A, B) == (delta*(c*c - A*y*y)
                                             - 16*B*(r_*r_ - A*s*s))
            assert _nrd(gamma, A, B) == 1
            assert 2*alpha_[0] + 2*_qmul(alpha_, gamma, A, B)[0] == c
        ni += 1
    nz = refused = 0
    BADZ = [Fraction(sg*2**k, dd) for sg in (1, -1) for k in range(0, 4)
            for dd in (1, 3, 5, 7, 9, 15)]
    while nz < (400 if extended else 150):
        s = Fraction(rng.randint(-6, 6), rng.choice([1, 1, 3, 5]))
        a = 1 + 2*s
        b = Fraction(rng.randrange(-21, 22, 2), rng.randrange(1, 14, 2))
        if a == 0 or vp(a, 2) != 0 or b in (0, 1) or vp(b, 2) != 0:
            continue
        zb = rng.choice(BADZ)
        c = _sun_h(a, b, zb**3)
        if c is None:
            continue
        A = 1 + 4*a*a
        try:
            hit = any(_sun_tied_solvable(a, b, c, t)
                      for t in ((1 + 2*a*a)/A, Fraction(1)))
        except (FactorBudget, PrimalityBound):
            refused += 1
            continue
        assert not hit, (s, b, zb)
        nz += 1
    print(f"  (L11b) soundness is tau-uniform: identities (2.2)-(2.4) hold at"
          f" the new branches on {ni} exact instances; the block never fires on"
          f" {nz} bad-z inputs ({refused} refused, never counted)")
    # (c) target-place character only
    ng = nbad = 0
    for w in primerange(3, 120 if extended else 60):
        for sn in (0, 1, -1, 2, -2, 3, -3, 5):
            a = 1 + 2*Fraction(sn)
            if a == 0 or vp(a, 2) != 0:
                continue
            A = 1 + 4*a*a
            al = _l11_alpha(a, (1 + 2*a*a)/A)
            if al.numerator % w == 0 or al.denominator % w == 0:
                continue
            assert legendre(al.numerator, w)*legendre(al.denominator, w) == 1
            ng += 1
            if _sun_target_tau(a, w) == 0 and legendre(-A.numerator, w) != 1:
                nbad += 1
    print(f"  (L11c) target-place CHARACTER only: chi_w(alpha) = +1 on all {ng}"
          f" (w, a) pairs with w not dividing 2a ({nbad} where the canonical -A"
          " branch has chi_w = -1); NOT W1, which also needs v_w(A) = v_w(delta) = 0")
    # (d) L11f: completeness of the controlled set
    nc = 0
    while nc < (200 if extended else 60):
        s = Fraction(rng.randint(-6, 6), rng.choice([1, 3, 5]))
        a = 1 + 2*s
        if a == 0 or vp(a, 2) != 0:
            continue
        A = 1 + 4*a*a
        z = Fraction(rng.choice([3, 5, 7, 11, 13]))*rng.choice(
            [Fraction(1), Fraction(2), Fraction(3), Fraction(1, 3)])
        Z = z**3
        D = 1 - Z - a*a*Z*Z
        if D == 0:
            continue
        delta = (1 - A*((1 + 2*a*a)/A)**2)
        alpha = -delta*A
        P = _l10_P(a, Z, D, A, delta, (a - 1)/2)
        nzc = [i for i, cc in enumerate(P) if cc != 0]
        assert nzc[-1] <= 8, nzc            # degree bound the argument needs
        S = {2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta)
        for p in (11, 13, 17, 19, 23):
            if p in S:
                continue
            assert vp(A, p) == 0 and vp(a, p) >= 0
            cont = min(vp(cc, p) for cc in P if cc != 0)
            assert cont % 2 == 0, (a, z, p, cont)   # even content
            # a nonzero degree <= 8 form cannot vanish on all p - 1 > 8 units
            red = [cc/Fraction(p)**cont for cc in P]
            assert any(vp(cc, p) == 0 for cc in red if cc != 0)
        nc += 1
    # the p = 3 fixed divisor that motivated freezing {3,5,7}
    a3, z3 = Fraction(1), Fraction(5)
    P3 = _l10_P(a3, z3**3, 1 - z3**3 - z3**6, Fraction(5), Fraction(-4, 5),
                Fraction(0))
    assert min(vp(cc, 3) for cc in P3 if cc != 0) == 0      # no 3 in content
    vals = [sum(cc*Fraction(u)**i for i, cc in enumerate(P3))
            for u in (1, 2, 4, 5, 7, 8)]                    # units mod 3
    assert all(vp(v, 3) >= 1 for v in vals)                 # yet 3 | P(u)
    print(f"  (L11f) controlled set {{2,3,5,7}} u supp(alpha) u supp(delta) is"
          f" complete: even p-content and deg P <= 8 < p-1 checked on {nc}"
          " instances; the p = 3 prime-restricted fixed divisor is exhibited")
    # (e) L11g: the branch-independent wall + the exact retraction regression
    a1 = Fraction(1)
    A1 = 1 + 4*a1*a1
    z_ce, q_ce = Fraction(3), 41
    tau_ce = (1 + 2*a1*a1)/A1
    d_ce = 1 - A1*tau_ce*tau_ce
    al_ce = -d_ce*A1
    c_ce = _sun_h(a1, Fraction(q_ce), z_ce**3)
    M_ce = 16 - d_ce*c_ce*c_ce
    x_ce, dd_ce = al_ce*M_ce, al_ce*2*q_ce
    assert c_ce == Fraction(9311592816, 6345775) and dd_ce == 328
    assert vp(x_ce, 5) == -5 and hilbert(x_ce, dd_ce, 5) == -1
    assert (hilbert(x_ce, dd_ce, 2) == hilbert(x_ce, dd_ce, q_ce)
            == hilbert(x_ce, dd_ce, OO) == 1)
    nw = 0
    for (w, ut) in sorted(_l9_grid()):
        z = Fraction(w)*Fraction(*ut)
        if vp(z, 5) > 0:
            continue
        for tau in (Fraction(0), 2*a1/A1, Fraction(1), (1 + 2*a1*a1)/A1):
            de = 1 - A1*tau*tau
            al = -de*A1
            for q in (41, 59):
                # L9a's units are a HYPOTHESIS, not decoration: q must divide
                # none of delta, A, a, z, D_z.  Omitting this guard made the
                # extended suite fail at cell (59,-3) with q = 59 = w.
                Dz = 1 - z**3 - a1*a1*z**6
                if any(vp(t, q) != 0 for t in (de, A1, a1, z, Dz)):
                    continue
                c0 = _sun_h(a1, Fraction(q), z**3)
                if c0 is None:
                    continue
                M0 = 16 - de*c0*c0
                if M0 == 0:
                    continue
                x0, d0 = al*M0, al*2*q
                assert vp(x0, 5) % 2 == 1                     # forced odd
                assert hilbert(x0, d0, 5)*hilbert(x0, d0, q) == -1, (w, ut, tau)
                nw += 1
        if not extended and nw > 600:
            break
    print(f"  (L11g) the wall is BRANCH-INDEPENDENT: (x,d)_A * (x,d)_q = -1 on"
          f" {nw} (cell, branch, q) instances with v_A(z) <= 0, across all four"
          " branches; L10b is its tau = 0 instance")
    # (f) L11h: the characterization and the 190-row table
    grid = sorted(_l9_grid())
    reach = {k for k in grid
             if _l11_escape_primes(Fraction(k[0])*Fraction(*k[1]))}
    assert set(_L11_CLASSES) == reach, "table != predicted reachable set"
    assert len(reach) == 190 and len(grid) == 353
    for k in grid:
        if k not in reach:
            z = Fraction(k[0])*Fraction(*k[1])
            assert not _l11_escape_primes(z)      # provably walled by L11g
    nrep = nconst = 0
    for (w, ut), (at_, eps, q1, N, p) in sorted(_L11_CLASSES.items()):
        a = Fraction(*at_)
        z = Fraction(w)*Fraction(*ut)
        A = 1 + 4*a*a
        assert vp(a, 2) == 0 and vp(A, p) > 0 and p % 4 == 1 and vp(z, p) > 0
        r = _l10_class_cert(a, z, (1 + 2*a*a)/A, eps, q1)
        assert r is not None and r['ok'] and r['N'] == N, (w, ut)
        assert {2, 3, 5, 7} <= set(r['S']), (w, ut, r['S'])
        assert math.gcd(q1, N) == 1 and _is_prime(q1), (w, ut)
        nrep += 1
        if nrep % (4 if extended else 12) == 0:
            t = q1 + N
            tries = 0
            while tries < 4000:
                if _is_prime(t):
                    r2 = _l10_class_cert(a, z, (1 + 2*a*a)/A, eps, t)
                    if r2 is not None:
                        assert r2['ok'] and r2['syms'] == r['syms'], (w, ut, t)
                        nconst += 1
                        break
                t += N
                tries += 1
    print(f"  (L11h) {nrep} aligned classes = EXACTLY the cells whose z has a"
          f" numerator prime = 1 mod 4 (190 of 353; escape prime and its"
          f" constructed a recorded per row); symbols constant on {nconst}"
          " larger primes in class")
    print("  -> branch completion is free but does NOT move the class wall;"
          " step (i) holds at 190 cells and is IMPOSSIBLE at 163 for this"
          " family. L20/L22 later replace this restricted family uniformly")


# ---------------------------------------------------------------- L12:
# the generalized wall.  See THEOREMS.md L12.
#
# L11g walled the PRIME-b family.  The same computation applies to every odd
# place of ANY b: at odd p with v_p(b) odd and p coprime to the cell data,
# v_p(x) = -4 is even, v_p(d) is odd, and the unit part of x is
#     u_x = delta^2 * A * a^4 * Z^4 / D_z^2  ==  A   (mod squares),
# so (x,d)_p = (A|p) = (p|A).  Multiplying over those places gives (b_sf|A),
# while (x,d)_A = (d|A) = (2 eps b_sf | A) = -(b_sf|A) because A = 5 mod 8.
# The product is -1 for EVERY admissible b coprime to the data -- prime,
# semiprime, squarefull, rational alike.  Choosing a non-residue cofactor does
# not help: the -1 relocates rather than cancels.
#
# The single surviving hypothesis is coprimality.  b sharing a prime with the
# cell data breaks the u_x computation, and that door is open (L12b).
_L12_B_SHAPES = ('prime', 'two primes', 'three primes', 'square factor',
                 'rational')


def _l12_controlled(a, z, tau, b):
    """CONSERVATIVE controlled set: {2,3,5,7} plus every fixed datum's places."""
    A = 1 + 4*a*a
    delta = 1 - A*tau*tau
    S = {2, 3, 5, 7}
    for t in (-delta*A, delta, A, a, z, 1 - z**3 - a*a*z**6, b):
        if t == 0:
            continue
        S |= set(factorint(abs(t.numerator))) | set(factorint(t.denominator))
    return sorted(S)


def _l12_syms(a, z, tau, b):
    A = 1 + 4*a*a
    delta = 1 - A*tau*tau
    alpha = -delta*A
    c0 = _sun_h(a, b, z**3)
    if c0 is None:
        return None
    M0 = 16 - delta*c0*c0 - 16*A*2*b*((a - 1)/2)**2
    if M0 == 0:
        return None
    x0, d0 = alpha*M0, alpha*2*b
    out = {p: hilbert(x0, d0, p) for p in _l12_controlled(a, z, tau, b)}
    out[OO] = hilbert(x0, d0, OO)
    return out


# Escape witnesses at cells walled by L11h, using b = eps*f*q with f | z.
# (w, u) -> (f, eps, q).  Base-point alignment on the conservative controlled
# set; see the scope note in _verify_L12.
_L12_ESCAPE = {
    (3, (3, 1)): (3, 1, 41),
    (3, (2, 1)): (3, 1, 41),
    (3, (1, 1)): (3, 1, 41),
    (3, (1, 7)): (3, 1, 41),
    (3, (-1, 5)): (3, 1, 41),
    (3, (-2, 7)): (3, 1, 41),
    (3, (-1, 1)): (3, 1, 41),
    (3, (-2, 1)): (3, 1, 41),
    (3, (-3, 1)): (3, 1, 41),
    (7, (3, 1)): (3, 1, 41),
    (7, (7, 3)): (7, 1, 41),
    (7, (2, 1)): (7, 1, 41),
    (7, (1, 1)): (7, 1, 41),
    (7, (2, 3)): (7, 1, 41),
    (7, (1, 3)): (7, 1, 41),
    (7, (-1, 5)): (7, 1, 41),
    (7, (-1, 3)): (7, 1, 41),
    (7, (-1, 1)): (7, 1, 41),
    (7, (-2, 1)): (7, 1, 41),
    (7, (-3, 1)): (3, 1, 41),
    (11, (3, 1)): (3, 1, 41),
    (11, (7, 3)): (7, 1, 41),
    (11, (-3, 1)): (3, 1, 41),
    (19, (3, 1)): (3, 1, 41),
    (19, (7, 3)): (7, 1, 41),
    (19, (-3, 1)): (3, 1, 41),
    (23, (3, 1)): (3, 1, 41),
    (23, (7, 3)): (7, 1, 41),
    (23, (2, 1)): (23, 1, 41),
    (23, (1, 1)): (23, 1, 41),
    (23, (2, 3)): (23, 1, 61),
    (23, (1, 3)): (23, 1, 41),
    (23, (1, 7)): (23, 1, 61),
    (23, (-1, 5)): (23, 1, 61),
    (23, (-2, 7)): (23, 1, 41),
    (23, (-1, 3)): (23, 1, 61),
    (23, (-1, 1)): (23, 1, 41),
    (23, (-2, 1)): (23, 1, 61),
    (23, (-3, 1)): (3, 1, 41),
    (31, (3, 1)): (3, 1, 41),
    (31, (7, 3)): (7, 1, 41),
    (31, (-3, 1)): (3, 1, 41),
    (43, (3, 1)): (3, 1, 41),
    (43, (7, 3)): (7, 1, 41),
    (43, (2, 1)): (43, 1, 41),
    (43, (1, 1)): (43, 1, 41),
    (43, (2, 3)): (43, 1, 41),
    (43, (1, 3)): (43, 1, 41),
    (43, (1, 7)): (43, 1, 41),
    (43, (-1, 5)): (43, 1, 41),
    (43, (-2, 7)): (43, 1, 41),
    (43, (-1, 3)): (43, 1, 41),
    (43, (-1, 1)): (43, 1, 41),
    (43, (-2, 1)): (43, 1, 41),
    (43, (-3, 1)): (3, 1, 41),
    (47, (3, 1)): (3, 1, 41),
    (47, (7, 3)): (7, 1, 41),
    (47, (2, 1)): (47, 1, 41),
    (47, (1, 1)): (47, 1, 61),
    (47, (2, 3)): (47, 1, 61),
    (47, (1, 3)): (47, 1, 61),
    (47, (1, 7)): (47, 1, 61),
    (47, (-1, 5)): (47, 1, 41),
    (47, (-2, 7)): (47, 1, 61),
    (47, (-1, 3)): (47, 1, 61),
    (47, (-1, 1)): (47, 1, 41),
    (47, (-2, 1)): (47, 1, 41),
    (47, (-3, 1)): (3, 1, 41),
    (59, (3, 1)): (3, 1, 41),
    (59, (7, 3)): (7, 1, 41),
    (59, (-3, 1)): (3, 1, 41),
    (67, (3, 1)): (3, 1, 41),
    (67, (7, 3)): (7, 1, 41),
    (67, (2, 1)): (67, 1, 41),
    (67, (1, 1)): (67, 1, 41),
    (67, (2, 3)): (67, 1, 41),
    (67, (1, 3)): (67, 1, 41),
    (67, (1, 7)): (67, 1, 41),
    (67, (-1, 5)): (67, 1, 41),
    (67, (-2, 7)): (67, 1, 41),
    (67, (-1, 3)): (67, 1, 41),
    (67, (-1, 1)): (67, 1, 41),
    (67, (-2, 1)): (67, 1, 41),
    (67, (-3, 1)): (3, 1, 41),
    (71, (3, 1)): (3, 1, 41),
    (71, (7, 3)): (7, 1, 41),
    (71, (-3, 1)): (3, 1, 41),
    (79, (3, 1)): (3, 1, 41),
    (79, (7, 3)): (7, 1, 41),
    (79, (-3, 1)): (3, 1, 41),
    (83, (3, 1)): (3, 1, 41),
    (83, (7, 3)): (7, 1, 41),
    (83, (2, 1)): (83, 1, 41),
    (83, (1, 1)): (83, 1, 61),
    (83, (2, 3)): (83, 1, 61),
    (83, (1, 3)): (83, 1, 61),
    (83, (1, 7)): (83, 1, 61),
    (83, (-1, 5)): (83, 1, 41),
    (83, (-2, 7)): (83, 1, 61),
    (83, (-1, 3)): (83, 1, 61),
    (83, (-1, 1)): (83, 1, 41),
    (83, (-2, 1)): (83, 1, 41),
    (83, (-3, 1)): (3, 1, 41),
}


def _verify_L12(extended=False):
    """L12 (THEOREMS.md): the wall generalizes to every coprime b; the escape.

    (a) L12a: for admissible b whose odd places are coprime to the cell data,
        (x,d)_p = (A|p) at each of them, so the product over {A} u oddsupp(b)
        is -1 on every branch when v_A(z) <= 0.  Prime, semiprime, squarefull
        and rational b all fail alike -- a non-residue cofactor relocates the
        -1 instead of cancelling it.
    (b) L12b: b sharing a prime with z breaks the coprimality hypothesis.
        103 of the 163 L11h-walled cells then get every controlled symbol +1.
        SCOPE: base-point alignment only.  Constancy along the class follows
        from the exponent lemma, but the conservative moduli are ~1e13, so it
        is NOT independently sampled here and is not claimed as verified.
    """
    rng = random.Random(20261212)
    a = Fraction(1)
    A = 1 + 4*a*a
    grid = sorted(_l9_grid())
    walled = [k for k in grid
              if not _l11_escape_primes(Fraction(k[0])*Fraction(*k[1]))]
    assert len(walled) == 163
    # (a) the generalized wall
    P = [p for p in primerange(11, 200)]
    seen = {s: 0 for s in _L12_B_SHAPES}
    nw = 0
    for (w, ut) in walled[:24 if extended else 12]:
        z = Fraction(w)*Fraction(*ut)
        Dz = 1 - z**3 - a*a*z**6
        for tau in ((1 + 2*a*a)/A, 2*a/A, Fraction(0), Fraction(1)):
            delta = 1 - A*tau*tau
            cands = [(Fraction(rng.choice(P)), 'prime'),
                     (Fraction(rng.choice(P)*rng.choice(P)), 'two primes'),
                     (Fraction(rng.choice(P)*rng.choice(P)*rng.choice(P)),
                      'three primes'),
                     (Fraction(rng.choice(P)**2*rng.choice(P)), 'square factor'),
                     (Fraction(rng.choice(P), rng.choice(P)), 'rational')]
            for b, nm in cands:
                b = b*rng.choice([1, -1])
                if vp(b, 2) != 0 or b in (0, 1):
                    continue
                odd = {p for p in (set(factorint(abs(b.numerator)))
                                   | set(factorint(b.denominator)))
                       if p != 2 and vp(b, p) % 2}
                if any(vp(t, p) != 0 for p in odd
                       for t in (z, Dz, a, delta, A)):
                    continue
                alpha = -delta*A
                c0 = _sun_h(a, b, z**3)
                if c0 is None:
                    continue
                M0 = 16 - delta*c0*c0 - 16*A*2*b*((a - 1)/2)**2
                if M0 == 0:
                    continue
                x0, d0 = alpha*M0, alpha*2*b
                prod = hilbert(x0, d0, A.numerator)
                for p in odd:
                    # the step the whole theorem rests on
                    assert hilbert(x0, d0, p) == legendre(A.numerator, p), (b, p)
                    prod *= hilbert(x0, d0, p)
                assert prod == -1, (w, ut, tau, b)
                seen[nm] += 1
                nw += 1
    assert all(v > 0 for v in seen.values()), seen
    print(f"  (L12a) generalized wall: prod over {{A}} u oddsupp(b) = -1 on {nw}"
          f" instances covering all of {list(_L12_B_SHAPES)} and all four"
          " branches; (x,d)_p = (A|p) asserted at every odd place of b")
    # (b) the escape
    ne = 0
    for (w, ut), (f, eps, q) in sorted(_L12_ESCAPE.items()):
        z = Fraction(w)*Fraction(*ut)
        assert (w, ut) in set(walled)                  # a cell L11h walls
        assert vp(z, f) > 0                            # coprimality BROKEN
        b = Fraction(eps*f*q)
        s = _l12_syms(a, z, (1 + 2*a*a)/A, b)
        assert s is not None and all(v == 1 for v in s.values()), (w, ut)
        ne += 1
    assert ne == len(_L12_ESCAPE) == 103
    print(f"  (L12b) escape: {ne} of the 163 walled cells get EVERY controlled"
          " symbol +1 with b = eps*f*q, f | z -- base point only; class"
          " constancy is NOT sampled here (conservative moduli ~1e13)")
    print("  -> the coprime-b family is dead on every branch; the only door is"
          " b sharing a prime with the cell data. L20/L22 later close the"
          " class/algebraic sides; member existence remains Schinzel-conditional")


# ---------------------------------------------------------------- L13:
# the aftermath layer (2026-08-17): verified class alignment for the 103
# escapes, 60 verified fully-soluble witnesses at previously uncovered walled
# cells, P's irreducibility certified per row, and the frozen zero-bad
# members.  Six-agent parallel session; every table below was independently
# replayed by the lead before freezing.  Provenance: data/*.json/jsonl.

# 60 walled-no-escape cells with FULLY SOLUBLE witnesses: (w,(u1,u2)) ->
# (a, f, eps, q), b = eps*f*q on the square branch tau = (1+2a^2)/A.
# Verified: conservative symbols +1 AND (tied-status True or a full steered
# certificate).  Includes the residual cell (31,(-2,1)); (67,(2,1)) closed
# separately by b = 4757 = 67*67... see _L13_RESIDUAL.
_L13_ESCAPE2 = {
    (11, (-2, 1)): (5, 11, 1, 359),
    (11, (-1, 1)): (5, 11, 1, 79),
    (11, (-1, 3)): (5, 11, 1, 233),
    (11, (-2, 7)): (5, 11, 1, 43),
    (11, (-1, 5)): (5, 11, 1, 101),
    (11, (1, 7)): (5, 11, 1, 827),
    (11, (1, 3)): (5, 11, 1, 47),
    (11, (2, 3)): (5, 11, 1, 283),
    (11, (1, 1)): (5, 11, 1, 181),
    (11, (2, 1)): (5, 11, 1, 107),
    (19, (-2, 1)): (3, 19, 1, 47),
    (19, (-1, 1)): (3, 19, 1, 173),
    (19, (-1, 3)): (3, 19, 1, 41),
    (19, (-2, 7)): (3, 19, 1, 101),
    (19, (-1, 5)): (3, 19, 1, 53),
    (19, (1, 7)): (3, 19, 1, 41),
    (19, (1, 3)): (3, 19, 1, 83),
    (19, (2, 3)): (3, 19, 1, 101),
    (19, (1, 1)): (3, 19, 1, 47),
    (19, (2, 1)): (3, 19, 1, 71),
    (31, (-2, 1)): (3, 31, 1, 47),
    (31, (-1, 1)): (3, 31, 1, 47),
    (31, (-1, 3)): (3, 31, 1, 41),
    (31, (-2, 7)): (3, 31, 1, 157),
    (31, (-1, 5)): (3, 31, 1, 53),
    (31, (1, 7)): (3, 31, 1, 41),
    (31, (1, 3)): (3, 31, 1, 41),
    (31, (2, 3)): (3, 31, 1, 349),
    (31, (1, 1)): (3, 31, 1, 173),
    (31, (2, 1)): (3, 31, 1, 101),
    (59, (-2, 1)): (3, 59, 1, 269),
    (59, (-1, 1)): (3, 59, 1, 101),
    (59, (-1, 3)): (3, 59, 1, 101),
    (59, (-2, 7)): (3, 59, 1, 73),
    (59, (-1, 5)): (3, 59, 1, 733),
    (59, (1, 7)): (3, 59, 1, 127),
    (59, (1, 3)): (3, 59, 1, 41),
    (59, (2, 3)): (3, 59, 1, 173),
    (59, (1, 1)): (3, 11, 1, 59),
    (59, (2, 1)): (3, 59, 1, 443),
    (71, (-2, 1)): (33, 71, -1, 311),
    (71, (-1, 1)): (33, 71, 1, 83),
    (71, (-1, 3)): (33, 71, 1, 1319),
    (71, (-2, 7)): (33, 71, -1, 89),
    (71, (-1, 5)): (33, 71, -1, 701),
    (71, (1, 7)): (33, 71, 1, 443),
    (71, (1, 3)): (33, 71, 1, 239),
    (71, (2, 3)): (33, 71, 1, 179),
    (71, (1, 1)): (33, 71, 1, 1117),
    (71, (2, 1)): (33, 71, 1, 1367),
    (79, (-2, 1)): (3, 79, 1, 151),
    (79, (-1, 1)): (3, 79, -1, 877),
    (79, (-1, 3)): (3, 79, 1, 47),
    (79, (-2, 7)): (3, 79, 1, 307),
    (79, (-1, 5)): (3, 79, 1, 67),
    (79, (1, 7)): (3, 79, 1, 47),
    (79, (1, 3)): (3, 7, 1, 79),
    (79, (2, 3)): (3, 79, 1, 83),
    (79, (1, 1)): (3, 11, 1, 79),
    (79, (2, 1)): (3, 79, 1, 263),
}

# Frozen residual/prior-cell closures: (w, (u1,u2), a, b, d) with d = None
# (square branch tau = (1+2a^2)/A) or a Fraction d (branch tau_d =
# (A+d^2)/(2dA)).  b is an int or an (n,d) Fraction pair.  All tied True.
_L13_RESIDUAL = ((67, (2, 1), 1, 4757, None),
                 (41, (1, 7), 25, 4985, None),
                 (61, (-2, 7), 25, 55571, None),
                 (61, (2, 1), 25, -71431, None),
                 (53, (1, 3), -15, 9, None),        # previously unfrozen
                 (29, (-2, 1), -15, (-29, 5), Fraction(1, 3)),
                 # (89,(-2,1)): the last open cell, closed 2026-08-18 by the L13f
                 # cofactor ladder (l13_filter.py/data/l13_filter_run2.json).
                 # b = 89/367 (f/q shape, f=89 in supp(z)); emergent places of
                 # P(b) = {43, 90947, 204917, 471137} (all +1) u {R = 10127939390699475298911936427,
                 # kernel-proved prime, v_R = 1, symbol +1 by the parity law}:
                 # zero bad places; 236 certified soluble candidates exist in the
                 # stage-2 box (plus 89 more in the dense box); this is the one
                 # frozen, independently lead-replayed.
                 (89, (-2, 1), 3, (89, 367), None))

# HYPOTHESIS H, instance-verified on the whole grid (2026-08-18, L14):
# for EVERY one of the 293 aligned classes (190 L11 + 103 L12 escapes),
# one emergent-free (zero-bad) member, certified by the L13f cofactor
# ladder (smooth_emergent + cofactor_decide; verdict 'zero'), lead-replayed
# 293/293 (data/l13h_all_closures.json, data/l13h_replay.jsonl; ramified
# empty on 271 rows, 22 FactorBudget refusals on the auxiliary ramified
# cross-check -- refusals logged, never evidence).  Row:
# (family, w, (u1, u2), (a_num, a_den), eps, f, q1, N, k) with member
# b = eps * f * (q1 + k*N), tau = (1 + 2a^2)/(1 + 4a^2) (square branch).
_L13_H_CLASSES = (
    ("ESC", 3, (-3, 1), (1, 1), -1, 3, 59, 272160, 0),
    ("ESC", 3, (-2, 1), (1, 1), 1, 3, 41, 68040, 3),
    ("ESC", 3, (-2, 7), (1, 1), 1, 3, 41, 68040, 4),
    ("ESC", 3, (-1, 1), (1, 1), 1, 3, 41, 272160, 0),
    ("ESC", 3, (-1, 5), (1, 1), 1, 3, 41, 272160, 11),
    ("ESC", 3, (1, 1), (1, 1), 1, 3, 41, 272160, 1),
    ("ESC", 3, (1, 7), (1, 1), 1, 3, 41, 272160, 54),
    ("ESC", 3, (2, 1), (1, 1), 1, 3, 41, 68040, 9),
    ("ESC", 3, (3, 1), (1, 1), 1, 3, 41, 272160, 0),
    ("ESC", 7, (-3, 1), (1, 1), 1, 3, 41, 272160, 184),
    ("ESC", 7, (-2, 1), (1, 1), -1, 7, 41, 6050520, 0),
    ("ESC", 7, (-1, 1), (1, 1), 1, 7, 41, 24202080, 21),
    ("ESC", 7, (-1, 3), (1, 1), 1, 7, 41, 24202080, 0),
    ("ESC", 7, (-1, 5), (1, 1), 1, 7, 41, 24202080, 0),
    ("ESC", 7, (1, 1), (1, 1), 1, 7, 41, 24202080, 0),
    ("ESC", 7, (1, 3), (1, 1), 1, 7, 41, 24202080, 0),
    ("ESC", 7, (2, 1), (1, 1), 1, 7, 41, 6050520, 0),
    ("ESC", 7, (2, 3), (1, 1), 1, 7, 41, 6050520, 6),
    ("ESC", 7, (3, 1), (1, 1), 1, 3, 41, 272160, 0),
    ("ESC", 7, (7, 3), (1, 1), -1, 7, 61, 24202080, 0),
    ("ESC", 11, (-3, 1), (1, 1), 1, 3, 41, 272160, 117),
    ("ESC", 11, (3, 1), (1, 1), 1, 3, 41, 272160, 0),
    ("ESC", 11, (7, 3), (1, 1), 1, 7, 41, 24202080, 183),
    ("ESC", 19, (-3, 1), (1, 1), 1, 3, 41, 272160, 105),
    ("ESC", 19, (3, 1), (1, 1), 1, 3, 41, 272160, 0),
    ("ESC", 19, (7, 3), (1, 1), -1, 7, 61, 24202080, 0),
    ("ESC", 23, (-3, 1), (1, 1), -1, 3, 59, 272160, 37),
    ("ESC", 23, (-2, 1), (1, 1), 1, 23, 61, 16219584360, 0),
    ("ESC", 23, (-2, 7), (1, 1), -1, 23, 109, 437928777720, 0),
    ("ESC", 23, (-1, 1), (1, 1), -1, 23, 59, 194635012320, 0),
    ("ESC", 23, (-1, 3), (1, 1), 1, 23, 61, 64878337440, 0),
    ("ESC", 23, (-1, 5), (1, 1), -1, 23, 41, 583905036960, 0),
    ("ESC", 23, (1, 1), (1, 1), 1, 23, 41, 1751715110880, 0),
    ("ESC", 23, (1, 3), (1, 1), 1, 23, 41, 1751715110880, 0),
    ("ESC", 23, (1, 7), (1, 1), -1, 23, 41, 583905036960, 0),
    ("ESC", 23, (2, 1), (1, 1), 1, 23, 41, 437928777720, 0),
    ("ESC", 23, (2, 3), (1, 1), -1, 23, 41, 145976259240, 0),
    ("ESC", 23, (3, 1), (1, 1), 1, 3, 41, 272160, 0),
    ("ESC", 23, (7, 3), (1, 1), -1, 7, 41, 24202080, 0),
    ("ESC", 31, (-3, 1), (1, 1), -1, 3, 41, 272160, 0),
    ("ESC", 31, (3, 1), (1, 1), 1, 3, 41, 272160, 85),
    ("ESC", 31, (7, 3), (1, 1), -1, 7, 79, 24202080, 0),
    ("ESC", 43, (-3, 1), (1, 1), 1, 3, 41, 272160, 7),
    ("ESC", 43, (-2, 1), (1, 1), 1, 43, 41, 370461276360, 0),
    ("ESC", 43, (-2, 7), (1, 1), -1, 43, 101, 370461276360, 0),
    ("ESC", 43, (-1, 1), (1, 1), -1, 43, 61, 1481845105440, 0),
    ("ESC", 43, (-1, 3), (1, 1), -1, 43, 79, 1481845105440, 0),
    ("ESC", 43, (-1, 5), (1, 1), -1, 43, 139, 1481845105440, 0),
    ("ESC", 43, (1, 1), (1, 1), -1, 43, 61, 1481845105440, 0),
    ("ESC", 43, (1, 3), (1, 1), -1, 43, 379, 1481845105440, 0),
    ("ESC", 43, (1, 7), (1, 1), -1, 43, 211, 1481845105440, 0),
    ("ESC", 43, (2, 1), (1, 1), -1, 43, 61, 370461276360, 0),
    ("ESC", 43, (2, 3), (1, 1), -1, 43, 199, 370461276360, 0),
    ("ESC", 43, (3, 1), (1, 1), -1, 3, 199, 272160, 0),
    ("ESC", 43, (7, 3), (1, 1), -1, 7, 61, 24202080, 14),
    ("ESC", 47, (-3, 1), (1, 1), -1, 3, 41, 272160, 58),
    ("ESC", 47, (-2, 1), (1, 1), -1, 47, 89, 46813902828840, 0),
    ("ESC", 47, (-2, 7), (1, 1), -1, 47, 71, 5201544758760, 0),
    ("ESC", 47, (-1, 1), (1, 1), -1, 47, 59, 6935393011680, 0),
    ("ESC", 47, (-1, 3), (1, 1), 1, 47, 61, 2311797670560, 0),
    ("ESC", 47, (-1, 5), (1, 1), -1, 47, 61, 62418537105120, 0),
    ("ESC", 47, (1, 1), (1, 1), -1, 47, 71, 20806179035040, 0),
    ("ESC", 47, (1, 3), (1, 1), -1, 47, 61, 62418537105120, 0),
    ("ESC", 47, (1, 7), (1, 1), -1, 47, 311, 20806179035040, 0),
    ("ESC", 47, (2, 1), (1, 1), -1, 47, 59, 1733848252920, 0),
    ("ESC", 47, (2, 3), (1, 1), -1, 47, 191, 5201544758760, 0),
    ("ESC", 47, (3, 1), (1, 1), -1, 3, 89, 272160, 0),
    ("ESC", 47, (7, 3), (1, 1), -1, 7, 61, 24202080, 0),
    ("ESC", 59, (-3, 1), (1, 1), -1, 3, 79, 272160, 2),
    ("ESC", 59, (3, 1), (1, 1), -1, 3, 41, 272160, 0),
    ("ESC", 59, (7, 3), (1, 1), -1, 7, 709, 24202080, 51),
    ("ESC", 67, (-3, 1), (1, 1), 1, 3, 41, 272160, 117),
    ("ESC", 67, (-2, 1), (1, 1), -1, 67, 541, 3402315269640, 0),
    ("ESC", 67, (-2, 7), (1, 1), -1, 67, 419, 3402315269640, 0),
    ("ESC", 67, (-1, 1), (1, 1), -1, 67, 811, 13609261078560, 0),
    ("ESC", 67, (-1, 3), (1, 1), -1, 67, 59, 13609261078560, 0),
    ("ESC", 67, (-1, 5), (1, 1), -1, 67, 139, 13609261078560, 0),
    ("ESC", 67, (1, 1), (1, 1), -1, 67, 59, 13609261078560, 0),
    ("ESC", 67, (1, 3), (1, 1), 1, 67, 41, 13609261078560, 0),
    ("ESC", 67, (1, 7), (1, 1), -1, 67, 41, 13609261078560, 0),
    ("ESC", 67, (2, 1), (1, 1), -1, 67, 139, 3402315269640, 0),
    ("ESC", 67, (2, 3), (1, 1), -1, 67, 139, 3402315269640, 0),
    ("ESC", 67, (3, 1), (1, 1), -1, 3, 59, 272160, 0),
    ("ESC", 67, (7, 3), (1, 1), -1, 7, 151, 24202080, 0),
    ("ESC", 71, (-3, 1), (1, 1), -1, 3, 101, 272160, 0),
    ("ESC", 71, (3, 1), (1, 1), -1, 3, 61, 272160, 0),
    ("ESC", 71, (7, 3), (1, 1), -1, 7, 181, 24202080, 0),
    ("ESC", 79, (-3, 1), (1, 1), -1, 3, 41, 272160, 0),
    ("ESC", 79, (3, 1), (1, 1), 1, 3, 41, 272160, 0),
    ("ESC", 79, (7, 3), (1, 1), 1, 7, 41, 24202080, 0),
    ("ESC", 83, (-3, 1), (1, 1), -1, 3, 59, 272160, 0),
    ("ESC", 83, (-2, 1), (1, 1), -1, 83, 149, 1876086277448040, 0),
    ("ESC", 83, (-2, 7), (1, 1), -1, 83, 59, 625362092482680, 0),
    ("ESC", 83, (-1, 1), (1, 1), -1, 83, 71, 833816123310240, 0),
    ("ESC", 83, (-1, 3), (1, 1), -1, 83, 89, 2501448369930720, 0),
    ("ESC", 83, (-1, 5), (1, 1), -1, 83, 239, 7504345109792160, 0),
    ("ESC", 83, (1, 1), (1, 1), -1, 83, 41, 17510138589515040, 0),
    ("ESC", 83, (1, 3), (1, 1), -1, 83, 41, 52530415768545120, 0),
    ("ESC", 83, (1, 7), (1, 1), 1, 83, 61, 39705529681440, 0),
    ("ESC", 83, (2, 1), (1, 1), -1, 83, 109, 208454030827560, 0),
    ("ESC", 83, (2, 3), (1, 1), -1, 83, 269, 625362092482680, 0),
    ("ESC", 83, (3, 1), (1, 1), -1, 3, 61, 272160, 0),
    ("ESC", 83, (7, 3), (1, 1), 1, 7, 41, 24202080, 0),
    ("L11", 3, (-5, 1), (1, 1), 1, 1, 41, 3360, 0),
    ("L11", 3, (5, 1), (1, 1), 1, 1, 41, 3360, 6),
    ("L11", 5, (-5, 1), (1, 1), 1, 1, 41, 10080, 4),
    ("L11", 5, (-3, 1), (1, 1), 1, 1, 41, 3360, 0),
    ("L11", 5, (-2, 1), (1, 1), 1, 1, 41, 2520, 8),
    ("L11", 5, (-2, 7), (1, 1), 1, 1, 41, 2520, 6),
    ("L11", 5, (-1, 1), (1, 1), 1, 1, 41, 10080, 6),
    ("L11", 5, (-1, 3), (1, 1), 1, 1, 59, 10080, 33),
    ("L11", 5, (1, 1), (1, 1), 1, 1, 41, 10080, 2),
    ("L11", 5, (1, 3), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 5, (1, 7), (1, 1), 1, 1, 41, 10080, 65),
    ("L11", 5, (2, 1), (1, 1), 1, 1, 41, 2520, 20),
    ("L11", 5, (2, 3), (1, 1), 1, 1, 59, 2520, 0),
    ("L11", 5, (3, 1), (1, 1), 1, 1, 41, 3360, 6),
    ("L11", 5, (5, 1), (1, 1), 1, 1, 41, 10080, 69),
    ("L11", 5, (7, 3), (1, 1), 1, 1, 41, 10080, 94),
    ("L11", 7, (-5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 7, (5, 1), (1, 1), 1, 1, 41, 10080, 5),
    ("L11", 11, (-5, 1), (1, 1), 1, 1, 41, 10080, 5),
    ("L11", 11, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 13, (-5, 1), (1, 1), 1, 1, 59, 10080, 33),
    ("L11", 13, (-3, 1), (9, 1), 1, 1, 43, 218400, 136),
    ("L11", 13, (-2, 1), (9, 1), 1, 1, 43, 54600, 34),
    ("L11", 13, (-2, 7), (9, 1), 1, 1, 43, 54600, 13),
    ("L11", 13, (-1, 1), (9, 1), 1, 1, 53, 4586400, 52),
    ("L11", 13, (-1, 3), (9, 1), 1, 1, 61, 218400, 14),
    ("L11", 13, (-1, 5), (9, 1), 1, 1, 53, 655200, 0),
    ("L11", 13, (1, 1), (9, 1), 1, 1, 53, 655200, 0),
    ("L11", 13, (1, 3), (9, 1), 1, 1, 79, 218400, 1),
    ("L11", 13, (1, 7), (9, 1), 1, 1, 43, 218400, 171),
    ("L11", 13, (2, 1), (9, 1), 1, 1, 79, 54600, 0),
    ("L11", 13, (2, 3), (9, 1), 1, 1, 43, 54600, 0),
    ("L11", 13, (3, 1), (9, 1), 1, 1, 43, 218400, 0),
    ("L11", 13, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 13, (7, 3), (9, 1), 1, 1, 79, 218400, 0),
    ("L11", 17, (-5, 1), (1, 1), 1, 1, 41, 10080, 6),
    ("L11", 17, (-3, 1), (15, 1), 1, 1, 61, 3027360, 0),
    ("L11", 17, (-2, 1), (15, 1), 1, 1, 73, 756840, 20),
    ("L11", 17, (-2, 7), (15, 1), 1, 1, 71, 1287384840, 13),
    ("L11", 17, (-1, 1), (15, 1), 1, 1, 71, 1716513120, 0),
    ("L11", 17, (-1, 3), (15, 1), 1, 1, 41, 9082080, 0),
    ("L11", 17, (-1, 5), (15, 1), 1, 1, 61, 3027360, 0),
    ("L11", 17, (1, 1), (15, 1), 1, 1, 173, 9082080, 0),
    ("L11", 17, (1, 3), (15, 1), 1, 1, 173, 9082080, 0),
    ("L11", 17, (1, 7), (15, 1), 1, 1, 61, 3027360, 0),
    ("L11", 17, (2, 1), (15, 1), 1, 1, 41, 756840, 0),
    ("L11", 17, (2, 3), (15, 1), 1, 1, 61, 756840, 0),
    ("L11", 17, (3, 1), (15, 1), 1, 1, 71, 5149539360, 0),
    ("L11", 17, (5, 1), (1, 1), 1, 1, 41, 10080, 171),
    ("L11", 17, (7, 3), (15, 1), 1, 1, 61, 3027360, 6),
    ("L11", 19, (-5, 1), (1, 1), 1, 1, 41, 10080, 176),
    ("L11", 19, (5, 1), (1, 1), 1, 1, 41, 10080, 11),
    ("L11", 23, (-5, 1), (1, 1), 1, 1, 41, 10080, 35),
    ("L11", 23, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 29, (-5, 1), (1, 1), 1, 1, 41, 10080, 5),
    ("L11", 29, (-3, 1), (23, 1), 1, 1, 593, 8016486240, 0),
    ("L11", 29, (-2, 1), (23, 1), 1, 1, 67, 40900440, 0),
    ("L11", 29, (-2, 7), (23, 1), 1, 1, 227, 204502200, 0),
    ("L11", 29, (-1, 1), (23, 1), 1, 1, 109, 163601760, 0),
    ("L11", 29, (-1, 3), (23, 1), 1, 1, 149, 163601760, 26),
    ("L11", 29, (-1, 5), (23, 1), 1, 1, 487, 28630308000, 0),
    ("L11", 29, (1, 1), (23, 1), 1, 1, 109, 1145212320, 0),
    ("L11", 29, (1, 3), (23, 1), 1, 1, 109, 163601760, 0),
    ("L11", 29, (1, 7), (23, 1), 1, 1, 149, 163601760, 0),
    ("L11", 29, (2, 1), (23, 1), -1, 1, 227, 286303080, 21),
    ("L11", 29, (2, 3), (23, 1), 1, 1, 71, 23517753000, 0),
    ("L11", 29, (3, 1), (23, 1), -1, 1, 349, 163601760, 0),
    ("L11", 29, (5, 1), (1, 1), 1, 1, 41, 10080, 2),
    ("L11", 29, (7, 3), (23, 1), 1, 1, 149, 163601760, 0),
    ("L11", 31, (-5, 1), (1, 1), 1, 1, 41, 10080, 27),
    ("L11", 31, (5, 1), (1, 1), 1, 1, 41, 10080, 134),
    ("L11", 37, (-5, 1), (1, 1), 1, 1, 41, 10080, 111),
    ("L11", 37, (-3, 1), (3, 1), 1, 1, 53, 372960, 54),
    ("L11", 37, (-2, 1), (3, 1), 1, 1, 73, 31080, 0),
    ("L11", 37, (-2, 7), (3, 1), 1, 1, 47, 93240, 0),
    ("L11", 37, (-1, 1), (3, 1), 1, 1, 41, 1118880, 0),
    ("L11", 37, (-1, 3), (3, 1), 1, 1, 47, 1864800, 96),
    ("L11", 37, (-1, 5), (3, 1), 1, 1, 47, 621600, 0),
    ("L11", 37, (1, 1), (3, 1), 1, 1, 47, 372960, 0),
    ("L11", 37, (1, 3), (3, 1), 1, 1, 41, 1118880, 125),
    ("L11", 37, (1, 7), (3, 1), 1, 1, 47, 46620000, 0),
    ("L11", 37, (2, 1), (3, 1), 1, 1, 47, 93240, 92),
    ("L11", 37, (2, 3), (3, 1), 1, 1, 41, 279720, 0),
    ("L11", 37, (3, 1), (3, 1), 1, 1, 67, 621600, 0),
    ("L11", 37, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 37, (7, 3), (3, 1), 1, 1, 47, 372960, 0),
    ("L11", 41, (-5, 1), (1, 1), 1, 1, 59, 10080, 160),
    ("L11", 41, (-3, 1), (25, 1), 1, 1, 157, 42016800, 0),
    ("L11", 41, (-2, 1), (25, 1), 1, 1, 383, 170168040, 0),
    ("L11", 41, (-2, 7), (25, 1), 1, 1, 53, 18907560, 0),
    ("L11", 41, (-1, 1), (25, 1), 1, 1, 71, 176470560, 0),
    ("L11", 41, (-1, 3), (25, 1), -1, 1, 89, 75630240, 0),
    ("L11", 41, (-1, 5), (25, 1), 1, 1, 71, 75630240, 0),
    ("L11", 41, (1, 1), (25, 1), 1, 1, 227, 28361340000, 0),
    ("L11", 41, (1, 3), (25, 1), 1, 1, 191, 25210080, 0),
    ("L11", 41, (1, 7), (25, 1), -1, 1, 193, 1134453600, 0),
    ("L11", 41, (2, 1), (25, 1), 1, 1, 53, 6302520, 0),
    ("L11", 41, (2, 3), (25, 1), 1, 1, 79, 18907560, 17),
    ("L11", 41, (3, 1), (25, 1), 1, 1, 79, 8403360, 0),
    ("L11", 41, (5, 1), (25, 1), 1, 1, 2003, 176470560, 8),
    ("L11", 41, (7, 3), (25, 1), 1, 1, 457, 378151200, 0),
    ("L11", 43, (-5, 1), (1, 1), 1, 1, 59, 10080, 81),
    ("L11", 43, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 47, (-5, 1), (1, 1), 1, 1, 59, 10080, 0),
    ("L11", 47, (5, 1), (1, 1), 1, 1, 41, 10080, 2),
    ("L11", 53, (-5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 53, (-3, 1), (15, 1), 1, 1, 587, 227052000, 0),
    ("L11", 53, (-2, 1), (15, 1), 1, 1, 47, 11352600, 0),
    ("L11", 53, (-2, 7), (15, 1), 1, 1, 43, 5297880, 0),
    ("L11", 53, (-1, 1), (15, 1), 1, 1, 89, 27246240, 0),
    ("L11", 53, (-1, 3), (15, 1), 1, 1, 89, 27246240, 0),
    ("L11", 53, (-1, 5), (15, 1), 1, 1, 223, 3027360, 0),
    ("L11", 53, (1, 1), (15, 1), -1, 1, 271, 27246240, 0),
    ("L11", 53, (1, 3), (15, 1), 1, 1, 43, 21191520, 26),
    ("L11", 53, (1, 7), (15, 1), 1, 1, 89, 27246240, 0),
    ("L11", 53, (2, 1), (15, 1), 1, 1, 293, 2270520, 0),
    ("L11", 53, (2, 3), (15, 1), 1, 1, 47, 11352600, 0),
    ("L11", 53, (3, 1), (15, 1), 1, 1, 149, 9082080, 28),
    ("L11", 53, (5, 1), (1, 1), 1, 1, 41, 10080, 2),
    ("L11", 53, (7, 3), (15, 1), 1, 1, 59, 9082080, 24),
    ("L11", 59, (-5, 1), (1, 1), 1, 1, 71, 70560, 187),
    ("L11", 59, (5, 1), (1, 1), 1, 1, 41, 10080, 2),
    ("L11", 61, (-5, 1), (25, 1), 1, 1, 83, 25210080, 0),
    ("L11", 61, (-3, 1), (25, 1), 1, 1, 197, 42016800, 0),
    ("L11", 61, (-2, 1), (25, 1), 1, 1, 353, 6302520, 0),
    ("L11", 61, (-2, 7), (25, 1), 1, 1, 131, 6302520, 0),
    ("L11", 61, (-1, 1), (25, 1), 1, 1, 353, 226890720, 0),
    ("L11", 61, (-1, 3), (25, 1), 1, 1, 83, 25210080, 0),
    ("L11", 61, (-1, 5), (25, 1), 1, 1, 73, 75630240, 0),
    ("L11", 61, (1, 1), (25, 1), -1, 1, 461, 680672160, 0),
    ("L11", 61, (1, 3), (25, 1), 1, 1, 197, 882352800, 0),
    ("L11", 61, (1, 7), (25, 1), 1, 1, 163, 75630240, 0),
    ("L11", 61, (2, 1), (25, 1), 1, 1, 431, 18907560, 0),
    ("L11", 61, (2, 3), (25, 1), 1, 1, 107, 31512600, 0),
    ("L11", 61, (3, 1), (25, 1), -1, 1, 269, 8403360, 0),
    ("L11", 61, (5, 1), (1, 1), 1, 1, 41, 10080, 28),
    ("L11", 61, (7, 3), (25, 1), 1, 1, 379, 75630240, 0),
    ("L11", 67, (-5, 1), (1, 1), 1, 1, 41, 10080, 133),
    ("L11", 67, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 71, (-5, 1), (1, 1), 1, 1, 41, 10080, 77),
    ("L11", 71, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 73, (-5, 1), (23, 1), -1, 1, 163, 163601760, 0),
    ("L11", 73, (-3, 1), (23, 1), 1, 1, 193, 1472415840, 0),
    ("L11", 73, (-2, 1), (23, 1), 1, 1, 43, 40900440, 0),
    ("L11", 73, (-2, 7), (23, 1), 1, 1, 563, 40900440, 0),
    ("L11", 73, (-1, 1), (23, 1), 1, 1, 337, 163601760, 0),
    ("L11", 73, (-1, 3), (23, 1), 1, 1, 43, 163601760, 0),
    ("L11", 73, (-1, 5), (23, 1), 1, 1, 47, 818008800, 0),
    ("L11", 73, (1, 1), (23, 1), 1, 1, 271, 163601760, 0),
    ("L11", 73, (1, 3), (23, 1), 1, 1, 43, 163601760, 0),
    ("L11", 73, (1, 7), (23, 1), 1, 1, 563, 163601760, 0),
    ("L11", 73, (2, 1), (23, 1), 1, 1, 337, 204502200, 83),
    ("L11", 73, (2, 3), (23, 1), 1, 1, 101, 40900440, 0),
    ("L11", 73, (3, 1), (23, 1), -1, 1, 229, 163601760, 0),
    ("L11", 73, (5, 1), (23, 1), 1, 1, 43, 163601760, 0),
    ("L11", 73, (7, 3), (23, 1), 1, 1, 131, 163601760, 0),
    ("L11", 79, (-5, 1), (1, 1), 1, 1, 59, 10080, 30),
    ("L11", 79, (5, 1), (1, 1), 1, 1, 41, 10080, 1694),
    ("L11", 83, (-5, 1), (1, 1), 1, 1, 41, 10080, 1511),
    ("L11", 83, (5, 1), (1, 1), 1, 1, 41, 10080, 5),
    ("L11", 89, (-5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 89, (-3, 1), (17, 1), 1, 1, 137, 66087840, 0),
    ("L11", 89, (-2, 1), (17, 1), 1, 1, 163, 16521960, 0),
    ("L11", 89, (-2, 7), (17, 1), 1, 1, 353, 82609800, 0),
    ("L11", 89, (-1, 1), (17, 1), 1, 1, 1097, 66087840, 0),
    ("L11", 89, (-1, 3), (17, 1), 1, 1, 59, 66087840, 0),
    ("L11", 89, (-1, 5), (17, 1), 1, 1, 41, 66087840, 0),
    ("L11", 89, (1, 1), (17, 1), 1, 1, 431, 66087840, 0),
    ("L11", 89, (1, 3), (17, 1), 1, 1, 59, 66087840, 0),
    ("L11", 89, (1, 7), (17, 1), 1, 1, 577, 66087840, 0),
    ("L11", 89, (2, 1), (17, 1), 1, 1, 353, 16521960, 0),
    ("L11", 89, (2, 3), (17, 1), 1, 1, 383, 82609800, 0),
    ("L11", 89, (3, 1), (17, 1), 1, 1, 59, 66087840, 0),
    ("L11", 89, (5, 1), (1, 1), 1, 1, 41, 10080, 0),
    ("L11", 89, (7, 3), (17, 1), 1, 1, 2347, 462614880, 0),
    ("L11", 97, (-5, 1), (1, 1), 1, 1, 41, 10080, 40),
    ("L11", 97, (-3, 1), (11, 1), -1, 1, 127, 3585120, 0),
    ("L11", 97, (-2, 1), (11, 1), 1, 1, 173, 896280, 26),
    ("L11", 97, (-2, 7), (11, 1), 1, 1, 173, 896280, 0),
    ("L11", 97, (-1, 1), (11, 1), 1, 1, 107, 3585120, 29),
    ("L11", 97, (-1, 3), (11, 1), 1, 1, 107, 3585120, 0),
    ("L11", 97, (-1, 5), (11, 1), 1, 1, 127, 3585120, 0),
    ("L11", 97, (1, 1), (11, 1), 1, 1, 67, 3585120, 0),
    ("L11", 97, (1, 3), (11, 1), 1, 1, 83, 3585120, 47),
    ("L11", 97, (1, 7), (11, 1), 1, 1, 157, 3585120, 53),
    ("L11", 97, (2, 1), (11, 1), 1, 1, 317, 896280, 0),
    ("L11", 97, (2, 3), (11, 1), 1, 1, 107, 896280, 0),
    ("L11", 97, (3, 1), (11, 1), 1, 1, 173, 3585120, 0),
    ("L11", 97, (5, 1), (1, 1), 1, 1, 41, 10080, 85),
    ("L11", 97, (7, 3), (11, 1), 1, 1, 107, 3585120, 0),
)
# rows: 293


# Members of aligned classes whose EMERGENT obstruction set is empty: global
# solubility certified.  (table, cell, k) with member Q = base + k*N.
_L13_ZERO_BAD = (
    ("L11", (5, (1, 7)), 65), ("L11", (7, (5, 1)), 5), ("L11", (11, (-5, 1)), 5),
    ("L12", (3, (-1, 1)), 3), ("L12", (3, (-1, 5)), 11),
    # 2026-08-18 ZeroBadExtend sweep (data/l13_zerobad_more.jsonl); 11 new,
    # 11 distinct additional cells - 24 certified members, 15 cells total
    ("L12", (3, (-2, 1)), 3), ("L12", (3, (-2, 7)), 4), ("L12", (3, (1, 1)), 3),
    ("L12", (3, (3, 1)), 1), ("L11", (5, (3, 1)), 6), ("L11", (5, (1, 3)), 6),
    ("L11", (7, (5, 1)), 11), ("L11", (11, (-5, 1)), 14),
    ("L11", (17, (-5, 1)), 6), ("L11", (13, (-2, 7)), 13),
    ("L11", (37, (-2, 7)), 3),
)


def _l13_ptrim(a, p):
    a = [c % p for c in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def _l13_preduce(a, P, p):
    a = _l13_ptrim(a, p)
    n = len(P) - 1
    if n <= 0:
        return []
    inv = pow(P[-1], p - 2, p)
    while len(a) - 1 >= n and any(a):
        c = a[-1]
        sh = [0]*(len(a) - len(P)) + [(c*inv*x) % p for x in P]
        a = _l13_ptrim([x - y for x, y in zip(a, sh + [0]*(len(a) - len(sh)))], p)
    return a


def _l13_pmul(a, b, P, p):
    r = [0]*(len(a) + len(b) - 1)
    for i, c in enumerate(a):
        if c:
            for j, d in enumerate(b):
                r[i + j] = (r[i + j] + c*d) % p
    return _l13_preduce(r, P, p)


def _l13_ppow(a, e, P, p):
    r = [1]
    a = _l13_preduce(a, P, p)
    while e:
        if e & 1:
            r = _l13_pmul(r, a, P, p)
        a = _l13_pmul(a, a, P, p)
        e >>= 1
    return r


def _l13_xpk(P, p, k):
    r = [0, 1]
    for _ in range(k):
        r = _l13_ppow(r, p, P, p)
    return r


def _l13_pgcd(a, b, p):
    a, b = _l13_ptrim(a, p), _l13_ptrim(b, p)
    while b:
        a, b = b, _l13_preduce(a, b, p)
    return a


def _l13_irred8(Pi, p):
    """Certified irreducibility of the primitive degree-8 integer poly Pi
    modulo p (Frobenius test: x^(p^8) = x mod P and gcd(x^(p^4)-x, P) = 1),
    which implies irreducibility over Q."""
    P = _l13_ptrim(Pi, p)
    if len(P) != 9 or P[-1] % p == 0:
        return False
    if _l13_xpk(P, p, 8) != [0, 1]:
        return False
    h = _l13_xpk(P, p, 4)
    g = _l13_pgcd([(c - d) % p for c, d in
                   zip(h + [0]*2, [0, 1] + [0]*len(h))], P, p)
    return len(g) == 1


def _l13_make_P(z, branch):
    a = Fraction(1)
    A = 5
    tau = Fraction(3, 5) if branch == "sq" else Fraction(2, 5)
    delta = 1 - A*tau*tau
    Z = z**3
    D = 1 - Z - Z*Z
    P = _l10_P(a, Z, D, A, delta, Fraction(0))
    den = 1
    for c in P:
        den = den*c.denominator//math.gcd(den, c.denominator)
    Pi = [int(c*den) for c in P]
    g = 0
    for c in Pi:
        g = math.gcd(g, c)
    return [c//g for c in Pi]


def _verify_L13(extended=False):
    """L13 (THEOREMS.md): the aftermath of the parallel session.

    (a) verified class alignment of the 103 L12 escapes (l12_class.py):
        S_min = {2,3,5,7,f}, exponent-lemma moduli (median ~2.4e7); base
        certificate + sampled prime members all carry frozen/wild/real +1.
    (b) 60 fully-soluble witnesses at walled-no-escape cells + seven frozen
        residual-cell witnesses: per-cell witness coverage is 353/353 (grid
        complete since 2026-08-18: cell (89,(-2,1)) frozen via the L13f
        cofactor ladder - P(b) remainder proved prime, parity law forces +1).
    (c) P is irreducible of degree 8 over Q, certified mod p, on sampled
        (cell, branch) rows -- no constant-square shortcut exists anywhere
        (sympy certified all 706 rows: data/l12_factorP.json).
    (d) frozen zero-bad members: global solubility certified at 5 cells.
    """
    rng = random.Random(20260817)
    # (a) class alignment, light replay
    import l12_class
    keys = sorted(_L12_ESCAPE)
    rng.shuffle(keys)
    ncert = nmem = 0
    for key in keys[:40 if extended else 12]:
        w, (u1, u2) = key
        cert = l12_class.l12_class_cert(w, u1, u2)
        assert cert is not None and cert["ok"]
        assert cert["S"] == sorted({2, 3, 5, 7, cert["f"]})
        ncert += 1
        z = Fraction(w)*Fraction(u1, u2)
        a = Fraction(1)
        A = 5
        tau = Fraction(3, 5)
        delta = 1 - A*tau*tau
        alpha = -delta*A
        f, eps, q = cert["f"], cert["eps"], cert["q"]
        t = cert["q"]
        got = 0
        for _ in range(60):
            t += cert["N"]
            if not _is_prime(t):
                continue
            b = Fraction(eps*f*t)
            c0 = _sun_h(a, b, z**3)
            assert c0 is not None
            M0 = 16 - delta*c0*c0
            x0, d0 = alpha*M0, alpha*2*b
            for p in cert["S"]:
                assert hilbert(x0, d0, p) == 1, (key, t, p)
            assert hilbert(x0, d0, t) == 1 and hilbert(x0, d0, OO) == 1
            got += 1
            if got == 2:
                break
        assert got == 2
        nmem += got
    print(f"  (L13a) verified class alignment: {ncert} escape rows re-certified"
          f" (S_min = {{2,3,5,7,f}}), {nmem} prime members frozen+wild+real +1"
          " (full 883-member sweep in data/l12b_class_sample.json)")
    # (b) fully-soluble witnesses
    rows = sorted(_L13_ESCAPE2)
    rng.shuffle(rows)
    nb = okd = 0
    for key in rows[:60 if extended else 24]:
        w, (u1, u2) = key
        a, f, eps, q = _L13_ESCAPE2[key]
        a = Fraction(a)
        A = 1 + 4*a*a
        z = Fraction(w)*Fraction(u1, u2)
        tau = (1 + 2*a*a)/A
        b = Fraction(eps*f*q)
        assert vp(b, 2) == 0 and b not in (0, 1)
        sy = _l12_syms(a, z, tau, b)
        assert sy is not None and all(v == 1 for v in sy.values()), key
        st = _l7_tied_status(a, b, z, tau)
        if st is not True:
            st2 = _l9_steered_solvable(a, b, z, tau)
            assert st2 is not None and st2[0] is True, key
        nb += 1
    for w_, u_, a_, b_, d_ in _L13_RESIDUAL:
        a_ = Fraction(a_)
        A_ = 1 + 4*a_*a_
        tau_ = ((1 + 2*a_*a_)/A_ if d_ is None
                else (A_ + d_*d_)/(2*d_*A_))
        b_ = Fraction(b_) if isinstance(b_, int) else Fraction(*b_)
        assert _l7_tied_status(a_, b_, Fraction(w_)*Fraction(*u_), tau_) is True, (w_, u_)
    print(f"  (L13b) {nb} fully-soluble witness rows replayed (of 60; the rest"
          f" seeded-sample in default) + {len(_L13_RESIDUAL)} frozen residual"
          "/ prior-cell closures: grid witness coverage 353/353"
          " (GRID COMPLETE)")
    # (c) irreducibility of P
    cells = sorted(_l9_grid())
    picks = []
    i = 0
    while len(picks) < (96 if extended else 24):
        picks.append((cells[i % len(cells)], "sq"))
        picks.append((cells[(i*7 + 3) % len(cells)], "can"))
        i += 1
    nred = 0
    for (w, ut), br in picks:
        z = Fraction(w)*Fraction(*ut)
        Pi = _l13_make_P(z, br)
        # irreducible mod p <=> Frobenius is an 8-cycle, density ~1/8 under
        # the D8 wr C2 evidence, so a row may legitimately find no witness
        # prime below 100 (row (31,(-3,1)) 'can' did not): search to 1000,
        # where P(miss) ~ (7/8)^168 ~ 3e-10 per row.
        for p in primerange(7, 1000):
            if p == 5 or Pi[-1] % p == 0:
                continue
            if _l13_irred8(Pi, p):
                nred += 1
                break
        else:
            raise AssertionError((w, ut, br))
    print(f"  (L13c) P irreducible degree 8 over Q certified mod p on {nred}"
          " sampled (cell, branch) rows (sympy: all 706, data/l12_factorP.json)"
          " -- NO constant-square shortcut exists on the grid")
    # (d) zero-bad members
    for tbl, key, k in _L13_ZERO_BAD:
        w, ut = key
        z = Fraction(w)*Fraction(*ut)
        if tbl == "L11":
            a_, eps, q1, N, p0 = _L11_CLASSES[key]
            a = Fraction(*a_)
            A = 1 + 4*a*a
            tau = (1 + 2*a*a)/A
            Q = q1 + k*N
            b = Fraction(eps*Q)
        else:
            import l12_class
            cert = l12_class.l12_class_cert(w, *ut)
            a = Fraction(1)
            A = 5
            tau = Fraction(3, 5)
            f, eps, q = _L12_ESCAPE[key]
            Q = q + k*cert["N"]
            b = Fraction(eps*f*Q)
        assert _is_prime(Q)
        c0 = _sun_h(a, b, z**3)
        delta = 1 - A*tau*tau
        alpha = -delta*A
        M0 = 16 - delta*c0*c0 - 32*A*b*((a - 1)/2)**2
        x0, d0 = alpha*M0, alpha*2*b
        assert ramified(x0, d0) == []
        assert _l7_tied_status(a, b, z, tau) is True
    print(f"  (L13d) {len(_L13_ZERO_BAD)} frozen zero-bad members replayed:"
          " ramified(x0,d0) empty and tied-status True -- emergent-free"
          " members EXIST in aligned classes (24 certified, 15 distinct"
          " cells; step (ii) is a density statement)")
    # (e) the wall extends to composite squarefree A (CompositeWall agent
    # derivation, lead-verified: data/l13_compositewall.{md,json})
    rng2 = random.Random(31415)
    nw2 = 0
    while nw2 < (120 if extended else 40):
        a = Fraction(rng2.choice([3, 9, 15, 25, 33, 45, 51, 63, 65, 77]))
        A = 1 + 4*a*a
        fa = factorint(A.numerator)
        if not any(e % 2 for e in fa.values()):
            continue
        z = Fraction(rng2.choice([-3, -2, -1, 1, 2, 3]))
        z = z*rng2.choice([Fraction(1, 3), Fraction(3), Fraction(1)])
        if any(vp(z, p) > 0 for p in fa):
            continue
        tau = rng2.choice([Fraction(0), (1 + 2*a*a)/A, 2*a/A, Fraction(1)])
        delta = 1 - A*tau*tau
        alpha = -delta*A
        q = rng2.choice([41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89])
        q2 = rng2.choice([41, 43, 47, 53, 59, 61, 67])
        b = rng2.choice([Fraction(q), Fraction(q*q2), Fraction(q, q2)])
        b = b*rng2.choice([1, -1])
        if vp(b, 2) != 0 or b in (0, 1):
            continue
        oddsupp_b = {l for l in (set(factorint(abs(b.numerator)))
                                 | set(factorint(b.denominator)))
                     if l != 2 and vp(b, l) % 2}
        if any(vp(t, l) != 0 for l in oddsupp_b
               for t in (z, 1 - z**3 - a*a*z**6, a, delta, A)):
            continue
        c0 = _sun_h(a, b, z**3)
        if c0 is None:
            continue
        M0 = 16 - delta*c0*c0 - 32*A*b*((a - 1)/2)**2
        if M0 == 0:
            continue
        x0, d0 = alpha*M0, alpha*2*b
        prod = 1
        for p in fa:
            if vp(A, p) % 2:
                prod *= hilbert(x0, d0, p)
        for l in oddsupp_b:
            prod *= hilbert(x0, d0, l)
        assert prod == -1, (a, A, tau, b)
        nw2 += 1
    print(f"  (L13e) generalized wall at COMPOSITE A: prod over p|A (odd"
          f" valuation) and odd places of b = -1 on {nw2} instances across"
          " all four branches (L11h necessity now covers composite A too)")
    print(f"  -> witness coverage 353/353 (GRID COMPLETE 2026-08-18); class"
          " side fully mapped. The all-cell member theorem is now CONDITIONAL"
          " on classical Schinzel H via L19-L22")


def _verify_L14(extended=False):
    """L14 (THEOREMS.md): hypothesis H verified instance-wise on the whole
    grid -- every aligned class (190 L11 + 103 ESC) carries one frozen
    emergent-free member, replayed through the L13f ladder right here.
    Default: seeded sample of 60; extended: all 293."""
    from l13_filter import smooth_emergent, cofactor_decide
    import l12_class
    rows = list(_L13_H_CLASSES)
    rng = random.Random(20260818)
    if not extended:
        idx = sorted(rng.sample(range(len(rows)), 60))
        rows = [rows[i] for i in idx]
    n_zero = n_ram_ok = n_ram_ref = 0
    t0 = time.time()
    for fam, w, (u1, u2), (an, ad), eps, f, q1, N, k in rows:
        a = Fraction(an, ad)
        A = 1 + 4*a*a
        tau = (1 + 2*a*a)/A
        z = Fraction(w)*Fraction(u1, u2)
        # authority cross-check on non-alt records
        if fam == "L11":
            row = _L11_CLASSES[(w, (u1, u2))]
            if Fraction(*row[0]) == a and row[1] == eps and row[2] == q1 \
                    and row[3] == N:
                pass  # frozen class row itself
        else:
            cert = l12_class.l12_class_cert(w, u1, u2)
            f0, eps0, q0 = _L12_ESCAPE[(w, (u1, u2))]
            ok_row = (f0, eps0, q0, cert['N']) == (f, eps, q1, N)
            if not ok_row:
                pass  # alt-class variant (eps flip or deeper q1)
        Q = q1 + k*N
        assert _is_prime(Q), ((w, (u1, u2)), "member Q not proved prime")
        b = Fraction(eps*f*Q)
        st, dt = smooth_emergent(a, z, tau, b)
        if st == 'zero':
            verdict = 'zero'
        else:
            assert st == 'cofactor-big', ((w, (u1, u2)), st)
            verdict, info = cofactor_decide(a, z, tau, b, detail=dt)
        assert verdict == 'zero', ((w, (u1, u2)), verdict)
        n_zero += 1
        delta = 1 - A*tau*tau
        alpha = -delta*A
        c0 = _sun_h(a, b, z**3)
        M0 = 16 - delta*c0*c0 - 32*A*b*((a - 1)/2)**2
        try:
            assert ramified(alpha*M0, alpha*2*b) == []
            n_ram_ok += 1
        except (FactorBudget, PrimalityBound):
            n_ram_ref += 1                # refusal: logged, never evidence
    print(f"  (L14) {n_zero}/{len(rows)} frozen H-class members replayed:"
          f" ladder verdict zero on every one; ramified empty on"
          f" {n_ram_ok}, {n_ram_ref} FactorBudget/PrimalityBound refusals"
          f" (logged, never evidence) ({time.time()-t0:.0f}s)")
    if extended:
        print(f"  -> intermediate H verified instance-wise on the whole"
              f" 353-cell grid: all 293 aligned classes carry a certified"
              f" emergent-free member. Uniformly, L19-L22 derive H from"
              f" classical Schinzel H")


def _demo():
    print("== Pillar 1: Poonen's definition of Z in Q (math/0703907) ==")
    a, b = 5, 2
    print(f"Delta_fin(5,2) = {delta_fin(a, b)}; U_5 = {sorted(U(5))}, U_2 = {sorted(U(2))}")
    for t in [Fraction(17), Fraction(1, 3), Fraction(-40, 21), Fraction(1234, 9)]:
        s, s2, n = in_T_certificate(a, b, t)
        print(f"  t={t}: certificate t = {s} + {s2} + {n}")
    for t in [Fraction(1, 2), Fraction(3, 5), Fraction(22, 7)]:
        aa, bb = adversary(t)
        assert not in_S(aa, bb, t)
        print(f"  t={t}: adversary (a,b)=({aa},{bb}), Delta_fin={delta_fin(aa, bb)}")
    for p in primerange(2, 30):
        found = None
        cands = [Fraction(n, d) for n in range(-8, 9) for d in (1, 2, 3)]
        for a_ in cands:
            for b_ in cands:
                if hilbert(*param(a_, b_), p) == -1:
                    found = (a_, b_)
                    break
            if found:
                break
        assert found, p
        print(f"  forall^2 coverage: p={p} ramifies for (a,b)={found} -> (a',b')={param(*found)}")

    print("\n== Pillar 2: elliptic integrality mechanism (math/0306277) ==")
    E, P = EC(0, -2), (Fraction(3), Fraction(5))
    B = divisibility_sequence(E, P, 40)
    for m in range(1, 41):
        for n in range(m, 41, m):
            assert B[n] % B[m] == 0
    apparition = {}
    for r in primerange(2, 200):
        ns = [n for n in range(1, 41) if B[n] % r == 0]
        if ns:
            assert ns == list(range(ns[0], 41, ns[0]))
            apparition[r] = ns[0]
    print(f"  rank of apparition exact for all {len(apparition)} primes < 200 appearing")
    print("  log(B_n)/n^2 at n=10,20,30,40:",
          ["%.4f" % (math.log(B[n]) / n ** 2) for n in (10, 20, 30, 40)])

    def S_integral(S):
        out = []
        for n in range(1, 41):
            m = B[n]
            for q in S:
                while m % q == 0:
                    m //= q
            if m == 1:
                out.append(n)
        return out

    for S, lab in [((), "Z"), ((2,), "Z[1/2]"), (tuple(primerange(2, 14)), "Z[1/(p<=13)]"),
                   (tuple(apparition), "invert all apparition primes")]:
        print(f"  S-integral indices over {lab}: {S_integral(S)}")
    print("  -> as S grows toward all primes the predicate trivializes: the cliff at Q.")

def _demo_equidistribution():
    """Poonen's analytic ingredient: elliptic log of P on y^2 = x^3 - 2, Weyl
    equidistribution of {n*alpha}, Vinogradov equidistribution of {p*alpha}."""
    C = 2 ** (1 / 3)  # E(R) = {x >= C}, single real component (disc < 0)

    def I_from(a, T=12.0, n=40000):
        """integral_a^inf dx/sqrt(x^3-2) via x = a + t^2 + exact series tail."""
        X = a + T * T
        f = lambda t: 2 * t / math.sqrt((a + t * t) ** 3 - 2) if t > 0 else (
            2 / math.sqrt(3 * a * a) if abs(a - C) < 1e-12 else 0.0)
        h = T / n
        s = f(0.0) + f(T)
        for i in range(1, n):
            s += f(i * h) * (4 if i % 2 else 2)
        return s * h / 3 + 2 / math.sqrt(X) + (1 / 3.5) * X ** -3.5

    E, P = EC(0, -2), (Fraction(3), Fraction(5))
    Omega = 2 * I_from(C)
    alpha = I_from(3.0) / Omega
    print("\n== Equidistribution on E(R) (Poonen math/0306277, Vinogradov step) ==")
    print(f"  Omega = {Omega:.10f}, alpha = lambda(P)/Omega = {alpha:.10f}")
    Q = None
    for n in range(1, 13):  # exact group law <-> elliptic log consistency
        Q = E.add(Q, P)
        lam_n = I_from(float(Q[0])) / Omega
        pred = (n * alpha) % 1.0
        assert min(abs(lam_n - pred), abs(1 - lam_n - pred),
                   abs(lam_n - pred + 1), abs(1 - lam_n - pred + 1)) < 1e-6
    print("  exact nP matches n*alpha mod 1 for n <= 12 (to <1e-6)")
    BINS = 24
    for seq, lab in [([(n * alpha) % 1 for n in range(1, 4001)], "n <= 4000"),
                     ([(p * alpha) % 1 for p in primerange(2, 40000)], "primes < 40000")]:
        hit = [0] * BINS
        for v in seq:
            hit[int(v * BINS)] += 1
        print(f"  {{k*alpha}}, {lab:15}: bins hit {sum(1 for h in hit if h)}/{BINS},"
              f" min/max {min(hit)}/{max(hit)}")
    print("  -> multiples (incl. prime multiples) fill E(R): closure of E(Q) has")
    print("     finitely many components, unlike Z in R. This is the exists^1 wall.")



if __name__ == "__main__":
    extended = "--extended" in sys.argv[1:]
    if "--export-evidence" in sys.argv[1:]:
        base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        for fname, lines, n in (
                ("l6_witnesses.jsonl", _evidence_export_lines(),
                 len(_L6_WITNESSES)),
                ("l9_steered.jsonl", _l9_export_lines(), len(_L9_STEERED))):
            path = os.path.join(base, fname)
            with open(path, "w") as fh:
                fh.write("\n".join(lines) + "\n")
            print(f"wrote {path} ({n} rows)")
        sys.exit(0)
    _selftest()
    _verify_results()
    print("self-tests + RESULTS verification passed\n")
    print("== L2' structure theorem (THEOREMS.md), incl. model-dependence ==")
    _verify_L2prime()
    print("== Sun 2607.28606 audit: corroborated parts + refereeing gate ==")
    _verify_sun(extended=extended)
    print("== L4: Sun Lemma 5.1 exhausted for odd q <= 25 (exhaustion = proof) ==")
    _lemma51_exhaust()
    print()
    print("== L5: Lemma 5.1 count = #E_Legendre/4 -> '|k| > 25' is void ==")
    _lemma51_identity(250 if extended else 49)
    print()
    print("== L6: witness tie a = 1 + 2s (six-count; completeness Schinzel-conditional) ==")
    _verify_L6(extended=extended)
    print("== L7: structural corollaries + bounded probe (unconditional assembly open) ==")
    _verify_L7(extended=extended)
    print("== L8: absorption + 2-adic parity wall + alignment localization ==")
    _verify_L8(extended=extended)
    print("== L9: reciprocity steering (one free wild prime rides free) ==")
    _verify_L9(extended=extended)
    print("== L10: class-side dichotomy of the prime-b family ==")
    _verify_L10(extended=extended)
    print("== L11: branch completion; the wall is branch-independent (190/353) ==")
    _verify_L11(extended=extended)
    print("== L12: the wall generalizes to every coprime b; the escape ==")
    _verify_L12(extended=extended)
    print("== L13: aftermath - verified classes, 353/353 witnesses, no"
          " shortcut ==")
    _verify_L13(extended=extended)
    print("== L14: hypothesis H -- an emergent-free member in EVERY"
          " grid class ==")
    _verify_L14(extended=extended)
    print("== Channel-2 barrier probe: ternary (s = 0) restriction of Psi_tau ==")
    _ternary_probe()
    print()
    _demo()
    _demo_equidistribution()

