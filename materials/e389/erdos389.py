"""Exact arithmetic for Erdős Problem #389.

The shifted problem asks whether, for every m >= 0, some k >= 1 satisfies

    binom(m + k, m) | binom(m + 2k, k).

All public predicates in this module use integer arithmetic.  The default
finite-factor path is deliberately capped at ``MAX_TRIAL_FACTORIZER_VALUE`` so
its retained prime sieve cannot grow without bound.  No claim about the open
universal quantifier is encoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb, isqrt
from typing import Mapping

MAX_TRIAL_FACTORIZER_VALUE = 100_000_000_000_000
_MILLER_RABIN_BASES_64 = (2, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022)


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _validate_pair(m: int, k: int) -> None:
    _require_int("m", m, 0)
    _require_int("k", k, 1)


def primes_up_to(limit: int) -> list[int]:
    """Return all primes <= limit by an exact sieve."""
    _require_int("limit", limit, 0)
    if limit < 2:
        return []
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            flags[start : limit + 1 : p] = b"\x00" * (((limit - start) // p) + 1)
    return [p for p in range(2, limit + 1) if flags[p]]


def _is_prime_64(n: int) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p

    d = n - 1
    shifts = 0
    while d % 2 == 0:
        d //= 2
        shifts += 1
    for base in _MILLER_RABIN_BASES_64:
        if base % n == 0:
            continue
        x = pow(base, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(shifts - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


class TrialFactorizer:
    """Reusable exact trial factorizer with an explicit resource bound.

    ``maximum`` may not exceed ``MAX_TRIAL_FACTORIZER_VALUE`` (10^14), which
    bounds the retained sieve at 10^7. Deterministic 64-bit Miller--Rabin
    short-circuits prime cofactors; trial division certifies composites.
    """

    def __init__(self, maximum: int):
        _require_int("maximum", maximum, 1)
        if maximum > MAX_TRIAL_FACTORIZER_VALUE:
            raise ValueError(
                f"maximum must not exceed {MAX_TRIAL_FACTORIZER_VALUE}"
            )
        self.maximum = maximum
        self._primes = tuple(primes_up_to(isqrt(maximum)))

    def factor(self, n: int) -> dict[int, int]:
        _require_int("n", n, 1)
        if n > self.maximum:
            raise ValueError(f"n={n} exceeds factorizer maximum {self.maximum}")
        if n == 1:
            return {}

        factors: dict[int, int] = {}
        remainder = n
        if _is_prime_64(remainder):
            return {remainder: 1}

        for p in self._primes:
            if p * p > remainder:
                break
            if remainder % p:
                continue
            exponent = 0
            while remainder % p == 0:
                remainder //= p
                exponent += 1
            factors[p] = exponent
            if remainder == 1:
                break
            if _is_prime_64(remainder):
                factors[remainder] = factors.get(remainder, 0) + 1
                remainder = 1
                break

        if remainder > 1:
            # The prime table reaches sqrt(maximum), hence any final cofactor is prime.
            factors[remainder] = factors.get(remainder, 0) + 1
        return factors


class SPFFactorizer:
    """Smallest-prime-factor table for bounded sequential searches."""

    def __init__(self, maximum: int):
        _require_int("maximum", maximum, 2)
        self.maximum = maximum
        spf = [0] * (maximum + 1)
        spf[1] = 1
        for p in range(2, maximum + 1):
            if spf[p] != 0:
                continue
            spf[p] = p
            if p > maximum // p:
                continue
            for multiple in range(p * p, maximum + 1, p):
                if spf[multiple] == 0:
                    spf[multiple] = p
        self._spf = spf

    def factor(self, n: int) -> dict[int, int]:
        _require_int("n", n, 1)
        if n > self.maximum:
            raise ValueError(f"n={n} exceeds factorizer maximum {self.maximum}")
        factors: dict[int, int] = {}
        while n > 1:
            p = self._spf[n]
            exponent = 0
            while n % p == 0:
                n //= p
                exponent += 1
            factors[p] = exponent
        return factors


def integer_valuation(n: int, p: int) -> int:
    """Return v_p(n); callers use prime p, but primality is not retested."""
    _require_int("n", n, 1)
    _require_int("p", p, 2)
    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    return exponent


def factorial_valuation(n: int, p: int) -> int:
    """Return v_p(n!) using Legendre's formula."""
    _require_int("n", n, 0)
    _require_int("p", p, 2)
    valuation = 0
    while n:
        n //= p
        valuation += n
    return valuation


def binomial_valuation(n: int, r: int, p: int) -> int:
    """Return v_p(binomial(n, r)) using Legendre's formula."""
    _require_int("n", n, 0)
    _require_int("r", r, 0)
    _require_int("p", p, 2)
    if r > n:
        raise ValueError("r must not exceed n")
    return (
        factorial_valuation(n, p)
        - factorial_valuation(r, p)
        - factorial_valuation(n - r, p)
    )


def carry_count(a: int, b: int, p: int) -> int:
    """Count carries when adding a and b in base p, least digit first."""
    _require_int("a", a, 0)
    _require_int("b", b, 0)
    _require_int("p", p, 2)
    carries = 0
    incoming = 0
    while a or b or incoming:
        total = a % p + b % p + incoming
        incoming = int(total >= p)
        carries += incoming
        a //= p
        b //= p
    return carries


def slack(m: int, k: int, p: int) -> int:
    """Return the p-adic valuation of the right binomial divided by the left."""
    _validate_pair(m, k)
    _require_int("p", p, 2)
    return binomial_valuation(m + 2 * k, k, p) - binomial_valuation(m + k, m, p)

def local_floor_contribution(m: int, x: int, q: int) -> int:
    """Return one floor-sum contribution for

    ``x! m! / (((x+m)/2)!)^2`` when ``x`` and ``m`` have the same parity.
    """
    _require_int("m", m, 0)
    _require_int("x", x, 0)
    _require_int("q", q, 2)
    if (x - m) % 2:
        raise ValueError("x and m must have the same parity")
    return x // q + m // q - 2 * ((x + m) // (2 * q))
def zone_contribution(m: int, k: int, q: int) -> int:
    """Classify one prime-power level from residue/carry zones.

    For ``q <= m`` this is ``c2 - 2*c1`` with the two residue carries.
    For ``q > m`` it is the bad/neutral/good interval classification of
    ``(m+k) mod q``.
    """
    _validate_pair(m, k)
    _require_int("q", q, 2)
    if q <= m:
        m_residue = m % q
        k_residue = k % q
        first_carry = (m_residue + k_residue) // q
        second_carry = (m_residue + 2 * k_residue) // q
        return second_carry - 2 * first_carry

    n_residue = (m + k) % q
    if 2 * n_residue < m:
        return -1
    if 2 * n_residue >= q + m:
        return 1
    return 0




def slack_floor_contributions(m: int, k: int, p: int) -> tuple[tuple[int, int], ...]:
    """Return ``(p^a, D_{p^a})`` terms whose sum is ``slack(m,k,p)``."""
    _validate_pair(m, k)
    _require_int("p", p, 2)
    x = m + 2 * k
    contributions: list[tuple[int, int]] = []
    q = p
    while q <= x:
        contributions.append((q, local_floor_contribution(m, x, q)))
        q *= p
    return tuple(contributions)

def _add_factors(target: dict[int, int], factors: Mapping[int, int], scale: int) -> None:
    for p, exponent in factors.items():
        updated = target.get(p, 0) + scale * exponent
        if updated:
            target[p] = updated
        else:
            target.pop(p, None)


def left_valuations(
    m: int, k: int, factorizer: TrialFactorizer | SPFFactorizer | None = None
) -> dict[int, int]:
    """Factor binomial(m+k, m) via its m numerator terms and m!.

    Without an injected factorizer, ``m+k`` must not exceed
    ``MAX_TRIAL_FACTORIZER_VALUE``.
    """
    _validate_pair(m, k)
    if m == 0:
        return {}
    if factorizer is None:
        factorizer = TrialFactorizer(m + k)

    valuations: dict[int, int] = {}
    for value in range(k + 1, k + m + 1):
        _add_factors(valuations, factorizer.factor(value), 1)
    for value in range(2, m + 1):
        _add_factors(valuations, factorizer.factor(value), -1)
    if any(exponent < 0 for exponent in valuations.values()):
        raise AssertionError("factor cancellation contradicted binomial integrality")
    return dict(sorted(valuations.items()))


@dataclass(frozen=True)
class WitnessCertificate:
    m: int
    k: int
    left_valuations: dict[int, int]
    right_valuations: dict[int, int]
    slacks: dict[int, int]

    @property
    def obstructions(self) -> dict[int, int]:
        return {p: value for p, value in self.slacks.items() if value < 0}

    @property
    def tight_primes(self) -> tuple[int, ...]:
        return tuple(p for p, value in self.slacks.items() if value == 0)

    @property
    def is_witness(self) -> bool:
        return not self.obstructions


def witness_certificate(
    m: int, k: int, factorizer: TrialFactorizer | SPFFactorizer | None = None
) -> WitnessCertificate:
    """Return every relevant prime valuation and slack.

    The default factorizer requires ``m+k <= MAX_TRIAL_FACTORIZER_VALUE``.
    """
    left = left_valuations(m, k, factorizer)
    right = {p: binomial_valuation(m + 2 * k, k, p) for p in left}
    slacks = {p: right[p] - left[p] for p in left}
    return WitnessCertificate(m, k, left, right, slacks)


def is_witness(
    m: int, k: int, factorizer: TrialFactorizer | SPFFactorizer | None = None
) -> bool:
    """Check divisibility exactly within the injected/default factorizer domain."""
    return witness_certificate(m, k, factorizer).is_witness


def is_witness_by_carries(
    m: int, k: int, factorizer: TrialFactorizer | SPFFactorizer | None = None
) -> bool:
    """Check the same divisibility through Kummer carry counts.

    The default factorizer requires ``m+k <= MAX_TRIAL_FACTORIZER_VALUE``.
    """
    left = left_valuations(m, k, factorizer)
    return all(carry_count(k, m + k, p) >= carry_count(k, m, p) for p in left)


def is_witness_direct(m: int, k: int, *, max_n: int = 20_000) -> bool:
    """Direct binomial check for bounded independent tests.

    The guard prevents accidental construction of a central binomial with
    billions of digits. The finite-factor checkers avoid that construction but
    enforce their factorizer's explicit maximum.
    """
    _validate_pair(m, k)
    _require_int("max_n", max_n, 1)
    if m + 2 * k > max_n:
        raise ValueError(f"direct check requires m + 2k <= {max_n}")
    return comb(m + 2 * k, k) % comb(m + k, m) == 0


@dataclass(frozen=True)
class ShiftObstruction:
    prime: int
    source_slack: int
    m_plus_one_valuation: int
    m_plus_two_k_valuation: int
    shifted_slack: int


def natural_shift_obstructions(
    m: int, k: int, factorizer: TrialFactorizer | SPFFactorizer | None = None
) -> tuple[ShiftObstruction, ...]:
    """Obstructions to sending a witness (m,k) to (m+1,k-1).

    The default factorizer requires ``m+2k <= MAX_TRIAL_FACTORIZER_VALUE``.
    """
    _validate_pair(m, k)
    if k < 2:
        raise ValueError("natural shift requires k >= 2")
    if factorizer is None:
        factorizer = TrialFactorizer(m + 2 * k)
    if not is_witness(m, k, factorizer):
        raise ValueError("natural shift analysis requires a source witness")

    obstructions: list[ShiftObstruction] = []
    for p, total_valuation in sorted(factorizer.factor(m + 2 * k).items()):
        source_slack = slack(m, k, p)
        added = integer_valuation(m + 1, p)
        shifted_slack = source_slack + added - total_valuation
        if shifted_slack < 0:
            obstructions.append(
                ShiftObstruction(p, source_slack, added, total_valuation, shifted_slack)
            )
    return tuple(obstructions)


class RatioSlackScanner:
    """Maintain every prime exponent of

        binom(m+2k,k) / binom(m+k,m)

    while k increases.  A state is a witness exactly when no exponent is
    negative.  The recurrence multiplies the ratio by
    (m+2k+1)(m+2k+2)/(m+k+1)^2 when advancing from k to k+1.
    """

    def __init__(self, m: int, factorizer: SPFFactorizer):
        _require_int("m", m, 0)
        self.m = m
        self.k = 1
        self._factorizer = factorizer
        self._exponents: dict[int, int] = {}
        self._negative_primes: set[int] = set()
        self._apply(m + 2, 1)
        self._apply(m + 1, -1)

    @property
    def exponents(self) -> dict[int, int]:
        return dict(sorted(self._exponents.items()))

    @property
    def negative_exponents(self) -> dict[int, int]:
        return {p: self._exponents[p] for p in sorted(self._negative_primes)}
    @property
    def negative_count(self) -> int:
        return len(self._negative_primes)


    @property
    def is_witness(self) -> bool:
        return not self._negative_primes

    def _apply(self, value: int, scale: int) -> None:
        for p, exponent in self._factorizer.factor(value).items():
            old = self._exponents.get(p, 0)
            new = old + scale * exponent
            if old < 0 <= new:
                self._negative_primes.remove(p)
            elif old >= 0 > new:
                self._negative_primes.add(p)
            if new:
                self._exponents[p] = new
            else:
                self._exponents.pop(p, None)

    def advance(self) -> None:
        old_k = self.k
        self._apply(self.m + 2 * old_k + 1, 1)
        self._apply(self.m + 2 * old_k + 2, 1)
        self._apply(self.m + old_k + 1, -2)
        self.k += 1


def first_witness(m: int, limit: int, factorizer: SPFFactorizer | None = None) -> int | None:
    """Return the first witness k <= limit, or None after exact exhaustion."""
    _require_int("m", m, 0)
    _require_int("limit", limit, 1)
    if factorizer is None:
        factorizer = SPFFactorizer(m + 2 * limit)
    scanner = RatioSlackScanner(m, factorizer)
    while True:
        if scanner.is_witness:
            return scanner.k
        if scanner.k == limit:
            return None
        scanner.advance()
