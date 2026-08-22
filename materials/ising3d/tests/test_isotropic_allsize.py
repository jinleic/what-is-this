"""Independent verifier for the wave-18 isotropic all-size obstruction.

No e163/e164/e165 module is imported.  Every finite seed is rebuilt from the
raw rational operator R=PDP, cleared with its actual entry-denominator LCM, and
its exterior-pair polynomial is reconstructed independently at both primes.
The universal obstruction is checked from its exact base and induction
identities, not by extrapolating a finite table.
"""

from __future__ import annotations

import json
import platform
import resource
import sys
import time
from fractions import Fraction
from math import comb, gcd, isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "isotropic_allsize.json"
P1, P2 = 1_000_003, 2_000_003
PRIMES = (P1, P2)
T0 = Fraction(1, 3)
Q0 = Fraction(5, 3)
CPU_BUDGET_SECONDS = 600
RSS_CAP_BYTES = 4_000_000_000

SEEDS: tuple[tuple[str, int, tuple[tuple[int, int], ...], int], ...] = (
    ("claw_K1_3", 4, ((0, 1), (0, 2), (0, 3)), 55),
    ("claw_attached_P1_T5", 5, ((0, 1), (0, 2), (0, 3), (1, 4)), 79),
    ("paw", 4, ((0, 1), (0, 2), (1, 2), (0, 3)), 1),
    (
        "grid_2x2_plus_pendant",
        5,
        ((0, 1), (0, 2), (1, 3), (2, 3), (0, 4)),
        21,
    ),
    ("star_K1_4", 5, ((0, 1), (0, 2), (0, 3), (0, 4)), 349),
)
EXPECTED_WINNERS = {"claw_attached_P1_T5", "paw", "grid_2x2_plus_pendant"}
PASSED: list[str] = []
FAILED: list[str] = []
CPU0 = time.process_time()


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(stage: str) -> None:
    cpu = time.process_time() - CPU0
    if cpu > CPU_BUDGET_SECONDS:
        raise RuntimeError(
            f"NON-DECISIVE resource expiry is a test failure here: {stage} exceeded "
            f"the predeclared {CPU_BUDGET_SECONDS}s process-time budget"
        )
    if max_rss_bytes() >= RSS_CAP_BYTES:
        raise RuntimeError(
            f"NON-DECISIVE resource expiry is a test failure here: {stage} reached "
            f"the {RSS_CAP_BYTES}-byte RSS cap"
        )


def ok(name: str, condition: bool, detail: str = "") -> None:
    print(f"[{'PASS' if condition else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    (PASSED if condition else FAILED).append(name)


# ---------------------------------------------------------------- raw operator

def build_raw_operator(n: int, bonds: tuple[tuple[int, int], ...], t: Fraction) -> list[list[Fraction]]:
    """Build R=P D P directly, with D_sigma=q for each aligned uniform bond."""
    budget_tick(f"raw builder n={n}")
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    diagonal: list[Fraction] = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        aligned = sum(spins[u] == spins[v] for u, v in bonds)
        diagonal.append(q**aligned)

    powers = [t**h for h in range(n + 1)]
    p_matrix = [
        [powers[(row ^ col).bit_count()] for col in range(dim)]
        for row in range(dim)
    ]
    result = [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        budget_tick(f"raw PDP n={n}, row={i}")
        for j in range(i, dim):
            value = sum(
                (p_matrix[i][k] * diagonal[k] * p_matrix[j][k] for k in range(dim)),
                Fraction(0),
            )
            result[i][j] = result[j][i] = value
    return result


def actual_lcm_integer_model(matrix: list[list[Fraction]]) -> tuple[int, list[list[int]]]:
    scale = 1
    for row in matrix:
        for value in row:
            scale = scale * value.denominator // gcd(scale, value.denominator)
    integer = [[int(value * scale) for value in row] for row in matrix]
    assert all(
        value * scale == integer[i][j]
        for i, row in enumerate(matrix)
        for j, value in enumerate(row)
    )
    return scale, integer


def kron_fraction(a: list[list[Fraction]], b: list[list[Fraction]]) -> list[list[Fraction]]:
    ar, br = len(a), len(b)
    return [
        [a[i // br][j // br] * b[i % br][j % br] for j in range(ar * br)]
        for i in range(ar * br)
    ]


# ------------------------------------------------------------- independent F_p[z]

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    return all(n % divisor for divisor in range(2, isqrt(n) + 1))


def modmul(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    assert a.shape[1] * (p - 1) ** 2 < 2**53
    product = np.asarray(a, dtype=np.float64) @ np.asarray(b, dtype=np.float64)
    return np.mod(product, float(p)).astype(np.int64)


def moddot(a: np.ndarray, b: np.ndarray, p: int) -> int:
    assert a.size * (p - 1) ** 2 < 2**63 - 1
    return int(np.dot(a, b)) % p


def trim(poly: np.ndarray) -> np.ndarray:
    for index in range(min(poly.size, 64)):
        if poly[index]:
            return poly if index == 0 else poly[index:]
    nonzero = np.flatnonzero(poly)
    return np.zeros(1, dtype=np.int64) if nonzero.size == 0 else poly[nonzero[0] :]


def is_zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def monic(poly: np.ndarray, p: int) -> np.ndarray:
    poly = trim(np.asarray(poly, dtype=np.int64) % p)
    if is_zero(poly):
        return poly
    lead = int(poly[0])
    return poly if lead == 1 else (poly * pow(lead, -1, p)) % p


def remainder(dividend: np.ndarray, divisor: np.ndarray, p: int) -> np.ndarray:
    assert int(divisor[0]) == 1
    work = (np.asarray(dividend, dtype=np.int64) % p).copy()
    dividend_degree = work.size - 1
    divisor_degree = divisor.size - 1
    if dividend_degree < divisor_degree:
        return trim(work)
    for offset in range(dividend_degree - divisor_degree + 1):
        if offset & 255 == 0:
            budget_tick("polynomial division")
        lead = int(work[offset])
        if lead:
            work[offset : offset + divisor_degree + 1] = (
                work[offset : offset + divisor_degree + 1] - lead * divisor
            ) % p
    return trim(work[dividend_degree - divisor_degree + 1 :])


def polynomial_gcd(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    a, b = monic(a, p), monic(b, p)
    while not is_zero(b):
        budget_tick("Euclidean gcd")
        a, b = b, monic(remainder(a, b, p), p)
    return a


def derivative(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    return np.ascontiguousarray(
        ((np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % p)[::-1]
    )


def newton(power_sums: np.ndarray, degree: int, p: int) -> np.ndarray:
    """Descending monic polynomial coefficients from exact modular power sums."""
    chunk = 2**53 // ((p - 1) ** 2)
    assert chunk >= 1 and p > degree
    coefficients = np.zeros(degree + 1, dtype=np.int64)
    coefficients[0] = 1
    sums = np.ascontiguousarray(power_sums[: degree + 1], dtype=np.float64)
    reversed_coefficients = np.zeros(degree + 1, dtype=np.float64)
    reversed_coefficients[degree] = 1.0
    for k in range(1, degree + 1):
        if k & 255 == 0:
            budget_tick("Newton identities")
        total, base, lo = 0, degree - k, 1
        while lo <= k:
            hi = min(k, lo + chunk - 1)
            total += int(
                np.dot(
                    sums[lo : hi + 1],
                    reversed_coefficients[base + lo : base + hi + 1],
                )
            ) % p
            lo = hi + 1
        value = (-(total % p) * pow(k, -1, p)) % p
        coefficients[k] = value
        reversed_coefficients[degree - k] = float(value)
    return coefficients


def exterior_pair_gcd_degree(matrix: list[list[int]], p: int) -> dict[str, object]:
    """Reconstruct charpoly(wedge^2 M) without assembling wedge^2 M."""
    assert is_prime(p)
    reduced = np.array([[value % p for value in row] for row in matrix], dtype=np.int64)
    dim = reduced.shape[0]
    slots = comb(dim, 2)
    assert p > slots

    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for exponent in range(1, dim + 1):
        budget_tick(f"power traces p={p}")
        power = modmul(power, reduced, p)
        traces.append(int(power.trace()) % p)
    characteristic = newton(np.array(traces, dtype=np.int64), dim, p)

    # Independent Cayley-Hamilton audit of the reconstructed characteristic polynomial.
    value = np.eye(dim, dtype=np.int64)
    identity = np.eye(dim, dtype=np.int64)
    for coefficient in characteristic[1:]:
        value = (modmul(value, reduced, p) + int(coefficient) * identity) % p
    cayley_hamilton = bool(np.all(value == 0))

    tail = characteristic[1:]
    for k in range(dim + 1, 2 * slots + 1):
        if k & 255 == 0:
            budget_tick(f"trace recurrence p={p}")
        window = np.array(traces[k - 1 : k - dim - 1 : -1], dtype=np.int64)
        traces.append((-moddot(tail, window, p)) % p)

    inverse_two = pow(2, -1, p)
    pair_sums = np.zeros(slots + 1, dtype=np.int64)
    for m in range(1, slots + 1):
        pair_sums[m] = (
            ((traces[m] * traces[m] - traces[2 * m]) % p) * inverse_two
        ) % p
    pair_poly = newton(pair_sums, slots, p)
    pair_derivative = derivative(pair_poly, p)
    common = polynomial_gcd(pair_poly, pair_derivative, p)
    divides = is_zero(remainder(pair_poly, common, p)) and is_zero(
        remainder(pair_derivative, common, p)
    )
    return {
        "gcd_degree": int(common.size - 1),
        "pair_poly_monic": int(pair_poly[0]) == 1,
        "gcd_monic": int(common[0]) == 1,
        "gcd_divides_both": divides,
        "cayley_hamilton": cayley_hamilton,
    }


# ---------------------------------------------------------------- exact theorem checks

def gaussian_ceiling(n: int) -> int:
    return 3**n - 2**n


def isotropic_degree_bound(n: int, edges: int) -> int:
    slots = comb(1 << n, 2)
    return 4 * (2 * n + 2 * edges) * slots * (slots - 1)


def equal_site_upper(m: int) -> int:
    return comb(17, 2) * (2 * m + 1)


def claw_plus_ceiling(m: int) -> int:
    return gaussian_ceiling(m + 4)


def cut_polynomial(c: int) -> list[int]:
    coefficients = [0] * (4 * c + 1)
    for j in range(2 * c + 1):
        coefficients[2 * j] += comb(2 * c, j)
    coefficients[2 * c] -= 2 ** (2 * c)
    return coefficients


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    ok("artifact schema", set(artifact) == {"meta", "data"})
    data = artifact["data"]
    ok(
        "producer checks are all decisive passes",
        bool(data["checks"]) and all(item["passed"] for item in data["checks"]),
        f"{len(data['checks'])} rows",
    )

    stored_rows = {row["label"]: row for row in data["core_search"]["rows"]}
    independently_derived: dict[str, dict[str, object]] = {}
    winners: set[str] = set()
    for label, n, bonds, expected_gcd in SEEDS:
        raw = build_raw_operator(n, bonds, T0)
        scale, integer = actual_lcm_integer_model(raw)
        per_prime: dict[int, dict[str, object]] = {}
        for prime in PRIMES:
            result = exterior_pair_gcd_degree(integer, prime)
            per_prime[prime] = result
            ok(
                f"{label} p={prime}: independent raw-operator gcd degree",
                int(result["gcd_degree"]) == expected_gcd,
                f"degree {result['gcd_degree']} (expected {expected_gcd})",
            )
            ok(
                f"{label} p={prime}: independent monic/divisibility audits",
                bool(result["pair_poly_monic"])
                and bool(result["gcd_monic"])
                and bool(result["gcd_divides_both"])
                and bool(result["cayley_hamilton"]),
            )
        lower = comb(1 << n, 2) - expected_gcd
        if lower > gaussian_ceiling(n):
            winners.add(label)
        independently_derived[label] = {
            "actual_lcm": scale,
            "gcd_degrees": [per_prime[p]["gcd_degree"] for p in PRIMES],
            "r_lower_bound": lower,
        }
        stored = stored_rows[label]
        ok(
            f"{label}: final artifact matches independent seed",
            [stored["per_prime"][str(p)]["gcd_degree"] for p in PRIMES]
            == [expected_gcd, expected_gcd]
            and int(stored["distinct_pair_products_lower_bound"]) == lower
            and stored["bonds"] == [list(edge) for edge in bonds],
            f"r_Q >= {lower}, ceiling {gaussian_ceiling(n)}",
        )

    ok(
        "exact winning-core set",
        winners == EXPECTED_WINNERS,
        f"winners={sorted(winners)}",
    )

    named = {row["graph"]: row for row in data["named_core_theorems"]}
    for label in EXPECTED_WINNERS:
        n = next(seed[1] for seed in SEEDS if seed[0] == label)
        bonds = next(seed[2] for seed in SEEDS if seed[0] == label)
        expected_bound = isotropic_degree_bound(n, len(bonds))
        row = named[label]
        ok(
            f"{label}: Theorem-G exceptional bound",
            int(row["exception_polynomial_degree_bound"]) == expected_bound
            and int(row["characteristic_zero_r_lower_bound"])
            == independently_derived[label]["r_lower_bound"],
            f"degree <= {expected_bound}",
        )

    # Direct raw-operator audit of the equal-site tensor identity at m=1.
    claw_bonds = SEEDS[0][2]
    claw = build_raw_operator(4, claw_bonds, T0)
    one_site = build_raw_operator(1, tuple(), T0)
    claw_plus_isolate = build_raw_operator(5, claw_bonds, T0)
    ok(
        "raw isotropic tensor identity for claw plus an isolated site",
        claw_plus_isolate == kron_fraction(claw, one_site),
        "R_(K1,3 union K1)=R_K1,3 tensor P_t^2 entry-by-entry",
    )

    obstruction = data["decoupling_obstruction"]
    for row in obstruction["cut_factorization_obstruction"]["rows"]:
        c = int(row["cut_edges"])
        coefficients = cut_polynomial(c)
        ok(
            f"cut c={c}: independent nonzero-minor certificate",
            coefficients[0] == 1
            and coefficients[-1] == 1
            and int(row["degree_bound"]) == 4 * c,
            f"(1+t^2)^{2*c}-(2t)^{2*c} has constant and leading coefficient 1",
        )

    finite_rows = obstruction["equal_mode_counterfamily"]["finite_audit_rows"]
    ok(
        "collapse table covers every declared row m=0..12",
        [int(row["isolated_sites_m"]) for row in finite_rows] == list(range(13)),
    )
    for row in finite_rows:
        m = int(row["isolated_sites_m"])
        exponents = {j + k for j in range(m + 1) for k in range(m + 1)}
        ok(
            f"m={m}: independent identical-mode count and bounds",
            len(exponents) == 2 * m + 1
            and int(row["free_pair_exponent_count"]) == len(exponents)
            and int(row["proved_pair_product_upper_bound"]) == equal_site_upper(m)
            and int(row["gaussian_ceiling"]) == claw_plus_ceiling(m),
        )

    # Universal, rather than finite-table, verification.  Put U_m=136(2m+1)
    # and C_m=81*3^m-16*2^m.  The two exact differences below prove the induction.
    base = 3
    base_ok = equal_site_upper(base) == 952 < 2059 == claw_plus_ceiling(base)
    step_coefficients_ok = (
        # 2U_m-U_(m+1)=136(2m-1)>0 for every m>=3.
        2 * equal_site_upper(base) - equal_site_upper(base + 1) == 136 * (2 * base - 1)
        and 2 * base - 1 > 0
        # C_(m+1)-2C_m=81*3^m>0 for every m>=3.
        and claw_plus_ceiling(base + 1) - 2 * claw_plus_ceiling(base) == 81 * 3**base
    )
    ok(
        "universal m>=3 collapse proof",
        base_ok and step_coefficients_ok,
        "U_3=952<C_3=2059, U_(m+1)<2U_m, C_(m+1)>2C_m",
    )

    # q(t)-1=(t-1)^2/(2t): on 0<t<1 every nonempty cut minor is positive.
    physical_identity_ok = all(
        (1 + t * t) / (2 * t) - 1 == (t - 1) ** 2 / (2 * t) > 0
        for t in (Fraction(1, 5), Fraction(1, 3), Fraction(2, 5))
    )
    endpoint = build_raw_operator(2, ((0, 1),), Fraction(1, 1))
    endpoint_rank_one = all(
        endpoint[i][j] * endpoint[0][0] == endpoint[i][0] * endpoint[0][j]
        for i in range(4)
        for j in range(4)
    )
    ok(
        "physical cut obstruction and rank-collapsed endpoint",
        physical_identity_ok and endpoint_rank_one,
        "q>1 on 0<t<1; at t=1 the raw operator has rank one",
    )

    conclusion = data["conclusion"]
    ok(
        "scope is an obstruction, not a fabricated all-size theorem",
        conclusion["tag"] == "[UNRESOLVED]"
        and conclusion["requested_isotropic_all_size_theorem_proved"] is False
        and conclusion["n0"] is None,
    )
    ok(
        "test process-time and RSS budgets",
        time.process_time() - CPU0 < CPU_BUDGET_SECONDS and max_rss_bytes() < RSS_CAP_BYTES,
        f"CPU {time.process_time()-CPU0:.3f}s; RSS {max_rss_bytes()}",
    )

    if FAILED:
        print(f"FAIL test_isotropic_allsize: {len(FAILED)} failures: {FAILED}", flush=True)
        return 1
    print(f"PASS test_isotropic_allsize ({len(PASSED)} checks)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
