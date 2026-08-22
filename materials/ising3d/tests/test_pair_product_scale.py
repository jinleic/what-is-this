"""Standalone independent verifier for the scaled pair-product certificate.

This test deliberately does not import ``e118_pair_product_scale`` or its
helpers.  It re-derives the rational transfer matrix, Newton reconstruction,
and finite-field polynomial gcd from raw definitions.  It checks one new
coupling (2x3 at t=1/5), the wave-11 t=1/3 regression, and an open-chain
Gaussian control before comparing only recomputed fields with e118's JSON.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from math import gcd as math_gcd
from math import isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ising.transfer_matrix import layer_bonds  # noqa: E402

INT64_MAX = 2**63 - 1
P = 1_000_003
DIM6 = 64
PAIRS6 = DIM6 * (DIM6 - 1) // 2
GAUSSIAN_MAX6 = 3**6 - 2**6
GAUSSIAN_FLOOR6 = PAIRS6 - GAUSSIAN_MAX6

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)
    return bool(ok)


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def matmul_mod(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    assert a.shape[0] * (p - 1) ** 2 < INT64_MAX
    return (a @ b) % p


def power_trace_binary(matrix: np.ndarray, exponent: int, p: int) -> int:
    result = np.eye(matrix.shape[0], dtype=np.int64)
    base = matrix.copy()
    while exponent:
        if exponent & 1:
            result = matmul_mod(result, base, p)
        exponent >>= 1
        if exponent:
            base = matmul_mod(base, base, p)
    return int(result.trace()) % p


def newton(power_sums: list[int], degree: int, p: int) -> list[int]:
    """Return [1,a1,...,a_degree] for z^degree+a1*z^(degree-1)+... ."""
    coeffs = [1]
    for k in range(1, degree + 1):
        total = sum(coeffs[k - i] * power_sums[i] for i in range(1, k + 1))
        coeffs.append((-total * pow(k, -1, p)) % p)
    return coeffs


def extend_traces(traces: list[int], charpoly_desc: list[int], upto: int, p: int) -> list[int]:
    degree = len(charpoly_desc) - 1
    for k in range(degree + 1, upto + 1):
        total = sum(charpoly_desc[i] * traces[k - i] for i in range(1, degree + 1))
        traces.append((-total) % p)
    return traces


def trim(poly: list[int]) -> list[int]:
    i = 0
    while i + 1 < len(poly) and poly[i] == 0:
        i += 1
    return poly[i:] if any(poly[i:]) else [0]


def remainder(dividend: list[int], divisor: list[int], p: int) -> list[int]:
    divisor = trim(divisor)
    assert divisor != [0]
    inverse = pow(divisor[0] % p, -1, p)
    divisor = [coefficient * inverse % p for coefficient in divisor]
    work = [coefficient % p for coefficient in dividend]
    degree = len(divisor) - 1
    while len(work) - 1 >= degree and any(work):
        lead = work[0]
        if lead:
            for j, coefficient in enumerate(divisor):
                work[j] = (work[j] - lead * coefficient) % p
        work = trim(work[1:])
    return work


def polynomial_gcd(a: list[int], b: list[int], p: int) -> list[int]:
    def monic(poly: list[int]) -> list[int]:
        poly = trim(poly)
        return [coefficient * pow(poly[0], -1, p) % p for coefficient in poly] if any(poly) else [0]

    a, b = monic(a), monic(b)
    while any(b):
        a, b = b, monic(remainder(a, b, p))
    return a


def raw_build_r(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Re-derive P_t diag(q^e) P_t without importing e38/e118."""
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    energies = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    epsilon = energies[0] % 2
    assert all((energy - epsilon) % 2 == 0 for energy in energies)
    weights = [q ** ((energy - epsilon) // 2) for energy in energies]
    powers = [t**hamming for hamming in range(n + 1)]
    kernel = [[powers[(left ^ right).bit_count()] for right in range(dim)] for left in range(dim)]
    matrix = [[Fraction(0) for _ in range(dim)] for _ in range(dim)]
    for left in range(dim):
        for right in range(left, dim):
            entry = sum(kernel[left][middle] * weights[middle] * kernel[middle][right] for middle in range(dim))
            matrix[left][right] = matrix[right][left] = entry
    return matrix, q, epsilon


def integralize(matrix: list[list[Fraction]]) -> tuple[int, list[list[int]]]:
    scale = 1
    for row in matrix:
        for entry in row:
            scale = scale * entry.denominator // math_gcd(scale, entry.denominator)
    output = []
    for row in matrix:
        integer_row = []
        for entry in row:
            scaled = scale * entry
            assert scaled.denominator == 1
            integer_row.append(scaled.numerator)
        output.append(integer_row)
    return scale, output


def pipeline(matrix_rows: list[list[int]], p: int) -> dict[str, object]:
    dim = len(matrix_rows)
    pairs = dim * (dim - 1) // 2
    assert is_prime(p) and p > pairs and p % 2 == 1
    matrix = np.array([[entry % p for entry in row] for row in matrix_rows], dtype=np.int64)
    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for _ in range(dim):
        power = matmul_mod(power, matrix, p)
        traces.append(int(power.trace()) % p)
    characteristic = newton(traces, dim, p)
    extended = extend_traces(traces, characteristic, 2 * pairs, p)
    pair_sums = [0] + [
        ((extended[m] * extended[m] - extended[2 * m]) % p) * pow(2, -1, p) % p
        for m in range(1, pairs + 1)
    ]
    pair_poly = newton(pair_sums, pairs, p)
    derivative = [((pairs - index) * pair_poly[index]) % p for index in range(pairs)]
    common = polynomial_gcd(pair_poly, derivative, p)
    return {
        "matrix": matrix,
        "traces": extended,
        "charpoly": characteristic,
        "pair_poly": pair_poly,
        "derivative": derivative,
        "gcd": common,
        "gcd_degree": len(trim(common)) - 1,
        "charpoly_digest": digest(characteristic[::-1]),
        "pair_poly_digest": digest(pair_poly[::-1]),
    }


# The 2x4 audit is deliberately a second, NumPy-vectorized implementation:
# it is kept in this standalone test rather than importing e118's producer
# helpers.  Pure-Python Newton summation is appropriate for degree 2016 but
# intentionally impractical for degree 32640.
def np_trim(poly: np.ndarray) -> np.ndarray:
    indices = np.flatnonzero(poly)
    return poly[indices[0] :] if indices.size else np.zeros(1, dtype=np.int64)


def np_zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def np_monic(poly: np.ndarray, p: int) -> np.ndarray:
    poly = np_trim(np.asarray(poly, dtype=np.int64) % p)
    return poly if np_zero(poly) else (poly * pow(int(poly[0]), -1, p)) % p


def np_remainder(dividend: np.ndarray, divisor: np.ndarray, p: int) -> np.ndarray:
    assert int(divisor[0]) == 1
    work = (np.asarray(dividend, dtype=np.int64) % p).copy()
    divisor_degree = divisor.size - 1
    dividend_degree = work.size - 1
    if dividend_degree < divisor_degree:
        return np_trim(work)
    for offset in range(dividend_degree - divisor_degree + 1):
        lead = int(work[offset])
        if lead:
            work[offset : offset + divisor_degree + 1] = (
                work[offset : offset + divisor_degree + 1] - lead * divisor
            ) % p
    return np_trim(work[dividend_degree - divisor_degree + 1 :])


def np_gcd(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    a, b = np_monic(a, p), np_monic(b, p)
    while not np_zero(b):
        a, b = b, np_monic(np_remainder(a, b, p), p)
    return a


def np_newton(power_sums: np.ndarray, degree: int, p: int) -> np.ndarray:
    assert degree * (p - 1) ** 2 < INT64_MAX
    coefficients = np.zeros(degree + 1, dtype=np.int64)
    coefficients[0] = 1
    for k in range(1, degree + 1):
        total = int(np.dot(power_sums[1 : k + 1], np.ascontiguousarray(coefficients[k - 1 :: -1]))) % p
        coefficients[k] = (-total * pow(k, -1, p)) % p
    return coefficients


def np_extend_traces(traces: list[int], charpoly_desc: np.ndarray, upto: int, p: int) -> list[int]:
    dim = charpoly_desc.size - 1
    coefficients = charpoly_desc[1:]
    for k in range(dim + 1, upto + 1):
        window = np.array(traces[k - 1 : k - dim - 1 : -1], dtype=np.int64)
        traces.append((-int(np.dot(coefficients, window))) % p)
    return traces


def np_horner_zero(charpoly_desc: np.ndarray, matrix: np.ndarray, p: int) -> bool:
    identity = np.eye(matrix.shape[0], dtype=np.int64)
    value = int(charpoly_desc[0]) * identity % p
    for coefficient in charpoly_desc[1:]:
        value = (matmul_mod(value, matrix, p) + int(coefficient) * identity) % p
    return bool(np.all(value == 0))


def pipeline_large(matrix_rows: list[list[int]], p: int) -> dict[str, object]:
    """Independent trace/Newton/recurrence/gcd path for the 32640-slot case."""
    dim = len(matrix_rows)
    pairs = dim * (dim - 1) // 2
    assert is_prime(p) and p > pairs and p % 2 == 1
    assert pairs * (p - 1) ** 2 < INT64_MAX
    matrix = np.array([[entry % p for entry in row] for row in matrix_rows], dtype=np.int64)
    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for _ in range(dim):
        power = matmul_mod(power, matrix, p)
        traces.append(int(power.trace()) % p)
    characteristic = np_newton(np.array(traces, dtype=np.int64), dim, p)
    horner_zero = np_horner_zero(characteristic, matrix, p)
    traces = np_extend_traces(traces, characteristic, 2 * pairs, p)
    inverse_two = pow(2, -1, p)
    pair_sums = np.zeros(pairs + 1, dtype=np.int64)
    for exponent in range(1, pairs + 1):
        pair_sums[exponent] = (
            (int(traces[exponent]) * int(traces[exponent]) - int(traces[2 * exponent])) % p
        ) * inverse_two % p
    pair_poly = np_newton(pair_sums, pairs, p)
    derivative = np.ascontiguousarray(
        (np.arange(pairs, 0, -1, dtype=np.int64) * pair_poly[:-1]) % p
    )
    common = np_gcd(pair_poly, derivative, p)
    recurrence_spots = {
        exponent: power_trace_binary(matrix, exponent, p) == traces[exponent]
        for exponent in (dim + 1, pairs, 2 * pairs)
    }
    return {
        "charpoly": characteristic,
        "pair_poly": pair_poly,
        "gcd": common,
        "gcd_degree": int(common.size) - 1,
        "charpoly_digest": digest([int(value) for value in characteristic[::-1] % p]),
        "pair_poly_digest": digest([int(value) for value in pair_poly[::-1] % p]),
        "horner_zero": horner_zero,
        "gcd_divides_both": np_zero(np_remainder(pair_poly, common, p))
        and np_zero(np_remainder(derivative, common, p)),
        "recurrence_spots": recurrence_spots,
    }


def digest(coefficients: list[int]) -> str:
    return hashlib.sha256(",".join(str(int(coefficient)) for coefficient in coefficients).encode()).hexdigest()


def main() -> int:
    artifact_path = ROOT / "results" / "spectral" / "pair_product_scale.json"
    # This is intentionally first: before e118 exists the contract is absent, so the
    # new test must fail rather than silently validate only the old wave-11 result.
    artifact = json.load(artifact_path.open())

    print(f"decisive prime p = {P}")
    check("p is prime, odd, and exceeds C(64,2)", is_prime(P) and P % 2 and P > PAIRS6)
    check("2x3 threshold C(64,2)=2016", PAIRS6 == 2016)
    check("2x3 Gaussian maximum 3^6-2^6=665", GAUSSIAN_MAX6 == 665)
    check("2x3 Gaussian gcd floor=1351", GAUSSIAN_FLOOR6 == 1351)

    # Catch a Newton-sign defect independently of transfer-matrix data.
    cubic = newton([0, 10, 38, 160], 3, P)
    check("Newton signs recover (z-2)(z-3)(z-5)", cubic == [1, -10 % P, 31 % P, -30 % P])

    bonds23 = [(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (4, 5)]
    check("hard-coded open 2x3 bonds match repository geometry", set(bonds23) == set(layer_bonds((2, 3), (False, False))))

    print("\nnew decisive coupling: open 2x3 at t=1/5")
    new_matrix, new_q, new_epsilon = raw_build_r(6, bonds23, Fraction(1, 5))
    new_scale, new_integer = integralize(new_matrix)
    new_result = pipeline(new_integer, P)
    check("new coupling has q=13/5 and odd bond parity", new_q == Fraction(13, 5) and new_epsilon == 1)
    check("new C_2 is monic of degree 2016", new_result["pair_poly"][0] == 1 and len(new_result["pair_poly"]) == 2017)
    check("new gcd divides C_2 and its formal derivative", remainder(new_result["pair_poly"], new_result["gcd"], P) == [0] and remainder(new_result["derivative"], new_result["gcd"], P) == [0])
    new_gcd = int(new_result["gcd_degree"])
    check("new t=1/5 layer is decisively below the Gaussian floor", new_gcd < GAUSSIAN_FLOOR6, f"gcd degree {new_gcd}")

    print("\nwave-11 regression: open 2x3 at t=1/3")
    old_matrix, old_q, old_epsilon = raw_build_r(6, bonds23, Fraction(1, 3))
    old_scale, old_integer = integralize(old_matrix)
    old_result = pipeline(old_integer, P)
    old_gcd = int(old_result["gcd_degree"])
    check("wave-11 input has q=5/3, parity 1, and its recorded exact scale", old_q == Fraction(5, 3) and old_epsilon == 1 and old_scale == 3**15 * 5**4)
    check("wave-11 modular gcd regression is 385", old_gcd == 385, f"got {old_gcd}")

    print("\nopen-chain n=6 Gaussian control at t=1/5")
    chain_matrix, _, _ = raw_build_r(6, [(site, site + 1) for site in range(5)], Fraction(1, 5))
    chain_scale, chain_integer = integralize(chain_matrix)
    chain_result = pipeline(chain_integer, P)
    chain_gcd = int(chain_result["gcd_degree"])
    check("same-site-count chain reaches exactly the Gaussian gcd floor", chain_gcd == GAUSSIAN_FLOOR6, f"got {chain_gcd}, D={chain_scale}")

    # Direct binary powers test the recurrence used in the pair formula at widely
    # separated exponents, including a scope-sensitive endpoint 2*C(64,2).
    recurrence_ok = all(
        power_trace_binary(new_result["matrix"], k, P) == new_result["traces"][k]
        for k in (65, 127, 2016, 4032)
    )
    check("new-case trace recurrence agrees with binary powering through 4032", recurrence_ok)

    print("\nnew decisive larger layer: open 2x4 at t=1/3")
    bonds24 = [
        (0, 4), (0, 1), (1, 5), (1, 2), (2, 6),
        (2, 3), (3, 7), (4, 5), (5, 6), (6, 7),
    ]
    check(
        "hard-coded open 2x4 bonds match repository geometry",
        set(bonds24) == set(layer_bonds((2, 4), (False, False))),
    )
    layer24_matrix, layer24_q, layer24_epsilon = raw_build_r(8, bonds24, Fraction(1, 3))
    layer24_scale, layer24_integer = integralize(layer24_matrix)
    del layer24_matrix
    layer24_result = pipeline_large(layer24_integer, P)
    pairs24 = 256 * 255 // 2
    gaussian_max24 = 3**8 - 2**8
    gaussian_floor24 = pairs24 - gaussian_max24
    layer24_gcd = int(layer24_result["gcd_degree"])
    layer24_distinct_lower = pairs24 - layer24_gcd
    check(
        "2x4 raw input has q=5/3, even parity, and actual scale 3^21*5^5",
        layer24_q == Fraction(5, 3)
        and layer24_epsilon == 0
        and layer24_scale == 3**21 * 5**5,
    )
    check(
        "2x4 defining threshold identities: C(256,2)=32640, max=6305, floor=26335",
        pairs24 == 32640 and gaussian_max24 == 6305 and gaussian_floor24 == 26335,
    )
    check(
        "independent 2x4 trace reconstruction has Cayley-Hamilton, recurrence, and gcd-division checks",
        bool(layer24_result["horner_zero"])
        and all(layer24_result["recurrence_spots"].values())
        and bool(layer24_result["gcd_divides_both"]),
    )
    check(
        "independent 2x4 modular gcd is 9329, giving 23311 > 6305 pair products",
        layer24_gcd == 9329
        and layer24_distinct_lower == 23311
        and layer24_distinct_lower > gaussian_max24,
        f"gcd={layer24_gcd}, lower={layer24_distinct_lower}",
    )

    print("\nopen-chain n=8 Gaussian control at t=1/3")
    chain8_matrix, chain8_q, chain8_epsilon = raw_build_r(
        8, [(site, site + 1) for site in range(7)], Fraction(1, 3)
    )
    chain8_scale, chain8_integer = integralize(chain8_matrix)
    del chain8_matrix
    chain8_result = pipeline_large(chain8_integer, P)
    chain8_gcd = int(chain8_result["gcd_degree"])
    check(
        "independent n=8 chain has q=5/3, odd parity, and exact Gaussian gcd floor",
        chain8_q == Fraction(5, 3)
        and chain8_epsilon == 1
        and chain8_gcd == gaussian_floor24,
        f"gcd={chain8_gcd}, D={chain8_scale}",
    )
    check(
        "independent n=8 chain recurrence and gcd checks hold",
        bool(chain8_result["horner_zero"])
        and all(chain8_result["recurrence_spots"].values())
        and bool(chain8_result["gcd_divides_both"]),
    )

    data = artifact["data"]
    cases = data["cases"]
    artifact_new = cases["layer_2x3_t_1_5"]
    artifact_old = cases["layer_2x3_t_1_3"]
    artifact_chain = cases["chain_n6_t_1_5"]
    check("artifact records the independently rebuilt new gcd degree", artifact_new["per_prime"][str(P)]["gcd_degree"] == new_gcd)
    check("artifact records the independently rebuilt new charpoly digest", artifact_new["per_prime"][str(P)]["charpoly_sha256"] == new_result["charpoly_digest"])
    check("artifact records the independently rebuilt new pair-poly digest", artifact_new["per_prime"][str(P)]["pair_poly_sha256"] == new_result["pair_poly_digest"])
    check("artifact preserves the wave-11 t=1/3 regression", artifact_old["per_prime"][str(P)]["gcd_degree"] == old_gcd)
    check("artifact preserves the raw chain-control value", artifact_chain["per_prime"][str(P)]["gcd_degree"] == chain_gcd)

    # The all-but-finitely-many result has separate algebraic content: its
    # numerical degree accounting must not silently degrade into a grid claim.
    generic = data["generic_coupling_theorem_2x3"]
    expected_resultant_bound = 2 * 52 * PAIRS6 * (PAIRS6 - 1)
    check(
        "generic polynomial representative uses degree-26 entries and degree-52 wedge entries",
        generic["polynomial_scalar"] == "(2t)^7 q(t)^4 = (2t)^3(1+t^2)^4"
        and generic["matrix_entry_t_degree_upper_bound"] == 26
        and generic["wedge2_entry_t_degree_upper_bound"] == 52,
    )
    check(
        "generic specialization bound is independently tied to the raw t=1/3 regression",
        generic["generic_gcd_degree_upper_bound"] == old_gcd == 385,
    )
    check(
        "generic resultant exceptional-root bound is 2*52*2016*2015",
        generic["resultant_exceptional_root_bound"] == expected_resultant_bound,
        f"got {generic['resultant_exceptional_root_bound']}",
    )
    check(
        "generic lower bound is 2016-385=1631 > 665",
        generic["generic_distinct_pair_products_lower_bound"] == PAIRS6 - old_gcd
        and generic["generic_distinct_pair_products_lower_bound"] > GAUSSIAN_MAX6,
    )
    grid = data["grid_cross_check"]["per_layer_case_modular_gcd_degrees"]
    check(
        "every stored deterministic 2x3 grid point has the recorded two-prime degree 385",
        all(values == [385, 385] for values in grid.values()) and len(grid) == 8,
    )

    layer24 = cases["layer_2x4_t_1_3"]
    chain8 = cases["chain_n8_t_1_3"]
    layer24_p = layer24["per_prime"][str(P)]
    chain8_p = chain8["per_prime"][str(P)]
    check(
        "artifact 2x4 threshold fields equal the independently recomputed identities",
        layer24["dim"] == 256
        and layer24["pair_count"] == pairs24
        and layer24["gaussian_max_distinct_pair_products"] == gaussian_max24
        and layer24["gaussian_gcd_floor"] == gaussian_floor24
        and int(layer24["D"]) == layer24_scale,
    )
    check(
        "artifact 2x4 p=1000003 gcd and both polynomial digests match the independent raw rebuild",
        layer24_p["gcd_degree"] == layer24_gcd
        and layer24_p["charpoly_sha256"] == layer24_result["charpoly_digest"]
        and layer24_p["pair_poly_sha256"] == layer24_result["pair_poly_digest"],
    )
    check(
        "artifact 2x4 inference fields implement the defining decisive comparison",
        layer24["modular_gcd_degree_upper_bound"] == layer24_gcd
        and layer24["distinct_pair_products_lower_bound"] == layer24_distinct_lower
        and layer24["decisive_non_gaussian"]
        == (layer24["modular_gcd_degree_upper_bound"] < layer24["gaussian_gcd_floor"])
        and layer24_distinct_lower > gaussian_max24,
    )
    check(
        "artifact 2x4 second-prime cross-check is also gcd degree 9329",
        layer24["per_prime"]["2000003"]["gcd_degree"] == 9329,
    )
    check(
        "artifact n=8 chain p=1000003 gcd and both polynomial digests match the independent raw rebuild",
        chain8_p["gcd_degree"] == chain8_gcd
        and chain8_p["charpoly_sha256"] == chain8_result["charpoly_digest"]
        and chain8_p["pair_poly_sha256"] == chain8_result["pair_poly_digest"]
        and int(chain8["D"]) == chain8_scale,
    )
    check(
        "same-site-count n=8 chain control equals its exact Gaussian threshold at both stored primes",
        chain8["modular_gcd_degree_upper_bound"] == gaussian_floor24
        and all(record["gcd_degree"] == gaussian_floor24 for record in chain8["per_prime"].values()),
    )
    check("all stored producer checks passed", all(record["passed"] for record in artifact["checks"]))

    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(
        "OK: raw t=1/5 rebuild gives gcd degree "
        f"{new_gcd}; wave-11 t=1/3 is {old_gcd}; chain control is {chain_gcd}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
