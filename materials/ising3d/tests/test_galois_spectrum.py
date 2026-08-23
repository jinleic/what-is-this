"""[COMPUTATION] Standalone clean-room replay of the Galois-spectrum certificate.

This verifier imports none of e203/e204/e205 and does not import build_R.  It
constructs an integer transfer representative directly from local tensor-entry
integers, recomputes every characteristic factorization, and audits the explicit
finite-field factors with coefficient-list arithmetic and Rabin's irreducibility
criterion.

Run from the repository root with:
    .venv/bin/python tests/test_galois_spectrum.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Iterable

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "galois_spectrum.json"
PROOF_NOTE = ROOT / "proofs" / "galois_spectrum.md"
X = sp.Symbol("x")
T_VALUE = Fraction(1, 3)
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2_000_000_000
C4_TO_SQUARE_SITE_MAP = (0, 1, 3, 2)
EXPECTED_PRIMES = {
    "ising_2d_open_width6": {
        "irreducible_sextic": 7,
        "five_cycle": 43,
        "transposition": 587,
    },
    "ising_3d_open_2x3": {
        "irreducible_sextic": 13,
        "five_cycle": 53,
        "transposition": 197,
    },
}
EXPECTED_CYCLES = {
    "irreducible_sextic": [6],
    "five_cycle": [1, 5],
    "transposition": [1, 1, 1, 1, 2],
}
CASES = (
    {
        "key": "ising_1d_single_spin",
        "n_sites": 1,
        "bonds": (),
    },
    {
        "key": "ising_2d_periodic_width4",
        "n_sites": 4,
        "bonds": ((0, 1), (1, 2), (2, 3), (3, 0)),
    },
    {
        "key": "ising_3d_open_2x2",
        "n_sites": 4,
        "bonds": ((0, 2), (0, 1), (1, 3), (2, 3)),
    },
    {
        "key": "ising_2d_open_width6",
        "n_sites": 6,
        "bonds": ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5)),
    },
    {
        "key": "ising_3d_open_2x3",
        "n_sites": 6,
        "bonds": ((0, 3), (0, 1), (1, 4), (1, 2), (2, 5), (3, 4), (4, 5)),
    },
)

PASSED: list[str] = []
FAILED: list[str] = []


class NonDecisiveBudget(RuntimeError):
    pass


class CpuBudget:
    def __init__(self, limit: float = CPU_BUDGET_SECONDS) -> None:
        self.limit = limit
        self.started = time.process_time()

    def check(self, stage: str) -> None:
        used = time.process_time() - self.started
        if used > self.limit:
            raise NonDecisiveBudget(
                f"[UNRESOLVED] process-CPU budget {self.limit:.1f}s expired "
                f"at {stage} after {used:.3f}s"
            )
        rss = max_rss_bytes()
        if rss > RSS_LIMIT_BYTES:
            raise NonDecisiveBudget(
                f"[UNRESOLVED] RSS limit {RSS_LIMIT_BYTES} exceeded at {stage}; "
                f"ru_maxrss={rss}"
            )


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def ok(name: str, condition: bool, detail: str = "") -> None:
    print(
        f"[{'PASS' if condition else 'FAIL'}] {name}"
        + (f": {detail}" if detail else ""),
        flush=True,
    )
    (PASSED if condition else FAILED).append(name)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_digest(rows: Iterable[Iterable[int]]) -> str:
    payload = ";".join(",".join(str(int(value)) for value in row) for row in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def polynomial_digest(poly: sp.Poly) -> str:
    exact = sp.Poly(poly, X, domain=sp.QQ).monic()
    payload = "\n".join(str(coefficient) for coefficient in exact.all_coeffs())
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


# ---------------------------------------------------------------- direct integer transfer construction


def direct_primitive_operator(
    n_sites: int,
    bonds: tuple[tuple[int, int], ...],
    t_value: Fraction,
    budget: CpuBudget,
) -> tuple[list[list[int]], dict[str, int | str]]:
    """Build an integer scalar multiple without Fraction matrix multiplication.

    If t=a/b and q=c/d, use
        P_int[i,j] = a^H(i,j) b^(n-H(i,j)),
        D_int[s]   = c^aligned(s) d^(m-aligned(s)).
    Then A=P_int D_int P_int is an integer scalar multiple of the canonical R.
    Dividing A by its entry content gives the unique positive primitive integer
    representative, independently of denominator-lcm clearing.
    """

    a, b = t_value.numerator, t_value.denominator
    q_value = (1 + t_value * t_value) / (2 * t_value)
    c, d = q_value.numerator, q_value.denominator
    dimension = 1 << n_sites
    edge_count = len(bonds)

    p_integer = [[0] * dimension for _ in range(dimension)]
    for left in range(dimension):
        if left % 16 == 0:
            budget.check(f"direct P row {left}/{dimension}")
        for right in range(dimension):
            distance = (left ^ right).bit_count()
            p_integer[left][right] = a**distance * b ** (n_sites - distance)

    diagonal: list[int] = []
    for state in range(dimension):
        spins = [
            1 - 2 * ((state >> (n_sites - 1 - site)) & 1)
            for site in range(n_sites)
        ]
        aligned = sum(spins[left] == spins[right] for left, right in bonds)
        diagonal.append(c**aligned * d ** (edge_count - aligned))

    integer = [[0] * dimension for _ in range(dimension)]
    for left in range(dimension):
        budget.check(f"direct integer PDP row {left}/{dimension}")
        p_left = p_integer[left]
        for right in range(left, dimension):
            value = sum(
                p_left[middle] * diagonal[middle] * p_integer[middle][right]
                for middle in range(dimension)
            )
            integer[left][right] = integer[right][left] = value

    content = math.gcd(*(abs(value) for row in integer for value in row))
    primitive = [[value // content for value in row] for row in integer]
    assert math.gcd(*(abs(value) for row in primitive for value in row)) == 1
    half_edges = (edge_count + 1) // 2
    natural_scale = b ** (2 * n_sites) * c**half_edges * d ** (edge_count - half_edges)
    return primitive, {
        "q": fraction_text(q_value),
        "natural_denominator_scale": natural_scale,
        "entry_content": content,
        "primitive_scale_B_over_R": fraction_text(Fraction(natural_scale, content)),
    }


def exact_factorization(
    matrix: list[list[int]], budget: CpuBudget, stage: str
) -> tuple[sp.Poly, list[tuple[sp.Poly, int]]]:
    budget.check(f"before clean-room charpoly {stage}")
    characteristic = sp.Poly(sp.Matrix(matrix).charpoly(X).as_expr(), X, domain=sp.ZZ)
    content, raw_factors = sp.factor_list(characteristic)
    assert content == 1
    factors = [
        (sp.Poly(factor, X, domain=sp.ZZ), int(exponent))
        for factor, exponent in raw_factors
    ]
    factors.sort(
        key=lambda item: (
            item[0].degree(),
            tuple(int(value) for value in item[0].all_coeffs()),
            item[1],
        )
    )
    reconstruction = sp.Poly(1, X, domain=sp.ZZ)
    for factor, exponent in factors:
        assert factor.LC() == 1 and factor.is_irreducible
        reconstruction *= factor**exponent
    assert reconstruction == characteristic
    budget.check(f"after clean-room factorization {stage}")
    return characteristic, factors


def compare_factor_record(
    artifact_case: dict[str, object],
    characteristic: sp.Poly,
    factors: list[tuple[sp.Poly, int]],
) -> bool:
    stored = artifact_case["characteristic_polynomial"]
    if int(stored["degree"]) != characteristic.degree():
        return False
    if stored["sha256"] != polynomial_digest(characteristic):
        return False
    rows = stored["factors"]
    if len(rows) != len(factors):
        return False
    for index, ((factor, exponent), row) in enumerate(zip(factors, rows)):
        coefficients = [str(value) for value in factor.all_coeffs()]
        if not (
            int(row["index"]) == index
            and int(row["degree"]) == factor.degree()
            and int(row["exponent"]) == exponent
            and list(row["coefficients_desc"]) == coefficients
            and row["sha256"] == polynomial_digest(factor)
            and row["irreducible_over_Q"] is True
        ):
            return False
    return bool(stored["factorization_reconstructs"]) and bool(
        stored["all_factors_irreducible_over_Q"]
    )


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def independent_root_scale(factor: sp.Poly) -> int:
    coefficients = [int(value) for value in factor.all_coeffs()]
    first_index = next(index for index, value in enumerate(coefficients[1:], 1) if value)
    first = abs(coefficients[first_index])
    scale = 1
    for prime_value, first_power in sp.factorint(first).items():
        prime = int(prime_value)
        bound = int(first_power) // first_index
        for index, coefficient in enumerate(coefficients[1:], 1):
            if coefficient:
                bound = min(bound, valuation(abs(coefficient), prime) // index)
        scale *= prime**bound
    assert all(
        coefficient % scale**index == 0
        for index, coefficient in enumerate(coefficients)
    )
    return scale


# ------------------------------------------------ independent finite-field polynomial arithmetic


def p_trim(poly: list[int], prime: int) -> list[int]:
    result = [value % prime for value in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def p_add(left: list[int], right: list[int], prime: int) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        ) % prime
    return p_trim(result, prime)


def p_sub(left: list[int], right: list[int], prime: int) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (
            (left[index] if index < len(left) else 0)
            - (right[index] if index < len(right) else 0)
        ) % prime
    return p_trim(result, prime)


def p_mul(left: list[int], right: list[int], prime: int) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = (
                result[left_index + right_index] + left_value * right_value
            ) % prime
    return p_trim(result, prime)


def p_divmod(
    dividend: list[int], divisor: list[int], prime: int
) -> tuple[list[int], list[int]]:
    remainder = p_trim(dividend, prime)
    divisor = p_trim(divisor, prime)
    if divisor == [0]:
        raise ZeroDivisionError("zero finite-field polynomial")
    if len(remainder) < len(divisor):
        return [0], remainder
    quotient = [0] * (len(remainder) - len(divisor) + 1)
    inverse_lead = pow(divisor[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse_lead % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] = (
                remainder[index + shift] - coefficient * value
            ) % prime
        remainder = p_trim(remainder, prime)
    return p_trim(quotient, prime), remainder


def p_mod(poly: list[int], modulus: list[int], prime: int) -> list[int]:
    return p_divmod(poly, modulus, prime)[1]


def p_monic(poly: list[int], prime: int) -> list[int]:
    poly = p_trim(poly, prime)
    inverse = pow(poly[-1], -1, prime)
    return [(value * inverse) % prime for value in poly]


def p_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    left = p_trim(left, prime)
    right = p_trim(right, prime)
    while right != [0]:
        left, right = right, p_mod(left, right, prime)
    return p_monic(left, prime)


def p_powmod(
    base: list[int], exponent: int, modulus: list[int], prime: int
) -> list[int]:
    result = [1]
    base = p_mod(base, modulus, prime)
    while exponent:
        if exponent & 1:
            result = p_mod(p_mul(result, base, prime), modulus, prime)
        base = p_mod(p_mul(base, base, prime), modulus, prime)
        exponent >>= 1
    return result


def integer_prime_divisors(value: int) -> list[int]:
    divisors: list[int] = []
    candidate = 2
    while candidate * candidate <= value:
        if value % candidate == 0:
            divisors.append(candidate)
            while value % candidate == 0:
                value //= candidate
        candidate += 1
    if value > 1:
        divisors.append(value)
    return divisors


def irreducible_over_fp(coefficients_desc: list[int], prime: int) -> bool:
    factor = p_trim(list(reversed(coefficients_desc)), prime)
    degree = len(factor) - 1
    if degree <= 0 or factor[-1] != 1:
        return False
    if degree == 1:
        return True
    x_poly = [0, 1]
    if p_powmod(x_poly, prime**degree, factor, prime) != p_mod(
        x_poly, factor, prime
    ):
        return False
    for divisor in integer_prime_divisors(degree):
        power = prime ** (degree // divisor)
        difference = p_sub(
            p_powmod(x_poly, power, factor, prime), x_poly, prime
        )
        if len(p_gcd(factor, difference, prime)) > 1:
            return False
    return True


def small_prime(value: int) -> bool:
    if value < 2:
        return False
    candidate = 2
    while candidate * candidate <= value:
        if value % candidate == 0:
            return value == candidate
        candidate += 1
    return True


def audit_modular_witness(
    normalized_coefficients_desc: list[int],
    discriminant: int,
    stored: dict[str, object],
) -> tuple[bool, list[int]]:
    prime = int(stored["prime"])
    if not small_prime(prime) or discriminant % prime == 0:
        return False, []
    reduced_desc = [value % prime for value in normalized_coefficients_desc]
    if list(stored["polynomial_coefficients_desc_mod_p"]) != reduced_desc:
        return False, []

    reconstructed = [1]
    cycle_type: list[int] = []
    all_irreducible = True
    for row in stored["factors"]:
        coefficients_desc = [int(value) % prime for value in row["coefficients_desc_mod_p"]]
        exponent = int(row["exponent"])
        irreducible = irreducible_over_fp(coefficients_desc, prime)
        all_irreducible = all_irreducible and irreducible and row["irreducible_over_Fp"] is True
        factor_ascending = list(reversed(coefficients_desc))
        for _ in range(exponent):
            reconstructed = p_mul(reconstructed, factor_ascending, prime)
            cycle_type.append(len(coefficients_desc) - 1)
    target_ascending = p_trim(list(reversed(reduced_desc)), prime)
    reconstructed = p_trim(reconstructed, prime)
    derivative = [
        index * target_ascending[index] % prime
        for index in range(1, len(target_ascending))
    ] or [0]
    squarefree = len(p_gcd(target_ascending, derivative, prime)) == 1
    cycle_type.sort()
    conditions = (
        reconstructed == target_ascending
        and all_irreducible
        and squarefree
        and int(stored["discriminant_mod_prime"]) == discriminant % prime
        and stored["prime_verified"] is True
        and stored["good_prime"] is True
        and stored["factor_product_reconstructs"] is True
        and stored["squarefree_mod_prime"] is True
        and stored["all_displayed_factors_irreducible_over_Fp"] is True
        and list(stored["frobenius_cycle_type"]) == cycle_type
    )
    return conditions, cycle_type


def main() -> None:
    budget = CpuBudget()
    artifact = json.loads(ARTIFACT.read_text())
    ok(
        "artifact top-level shape",
        set(artifact) == {"meta", "data", "checks"},
        str(sorted(artifact)),
    )
    ok(
        "producer stored checks pass",
        bool(artifact["checks"])
        and all(row["passed"] is True for row in artifact["checks"]),
        f"{len(artifact['checks'])} checks",
    )

    source_hashes = artifact["meta"]["source_sha256"]
    sources_match = all(
        sha256_file(ROOT / relative_path) == expected_hash
        for relative_path, expected_hash in source_hashes.items()
    )
    ok("producer source hashes", sources_match)
    ok("proof note exists", PROOF_NOTE.is_file())

    stored_cases = {
        row["key"]: row for row in artifact["data"]["exact_cases"]
    }
    rebuilt: dict[str, dict[str, object]] = {}
    for spec in CASES:
        key = str(spec["key"])
        budget.check(f"start clean-room case {key}")
        matrix, construction = direct_primitive_operator(
            int(spec["n_sites"]), spec["bonds"], T_VALUE, budget
        )
        characteristic, factors = exact_factorization(matrix, budget, key)
        stored = stored_cases[key]
        construction_matches = (
            stored["q"] == construction["q"]
            and int(stored["entry_denominator_lcm"])
            == int(construction["natural_denominator_scale"])
            and int(stored["cleared_matrix_entry_content"])
            == int(construction["entry_content"])
            and stored["primitive_scale_B_over_R"]
            == construction["primitive_scale_B_over_R"]
            and stored["primitive_matrix_sha256"] == matrix_digest(matrix)
            and stored["primitive_matrix_symmetric"] is True
        )
        factor_matches = compare_factor_record(stored, characteristic, factors)
        ok(
            f"clean-room matrix and exact factors for {key}",
            construction_matches and factor_matches,
            str(stored["characteristic_polynomial"]["factor_degree_exponent_pattern"]),
        )
        rebuilt[key] = {
            "matrix": matrix,
            "characteristic": characteristic,
            "factors": factors,
        }

    one_factors = [
        [int(value) for value in factor.all_coeffs()]
        for factor, _exponent in rebuilt["ising_1d_single_spin"]["factors"]
    ]
    ok("one-dimensional factors", one_factors == [[1, -8], [1, -2]])

    c4_factors = rebuilt["ising_2d_periodic_width4"]["factors"]
    square_factors = rebuilt["ising_3d_open_2x2"]["factors"]
    mapped_edges = {
        tuple(sorted((C4_TO_SQUARE_SITE_MAP[left], C4_TO_SQUARE_SITE_MAP[right])))
        for left, right in CASES[1]["bonds"]
    }
    square_edges = {tuple(sorted(edge)) for edge in CASES[2]["bonds"]}
    exact_factor_equality = [
        ([int(value) for value in factor.all_coeffs()], exponent)
        for factor, exponent in c4_factors
    ] == [
        ([int(value) for value in factor.all_coeffs()], exponent)
        for factor, exponent in square_factors
    ]
    ok(
        "C4/open-square graph and factor identity",
        mapped_edges == square_edges and exact_factor_equality,
    )
    ok(
        "C4 solvable-degree control",
        max(factor.degree() for factor, _exponent in c4_factors) <= 4,
    )

    stored_certificates = artifact["data"]["galois_certificates"]
    for key, expected_primes in EXPECTED_PRIMES.items():
        factors = rebuilt[key]["factors"]
        sextics = [factor for factor, _exponent in factors if factor.degree() == 6]
        sextics.sort(key=lambda factor: tuple(int(value) for value in factor.all_coeffs()))
        selected = sextics[0]
        scale = independent_root_scale(selected)
        raw_coefficients = [int(value) for value in selected.all_coeffs()]
        normalized_coefficients = [
            value // scale**index for index, value in enumerate(raw_coefficients)
        ]
        normalized = sp.Poly.from_list(normalized_coefficients, gens=X, domain=sp.ZZ)
        discriminant = int(sp.discriminant(normalized, X))
        factorization = {
            str(int(prime)): int(exponent)
            for prime, exponent in sorted(sp.factorint(abs(discriminant)).items())
        }
        selected_stored = stored_cases[key]["selected_sextic"]
        exact_selected = (
            int(selected_stored["rational_root_scale_s"]) == scale
            and [int(value) for value in selected_stored["raw_factor_coefficients_desc"]]
            == raw_coefficients
            and [int(value) for value in selected_stored["normalized_factor_coefficients_desc"]]
            == normalized_coefficients
            and int(selected_stored["discriminant"]) == discriminant
            and selected_stored["discriminant_factorization"] == factorization
            and math.isqrt(abs(discriminant)) ** 2 != abs(discriminant)
            and all(sp.isprime(int(prime)) for prime in factorization)
        )
        ok(
            f"clean-room normalized sextic and discriminant for {key}",
            exact_selected,
            f"scale={scale}, discriminant={discriminant}",
        )

        certificate = stored_certificates[key]
        observed_primes: dict[str, int] = {}
        observed_cycles: dict[str, list[int]] = {}
        modular_ok = True
        for witness_name, stored_witness in certificate["witnesses"].items():
            witness_ok, cycle_type = audit_modular_witness(
                normalized_coefficients, discriminant, stored_witness
            )
            modular_ok = modular_ok and witness_ok
            observed_primes[witness_name] = int(stored_witness["prime"])
            observed_cycles[witness_name] = cycle_type
        s6_logic = (
            observed_primes == expected_primes
            and observed_cycles == EXPECTED_CYCLES
            and certificate["galois_group_over_Q"] == "S6"
            and int(certificate["galois_group_order"]) == 720
            and certificate["solvable_group"] is False
        )
        ok(
            f"independent finite-field factors and S6 logic for {key}",
            modular_ok and s6_logic,
            f"primes={observed_primes}, cycles={observed_cycles}",
        )

    counterexample = (
        stored_certificates["ising_2d_open_width6"]["galois_group_over_Q"] == "S6"
        and stored_certificates["ising_3d_open_2x3"]["galois_group_over_Q"] == "S6"
        and tuple(CASES[3]["bonds"]) == tuple((index, index + 1) for index in range(5))
    )
    ok(
        "open-chain counterexample to Galois-solvability contrast",
        counterexample,
    )

    budget.check("final verifier resource audit")
    ok(
        "verifier resource limits",
        time.process_time() - budget.started <= CPU_BUDGET_SECONDS
        and max_rss_bytes() <= RSS_LIMIT_BYTES,
        f"cpu={time.process_time() - budget.started:.6f}s rss={max_rss_bytes()}",
    )
    if FAILED:
        raise SystemExit(f"FAIL test_galois_spectrum: {FAILED}")
    print(
        "PASS test_galois_spectrum: independent exact factors, discriminants, "
        "modular cycle types, and both S6 conclusions replayed",
        flush=True,
    )


if __name__ == "__main__":
    main()
