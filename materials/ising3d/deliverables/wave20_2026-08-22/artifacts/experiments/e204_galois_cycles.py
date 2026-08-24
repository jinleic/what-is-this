"""[COMPUTATION] Modular Frobenius cycle certificates for selected sextics.

For a monic separable polynomial over ``Z``, factor degrees modulo a prime not
dividing the discriminant give the cycle type of a Frobenius element in the
characteristic-zero Galois group.  This module finds the least good prime below
1000 for each of the requested types ``(6)``, ``(1,5)``, and
``(1,1,1,1,2)`` and stores the explicit modular factors.

Run from the repository root with:
    .venv/bin/python experiments/e204_galois_cycles.py
"""

from __future__ import annotations

import math
import warnings
from typing import Iterable

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

X = sp.Symbol("x")
SEARCH_BOUND_EXCLUSIVE = 1_000
TARGET_CYCLE_TYPES: tuple[tuple[str, tuple[int, ...]], ...] = (
    ("irreducible_sextic", (6,)),
    ("five_cycle", (1, 5)),
    ("transposition", (1, 1, 1, 1, 2)),
)


def _poly_from_coefficients(coefficients: Iterable[int | str]) -> sp.Poly:
    return sp.Poly.from_list([int(value) for value in coefficients], gens=X, domain=sp.ZZ)


def _canonical_modular_coefficients(factor: sp.Poly, prime: int) -> list[int]:
    monic = factor.monic()
    return [int(coefficient) % prime for coefficient in monic.all_coeffs()]


def modular_factorization_record(
    polynomial: sp.Poly,
    discriminant: int,
    prime: int,
) -> dict[str, object]:
    """[COMPUTATION] Factor exactly in F_p[x] and audit every stored factor."""

    if not bool(sp.isprime(prime)):
        raise ValueError(f"nonprime modulus {prime}")
    exact = sp.Poly(polynomial, X, domain=sp.ZZ)
    good_prime = discriminant % prime != 0
    reduced = sp.Poly(exact, X, modulus=prime).monic()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SymPyDeprecationWarning)
        content, raw_factors = sp.factor_list(reduced, modulus=prime)
    if int(content) % prime != 1:
        raise ArithmeticError("unexpected nonunit modular factorization content")

    factors: list[tuple[sp.Poly, int, list[int]]] = []
    for raw_factor, exponent in raw_factors:
        factor = sp.Poly(raw_factor, X, modulus=prime).monic()
        coefficients = _canonical_modular_coefficients(factor, prime)
        factors.append((factor, int(exponent), coefficients))
    factors.sort(key=lambda item: (item[0].degree(), tuple(item[2]), item[1]))

    reconstruction = sp.Poly(1, X, modulus=prime)
    factor_rows: list[dict[str, object]] = []
    cycle_type: list[int] = []
    for factor, exponent, coefficients in factors:
        reconstruction *= factor**exponent
        irreducible = bool(factor.is_irreducible)
        factor_rows.append(
            {
                "claim_tag": "[COMPUTATION]",
                "degree": factor.degree(),
                "exponent": exponent,
                "coefficients_desc_mod_p": coefficients,
                "irreducible_over_Fp": irreducible,
            }
        )
        cycle_type.extend([factor.degree()] * exponent)
    cycle_type.sort()
    reconstruction = sp.Poly(reconstruction, X, modulus=prime).monic()
    reconstruction_ok = reconstruction == reduced
    squarefree = sp.gcd(reduced, reduced.diff()).degree() == 0
    all_irreducible = all(bool(row["irreducible_over_Fp"]) for row in factor_rows)
    if good_prime and not squarefree:
        raise ArithmeticError("a good discriminant prime produced a repeated factor")
    if not reconstruction_ok or not all_irreducible:
        raise ArithmeticError("modular factor audit failed")

    return {
        "claim_tag": "[COMPUTATION]",
        "prime": prime,
        "prime_verified": bool(sp.isprime(prime)),
        "discriminant_mod_prime": discriminant % prime,
        "good_prime": good_prime,
        "polynomial_coefficients_desc_mod_p": [
            int(coefficient) % prime for coefficient in exact.all_coeffs()
        ],
        "factors": factor_rows,
        "factor_product_reconstructs": reconstruction_ok,
        "squarefree_mod_prime": squarefree,
        "all_displayed_factors_irreducible_over_Fp": all_irreducible,
        "frobenius_cycle_type": cycle_type,
    }


def find_cycle_witnesses(
    normalized_coefficients: Iterable[int | str],
    discriminant: int,
    search_bound_exclusive: int = SEARCH_BOUND_EXCLUSIVE,
) -> dict[str, dict[str, object]]:
    """[COMPUTATION] Find the least good prime for each requested cycle type."""

    polynomial = _poly_from_coefficients(normalized_coefficients)
    if polynomial.degree() != 6 or polynomial.LC() != 1:
        raise ValueError("the S6 certificate requires a monic sextic")
    wanted = {cycle_type: name for name, cycle_type in TARGET_CYCLE_TYPES}
    found: dict[str, dict[str, object]] = {}
    for prime_value in sp.primerange(2, search_bound_exclusive):
        prime = int(prime_value)
        if discriminant % prime == 0:
            continue
        row = modular_factorization_record(polynomial, discriminant, prime)
        cycle_type = tuple(int(value) for value in row["frobenius_cycle_type"])
        name = wanted.get(cycle_type)
        if name is not None and name not in found:
            row["selection"] = (
                f"[COMPUTATION] least good prime below {search_bound_exclusive} "
                f"with cycle type {cycle_type}"
            )
            found[name] = row
        if len(found) == len(TARGET_CYCLE_TYPES):
            break
    missing = [name for name, _cycle_type in TARGET_CYCLE_TYPES if name not in found]
    if missing:
        raise ArithmeticError(
            f"no good-prime witness below {search_bound_exclusive} for {missing}"
        )
    return {name: found[name] for name, _cycle_type in TARGET_CYCLE_TYPES}


def certify_s6(selected_factor: dict[str, object]) -> dict[str, object]:
    """[THEOREM] Combine exact and Frobenius data into an S6 certificate."""

    coefficients = selected_factor["normalized_factor_coefficients_desc"]
    polynomial = _poly_from_coefficients(coefficients)
    discriminant = int(selected_factor["discriminant"])
    witnesses = find_cycle_witnesses(coefficients, discriminant)

    expected = {name: list(cycle_type) for name, cycle_type in TARGET_CYCLE_TYPES}
    witnessed = {
        name: list(row["frobenius_cycle_type"]) for name, row in witnesses.items()
    }
    all_good = all(
        bool(row["prime_verified"])
        and bool(row["good_prime"])
        and bool(row["factor_product_reconstructs"])
        and bool(row["squarefree_mod_prime"])
        and bool(row["all_displayed_factors_irreducible_over_Fp"])
        for row in witnesses.values()
    )
    irreducible_reduction = witnessed["irreducible_sextic"] == [6]
    has_five_cycle = witnessed["five_cycle"] == [1, 5]
    has_transposition = witnessed["transposition"] == [1, 1, 1, 1, 2]
    conditions = {
        "degree_six": polynomial.degree() == 6,
        "irreducible_mod_prime_hence_irreducible_over_Q": irreducible_reduction,
        "transitive_over_Q": irreducible_reduction,
        "five_cycle_witnessed": has_five_cycle,
        "primitivity_from_transitivity_and_five_cycle": (
            irreducible_reduction and has_five_cycle
        ),
        "transposition_witnessed": has_transposition,
        "all_witness_primes_good": all_good,
        "discriminant_nonzero": discriminant != 0,
        "discriminant_nonsquare": not bool(selected_factor["discriminant_is_square"]),
    }
    is_s6 = all(
        conditions[key]
        for key in (
            "degree_six",
            "irreducible_mod_prime_hence_irreducible_over_Q",
            "five_cycle_witnessed",
            "primitivity_from_transitivity_and_five_cycle",
            "transposition_witnessed",
            "all_witness_primes_good",
            "discriminant_nonzero",
        )
    )
    return {
        "claim_tag": "[THEOREM]",
        "normalized_factor_coefficients_desc": [str(value) for value in coefficients],
        "discriminant": str(discriminant),
        "witness_search_bound_exclusive": SEARCH_BOUND_EXCLUSIVE,
        "requested_cycle_types": expected,
        "witnesses": witnesses,
        "group_theorem_conditions": conditions,
        "galois_group_over_Q": "S6" if is_s6 else "[UNRESOLVED]",
        "galois_group_order": math.factorial(6) if is_s6 else None,
        "solvable_group": False if is_s6 else None,
        "proof": (
            "[THEOREM] The irreducible reduction makes G transitive. A transitive "
            "degree-six group containing a 5-cycle is primitive: a nontrivial block "
            "would have size 2 or 3, while the 5-cycle and its fixed point force a "
            "block of size 1 or 6. Conjugates of a transposition are edges of a "
            "G-invariant graph; primitivity makes that graph connected, and edge "
            "transpositions of a connected graph generate S6. Hence G=S6."
        ),
        "frobenius_direction": (
            "[LEMMA] At p not dividing the discriminant, the displayed factor "
            "degrees are the cycle lengths of an element of Gal(g/Q); modular "
            "factorization supplies elements of the characteristic-zero group, "
            "not an equality of groups by itself."
        ),
    }


def main() -> None:
    from e203_galois_charpoly import CpuBudget, TARGET_CASES, compute_all_cases

    rows, _factors = compute_all_cases(CpuBudget())
    by_key = {str(row["key"]): row for row in rows}
    certificates = {
        key: certify_s6(by_key[key]["selected_sextic"]) for key in TARGET_CASES
    }
    if not all(
        certificate["galois_group_over_Q"] == "S6"
        and certificate["galois_group_order"] == 720
        and certificate["solvable_group"] is False
        for certificate in certificates.values()
    ):
        raise SystemExit("FAIL e204_galois_cycles")
    primes = {
        key: {
            name: row["prime"]
            for name, row in certificate["witnesses"].items()
        }
        for key, certificate in certificates.items()
    }
    print(f"PASS e204_galois_cycles: exact S6 witnesses {primes}", flush=True)


if __name__ == "__main__":
    main()
