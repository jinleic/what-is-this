"""[THEOREM] Assemble the exact small-layer Galois-spectrum certificate.

The artifact has top-level ``meta``, ``data``, and ``checks`` entries.  Exact
characteristic-zero factorization is delegated to e203 and the explicit
finite-field/Frobenius certificates to e204.

Run from the repository root with:
    .venv/bin/python experiments/e205_galois_certificate.py
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import sympy as sp

from e203_galois_charpoly import (
    CASES,
    CPU_BUDGET_SECONDS,
    Q_SPECIALIZATION,
    ROOT,
    RSS_LIMIT_BYTES,
    TARGET_CASES,
    T_SPECIALIZATION,
    CpuBudget,
    compute_all_cases,
    fraction_text,
    max_rss_bytes,
)
from e204_galois_cycles import certify_s6

OUTPUT = ROOT / "results" / "spectral" / "galois_spectrum.json"
PROOF_NOTE = ROOT / "proofs" / "galois_spectrum.md"
EXPECTED_LAYER_SEXTIC = [
    1,
    -6613,
    6886615,
    -2793624000,
    498389850000,
    -34838100000000,
    437400000000000,
]
EXPECTED_CHAIN_SEXTIC = [
    1,
    -782,
    204493,
    -20321664,
    615904560,
    -7138368000,
    27993600000,
]
C4_TO_SQUARE_SITE_MAP = (0, 1, 3, 2)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_check(
    checks: list[dict[str, object]],
    name: str,
    passed: bool,
    detail: str,
    claim_tag: str = "[COMPUTATION]",
) -> None:
    row = {
        "claim_tag": claim_tag,
        "name": name,
        "passed": bool(passed),
        "detail": detail,
    }
    checks.append(row)
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def _undirected_edges(edges: list[list[int]] | tuple[tuple[int, int], ...]) -> set[tuple[int, int]]:
    return {tuple(sorted((int(left), int(right)))) for left, right in edges}


def _mapped_edges(
    edges: list[list[int]], mapping: tuple[int, ...]
) -> set[tuple[int, int]]:
    return {
        tuple(sorted((mapping[int(left)], mapping[int(right)])))
        for left, right in edges
    }


def _factor_coefficients(case: dict[str, object]) -> list[list[str]]:
    characteristic = case["characteristic_polynomial"]
    return [list(row["coefficients_desc"]) for row in characteristic["factors"]]


def _cycle_prime_map(certificate: dict[str, object]) -> dict[str, int]:
    return {
        name: int(row["prime"])
        for name, row in certificate["witnesses"].items()
    }


def main() -> None:
    started_at = utc_now()
    process_started = time.process_time()
    budget = CpuBudget(limit_seconds=CPU_BUDGET_SECONDS)
    checks: list[dict[str, object]] = []

    cases, _factors = compute_all_cases(budget)
    by_key = {str(row["key"]): row for row in cases}
    certificates = {
        key: certify_s6(by_key[key]["selected_sextic"]) for key in TARGET_CASES
    }
    budget.check("after Galois certificates")

    expected_keys = [spec.key for spec in CASES]
    add_check(
        checks,
        "declared exact case coverage",
        [row["key"] for row in cases] == expected_keys,
        f"computed cases {', '.join(str(row['key']) for row in cases)}",
    )

    exact_factor_audits = all(
        int(row["operator_dimension"])
        == int(row["characteristic_polynomial"]["degree"])
        and row["primitive_matrix_symmetric"] is True
        and row["characteristic_polynomial"]["factorization_reconstructs"] is True
        and row["characteristic_polynomial"]["all_factors_irreducible_over_Q"] is True
        and sum(
            int(factor["degree"]) * int(factor["exponent"])
            for factor in row["characteristic_polynomial"]["factors"]
        )
        == int(row["operator_dimension"])
        for row in cases
    )
    add_check(
        checks,
        "all exact characteristic factorizations reconstruct",
        exact_factor_audits,
        "every monic irreducible factor list reconstructs a degree-2^n charpoly",
    )

    one_dimensional = by_key["ising_1d_single_spin"]
    one_dimensional_factors = _factor_coefficients(one_dimensional)
    one_dimensional_exact = one_dimensional_factors == [["1", "-8"], ["1", "-2"]]
    add_check(
        checks,
        "one-dimensional transfer control",
        one_dimensional_exact,
        f"primitive charpoly factors are {one_dimensional_factors}",
        "[THEOREM]",
    )

    c4 = by_key["ising_2d_periodic_width4"]
    square = by_key["ising_3d_open_2x2"]
    mapped_c4 = _mapped_edges(c4["bonds"], C4_TO_SQUARE_SITE_MAP)
    square_edges = _undirected_edges(square["bonds"])
    graph_isomorphism = mapped_c4 == square_edges
    equal_factorization = (
        c4["characteristic_polynomial"]["factors"]
        == square["characteristic_polynomial"]["factors"]
    )
    add_check(
        checks,
        "open 2x2 square equals the periodic C4 control",
        graph_isomorphism and equal_factorization,
        f"site map {list(C4_TO_SQUARE_SITE_MAP)} maps all bonds and exact factors agree",
        "[THEOREM]",
    )
    c4_degrees = [
        int(row["degree"]) for row in c4["characteristic_polynomial"]["factors"]
    ]
    c4_solvable = max(c4_degrees) <= 4
    add_check(
        checks,
        "C4 full splitting field is solvable",
        c4_solvable,
        f"irreducible factor degrees are {c4_degrees}, all at most four",
        "[THEOREM]",
    )

    layer_selected = by_key["ising_3d_open_2x3"]["selected_sextic"]
    chain_selected = by_key["ising_2d_open_width6"]["selected_sextic"]
    layer_coefficients = [
        int(value) for value in layer_selected["normalized_factor_coefficients_desc"]
    ]
    chain_coefficients = [
        int(value) for value in chain_selected["normalized_factor_coefficients_desc"]
    ]
    add_check(
        checks,
        "explicit 2x3 layer sextic",
        layer_coefficients == EXPECTED_LAYER_SEXTIC
        and int(layer_selected["rational_root_scale_s"]) == 57_600,
        f"g_3D coefficients {layer_coefficients}; rational root scale 57600",
    )
    add_check(
        checks,
        "explicit open-chain control sextic",
        chain_coefficients == EXPECTED_CHAIN_SEXTIC
        and int(chain_selected["rational_root_scale_s"]) == 4_800,
        f"g_open coefficients {chain_coefficients}; rational root scale 4800",
    )

    for key, expected_primes in (
        (
            "ising_3d_open_2x3",
            {"irreducible_sextic": 13, "five_cycle": 53, "transposition": 197},
        ),
        (
            "ising_2d_open_width6",
            {"irreducible_sextic": 7, "five_cycle": 43, "transposition": 587},
        ),
    ):
        certificate = certificates[key]
        prime_map = _cycle_prime_map(certificate)
        witness_checks = all(
            bool(row["prime_verified"])
            and bool(row["good_prime"])
            and bool(row["factor_product_reconstructs"])
            and bool(row["squarefree_mod_prime"])
            and bool(row["all_displayed_factors_irreducible_over_Fp"])
            for row in certificate["witnesses"].values()
        )
        expected_cycles = {
            "irreducible_sextic": [6],
            "five_cycle": [1, 5],
            "transposition": [1, 1, 1, 1, 2],
        }
        observed_cycles = {
            name: list(row["frobenius_cycle_type"])
            for name, row in certificate["witnesses"].items()
        }
        add_check(
            checks,
            f"explicit good-prime Frobenius witnesses for {key}",
            witness_checks
            and prime_map == expected_primes
            and observed_cycles == expected_cycles,
            f"primes {prime_map}; cycles {observed_cycles}",
        )
        add_check(
            checks,
            f"Galois group S6 for {key}",
            certificate["galois_group_over_Q"] == "S6"
            and int(certificate["galois_group_order"]) == math.factorial(6)
            and certificate["solvable_group"] is False
            and all(bool(value) for value in certificate["group_theorem_conditions"].values()),
            "irreducible + 5-cycle + transposition at good primes gives S6",
            "[THEOREM]",
        )

    discriminants_exact = all(
        selected["discriminant_factorization_reconstructs"] is True
        and selected["discriminant_factors_all_prime"] is True
        and selected["discriminant_is_square"] is False
        and bool(selected["odd_exponent_discriminant_primes"])
        for selected in (layer_selected, chain_selected)
    )
    add_check(
        checks,
        "exact nonsquare discriminants",
        discriminants_exact,
        (
            f"2x3 discriminant {layer_selected['discriminant']}; "
            f"open-chain discriminant {chain_selected['discriminant']}"
        ),
    )

    contrast_counterexample = (
        certificates["ising_3d_open_2x3"]["galois_group_over_Q"] == "S6"
        and certificates["ising_2d_open_width6"]["galois_group_over_Q"] == "S6"
        and by_key["ising_2d_open_width6"]["bonds"]
        == [[index, index + 1] for index in range(5)]
    )
    add_check(
        checks,
        "non-solvability contrast fails on the open 2D control",
        contrast_counterexample,
        "both the open P6 control and open 2x3 layer have an exact S6 sextic",
        "[THEOREM][EXTERNAL]",
    )

    process_cpu_seconds = time.process_time() - process_started
    peak_rss = max_rss_bytes()
    resource_ok = process_cpu_seconds <= CPU_BUDGET_SECONDS and peak_rss <= RSS_LIMIT_BYTES
    add_check(
        checks,
        "declared resource limits",
        resource_ok,
        f"process CPU {process_cpu_seconds:.6f}s; peak RSS {peak_rss} bytes",
    )

    data: dict[str, object] = {
        "definitions": {
            "claim_tag": "[LEMMA]",
            "specialization": (
                f"[COMPUTATION] t={fraction_text(T_SPECIALIZATION)} and "
                f"q=(1+t^2)/(2t)={fraction_text(Q_SPECIALIZATION)} exactly"
            ),
            "operator": (
                "[LEMMA] R=P_t diag(q^((b_sigma-epsilon)/2)) P_t is the "
                "canonical rational representative from e38; each primitive B and "
                "each displayed root normalization differ from R by a nonzero "
                "rational scalar, so their splitting fields over Q are identical."
            ),
            "factor_convention": (
                "[COMPUTATION] coefficients_desc lists monic coefficients from the "
                "highest power to the constant term; exponents retain charpoly "
                "multiplicity."
            ),
        },
        "exact_cases": cases,
        "galois_certificates": certificates,
        "theorems": {
            "open_2x3_layer": {
                "claim_tag": "[THEOREM]",
                "statement": (
                    "[THEOREM] At the exact rational specialization t=1/3, the "
                    "normalized open 2x3 layer characteristic polynomial has the "
                    "displayed irreducible sextic factor with Galois group S6 over Q; "
                    "hence that factor is not solvable by radicals over Q."
                ),
                "normalized_factor_coefficients_desc": [
                    str(value) for value in layer_coefficients
                ],
                "good_primes": _cycle_prime_map(
                    certificates["ising_3d_open_2x3"]
                ),
            },
            "open_chain_counterexample": {
                "claim_tag": "[THEOREM][EXTERNAL]",
                "statement": (
                    "[THEOREM] The exact open P6 layer control at the same t has its "
                    "own irreducible S6 sextic. [EXTERNAL] This is the standard "
                    "finite-width 2D Ising/free-fermion control. Therefore an S6 "
                    "factor is not a discriminator between this control and the 3D "
                    "layer and, by itself, proves neither nonintegrability nor "
                    "non-Gaussianity."
                ),
                "normalized_factor_coefficients_desc": [
                    str(value) for value in chain_coefficients
                ],
                "good_primes": _cycle_prime_map(
                    certificates["ising_2d_open_width6"]
                ),
            },
            "square_cycle_identity": {
                "claim_tag": "[THEOREM]",
                "statement": (
                    "[THEOREM] The open 2x2 square layer graph is isomorphic to C4 "
                    "under 0,1,2,3 -> 0,1,3,2; the two normalized transfer matrices "
                    "are permutation-similar and have identical exact factors."
                ),
                "site_map_C4_to_square": list(C4_TO_SQUARE_SITE_MAP),
                "all_factor_degrees_at_most_four": c4_solvable,
                "full_splitting_group_solvable": c4_solvable,
            },
        },
        "specialization_direction": {
            "claim_tag": "[LEMMA][UNRESOLVED]",
            "finite_field_to_characteristic_zero": (
                "[LEMMA] For each fixed specialized sextic g in Z[x] and each "
                "displayed p with p not dividing disc(g), modular factor degrees "
                "exhibit a Frobenius element of Gal(g/Q) with that cycle type."
            ),
            "rational_parameter_specialization": (
                "[LEMMA] For a separately constructed separable family F(t,x), a "
                "good specialized Galois group embeds into the generic group. Thus "
                "an S6 specialization would force a degree-six generic group to be "
                "S6 only after identifying this factor in Q(t)[x] and proving the "
                "specialization is good."
            ),
            "not_claimed": (
                "[UNRESOLVED] No degree-six factor family over Q(t) is constructed "
                "here, so the exact t=1/3 theorem is not promoted to generic t, all "
                "physical t, critical coupling, or any all-size statement."
            ),
        },
        "scope_limits": [
            "[UNRESOLVED] The computations concern only the five finite layer graphs explicitly stored in exact_cases.",
            "[UNRESOLVED] K_c=0.221654626 was not used to select, fit, or validate t=1/3 or any prime.",
            "[UNRESOLVED] Non-solvability by radicals of a characteristic factor is not identified with nonintegrability or non-Gaussianity.",
            "[UNRESOLVED] The theorem is for the rationally normalized representative R; no claim about an unrelated normalization field is needed.",
        ],
    }

    source_paths = [
        ROOT / "experiments" / "e203_galois_charpoly.py",
        ROOT / "experiments" / "e204_galois_cycles.py",
        Path(__file__).resolve(),
    ]
    meta: dict[str, object] = {
        "claim_tag": "[COMPUTATION]",
        "producer": "experiments/e205_galois_certificate.py",
        "producer_sha256": sha256_file(Path(__file__).resolve()),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256_file(path) for path in source_paths
        },
        "proof_note": str(PROOF_NOTE.relative_to(ROOT)),
        "started_at_utc": started_at,
        "generated_at_utc": utc_now(),
        "cwd": str(Path.cwd()),
        "argv": sys.argv,
        "environment": {
            "python_executable": sys.executable,
            "python_version": sys.version,
            "sympy_version": sp.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        "arithmetic": (
            "[COMPUTATION] exact integers/rationals, exact factorization over Q, "
            "and exact polynomial arithmetic over finite prime fields; no floating "
            "point enters a mathematical decision"
        ),
        "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
        "process_cpu_seconds": round(process_cpu_seconds, 6),
        "peak_rss_bytes": peak_rss,
        "rss_limit_bytes": RSS_LIMIT_BYTES,
    }

    artifact = {"meta": meta, "data": data, "checks": checks}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUTPUT)

    if not all(bool(row["passed"]) for row in checks):
        raise SystemExit("FAIL e205_galois_certificate: one or more checks failed")
    print(
        "PASS e205_galois_certificate: exact S6 factors in the open 2x3 layer "
        "and the open P6 integrable control; contrast refuted",
        flush=True,
    )


if __name__ == "__main__":
    main()
