"""Produce the exact characteristic-zero rank/count certificate for Theorem U2.

The output is results/spectral/rank_char0.json.  The characteristic-zero work is
performed by e166/e167.  Previously certified modular lower bounds for T disjoint
union P_2 and P_4 are then squeezed against the now-exact Lemma 5' upper bounds.

Run from the repository root with:
    .venv/bin/python experiments/e168_rank_char0.py
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
from pathlib import Path

import sympy as sp

from e166_rank_char0_spectrum import (
    CPU_BUDGET_SECONDS,
    NAMED_T,
    ROOT,
    CpuBudget,
    exact_spectrum_record,
    fraction_text,
    max_rss_bytes,
    physical_weight,
)
from e167_rank_char0_products import exact_pair_product_counts

OUTPUT = ROOT / "results" / "spectral" / "rank_char0.json"
PRIOR_ARTIFACT = ROOT / "results" / "spectral" / "allsize_gaussian.json"
PROOF_NOTE = ROOT / "proofs" / "rank_char0.md"
SOURCE_PRIME = 1_000_003


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(
    checks: list[dict[str, object]],
    name: str,
    passed: bool,
    detail: str,
) -> None:
    row = {"name": name, "passed": bool(passed), "detail": detail}
    checks.append(row)
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def exact_tensor_upper(m: int, r_core: int, r_all_core: int) -> int:
    return (3**m - 2**m) * r_all_core + 2**m * r_core


def gaussian_ceiling(n: int) -> int:
    return 3**n - 2**n


def main() -> None:
    started_at = utc_now()
    process_started = time.process_time()
    budget = CpuBudget(limit_seconds=CPU_BUDGET_SECONDS)
    checks: list[dict[str, object]] = []

    prior_bytes = PRIOR_ARTIFACT.read_bytes()
    prior = json.loads(prior_bytes)
    prior_sha = hashlib.sha256(prior_bytes).hexdigest()
    check(
        checks,
        "prior certificate checks passed",
        all(bool(row.get("ok")) for row in prior.get("checks", [])),
        f"results/spectral/allsize_gaussian.json sha256={prior_sha}",
    )

    point_rows: list[dict[str, object]] = []
    exact_by_t: dict[str, dict[str, object]] = {}
    expected_pattern = [
        {"degree": 1, "exponent": 2},
        {"degree": 1, "exponent": 2},
        {"degree": 3, "exponent": 1},
        {"degree": 3, "exponent": 1},
        {"degree": 11, "exponent": 1},
        {"degree": 11, "exponent": 1},
    ]

    for t in NAMED_T:
        budget.check(f"start named point {fraction_text(t)}")
        spectrum, factors = exact_spectrum_record(t, budget)
        counts = exact_pair_product_counts(factors, budget)
        row = {**spectrum, "claim_tag": "[COMPUTATION]", "pair_product_certificate": counts}
        point_rows.append(row)
        exact_by_t[fraction_text(t)] = row
        characteristic = spectrum["characteristic_polynomial"]
        check(
            checks,
            f"exact factorization at t={fraction_text(t)}",
            characteristic["factorization_reconstructs"] is True
            and characteristic["factor_degree_exponent_pattern"] == expected_pattern
            and characteristic["distinct_eigenvalue_values"] == 30,
            "charpoly degree 32; factors (1^2,1^2,3,3,11,11) over Q",
        )
        check(
            checks,
            f"two-sided characteristic-zero counts at t={fraction_text(t)}",
            counts["r_two_sided_sandwich"] == {"lower": 417, "upper": 417}
            and counts["r_all_two_sided_sandwich"] == {"lower": 445, "upper": 445},
            "417 <= r(T) <= 417 and 445 <= r_all(T) <= 445",
        )
        check(
            checks,
            f"characteristic-zero derivative gcd degrees at t={fraction_text(t)}",
            counts["pair_derivative_gcd_degree_over_Q"] == 79
            and counts["symmetric_derivative_gcd_degree_over_Q"] == 83,
            "deg gcd_Q(C_2,C_2')=496-417=79; deg gcd_Q(C_all,C_all')=528-445=83",
        )

    check(
        checks,
        "all four declared characteristic-zero rows covered",
        set(exact_by_t) == {fraction_text(value) for value in NAMED_T},
        ", ".join(sorted(exact_by_t)),
    )
    check(
        checks,
        "core excess and square gap",
        445 - 3**5 == 202 and 2**5 - (445 - 417) == 4,
        "E_T=445-243=202 and 32-(445-417)=4",
    )

    prior_chain = {str(block["t"]): block for block in prior["uniform_field_chain_route"]}
    u2_rows: list[dict[str, object]] = []
    for t in NAMED_T:
        t_text = fraction_text(t)
        source = prior_chain[t_text]
        source_m2 = next(row for row in source["rows"] if int(row["m"]) == 2)
        exact_core = exact_by_t[t_text]["pair_product_certificate"]
        upper = exact_tensor_upper(2, int(exact_core["r"]), int(exact_core["r_all"]))
        lower = int(source_m2["certified_r_lower_bound"])
        all_upper = 3**2 * int(exact_core["r_all"])
        all_lower = int(source_m2["certified_r_all_lower_bound"])
        excess = lower - gaussian_ceiling(7)
        formula_excess = 202 * 3**2 + 4 * 2**2
        theorem_row = {
            "claim_tag": "[THEOREM]",
            "m": 2,
            "n": 7,
            "t": t_text,
            "w": fraction_text(physical_weight(t)),
            "modular_lower_bound_direction": f"r >= {lower} over characteristic zero",
            "modular_certificate_prime": SOURCE_PRIME,
            "lemma5prime_upper_bound_direction": f"r <= {upper} using exact core counts",
            "two_sided_sandwich": {"lower": lower, "upper": upper},
            "exact_r": lower if lower == upper else None,
            "r_all_two_sided_sandwich": {"lower": all_lower, "upper": all_upper},
            "exact_r_all": all_lower if all_lower == all_upper else None,
            "injectivity_INJ_m_certified": lower == upper,
            "gaussian_ceiling": gaussian_ceiling(7),
            "exact_excess": excess if lower == upper else None,
            "excess_formula_value": formula_excess,
            "source_artifact": "results/spectral/allsize_gaussian.json",
        }
        u2_rows.append(theorem_row)
        check(
            checks,
            f"U2 squeeze m=2 at t={t_text}",
            lower == upper == 3893
            and all_lower == all_upper == 4005
            and excess == formula_excess == 1834,
            f"3893 <= r <= 3893; excess {excess}; INJ(2) follows from equality",
        )

    source_p4 = prior["uniform_T_plus_P4"]
    exact_core = exact_by_t["1/3"]["pair_product_certificate"]
    upper_p4 = exact_tensor_upper(4, int(exact_core["r"]), int(exact_core["r_all"]))
    lower_p4 = int(source_p4["certified_r_lower_bound"])
    excess_p4 = lower_p4 - gaussian_ceiling(9)
    formula_p4 = 202 * 3**4 + 4 * 2**4
    u2_p4 = {
        "claim_tag": "[THEOREM]",
        "m": 4,
        "n": 9,
        "t": "1/3",
        "w": "5/3",
        "modular_lower_bound_direction": f"r >= {lower_p4} over characteristic zero",
        "modular_certificate_prime": int(source_p4["prime"]),
        "pair_polynomial_sha256_mod_prime": source_p4["pair_poly_sha256"],
        "lemma5prime_upper_bound_direction": f"r <= {upper_p4} using exact core counts",
        "two_sided_sandwich": {"lower": lower_p4, "upper": upper_p4},
        "exact_r": lower_p4 if lower_p4 == upper_p4 else None,
        "injectivity_INJ_m_certified": lower_p4 == upper_p4,
        "gaussian_ceiling": gaussian_ceiling(9),
        "exact_excess": excess_p4 if lower_p4 == upper_p4 else None,
        "excess_formula_value": formula_p4,
        "source_artifact": "results/spectral/allsize_gaussian.json",
    }
    u2_rows.append(u2_p4)
    check(
        checks,
        "U2 squeeze m=4 at t=1/3",
        lower_p4 == upper_p4 == 35597 and excess_p4 == formula_p4 == 16426,
        f"35597 <= r <= 35597; excess {excess_p4}; INJ(4) follows from equality",
    )

    generic_rows: list[dict[str, object]] = []
    for m in (2, 4):
        n = 5 + m
        excess = 202 * 3**m + 4 * 2**m
        generic_rows.append({
            "claim_tag": "[THEOREM]",
            "m": m,
            "n": n,
            "pair_slot_degree_S": comb(1 << n, 2),
            "generic_lower_bound": gaussian_ceiling(n) + excess,
            "gaussian_ceiling": gaussian_ceiling(n),
            "transported_excess_lower_bound": excess,
            "exceptional_hypersurface_degree_bound": "4*(2*n+|E|)*S*(S-1)",
            "direction": (
                "[THEOREM] Outside the nonzero resultant hypersurface, generic r is "
                "at least the exact decoupled count; equality at the generic point "
                "is not claimed."
            ),
            "all_F2_hypotheses_discharged_for_this_m": True,
        })
    check(
        checks,
        "strengthened F2 specializations",
        [row["generic_lower_bound"] for row in generic_rows] == [3893, 35597],
        "m=2 and m=4 have no remaining CZ/INJ hypothesis at the named witness",
    )

    data: dict[str, object] = {
        "definitions": {
            "claim_tag": "[LEMMA]",
            "graph_T": {
                "vertices": [0, 1, 2, 3, 4],
                "edges": [[0, 1], [0, 2], [0, 3], [3, 4]],
                "description": "degree-3 root 0 with branch lengths 1,1,2",
            },
            "operator": (
                "R_T(t,w)=P(t) D_T(w) P(t), P(t)=tensor_v [[1,t],[t,1]], "
                "D_T(w)[sigma,sigma]=w^(number of aligned T edges)."
            ),
            "r": "number of distinct lambda_i*lambda_j over eigenvalue slots i<j",
            "r_all": "number of distinct lambda_i*lambda_j over eigenvalue slots i<=j",
            "slot_multiplicity_rule": (
                "Repeated eigenvalue values retain separate slots; equality is imposed only "
                "after the allowed unordered slot pairs are formed."
            ),
        },
        "exact_characteristic_zero_points": point_rows,
        "core_theorem": {
            "claim_tag": "[THEOREM]",
            "scope": (
                "[THEOREM] At each named physical point t in {1/3,1/4,2/5,1/5}, "
                "w=(1+t^2)/(2t), exactly r(T)=417 and r_all(T)=445."
            ),
            "r": 417,
            "r_all": 445,
            "r_all_minus_r": 28,
            "core_excess_E_T": 202,
            "exact_gcd_degrees": {"pair": 79, "symmetric": 83},
        },
        "restored_U2_named_points": u2_rows,
        "excess_identity": {
            "claim_tag": "[THEOREM]",
            "under_certified_INJ_m": (
                "[THEOREM] r(T disjoint_union P_m)-(3^(5+m)-2^(5+m)) "
                "=202*3^m+4*2^m"
            ),
            "certified_m": [2, 4],
            "values": {"m=2": 1834, "m=4": 16426},
        },
        "F2_status": {
            "claim_tag": "[THEOREM]",
            "CZ": (
                "[THEOREM] DISCHARGED by the exact t=1/3,w=5/3 core witness "
                "(and independently at the other three named physical points)."
            ),
            "INJ": (
                "[THEOREM] DISCHARGED for m=2 and m=4 by the two-sided squeezes; "
                "for the all-n theorem, INJ(n-5) remains the sole hypothesis at "
                "other chain lengths."
            ),
            "unconditional_generic_rows": generic_rows,
            "fully_isotropic_y_equals_w": (
                "[UNRESOLVED] The generic resultant theorem does not decide whether "
                "y=w lies on the exceptional hypersurface."
            ),
        },
        "scope_limits": [
            "[UNRESOLVED] No exact core count is claimed for arbitrary t or w beyond the four named points.",
            "[UNRESOLVED] No INJ(m) claim is made for chain lengths other than m=2 and m=4.",
            "[UNRESOLVED] The generic three-parameter conclusion is a lower bound outside a hypersurface, not an exact generic count.",
            "[UNRESOLVED] The physical fully isotropic specialization y=w is not located relative to that hypersurface.",
        ],
        "checks": checks,
    }

    source_files = [
        Path(__file__).resolve(),
        ROOT / "experiments" / "e166_rank_char0_spectrum.py",
        ROOT / "experiments" / "e167_rank_char0_products.py",
    ]
    meta: dict[str, object] = {
        "producer": "experiments/e168_rank_char0.py",
        "producer_sha256": sha256_file(Path(__file__).resolve()),
        "source_sha256": {str(path.relative_to(ROOT)): sha256_file(path) for path in source_files},
        "prior_certificate": {
            "path": "results/spectral/allsize_gaussian.json",
            "sha256": prior_sha,
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
        "arithmetic": "exact integers/rationals and exact polynomial arithmetic over Q only",
        "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
        "process_cpu_seconds": round(time.process_time() - process_started, 6),
        "peak_rss_bytes": max_rss_bytes(),
        "rss_limit_bytes": 4_000_000_000,
    }

    all_passed = all(bool(row["passed"]) for row in checks)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps({"meta": meta, "data": data}, indent=2, sort_keys=True) + "\n")
    temporary.replace(OUTPUT)
    if not all_passed:
        raise SystemExit("FAIL e168_rank_char0: one or more checks failed")
    print(
        f"PASS e168_rank_char0: r(T)=417, r_all(T)=445; artifact {OUTPUT.relative_to(ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
