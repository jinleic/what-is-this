"""Integrate the isotropic-core success and the exact all-size decoupling obstruction.

This producer deliberately does not promote the surviving finite core to an
all-graph theorem.  The final artifact records the strongest valid conclusion:
three named isotropic core theorems (generic in the common coupling) and an
all-size counterfamily showing why Lemma-5 localization cannot close the
requested every-graph quantifier.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT, ROOT / "src", ROOT / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np  # noqa: E402
from e163_isotropic_core_search import PRIMES, run_core_search  # noqa: E402
from e164_isotropic_decoupling import run_decoupling_obstruction  # noqa: E402

ARTIFACT = ROOT / "results" / "spectral" / "isotropic_allsize.json"
CPU_BUDGET_SECONDS = 1_300
RSS_CAP_BYTES = 4_000_000_000
EXPECTED_GCD_DEGREES = {
    "claw_K1_3": 55,
    "claw_attached_P1_T5": 79,
    "paw": 1,
    "grid_2x2_plus_pendant": 21,
    "star_K1_4": 349,
}
WINNERS = ("claw_attached_P1_T5", "paw", "grid_2x2_plus_pendant")


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def isotropic_exception_degree_bound(n: int, edges: int) -> int:
    slots = comb(1 << n, 2)
    return 4 * (2 * n + 2 * edges) * slots * (slots - 1)


def run_all() -> dict[str, object]:
    cpu0 = time.process_time()
    core = run_core_search()
    decoupling = run_decoupling_obstruction()
    checks: list[dict[str, object]] = []

    rows_by_label = {str(row["label"]): row for row in core["rows"]}
    complete_labels = set(rows_by_label) == set(EXPECTED_GCD_DEGREES)
    check(
        checks,
        "declared core search is completely covered",
        complete_labels,
        f"labels={sorted(rows_by_label)}",
    )

    for label, expected_degree in EXPECTED_GCD_DEGREES.items():
        row = rows_by_label.get(label)
        covered = row is not None and row.get("blocked_stage") is None
        check(checks, f"{label}: decisive row present", covered, "no resource wall" if covered else "missing/blocked")
        if not covered:
            continue
        degrees = [int(row["per_prime"][str(p)]["gcd_degree"]) for p in PRIMES]
        check(
            checks,
            f"{label}: exact modular gcd degrees",
            degrees == [expected_degree, expected_degree],
            f"{degrees} at {list(PRIMES)}",
        )
        expected_r = int(row["pair_count"]) - expected_degree
        check(
            checks,
            f"{label}: characteristic-zero lower bound direction",
            int(row["distinct_pair_products_lower_bound"]) == expected_r,
            f"r_Q >= {row['pair_count']}-{expected_degree}={expected_r}",
        )

    named_core_theorems: list[dict[str, object]] = []
    for label in WINNERS:
        row = rows_by_label[label]
        n = int(row["sites"])
        edges = int(row["bond_count"])
        r0 = int(row["distinct_pair_products_lower_bound"])
        ceiling = int(row["gaussian_ceiling"])
        degree_bound = isotropic_exception_degree_bound(n, edges)
        theorem = {
            "tag": "[THEOREM]",
            "graph": label,
            "sites": n,
            "edges": edges,
            "bonds": row["bonds"],
            "named_point_t": "1/3",
            "named_point_q": "5/3",
            "modular_gcd_degrees": [
                int(row["per_prime"][str(p)]["gcd_degree"]) for p in PRIMES
            ],
            "pair_slots": int(row["pair_count"]),
            "characteristic_zero_r_lower_bound": r0,
            "gaussian_ceiling": ceiling,
            "certified_excess": r0 - ceiling,
            "exception_polynomial_degree_bound": degree_bound,
            "statement": (
                f"Along the isotropic curve, outside the roots of a nonzero polynomial "
                f"rho_{label}(t) of degree at most {degree_bound}, the graph has at least "
                f"{r0}>{ceiling} distinct unordered-slot pair products and is not a full "
                f"{n}-mode Gaussian subset-product multiset."
            ),
            "scope": "this named finite graph only; no supergraph or all-size monotonicity is asserted",
        }
        named_core_theorems.append(theorem)
        check(
            checks,
            f"{label}: named core is decisively above ceiling",
            r0 > ceiling,
            f"{r0}>{ceiling}; exceptional degree <= {degree_bound}",
        )
        check(
            checks,
            f"{label}: exceptional degree formula",
            degree_bound == 4 * (2 * n + 2 * edges) * comb(1 << n, 2) * (comb(1 << n, 2) - 1),
            str(degree_bound),
        )

    failed_piece_checks = all(bool(item["passed"]) for item in decoupling["checks"])
    check(
        checks,
        "decoupling obstruction certificate is complete",
        failed_piece_checks,
        f"{len(decoupling['checks'])} exact checks",
    )
    obstruction = decoupling["equal_mode_counterfamily"]["universal_induction_certificate"]
    check(
        checks,
        "all-size counterfamily begins at n=7",
        int(obstruction["base_m"]) == 3
        and int(obstruction["base_upper"]) == 952
        and int(obstruction["base_ceiling"]) == 2059,
        "G_m has n=m+4; 952<2059 at m=3 and the exact induction propagates",
    )

    process_used = time.process_time() - cpu0
    check(
        checks,
        "process-time and RSS budgets",
        process_used < CPU_BUDGET_SECONDS and max_rss_bytes() < RSS_CAP_BYTES,
        f"CPU {process_used:.3f}s < {CPU_BUDGET_SECONDS}; RSS {max_rss_bytes()} < {RSS_CAP_BYTES}",
    )

    conclusion = {
        "tag": "[UNRESOLVED]",
        "requested_isotropic_all_size_theorem_proved": False,
        "n0": None,
        "surviving_piece": (
            "The isotropic core search succeeds: paw, T5, and 2x2-grid-plus-pendant "
            "have certified pair-product excess at t=1/3 and therefore all but finitely "
            "many uniform couplings for each named graph."
        ),
        "failed_piece": (
            "The required every-graph isotropic decoupling identity does not exist in "
            "the Lemma-5 edge-cut form. A nonempty uniform cut has the nonzero minor "
            "q(t)^(2c)-1; the valid family K1,3 disjoint-union m K1 has only identical "
            "free modes and pair count at most 136(2m+1), below the Gaussian ceiling "
            "for every m>=3."
        ),
        "precise_scope": (
            "This is an obstruction to the proposed pair-product/tensor-localization proof, "
            "not a proof that any G_m spectrum is Gaussian and not a disproof of a future "
            "isotropic all-size spectral no-go by another invariant."
        ),
        "what_remains": (
            "A non-localizing isotropic argument, or one certified named-point obstruction "
            "for every graph/size, is still required."
        ),
    }

    return {
        "core_search": core,
        "named_core_theorems": named_core_theorems,
        "decoupling_obstruction": decoupling,
        "conclusion": conclusion,
        "checks": checks,
        "process_time_seconds": round(process_used, 3),
        "peak_rss_bytes": max_rss_bytes(),
    }


def make_artifact(data: dict[str, object]) -> dict[str, object]:
    sources = [
        Path(__file__),
        ROOT / "experiments" / "e163_isotropic_core_search.py",
        ROOT / "experiments" / "e164_isotropic_decoupling.py",
        ROOT / "experiments" / "e136_pair_product_2x5.py",
        ROOT / "experiments" / "e38_gaussianity_certificate.py",
    ]
    return {
        "meta": {
            "experiment": "e165_isotropic_allsize",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "command": ".venv/bin/python experiments/e165_isotropic_allsize.py",
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "numpy_version": np.__version__,
            "primes": list(PRIMES),
            "source_sha256": {str(path.relative_to(ROOT)): sha256_file(path) for path in sources},
            "arithmetic": (
                "raw rational Ising operator -> actual-LCM monic integer matrix -> exact F_p "
                "power traces/Newton reconstruction/Euclidean gcd; symbolic obstruction uses "
                "exact integers and fractions"
            ),
            "modular_direction": (
                "deg gcd_Q <= deg gcd_Fp; every modular row is used only as a lower bound "
                "r_Q >= pair_slots-deg gcd_Fp"
            ),
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
            },
        },
        "data": data,
    }


def main() -> int:
    data = run_all()
    artifact = make_artifact(data)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    all_checks = (
        list(data["checks"])
        + list(data["core_search"]["checks"])
        + list(data["decoupling_obstruction"]["checks"])
    )
    passed = all(bool(item["passed"]) for item in all_checks)
    print(f"wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if passed:
        print("PASS e165 isotropic all-size obstruction", flush=True)
        return 0
    print("FAIL e165 isotropic all-size obstruction", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
