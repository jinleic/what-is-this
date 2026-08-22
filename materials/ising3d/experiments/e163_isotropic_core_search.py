"""Exact isotropic core search for the all-size spectral no-go (wave 18).

Every row uses the actual-denominator monic integral model from e136 and
reconstructs the exterior-square characteristic polynomial over two prime
fields.  A modular derivative-gcd degree is an upper bound for the
characteristic-zero gcd degree, hence gives the recorded lower bound on the
number of distinct pair products.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT, ROOT / "src", ROOT / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

import numpy as np  # noqa: E402
import e136_pair_product_2x5 as engine  # noqa: E402

ARTIFACT = ROOT / "results" / "spectral" / "isotropic_core_search.json"
PRIMES = tuple(engine.PRIMES)
T0 = Fraction(1, 3)
CPU_BUDGET_SECONDS = 1_200
RSS_CAP_BYTES = 4_000_000_000

# The labels are mathematical names; the listed edge sets are the raw inputs.
CORE_CANDIDATES: tuple[tuple[str, int, tuple[tuple[int, int], ...], str], ...] = (
    ("claw_K1_3", 4, ((0, 1), (0, 2), (0, 3)), "star K_{1,3}"),
    (
        "claw_attached_P1_T5",
        5,
        ((0, 1), (0, 2), (0, 3), (1, 4)),
        "degree-3 tree with branch lengths 2,1,1",
    ),
    ("paw", 4, ((0, 1), (0, 2), (1, 2), (0, 3)), "triangle with one pendant vertex"),
    (
        "grid_2x2_plus_pendant",
        5,
        ((0, 1), (0, 2), (1, 3), (2, 3), (0, 4)),
        "open 2x2 grid with a pendant vertex at site 0",
    ),
    ("star_K1_4", 5, ((0, 1), (0, 2), (0, 3), (0, 4)), "star K_{1,4} control"),
)


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gaussian_ceiling(n: int) -> int:
    return 3**n - 2**n


def maximum_degree(n: int, bonds: tuple[tuple[int, int], ...]) -> int:
    degrees = [0] * n
    for u, v in bonds:
        degrees[u] += 1
        degrees[v] += 1
    return max(degrees)


def check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def run_core_search() -> dict[str, object]:
    """Run every declared core row; a budget expiry remains explicitly non-decisive."""
    cpu0 = time.process_time()
    rows: list[dict[str, object]] = []
    checks: list[dict[str, object]] = []

    for label, n, bonds_tuple, geometry in CORE_CANDIDATES:
        elapsed = time.process_time() - cpu0
        if elapsed >= CPU_BUDGET_SECONDS:
            rows.append(
                {
                    "label": label,
                    "sites": n,
                    "bonds": [list(edge) for edge in bonds_tuple],
                    "status": "NON-DECISIVE",
                    "blocking_reason": "predeclared total process-time budget expired",
                }
            )
            continue
        record = engine.build_and_run(
            name=label,
            geometry=geometry,
            n_sites=n,
            bonds=list(bonds_tuple),
            t=T0,
        )
        record["label"] = label
        record["sites"] = n
        record["bonds"] = [list(edge) for edge in bonds_tuple]
        record["max_degree"] = maximum_degree(n, bonds_tuple)
        record["gaussian_ceiling"] = gaussian_ceiling(n)
        record["status"] = "decisive" if record.get("blocked_stage") is None else "NON-DECISIVE"
        rows.append(record)

    for row in rows:
        covered = row.get("blocked_stage") is None and row.get("status") == "decisive"
        check(
            checks,
            f"{row['label']}: both modular computations completed",
            covered,
            "blocked=" + str(row.get("blocked_stage") or "none"),
        )
        if not covered:
            continue
        per_prime = row["per_prime"]
        degree_pair = [int(per_prime[str(p)]["gcd_degree"]) for p in PRIMES]
        audits = all(
            bool(per_prime[str(p)]["pair_poly_monic"])
            and bool(per_prime[str(p)]["gcd_monic"])
            and bool(per_prime[str(p)]["gcd_divides_both"])
            and bool(per_prime[str(p)]["horner_charpoly_zero_matrix"])
            and bool(per_prime[str(p)]["spot_trace_checks_all_match"])
            and bool(per_prime[str(p)]["float64_vs_int64_first_product_equal"])
            for p in PRIMES
        )
        check(
            checks,
            f"{row['label']}: monic reconstruction audits",
            audits,
            f"gcd degrees {degree_pair} at primes {list(PRIMES)}",
        )
        expected_lower = int(row["pair_count"]) - min(degree_pair)
        check(
            checks,
            f"{row['label']}: one-sided direction recorded correctly",
            int(row["distinct_pair_products_lower_bound"]) == expected_lower,
            f"r_Q >= {row['pair_count']} - min{degree_pair} = {expected_lower}",
        )

    winners = [
        row["label"]
        for row in rows
        if row.get("blocked_stage") is None
        and int(row["distinct_pair_products_lower_bound"]) > int(row["gaussian_ceiling"])
    ]
    check(
        checks,
        "at least one uniform isotropic core exceeds its Gaussian ceiling",
        bool(winners),
        ", ".join(winners) if winners else "no decisive winner",
    )
    check(
        checks,
        "resource cap",
        max_rss_bytes() < RSS_CAP_BYTES,
        f"peak RSS {max_rss_bytes()} < {RSS_CAP_BYTES}",
    )

    return {
        "parameters": {
            "t": f"{T0.numerator}/{T0.denominator}",
            "q": "5/3",
            "primes": list(PRIMES),
            "total_process_time_budget_seconds": CPU_BUDGET_SECONDS,
            "rss_cap_bytes": RSS_CAP_BYTES,
        },
        "rows": rows,
        "winning_cores": winners,
        "checks": checks,
        "process_time_seconds": round(time.process_time() - cpu0, 3),
        "peak_rss_bytes": max_rss_bytes(),
    }


def make_artifact(data: dict[str, object]) -> dict[str, object]:
    engine_path = ROOT / "experiments" / "e136_pair_product_2x5.py"
    builder_path = ROOT / "experiments" / "e38_gaussianity_certificate.py"
    return {
        "meta": {
            "experiment": "e163_isotropic_core_search",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "command": ".venv/bin/python experiments/e163_isotropic_core_search.py",
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "numpy_version": np.__version__,
            "producer_sha256": sha256_file(Path(__file__)),
            "engine_sha256": sha256_file(engine_path),
            "canonical_builder_sha256": sha256_file(builder_path),
            "arithmetic": "exact integer/F_p; float64 BLAS only under explicit <2^53 integer guards",
            "claim_direction": "deg gcd_Q <= deg gcd_Fp, hence r_Q >= pair_slots - deg gcd_Fp",
        },
        "data": data,
    }


def main() -> int:
    data = run_core_search()
    artifact = make_artifact(data)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    passed = all(bool(item["passed"]) for item in data["checks"])
    print(f"wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if passed:
        print("PASS e163 isotropic core search", flush=True)
        return 0
    print("FAIL e163 isotropic core search", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
