"""Integrate the exact isotropic trace-invariant certificate.

This producer writes results/spectral/isotropic_invariant.json.  It uses only
exact integer/rational graph enumerations and symbolic polynomial arithmetic;
there is no modular-rank inference and no finite-size extrapolation.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _path in (ROOT, ROOT / "experiments"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from e200_isotropic_puiseux import run_puiseux  # noqa: E402
from e201_isotropic_invariant import run_invariant  # noqa: E402

ARTIFACT = ROOT / "results" / "spectral" / "isotropic_invariant.json"
CPU_BUDGET_SECONDS = 30.0
RSS_CAP_BYTES = 2_000_000_000


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


class Budget:
    def __init__(self) -> None:
        self.started = time.process_time()

    def used(self) -> float:
        return time.process_time() - self.started

    def check(self, stage: str) -> None:
        used = self.used()
        rss = max_rss_bytes()
        if used >= CPU_BUDGET_SECONDS:
            raise RuntimeError(f"process-time budget exceeded at {stage}: {used:.3f}s")
        if rss >= RSS_CAP_BYTES:
            raise RuntimeError(f"RSS cap exceeded at {stage}: {rss} bytes")


def add_check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}", flush=True)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_all() -> tuple[dict[str, object], list[dict[str, object]], Budget]:
    budget = Budget()
    checks: list[dict[str, object]] = []

    puiseux = run_puiseux()
    budget.check("high-temperature component")
    invariant = run_invariant()
    budget.check("spectral-invariant component")

    for item in puiseux["checks"]:
        add_check(checks, f"e200: {item['name']}", bool(item["passed"]), str(item["detail"]))
    for item in invariant["checks"]:
        add_check(checks, f"e201: {item['name']}", bool(item["passed"]), str(item["detail"]))

    puiseux_rows = puiseux["data"]["graph_rows"]
    invariant_rows = invariant["data"]["physical_graph_reduction"]["graph_rows"]
    puiseux_labels = {str(row["label"]) for row in puiseux_rows}
    invariant_labels = {str(row["label"]) for row in invariant_rows}
    add_check(
        checks,
        "component graph labels agree",
        puiseux_labels == invariant_labels and len(puiseux_labels) > 0,
        f"labels={sorted(puiseux_labels)}",
    )

    bipartite_by_label = {str(row["label"]): bool(row["bipartite"]) for row in puiseux_rows}
    psi_zero_by_label = {str(row["label"]): bool(row["psi_zero"]) for row in invariant_rows}
    add_check(
        checks,
        "trace route vanishes exactly on bipartite audit rows",
        bipartite_by_label.keys() == psi_zero_by_label.keys()
        and all(bipartite_by_label[label] == psi_zero_by_label[label] for label in bipartite_by_label),
        f"bipartite={sum(bipartite_by_label.values())}, nonbipartite={len(bipartite_by_label)-sum(bipartite_by_label.values())}",
    )

    nonbipartite_rows = [row for row in invariant_rows if not bool(row["psi_zero"])]
    add_check(
        checks,
        "every nonbipartite audit row is nonzero on the physical curve",
        len(nonbipartite_rows) > 0
        and all(bool(row["nonzero_modulo_physical_curve"]) for row in nonbipartite_rows)
        and all(int(row["psi_at_y_5_over_3"]["numerator"]) > 0 for row in nonbipartite_rows),
        f"certified audit rows={len(nonbipartite_rows)} at exact y=5/3",
    )

    lollipop_rows = invariant["data"]["leaf_extension"]["finite_lemma_tests"]
    add_check(
        checks,
        "connected branching family recurrence audits",
        len(lollipop_rows) >= 2
        and all(bool(row["checks_passed"]) for row in lollipop_rows)
        and all(int(row["sites"]) == int(row["tail_edges_r"]) + 3 for row in lollipop_rows),
        f"tail lengths={[row['tail_edges_r'] for row in lollipop_rows]}",
    )

    budget.check("final resource audit")
    add_check(
        checks,
        "process-time and RSS budgets",
        budget.used() < CPU_BUDGET_SECONDS and max_rss_bytes() < RSS_CAP_BYTES,
        (
            f"CPU={budget.used():.3f}s<{CPU_BUDGET_SECONDS}; "
            f"RSS={max_rss_bytes()}<{RSS_CAP_BYTES}"
        ),
    )

    data = {
        "high_temperature_leading_terms": puiseux["data"],
        "spectral_invariant": invariant["data"],
        "conclusion": {
            "tag": "[THEOREM]",
            "headline": (
                "Every finite simple non-bipartite layer graph is spectrally non-Gaussian at "
                "every real physical isotropic coupling 0<t<1."
            ),
            "all_size_quantifier": "all finite n; no finite graph census enters the proof",
            "nonzero_at_one_point": (
                "at t=1/3, (x,y)=(2,5/3) lies on F_phys and the odd-Eulerian sum makes "
                "Psi_G(5/3)>0 for every non-bipartite G"
            ),
            "finite_exception_conclusion": (
                "for fixed non-bipartite G, subset-product points on the nonsingular complex "
                "physical curve lie among at most 2m*2^(n-1) roots of a nonzero cleared polynomial"
            ),
            "every_coupling_conclusion": (
                "for real 0<t<1, q(t)>1 and positivity excludes every root, so the exceptional "
                "set in the physical interval is empty"
            ),
            "new_infinite_family": (
                "triangle-with-tail graphs of every tail length r>=1 are connected and branching; "
                "their trace skew is 2(y-1)^3(y+1)^r"
            ),
        },
        "route_boundary": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "Bipartite graphs have Q_G=Q_G^*, so this determinant/trace-reciprocity invariant "
                "vanishes identically. Open rectangular grids and the universal every-branching-graph "
                "isotropic theorem remain unresolved by this route."
            ),
            "not_claimed": [
                "no solution of the three-dimensional Ising model",
                "no conclusion for bipartite branching layers",
                "no promotion of the finite graph rows to an all-size theorem",
            ],
        },
    }
    return data, checks, budget


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]], budget: Budget
) -> dict[str, object]:
    sources = [
        ROOT / "experiments" / "e200_isotropic_puiseux.py",
        ROOT / "experiments" / "e201_isotropic_invariant.py",
        Path(__file__),
    ]
    return {
        "meta": {
            "experiment": "e202_isotropic_certificate",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e202_isotropic_certificate.py",
            "repository_root": str(ROOT),
            "working_directory": str(Path.cwd()),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": (
                "exact integer spin enumeration, exact Eulerian-subgraph enumeration, integer "
                "polynomial arithmetic, and Fraction evaluation; no floating point or modular rank"
            ),
            "proof_scope": (
                "finite computations audit identities only; all-size claims follow from exact "
                "complement pairing and the Eulerian expansion"
            ),
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
                "observed_process_time_seconds": round(budget.used(), 6),
                "observed_peak_rss_bytes": max_rss_bytes(),
            },
            "source_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in sources
            },
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks, budget = run_all()
    artifact = make_artifact(data, checks, budget)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {ARTIFACT.relative_to(ROOT)}", flush=True)
    if all(bool(item["passed"]) for item in checks):
        print("PASS e202 isotropic invariant certificate", flush=True)
        return 0
    print("FAIL e202 isotropic invariant certificate", flush=True)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
