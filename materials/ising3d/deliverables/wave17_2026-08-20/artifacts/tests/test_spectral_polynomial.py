"""Standalone exact regression for e69's falsified generalization."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


path = ROOT / "results" / "spectral" / "polynomial_identity.json"
with path.open() as handle:
    stored = json.load(handle)

check("artifact envelope", set(stored) == {"provenance", "data", "checks"})
check("all embedded checks pass", all(item["passed"] for item in stored["checks"]))
check("status is honestly unresolved", stored["data"]["status"] == "UNRESOLVED")
check("invalid theorem marked falsified", stored["data"]["claim_tag"] == "[FALSIFIED]")
check("invalidated flag present", stored["data"]["invalidated"] is True)
check("headline explicitly not achieved", stored["data"]["requested_headline_achieved"] is False)
check("critical benchmark absent", stored["provenance"]["critical_benchmark_used"] is False)

Q2, Q3, Q4, Q5 = sp.symbols("Q2 Q3 Q4 Q5")
candidate = sp.sympify(
    stored["data"]["invalidated_candidate_relation"],
    locals={"Q2": Q2, "Q3": Q3, "Q4": Q4, "Q5": Q5},
)

# Independently derive a genuine one-mode identity using correct squared
# centered traces.  With z=r+r^-1+2, Q1=z and Q2=(z-2)^2, hence
# Q2-(Q1-2)^2=0.  This is deliberately smaller than the failed target.
z, q1, q2 = sp.symbols("z q1 q2")
one_mode_relation = sp.resultant(q1 - z, q2 - (z - 2) ** 2, z)
for mode in (4, 9, 25):
    root = sp.sqrt(mode)
    s1 = root + 1 / root
    s2 = root**2 + root ** (-2)
    residual = sp.factor(one_mode_relation.subs({q1: s1**2, q2: s2**2}))
    check(f"independent correct one-mode identity at r={mode}", residual == 0)

# Independently construct the exact positive six-mode Gaussian control.  Each
# r_i=1, so every centered trace S_k is 2^6=64.
six_modes = (1, 1, 1, 1, 1, 1)
roots = [sp.sqrt(mode) for mode in six_modes]
centered = [
    sp.prod(root**order + root ** (-order) for root in roots)
    for order in range(1, 6)
]
check("positive six-mode control centered traces", centered == [64] * 5)

ratio_residual = sp.factor(candidate.subs({Q2: 64, Q3: 1, Q4: 64, Q5: 1}))
ratio_expected = 64 * (64 - 8) * (64 - 1) * (64 + 1)
check("invalid ratio normalization caught", ratio_residual == ratio_expected != 0)

norm_residual = sp.factor(candidate.subs({Q2: 4096, Q3: 4096, Q4: 4096, Q5: 4096}))
check("six-mode generalization independently falsified", norm_residual != 0)
counterexample = stored["data"]["positive_six_mode_counterexample"]
check("stored ratio residual reproduced", str(ratio_residual) == counterexample["invalid_ratio_residual"])
check("stored norm residual reproduced", str(norm_residual) == counterexample["norm_residual"])
check("control explicitly positive Gaussian", counterexample["is_positive_gaussian_subset_product_spectrum"])

# Independently catch the single-mode normalization error at k=2.
correct_g2 = (z - 2) ** 2
proposed_g2 = z - 2
discrepancy = sp.expand(correct_g2 - proposed_g2)
check("wrong Chebyshev factor caught", discrepancy != 0)
check("stored normalization discrepancy reproduced", str(discrepancy) == stored["data"]["normalization_error"]["difference"])
check(
    "no obstruction or Sturm theorem retained",
    "obstruction" not in stored["data"] and "sturm_certificate" not in stored["data"],
)

if FAILURES:
    print(f"FAIL: {FAILURES}")
    raise SystemExit(1)
print("PASS")
