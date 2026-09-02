#!/usr/bin/env python3
"""Additive post-run audit: compare every completed parent prefix cell strictly.

The frozen C26 control omitted ``complete`` together with the two intended
provenance fields. This audit removes only ``header`` and ``t_utc``, as the
pre-statement required, and covers all 14 completed parent cells rather than
only beta 8. It does not alter or retroactively relabel the frozen control.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "2026-08-31T08-44Z_DA168C79"
EXPECTED_PARENT_LEDGER_SHA256 = (
    "6cf79b7adc2c1513d134b7e5f6b5d12080c026af6a92e8a9a15cd4242dbbb6a6"
)
EXPECTED_SUCCESSOR_LEDGER_SHA256 = (
    "0d3abbcb7d8795c6034cd20f9fd8125fe338a6a784305e2e3af69b52c3993937"
)
EXPECTED_BETAS = list(range(8, 22))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_ledger(root: Path, ledger: Path, expected_hash: str) -> int:
    if sha256(ledger) != expected_hash:
        raise RuntimeError(f"ledger hash mismatch: {ledger}")
    count = 0
    for line in ledger.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        expected, relative = line.split("  ", 1)
        artifact = root / relative
        if sha256(artifact) != expected:
            raise RuntimeError(f"artifact hash mismatch: {artifact}")
        count += 1
    return count


def strict_semantic_payload(record: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in record.items() if key not in {"header", "t_utc"}}


def main() -> None:
    parent_rows = verify_ledger(
        PARENT, PARENT / "checksums_partial.sha256", EXPECTED_PARENT_LEDGER_SHA256
    )
    successor_rows = verify_ledger(
        HERE, HERE / "state" / "checksums.sha256", EXPECTED_SUCCESSOR_LEDGER_SHA256
    )

    comparisons = []
    for beta in EXPECTED_BETAS:
        parent_path = PARENT / "state" / f"cell_{beta}.json"
        successor_path = HERE / "state" / f"cell_{beta}.json"
        parent = json.loads(parent_path.read_text(encoding="utf-8"))
        successor = json.loads(successor_path.read_text(encoding="utf-8"))
        equal = strict_semantic_payload(parent) == strict_semantic_payload(successor)
        comparisons.append(
            {
                "beta": beta,
                "strict_semantic_equal": equal,
                "parent_complete": parent.get("complete"),
                "successor_complete": successor.get("complete"),
                "parent_sha256": sha256(parent_path),
                "successor_sha256": sha256(successor_path),
            }
        )

    all_equal = all(row["strict_semantic_equal"] for row in comparisons)
    complete_equal_true = all(
        row["parent_complete"] is True and row["successor_complete"] is True
        for row in comparisons
    )
    result = {
        "schema": "mceliece-postrun-strict-parent-replay-v1",
        "status": "PASS" if all_equal and complete_equal_true else "FAIL",
        "timing": "POST-RUN CORRECTIVE AUDIT; not a pre-production control",
        "removed_fields": ["header", "t_utc"],
        "parent_ledger_sha256": EXPECTED_PARENT_LEDGER_SHA256,
        "parent_ledger_rows_verified": parent_rows,
        "successor_ledger_sha256": EXPECTED_SUCCESSOR_LEDGER_SHA256,
        "successor_ledger_rows_verified": successor_rows,
        "betas": EXPECTED_BETAS,
        "all_strict_semantic_equal": all_equal,
        "all_complete_equal_true": complete_equal_true,
        "comparisons": comparisons,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
