"""Independent exact verifier for frozen RS witness records.

Usage:
    python -m src.verify campaigns/<run>/witnesses.json

The verifier reconstructs the instance, evaluates every supplied line
polynomial in F_q, checks the sum and weight, and optionally exhausts every
line-support pattern cheaper than the claimed delta.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from .delta import DeltaEngine
from .rs import Inst


def _instance(record: dict) -> Inst:
    return Inst(
        int(record["q"]),
        tuple(int(x) for x in record["s"]),
        tuple(int(x) for x in record["t"]),
        record.get("lam"),
    )


def evaluate_components(inst: Inst, coefficients) -> tuple[list[int], list[list[int]]]:
    """Return sum_i M_i and each evaluated component as flat F_q vectors."""
    q = inst.q
    total = [0] * inst.N
    evaluated: list[list[int]] = []
    for i in range(3):
        rows = coefficients[i]
        if len(rows) != inst.num_lines(i):
            raise ValueError(
                f"direction {i}: coefficient rows {len(rows)} != "
                f"{inst.num_lines(i)}"
            )
        comp = [0] * inst.N
        for line, coeff in enumerate(rows):
            if len(coeff) != inst.t[i]:
                raise ValueError(
                    f"direction {i}, line {line}: coefficient length "
                    f"{len(coeff)} != {inst.t[i]}"
                )
            idxs = inst.line_indices(i, line)
            for ai, p in enumerate(idxs):
                value = sum(
                    int(coeff[d]) * inst.lam[i][ai] * pow(inst.S[i][ai], d, q)
                    for d in range(inst.t[i])
                ) % q
                comp[p] = value
                total[p] = (total[p] + value) % q
        evaluated.append(comp)
    return total, evaluated


def verify_record(record: dict, check_lower: bool = True) -> dict:
    inst = _instance(record)
    q = inst.q
    M = [int(x) % q for x in record["M"]]
    if len(M) != inst.N:
        raise ValueError(f"M length {len(M)} != N={inst.N}")
    coefficients = record["coefficients"]
    components = tuple(tuple(int(x) for x in xs) for xs in record["components"])
    if len(components) != 3:
        raise ValueError("components must have three directions")

    total, evaluated = evaluate_components(inst, coefficients)
    if total != M:
        bad = next(i for i, (a, b) in enumerate(zip(total, M)) if a != b)
        raise ValueError(f"F_q sum mismatch at flat point {bad}: {total[bad]} != {M[bad]}")

    actual = []
    for i in range(3):
        actual_i = tuple(
            line for line in range(inst.num_lines(i))
            if any(evaluated[i][p] for p in inst.line_indices(i, line))
        )
        actual.append(actual_i)
        if actual_i != components[i]:
            raise ValueError(
                f"direction {i}: supplied lines {components[i]} != "
                f"evaluated nonzero lines {actual_i}"
            )

    wt = sum(x != 0 for x in M)
    expected_wt = int(record.get("weight", wt))
    if wt != expected_wt:
        raise ValueError(f"weight {wt} != recorded {expected_wt}")
    delta = int(record["delta"])
    recomputed_cost = sum(inst.s[i] * len(actual[i]) for i in range(3))
    if recomputed_cost != delta:
        raise ValueError(f"delta {delta} != support cost {recomputed_cost}")

    lower_checked = False
    lower_result = None
    if check_lower and delta > 0:
        engine = DeltaEngine(inst)
        lower_result = engine.delta_exact(M, max_cost=delta - 1)
        lower_checked = True
        if lower_result[0] is not None:
            raise ValueError(
                f"cheaper exact decomposition found: cost={lower_result[0]} "
                f"below claimed delta={delta}"
            )

    return {
        "ok": True,
        "q": q,
        "s": list(inst.s),
        "t": list(inst.t),
        "N": inst.N,
        "weight": wt,
        "delta": delta,
        "ratio": f"{wt}/{delta}",
        "components": [list(x) for x in actual],
        "lower_support_exhaustive": lower_checked,
        "lower_support_result": None if lower_result is None else lower_result[0],
    }


def verify_path(path: str | Path, check_lower: bool = True) -> dict:
    payload = json.loads(Path(path).read_text())
    if isinstance(payload, dict):
        records = payload.get("records", [payload])
    else:
        records = payload
    if not isinstance(records, list):
        raise ValueError("expected a record or a JSON object with records")
    results = [verify_record(record, check_lower=check_lower) for record in records]
    return {"ok": True, "count": len(results), "records": results}


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) not in (1, 2) or (len(argv) == 2 and argv[1] != "--no-lower-check"):
        print("usage: python -m src.verify FILE [--no-lower-check]", file=sys.stderr)
        return 2
    try:
        out = verify_path(argv[0], check_lower=(len(argv) == 1))
    except Exception as exc:  # CLI emits a machine-readable failure line.
        print(json.dumps({"ok": False, "error": str(exc)}, sort_keys=True))
        return 1
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
