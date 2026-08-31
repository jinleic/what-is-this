"""Pre-campaign unit anchors for the RS line-sum implementation.

Anchor 2 is exhaustive on q=5,s=(2,2,2),t=(1,1,1): all 5^7 nonzero
vectors in V are enumerated by a direct component-array product, and the
minimum wt/delta is computed exactly. A deterministic sample is cross-checked
against DeltaEngine and the HiGHS modular model.

Anchor 1 is run on the nearest feasible subgroup toy q=13,s=(2,3,4),t=s.
It records a line-sum equality that I derived as a unit hypothesis and
fails it with an independently F_q-verified split certificate; this is a
scratch/unit premise failure, not a campaign or paper claim.
"""
from __future__ import annotations

import datetime as dt
import itertools
import json
import random
from fractions import Fraction
from pathlib import Path

import numpy as np

from .delta import DeltaEngine
from .milp import delta_milp
from .rs import Inst
from .verify import verify_record

ROOT = Path(__file__).resolve().parents[1]


def _component_arrays(inst: Inst):
    q = inst.q
    out = []
    for i in range(3):
        rows = []
        for vals in itertools.product(range(q), repeat=inst.num_lines(i)):
            vec = np.zeros(inst.N, dtype=np.int16)
            for line, value in enumerate(vals):
                if value:
                    for p in inst.line_indices(i, line):
                        vec[p] = value
            rows.append((vec, inst.s[i] * sum(x != 0 for x in vals)))
        out.append(rows)
    return out


def _decode(code: int, q: int, n: int) -> list[int]:
    return [(code // (q**j)) % q for j in range(n)]


def anchor2() -> dict:
    inst = Inst(5, (2, 2, 2), (1, 1, 1))
    arrays = _component_arrays(inst)
    A = np.stack([x[0] for x in arrays[0]])
    B = np.stack([x[0] for x in arrays[1]])
    C = arrays[2]
    AB = ((A[:, None, :] + B[None, :, :]) % inst.q).reshape(-1, inst.N)
    cost_ab = np.array(
        [a[1] + b[1] for a in arrays[0] for b in arrays[1]], dtype=np.int16
    )
    powers = (inst.q ** np.arange(inst.N)).astype(np.int64)
    best = np.full(inst.q**inst.N, 10**9, dtype=np.int32)
    for cvec, ccost in C:
        codes = ((AB + cvec) % inst.q).astype(np.int64) @ powers
        np.minimum.at(best, codes, cost_ab + ccost)

    reachable = np.flatnonzero(best < 10**9)
    expected = inst.q**7
    if len(reachable) != expected:
        raise AssertionError(f"reachable={len(reachable)} expected={expected}")
    ratios: list[tuple[Fraction, int, int]] = []
    for code in reachable:
        if code == 0:
            continue
        digits = _decode(int(code), inst.q, inst.N)
        weight = sum(x != 0 for x in digits)
        ratios.append((Fraction(weight, int(best[code])), weight, int(code)))
    min_ratio, min_weight, min_code = min(ratios, key=lambda x: x[0])
    witness = _decode(min_code, inst.q, inst.N)

    # Independent line-support enumeration and the modular MILP are sampled
    # at deterministic codes; the full exact rho remains the DP result above.
    rng = random.Random(0)
    sample_codes = [int(x) for x in rng.sample(list(reachable[1:]), 24)]
    sample_checks = []
    engine = DeltaEngine(inst)
    for code in sample_codes:
        M = _decode(code, inst.q, inst.N)
        brute_delta = int(best[code])
        exact_delta, _exact_lines = engine.delta_exact(M)
        if exact_delta != brute_delta:
            raise AssertionError((code, brute_delta, exact_delta))
        milp = delta_milp(inst, M, time_limit=20)
        if milp.objective != brute_delta or milp.components is None:
            raise AssertionError((code, brute_delta, milp.status, milp.objective))
        rec = {
            "q": inst.q,
            "s": list(inst.s),
            "t": list(inst.t),
            "lam": [list(x) for x in inst.lam],
            "M": M,
            "weight": sum(x != 0 for x in M),
            "delta": brute_delta,
            "components": [list(x) for x in milp.components],
            "coefficients": [
                [list(c) for c in rows] for rows in milp.coefficients
            ],
        }
        verify_record(rec, check_lower=True)
        sample_checks.append({"code": code, "delta": brute_delta})

    return {
        "anchor": "2",
        "status": "MACHINE-VERIFIED",
        "q": inst.q,
        "s": list(inst.s),
        "t": list(inst.t),
        "V_nonzero_count": len(reachable) - 1,
        "V_expected_nonzero_count": expected - 1,
        "rho_exact": f"{min_ratio.numerator}/{min_ratio.denominator}",
        "witness": witness,
        "witness_weight": min_weight,
        "witness_delta": int(best[min_code]),
        "sample_cross_checks": len(sample_checks),
    }


def anchor1() -> dict:
    # q=7,s=(2,3,4) cannot have S_3 of order 4; q=13 is the fixed feasible
    # substitution recorded in pre_statement.md.
    inst = Inst(13, (2, 3, 4), (2, 3, 4))
    rng = random.Random(0)
    mismatch = None
    for sample in range(20):
        M = [rng.randrange(inst.q) for _ in range(inst.N)]
        direct = min(
            inst.s[i]
            * sum(
                any(M[p] % inst.q for p in inst.line_indices(i, line))
                for line in range(inst.num_lines(i))
            )
            for i in range(3)
        )
        result = delta_milp(inst, M, time_limit=20)
        if result.status not in {"kOptimal", "kModelOptimal"}:
            raise AssertionError(result)
        if result.objective < direct:
            rec = {
                "q": inst.q,
                "s": list(inst.s),
                "t": list(inst.t),
                "lam": [list(x) for x in inst.lam],
                "M": M,
                "weight": sum(x != 0 for x in M),
                "delta": result.objective,
                "components": [list(x) for x in result.components],
                "coefficients": [
                    [list(c) for c in rows] for rows in result.coefficients
                ],
            }
            checked = verify_record(rec, check_lower=False)
            mismatch = {
                "sample": sample,
                "direct_one_direction_cost": direct,
                "split_cost": result.objective,
                "verified": checked,
                "record": rec,
            }
            break
    if mismatch is None:
        raise AssertionError("anchor-1 mismatch was not reproduced")
    return {
        "anchor": "1",
        "status": "FAILED-ASSERTION",
        "reason": "the asserted one-direction equality is false",
        "feasible_substitution": "q=13,s=(2,3,4),t=s",
        "counterexample": mismatch,
    }


def run() -> Path:
    result = {
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "anchor1": anchor1(),
        "anchor2": anchor2(),
        "evidence_note": "Anchor 2 exhaustive; Anchor 1 failure is a unit-premise finding",
    }
    out = ROOT / "campaigns-smoke" / "unit-anchors.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return out


if __name__ == "__main__":
    print(run())
