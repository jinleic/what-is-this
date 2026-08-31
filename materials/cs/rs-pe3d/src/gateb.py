"""Gate B: fixed-(q, eta) balance sweep over subgroup triples.

Pre-registration: campaigns/2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB/
pre_statement.md, committed BEFORE this code ran. beta(s) = max(s)/min(s).

Phase 1 (anchors, must pass before anything else):
  - unit anchor (5,(2,2,2),(1,1,1)) exhaustive -> exact rho = 1/3;
  - B01 (61,(2,3,5),(1,1,1)) -> 3/5; B02 (61,(3,4,5),(1,1,1)) -> 1/1;
Phase 2: B03-B09 frozen gate-A value reproduction (second route).
Phase 3: B10-B13, the 241 wide-beta tier.

Window per pre-statement GB3/GB4: weight<=3 point-support census over the
complete support range, candidates = normalized reduced basis of each nonzero
V n F^S, per-witness delta closed by the exhaustive cheaper-line-tuple sweep
(frozen optcheck._proof_delta procedure). Checkpoint after every instance.
"""
from __future__ import annotations

import datetime as dt
import itertools
import json
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

from src.campaign import _code_hash, _normalize, _profile_census
from src.delta import DeltaEngine
from src.milp import delta_milp
from src.rs import Inst
from src.verify import verify_record

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "campaigns" / "2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB"
CHECKPOINT = ROOT / "scratch" / "gateB_checkpoint.json"
MILP_LIMIT = 120.0

INSTANCES = [
    # id, q, s, t, frozen expected ratio or None, tier
    {"id": "B01", "q": 61, "s": (2, 3, 5), "t": (1, 1, 1), "expect": "3/5", "tier": 1},
    {"id": "B02", "q": 61, "s": (3, 4, 5), "t": (1, 1, 1), "expect": "1/1", "tier": 1},
    {"id": "B03", "q": 31, "s": (2, 3, 5), "t": (1, 1, 1), "expect": "3/5", "tier": 2},
    {"id": "B04", "q": 31, "s": (2, 3, 10), "t": (1, 1, 1), "expect": "3/5", "tier": 2},
    {"id": "B05", "q": 31, "s": (2, 6, 5), "t": (1, 1, 1), "expect": "1/1", "tier": 2},
    {"id": "B06", "q": 31, "s": (3, 10, 2), "t": (1, 1, 1), "expect": "3/5", "tier": 2},
    {"id": "B07", "q": 13, "s": (2, 2, 4), "t": (1, 1, 1), "expect": "1/2", "tier": 2},
    {"id": "B08", "q": 13, "s": (2, 2, 4), "t": (1, 1, 2), "expect": "3/8", "tier": 2},
    {"id": "B09", "q": 13, "s": (2, 2, 4), "t": (1, 1, 4), "expect": "1/4", "tier": 2},
    {"id": "B10", "q": 241, "s": (3, 5, 8), "t": (1, 1, 1), "expect": None, "tier": 3},
    {"id": "B11", "q": 241, "s": (3, 5, 16), "t": (1, 1, 1), "expect": None, "tier": 3},
    {"id": "B12", "q": 241, "s": (2, 3, 5), "t": (1, 1, 1), "expect": None, "tier": 3},
    {"id": "B13", "q": 241, "s": (3, 4, 5), "t": (1, 1, 1), "expect": None, "tier": 3},
]


def beta_of(s: tuple[int, int, int]) -> Fraction:
    return Fraction(max(s), min(s))


def _checkpoint(state: dict) -> None:
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    tmp = CHECKPOINT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
    tmp.replace(CHECKPOINT)


def unit_anchor() -> dict:
    """Exhaustive (5,(2,2,2),(1,1,1)) -> exact rho. Reuses units anchor-2 DP."""
    from src.units import anchor2

    t0 = time.monotonic()
    res = anchor2()
    res["seconds"] = round(time.monotonic() - t0, 3)
    ok = res["rho_exact"] == "1/3" and res["V_nonzero_count"] == 5 ** 7 - 1
    res["anchor_ok"] = ok
    return res


def proof_delta(engine: DeltaEngine, M: list[int], best_cost: int) -> dict:
    """Exact per-witness closure: all tuples cheaper than best_cost absent,
    >= 1 tuple at best_cost present. Identical logic to frozen optcheck."""
    inst = engine.inst
    n_lines = [inst.num_lines(i) for i in range(3)]
    cheaper_hit = None
    equal_hits = 0
    for cost, n0, n1, n2 in engine.comps:
        if cost >= best_cost:
            break
        hit = _sweep_cost(engine, M, n0, n1, n2, n_lines)
        if hit:
            cheaper_hit = cost
            break
    if cheaper_hit is None:
        for cost, n0, n1, n2 in engine.comps:
            if cost != best_cost:
                continue
            if _sweep_cost(engine, M, n0, n1, n2, n_lines):
                equal_hits += 1
    return {
        "cheaper_hit": cheaper_hit,
        "equal_hits": equal_hits,
        "proved": cheaper_hit is None and equal_hits > 0,
    }


def _sweep_cost(engine: DeltaEngine, M, n0, n1, n2, n_lines) -> bool:
    inst = engine.inst
    for B0 in itertools.combinations(range(n_lines[0]), n0):
        rows0 = [r for l in B0 for r in engine.line_rows[0][l]]
        for B1 in itertools.combinations(range(n_lines[1]), n1):
            rows01 = rows0 + [r for l in B1 for r in engine.line_rows[1][l]] if n1 else rows0
            for B2 in itertools.combinations(range(n_lines[2]), n2):
                rows = rows01 + [r for l in B2 for r in engine.line_rows[2][l]] if n2 else rows01
                if engine._membership(rows, M):
                    return True
    return False


def run_instance(spec: dict) -> dict:
    inst = Inst(spec["q"], spec["s"], spec["t"])
    start = time.monotonic()
    engine = DeltaEngine(inst)
    census, candidates = _profile_census(inst, engine, 3)
    best = None
    witnesses = []
    for M, info in sorted(candidates.items()):
        result = delta_milp(inst, list(M), time_limit=MILP_LIMIT)
        if result.status not in {"kOptimal", "kModelOptimal"}:
            continue
        wt = sum(x != 0 for x in M)
        ratio = Fraction(wt, result.objective)
        if best is None or ratio < best[0]:
            record = {
                "id": spec["id"],
                "q": inst.q,
                "s": list(inst.s),
                "t": list(inst.t),
                "lam": [list(x) for x in inst.lam],
                "M": list(M),
                "weight": wt,
                "delta": result.objective,
                "components": [list(x) for x in result.components],
                "coefficients": [
                    [list(c) for c in rows_] for rows_ in result.coefficients
                ],
            }
            best = (ratio, wt, result.objective, record)
    if best is None:
        raise AssertionError(f"no optimal candidate rows for {spec['id']}")
    ratio, wt, delta, record = best
    checked = verify_record(record, check_lower=False)
    if not checked["ok"]:
        raise AssertionError((spec["id"], checked))
    # Exact per-witness closure at the best witness (this is the more precise
    # gate-A optimality step; verify_record's check_lower cannot bound above).
    engine2 = DeltaEngine(inst)
    if not hasattr(engine2, "rank_cache"):
        engine2.rank_cache = {}
    proof = proof_delta(engine2, record["M"], delta)
    row = {
        "id": spec["id"],
        "q": inst.q,
        "s": list(inst.s),
        "t": list(inst.t),
        "N": inst.N,
        "V_dim": engine.V_dim(),
        "beta": f"{beta_of(spec['s']).numerator}/{beta_of(spec['s']).denominator}",
        "support_window": "all point supports of weight <= 3",
        "candidate_class": "normalized reduced basis of each V n F^S",
        "candidate_count_total": len(candidates),
        "ratio": f"{ratio.numerator}/{ratio.denominator}",
        "weight": wt,
        "delta": delta,
        "attained_in_window": len(record["M"]) <= 3 and wt <= 3,
        "delta_proven_exact_for_best_witness": proof["proved"],
        "proof_cheaper_hit": proof["cheaper_hit"],
        "proof_equal_hits": proof["equal_hits"],
        "evidence_label": (
            "MACHINE-VERIFIED (finite window)" if proof["proved"]
            else "COMPUTATIONAL-EVIDENCE"
        ),
        "expect_match": None,
        "seconds": round(time.monotonic() - start, 3),
    }
    if spec["expect"] is not None:
        row["expect_match"] = row["ratio"] == spec["expect"]
    return row, record


def main() -> int:
    started = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    state = {
        "campaign": CAMPAIGN.name,
        "started_utc": started.isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "code_sha256": _code_hash(),
        "phase": "unit-anchor",
        "rows": [],
        "unit_anchor": None,
    }
    _checkpoint(state)
    print("checkpoint written (pre-anchor)", flush=True)

    # ---- Phase 1a: unit anchor ----
    unit = unit_anchor()
    state["unit_anchor"] = unit
    state["phase"] = "unit-anchor-done"
    _checkpoint(state)
    print(f"unit anchor rho={unit['rho_exact']} ok={unit['anchor_ok']}", flush=True)
    if not unit["anchor_ok"]:
        state["phase"] = "UNIT-ANCHOR-FAILED"
        _checkpoint(state)
        print("UNIT ANCHOR FAILED - stopping per pre-statement GB5", flush=True)
        return 1

    # ---- Phase 1b + 2 + 3: instance rows in committed order ----
    for spec in INSTANCES:
        t0 = time.monotonic()
        try:
            row, record = run_instance(spec)
            state["rows"].append(row)
            state["last_records"] = state.get("last_records", [])[-2:] + [record]
            state["phase"] = f"done-{spec['id']}"
            _checkpoint(state)
            print(
                f"{spec['id']} q={spec['q']} s={spec['s']} beta={row['beta']} "
                f"ratio={row['ratio']} wt={row['weight']} delta={row['delta']} "
                f"proved={row['delta_proven_exact_for_best_witness']} "
                f"expect={row['expect_match']} ({row['seconds']}s)",
                flush=True,
            )
            if row["expect_match"] is False:
                state["phase"] = f"ANCHOR-MISMATCH-{spec['id']}"
                state["anchor_mismatch"] = {
                    "id": spec["id"],
                    "expected": spec["expect"],
                    "got": row["ratio"],
                }
                _checkpoint(state)
                print(
                    f"ANCHOR MISMATCH at {spec['id']}: expected {spec['expect']}, "
                    f"got {row['ratio']} - STOPPING, escalation to Main required",
                    flush=True,
                )
                return 2
        except Exception as exc:  # checkpoint the failure, then re-raise
            state["phase"] = f"FAILED-{spec['id']}"
            state.setdefault("failures", []).append(
                {"id": spec.get("id"), "error": repr(exc), "at": dt.datetime.now(dt.timezone.utc).isoformat()}
            )
            _checkpoint(state)
            raise

    state["phase"] = "complete"
    state["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    _checkpoint(state)
    print("all rows complete", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
