"""Turn best-in-window rows into exact per-witness optimality proofs.

For each frozen best witness this exhaustively enumerates EVERY line-support
tuple B=(B1,B2,B3) whose weighted cost is <= best_delta and proves exact
membership (F_q row reduction) fails everywhere cheaper and succeeds only at
the claimed minimum cost. This closes delta(M) exactly for that witness.
Combined with the support census the row is then:
  - delta(witness) EXACT (complete line-tuple proof);
  - ratio best among wt<=3 supports (complete in that window).
Heavy-support global minimization remains out of scope and is labeled OPEN
per row.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

from .delta import DeltaEngine
from .rs import Inst
from .verify import verify_record

ROOT = Path(__file__).resolve().parents[1]


def _enumerate_line_tuples(engine, best_cost: int, max_cost: int):
    """All (B0,B1,B2) with cost <= max_cost, ascending, caches membership."""
    inst = engine.inst
    s = inst.s
    n_lines = [inst.num_lines(i) for i in range(3)]
    for n0 in range(n_lines[0] + 1):
        base0 = s[0] * n0
        if base0 > max_cost:
            break
        for n1 in range(n_lines[1] + 1):
            base1 = base0 + s[1] * n1
            if base1 > max_cost:
                break
            for n2 in range(n_lines[2] + 1):
                cost = base1 + s[2] * n2
                if cost > max_cost:
                    break
                for B0 in itertools.combinations(range(n_lines[0]), n0):
                    rows0 = [r for l in B0 for r in engine.line_rows[0][l]]
                    rank0 = engine.rank_cache.get(tuple(B0))
                    if rank0 is None:
                        from .field import rank_mod

                        rank0 = rank_mod(rows0, inst.q)
                        engine.rank_cache[tuple(B0)] = rank0
                    if rank0 == n0 * inst.t[0]:
                        # full, skip
                        pass  # not pruning; keep semantics simple
                    yield cost, (B0, n1, n2), rows0


def _proof_delta(engine, M: list[int], best_cost: int) -> dict:
    """Exactness proof for delta(M)==best_cost via complete sweep to best_cost."""
    inst = engine.inst
    from .field import rref, rank_mod

    n_lines = [inst.num_lines(i) for i in range(3)]
    s = inst.s
    cheaper_hit = None
    equal_hits = 0
    # all tuples with cost < best_cost must fail
    for cost, n0, n1, n2 in engine.comps:
        if cost >= best_cost:
            break
        for B0 in itertools.combinations(range(n_lines[0]), n0):
            rows0 = [r for l in B0 for r in engine.line_rows[0][l]]
            for B1 in itertools.combinations(range(n_lines[1]), n1):
                rows01 = (
                    rows0 + [r for l in B1 for r in engine.line_rows[1][l]]
                    if n1
                    else rows0
                )
                for B2 in itertools.combinations(range(n_lines[2]), n2):
                    rows = (
                        rows01 + [r for l in B2 for r in engine.line_rows[2][l]]
                        if n2
                        else rows01
                    )
                    if engine._membership(rows, M):
                        cheaper_hit = (cost, B0, B1, B2)
                        break
                if cheaper_hit:
                    break
            if cheaper_hit:
                break
        if cheaper_hit:
            break
    # all tuples with cost == best_cost: count successes (>=1 expected)
    if cheaper_hit is None:
        for cost, n0, n1, n2 in engine.comps:
            if cost != best_cost:
                continue
            for B0 in itertools.combinations(range(n_lines[0]), n0):
                rows0 = [r for l in B0 for r in engine.line_rows[0][l]]
                for B1 in itertools.combinations(range(n_lines[1]), n1):
                    rows01 = (
                        rows0 + [r for l in B1 for r in engine.line_rows[1][l]]
                        if n1
                        else rows0
                    )
                    for B2 in itertools.combinations(range(n_lines[2]), n2):
                        rows = (
                            rows01 + [r for l in B2 for r in engine.line_rows[2][l]]
                            if n2
                            else rows01
                        )
                        if engine._membership(rows, M):
                            equal_hits += 1
    return {
        "cheaper_hit": cheaper_hit,
        "equal_hits": equal_hits,
        "proved": cheaper_hit is None and equal_hits > 0,
    }


def run(path: str | None = None) -> Path:
    import sys

    payload_path = Path(
        path
        or (ROOT / "campaigns" / "2026-08-30T01-04-49Z_exacttable" / "payload.json")
    )
    payload = json.loads(payload_path.read_text())
    records = json.loads(
        (payload_path.parent / "witnesses.json").read_text()
    )["records"]

    out_rows = []
    t0 = time.monotonic()
    for row in payload["rows"]:
        inst = Inst(row["q"], tuple(row["s"]), tuple(row["t"]))
        engine = DeltaEngine(inst)
        if not hasattr(engine, "rank_cache"):
            engine.rank_cache = {}
        best_witnesses = [
            r
            for r in records
            if r["q"] == row["q"]
            and tuple(r["s"]) == tuple(row["s"])
            and tuple(r["t"]) == tuple(row["t"])
            and r["weight"] == row["weight"]
            and r["delta"] == row["delta"]
        ]
        witness = best_witnesses[0]
        proof = _proof_delta(engine, witness["M"], witness["delta"])
        verdict = {
            **{k: row[k] for k in ("q", "s", "t", "N", "V_dim", "ratio", "weight", "delta")},
            "delta_proven_exact_for_best_witness": proof["proved"],
            "proof_cheaper_hit": proof["cheaper_hit"],
            "proof_equal_hits": proof["equal_hits"],
            "scope": "per-witness exact delta; global rho over heavy supports OPEN",
            "evidence_label": (
                "MACHINE-VERIFIED" if proof["proved"] else row["evidence_label"]
            ),
        }
        out_rows.append(verdict)
        print(json.dumps(verdict), flush=True)
    out = payload_path.parent / "optimality.json"
    out.write_text(json.dumps(out_rows, indent=2, sort_keys=True) + "\n")
    return out


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    print(run(arg))
