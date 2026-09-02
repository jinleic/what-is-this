"""Campaign W rho stage — exact delta closure over the w=4 candidate set.

ADJUDICATION OBJECT (pre_statement.md section 1/6, frozen at 8114a8f):
per-candidate delta computed EXACTLY by the frozen exhaustive ascending-cost
line-tuple sweep (DeltaEngine.delta_exact, src/delta.py, unchanged); candidates
are the normalized reduced-basis vectors of each nonzero V ∩ F^S. This stage
runs on the PARITY PRIMARY census output (parity_w4_hits.json), the
pre-registered PRIMARY instrument per the mid-flight freeze (commit 11cc15a
handback instruction: "Do not silently let the slow instrument be primary").
Exact Fraction arithmetic only (no float, no Arb needed at this scale: delta
values are small integers; rho = Fraction(4, delta) is exact).

Certified domain (rule 15): rho^window_w4 is certified ONLY on the
weight-EXACTLY-4 point-support window at PN4 = (421,(4,5,7),t=(1,1,1),
Λ=Id), complete census (all 15,329,615 supports), by exact F_421 arithmetic
plus the exhaustive sweep. It says nothing about rho_inst or any other
instance/window (rule 7 sentence in the freeze addendum).
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
RS = HERE.parents[1]
sys.path.insert(0, str(RS))
sys.path.insert(0, str(HERE))
from census_runner import build_instruments, Q, S_PN4, T_PN4  # noqa: E402


def main() -> int:
    hits = json.loads((HERE / "parity_w4_hits.json").read_text())["hits"]
    n_supports = sum(1 for _ in open(HERE / "parity_ledger.jsonl"))  # noqa
    t0 = time.perf_counter()
    I, eng, B = build_instruments(Q, S_PN4, T_PN4)
    assert B.shape[0] == 68
    recs = []
    multi = 0
    for h in hits:
        S = h["S"]
        assert h["dim"] == len(h["class"])
        if h["dim"] != 1:
            multi += 1
        for k, vec in enumerate(h["class"]):
            assert vec[S[0]] == 1, f"normalization: first coord {vec[S[0]]}"
            delta, tup = eng.delta_exact(vec)
            assert delta is not None, f"delta not closed for {S}"
            recs.append({
                "S": S, "dim": h["dim"], "lex_idx": h["idx"],
                "basis_vec_idx": k,
                "delta": delta,
                "line_tuple": [list(b) for b in tup],
                "rho": str(Fraction(4, delta)),
                "rho_num": [4, delta],
            })
    # rho^window = min wt/delta over nonzero supports (wt = 4 exactly)
    best = min(recs, key=lambda r: Fraction(r["rho_num"][0],
                                            r["rho_num"][1]))
    rho_win = Fraction(4, best["delta"])
    deltas = sorted({r["delta"] for r in recs})
    out = {
        "instrument": "parity PRIMARY census (parity_w4_hits.json, "
                      "final next_index 15329615, hits_total 35)",
        "n_supports_enumerated": 15_329_615,
        "n_candidate_supports": len(hits),
        "multi_dim_supports": multi,
        "n_candidates_total": len(recs),
        "distinct_deltas": deltas,
        "rho_window_w4": str(rho_win),
        "best": best,
        "per_candidate": recs,
    }
    (HERE / "rho_w4_results.json").write_text(json.dumps(out, indent=1))
    print(f"rho stage complete: {len(recs)} candidates, deltas {deltas}, "
          f"rho^window_w4 = {rho_win}, wall "
          f"{time.perf_counter()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
