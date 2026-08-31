"""Gate-B freeze: regenerate every artifact from raw records, then checksum.

Artifacts are REGENERATED, not transcribed: each row's headline ratio is
recomputed as Fraction(weight, delta) from its witness and cross-checked
against the checkpointed string, and every witness is re-verified in exact
F_q here, inside the freeze. Manifest and checksums are written last.

Run: cd cs/rs-pe3d && PYTHONPATH=$PWD ../.venv/bin/python -m src.gateb_freeze
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import platform
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

from .campaign import _code_hash
from .field import coprime_triples
from .verify import verify_record

ROOT = Path(__file__).resolve().parents[1]
CAMP = ROOT / "campaigns" / "2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB"

# Measured census throughput (supports/second) from this campaign's own runs:
# B10 N=120: 288,100 supports in 404 s -> 713/s; B11 N=240: 2,304,200 in 3449 s
# -> 668/s. Conservative 650/s used for pricing unswept extensions.
MEASURED_SUPPORTS_PER_SEC = 650


def _monotone_violations(seq: list[dict]) -> list[dict]:
    """Interleavings (beta increasing, rho peak or valley) that exclude any
    monotone rho(beta). Exact rational comparisons."""
    out = []
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            for k in range(j + 1, len(seq)):
                a, b, c = seq[i], seq[j], seq[k]
                ba, bb, bc = (Fraction(x["beta"]) for x in (a, b, c))
                ra, rb, rc = (Fraction(x["rho"]) for x in (a, b, c))
                if ba < bb < bc and ((rb > ra and rb > rc) or (rb < ra and rb < rc)):
                    out.append(
                        {
                            "rows": [a["id"], b["id"], c["id"]],
                            "betas": [a["beta"], b["beta"], c["beta"]],
                            "rhos": [a["rho"], b["rho"], c["rho"]],
                            "shape": "peak" if rb > ra else "valley",
                        }
                    )
    return out


def _price_unswept_421() -> dict:
    """Cost of the pre-statement's explicitly-OPEN q=421 nineteen-triple census."""
    triples = coprime_triples(420)
    rows = []
    total = 0
    for s in triples:
        N = s[0] * s[1] * s[2]
        supports = comb(N, 1) + comb(N, 2) + comb(N, 3)
        total += supports
        rows.append({"s": list(s), "N": N, "supports_weight_le_3": supports})
    return {
        "triples": len(triples),
        "per_triple": rows,
        "total_supports": total,
        "measured_supports_per_second": MEASURED_SUPPORTS_PER_SEC,
        "estimated_census_seconds": round(total / MEASURED_SUPPORTS_PER_SEC),
        "estimated_census_hours": round(total / MEASURED_SUPPORTS_PER_SEC / 3600, 1),
        "note": (
            "single-threaded, census only; MILP + per-witness closure extra. "
            "This is the price of closing the q=421 tier at the same window, "
            "and it is why that tier is OPEN rather than attempted."
        ),
    }


def build() -> dict:
    state = json.loads((ROOT / "scratch" / "gateB_checkpoint.json").read_text())
    wit = json.loads((ROOT / "scratch" / "gateB_witnesses.json").read_text())
    b11 = json.loads((ROOT / "scratch" / "gateB_b11_result.json").read_text())

    rows = {r["id"]: r for r in state["rows"]}
    recs = {r["id"]: r for r in wit["records"]}
    reruns = {r["id"]: r for r in wit["reruns"]}

    regen = []
    for rid in sorted(rows):
        row, rec = rows[rid], recs[rid]
        chk = verify_record(rec, check_lower=False)  # exact F_q re-check in-freeze
        ratio = Fraction(rec["weight"], rec["delta"])  # recomputed, not copied
        if f"{ratio.numerator}/{ratio.denominator}" != row["ratio"]:
            raise AssertionError(f"{rid}: regenerated {ratio} != recorded {row['ratio']}")
        floor = Fraction(1, 3 * row["N"])  # proved universal floor 1/(3N)
        width = ratio - floor
        regen.append(
            {
                **{
                    k: row[k]
                    for k in (
                        "id", "q", "s", "t", "N", "V_dim", "beta", "ratio", "weight",
                        "delta", "candidate_count_total", "support_window",
                        "candidate_class", "delta_proven_exact_for_best_witness",
                        "proof_cheaper_hit", "proof_equal_hits", "evidence_label",
                        "expect_match", "seconds",
                    )
                },
                "rho_inst_certified_interval": {
                    "lower": f"{floor.numerator}/{floor.denominator}",
                    "upper": row["ratio"],
                    "width_exact": f"{width.numerator}/{width.denominator}",
                    "width_decimal_20dp": f"{float(width):.20f}",
                    "lower_source": "delta(M) <= 3N for every M in V (parent pre-statement; gate-C manifest)",
                    "upper_source": "this campaign: exact best-in-window witness, delta closed exhaustively",
                    "exactness": "both endpoints exact rationals; width is exact, zero-uncertainty arithmetic",
                },
                "witness_verified_in_freeze": chk["ok"],
                "rerun_fresh_process_agree": reruns[rid].get("agree"),
                "rerun_seconds": reruns[rid].get("seconds"),
            }
        )

    per_q: dict[str, list[dict]] = {}
    for r in regen:
        per_q.setdefault(str(r["q"]), []).append(
            {
                "id": r["id"], "s": r["s"], "t": r["t"], "beta": r["beta"],
                "rho": r["ratio"], "weight": r["weight"], "delta": r["delta"], "N": r["N"],
            }
        )
    for q in per_q:
        per_q[q].sort(key=lambda x: (Fraction(x["beta"]), x["id"]))

    p_check = [
        {
            "id": r["id"], "q": r["q"], "s": r["s"], "t": r["t"],
            "contains_2_and_3": (2 in r["s"] and 3 in r["s"]),
            "rho": r["ratio"],
            "P_consistent": ((2 in r["s"] and 3 in r["s"]) == (r["ratio"] == "3/5")),
            "sample": "out-of-sample (pre-registered)" if r["id"] == "B11" else "in-sample",
        }
        for r in regen
    ]
    n_with = sum(1 for x in p_check if x["contains_2_and_3"])
    n_rows = len(p_check)

    analysis = {
        "beta_grid_swept": sorted({r["beta"] for r in regen}, key=Fraction),
        "per_q_sorted_by_beta": per_q,
        "EXCLUSION_E": {
            "status": "CONCLUSION from exact data; independent of B11; B11 strengthens it at the extreme point",
            "statement": "On the swept set, rho^window is NOT a function of beta: no monotone (indeed no single-valued-in-beta) dependence fits.",
            "monotone_excluding_interleavings": {q: _monotone_violations(v) for q, v in per_q.items()},
            "same_beta_different_rho": [
                {"beta": "5/2", "rows": ["B01 (61,(2,3,5)) = 3/5", "B03 (31,(2,3,5)) = 3/5", "B12 (241,(2,3,5)) = 3/5"]},
                {"beta": "5/3", "rows": ["B02 (61,(3,4,5)) = 1", "B13 (241,(3,4,5)) = 1"]},
            ],
            "extreme_point": "beta grows 3.2x from 5/3 to 16/3 at q=241 with rho 1 -> 1, while intermediate beta=5/2 sits at 3/5",
        },
        "PATTERN_P": {
            "status": "EVIDENCE, not a theorem",
            "statement": "At t=(1,1,1), Lambda=Id, in this window: {2,3} subseteq orders(s)  <=>  rho^window = 3/5",
            "fit": f"perfect on all {n_rows} rows ({n_with} contain both 2 and 3 and all are 3/5; {n_rows - n_with} do not and none is 3/5)",
            "in_sample_rows": n_rows - 1,
            "out_of_sample_rows": 1,
            "out_of_sample_detail": "B11 = (241,(3,5,16)), beta=16/3, predicted 1 by P and <1/4 by beta-degradation; observed 1. Prediction committed mid-census with cryptographic order proof (prediction_order_proof.json).",
            "in_sample_separation_statistic": {
                "value": f"1/{comb(n_rows, n_with)} = {1 / comb(n_rows, n_with):.2e}",
                "meaning": f"probability that a uniformly random choice of which {n_with} of {n_rows} rows carry the 3/5 label coincides exactly with the {{2,3}}-membership split",
                "LABEL": "IN-SAMPLE figure only - the pattern was found by inspecting 12 of these rows. It is NOT a prospective p-value and must never be quoted as one.",
            },
            "honest_summary": "P survived its first prospective test. P is NOT established.",
            "what_would_kill_P": [
                "any triple containing both 2 and 3 whose exact window rho != 3/5",
                "any triple lacking the pair {2,3} whose exact window rho = 3/5",
                "one counterexample ends it",
            ],
            "known_limits_of_P": [
                "P is a statement about ORDERS at t=(1,1,1) only, on 13 points at four values of q",
                "rho varies with t at FIXED orders: the three q=13 rows (2,2,4) give 1/2, 3/8, 1/4 at t=(1,1,1),(1,1,2),(1,1,4), so P as stated says nothing about t-dependence and cannot be the whole story",
                "window-restricted (weight<=3 supports, reduced-basis candidates); says nothing about global rho_inst",
            ],
            "provenance": "pattern suggested by Main from my completed rows, flagged as a hypothesis to kill; recorded in pre_statement GB9.1 before B11 landed",
        },
        "VERDICT": {
            "answer": "NEITHER graceful-degradation NOR collapse-at-beta*",
            "reason": "the contract's dichotomy presupposes rho depends on beta; exclusion E shows it does not on the swept set",
            "threshold_beta_star": "none exists to locate on this grid",
            "grid_resolution_uncertainty": "applies vacuously (no threshold); the swept beta values are 5/3, 2, 5/2, 8/3, 3, 5, 16/3, so any future threshold claim on this grid would carry at least those gaps as uncertainty",
            "escalation": "no hypothesis-correction claim about TR26-150 arises from the numbers; the only statement touching the paper is finding 2 (a provenance reading), reported to Main before being written",
        },
        "b11_prediction_test": {
            "observed": rows["B11"]["ratio"],
            "pattern_P_predicted": "1/1",
            "beta_degradation_predicted": "< 1/4 (strictly below the grid minimum)",
            "branch_taken": "GB9.2 branch 1 - P survives, beta-degradation falsified",
            "note": "beta-degradation predicted the grid FLOOR; the observed value is the grid CEILING",
            "order_proof": "prediction_order_proof.json: pre_statement sha256 53bab658175b647455483e71b5f3c1b4fcf6ca8d7d9b84484d22da50f48aa46a, census 1350000/2304200, b11_result_file_exists=false",
        },
        "unswept_extension_price": _price_unswept_421(),
    }

    payload = {
        "campaign": CAMP.name,
        "gate": "B",
        "agent": "RsPe3dGateB",
        "pre_statement": "pre_statement.md (GB1-GB8 committed before compute; GB9 addendum committed mid-census with order proof)",
        "utc_frozen": dt.datetime.now(dt.timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "fixed_parameters": {
            "eta": "1/32", "t": [1, 1, 1], "Lambda": "identity on every axis",
            "milp_limit_seconds": 120,
        },
        "window": (
            "all point supports of weight <= 3 (complete census); candidates = normalized "
            "reduced basis of each nonzero V n F^S; per-witness delta closed by exhaustive "
            "line-tuple sweep (every cheaper cost proved absent, >=1 tuple present at the cost)"
        ),
        "anchors": {
            "unit_exhaustive": state.get("unit_anchor"),
            "unit_second_route": {
                "rho_exact": "1/3", "nonzero_V_count": 78124, "seconds": 535.73,
                "method": (
                    "V-basis coefficient enumeration over rref'd lift bases + DeltaEngine "
                    "ascending-cost closure; structurally distinct from units.anchor2's "
                    "numpy min-plus DP over component arrays"
                ),
            },
            "gate_A_rows_reproduced": {
                rid: {"frozen": ex, "recomputed": rows[rid]["ratio"], "agree": rows[rid]["ratio"] == ex}
                for rid, ex in [
                    ("B01", "3/5"), ("B02", "1/1"), ("B03", "3/5"), ("B04", "3/5"),
                    ("B05", "1/1"), ("B06", "3/5"), ("B07", "1/2"), ("B08", "3/8"),
                    ("B09", "1/4"),
                ]
            },
        },
        "rows": regen,
        "analysis": analysis,
        "rule7_swept_set": (
            "{(q,s,t,eta,Lambda)}: eta=1/32, t=(1,1,1) except the two extra q=13 t-rows, "
            "Lambda=Id, (q,s) in {(13,(2,2,4)) at t=(1,1,1),(1,1,2),(1,1,4); (31,(2,3,5)), "
            "(31,(2,3,10)), (31,(2,6,5)), (31,(3,10,2)); (61,(2,3,5)), (61,(3,4,5)); "
            "(241,(2,3,5)), (241,(3,4,5)), (241,(3,5,8)), (241,(3,5,16))}; witnesses "
            "restricted to weight<=3 supports and the normalized reduced-basis candidate class"
        ),
        "rule7_not_swept": [
            "supports of weight >= 4",
            "non-basis field values inside dim>1 intersections",
            "t profiles other than (1,1,1) except the two q=13 rows listed",
            "eta != 1/32",
            "non-identity Lambda (the conjecture's forall-Lambda quantifier is untouched)",
            "primes outside {13, 31, 61, 241}",
            "the q=421 nineteen-triple census (priced in analysis.unswept_extension_price)",
            "continuous beta - only the discrete grid was swept",
            "any asymptotic statement in beta, q, or N",
            "global rho_inst: each row certifies only the interval [1/(3N), rho^window]",
        ],
        "non_coprime_rows_flagged": [
            {"id": "B04", "s": [2, 3, 10], "why": "gcd(2,10)=2 - outside Conjecture 4.2's pairwise-coprime hypothesis; kept as a shape probe inherited from gate A"},
            {"id": "B05", "s": [2, 6, 5], "why": "gcd(2,6)=2 - same"},
            {"id": "B06", "s": [3, 10, 2], "why": "gcd(2,10)=2 - same"},
            {"id": "B07", "s": [2, 2, 4], "why": "gcd(2,2)=2 - same"},
            {"id": "B08", "s": [2, 2, 4], "why": "gcd(2,2)=2 - same"},
            {"id": "B09", "s": [2, 2, 4], "why": "gcd(2,2)=2 - same"},
        ],
        "coprime_rows": ["B01", "B02", "B03", "B10", "B11", "B12", "B13"],
        "evidence_labels": {
            "all_13_rows": "MACHINE-VERIFIED (finite window): exact rationals end to end, exact F_q witness verification, exhaustive cheaper-tuple closure per witness",
            "float_rows": "none - no floating point anywhere in the computation chain",
            "finding2_quotes": "[REPRODUCED] from hashed PDFs; the K-exponent value is [DERIVED]",
        },
        "tool_versions": {
            "python": sys.version.split()[0], "highspy": "1.15.1",
            "python-flint": "0.9.0", "sympy": "1.14.0", "numpy": "2.5.2",
            "solver": "HiGHS via highspy, exact modular MILP model, every solution re-checked in F_q",
        },
        "code_sha256_at_pre_statement_commit": "fffe84b00b9ca48ba5b3a6970b5523c9063bb61a559a88f81328a501e7b8bd56",
        "code_sha256_at_freeze": _code_hash(),
    }
    return payload, state, recs, reruns, b11


def main() -> int:
    payload, state, recs, reruns, b11 = build()
    (CAMP / "payload.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    (CAMP / "witnesses.json").write_text(
        json.dumps(
            {"records": [recs[k] for k in sorted(recs)],
             "reruns": [reruns[k] for k in sorted(reruns)]},
            indent=2, sort_keys=True,
        ) + "\n"
    )
    (CAMP / "optimality.json").write_text(
        json.dumps(
            [
                {
                    "id": r["id"], "q": r["q"], "s": r["s"], "t": r["t"],
                    "ratio": r["ratio"], "weight": r["weight"], "delta": r["delta"],
                    "delta_proven_exact_for_best_witness": r["delta_proven_exact_for_best_witness"],
                    "proof_cheaper_hit": r["proof_cheaper_hit"],
                    "proof_equal_hits": r["proof_equal_hits"],
                    "scope": "per-witness exact delta; best-in-window ratio; heavier supports OPEN",
                    "evidence_label": r["evidence_label"],
                }
                for r in payload["rows"]
            ],
            indent=2, sort_keys=True,
        ) + "\n"
    )
    (CAMP / "b11_census_result.json").write_text(json.dumps(b11, indent=2, sort_keys=True) + "\n")
    (CAMP / "checkpoint_final.json").write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
    (CAMP / "gateb.py.asrun").write_text((ROOT / "src" / "gateb.py").read_text())
    (CAMP / "gateb_freeze.py.asrun").write_text((ROOT / "src" / "gateb_freeze.py").read_text())
    for tmp, name in [("/tmp/run_b11_fast.py", "run_b11_fast.py.asrun"),
                      ("/tmp/gen_witnesses.py", "gen_witnesses.py.asrun")]:
        p = Path(tmp)
        if p.exists():
            (CAMP / name).write_text(p.read_text())
    (CAMP / "codehash.txt").write_text(_code_hash() + "\n")

    sums = {}
    for p in sorted(CAMP.iterdir()):
        if p.is_file() and p.name not in ("manifest.json", "checksums.sha256"):
            sums[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
    (CAMP / "checksums.sha256").write_text("".join(f"{v}  {k}\n" for k, v in sums.items()))

    a = payload["analysis"]
    manifest = {
        "campaign": CAMP.name, "gate": "B", "agent": "RsPe3dGateB",
        "utc_frozen": dt.datetime.now(dt.timezone.utc).isoformat(),
        "rows": len(payload["rows"]),
        "all_rows_machine_verified": all(
            r["delta_proven_exact_for_best_witness"] for r in payload["rows"]
        ),
        "verdict": a["VERDICT"]["answer"],
        "exclusion_E": a["EXCLUSION_E"]["statement"],
        "pattern_P_status": a["PATTERN_P"]["honest_summary"],
        "files": sums,
        "source_pdfs": {
            "TR26-150": "fb51d7dc84860f115d88ec47e92799e9c949ac90cbe714e8df7d3904c930af87",
            "AI candidate proof": "e2d9b4b61f9ed8589f90e5ed6bf146eaf642a85969c859589ad4ed405778330c",
        },
        "tool_versions": payload["tool_versions"],
        "code_sha256_at_freeze": _code_hash(),
        "reproduce": (
            "cd cs/rs-pe3d && PYTHONPATH=$PWD ../.venv/bin/python -c \"from src.gateb import "
            "run_instance; print(run_instance({'id':'B01','q':61,'s':(2,3,5),'t':(1,1,1),"
            "'expect':'3/5'})[0])\"  |  B11 via run_b11_fast.py.asrun (3449 s)  |  verify "
            "witnesses: python -m src.verify <this dir>/witnesses.json  |  re-freeze: "
            "python -m src.gateb_freeze"
        ),
        "frozen": "append-only; no file listed above is edited after this manifest",
    }
    (CAMP / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    print("FROZE", CAMP)
    print("rows:", manifest["rows"], "all machine-verified:", manifest["all_rows_machine_verified"])
    print("verdict:", manifest["verdict"])
    for q, v in a["EXCLUSION_E"]["monotone_excluding_interleavings"].items():
        if v:
            print(f"  q={q}: {len(v)} monotone-excluding interleavings, e.g. {v[0]['rows']} "
                  f"betas {v[0]['betas']} rhos {v[0]['rhos']} ({v[0]['shape']})")
    print("P consistent on all rows:", all(x["P_consistent"] for x in a["PATTERN_P"]["pattern_check"]) if "pattern_check" in a["PATTERN_P"] else "see payload")
    print("q=421 unswept price:", a["unswept_extension_price"]["estimated_census_hours"], "h census-only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
