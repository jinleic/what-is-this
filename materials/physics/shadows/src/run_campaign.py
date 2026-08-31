"""Campaign runner for physics/shadows (owner: shadows).

Runs the frozen gates from ../pre_statement.md (2026-08-29) with EXACT
arithmetic, and writes an immutable campaign snapshot under ../campaigns/.

Usage:  nice -n 10 physics/.venv/bin/python physics/shadows/src/run_campaign.py

No RNG anywhere. Floats only in fields explicitly keyed "numerical_*".
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import platform
import sys
import time
import uuid
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import shadows_exact as se  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def frac_str(x: Fraction | int) -> str:
    x = Fraction(x)
    return f"{x.numerator}/{x.denominator}"


def gate_b_exact_extra() -> dict:
    """Exact Eq.(5)-summed weight == Eq.(4) for ALL k <= 16 (integer counts,
    Fractions). Closes the k in 11..16 gap that the exhaustive-token path
    (k <= 10) leaves open."""
    mism, rows = [], []
    for k in range(1, 17):
        counts = se.eq5_hist_counts(k)
        wsum = sum(Fraction(c, 3 ** k) * Fraction(1, 3 ** m)
                   for m, c in counts.items() if c)
        w4 = se.eq4_weight(k)
        ok = (wsum == w4)
        if not ok:
            mism.append(k)
        rows.append({"k": k, "eq5_sum_weight": frac_str(wsum),
                     "eq4_weight": frac_str(w4), "match": ok})
    return {"eq5_vs_eq4_exact_k16": rows, "mismatches": mism}


def main() -> int:
    t_start = time.time()
    ts = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    uid = uuid.uuid4().hex[:8]
    camp = None  # set after content known? Contract wants dir first; we
    # create dir now and finalize the inner hash at the end.
    camp = ROOT / "campaigns" / f"{ts}_{uid}"
    camp.mkdir(parents=True, exist_ok=False)
    scratch_dir = camp  # artifacts live flat: gates_a.json, ...
    print(f"[campaign] {camp}")

    report: dict = {"frozen_pre_statement": "../pre_statement.md (2026-08-29)",
                    "started_utc": ts, "campaign_dir": str(camp)}

    # ---------------- selftest -------------------------------------------
    t0 = time.time()
    se.selftest()
    report["selftest_ok"] = True
    report["selftest_elapsed_s"] = round(time.time() - t0, 3)
    report["python"] = sys.version.split()[0]
    report["platform"] = platform.platform()

    # ---------------- gate A ---------------------------------------------
    t0 = time.time()
    a = se.gate_a_census()
    a_elapsed = round(time.time() - t0, 3)
    a_pass = (a["sp42_order"] == 720 and a["mod_phase_labels"] == 11520
              and a["mod_phase_labels_crosscheck"] == 11520
              and a["max_contracted"] == 4 and a["anchor_contracted"] == 4)
    gate_a = {
        "gate": "A", "quantity": "exhaustive two-qubit contraction census",
        "result": a, "elapsed_s": a_elapsed,
        "pass_fail": "PASS" if a_pass else "FAIL",
        "verdict": "CONFIRMED" if a_pass else "REFUTED",
        "evidence_label": "PROVED (exhaustive, exact GF(2) integer arithmetic)",
        "paper_claim": "Methods lemma: at most 4 of the 9 size-2 two-qubit "
                       "Paulis can be contracted to size 1 by any Clifford "
                       "(npj QI 12, 86 (2026), DOI 10.1038/s41534-026-"
                       "01227-w = arXiv:2412.01850).",
    }
    print(f"[gate A] {gate_a['verdict']} max={a['max_contracted']} "
          f"hist={a['histogram_by_contracted']} ({a_elapsed}s)")

    # ---------------- gate B ---------------------------------------------
    t0 = time.time()
    b = se.gate_b_run()
    b_elapsed = round(time.time() - t0, 3)
    extra = gate_b_exact_extra()
    b["eq5_vs_eq4_exact_k16"] = extra
    checks = dict(b["checks"])
    checks["eq5_sum_eq4_exact_k16"] = not extra["mismatches"]
    ratios = [row["numerical_ratio_norm2_over_2x1.8^k"]
              for row in b["per_k"]]
    trend_ok = ratios[-1] > ratios[len(ratios) // 2] > ratios[0] and \
        ratios[-1] > 0.999
    b_pass = all(bool(v) for v in checks.values()) and trend_ok
    gate_b = {
        "gate": "B", "quantity": "exact re-derivation of Eq. (4)/(5)/(7) "
        "weights and random-Clifford baseline",
        "result": b, "elapsed_s": b_elapsed,
        "checks": checks,
        "asymptote_trend_note": "numerical ratios ||O||^2/(2*1.8^k) rise "
                                "monotonically toward 1 (k=1..16); NUMERICAL "
                                "tag, not evidence",
        "pass_fail": "PASS" if b_pass else "FAIL",
        "evidence_label": "PROVED for k<=10 (exhaustive tokens) and "
                          "integer-exact for k<=16 (Eq.(5)-sum vs Eq.(4)); "
                          "ratios NUMERICAL",
    }
    print(f"[gate B] {'PASS' if b_pass else 'FAIL'} checks={checks} "
          f"({b_elapsed}s)")

    # ---------------- gate C ---------------------------------------------
    t0 = time.time()
    c = se.gate_c_run(kmax=8)
    c_elapsed = round(time.time() - t0, 3)
    n_beats = len(c["beats_uct"])
    gate_c = {
        "gate": "C", "quantity": "bounded alternative-family search "
        "(pre_statement pinned families a/b/c, k <= 8)",
        "n_composites_family_b": c["n_composites_family_b"],
        "notes": c["notes"],
        "best_family": (
            {"family": c["best_family"]["family"],
             "row": c["best_family"]["row"],
             "weight_ratio_vs_uct":
                 c["best_family"]["weight_ratio_vs_uct"]}
            if c["best_family"] else None),
        "max_w_over_uct_exact": c["max_w_over_uct_exact"],
        "families": c["families"],
        "beats_uct": c["beats_uct"],
        "elapsed_s": c_elapsed,
        "pass_fail": "PASS" if n_beats else "NOT-FOUND",
        "verdict": ("IMPROVEMENT-FOUND (smaller shadow norm than U_ct "
                    "within pinned scope; paper's optimality framing is "
                    "family-limited)" if n_beats
                    else "NOT-FOUND-in-pinned-scope (NOT an obstruction "
                         "proof; non-goal: full family space)"),
        "evidence_label": "PROVED within pinned scope (exact per composite)",
    }
    print(f"[gate C] {gate_c['verdict']} beats={n_beats} ({c_elapsed}s)")

    report["gates"] = {"A": gate_a, "B": gate_b, "C": gate_c}
    elapsed = round(time.time() - t_start, 3)
    report["total_elapsed_s"] = elapsed

    # ---------------- serialize + hash ------------------------------------
    def dump(name: str, obj: dict) -> str:
        p = camp / name
        p.write_text(json.dumps(obj, indent=1, sort_keys=False) + "\n")
        return p.read_bytes()

    bytes_a = dump("gates_a.json", gate_a)
    bytes_b = dump("gates_b.json", gate_b)
    bytes_c = dump("gates_c.json", gate_c)
    summary = {
        "campaign": camp.name,
        "frozen_statement": "pre_statement.md (2026-08-29), run may only "
                            "confirm/refute, never redefine",
        "paper": "npj QI 12, 86 (2026), DOI 10.1038/s41534-026-01227-w "
                 "(= arXiv:2412.01850)",
        "gate_A": {"verdict": gate_a["verdict"],
                   "pass_fail": gate_a["pass_fail"],
                   "label": gate_a["evidence_label"],
                   "max_contracted": a["max_contracted"],
                   "histogram": a["histogram_by_contracted"],
                   "n_achieving_actions": a["n_achieving_actions"],
                   "note_derived": "all 9 achieving class-4 sets equal "
                                   "(row∪column)\\{cell} of paper Table 2 "
                                   "(hand-checked against the published "
                                   "table layout) [DERIVED]",
                   "elapsed_s": a_elapsed},
        "gate_B": {"pass_fail": gate_b["pass_fail"],
                   "checks": checks,
                   "label": gate_b["evidence_label"],
                   "eq7_mismatches": sum(
                       1 for r in b["identity_insertion"] if not r["match"]),
                   "ratio_trend_last": ratios[-1],
                   "elapsed_s": b_elapsed},
        "gate_C": {"pass_fail": gate_c["pass_fail"],
                   "verdict": gate_c["verdict"],
                   "label": gate_c["evidence_label"],
                   "elapsed_s": c_elapsed},
        "exactness_argument": "all census/weight quantities computed on "
                              "Python ints as GF(2) bit labels and "
                              "fractions.Fraction; no RNG anywhere; floats "
                              "only in keys tagged numerical_*; gate A "
                              "enumerates 2^16 maps -> asserts |Sp(4,2)|=720 "
                              "-> 720 x 16 = 11520 mod-phase labels; "
                              "exhaustive therefore universal, not sampled.",
        "single_core_low_priority": "nice -n 10, one process",
        "total_elapsed_s": elapsed,
        "python": report["python"],
        "platform": report["platform"],
    }
    h = hashlib.sha256(bytes_a + bytes_b + bytes_c).hexdigest()[:12]
    final_name = f"{ts}_{uid}_{h}"
    summary["campaign"] = final_name
    (camp / "SUMMARY.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    (camp / "MANIFEST.txt").write_text(
        "campaign " + final_name + "\n"
        "artifact_sha256_12(gates_a+b+c): " + h + "\n"
        "frozen pre-statement: pre_statement.md (2026-08-29)\n"
        "paper: npj QI 12, 86 (2026), DOI 10.1038/s41534-026-01227-w "
        "(= arXiv:2412.01850); follow-up: arXiv:2608.18935\n"
        "exactness: GF(2) bit-label integers + fractions.Fraction only; "
        "no RNG; floats only under keys containing 'numerical'\n"
        "evidence labels: gate A PROVED; gate B PROVED(k<=10)/exact-integer"
        "(k<=16)+NUMERICAL ratio; gate C PROVED-in-scope\n"
        f"total_elapsed_s: {elapsed}\n")
    # FINAL dir name per contract <UTC>_<uuid>_<hash>:
    final = camp.parent / final_name
    camp.rename(final)
    print(f"[campaign] finalized {final.name}")
    print(json.dumps({k: summary[k] for k in ("gate_A", "gate_B", "gate_C")},
                     indent=1)[:1600])
    return 0 if (gate_a["pass_fail"] == "PASS" and gate_b["pass_fail"]
                 == "PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
