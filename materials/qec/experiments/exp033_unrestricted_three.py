"""EXP-033: unrestricted-class depth decision for the three TI-refuted PBB members.

EXP-031 (rerun, verdict NEGATIVE) certified that pbb144-01, pbb144-03 and
pbb144-11 (all w_max = 9) admit NO translation-invariant schedule at any depth
9..13: the depth criterion's prediction (10) and the two-value law are refuted
IN THE TRANSLATION-INVARIANT CLASS for these members.

Decisive question: does the two-value law survive in the UNRESTRICTED class?
  - feasible at 9 or 10 unrestricted  -> law survives unrestricted; the TI
    failure is a structural cost of translation invariance;
  - certified INFEASIBLE at 9 and 10  -> universal counterexample to the law;
  - UNKNOWN                            -> honestly undecided at this budget.

v2: the first version stalled inside exp031's PySAT criterion stage, which has
no time limit and is unnecessary here - only the CP-SAT ground truth matters.
This version calls _certify_minimum_depth directly (every solve time-limited).
Canonical write only after all three instances complete (any exception aborts).
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(var, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

EXPERIMENT = "EXP-033"
MASTER_SEED = 20260812
TARGET_INDICES = (1, 3, 11)          # pbb144-01 / -03 / -11 in exp031 order
TIME_LIMIT_S = 3600.0                # per depth (CP-SAT, certified statuses only)
MAX_EXTRA_DEPTH = 1                  # depths w_max, w_max+1: exactly the law's range
WORKERS = 8
OUT = ROOT / "results" / "processed" / "exp033_unrestricted_three.json"

_SPEC = importlib.util.spec_from_file_location(
    "exp031_depth_criterion", ROOT / "experiments" / "exp031_depth_criterion.py")
_E31 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _E31
_SPEC.loader.exec_module(_E31)

_E31.CP_SAT_TIME_LIMIT_S = TIME_LIMIT_S
_E31.CP_SAT_MAX_EXTRA_DEPTH = MAX_EXTRA_DEPTH
_E31.CP_SAT_WORKERS = WORKERS


def main() -> None:
    t0 = time.time()
    from qec_research.circuits.bicycle_schedule import pbb_supports_and_orbits
    from qec_research.circuits.scheduling import generator_supports

    items = _E31._pbb_instances()
    results = []
    for idx in TARGET_INDICES:
        item = items[idx]
        code, _, _orbits = pbb_supports_and_orbits(item["spec"])
        n = code.H.shape[1] // 2
        supports = generator_supports(code.H)
        w_max = max(s.weight for s in supports)
        seed = (MASTER_SEED + 7919 * idx) % 2_147_483_647
        cert = _E31._certify_minimum_depth(
            supports, n, w_max, seed,
            symmetry_orbits=None,        # UNRESTRICTED class
            criterion_witness=None,
        )
        trace = [(e["depth"], e["status"]) for e in cert["trace"]]
        entry = {
            "exp031_instance": f"pbb144-{idx:02d}",
            "label": item["label"],
            "catalogue_row_number": item["row_number"],
            "w_max": w_max,
            "schedule_class": "unrestricted",
            "decided": cert["decided"],
            "certified_depth": cert["certified_depth"],
            "trace": cert["trace"],
            # Persist the independently verified witness, not only the solver
            # status, so the exact-depth claim is machine-checkable offline.
            "slot": cert.get("slot"),
            "verification": cert.get("verification"),
            "two_value_survives": (cert["decided"]
                                   and cert["certified_depth"] in {w_max, w_max + 1}),
            "two_value_refuted": (_E31.certified_infeasible_at(cert["trace"], w_max)
                                  and _E31.certified_infeasible_at(cert["trace"], w_max + 1)),
        }
        print(f"pbb144-{idx:02d} ({item['label']}): w_max={w_max} "
              f"cert={cert['certified_depth']} trace={trace}", flush=True)
        results.append(entry)

    n_survive = sum(1 for r in results if r["two_value_survives"])
    n_refuted = sum(1 for r in results if r["two_value_refuted"])
    n_unknown = len(results) - n_survive - n_refuted
    if n_refuted:
        verdict, reason = "BREAKTHROUGH_CANDIDATE", (
            f"{n_refuted}/3 members certifiably violate the two-value law in the "
            "UNRESTRICTED class - universal counterexample")
    elif n_survive == 3:
        verdict, reason = "POSITIVE", (
            "two-value law survives unrestricted for all three members; the TI-class "
            "failure is a structural cost of translation invariance")
    elif n_survive:
        verdict, reason = "POSITIVE", (
            f"law survives unrestricted for {n_survive}/3; {n_unknown} undecided at budget")
    else:
        verdict, reason = "INCONCLUSIVE", (
            f"all three undecided at {TIME_LIMIT_S:.0f}s/depth unrestricted budget")

    payload = {
        "experiment": EXPERIMENT,
        "protocol_constants": {
            "master_seed": MASTER_SEED,
            "target_indices": list(TARGET_INDICES),
            "time_limit_s_per_depth": TIME_LIMIT_S,
            "max_extra_depth": MAX_EXTRA_DEPTH,
            "workers": WORKERS,
            "schedule_class": "unrestricted",
            "provenance": "exp031 rerun mismatches (TI-class certified refutations)",
            "note": ("v3 persists slot + verification witnesses and gates the "
                     "canonical write on all three rows carrying them"),
        },
        "results": results,
        "n_survive_two_value": n_survive,
        "n_refuted_unrestricted": n_refuted,
        "n_undecided": n_unknown,
        "verdict": verdict,
        "verdict_reason": reason,
        "wall_s": round(time.time() - t0, 1),
    }
    # FR-012 discipline: the canonical exact-status artifact may only be
    # replaced by a strictly stronger one - every row decided, with a persisted
    # schedule witness that passed independent verification.  Anything weaker
    # (UNKNOWN rows, missing witnesses) is routed to results/partial_runs/.
    fully_witnessed = all(
        r["decided"]
        and r.get("slot") is not None
        and bool((r.get("verification") or {}).get("valid"))
        for r in results
    )
    if fully_witnessed:
        destination = OUT
    else:
        stamp = int(time.time())
        destination = (ROOT / "results" / "partial_runs"
                       / f"exp033_unrestricted_three-UNGATED-{stamp}.json")
        payload["not_canonical_reason"] = (
            "at least one row lacks a decided status or a valid persisted witness"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=str))
    tmp.replace(destination)
    print(f"wrote {destination.relative_to(ROOT)}: {verdict} - {reason}", flush=True)


if __name__ == "__main__":
    main()
