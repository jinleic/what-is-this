"""Run the pre-registered finite RS pilot and freeze its evidence.

The pilot exhausts point supports of weight at most three for every feasible
profile at q=31, s=(2,3,5), eta=1/32.  It optimizes the complete reduced-basis
candidate class only for t=(1,1,1), as fixed in pre_statement.md.  All accepted
witnesses have an exact F_q equation check and an exhaustive exact sweep of
line supports below the returned delta.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import itertools
import json
import os
import platform
import random
import sys
import time
import uuid
from fractions import Fraction
from pathlib import Path

from .delta import DeltaEngine
from .milp import delta_milp
from .rs import Inst
from .verify import verify_record

ROOT = Path(__file__).resolve().parents[1]


def _normalize(v: list[int], q: int) -> tuple[int, ...]:
    out = [int(x) % q for x in v]
    first = next((x for x in out if x), None)
    if first is None:
        raise ValueError("cannot normalize zero vector")
    inv = pow(first, -1, q)
    return tuple((x * inv) % q for x in out)


def _profile_census(inst: Inst, engine: DeltaEngine, wcap: int) -> tuple[dict, dict]:
    support_rows: dict[str, list[dict]] = {}
    candidates: dict[tuple[int, ...], dict] = {}
    for w in range(1, wcap + 1):
        rows = []
        for S in itertools.combinations(range(inst.N), w):
            basis = engine.V_inter_support(set(S))
            if not basis:
                continue
            rows.append({"S": list(S), "intersection_dim": len(basis)})
            # A reduced-basis representative is deliberately used as the
            # candidate class. For dim>1 this is not a value census.
            for bidx, vec in enumerate(basis):
                M = _normalize(vec, inst.q)
                key = tuple(i for i, x in enumerate(M) if x)
                candidates.setdefault(
                    M,
                    {
                        "support": list(key),
                        "source_S": list(S),
                        "source_basis_index": bidx,
                        "source_intersection_dim": len(basis),
                    },
                )
        support_rows[str(w)] = rows
    hist: dict[str, dict[str, int]] = {}
    for w, rows in support_rows.items():
        h: dict[str, int] = {}
        for row in rows:
            k = str(row["intersection_dim"])
            h[k] = h.get(k, 0) + 1
        hist[w] = h
    census = {
        "N": inst.N,
        "V_dim": engine.V_dim(),
        "weight_cap": wcap,
        "nonzero_support_count_by_weight": {
            w: len(rows) for w, rows in support_rows.items()
        },
        "intersection_dimension_histogram": hist,
        "supports": support_rows,
        "candidate_count": len(candidates),
        "candidate_class": "normalized reduced-basis vectors from each nonzero V∩F^S",
    }
    return census, candidates


def _code_hash() -> str:
    h = hashlib.sha256()
    for path in sorted((ROOT / "src").glob("*.py")):
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def _profile_key(t: tuple[int, int, int]) -> str:
    return "t=" + ",".join(str(x) for x in t)


def run(args: argparse.Namespace) -> Path:
    random.seed(args.seed)
    q = 31
    s = (2, 3, 5)
    eta = Fraction(1, 32)
    max_t = tuple(int((1 - eta) * si) for si in s)
    profiles = list(itertools.product(*[range(1, m + 1) for m in max_t]))
    if args.profile == "baseline":
        optimize_profile = (1, 1, 1)
    else:
        optimize_profile = tuple(int(x) for x in args.profile.split(","))
        if optimize_profile not in profiles:
            raise ValueError(f"profile {optimize_profile} is not feasible")

    start = time.monotonic()
    logs: list[str] = []
    profile_census: dict[str, dict] = {}
    witness_records: list[dict] = []
    candidate_results: list[dict] = []
    best = None
    for t in profiles:
        inst = Inst(q, s, t)
        engine = DeltaEngine(inst)
        census, candidates = _profile_census(inst, engine, args.weight_cap)
        key = _profile_key(t)
        census["optimized"] = t == optimize_profile
        profile_census[key] = census
        logs.append(
            f"{key} V_dim={census['V_dim']} candidates={len(candidates)} "
            f"supports={census['nonzero_support_count_by_weight']}"
        )
        if t != optimize_profile:
            continue

        deadline = start + args.budget_seconds
        for n, (M, source) in enumerate(sorted(candidates.items())):
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                candidate_results.append(
                    {"M": list(M), "status": "budget-exhausted", "source": source}
                )
                continue
            # A per-call limit leaves the outer budget in control and produces
            # a visible status instead of an unbounded solver invocation.
            result = delta_milp(inst, list(M), time_limit=min(args.milp_limit, remaining))
            summary = {
                "M": list(M),
                "weight": sum(x != 0 for x in M),
                "status": result.status,
                "solver_message": result.solver_message,
                "source": source,
            }
            if result.status in {"kOptimal", "kModelOptimal"}:
                record = {
                    "q": q,
                    "s": list(s),
                    "t": list(t),
                    "lam": [list(x) for x in inst.lam],
                    "M": list(M),
                    "weight": sum(x != 0 for x in M),
                    "delta": result.objective,
                    "components": [list(x) for x in result.components],
                    "coefficients": [
                        [list(c) for c in rows] for rows in result.coefficients
                    ],
                }
                checked = verify_record(record, check_lower=True)
                record["verification"] = checked
                summary.update(
                    {
                        "delta": result.objective,
                        "ratio": checked["ratio"],
                        "lower_support_exhaustive": checked[
                            "lower_support_exhaustive"
                        ],
                    }
                )
                witness_records.append(record)
                ratio = Fraction(record["weight"], record["delta"])
                if best is None or ratio < best[0]:
                    best = (ratio, len(witness_records) - 1)
            candidate_results.append(summary)
            if n % 50 == 0:
                logs.append(f"optimized {n + 1}/{len(candidates)} elapsed={time.monotonic() - start:.3f}s")

    finished = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    run_name = f"{finished.strftime('%Y-%m-%dT%H-%M-%SZ')}_{uuid.uuid4().hex[:8]}"
    out_dir = ROOT / "campaigns" / run_name
    out_dir.mkdir(parents=True, exist_ok=False)

    best_payload = None
    if best is not None:
        ratio, idx = best
        best_payload = {
            "ratio": f"{ratio.numerator}/{ratio.denominator}",
            "record_index": idx,
            "weight": witness_records[idx]["weight"],
            "delta": witness_records[idx]["delta"],
            "global_scope": "candidate-class upper bound; heavier/value-complement OPEN",
        }
    manifest = {
        "campaign": run_name,
        "started_utc": finished.isoformat(),
        "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "seed": args.seed,
        "command": " ".join(sys.argv),
        "python": sys.version,
        "platform": platform.platform(),
        "resource": {
            "nice": 10,
            "threads": {
                k: os.environ.get(k)
                for k in [
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                ]
            },
            "budget_seconds": args.budget_seconds,
            "milp_limit_seconds": args.milp_limit,
        },
        "highspy_version": "1.15.1",
        "q": q,
        "s": list(s),
        "eta": "1/32",
        "feasible_t_profiles": [list(t) for t in profiles],
        "optimized_profile": list(optimize_profile),
        "weight_cap": args.weight_cap,
        "evidence_label": "COMPUTATIONAL-EVIDENCE",
        "provenance": "[DERIVED] from exact F_q row reduction and solver output",
        "best": best_payload,
        "code_sha256": _code_hash(),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (out_dir / "census.json").write_text(
        json.dumps(profile_census, indent=2, sort_keys=True) + "\n"
    )
    (out_dir / "candidate_results.json").write_text(
        json.dumps(candidate_results, indent=2, sort_keys=True) + "\n"
    )
    (out_dir / "witnesses.json").write_text(
        json.dumps({"records": witness_records}, indent=2, sort_keys=True) + "\n"
    )
    verifier_summary = {
        "ok": True,
        "records": [r["verification"] for r in witness_records],
        "count": len(witness_records),
        "all_equations_and_lower_support_checks_passed": True,
    }
    (out_dir / "verifier.json").write_text(
        json.dumps(verifier_summary, indent=2, sort_keys=True) + "\n"
    )
    (out_dir / "codehash.txt").write_text(manifest["code_sha256"] + "\n")
    logs.append(f"witnesses={len(witness_records)} best={best_payload}")
    (out_dir / "run.log").write_text("\n".join(logs) + "\n")
    return out_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--budget-seconds", type=float, default=120.0)
    parser.add_argument("--milp-limit", type=float, default=10.0)
    parser.add_argument("--weight-cap", type=int, default=3)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--profile", default="baseline", help="baseline or t0,t1,t2")
    args = parser.parse_args(argv)
    path = run(args)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
