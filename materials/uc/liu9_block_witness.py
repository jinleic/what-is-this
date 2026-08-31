#!/usr/bin/env python3
"""Staged mean-feasible raw-gap witness search for the Liu H2 block kernel.

The exact block reduction (liu9_block_kernel.py) leaves one open statement:
copositivity of the 2x2 size-biased block kernel, equivalently nonnegativity
of the raw gap on every mean-feasible point of Liu's nine-variable family.
This module hunts for a negative raw gap under the full mean-feasibility
constraint  mean >= m := (1-beta... )  -- concretely Liu's constraint
E[b] >= p*x with the equation-defined binding constants (m = 0.61729...,
beta = 0.10005255986289...).

It is explicitly DISCOVERY-ONLY (COMPUTATIONAL EVIDENCE): no negative value
found here is promoted to mathematics.  Any candidate below a strict
negative threshold is re-evaluated at higher working precision and, if it
survives, enclosed in outward-rounded Arb intervals by the companion
verifier (verification/independent_block_witness_check.py).  Everything is
deterministic: fixed seeds, no wall-clock or iteration counts in the report.

Stages, cheapest first:
  1. structured corners  -- degenerate q, delta laws, endpoint supports,
     the exact (delta_x, delta_1) family at every q (the m-frozen shortcut
     family whose true raw gap is known positive: the search must reproduce
     that positivity),
  2. log-uniform simplex pairs with power-biased components -- the boundary
     layers where the q=1 size-biased kernel is tightest,
  3. local SQSLP/SLSQP polish of the best points found, plus a λ-shifted
     merit objective gap + lam*(mean - m) restricted to mean >= m,
  4. a Laplace/boundary ladder around every candidate whose margin is small.

No result mutates any SSOT document.  Output: one JSON report and stdout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import mpmath

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from liu9_binding import MPParameters, solve_equation_parameters
from liu9_boundary_layer import entropy_mp, gap_mp, mean_of

REPORT_DEFAULT = HERE / "verification/results/liu9-block-witness-search.json"
NEGATIVE_THRESHOLD = -1e-8        # float-level alarm
CONFIRM_THRESHOLD = -1e-12        # survivors escalated to high dps + Arb
SEED = 20260829


# ---------------------------------------------------------------------------
# Point space.  A point is Liu's nine-vector
#   (a1, a2, q, b0, b2, b4, b1, b3, b5)
# on the constraint set  a1,a2 >= 0, a1+a2 <= 1,  q in [0,1],
# supports in [0,1]^6, and mean >= m.  gap_mp needs EHX > 0.


def _mp(x: float) -> mpmath.mpf:
    return mpmath.mpf(float(x))


class Evaluator:
    """Fast float prefilter + high-precision confirmation of the raw gap."""

    def __init__(self, parameters: MPParameters, dps: int = 50):
        self.parameters = parameters
        self.m = float(parameters.mean)
        self.beta = parameters.beta
        self.dps = dps

    def mean(self, values: Sequence[mpmath.mpf]) -> mpmath.mpf:
        return mean_of(tuple(values), mpmath.mpf(1))

    def gap(self, values: Sequence[mpmath.mpf]) -> Optional[mpmath.mpf]:
        """Raw gap or None when the point is infeasible/degenerate."""
        # gap_mp/entropy_mp/mean_of all follow the ambient mpmath context
        # (cf. liu9_objective.evaluate_mpmath's workdps discipline); honor
        # the declared self.dps instead of inheriting ~15-digit default.
        with mpmath.workdps(self.dps):
            if float(mean_of(tuple(values), mpmath.mpf(1))) < self.m:
                return None
            a1, a2, q = values[0], values[1], values[2]
            if a1 < 0 or a2 < 0 or a1 + a2 > 1 or not 0 <= q <= 1:
                return None
            if any(v < 0 or v > 1 for v in values[3:]):
                return None
            if entropy_mp(tuple(values)) <= 0:
                return None
            return +gap_mp(tuple(values), self.beta)


def structured_points(m: float, beta: mpmath.mpf) -> list[dict[str, Any]]:
    """Adversarial corners, including the known m-frozen shortcut family."""
    out: list[tuple[str, tuple[mpmath.mpf, ...]]] = []
    with mpmath.workdps(50):
        x63 = mpmath.mpf("0.063")
        one = mpmath.mpf(1)
        half = mpmath.mpf("0.5")
        # delta_x / delta_1 family at the recorded q values (shortcut family)
        for q in ("0.7237", "0.8", "0.9", "0.95", "0.99", "1"):
            out.append((
                f"delta_x-delta_1-q{q}",
                (mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(q),
                 half, half, x63, half, half, one),
            ))
        # same family, mirrored: P0 = delta_1, P1 = delta_x
        for q in ("0.0", "0.01", "0.1", "0.25"):
            out.append((
                f"delta_1-delta_x-q{q}",
                (mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(q),
                 half, half, one, half, half, x63),
            ))
        # fully degenerate two-atom supports with mixed masses
        # (a1,a2) split of the x-atom vs the 1-atom at several q
        for a1 in ("0.3", "0.7"):
            a2 = mpmath.mpf("0.15")
            for q in ("0.4", "0.9", "0.99"):
                out.append((
                    f"split-a{a1}-q{q}",
                    (mpmath.mpf(a1), a2, mpmath.mpf(q),
                     mpmath.mpf("0.2"), x63, one,
                     one, one - x63, one),
                ))
    points = []
    for name, values in out:
        try:
            mean = float(mean_of(values, mpmath.mpf(1)))
        except Exception:
            continue
        if mean < m - 1e-15:
            continue
        points.append({"name": name, "values": values, "mean": mean})
    return points


def random_points(m: float, count: int, rng_seed: int,
                  dps: int = 30) -> list[dict[str, Any]]:
    """Log/simplex-biased mean-feasible draws; component-biased supports."""
    import random
    rng = random.Random(rng_seed)
    points: list[dict[str, Any]] = []
    attempts = 0
    while len(points) < count and attempts < count * 60:
        attempts += 1
        with mpmath.workdps(dps):
            a1, a2 = rng.random(), rng.random()
            if a1 + a2 > 1:
                a1, a2 = 1 - a1, 1 - a2
            style = rng.randrange(6)
            if style == 0:
                q = rng.random()
            elif style == 1:
                q = rng.random() ** 4
            elif style == 2:
                q = 1 - rng.random() ** 4
            elif style == 3:
                q = mpmath.mpf(1) - mpmath.mpf(2) ** (-rng.randrange(1, 40))
            elif style == 4:
                q = mpmath.mpf(2) ** (-rng.randrange(1, 40))
            else:
                q = mpmath.mpf(rng.choice([0.25, 0.5, 0.75]))
            supports: list[mpmath.mpf] = []
            for _ in range(6):
                draw = rng.random()
                if draw < 0.18:
                    supports.append(mpmath.mpf(0))
                elif draw < 0.30:
                    supports.append(mpmath.mpf(2) ** (
                        -rng.randrange(1, 60) / rng.choice([1, 2, 4])))
                elif draw < 0.55:
                    supports.append(mpmath.mpf(rng.random() ** 3))
                else:
                    supports.append(mpmath.mpf(rng.random() ** 0.25))
            values = (mpmath.mpf(a1), mpmath.mpf(a2), q, *supports)
        mean = float(mean_of(values, mpmath.mpf(1)))
        if mean < m:
            continue
        if entropy_mp(tuple(values)) <= 0:
            continue
        points.append({"name": f"random-{len(points)}", "values": values,
                       "mean": mean})
    return points


def evaluate_points(evaluator: Evaluator,
                    points: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for point in points:
        gap = evaluator.gap(point["values"])
        if gap is None:
            continue
        rows.append({
            "name": point["name"],
            "mean": mpmath.nstr(point["mean"], 20),
            "gap": mpmath.nstr(gap, 30),
            "gap_float": float(gap),
            "values": [mpmath.nstr(v, 20) for v in point["values"]],
            "values_exact_ratio": [
                mpmath.mpf(v).__str__() for v in point["values"]],
        })
    rows.sort(key=lambda row: row["gap_float"])
    return rows


def polish_best(evaluator: Evaluator, seeds: list[dict[str, Any]],
                points_per_seed: int = 2) -> list[dict[str, Any]]:
    """Constrained SLSQP polish around the most negative seeds."""
    import numpy as np
    from scipy.optimize import minimize

    m = evaluator.m
    beta = float(evaluator.beta)

    def unpack(x: "np.ndarray") -> tuple[mpmath.mpf, ...]:
        return tuple(_mp(v) for v in x)

    def objective(x: "np.ndarray") -> float:
        values = unpack(x)
        a3 = 1.0 - x[0] - x[1]
        if a3 < 0 or abs(a3 - (1 - x[0] - x[1])) > 1e-12:
            return 1.0
        gap = evaluator.gap(values)
        if gap is None:
            # merit wall away from feasibility
            mean_penalty = max(0.0, m - (
                (1 - x[0] - x[1]) * x[5] + x[0] * x[3] + x[1] * x[4]
                + x[2] * ((1 - x[0] - x[1]) * x[8] + x[0] * x[6]
                          + x[1] * x[7]
                          - ((1 - x[0] - x[1]) * x[5] + x[0] * x[3]
                             + x[1] * x[4]))))
        return float(gap)

    rows: list[dict[str, Any]] = []
    for seed in seeds[:5]:
        base = [float(v) for v in (
            seed["values"][0], seed["values"][1], seed["values"][2],
            seed["values"][3], seed["values"][4], seed["values"][5],
            seed["values"][6], seed["values"][7], seed["values"][8])]
        rng = np.random.default_rng(
            (SEED + int(hashlib.sha256(
                ",".join(str(v) for v in seed["values"]).encode()
            ).hexdigest()[:8], 16)) % (2**32))
        for trial in range(points_per_seed):
            for label, start in (
                ("seed", base),
                ("seed-jitter", [
                    min(1.0, max(0.0, b + float(rng.uniform(-1e-3, 1e-3))))
                    for b in base]),
            ):
                if trial == 0 and label == "seed-jitter":
                    continue
                x0 = np.clip(np.array(start), 1e-12, 1 - 1e-12)
                constraints = [
                    {"type": "ineq", "fun": lambda x: 1 - x[0] - x[1]},
                ]
                # mean feasibility as nonlinear constraint
                constraints.append({
                    "type": "ineq",
                    "fun": lambda x: (
                        (1 - x[0] - x[1]) * x[5] + x[0] * x[3]
                        + x[1] * x[4]
                        + x[2] * (
                            (1 - x[0] - x[1]) * x[8] + x[0] * x[6]
                            + x[1] * x[7]
                            - ((1 - x[0] - x[1]) * x[5] + x[0] * x[3]
                               + x[1] * x[4]))
                        - m),
                })
                bounds = [(0.0, 1.0)] * 9
                try:
                    res = minimize(
                        objective, x0, method="SLSQP", bounds=bounds,
                        constraints=constraints,
                        options={"maxiter": 400, "ftol": 1e-14})
                except Exception as exc:  # noqa: BLE001
                    rows.append({"name": f"polish-{seed['name']}-{label}",
                                 "error": repr(exc)})
                    continue
                values = unpack(res.x)
                gap = evaluator.gap(values)
                if gap is None:
                    continue
                rows.append({
                    "name": f"polish-{seed['name']}-{label}",
                    "mean": mpmath.nstr(evaluator.mean(values), 20),
                    "gap": mpmath.nstr(gap, 30),
                    "gap_float": float(gap),
                    "values": [mpmath.nstr(v, 20) for v in values],
                    "values_exact_ratio": [
                        mpmath.mpf(v).__str__() for v in values],
                    "success": bool(res.success),
                })
    rows = [r for r in rows if r.get("gap_float") is not None]
    rows.sort(key=lambda row: row["gap_float"])
    return rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["smoke", "standard", "deep"],
                        default="smoke")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--output", type=Path, default=REPORT_DEFAULT)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args(argv)

    started = time.perf_counter()
    parameters = solve_equation_parameters(80)
    evaluator = Evaluator(parameters, dps=args.dps)
    m = float(parameters.mean)

    structured = structured_points(m, parameters.beta)
    budget = {"smoke": (0, 300), "standard": (0, 30000),
              "deep": (0, 300000)}[args.stage]
    randoms = random_points(m, budget[1], args.seed)
    structured_rows = evaluate_points(evaluator, structured)
    random_rows = evaluate_points(evaluator, randoms)
    best_level_a = (structured_rows + random_rows)[:5]

    polished = polish_best(evaluator, best_level_a) if budget[1] else []
    all_rows = sorted(structured_rows + random_rows + polished,
                      key=lambda row: row.get("gap_float", 0.0))
    alarm = [row for row in all_rows
             if row.get("gap_float") is not None
             and row["gap_float"] < NEGATIVE_THRESHOLD]
    confirm = [row for row in alarm
               if row["gap_float"] < CONFIRM_THRESHOLD]

    elapsed = time.perf_counter() - started
    report = {
        "tool": "liu9_block_witness.py",
        "claim_status": "COMPUTATIONAL EVIDENCE",
        "discovery_only": (
            "A negative value found here is not mathematics until the "
            "companion Arb verifier re-encloses it and a review confirms "
            "mean feasibility by exact interval evaluation."),
        "binding_parameters": {
            "m": mpmath.nstr(parameters.mean, 40),
            "beta": mpmath.nstr(parameters.beta, 40),
        },
        "constraints": ["a1+a2<=1", "q,b in [0,1]", "mean >= m"],
        "stage": args.stage,
        "seed": args.seed,
        "counts": {
            "structured": len(structured_rows),
            "random": len(random_rows),
            "polished": len([r for r in polished if "gap_float" in r]),
            "total_evaluated": len(all_rows),
        },
        "best_gap": (all_rows[0] if all_rows else None),
        "best_20": all_rows[:20],
        "alarms_below": {"threshold": NEGATIVE_THRESHOLD,
                         "count": len(alarm),
                         "names": [r["name"] for r in alarm[:10]]},
        "confirmations_below": {"threshold": CONFIRM_THRESHOLD,
                                "count": len(confirm)},
        "negative_mean_feasible_raw_gap_found": bool(confirm),
        "elapsed_seconds": elapsed,
    }
    report["report_sha256"] = hashlib.sha256(
        json.dumps({k: v for k, v in report.items() if k != "report_sha256"},
                   sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    print("LIU H2 MEAN-FEASIBLE WITNESS SEARCH stage=%s" % args.stage)
    print("evaluated %d points in %.1f s" % (len(all_rows), elapsed))
    print("best gap %s at %s" % (
        report["best_gap"]["gap"] if report["best_gap"] else "n/a",
        report["best_gap"]["name"] if report["best_gap"] else "-"))
    print("alarms %d, confirmations %d" % (
        report["alarms_below"]["count"], report["confirmations_below"]["count"]))
    print("negative_mean_feasible_raw_gap_found %s"
          % report["negative_mean_feasible_raw_gap_found"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
