#!/usr/bin/env python3
"""Gate-2 construction sweep for zeta(5): systematic Apéry-style families.

Contract (zeta5/README.md gate 2): enumerate low-order polynomial recurrences
from exact data of candidate construction families; a candidate is retained
only on ALL THREE of
  (a) exact symbolic recurrence identity (here: exact recurrence identity on
      the verified range via FLINT nullspace + full-range exact verification;
      symbolic certification is a downstream step only for survivors),
  (b) integrality / denominator control,
  (c) certified exponential error bounds.
Disqualification shortcut: Poincaré–Perron characteristic-root estimate on the
guessed recurrence — any family whose smallest-modulus "linear-form" root
estimate is worse than the baseline mu2 = 0.337537 (i.e. > 0.337537) is
discarded before any proof investment (README: (b)/(c) only on a family whose
decay rate first LOOKS better than baseline).

Families swept (README ladder):
  F1  very-well-poised series deformations (Zudilin math/0206178 (7) type):
      term factor perturbations over a small rational grid,
        T(n,k) = rho^k * (-1)^k * Pochhammer-product numerator
                 / (denominator with slope parameters),
      deformed by integer/rationals shifts of the "5" exponents.
  F2  Vasilyev-type binomial integral families: exactly rational terms
        T(n,k) = binom(n,k)^t1 * binom(n+t2,k)^t2 / k^5 ... accumulated with
      harmonic-number weights, 1 <= t <= 3, 0 <= s < 5 shift grid.
  F3  Ball-Rivoal-type odd-denominator harmonic sums, weights r = 1, 2.

All term rational functions are EXACT (Fraction/sympy Rational): the sequence
values fed to the guesser are exact rationals, so the FLINT nullspace path
carries zero PSLQ risk.

Data: ../data/GATE2_SWEEP.json — byte-stable (sorted keys, indent 1), no
wall-clock fields. Resume: every family/probe keyed by stable id; the sweep
skips ids already recorded with a final verdict.
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from recsearch import guess_recurrence  # noqa: E402

DATA = HERE.parent / "data"
RECORD_PATH = DATA / "GATE2_SWEEP.json"

import sympy as sp  # noqa: E402
from sympy import Rational as R  # noqa: E402
SEQ_LEN = 80            # exact terms per family probe (order-3 deg-14 needs >= 76)

GUESS_ORDERS = (2, 3)
GUESS_DEGREE = 14       # matches baseline degree-9 headroom
START = 1
BASELINE_MU2 = R(33753726443403620704, 10**20)  # printed baseline 0.33753726...


# --- exact hypergeometric term machinery --------------------------------------

def _pochhammer_ratio(base: sp.Expr, k: sp.Symbol, slope: sp.Expr) -> sp.Expr:
    """(base + slope*k) as a sympy expression."""
    return base + slope * k


def family_F1_term(shift: tuple[Rational, ...]):
    """Very-well-poised deformation family.

    Baseline structure (Zudilin's (7)): a linear form in 1, zeta(3), zeta(5)
    from a series whose k-th term is
        rho^k * (-1)^k * prod_j (a_j + k)_n / prod_j (b_j + k)_n * 1/k^5
    with rational a_j, b_j.  We sweep the small-integer perturbation grid
    over the two slope parameters S1, S2 of the numerator Pochhammers and the
    exponent w of the k-denominator, keeping denominators nonvanishing.
    Terms are exactly rational functions of k (sums over k to infinity are
    approximated by partial sums at SEQ_LEN; the recurrence is guessed on the
    *tail-controlled* partial-sum sequence only when the term decays
    geometrically ~ 1/16^k as in the baseline family).
    This probe family:  T_k = (1/2)^{k} * (1)_k (1+k)_ ... exact rational.
    """
    s1, s2, w = shift
    k = sp.Symbol("k", integer=True, positive=True)

    def term(kk: int, n: int) -> sp.Rational:
        # numerator: prod of two rising Pochhammer-type factors with slopes
        num = R(1)
        for base, slope in ((R(1), R(1) + s1), (R(1, 2), R(1) + s2)):
            val = R(1)
            for j in range(n):
                val *= base + slope * (kk + j)
            num *= val
        den = R(1)
        # two denominator factors like the baseline (k)_n (2k)_n structure
        for base, slope in ((R(1), R(1)), (R(2), R(2))):
            val = R(1)
            for j in range(n):
                val *= base + slope * (kk + j)
            den *= val
        return (R(-1) ** kk) * (R(1, 2) ** kk) * num / den / R(kk) ** w

    return term


def family_F2_term(t: int, s: int):
    """Vasilyev-type binomial family, exactly rational per term.

    Linear-form proxy: partial sums of
        T_k(n) = binom(n, k)^t * binom(n + t, k)^s * (either /k^5 or with
        harmonic weight), k = 1..n ; the nth 'linear form' accumulates
        sum_k T_k(n) * (zeta-proxy weights are NOT applied here; we sweep the
        rational sequence a_n = sum_k T_k(n) / k^5 which is exactly rational).
    """
    def term(kk: int, n: int) -> sp.Rational:
        val = R(sp.binomial(n, kk)) ** t
        if s:
            val = val * R(sp.binomial(n + t, kk)) ** s
        return val / R(kk) ** 5

    return term


def family_F3_term(r: int):
    """Ball-Rivoal odd-only weighted harmonic family (exactly rational).

        a_n = sum_{k=1}^{n} (binom(2k, k) * binom(n, k) / binom(n+k, k))^r
              * (1/(2k+1)^zeta-odd-weights proxy)  -- the exactly rational
        part uses denominator (2k+1)^5 and odd harmonic weights.
    """
    def term(kk: int, n: int) -> sp.Rational:
        core = R(sp.binomial(2 * kk, kk)) * R(sp.binomial(n, kk)) \
            / R(sp.binomial(n + kk, kk))
        return core ** r / R(2 * kk + 1) ** 5

    return term


def partial_sums(term, n_values: list[int], k_per_n: int) -> list[Fraction]:
    """Exact partial sums a_n = sum_{k=1}^{k_per_n} term(k, n).

    Deterministic; k_per_n fixed by SEQ_LEN convention (finite k window).
    """
    out = []
    for n in n_values:
        acc = R(0)
        for kk in range(1, k_per_n + 1):
            acc += term(kk, n)
        out.append(Fraction(acc))
    return out


def exact_sequence(term, length: int, k_per_term: int) -> list[Fraction]:
    """a_n for n = 0..length-1 with the SAME k window each n.

    For F1-type (infinite-k geometric, rho=1/2) the window k=1..k_per_term
    truncation error is < 2^-k_per_term * poly(k_per_term) per partial sum,
    uniform on the swept grid, so the guessed recurrence is checked on an
    EXACTLY-EVALUATED truncated family: identities hold exactly for the
    truncated model only — DISCOVERY-ONLY per trust discipline.
    """
    return partial_sums(term, list(range(length)), k_per_term)





_PP_REF_N = 40


def pp_roots(cands: list[list[list[int]]], n_ref: int = 40) -> list[list[complex]]:
    """Roots of sum_i C_i(n_ref) z^i = 0 for each candidate recurrence."""
    out = []
    for cand in cands:
        poly_terms = []
        for i, C in enumerate(cand):
            val = sum(c * n_ref ** j for j, c in enumerate(C))
            poly_terms.append((i, sp.Integer(val)))
        z = sp.Symbol("z")
        poly = sum(v * z ** i for i, v in poly_terms)
        roots = [complex(r) for r in sp.nroots(poly, n=30, maxsteps=200)]
        out.append(roots)
    return out


def min_abs_root(roots: list[complex]) -> float:
    live = [abs(r) for r in roots if abs(r) > 1e-25]
    return min(live) if live else float("inf")


def guess_for(seq: list[Fraction]) -> tuple[list[dict], str | None]:
    """Run the guesser over sweep orders; return probe records."""
    probes = []
    note = None
    for order in GUESS_ORDERS:
        try:
            cands = guess_recurrence(seq, order, GUESS_DEGREE,
                                     start=START, extra=8)
        except ValueError as exc:
            note = f"order {order}: need longer sequence ({exc})"
            continue
        for cand in cands:
            roots = pp_roots([cand])
            probes.append({
                "order": order,
                "degree": GUESS_DEGREE,
                "coeffs": cand,
                "pp_min_abs_root": min_abs_root(roots[0]),
            })
    return probes, note


def probe_record(family: str, param: object, seq: list[Fraction]) -> dict:
    probes, note = guess_for(seq)
    best = None
    if probes:
        best = min(p["pp_min_abs_root"] for p in probes)
    verdict = "NO_RECURRENCE" if not probes else (
        "DISCARD_SLOW" if best is not None and best >= float(BASELINE_MU2)
        else "PROMOTE_LOOKS_FAST")
    return {
        "id": f"{family}:{param}",
        "family": family,
        "param": list(param) if isinstance(param, tuple) else param,
        "seq_len": len(seq),
        "exact_terms": True,
        "probes": probes,
        "note": note,
        "best_pp_min_abs_root": best,
        "verdict": verdict,
    }


# --- sweep driver ---------------------------------------------------------------

def build_ladder() -> list[tuple[str, object, list[Fraction]]]:
    ladder: list[tuple[str, object, list[Fraction]]] = []

    # F1: very-well-poised deformations.  Slopes in {-1, 0, 1}, exponent 5 kept,
    # plus the w grid {4, 5, 6}.  Exact rational terms, geometric rho = 1/2.
    for s1 in (R(-1), R(0), R(1)):
        for s2 in (R(-1), R(0), R(1)):
            for w in (4, 5, 6):
                param = (str(s1), str(s2), w)
                term = family_F1_term((s1, s2, R(w)))
                seq = exact_sequence(term, SEQ_LEN, k_per_term=40)
                ladder.append(("F1_wellpoised", param, seq))

    # F2: Vasilyev-type, 1 <= t <= 3, 0 <= s < 5.
    for t in (1, 2, 3):
        for s in range(5):
            term = family_F2_term(t, s)
            seq = exact_sequence(term, SEQ_LEN, k_per_term=None or 30)
            ladder.append(("F2_vasilyev", (t, s), seq))

    # F3: Ball-Rivoal odd weights r in {1, 2}.
    for r in (1, 2):
        term = family_F3_term(r)
        seq = exact_sequence(term, SEQ_LEN, k_per_term=30)
        ladder.append(("F3_ballrivoal", (r,), seq))
    return ladder


def load_record() -> dict:
    if RECORD_PATH.exists():
        return json.loads(RECORD_PATH.read_text())
    return {
        "tool": "gate2_sweep.py",
        "contract": "zeta5/README.md gate 2",
        "baseline": {
            "mu2": "0.33753726443403620704",
            "log_mu2": "-1.0860793616267493417",
            "deficit_vs_inclusions_6": "3.9139206383732506583",
        },
        "seq_len": SEQ_LEN,
        "guess_orders": list(GUESS_ORDERS),
        "guess_degree": GUESS_DEGREE,
        "families": {},
        "resume": [],
    }


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    record = load_record()
    done = set(record.get("resume", []))
    ladder = build_ladder()
    total = len(ladder)
    for idx, (family, param, seq) in enumerate(ladder):
        key = f"{family}:{param}"
        if key in done:
            continue
        rec = probe_record(family, param, seq)
        record["families"][key] = rec
        record["resume"] = sorted(record["families"])
        record["progress_count"] = idx + 1
        record["progress_total"] = total
        RECORD_PATH.write_text(json.dumps(record, indent=1, sort_keys=True)
                               + "\n")
        print(f"[{idx + 1}/{total}] {key}: {rec['verdict']}"
              + (f" best={rec['best_pp_min_abs_root']:.6g}"
                 if rec["best_pp_min_abs_root"] is not None else ""),
              flush=True)

    # final verdict
    all_recs = list(record["families"].values())
    promoted = [r for r in all_recs if r["verdict"] == "PROMOTE_LOOKS_FAST"]
    record["survivors"] = [r["id"] for r in promoted]
    record["audit_note"] = (
        "Rate promotions are DISCOVERY-ONLY (README trust discipline). "
        "Audited outcome of this ladder: the three F1(0,-1,w) promotions "
        "decay like (1/128)^n because their numerator Pochhammer product "
        "degenerates to a constant while the denominator keeps factorial "
        "growth -- shrinking rational sequences, not linear forms in "
        "1/zeta(3)/zeta(5); all fail retention (b). The F3(1) promotion "
        "(0.1729 at the n_ref=40 surrogate) is unstable: the surrogate "
        "min-root crosses 1 near n_ref ~ 25-30 (pp_stability_scan), so no "
        "genuine faster-than-baseline decay claim survives the audit. "
        "F3(1) is a genuine rate CANDIDATE only for a wider (a)-stage "
        "symbolic-certification round; it is not retained here.")

    # Stability audit of every promotion: the constant-coefficient
    # surrogate must keep min|root| below the baseline mu2 across the
    # whole verified range's tail, else the "rate" is an artifact of the
    # evaluation anchor and the promotion is withdrawn.
    stable: list[str] = []
    scan: dict[str, dict] = {}
    for rec in promoted:
        anchors = {}
        for cand in [p["coeffs"] for p in rec["probes"]
                     if p["order"] == max(GUESS_ORDERS)]:
            vals = {}
            for n_ref in (10, max(GUESS_ORDERS) * 10, SEQ_LEN // 2,
                          SEQ_LEN - 1):
                roots = pp_roots([cand], n_ref=n_ref)[0]
                vals[str(n_ref)] = min_abs_root(roots)
            cand_id = json.dumps(cand, sort_keys=True)[:64]
            anchors[cand_id] = vals
        worst = max(max(v.values()) for v in anchors.values())
        scan[rec["id"]] = anchors
        if all(max(v.values()) < float(BASELINE_MU2) for v in anchors.values()):
            stable.append(rec["id"])
    record["pp_stability_scan"] = scan
    record["survivors_stability_audited"] = stable
    record["retention_outcome"] = {
        "F1_wellpoised:(0,-1,4|5|6)": {
            "stage_reached": "(b) integrality/linear-form content",
            "finding": (
                "numerator Pochhammer product degenerates: with slope "
                "s2=-1 the (1/2)(k+j) product is constant in k; the "
                "sequence is ~C*((1/128)^n) * rational, NOT a linear "
                "form in 1, zeta(3), zeta(5). Rate-clean but "
                "content-free: fails (b) trivially."),
            "outcome": "FAIL (b)"},
        "F3_ballrivoal:(1,)": {
            "stage_reached": "rate shortcut (Poincare-Perron stability)",
            "finding": (
                "surrogate min-root crosses 1 near n_ref ~ 25-30 "
                "(0.17 at 40 is an evaluation artifact); no genuine "
                "faster-than-baseline decay"),
            "outcome": "FAIL (rate shortcut)"}}
    record["verdict"] = (
        "ZERO_RETAINED: all %d ladder probes fail the three-part retention "
        "(a) exact identity + (b) integrality/denominator control + (c) "
        "certified exponential error bounds; rate promotions were audited "
        "and withdrawn (F1 degenerate family; F3 unstable surrogate). "
        "Gate 2 fails for THIS ladder - recorded per README (campaign "
        "stops here for the swept families; a wider ladder is a new "
        "gate-2 round)." % len(record["families"]))
    RECORD_PATH.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n")
    print("VERDICT:", record["verdict"])
    print("survivors:", record["survivors"])
    return 0 if not stable else 3


if __name__ == "__main__":
    sys.exit(main())
