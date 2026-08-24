r"""e220 -- Directed-bound audit for the susceptibility-coefficient route.

Audits every Fekete / radius-of-convergence direction for the exact
simple-cubic susceptibility coefficients a_n (theorem-grade n <= 10 from
e218, [EXTERNAL] published values through n = 32), and proves the
finite-verification no-go by exhibiting two explicit completions of the
verified data with contradictory endpoint conclusions.

Direction table (proofs/susceptibility_coeffs.md, Theorem 2):

  * positivity a_n > 0 (all n) together with submultiplicativity (all m, n)
    ==> Fekete: lim a_n^{1/n} = inf_n a_n^{1/n}, so every exact a_n gives
    an UPPER bound on the growth constant and hence a candidate LOWER
    endpoint K_c >= atanh(a_n^{-1/n}), conditional also on the bridge
    premise P1' (|[v^n] chi_box| <= a_n uniformly, tested by e219, open).
    Status: all three premises [CONJECTURE]; no bound is claimed here.
  * supermultiplicativity (all m, n) ==> Fekete would give candidate UPPER
    endpoints.  Status: REFUTED at the smallest pair (1, 1): a_2 = 30 < 36
    (e219, theorem-grade), and at every anchored variant.
  * an all-n ratio floor a_{n+1} >= rho a_n would give a candidate upper
    endpoint atanh(1/rho) conditional on complex-disk analyticity; no
    finite prefix certifies any floor because the observed ratios are
    nonincreasing, with one early plateau, throughout the checked range.
    Status: [UNRESOLVED], no finite certificate possible from prefix data.
  * the proved ratio ceiling direction a_{n+1} <= 6 a_n (checked range)
    would, even if proved for all n, give only radius >= 1/6 and
    K_c >= atanh(1/6) = 0.1682..., strictly weaker than the certified
    endpoint: that route is valueless.

All enclosures below are computed in exact rational arithmetic: integer
n-th roots bracket a_n^{-1/n} between adjacent rationals with denominator
10^40, and atanh is enclosed by its Taylor partial sum plus an exact tail
bound (atanh(x) - sum_{k<=K} x^{2k+1}/(2k+1) <= x^{2K+3}/((2K+3)(1-x^2))).
No floating point enters any stored value.

The incumbent certified lower endpoint used for comparisons is the wave-14
repository constant 0.2122159753270231627267517174278577806020 (README,
proofs/saw_union4.md).  The benchmark value 0.221654626 is never used.
"""

from __future__ import annotations

import json
import resource
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

SCRIPT = "experiments/e220_coeff_bound.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = ROOT / "results" / "bounds" / "susceptibility_coeffs.json"

SCALE_DIGITS = 40
SCALE = 10**SCALE_DIGITS
ATANH_TERMS = 70  # x <= 0.25 => tail < 0.25^143/143 < 10^-87
INCUMBENT_LOWER = Fraction(
    2122159753270231627267517174278577806020, 10**40
)  # repo-certified wave-14 lower endpoint (README; proofs/saw_union4.md)
CLOSURE_HORIZON = 64


def introot(x: int, n: int) -> int:
    """floor(x**(1/n)) for x >= 0, n >= 1, exact."""
    if x < 0 or n < 1:
        raise ValueError
    if x in (0, 1) or n == 1:
        return x
    r = 1 << (-(-x.bit_length() // n))  # upper seed: 2^ceil(bits/n) >= x^(1/n)
    while True:
        nr = ((n - 1) * r + x // r ** (n - 1)) // n
        if nr >= r:
            break
        r = nr
    while r**n > x:
        r -= 1
    while (r + 1) ** n <= x:
        r += 1
    return r


def bracket_inverse_root(a_n: int, n: int) -> tuple[Fraction, Fraction]:
    """Adjacent rationals lo < a_n^{-1/n} <= hi with denominator SCALE."""
    q = SCALE**n // a_n
    r = introot(q, n)
    lo, hi = Fraction(r, SCALE), Fraction(r + 1, SCALE)
    assert lo**n * a_n <= 1 and hi**n * a_n > 1, "bracket sanity"
    return lo, hi


def atanh_enclosure(x: Fraction) -> tuple[Fraction, Fraction]:
    """Exact rational lo <= atanh(x) <= hi for 0 <= x < 1/2."""
    assert 0 <= x < Fraction(1, 2)
    total = Fraction(0)
    xsq = x * x
    power = x
    for k in range(ATANH_TERMS + 1):
        total += power / (2 * k + 1)
        power *= xsq
    # power is now x^(2*ATANH_TERMS+3); exact tail bound (all terms positive)
    tail = power / ((2 * ATANH_TERMS + 3) * (1 - xsq))
    return total, total + tail


def dec(f: Fraction, places: int, direction: str) -> str:
    """Directed decimal rendering of a rational."""
    s = 10**places
    n = f.numerator * s
    d = f.denominator
    q, rem = divmod(n, d)
    if direction == "up" and rem:
        q += 1
    if direction not in ("up", "down"):
        raise ValueError
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q // s}.{q % s:0{places}d}"


def main() -> None:
    checks: list[dict[str, object]] = []
    t_start = time.process_time()

    def record(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        print(("PASS" if passed else "FAIL"), name + ":", detail)
        if not passed:
            raise AssertionError(f"check failed: {name}")

    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    a_own = [int(s) for s in payload["data"]["coefficients"]["a"]]
    a_ext = [int(s) for s in payload["data"]["external"]["oeis_A002913"]]
    N_own = len(a_own) - 1
    N_ext = len(a_ext) - 1

    # -- self-validation of the enclosure machinery -------------------------
    lo, hi = atanh_enclosure(Fraction(1, 5))
    # atanh(1/5) = 0.20273255405408219098900655773217... (proofs/upper_infrared.md)
    ref_lo = Fraction(20273255405408219098900655773217, 10**32)
    ref_hi = ref_lo + Fraction(1, 10**31)
    record(
        "atanh_enclosure_control",
        ref_lo <= lo and hi <= ref_hi and hi - lo < Fraction(1, 10**60),
        "rational atanh enclosure of atanh(1/5) lies inside the "
        "independently recorded repository bracket "
        "[0.20273255405408219098900655773217, +10^-31] "
        "(results/bounds/upper_infrared.json) and is tighter than 10^-60",
    )
    r = introot(10**80 // 7, 3)
    record(
        "introot_control",
        r**3 <= 10**80 // 7 < (r + 1) ** 3 and introot(3883554493872938687622, 1)
        == 3883554493872938687622,
        "integer n-th root: floor bracket verified on a 80-digit cube and "
        "the identity root",
    )

    # -- conditional lower-endpoint targets ---------------------------------
    targets = {}
    beats = {}
    for n in range(1, N_ext + 1):
        blo, bhi = bracket_inverse_root(a_ext[n], n)
        alo, _ = atanh_enclosure(blo)
        _, ahi = atanh_enclosure(bhi)
        beats[n] = alo > INCUMBENT_LOWER
        if not beats[n]:
            assert ahi < INCUMBENT_LOWER or alo <= INCUMBENT_LOWER <= ahi
        if n in (N_own, 25, N_ext):
            targets[n] = {
                "bracket_v": [str(blo), str(bhi)],
                "atanh_enclosure_40dp": [dec(alo, 40, "down"), dec(ahi, 40, "up")],
                "grade": (
                    "theorem-grade coefficients (finite-graph derivation e218)"
                    if n <= N_own
                    else "[EXTERNAL] peer-reviewed table values"
                    if n <= 25
                    else "[EXTERNAL] workshop-slide values (weaker grade)"
                ),
            }
    crossing = min(n for n in range(1, N_ext + 1) if beats[n])
    record(
        "conditional_targets_certified_brackets",
        all(
            Fraction(targets[n]["atanh_enclosure_40dp"][0])
            < Fraction(targets[n]["atanh_enclosure_40dp"][1])
            for n in targets
        ),
        "directed 40-decimal enclosures of atanh(a_n^(-1/n)) for "
        f"n in {sorted(targets)} (hypothetical targets, not claimed bounds)",
    )
    record(
        "own_range_conditional_value_insufficient",
        not beats[N_own],
        f"even if positivity, submultiplicativity, and the bridge premise were theorems, "
        f"the theorem-grade range n <= {N_own} yields at most "
        f"atanh(a_{N_own}^(-1/{N_own})) <= "
        f"{targets[N_own]['atanh_enclosure_40dp'][1][:14]}... < incumbent "
        "0.2122159753...: the route cannot beat the certified endpoint from "
        "own-derived coefficients",
    )
    record(
        "external_range_conditional_value",
        beats[25] and beats[N_ext] and crossing == 21,
        "the conditional route would beat the incumbent from order n = 21 "
        "on ([EXTERNAL] coefficients); at n = 25 (peer-reviewed grade) the "
        f"target is {targets[25]['atanh_enclosure_40dp'][0][:13]}..., at "
        f"n = 32 {targets[N_ext]['atanh_enclosure_40dp'][0][:13]}...; both "
        "remain conditional on three unproved premises and are NOT claimed",
    )

    # -- ratio-ceiling route is valueless ------------------------------------
    alo6, ahi6 = atanh_enclosure(Fraction(1, 6))
    record(
        "ratio_ceiling_route_valueless",
        ahi6 < INCUMBENT_LOWER,
        f"a proved all-n ceiling a_(n+1) <= 6 a_n would give only "
        f"K_c >= atanh(1/6) <= {dec(ahi6, 20, 'up')} < 0.2122159753...: "
        "strictly weaker than the certified endpoint",
    )

    # -- no-go completions ----------------------------------------------------
    s = list(a_ext)
    for n in range(N_ext + 1, CLOSURE_HORIZON + 1):
        s.append(min(s[m] * s[n - m] for m in range(1, n)))
    closure_ok = all(
        s[m + n] <= s[m] * s[n]
        for m in range(1, CLOSURE_HORIZON)
        for n in range(m, CLOSURE_HORIZON + 1 - m)
    )
    record(
        "closure_completion_submultiplicative",
        closure_ok and s[: N_ext + 1] == a_ext,
        f"the min-plus closure extension of the verified prefix is fully "
        f"submultiplicative on every pair with m+n <= {CLOSURE_HORIZON} and "
        "by construction for all larger n (Theorem 4 of the proof note): "
        "the verified data can never refute submultiplicativity",
    )
    d = list(a_ext) + [6**n for n in range(N_ext + 1, CLOSURE_HORIZON + 1)]
    d_violates_at = next(
        (m, n)
        for m in range(1, CLOSURE_HORIZON)
        for n in range(m, CLOSURE_HORIZON + 1 - m)
        if d[m + n] > d[m] * d[n]
    )
    record(
        "divergent_completion_consistent_with_checks",
        d[: N_ext + 1] == a_ext and d_violates_at == (1, 32)
        and d[33] > d[1] * d[32],
        "the completion d_n = a_n (n <= 32), d_n = 6^n (n > 32) agrees with "
        "every verified value and every verified pair inequality (all have "
        "m+n <= 32), yet has limsup d_n^(1/n) = 6 and radius 1/6; its first "
        "submultiplicativity violation is (m,n)=(1,32), at total order 33, "
        "and is invisible to any finite check through order 32: no endpoint "
        "conclusion follows from the verified constraints alone",
    )
    # growth envelope (Lemma 6 of the note): E_n = (1000 (n+1)^3)^(n+2),
    # proved from the ball(2n) stabilization plus the truncated-division
    # composition bound |[v^j] 1/Z| <= 2^(j-1) E^j
    def envelope(n: int) -> int:
        return (1000 * (n + 1) ** 3) ** (n + 2)

    env_ok = all(
        a_ext[n] <= envelope(n) for n in range(1, N_ext + 1)
    ) and all(6**n <= envelope(n) for n in range(1, CLOSURE_HORIZON + 1)) and all(
        s[n] <= envelope(n) for n in range(1, CLOSURE_HORIZON + 1)
    )
    record(
        "envelope_consistency",
        env_ok,
        "both completions and all exact values respect the provable "
        "all-n envelope |a_n| <= (1000 (n+1)^3)^(n+2) (Lemma 6), so the "
        "no-go survives adding every rigorously known all-n constraint",
    )

    elapsed = time.process_time() - t_start
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    record(
        "resource_budget",
        elapsed < 600.0 and peak_rss < 2 * 1024**3,
        f"process time {elapsed:.1f}s < 600s, peak RSS "
        f"{peak_rss/1024**2:.0f} MiB < 2048 MiB",
    )

    payload["meta"]["e220"] = {
        "script": SCRIPT,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "interpreter": ".venv/bin/python",
        "arithmetic": "[COMPUTATION] exact integer n-th-root brackets and "
        "rational atanh Taylor enclosures with exact tail bounds; no floats",
        "elapsed_process_seconds": elapsed,
        "peak_rss_bytes": peak_rss,
    }
    payload["data"]["directed"] = {
        "status": "[THEOREM] direction table and no-go completions; "
        "[CONJECTURE] all three premises of the conditional lower-endpoint "
        "route; [COMPUTATION] certified enclosures. No new K_c bound is "
        "claimed: no all-n inequality is proved.",
        "claimed_new_bound": None,
        "reason_no_bound": "all-n positivity, submultiplicativity, and the "
        "bridge premise P1' are finite evidence only; the task rule 'only "
        "compute a new Kc bound if an all-n inequality is actually proved' "
        "applies",
        "incumbent_lower_endpoint": dec(INCUMBENT_LOWER, 40, "down"),
        "incumbent_source": "README.md / proofs/saw_union4.md (wave-14)",
        "conditional_route": {
            "premises": [
                "(G0) [CONJECTURE] a_n > 0 for every n",
                "(G1) [CONJECTURE] a_(m+n) <= a_m a_n for all m, n >= 1",
                "(P1') [CONJECTURE] |[v^n] chi_box| <= a_n for every finite "
                "box, origin, and n",
            ],
            "conclusion_if_all_held": "K_c >= atanh(a_n^(-1/n)) for every "
            "exactly known a_n (proof note, Theorem 5)",
            "targets": {str(n): targets[n] for n in sorted(targets)},
            "first_order_beating_incumbent": crossing,
            "beats_incumbent_by_order": {
                str(n): bool(beats[n]) for n in range(1, N_ext + 1)
            },
        },
        "supermultiplicativity_route": {
            "status": "[THEOREM] refuted",
            "smallest_counterexample": "(m,n)=(1,1): a_2 = 30 < 36 = a_1^2",
            "consequence": "no Fekete sup-direction upper endpoint exists "
            "for this sequence",
        },
        "ratio_floor_route": {
            "status": "[UNRESOLVED] no finite prefix can certify an all-n "
            "ratio floor (ratios are nonincreasing on the checked range, with "
            "a_2/a_1 = a_3/a_2 = 5); would additionally need complex-disk "
            "analyticity of chi up to the radius (open)",
        },
        "ratio_ceiling_route": {
            "status": "[COMPUTATION] valueless even if proved",
            "atanh_one_sixth_upper": dec(ahi6, 40, "up"),
        },
        "completions": {
            "closure_tail_33_to_40": [str(x) for x in s[33:41]],
            "closure_horizon": CLOSURE_HORIZON,
            "divergent_tail_rule": "d_n = 6^n for n > 32",
            "divergent_first_violation_pair": list(d_violates_at),
            "envelope": "E_n = (1000 (n+1)^3)^(n+2) (Lemma 6)",
        },
    }
    fresh = {c["name"] for c in checks}
    payload["checks"] = [
        c for c in payload["checks"] if c["name"] not in fresh
    ] + checks
    RESULT_PATH.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"extended {RESULT_PATH.relative_to(ROOT)}; {elapsed:.1f} CPU s")


if __name__ == "__main__":
    main()
