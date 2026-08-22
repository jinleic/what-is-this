#!/usr/bin/env python3
"""Replay the fixed-class implication Schinzel H => operational ladder-zero H.

This is an exact, lightweight audit.  It checks the 293 ledgered class pairs,
then replays 24 prime-rung closure members through the kernel ladder.  It does
not run a search and does not treat a primality refusal as evidence.
"""
from __future__ import annotations

import json
import math
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import h10q
from l13_filter import cofactor_decide, smooth_emergent


HERE = Path(__file__).resolve().parent
CLOSURES = HERE / "data" / "l13h_all_closures.json"
SCHINZEL = HERE / "data" / "l17_schinzel_audit.jsonl"
REMAINDERS = HERE / "data" / "l15_remainders.jsonl"
COFACTOR = HERE / "data" / "l17_cofactor.jsonl"
OUT = HERE / "data" / "l18_schinzel_implies_h.jsonl"
REPORT = Path("/tmp/l18_schinzel_implies_h.md")
SAMPLE_PER_FAMILY = 12


def gcd_many(values: list[int]) -> int:
    g = 0
    for value in values:
        g = math.gcd(g, abs(value))
    return g


def lcm(a: int, b: int) -> int:
    return a // math.gcd(a, b) * b


def primitive_integer_poly(coeffs: list[F]) -> tuple[list[int], F]:
    """pi, scale with coeff-polynomial = scale*pi and pi primitive."""
    den = 1
    for coeff in coeffs:
        den = lcm(den, coeff.denominator)
    ints = [int(coeff * den) for coeff in coeffs]
    content = gcd_many(ints)
    assert content
    return [value // content for value in ints], F(content, den)


def compose(P: list[F], b0: F, step: F) -> list[F]:
    out = [F(0)] * len(P)
    for i, coeff in enumerate(P):
        if not coeff:
            continue
        for j in range(i + 1):
            out[j] += coeff * math.comb(i, j) * b0 ** (i - j) * step ** j
    return out


def evaluate(coeffs: list[int] | list[F], t: int | F) -> int | F:
    value = 0
    for coeff in reversed(coeffs):
        value = value * t + coeff
    return value


def remove_support(value: F, support: set[int] | list[int]) -> F:
    out = F(value)
    for p in support:
        out /= F(p) ** h10q.vp(out, p)
    return out


def strip_integer_parts(value: F, support: set[int] | list[int]) -> tuple[int, int]:
    rn, rd = abs(value.numerator), value.denominator
    for p in support:
        while rn % p == 0:
            rn //= p
        while rd % p == 0:
            rd //= p
    return rn, rd


def setup(rec: dict) -> dict:
    a = F(rec["a"])
    w, ut = rec["cell"]
    z = F(w) * F(*ut)
    A = 1 + 4 * a * a
    if rec["family"] == "L11":
        tau = (1 + 2 * a * a) / A
    elif rec["family"] == "ESC":
        tau = F(3, 5)
    else:
        raise AssertionError(rec["family"])
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z = z ** 3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    eps, f = int(rec["eps"]), int(rec["f"])
    q1, N = int(rec["q1"]), int(rec["N"])
    P = h10q._l10_P(a, Z, D, A, delta, s)
    S = sorted(
        {2, 3, 5, 7}
        | h10q._l10_supp(alpha)
        | h10q._l10_supp(delta)
        | h10q._l10_supp(F(f))
    )
    return {
        "a": a,
        "z": z,
        "A": A,
        "tau": tau,
        "delta": delta,
        "alpha": alpha,
        "Z": Z,
        "D": D,
        "s": s,
        "eps": eps,
        "f": f,
        "q1": q1,
        "N": N,
        "P": P,
        "S": S,
    }


def fixed_class_audit(records: list[dict], sch_rows: list[dict]) -> tuple[dict, list[dict]]:
    assert len(records) == len(sch_rows) == 293
    details: list[dict] = []
    fallback: list[dict] = []
    family_counts = Counter()

    for index, (rec, prior) in enumerate(zip(records, sch_rows)):
        d = setup(rec)
        family_counts[rec["family"]] += 1
        assert prior["record_index"] == index
        assert prior["family"] == rec["family"] and prior["cell"] == rec["cell"]
        assert prior["q1_0"] == d["q1"] and int(prior["N"]) == d["N"]
        assert d["delta"] and d["D"] and d["a"] and d["z"]
        assert d["eps"] in (-1, 1) and d["f"] > 0
        assert h10q.vp(d["a"], 2) == 0 and h10q.vp(F(d["f"]), 2) == 0
        assert d["N"] > 0 and d["N"] % 8 == 0
        assert math.gcd(d["q1"], d["N"]) == 1
        assert h10q._is_prime(d["q1"])
        if d["f"] != 1:
            assert h10q._is_prime(d["f"])

        b0 = F(d["eps"] * d["f"] * d["q1"])
        step = F(d["eps"] * d["f"] * d["N"])
        F_rat = compose(d["P"], b0, step)
        G, scale = primitive_integer_poly(F_rat)
        assert G == [int(x) for x in prior["F_coeffs_t"]]
        assert scale == F(prior["F_scale"])
        assert len(G) == 9 and G[-1] > 0 and gcd_many(G) == 1
        assert d["P"][-1] > 0

        # Exact polynomial meaning of "divide by the S-part": the rational
        # scale is entirely S-supported, and G(t) is an S-unit for every t.
        assert remove_support(scale, d["S"]) == 1
        for p in d["S"]:
            assert d["N"] % p == 0
            assert G[0] % p != 0
            assert all(coeff % p == 0 for coeff in G[1:])
        # At 2 the stronger congruence is needed for square classes.
        assert G[0] % 2 == 1 and all(coeff % 8 == 0 for coeff in G[1:])

        # Base frozen signs; the displayed polynomial congruences prove these
        # square classes persist for every integer t, not merely samples.
        c0 = h10q._sun_h(d["a"], b0, d["Z"])
        assert c0 is not None
        M0 = 16 - d["delta"] * c0 * c0 - 32 * d["A"] * b0 * d["s"] * d["s"]
        Pb0 = evaluate(d["P"], b0)
        assert M0 == Pb0 / (b0 ** 4 * d["D"] ** 2 * d["A"] ** 2)
        x0, d0 = d["alpha"] * M0, d["alpha"] * 2 * b0
        assert all(h10q.hilbert(x0, d0, p) == 1 for p in d["S"])

        # The moving symbol is (A|Q), fixed by Q == q1 (mod N); the unit
        # hypotheses exclude only finitely many Q dividing fixed data.
        assert h10q.legendre(d["A"], d["q1"]) == 1
        modulus_A = 4 * d["A"].numerator * d["A"].denominator
        assert d["N"] % modulus_A == 0

        nz = [i for i, coeff in enumerate(d["P"]) if coeff]
        degree = nz[-1]
        lead = d["P"][degree]
        sx = (1 if d["alpha"] > 0 else -1) * (1 if lead > 0 else -1) * d["eps"] ** degree
        sd = 1 if d["alpha"] * d["eps"] > 0 else -1
        infinity_symbol = -1 if sx < 0 and sd < 0 else 1
        assert infinity_symbol == 1

        # Recheck the standard Schinzel-H hypotheses.  Degree nine means ten
        # consecutive product values are already an exact fixed-divisor test.
        cert_p = int(prior["F_irreducibility_certificate_prime"])
        assert h10q._l13_irred8(G, cert_p)
        product_value_gcd = 0
        for t in range(10):
            product_value_gcd = math.gcd(
                product_value_gcd,
                abs((d["q1"] + d["N"] * t) * int(evaluate(G, t))),
            )
        assert product_value_gcd == prior["pair_product_gcd"] == 1
        root_resultant = sum(
            G[j] * (-d["q1"]) ** j * d["N"] ** (8 - j)
            for j in range(9)
        )
        assert root_resultant != 0  # q1(t) is not a factor of G(t)

        # Record where the generic Taylor exponent alone is not sufficient.
        # The exact displayed congruence of G above still proves constancy:
        # G(t)/G(0) is in 1+pZ_p (1+8Z_2), and Q(t)/q1 likewise.
        ks = {p: h10q._l10_exponent(d["P"], b0, p) for p in d["S"]}
        short = [
            {"p": p, "required_k": k, "step_vp": h10q.vp(step, p)}
            for p, k in ks.items()
            if h10q.vp(step, p) < k
        ]
        if short:
            assert all(item["p"] != 2 for item in short)
            fallback.append({"record_index": index, "cell": rec["cell"], "places": short})

        details.append(
            {
                "record_index": index,
                "cell": rec["cell"],
                "family": rec["family"],
                "S": d["S"],
                "scale": str(scale),
                "fallback_places": short,
            }
        )
        time.sleep(0.1)

    summary = {
        "canonical_records": len(records),
        "family_counts": dict(sorted(family_counts.items())),
        "q1_proven_prime_and_gcd_one": len(records),
        "G_degree_8_primitive_irreducible": len(records),
        "G_positive_leading_coefficient": len(records),
        "pair_product_fixed_divisor_one": len(records),
        "pair_polynomials_coprime": len(records),
        "scale_entirely_S_supported": len(records),
        "G_S_unit_for_every_integer_t": len(records),
        "frozen_base_symbols_all_plus_one": len(records),
        "moving_legendre_symbol_plus_one": len(records),
        "asymptotic_infinity_symbol_plus_one": len(records),
        "class_constancy_exact": len(records),
        "taylor_exponent_direct_classes": len(records) - len(fallback),
        "polynomial_congruence_fallback_classes": len(fallback),
        "polynomial_congruence_fallbacks": fallback,
        "fixed_divisor_certificate": "gcd of t=0..9 product values; exact by degree-9 finite differences",
    }
    return summary, details


def prime_rung_sample(remainder_rows: list[dict], records: list[dict]) -> tuple[list[dict], dict]:
    prime_rows = [row for row in remainder_rows if row["R_class"] == "prime"]
    assert len(prime_rows) == 197
    by_key = {(rec["family"], tuple([rec["cell"][0], *rec["cell"][1]])): rec for rec in records}
    selected: list[dict] = []
    for family in ("ESC", "L11"):
        candidates = sorted(
            (row for row in prime_rows if row["family"] == family),
            key=lambda row: (row["rn_digits"], row["cell"]),
        )
        selected.extend(candidates[:SAMPLE_PER_FAMILY])

    out: list[dict] = []
    for sample_index, row in enumerate(selected):
        key = (row["family"], tuple([row["cell"][0], *row["cell"][1]]))
        authority = by_key[key]
        for field in ("a", "eps", "f", "q1", "N", "k_zero", "Q"):
            assert str(row[field]) == str(authority[field])
        d = setup(row)
        Q = int(row["Q"])
        assert Q == d["q1"] + int(row["k_zero"]) * d["N"]
        assert h10q._is_prime(Q)
        b = F(d["eps"] * d["f"] * Q)

        status, detail = smooth_emergent(d["a"], d["z"], d["tau"], b)
        assert status == "cofactor-big"
        small_bad, (rn, rd), smooth = detail
        assert small_bad == [] and smooth is False
        assert rn == int(row["rn"]) and rd == int(row["rd"])
        verdict, info = cofactor_decide(d["a"], d["z"], d["tau"], b, detail=detail)
        assert verdict == "zero"
        assert info["rung"] == "prime"
        assert len(info["places"]) == 1
        place, sign, tier, rung = info["places"][0]
        assert int(place) == rn and sign == 1 and tier == "proved" and rung == "prime"
        kinds = sorted(part["kind"] for part in info["parts"])
        assert kinds == ["prime", "unit"]

        P_value = F(evaluate(d["P"], b))
        raw_rn, raw_rd = strip_integer_parts(P_value, set(d["S"]) | {Q})
        assert raw_rn % rn == 0 and raw_rd % rd == 0

        out.append(
            {
                "type": "verified_prime_rung",
                "label": "PROVED_PER_ROW",
                "sample_index": sample_index,
                "source_record": "data/l15_remainders.jsonl",
                "family": row["family"],
                "cell": row["cell"],
                "a": row["a"],
                "eps": int(row["eps"]),
                "f": int(row["f"]),
                "q1": int(row["q1"]),
                "N": str(row["N"]),
                "k_zero": int(row["k_zero"]),
                "Q": str(Q),
                "Q_prime_proved": True,
                "smooth_status": status,
                "small_bad_places": small_bad,
                "ladder_residual_rn": str(rn),
                "ladder_residual_rd": str(rd),
                "ladder_residual_single_proved_prime": True,
                "cofactor_place_hilbert": sign,
                "cofactor_tier": tier,
                "cofactor_verdict": verdict,
                "operational_emergent_free": True,
                "event_equality_on_row": True,
                "raw_frozen_strip_rn": str(raw_rn),
                "raw_frozen_strip_rd": str(raw_rd),
                "raw_frozen_strip_is_same_single_prime": raw_rn == rn and raw_rd == rd,
                "smooth_factor_removed_num": str(raw_rn // rn),
                "smooth_factor_removed_den": str(raw_rd // rd),
                "strip_note": "ladder residual additionally removes every prime <=1e6 after the frozen S union supp(b) strip",
            }
        )
        time.sleep(0.1)

    summary = {
        "prime_rung_authority_rows": len(prime_rows),
        "prime_rung_family_counts": dict(Counter(row["family"] for row in prime_rows)),
        "sample_rows": len(out),
        "sample_family_counts": dict(Counter(row["family"] for row in out)),
        "Q_prime_proved": sum(row["Q_prime_proved"] for row in out),
        "ladder_residual_single_proved_prime": sum(row["ladder_residual_single_proved_prime"] for row in out),
        "cofactor_hilbert_plus_one": sum(row["cofactor_place_hilbert"] == 1 for row in out),
        "cofactor_verdict_zero": sum(row["cofactor_verdict"] == "zero" for row in out),
        "event_equalities": sum(row["event_equality_on_row"] for row in out),
        "raw_frozen_strip_already_single_prime": sum(row["raw_frozen_strip_is_same_single_prime"] for row in out),
        "strip_scope_warning": "The 197-way prime-rung classification concerns the post-1e6 ladder residual. The raw S union supp(b) strip can retain additional <=1e6 +1-symbol factors.",
        "global_converse": False,
        "global_converse_reason": "96 factorint closures are also zero; singleton prime is sufficient, not necessary",
    }
    return out, summary


def write_report(class_summary: dict, sample_summary: dict, cofactor_summary: dict) -> None:
    fallback = class_summary["polynomial_congruence_fallback_classes"]
    raw_single = sample_summary["raw_frozen_strip_already_single_prime"]
    n_sample = sample_summary["sample_rows"]
    lines = [
        "# L18 — Schinzel H implies operational H for a fixed aligned class",
        "",
        "**Verdict — PROVED implication; CONDITIONAL conclusion.** Yes.  For a fixed verified aligned class, once the stripped class polynomial is normalized as a primitive integral polynomial `G`, Schinzel's Hypothesis H for `{q1+Nt, G(t)}` gives infinitely many prime members with ladder verdict `zero`.  The implication is a theorem; existence of the prime values remains CONDITIONAL on the named classical conjecture.",
        "",
        "This is a per-class result.  Schinzel H does not construct an aligned class for a cell, and it does not by itself prove that every future class polynomial is irreducible/admissible.  Those are separate finite hypotheses.",
        "",
        "## 1. Operational meaning and the terminology trap",
        "",
        "`smooth_emergent` first checks every place in `S ∪ supp(b)` and infinity, strips those primes from numerator and denominator, then removes primes up to `10^6`; it appends a place only when its valuation is odd **and** its Hilbert symbol is `-1` (`l13_filter.py:90-144`).  `cofactor_decide` proves the remaining prime factors and computes their Hilbert symbols; it returns `zero` only when every remaining place is proved and every sign is `+1` (a proved `-1` is `bad-*`, while an unproved cofactor remains Jacobi-tier; `l13_filter.py:281-366`).  Thus operational **emergent-free/zero-bad means no bad (-1) emergent place**, not “no odd valuation outside the frozen set.”",
        "",
        "This distinction is load-bearing.  A prime-rung row has one odd-valuation cofactor prime, but its sign is `+1`; L13f explicitly calls that zero-bad (`THEOREMS.md:1433-1444`), L14 contains 197 such closures (`THEOREMS.md:1481-1503`), and L17 records all 197 singleton signs as reciprocity-forced (`THEOREMS.md:1605-1615`).  The literal wording in `CONDITIONAL.md:103-109` (“every outside valuation even”) is therefore stronger than the kernel/L13/L14 contract and is contradicted by the 197 accepted prime-rung rows.  The theorem below uses the operational contract required by this assignment.  If the literal no-odd-valuation wording were retained, Schinzel prime values would **not** imply it: they deliberately create one odd prime.",
        "",
        "The polynomial is the kernel polynomial",
        "`P(b)=16 D^2 A^2 b^4 - delta a^4 Z^4 Ng(b)^2 - 32 A^3 s^2 D^2 b^5` (`h10q.py:2579-2588`).  At replay points `x=alpha P(b)/(D^2 A^2 b^4)` and `d=2 alpha b` (`l12_class.py:129-138`; normalization also `THEOREMS.md:1475-1479`).  Consequently, outside `S ∪ supp(b)`, parity of `v_p(x)` is parity of `v_p(P(b))`.",
        "",
        "## 2. Fixed-class theorem",
        "",
        "Let `Q(t)=q1+Nt`, `b(t)=eps*f*Q(t)`, and let",
        "`S={2,3,5,7} ∪ supp(alpha) ∪ supp(delta) ∪ supp(f)`.  Assume:",
        "",
        "1. the fixed data are nondegenerate, `a` is a 2-adic odd unit, `f,q1` are odd, `N>0`, `8|N`, and `gcd(q1,N)=1`;",
        "2. the class certificate makes every symbol at `S`, the moving prime `Q`, and infinity equal to `+1` for every sufficiently large prime member outside a finite excluded set;",
        "3. `F(t)=P(eps*f*(q1+Nt))=c G(t)`, where `c>0` is supported on `S`, `G∈Z[t]` is primitive, and `G(t)` is an `S`-unit for every integer `t` (this is the precise polynomial meaning of dividing by the frozen S-part); and",
        "4. `Q` and `G` are irreducible with positive leading coefficients and their product has no fixed prime divisor.",
        "",
        "Then Schinzel's Hypothesis H for `{Q,G}` implies that the class contains infinitely many operationally emergent-free members.",
        "",
        "### Proof",
        "",
        "Schinzel H supplies infinitely many positive integers `t` for which both `Q(t)` and `G(t)` are positive primes.  The side conditions remove only finitely many `t`: the Cauchy/sign bound, `b∉{0,1}`, and each forbidden equality `Q(t)=p` for a prime dividing fixed data.  Hence infinitely many Schinzel values remain eligible.",
        "",
        "For such a value, stripping `S` from `F(t)=cG(t)` removes `c` and removes nothing from the prime `G(t)`.  At the moving prime, `b≡0 (mod Q)` and the displayed kernel formula gives `P(b)≡-delta*a^4*Z^4*A^2 (mod Q)`, a unit under the finite-exclusion hypotheses; hence `v_Q(P(b))=0` and `G(t)≠Q(t)`.  Thus there is exactly one outside odd-valuation place `R=G(t)`.  All other places outside `S∪{Q}` have even `v_p(x)` and `v_p(d)=0`, hence symbol `+1` by the odd-prime Hilbert formula (`h10q.py:359-380`).",
        "",
        "Every place in `S`, the moving place, and infinity is `+1` by the aligned certificate.  Hilbert reciprocity now forces `(x,d)_R=+1`: L9 states precisely that if all places in the finite support set but one are `+1`, the last is `+1` (`THEOREMS.md:828-853`).  Its exact hypotheses here are `x,d≠0`, `R∉supp(d)`, one singleton odd place, and all frozen/moving/infinite signs `+1`.  Therefore every odd-valuation outside place has sign `+1`, so the ladder verdict is `zero`.",
        "",
        "The associated tied conic is mathematically soluble by W2/Hasse–Minkowski.  `_l7_tied_status` is only the computational wrapper and can return `budget` on factorization refusal (`h10q.py:1514-1522`); the exact conic it calls is in `h10q.py:833-844`.  A refusal is irrelevant to this proof because primality and reciprocity are hypotheses/theorems, not search evidence.  This proves the implication.  □",
        "",
        "## 3. Why the frozen signs hold for every t",
        "",
        "The generic certificate is genuinely uniform: `_l10_exponent` states that `b≡b0 (mod p^k)` fixes both `v_p(P(b))` and the local square class of `P(b)/P(b0)` (`h10q.py:2601-2608`); `_l10_class_cert` puts every such prime power, `8`, and `4*num(A)*den(A)` into `N`, rejects `gcd(q1,N)≠1`, and records the asymptotic sign (`h10q.py:2617-2688`).  L9a proves `(x,d)_Q=(A|Q)` under the factor-by-factor unit hypotheses (`THEOREMS.md:887-903`).  For escape classes the same every-member statement, the finite exclusions, moving symbol, and real bound are explicit in `l12_class.py:1-49,70-151`.",
        "",
        f"The present replay gives an independent exact check on all {class_summary['canonical_records']} canonical pairs.  In every row `F_scale` is entirely S-supported; the primitive `G` has nonzero constant coefficient modulo every `p∈S` and every higher coefficient is `0 mod p` (at `p=2`, `0 mod 8`).  Hence for **every integer t**, not a sample, `G(t)/G(0)∈1+pZ_p` (or `1+8Z_2`) and `Q(t)/q1` has the same property; both are local squares.  All {class_summary['canonical_records']} base frozen signs are `+1`, so all remain `+1` throughout the class.  The Taylor exponent alone directly covers {class_summary['taylor_exponent_direct_classes']} rows; {fallback} alternate ESC rows need the displayed normalized-polynomial congruence, which closes them exactly rather than by sampling.",
        "",
        "## 4. Schinzel hypotheses on the 293 audited pairs",
        "",
        f"**PROVED per row.** {class_summary['G_degree_8_primitive_irreducible']}/293 primitive degree-8 `G` polynomials replay their Frobenius irreducibility certificate (`h10q.py:4083-4095`); {class_summary['G_positive_leading_coefficient']}/293 have positive leading coefficient; all linear `Q` have positive lead and `gcd(q1,N)=1`; and {class_summary['pair_product_fixed_divisor_one']}/293 pair products have fixed divisor one.  The producer constructs the exact composed polynomial and primitive normalization at `l17_schinzel_audit.py:118-142` and computes the pair-product gcd at `l17_schinzel_audit.py:148-158`.  Because the product has degree nine, gcd one already on values `t=0,…,9` is an exact finite-difference certificate, not statistical evidence; this replay recomputed those ten values.  The two polynomials are coprime in `Q[t]` on 293/293 rows (the degree-8 irreducible `G` does not vanish at the root of the linear `Q`).",
        "",
        "All leading coefficients are already positive on the canonical branches.  For a future class with negative leading coefficient one must apply Schinzel H to `-G`; this sign change does not affect roots, irreducibility, fixed divisors, valuation parity, or Hilbert symbols.",
        "",
        "## 5. Admissibility and density of the side conditions",
        "",
        "`N` is divisible by 8 and `q1,f` are odd, so every member has `v2(b)=0`; the inadmissible `v2(b)=1` wall boundary is never used (the exact boundary is recorded at `THEOREMS.md:1423-1426`).  Positivity, the real Cauchy bound, `b≠0,1`, and the factor-by-factor moving-prime unit conditions fail at only finitely many `t`.  Thus the deterministic eligibility set is cofinite and has natural density 1.  This is not a claim that prime values have positive density: their infinitude is exactly the Schinzel input.",
        "",
        "## 6. Prime-rung cross-check",
        "",
        f"**PROVED per sampled row.** The existing authority has 197 prime-rung closures; its cofactor summary records 197 singleton cofactor sets, 197 frozen/infinity-aligned rows, and 197 reciprocity-forced `+1` signs (`data/l17_cofactor.jsonl:294`).  This replay selected {n_sample} deterministic low-digit rows ({sample_summary['sample_family_counts'].get('ESC', 0)} ESC, {sample_summary['sample_family_counts'].get('L11', 0)} L11) and reran the proven-primality kernel, `smooth_emergent`, and `cofactor_decide`: {n_sample}/{n_sample} have a single proved-prime ladder residual, Hilbert sign `+1`, and verdict `zero`.",
        "",
        f"Strip convention matters: `cofactor_decide` receives the residual **after** `smooth_emergent` has additionally removed every prime `≤10^6`.  Only {raw_single}/{n_sample} sampled rows already have a literal single prime immediately after the raw `S∪supp(b)` strip; the others have extra small factors, all checked `+1`, before the singleton cofactor.  The 197/197 statement in `data/l15_remainders.jsonl:294` is about this post-smooth ladder residual.  By contrast, a Schinzel value of the normalized `G` is stronger: its raw S-stripped value itself is prime.",
        "The equality checked in the 24 records is the exact prime-rung mechanism, not a global biconditional: the 96 factorint closures also have verdict `zero`, so a singleton prime is sufficient but not necessary.",
        "",
        "## 7. Exact scope of the upgrade",
        "",
        "**PROVED:** for each fixed verified aligned class whose pair passes the standard polynomial hypotheses, classical Schinzel H implies the operational per-class H clause, indeed infinitely often.  **CONDITIONAL:** the actual members exist only under Schinzel H.  **Not supplied by Schinzel:** existence of an aligned class for every cell.  Therefore the analytic member-existence clause can be named Schinzel H rather than a bespoke conjecture, while any global all-cell statement still needs an independent class-construction theorem.",
        "",
        "Artifacts: `data/l18_schinzel_implies_h.jsonl` (24 member rows plus one summary) and generator `l18_schinzel_implies_h.py`.",
    ]
    REPORT.write_text("\n".join(lines) + "\n")


def main() -> None:
    closure_payload = json.loads(CLOSURES.read_text())
    records = closure_payload["records"]
    sch_lines = [json.loads(line) for line in SCHINZEL.read_text().splitlines()]
    sch_meta, sch_rows = sch_lines[0], sch_lines[1:294]
    assert sch_meta["canonical_records"] == 293
    remainder_rows = [
        row
        for row in (json.loads(line) for line in REMAINDERS.read_text().splitlines())
        if row.get("type") == "record"
    ]
    cofactor_summary = json.loads(COFACTOR.read_text().splitlines()[-1])
    assert cofactor_summary["counts"]["prime_rung_closures"] == 197
    assert cofactor_summary["parity"]["prime_rung_auto_plus1_count"] == 197

    class_summary, _ = fixed_class_audit(records, sch_rows)
    sample_rows, sample_summary = prime_rung_sample(remainder_rows, records)
    summary = {
        "type": "summary",
        "label": "PROVED_IMPLICATION_CONDITIONAL_ON_SCHINZEL_H",
        "verdict": "Schinzel H for an admissible fixed-class pair implies infinitely many operational ladder-zero members",
        "operational_definition": "zero means every outside odd-valuation place has Hilbert symbol +1; it does not mean there are no odd valuations",
        "class_audit": class_summary,
        "prime_rung_crosscheck": sample_summary,
        "cofactor_authority": {
            "prime_rung_rows": cofactor_summary["counts"]["prime_rung_closures"],
            "singleton_rows": cofactor_summary["parity"]["prime_rung_singleton_E_rows"],
            "auto_plus_one": cofactor_summary["parity"]["prime_rung_auto_plus1_count"],
        },
        "proof_dependencies": {
            "parity_law": "THEOREMS.md:828-853",
            "moving_symbol": "THEOREMS.md:887-903 and l12_class.py:21-28",
            "frozen_constancy": "h10q.py:2601-2608,2617-2688",
            "ladder": "l13_filter.py:90-144,281-366",
        },
        "scope": "per fixed verified aligned class; class existence is separate",
    }
    with OUT.open("w") as fh:
        for row in sample_rows:
            fh.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        fh.write(json.dumps(summary, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(class_summary, sample_summary, cofactor_summary)
    print(json.dumps(summary, sort_keys=True))
    print(f"wrote {OUT} ({len(sample_rows)} verified rows + summary)")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
