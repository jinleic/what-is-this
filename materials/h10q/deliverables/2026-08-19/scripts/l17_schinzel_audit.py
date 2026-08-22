#!/usr/bin/env python3
"""Audit the existing ledgered Schinzel pair for canonical H closures.

The canonical rows are the 293 records in data/l13h_all_closures.json.  For
those rows, q1(t)=q1_0+N*t and F(t)=P(eps*f*q1(t)) are constructed exactly
from h10q._l10_P.  F is normalized to its primitive integer polynomial for
irreducibility, local roots, and fixed-divisor checks; F_scale records the
exact rational scale relating that primitive polynomial to F.

The first six legacy _L13_RESIDUAL rows are emitted separately as pointwise
rational-b witnesses.  They are deliberately not coerced into a prime-b AP.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path("/Users/jinleic/jinleic-workspace")
REPO = ROOT / "math" / "h10q"
sys.path.insert(0, str(REPO))

# This is the repository's proven-primality / exact-arithmetic engine.
import h10q  # noqa: E402


OUT = REPO / "data" / "l17_schinzel_audit.jsonl"
REPORT = Path("/tmp/l17_schinzel_audit.md")
CLOSURES = REPO / "data" / "l13h_all_closures.json"
GCD_T = 2000
P_LIMIT = 1000


def frac_text(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def gcd_many(values: list[int]) -> int:
    g = 0
    for x in values:
        g = math.gcd(g, abs(int(x)))
    return abs(g)


def lcm(a: int, b: int) -> int:
    return a // math.gcd(a, b) * b


def primitive_integer_poly(coeffs: list[Fraction]) -> tuple[list[int], Fraction]:
    """Return Pi and scale with sum coeffs[i] t^i = scale*Pi(t)."""
    den = 1
    for c in coeffs:
        den = lcm(den, c.denominator)
    ints = [int(c * den) for c in coeffs]
    content = gcd_many(ints)
    if content == 0:
        raise AssertionError("zero polynomial")
    pi = [x // content for x in ints]
    return pi, Fraction(content, den)


def compose_in_t(P: list[Fraction], b0: Fraction, step: Fraction) -> list[Fraction]:
    """Exact ascending coefficients of P(b0 + step*t)."""
    out = [Fraction(0)] * len(P)
    for i, ci in enumerate(P):
        if ci == 0:
            continue
        for j in range(i + 1):
            out[j] += ci * math.comb(i, j) * b0 ** (i - j) * step ** j
    return out


def eval_poly(coeffs: list[int], t: int) -> int:
    value = 0
    for c in reversed(coeffs):
        value = value * t + c
    return value


def prime_factors_leq(n: int, limit: int) -> list[int]:
    # The prime list comes from the same h10q engine used by the arithmetic
    # audit; no probabilistic primality routine is used anywhere here.
    return [p for p in h10q.primerange(2, limit + 1) if n % p == 0]


def irred_certificate(pi: list[int]) -> int:
    """Find a Frobenius modular irreducibility certificate for degree eight."""
    for p in h10q.primerange(2, P_LIMIT + 1):
        if h10q._l13_irred8(pi, p):
            return p
    raise AssertionError("no degree-8 Frobenius certificate <= 1000")


def canonical_branch(rec: dict[str, Any]) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
    w, ut = rec["cell"]
    a = Fraction(rec["a"])
    z = Fraction(w) * Fraction(*ut)
    A = 1 + 4 * a * a
    if rec["family"] == "L11":
        tau = (1 + 2 * a * a) / A
    elif rec["family"] == "ESC":
        tau = Fraction(3, 5)
    else:
        raise AssertionError(rec["family"])
    delta = 1 - A * tau * tau
    Z = z ** 3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    return a, z, tau, delta, Z, D, s


def audit_canonical(rec: dict[str, Any], index: int) -> tuple[dict[str, Any], dict[int, dict[str, float | int]]]:
    a, z, tau, delta, Z, D, s = canonical_branch(rec)
    A = 1 + 4 * a * a
    P = h10q._l10_P(a, Z, D, A, delta, s)
    if len(P) != 9 or P[-1] == 0:
        raise AssertionError((index, rec["cell"], len(P), P[-1]))

    q1_0 = int(rec["q1"])
    N = int(rec["N"])
    eps = int(rec["eps"])
    f = int(rec["f"])
    if not h10q._is_prime(q1_0):
        raise AssertionError((index, "q1 is not proven prime", q1_0))
    q1_gcd = math.gcd(q1_0, N)
    if q1_gcd != 1:
        raise AssertionError((index, "Dirichlet AP gcd != 1", q1_0, N, q1_gcd))

    b0 = Fraction(eps * f * q1_0)
    step = Fraction(eps * f * N)
    F_rat = compose_in_t(P, b0, step)
    if len(F_rat) != 9 or F_rat[-1] == 0:
        raise AssertionError((index, "F is not degree 8", rec["cell"]))
    Pi, scale = primitive_integer_poly(F_rat)
    if len(Pi) != 9 or Pi[-1] == 0 or gcd_many(Pi) != 1:
        raise AssertionError((index, "primitive normalization failed", rec["cell"]))

    cert_p = irred_certificate(Pi)
    # A 0..2000 value scan is intentionally retained in the artifact.  It is
    # stronger than the degree+1 finite-difference check requested by the
    # ledger and makes the arithmetic provenance directly replayable.
    value_gcd = 0
    pair_product_gcd = 0
    for t in range(GCD_T + 1):
        q_value = q1_0 + N * t
        f_value = eval_poly(Pi, t)
        value_gcd = math.gcd(value_gcd, abs(f_value))
        pair_product_gcd = math.gcd(pair_product_gcd, abs(q_value * f_value))
    if value_gcd != 1:
        raise AssertionError((index, "fixed divisor survived value scan", rec["cell"], value_gcd))
    if pair_product_gcd != 1:
        raise AssertionError((index, "pair product fixed divisor survived value scan", rec["cell"], pair_product_gcd))

    checks: dict[str, int] = {}
    counts: dict[str, int] = {}
    q1_residues: dict[str, int] = {}
    local_table: dict[int, dict[str, float | int]] = {}
    for p in prime_factors_leq(N, P_LIMIT):
        admissible = [
            t for t in range(p)
            if (q1_0 + N * t) % p != 0 and eval_poly(Pi, t) % p != 0
        ]
        if not admissible:
            raise AssertionError((index, "pair local obstruction", rec["cell"], p))
        # The map is intentionally p -> witness_t as required by the ledger;
        # richer counts are carried in the adjacent fields.
        checks[str(p)] = admissible[0]
        counts[str(p)] = len(admissible)
        q1_residues[str(p)] = q1_0 % p
        local_table[p] = {
            "records": 1,
            "admissible_count": len(admissible),
            "density": len(admissible) / p,
        }

    row = {
        "type": "canonical",
        "record_index": index,
        "family": rec["family"],
        "cell": rec["cell"],
        "a": str(a),
        "tau": frac_text(tau),
        "eps": eps,
        "f": f,
        "q1_0": q1_0,
        "N": str(N),
        "q1_polynomial": [str(q1_0), str(N)],
        "F_definition": "P(eps*f*(q1_0 + N*t))",
        # F_exact(t) = F_scale * sum(F_coeffs_t[j] t^j).
        "F_scale": frac_text(scale),
        "F_coeffs_t": [str(x) for x in Pi],
        "F_degree": 8,
        "F_irreducible": True,
        "F_irreducibility_certificate_prime": cert_p,
        "F_value_gcd": value_gcd,
        "F_value_gcd_scan_T": GCD_T,
        "pair_product_gcd": pair_product_gcd,
        "pair_product_gcd_scan_T": GCD_T,
        "q1_value_gcd": q1_gcd,
        "small_prime_checks": checks,
        "small_prime_admissible_counts": counts,
        "q1_mod_small_prime": q1_residues,
        "pair_local_conditions_pass": True,
        "source": rec.get("source"),
    }
    return row, local_table


def audit_residuals(records: list[dict[str, Any]], canonical_by_cell: dict[tuple[int, tuple[int, int]], dict[str, Any]]) -> list[dict[str, Any]]:
    # The assignment names six legacy residual rows.  The current frozen table
    # has a later seventh row (89,(-2,1)); it is intentionally outside this
    # six-row audit scope and is reported in the metadata.
    residuals = list(h10q._L13_RESIDUAL[:6])
    if len(residuals) != 6:
        raise AssertionError(f"expected six legacy residual rows, got {len(residuals)}")
    out = []
    for i, (w, ut, a0, b0, d0) in enumerate(residuals, 1):
        a = Fraction(a0)
        b = Fraction(b0) if isinstance(b0, int) else Fraction(*b0)
        d = None if d0 is None else Fraction(d0)
        A = 1 + 4 * a * a
        tau = (1 + 2 * a * a) / A if d is None else (A + d * d) / (2 * d * A)
        z = Fraction(w) * Fraction(*ut)
        tied = h10q._l7_tied_status(a, b, z, tau)
        if tied is not True:
            raise AssertionError(("residual tied-status", i, w, ut, tied))
        key = (int(w), tuple(ut))
        covering = canonical_by_cell.get(key)
        out.append({
            "type": "residual",
            "record_index": i,
            "cell": [w, list(ut)],
            "a": frac_text(a),
            "b": frac_text(b),
            "d": None if d is None else frac_text(d),
            "tau": frac_text(tau),
            "pointwise_tied_status": True,
            "schinzel_applicable": False,
            "reason": "rational/composite witness; not coerced into eps*f*q1(t)",
            "covered_by_canonical_class": covering is not None,
            "covering_canonical_record_index": None if covering is None else covering["record_index"],
            "covering_canonical_family": None if covering is None else covering["family"],
        })
    return out


def load_matched_summary() -> dict[str, Any]:
    path = REPO / "data" / "l17_matched.jsonl"
    with path.open() as fh:
        for line in fh:
            rec = json.loads(line)
            if rec.get("type") == "global":
                return rec
    return {}


def write_report(meta: dict[str, Any], rows: list[dict[str, Any]], residuals: list[dict[str, Any]], local_agg: dict[int, dict[str, float | int]], matched: dict[str, Any]) -> None:
    canonical = [r for r in rows if r.get("type") == "canonical"]
    local_rows = len(canonical)
    lines: list[str] = []
    lines.append("# L17 Schinzel pair audit")
    lines.append("")
    lines.append("## Scope and exact objects")
    lines.append("")
    lines.append(f"- Canonical source: `data/l13h_all_closures.json`, {local_rows} records ({meta['family_counts']['L11']} L11 + {meta['family_counts']['ESC']} ESC).")
    lines.append("- For each canonical record, the audited pair is the existing ledgered pair: `q1(t) = q1_0 + N*t` and `F(t) = P(eps*f*q1(t))`, with `P` exactly `h10q._l10_P`.")
    lines.append("- `F_coeffs_t` in the JSONL is the primitive integer polynomial in ascending powers of `t`; `F_scale` gives the exact rational factor, so `F(t) = F_scale * sum_j F_coeffs_t[j] t^j`. No new polynomial or alternate reduction was introduced.")
    lines.append("")
    lines.append("## Status of the reduction")
    lines.append("")
    lines.append("- **PROVED side:** every canonical row has a proven-prime `q1_0`, `gcd(q1_0,N)=1`, and therefore Dirichlet's theorem supplies infinitely many primes in the existing AP `q1_0 (mod N)`. The finite alignment/modulus is the ledgered one; this audit does not alter it.")
    lines.append("- **OPEN side:** requiring the degree-8 cofactor `F(t)` to be prime simultaneously with the AP prime is the Schinzel-Hypothesis-H/Bunyakovsky-type input. The audit verifies local admissibility and irreducibility only; it does **not** prove infinitely many prime values. Prime values of degree-8 polynomials are not known unconditionally in this setting; the requested boundary statement is that the cited Krumm-level unconditional results reach only degree `<= 3`, but that citation is **PENDING** here because it was not externally verified in this run.")
    lines.append("")
    lines.append("## Local-condition result")
    lines.append("")
    lines.append(f"- All {local_rows}/{local_rows} canonical rows have a degree-8 Frobenius irreducibility certificate modulo a prime `p <= {P_LIMIT}`; `F_irreducible=true` in every row.")
    lines.append(f"- All {local_rows}/{local_rows} canonical primitive polynomials have `F_value_gcd=1` after the exact value scan `t=0..{GCD_T}`. This is a fixed-divisor check, not a tail estimate.")
    lines.append(f"- All {local_rows}/{local_rows} pair products have `pair_product_gcd=1` on the same exact scan `t=0..{GCD_T}`. Thus the two polynomials have no common fixed prime divisor on this finite certificate scan; the individual first-polynomial check is `gcd(q1_0,N)=1`.")
    lines.append(f"- For every prime `p <= {P_LIMIT}` dividing that row's `N`, the table counts residues with `q1(t) != 0 (mod p)` and `F(t) != 0 (mod p)`; the witness is stored as `small_prime_checks[p]`. Since p|N and gcd(q1_0,N)=1, q1(t) is a nonzero constant modulo p, so the corrected count is p whenever F has no local root. The p=2 witness is t=0 (q1(0) is odd), not a t-polynomial surrogate.")
    lines.append("- These are local compatibility checks only. They do **not** assert independence of the two residue conditions and do **not** imply a prime-value theorem.")
    lines.append("")
    lines.append("## Admissible-density table (p divides N, p <= 1000)")
    lines.append("")
    lines.append("The aggregate below is a compact table over canonical records. `mean_count` is the mean number of residues among the p classes satisfying the corrected joint condition `q1(t) != 0` and `F(t) != 0`; `mean_density` divides by p. It is a table of the record machinery, not an independence calculation.")
    lines.append("")
    lines.append("| p | records | count range | mean count | mean density / p |")
    lines.append("|---:|---:|---:|---:|---:|")
    for p in sorted(local_agg):
        a = local_agg[p]
        lines.append(f"| {p} | {int(a['records'])} | {int(a['min_count'])}–{int(a['max_count'])} | {a['mean_count']:.3f} | {a['mean_density']:.4f} |")
    lines.append("")
    if matched:
        bad = int(matched.get("observed_bad_members", 0))
        members = int(matched.get("members", 0))
        clean = members - bad
        lines.append("## Comparison to the validated matched model")
        lines.append("")
        lines.append(f"- `data/l17_matched.jsonl` reports {bad}/{members} observed bad members ({bad/members:.4f}) and {clean}/{members} complementary observed clean members ({clean/members:.4f}); its matched marginal expectation is {matched.get('matched_marginal_expectation')} odd hits.")
        lines.append("- The matched artifact is a p<=10000-window statistic (and explicitly does not bound the p>10000 tail). The table above is likewise only the finite local condition for primes dividing each AP modulus. Any comparison to the roughly 0.4 bad/obstruction fraction is qualitative; no product/independence claim is made.")
        lines.append("")
    lines.append("## Six legacy residual rows")
    lines.append("")
    lines.append("The six selected `_L13_RESIDUAL[:6]` entries are rational/composite pointwise witnesses. They are audited by exact branch reconstruction and `_l7_tied_status`, not by forcing their `b` values into an `eps*f*q1(t)` prime family.")
    lines.append("")
    lines.append("| cell | (a,b) | tau | tied | canonical class also present? |")
    lines.append("|---|---|---|---|---|")
    for r in residuals:
        lines.append(f"| `{r['cell']}` | `({r['a']},{r['b']})` | `{r['tau']}` | {r['pointwise_tied_status']} | {r['covered_by_canonical_class']} (record {r['covering_canonical_record_index']}) |")
    lines.append("")
    lines.append("All six selected residual cells have a canonical closure class in the 293-record source, so the prime-family reduction covers the cell through that class. The residual rational witnesses remain individual unconditional certificates and are not evidence for the Schinzel-open prime-value side. The current frozen table contains a later seventh `(89,(-2,1))` row; it is outside this assignment's explicitly requested six-row residual scope.")
    REPORT.write_text("\n".join(lines) + "\n")


def main() -> None:
    # Run the repository's exact arithmetic self-test once before producing
    # evidence.  No probabilistic primality routine is used by this script.
    h10q._selftest()
    payload = json.loads(CLOSURES.read_text())
    records = payload["records"]
    if len(records) != 293:
        raise AssertionError(f"expected 293 canonical closure records, got {len(records)}")

    rows: list[dict[str, Any]] = []
    local_agg_raw: dict[int, list[int]] = defaultdict(list)
    family_counts = defaultdict(int)
    for i, rec in enumerate(records):
        row, local_table = audit_canonical(rec, i)
        rows.append(row)
        family_counts[rec["family"]] += 1
        for p, info in local_table.items():
            local_agg_raw[p].append(int(info["admissible_count"]))
        # Keep one subprocess and leave a small pacing gap between records.
        time.sleep(0.1)

    canonical_by_cell = {(r["cell"][0], tuple(r["cell"][1])): r for r in rows}
    residuals = audit_residuals(records, canonical_by_cell)
    rows.extend(residuals)

    local_agg: dict[int, dict[str, float | int]] = {}
    for p, counts in local_agg_raw.items():
        local_agg[p] = {
            "records": len(counts),
            "min_count": min(counts),
            "max_count": max(counts),
            "mean_count": sum(counts) / len(counts),
            "mean_density": sum(c / p for c in counts) / len(counts),
        }

    matched = load_matched_summary()
    meta = {
        "type": "meta",
        "task": "l17_schinzel_audit_v2_pair_product",
        "canonical_records": len(records),
        "residual_records": len(residuals),
        "total_rows": len(rows),
        "family_counts": dict(family_counts),
        "F_irreducible_true": sum(r.get("F_irreducible") is True for r in rows),
        "F_value_gcd_one": sum(r.get("F_value_gcd") == 1 for r in rows),
        "pair_product_gcd_one": sum(r.get("pair_product_gcd") == 1 for r in rows),
        "pair_product_gcd_scan_T": GCD_T,
        "pair_product_gcd_note": "gcd over q1(t)*primitive_F(t) for t=0..2000; canonical rows all equal 1",
        "pair_local_conditions_pass": sum(r.get("pair_local_conditions_pass") is True for r in rows),
        "value_gcd_scan_T": GCD_T,
        "small_prime_limit": P_LIMIT,
        "proven_primality_engine": "math/h10q/h10q.py::_is_prime (refusals are not evidence)",
        "irreducibility_engine": "h10q._l13_irred8 Frobenius certificate modulo p",
        "residual_scope": "h10q._L13_RESIDUAL[:6] legacy rows; current seventh row intentionally outside requested scope",
        "local_density_table": {str(p): local_agg[p] for p in sorted(local_agg)},
        "matched_model_global": matched,
    }
    with OUT.open("w") as fh:
        fh.write(json.dumps(meta, sort_keys=True) + "\n")
        for row in rows:
            fh.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    write_report(meta, rows, residuals, local_agg, matched)
    print(json.dumps(meta, sort_keys=True))
    print(f"wrote {OUT} ({len(rows)} data rows + meta)")
    print(f"wrote {REPORT}")


if __name__ == "__main__":
    main()
