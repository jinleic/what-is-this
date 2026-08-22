#!/usr/bin/env python3
"""Classify the canonical frozen-place-5 ESC wall for z=-w, 101<=w<=400.

The universal statement in this file is symbolic.  The bounded q1<=200 scan
only audits that identity and finds certified L11/ESC base classes; a bounded
absence is always labelled EVIDENCE, and kernel refusals are never evidence.
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


def find_root() -> Path:
    candidates = []
    if os.environ.get("H10Q_ROOT"):
        candidates.append(Path(os.environ["H10Q_ROOT"]).expanduser())
    here = Path(__file__).resolve().parent
    cwd = Path.cwd()
    candidates.extend((cwd, cwd / "math" / "h10q", here))
    for candidate in candidates:
        if (candidate / "h10q.py").is_file():
            return candidate.resolve()
    raise RuntimeError("set H10Q_ROOT to the math/h10q directory")


ROOT = find_root()
sys.path.insert(0, str(ROOT))

import h10q  # noqa: E402  (the proven kernel is the authority)
import l12_class  # noqa: E402

W_MIN = 101
W_MAX = 400
Q1_MAX = 200
UT = (-1, 1)
OUT = ROOT / "data" / "l18_fivewall.jsonl"
REPORT = Path(os.environ.get("L18_FIVEWALL_REPORT", "/tmp/l18_fivewall.md"))


def refusal(exc: BaseException) -> str:
    return "refused:" + type(exc).__name__


def place_name(place) -> str:
    return str(place)


def symbol_map(syms: dict) -> dict[str, int]:
    return {place_name(place): int(value) for place, value in syms.items()}


def cert_summary(cert: dict | None) -> dict:
    if cert is None:
        return {}
    return {
        "ok": bool(cert["ok"]),
        "N": str(cert["N"]),
        "S": list(cert["S"]),
        "ks": {str(p): int(k) for p, k in cert["ks"].items()},
        "syms": symbol_map(cert["syms"]),
        "Q0": str(cert["Q0"]),
        "excluded": list(cert.get("excluded", [])),
    }


def chi5(value: int) -> int:
    """Quadratic character modulo 5, for a 5-adic unit."""
    assert value % 5
    return h10q.legendre(F(value), 5)


def fresh_l12_cert(w: int, eps: int, q1: int):
    """Call the canonical ESC certificate with one temporary off-grid row."""
    key = (w, UT)
    sentinel = object()
    previous = l12_class._L12_ESCAPE.get(key, sentinel)
    l12_class._L12_ESCAPE[key] = (w, eps, q1)
    try:
        return l12_class.l12_class_cert(w, *UT)
    finally:
        if previous is sentinel:
            l12_class._L12_ESCAPE.pop(key, None)
        else:
            l12_class._L12_ESCAPE[key] = previous


def negative_profile(cert: dict) -> tuple[str, ...]:
    return tuple(
        sorted(
            (place_name(place) for place, value in cert["syms"].items() if value == -1),
            key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value),
        )
    )


def l11_candidates(w: int) -> tuple[list[int], list[int]]:
    roots = [r for r in range(1, w) if (1 + 4 * r * r) % w == 0]
    a_pool = sorted(
        {
            abs(r + sign * w * lift)
            for r in roots
            for sign in (1, -1)
            for lift in range(3)
            if r + sign * w * lift != 0
        }
    )
    return roots, [a for a in a_pool if a % 2 == 1]


def scan_l11(w: int, q1_primes: list[int]) -> dict:
    """Find one proven aligned L11 class, if q1<=Q1_MAX supplies one."""
    roots, a_pool = l11_candidates(w)
    minus_one = h10q.legendre(F(-1), w)
    assert bool(roots) == (minus_one == 1)
    status_counts: Counter[str] = Counter()
    refusals: Counter[str] = Counter()
    witness = None

    for a_int in a_pool:
        a = F(a_int)
        A = 1 + 4 * a * a
        tau = (1 + 2 * a * a) / A
        for eps in (1, -1):
            for q1 in q1_primes:
                if q1 in {2, 3, 5, 7}:
                    continue
                time.sleep(0.1)
                try:
                    cert = h10q._l10_class_cert(a, F(-w), tau, eps, q1)
                except Exception as exc:  # proven-engine refusals are recorded only
                    label = refusal(exc)
                    refusals[label] += 1
                    status_counts[label] += 1
                    continue
                if cert is None:
                    status_counts["refused:none"] += 1
                    refusals["refused:none"] += 1
                    continue
                if not cert["ok"]:
                    status_counts["cert-not-ok"] += 1
                    continue
                status_counts["aligned"] += 1
                witness = {
                    "a": a_int,
                    "eps": eps,
                    "q1": q1,
                    "tau": str(tau),
                    "cert": cert_summary(cert),
                }
                break
            if witness is not None:
                break
        if witness is not None:
            break

    return {
        "minus_one_legendre_w": minus_one,
        "roots": roots,
        "odd_a_pool": a_pool,
        "root_criterion": "available" if roots else "unavailable",
        "root_criterion_label": "PROVED",
        "bounded_q1_max": Q1_MAX,
        "attempt_status_counts": dict(sorted(status_counts.items())),
        "refusal_counts": dict(sorted(refusals.items())),
        "witness": witness,
    }


def blank_branch(q1_values: list[int]) -> dict:
    return {
        "q1_values": q1_values,
        "q1_count": len(q1_values),
        "attempt_count": 0,
        "verified_certificate_count": 0,
        "aligned_count": 0,
        "hilbert5_wild_pairs": Counter(),
        "negative_profiles": Counter(),
        "refusal_counts": Counter(),
    }


def scan_esc(w: int, q1_primes: list[int]) -> dict:
    """Exhaust canonical f=w, eps=+-1 base classes through q1<=Q1_MAX."""
    controlled = {2, 3, 5, 7, w}
    eligible = [q1 for q1 in q1_primes if q1 not in controlled]
    branches = {
        str(character): blank_branch([q1 for q1 in eligible if chi5(q1) == character])
        for character in (1, -1)
    }
    ineligible = Counter()
    refusals: list[dict] = []
    mismatches: list[dict] = []
    first_aligned = None
    aligned_total = 0
    z = F(-w)
    Z = z**3
    D = 1 - Z - Z * Z
    k_D = h10q.vp(D, 5)
    assert h10q.vp(Z, 5) == 0 and k_D >= 0
    w_character = chi5(w)

    for eps in (1, -1):
        for q1 in q1_primes:
            time.sleep(0.1)
            try:
                cert = fresh_l12_cert(w, eps, q1)
            except Exception as exc:  # proven-engine refusals are labels, never evidence
                label = refusal(exc)
                record = {"eps": eps, "q1": q1, "label": label}
                refusals.append(record)
                if q1 not in controlled:
                    branch = branches[str(chi5(q1))]
                    branch["attempt_count"] += 1
                    branch["refusal_counts"][label] += 1
                else:
                    ineligible[label] += 1
                continue

            if q1 in controlled:
                status = "ineligible:q1-in-controlled-set"
                ineligible[status] += 1
                if cert is not None:
                    mismatches.append(
                        {"eps": eps, "q1": q1, "kind": "controlled-q1-not-refused"}
                    )
                continue

            character = chi5(q1)
            branch = branches[str(character)]
            branch["attempt_count"] += 1
            if cert is None:
                label = "refused:none"
                branch["refusal_counts"][label] += 1
                refusals.append({"eps": eps, "q1": q1, "label": label})
                continue

            branch["verified_certificate_count"] += 1
            b = F(eps * w * q1)
            c0 = h10q._sun_h(F(1), b, Z)
            assert c0 is not None
            delta = F(-4, 5)
            x0 = F(4) * (16 - delta * c0 * c0)
            d0 = F(8) * b
            expected_v5_c = -1 - k_D
            expected_v5_x = -3 - 2 * k_D
            hilbert_5 = h10q.hilbert(x0, d0, 5)
            wild_q1 = h10q.legendre(F(5), q1)
            predicted_hilbert_5 = -w_character * character
            checks = {
                "v5_c": h10q.vp(c0, 5) == expected_v5_c,
                "v5_x": h10q.vp(x0, 5) == expected_v5_x,
                "v5_d": h10q.vp(d0, 5) == 0,
                "cert_hilbert_5": cert["syms"][5] == hilbert_5,
                "cert_wild_q1": cert["syms"]["q1"] == wild_q1,
                "hilbert_formula": hilbert_5 == predicted_hilbert_5,
                "reciprocity": wild_q1 == character,
                "product_formula": hilbert_5 * wild_q1 == -w_character,
            }
            failed = [name for name, ok in checks.items() if not ok]
            if failed:
                mismatches.append(
                    {
                        "eps": eps,
                        "q1": q1,
                        "failed_checks": failed,
                        "actual": {"hilbert_5": hilbert_5, "wild_q1": wild_q1},
                        "predicted": {
                            "hilbert_5": predicted_hilbert_5,
                            "product": -w_character,
                        },
                    }
                )

            pair = f"{hilbert_5:+d}/{wild_q1:+d}"
            branch["hilbert5_wild_pairs"][pair] += 1
            profile = negative_profile(cert)
            branch["negative_profiles"][",".join(profile) if profile else "(none)"] += 1
            if cert["ok"]:
                branch["aligned_count"] += 1
                aligned_total += 1
                if first_aligned is None:
                    first_aligned = {
                        "eps": eps,
                        "q1": q1,
                        "q1_character_mod5": character,
                        "cert": cert_summary(cert),
                    }

    for branch in branches.values():
        branch["hilbert5_wild_pairs"] = dict(sorted(branch["hilbert5_wild_pairs"].items()))
        branch["negative_profiles"] = dict(sorted(branch["negative_profiles"].items()))
        branch["refusal_counts"] = dict(sorted(branch["refusal_counts"].items()))

    predicted_wall = w_character == 1
    prediction = (
        "5-obstructed-for-every-admissible-q1"
        if predicted_wall
        else "5-compatible-exactly-when-(q1|5)=+1"
    )
    return {
        "f": w,
        "eps": [1, -1],
        "a": "1",
        "tau": "3/5",
        "bounded_q1_max": Q1_MAX,
        "controlled_set": sorted(controlled),
        "ineligible_status_counts": dict(sorted(ineligible.items())),
        "refusals": refusals,
        "v5_D": k_D,
        "w_character_mod5": w_character,
        "criterion_prediction": prediction,
        "criterion_label": "PROVED",
        "branches_by_q1_character_mod5": branches,
        "aligned_count": aligned_total,
        "first_aligned": first_aligned,
        "identity_mismatches": mismatches,
    }


def common_escape_branch_negatives(esc: dict) -> list[str]:
    """Common negatives among verified (q1|5)=+1 non-aligned attempts."""
    profiles = esc["branches_by_q1_character_mod5"]["1"]["negative_profiles"]
    expanded: list[set[str]] = []
    for profile, count in profiles.items():
        if profile == "(none)":
            expanded.extend([set()] * count)
        else:
            expanded.extend([set(profile.split(","))] * count)
    if not expanded:
        return []
    return sorted(set.intersection(*expanded), key=lambda x: (not x.isdigit(), x))


def classify(w: int, l11: dict, esc: dict) -> tuple[str, str, str]:
    if l11["witness"] is not None:
        return "L11 route available", "PROVED", "an aligned L11 class is certified"
    if esc["first_aligned"] is not None:
        return "ESC route available", "PROVED", "an aligned canonical ESC class is certified"
    if esc["w_character_mod5"] == 1:
        return (
            "5-obstructed",
            "PROVED",
            "Hilb_5 * wild_q1 = -1 for every admissible canonical ESC q1",
        )
    common = common_escape_branch_negatives(esc)
    detail = (
        "bounded q1<=200 scan: common negative place(s) on the 5-compatible branch: "
        + ",".join(common)
        if common
        else "bounded q1<=200 scan found no aligned class; no universal claim"
    )
    return "other-obstructed", "EVIDENCE", detail


def per_w_record(w: int, q1_primes: list[int]) -> dict:
    time.sleep(0.1)
    assert h10q._is_prime(w)
    l11 = scan_l11(w, q1_primes)
    esc = scan_esc(w, q1_primes)
    category, label, reason = classify(w, l11, esc)
    full_wall_prediction = (
        l11["minus_one_legendre_w"] == -1 and esc["w_character_mod5"] == 1
    )
    assert full_wall_prediction == (w % 20 in (11, 19))
    return {
        "type": "per-w",
        "label": label,
        "w": w,
        "cell": [w, [-1, 1]],
        "z": str(-w),
        "w_mod_20": w % 20,
        "criterion": {
            "esc_fivewall": esc["criterion_prediction"],
            "esc_fivewall_label": "PROVED",
            "w131_type_two_route_wall": full_wall_prediction,
            "w131_type_congruence": "w mod 20 in {11,19}",
            "w131_type_label": "PROVED",
        },
        "l11": l11,
        "esc": esc,
        "classification": category,
        "classification_label": label,
        "classification_reason": reason,
        "classification_scope": (
            "route availability is a positive certificate with q1<=200; "
            "5-obstructed is universal; other-obstructed is bounded evidence only"
        ),
    }


def branch_cell(record: dict, character: int) -> str:
    branch = record["esc"]["branches_by_q1_character_mod5"][str(character)]
    return (
        f"{branch['q1_count']}q/{branch['verified_certificate_count']}v/"
        f"{branch['aligned_count']}a"
    )


def obstruction_note(record: dict) -> str:
    if record["classification"] != "other-obstructed":
        return "-"
    common = common_escape_branch_negatives(record["esc"])
    if common:
        return "common -1: " + ",".join(common)
    refusals = len(record["esc"]["refusals"])
    return f"bounded no-hit; refusals={refusals}"


def write_jsonl(records: list[dict], summary: dict) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "type": "meta",
        "label": "EVIDENCE",
        "schema": "l18-fivewall-v1",
        "source_script": "math/h10q/l18_fivewall.py",
        "engine": "math/h10q/h10q.py proven-primality only",
        "w_range": [W_MIN, W_MAX],
        "q1_bound": Q1_MAX,
        "cell_shape": "[w,[-1,1]] (z=-w)",
        "canonical_ESC": {"a": "1", "tau": "3/5", "f": "w", "eps": [1, -1]},
        "scope_note": "canonical L11 root route plus canonical ESC a=1,tau=3/5,f=w only",
        "outside_protocol_counterexample": {
            "label": "PROVED",
            "cell": [131, [-1, 1]],
            "a": 7,
            "tau": "0",
            "f": 131,
            "eps": 1,
            "q1": 41,
            "result": "verified closure; the canonical 5-wall is not intrinsic",
        },
        "refusals_are_not_evidence": True,
    }
    with OUT.open("w", encoding="utf-8") as handle:
        for item in (meta, *records, summary):
            handle.write(json.dumps(item) + "\n")


def write_report(records: list[dict], summary: dict) -> None:
    lines = [
        "# L18 frozen-place-5 wall for the canonical ESC route",
        "",
        "## PROVED theorem",
        "",
        "Let `w != 5` be an odd prime, `z=-w`, `f=w`, `a=1`,",
        "`tau=3/5`, `eps in {+1,-1}`, and `q1` an admissible canonical ESC",
        "prime.  Put `b=eps*w*q1`, `Z=z^3`, and `D=1-Z-Z^2`.  Then",
        "",
        "    Hilb_5 = -(w|5)(q1|5),        wild_q1 = (5|q1) = (q1|5),",
        "    Hilb_5 * wild_q1 = -(w|5).",
        "",
        "Consequently the canonical ESC route is 5-obstructed for every",
        "admissible `q1` iff `(w|5)=+1`, equivalently `w == 1 or 4 (mod 5)`.",
        "If `(w|5)=-1`, place 5 and the wild condition are simultaneously +1",
        "exactly when `(q1|5)=+1`; the other residue branch has both symbols -1.",
        "",
        "| `(w|5)` | `(q1|5)` | `(Hilb_5, wild_q1)` | canonical consequence |",
        "|---:|---:|:---:|:---|",
        "| +1 | +1 | `(-1,+1)` | reciprocity wall |",
        "| +1 | -1 | `(+1,-1)` | reciprocity wall |",
        "| -1 | +1 | `(+1,+1)` | 5-compatible |",
        "| -1 | -1 | `(-1,-1)` | rejected at both places |",
        "",
        "By the proved L11 root construction, that route is available iff",
        "`1+4r^2 == 0 (mod w)` has a root, equivalently `(-1|w)=+1` or",
        "`w == 1 (mod 4)`.  Hence the w=131-type pair of canonical obstructions",
        "(no L11 root and a universal canonical-ESC 5-wall) occurs exactly for",
        "primes `w == 11 or 19 (mod 20)`.",
        "",
        "### Proof",
        "",
        "For `t=(b-1)^2/b`, the canonical Sun quantities are",
        "",
        "    g=(16-5*t^2)/5,  c=Z^2*g/D,  delta=-4/5,",
        "    x=4*(16-delta*c^2),  d=8*b.",
        "",
        "All of `b`, `w`, `q1`, and `Z` are 5-adic units.  The numerator of `g`",
        "is `16` modulo 5, so `v5(g)=-1`.  Write `k=v5(D)>=0`.  It follows that",
        "`v5(c)=-1-k`, `v5(x)=-3-2k` (odd), and `v5(d)=0`; the negative-valuation",
        "term cannot cancel against 16.  The odd-prime Hilbert formula with a unit",
        "second entry therefore gives",
        "",
        "    (x,d)_5 = (d|5) = (8*eps*w*q1|5) = -(w|5)(q1|5),",
        "",
        "because `(8|5)=-1` and both `+1` and `-1` are squares modulo 5.  The",
        "canonical wild class symbol is `(A|q1)` with `A=5`; quadratic reciprocity",
        "gives `(5|q1)=(q1|5)`.  Multiplication proves the criterion.  Finally,",
        "`1+4r^2 == 0 (mod w)` has a root exactly when `(-1|w)=+1`; combining the",
        "mod-5 and mod-4 conditions by CRT gives residues 11 and 19 modulo 20.",
        "",
        "**Scope (PROVED).**  This is a wall of the frozen canonical protocol, not",
        "an intrinsic obstruction at the cell.  The separately verified Route131",
        "construction closes `w=131` outside this protocol with `a=7`, `tau=0`,",
        "`f=131`, `eps=+1`, and `q1=41`; changing `a,tau` changes the frozen-symbol",
        "identity.  The congruence theorem above applies only to the L11 root route",
        "and the canonical `a=1,tau=3/5,f=w` ESC route.",
        "",
        "## EVIDENCE bounded discriminating experiment",
        "",
        f"The proven kernel scanned every prime `w` in [{W_MIN},{W_MAX}], both signs",
        f"of `eps`, and every proven prime `q1 <= {Q1_MAX}`.  Controlled-set q1",
        "values were called through the certificate and labelled ineligible; exceptions",
        "and `None` on eligible inputs are labelled refusals and never counted as",
        "evidence.  The two `(q1|5)` branches are counted separately below.  A branch",
        "entry `nq/nv/na` means distinct admissible q1 values / verified eps-q1",
        "certificates / aligned certificates.",
        "",
        "| w | w mod 20 | (-1|w) | (w|5) | L11 | ESC | qchar +1 (q/v/a) | qchar -1 (q/v/a) | criterion | classification | other note |",
        "|---:|---:|---:|---:|:---|:---|:---|:---|:---|:---|:---|",
    ]
    for record in records:
        l11_witness = record["l11"]["witness"]
        esc_witness = record["esc"]["first_aligned"]
        l11_cell = (
            f"yes a={l11_witness['a']},q={l11_witness['q1']}"
            if l11_witness
            else "no root"
            if not record["l11"]["roots"]
            else "bounded no-hit"
        )
        esc_cell = (
            f"yes eps={esc_witness['eps']},q={esc_witness['q1']}"
            if esc_witness
            else "no"
        )
        criterion = (
            "wall"
            if record["esc"]["w_character_mod5"] == 1
            else "qchar +1 escapes 5"
        )
        lines.append(
            "| {w} | {mod} | {minus} | {chi} | {l11} | {esc} | {plus} | {minus_branch} | {criterion} | {classification} ({label}) | {note} |".format(
                w=record["w"],
                mod=record["w_mod_20"],
                minus=record["l11"]["minus_one_legendre_w"],
                chi=record["esc"]["w_character_mod5"],
                l11=l11_cell,
                esc=esc_cell,
                plus=branch_cell(record, 1),
                minus_branch=branch_cell(record, -1),
                criterion=criterion,
                classification=record["classification"],
                label=record["classification_label"],
                note=obstruction_note(record),
            )
        )

    lines.extend(
        [
            "",
            "## Verdict",
            "",
            f"- PROVED identity mismatches: {summary['identity_mismatch_count']}.",
            f"- Kernel refusals on admissible ESC attempts: {summary['eligible_esc_refusal_count']} (never evidence).",
            f"- Both nonzero q1-character branches had verified attempts for all {summary['prime_w_count']} primes: {summary['both_branches_covered_for_all_w']}.",
            "- Classification counts: " + json.dumps(summary["classification_counts"], sort_keys=True) + ".",
            f"- w=131-type criterion counterexamples: {summary['criterion_counterexamples']}.",
            f"- Bounded L11 iff `w == 1 (mod 4)` mismatches: {summary['l11_scan_mismatches']}.",
            f"- Bounded canonical ESC iff `(w|5) == -1` mismatches: {summary['esc_scan_mismatches']}.",
            f"- Full bounded classification versus the mod-20 criterion mismatches: {summary['classification_mismatches']}.",
            "- PROVED scope warning: the alternate `a=7,tau=0,f=131,q1=41` Route131",
            "  closes `w=131`; residues 11 and 19 classify only the frozen canonical protocol.",
            "- `other-obstructed` means only that the bounded 5-compatible branch had no aligned",
            "  base class; the table names any common negative frozen place and makes no universal claim.",
            "",
            f"Data: `{OUT}`",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    h10q._selftest()
    q1_primes = list(h10q.primerange(2, Q1_MAX + 1))
    w_primes = list(h10q.primerange(W_MIN, W_MAX + 1))
    records = [per_w_record(w, q1_primes) for w in w_primes]

    mismatches = [
        {"w": record["w"], **mismatch}
        for record in records
        for mismatch in record["esc"]["identity_mismatches"]
    ]
    l11_scan_mismatches = [
        record["w"]
        for record in records
        if (record["l11"]["witness"] is not None)
        != (record["l11"]["minus_one_legendre_w"] == 1)
    ]
    esc_scan_mismatches = [
        record["w"]
        for record in records
        if (record["esc"]["first_aligned"] is not None)
        != (record["esc"]["w_character_mod5"] == -1)
    ]
    classification_mismatches = []
    for record in records:
        expected = (
            "L11 route available"
            if record["l11"]["minus_one_legendre_w"] == 1
            else "ESC route available"
            if record["esc"]["w_character_mod5"] == -1
            else "5-obstructed"
        )
        if record["classification"] != expected:
            classification_mismatches.append(
                {
                    "w": record["w"],
                    "expected": expected,
                    "observed": record["classification"],
                }
            )
    both_covered = all(
        all(
            record["esc"]["branches_by_q1_character_mod5"][str(character)][
                "verified_certificate_count"
            ]
            > 0
            for character in (1, -1)
        )
        for record in records
    )
    eligible_refusals = sum(
        sum(
            sum(branch["refusal_counts"].values())
            for branch in record["esc"]["branches_by_q1_character_mod5"].values()
        )
        for record in records
    )
    classification_counts = dict(
        sorted(Counter(record["classification"] for record in records).items())
    )
    summary = {
        "type": "summary",
        "label": "EVIDENCE",
        "prime_w_count": len(records),
        "classification_counts": classification_counts,
        "identity_mismatch_count": len(mismatches),
        "identity_mismatches": mismatches,
        "eligible_esc_refusal_count": eligible_refusals,
        "both_branches_covered_for_all_w": both_covered,
        "l11_scan_mismatches": l11_scan_mismatches,
        "esc_scan_mismatches": esc_scan_mismatches,
        "classification_mismatches": classification_mismatches,
        "criterion_counterexamples": classification_mismatches,
        "criterion_verdict": (
            "PROVED formula; bounded audit agrees"
            if not mismatches and not classification_mismatches
            else "counterexample-listed"
        ),
    }
    write_jsonl(records, summary)
    write_report(records, summary)
    print("WROTE", OUT)
    print("WROTE", REPORT)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
