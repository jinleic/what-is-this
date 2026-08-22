#!/usr/bin/env python3
"""L17 two-layer rate model on the predeclared censored cohort.

For each class and each p<=10000, evaluate P(b(k)) modulo p at every k mod p.
A simple root with symbol -1 contributes exact odd-valuation p-adic mass
1/(p+1), not 1/p. Singular roots are recursively lifted with a Taylor
valuation tree; unresolved roots exclude that (class,p) pair and log a bound.

The unknown refusal/zero-jacobi rows are handled as a two-sided band: c_lo
counts every unknown trial as non-success, c_hi counts every unknown trial as
success. Both are family-stratified fits on prime-Q exposure; classes without
zero are right-censored at the evaluated exposure. No final-rung conditioning
is used because open classes have null rung.
"""
import json
import math
import os
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from fractions import Fraction
from math import comb

ROOT = os.environ.get("H10Q_ROOT", str(Path(__file__).resolve().parent))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from h10q import _is_prime, _l10_P, _l10_supp, primerange, vp  # noqa: E402

LIMIT = int(os.environ.get("L17_P_LIMIT", "10000"))
MAX_K = int(os.environ.get("L17_MAX_K", "200"))
AUTHORITY_FILES = (
    os.path.join(ROOT, "data", "l13h_l11_a.jsonl"),
    os.path.join(ROOT, "data", "l13h_l11_b.jsonl"),
    os.path.join(ROOT, "data", "l13h_esc.jsonl"),
)
L15_DENSITY = os.path.join(ROOT, "data", "l15_density.json")
DEFAULT_OUT = os.path.join(ROOT, "data", "l17_ratemodel_censored.jsonl")
PRIMES = tuple(primerange(2, LIMIT + 1))
SINGULAR_MAX_DEPTH = 12
SINGULAR_MAX_NODES = 200_000


def _mod_fraction(x, p):
    return (x.numerator % p) * pow(x.denominator % p, -1, p) % p


def _compose_affine(P, B, S):
    K = [Fraction(0)] * len(P)
    for i, ci in enumerate(P):
        for j in range(i + 1):
            K[j] += ci * comb(i, j) * B ** (i - j) * S ** j
    return K


def _taylor_coeffs(K, a):
    return [
        sum(K[i] * comb(i, j) * a ** (i - j) for i in range(j, len(K)))
        for j in range(len(K))
    ]


def _valuation_or_none(x, p):
    return None if x == 0 else vp(x, p)


def _singular_bad_mass(K, p, root):
    """Return exact odd-v_p mass and resolved prefix, or exclude residual."""
    stack = [(int(root), 1)]
    mass = Fraction(0)
    nodes = 0
    max_depth = 1
    unresolved = []
    while stack:
        if nodes >= SINGULAR_MAX_NODES:
            unresolved.extend(stack)
            break
        a, d = stack.pop()
        nodes += 1
        max_depth = max(max_depth, d)
        T = _taylor_coeffs(K, a)
        v0 = _valuation_or_none(T[0], p)
        tails = []
        for j in range(1, len(T)):
            vj = _valuation_or_none(T[j], p)
            if vj is not None:
                tails.append(vj + d * j)
        vtail = min(tails) if tails else None
        if (v0 is None and vtail is None) or (
            v0 is not None and (vtail is None or v0 < vtail)
        ):
            if v0 is not None and v0 % 2:
                mass += Fraction(1, p ** d)
            continue
        if d >= SINGULAR_MAX_DEPTH:
            unresolved.append((a, d))
            continue
        pd = p ** d
        for t in range(p):
            stack.append((a + t * pd, d + 1))
    residual = sum((Fraction(1, p ** d) for _, d in unresolved), Fraction(0))
    info = {
        "nodes": nodes,
        "max_depth": max_depth,
        "unresolved_nodes": len(unresolved),
        "resolved_mass": str(mass),
        "resolved_mass_num": mass.numerator,
        "resolved_mass_den": mass.denominator,
        "residual_mass_upper_bound": str(residual),
        "residual_mass_upper_bound_float": float(residual),
        "method": "recursive Taylor valuation tree",
    }
    if unresolved:
        info["status"] = "excluded-unresolved-singular-branch"
        return None, info
    info["status"] = "certified"
    return mass, info


def _class_small_factor(r):
    w, (u1, u2) = r["cell"]
    a = Fraction(r["a"])
    A = 1 + 4 * a * a
    tau = (1 + 2 * a * a) / A
    z = Fraction(w) * Fraction(u1, u2)
    delta = 1 - A * tau * tau
    alpha = -delta * A
    Z = z ** 3
    D = 1 - Z - a * a * Z * Z
    s = (a - 1) / 2
    P = _l10_P(a, Z, D, A, delta, s)
    eps, f, q1, N = int(r["eps"]), int(r["f"]), int(r["q1"]), int(r["N"])
    B, S = eps * f * q1, eps * f * N
    K = _compose_affine(P, B, S)

    fixed = {2, 3, 5, 7} | _l10_supp(alpha) | _l10_supp(delta)
    denominator_support = set()
    for c in P:
        denominator_support |= _l10_supp(Fraction(c.denominator))
    missing = denominator_support - fixed
    if missing:
        raise RuntimeError((r["family"], r["cell"], sorted(missing)))

    table = []
    root_stats = []
    bad_root_lists = []
    product = 1.0
    n_scanned = 0
    n_simple = n_step = n_bad_simple = n_bad_step = 0
    excluded_pairs = []
    for p in PRIMES:
        local_mass = Fraction(0)
        roots_total = simple = step_roots = bad_simple = bad_step = 0
        step_details = []
        bad_entries = []
        excluded = False
        if p not in fixed:
            n_scanned += 1
            cm = [_mod_fraction(c, p) for c in P]
            dP = [i * cm[i] for i in range(1, len(cm))]
            am = _mod_fraction(alpha, p)
            bmod = B % p
            step_mod = S % p
            for residue in range(p):
                val = cm[-1]
                for c in reversed(cm[:-1]):
                    val = (val * bmod + c) % p
                if val == 0:
                    roots_total += 1
                    dval = 0
                    if dP:
                        dval = dP[-1]
                        for c in reversed(dP[:-1]):
                            dval = (dval * bmod + c) % p
                    is_simple = bool(step_mod and dval)
                    if is_simple:
                        simple += 1
                    else:
                        step_roots += 1
                    arg = (2 * am * bmod) % p
                    bad_sign = bool(arg and pow(arg, (p - 1) // 2, p) == p - 1)
                    if bad_sign and is_simple:
                        bad_simple += 1
                        bad_entries.append({
                            "r": residue,
                            "kind": "simple",
                            "status": "certified",
                            "mass_num": 1,
                            "mass_den": p + 1,
                            "resolved_mass_num": 1,
                            "resolved_mass_den": p + 1,
                        })
                    elif bad_sign:
                        bad_step += 1
                        smass, sinfo = _singular_bad_mass(K, p, residue)
                        entry = {"r": residue, "kind": "step", **sinfo}
                        if smass is None:
                            excluded = True
                            entry["mass_num"] = None
                            entry["mass_den"] = None
                        else:
                            local_mass += smass
                            entry["mass_num"] = smass.numerator
                            entry["mass_den"] = smass.denominator
                        bad_entries.append(entry)
                bmod += step_mod
                if bmod >= p:
                    bmod -= p
            if not excluded:
                local_mass += Fraction(bad_simple, p + 1)
                if not (0 <= local_mass <= 1):
                    raise RuntimeError((r["cell"], p, local_mass))
                product *= float(1 - local_mass)
            else:
                excluded_pairs.append({
                    "p": p,
                    "roots_total": roots_total,
                    "step_bad_roots": bad_step,
                    "step_details": step_details,
                    "disposition": "excluded from product; unresolved branch",
                })
        if bad_entries:
            bad_root_lists.append({
                "p": p,
                "roots": bad_entries,
                "local_mass_num": local_mass.numerator if not excluded else None,
                "local_mass_den": local_mass.denominator if not excluded else None,
                "excluded": excluded,
            })
        table.append(
            [p, local_mass.numerator, local_mass.denominator]
            if not excluded else [p, 0, 1]
        )
        if roots_total or excluded:
            root_stats.append({
                "p": p,
                "roots_total": roots_total,
                "simple_roots": simple,
                "step_roots": step_roots,
                "bad_simple_roots": bad_simple,
                "bad_step_roots": bad_step,
                "bad_mass_num": local_mass.numerator if not excluded else 0,
                "bad_mass_den": local_mass.denominator if not excluded else 1,
                "excluded": excluded,
                "step_details": step_details if step_details else None,
            })
        n_simple += simple
        n_step += step_roots
        n_bad_simple += bad_simple
        n_bad_step += bad_step
    return product, table, {
        "prime_limit": LIMIT,
        "root_stats": root_stats,
        "bad_root_lists": bad_root_lists,
        "n_primes_total": len(PRIMES),
        "n_primes_scanned": n_scanned,
        "fixed_primes": sorted(p for p in PRIMES if p in fixed),
        "table_encoding": "[p,mass_num,mass_den] for every p; exact odd-v_p bad mass in Z_p",
        "n_simple_roots": n_simple,
        "n_step_roots": n_step,
        "n_bad_simple_roots": n_bad_simple,
        "n_bad_step_roots": n_bad_step,
        "n_excluded_pairs": len(excluded_pairs),
        "excluded_pairs": excluded_pairs,
        "partial": bool(excluded_pairs),
        "capped": True,
        "singular_depth_cap": SINGULAR_MAX_DEPTH,
        "singular_node_cap": SINGULAR_MAX_NODES,
        "simple_root_mass": "bad_simple_roots/(p+1)",
        "step_root_method": "recursive Taylor valuation tree; unresolved pairs excluded with residual bound",
    }


def _rank_average(values):
    pairs = sorted((v, i) for i, v in enumerate(values))
    ranks = [0.0] * len(values)
    j = 0
    while j < len(pairs):
        k = j + 1
        while k < len(pairs) and pairs[k][0] == pairs[j][0]:
            k += 1
        rank = (j + 1 + k) / 2.0
        for _, i in pairs[j:k]:
            ranks[i] = rank
        j = k
    return ranks


def _pearson(x, y):
    if not x or not y:
        return None
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denx = math.sqrt(sum((a - mx) ** 2 for a in x))
    deny = math.sqrt(sum((b - my) ** 2 for b in y))
    return num / (denx * deny) if denx and deny else None


def _spearman(x, y):
    return _pearson(_rank_average(x), _rank_average(y))


def _load_classes():
    rows = []
    for path in AUTHORITY_FILES:
        with open(path) as fh:
            for line in fh:
                r = json.loads(line)
                if r.get("type") == "class":
                    r["_cohort_file"] = os.path.basename(path)
                    rows.append(r)
    assert len(rows) == 293
    assert len({(r["family"], r["class_idx"]) for r in rows}) == 293
    return rows


def _known_counts(r):
    v = r["verdicts"]
    event = int(v.get("zero", 0))
    bad = sum(int(n) for k, n in v.items() if k.startswith("bad"))
    unknown = sum(int(n) for k, n in v.items() if k.startswith("refused") or k == "zero-jacobi")
    evaluated = sum(int(n) for n in v.values())
    return event, bad, unknown, evaluated


def _fit_c(rows, optimistic):
    """Family MLE; unknown rows are failures (lo) or successes (hi)."""
    groups = defaultdict(list)
    for row in rows:
        groups[row["family"]].append(row)
    out = {}
    for family in sorted(groups):
        rs = groups[family]
        successes = sum(r["n_known_events"] + (r["n_unknown"] if optimistic else 0) for r in rs)
        max_s = max(r["small_factor"] for r in rs)
        hi = min(1.0, (1.0 / max_s) * (1.0 - 1e-14))

        def score(c):
            ans = successes / c
            for r in rs:
                failures = r["n_evaluated_prime_Q_trials"] - r["n_known_events"]
                if optimistic:
                    failures -= r["n_unknown"]
                ans -= failures * r["small_factor"] / (1.0 - c * r["small_factor"])
            return ans

        if score(hi) >= 0:
            c = hi
            boundary = True
        else:
            lo = 1e-15
            for _ in range(120):
                mid = (lo + hi) / 2
                if score(mid) > 0:
                    lo = mid
                else:
                    hi = mid
            c = (lo + hi) / 2
            boundary = False
        ll = 0.0
        for r in rs:
            p = c * r["small_factor"]
            succ = r["n_known_events"] + (r["n_unknown"] if optimistic else 0)
            fail = r["n_evaluated_prime_Q_trials"] - succ
            if succ:
                ll += succ * math.log(p)
            if fail:
                ll += fail * math.log1p(-p)
        out[family] = {
            "family": family,
            "c": c,
            "n_classes": len(rs),
            "n_success_trials": successes,
            "n_failure_trials": sum(r["n_evaluated_prime_Q_trials"] for r in rs) - successes,
            "n_unknown_rows": sum(r["n_unknown"] for r in rs),
            "log_likelihood": ll,
            "boundary_mle": boundary,
            "calibration": "family-stratified trial MLE with right-censored no-zero rows; unknown as " + ("success" if optimistic else "non-success") + "; empirical, not theorem",
        }
    return out


def _class_label(row):
    u = row["cell"][1]
    return f"{row['family']}|{row['cell'][0]}|({u[0]},{u[1]})"

def _prime_counts_by_k(row, max_k=None):
    if max_k is None:
        max_k = MAX_K
    q1, N = int(row["q1"]), int(row["N"])
    n = 0
    out = []
    refusals = []
    for k in range(max_k + 1):
        Q = q1 + k * N
        try:
            prime = _is_prime(Q)
        except Exception as exc:
            prime = False
            refusals.append({"k": k, "Q": str(Q), "reason": type(exc).__name__})
        if prime:
            n += 1
        out.append(n)
    return out, refusals


def _shape_mix(rows):
    all_shapes = Counter()
    by_family = defaultdict(Counter)
    for r in rows:
        if r["rung_shape_closed"] is not None:
            all_shapes[r["rung_shape_closed"]] += 1
            by_family[r["family"]][r["rung_shape_closed"]] += 1
    return {"all_closed_rows": dict(all_shapes), "by_family_closed_rows": {f: dict(c) for f, c in sorted(by_family.items())}, "R=1": 0, "note": "descriptive only; no final-rung calibration; L15 R=1 absent"}


def _top(rows, tag, reverse):
    eligible = [r for r in rows if r["event_observed"] and r["observed_k_zero"] is not None and r["observed_k_zero"] > 0]
    eligible.sort(key=lambda r: r[f"ratio_member_trials_{tag}"], reverse=reverse)
    out = []
    for r in eligible[:10]:
        out.append({"class": _class_label(r), "cell": r["cell"], "family": r["family"], "k_zero": r["observed_k_zero"], "observed_member_trials": r["observed_member_trials"], "small_factor": r["small_factor"], "cofactor_c": r[f"cofactor_c_{tag}"], "model_wait_member_trials": r[f"model_wait_{tag}"], "ratio_member_trials": r[f"ratio_member_trials_{tag}"], "ratio_k_zero": r[f"ratio_k_zero_{tag}"]})
    return out


def build(out_path=DEFAULT_OUT, limit=None):
    source = _load_classes()
    if limit is not None:
        source = source[:limit]
    t0 = time.time()
    rows = []
    for i, src in enumerate(source, 1):
        sf, table, smeta = _class_small_factor(src)
        n_event, n_bad, n_unknown, n_eval = _known_counts(src)
        assert n_event == (1 if src["status"] == "closed" else 0)
        rows.append({"type": "class", "cell": src["cell"], "family": src["family"], "class_idx": src["class_idx"], "small_factor": sf, "small_log10_factor": math.log10(sf) if sf > 0 else None, "f_p_table": table, "f_p_meta": smeta, "observed_k_zero": src["k_zero"], "rung_shape_closed": src["rung"], "event_observed": bool(n_event), "right_censored": not bool(n_event), "n_evaluated_prime_Q_trials": n_eval, "observed_member_trials": n_eval if n_event else None, "censor_member_trials": n_eval if not n_event else None, "n_known_events": n_event, "n_known_bad": n_bad, "n_unknown": n_unknown, "n_nonprime_skips": int(src["nonprime_skips"]), "n_excluded_skips": int(src.get("excluded_skips", 0)), "cohort_file": src["_cohort_file"], "q1": int(src["q1"]), "N": int(src["N"])})
        rows[-1]["bad_root_lists"] = smeta["bad_root_lists"]
        rows[-1]["excluded_mass_upper_bound"] = sum(
            d.get("residual_mass_upper_bound_float", 0.0)
            for ex in smeta["excluded_pairs"]
            for d in ex["step_details"]
        )
        print(f"small {i}/{len(source)} {src['family']} {src['cell']} sf={sf:.9g}", flush=True)
    fits = {"lo": _fit_c(rows, False), "hi": _fit_c(rows, True)}
    for row in rows:
        for tag in ("lo", "hi"):
            c = fits[tag][row["family"]]["c"]
            row[f"cofactor_c_{tag}"] = c
            row[f"model_rate_{tag}"] = row["small_factor"] * c
            row[f"model_wait_{tag}"] = 1.0 / row[f"model_rate_{tag}"]
            if row["event_observed"]:
                row[f"ratio_member_trials_{tag}"] = row[f"model_wait_{tag}"] / row["observed_member_trials"]
                row[f"ratio_k_zero_{tag}"] = row[f"model_wait_{tag}"] / row["observed_k_zero"] if row["observed_k_zero"] > 0 else None
            else:
                row[f"ratio_member_trials_{tag}"] = None
                row[f"ratio_k_zero_{tag}"] = None
        row["cofactor_c_band"] = [row["cofactor_c_lo"], row["cofactor_c_hi"]]
        row["model_wait_band"] = [row["model_wait_hi"], row["model_wait_lo"]]
        row["model_wait"] = row["model_wait_lo"]
        row["ratio"] = row["ratio_k_zero_lo"]
    prime_refusals = []
    for row, src in zip(rows, source):
        prefix, refusals = _prime_counts_by_k(src)
        row["prime_Q_prefix_counts_k0_to_200"] = prefix
        prime_refusals.extend(refusals)

    events = [r for r in rows if r["event_observed"]]
    tail = [r for r in events if r["observed_k_zero"] > 0]
    fit_summary = {"n_tail": len(tail), "ratio_definition": "model_wait / event prime-member trial count; raw ratio field uses k_zero"}
    for tag in ("lo", "hi"):
        pm = [r[f"model_wait_{tag}"] for r in tail]
        om = [float(r["observed_member_trials"]) for r in tail]
        pk = pm
        ok = [float(r["observed_k_zero"]) for r in tail]
        lm = [math.log10(a / b) for a, b in zip(pm, om)]
        lk = [math.log10(a / b) for a, b in zip(pk, ok)]
        fit_summary.update({f"{tag}_tail_rank_corr_member_trials": _spearman(pm, om), f"{tag}_tail_log10_member_wait_over_observed_median": statistics.median(lm) if lm else None, f"{tag}_tail_rank_corr_raw_k_zero": _spearman(pk, ok), f"{tag}_tail_log10_wait_over_k_zero_median": statistics.median(lk) if lk else None})

    curve = []
    for K in range(MAX_K + 1):
        pred_lo = pred_hi = 0.0
        actual = 0
        for r in rows:
            n = r["prime_Q_prefix_counts_k0_to_200"][K]
            pred_lo += 1.0 - (1.0 - r["model_rate_lo"]) ** n
            pred_hi += 1.0 - (1.0 - r["model_rate_hi"]) ** n
            if r["event_observed"] and r["observed_k_zero"] <= K:
                actual += 1
        curve.append({"k": K, "predicted_cumulative_closures_lo": pred_lo, "predicted_cumulative_closures_hi": pred_hi, "observed_wave1_cumulative_closures": actual, "wave12_target_at_k200": 117 if K == 200 else None})
    with open(L15_DENSITY) as fh:
        measured_p0 = float(json.load(fh)["ii_fits"]["p0_mle"]["global"])
    mean_small = sum(r["small_factor"] for r in rows) / len(rows)
    global_rates = {tag: sum(r[f"model_rate_{tag}"] for r in rows) / len(rows) for tag in ("lo", "hi")}
    summary = {"type": "summary", "n_classes": len(rows), "n_event_classes_wave1": len(events), "n_tail_event_classes_k_gt_0": len(tail), "n_right_censored_classes": sum(r["right_censored"] for r in rows), "calibration_constants_band": fits, "closed_rung_shape_mix": _shape_mix(rows), "fit": fit_summary, "aggregate": {"curve_k0_to_200": curve, "predicted_total_closures_k_le_200_lo": curve[-1]["predicted_cumulative_closures_lo"], "predicted_total_closures_k_le_200_hi": curve[-1]["predicted_cumulative_closures_hi"], "observed_wave1_total_k_le_200": curve[-1]["observed_wave1_cumulative_closures"], "actual_wave1_plus_wave2_target_k_le_200": 117, "comparison_to_117_band": [curve[-1]["predicted_cumulative_closures_lo"] - 117, curve[-1]["predicted_cumulative_closures_hi"] - 117], "curve_units": "progression k; prime-Q prefixes via kernel _is_prime", "right_censoring": "all no-zero classes censored after evaluated prime-Q exposure through k<=200"}, "global_law_validation_band": {"measured_p0_l15": measured_p0, "predicted_mean_member_rate_band": [global_rates["lo"], global_rates["hi"]], "predicted_over_measured_band": [global_rates["lo"] / measured_p0, global_rates["hi"] / measured_p0], "mean_small_factor": mean_small, "effective_cofactor_factor_band": [global_rates["lo"] / mean_small, global_rates["hi"] / mean_small], "measured_bateman_horn_gap": 15.7, "small_layer_suppression": 1.0 / mean_small, "cofactor_layer_suppression_after_small_band": [mean_small / global_rates["hi"], mean_small / global_rates["lo"]], "dominance_note": "cofactor layer dominates the low/identified side; the unknown band must be reported and may be too wide to resolve"}, "top10_overpredicted_tail_member_trials_lo": _top(rows, "lo", True), "top10_underpredicted_tail_member_trials_lo": _top(rows, "lo", False), "top10_overpredicted_tail_member_trials_hi": _top(rows, "hi", True), "top10_underpredicted_tail_member_trials_hi": _top(rows, "hi", False), "prime_count_refusals": prime_refusals, "provenance": {"cohort_authority": [os.path.relpath(p, ROOT) for p in AUTHORITY_FILES], "small_kernel": "h10q._l10_P, _l10_supp, primerange; direct all-residue modular evaluation; no sympy", "strip_convention": "{2,3,5,7}|supp(alpha)|supp(delta)|supp(b); p|b has Legendre 0", "prime_claim_policy": "all p and Q claims use h10q kernel; no probabilistic claims", "unknown_policy": "c_lo treats every refused/zero-jacobi trial as non-success; c_hi treats every such trial as success; both bands are reported", "padic_note": "simple root bad mass=1/(p+1); singular roots recursively Taylor-lifted; unresolved pairs excluded with residual bound", "rung_note": "closed prime/factorint mix descriptive only; no final-rung conditioning; L15 R=1 absent"}, "runtime_seconds": time.time() - t0}
    with open(out_path, "w") as out:
        out.write(json.dumps({"type": "meta", "artifact": "l17_ratemodel_censored", "cohort_files": [os.path.relpath(p, ROOT) for p in AUTHORITY_FILES], "n_classes": len(rows), "small_prime_cap": LIMIT, "prime_count": len(PRIMES), "small_factor_definition": "product_{p<=10000}(1-f_p), f_p exact odd-v_p bad mass; simple f_p=bad_simple/(p+1)", "f_p_table_encoding": "[p,mass_num,mass_den] all p; roots and singular dispositions in f_p_meta", "cofactor_definition": "family c_lo/c_hi missing-outcome band; empirical not theorem", "censoring_definition": "no-zero class right-censored after all evaluated prime-Q trials through k<=200", "rung_shape_note": "closed rung mix descriptive; no final-rung fit; R=1 absent"}, separators=(",", ":")) + "\n")
        for row in rows:
            out.write(json.dumps(row, separators=(",", ":")) + "\n")
        out.write(json.dumps(summary, separators=(",", ":")) + "\n")
    print(json.dumps({"out": out_path, "n": len(rows), "events": len(events), "tail": len(tail), "runtime_seconds": time.time() - t0, "summary": summary}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    out = DEFAULT_OUT
    limit = None
    args = list(sys.argv[1:])
    if args and args[0] == "--limit":
        limit = int(args[1]); args = args[2:]
    if args and args[0] == "--out":
        out = args[1]; args = args[2:]
    if args:
        raise SystemExit(f"usage: {sys.argv[0]} [--limit N] [--out PATH]")
    build(out, limit)
