"""Certificate-grade correction replay for the frozen q=3 orbit ladder.

The adjudication is frozen in the correction campaign pre_statement.md.  This
runner uses total orbit masses for both the invariant input p and output
reference D, and directly checks every input word before accepting compression.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.machinery
import importlib.util
import itertools
import json
import math
import os
import signal
import sys
import time
import traceback
from pathlib import Path

from flint import arb, ctx, fmpq

ROOT = Path("/Users/jinleic/jinleic-workspace/cs/delcap")
CORRECTION = ROOT / "campaigns" / (
    "2026-08-31T09:13:25Z_7dc5babe-5e1b-44fa-a9e5-03d2f10183b6_q3-invariance-correction"
)
FROZEN = ROOT / "campaigns" / (
    "2026-08-30T17:50:25Z_ddd57c90-d088-47c1-acda-dcb54bcf3555"
)


def load_dependency(name: str, path: Path):
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


q4 = load_dependency(
    "delcap_q4_v3_dependency",
    CORRECTION / "gate_c_q4_v3_dependency.py.asrun",
)
frozen = q4.frozen

PREC = 400
SNAP_IN = 1 << 30
SNAP_OUT = 1 << 30
BA_ITERS = 4000
D_GRID = (fmpq(1, 2), fmpq(1, 5), fmpq(1, 10), fmpq(1, 20))
ROW_TIMEOUT_S = 90 * 60
STAGE_TIMEOUT_S = 4 * 60 * 60
ROWS_PATH = CORRECTION / "corrected_q3_rows.jsonl"
LINE_HASH_PATH = CORRECTION / "corrected_q3_rows.line_checksums.json"
ANCHOR_PATH = CORRECTION / "anchor_check.json"
SUMMARY_PATH = CORRECTION / "summary.json"
FAILURE_PATH = CORRECTION / "FAILURE.json"
LOG_PATH = CORRECTION / "run.log"
FROZEN_ROWS_PATH = FROZEN / "orbit_rows.jsonl"


def utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def log(message: str) -> None:
    line = f"[{utc()}] {message}"
    print(line, flush=True)
    with LOG_PATH.open("a") as handle:
        handle.write(line + "\n")


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def row_key(row: dict) -> tuple[int, int, str]:
    return int(row["q"]), int(row["n"]), str(row["d"])


def row_key_text(key: tuple[int, int, str]) -> str:
    return f"{key[0]}:{key[1]}:{key[2]}"


def verified_done_keys() -> set[tuple[int, int, str]]:
    if not ROWS_PATH.exists() or ROWS_PATH.stat().st_size == 0:
        return set()
    assert LINE_HASH_PATH.exists()
    ledger = json.loads(LINE_HASH_PATH.read_text())
    assert ledger["algorithm"] == "sha256"
    expected = ledger["rows"]
    done = set()
    seen = set()
    for raw in ROWS_PATH.read_bytes().splitlines(keepends=True):
        if not raw.strip():
            continue
        row = json.loads(raw)
        key = row_key(row)
        text = row_key_text(key)
        assert key not in done
        assert hashlib.sha256(raw).hexdigest() == expected[text]
        done.add(key)
        seen.add(text)
    assert seen == set(expected)
    return done


def append_row(row: dict) -> None:
    done = verified_done_keys()
    key = row_key(row)
    if key in done:
        return
    ledger = (
        json.loads(LINE_HASH_PATH.read_text())
        if LINE_HASH_PATH.exists()
        else {"algorithm": "sha256", "rows": {}}
    )
    raw = (json.dumps(row, sort_keys=True) + "\n").encode()
    with ROWS_PATH.open("ab") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())
    ledger["rows"][row_key_text(key)] = hashlib.sha256(raw).hexdigest()
    write_json(LINE_HASH_PATH, ledger)


def grid() -> list[tuple[int, int, fmpq]]:
    return [(3, n, d) for n in range(6, 11) for d in D_GRID]


def generators(q: int, word: tuple[int, ...]):
    swap = list(range(q))
    swap[0], swap[1] = swap[1], swap[0]
    cycle = list(range(1, q)) + [0]
    return (
        tuple(reversed(word)),
        tuple(swap[symbol] for symbol in word),
        tuple(cycle[symbol] for symbol in word),
    )


def assert_input_distribution(structure, masses: list[int]) -> dict:
    total = sum(masses)
    assert total == SNAP_IN
    per_word = [
        fmpq(masses[j], total * structure.od.gin.sizes[j])
        for j in range(len(masses))
    ]
    assert sum(
        structure.od.gin.sizes[j] * per_word[j]
        for j in range(len(masses))
    ) == 1
    checks = 0
    positive_words = 0
    for word in itertools.product(range(structure.q), repeat=structure.n):
        orbit = structure.od.gin.word_orbit[word]
        value = per_word[orbit]
        positive_words += int(value > 0)
        for transformed in generators(structure.q, word):
            assert value == per_word[structure.od.gin.word_orbit[transformed]]
            checks += 1
    return {
        "orbit_total_mass_numerators": masses,
        "total_mass_denominator": total,
        "per_word_formula": "mass(O)/(S*|O|)",
        "orbit_sizes": list(structure.od.gin.sizes),
        "expanded_g_generator_checks": checks,
        "positive_input_words": positive_words,
        "zero_input_words": structure.q**structure.n - positive_words,
        "exact_sum": "1",
        "status": "PASS",
        "evidence_label": "MACHINE-VERIFIED",
    }


def kl_evaluation(
    structure,
    dp: list[fmpq],
    d: fmpq,
    *,
    input_masses: list[int] | None = None,
) -> dict:
    """Direct every-input KL and representative overlap; optional primal sum."""
    ctx.prec = PREC
    log2 = arb(2).log()
    a, b, denominator = q4.d_parts(d, structure.n)
    cache: dict[tuple[int, int], arb] = {}
    entries = 0
    positive_checks = 0

    def value(word: tuple[int, ...]) -> arb:
        nonlocal entries, positive_checks
        subtotal = arb(0)
        for output, count in frozen.dpdict(word, structure.n).items():
            entries += 1
            t = structure.od.flat_of_word(output)
            dy = dp[t]
            assert dy > 0
            positive_checks += 1
            key = (t, count)
            term = cache.get(key)
            if term is None:
                k = len(output)
                weight = fmpq(count * a ** (structure.n - k) * b**k, denominator)
                term = arb(weight) * (arb(weight / dy).log() / log2)
                cache[key] = term
            subtotal += term
        return subtotal

    if input_masses is None:
        representative_indices = list(range(len(structure.od.gin.reps)))
    else:
        representative_indices = [
            j for j, mass in enumerate(input_masses) if mass > 0
        ]
    assert representative_indices
    rep_values = [None] * len(structure.od.gin.reps)
    for j in representative_indices:
        rep_values[j] = value(structure.od.gin.reps[j])
    rep_argmax = max(
        representative_indices, key=lambda j: rep_values[j].upper()
    )
    rep_best = rep_values[rep_argmax]
    full_best = None
    full_argmax = None
    overlap_failures = 0
    full_mutual = arb(0) if input_masses is not None else None
    input_total = sum(input_masses) if input_masses is not None else None
    positive_input_words = 0
    for index, word in enumerate(
        itertools.product(range(structure.q), repeat=structure.n)
    ):
        orbit = structure.od.gin.word_orbit[word]
        if input_masses is not None and input_masses[orbit] == 0:
            continue
        current = value(word)
        positive_input_words += 1
        representative = rep_values[orbit]
        assert representative is not None
        if (
            current.upper() < representative.lower()
            or representative.upper() < current.lower()
        ):
            overlap_failures += 1
        if full_best is None or current.upper() > full_best.upper():
            full_best = current
            full_argmax = word
        if input_masses is not None:
            px = fmpq(
                input_masses[orbit],
                input_total * structure.od.gin.sizes[orbit],
            )
            full_mutual += arb(px) * current
        if index and index % 10000 == 0:
            log(
                f"  full KL q={structure.q} n={structure.n}: "
                f"{index}/{structure.q**structure.n} inputs"
            )
    assert full_best is not None
    assert overlap_failures == 0
    assert not (
        full_best.upper() < rep_best.lower()
        or rep_best.upper() < full_best.lower()
    )
    result = {
        "representative_max_ball": str(rep_best),
        "full_alphabet_max_ball": str(full_best),
        "representative_max_lower": float(rep_best.lower()),
        "representative_max_upper": float(rep_best.upper()),
        "full_alphabet_max_lower": float(full_best.lower()),
        "full_alphabet_max_upper": float(full_best.upper()),
        "representative_argmax_orbit": rep_argmax,
        "full_argmax_word": list(full_argmax),
        "full_argmax_orbit": structure.od.gin.word_orbit[full_argmax],
        "representatives_evaluated": len(representative_indices),
        "positive_full_input_words_evaluated": positive_input_words,
        "zero_mass_input_words_skipped": (
            structure.q**structure.n - positive_input_words
            if input_masses is not None
            else 0
        ),
        "full_word_vs_representative_overlap_failures": overlap_failures,
        "conditional_entries_evaluated": entries,
        "positive_W_terms_with_positive_D_checked": positive_checks,
        "cached_exact_terms": len(cache),
        "status": "PASS",
        "evidence_label": "MACHINE-VERIFIED",
        "_representative_max_arb": rep_best,
        "_full_alphabet_max_arb": full_best,
    }
    if full_mutual is not None:
        orbit_mutual = arb(0)
        for j in representative_indices:
            representative = rep_values[j]
            assert representative is not None
            orbit_mutual += (
                arb(fmpq(input_masses[j], input_total)) * representative
            )
        assert not (
            full_mutual.upper() < orbit_mutual.lower()
            or orbit_mutual.upper() < full_mutual.lower()
        )
        result.update(
            {
                "orbit_weighted_primal_ball": str(orbit_mutual),
                "full_alphabet_primal_ball": str(full_mutual),
                "orbit_weighted_primal_lower": float(orbit_mutual.lower()),
                "orbit_weighted_primal_upper": float(orbit_mutual.upper()),
                "full_alphabet_primal_lower": float(full_mutual.lower()),
                "full_alphabet_primal_upper": float(full_mutual.upper()),
                "full_vs_orbit_primal_overlap": "PASS",
                "_orbit_weighted_primal_arb": orbit_mutual,
                "_full_alphabet_primal_arb": full_mutual,
            }
        )
    return result


def output_candidate(structure, name: str, masses: list[int], d: fmpq) -> dict:
    total = sum(masses)
    assert len(masses) == structure.od.n_out
    assert all(value > 0 for value in masses)
    dp = [
        fmpq(masses[t], total * int(structure.flat_sizes[t]))
        for t in range(len(masses))
    ]
    assert sum(int(structure.flat_sizes[t]) * dp[t] for t in range(len(dp))) == 1
    invariance_checks = q4.assert_expanded_g_invariance(structure, dp)
    evaluation = kl_evaluation(structure, dp, d)
    representative_arb = evaluation.pop("_representative_max_arb")
    full_arb = evaluation.pop("_full_alphabet_max_arb")
    return {
        "name": name,
        "status": "ADMITTED_EXACT_OUTPUT_ORBIT_INVARIANT",
        "output_orbit_total_mass_numerators": masses,
        "total_mass_denominator": total,
        "output_orbit_sizes": [int(value) for value in structure.flat_sizes],
        "per_word_formula": "mass(T)/(S*|T|)",
        "exact_sum": "1",
        "minimum_orbit_mass_numerator": min(masses),
        "expanded_g_generator_checks": invariance_checks,
        "evaluation": evaluation,
        "_representative_max_arb": representative_arb,
        "_full_alphabet_max_arb": full_arb,
        "evidence_label": "MACHINE-VERIFIED",
    }


def legacy_word_snap_audit(q: int, n: int, d: fmpq) -> dict:
    od = frozen.OrbitData(q, n, d, log=lambda *_: None)
    os_rows = frozen.merged_OS(od)
    _, dword, rate, dual = frozen.ba_orbit_fine(od, os_rows, iters=BA_ITERS)
    words = sorted(q4.word_iter(q, n))
    orbit_ids = []
    floats = []
    for output in words:
        gout = od.gout[len(output)]
        orbit = gout.word_orbit[output]
        value = None
        for key, candidate in dword.items():
            if len(key) == len(output) and gout.word_orbit[key] == orbit:
                value = candidate
                break
        assert value is not None
        orbit_ids.append(od.base[len(output)] + orbit)
        floats.append(float(value))
    counts = frozen.snap(floats, SNAP_OUT)
    if min(counts) == 0:
        counts = [value + 1 for value in counts]
    values: dict[int, set[int]] = {}
    digest = hashlib.sha256()
    for orbit, count in zip(orbit_ids, counts):
        values.setdefault(orbit, set()).add(count)
        digest.update(int(count).to_bytes(8, "little", signed=False))
    flat_sizes = [
        size for k in range(n + 1) for size in od.gout[k].sizes
    ]
    mixed = [
        {
            "flat_output_orbit": orbit,
            "orbit_size": flat_sizes[orbit],
            "distinct_word_counts": sorted(items),
        }
        for orbit, items in sorted(values.items())
        if len(items) != 1
    ]
    return {
        "candidate": "legacy_ba_word",
        "disposition": (
            "REJECTED_NONINVARIANT_BEFORE_BOUND"
            if mixed
            else "LEGACY_VECTOR_EXACTLY_INVARIANT_NOT_ADMITTED_BY_NEW_RULE"
        ),
        "mixed_output_orbits": len(mixed),
        "mixed_output_orbit_details": mixed,
        "word_count": len(words),
        "word_count_denominator": sum(counts),
        "word_counts_sha256_le_u64": digest.hexdigest(),
        "ba_rate_float_block": rate,
        "ba_dual_float_block": dual,
        "integer_audit_evidence_label": "MACHINE-VERIFIED",
        "float_locator_evidence_label": "COMPUTATIONAL-EVIDENCE",
    }


def correction_row(structure, d: fmpq, old: dict) -> dict:
    start_wall = time.perf_counter()
    start_cpu = time.process_time()
    ctx.prec = PREC
    locator = q4.locator(structure, d, iters=BA_ITERS)

    input_masses = frozen.snap(locator["m"], SNAP_IN)
    input_archive = assert_input_distribution(structure, input_masses)
    orbit_primal, induced = q4.cert_primal_structure(structure, input_masses, d)
    assert sum(
        int(structure.flat_sizes[t]) * induced[t]
        for t in range(len(induced))
    ) == 1
    induced_invariance_checks = q4.assert_expanded_g_invariance(structure, induced)
    primal_check = kl_evaluation(
        structure, induced, d, input_masses=input_masses
    )
    primal_check.pop("_representative_max_arb")
    primal_check.pop("_full_alphabet_max_arb")
    orbit_weighted_check = primal_check.pop("_orbit_weighted_primal_arb")
    full_primal_check = primal_check.pop("_full_alphabet_primal_arb")
    assert not (
        orbit_primal.upper() < full_primal_check.lower()
        or full_primal_check.upper() < orbit_primal.lower()
    )
    assert not (
        orbit_primal.upper() < orbit_weighted_check.lower()
        or orbit_weighted_check.upper() < orbit_primal.lower()
    )

    total_float_masses = [
        float(locator["dm_orbit"][t]) * int(structure.flat_sizes[t])
        for t in range(structure.od.n_out)
    ]
    assert abs(sum(total_float_masses) - 1.0) < 1e-9
    ba_masses = frozen.snap(total_float_masses, SNAP_OUT)
    ba_zero_before_bump = sum(value == 0 for value in ba_masses)
    if ba_zero_before_bump:
        ba_masses = [value + 1 for value in ba_masses]
    uniform_masses = [int(value) for value in structure.flat_sizes]
    candidates = [
        output_candidate(structure, "ba_total_orbit_mass", ba_masses, d),
        output_candidate(structure, "uniform", uniform_masses, d),
    ]
    for candidate in candidates:
        candidate["zero_orbit_masses_before_full_support_bump"] = (
            ba_zero_before_bump if candidate["name"] == "ba_total_orbit_mass" else 0
        )
    finite = [
        candidate
        for candidate in candidates
        if candidate["_full_alphabet_max_arb"].upper() != float("inf")
    ]
    assert finite
    chosen = min(
        finite,
        key=lambda candidate: candidate["_full_alphabet_max_arb"].upper(),
    )
    dual = chosen["_full_alphabet_max_arb"]

    sandwich, delta = q4.sandwich_structure(structure, d)
    primal = orbit_primal
    lo = float((primal / structure.n).lower())
    hi = float((dual / structure.n).upper())
    assert hi >= lo
    lower_margin_ball = (
        (primal / structure.n).lower() - sandwich["lbplus"].upper()
    )
    upper_margin_ball = (
        sandwich["ub"].lower() - (dual / structure.n).upper()
    )
    lower_margin = float(lower_margin_ball)
    upper_margin = float(upper_margin_ball)
    verdicts = []
    if upper_margin_ball > 0:
        verdicts.append("CERT_UPPER_BEATS_UB")
    if lower_margin_ball > 0:
        verdicts.append("CERT_LOWER_BEATS_LBplus")
    if not verdicts:
        verdicts.append("NO_STRICT_IMPROVEMENT")
    lbplus = float(sandwich["lbplus"])
    ub = float(sandwich["ub"])
    if (
        (dual / structure.n).upper() < sandwich["lbplus"].lower()
        or (primal / structure.n).lower() > sandwich["ub"].upper()
    ):
        raise AssertionError("corrected interval excludes published sandwich")
    for candidate in candidates:
        candidate.pop("_representative_max_arb")
        candidate.pop("_full_alphabet_max_arb")

    legacy = legacy_word_snap_audit(3, structure.n, d)
    row = {
        "status": "CORRECTED_CERTIFIED",
        "q": 3,
        "n": structure.n,
        "d": str(d),
        "old_source": str(FROZEN_ROWS_PATH),
        "old_interval": [
            old["cert_lo_per_symbol"],
            old["cert_hi_per_symbol"],
        ],
        "old_dual_from": old["dual_from"],
        "legacy_output_reference_audit": legacy,
        "replacement_source": "gate_c_q3_invariance_correction.py.asrun",
        "input_distribution": input_archive,
        "induced_output_exact_sum": "1",
        "induced_output_expanded_g_generator_checks": induced_invariance_checks,
        "primal_orbit_ball": str(primal),
        "primal_full_check": primal_check,
        "dual_candidates": candidates,
        "dual_selection_rule": (
            "minimum finite direct-full-alphabet Arb upper endpoint; "
            "fixed tie order ba_total_orbit_mass then uniform"
        ),
        "dual_from": chosen["name"],
        "dual_full_alphabet_ball": chosen["evaluation"]["full_alphabet_max_ball"],
        "lbplus": lbplus,
        "ub": ub,
        "delta_n": float(delta),
        "cert_lo_per_symbol": lo,
        "cert_hi_per_symbol": hi,
        "cert_width_per_symbol": hi - lo,
        "lower_gain_over_lbplus": lower_margin,
        "upper_gain_below_ub": upper_margin,
        "width_sign_check": "PASS",
        "published_sandwich_containment": "PASS",
        "verdict": "+".join(verdicts),
        "locator": {
            "rate_float_block": locator["rate"],
            "dual_float_block": locator["dual"],
            "weighted_mass_float": locator["weighted_mass"],
            "wall_s": locator["wall_s"],
            "evidence_label": "COMPUTATIONAL-EVIDENCE",
        },
        "n_input_words": 3**structure.n,
        "n_input_orbits": len(structure.od.gin.reps),
        "n_output_words": sum(3**k for k in range(structure.n + 1)),
        "n_output_orbits": structure.od.n_out,
        "structure_build_wall_s": structure.build_wall_s,
        "wall_s": round(time.perf_counter() - start_wall, 3),
        "cpu_s": round(time.process_time() - start_cpu, 3),
        "evidence_label": "MACHINE-VERIFIED",
    }
    return row


class ResourceStop(Exception):
    pass


ALARM_REASON = "resource stop"


def alarm_handler(signum, frame):
    del signum, frame
    raise ResourceStop(ALARM_REASON)


def startup_anchors(old_by_key: dict) -> dict:
    start_wall = time.perf_counter()
    start_cpu = time.process_time()
    result = {"status": "IN_PROGRESS", "started_utc": utc()}
    write_json(ANCHOR_PATH, result)
    frozen._selftest(log=log)
    result["frozen_A1_A6"] = "PASS"
    raw = frozen.raw_snap_dual_check(2, 10, fmpq(1, 20))
    assert raw["n_zeroed"] == 1 and raw["dual_raw_is_inf"]
    result["raw_inf_plant"] = raw
    actual = frozen.row(3, 6, fmpq(1, 2), log=lambda *_: None)
    expected = old_by_key[(3, 6, "1/2")]
    q4.assert_frozen_row(actual, expected)
    result["frozen_3_6_half_bit_anchor"] = "PASS"
    split = q4.split_snap_full_alphabet_counterfactual()
    assert split["full_lower_minus_representative_upper"] > 0
    result["out_of_grid_split_plant"] = split
    result["status"] = "PASS"
    result["wall_s"] = round(time.perf_counter() - start_wall, 3)
    result["cpu_s"] = round(time.process_time() - start_cpu, 3)
    result["completed_utc"] = utc()
    write_json(ANCHOR_PATH, result)
    return result


def run_stage() -> None:
    global ALARM_REASON
    start_wall = time.perf_counter()
    start_cpu = time.process_time()
    deadline = time.monotonic() + STAGE_TIMEOUT_S
    old_rows = [
        json.loads(line)
        for line in FROZEN_ROWS_PATH.read_text().splitlines()
        if line.strip()
    ]
    old_by_key = {row_key(row): row for row in old_rows}
    assert len(old_by_key) == 20
    startup_anchors(old_by_key)
    done = verified_done_keys()
    old_handler = signal.signal(signal.SIGALRM, alarm_handler)
    try:
        for n in range(6, 11):
            remaining_ds = [
                d for d in D_GRID if (3, n, str(d)) not in done
            ]
            if not remaining_ds:
                continue
            remaining_stage = deadline - time.monotonic()
            assert remaining_stage > 0
            ALARM_REASON = "4-hour correction-stage wall reached during structure build"
            signal.alarm(max(1, math.ceil(remaining_stage)))
            structure = q4.build_structure(3, n, logger=log)
            signal.alarm(0)
            for d in remaining_ds:
                remaining_stage = deadline - time.monotonic()
                assert remaining_stage > 0
                alarm_s = min(ROW_TIMEOUT_S, max(1, math.ceil(remaining_stage)))
                ALARM_REASON = (
                    "90-minute correction-row wall reached"
                    if alarm_s == ROW_TIMEOUT_S
                    else "4-hour correction-stage wall reached"
                )
                signal.alarm(alarm_s)
                try:
                    row = correction_row(
                        structure, d, old_by_key[(3, n, str(d))]
                    )
                except Exception as error:
                    failure = {
                        "status": "FAIL",
                        "q": 3,
                        "n": n,
                        "d": str(d),
                        "error": repr(error),
                        "traceback": traceback.format_exc(),
                        "completed_rows": len(verified_done_keys()),
                        "utc": utc(),
                    }
                    write_json(FAILURE_PATH, failure)
                    raise
                finally:
                    signal.alarm(0)
                append_row(row)
                done.add((3, n, str(d)))
                log(
                    f"corrected q=3 n={n} d={d}: old={row['old_interval']} "
                    f"new=[{row['cert_lo_per_symbol']:.12f},"
                    f"{row['cert_hi_per_symbol']:.12f}] "
                    f"legacy_mixed={row['legacy_output_reference_audit']['mixed_output_orbits']} "
                    f"{row['verdict']} [{row['wall_s']}s wall/{row['cpu_s']}s cpu]"
                )
            del structure
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
    rows = [
        json.loads(line)
        for line in ROWS_PATH.read_text().splitlines()
        if line.strip()
    ]
    assert len(rows) == 20
    assert len(verified_done_keys()) == 20
    summary = {
        "status": "PASS",
        "rows": len(rows),
        "legacy_split_rows": sum(
            row["legacy_output_reference_audit"]["mixed_output_orbits"] > 0
            for row in rows
        ),
        "both_end_improvements": sum(
            "CERT_LOWER_BEATS_LBplus" in row["verdict"]
            and "CERT_UPPER_BEATS_UB" in row["verdict"]
            for row in rows
        ),
        "max_n": max(row["n"] for row in rows),
        "wall_s": round(time.perf_counter() - start_wall, 3),
        "cpu_s": round(time.process_time() - start_cpu, 3),
        "completed_utc": utc(),
        "evidence_label": "MACHINE-VERIFIED",
        "timing_evidence_label": "COMPUTATIONAL-EVIDENCE",
    }
    write_json(SUMMARY_PATH, summary)
    log(f"correction stage PASS {summary}")


def table_stage() -> None:
    rows = [
        json.loads(line)
        for line in ROWS_PATH.read_text().splitlines()
        if line.strip()
    ]
    assert len(rows) == len(grid()) == len(verified_done_keys())
    fields = [
        "n",
        "d",
        "old_lo",
        "old_hi",
        "new_lo",
        "new_hi",
        "width",
        "lbplus",
        "ub",
        "lower_gain",
        "upper_gain",
        "legacy_mixed_output_orbits",
        "dual_from",
        "verdict",
        "wall_s",
        "cpu_s",
    ]
    with (CORRECTION / "TABLE_corrections.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in sorted(rows, key=lambda item: (item["n"], D_GRID.index(fmpq(item["d"])))):
            writer.writerow(
                {
                    "n": row["n"],
                    "d": row["d"],
                    "old_lo": row["old_interval"][0],
                    "old_hi": row["old_interval"][1],
                    "new_lo": row["cert_lo_per_symbol"],
                    "new_hi": row["cert_hi_per_symbol"],
                    "width": row["cert_width_per_symbol"],
                    "lbplus": row["lbplus"],
                    "ub": row["ub"],
                    "lower_gain": row["lower_gain_over_lbplus"],
                    "upper_gain": row["upper_gain_below_ub"],
                    "legacy_mixed_output_orbits": row["legacy_output_reference_audit"]["mixed_output_orbits"],
                    "dual_from": row["dual_from"],
                    "verdict": row["verdict"],
                    "wall_s": row["wall_s"],
                    "cpu_s": row["cpu_s"],
                }
            )
    log("correction table PASS")


def selfcheck() -> None:
    old_rows = [
        json.loads(line)
        for line in FROZEN_ROWS_PATH.read_text().splitlines()
        if line.strip()
    ]
    old_by_key = {row_key(row): row for row in old_rows}
    structure = q4.build_structure(3, 3, logger=lambda *_: None)
    locator = q4.locator(structure, fmpq(1, 2), iters=32)
    input_masses = frozen.snap(locator["m"], SNAP_IN)
    assert_input_distribution(structure, input_masses)
    _, induced = q4.cert_primal_structure(structure, input_masses, fmpq(1, 2))
    kl_evaluation(structure, induced, fmpq(1, 2), input_masses=input_masses)
    total_float = [
        float(locator["dm_orbit"][t]) * int(structure.flat_sizes[t])
        for t in range(structure.od.n_out)
    ]
    masses = frozen.snap(total_float, SNAP_OUT)
    if min(masses) == 0:
        masses = [value + 1 for value in masses]
    output_candidate(structure, "selfcheck", masses, fmpq(1, 2))
    split = q4.split_snap_full_alphabet_counterfactual()
    assert split["full_lower_minus_representative_upper"] > 0
    assert old_by_key[(3, 6, "1/2")]["dual_from"] == "ba_word"
    print("correction selfcheck PASS", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("selfcheck", "run", "table"))
    args = parser.parse_args()
    ctx.prec = PREC
    code_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    log(f"stage={args.stage}; code_sha256={code_hash}; Arb={PREC} bits")
    if args.stage == "selfcheck":
        selfcheck()
    elif args.stage == "run":
        run_stage()
    else:
        table_stage()


if __name__ == "__main__":
    main()
