#!/usr/bin/env python3
"""Independent bounded verifier for the H681 two-column trace-nine canary.

Rebuilds the exact QQ(q) lift through the frozen Wave-27 verification stack
(tests/test_trace_nine_norm_envelope.py and its e251 base), pins the full
35-rule digest, checks that F4..F8 reduce to zero, releases every lift-only
object, and then computes only the constant and the maximum-total-degree
multiplication columns one at a time with a streaming canonical digest.

Scope: two exact columns of the trace-nine multiplication matrix.  No full
matrix, norm, height, exceptional-root, or thermodynamic claim follows.
Live memory and timings go to stdout only; they are never part of the data.
"""
from __future__ import annotations

import argparse
import copy
import gc
import hashlib
import importlib.util
import json
import lzma
import math
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from flint import fmpz_poly
from sympy.core.cache import clear_cache

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent.parent
SELF = ROOT / "tests" / "test_trace_nine_bounded_canary.py"
PRODUCER = ROOT / "experiments" / "e254_trace_nine_bounded_canary.py"
GUARD = ROOT / "tools" / "resource_guard.py"
ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_bounded_canary.json"
TEMPLATE = ROOT / "results" / "spectral" / "trace_nine_lift_template.json.xz"
WAVE27_PRODUCER = ROOT / "experiments" / "e252_trace_nine_norm_envelope.py"
WAVE27_VERIFIER = ROOT / "tests" / "test_trace_nine_norm_envelope.py"
WAVE27_PROOF = ROOT / "proofs" / "trace_nine_norm_envelope.md"
WAVE27_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_norm_envelope.json"
E251_PRODUCER = ROOT / "experiments" / "e251_trace_nine_projection.py"
E251_VERIFIER = ROOT / "tests" / "test_trace_nine_projection.py"
E251_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_projection.json"
E248_PRODUCER = ROOT / "experiments" / "e248_replica_trace_nine.py"
E248_VERIFIER = ROOT / "tests" / "test_replica_trace_nine.py"
TRACE_ARTIFACT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
LOCKFILE = ROOT / "uv.lock"
RUNTIME_FREEZE = ROOT.parent / "requirements-freeze.txt"
SOURCE_PATHS = (
    PRODUCER,
    SELF,
    GUARD,
    TEMPLATE,
    WAVE27_PRODUCER,
    WAVE27_VERIFIER,
    WAVE27_PROOF,
    WAVE27_ARTIFACT,
    E251_PRODUCER,
    E251_VERIFIER,
    E251_ARTIFACT,
    E248_PRODUCER,
    E248_VERIFIER,
    TRACE_ARTIFACT,
    LOCKFILE,
    RUNTIME_FREEZE,
)

SCOPE = "BOUNDED_TWO_COLUMN_CANARY_ONLY"
QUOTIENT_RANK = 96
RULE_COUNT = 35
PURE_POWER_BOUNDS = [2, 2, 3, 5, 9]
HILBERT_VECTOR = [1, 5, 12, 19, 22, 19, 12, 5, 1]
RULE_SHA256 = "f4cfe038d7625c6a411f2e01b9af12d9ae9ca06d980e64ad38d8eff030fe22cc"
CONSTANT_MONOMIAL = (0, 0, 0, 0, 0)
LARGEST_MONOMIAL = (0, 0, 0, 0, 8)
EXPECTED_COLUMNS = [
    {
        "basis_monomial": [0, 0, 0, 0, 0],
        "column_scalar_lcm_bits": 123,
        "column_sha256": "69065ff900271f44f3f34856b1bad85fd148182b8a045aefd1cc66cb1b8bdfb9",
        "entries": 37,
        "maximum_denominator_degree": 145,
        "maximum_numerator_degree": 181,
        "maximum_numerator_l1_log2": 224,
        "ordered_reduction_steps": 19,
    },
    {
        "basis_monomial": [0, 0, 0, 0, 8],
        "column_scalar_lcm_bits": 1321,
        "column_sha256": "63b8a477b5d9dd25a354d1bd2f7d18666cbcf1bf21cb6060d0c2b1b51441f47b",
        "entries": 96,
        "maximum_denominator_degree": 513,
        "maximum_numerator_degree": 629,
        "maximum_numerator_l1_log2": 1771,
        "ordered_reduction_steps": 578,
    },
]
META_KEYS = frozenset({"schema_version", "source_sha256", "data_sha256"})
DATA_KEYS = frozenset(
    {
        "scope",
        "quotient_rank",
        "final_rational_rule_sha256",
        "columns",
        "full_matrix_status",
        "physical_root_status",
    }
)
COLUMN_KEYS = frozenset(EXPECTED_COLUMNS[0])
COLUMN_COUNT_KEYS = (
    "entries",
    "ordered_reduction_steps",
    "maximum_numerator_degree",
    "maximum_denominator_degree",
    "maximum_numerator_l1_log2",
    "column_scalar_lcm_bits",
)
MUTATION_IDS = [
    "M1_largest_column_digest_last_nibble",
    "M2_constant_column_entry_count_plus_one",
    "M3_column_order_swapped",
    "M4_rule_digest_last_nibble",
    "M5_quotient_rank_minus_one",
    "M6_producer_source_hash_last_nibble",
    "M7_full_matrix_status_escalated",
]
HEX_DIGITS = frozenset("0123456789abcdef")
Monomial = tuple[int, int, int, int, int]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


wave27 = load_module("bounded_canary_wave27_verifier", WAVE27_VERIFIER)
e251 = wave27.base
guard = load_module("bounded_canary_resource_guard", GUARD)
Rat = wave27.Rat


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expected_sources() -> dict[str, str]:
    return {str(path.relative_to(WORKSPACE)): file_sha256(path) for path in SOURCE_PATHS}


def memory_stage(stage: str, started: float, wall: float) -> dict[str, object]:
    """Print CURRENT residency/footprint next to the cumulative peak; stdout only."""
    gc.collect()
    sample = guard.sample_process(os.getpid())
    record = {
        "stage": stage,
        "cpu_seconds": round(time.process_time() - started, 3),
        "wall_seconds": round(time.monotonic() - wall, 3),
        "current_resident_bytes": int(sample["resident_bytes"]),
        "current_footprint_bytes": int(sample["footprint_bytes"]),
        "peak_rss_bytes": wave27.peak_rss_bytes(),
    }
    print(
        "TRACE_NINE_BOUNDED_CANARY_VERIFIER STAGE "
        + json.dumps(record, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    return record


def stream_integer_list(digest, polynomial: fmpz_poly) -> None:
    digest.update(b"[")
    separator = b""
    for coefficient in polynomial:
        digest.update(separator)
        separator = b","
        digest.update(str(int(coefficient)).encode())
    digest.update(b"]")


def stream_column_sha256(column: dict[Monomial, Rat]) -> str:
    """Streamed canonical digest of sorted [[monomial],[numerator,denominator]] entries.

    Feeds one decimal coefficient at a time; never materializes the JSON text.
    Byte-identical to json.dumps(sort_keys=True, separators=(",", ":")) of the
    H679 column encoding, which reference_column_sha256 spells out.
    """
    digest = hashlib.sha256()
    digest.update(b"[")
    separator = b""
    for monomial in sorted(column):
        value = column[monomial]
        digest.update(separator)
        separator = b","
        digest.update(b"[[")
        digest.update(",".join(str(int(power)) for power in monomial).encode())
        digest.update(b"],[")
        stream_integer_list(digest, value.numerator)
        digest.update(b",")
        stream_integer_list(digest, value.denominator)
        digest.update(b"]]")
    digest.update(b"]")
    return digest.hexdigest()


def stream_rule_sha256(rules: list[tuple[Monomial, dict[Monomial, Rat]]]) -> str:
    return canonical_sha256([stream_column_sha256(tail) for _, tail in rules])


def reference_column_sha256(column: dict[Monomial, Rat]) -> str:
    """The H679/H680 column digest definition, materialized; smoke reference only."""
    return canonical_sha256(
        [
            [
                list(monomial),
                [list(map(int, value.numerator)), list(map(int, value.denominator))],
            ]
            for monomial, value in sorted(column.items())
        ]
    )


def reference_rule_sha256(rules: list[tuple[Monomial, dict[Monomial, Rat]]]) -> str:
    return canonical_sha256([reference_column_sha256(tail) for _, tail in rules])


def load_template() -> tuple[dict[str, object], str]:
    """Decompress and parse the frozen lift template without recompressing it.

    The Wave-27 verifier certifies canonical XZ by recompressing at preset 9,
    a transient of several hundred MiB; this bounded verifier pins the
    compressed file digest through provenance and the payload digest through
    the inherited Wave-27 artifact instead.
    """
    raw = lzma.decompress(TEMPLATE.read_bytes(), format=lzma.FORMAT_XZ)
    payload = json.loads(raw)
    if raw != canonical_bytes(payload):
        raise AssertionError("lift template is not canonical JSON")
    if payload.get("schema") != "ising3d.trace-nine-lift-template/v1":
        raise AssertionError("unexpected lift-template schema")
    return payload, hashlib.sha256(raw).hexdigest()


def inherited_wave27() -> dict[str, object]:
    artifact = json.loads(WAVE27_ARTIFACT.read_text())
    if artifact["meta"]["data_sha256"] != canonical_sha256(artifact["data"]):
        raise AssertionError("inherited Wave-27 data digest mismatch")
    if not all(check["passed"] for check in artifact["checks"]):
        raise AssertionError("inherited Wave-27 artifact contains a failed check")
    data = artifact["data"]
    if data["physical_lift"]["final_rule_sha256"] != RULE_SHA256:
        raise AssertionError("inherited Wave-27 rule digest is not the pinned digest")
    if data["template"]["compressed_sha256"] != file_sha256(TEMPLATE):
        raise AssertionError("inherited Wave-27 template file digest mismatch")
    if data["template"]["standard_monomial_count"] != QUOTIENT_RANK:
        raise AssertionError("inherited Wave-27 quotient rank is not 96")
    if data["nonzero_witness"]["source_artifact_sha256"] != file_sha256(E251_ARTIFACT):
        raise AssertionError("inherited Wave-27 and e251 artifacts are not cross-linked")
    return data


def inherited_trace_rows() -> list[list[int]]:
    trace = json.loads(TRACE_ARTIFACT.read_text())
    rows = [[int(value) for value in row] for row in trace["data"]["trace_coefficients_ascending"]]
    digest = canonical_sha256([[str(value) for value in row] for row in rows])
    if digest != trace["data"]["trace_coefficients_sha256"]:
        raise AssertionError("inherited trace coefficient digest mismatch")
    projection = json.loads(E251_ARTIFACT.read_text())
    if projection["meta"]["data_sha256"] != canonical_sha256(projection["data"]):
        raise AssertionError("inherited e251 data digest mismatch")
    if not all(check["passed"] for check in projection["checks"]):
        raise AssertionError("inherited e251 artifact contains a failed check")
    provenance = next(
        check
        for check in projection["checks"]
        if check["name"] == "C0_inherited_e248_provenance"
    )
    if provenance["detail"] != digest:
        raise AssertionError("e251 witness and trace artifact are not cross-linked")
    return rows


def template_identities(
    records: list[dict[str, object]],
    leaders: list[Monomial],
    leading_forms: list[dict[Monomial, Fraction]],
    degrees: list[int],
) -> None:
    """Exact QQ check that every template basis is sum_i H_ji L_i, graded and reduced."""
    for index, (record, leading) in enumerate(zip(records, leaders)):
        basis = wave27.parse_q_polynomial(record["basis"])
        basis_degree = sum(leading)
        if basis.get(leading) != 1 or max(basis, key=wave27.grevlex_order) != leading:
            raise AssertionError(f"lift-template record {index} is not monic at its leader")
        if any(sum(monomial) != basis_degree for monomial in basis):
            raise AssertionError(f"lift-template basis {index} is not homogeneous")
        if any(
            wave27.divides(divisor, monomial)
            for monomial in basis
            if monomial != leading
            for divisor in leaders
        ):
            raise AssertionError(f"lift-template basis {index} is not reduced below its leader")
        if len(record["transform"]) != len(leading_forms):
            raise AssertionError(f"lift-template record {index} transformation arity changed")
        reconstructed: dict[Monomial, Fraction] = {}
        for entries, form, degree in zip(record["transform"], leading_forms, degrees):
            multiplier = wave27.parse_q_polynomial(entries)
            if any(sum(monomial) != basis_degree - degree for monomial in multiplier):
                raise AssertionError(f"lift-template transformation {index} is not graded")
            for left, scalar in multiplier.items():
                for right, coefficient in form.items():
                    monomial = wave27.add_monomials(left, right)
                    reconstructed[monomial] = reconstructed.get(monomial, 0) + scalar * coefficient
        if {monomial: value for monomial, value in reconstructed.items() if value} != basis:
            raise AssertionError(f"lift-template transformation identity {index} failed")


def evaluate_rows(compiled, targets: tuple[Rat, ...]) -> list[dict[Monomial, Rat]]:
    """Specialize every compiled equation row at the exact physical targets.

    Same arithmetic as the Wave-27 evaluate_target_polynomial; target powers
    are computed once and shared across rows instead of per term.
    """
    zero = Rat.scalar(0)
    powers: list[dict[int, Rat]] = [{} for _ in targets]
    rows: list[dict[Monomial, Rat]] = []
    for row in compiled:
        result: dict[Monomial, Rat] = {}
        for monomial, terms in row:
            value = zero
            for exponents, coefficient in terms:
                term = Rat.scalar(coefficient)
                for index, exponent in enumerate(exponents):
                    if exponent:
                        cache = powers[index]
                        power = cache.get(exponent)
                        if power is None:
                            power = cache[exponent] = targets[index] ** exponent
                        term = term * power
                value = value + term
            if value:
                result[monomial] = value
        rows.append(result)
    return rows


def rebuild_rules(
    records: list[dict[str, object]],
    leaders: list[Monomial],
    remainders: list[dict[Monomial, Rat]],
    started: float,
) -> tuple[list[tuple[Monomial, dict[Monomial, Rat]]], list[int]]:
    """Lift each template record on demand and interreduce in (degree, leader) order.

    Only already-interreduced tails are retained; no unreduced lift survives.
    """
    zero = Rat.scalar(0)
    one = Rat.scalar(1)
    processed: list[tuple[Monomial, dict[Monomial, Rat]]] = []
    by_index: dict[int, tuple[Monomial, dict[Monomial, Rat]]] = {}
    steps_by_index: dict[int, int] = {}
    for index in sorted(range(len(leaders)), key=lambda i: (sum(leaders[i]), leaders[i])):
        record = records[index]
        leading = leaders[index]
        polynomial = {
            monomial: Rat.scalar(coefficient)
            for monomial, coefficient in wave27.parse_q_polynomial(record["basis"]).items()
        }
        for source, entries in zip(remainders, record["transform"]):
            for left, scalar in wave27.parse_q_polynomial(entries).items():
                for right, coefficient in source.items():
                    monomial = wave27.add_monomials(left, right)
                    value = polynomial.get(monomial, zero) + coefficient * scalar
                    if value:
                        polynomial[monomial] = value
                    else:
                        polynomial.pop(monomial, None)
        if polynomial.pop(leading, None) != one:
            raise AssertionError(f"lift {index} lost its constant leader")
        tail, steps = wave27.reduce_exact(polynomial, processed)
        del polynomial
        by_index[index] = (leading, tail)
        processed.append((leading, tail))
        steps_by_index[index] = steps
        wave27.guard_resources(started, f"lift {index}")
    rules = [by_index[index] for index in range(len(leaders))]
    return rules, [steps_by_index[index] for index in range(len(leaders))]


def exact_column(
    basis_monomial: Monomial,
    trace_row: dict[Monomial, Rat],
    rules: list[tuple[Monomial, dict[Monomial, Rat]]],
    standard: frozenset[Monomial],
) -> tuple[dict[Monomial, Rat], int]:
    initial = {
        wave27.add_monomials(monomial, basis_monomial): value
        for monomial, value in trace_row.items()
    }
    column, steps = wave27.reduce_exact(initial, rules)
    del initial
    if not standard.issuperset(column):
        raise AssertionError("exact column contains a nonstandard monomial")
    return column, steps


def column_record(
    basis_monomial: Monomial, column: dict[Monomial, Rat], steps: int
) -> dict[str, object]:
    if not column:
        raise AssertionError("exact column is identically zero")
    scalar_lcm = 1
    numerator_degree = -1
    denominator_degree = -1
    numerator_l1_log2 = -1
    for value in column.values():
        scalar_lcm = math.lcm(scalar_lcm, abs(int(value.denominator.content())))
        numerator_degree = max(numerator_degree, value.numerator.degree())
        denominator_degree = max(denominator_degree, value.denominator.degree())
        numerator_l1_log2 = max(
            numerator_l1_log2,
            wave27.ceil_log2(sum(abs(int(coefficient)) for coefficient in value.numerator)),
        )
    return {
        "basis_monomial": list(basis_monomial),
        "entries": len(column),
        "ordered_reduction_steps": steps,
        "maximum_numerator_degree": numerator_degree,
        "maximum_denominator_degree": denominator_degree,
        "maximum_numerator_l1_log2": numerator_l1_log2,
        "column_scalar_lcm_bits": scalar_lcm.bit_length(),
        "column_sha256": stream_column_sha256(column),
    }


def recompute(started: float, wall: float) -> dict[str, object]:
    runtime = wave27.runtime_manifest()
    print(
        "TRACE_NINE_BOUNDED_CANARY_VERIFIER RUNTIME "
        + json.dumps(runtime, sort_keys=True, separators=(",", ":")),
        flush=True,
    )
    inherited = inherited_wave27()
    rows = inherited_trace_rows()
    template, payload_sha256 = load_template()
    if payload_sha256 != inherited["template"]["payload_sha256"]:
        raise AssertionError("lift-template payload digest differs from the Wave-27 artifact")
    records = template["records"]
    del template
    if len(records) != RULE_COUNT:
        raise AssertionError("lift-template record count changed")
    leaders = [tuple(record["leading_monomial"]) for record in records]
    if len(set(leaders)) != RULE_COUNT or any(
        len(leader) != 5 or any(type(power) is not int or power < 0 for power in leader)
        for leader in leaders
    ):
        raise AssertionError("lift-template leaders are not distinct 5-exponent monomials")
    pure_power_bounds, standard_list = wave27.standard_monomials(leaders)
    hilbert = [sum(sum(monomial) == degree for monomial in standard_list) for degree in range(9)]
    if (
        pure_power_bounds != PURE_POWER_BOUNDS
        or len(standard_list) != QUOTIENT_RANK
        or hilbert != HILBERT_VECTOR
    ):
        raise AssertionError("quotient dimensions changed")
    standard = frozenset(standard_list)
    memory_stage("inherited inputs loaded", started, wall)

    equations, _ = e251.independent_incidence()
    if len(equations) != 6:
        raise AssertionError("independent incidence must have six equations F4..F9")
    leading_forms = [
        {
            monomial: Fraction(int(sp.numer(coefficient)), int(sp.denom(coefficient)))
            for monomial, coefficient in wave27.homogeneous_leading(equation).terms()
        }
        for equation in equations[:5]
    ]
    degrees = [max(map(sum, form)) for form in leading_forms]
    template_identities(records, leaders, leading_forms, degrees)
    wave27.guard_resources(started, "template identities")
    memory_stage("template identities", started, wall)

    targets = tuple(Rat.from_sympy(value) for value in e251.target_expressions(rows).values())
    physical = evaluate_rows(e251.compile_rows(equations), targets)
    del equations, targets
    for index, form in enumerate(leading_forms):
        degree = degrees[index]
        observed = {
            monomial: value
            for monomial, value in physical[index].items()
            if sum(monomial) == degree
        }
        expected = {monomial: Rat.scalar(value) for monomial, value in form.items()}
        if observed != expected or any(sum(monomial) > degree for monomial in physical[index]):
            raise AssertionError(f"physical equation F{index + 4} changed its leading form")
    del leading_forms
    clear_cache()
    wave27.guard_resources(started, "physical targets")
    memory_stage("physical targets; sympy released", started, wall)

    remainders = [
        {monomial: value for monomial, value in physical[index].items() if sum(monomial) < degrees[index]}
        for index in range(5)
    ]
    rules, interreduction_steps = rebuild_rules(records, leaders, remainders, started)
    del records, remainders
    lift = inherited["physical_lift"]
    if interreduction_steps != lift["interreduction_steps_by_record"]:
        raise AssertionError("interreduction step transcript differs from the Wave-27 artifact")
    if [1 + len(tail) for _, tail in rules] != lift["final_rule_term_counts"]:
        raise AssertionError("final rule term counts differ from the Wave-27 artifact")
    if any(monomial not in standard for _, tail in rules for monomial in tail):
        raise AssertionError("interreduced tail contains a nonstandard monomial")
    rule_sha256 = stream_rule_sha256(rules)
    if rule_sha256 != RULE_SHA256:
        raise AssertionError(f"final rational rule digest {rule_sha256} is not the pinned digest")
    for index in range(5):
        remainder, _ = wave27.reduce_exact(physical[index], rules)
        if remainder:
            raise AssertionError(f"physical equation F{index + 4} does not reduce to zero")
    trace_row = physical[5]
    del physical
    scalar = 1
    for _, tail in rules:
        for value in tail.values():
            scalar = math.lcm(scalar, abs(int(value.denominator.content())))
    for value in trace_row.values():
        scalar = math.lcm(scalar, abs(int(value.denominator.content())))
    if scalar.bit_length() != lift["scalar_denominator_lcm_bits"]:
        raise AssertionError("lift scalar denominator LCM differs from the Wave-27 artifact")
    wave27.guard_resources(started, "validated lift")
    memory_stage("validated lift; template and F4..F8 released", started, wall)

    if CONSTANT_MONOMIAL not in standard:
        raise AssertionError("constant monomial is not standard")
    largest = max(standard_list, key=lambda monomial: (sum(monomial), monomial))
    if largest != LARGEST_MONOMIAL or sum(sum(monomial) == sum(largest) for monomial in standard_list) != 1:
        raise AssertionError("maximum-total-degree standard monomial is not the unique (0,0,0,0,8)")
    columns = []
    for label, basis_monomial in (("constant", CONSTANT_MONOMIAL), ("largest", largest)):
        column, steps = exact_column(basis_monomial, trace_row, rules, standard)
        columns.append(column_record(basis_monomial, column, steps))
        del column
        wave27.guard_resources(started, f"{label} column")
        memory_stage(f"{label} column digested and released", started, wall)
    if columns != EXPECTED_COLUMNS:
        raise AssertionError(f"reconstruction drifted from the pinned H679/H680 columns: {columns}")
    return {
        "quotient_rank": len(standard_list),
        "final_rational_rule_sha256": rule_sha256,
        "columns": columns,
    }


def validate_column_record(record: object) -> None:
    if not isinstance(record, dict) or set(record) != COLUMN_KEYS:
        raise AssertionError("column record keys mismatch")
    monomial = record["basis_monomial"]
    if (
        not isinstance(monomial, list)
        or len(monomial) != 5
        or any(type(power) is not int or power < 0 for power in monomial)
    ):
        raise AssertionError("column basis monomial is not a 5-exponent list")
    for key in COLUMN_COUNT_KEYS:
        if type(record[key]) is not int or record[key] < 0:
            raise AssertionError(f"column {key} is not a non-negative integer")
    if not 1 <= record["entries"] <= QUOTIENT_RANK:
        raise AssertionError("column entry count is outside 1..96")
    digest = record["column_sha256"]
    if not isinstance(digest, str) or len(digest) != 64 or not set(digest) <= HEX_DIGITS:
        raise AssertionError("column digest is not lowercase sha256 hex")


def validate_artifact(
    artifact: object,
    sources: dict[str, str],
    recomputed: dict[str, object] | None = None,
) -> None:
    if not isinstance(artifact, dict) or set(artifact) != {"data", "meta"}:
        raise AssertionError("artifact must contain exactly data and meta")
    meta = artifact["meta"]
    data = artifact["data"]
    if not isinstance(meta, dict) or set(meta) != META_KEYS:
        raise AssertionError("meta keys mismatch")
    if type(meta["schema_version"]) is not int or meta["schema_version"] != 1:
        raise AssertionError("unexpected schema version")
    if not isinstance(data, dict) or set(data) != DATA_KEYS:
        raise AssertionError("data keys mismatch")
    if meta["data_sha256"] != canonical_sha256(data):
        raise AssertionError("artifact data digest mismatch")
    if meta["source_sha256"] != sources:
        raise AssertionError("source provenance inventory or hash mismatch")
    if data["scope"] != SCOPE:
        raise AssertionError("artifact scope is not the bounded two-column canary")
    if data["full_matrix_status"] != "NOT_COMPUTED":
        raise AssertionError("artifact claims a full matrix status beyond the canary")
    if data["physical_root_status"] != "UNRESOLVED":
        raise AssertionError("artifact claims a physical root status beyond the canary")
    if type(data["quotient_rank"]) is not int or data["quotient_rank"] != QUOTIENT_RANK:
        raise AssertionError("quotient rank is not 96")
    if data["final_rational_rule_sha256"] != RULE_SHA256:
        raise AssertionError("final rational rule digest is not the pinned digest")
    columns = data["columns"]
    if not isinstance(columns, list) or len(columns) != 2:
        raise AssertionError("artifact must carry exactly two column records")
    for record in columns:
        validate_column_record(record)
    if columns != EXPECTED_COLUMNS:
        raise AssertionError("column records differ from the pinned H679/H680 observations")
    if recomputed is None:
        return
    if data["quotient_rank"] != recomputed["quotient_rank"]:
        raise AssertionError("quotient rank differs from the independent reconstruction")
    if data["final_rational_rule_sha256"] != recomputed["final_rational_rule_sha256"]:
        raise AssertionError("rule digest differs from the independent reconstruction")
    if columns != recomputed["columns"]:
        raise AssertionError("column records differ from the independent reconstruction")


def rejects(artifact: dict[str, object], sources: dict[str, str], recomputed: dict[str, object]) -> bool:
    try:
        validate_artifact(artifact, sources, recomputed)
    except AssertionError:
        return True
    return False


def flip_last_nibble(container: dict, key: str) -> None:
    digest = container[key]
    container[key] = digest[:-1] + ("0" if digest[-1] != "0" else "1")


def mutation_rejections(
    artifact: dict[str, object], sources: dict[str, str], recomputed: dict[str, object]
) -> list[str]:
    producer_key = str(PRODUCER.relative_to(WORKSPACE))
    outcomes = []

    def attempt(label: str, mutate, reseal: bool = True) -> None:
        mutated = copy.deepcopy(artifact)
        mutate(mutated)
        if reseal:
            mutated["meta"]["data_sha256"] = canonical_sha256(mutated["data"])
        outcomes.append(label if rejects(mutated, sources, recomputed) else "")

    def bump_entries(mutated: dict) -> None:
        mutated["data"]["columns"][0]["entries"] += 1

    def escalate(mutated: dict) -> None:
        mutated["data"]["full_matrix_status"] = "COMPUTED"

    attempt(MUTATION_IDS[0], lambda m: flip_last_nibble(m["data"]["columns"][1], "column_sha256"))
    attempt(MUTATION_IDS[1], bump_entries)
    attempt(MUTATION_IDS[2], lambda m: m["data"]["columns"].reverse())
    attempt(MUTATION_IDS[3], lambda m: flip_last_nibble(m["data"], "final_rational_rule_sha256"))
    attempt(MUTATION_IDS[4], lambda m: m["data"].__setitem__("quotient_rank", QUOTIENT_RANK - 1))
    attempt(
        MUTATION_IDS[5],
        lambda m: flip_last_nibble(m["meta"]["source_sha256"], producer_key),
        reseal=False,
    )
    attempt(MUTATION_IDS[6], escalate)
    return outcomes


def verify(artifact_path: Path) -> int:
    started = time.process_time()
    wall = time.monotonic()
    artifact_bytes = artifact_path.read_bytes()
    artifact = json.loads(artifact_bytes)
    if artifact_bytes != canonical_bytes(artifact):
        raise AssertionError("canary artifact is not canonical JSON")
    sources = expected_sources()
    validate_artifact(artifact, sources)
    memory_stage("artifact shape and provenance validated", started, wall)

    recomputed = recompute(started, wall)
    validate_artifact(artifact, sources, recomputed)
    rejected = mutation_rejections(artifact, sources, recomputed)
    if rejected != MUTATION_IDS:
        raise AssertionError(f"mutation rejection mismatch: {rejected}")
    if artifact_path.read_bytes() != artifact_bytes:
        raise AssertionError("artifact bytes changed during verification")
    final = memory_stage("complete", started, wall)
    columns = recomputed["columns"]
    print(
        "TRACE_NINE_BOUNDED_CANARY_VERIFIER PASS "
        f"scope={SCOPE} columns={len(columns)} mutations={len(rejected)} "
        f"rules_sha256={recomputed['final_rational_rule_sha256']} "
        f"constant_sha256={columns[0]['column_sha256']} "
        f"largest_sha256={columns[1]['column_sha256']} "
        f"cpu={final['cpu_seconds']}s wall={final['wall_seconds']}s "
        f"current_resident={final['current_resident_bytes']} "
        f"current_footprint={final['current_footprint_bytes']} "
        f"peak_rss={final['peak_rss_bytes']}",
        flush=True,
    )
    return 0


def smoke_columns() -> list[tuple[str, dict[Monomial, Rat]]]:
    def rat(numerator: list[int], denominator: list[int]) -> Rat:
        return Rat(fmpz_poly(numerator), fmpz_poly(denominator))

    huge = 2**1771
    big = 2**1321
    mixed = {
        (0, 1, 0, 12, 3): rat([-7, 0, 123456789], [1, -1]),
        (0, 0, 0, 0, 10): rat([0, 1], [3]),
        (0, 0, 0, 1, 0): rat([huge - 1, -huge, 0, 1], [big + 1, 5]),
        (2, 0, 0, 0, 0): rat([], [1]),
    }
    return [
        ("empty", {}),
        ("zero_numerator", {CONSTANT_MONOMIAL: rat([], [1])}),
        ("unit", {CONSTANT_MONOMIAL: rat([1], [1])}),
        ("mixed_signs_sizes_exponents", mixed),
        ("reverse_insertion_order", dict(reversed(list(mixed.items())))),
    ]


def smoke() -> int:
    started = time.process_time()
    wall = time.monotonic()
    cases = smoke_columns()
    digests = {}
    for name, column in cases:
        expected = reference_column_sha256(column)
        observed = stream_column_sha256(column)
        if observed != expected:
            raise AssertionError(f"streaming column digest differs from the reference: {name}")
        digests[name] = observed
    if digests["mixed_signs_sizes_exponents"] != digests["reverse_insertion_order"]:
        raise AssertionError("column digest depends on insertion order")
    if len(set(digests.values())) != len(digests) - 1:
        raise AssertionError("distinct smoke columns collided or the empty digest failed")
    mixed = dict(cases[3][1])
    mutated_key = (0, 0, 0, 0, 10)
    mixed[mutated_key] = Rat(fmpz_poly([1, 1]), mixed[mutated_key].denominator)
    if stream_column_sha256(mixed) == digests["mixed_signs_sizes_exponents"]:
        raise AssertionError("a coefficient mutation left the streamed digest unchanged")
    rules = [((0, 0, 0, 0, 9), cases[3][1]), ((0, 0, 0, 5, 0), {}), ((2, 0, 0, 0, 0), cases[2][1])]
    if stream_rule_sha256(rules) != reference_rule_sha256(rules):
        raise AssertionError("streaming rule digest differs from the reference")
    sample = guard.sample_process(os.getpid())
    for key in ("resident_bytes", "footprint_bytes"):
        if type(sample.get(key)) is not int or sample[key] <= 0:
            raise AssertionError(f"resource guard sample lacks a positive {key}")
    final = memory_stage("smoke complete", started, wall)
    print(
        "TRACE_NINE_BOUNDED_CANARY_SMOKE PASS "
        f"cases={len(cases)} rule_lists=1 coefficient_mutations=1 "
        f"current_resident={final['current_resident_bytes']} "
        f"current_footprint={final['current_footprint_bytes']} "
        f"peak_rss={final['peak_rss_bytes']}",
        flush=True,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="streaming-digest and guard-sample checks only; no lift, no artifact",
    )
    args = parser.parse_args(argv)
    if args.smoke:
        return smoke()
    return verify(args.artifact.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
