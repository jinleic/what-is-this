#!/usr/bin/env python3
"""Independent verifier for the H682/H684 exact rational-point sign census (e255).

Target: results/spectral/trace_nine_exact_sign_census.json written by
experiments/e255_trace_nine_exact_sign_census.py.  The producer is never
imported and none of its reduction, clearing, or determinant code is reused.

Independent route
  * Provenance: canonical bytes, data digest, and an 18-path source inventory
    recomputed from the working tree (e255, e254, the guard, the e253 and e251
    artifacts, and every e254/e252 source_hashes path).
  * Lift: the pinned Wave-27 QQ(q) lift is rebuilt through the frozen
    verification stack (tests/test_trace_nine_bounded_canary.py and its
    Wave-27/e251 verifier base) with the transcript, term-count, rule-digest,
    and F4..F8 -> 0 checks that H681 passed.
  * Specialization: every rule/F9 coefficient is evaluated at each census
    point by integer Horner into reduced Fractions (no fmpq), and every
    multiplication column is normal-formed by a separately written reducer
    that picks the LAST divisible leader (the producer picks the first).
    proofs/trace_nine_norm_envelope.md section 3 makes the lifted rules a
    Groebner basis, so both strategies must agree wherever no denominator
    vanishes.
  * Determinant: at the exact point(s) the column-cleared integer matrix is
    digested and its determinant recomputed by a Python-driven fraction-free
    Bareiss elimination on flint fmpz scalars (not fmpz_mat.det); the reduced
    rational must equal the recorded decimal strings digit for digit.
  * Modular cross-check: at all six points the specialized matrix is reduced
    modulo 2147483647 and 2147483629 and its determinant is computed by
    Gaussian elimination mod p; both recorded residues must match.
  * Landed witnesses: the four pinned residues and two pinned step totals are
    read from the e253 and e251 artifacts; the e253 reduced modular norm is
    Horner-evaluated at q mod 2147483647 and must agree at every point.
  * Denominator atoms: every e252 physical_lift.denominator_atoms polynomial
    is evaluated exactly at all six points and must be positive.
  * Mutations: nine claim-falsifying edits of the artifact must be rejected.

Scope of a PASS: the recorded exact rational determinants, signs, residues,
and bracket bookkeeping at six rational points of the fixed rank-96 quotient.
No polynomial reconstruction, root isolation, multiplicity, emptiness, shifted
H0/H1 branch, height, or thermodynamic claim is verified or implied.
Read-only: this script never writes a file.  Live memory and timings go to
stdout only.
"""
from __future__ import annotations

import argparse
import copy
import gc
import hashlib
import importlib.util
import json
import math
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from flint import fmpz, fmpz_mat
from sympy.core.cache import clear_cache

sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent.parent
SELF = ROOT / "tests" / "test_trace_nine_exact_sign_census.py"
PRODUCER = ROOT / "experiments" / "e255_trace_nine_exact_sign_census.py"
CANARY_PRODUCER = ROOT / "experiments" / "e254_trace_nine_bounded_canary.py"
CANARY_VERIFIER = ROOT / "tests" / "test_trace_nine_bounded_canary.py"
GUARD = ROOT / "tools" / "resource_guard.py"
ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_exact_sign_census.json"
TEMPLATE = ROOT / "results" / "spectral" / "trace_nine_lift_template.json.xz"
WAVE27_PRODUCER = ROOT / "experiments" / "e252_trace_nine_norm_envelope.py"
WAVE27_VERIFIER = ROOT / "tests" / "test_trace_nine_norm_envelope.py"
WAVE27_PROOF = ROOT / "proofs" / "trace_nine_norm_envelope.md"
WAVE27_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_norm_envelope.json"
E251_PRODUCER = ROOT / "experiments" / "e251_trace_nine_projection.py"
E251_VERIFIER = ROOT / "tests" / "test_trace_nine_projection.py"
E251_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_projection.json"
E253_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_structured_determinant.json"
E248_PRODUCER = ROOT / "experiments" / "e248_replica_trace_nine.py"
E248_VERIFIER = ROOT / "tests" / "test_replica_trace_nine.py"
TRACE_ARTIFACT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
LOCKFILE = ROOT / "uv.lock"
RUNTIME_FREEZE = ROOT.parent / "requirements-freeze.txt"
# e252.source_hashes (9) + e254 additions (7) + e255 additions (e255 script,
# e253 artifact; the e251 artifact is already present) = 18 distinct paths.
SOURCE_PATHS = (
    PRODUCER,
    CANARY_PRODUCER,
    CANARY_VERIFIER,
    GUARD,
    TEMPLATE,
    WAVE27_PRODUCER,
    WAVE27_VERIFIER,
    WAVE27_PROOF,
    WAVE27_ARTIFACT,
    E251_PRODUCER,
    E251_VERIFIER,
    E251_ARTIFACT,
    E253_ARTIFACT,
    E248_PRODUCER,
    E248_VERIFIER,
    TRACE_ARTIFACT,
    LOCKFILE,
    RUNTIME_FREEZE,
)

PREFIX = "TRACE_NINE_EXACT_SIGN_CENSUS_VERIFIER"
SCOPE = "EXACT_RATIONAL_POINT_SIGN_CENSUS_ONLY"
QUOTIENT_RANK = 96
RULE_COUNT = 35
PURE_POWER_BOUNDS = [2, 2, 3, 5, 9]
HILBERT_VECTOR = [1, 5, 12, 19, 22, 19, 12, 5, 1]
RULE_SHA256 = "f4cfe038d7625c6a411f2e01b9af12d9ae9ca06d980e64ad38d8eff030fe22cc"
PRIMARY_PRIME = 2147483647
SECONDARY_PRIME = 2147483629
POINT_LABELS = ("5/4", "3/2", "5/3", "2", "3", "5")
POINT_VALUES = tuple(Fraction(label) for label in POINT_LABELS)
# Planner-pinned landed witnesses; the same values are re-read from the e253
# and e251 artifacts at run time and both sources must agree.
PINNED_RESIDUES = {
    ("2", PRIMARY_PRIME): 1660951362,
    ("3", PRIMARY_PRIME): 2143576240,
    ("5/3", PRIMARY_PRIME): 843620109,
    ("5/3", SECONDARY_PRIME): 1100728375,
}
PINNED_STEPS = {"2": 26131, "5/3": 26131}
E253_GENERIC_STEPS = 26131
E253_GENERIC_NONZERO = 8557
NUMERATOR_SHA256 = "45834f272f6e171ae6072008bfd55556d2b3f39f0621c2f57684872999e61686"
DENOMINATOR_SHA256 = "1b1bfabcff271a3ac3ec937a509b66f483f2c1b53fc390e28f99a6a68ea4a77e"
NUMERATOR_DEGREE = 19846
DENOMINATOR_DEGREE = 16048
ROOT_STATUS = {
    "zero": "EXACT_ROOT_AT_SAMPLED_POINT",
    "bracket": "REAL_ROOT_BRACKET_PROVED_BY_SIGN_CHANGE",
    "none": "NO_SIGN_CHANGE_AMONG_SAMPLED_POINTS_NOT_EVIDENCE_OF_EMPTINESS",
}
META_KEYS = frozenset({"schema_version", "source_sha256", "data_sha256"})
DATA_KEYS = frozenset(
    {
        "scope",
        "quotient_rank",
        "final_rational_rule_sha256",
        "point_count",
        "points",
        "summary",
        "physical_root_status",
        "full_polynomial_status",
        "interpretation",
    }
)
SUMMARY_KEYS = frozenset(
    {
        "signs",
        "zero_count",
        "sign_change_brackets",
        "reduced_modular_norm_agreement_count",
        "landed_residue_matches",
    }
)
POINT_KEYS = frozenset(
    {
        "q",
        "q_numerator",
        "q_denominator",
        "ordered_reduction_steps_total",
        "ordered_reduction_steps_maximum",
        "matrix_nonzero_entries",
        "column_lcm_bits_maximum",
        "column_lcm_product_bits",
        "cleared_integer_matrix_sha256",
        "cleared_integer_determinant_sign",
        "cleared_integer_determinant_bits",
        "determinant_sign",
        "determinant_numerator_bits",
        "determinant_denominator_bits",
        "determinant_numerator_decimal",
        "determinant_denominator_decimal",
        "determinant_sha256",
        "residue_mod_2147483647",
        "residue_mod_2147483629",
        "reduced_modular_norm_residue_mod_2147483647",
        "matches_reduced_modular_norm",
    }
)
EXACT_COMPARED_KEYS = (
    "matrix_nonzero_entries",
    "column_lcm_bits_maximum",
    "column_lcm_product_bits",
    "cleared_integer_matrix_sha256",
    "cleared_integer_determinant_sign",
    "cleared_integer_determinant_bits",
    "determinant_sign",
    "determinant_numerator_bits",
    "determinant_denominator_bits",
    "determinant_numerator_decimal",
    "determinant_denominator_decimal",
    "determinant_sha256",
    "residue_mod_2147483647",
    "residue_mod_2147483629",
)
MUTATION_IDS = [
    "M1_sign_flipped_at_first_point",
    "M2_primary_residue_plus_one_at_second_point",
    "M3_step_total_plus_one_at_five_thirds",
    "M4_bracket_removed_or_spurious_bracket_inserted",
    "M5_producer_source_hash_last_nibble",
    "M6_physical_root_status_swapped",
    "M7_exact_numerator_negated_at_q_two",
    "M8_full_polynomial_status_escalated",
    "M9_secondary_residue_plus_one_at_five_thirds",
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


# Frozen verification stack only: H681 canary verifier -> Wave-27 verifier ->
# e251 verifier.  No producer module is loaded.
canary = load_module("sign_census_canary_verifier", CANARY_VERIFIER)
wave27 = canary.wave27
e251 = canary.e251
guard = canary.guard
Rat = canary.Rat


def fail(message: str) -> None:
    raise AssertionError(message)


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


def stage(label: str, started: float, wall: float) -> dict[str, object]:
    gc.collect()
    sample = guard.sample_process(os.getpid())
    record = {
        "stage": label,
        "cpu_seconds": round(time.process_time() - started, 3),
        "wall_seconds": round(time.monotonic() - wall, 3),
        "current_resident_bytes": int(sample["resident_bytes"]),
        "current_footprint_bytes": int(sample["footprint_bytes"]),
        "peak_rss_bytes": wave27.peak_rss_bytes(),
    }
    print(PREFIX + " STAGE " + json.dumps(record, sort_keys=True, separators=(",", ":")), flush=True)
    return record


def is_int(value: object) -> bool:
    return type(value) is int


def sign_of(value) -> int:
    return (value > 0) - (value < 0)


def residue(value: Fraction, prime: int) -> int | None:
    if value.denominator % prime == 0:
        return None
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def horner_mod(coefficients: list[int], point: int, prime: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * point + coefficient) % prime
    return value


# ---------------------------------------------------------------------------
# Exact scalar specialization and reduction (independent of e255.reduce_point)
# ---------------------------------------------------------------------------

_KEYS: dict[Monomial, tuple[int, tuple[int, ...]]] = {}


def grevlex_key(monomial: Monomial) -> tuple[int, tuple[int, ...]]:
    key = _KEYS.get(monomial)
    if key is None:
        key = _KEYS[monomial] = (sum(monomial), tuple(-value for value in reversed(monomial)))
    return key


def add(left: Monomial, right: Monomial) -> Monomial:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2],
            left[3] + right[3], left[4] + right[4])


def scaled_polynomial_value(coefficients: list[int], a: int, bpows: list[int]) -> int:
    """Return b^d * P(a/b) for ascending integer coefficients (0 for the zero polynomial)."""
    degree = len(coefficients) - 1
    accumulator = 0
    for index in range(degree, -1, -1):
        accumulator = accumulator * a + coefficients[index] * bpows[degree - index]
    return accumulator


def specialize_rat(value, a: int, bpows: list[int]) -> Fraction:
    numerator = [int(c) for c in value.numerator]
    denominator = [int(c) for c in value.denominator]
    if not numerator:
        return Fraction(0)
    top = scaled_polynomial_value(numerator, a, bpows)
    bottom = scaled_polynomial_value(denominator, a, bpows)
    if bottom == 0:
        fail("a lifted coefficient denominator vanishes at a census point")
    return Fraction(top * bpows[len(denominator) - 1], bottom * bpows[len(numerator) - 1])


def specialize_inputs(rules, trace_row, point: Fraction):
    a, b = point.numerator, point.denominator
    maximum_degree = 0
    for _, tail in rules:
        for value in tail.values():
            maximum_degree = max(maximum_degree, value.numerator.degree(), value.denominator.degree())
    for value in trace_row.values():
        maximum_degree = max(maximum_degree, value.numerator.degree(), value.denominator.degree())
    bpows = [b ** k for k in range(maximum_degree + 1)]
    rules_point = []
    for leading, tail in rules:
        specialized = {}
        for monomial, value in tail.items():
            scalar = specialize_rat(value, a, bpows)
            if scalar:
                specialized[monomial] = scalar
        rules_point.append((leading, specialized))
    trace_point = {}
    for monomial, value in trace_row.items():
        scalar = specialize_rat(value, a, bpows)
        if scalar:
            trace_point[monomial] = scalar
    return rules_point, trace_point


def to_mod(value: Fraction, prime: int) -> int:
    if value.denominator % prime == 0:
        fail(f"mod-{prime} route unavailable: a specialized coefficient denominator vanishes")
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def reduce_inputs_mod(rules_point, trace_point, prime: int):
    rules_mod = []
    for leading, tail in rules_point:
        specialized = {}
        for monomial, value in tail.items():
            scalar = to_mod(value, prime)
            if scalar:
                specialized[monomial] = scalar
        rules_mod.append((leading, specialized))
    trace_mod = {}
    for monomial, value in trace_point.items():
        scalar = to_mod(value, prime)
        if scalar:
            trace_mod[monomial] = scalar
    return rules_mod, trace_mod


def normal_form(polynomial, rules, leaders, lookup, normalize):
    """Largest-grevlex-first reduction choosing the LAST divisible leader.

    Any strategy yields the same remainder for a Groebner basis; the producer
    uses a priority heap and the first divisible leader, so agreement is a
    genuine cross-check of both the columns and the Groebner claim.
    """
    work = {monomial: coefficient for monomial, coefficient in polynomial.items() if coefficient}
    remainder = {}
    steps = 0
    while work:
        monomial = max(work, key=grevlex_key)
        coefficient = work.pop(monomial)
        index = lookup.get(monomial, -1)
        if index == -1:
            index = None
            for candidate, leader in enumerate(leaders):
                if (leader[0] <= monomial[0] and leader[1] <= monomial[1] and leader[2] <= monomial[2]
                        and leader[3] <= monomial[3] and leader[4] <= monomial[4]):
                    index = candidate
            lookup[monomial] = index
        if index is None:
            updated = normalize(remainder.get(monomial, 0) + coefficient)
            if updated:
                remainder[monomial] = updated
            else:
                remainder.pop(monomial, None)
            continue
        leader, tail = rules[index]
        shift = (monomial[0] - leader[0], monomial[1] - leader[1], monomial[2] - leader[2],
                 monomial[3] - leader[3], monomial[4] - leader[4])
        for tail_monomial, tail_coefficient in tail.items():
            output = add(shift, tail_monomial)
            updated = normalize(work.get(output, 0) - coefficient * tail_coefficient)
            if updated:
                work[output] = updated
            else:
                work.pop(output, None)
        steps += 1
    return remainder, steps


def identity(value):
    return value


def census_columns(rules, trace, standard_list, normalize):
    leaders = [leading for leading, _ in rules]
    row_index = {monomial: index for index, monomial in enumerate(standard_list)}
    lookup: dict[Monomial, int | None] = {}
    columns = []
    total = 0
    maximum = 0
    for basis in standard_list:
        initial = {add(term, basis): value for term, value in trace.items()}
        remainder, steps = normal_form(initial, rules, leaders, lookup, normalize)
        column = {}
        for monomial, value in remainder.items():
            if monomial not in row_index:
                fail("reduced column contains a nonstandard monomial")
            column[row_index[monomial]] = value
        columns.append(column)
        total += steps
        maximum = max(maximum, steps)
    return columns, total, maximum


# ---------------------------------------------------------------------------
# Determinants: fraction-free Bareiss on fmpz scalars, and Gaussian mod p
# ---------------------------------------------------------------------------

def bareiss_determinant(rows: list[list[int]]) -> int:
    size = len(rows)
    if size == 0:
        return 1
    matrix = [[fmpz(entry) for entry in row] for row in rows]
    sign = 1
    previous = fmpz(1)
    zero = fmpz(0)
    for k in range(size - 1):
        if matrix[k][k] == zero:
            swap = None
            for i in range(k + 1, size):
                if matrix[i][k] != zero:
                    swap = i
                    break
            if swap is None:
                return 0
            matrix[k], matrix[swap] = matrix[swap], matrix[k]
            sign = -sign
        pivot_row = matrix[k]
        pivot = pivot_row[k]
        for i in range(k + 1, size):
            row = matrix[i]
            factor = row[k]
            if factor == zero:
                for j in range(k + 1, size):
                    row[j] = (row[j] * pivot) // previous
            else:
                for j in range(k + 1, size):
                    row[j] = (row[j] * pivot - factor * pivot_row[j]) // previous
            row[k] = zero
        previous = pivot
        matrix[k] = None  # release processed row
    return sign * int(matrix[size - 1][size - 1])


def determinant_mod(rows: list[list[int]], prime: int) -> int:
    size = len(rows)
    matrix = [[entry % prime for entry in row] for row in rows]
    determinant = 1
    for k in range(size):
        pivot_index = None
        for i in range(k, size):
            if matrix[i][k]:
                pivot_index = i
                break
        if pivot_index is None:
            return 0
        if pivot_index != k:
            matrix[k], matrix[pivot_index] = matrix[pivot_index], matrix[k]
            determinant = -determinant
        pivot_row = matrix[k]
        pivot = pivot_row[k]
        determinant = determinant * pivot % prime
        inverse = pow(pivot, -1, prime)
        for i in range(k + 1, size):
            row = matrix[i]
            factor = row[k] * inverse % prime
            if factor:
                for j in range(k, size):
                    row[j] = (row[j] - factor * pivot_row[j]) % prime
    return determinant % prime


def dense_from_columns(columns, size: int) -> list[list[int]]:
    dense = [[0] * size for _ in range(size)]
    for column_index, column in enumerate(columns):
        for row, value in column.items():
            dense[row][column_index] = value
    return dense


def exact_point_record(label: str, rules_point, trace_point, standard_list) -> dict[str, object]:
    columns, my_total, my_maximum = census_columns(rules_point, trace_point, standard_list, identity)
    size = len(standard_list)
    column_lcms = []
    for column in columns:
        scale = 1
        for value in column.values():
            scale = math.lcm(scale, value.denominator)
        column_lcms.append(scale)
    cleared = [[0] * size for _ in range(size)]
    for column_index, column in enumerate(columns):
        scale = column_lcms[column_index]
        for row, value in column.items():
            cleared[row][column_index] = value.numerator * (scale // value.denominator)
    nonzero_entries = sum(1 for row in cleared for entry in row if entry)
    entry_bits_maximum = max(abs(entry).bit_length() for row in cleared for entry in row)
    cleared_digest = canonical_sha256([[str(entry) for entry in row] for row in cleared])
    del columns
    bareiss_started = time.process_time()
    cleared_determinant = bareiss_determinant(cleared)
    bareiss_cpu = time.process_time() - bareiss_started
    del cleared
    scale_product = math.prod(column_lcms)
    determinant = Fraction(cleared_determinant, scale_product)
    record = {
        "q": label,
        "independent_reduction_steps_total": my_total,
        "independent_reduction_steps_maximum": my_maximum,
        "cleared_entry_bits_maximum": entry_bits_maximum,
        "bareiss_cpu_seconds": round(bareiss_cpu, 3),
        "matrix_nonzero_entries": nonzero_entries,
        "column_lcm_bits_maximum": max(value.bit_length() for value in column_lcms),
        "column_lcm_product_bits": scale_product.bit_length(),
        "cleared_integer_matrix_sha256": cleared_digest,
        "cleared_integer_determinant_sign": sign_of(cleared_determinant),
        "cleared_integer_determinant_bits": abs(cleared_determinant).bit_length(),
        "determinant_sign": sign_of(determinant),
        "determinant_numerator_bits": abs(determinant.numerator).bit_length(),
        "determinant_denominator_bits": determinant.denominator.bit_length(),
        "determinant_numerator_decimal": str(determinant.numerator),
        "determinant_denominator_decimal": str(determinant.denominator),
        "determinant_sha256": hashlib.sha256(
            f"{determinant.numerator}/{determinant.denominator}".encode("ascii")
        ).hexdigest(),
        "residue_mod_2147483647": residue(determinant, PRIMARY_PRIME),
        "residue_mod_2147483629": residue(determinant, SECONDARY_PRIME),
    }
    print(PREFIX + " EXACT_POINT " + json.dumps(
        {key: value for key, value in record.items() if not key.endswith("_decimal")},
        sort_keys=True, separators=(",", ":")), flush=True)
    return record


# ---------------------------------------------------------------------------
# Inherited artifacts: e253 reduced modular norm, e251 witnesses, e252 atoms
# ---------------------------------------------------------------------------

def load_e253() -> tuple[list[int], list[int], dict[tuple[str, int], int]]:
    artifact = json.loads(E253_ARTIFACT.read_text())
    data = artifact["data"]
    if artifact["meta"]["data_sha256"] != canonical_sha256(data):
        fail("e253 artifact data digest mismatch")
    if data["field"] != {"parameter": "q", "prime": PRIMARY_PRIME}:
        fail("e253 field is not F_2147483647(q)")
    if data["matrix"]["ordered_reduction_steps"] != E253_GENERIC_STEPS:
        fail("e253 generic ordered reduction count changed")
    if data["matrix"]["exact_nonzero_entries"] != E253_GENERIC_NONZERO:
        fail("e253 generic nonzero entry count changed")
    norm = data["reduced_modular_norm"]
    numerator = norm["numerator_coefficients_ascending"]
    denominator = norm["denominator_coefficients_ascending"]
    for name, coefficients, degree, digest in (
        ("numerator", numerator, NUMERATOR_DEGREE, NUMERATOR_SHA256),
        ("denominator", denominator, DENOMINATOR_DEGREE, DENOMINATOR_SHA256),
    ):
        if not isinstance(coefficients, list) or len(coefficients) != degree + 1:
            fail(f"e253 reduced modular norm {name} degree is not {degree}")
        if any(not is_int(value) or abs(value) >= PRIMARY_PRIME for value in coefficients):
            fail(f"e253 reduced modular norm {name} has a non-residue coefficient")
        if coefficients[-1] == 0:
            fail(f"e253 reduced modular norm {name} leading coefficient vanishes")
        if norm[f"{name}_degree"] != degree:
            fail(f"e253 recorded {name} degree mismatch")
        if canonical_sha256(coefficients) != digest or norm[f"{name}_coefficients_sha256"] != digest:
            fail(f"e253 reduced modular norm {name} digest mismatch")
    pinned: dict[tuple[str, int], int] = {}
    for check in data["direct_point_checks"]:
        q_value = Fraction(check["q"])
        if residue(q_value, PRIMARY_PRIME) != check["q_residue"]:
            fail("e253 direct point check q_residue inconsistent")
        if not is_int(check["determinant_residue"]):
            fail("e253 direct point residue is not an integer")
        pinned[(check["q"], PRIMARY_PRIME)] = check["determinant_residue"]
    return numerator, denominator, pinned


def load_e251() -> tuple[dict[tuple[str, int], int], dict[str, int]]:
    artifact = json.loads(E251_ARTIFACT.read_text())
    data = artifact["data"]
    if artifact["meta"]["data_sha256"] != canonical_sha256(data):
        fail("e251 artifact data digest mismatch")
    if not all(check["passed"] for check in artifact["checks"]):
        fail("e251 artifact contains a failed check")
    residues: dict[tuple[str, int], int] = {}
    steps: dict[str, int] = {}
    for witness in data["modular_certificates"]["physical_witnesses"]:
        prime = witness["prime"]
        q_value = Fraction(witness["q"])
        if residue(q_value, prime) != witness["q_residue"]:
            fail("e251 witness q_residue inconsistent")
        if witness["matrix_rank"] != QUOTIENT_RANK or witness["quotient_basis_rank"] != QUOTIENT_RANK:
            fail("e251 witness rank is not 96")
        residues[(witness["q"], prime)] = witness["determinant_residue"]
        total = witness["reduction_steps"]["sum"]
        if steps.get(witness["q"], total) != total:
            fail("e251 witnesses disagree on the reduction step total")
        steps[witness["q"]] = total
    return residues, steps


def atom_polynomials(wave27_data: dict[str, object]) -> list[list[int]]:
    atoms = wave27_data["physical_lift"]["denominator_atoms"]
    if not isinstance(atoms, list) or not atoms:
        fail("e252 denominator atom list is missing or empty")
    polynomials = []
    for index, atom in enumerate(atoms):
        if atom["id"] != index:
            fail("e252 denominator atom ids are not sequential")
        coefficients = [int(value) for value in atom["coefficients_ascending"]]
        if [str(value) for value in coefficients] != atom["coefficients_ascending"]:
            fail("e252 denominator atom coefficients are not canonical decimals")
        if not coefficients or coefficients[-1] == 0 or atom["degree"] != len(coefficients) - 1:
            fail("e252 denominator atom degree mismatch")
        l1 = sum(abs(value) for value in coefficients)
        if l1 <= 0 or (l1 - 1).bit_length() != atom["l1_log2_bound"]:
            fail("e252 denominator atom l1 bound mismatch")
        polynomials.append(coefficients)
    return polynomials


def atoms_positive(polynomials: list[list[int]]) -> dict[str, list[int]]:
    """Exact sign of every denominator atom at every census point (must be +1)."""
    signs: dict[str, list[int]] = {}
    for label, point in zip(POINT_LABELS, POINT_VALUES):
        a, b = point.numerator, point.denominator
        maximum_degree = max(len(polynomial) - 1 for polynomial in polynomials)
        bpows = [b ** k for k in range(maximum_degree + 1)]
        values = [sign_of(scaled_polynomial_value(polynomial, a, bpows)) for polynomial in polynomials]
        if any(value != 1 for value in values):
            fail(f"a denominator atom is not positive at q={label}")
        signs[label] = values
    return signs


# ---------------------------------------------------------------------------
# Lift rebuild through the frozen verification stack (no producer import)
# ---------------------------------------------------------------------------

def rebuild_lift(started: float, wall: float):
    runtime = wave27.runtime_manifest()
    print(PREFIX + " RUNTIME " + json.dumps(runtime, sort_keys=True, separators=(",", ":")), flush=True)
    inherited = canary.inherited_wave27()
    rows = canary.inherited_trace_rows()
    template, payload_sha256 = canary.load_template()
    if payload_sha256 != inherited["template"]["payload_sha256"]:
        fail("lift-template payload digest differs from the Wave-27 artifact")
    records = template["records"]
    del template
    if len(records) != RULE_COUNT:
        fail("lift-template record count changed")
    leaders = [tuple(record["leading_monomial"]) for record in records]
    if len(set(leaders)) != RULE_COUNT or any(
        len(leader) != 5 or any(not is_int(power) or power < 0 for power in leader) for leader in leaders
    ):
        fail("lift-template leaders are not distinct 5-exponent monomials")
    pure_power_bounds, standard_list = wave27.standard_monomials(leaders)
    hilbert = [sum(sum(monomial) == degree for monomial in standard_list) for degree in range(9)]
    if pure_power_bounds != PURE_POWER_BOUNDS or len(standard_list) != QUOTIENT_RANK or hilbert != HILBERT_VECTOR:
        fail("quotient dimensions changed")
    standard = frozenset(standard_list)
    stage("inherited inputs loaded", started, wall)

    equations, _ = e251.independent_incidence()
    if len(equations) != 6:
        fail("independent incidence must have six equations F4..F9")
    leading_forms = [
        {
            monomial: Fraction(int(sp.numer(coefficient)), int(sp.denom(coefficient)))
            for monomial, coefficient in wave27.homogeneous_leading(equation).terms()
        }
        for equation in equations[:5]
    ]
    degrees = [max(map(sum, form)) for form in leading_forms]
    canary.template_identities(records, leaders, leading_forms, degrees)
    wave27.guard_resources(started, "template identities")
    stage("template identities", started, wall)

    targets = tuple(Rat.from_sympy(value) for value in e251.target_expressions(rows).values())
    physical = canary.evaluate_rows(e251.compile_rows(equations), targets)
    del equations, targets
    for index, form in enumerate(leading_forms):
        degree = degrees[index]
        observed = {monomial: value for monomial, value in physical[index].items() if sum(monomial) == degree}
        expected = {monomial: Rat.scalar(value) for monomial, value in form.items()}
        if observed != expected or any(sum(monomial) > degree for monomial in physical[index]):
            fail(f"physical equation F{index + 4} changed its leading form")
    del leading_forms
    clear_cache()
    wave27.guard_resources(started, "physical targets")
    stage("physical targets; sympy released", started, wall)

    remainders = [
        {monomial: value for monomial, value in physical[index].items() if sum(monomial) < degrees[index]}
        for index in range(5)
    ]
    rules, interreduction_steps = canary.rebuild_rules(records, leaders, remainders, started)
    del records, remainders
    lift = inherited["physical_lift"]
    if interreduction_steps != lift["interreduction_steps_by_record"]:
        fail("interreduction step transcript differs from the Wave-27 artifact")
    if [1 + len(tail) for _, tail in rules] != lift["final_rule_term_counts"]:
        fail("final rule term counts differ from the Wave-27 artifact")
    if any(monomial not in standard for _, tail in rules for monomial in tail):
        fail("interreduced tail contains a nonstandard monomial")
    if [leading for leading, _ in rules] != leaders:
        fail("rule leaders are not in template order")
    rule_sha256 = canary.stream_rule_sha256(rules)
    if rule_sha256 != RULE_SHA256:
        fail(f"final rational rule digest {rule_sha256} is not the pinned digest")
    for index in range(5):
        remainder, _ = wave27.reduce_exact(physical[index], rules)
        if remainder:
            fail(f"physical equation F{index + 4} does not reduce to zero")
    trace_row = physical[5]
    del physical
    if not trace_row:
        fail("F9 row is identically zero")
    wave27.guard_resources(started, "validated lift")
    stage("validated lift; template and F4..F8 released", started, wall)
    return inherited, rules, trace_row, standard_list


# ---------------------------------------------------------------------------
# Artifact validation (pure; all heavy results arrive through ctx)
# ---------------------------------------------------------------------------

def parse_decimal(text: object, name: str) -> int:
    if not isinstance(text, str) or not text or text.lstrip("-") == "" or not text.lstrip("-").isdigit():
        fail(f"{name} is not a decimal integer string")
    value = int(text)
    if str(value) != text:
        fail(f"{name} is not a canonical decimal string")
    return value


def validate_point(index: int, record: object, ctx: dict[str, object]) -> Fraction:
    label = POINT_LABELS[index]
    point = POINT_VALUES[index]
    if not isinstance(record, dict) or set(record) != POINT_KEYS:
        fail(f"point {index} keys mismatch")
    if record["q"] != label:
        fail(f"point {index} is not q={label}")
    if not is_int(record["q_numerator"]) or not is_int(record["q_denominator"]):
        fail(f"point {label} q_numerator/q_denominator are not integers")
    if Fraction(record["q_numerator"], record["q_denominator"]) != point or record["q_denominator"] <= 0:
        fail(f"point {label} q_numerator/q_denominator inconsistent")
    if math.gcd(record["q_numerator"], record["q_denominator"]) != 1 or point <= 1:
        fail(f"point {label} is not a reduced rational greater than one")
    for key in (
        "ordered_reduction_steps_total",
        "ordered_reduction_steps_maximum",
        "matrix_nonzero_entries",
        "column_lcm_bits_maximum",
        "column_lcm_product_bits",
        "cleared_integer_determinant_bits",
        "determinant_numerator_bits",
        "determinant_denominator_bits",
        "residue_mod_2147483647",
        "residue_mod_2147483629",
        "reduced_modular_norm_residue_mod_2147483647",
    ):
        if not is_int(record[key]) or record[key] < 0:
            fail(f"point {label} {key} is not a non-negative integer")
    total = record["ordered_reduction_steps_total"]
    if not 0 < record["ordered_reduction_steps_maximum"] <= total <= ctx["step_bound"]:
        fail(f"point {label} ordered reduction totals violate the generic bound")
    if label in ctx["pinned_steps"] and total != ctx["pinned_steps"][label]:
        fail(f"point {label} ordered reduction total differs from the landed e251 witness")
    if not QUOTIENT_RANK <= record["matrix_nonzero_entries"] <= QUOTIENT_RANK * QUOTIENT_RANK:
        fail(f"point {label} nonzero entry count outside 96..9216")
    if record["column_lcm_bits_maximum"] < 1 or record["column_lcm_product_bits"] < record["column_lcm_bits_maximum"]:
        fail(f"point {label} column LCM bit lengths inconsistent")
    numerator = parse_decimal(record["determinant_numerator_decimal"], f"point {label} numerator")
    denominator = parse_decimal(record["determinant_denominator_decimal"], f"point {label} denominator")
    if denominator <= 0 or math.gcd(numerator, denominator) != 1:
        fail(f"point {label} determinant is not a reduced fraction with positive denominator")
    determinant = Fraction(numerator, denominator)
    if not is_int(record["determinant_sign"]) or not is_int(record["cleared_integer_determinant_sign"]):
        fail(f"point {label} determinant signs are not integers")
    if record["determinant_sign"] != sign_of(numerator):
        fail(f"point {label} determinant_sign disagrees with the decimal numerator")
    if record["cleared_integer_determinant_sign"] != record["determinant_sign"]:
        fail(f"point {label} cleared determinant sign disagrees with the rational sign")
    if record["determinant_numerator_bits"] != abs(numerator).bit_length():
        fail(f"point {label} numerator bit length mismatch")
    if record["determinant_denominator_bits"] != denominator.bit_length():
        fail(f"point {label} denominator bit length mismatch")
    if numerator == 0:
        if record["cleared_integer_determinant_bits"] != 0:
            fail(f"point {label} zero determinant with nonzero cleared bits")
    else:
        drift = (record["cleared_integer_determinant_bits"] - record["column_lcm_product_bits"]
                 - record["determinant_numerator_bits"] + record["determinant_denominator_bits"])
        if record["cleared_integer_determinant_bits"] < record["determinant_numerator_bits"] or abs(drift) > 1:
            fail(f"point {label} cleared determinant bit length inconsistent with the LCM product")
    expected_digest = hashlib.sha256(f"{numerator}/{denominator}".encode("ascii")).hexdigest()
    if record["determinant_sha256"] != expected_digest:
        fail(f"point {label} determinant digest mismatch")
    if record["cleared_integer_matrix_sha256"] is None or not isinstance(record["cleared_integer_matrix_sha256"], str) \
            or len(record["cleared_integer_matrix_sha256"]) != 64 \
            or not set(record["cleared_integer_matrix_sha256"]) <= HEX_DIGITS:
        fail(f"point {label} cleared matrix digest is not lowercase sha256 hex")
    for prime, key in ((PRIMARY_PRIME, "residue_mod_2147483647"), (SECONDARY_PRIME, "residue_mod_2147483629")):
        expected = residue(determinant, prime)
        if expected is None or record[key] != expected:
            fail(f"point {label} {key} disagrees with the decimal determinant")
        modular_determinant, modular_nonzero = ctx["modp"][(label, prime)]
        if modular_determinant != expected:
            fail(f"point {label} {key} disagrees with the independent mod-{prime} determinant")
        if modular_nonzero > record["matrix_nonzero_entries"]:
            fail(f"point {label} records fewer nonzero entries than the mod-{prime} matrix")
        pinned = ctx["pinned_residues"].get((label, prime))
        if pinned is not None and record[key] != pinned:
            fail(f"point {label} {key} misses the landed residue {pinned}")
    norm_residue = ctx["norm_residue"][label]
    if record["reduced_modular_norm_residue_mod_2147483647"] != norm_residue:
        fail(f"point {label} reduced modular norm residue differs from the independent Horner value")
    if record["matches_reduced_modular_norm"] is not True or record["residue_mod_2147483647"] != norm_residue:
        fail(f"point {label} exact determinant disagrees with the e253 reduced modular norm")
    exact = ctx["exact"].get(label)
    if exact is not None:
        for key in EXACT_COMPARED_KEYS:
            if record[key] != exact[key]:
                fail(f"point {label} {key} differs from the independent exact recomputation")
    return determinant


def validate_artifact(artifact: object, sources: dict[str, str], ctx: dict[str, object]) -> dict[str, object]:
    if not isinstance(artifact, dict) or set(artifact) != {"data", "meta"}:
        fail("artifact must contain exactly data and meta")
    meta = artifact["meta"]
    data = artifact["data"]
    if not isinstance(meta, dict) or set(meta) != META_KEYS:
        fail("meta keys mismatch")
    if not is_int(meta["schema_version"]) or meta["schema_version"] != 1:
        fail("unexpected schema version")
    if not isinstance(data, dict) or set(data) != DATA_KEYS:
        fail("data keys mismatch")
    if meta["data_sha256"] != canonical_sha256(data):
        fail("artifact data digest mismatch")
    if meta["source_sha256"] != sources:
        fail("source provenance inventory or hash mismatch")
    if data["scope"] != SCOPE:
        fail("artifact scope is not the exact rational-point sign census")
    if data["full_polynomial_status"] != "NOT_RECONSTRUCTED":
        fail("artifact claims a polynomial status beyond the census")
    if not is_int(data["quotient_rank"]) or data["quotient_rank"] != QUOTIENT_RANK:
        fail("quotient rank is not 96")
    if data["final_rational_rule_sha256"] != RULE_SHA256:
        fail("final rational rule digest is not the pinned digest")
    if not isinstance(data["interpretation"], str) or not data["interpretation"]:
        fail("interpretation is not a nonempty string")
    points = data["points"]
    if not isinstance(points, list) or len(points) != len(POINT_LABELS):
        fail("artifact must carry exactly six point records")
    if not is_int(data["point_count"]) or data["point_count"] != len(points):
        fail("point_count disagrees with the point list")
    determinants = [validate_point(index, record, ctx) for index, record in enumerate(points)]
    signs = [sign_of(value) for value in determinants]
    zero_count = sum(1 for value in signs if value == 0)
    brackets = [
        {"lower": POINT_LABELS[index], "upper": POINT_LABELS[index + 1]}
        for index in range(len(signs) - 1)
        if signs[index] * signs[index + 1] < 0
    ]
    summary = data["summary"]
    if not isinstance(summary, dict) or set(summary) != SUMMARY_KEYS:
        fail("summary keys mismatch")
    if summary["signs"] != signs or summary["signs"] != [record["determinant_sign"] for record in points]:
        fail("summary signs differ from the recorded determinants")
    if not is_int(summary["zero_count"]) or summary["zero_count"] != zero_count:
        fail("summary zero_count mismatch")
    if summary["sign_change_brackets"] != brackets:
        fail("summary sign_change_brackets differ from the recomputed brackets")
    if not is_int(summary["reduced_modular_norm_agreement_count"]) or \
            summary["reduced_modular_norm_agreement_count"] != len(points):
        fail("reduced modular norm agreement count is not six")
    if not is_int(summary["landed_residue_matches"]) or \
            summary["landed_residue_matches"] != len(ctx["pinned_residues"]):
        fail("landed residue match count differs from the pinned witness count")
    expected_status = ROOT_STATUS["zero"] if zero_count else ROOT_STATUS["bracket"] if brackets else ROOT_STATUS["none"]
    if data["physical_root_status"] != expected_status:
        fail("physical_root_status violates the three-way rule")
    return {"signs": signs, "zero_count": zero_count, "brackets": brackets}


def rejects(artifact: dict[str, object], sources: dict[str, str], ctx: dict[str, object]) -> bool:
    try:
        validate_artifact(artifact, sources, ctx)
    except (AssertionError, KeyError, TypeError, ValueError):
        return True
    return False


def flip_last_nibble(container: dict, key: str) -> None:
    digest = container[key]
    container[key] = digest[:-1] + ("0" if digest[-1] != "0" else "1")


def mutation_rejections(artifact: dict[str, object], sources: dict[str, str], ctx: dict[str, object]) -> list[str]:
    producer_key = str(PRODUCER.relative_to(WORKSPACE))
    outcomes = []

    def attempt(label: str, mutate, reseal: bool = True) -> None:
        mutated = copy.deepcopy(artifact)
        mutate(mutated)
        if reseal:
            mutated["meta"]["data_sha256"] = canonical_sha256(mutated["data"])
        outcomes.append(label if rejects(mutated, sources, ctx) else "")

    def flip_sign(mutated: dict) -> None:
        record = mutated["data"]["points"][0]
        record["determinant_sign"] = -record["determinant_sign"] if record["determinant_sign"] else 1

    def bump_primary(mutated: dict) -> None:
        record = mutated["data"]["points"][1]
        record["residue_mod_2147483647"] = (record["residue_mod_2147483647"] + 1) % PRIMARY_PRIME

    def bump_steps(mutated: dict) -> None:
        mutated["data"]["points"][2]["ordered_reduction_steps_total"] += 1

    def alter_brackets(mutated: dict) -> None:
        brackets = mutated["data"]["summary"]["sign_change_brackets"]
        if brackets:
            brackets.pop()
        else:
            brackets.append({"lower": POINT_LABELS[0], "upper": POINT_LABELS[1]})

    def swap_status(mutated: dict) -> None:
        current = mutated["data"]["physical_root_status"]
        mutated["data"]["physical_root_status"] = (
            ROOT_STATUS["none"] if current == ROOT_STATUS["bracket"] else ROOT_STATUS["bracket"]
        )

    def negate_numerator(mutated: dict) -> None:
        record = mutated["data"]["points"][3]
        text = record["determinant_numerator_decimal"]
        record["determinant_numerator_decimal"] = text[1:] if text.startswith("-") else "-" + text

    def escalate(mutated: dict) -> None:
        mutated["data"]["full_polynomial_status"] = "RECONSTRUCTED"

    def bump_secondary(mutated: dict) -> None:
        record = mutated["data"]["points"][2]
        record["residue_mod_2147483629"] = (record["residue_mod_2147483629"] + 1) % SECONDARY_PRIME

    attempt(MUTATION_IDS[0], flip_sign)
    attempt(MUTATION_IDS[1], bump_primary)
    attempt(MUTATION_IDS[2], bump_steps)
    attempt(MUTATION_IDS[3], alter_brackets)
    attempt(MUTATION_IDS[4], lambda m: flip_last_nibble(m["meta"]["source_sha256"], producer_key), reseal=False)
    attempt(MUTATION_IDS[5], swap_status)
    attempt(MUTATION_IDS[6], negate_numerator)
    attempt(MUTATION_IDS[7], escalate)
    attempt(MUTATION_IDS[8], bump_secondary)
    return outcomes


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def verify(artifact_path: Path, exact_labels: tuple[str, ...]) -> int:
    started = time.process_time()
    wall = time.monotonic()
    artifact_bytes = artifact_path.read_bytes()
    artifact = json.loads(artifact_bytes)
    if artifact_bytes != canonical_bytes(artifact):
        fail("census artifact is not canonical JSON")
    sources = expected_sources()
    if len(sources) != len(SOURCE_PATHS):
        fail("source inventory has duplicate keys")
    if not isinstance(artifact, dict) or artifact.get("meta", {}).get("source_sha256") != sources:
        fail("source provenance inventory or hash mismatch")
    if artifact["meta"].get("data_sha256") != canonical_sha256(artifact["data"]):
        fail("artifact data digest mismatch")
    stage("artifact bytes and provenance validated", started, wall)

    numerator, denominator, pinned_e253 = load_e253()
    residues_e251, steps_e251 = load_e251()
    pinned_residues = dict(pinned_e253)
    for key, value in residues_e251.items():
        if pinned_residues.get(key, value) != value:
            fail("e251 and e253 disagree on a landed residue")
        pinned_residues[key] = value
    if pinned_residues != PINNED_RESIDUES:
        fail(f"landed residues read from the artifacts differ from the planner pins: {pinned_residues}")
    pinned_steps = {label: total for label, total in steps_e251.items() if label in POINT_LABELS}
    if pinned_steps != PINNED_STEPS:
        fail(f"landed step totals read from e251 differ from the planner pins: {pinned_steps}")
    norm_residue: dict[str, int] = {}
    for label, point in zip(POINT_LABELS, POINT_VALUES):
        q_residue = residue(point, PRIMARY_PRIME)
        if q_residue is None:
            fail("census point is a pole modulo the primary prime")
        bottom = horner_mod(denominator, q_residue, PRIMARY_PRIME)
        if bottom == 0:
            fail(f"e253 reduced modular denominator vanishes at q={label}")
        norm_residue[label] = horner_mod(numerator, q_residue, PRIMARY_PRIME) * pow(bottom, -1, PRIMARY_PRIME) % PRIMARY_PRIME
    del numerator, denominator
    stage("e253 norm residues and e251/e253 pins loaded", started, wall)

    inherited, rules, trace_row, standard_list = rebuild_lift(started, wall)
    atoms = atom_polynomials(inherited)
    atom_signs = atoms_positive(atoms)
    step_bound = inherited["matrix_envelope"]["column_reduction_steps"]
    if not isinstance(step_bound, list) or len(step_bound) != QUOTIENT_RANK or any(not is_int(v) or v < 0 for v in step_bound):
        fail("Wave-27 column reduction step envelope malformed")
    step_bound = sum(step_bound)
    if step_bound < max(PINNED_STEPS.values()):
        fail("Wave-27 step envelope is below the landed witness count")
    del inherited
    print(PREFIX + " ATOMS " + json.dumps({"count": len(atoms), "signs": atom_signs, "step_bound": step_bound},
                                          sort_keys=True, separators=(",", ":")), flush=True)
    stage("denominator atoms positive at all six points", started, wall)

    exact: dict[str, dict[str, object]] = {}
    modp: dict[tuple[str, int], tuple[int, int]] = {}
    for label, point in zip(POINT_LABELS, POINT_VALUES):
        rules_point, trace_point = specialize_inputs(rules, trace_row, point)
        for prime in (PRIMARY_PRIME, SECONDARY_PRIME):
            rules_mod, trace_mod = reduce_inputs_mod(rules_point, trace_point, prime)
            columns, total, maximum = census_columns(rules_mod, trace_mod, standard_list, lambda v, p=prime: v % p)
            dense = dense_from_columns(columns, QUOTIENT_RANK)
            nonzero = sum(1 for row in dense for entry in row if entry)
            modp[(label, prime)] = (determinant_mod(dense, prime), nonzero)
            print(PREFIX + " MODP " + json.dumps(
                {"q": label, "prime": prime, "determinant_residue": modp[(label, prime)][0],
                 "nonzero_entries": nonzero, "independent_steps_total": total,
                 "independent_steps_maximum": maximum},
                sort_keys=True, separators=(",", ":")), flush=True)
            del columns, dense, rules_mod, trace_mod
        if label in exact_labels:
            exact[label] = exact_point_record(label, rules_point, trace_point, standard_list)
        del rules_point, trace_point
        gc.collect()
        wave27.guard_resources(started, f"census point q={label}")
        stage(f"point q={label} specialized and cross-checked", started, wall)
    del rules, trace_row
    for label in exact_labels:
        pinned = pinned_residues.get((label, PRIMARY_PRIME))
        if pinned is not None and exact[label]["residue_mod_2147483647"] != pinned:
            fail(f"independent exact determinant at q={label} misses the landed residue")
        for prime in (PRIMARY_PRIME, SECONDARY_PRIME):
            key = "residue_mod_2147483647" if prime == PRIMARY_PRIME else "residue_mod_2147483629"
            if exact[label][key] != modp[(label, prime)][0]:
                fail(f"independent exact and mod-{prime} determinants disagree at q={label}")
        if exact[label]["residue_mod_2147483647"] != norm_residue[label]:
            fail(f"independent exact determinant at q={label} disagrees with the e253 reduced modular norm")

    ctx = {
        "pinned_residues": pinned_residues,
        "pinned_steps": pinned_steps,
        "norm_residue": norm_residue,
        "modp": modp,
        "exact": exact,
        "step_bound": step_bound,
    }
    summary = validate_artifact(artifact, sources, ctx)
    rejected = mutation_rejections(artifact, sources, ctx)
    if rejected != MUTATION_IDS:
        fail(f"mutation rejection mismatch: {rejected}")
    if artifact_path.read_bytes() != artifact_bytes:
        fail("artifact bytes changed during verification")
    final = stage("complete", started, wall)
    print(
        PREFIX + " PASS "
        f"scope={SCOPE} points={len(POINT_LABELS)} exact_points={','.join(exact_labels)} "
        f"signs={''.join('+' if s > 0 else '-' if s < 0 else '0' for s in summary['signs'])} "
        f"zero_count={summary['zero_count']} brackets={len(summary['brackets'])} "
        f"mutations={len(rejected)} rules_sha256={RULE_SHA256} "
        f"cpu={final['cpu_seconds']}s wall={final['wall_seconds']}s "
        f"current_resident={final['current_resident_bytes']} "
        f"current_footprint={final['current_footprint_bytes']} "
        f"peak_rss={final['peak_rss_bytes']}",
        flush=True,
    )
    return 0


def smoke() -> int:
    """Self-tests of the independent scalar routines; no lift, no artifact."""
    # Bareiss versus FLINT on small deterministic matrices, including swaps and singular cases.
    seed = 12345

    def next_value() -> int:
        nonlocal seed
        seed = (seed * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        return (seed >> 20) % 2001 - 1000

    cases = 0
    for size in (1, 2, 3, 5, 8, 11):
        for trial in range(4):
            rows = [[next_value() * (3 ** (trial * i)) for j in range(size)] for i in range(size)]
            if trial == 1 and size > 1:
                rows[0][0] = 0
            if trial == 2 and size > 1:
                rows[1] = list(rows[0])
            expected = int(fmpz_mat(rows).det())
            if bareiss_determinant(rows) != expected:
                fail(f"bareiss determinant mismatch at size {size} trial {trial}")
            prime = 1000003
            if determinant_mod(rows, prime) != expected % prime:
                fail(f"modular determinant mismatch at size {size} trial {trial}")
            cases += 1
    huge = [[(-1) ** ((i * j) % 3) * (7 ** (60 * i + 11 * j) + i - j) for j in range(6)] for i in range(6)]
    if bareiss_determinant(huge) != int(fmpz_mat(huge).det()):
        fail("bareiss determinant mismatch on a wide-entry matrix")
    cases += 1
    # Scaled polynomial evaluation versus Fraction arithmetic.
    coefficients = [3, -2, 0, 5, -(1 << 70)]
    bpows = [3 ** k for k in range(5)]
    expected = sum(Fraction(c) * Fraction(5, 3) ** i for i, c in enumerate(coefficients))
    if Fraction(scaled_polynomial_value(coefficients, 5, bpows), 3 ** 4) != expected:
        fail("scaled polynomial evaluation mismatch")
    if scaled_polynomial_value([], 5, bpows) != 0:
        fail("empty polynomial must evaluate to zero")
    cases += 1
    # Horner mod p versus direct evaluation.
    prime = PRIMARY_PRIME
    direct = sum(c * pow(715827884, i, prime) for i, c in enumerate(coefficients)) % prime
    if horner_mod(coefficients, 715827884, prime) != direct:
        fail("horner_mod mismatch")
    cases += 1
    # Reducer on a Groebner basis with two competing leaders: I = (x0^2 - x0/2, x1^2 - 3).
    # Lift convention: a rule (leader, tail) encodes leader + sum(tail) in the ideal,
    # so x0^2 -> x0/2 is stored with tail coefficient -1/2 and x1^2 -> 3 with -3.
    rules = [
        ((2, 0, 0, 0, 0), {(1, 0, 0, 0, 0): Fraction(-1, 2)}),
        ((0, 2, 0, 0, 0), {(0, 0, 0, 0, 0): Fraction(-3)}),
    ]
    leaders = [leading for leading, _ in rules]
    lookup: dict[Monomial, int | None] = {}
    remainder, steps = normal_form({(3, 1, 0, 0, 0): Fraction(1), (2, 2, 0, 0, 0): Fraction(2)},
                                   rules, leaders, lookup, identity)
    # x0^3 x1 -> x0 x1 / 4 ; 2 x0^2 x1^2 -> 2 * (x0/2) * 3 = 3 x0
    if remainder != {(1, 1, 0, 0, 0): Fraction(1, 4), (1, 0, 0, 0, 0): Fraction(3)}:
        fail(f"exact reducer mismatch: {remainder}")
    if lookup[(2, 2, 0, 0, 0)] != 1:
        fail("reducer must choose the last divisible leader")
    remainder_mod, _ = normal_form({(3, 1, 0, 0, 0): 1, (2, 2, 0, 0, 0): 2}, [
        (leading, {m: to_mod(v, 101) for m, v in tail.items()}) for leading, tail in rules
    ], leaders, {}, lambda v: v % 101)
    if remainder_mod != {(1, 1, 0, 0, 0): to_mod(Fraction(1, 4), 101), (1, 0, 0, 0, 0): 3}:
        fail("modular reducer mismatch")
    if steps <= 0:
        fail("reducer step count must be positive")
    cases += 1
    # Residue helper and pole detection.
    if residue(Fraction(-7, 12), 101) != (-7 * pow(12, -1, 101)) % 101 or residue(Fraction(5, 101), 101) is not None:
        fail("residue helper mismatch")
    cases += 1
    sample = guard.sample_process(os.getpid())
    for key in ("resident_bytes", "footprint_bytes"):
        if not is_int(sample.get(key)) or sample[key] <= 0:
            fail(f"resource guard sample lacks a positive {key}")
    print(PREFIX + f" SMOKE PASS cases={cases}", flush=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--exact-points", default="2",
                        help="comma-separated census points recomputed exactly (default: 2)")
    parser.add_argument("--smoke", action="store_true", help="scalar self-tests only; no lift, no artifact")
    args = parser.parse_args(argv)
    try:
        smoke()
        if args.smoke:
            return 0
        exact_labels = tuple(label.strip() for label in args.exact_points.split(",") if label.strip())
        if not exact_labels or any(label not in POINT_LABELS for label in exact_labels):
            fail(f"--exact-points must name census points among {POINT_LABELS}")
        return verify(args.artifact.resolve(), exact_labels)
    except Exception as error:  # noqa: BLE001 - any failure is a verification failure
        print(PREFIX + f" FAIL {type(error).__name__}: {error}", flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
