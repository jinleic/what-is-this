#!/usr/bin/env python3
"""Exhaustive finite check of the entropy-to-union-closed bridge.

Scope
-----
The union-closed corollary in ``uc/paper/main.tex`` (Corollary ``cor:uc``) is a
human proof.  Its load-bearing finite content is:

1. the sequential coupling of ``A`` and ``C`` exists, with nonnegative masses,
   Bernoulli marginals ``p`` and ``r``, and ``P(A_i or C_i)=s*(p,r)``;
2. both ``A`` and ``C`` stay exactly uniform on the family, prefix by prefix;
3. ``(P_i,R_i)`` is a self-coupling: ``P_i`` and ``R_i`` are identically
   distributed, and their common mean is the frequency of element ``i``;
4. the chain-rule identity ``H(A)=sum_i E h(P_i)``;
5. the two data-processing inequalities
   ``H(A u B) >= sum_i E h(P_i+Q_i-P_iQ_i)`` and
   ``H(A u C) >= sum_i E h(s*(P_i,R_i))``.

This module rebuilds all of that from the manuscript formulas alone and
verifies it on **every** nonempty family of subsets of ``[n]`` for the
requested ``n``.  Items 1-4 are checked in exact rational arithmetic; item 5 is
checked with 256-bit Arb balls, so a reported ``certified_nonnegative`` slack is
a rigorous sign statement about that family, not a floating-point impression.

Independence
------------
Nothing here imports ``uc/bridge_uc.py``, ``uc/cert3*.py``, or any frozen
campaign snapshot.  ``s*``, the coupling masses, the entropies, and the family
enumeration are written directly from the manuscript.

Limitations
-----------
This is an exhaustive *finite* control.  It proves the listed statements for
the enumerated families only; the universal statement is the induction printed
in the manuscript.  Arb (python-flint) is a trusted dependency.

Usage:
    python -I -B uc/verification/entropy_bridge_exhaustive.py \
        --max-coordinates 3 [--output PATH]
"""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
import os
import random
import sys
import time
import uuid
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx, fmpq

SCHEMA = "uc-entropy-bridge-exhaustive-report-v1"
PRECISION_BITS = 256
EXPECTED_FLINT_VERSION = "0.9.0"
HERE = Path(__file__).resolve().parent
SCHEMA_PATH = HERE / "report.schema.json"
ZERO = Fraction(0)
ONE = Fraction(1)
HALF = Fraction(1, 2)
T_CERT = Fraction(955165028125263, 2500000000000000)
# Exactly representable 2^-200; the structural ties below are true equalities,
# so a straddling enclosure must still be this small to be accepted.
EQUALITY_TOLERANCE = arb(2) ** -200
SAMPLE_SEED = 20260827


class CheckError(RuntimeError):
    pass


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def object_hash(value: dict[str, object], omitted: str) -> str:
    reduced = {key: item for key, item in value.items() if key != omitted}
    return hashlib.sha256(canonical_bytes(reduced)).hexdigest()


def file_hash(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1 << 20)
            if not chunk:
                break
            size += len(chunk)
            digest.update(chunk)
    return {"sha256": digest.hexdigest(), "size": size}


def sstar(p: Fraction, r: Fraction) -> Fraction:
    """Cambie's ``max(p,r,min(p+r,1/2))`` on exact rationals."""
    return max(p, r, min(p + r, HALF))


_ENTROPY_CACHE: dict[Fraction, arb] = {}


def binary_entropy(value: Fraction) -> arb:
    """Base-two binary entropy of an exact rational, as a 256-bit Arb ball."""
    cached = _ENTROPY_CACHE.get(value)
    if cached is not None:
        return cached
    if value < 0 or value > 1:
        raise CheckError(f"entropy argument {value} left [0,1]")
    if value in (ZERO, ONE):
        answer = arb(0)
    else:
        ball = arb(fmpq(value.numerator, value.denominator))
        one = arb(1)
        answer = -(ball * ball.log() + (one - ball) * (one - ball).log()) / arb(2).log()
    _ENTROPY_CACHE[value] = answer
    return answer


def distribution_entropy(masses: dict[object, Fraction]) -> arb:
    """Shannon entropy in bits of an exact rational distribution."""
    total = sum(masses.values(), ZERO)
    if total != ONE:
        raise CheckError(f"distribution mass {total} is not one")
    accumulator = arb(0)
    for mass in masses.values():
        if mass < 0:
            raise CheckError("negative mass in a distribution")
        if mass == ZERO:
            continue
        ball = arb(fmpq(mass.numerator, mass.denominator))
        accumulator -= ball * ball.log()
    return accumulator / arb(2).log()


def add_mass(mapping: dict, key: object, mass: Fraction) -> None:
    mapping[key] = mapping.get(key, ZERO) + mass


def prefix_parameter(family: tuple[int, ...], prefix: int, coordinate: int) -> Fraction:
    """``p_i(a_{<i})`` for the uniform law on ``family``."""
    mask = (1 << coordinate) - 1
    eligible = [member for member in family if member & mask == prefix]
    if not eligible:
        raise CheckError("prefix parameter requested at a null prefix")
    bit = 1 << coordinate
    return Fraction(sum(1 for member in eligible if member & bit), len(eligible))


def enumerate_families(coordinate_count: int):
    universe_size = 1 << coordinate_count
    for family_mask in range(1, 1 << universe_size):
        yield tuple(
            member for member in range(universe_size) if family_mask & (1 << member)
        )


def random_union_closed(
    generator: random.Random, coordinate_count: int, generator_count: int
) -> tuple[int, ...]:
    """A random union-closed family: close random generators under union."""
    universe_size = 1 << coordinate_count
    while True:
        members = {
            generator.randrange(universe_size) for _ in range(generator_count)
        }
        frontier = list(members)
        while frontier:
            current = frontier.pop()
            for other in tuple(members):
                union = current | other
                if union not in members:
                    members.add(union)
                    frontier.append(union)
        if members != {0}:
            return tuple(sorted(members))


def classify_slack(slack: arb, label: str) -> tuple[str, float]:
    """Certify the sign of one entropy slack, or place it at exact equality.

    ``nonnegative`` means Arb proved ``slack >= 0``.  ``equality`` means the
    enclosure straddles zero but lies inside ``+/-1e-60``, which is the
    structural tie the chain rule produces when conditioning adds nothing.
    Anything else is a failure.
    """
    if slack >= 0:
        return ("positive" if slack > 0 else "nonnegative"), float(slack.mid())
    if slack < -EQUALITY_TOLERANCE or slack > EQUALITY_TOLERANCE:
        raise CheckError(f"{label} data-processing inequality failed")
    return "equality", float(slack.mid())


def check_family(family: tuple[int, ...], coordinate_count: int) -> dict[str, object]:
    """Verify every finite ingredient of the bridge for one family."""
    size = len(family)
    uniform = Fraction(1, size)
    members = set(family)

    # A and B: two independent uniform copies; the pair law is a product.
    # A and C: the sequential coupling of the manuscript.
    joint = {(0, 0): ONE}
    prefix_law = {0: ONE}

    iid_or_terms = arb(0)
    coupled_or_terms = arb(0)
    single_terms = arb(0)
    coupled_prefix_states = 0

    for coordinate in range(coordinate_count):
        frequency = Fraction(
            sum(1 for member in family if member & (1 << coordinate)), size
        )

        # --- independent pair: laws of P_i and Q_i come from prefix_law ---
        law_p: dict[Fraction, Fraction] = {}
        for prefix, mass in prefix_law.items():
            add_mass(law_p, prefix_parameter(family, prefix, coordinate), mass)
        if sum(law_p.values(), ZERO) != ONE:
            raise CheckError("prefix law lost mass")
        if sum(value * mass for value, mass in law_p.items()) != frequency:
            raise CheckError("E[P_i] differs from the frequency of element i")

        for value, mass in law_p.items():
            single_terms += arb(fmpq(mass.numerator, mass.denominator)) * binary_entropy(value)
        for value_p, mass_p in law_p.items():
            for value_q, mass_q in law_p.items():
                weight = mass_p * mass_q
                union = value_p + value_q - value_p * value_q
                iid_or_terms += (
                    arb(fmpq(weight.numerator, weight.denominator))
                    * binary_entropy(union)
                )

        # --- sequential coupling of A and C ---
        law_r: dict[Fraction, Fraction] = {}
        law_p_joint: dict[Fraction, Fraction] = {}
        next_joint: dict[tuple[int, int], Fraction] = {}
        for (a_prefix, c_prefix), state_mass in joint.items():
            coupled_prefix_states += 1
            p = prefix_parameter(family, a_prefix, coordinate)
            r = prefix_parameter(family, c_prefix, coordinate)
            add_mass(law_p_joint, p, state_mass)
            add_mass(law_r, r, state_mass)
            s = sstar(p, r)
            masses = {
                (1, 1): p + r - s,
                (1, 0): s - r,
                (0, 1): s - p,
                (0, 0): ONE - s,
            }
            if any(mass < 0 for mass in masses.values()):
                raise CheckError("sequential coupling produced a negative mass")
            if sum(masses.values(), ZERO) != ONE:
                raise CheckError("sequential coupling masses do not sum to one")
            if sum(m for (a_bit, _), m in masses.items() if a_bit) != p:
                raise CheckError("first Bernoulli marginal is not p")
            if sum(m for (_, c_bit), m in masses.items() if c_bit) != r:
                raise CheckError("second Bernoulli marginal is not r")
            if sum(m for bits, m in masses.items() if any(bits)) != s:
                raise CheckError("coupled OR probability is not s*")
            coupled_or_terms += (
                arb(fmpq(state_mass.numerator, state_mass.denominator))
                * binary_entropy(s)
            )
            for (a_bit, c_bit), mass in masses.items():
                if mass == ZERO:
                    continue
                add_mass(
                    next_joint,
                    (a_prefix | (a_bit << coordinate), c_prefix | (c_bit << coordinate)),
                    state_mass * mass,
                )

        if law_p_joint != law_p or law_r != law_p:
            raise CheckError("(P_i,R_i) is not an equal-marginal self-coupling")
        joint = next_joint

        # --- prefix uniformity of both sequential copies ---
        prefix_mask = (1 << (coordinate + 1)) - 1
        expected_prefix: dict[int, Fraction] = {}
        for member in family:
            add_mass(expected_prefix, member & prefix_mask, uniform)
        first: dict[int, Fraction] = {}
        second: dict[int, Fraction] = {}
        for (a_prefix, c_prefix), mass in joint.items():
            add_mass(first, a_prefix, mass)
            add_mass(second, c_prefix, mass)
        if first != expected_prefix or second != expected_prefix:
            raise CheckError("a sequential copy left the uniform prefix law")
        prefix_law = expected_prefix

    # --- terminal laws ---
    marginal_a: dict[int, Fraction] = {}
    marginal_c: dict[int, Fraction] = {}
    union_ac: dict[int, Fraction] = {}
    for (a_value, c_value), mass in joint.items():
        add_mass(marginal_a, a_value, mass)
        add_mass(marginal_c, c_value, mass)
        add_mass(union_ac, a_value | c_value, mass)
    expected_uniform = {member: uniform for member in family}
    if marginal_a != expected_uniform or marginal_c != expected_uniform:
        raise CheckError("a sequential copy is not uniform on the family")

    union_ab: dict[int, Fraction] = {}
    pair_mass = uniform * uniform
    for first_member in family:
        for second_member in family:
            add_mass(union_ab, first_member | second_member, pair_mass)

    entropy_a = distribution_entropy(expected_uniform)
    entropy_ab = distribution_entropy(union_ab)
    entropy_ac = distribution_entropy(union_ac)

    chain_gap = entropy_a - single_terms
    iid_slack = entropy_ab - iid_or_terms
    coupled_slack = entropy_ac - coupled_or_terms

    if not (chain_gap < EQUALITY_TOLERANCE and chain_gap > -EQUALITY_TOLERANCE):
        raise CheckError("chain rule H(A)=sum_i E h(P_i) was not reproduced")
    iid_class, iid_value = classify_slack(iid_slack, "iid")
    coupled_class, coupled_value = classify_slack(coupled_slack, "coupled")

    union_closed = all(
        (first_member | second_member) in members
        for first_member in family
        for second_member in family
    )
    maximum_frequency = max(
        Fraction(sum(1 for member in family if member & (1 << coordinate)), size)
        for coordinate in range(coordinate_count)
    )
    if union_closed and family != (0,) and maximum_frequency < T_CERT:
        raise CheckError(
            "a union-closed family violated the certified frequency conclusion"
        )
    return {
        "coupled_class": coupled_class,
        "coupled_prefix_states": coupled_prefix_states,
        "coupled_value": coupled_value,
        "iid_class": iid_class,
        "iid_value": iid_value,
        "maximum_frequency": maximum_frequency,
        "union_closed": union_closed,
    }


def run(
    max_coordinates: int,
    deep_samples: dict[int, int],
    closed_samples: dict[int, int],
    output: Path | None,
) -> int:
    if not (
        sys.flags.isolated == 1
        and sys.flags.no_user_site == 1
        and sys.flags.ignore_environment == 1
        and sys.flags.dont_write_bytecode == 1
    ):
        raise CheckError("invoke with an isolated no-bytecode interpreter: python -I -B")
    import flint

    version = getattr(flint, "__version__", "unknown")
    if version != EXPECTED_FLINT_VERSION:
        raise CheckError(
            f"python-flint version {version!r} differs from pinned {EXPECTED_FLINT_VERSION!r}"
        )
    ctx.prec = PRECISION_BITS

    started = time.monotonic()
    levels: dict[str, object] = {}
    for coordinate_count in range(1, max_coordinates + 1):
        families = 0
        union_closed_families = 0
        prefix_states = 0
        classes = {
            "iid": {"positive": 0, "nonnegative": 0, "equality": 0},
            "coupled": {"positive": 0, "nonnegative": 0, "equality": 0},
        }
        smallest = {"iid": None, "coupled": None}
        tightest_union_closed = None
        for family in enumerate_families(coordinate_count):
            summary = check_family(family, coordinate_count)
            families += 1
            prefix_states += summary["coupled_prefix_states"]
            if summary["union_closed"]:
                union_closed_families += 1
                if family != (0,):
                    frequency = summary["maximum_frequency"]
                    if tightest_union_closed is None or frequency < tightest_union_closed:
                        tightest_union_closed = frequency
            for label in ("iid", "coupled"):
                classification = summary[f"{label}_class"]
                classes[label][classification] += 1
                if classification == "positive":
                    value = summary[f"{label}_value"]
                    if smallest[label] is None or value < smallest[label]:
                        smallest[label] = value
        levels[str(coordinate_count)] = {
            "coupled_prefix_states": prefix_states,
            "coupled_slack_classes": classes["coupled"],
            "families": families,
            "iid_slack_classes": classes["iid"],
            "smallest_certified_positive_coupled_slack": smallest["coupled"],
            "smallest_certified_positive_iid_slack": smallest["iid"],
            "smallest_union_closed_max_frequency": (
                None if tightest_union_closed is None
                else f"{tightest_union_closed.numerator}/{tightest_union_closed.denominator}"
            ),
            "union_closed_families": union_closed_families,
        }

    sampled: dict[str, object] = {}
    for coordinate_count, sample_size in sorted(deep_samples.items()):
        generator = random.Random(SAMPLE_SEED + coordinate_count)
        universe_size = 1 << coordinate_count
        families = 0
        union_closed_families = 0
        prefix_states = 0
        classes = {
            "iid": {"positive": 0, "nonnegative": 0, "equality": 0},
            "coupled": {"positive": 0, "nonnegative": 0, "equality": 0},
        }
        while families < sample_size:
            mask = generator.randrange(1, 1 << universe_size)
            family = tuple(
                member for member in range(universe_size) if mask & (1 << member)
            )
            summary = check_family(family, coordinate_count)
            families += 1
            union_closed_families += int(summary["union_closed"])
            prefix_states += summary["coupled_prefix_states"]
            for label in ("iid", "coupled"):
                classes[label][summary[f"{label}_class"]] += 1
        sampled[str(coordinate_count)] = {
            "coupled_prefix_states": prefix_states,
            "coupled_slack_classes": classes["coupled"],
            "families": families,
            "iid_slack_classes": classes["iid"],
            "seed": SAMPLE_SEED + coordinate_count,
            "union_closed_families": union_closed_families,
        }

    closed_sampled: dict[str, object] = {}
    for coordinate_count, sample_size in sorted(closed_samples.items()):
        generator = random.Random(SAMPLE_SEED + 1000 + coordinate_count)
        families = 0
        prefix_states = 0
        classes = {
            "iid": {"positive": 0, "nonnegative": 0, "equality": 0},
            "coupled": {"positive": 0, "nonnegative": 0, "equality": 0},
        }
        sizes = []
        tightest = None
        while families < sample_size:
            family = random_union_closed(
                generator,
                coordinate_count,
                2 + generator.randrange(4 * coordinate_count),
            )
            summary = check_family(family, coordinate_count)
            if not summary["union_closed"]:
                raise CheckError("union-closure generator produced a non-closed family")
            families += 1
            sizes.append(len(family))
            prefix_states += summary["coupled_prefix_states"]
            for label in ("iid", "coupled"):
                classes[label][summary[f"{label}_class"]] += 1
            frequency = summary["maximum_frequency"]
            if tightest is None or frequency < tightest:
                tightest = frequency
        closed_sampled[str(coordinate_count)] = {
            "coupled_prefix_states": prefix_states,
            "coupled_slack_classes": classes["coupled"],
            "families": families,
            "iid_slack_classes": classes["iid"],
            "largest_family_size": max(sizes),
            "seed": SAMPLE_SEED + 1000 + coordinate_count,
            "smallest_union_closed_max_frequency": (
                None if tightest is None
                else f"{tightest.numerator}/{tightest.denominator}"
            ),
        }
    elapsed = time.monotonic() - started

    report: dict[str, object] = {
        "$schema": str(SCHEMA_PATH),
        "checked_statements": [
            "sequential coupling masses are nonnegative and sum to one",
            "coupled Bernoulli marginals equal p and r",
            "P(A_i or C_i | prefixes) equals s*(p,r)",
            "A and C stay exactly uniform at every prefix and at the end",
            "law(P_i)=law(R_i) and E[P_i] equals the frequency of element i",
            "H(A)=sum_i E h(P_i) to within 1e-60",
            "H(A u B)>=sum_i E h(P_i+Q_i-P_iQ_i) with a certified Arb sign",
            "H(A u C)>=sum_i E h(s*(P_i,R_i)) with a certified Arb sign",
            "every enumerated union-closed family other than {emptyset} has an element of frequency at least t_cert",
        ],
        "claim_status": "MACHINE-VERIFIED (exhaustive finite control)",
        "elapsed_seconds": round(elapsed, 3),
        "finished_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "levels": levels,
        "sampled_levels": sampled,
        "sampled_union_closed_levels": closed_sampled,
        "limitations": [
            "Exhaustive only for the enumerated coordinate counts; deeper coordinate counts are a seeded uniform sample, not a proof.",
            "Arb (python-flint) is a trusted dependency for the entropy inequalities.",
            "This control does not verify Theorem thm:main, the certificate, or the strictness estimate at t_cert.",
        ],
        "max_coordinates": max_coordinates,
        "outcome": "PASS",
        "precision_bits": PRECISION_BITS,
        "report_type": "entropy_bridge_exhaustive",
        "runtime": {
            "python": sys.version,
            "python_executable": str(Path(sys.executable).resolve()),
            "python_flint": version,
        },
        "schema": SCHEMA,
        "verifier": file_hash(Path(__file__).resolve()),
    }
    report["report_sha256"] = object_hash(report, "report_sha256")
    payload = json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
    if output is None:
        print(payload, end="")
        return 0
    output = output.expanduser().resolve()
    output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{uuid.uuid4().hex}.tmp")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        view = memoryview(payload.encode("ascii"))
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise CheckError(f"short report write: {output}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(temporary, output)
    print(
        canonical_bytes(
            {
                "outcome": "PASS",
                "report": str(output),
                "report_sha256": report["report_sha256"],
            }
        ).decode("ascii")
    )
    return 0


def parse_sample(parser, specifications, max_coordinates, flag):
    parsed: dict[int, int] = {}
    for specification in specifications or ():
        try:
            coordinate_text, count_text = specification.split(":", 1)
            coordinate_count = int(coordinate_text)
            count = int(count_text)
        except ValueError:
            parser.error(f"{flag} expects N:COUNT with integer N and COUNT")
        if not max_coordinates < coordinate_count <= 8:
            parser.error(f"{flag} N must exceed --max-coordinates and be at most 8")
        if count < 1:
            parser.error(f"{flag} COUNT must be positive")
        parsed[coordinate_count] = count
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-coordinates", type=int, default=3)
    parser.add_argument(
        "--sample",
        action="append",
        default=None,
        metavar="N:COUNT",
        help="additionally check COUNT seeded random families on N coordinates",
    )
    parser.add_argument(
        "--sample-closed",
        action="append",
        default=None,
        metavar="N:COUNT",
        help="additionally check COUNT seeded random union-closed families on N coordinates",
    )
    parser.add_argument("--output", type=Path, default=None)
    arguments = parser.parse_args()
    if not 1 <= arguments.max_coordinates <= 4:
        parser.error("--max-coordinates must lie in 1..4")
    deep_samples = parse_sample(
        parser, arguments.sample, arguments.max_coordinates, "--sample"
    )
    closed_samples = parse_sample(
        parser, arguments.sample_closed, arguments.max_coordinates, "--sample-closed"
    )
    try:
        return run(
            arguments.max_coordinates, deep_samples, closed_samples, arguments.output
        )
    except CheckError as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
