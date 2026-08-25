#!/usr/bin/env python3
"""Exact model for fixed-five involution signed completion.

The module owns the deterministic four-state pseudo-Boolean encoding used after
the exact 705-orbit Ramsey support filter.  It can enumerate every linear W row
domain, emit the W^2/T^2 completion formula, reconstruct a 45-vertex graph, and
spell out the D5 and residual identical-column actions.  A solver's UNSAT answer
is not a result by itself: negative claims require a RoundingSat proof checked
by the pinned VeriPB/CakePB pipeline declared below.
"""

from __future__ import annotations

import argparse
import collections
import functools
import hashlib
import itertools
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_signed_completion"
RELATION_STATES = ("00", "01", "10", "11")
W_NONZERO_STATES = ("00", "11")
T_NONZERO_STATES = ("01", "10")
VERIPB_COMMIT = "6d38dab246af9c321b8f17cb5a187f2fbb9e491d"
CAKEPB_COMMIT = "a7593ef22de2fc0b47a688f2d4f08e6b742735af"
ROUNDINGSAT_COMMIT = "d4edbf7908a9bb951fd181940919e0f3ac7ab1ee"
CHECKED_DELETION_REQUIRED = True
PROJECTED_COUNTS_ARE_ORBIT_COUNTS = False

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FRONTIER = ROOT / "data" / "involution_f5_ramsey_filter.json"
DEFAULT_SUPPORT = ROOT / "data" / "involution_f5_support_census.json"
FRONTIER_BYTES = 513_707
FRONTIER_SHA256 = "9dfc89b79b36041a71c113d1f1d4c4a76478b95a4131fb5b12666c16be681f02"
SUPPORT_BYTES = 503_066
SUPPORT_SHA256 = "520e2cb453cc2efee7925636af018316acca3e375a8107886c5e9c90c3727aaa"

Q_FIXED = tuple(tuple(
    0 if left == right else
    -1 if (left - right) % 5 in (1, 4) else 1
    for right in range(5)
) for left in range(5))
DIHEDRAL_PERMUTATIONS = tuple(
    [tuple((vertex + shift) % 5 for vertex in range(5))
     for shift in range(5)]
    + [tuple((shift - vertex) % 5 for vertex in range(5))
       for shift in range(5)]
)

SemanticVariable = tuple


class CompletionViolation(RuntimeError):
    """A pinned input, exact model, witness, or group action is invalid."""


@dataclass(frozen=True)
class CompletionInstance:
    """One canonical Ramsey-support survivor in deterministic orbit order."""

    source_index: int
    orbit_size: int
    masks: tuple[int, ...]
    internal: tuple[int, ...]
    restrictions: tuple[tuple[int, int, tuple[int, ...]], ...]
    edge_mask_order: tuple[int, ...]
    nonedge_mask_order: tuple[int, ...]

    @property
    def r_columns(self) -> tuple[tuple[int, ...], ...]:
        return tuple(tuple(1 - 2 * ((mask >> vertex) & 1)
                           for vertex in range(5))
                     for mask in self.masks)


@dataclass(frozen=True)
class PBConstraint:
    terms: tuple[tuple[int, int], ...]
    comparator: str
    rhs: int
    label: str

    def is_satisfied(self, assignment: Mapping[int, int]) -> bool:
        value = sum(coefficient * assignment.get(variable, 0)
                    for coefficient, variable in self.terms)
        if self.comparator == "=":
            return value == self.rhs
        if self.comparator == ">=":
            return value >= self.rhs
        raise CompletionViolation(f"unsupported comparator {self.comparator}")


class VariableRegistry:
    """Deterministic bijection between semantic auxiliaries and OPB x-IDs."""

    def __init__(self) -> None:
        self._semantic_to_id: dict[SemanticVariable, int] = {}
        self._id_to_semantic: list[SemanticVariable | None] = [None]

    def id(self, semantic: SemanticVariable) -> int:
        if semantic not in self._semantic_to_id:
            variable = len(self._id_to_semantic)
            self._semantic_to_id[semantic] = variable
            self._id_to_semantic.append(semantic)
        return self._semantic_to_id[semantic]

    def lookup(self, semantic: SemanticVariable) -> int:
        try:
            return self._semantic_to_id[semantic]
        except KeyError as error:
            raise CompletionViolation(f"unknown semantic variable {semantic}") from error

    def semantic(self, variable: int) -> SemanticVariable:
        if type(variable) is not int or not 0 < variable < len(self._id_to_semantic):
            raise CompletionViolation(f"invalid variable ID {variable}")
        semantic = self._id_to_semantic[variable]
        assert semantic is not None
        return semantic

    def semantic_names(self) -> tuple[SemanticVariable, ...]:
        return tuple(self._id_to_semantic[1:])  # type: ignore[return-value]

    @property
    def variable_count(self) -> int:
        return len(self._id_to_semantic) - 1


class PBFormula:
    """Small exact OPB formula with normalized integer terms."""

    def __init__(self) -> None:
        self.constraints: list[PBConstraint] = []

    @staticmethod
    def _normalize(terms: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
        coefficients: collections.Counter[int] = collections.Counter()
        for coefficient, variable in terms:
            if type(coefficient) is not int or type(variable) is not int:
                raise CompletionViolation("PB terms must be exact integers")
            if variable <= 0:
                raise CompletionViolation("PB variable IDs must be positive")
            coefficients[variable] += coefficient
        return tuple((coefficient, variable)
                     for variable, coefficient in sorted(coefficients.items())
                     if coefficient)

    def add_equality(self, terms: Iterable[tuple[int, int]], rhs: int,
                     label: str) -> None:
        if type(rhs) is not int:
            raise CompletionViolation("PB right-hand sides must be exact integers")
        self.constraints.append(PBConstraint(
            self._normalize(terms), "=", rhs, label))

    def add_at_least(self, terms: Iterable[tuple[int, int]], rhs: int,
                     label: str) -> None:
        if type(rhs) is not int:
            raise CompletionViolation("PB right-hand sides must be exact integers")
        self.constraints.append(PBConstraint(
            self._normalize(terms), ">=", rhs, label))

    @property
    def constraint_count(self) -> int:
        return len(self.constraints)

    def to_opb(self, registry: VariableRegistry) -> str:
        lines = [
            f"* #variable= {registry.variable_count} "
            f"#constraint= {self.constraint_count} "
            f"#equal= {sum(constraint.comparator == '=' for constraint in self.constraints)} "
            "intsize= 64",
            f"* campaign={CAMPAIGN_ID}",
            f"* roundingsat_base_commit={ROUNDINGSAT_COMMIT}",
            f"* veripb_commit={VERIPB_COMMIT}",
            f"* cakepb_commit={CAKEPB_COMMIT}",
            "* checked_deletion_required=true",
        ]
        for constraint in self.constraints:
            if constraint.label:
                lines.append(f"* {constraint.label}")
            body = " ".join(
                f"{coefficient:+d} x{variable}"
                for coefficient, variable in constraint.terms
            )
            if not body:
                body = "0"
            lines.append(
                f"{body} {constraint.comparator} {constraint.rhs};")
        return "\n".join(lines) + "\n"


def _load_pinned(path: Path, canonical: Path, expected_bytes: int,
                 expected_sha256: str, label: str) -> dict:
    if path.resolve() != canonical.resolve():
        raise CompletionViolation(f"{label} path is not canonical")
    raw = path.read_bytes()
    if len(raw) != expected_bytes:
        raise CompletionViolation(
            f"{label} byte count {len(raw)} != {expected_bytes}")
    digest = hashlib.sha256(raw).hexdigest()
    if digest != expected_sha256:
        raise CompletionViolation(f"{label} SHA-256 {digest} is not pinned")
    try:
        document = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise CompletionViolation(f"cannot parse {label}") from error
    if type(document) is not dict:
        raise CompletionViolation(f"{label} root must be an object")
    return document


def _expand_representative(census: dict, representative: dict
                           ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    masks: list[int] = []
    internal: list[int] = []
    for mask, count in zip(census["edge_mask_order"],
                           representative["edge_even_counts"], strict=True):
        masks.extend([mask] * count)
        internal.extend([1] * count)
    for mask, count in zip(census["nonedge_mask_order"],
                           representative["nonedge_odd_counts"], strict=True):
        masks.extend([mask] * count)
        internal.extend([0] * count)
    if len(masks) != 20 or len(internal) != 20:
        raise CompletionViolation("support representative does not expand to 20 columns")
    return tuple(masks), tuple(internal)


@functools.cache
def load_all_instances() -> tuple[CompletionInstance, ...]:
    """Load all 705 canonical survivors from the two hash-pinned artifacts."""
    frontier = _load_pinned(
        DEFAULT_FRONTIER, DEFAULT_FRONTIER, FRONTIER_BYTES, FRONTIER_SHA256,
        "Ramsey frontier")
    support = _load_pinned(
        DEFAULT_SUPPORT, DEFAULT_SUPPORT, SUPPORT_BYTES, SUPPORT_SHA256,
        "support census")
    candidates = frontier.get("census", {}).get("candidates")
    representatives = support.get("support_census", {}).get("representatives")
    if type(candidates) is not list or type(representatives) is not list:
        raise CompletionViolation("pinned artifacts have invalid candidate arrays")
    census = support["support_census"]
    instances = []
    for candidate in candidates:
        if candidate.get("status") != "SURVIVES_R33_SUPPORT_FILTER":
            continue
        source_index = candidate.get("source_index")
        orbit_size = candidate.get("orbit_size")
        if type(source_index) is not int or not 0 <= source_index < len(representatives):
            raise CompletionViolation("invalid survivor source index")
        if orbit_size not in (1, 5, 10):
            raise CompletionViolation("invalid survivor orbit size")
        masks, internal = _expand_representative(
            census, representatives[source_index])
        restrictions = []
        for record in candidate.get("relation_restrictions", []):
            left = record.get("left_orbit")
            right = record.get("right_orbit")
            values = record.get("w_ij_over_2")
            if (type(left) is not int or type(right) is not int
                    or type(values) is not list
                    or any(type(value) is not int or value not in (-1, 0, 1)
                           for value in values)):
                raise CompletionViolation("invalid R(3,3) relation restriction")
            restrictions.append((min(left, right), max(left, right), tuple(values)))
        instances.append(CompletionInstance(
            source_index=source_index,
            orbit_size=orbit_size,
            masks=masks,
            internal=internal,
            restrictions=tuple(restrictions),
            edge_mask_order=tuple(census["edge_mask_order"]),
            nonedge_mask_order=tuple(census["nonedge_mask_order"]),
        ))
    if len(instances) != 705 or sum(item.orbit_size for item in instances) != 6627:
        raise CompletionViolation("pinned frontier does not contain the exact 705/6627 census")
    return tuple(instances)


@functools.cache
def _instances_by_source() -> dict[int, CompletionInstance]:
    return {instance.source_index: instance for instance in load_all_instances()}


def load_instance(source_index: int) -> CompletionInstance:
    if type(source_index) is not int:
        raise CompletionViolation("source index must be an exact integer")
    try:
        return _instances_by_source()[source_index]
    except KeyError as error:
        raise CompletionViolation(
            f"source index {source_index} is not one of the 705 survivors") from error


def relation_values(state: str) -> tuple[int, int, int, int]:
    """Return ``(p,q,W/2,T/2)`` for one exact relation state."""
    if state not in RELATION_STATES:
        raise CompletionViolation(f"invalid relation state {state!r}")
    p, q = (int(bit) for bit in state)
    return p, q, 1 - p - q, q - p
def w_domain_code(domain: tuple[int, ...]) -> int:
    """Encode a twenty-trit W row as a canonical unsigned base-three word."""
    if (type(domain) is not tuple or len(domain) != 20
            or any(type(value) is not int or value not in (-1, 0, 1)
                   for value in domain)):
        raise CompletionViolation("invalid W domain for base-three encoding")
    code = 0
    for value in domain:
        code = 3 * code + value + 1
    return code


def decode_w_domain(code: int) -> tuple[int, ...]:
    if type(code) is not int or not 0 <= code < 3 ** 20:
        raise CompletionViolation("invalid base-three W domain code")
    values = [0] * 20
    for index in range(19, -1, -1):
        values[index] = code % 3 - 1
        code //= 3
    return tuple(values)


def w_row_semantic(row: int, domain: tuple[int, ...]) -> SemanticVariable:
    if type(row) is not int or not 0 <= row < 20 or domain[row] != 0:
        raise CompletionViolation("invalid W row selector")
    return ("wrow", row, w_domain_code(domain))




def _pair(left: int, right: int) -> tuple[int, int]:
    if (type(left) is not int or type(right) is not int
            or not 0 <= left < 20 or not 0 <= right < 20 or left == right):
        raise CompletionViolation("relation endpoints must be distinct orbit indices")
    return (left, right) if left < right else (right, left)


def state_semantic(left: int, right: int, state: str) -> SemanticVariable:
    if state not in RELATION_STATES:
        raise CompletionViolation(f"invalid relation state {state!r}")
    left, right = _pair(left, right)
    return ("state", left, right, state)


def state_variable(registry: VariableRegistry, left: int, right: int,
                   state: str) -> int:
    return registry.id(state_semantic(left, right, state))


def _product_semantic(kind: str, left: int, right: int, centre: int,
                      left_state: str, right_state: str) -> SemanticVariable:
    if kind not in ("wprod", "tprod"):
        raise CompletionViolation(f"invalid product family {kind}")
    left, right = _pair(left, right)
    allowed = W_NONZERO_STATES if kind == "wprod" else T_NONZERO_STATES
    if left_state not in allowed or right_state not in allowed:
        raise CompletionViolation("invalid signed product state")
    if type(centre) is not int or not 0 <= centre < 20 or centre in (left, right):
        raise CompletionViolation("invalid product centre")
    return (kind, left, right, centre, left_state, right_state)


def build_variable_registry(
        instance: CompletionInstance, *, include_t_products: bool = True,
        w_product_pairs: tuple[tuple[int, int], ...] | None = None
        ) -> VariableRegistry:
    del instance  # The complete semantic family is support-independent.
    registry = VariableRegistry()
    all_pairs = tuple((left, right) for left in range(20)
                      for right in range(left + 1, 20))
    w_pairs = all_pairs if w_product_pairs is None else w_product_pairs
    if (type(w_pairs) is not tuple or len(set(w_pairs)) != len(w_pairs)
            or any(type(pair) is not tuple or len(pair) != 2
                   or pair[0] >= pair[1] or not 0 <= pair[0] < pair[1] < 20
                   for pair in w_pairs)):
        raise CompletionViolation("invalid W product-pair set")
    for left, right in all_pairs:
        for state in RELATION_STATES:
            registry.id(state_semantic(left, right, state))
    product_families = [("wprod", W_NONZERO_STATES, w_pairs)]
    if include_t_products:
        product_families.append(("tprod", T_NONZERO_STATES, all_pairs))
    for kind, states, pairs in product_families:
        for left, right in pairs:
            for centre in range(20):
                if centre in (left, right):
                    continue
                for left_state in states:
                    for right_state in states:
                        registry.id(_product_semantic(
                            kind, left, right, centre,
                            left_state, right_state))
    return registry


def _permuted_product(name: SemanticVariable,
                      permutation: tuple[int, ...]) -> SemanticVariable:
    kind, left, right, centre, left_state, right_state = name
    image_left = permutation[left]
    image_right = permutation[right]
    image_centre = permutation[centre]
    if image_left < image_right:
        return _product_semantic(
            kind, image_left, image_right, image_centre,
            left_state, right_state)
    return _product_semantic(
        kind, image_right, image_left, image_centre,
        right_state, left_state)


def permute_semantic_variable(name: SemanticVariable,
                              permutation: tuple[int, ...]) -> SemanticVariable:
    if (type(permutation) is not tuple or len(permutation) != 20
            or set(permutation) != set(range(20))):
        raise CompletionViolation("orbit action is not a permutation of 0..19")
    if type(name) is not tuple or not name:
        raise CompletionViolation("invalid semantic variable")
    if name[0] == "state" and len(name) == 4:
        _, left, right, state = name
        return state_semantic(permutation[left], permutation[right], state)
    if name[0] in ("wprod", "tprod") and len(name) == 6:
        return _permuted_product(name, permutation)
    if name[0] == "wrow" and len(name) == 3:
        _, row, code = name
        domain = decode_w_domain(code)
        image = [0] * 20
        for index, value in enumerate(domain):
            image[permutation[index]] = value
        return w_row_semantic(permutation[row], tuple(image))
    raise CompletionViolation(f"unknown semantic variable {name}")


def inverse_permutation(permutation: tuple[int, ...]) -> tuple[int, ...]:
    if (type(permutation) is not tuple or len(permutation) != 20
            or set(permutation) != set(range(20))):
        raise CompletionViolation("not a permutation of 0..19")
    inverse = [0] * 20
    for source, target in enumerate(permutation):
        inverse[target] = source
    return tuple(inverse)
def permute_constraint_signature(
        constraint: PBConstraint, registry: VariableRegistry,
        permutation: tuple[int, ...]
        ) -> tuple[tuple[tuple[int, int], ...], str, int]:
    """Apply an orbit permutation to every semantic variable in a PB row."""
    terms = []
    for coefficient, variable in constraint.terms:
        image = permute_semantic_variable(
            registry.semantic(variable), permutation)
        terms.append((coefficient, registry.lookup(image)))
    return PBFormula._normalize(terms), constraint.comparator, constraint.rhs




def residual_generators(instance: CompletionInstance) -> tuple[tuple[int, ...], ...]:
    """Return one exact swap generator per identical ``(support,a)`` pair."""
    classes: dict[tuple[int, int], list[int]] = collections.defaultdict(list)
    for index, key in enumerate(zip(instance.masks, instance.internal, strict=True)):
        classes[key].append(index)
    generators = []
    for indices in sorted(classes.values()):
        for left, right in zip(indices, indices[1:], strict=False):
            permutation = list(range(20))
            permutation[left], permutation[right] = right, left
            generators.append(tuple(permutation))
    return tuple(generators)


def _permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(1 << permutation[vertex] for vertex in range(5)
               if mask >> vertex & 1)


def _support_signature(masks: tuple[int, ...], internal: tuple[int, ...],
                       edge_order: tuple[int, ...],
                       nonedge_order: tuple[int, ...]
                       ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    edge = tuple(sum(bit == 1 and mask == target
                     for mask, bit in zip(masks, internal, strict=True))
                 for target in edge_order)
    nonedge = tuple(sum(bit == 0 and mask == target
                        for mask, bit in zip(masks, internal, strict=True))
                    for target in nonedge_order)
    return edge, nonedge


def explicit_d5_images(instance: CompletionInstance
                       ) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    """Expand one canonical support through all ten fixed-cycle relabellings."""
    images = {
        _support_signature(
            tuple(_permute_mask(mask, permutation) for mask in instance.masks),
            instance.internal,
            instance.edge_mask_order,
            instance.nonedge_mask_order,
        )
        for permutation in DIHEDRAL_PERMUTATIONS
    }
    if len(images) != instance.orbit_size:
        raise CompletionViolation("explicit D5 image count disagrees with orbit size")
    return tuple(sorted(images))


def _restriction_map(instance: CompletionInstance, row: int) -> dict[int, tuple[int, ...]]:
    result = {}
    for left, right, values in instance.restrictions:
        if left == row:
            result[right] = values
        elif right == row:
            result[left] = values
    return result


def _row_target(instance: CompletionInstance, row: int) -> tuple[int, ...]:
    column = instance.r_columns[row]
    diagonal = 1 - 2 * instance.internal[row]
    return tuple((
        -1
        - sum(Q_FIXED[vertex][other] * column[other]
              for other in range(5))
        - diagonal * column[vertex]
    ) // 2 for vertex in range(5))


def _half_assignments(instance: CompletionInstance, indices: tuple[int, ...],
                      domains: Mapping[int, tuple[int, ...]], max_plus: int,
                      max_minus: int
                      ) -> dict[tuple[int, ...], list[tuple[int, int]]]:
    columns = instance.r_columns
    records: dict[tuple[int, ...], list[tuple[int, int]]] = collections.defaultdict(list)

    def visit(position: int, sums: tuple[int, ...], plus: int, minus: int,
              plus_mask: int, minus_mask: int) -> None:
        if plus > max_plus or minus > max_minus:
            return
        if position == len(indices):
            records[sums + (plus, minus)].append((plus_mask, minus_mask))
            return
        index = indices[position]
        for value in domains.get(index, (-1, 0, 1)):
            visit(
                position + 1,
                tuple(sums[vertex] + value * columns[index][vertex]
                      for vertex in range(5)),
                plus + (value == 1),
                minus + (value == -1),
                plus_mask | ((value == 1) << index),
                minus_mask | ((value == -1) << index),
            )

    visit(0, (0, 0, 0, 0, 0), 0, 0, 0, 0)
    return records


def enumerate_w_row_domains(instance: CompletionInstance, row: int
                            ) -> tuple[tuple[int, ...], ...]:
    """Enumerate exactly the linear cross/degree/R(3,3) domain for one W row."""
    if type(row) is not int or not 0 <= row < 20:
        raise CompletionViolation("W row index must lie in 0..19")
    plus = 3 + instance.masks[row].bit_count() // 2
    minus = 8 - plus
    restrictions = _restriction_map(instance, row)
    indices = tuple(index for index in range(20) if index != row)
    left_indices, right_indices = indices[:9], indices[9:]
    left = _half_assignments(
        instance, left_indices, restrictions, plus, minus)
    right = _half_assignments(
        instance, right_indices, restrictions, plus, minus)
    target = _row_target(instance, row)
    domains = []
    for key in sorted(left):
        need = tuple(target[vertex] - key[vertex] for vertex in range(5)) + (
            plus - key[5], minus - key[6])
        if need not in right:
            continue
        for left_plus, left_minus in left[key]:
            for right_plus, right_minus in right[need]:
                plus_mask = left_plus | right_plus
                minus_mask = left_minus | right_minus
                domain = tuple(
                    1 if plus_mask >> index & 1 else
                    -1 if minus_mask >> index & 1 else 0
                    for index in range(20)
                )
                domains.append(domain)
    result = tuple(sorted(domains))
    if len(result) != len(set(result)):
        raise CompletionViolation("duplicate W row domain")
    return result


def verify_w_row_domain(instance: CompletionInstance, row: int,
                        domain: tuple[int, ...]) -> None:
    if (type(domain) is not tuple or len(domain) != 20
            or any(type(value) is not int or value not in (-1, 0, 1)
                   for value in domain)
            or domain[row] != 0):
        raise CompletionViolation("invalid W row domain shape")
    plus = 3 + instance.masks[row].bit_count() // 2
    if domain.count(1) != plus or domain.count(-1) != 8 - plus:
        raise CompletionViolation("W row sign degree failed")
    for other, allowed in _restriction_map(instance, row).items():
        if domain[other] not in allowed:
            raise CompletionViolation("W row violates an R(3,3) restriction")
    target = _row_target(instance, row)
    for vertex in range(5):
        value = sum(instance.r_columns[other][vertex] * domain[other]
                    for other in range(20))
        if value != target[vertex]:
            raise CompletionViolation("W row cross equation failed")


def _add_and(formula: PBFormula, output: int, left: int, right: int,
             label: str) -> None:
    formula.add_at_least(((1, left), (-1, output)), 0, f"{label}:z<=left")
    formula.add_at_least(((1, right), (-1, output)), 0, f"{label}:z<=right")
    formula.add_at_least(
        ((1, output), (-1, left), (-1, right)), -1,
        f"{label}:z>=left+right-1")


def _signed_state_terms(registry: VariableRegistry, left: int, right: int,
                        kind: str, multiplier: int = 1
                        ) -> list[tuple[int, int]]:
    if kind == "w":
        return [
            (multiplier, registry.lookup(state_semantic(left, right, "00"))),
            (-multiplier, registry.lookup(state_semantic(left, right, "11"))),
        ]
    if kind == "t":
        return [
            (multiplier, registry.lookup(state_semantic(left, right, "01"))),
            (-multiplier, registry.lookup(state_semantic(left, right, "10"))),
        ]
    raise CompletionViolation(f"unknown signed relation kind {kind}")


def build_signed_formula(
        instance: CompletionInstance, *, include_t_square: bool = True,
        w_square_rows: tuple[int, ...] | None = None,
        include_w_domain_tables: bool = False,
        w_domain_table_rows: tuple[int, ...] | None = None,
        include_ramsey: bool = False
        ) -> tuple[PBFormula, VariableRegistry]:
    """Build the exact shared-pair W² and optional T² PB formula."""
    if include_ramsey and not include_t_square:
        raise CompletionViolation("Ramsey graph constraints require the T layer")
    all_pairs = tuple((left, right) for left in range(20)
                      for right in range(left + 1, 20))
    if w_square_rows is None:
        w_pairs = all_pairs
    else:
        if (type(w_square_rows) is not tuple
                or len(set(w_square_rows)) != len(w_square_rows)
                or len(w_square_rows) < 2
                or any(type(row) is not int or not 0 <= row < 20
                       for row in w_square_rows)):
            raise CompletionViolation("invalid W square row set")
        w_pairs = tuple(
            (left, right)
            for left, right in itertools.combinations(sorted(w_square_rows), 2)
        )
    w_pair_set = set(w_pairs)
    registry = build_variable_registry(
        instance, include_t_products=include_t_square,
        w_product_pairs=w_pairs)
    w_domains = None
    if include_w_domain_tables:
        rows = tuple(range(20)) if w_domain_table_rows is None else (
            w_domain_table_rows)
        if (type(rows) is not tuple or len(set(rows)) != len(rows)
                or any(type(row) is not int or not 0 <= row < 20
                       for row in rows)):
            raise CompletionViolation("invalid W domain-table row set")
        w_domains = {
            row: enumerate_w_row_domains(instance, row) for row in rows
        }
        for row, domains in w_domains.items():
            for domain in domains:
                registry.id(w_row_semantic(row, domain))
    elif w_domain_table_rows is not None:
        raise CompletionViolation(
            "W domain-table rows require include_w_domain_tables")
    formula = PBFormula()

    for left in range(20):
        for right in range(left + 1, 20):
            formula.add_equality(
                ((1, registry.lookup(state_semantic(left, right, state)))
                 for state in RELATION_STATES),
                1,
                f"onehot:{left}:{right}",
            )

    for left, right, allowed in instance.restrictions:
        allowed_set = set(allowed)
        for state in RELATION_STATES:
            w_value = relation_values(state)[2]
            if w_value not in allowed_set:
                formula.add_equality(
                    ((1, registry.lookup(state_semantic(left, right, state))),),
                    0,
                    f"r33:{left}:{right}:{state}",
                )

    columns = instance.r_columns
    for row in range(20):
        plus = 3 + instance.masks[row].bit_count() // 2
        minus = 8 - plus
        formula.add_equality(
            ((1, registry.lookup(state_semantic(row, other, "00")))
             for other in range(20) if other != row),
            plus,
            f"w-plus-degree:{row}",
        )
        formula.add_equality(
            ((1, registry.lookup(state_semantic(row, other, "11")))
             for other in range(20) if other != row),
            minus,
            f"w-minus-degree:{row}",
        )
        target = _row_target(instance, row)
        for vertex in range(5):
            terms = []
            for other in range(20):
                if other == row:
                    continue
                terms.extend(_signed_state_terms(
                    registry, row, other, "w", columns[other][vertex]))
            formula.add_equality(
                terms, target[vertex], f"cross:{vertex}:{row}")
    if w_domains is not None:
        for row, domains in w_domains.items():
            selectors = [
                registry.lookup(w_row_semantic(row, domain))
                for domain in domains
            ]
            formula.add_equality(
                ((1, variable) for variable in selectors), 1,
                f"wrow-onehot:{row}")
            for other in range(20):
                if other == row:
                    continue
                for value, states in (
                        (1, ("00",)), (0, ("01", "10")), (-1, ("11",))):
                    terms = [
                        (1, selector)
                        for selector, domain in zip(
                            selectors, domains, strict=True)
                        if domain[other] == value
                    ]
                    terms.extend(
                        (-1, registry.lookup(
                            state_semantic(row, other, state)))
                        for state in states
                    )
                    formula.add_equality(
                        terms, 0,
                        f"wrow-link:{row}:{other}:{value}")



    product_families = [("wprod", W_NONZERO_STATES, w_pairs)]
    if include_t_square:
        product_families.append(("tprod", T_NONZERO_STATES, all_pairs))
    for kind, states, pairs in product_families:
        for left, right in pairs:
            for centre in range(20):
                if centre in (left, right):
                    continue
                for left_state in states:
                    for right_state in states:
                        output = registry.lookup(_product_semantic(
                            kind, left, right, centre,
                            left_state, right_state))
                        left_input = registry.lookup(state_semantic(
                            left, centre, left_state))
                        right_input = registry.lookup(state_semantic(
                            right, centre, right_state))
                        _add_and(
                            formula, output, left_input, right_input,
                            f"{kind}:{left}:{right}:{centre}:"
                            f"{left_state}:{right_state}",
                        )

    w_sign = {"00": 1, "11": -1}
    t_sign = {"01": 1, "10": -1}
    for left in range(20):
        for right in range(left + 1, 20):
            include_w_square = (left, right) in w_pair_set
            w_terms = []
            if include_w_square:
                w_diagonal_factor = (
                    (1 - 2 * instance.internal[left])
                    + (1 - 2 * instance.internal[right])
                ) // 2
                if w_diagonal_factor:
                    w_terms.extend(_signed_state_terms(
                        registry, left, right, "w", w_diagonal_factor))
            t_terms = []
            if include_t_square:
                t_diagonal_factor = (
                    (2 * instance.internal[left] - 1)
                    + (2 * instance.internal[right] - 1)
                ) // 2
                if t_diagonal_factor:
                    t_terms.extend(_signed_state_terms(
                        registry, left, right, "t", t_diagonal_factor))
            for centre in range(20):
                if centre in (left, right):
                    continue
                if include_w_square:
                    for left_state in W_NONZERO_STATES:
                        for right_state in W_NONZERO_STATES:
                            w_terms.append((
                                w_sign[left_state] * w_sign[right_state],
                                registry.lookup(_product_semantic(
                                    "wprod", left, right, centre,
                                    left_state, right_state)),
                            ))
                if include_t_square:
                    for left_state in T_NONZERO_STATES:
                        for right_state in T_NONZERO_STATES:
                            t_terms.append((
                                t_sign[left_state] * t_sign[right_state],
                                registry.lookup(_product_semantic(
                                    "tprod", left, right, centre,
                                    left_state, right_state)),
                            ))
            gram = sum(columns[left][vertex] * columns[right][vertex]
                       for vertex in range(5))
            if include_w_square:
                formula.add_equality(
                    w_terms, (-1 - gram) // 2,
                    f"w-square:{left}:{right}")
            if include_t_square:
                formula.add_equality(
                    t_terms, 0, f"t-square:{left}:{right}")

    if include_ramsey:
        add_all_ramsey_constraints(formula, registry, instance)
    return formula, registry


def orbit_vertices(orbit: int) -> tuple[int, int]:
    if type(orbit) is not int or not 0 <= orbit < 20:
        raise CompletionViolation("orbit index must lie in 0..19")
    return 5 + 2 * orbit, 6 + 2 * orbit


def _vertex_orbit(vertex: int) -> tuple[int, int]:
    if not 5 <= vertex < 45:
        raise CompletionViolation("nonfixed vertex must lie in 5..44")
    return (vertex - 5) // 2, (vertex - 5) % 2


def reconstruct_graph(instance: CompletionInstance,
                      states: Mapping[tuple[int, int], str]) -> tuple[int, ...]:
    """Recover the labelled 45-vertex graph from all 190 exact p/q states."""
    expected_pairs = {(left, right) for left in range(20)
                      for right in range(left + 1, 20)}
    if set(states) != expected_pairs:
        raise CompletionViolation("relation assignment must cover all 190 pairs")
    if any(state not in RELATION_STATES for state in states.values()):
        raise CompletionViolation("relation assignment contains an invalid state")
    adjacency = [0] * 45

    def add_edge(left: int, right: int) -> None:
        adjacency[left] |= 1 << right
        adjacency[right] |= 1 << left

    for left in range(5):
        for right in range(left + 1, 5):
            if (left - right) % 5 in (1, 4):
                add_edge(left, right)
    for orbit, mask in enumerate(instance.masks):
        first, second = orbit_vertices(orbit)
        for fixed in range(5):
            if mask >> fixed & 1:
                add_edge(fixed, first)
                add_edge(fixed, second)
        if instance.internal[orbit]:
            add_edge(first, second)
    for (left, right), state in states.items():
        p, q, _, _ = relation_values(state)
        left_first, left_second = orbit_vertices(left)
        right_first, right_second = orbit_vertices(right)
        if p:
            add_edge(left_first, right_first)
            add_edge(left_second, right_second)
        if q:
            add_edge(left_first, right_second)
            add_edge(left_second, right_first)
    return tuple(adjacency)


def has_edge(adjacency: Sequence[int], left: int, right: int) -> bool:
    if not 0 <= left < len(adjacency) or not 0 <= right < len(adjacency):
        raise CompletionViolation("graph vertex out of range")
    return bool(adjacency[left] >> right & 1)


def is_srg_45_22_10_11(adjacency: Sequence[int]) -> bool:
    if len(adjacency) != 45:
        return False
    full = (1 << 45) - 1
    for vertex, row in enumerate(adjacency):
        if type(row) is not int or row & ~full or row >> vertex & 1:
            return False
        if row.bit_count() != 22:
            return False
        for other in range(vertex + 1, 45):
            if bool(row >> other & 1) != bool(adjacency[other] >> vertex & 1):
                return False
            common = (row & adjacency[other]).bit_count()
            if common != (10 if row >> other & 1 else 11):
                return False
    return True


def first_clique(adjacency: Sequence[int], order: int) -> tuple[int, ...] | None:
    def visit(candidates: int, need: int,
              chosen: tuple[int, ...]) -> tuple[int, ...] | None:
        if need == 0:
            return chosen
        if candidates.bit_count() < need:
            return None
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            witness = visit(candidates & adjacency[vertex], need - 1,
                            chosen + (vertex,))
            if witness is not None:
                return witness
        return None

    if type(order) is not int or order < 0:
        raise CompletionViolation("clique order must be a nonnegative integer")
    return visit((1 << len(adjacency)) - 1, order, ())


def first_independent_set(adjacency: Sequence[int], order: int
                          ) -> tuple[int, ...] | None:
    full = (1 << len(adjacency)) - 1
    complement = tuple(full ^ (1 << vertex) ^ row
                       for vertex, row in enumerate(adjacency))
    return first_clique(complement, order)


def _edge_expression(instance: CompletionInstance, left: int, right: int,
                     registry: VariableRegistry
                     ) -> tuple[int, tuple[tuple[int, int], ...]]:
    if left > right:
        left, right = right, left
    if left == right:
        raise CompletionViolation("edge expression requires two vertices")
    if right < 5:
        return (int((left - right) % 5 in (1, 4)), ())
    if left < 5:
        orbit, _ = _vertex_orbit(right)
        return ((instance.masks[orbit] >> left) & 1, ())
    left_orbit, left_side = _vertex_orbit(left)
    right_orbit, right_side = _vertex_orbit(right)
    if left_orbit == right_orbit:
        return (instance.internal[left_orbit], ())
    parallel = left_side == right_side
    states = ("10", "11") if parallel else ("01", "11")
    return (0, tuple(
        (1, registry.lookup(state_semantic(left_orbit, right_orbit, state)))
        for state in states
    ))


def ramsey_constraint(instance: CompletionInstance, registry: VariableRegistry,
                      vertices: tuple[int, ...], kind: str
                      ) -> PBConstraint | None:
    if (type(vertices) is not tuple or len(vertices) != 5
            or tuple(sorted(vertices)) != vertices
            or len(set(vertices)) != 5
            or not all(type(vertex) is int and 0 <= vertex < 45
                       for vertex in vertices)):
        raise CompletionViolation("Ramsey witness must be a sorted five-set")
    constant = 0
    terms = []
    for left, right in itertools.combinations(vertices, 2):
        edge_constant, edge_terms = _edge_expression(
            instance, left, right, registry)
        constant += edge_constant
        terms.extend(edge_terms)
    normalized = PBFormula._normalize(terms)
    if kind == "clique":
        maximum = constant + sum(max(0, coefficient)
                                 for coefficient, _ in normalized)
        if maximum <= 9:
            return None
        return PBConstraint(
            tuple((-coefficient, variable)
                  for coefficient, variable in normalized),
            ">=", constant - 9,
            "ramsey-clique:" + ",".join(map(str, vertices)),
        )
    if kind == "independent":
        if constant >= 1:
            return None
        return PBConstraint(
            normalized, ">=", 1 - constant,
            "ramsey-independent:" + ",".join(map(str, vertices)),
        )
    raise CompletionViolation(f"unknown Ramsey obstruction kind {kind}")


def add_ramsey_obstruction(formula: PBFormula, registry: VariableRegistry,
                           instance: CompletionInstance,
                           vertices: tuple[int, ...], kind: str) -> bool:
    constraint = ramsey_constraint(instance, registry, vertices, kind)
    if constraint is None:
        return False
    formula.constraints.append(constraint)
    return True


def add_all_ramsey_constraints(formula: PBFormula, registry: VariableRegistry,
                               instance: CompletionInstance) -> int:
    """Add every distinct nontrivial K5/I5 exclusion; intended for auditing."""
    existing = {(constraint.terms, constraint.comparator, constraint.rhs)
                for constraint in formula.constraints}
    added = 0
    for vertices in itertools.combinations(range(45), 5):
        for kind in ("clique", "independent"):
            constraint = ramsey_constraint(
                instance, registry, vertices, kind)
            if constraint is None:
                continue
            key = (constraint.terms, constraint.comparator, constraint.rhs)
            if key in existing:
                continue
            existing.add(key)
            formula.constraints.append(constraint)
            added += 1
    return added


def decode_states(registry: VariableRegistry,
                  true_variables: Iterable[int]) -> dict[tuple[int, int], str]:
    true_set = set(true_variables)
    states = {}
    for left in range(20):
        for right in range(left + 1, 20):
            selected = [state for state in RELATION_STATES
                        if registry.lookup(state_semantic(left, right, state))
                        in true_set]
            if len(selected) != 1:
                raise CompletionViolation(
                    f"pair {(left, right)} has {len(selected)} true states")
            states[left, right] = selected[0]
    return states


def verify_state_assignment(instance: CompletionInstance,
                            states: Mapping[tuple[int, int], str]) -> dict:
    """Independently replay W/T equations and graph properties for one model."""
    adjacency = reconstruct_graph(instance, states)
    columns = instance.r_columns
    w = [[0] * 20 for _ in range(20)]
    t = [[0] * 20 for _ in range(20)]
    for index in range(20):
        w[index][index] = 1 - 2 * instance.internal[index]
        t[index][index] = -w[index][index]
    for (left, right), state in states.items():
        _, _, w_half, t_half = relation_values(state)
        w[left][right] = w[right][left] = 2 * w_half
        t[left][right] = t[right][left] = 2 * t_half
    for vertex in range(5):
        for orbit in range(20):
            value = sum(Q_FIXED[vertex][other] * columns[orbit][other]
                        for other in range(5))
            value += sum(columns[other][vertex] * w[other][orbit]
                         for other in range(20))
            if value != -1:
                raise CompletionViolation("decoded model fails the cross equation")
    for left in range(20):
        for right in range(20):
            w_square = sum(w[left][index] * w[index][right]
                           for index in range(20))
            gram = sum(columns[left][vertex] * columns[right][vertex]
                       for vertex in range(5))
            if w_square + 2 * gram != 45 * (left == right) - 2:
                raise CompletionViolation("decoded model fails the W square equation")
            t_square = sum(t[left][index] * t[index][right]
                           for index in range(20))
            if t_square != 45 * (left == right):
                raise CompletionViolation("decoded model fails the T square equation")
    clique = first_clique(adjacency, 5)
    independent = first_independent_set(adjacency, 5)
    return {
        "srg_45_22_10_11": is_srg_45_22_10_11(adjacency),
        "first_K5": list(clique) if clique is not None else None,
        "first_I5": list(independent) if independent is not None else None,
        "ramsey_good": clique is None and independent is None,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--include-all-ramsey", action="store_true")
    parser.add_argument("--w-only", action="store_true")
    parser.add_argument("--w-domain-tables", action="store_true")
    args = parser.parse_args(argv)
    instance = load_instance(args.source_index)
    formula, registry = build_signed_formula(
        instance,
        include_t_square=not args.w_only,
        include_w_domain_tables=args.w_domain_tables,
        include_ramsey=args.include_all_ramsey,
    )
    args.output.write_text(formula.to_opb(registry))
    print(json.dumps({
        "source_index": instance.source_index,
        "variables": registry.variable_count,
        "constraints": formula.constraint_count,
        "include_all_ramsey": args.include_all_ramsey,
        "t_square": not args.w_only,
        "w_domain_tables": args.w_domain_tables,
        "output": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
