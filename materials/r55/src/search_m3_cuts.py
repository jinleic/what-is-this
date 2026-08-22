#!/usr/bin/env python3
"""Exact affine cuts and exact obstruction witnesses over the m=3 cone.

Setting
-------
`m3_deficiency_cone` writes a finite, sorted set of states (d, a, b), each
carrying the exact integers

    deficiency(s)     = a + b,
    excess_balance(s) = 2 (a + b) - budget2(d),
    g_lo(s) <= g(s) <= g_hi(s).

Every vertex of a hypothetical R(5,5,45) graph occupies exactly one state,
so such a graph induces nonnegative state weights w with

    sum_s w_s = 45,
    sum_s w_s excess_balance(s) = 0,
    sum_s w_s g_lo(s) <= 0 <= sum_s w_s g_hi(s),

because the m=2 excess and the m=3 contribution both sum to zero over the
45 vertices while every vertex g stays inside its own state window. This
module reasons about that finite relaxation and about nothing else.

Upper certificates
------------------
A certificate (alpha, beta, gamma, delta) with gamma, delta >= 0 that
satisfies, at every state s of the cone,

    objective(s) <= alpha + beta excess_balance(s)
                    + gamma g_lo(s) - delta g_hi(s),

gives, after weighting by any admissible w and using the three global
relations,

    sum_v objective <= 45 alpha + beta * 0 + gamma * (<= 0) - delta * (>= 0)
                    <= 45 alpha.

So `45 * alpha` is an exact upper bound for every R(5,5,45) graph, and a
route is accepted only when that exact bound clears the frozen acceptance
edge of its objective.

Lower witnesses
---------------
A witness is an exactly rational nonnegative weighting that satisfies the
same three relations and whose exact objective value falls beyond that
acceptance edge. Since value <= 45 alpha holds for every certificate and
every witness, one such witness proves that no certificate of this frozen
cone can accept the route. It refutes the relaxation only; it says nothing
about R(5,5) itself, and eliminating a route is not a positive result.

Trust boundary
--------------
HiGHS proposes dual coefficients and a sparse primal support; no numerical
value is ever believed. Coefficients are rationalized on a frozen
denominator ladder, alpha is rederived exactly, every proposed support is
re-solved by `Fraction` Gaussian elimination, and every certificate,
witness, and persisted record is re-checked by the exact verifiers here. A
solver status is discovery metadata and never a route verdict.

Frozen registry
---------------
Three routes are searched, in this order: `total_deficiency` (accepted at
an exact bound <= 315), then `degree20_count` and `deficiency_ge8_count`
(accepted at an exact bound < 1, which forces a nonnegative integer count
to be zero). A fourth route, `required_local_family`, is not instantiable:
this repository holds no hash-pinned complete-gluing-cover manifest, so it
is persisted as UNAVAILABLE_NO_COVER_CERTIFICATE with no evidence at all.
No caller may supply an objective, a threshold, or an acceptance flag.

Search CLI
----------
`main` revalidates the Task 2 artifact and rehashes every recorded input,
searches the three frozen routes, records the unavailable fourth, verifies
every record from its own stored fields, adds this source file to the
provenance, and atomically replaces the artifact only after the serialized
bytes parse back and re-verify. `CERTIFICATION_UNRESOLVED` exits nonzero
and can never carry a scientific disposition.
"""

import argparse
import json
import os
import sys
import tempfile
from dataclasses import asdict, dataclass, fields
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Task 2 owns the artifact schema, the exact-integer guard, and the
# provenance/hashing helpers; importing them keeps one owner per fact.
from m3_deficiency_cone import (  # noqa: E402
    CAMPAIGN_ID, DATA_DIRNAME, DISPOSITION, N_TARGET, SCHEMA_VERSION, State,
    _INPUT_KEYS, _STATE_KEYS, _TOP_LEVEL_KEYS, _exact_int, _input_record,
    _is_sha256, _load_json, _relative_path, _require_keys, _sha256,
    _sorted_inputs,
)

VERTICES = N_TARGET               # every vertex occupies exactly one state
EXPECTED_STATES = 3215            # the Task 2 cone, frozen before searching
ANALYSIS_FILE = "higher_identity_m3.json"   # the only artifact this updates

# Numerical discovery knobs, frozen before the first full run.
DENOMINATOR_LIMITS = (100, 1000, 10000, 100000, 1000000)
SUPPORT_TOLERANCE = 1e-9          # a primal weight below this is not support
FALLBACK_SUPPORT = 32             # dual-slack neighbours of the one retry
MAX_SUPPORT = 4                   # 2 equalities + 2 inequalities span a basis

ACCEPTED = "ACCEPTED_EXACT_CUT"
REJECTED = "REJECTED_BY_EXACT_WITNESS"
UNRESOLVED = "CERTIFICATION_UNRESOLVED"
UNAVAILABLE = "UNAVAILABLE_NO_COVER_CERTIFICATE"
ROUTE_STATUSES = (ACCEPTED, REJECTED, UNRESOLVED)

ACCEPTED_CUT = "ACCEPTED_CUT"
NO_CUT = "NO_CUT_IN_FROZEN_CONE"
_DISPOSITIONS = {
    ACCEPTED_CUT: "M3_ACCEPTED_CUT",
    NO_CUT: "M3_NO_CUT_IN_FROZEN_CONE",
    UNRESOLVED: "M3_CERTIFICATION_UNRESOLVED",
}

UNAVAILABLE_ROUTE = "required_local_family"

_SELF_PATH = Path(__file__).resolve()
_SRC_DIR = _SELF_PATH.parent
_MATH_ROOT = _SRC_DIR.parents[1]
_ANALYSIS_PATH = _SRC_DIR.parent / DATA_DIRNAME / ANALYSIS_FILE
_STATE_FIELDS = tuple(field.name for field in fields(State))
_ZERO = Fraction(0)


# -------------------------------------------------------- frozen registry ---


def _total_deficiency(state):
    return state.deficiency


def _degree20_count(state):
    return int(state.d == 20)


def _deficiency_ge8_count(state):
    return int(state.deficiency >= 8)


@dataclass(frozen=True)
class _Objective:
    """One frozen route: its exact per-state value and its acceptance edge.

    `limit` is the single owner of the route threshold. An exact upper
    bound at or under the edge accepts; an exactly achieved value beyond
    the same edge rejects. The two predicates are complementary by
    construction, so no route can be tuned on one side only.
    """

    objective_id: str
    value: object
    limit: Fraction
    inclusive: bool

    def accepts(self, bound):
        """True when this exact upper bound clears the frozen edge."""
        return bound <= self.limit if self.inclusive else bound < self.limit

    def rejects(self, value):
        """True when this exactly achieved value falls beyond the edge."""
        return not self.accepts(value)


_REGISTRY = {
    "total_deficiency": _Objective(
        "total_deficiency", _total_deficiency, Fraction(315), True,
    ),
    "degree20_count": _Objective(
        "degree20_count", _degree20_count, Fraction(1), False,
    ),
    "deficiency_ge8_count": _Objective(
        "deficiency_ge8_count", _deficiency_ge8_count, Fraction(1), False,
    ),
}
ROUTE_IDS = tuple(_REGISTRY)


def _objective(objective_id):
    """The frozen objective of a route; an unregistered ID is rejected."""
    if type(objective_id) is not str or objective_id not in _REGISTRY:
        raise ValueError(f"objective_id {objective_id!r} is not registered")
    return _REGISTRY[objective_id]


# ------------------------------------------------- exact rationals, states ---


def _canonical_fraction(text, label):
    """The value of an already reduced canonical fraction string.

    One round trip rejects every non-canonical spelling at once: padded or
    signed text, decimal strings, unreduced fractions, and `n/1`.
    """
    if type(text) is not str:
        raise ValueError(f"{label} must be a string, got {text!r}")
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        raise ValueError(f"{label} is not a rational: {text!r}") from None
    if str(value) != text:
        raise ValueError(f"{label} {text!r} is not canonical, want {value}")
    return value


def _exact_state(state, index):
    """One immutable cell whose seven coordinates are exact integers."""
    if not isinstance(state, State):
        raise ValueError(f"state {index} is not a State: {state!r}")
    for field in _STATE_FIELDS:
        _exact_int(getattr(state, field), f"state {index} {field}")
    return state


def _require_states(states):
    """A non-empty tuple of exact State cells."""
    cells = tuple(states)
    if not cells:
        raise ValueError("no states to certify")
    for index, state in enumerate(cells):
        _exact_state(state, index)
    return cells


# ------------------------------------------------------ exact certificates ---


@dataclass(frozen=True)
class Certificate:
    """An affine upper certificate; every coefficient is canonical text."""

    objective_id: str
    alpha: str
    beta: str
    gamma: str
    delta: str


@dataclass(frozen=True)
class PrimalWitness:
    """A rational state weighting as sorted (state_index, weight) pairs."""

    objective_id: str
    weights: tuple[tuple[int, str], ...]


_CERTIFICATE_KEYS = tuple(sorted(field.name for field in fields(Certificate)))
_WITNESS_KEYS = tuple(sorted(field.name for field in fields(PrimalWitness)))
_WEIGHT_KEYS = ("state_index", "weight")
_RECORD_KEYS = (
    "certificate", "exact_upper_bound", "exact_witness_value", "objective_id",
    "primal_witness", "route_status",
)


def verify_certificate(states, certificate):
    """Check one affine certificate exactly; return the exact 45*alpha.

    Resolves the objective through the frozen registry, requires canonical
    coefficients with gamma, delta >= 0, and checks the affine inequality
    at every single state. Nothing here trusts how the coefficients were
    found.
    """
    if not isinstance(certificate, Certificate):
        raise ValueError(f"not a Certificate: {certificate!r}")
    objective = _objective(certificate.objective_id)
    alpha = _canonical_fraction(certificate.alpha, "alpha")
    beta = _canonical_fraction(certificate.beta, "beta")
    gamma = _canonical_fraction(certificate.gamma, "gamma")
    delta = _canonical_fraction(certificate.delta, "delta")
    if gamma < 0:
        raise ValueError(f"gamma {gamma} is negative")
    if delta < 0:
        raise ValueError(f"delta {delta} is negative")
    for index, state in enumerate(_require_states(states)):
        slack = (
            alpha
            + beta * state.excess_balance
            + gamma * state.g_lo
            - delta * state.g_hi
            - objective.value(state)
        )
        if slack < 0:
            raise ValueError(
                f"state {index} {state} violates the certificate by {-slack}"
            )
    return VERTICES * alpha


def verify_primal_witness(states, witness):
    """Check one witness exactly; return its exact objective value.

    Requires strictly increasing in-range exact integer indices, canonical
    nonnegative weights, total weight 45, zero aggregate excess balance,
    and aggregate g_lo <= 0 <= g_hi: exactly the relations a real 45-vertex
    graph satisfies.
    """
    if not isinstance(witness, PrimalWitness):
        raise ValueError(f"not a PrimalWitness: {witness!r}")
    objective = _objective(witness.objective_id)
    # Indexed access only: a witness never scans the whole cone.
    cells = states if isinstance(states, (list, tuple)) else tuple(states)
    if not cells:
        raise ValueError("no states to weight")
    if not isinstance(witness.weights, tuple) or not witness.weights:
        raise ValueError("witness carries no weights")
    total = _ZERO
    balance = _ZERO
    low = _ZERO
    high = _ZERO
    value = _ZERO
    previous = None
    for entry in witness.weights:
        if not isinstance(entry, tuple) or len(entry) != 2:
            raise ValueError(f"weight entry {entry!r} is not a pair")
        index, text = entry
        _exact_int(index, "state index")
        if not 0 <= index < len(cells):
            raise ValueError(f"state index {index} is out of range")
        if previous is not None and index <= previous:
            raise ValueError(
                f"state index {index} does not increase past {previous}"
            )
        previous = index
        weight = _canonical_fraction(text, f"weight of state {index}")
        if weight < 0:
            raise ValueError(f"weight {text} of state {index} is negative")
        state = _exact_state(cells[index], index)
        total += weight
        balance += weight * state.excess_balance
        low += weight * state.g_lo
        high += weight * state.g_hi
        value += weight * objective.value(state)
    if total != VERTICES:
        raise ValueError(f"weights sum to {total}, not {VERTICES}")
    if balance != 0:
        raise ValueError(f"aggregate excess balance is {balance}, not zero")
    if low > 0:
        raise ValueError(f"aggregate g_lo is {low}, which exceeds zero")
    if high < 0:
        raise ValueError(f"aggregate g_hi is {high}, which is below zero")
    return value


def _certificate(objective_id, alpha, beta, gamma, delta):
    """One certificate with canonical reduced coefficient strings."""
    return Certificate(
        objective_id, str(alpha), str(beta), str(gamma), str(delta),
    )


def _exact_alpha(states, objective, beta, gamma, delta):
    """The least alpha satisfying the affine inequality at every state."""
    return max(
        Fraction(objective.value(state))
        - beta * state.excess_balance
        - gamma * state.g_lo
        + delta * state.g_hi
        for state in states
    )


# ---------------------------------------------------------- dual discovery ---


def _highs(**problem):
    """One HiGHS dual-simplex proposal; a solver failure is never fatal."""
    try:
        return linprog(method="highs-ds", **problem)
    except (ValueError, TypeError) as exc:
        print(f"linprog proposal failed: {exc}", file=sys.stderr, flush=True)
        return None


def _dual_proposal(states, objective):
    """HiGHS floats for (beta, gamma, delta), or None when it fails.

    Minimizes alpha over `[-1, -balance, -g_lo, g_hi] . x <= -objective(s)`,
    one row per state, with gamma, delta >= 0 and alpha, beta free.
    """
    rows = np.empty((len(states), 4))
    limits = np.empty(len(states))
    for index, state in enumerate(states):
        rows[index, 0] = -1.0
        rows[index, 1] = -state.excess_balance
        rows[index, 2] = -state.g_lo
        rows[index, 3] = state.g_hi
        limits[index] = -objective.value(state)
    result = _highs(
        c=np.array([1.0, 0.0, 0.0, 0.0]),
        A_ub=rows,
        b_ub=limits,
        bounds=[(None, None), (None, None), (0, None), (0, None)],
    )
    if result is None or not result.success or result.x is None:
        return None
    if not np.isfinite(result.x).all():
        return None                   # never rationalize a nonfinite float
    return tuple(float(value) for value in result.x[1:4])


def _dual_candidates(states, objective):
    """Rationalized (beta, gamma, delta) proposals on the frozen ladder."""
    proposal = _dual_proposal(states, objective)
    if proposal is None:
        return ()
    candidates = []
    for limit in DENOMINATOR_LIMITS:
        triple = tuple(
            Fraction(value).limit_denominator(limit) for value in proposal
        )
        if triple[1] < 0 or triple[2] < 0:
            continue                      # gamma and delta must stay dual
        if triple not in candidates:
            candidates.append(triple)
    return tuple(candidates)


def exact_affine_bound(states, objective_id):
    """The best exactly verified affine certificate of one frozen route.

    HiGHS only proposes (beta, gamma, delta): alpha is rederived exactly as
    the largest slack the inequality demands, and every candidate is passed
    through `verify_certificate` before it can win. The trivial (0, 0, 0)
    triple always competes, so a poor rationalization can never persist a
    bound above the constant one; winning on that triple is still never
    evidence that no smaller alpha exists. The exact minimum is taken over
    the rational coefficients, so the result never depends on solver order.
    """
    objective = _objective(objective_id)
    cells = _require_states(states)
    triples = [(_ZERO, _ZERO, _ZERO)]
    for triple in _dual_candidates(cells, objective):
        if triple not in triples:
            triples.append(triple)
    candidates = []
    for beta, gamma, delta in triples:
        alpha = _exact_alpha(cells, objective, beta, gamma, delta)
        certificate = _certificate(objective_id, alpha, beta, gamma, delta)
        verify_certificate(cells, certificate)
        candidates.append(((alpha, beta, gamma, delta), certificate))
    return min(candidates, key=lambda candidate: candidate[0])[1]


# -------------------------------------------------------- primal discovery ---


_ACTIVE_CHOICES = ((False, False), (False, True), (True, False), (True, True))


def _primal_proposal(states, objective):
    """HiGHS positive support and lower-bound reduced costs, or None.

    Maximizes `sum objective(s) x_s` over the relaxation itself, so the
    optimal vertex exposes a sparse support to reconstruct exactly.
    """
    size = len(states)
    values = np.empty(size)
    balance = np.empty(size)
    low = np.empty(size)
    high = np.empty(size)
    for index, state in enumerate(states):
        values[index] = objective.value(state)
        balance[index] = state.excess_balance
        low[index] = state.g_lo
        high[index] = state.g_hi
    result = _highs(
        c=-values,
        A_ub=np.stack((low, -high)),
        b_ub=np.zeros(2),
        A_eq=np.stack((np.ones(size), balance)),
        b_eq=np.array([float(VERTICES), 0.0]),
        bounds=(0, None),
    )
    if result is None or not result.success or result.x is None:
        return None
    support = tuple(np.flatnonzero(result.x > SUPPORT_TOLERANCE).tolist())
    lower = getattr(result, "lower", None)
    marginals = None if lower is None else getattr(lower, "marginals", None)
    slacks = None
    if marginals is not None:
        slacks = np.abs(np.asarray(marginals, dtype=float))
        if slacks.shape != (size,):
            slacks = None
    return support, slacks


def _fallback_support(support, slacks):
    """The one frozen retry set: the support and 32 dual-slack neighbours.

    Ordering by `(absolute_dual_slack, state_index)` is fixed in advance;
    this set is never widened after seeing a result.
    """
    order = sorted(
        range(len(slacks)), key=lambda index: (float(slacks[index]), index),
    )
    return tuple(sorted(set(support) | set(order[:FALLBACK_SUPPORT])))


def _unique_solution(rows, width):
    """The unique consistent nonnegative solution of an exact system.

    `rows` is an augmented matrix over `Fraction`. A free column, an
    inconsistent surplus equation, or a negative coordinate all reject the
    proposed basis outright.
    """
    row = 0
    for column in range(width):
        pivot = next(
            (index for index in range(row, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            return None                  # rank deficit: no unique solution
        rows[row], rows[pivot] = rows[pivot], rows[row]
        scale = rows[row][column]
        rows[row] = [value / scale for value in rows[row]]
        for other in range(len(rows)):
            if other != row and rows[other][column]:
                factor = rows[other][column]
                rows[other] = [
                    value - factor * base
                    for value, base in zip(rows[other], rows[row])
                ]
        row += 1
    for surplus in rows[row:]:
        if surplus[width]:
            return None                  # an active choice contradicts itself
    solution = tuple(rows[index][width] for index in range(width))
    if any(value < 0 for value in solution):
        return None
    return solution


def _solve_support(states, support, active):
    """Exact weights of one proposed support and active set, or None."""
    rows = [
        [Fraction(1)] * len(support) + [Fraction(VERTICES)],
        [Fraction(states[index].excess_balance) for index in support]
        + [_ZERO],
    ]
    if active[0]:
        rows.append(
            [Fraction(states[index].g_lo) for index in support] + [_ZERO]
        )
    if active[1]:
        rows.append(
            [Fraction(states[index].g_hi) for index in support] + [_ZERO]
        )
    return _unique_solution(rows, len(support))


def _better(candidate, best):
    """Frozen preference: larger exact value, then sparser and smaller."""
    if best is None:
        return True
    if candidate[0] != best[0]:
        return candidate[0] > best[0]
    return candidate[1:] < best[1:]


def _best_witness(states, objective, indices):
    """The largest exactly reconstructed witness over one candidate set.

    Enumerates every support of at most four states and every active or
    inactive choice for the two g inequalities. A unique solution needs at
    least as many equations as unknowns, so a support wider than the two
    mandatory equations plus the active choices is skipped unsolved. The
    preference key is built from the emitted nonzero entries alone, so two
    proposed bases that reduce to one witness can never disagree on it.
    """
    best = None
    winner = None
    for size in range(1, MAX_SUPPORT + 1):
        for support in combinations(indices, size):
            for active in _ACTIVE_CHOICES:
                if 2 + active[0] + active[1] < size:
                    continue
                weights = _solve_support(states, support, active)
                if weights is None:
                    continue
                pairs = tuple(
                    (index, weight)
                    for index, weight in zip(support, weights)
                    if weight
                )
                entries = tuple(
                    (index, str(weight)) for index, weight in pairs
                )
                witness = PrimalWitness(objective.objective_id, entries)
                try:
                    value = verify_primal_witness(states, witness)
                except ValueError:
                    continue
                candidate = (
                    value,
                    len(pairs),
                    tuple(index for index, _ in pairs),
                    tuple(weight for _, weight in pairs),
                )
                if _better(candidate, best):
                    best = candidate
                    winner = witness
    return winner


def _primal_witness(states, objective):
    """The best exact witness for one route, or None when none reconstructs.

    The positive support of the HiGHS vertex is tried first. Only if that
    reconstructs nothing is the single frozen fallback set enumerated once;
    without reduced costs that set cannot be formed, and the route stays
    unresolved rather than guessing one.
    """
    proposal = _primal_proposal(states, objective)
    if proposal is None:
        return None
    support, slacks = proposal
    witness = _best_witness(states, objective, support)
    if witness is not None:
        return witness
    if slacks is None:
        return None
    return _best_witness(
        states, objective, _fallback_support(support, slacks),
    )


# ---------------------------------------------------------- search records ---


def _witness_record(witness):
    """The canonical record of one witness; weights stay index sorted."""
    return {
        "objective_id": witness.objective_id,
        "weights": [
            {"state_index": index, "weight": weight}
            for index, weight in witness.weights
        ],
    }


def _certificate_from_record(entry, objective_id):
    """Rebuild one certificate from its stored fields."""
    _require_keys(entry, _CERTIFICATE_KEYS, "certificate")
    if entry["objective_id"] != objective_id:
        raise ValueError(
            f"certificate objective {entry['objective_id']!r} is not the "
            f"record objective {objective_id!r}"
        )
    return Certificate(
        entry["objective_id"], entry["alpha"], entry["beta"], entry["gamma"],
        entry["delta"],
    )


def _weight_pair(entry):
    """One (state_index, weight) pair from its persisted spelling.

    A persisted weight is spelled exactly one way: an object carrying
    `state_index` and `weight`. Every other shape, a two-element array
    included, is rejected before the exact verifier ever sees it.
    """
    _require_keys(entry, _WEIGHT_KEYS, "weight")
    return (entry["state_index"], entry["weight"])


def _witness_from_record(entry, objective_id):
    """Rebuild one witness from its stored fields."""
    _require_keys(entry, _WITNESS_KEYS, "primal witness")
    if entry["objective_id"] != objective_id:
        raise ValueError(
            f"witness objective {entry['objective_id']!r} is not the record "
            f"objective {objective_id!r}"
        )
    weights = entry["weights"]
    if not isinstance(weights, (list, tuple)) or not weights:
        raise ValueError("primal witness carries no weights")
    return PrimalWitness(
        entry["objective_id"], tuple(_weight_pair(item) for item in weights),
    )


def _require_stored_value(text, expected, label):
    """A stored exact value exists exactly when its evidence does."""
    if expected is None:
        if text is not None:
            raise ValueError(f"{label} {text!r} has no evidence")
        return
    if text is None:
        raise ValueError(f"{label} is null but its evidence gives {expected}")
    if _canonical_fraction(text, label) != expected:
        raise ValueError(f"{label} {text!r} is not the verified {expected}")


def _derived_status(objective, bound, value):
    """The route status implied by verified evidence and the frozen edge.

    Every verified witness value is at most `45*alpha` of every verified
    certificate, so an accepted bound and a rejecting witness cannot both
    hold. Contradictory evidence therefore means a broken verifier, and it
    fails loudly instead of choosing a side.
    """
    accepted = bound is not None and objective.accepts(bound)
    rejected = value is not None and objective.rejects(value)
    if accepted and rejected:
        raise ValueError(
            f"{objective.objective_id}: accepted bound {bound} contradicts "
            f"witness value {value}"
        )
    if accepted:
        return ACCEPTED
    if rejected:
        return REJECTED
    return UNRESOLVED


def verify_search_record(states, record):
    """Verify one search record independently; return its route status.

    The stored status decides nothing: every piece of evidence is rebuilt
    and re-verified, each stored exact value must equal the verifier's own
    output, and the status is rederived through the frozen threshold.
    """
    _require_keys(record, _RECORD_KEYS, "search record")
    objective = _objective(record["objective_id"])
    bound = None
    if record["certificate"] is not None:
        bound = verify_certificate(
            states,
            _certificate_from_record(
                record["certificate"], objective.objective_id,
            ),
        )
    _require_stored_value(
        record["exact_upper_bound"], bound, "exact_upper_bound",
    )
    value = None
    if record["primal_witness"] is not None:
        value = verify_primal_witness(
            states,
            _witness_from_record(
                record["primal_witness"], objective.objective_id,
            ),
        )
    _require_stored_value(
        record["exact_witness_value"], value, "exact_witness_value",
    )
    status = record["route_status"]
    if status not in ROUTE_STATUSES:
        raise ValueError(f"route_status {status!r} is not a route status")
    derived = _derived_status(objective, bound, value)
    if status != derived:
        raise ValueError(
            f"route_status {status!r} is not the derived {derived}"
        )
    return derived


def _unavailable_record():
    """Route 4: no hash-pinned complete-gluing-cover manifest exists."""
    return {
        "objective_id": UNAVAILABLE_ROUTE,
        "certificate": None,
        "exact_upper_bound": None,
        "primal_witness": None,
        "exact_witness_value": None,
        "route_status": UNAVAILABLE,
    }


def _verify_unavailable_record(record):
    """The fourth record carries the same keys and no evidence at all."""
    if record != _unavailable_record():
        raise ValueError(
            f"the {UNAVAILABLE_ROUTE} record is not the frozen unavailable "
            f"record: {record!r}"
        )
    return UNAVAILABLE


def _search_route(states, objective_id):
    """Search one frozen route and return its fully verified record.

    The exact upper certificate always exists and is always recorded; the
    exact witness is searched only when that bound fails to accept, and a
    route with no reconstructed witness stays unresolved.
    """
    objective = _objective(objective_id)
    certificate = exact_affine_bound(states, objective_id)
    bound = verify_certificate(states, certificate)
    witness = None
    value = None
    if not objective.accepts(bound):
        witness = _primal_witness(states, objective)
        if witness is not None:
            value = verify_primal_witness(states, witness)
    record = {
        "objective_id": objective_id,
        "certificate": asdict(certificate),
        "exact_upper_bound": str(bound),
        "primal_witness": (
            None if witness is None else _witness_record(witness)
        ),
        "exact_witness_value": None if value is None else str(value),
        "route_status": _derived_status(objective, bound, value),
    }
    verify_search_record(states, record)
    return record


def _verify_searches(states, searches):
    """Verify all four persisted records; return the three route statuses."""
    expected = len(ROUTE_IDS) + 1
    if not isinstance(searches, list) or len(searches) != expected:
        raise ValueError(f"searches must hold exactly {expected} records")
    statuses = []
    for record, objective_id in zip(searches, ROUTE_IDS):
        statuses.append(verify_search_record(states, record))
        if record["objective_id"] != objective_id:
            raise ValueError(
                f"search record {record['objective_id']!r} is out of the "
                f"frozen route order, expected {objective_id!r}"
            )
    _verify_unavailable_record(searches[-1])
    return tuple(statuses)


def _terminal_status(statuses):
    """The aggregate verdict and the number of accepted routes."""
    accepted = sum(status == ACCEPTED for status in statuses)
    if accepted:
        return ACCEPTED_CUT, accepted
    if all(status == REJECTED for status in statuses):
        return NO_CUT, 0
    return UNRESOLVED, 0


# --------------------------------------------------- artifact and full CLI ---


def _verify_inputs(inputs, math_root):
    """Every recorded input is still exactly the bytes Task 2 consumed."""
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("the analysis records no inputs")
    for record in inputs:
        _require_keys(record, _INPUT_KEYS, "input record")
        label = record["relative_path"]
        if type(label) is not str or not label:
            raise ValueError("input record without a relative path")
        if not _is_sha256(record["sha256"]):
            raise ValueError(f"{label}: malformed sha256")
        path = math_root / label
        if _relative_path(path, math_root) != label:
            raise ValueError(f"{label}: not a canonical math-relative path")
        digest = _sha256(path)
        if digest != record["sha256"]:
            raise ValueError(
                f"{label}: sha256 {digest} is not the recorded "
                f"{record['sha256']}"
            )
        size = path.stat().st_size
        if size != _exact_int(record["bytes"], f"{label} bytes"):
            raise ValueError(
                f"{label}: {size} bytes, recorded {record['bytes']}"
            )


def _load_states(records):
    """Immutable State cells; the source annotations are never read."""
    if not isinstance(records, list):
        raise ValueError("the analysis states are not a list")
    if len(records) != EXPECTED_STATES:
        raise ValueError(
            f"{len(records)} state records, expected {EXPECTED_STATES}"
        )
    states = []
    cells = []
    for index, record in enumerate(records):
        _require_keys(record, _STATE_KEYS, f"state {index}")
        state = State(*(
            _exact_int(record[field], f"state {index} {field}")
            for field in _STATE_FIELDS
        ))
        if state.deficiency != state.a + state.b:
            raise ValueError(f"state {index} deficiency disagrees with a + b")
        if state.g_lo > state.g_hi:
            raise ValueError(f"state {index} has an empty g interval")
        states.append(state)
        cells.append((state.d, state.a, state.b))
    if cells != sorted(set(cells)):
        raise ValueError("state records are not sorted and unique")
    return tuple(states)


def _load_analysis(target):
    """The Task 2 artifact, fully revalidated before any search runs."""
    document = _load_json(target)
    _require_keys(document, _TOP_LEVEL_KEYS, "analysis document")
    if document["schema_version"] != SCHEMA_VERSION:
        raise ValueError(
            f"schema_version {document['schema_version']!r}, expected "
            f"{SCHEMA_VERSION}"
        )
    if document["campaign_id"] != CAMPAIGN_ID:
        raise ValueError(f"campaign_id {document['campaign_id']!r}")
    if document["disposition"] != DISPOSITION:
        raise ValueError(
            f"disposition {document['disposition']!r}, expected "
            f"{DISPOSITION}"
        )
    if document["searches"] != []:
        raise ValueError("searches is not empty; refusing to search twice")
    _verify_inputs(document["inputs"], _MATH_ROOT)
    return document, _load_states(document["states"])


def _updated(document, records, terminal):
    """Task 2's document plus the search records and this run's source.

    Every other field is carried through exactly as parsed, so the artifact
    stays byte-semantically identical outside the three touched fields.
    """
    document["inputs"] = _sorted_inputs(
        list(document["inputs"])
        + [_input_record(_SELF_PATH, _MATH_ROOT, None, None)]
    )
    document["searches"] = records
    document["disposition"] = _DISPOSITIONS[terminal]
    return document


def _write_analysis(target, document, states):
    """Serialize canonically, re-verify the parsed bytes, then replace.

    Task 2's writer pins `searches == []` and its own disposition, so this
    run needs its own atomic write with the search-record checks. A failure
    anywhere leaves the original artifact untouched.
    """
    text = json.dumps(document, indent=2, sort_keys=True, ensure_ascii=True)
    parent = target.parent
    if not parent.is_dir():
        raise ValueError(f"{target.name}: no output directory")
    handle, name = tempfile.mkstemp(
        dir=str(parent), prefix=target.name + ".", suffix=".tmp"
    )
    temporary = Path(name)
    try:
        with os.fdopen(handle, "w", encoding="ascii") as stream:
            stream.write(text)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(name, 0o644)   # mkstemp is 0600; artifacts are readable
        parsed = _load_json(temporary)
        if parsed != document:
            raise ValueError("the serialized search does not parse back")
        _verify_searches(states, parsed["searches"])
        os.replace(name, str(target))
    except BaseException:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise


def _analysis_target(value):
    """The single artifact this CLI is allowed to replace."""
    target = Path(value).resolve()
    if target != _ANALYSIS_PATH:
        raise ValueError(
            f"--analysis must be "
            f"{_relative_path(_ANALYSIS_PATH, _MATH_ROOT)}, got {value}"
        )
    return target


def _run(target):
    """Load, certify, witness, verify, and persist the frozen search."""
    document, states = _load_analysis(target)
    records = []
    for objective_id in ROUTE_IDS:
        record = _search_route(states, objective_id)
        records.append(record)
        print(
            f"route {objective_id}: bound={record['exact_upper_bound']} "
            f"witness={record['exact_witness_value']} "
            f"status={record['route_status']}",
            file=sys.stderr, flush=True,
        )
    records.append(_unavailable_record())
    terminal, accepted = _terminal_status(_verify_searches(states, records))
    _write_analysis(target, _updated(document, records, terminal), states)
    return terminal, accepted


def main(argv=None):
    """Search the frozen m=3 cut routes and update the analysis artifact."""
    parser = argparse.ArgumentParser(
        description=(
            "Search the frozen m=3 cut routes with exact affine "
            "certificates and exact primal witnesses, then record every "
            "verified decision in the canonical analysis JSON."
        )
    )
    parser.add_argument(
        "--analysis", required=True,
        help="canonical Task 2 analysis JSON, replaced atomically after "
             "every check",
    )
    args = parser.parse_args(argv)
    try:
        terminal, accepted = _run(_analysis_target(args.analysis))
    except (ValueError, OSError) as exc:
        print(f"m3 cut search failed: {exc}", file=sys.stderr)
        return 1
    print(f"M3 CUT SEARCH COMPLETE status={terminal} accepted={accepted}")
    return 0 if terminal in (ACCEPTED_CUT, NO_CUT) else 1


if __name__ == "__main__":
    sys.exit(main())
