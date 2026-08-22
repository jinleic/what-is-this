#!/usr/bin/env python3
"""Exact affine cuts and exact obstruction witnesses over the mixed cone.

Setting
-------
`mixed_deficiency_cone` writes a finite, sorted set of states (d, a, b),
each carrying the exact integers

    deficiency(s)     = a + b,
    excess_balance(s) = 2 (a + b) - budget2(d),
    g_lo(s) <= g(s) <= g_hi(s),      (frozen m=2/m=3 rows)
    h_lo(s) <= h(s) <= h_hi(s),      (frozen m=4 row)
    f_lo(s) <= F(s) <= f_hi(s),      (Engstrom mixed K=4 row)

Every vertex of a hypothetical R(5,5,45) graph occupies exactly one state,
so such a graph induces nonnegative state weights w with

    sum_s w_s = 45,
    sum_s w_s excess_balance(s) = 0,
    sum_s w_s g_lo(s) <= 0 <= sum_s w_s g_hi(s),
    sum_s w_s h_lo(s) <= 0 <= sum_s w_s h_hi(s),
    sum_s w_s f_lo(s) <= 0 <= sum_s w_s f_hi(s),

because the m=2 excess, the m=3 contribution, the summed m=4 vertex row
sum_v h_v = 0, and the summed Engstrom mixed row sum_v F_v = 0 all vanish
over the 45 vertices while every vertex g, h and F stays inside its own
state window.  This module reasons about that finite relaxation and about
nothing else.  The mixed identity enters ONLY through the per-state
interval endpoints and these two extra aggregate sign conditions; no motif
coordinate is ever a search feature.

Upper certificates
------------------
A certificate (alpha, beta, gamma, delta, epsilon, zeta, eta, theta) with
gamma, delta, epsilon, zeta, eta, theta >= 0 that satisfies, at every
state s of the cone,

    objective(s) <= alpha + beta excess_balance(s)
                    + gamma g_lo(s) - delta g_hi(s)
                    + epsilon h_lo(s) - zeta h_hi(s)
                    + eta f_lo(s) - theta f_hi(s),

gives, after weighting by any admissible w and using the seven global
relations, `sum_v objective <= 45 alpha`.  So `45 * alpha` is an exact
upper bound for every R(5,5,45) graph, and a route is accepted only when
that exact bound clears the frozen acceptance edge of its objective.  The
two count objectives are zero or one per vertex of an actual graph, so an
exact bound strictly below one forces the count to zero; no integrality
enters the relaxation itself.

Every certificate this module accepts is also canonical: `alpha` must be
the least constant admissible for its own coefficients, so a padded alpha
is rejected exactly like a violated one.

Lower witnesses
---------------
A witness is an exactly rational nonnegative weighting that satisfies the
same seven relations and whose exact objective value falls beyond that
acceptance edge.  Since value <= 45 alpha holds for every certificate and
every witness, one such witness proves that no certificate of this frozen
cone -- the mixed row included -- can accept the route.  It refutes the
relaxation only; it says nothing about R(5,5) itself, and eliminating a
route is not a positive result.

Trust boundary
--------------
HiGHS proposes dual coefficients and a sparse primal support; no numerical
value is ever believed.  Coefficients are rationalized on a frozen
denominator ladder, alpha is rederived exactly, every proposed support is
re-solved by `Fraction` Gaussian elimination, and every certificate,
witness, and persisted record is re-checked by the exact verifiers here.
A solver status is discovery metadata and never a route verdict.

Frozen registry
---------------
The June 2026 sweep registry, identical routes and edges to m=3 and m=4.
Three routes are searched, in this order: `total_deficiency` (accepted at
an exact bound <= 315), then `degree20_count` and `deficiency_ge8_count`
(accepted at an exact bound < 1, which forces a nonnegative integer count
to be zero).  A fourth route, `required_local_family`, is not
instantiable: this repository holds no hash-pinned complete-gluing-cover
manifest, so it is persisted as UNAVAILABLE_NO_COVER_CERTIFICATE with no
evidence at all.  No caller may supply an objective, a threshold, or an
acceptance flag.

Search CLI
----------
`main` revalidates the Task 2 artifact and rehashes every recorded input,
searches the three frozen routes over the eight-coefficient basis, records
the unavailable fourth, verifies every record from its own stored fields,
adds this source file to the provenance, and atomically replaces the
artifact only after the serialized bytes parse back and re-verify.
`MIXED_CERTIFICATION_UNRESOLVED` exits nonzero and can never carry a
scientific disposition.
"""

import argparse
import json
import os
import sys
import tempfile
from dataclasses import asdict, dataclass, fields
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Task 2 owns the artifact schema, the exact-integer guard, and the
# provenance/hashing helpers; importing them keeps one owner per fact.
from mixed_deficiency_cone import (  # noqa: E402
    CAMPAIGN_ID, DATA_DIRNAME, DISPOSITION, N_TARGET, SCHEMA_VERSION, State,
    _TOP_LEVEL_KEYS, _exact_int, _input_record, _load_json, _relative_path,
    _sorted_inputs,
)
from m3_deficiency_cone import (  # noqa: E402
    _INPUT_KEYS, _is_sha256, _require_keys, _sha256,
)

VERTICES = N_TARGET               # every vertex occupies exactly one state
EXPECTED_STATES = 3215            # the Task 2 cone, frozen before searching
ANALYSIS_FILE = "engstrom_identity.json"   # the only artifact this updates

# Numerical discovery knobs, frozen before the first full run.
DENOMINATOR_LIMITS = (100, 1000, 10000, 100000, 1000000)
SUPPORT_TOLERANCE = 1e-9          # a primal weight below this is not support
FALLBACK_SUPPORT = 24             # dual-slack neighbours of the one retry
MAX_SUPPORT = 8                   # 2 equalities + 6 inequalities span a basis

ACCEPTED = "ACCEPTED_EXACT_CUT"
REJECTED = "REJECTED_BY_EXACT_WITNESS"
UNRESOLVED = "CERTIFICATION_UNRESOLVED"
UNAVAILABLE = "UNAVAILABLE_NO_COVER_CERTIFICATE"
ROUTE_STATUSES = (ACCEPTED, REJECTED, UNRESOLVED)

ACCEPTED_CUT = "MIXED_ACCEPTED_CUT"
NO_CUT = "MIXED_NO_CUT_IN_FROZEN_BASIS"
_DISPOSITIONS = {
    ACCEPTED_CUT: "MIXED_ACCEPTED_CUT",
    NO_CUT: "MIXED_NO_CUT_IN_FROZEN_BASIS",
    UNRESOLVED: "MIXED_CERTIFICATION_UNRESOLVED",
}

UNAVAILABLE_ROUTE = "required_local_family"

_SELF_PATH = Path(__file__).resolve()
_SRC_DIR = _SELF_PATH.parent
_MATH_ROOT = _SRC_DIR.parents[1]
_ANALYSIS_PATH = _SRC_DIR.parent / DATA_DIRNAME / ANALYSIS_FILE
_STATE_FIELDS = tuple(field.name for field in fields(State))
# The v3 state record: eleven exact fields plus the four source
# annotations, which are carried through untouched and never read here.
_STATE_RECORD_KEYS = (
    "a", "b", "d", "deficiency", "excess_balance", "f_hi", "f_lo",
    "g_hi", "g_lo", "h_hi", "h_lo", "x_interval_source", "x_motif_source",
    "y_interval_source", "y_motif_source",
)
_DOCUMENT_KEYS = tuple(sorted(_TOP_LEVEL_KEYS))
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
    """One frozen route: its exact per-state value and its acceptance edge."""

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
    if isinstance(objective_id, _Objective):
        return objective_id
    if type(objective_id) is not str or objective_id not in _REGISTRY:
        raise ValueError(f"objective_id {objective_id!r} is not registered")
    return _REGISTRY[objective_id]


# ------------------------------------------------- exact rationals, states ---


def _canonical_fraction(text, label):
    """The value of an already reduced canonical fraction string."""
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
    """One immutable cell whose eleven coordinates are exact integers."""
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
    epsilon: str
    zeta: str
    eta: str
    theta: str


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


def _exact_alpha(states, objective, beta, gamma, delta, epsilon, zeta, eta,
                 theta):
    """The least alpha satisfying the affine inequality at every state."""
    resolved = _objective(objective)
    return max(
        Fraction(resolved.value(state))
        - beta * state.excess_balance
        - gamma * state.g_lo
        + delta * state.g_hi
        - epsilon * state.h_lo
        + zeta * state.h_hi
        - eta * state.f_lo
        + theta * state.f_hi
        for state in _require_states(states)
    )


def verify_certificate(states, certificate):
    """Check one affine certificate exactly; return the exact 45*alpha.

    Resolves the objective through the frozen registry, requires canonical
    coefficients with gamma, delta, epsilon, zeta, eta, theta >= 0, checks
    the affine inequality at every single state, and requires alpha to be
    the least constant its own coefficients admit.  Nothing here trusts
    how the coefficients were found.
    """
    if not isinstance(certificate, Certificate):
        raise ValueError(f"not a Certificate: {certificate!r}")
    objective = _objective(certificate.objective_id)
    alpha = _canonical_fraction(certificate.alpha, "alpha")
    beta = _canonical_fraction(certificate.beta, "beta")
    gamma = _canonical_fraction(certificate.gamma, "gamma")
    delta = _canonical_fraction(certificate.delta, "delta")
    epsilon = _canonical_fraction(certificate.epsilon, "epsilon")
    zeta = _canonical_fraction(certificate.zeta, "zeta")
    eta = _canonical_fraction(certificate.eta, "eta")
    theta = _canonical_fraction(certificate.theta, "theta")
    for name, value in (
        ("gamma", gamma), ("delta", delta), ("epsilon", epsilon),
        ("zeta", zeta), ("eta", eta), ("theta", theta),
    ):
        if value < 0:
            raise ValueError(f"{name} {value} is negative")
    cells = _require_states(states)
    for index, state in enumerate(cells):
        slack = (
            alpha
            + beta * state.excess_balance
            + gamma * state.g_lo
            - delta * state.g_hi
            + epsilon * state.h_lo
            - zeta * state.h_hi
            + eta * state.f_lo
            - theta * state.f_hi
            - objective.value(state)
        )
        if slack < 0:
            raise ValueError(
                f"state {index} {state} violates the certificate by {-slack}"
            )
    least = _exact_alpha(
        cells, objective, beta, gamma, delta, epsilon, zeta, eta, theta)
    if alpha != least:
        raise ValueError(
            f"alpha {alpha} is not the least admissible constant {least}"
        )
    return VERTICES * alpha


def verify_primal_witness(states, witness):
    """Check one witness exactly; return its exact objective value.

    Requires strictly increasing in-range exact integer indices, canonical
    nonnegative weights, total weight 45, zero aggregate excess balance,
    and the six interval aggregates g_lo <= 0 <= g_hi, h_lo <= 0 <= h_hi
    and f_lo <= 0 <= f_hi: exactly the relations a real 45-vertex graph
    satisfies, the Engstrom mixed row included.
    """
    if not isinstance(witness, PrimalWitness):
        raise ValueError(f"not a PrimalWitness: {witness!r}")
    objective = _objective(witness.objective_id)
    cells = states if isinstance(states, (list, tuple)) else tuple(states)
    if not cells:
        raise ValueError("no states to weight")
    if not isinstance(witness.weights, tuple) or not witness.weights:
        raise ValueError("witness carries no weights")
    total = balance = _ZERO
    low = high = h_low = h_high = f_low = f_high = _ZERO
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
        h_low += weight * state.h_lo
        h_high += weight * state.h_hi
        f_low += weight * state.f_lo
        f_high += weight * state.f_hi
        value += weight * objective.value(state)
    if total != VERTICES:
        raise ValueError(f"weights sum to {total}, not {VERTICES}")
    if balance != 0:
        raise ValueError(f"aggregate excess balance is {balance}, not zero")
    if low > 0:
        raise ValueError(f"aggregate g_lo is {low}, which exceeds zero")
    if high < 0:
        raise ValueError(f"aggregate g_hi is {high}, which is below zero")
    if h_low > 0:
        raise ValueError(f"aggregate h_lo is {h_low}, which exceeds zero")
    if h_high < 0:
        raise ValueError(f"aggregate h_hi is {h_high}, which is below zero")
    if f_low > 0:
        raise ValueError(f"aggregate f_lo is {f_low}, which exceeds zero")
    if f_high < 0:
        raise ValueError(f"aggregate f_hi is {f_high}, which is below zero")
    return value


def _certificate(objective_id, alpha, beta, gamma, delta, epsilon, zeta, eta,
                 theta):
    """One certificate with canonical reduced coefficient strings."""
    return Certificate(
        objective_id, str(alpha), str(beta), str(gamma), str(delta),
        str(epsilon), str(zeta), str(eta), str(theta),
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
    """HiGHS floats for the seven free coefficients, or None.

    Minimizes alpha over `[-1, -balance, -g_lo, g_hi, -h_lo, h_hi, -f_lo,
    f_hi] . x <= -objective(s)`, one row per state, with the six interval
    coefficients nonnegative and alpha, beta free.
    """
    rows = np.empty((len(states), 8))
    limits = np.empty(len(states))
    for index, state in enumerate(states):
        rows[index, 0] = -1.0
        rows[index, 1] = -state.excess_balance
        rows[index, 2] = -state.g_lo
        rows[index, 3] = state.g_hi
        rows[index, 4] = -state.h_lo
        rows[index, 5] = state.h_hi
        rows[index, 6] = -state.f_lo
        rows[index, 7] = state.f_hi
        limits[index] = -objective.value(state)
    result = _highs(
        c=np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        A_ub=rows,
        b_ub=limits,
        bounds=[(None, None), (None, None), (0, None), (0, None), (0, None),
                (0, None), (0, None), (0, None)],
    )
    if result is None or not result.success or result.x is None:
        return None
    if not np.isfinite(result.x).all():
        return None                   # never rationalize a nonfinite float
    return tuple(float(value) for value in result.x[1:8])


def _dual_candidates(states, objective):
    """Rationalized seven-coefficient proposals."""
    proposal = _dual_proposal(states, objective)
    if proposal is None:
        return ()
    candidates = []
    for limit in DENOMINATOR_LIMITS:
        septet = tuple(
            Fraction(value).limit_denominator(limit) for value in proposal
        )
        if any(value < 0 for value in septet[1:]):
            continue          # the six interval coefficients stay feasible
        if septet not in candidates:
            candidates.append(septet)
    return tuple(candidates)


def exact_affine_bound(states, objective_id):
    """The best exactly verified affine certificate of one frozen route.

    HiGHS only proposes the seven free coefficients: alpha is rederived
    exactly as the largest slack the inequality demands, and every
    candidate is passed through `verify_certificate` before it can win.
    The trivial all-zero septet always competes, so a poor rationalization
    can never persist a bound above the constant one; winning on that
    septet is still never evidence that no smaller alpha exists.
    """
    objective = _objective(objective_id)
    cells = _require_states(states)
    septets = [(_ZERO,) * 7]
    for septet in _dual_candidates(cells, objective):
        if septet not in septets:
            septets.append(septet)
    candidates = []
    for beta, gamma, delta, epsilon, zeta, eta, theta in septets:
        alpha = _exact_alpha(
            cells, objective, beta, gamma, delta, epsilon, zeta, eta, theta)
        certificate = _certificate(
            objective.objective_id, alpha, beta, gamma, delta, epsilon, zeta,
            eta, theta)
        verify_certificate(cells, certificate)
        candidates.append(
            ((alpha, beta, gamma, delta, epsilon, zeta, eta, theta),
             certificate)
        )
    return min(candidates, key=lambda candidate: candidate[0])[1]


# -------------------------------------------------------- primal discovery ---

# Every active/inactive choice for the six inequality rows
# (g_lo <= 0, g_hi >= 0, h_lo <= 0, h_hi >= 0, f_lo <= 0, f_hi >= 0).
_ACTIVE_CHOICES = tuple(product((False, True), repeat=6))


def _primal_proposal(states, objective):
    """HiGHS positive support and lower-bound reduced costs, or None.

    Maximizes `sum objective(s) x_s` over the relaxation itself -- two
    equalities and six inequalities -- so the optimal vertex exposes a
    sparse support to reconstruct exactly.
    """
    size = len(states)
    values = np.empty(size)
    balance = np.empty(size)
    low = np.empty(size)
    high = np.empty(size)
    h_low = np.empty(size)
    h_high = np.empty(size)
    f_low = np.empty(size)
    f_high = np.empty(size)
    for index, state in enumerate(states):
        values[index] = objective.value(state)
        balance[index] = state.excess_balance
        low[index] = state.g_lo
        high[index] = state.g_hi
        h_low[index] = state.h_lo
        h_high[index] = state.h_hi
        f_low[index] = state.f_lo
        f_high[index] = state.f_hi
    result = _highs(
        c=-values,
        A_ub=np.stack((low, -high, h_low, -h_high, f_low, -f_high)),
        b_ub=np.zeros(6),
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
    """The one frozen retry set: the support and 24 dual-slack neighbours."""
    order = sorted(
        range(len(slacks)), key=lambda index: (float(slacks[index]), index),
    )
    return tuple(sorted(set(support) | set(order[:FALLBACK_SUPPORT])))


def _unique_solution(rows, width):
    """The unique consistent nonnegative solution of an exact system."""
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
    for flag, field in zip(
        active, ("g_lo", "g_hi", "h_lo", "h_hi", "f_lo", "f_hi")
    ):
        if flag:
            rows.append(
                [Fraction(getattr(states[index], field)) for index in support]
                + [_ZERO]
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

    Enumerates every support of at most eight states and every active or
    inactive choice for the six interval inequalities.  A unique solution
    needs at least as many equations as unknowns, so a support wider than
    the two mandatory equations plus the active choices is skipped
    unsolved.
    """
    best = None
    winner = None
    for size in range(1, MAX_SUPPORT + 1):
        for support in combinations(indices, size):
            for active in _ACTIVE_CHOICES:
                if 2 + sum(active) < size:
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
    """The best exact witness for one route, or None when none reconstructs."""
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
        entry["delta"], entry["epsilon"], entry["zeta"], entry["eta"],
        entry["theta"],
    )


def _weight_pair(entry):
    """One (state_index, weight) pair from its persisted spelling."""
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
    """The route status implied by verified evidence and the frozen edge."""
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
    """Verify one search record independently; return its route status."""
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
    """Search one frozen route and return its fully verified record."""
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
        _require_keys(record, _STATE_RECORD_KEYS, f"state {index}")
        state = State(*(
            _exact_int(record[field], f"state {index} {field}")
            for field in _STATE_FIELDS
        ))
        if state.deficiency != state.a + state.b:
            raise ValueError(f"state {index} deficiency disagrees with a + b")
        if state.g_lo > state.g_hi:
            raise ValueError(f"state {index} has an empty g interval")
        if state.h_lo > state.h_hi:
            raise ValueError(f"state {index} has an empty h interval")
        if state.f_lo > state.f_hi:
            raise ValueError(f"state {index} has an empty f interval")
        states.append(state)
        cells.append((state.d, state.a, state.b))
    if cells != sorted(set(cells)):
        raise ValueError("state records are not sorted and unique")
    return tuple(states)


def _load_analysis(target):
    """The Task 2 artifact, fully revalidated before any search runs."""
    document = _load_json(target)
    _require_keys(document, _DOCUMENT_KEYS, "analysis document")
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
    """Task 2's document plus the search records and this run's source."""
    document["inputs"] = _sorted_inputs(
        list(document["inputs"])
        + [_input_record(_SELF_PATH, _MATH_ROOT, None, None)]
    )
    document["searches"] = records
    document["disposition"] = _DISPOSITIONS[terminal]
    return document


def _write_analysis(target, document, states):
    """Serialize canonically, re-verify the parsed bytes, then replace."""
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
    """Search the frozen mixed cut routes and update the analysis artifact."""
    parser = argparse.ArgumentParser(
        description=(
            "Search the frozen mixed cut routes with exact eight-coefficient "
            "affine certificates and exact primal witnesses, then record "
            "every verified decision in the canonical analysis JSON."
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
        print(f"mixed cut search failed: {exc}", file=sys.stderr)
        return 1
    print(f"MIXED CUT SEARCH COMPLETE status={terminal} accepted={accepted}")
    return 0 if terminal in (ACCEPTED_CUT, NO_CUT) else 1


if __name__ == "__main__":
    sys.exit(main())
