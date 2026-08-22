"""Level-2 (pair-coupled) lift of the n=45 mixed relaxation, and its verdict.

This module answers the falsification criterion recorded at the end of the
mixed-tier campaign (see notes/mixed_campaign_result_2026-08-21.md, section
6.5, and paper/engstrom_mixed_n45.tex, section "Window quality is not the
bottleneck"):

    "if the level-2 lift's exact optimum on `degree20_count` stays above 1,
     the pair-coupled relaxation of this basis is closed as well."

The answer is that it stays at 41.  The pair-coupled relaxation is closed.

Two independent findings are certified here.


1.  The neighbour-degree identity is ALREADY the frozen m=2 row
    ----------------------------------------------------------
For any finite simple graph G on N vertices with e edges, writing
e_x(v) = e(G[N(v)]) and e_D(v) = e(G[V - {v} - N(v)]), pure double counting
gives, for every vertex v,

    (P)   sum_{u in N(v)} deg(u) = e + e_x(v) - e_D(v),

because the edges incident to N(v) split into the d_v edges at v, the
2 e_x(v) endpoints inside N(v), and the e_xz(v) crossing edges, while
e = d_v + e_x(v) + e_D(v) + e_xz(v).  Summing (P) over v and using
sum_v sum_{u in N(v)} deg(u) = sum_u deg(u)^2 yields

    (S)   sum_v deg(v)^2 = N e + sum_v e_x(v) - sum_v e_D(v).

`verify_identity_exhaustive` confirms (P) and (S) on every labelled graph up
to n = 7 and on random graphs at larger orders.

It is tempting to adjoin (S) to the relaxation as a new aggregate row.  It is
not new.  In this repository the z stratum is the COMPLEMENT of the
anti-neighbourhood, so the recorded `e_z = E45[m] - b` is e(Gbar[D(v)]) and the
actual edge count inside D(v) is

    e_D(v) = e_y = C(m, 2) - e_z,      m = 44 - d.

With that substitution the state functional of (S) collapses exactly onto the
frozen m=2 row:

    2 d^2 - 45 d - 2 e_x + 2 e_y  ==  excess_balance(s)      for every state s.

`neighbour_degree_row` computes the left side and `assert_row_is_excess_balance`
checks the equality on all 3215 states.  The exact rank of the frozen basis
{1, excess_balance, g_lo, g_hi, h_lo, h_hi, f_lo, f_hi} over the states is 8
and stays 8 when this row is adjoined.

This closes an entire family of proposed repairs: no handshake-, assortativity-
or degree-squared-style aggregate row can strengthen the relaxation, because
`excess_balance` already IS that row.

    TRAP, recorded deliberately.  Substituting the recorded `e_z` for e_D(v)
    produces `invalid_ez_row`, which IS linearly independent of the frozen
    basis (rank 8 -> 9) and appears to cut all three routes hard.  It is not a
    valid constraint and its "bounds" are meaningless.  `test_pair_lift_mixed`
    pins both facts so the mistake cannot be reintroduced silently.


2.  The level-2 pair lift does not open the routes
    ----------------------------------------------
The lift adds, for a fixed integer edge count e, ordered degree-pair variables
M[d, d'] (adjacent) and Nn[d, d'] (non-adjacent) coupled to the state weights:

    sum_{d'} M[d, d']      = d n_d,          n_d = sum_{s: d_s = d} w_s
    sum_{d'} Nn[d, d']     = (44 - d) n_d
    sum_{d'} M[d, d'] d'   = e n_d + sum_{s: d_s = d} w_s phi(s)
    M, Nn symmetric and nonnegative,         phi(s) = e_x(s) - e_y(s)

The third family is the per-degree-class form of (P); it is strictly finer than
its aggregate, which is the m=2 row.  Each occupied state must also satisfy the
neighbour-degree window 20 d_s <= e + phi(s) <= 24 d_s, since every degree lies
in {20, ..., 24}.

Three codegree facts are adjoined.  Write lam(u, v) = |N(u) cap N(v)|.

  * u ~ v:  N(u) cap N(v) is K3-free (a triangle there plus u and v would be a
    K5) and has independence number <= 4, so lam <= R(3,5) - 1 = 13.  Also
    lam <= min(d_u, d_v) - 1 and, since D(u) cap D(v) has independence number
    <= 3 and is K5-free, |D(u) cap D(v)| = 45 - d_u - d_v + lam <= R(5,4) - 1
    = 24, i.e. lam <= d_u + d_v - 21.
  * u !~ v:  D(u) cap D(v) has independence number <= 2 (three independent
    vertices there plus u and v would be an independent 5-set) and is K5-free,
    so |D(u) cap D(v)| = 43 - d_u - d_v + lam <= R(5,3) - 1 = 13, giving the
    UPPER bound lam <= d_u + d_v - 30.  Also lam <= min(d_u, d_v) and, since
    N(u) cup N(v) avoids u and v, lam >= d_u + d_v - 43.
  * Counting: sum_v e_x(v) = 3T, sum over ordered adjacent pairs of lam = 6T,
    and sum over ordered non-adjacent pairs of lam = 2 (sum_v C(d_v, 2) - 3T).

    The non-adjacent bound is an UPPER bound on lam.  Reading it as a lower
    bound makes every e infeasible, i.e. manufactures a false refutation.

Verdict, exact, maximised over every integer e for which the programme is
feasible (e in [454, 536]):

    total_deficiency       360          (frozen 360)          unchanged
    degree20_count          41          (frozen 14310/349)    -1/349
    deficiency_ge8_count    45          (frozen 45)           unchanged

No route reaches its acceptance edge (315, 1, 1).  NO BOUND ON R(5,5) IS
CLAIMED OR IMPLIED.

Trust discipline, as in the rest of this tier: SciPy/HiGHS only PROPOSES dual
multipliers.  Every accepted number is rationalised and then re-verified in
exact `Fraction` arithmetic by `exact_dual_bound`, which checks dual
feasibility A_eq^T y_eq + A_ub^T y_ub >= c componentwise; any rational
multiplier vector passing that check yields a soundly certified upper bound,
whether or not it is optimal.

Usage:
    python r55/src/pair_lift_mixed.py r55/data/engstrom_identity.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

VERTICES = 45
STATE_DEGREES = (20, 21, 22, 23, 24)

# Published R(4,5) edge maxima per order; owned by data/structural_tables.json
# and re-derived there.  Repeated here only so this module stays importable
# without the producer chain.
E45 = {20: 100, 21: 107, 22: 114, 23: 122, 24: 132}

RAMSEY_3_5 = 14
RAMSEY_4_5 = 25
RAMSEY_5_3 = 14
RAMSEY_5_4 = 25

ROUTES = ("total_deficiency", "degree20_count", "deficiency_ge8_count")
ACCEPTANCE = {
    "total_deficiency": (Fraction(315), "<="),
    "degree20_count": (Fraction(1), "<"),
    "deficiency_ge8_count": (Fraction(1), "<"),
}
FROZEN = {
    "total_deficiency": Fraction(360),
    "degree20_count": Fraction(14310, 349),
    "deficiency_ge8_count": Fraction(45),
}
EXPECTED_LIFT = {
    "total_deficiency": Fraction(360),
    "degree20_count": Fraction(41),
    "deficiency_ge8_count": Fraction(45),
}
EXPECTED_E_WINDOW = (454, 536)

_ROW_FIELDS = (("g_lo", "g_hi"), ("h_lo", "h_hi"), ("f_lo", "f_hi"))


class LiftViolation(Exception):
    """Raised when a claimed fact fails exact re-verification."""


# --------------------------------------------------------------------------
# state algebra
# --------------------------------------------------------------------------

def c2(k: int) -> int:
    return k * (k - 1) // 2


def objective_value(route: str, state: dict) -> int:
    if route == "total_deficiency":
        return state["deficiency"]
    if route == "degree20_count":
        return 1 if state["d"] == 20 else 0
    if route == "deficiency_ge8_count":
        return 1 if state["deficiency"] >= 8 else 0
    raise LiftViolation(f"unknown route {route!r}")


def e_x(state: dict) -> int:
    """Edges inside N(v): the x stratum is the neighbourhood itself."""
    return E45[state["d"]] - state["a"]


def e_y(state: dict) -> int:
    """Edges inside D(v).  The recorded e_z counts the COMPLEMENT on D(v)."""
    m = 44 - state["d"]
    return c2(m) - (E45[m] - state["b"])


def phi(state: dict) -> int:
    """Per-vertex neighbour-degree offset: sum_{u~v} d_u = e + phi(s)."""
    return e_x(state) - e_y(state)


def neighbour_degree_row(state: dict) -> int:
    """State functional of identity (S), doubled to stay integral.

    Equals `excess_balance` on every state; see the module docstring.
    """
    d = state["d"]
    return 2 * d * d - 45 * d - 2 * e_x(state) + 2 * e_y(state)


def invalid_ez_row(state: dict) -> int:
    """The TRAP row: identity (S) with the recorded e_z misread as e(G[D(v)]).

    Linearly independent of the frozen basis but NOT a valid constraint.
    Retained so the regression test can pin the trap.
    """
    d = state["d"]
    m = 44 - d
    return 2 * d * d - 45 * d - 2 * (E45[d] - E45[m]) + 2 * state["a"] - 2 * state["b"]


# --------------------------------------------------------------------------
# codegree facts
# --------------------------------------------------------------------------

def lam_adjacent_max(d: int, dp: int) -> int:
    """u ~ v: K3-free common neighbourhood, plus the D-side order bound."""
    return min(RAMSEY_3_5 - 1, d - 1, dp - 1, d + dp - (VERTICES - RAMSEY_5_4 + 1))


def lam_nonadjacent_max(d: int, dp: int) -> int:
    """u !~ v: |D(u) cap D(v)| = 43 - d - d' + lam <= R(5,3) - 1 = 13."""
    return min(d, dp, d + dp - (VERTICES - 2 - RAMSEY_5_3 + 1))


def lam_nonadjacent_min(d: int, dp: int) -> int:
    """u !~ v: N(u) cup N(v) avoids both u and v, so it has at most 43 members."""
    return max(0, d + dp - (VERTICES - 2))


# --------------------------------------------------------------------------
# identity verification
# --------------------------------------------------------------------------

def _edges_within(adj: list[int], members: int) -> int:
    total = 0
    rest = members
    while rest:
        v = (rest & -rest).bit_length() - 1
        rest &= rest - 1
        total += bin(adj[v] & members).count("1")
    return total // 2


def _identity_holds(adj: list[int], n: int) -> bool:
    full = (1 << n) - 1
    deg = [bin(a).count("1") for a in adj]
    e = sum(deg) // 2
    sum_x = sum_d = 0
    for v in range(n):
        nv = adj[v]
        dv = full & ~nv & ~(1 << v)
        ex = _edges_within(adj, nv)
        ed = _edges_within(adj, dv)
        sum_x += ex
        sum_d += ed
        walk = nv
        nbr_deg = 0
        while walk:
            u = (walk & -walk).bit_length() - 1
            walk &= walk - 1
            nbr_deg += deg[u]
        if nbr_deg != e + ex - ed:                      # identity (P)
            return False
    return sum(x * x for x in deg) == n * e + sum_x - sum_d      # identity (S)


def verify_identity_exhaustive(max_order: int = 7, random_orders=(12, 20, 45),
                               per_order: int = 40, seed: int = 20260821) -> dict:
    """Confirm (P) and (S) exhaustively to `max_order`, then on random graphs."""
    checked = 0
    for n in range(2, max_order + 1):
        pairs = list(itertools.combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            adj = [0] * n
            for bit, (i, j) in enumerate(pairs):
                if mask >> bit & 1:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
            if not _identity_holds(adj, n):
                raise LiftViolation(f"identity fails on n={n} mask={mask}")
            checked += 1
    rng = random.Random(seed)
    randoms = 0
    for n in random_orders:
        for _ in range(per_order):
            p = rng.choice((0.25, 0.5, 0.75))
            adj = [0] * n
            for i, j in itertools.combinations(range(n), 2):
                if rng.random() < p:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
            if not _identity_holds(adj, n):
                raise LiftViolation(f"identity fails on a random graph at n={n}")
            randoms += 1
    return {"exhaustive_graphs": checked, "max_order": max_order,
            "random_graphs": randoms}


def assert_row_is_excess_balance(states: list[dict]) -> int:
    """The neighbour-degree row equals the frozen m=2 row on every state."""
    for state in states:
        if neighbour_degree_row(state) != state["excess_balance"]:
            raise LiftViolation(
                f"neighbour-degree row != excess_balance at "
                f"(d,a,b)=({state['d']},{state['a']},{state['b']})")
    return len(states)


def exact_rank(states: list[dict], extra=None) -> int:
    """Exact rank of the frozen basis over the states, optionally plus a row."""
    cols = [
        lambda s: Fraction(1),
        lambda s: Fraction(s["excess_balance"]),
        lambda s: Fraction(s["g_lo"]),
        lambda s: Fraction(s["g_hi"]),
        lambda s: Fraction(s["h_lo"]),
        lambda s: Fraction(s["h_hi"]),
        lambda s: Fraction(s["f_lo"]),
        lambda s: Fraction(s["f_hi"]),
    ]
    if extra is not None:
        cols = cols + [lambda s: Fraction(extra(s))]
    rows = [[col(s) for col in cols] for s in states]
    width = len(cols)
    pivot = 0
    for c in range(width):
        target = None
        for r in range(pivot, len(rows)):
            if rows[r][c] != 0:
                target = r
                break
        if target is None:
            continue
        rows[pivot], rows[target] = rows[target], rows[pivot]
        head = rows[pivot][c]
        for r in range(pivot + 1, len(rows)):
            if rows[r][c] != 0:
                factor = rows[r][c] / head
                for k in range(c, width):
                    rows[r][k] -= factor * rows[pivot][k]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


# --------------------------------------------------------------------------
# level-2 programme
# --------------------------------------------------------------------------

def build_programme(states: list[dict], route: str, e: int):
    """Exact-Fraction data for the lift: max c.x s.t. Aeq x = beq, Aub x <= bub, x >= 0."""
    keep = [i for i, s in enumerate(states)
            if 20 * s["d"] <= e + phi(s) <= 24 * s["d"]]
    if not keep:
        return None
    degrees = list(STATE_DEGREES)
    pairs = [(i, j) for i in range(5) for j in range(5)]
    nw = len(keep)
    n_pair = len(pairs)
    size = nw + 2 * n_pair
    m_at = lambda p: nw + p                                   # noqa: E731
    n_at = lambda p: nw + n_pair + p                          # noqa: E731
    blank = lambda: [Fraction(0)] * size                      # noqa: E731
    a_eq, b_eq, a_ub, b_ub = [], [], [], []

    row = blank()
    for k in range(nw):
        row[k] = Fraction(1)
    a_eq.append(row)
    b_eq.append(Fraction(VERTICES))

    row = blank()
    for k, i in enumerate(keep):
        row[k] = Fraction(states[i]["excess_balance"])
    a_eq.append(row)
    b_eq.append(Fraction(0))

    row = blank()
    for k, i in enumerate(keep):
        row[k] = Fraction(states[i]["d"])
    a_eq.append(row)
    b_eq.append(Fraction(2 * e))

    for slot, d in enumerate(degrees):
        row = blank()
        for p, (i, _j) in enumerate(pairs):
            if i == slot:
                row[m_at(p)] = Fraction(1)
        for k, i in enumerate(keep):
            if states[i]["d"] == d:
                row[k] = Fraction(-d)
        a_eq.append(row)
        b_eq.append(Fraction(0))

        row = blank()
        for p, (i, _j) in enumerate(pairs):
            if i == slot:
                row[n_at(p)] = Fraction(1)
        for k, i in enumerate(keep):
            if states[i]["d"] == d:
                row[k] = Fraction(-(44 - d))
        a_eq.append(row)
        b_eq.append(Fraction(0))

        row = blank()
        for p, (i, j) in enumerate(pairs):
            if i == slot:
                row[m_at(p)] = Fraction(degrees[j])
        for k, i in enumerate(keep):
            if states[i]["d"] == d:
                row[k] = -(Fraction(e) + Fraction(phi(states[i])))
        a_eq.append(row)
        b_eq.append(Fraction(0))

    for p, (i, j) in enumerate(pairs):
        if i < j:
            q = pairs.index((j, i))
            for at in (m_at, n_at):
                row = blank()
                row[at(p)] = Fraction(1)
                row[at(q)] = Fraction(-1)
                a_eq.append(row)
                b_eq.append(Fraction(0))

    for lo, hi in _ROW_FIELDS:
        row = blank()
        for k, i in enumerate(keep):
            row[k] = Fraction(states[i][lo])
        a_ub.append(row)
        b_ub.append(Fraction(0))
        row = blank()
        for k, i in enumerate(keep):
            row[k] = Fraction(-states[i][hi])
        a_ub.append(row)
        b_ub.append(Fraction(0))

    row = blank()
    for k, i in enumerate(keep):
        row[k] = Fraction(e_x(states[i]))
    a_ub.append(row)
    b_ub.append(Fraction((RAMSEY_3_5 - 1) * e))               # 3T <= 13 e

    for sign, lam in ((Fraction(1), lam_nonadjacent_max),
                      (Fraction(-1), lam_nonadjacent_min)):
        row = blank()
        for k, i in enumerate(keep):
            row[k] = sign * (Fraction(2 * c2(states[i]["d"]))
                             - Fraction(2 * e_x(states[i])))
        for p, (i, j) in enumerate(pairs):
            row[n_at(p)] = -sign * Fraction(lam(degrees[i], degrees[j]))
        a_ub.append(row)
        b_ub.append(Fraction(0))

    cost = blank()
    for k, i in enumerate(keep):
        cost[k] = Fraction(objective_value(route, states[i]))
    return cost, a_eq, b_eq, a_ub, b_ub, size


def exact_dual_bound(states: list[dict], route: str, e: int,
                     denominator: int = 10 ** 7):
    """Propose duals with HiGHS, then certify the bound in exact arithmetic."""
    from scipy.optimize import linprog                        # noqa: PLC0415

    built = build_programme(states, route, e)
    if built is None:
        return None
    cost, a_eq, b_eq, a_ub, b_ub, size = built
    result = linprog(
        [-float(v) for v in cost],
        A_ub=[[float(v) for v in r] for r in a_ub],
        b_ub=[float(v) for v in b_ub],
        A_eq=[[float(v) for v in r] for r in a_eq],
        b_eq=[float(v) for v in b_eq],
        bounds=(0, None), method="highs")
    if not result.success:
        return None
    y_eq = [Fraction(v).limit_denominator(denominator)
            for v in -result.eqlin.marginals]
    y_ub = [max(Fraction(0), Fraction(v).limit_denominator(denominator))
            for v in -result.ineqlin.marginals]
    for col in range(size):
        value = sum(a_eq[k][col] * y_eq[k] for k in range(len(a_eq)))
        value += sum(a_ub[k][col] * y_ub[k] for k in range(len(a_ub)))
        if value - cost[col] < 0:
            raise LiftViolation(
                f"rationalised duals infeasible for {route} at e={e}")
    return (sum(b_eq[k] * y_eq[k] for k in range(len(b_eq)))
            + sum(b_ub[k] * y_ub[k] for k in range(len(b_ub))))


def lift_bounds(states: list[dict], verbose: bool = False) -> dict:
    """Exact level-2 bound per route, maximised over every feasible integer e."""
    lo = VERTICES * min(STATE_DEGREES) // 2
    hi = VERTICES * max(STATE_DEGREES) // 2
    out = {}
    feasible = None
    for route in ROUTES:
        best = None
        seen = []
        for e in range(lo, hi + 1):
            value = exact_dual_bound(states, route, e)
            if value is None:
                continue
            seen.append(e)
            if best is None or value > best:
                best = value
        if best is None:
            raise LiftViolation(f"programme infeasible for every e on {route}")
        out[route] = best
        if feasible is None:
            feasible = (min(seen), max(seen))
        if verbose:
            print(f"  {route:24s} exact level-2 bound = {best}")
    out["_e_window"] = feasible
    return out


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------

def load_states(path: Path) -> list[dict]:
    document = json.loads(path.read_text())
    if document.get("schema_version") != 3:
        raise LiftViolation("expected schema version 3")
    states = document["states"]
    if len(states) != 3215:
        raise LiftViolation(f"expected 3215 states, found {len(states)}")
    return states


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("analysis", type=Path)
    parser.add_argument("--skip-exhaustive", action="store_true",
                        help="skip the n<=7 identity sweep (keeps the state checks)")
    args = parser.parse_args(argv)

    states = load_states(args.analysis)
    print(f"LIFT INPUT: {len(states)} states from {args.analysis}")

    if args.skip_exhaustive:
        print("LIFT IDENTITY: exhaustive sweep SKIPPED by flag")
    else:
        stats = verify_identity_exhaustive()
        print(f"LIFT IDENTITY: (P) and (S) hold on all "
              f"{stats['exhaustive_graphs']} labelled graphs up to n="
              f"{stats['max_order']} and {stats['random_graphs']} random graphs")

    count = assert_row_is_excess_balance(states)
    print(f"LIFT ROW: neighbour-degree row == excess_balance on all {count} "
          f"states -> identity (S) is ALREADY the frozen m=2 row")

    base = exact_rank(states)
    same = exact_rank(states, neighbour_degree_row)
    trap = exact_rank(states, invalid_ez_row)
    if base != 8 or same != 8:
        raise LiftViolation(f"unexpected ranks: base={base} with_row={same}")
    if trap != 9:
        raise LiftViolation(f"trap row rank changed: {trap}")
    print(f"LIFT RANK: frozen basis {base}; + neighbour-degree row {same} "
          f"(no gain); + invalid e_z row {trap} (independent but UNSOUND)")

    bounds = lift_bounds(states, verbose=True)
    window = bounds.pop("_e_window")
    if window != EXPECTED_E_WINDOW:
        raise LiftViolation(f"feasible e window drifted: {window}")
    print(f"LIFT E WINDOW: programme feasible exactly for e in "
          f"[{window[0]}, {window[1]}]")

    opened = []
    for route in ROUTES:
        value = bounds[route]
        if value != EXPECTED_LIFT[route]:
            raise LiftViolation(
                f"{route}: expected {EXPECTED_LIFT[route]}, got {value}")
        edge, op = ACCEPTANCE[route]
        accepted = value <= edge if op == "<=" else value < edge
        if accepted:
            opened.append(route)
        delta = value - FROZEN[route]
        print(f"LIFT ROUTE {route:24s} frozen={FROZEN[route]} "
              f"level-2={value} delta={delta} edge {op}{edge} "
              f"{'ACCEPTED' if accepted else 'REJECTED'}")

    if opened:
        raise LiftViolation(
            "a route was accepted; this needs manual review before any claim: "
            + ", ".join(opened))
    print("LIFT VERDICT: MIXED_PAIR_LIFT_CLOSED "
          "(no route reaches its acceptance edge; no R(5,5) bound claimed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
