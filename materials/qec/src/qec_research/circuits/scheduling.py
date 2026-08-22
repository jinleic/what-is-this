"""Correct scheduling of ancilla-mediated stabilizer measurements.

THE CORRECTNESS CRITERION
-------------------------
Ancilla a measures generator g_a by applying controlled-P^{a}_j to each data
qubit j in its support.  For two ancillas a,b sharing data qubit j,

    C-P^{a}_j  C-P^{b}_j  =  C-P^{b}_j  C-P^{a}_j  ,     if [P^a_j, P^b_j] = 0
    C-P^{a}_j  C-P^{b}_j  =  C-P^{b}_j  C-P^{a}_j . CZ_{a,b} ,   if they anticommute

Let  J(a,b) = { j : P^a_j and P^b_j anticommute }.  Since g_a and g_b commute
as Pauli operators, |J(a,b)| is even.  Reordering the schedule to the
canonical "all of a before all of b" order emits one CZ_{a,b} for every
j in J(a,b) at which b acts first.  Because CZ^2 = I the circuit measures the
intended operators **iff**

    #{ j in J(a,b) : a acts before b }  is even,   for every pair (a,b).   (*)

An arbitrary edge colouring does NOT satisfy (*).  Empirically, a Konig
colouring of the Steane code violates (*) on 5 pairs and the resulting
"noiseless" circuit fires 1402 detectors in 200 shots.

LOWER BOUND ON DEPTH
--------------------
Any schedule needs at least
    T >= max( max_a |supp(g_a)| , max_j deg(j) )
two-qubit layers, because a single ancilla (resp. data qubit) can take part in
only one two-qubit gate per layer.  Constraint (*) can push the true optimum
higher; we search for the minimum T that admits a valid schedule.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
from ortools.sat.python import cp_model

from .mixed_stabilizer import PauliSupport, generator_supports

__all__ = ["anticommuting_overlaps", "parity_defects", "depth_lower_bound",
           "cpsat_schedule", "ScheduleResult", "slots_to_layers", "verify_schedule"]


def _anticommute(p: str, q: str) -> bool:
    return p != q


def anticommuting_overlaps(supports: list[PauliSupport]) -> dict[tuple[int, int], list[int]]:
    """J(a,b) for every pair with a nonempty anticommuting overlap."""
    by_qubit: dict[int, list[int]] = {}
    for s in supports:
        for j in s.paulis:
            by_qubit.setdefault(j, []).append(s.index)
    out: dict[tuple[int, int], list[int]] = {}
    for j, cs in by_qubit.items():
        for a, b in itertools.combinations(sorted(cs), 2):
            if _anticommute(supports[a].paulis[j], supports[b].paulis[j]):
                out.setdefault((a, b), []).append(j)
    return out


def parity_defects(supports: list[PauliSupport],
                   slot: dict[tuple[int, int], int]) -> list[tuple[int, int, int, int]]:
    """Pairs violating (*).  Entries are (a, b, |J|, #{a before b})."""
    bad = []
    for (a, b), J in anticommuting_overlaps(supports).items():
        cnt = sum(1 for j in J if slot[(a, j)] < slot[(b, j)])
        if cnt % 2:
            bad.append((a, b, len(J), cnt))
    return bad


def depth_lower_bound(supports: list[PauliSupport], n: int) -> int:
    deg = np.zeros(n, dtype=int)
    for s in supports:
        for j in s.paulis:
            deg[j] += 1
    return max(max((s.weight for s in supports), default=0), int(deg.max()) if n else 0)


def slots_to_layers(slot: dict[tuple[int, int], int], T: int) -> list[list[tuple[int, int]]]:
    layers: list[list[tuple[int, int]]] = [[] for _ in range(T)]
    for (c, j), t in slot.items():
        layers[t].append((c, j))
    return layers


def verify_schedule(supports: list[PauliSupport], slot: dict[tuple[int, int], int]) -> dict:
    """All three structural requirements, checked independently of the solver."""
    per_check: dict[int, set[int]] = {}
    per_qubit: dict[int, set[int]] = {}
    ok_check = ok_qubit = True
    for (c, j), t in slot.items():
        if t in per_check.setdefault(c, set()):
            ok_check = False
        per_check[c].add(t)
        if t in per_qubit.setdefault(j, set()):
            ok_qubit = False
        per_qubit[j].add(t)
    edges = {(s.index, j) for s in supports for j in s.paulis}
    defects = parity_defects(supports, slot)
    return {
        "covers_all_edges": set(slot) == edges,
        "ancilla_conflict_free": ok_check,
        "qubit_conflict_free": ok_qubit,
        "parity_defects": len(defects),
        "valid": (set(slot) == edges and ok_check and ok_qubit and not defects),
        "defect_examples": defects[:8],
    }


@dataclass
class ScheduleResult:
    slot: dict[tuple[int, int], int] | None
    depth: int | None
    status: str
    lower_bound: int
    wall_time_s: float
    verification: dict | None = None


def cpsat_schedule(supports: list[PauliSupport], n: int, *, T: int | None = None,
                   time_limit_s: float = 300.0, workers: int = 12,
                   seed_slot: dict[tuple[int, int], int] | None = None,
                   symmetry_orbits: dict[tuple[int, int], int] | None = None,
                   random_seed: int = 0, deterministic: bool = True,
                   ) -> ScheduleResult:
    """Find a depth-T schedule satisfying ancilla/qubit disjointness AND (*).

    ``random_seed`` selects among the (generally many) valid schedules; with
    ``deterministic=True`` the same seed always returns the same schedule, so a
    benchmark can pin one circuit and reuse it across noise strengths.

    ``symmetry_orbits`` optionally maps each edge to an orbit id; all edges in
    an orbit are forced to share a time slot.  For group-structured codes
    (BB / PBB) the orbits are the monomial 'directions', which collapses the
    model from thousands of variables to a handful.
    """
    import time

    t0 = time.time()
    lb = depth_lower_bound(supports, n)
    T = T or lb
    if T < lb:
        return ScheduleResult(None, None, "BELOW_LOWER_BOUND", lb, 0.0)

    model = cp_model.CpModel()
    edges = [(s.index, j) for s in supports for j in s.paulis]

    if symmetry_orbits is not None:
        orbit_var: dict[int, cp_model.IntVar] = {}
        for e in edges:
            o = symmetry_orbits[e]
            if o not in orbit_var:
                orbit_var[o] = model.new_int_var(0, T - 1, f"o{o}")
        t = {e: orbit_var[symmetry_orbits[e]] for e in edges}
    else:
        t = {e: model.new_int_var(0, T - 1, f"t{e[0]}_{e[1]}") for e in edges}

    by_check: dict[int, list] = {}
    by_qubit: dict[int, list] = {}
    for (c, j) in edges:
        by_check.setdefault(c, []).append(t[(c, j)])
        by_qubit.setdefault(j, []).append(t[(c, j)])
    for v in by_check.values():
        if len(v) > 1:
            model.add_all_different(v)
    for v in by_qubit.values():
        if len(v) > 1:
            model.add_all_different(v)

    for (a, b), J in anticommuting_overlaps(supports).items():
        bs = []
        for j in J:
            bv = model.new_bool_var(f"b{a}_{b}_{j}")
            model.add(t[(a, j)] < t[(b, j)]).only_enforce_if(bv)
            model.add(t[(a, j)] > t[(b, j)]).only_enforce_if(bv.negated())
            bs.append(bv)
        z = model.new_int_var(0, len(bs) // 2, f"z{a}_{b}")
        model.add(sum(bs) == 2 * z)

    if seed_slot:
        for e, v in seed_slot.items():
            if e in t:
                model.add_hint(t[e], v)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = workers
    solver.parameters.random_seed = random_seed
    if deterministic:
        # Without this, CP-SAT returns whichever feasible schedule its worker
        # threads happen to find first.  Different schedules for the same code
        # give materially different logical error rates (a factor of ~2 was
        # measured), so a benchmark that re-solves per noise point is
        # confounded.  Pin the search.
        solver.parameters.num_workers = 1
        solver.parameters.interleave_search = False
    st = solver.solve(model)
    name = solver.status_name(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        slot = {e: int(solver.value(t[e])) for e in edges}
        ver = verify_schedule(supports, slot)
        return ScheduleResult(slot, T, name, lb, time.time() - t0, ver)
    return ScheduleResult(None, None, name, lb, time.time() - t0)


def enumerate_orbit_schedules(supports: list[PauliSupport], n: int, T: int,
                              symmetry_orbits: dict[tuple[int, int], int],
                              *, limit: int = 200000,
                              time_limit_s: float = 600.0) -> tuple[list[dict], bool]:
    """Enumerate ALL valid translation-invariant schedules at depth T.

    Returns ``(schedules, complete)``.  ``complete`` is True only when CP-SAT
    proved the enumeration exhaustive (status OPTIMAL) without hitting ``limit``.
    A uniform sample can then be drawn from the full solution set, which a
    seed-varied search cannot provide.
    """
    model = cp_model.CpModel()
    edges = [(s.index, j) for s in supports for j in s.paulis]
    orbit_var: dict[int, cp_model.IntVar] = {}
    for e in edges:
        o = symmetry_orbits[e]
        if o not in orbit_var:
            orbit_var[o] = model.new_int_var(0, T - 1, f"o{o}")
    t = {e: orbit_var[symmetry_orbits[e]] for e in edges}

    by_check: dict[int, list] = {}
    by_qubit: dict[int, list] = {}
    for (c, j) in edges:
        by_check.setdefault(c, []).append(t[(c, j)])
        by_qubit.setdefault(j, []).append(t[(c, j)])
    for v in by_check.values():
        if len(v) > 1:
            model.add_all_different(v)
    for v in by_qubit.values():
        if len(v) > 1:
            model.add_all_different(v)
    for (a, b), J in anticommuting_overlaps(supports).items():
        bs = []
        for j in J:
            bv = model.new_bool_var(f"b{a}_{b}_{j}")
            model.add(t[(a, j)] < t[(b, j)]).only_enforce_if(bv)
            model.add(t[(a, j)] > t[(b, j)]).only_enforce_if(bv.negated())
            bs.append(bv)
        z = model.new_int_var(0, len(bs) // 2, f"z{a}_{b}")
        model.add(sum(bs) == 2 * z)

    keys = sorted(orbit_var)
    found: list[dict] = []

    class _Collect(cp_model.CpSolverSolutionCallback):
        def __init__(self) -> None:
            super().__init__()
            self.n = 0

        def on_solution_callback(self) -> None:
            found.append({o: self.value(orbit_var[o]) for o in keys})
            self.n += 1
            if self.n >= limit:
                self.stop_search()

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = 1
    cb = _Collect()
    st = solver.solve(model, cb)
    complete = (st == cp_model.OPTIMAL) and (len(found) < limit)
    out = [{e: assign[symmetry_orbits[e]] for e in edges} for assign in found]
    return out, complete
