"""Gate A1 certificate engine: exact min-batch and min-duration solvers.

Certified formulations (pre_statement.md section 5):

* Min-batches = chromatic number chi(G) of the strict conflict graph.
  Z3 decides k-colorability; the k-1 UNSAT + k SAT pair (or the last UNSAT at
  k=chi) is archived verbatim as the optimality proof.
* Min-duration = 0/1 ILP with batch *distance levels*:
    y[b][v] selects level (distance value) v for batch b
    x[i][b] assigns move i to batch b
    sum_v y[b][v] <= 1                    (>= 0: empty batch)
    sum_b x[i][b] = 1                     (every move assigned)
    x[i][b] <= sum_{v: uniq[v] >= d_i} y[b][v]   (level covers distance)
    x[i][b] + x[j][b] <= 1 for (i,j) in E        (compatibility)
  Objective: sum_b sum_v t(uniq[v]) y[b][v], linear in binaries. Every legal
  partition has the feasible y (level = its true max distance), and any
  feasible y costs >= the induced partition's true cost, so the ILP optimum
  equals the true minimum transport duration over ALL legal partitions under
  the pinned duration model.
  Per-moved-atom load/store (2 x 15 us each) is constant and added by the
  caller.

Certificates: solver status and dual bounds are returned in ``log`` /
``dual_bound`` and archived verbatim by the campaign runner.
"""

from __future__ import annotations

from dataclasses import dataclass

from .conflict import Move, conflict_graph


@dataclass
class TransitionCertificate:
    n_moves: int
    chi: int | None
    chi_witness: list[list[int]] | None
    duration_opt: float | None
    duration_batches: list[list[int]] | None
    engine: str
    log: str


def move_len(m: Move) -> float:
    import math

    return math.dist(m.start, m.target)


# ---------------------------------------------------------------- Z3: chi ---
def cert_min_batches(moves: list[Move], ghost_spots: bool = False):
    """Exact chromatic number of the conflict graph via Z3.

    Returns (chi, witness_batches, log): UNSAT logs retained for every k
    below chi, plus the SAT witness at k = chi (pre_statement section 5
    requires the (chi-1)-UNSAT + chi-SAT pair).
    """
    import networkx as nx

    cg = conflict_graph(moves, ghost_spots)
    n = len(moves)
    if n == 0:
        return 0, [], "empty transition: chi=0"
    if not cg:
        return 1, [[i for i in range(n)]], "edgeless conflict graph: chi=1"
    g = nx.Graph()
    g.add_nodes_from(range(n))
    g.add_edges_from(cg)
    ub = max(nx.greedy_color(g, strategy="smallest_last").values()) + 1

    import z3

    log_lines = []
    for k in range(1, ub + 1):
        s = z3.Solver()
        xs = [z3.Int(f"c{i}") for i in range(n)]
        s.add(*[x >= 0 for x in xs], *[x < k for x in xs])
        for i, j in cg:
            s.add(xs[i] != xs[j])
        res = s.check()
        log_lines.append(f"k={k} check={res}")
        if res == z3.sat:
            m = s.model()
            colors = {i: m.eval(xs[i]).as_long() for i in range(n)}
            batches: dict[int, list[int]] = {}
            for i, c in colors.items():
                batches.setdefault(c, []).append(i)
            witness = [batches[c] for c in sorted(batches)]
            log_lines.append(f"chi={k} (sat), witness={witness}")
            log_lines.append(
                "certificate: unsat for every k<chi (see k lines above), sat at chi"
            )
            return k, witness, "\n".join(log_lines)
        assert res == z3.unsat, f"unexpected z3 result {res} at k={k}"
    raise AssertionError("unreachable: DP coloring upper bound must be colorable")


# -------------------------------------------------------- MIP: min duration --
def cert_min_duration(
    moves: list[Move],
    max_batches: int | None = None,
    ghost_spots: bool = False,
    time_limit_s: float = 300.0,
):
    """Exact min transport-duration (batch-duration part) over all legal
    partitions (optionally restricted to <= max_batches batches). Returns
    dict with duration_opt, batches, log, dual_bound, mip_gap.
    """
    import highspy

    cg = conflict_graph(moves, ghost_spots)
    n = len(moves)
    if n == 0:
        return {"duration_opt": 0.0, "batches": [], "log": "empty transition"}
    B = max_batches if max_batches is not None else n
    import math

    uniq = sorted({round(math.dist(m.start, m.target), 9) for m in moves})
    V = len(uniq)

    def yidx(b, v):
        return b * V + v

    xoff = B * V

    def xidx(i, b):
        return xoff + i * B + b

    nvar = B * V + n * B
    from .config import move_batch_duration_us

    cost_v = [move_batch_duration_us(d) for d in uniq]

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("mip_rel_gap", 0.0)
    h.setOptionValue("mip_abs_gap", 0.0)
    h.setOptionValue("time_limit", float(time_limit_s))
    h.addVars(nvar, [0.0] * nvar, [1.0] * nvar)
    h.changeColsIntegrality(
        nvar, list(range(nvar)), [highspy.HighsVarType.kInteger] * nvar
    )
    obj = [0.0] * nvar
    for b in range(B):
        for v in range(V):
            obj[yidx(b, v)] = cost_v[v]
    h.changeColsCost(nvar, list(range(nvar)), obj)
    h.changeObjectiveSense(highspy.ObjSense.kMinimize)
    inf = highspy.kHighsInf
    # 1. sum_v y[b][v] <= 1
    for b in range(B):
        idx = [yidx(b, v) for v in range(V)]
        h.addRow(0.0, 1.0, len(idx), idx, [1.0] * len(idx))
    # 2. sum_b x[i][b] = 1
    for i in range(n):
        idx = [xidx(i, b) for b in range(B)]
        h.addRow(1.0, 1.0, len(idx), idx, [1.0] * len(idx))
    # 3. x[i][b] <= sum_{v: uniq[v] >= d_i} y[b][v]
    for i in range(n):
        di = move_len(moves[i])
        for b in range(B):
            cover = [yidx(b, v) for v in range(V) if uniq[v] >= di - 1e-9]
            assert cover, "every distance must be covered by the largest level"
            idx = [xidx(i, b)] + cover
            coef = [1.0] + [-1.0] * len(cover)
            h.addRow(-inf, 0.0, len(idx), idx, coef)
    # 4. conflicts
    for i, j in cg:
        for b in range(B):
            h.addRow(-inf, 1.0, 2, [xidx(i, b), xidx(j, b)], [1.0, 1.0])
    h.run()
    status = h.getModelStatus()
    info = h.getInfo()
    log = (
        f"highs status={h.modelStatusToString(status)} "
        f"mip_gap={info.mip_gap:.3e} dual_bound={info.mip_dual_bound:.9g} "
        f"obj={info.objective_function_value:.9g} nodes={info.mip_node_count}"
    )
    if h.modelStatusToString(status) != "Optimal":
        return {
            "duration_opt": None,
            "batches": None,
            "log": log,
            "dual_bound": float(info.mip_dual_bound),
            "mip_gap": float(info.mip_gap),
        }
    sol = h.getSolution()
    col = list(sol.col_value)
    batches = []
    for b in range(B):
        batch = [i for i in range(n) if col[xidx(i, b)] > 0.5]
        if batch:
            batches.append(batch)
    assert len(batches) <= B
    # independent re-check: the model's objective must equal the true cost of
    # the reported partition under the pinned duration model
    from .config import move_batch_duration_us

    true_cost = 0.0
    for batch in batches:
        dmax = max(move_len(moves[i]) for i in batch)
        true_cost += move_batch_duration_us(dmax)
    opt = float(info.objective_function_value)
    assert abs(true_cost - opt) < 1e-6 * max(1.0, abs(opt)), (
        f"witness cost {true_cost} != model optimum {opt}"
    )
    return {
        "duration_opt": opt,
        "batches": batches,
        "log": log,
        "dual_bound": float(info.mip_dual_bound),
        "mip_gap": float(info.mip_gap),
    }
