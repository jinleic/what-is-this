"""Integer optimization for a fixed RS-sum witness.

The field equations are represented over the integers with one unrestricted
integer multiple of q per point. HiGHS searches the integer feasible set;
every returned vector is rechecked in F_q before it is accepted.
"""
from __future__ import annotations

from dataclasses import dataclass

import highspy
import numpy as np


@dataclass(frozen=True)
class DeltaMILPResult:
    status: str
    objective: int | None
    components: tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]] | None
    coefficients: tuple[tuple[tuple[int, ...], ...], ...] | None
    solver_message: str


def _status_name(status) -> str:
    return str(status).split(".")[-1]


def delta_milp(inst, M: list[int], time_limit: float = 60.0) -> DeltaMILPResult:
    """Minimize ``sum_i s_i |B_i|`` for a fixed M.

    Variables are integer field representatives ``c[i,l,d]`` in
    ``[0,q-1]``, binary line indicators ``u[i,l]``, and integer equation
    slacks ``k[p]``. For each point p the model has

        sum(i,d) T[i,p,d] c[i,line_i(p),d] - q*k[p] = M[p],

    where T is the exact Lambda-scaled monomial evaluation. Linking
    ``c <= (q-1)u`` makes every used line explicit; the extraction below also
    rejects redundant indicators. The result is solver evidence only until
    the caller performs exact F_q verification and any cheaper-support check.
    """
    q = int(inst.q)
    if len(M) != inst.N:
        raise ValueError(f"M has length {len(M)}, expected {inst.N}")
    M = [int(x) % q for x in M]

    # Variable layout: c variables, then u variables, then one k per point.
    cvar: dict[tuple[int, int, int], int] = {}
    uvar: dict[tuple[int, int], int] = {}
    nvar = 0
    for i in range(3):
        for line in range(inst.num_lines(i)):
            for d in range(inst.t[i]):
                cvar[(i, line, d)] = nvar
                nvar += 1
    for i in range(3):
        for line in range(inst.num_lines(i)):
            uvar[(i, line)] = nvar
            nvar += 1
    kvar = []
    for _ in range(inst.N):
        kvar.append(nvar)
        nvar += 1

    lower = np.zeros(nvar, dtype=np.float64)
    upper = np.zeros(nvar, dtype=np.float64)
    cost = np.zeros(nvar, dtype=np.float64)
    integer = np.zeros(nvar, dtype=np.uint8)
    for v in cvar.values():
        upper[v] = q - 1
        integer[v] = 1
    for (i, _line), v in uvar.items():
        upper[v] = 1
        cost[v] = inst.s[i]
        integer[v] = 1

    # A positive equation contribution is bounded by term_count*(q-1)^2.
    # This finite slack bound is generous and avoids solver-dependent infinity.
    term_count = sum(inst.t)
    K = 1 + (term_count * (q - 1) * (q - 1) + q - 1) // q
    for v in kvar:
        lower[v] = -K
        upper[v] = K
        integer[v] = 1

    rows: list[dict[int, int]] = []
    row_lower: list[float] = []
    row_upper: list[float] = []

    # Exact modular equation rows, one per flat point.
    for p in range(inst.N):
        coord = inst.coord_of(p)
        row: dict[int, int] = {kvar[p]: -q}
        for i in range(3):
            line = inst.line_of_point(i, p)
            ai = coord[i]
            for d in range(inst.t[i]):
                coeff = int(inst.lam[i][ai]) * pow(int(inst.S[i][ai]), d, q) % q
                if coeff:
                    v = cvar[(i, line, d)]
                    row[v] = row.get(v, 0) + coeff
        rows.append(row)
        row_lower.append(float(M[p]))
        row_upper.append(float(M[p]))

    # c[i,l,d] <= (q-1) u[i,l].
    for (i, line, d), cv in cvar.items():
        uv = uvar[(i, line)]
        rows.append({cv: 1, uv: -(q - 1)})
        row_lower.append(-highspy.kHighsInf)
        row_upper.append(0.0)

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("mip_rel_gap", 0.0)
    h.addVars(nvar, lower, upper)
    indices = np.arange(nvar, dtype=np.int32)
    h.changeColsCost(nvar, indices, cost)
    h.changeColsIntegrality(nvar, indices, integer)

    # addRows uses row-wise sparse starts/indices/values.
    starts = np.zeros(len(rows) + 1, dtype=np.int32)
    row_indices: list[int] = []
    row_values: list[float] = []
    for r, row in enumerate(rows):
        for c, a in sorted(row.items()):
            row_indices.append(c)
            row_values.append(float(a))
        starts[r + 1] = len(row_indices)
    h.addRows(
        len(rows),
        np.asarray(row_lower, dtype=np.float64),
        np.asarray(row_upper, dtype=np.float64),
        len(row_indices),
        starts,
        np.asarray(row_indices, dtype=np.int32),
        np.asarray(row_values, dtype=np.float64),
    )
    h.setMinimize()
    run_status = h.run()
    model_status = h.getModelStatus()
    status = _status_name(model_status)
    if status not in {"kOptimal", "kModelOptimal"}:
        return DeltaMILPResult(
            status=status,
            objective=None,
            components=None,
            coefficients=None,
            solver_message=f"run={_status_name(run_status)} model={status}",
        )

    sol = h.getSolution().col_value
    if len(sol) != nvar:
        return DeltaMILPResult(
            status="kNoSolutionVector",
            objective=None,
            components=None,
            coefficients=None,
            solver_message=f"run={_status_name(run_status)} model={status}",
        )

    def rounded(v: float) -> int:
        r = int(round(float(v)))
        if abs(float(v) - r) > 1e-6:
            raise ArithmeticError(f"fractional integer variable {v}")
        return r

    try:
        coeffs: list[list[tuple[int, ...]]] = []
        comps: list[tuple[int, ...]] = []
        for i in range(3):
            ci: list[tuple[int, ...]] = []
            bi: list[int] = []
            for line in range(inst.num_lines(i)):
                vals = tuple(
                    rounded(sol[cvar[(i, line, d)]]) % q
                    for d in range(inst.t[i])
                )
                uv = rounded(sol[uvar[(i, line)]])
                if uv not in (0, 1):
                    raise ArithmeticError(f"non-binary line variable {uv}")
                if any(vals) != (uv == 1):
                    raise ArithmeticError(
                        "line indicator disagrees with coefficient support"
                    )
                if uv:
                    bi.append(line)
                ci.append(vals)
            coeffs.append(ci)
            comps.append(tuple(bi))
        objective = sum(inst.s[i] * len(comps[i]) for i in range(3))
    except ArithmeticError as exc:
        return DeltaMILPResult(
            status="kInvalidSolution",
            objective=None,
            components=None,
            coefficients=None,
            solver_message=str(exc),
        )

    return DeltaMILPResult(
        status=status,
        objective=int(objective),
        components=(comps[0], comps[1], comps[2]),
        coefficients=(tuple(coeffs[0]), tuple(coeffs[1]), tuple(coeffs[2])),
        solver_message=f"run={_status_name(run_status)} model={status}",
    )
