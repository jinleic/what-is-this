"""CDCL decision procedure for weight-bounded logical existence.

Frozen decision semantics (2026-08-15, before any EXP-036 record existed)
-------------------------------------------------------------------------
``decide_weight_bounded(instance, cap)`` asks: does a *nontrivial logical*
vector of weight <= cap exist?

* ``SAT``   => a witness exists; it is decoded and MUST be re-verified through
  two independent GF(2) paths before any bound may move.  A verified witness
  at weight w certifies the upper bound d <= w.
* ``UNSAT`` => no such vector exists; certifies the lower bound d > cap.
* ``UNDECIDED_BUDGET`` => the conflict budget expired; certifies nothing.

Comparison protocol built on these semantics (EXP-036):

* Weak parent domination (d_parent >= d_pbb) is certified by a verified PBB
  witness at U_b plus parent ``UNSAT`` at cap = U_b - 1.
* Strict reversal (d_pbb > d_parent) is certified by a verified parent
  witness at U_p plus PBB ``UNSAT`` at cap = U_p  (NOT U_p - 1: a tie is
  weak domination, never a reversal, because the parent carries more logical
  qubits at the same length).

Encodings
---------
The unknown vector v has V binary variables (ids 1..V).  Constraints:

* every parity row r: XOR_{j in supp(r)} v_j = 0    (Tseitin XOR chain);
* nontriviality: at least one of the *pairing rows* has odd inner product
  with v — encoded as a plain clause over the XOR-chain output literals;
* weight: one weight variable per qubit group (OR of the group's bits),
  bounded by a sequential-counter cardinality constraint.

For a general symplectic instance over (x|z), parity rows are
``lambda_swap(H)``, pairing rows are ``lambda_swap(logical_basis)`` and the
groups are {j, n+j}.  For one CSS side (classical picture) parity rows are
the opposite-type check matrix, pairing rows are the opposite-type logical
representatives, and every group is the singleton {j}.
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field
from threading import Lock
from typing import Any

import numpy as np

from ..gf2.linalg import (
    matmul,
    nullspace_np,
    rank_bitset,
    rows_to_bitsets,
)
from ..symplectic.core import StabilizerCode, lambda_swap

ENCODING_VERSION = "sat-decide-v2"
# PySAT mutates a process-global Formula variable pool while CNFs are built.
# Residual screen workers run in threads; serialize only construction, not solve.
_CNF_BUILD_LOCK = Lock()

__all__ = [
    "ENCODING_VERSION",
    "DecisionInstance",
    "symplectic_instance",
    "css_side_instance",
    "css_logical_bases",
    "build_decision_cnf",
    "cnf_sha256",
    "decide_weight_bounded",
    "decide_by_sectors",
    "sector_instance",
    "verify_witness_two_paths",
]


def cnf_sha256(cnf) -> str:
    """Canonical hash of the exact clause sequence the solver received."""

    payload = json.dumps(cnf.clauses, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


@dataclass
class DecisionInstance:
    """One weight-bounded logical-existence question, encoding-ready."""

    parity_rows: np.ndarray      # (r, V) — every row must pair to 0 with v
    pairing_rows: np.ndarray     # (m, V) — at least one must pair to 1
    groups: list[list[int]]      # weight groups of 0-based variable indexes
    kind: str                    # "symplectic" | "css_x" | "css_z"
    meta: dict[str, Any] = field(default_factory=dict)
    # Optional sound symmetry break: 0-based coordinates, at least one of
    # which may be assumed set.  For bivariate-bicycle-derived codes the
    # translation group x^a y^b acts as qubit permutations that preserve the
    # stabilizer group, the centralizer, nontriviality and weight, and acts
    # transitively on the lm positions of each block.  Hence any solution can
    # be translated so that block A's index 0 or block B's index 0 carries
    # support, making the single clause satisfiability-preserving.
    symmetry_clause: list[int] | None = None

    @property
    def num_vars(self) -> int:
        return int(self.parity_rows.shape[1])


def symplectic_instance(
    code: StabilizerCode,
    logicals: np.ndarray,
    *,
    block_length: int | None = None,
) -> DecisionInstance:
    """General mixed-stabilizer instance over (x|z) in GF(2)^{2n}.

    ``block_length`` = lm for a bivariate-bicycle-derived code enables the
    sound translation symmetry break (qubit 0 of block A or of block B is
    supported).
    """

    logicals = np.asarray(logicals, dtype=np.uint8) & 1
    n = code.n
    if logicals.shape[1] != 2 * n:
        raise ValueError(f"logical basis width {logicals.shape} != 2n={2 * n}")
    symmetry = None
    if block_length is not None:
        if 2 * int(block_length) != n:
            raise ValueError(f"block_length {block_length} incompatible with n={n}")
        b = int(block_length)
        # qubit 0 supported (x_0 or z_0) OR qubit b supported (x_b or z_b)
        symmetry = [0, n, b, n + b]
    return DecisionInstance(
        parity_rows=lambda_swap(code.H),
        pairing_rows=lambda_swap(logicals),
        groups=[[j, n + j] for j in range(n)],
        kind="symplectic",
        meta={"n": n, "block_length": block_length},
        symmetry_clause=symmetry,
    )


def css_logical_bases(HX: np.ndarray, HZ: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (LX, LZ): X-type and Z-type logical representatives.

    LX is a basis of ker(HZ) modulo rowspace(HX); LZ is a basis of ker(HX)
    modulo rowspace(HZ).  Nontriviality of x in ker(HZ) is equivalent to
    ``exists z in LZ with x . z = 1`` (and symmetrically), because the dot
    pairing between ker(HZ)/row(HX) and ker(HX)/row(HZ) is nondegenerate.
    """

    HX = np.asarray(HX, dtype=np.uint8) & 1
    HZ = np.asarray(HZ, dtype=np.uint8) & 1
    n = HX.shape[1]
    if HZ.shape[1] != n:
        raise ValueError("HX/HZ width mismatch")

    def quotient_basis(kernel_of: np.ndarray, modulo: np.ndarray) -> np.ndarray:
        kernel = nullspace_np(kernel_of)
        base_rows = rows_to_bitsets(modulo)
        base_rank = rank_bitset(base_rows, n)
        chosen: list[np.ndarray] = []
        current = list(base_rows)
        current_rank = base_rank
        for row in kernel:
            candidate = rows_to_bitsets(row[None, :])[0]
            trial = current + [candidate]
            rank = rank_bitset(trial, n)
            if rank > current_rank:
                chosen.append(row)
                current = trial
                current_rank = rank
        if not chosen:
            return np.zeros((0, n), dtype=np.uint8)
        return np.asarray(chosen, dtype=np.uint8)

    LX = quotient_basis(HZ, HX)
    LZ = quotient_basis(HX, HZ)
    if LX.shape[0] != LZ.shape[0]:
        raise RuntimeError(
            f"CSS logical bases disagree: |LX|={LX.shape[0]} |LZ|={LZ.shape[0]}"
        )
    pairing = matmul(LX, LZ.T)
    if rank_bitset(rows_to_bitsets(pairing), pairing.shape[1]) != LX.shape[0]:
        raise RuntimeError("CSS logical pairing is degenerate; bases invalid")
    return LX, LZ


def css_side_instance(
    HX: np.ndarray,
    HZ: np.ndarray,
    side: str,
    *,
    block_length: int | None = None,
) -> DecisionInstance:
    """Instance for one CSS side: d_X (side='x') or d_Z (side='z').

    d(CSS) = min(d_X, d_Z), so ``UNSAT`` on BOTH sides at cap c certifies
    d > c, and a verified witness on either side at weight w certifies
    d <= w.
    """

    LX, LZ = css_logical_bases(HX, HZ)
    n = int(np.asarray(HX).shape[1])
    if side == "x":
        parity, pairing = np.asarray(HZ, dtype=np.uint8) & 1, LZ
    elif side == "z":
        parity, pairing = np.asarray(HX, dtype=np.uint8) & 1, LX
    else:
        raise ValueError(f"side must be 'x' or 'z', got {side!r}")
    symmetry = None
    if block_length is not None:
        if 2 * int(block_length) != n:
            raise ValueError(f"block_length {block_length} incompatible with n={n}")
        symmetry = [0, int(block_length)]
    return DecisionInstance(
        parity_rows=parity,
        pairing_rows=pairing,
        groups=[[j] for j in range(n)],
        kind=f"css_{side}",
        meta={"n": n, "block_length": block_length},
        symmetry_clause=symmetry,
    )


def _build_decision_cnf_unlocked(instance: DecisionInstance, weight_cap: int):
    """Equisatisfiable CNF for ``exists nontrivial logical, weight <= cap``."""

    from pysat.card import CardEnc, EncType
    from pysat.formula import CNF, IDPool

    V = instance.num_vars
    groups = instance.groups
    singleton_groups = all(len(group) == 1 for group in groups)
    weight_var_base = V if singleton_groups else V + len(groups)
    cnf = CNF()
    vpool = IDPool(start_from=weight_var_base + 1)

    def xor_output(literals: list[int]) -> int | None:
        """Tseitin XOR chain; returns the output literal, or None if empty."""

        accumulator: int | None = None
        for literal in literals:
            if accumulator is None:
                accumulator = literal
                continue
            auxiliary = vpool.id()
            cnf.extend(
                [
                    [-accumulator, -literal, -auxiliary],
                    [accumulator, literal, -auxiliary],
                    [accumulator, -literal, auxiliary],
                    [-accumulator, literal, auxiliary],
                ]
            )
            accumulator = auxiliary
        return accumulator

    if instance.symmetry_clause:
        cnf.append([int(j) + 1 for j in instance.symmetry_clause])
    for row in instance.parity_rows:
        output = xor_output([int(j) + 1 for j in np.flatnonzero(row)])
        if output is not None:
            cnf.append([-output])
    pairing_outputs: list[int] = []
    for row in instance.pairing_rows:
        output = xor_output([int(j) + 1 for j in np.flatnonzero(row)])
        if output is not None:
            pairing_outputs.append(output)
    if not pairing_outputs:
        # No pairing row can ever be odd: the instance is trivially UNSAT.
        cnf.append([])
        return cnf
    cnf.append(pairing_outputs)

    if singleton_groups:
        weight_literals = [group[0] + 1 for group in groups]
    else:
        weight_literals = []
        for index, group in enumerate(groups):
            w_literal = V + index + 1
            for member in group:
                cnf.append([-(member + 1), w_literal])
            cnf.append([member + 1 for member in group] + [-w_literal])
            weight_literals.append(w_literal)
    cnf.extend(
        CardEnc.atmost(
            lits=weight_literals,
            bound=weight_cap,
            vpool=vpool,
            encoding=EncType.seqcounter,
        ).clauses
    )
    return cnf


def build_decision_cnf(instance: DecisionInstance, weight_cap: int):
    """Thread-safe construction; preserves the exact canonical clause order."""
    with _CNF_BUILD_LOCK:
        return _build_decision_cnf_unlocked(instance, weight_cap)


def sector_bundle_digest(selected: list[int], sector_hashes: list[str]) -> str:
    """Canonical hash of an ordered sector bundle."""

    return hashlib.sha256(
        json.dumps(
            {
                "mode": "sectors",
                "sectors": list(selected),
                "hashes": list(sector_hashes),
            },
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def decision_cnf_digest(
    instance: DecisionInstance,
    weight_cap: int,
    *,
    decomposition: str = "monolithic",
    sectors: list[int] | None = None,
) -> str:
    """Pure canonical digest of the CNF(s) a decision would receive.

    Shared by the solver path and by stamp validation so a replayed
    decisive UNSAT hash-binds under the SAME encoding mode it was proved in.
    """

    if decomposition == "monolithic":
        return cnf_sha256(build_decision_cnf(instance, weight_cap))
    if decomposition == "sectors":
        total = int(instance.pairing_rows.shape[0])
        selected = list(range(total)) if sectors is None else list(sectors)
        hashes = [
            cnf_sha256(build_decision_cnf(sector_instance(instance, s), weight_cap))
            for s in selected
        ]
        return sector_bundle_digest(selected, hashes)
    raise ValueError(f"unknown decomposition {decomposition!r}")


def sector_instance(instance: DecisionInstance, sector: int) -> DecisionInstance:
    """Restrict nontriviality to one logical functional.

    The full question fixes ``at least one pairing row is odd`` — a
    disjunction the solver must carry through the whole search.  Sector
    ``j`` replaces it with the single parity ``<v, L_j> = 1``.  Since

        (exists v: parity & OR_j <v,L_j>=1 & wt<=c)
          <=>  (exists j: exists v: parity & <v,L_j>=1 & wt<=c),

    the full instance is UNSAT iff every sector is UNSAT, and any sector
    witness is a witness for the full instance.

    A symmetry clause sound for the full disjunction need not preserve a fixed
    sector: translation can move ``L_j`` to a linear combination of pairing
    rows.  Therefore sector instances deliberately drop the parent clause.
    Sector-preserving anchors must be proved and attached to the returned
    instance by the caller.
    """
    rows = instance.pairing_rows
    if not 0 <= sector < rows.shape[0]:
        raise ValueError(f"sector {sector} out of range for {rows.shape[0]} rows")
    return DecisionInstance(
        parity_rows=instance.parity_rows,
        pairing_rows=rows[sector : sector + 1],
        groups=instance.groups,
        kind=f"{instance.kind}:sector{sector}",
        meta={**instance.meta, "sector": sector, "parent_kind": instance.kind},
        symmetry_clause=None,
    )


def decide_by_sectors(
    instance: DecisionInstance,
    weight_cap: int,
    *,
    solver_name: str = "cadical195",
    conflict_budget: int = 0,
    sectors: list[int] | None = None,
) -> dict[str, Any]:
    """Decide the full question by exact decomposition over logical sectors.

    SAT as soon as any sector is SAT (that witness is verified against the
    FULL instance); UNSAT only when every sector is UNSAT; otherwise
    UNDECIDED_BUDGET.  The record's ``cnf_sha256`` binds the ordered list of
    per-sector CNF hashes, so a replay reproduces it exactly.
    """

    total = int(instance.pairing_rows.shape[0])
    selected = list(range(total)) if sectors is None else list(sectors)
    per_sector: list[dict[str, Any]] = []
    sector_hashes: list[str] = []
    witness: dict[str, Any] | None = None
    undecided = False
    started = time.perf_counter()
    for sector in selected:
        sub = sector_instance(instance, sector)
        record = decide_weight_bounded(
            sub,
            weight_cap,
            solver_name=solver_name,
            conflict_budget=conflict_budget,
        )
        sector_hashes.append(record["cnf_sha256"])
        per_sector.append(
            {
                "sector": sector,
                "status": record["status"],
                "cnf_sha256": record["cnf_sha256"],
                "wall_time_s": record["solver"]["wall_time_s"],
                "clauses": record["solver"]["clauses"],
            }
        )
        if record["status"] == "SAT":
            witness = record
            break
        if record["status"] == "UNDECIDED_BUDGET":
            undecided = True
    wall = time.perf_counter() - started
    full_coverage = selected == list(range(total))
    if witness is not None:
        status = "SAT"
    elif undecided:
        status = "UNDECIDED_BUDGET"
    elif full_coverage:
        status = "UNSAT"
    else:
        # Exhausting a strict subset proves nothing about the full question:
        # an unexamined sector may still be satisfiable.  Never call this
        # UNSAT; the distinct status can never pass a proof gate.
        status = "UNSAT_SUBSET"
    digest = sector_bundle_digest(selected, sector_hashes)
    out: dict[str, Any] = {
        "kind": instance.kind,
        "weight_cap": int(weight_cap),
        "status": status,
        "cnf_sha256": digest,
        "encoding_version": ENCODING_VERSION,
        "symmetry_break": False,
        "source_symmetry_break_dropped": bool(instance.symmetry_clause),
        "decomposition": "sectors",
        "solver": {
            "name": f"PySAT {solver_name}",
            "backend": "pysat",
            "conflict_budget": int(conflict_budget),
            "num_sectors_total": total,
            "num_sectors_solved": len(per_sector),
            "per_sector": per_sector,
            "wall_time_s": wall,
        },
    }
    if witness is not None:
        # Re-verify the sector witness against the FULL instance before use.
        vector = np.asarray(witness["vector"], dtype=np.uint8)
        verification = verify_witness_two_paths(instance, vector)
        weight = group_weight(instance, vector)
        if not verification["valid"] or weight > weight_cap:
            raise RuntimeError(
                "sector witness failed full-instance verification: "
                f"{verification!r} weight={weight} cap={weight_cap}"
            )
        out["vector"] = witness["vector"]
        out["weight"] = int(weight)
        out["verification"] = verification
    return out


def decide_weight_bounded(
    instance: DecisionInstance,
    weight_cap: int,
    *,
    solver_name: str = "cadical195",
    conflict_budget: int = 0,
) -> dict[str, Any]:
    """Decide ``exists nontrivial logical of weight <= cap``; verify SAT models.

    Returns a record with ``status`` in {"SAT", "UNSAT", "UNDECIDED_BUDGET"};
    a SAT record carries the decoded ``vector`` plus the two-path
    ``verification`` block and fail-stops if verification fails.
    """

    import pysat
    from pysat.solvers import Solver

    cnf = build_decision_cnf(instance, weight_cap)
    cnf_hash = cnf_sha256(cnf)
    started = time.perf_counter()
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as engine:
        if conflict_budget > 0:
            engine.conf_budget(conflict_budget)
            answer = engine.solve_limited(expect_interrupt=False)
        else:
            answer = engine.solve()
        stats = dict(engine.accum_stats())
        model = engine.get_model() if answer else None
    wall = time.perf_counter() - started
    status = {True: "SAT", False: "UNSAT", None: "UNDECIDED_BUDGET"}[answer]
    record: dict[str, Any] = {
        "kind": instance.kind,
        "weight_cap": int(weight_cap),
        "status": status,
        "cnf_sha256": cnf_hash,
        "encoding_version": ENCODING_VERSION,
        "symmetry_break": bool(instance.symmetry_clause),
        "solver": {
            "name": f"PySAT {solver_name}",
            "backend": "pysat",
            "pysat_version": pysat.__version__,
            "conflict_budget": int(conflict_budget),
            "clauses": len(cnf.clauses),
            "cnf_variables": cnf.nv,
            "stats": {key: int(value) for key, value in stats.items()},
            "wall_time_s": wall,
        },
    }
    if model is not None:
        assignment = {abs(literal): literal > 0 for literal in model}
        vector = np.asarray(
            [1 if assignment.get(j + 1, False) else 0 for j in range(instance.num_vars)],
            dtype=np.uint8,
        )
        verification = verify_witness_two_paths(instance, vector)
        weight = group_weight(instance, vector)
        if not verification["valid"] or weight > weight_cap:
            raise RuntimeError(
                f"SAT model failed independent verification: {verification!r} "
                f"weight={weight} cap={weight_cap}"
            )
        record["vector"] = vector.tolist()
        record["weight"] = int(weight)
        record["verification"] = verification
    return record


def group_weight(instance: DecisionInstance, vector: np.ndarray) -> int:
    vector = np.asarray(vector, dtype=np.uint8) & 1
    return sum(
        1 for group in instance.groups if any(vector[member] for member in group)
    )


def verify_witness_two_paths(
    instance: DecisionInstance, vector: np.ndarray
) -> dict[str, Any]:
    """Re-derive parity and pairing claims through numpy AND bitset paths."""

    vector = np.asarray(vector, dtype=np.uint8) & 1
    V = instance.num_vars
    parity_np = matmul(instance.parity_rows, vector[:, None])[:, 0]
    pairing_np = matmul(instance.pairing_rows, vector[:, None])[:, 0]
    vector_bits = rows_to_bitsets(vector[None, :])[0]
    parity_bits = [
        (row & vector_bits).bit_count() & 1
        for row in rows_to_bitsets(instance.parity_rows)
    ]
    pairing_bits = [
        (row & vector_bits).bit_count() & 1
        for row in rows_to_bitsets(instance.pairing_rows)
    ]
    parity_clean_np = not parity_np.any()
    parity_clean_bits = not any(parity_bits)
    pairing_hit_np = bool(pairing_np.any())
    pairing_hit_bits = bool(any(pairing_bits))
    paths_agree = (
        parity_clean_np == parity_clean_bits
        and pairing_hit_np == pairing_hit_bits
        and list(parity_np) == parity_bits
        and list(pairing_np) == pairing_bits
    )
    return {
        "parity_clean_numpy": parity_clean_np,
        "parity_clean_bitset": parity_clean_bits,
        "pairing_hit_numpy": pairing_hit_np,
        "pairing_hit_bitset": pairing_hit_bits,
        "independent_paths_agree": paths_agree,
        "valid": bool(
            parity_clean_np
            and parity_clean_bits
            and pairing_hit_np
            and pairing_hit_bits
            and paths_agree
        ),
        "num_vars": V,
    }
