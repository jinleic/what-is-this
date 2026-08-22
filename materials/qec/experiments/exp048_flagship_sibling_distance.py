"""EXP-048 — exact distance of the flagship's EXP-046 collapsed sibling.

The EXP-046 constructor for ``12_6_0193`` is the valid perturbation
``C = witness.c_terms, D = 0`` (``build_CD()``: C = poly(cvec), D = zeros,
M = A C^T symmetric; see ``notes/theorem_jc_syzygy_constructor.md``),
whose demoted pure-X witness is the parent row ``x0 = H_X[0]`` of weight 6:
d_X(Q) <= 6.  This script:
  1. builds that exact sibling (D empty) and checks k;
  2. independently verifies the weight-6 pure-X witness lives in
     ker(H_Z(Q)) \\ rowspace(H_X);
  3. runs ``decide_weight_bounded(Q, cap=5)`` on the full-symplectic
     instance: UNSAT certifies d(Q) >= 6 (hence d(Q) = d_X(Q) = 6 exactly —
     a certified 12 -> 6 collapse); SAT instead returns a lower witness and
     deepens the collapse.
Output: atomic JSON record.
"""
from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "processed" / "exp048_flagship_sibling_distance.json"

_spec = importlib.util.spec_from_file_location("e27", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E27)

from qec_research.distance.sat_decide import decide_weight_bounded, symplectic_instance  # noqa: E402


def main() -> None:
    rows = E27.load_catalogue()
    index = next(
        i for i, r in enumerate(rows)
        if r["ell"] == 12 and r["m"] == 6 and E27.catalogue_label(r, i).endswith("0193")
    )
    row = rows[index]
    label = E27.catalogue_label(row, index)

    exp046 = json.loads((ROOT / "results" / "processed" / "exp046_generic_syzygy.json").read_text())
    witness = next(r for r in exp046["rows"] if r["label"] == "12_6_0193")["witness"]

    # The EXP-046 sibling: C = witness.c_terms, D = 0 (build_CD replacement).
    perturbed = dict(row, C_terms=[list(t) for t in witness["c_terms"]], D_terms=[])
    spec, code = E27.pbb_code(perturbed)
    k_q = int(code.n - E27.rank_np(code.H))
    n = code.n

    # Independent verification of the weight-6 pure-X witness: x0 = parent H_X[0]
    # as a full-symplectic vector must be a nontrivial logical of Q (centralizer
    # minus Q's own stabilizer row space, in full (x|z) space).
    _, HX_parent, _ = E27.parent_matrices(row)
    x0 = HX_parent[0].astype(np.uint8)  # (n,) pure-X support
    assert int(x0.sum()) == 6
    x0_full = np.concatenate([x0, np.zeros(n, dtype=np.uint8)])
    witness_valid = bool(code.is_logical(x0_full))
    assert k_q == 12, f"k changed under perturbation: {k_q} != 12"
    assert witness_valid, "independent weight-6 witness verification failed"

    instance = symplectic_instance(code, code.logical_basis(), block_length=72)
    started = time.perf_counter()
    record = decide_weight_bounded(instance, 5, solver_name="cadical195", conflict_budget=0)
    wall = time.perf_counter() - started
    sat_fields = {}
    if record["status"] == "SAT":
        sat_fields = {
            "sat_vector": [int(b) for b in record["vector"]],
            "sat_weight": int(sum((record["vector"][j] or record["vector"][n + j]) for j in range(n))),
            "sat_verification": record.get("verification"),
        }
    payload = {
        "experiment": "exp048_flagship_sibling_distance",
        "parent": label,
        "catalogue_index": index,
        "parent_d_certified": 12,
        "perturbation": {"c_terms": [list(t) for t in witness["c_terms"]], "d_terms": []},
        "n": n,
        "k": k_q,
        "k_parent": 12,
        "witness_x0": {
            "weight": 6,
            "is_logical_full_symplectic": witness_valid,
        },
        "upper_bound": 6 if witness_valid else None,
        "cap": 5,
        "decision": record["status"],
        **sat_fields,
        "cnf_sha256": record["cnf_sha256"],
        "encoding_version": record.get("encoding_version"),
        "solver_wall_s": record["solver"]["wall_time_s"],
        "wall_s": wall,
        "conclusion": (
            "d(Q) = d_X(Q) = 6 exact: certified 12 -> 6 collapse at full k"
            if record["status"] == "UNSAT" and witness_valid
            else f"decision {record['status']} (witness_valid={witness_valid})"
        ),
    }
    tmp = OUT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2))
    tmp.replace(OUT)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
