"""EXP-001: certified exact distances for the BB reference instances.

A distance may be called EXACT only when both halves are established:

  * lower bound -- every logical sector proved `INFEASIBLE` below the value;
  * upper bound -- an explicit minimum-weight nontrivial logical operator is
    produced, independently verified, and persisted.

An earlier version of this script seeded the solver with the catalogue value
and only ever searched strictly below it.  When every sector came back
infeasible it returned that seed as `exact=True` with no witness, which proves
only `d >= D`.  The solver now treats `upper_bound` purely as a search cap and
refuses to report `exact` without a witness of matching weight; this script
additionally re-verifies each witness through a second, independent code path
(pure-Python GF(2) bitsets) before writing the certificate.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import BRAVYI_BB, bb_stabilizer, build_bb  # noqa: E402
from qec_research.distance.exact import exact_distance_css  # noqa: E402
from qec_research.gf2.linalg import (  # noqa: E402
    matmul, rows_to_bitsets, row_space_contains)
from qec_research.symplectic.core import StabilizerCode, symplectic_weight  # noqa: E402

OUT = ROOT / "results" / "certificates"
OUT.mkdir(parents=True, exist_ok=True)


def verify_witness(support: list[int], Hc: np.ndarray, Hs: np.ndarray,
                   n: int, kind: str, code: StabilizerCode) -> dict:
    """Independently confirm a claimed minimum-weight logical operator.

    Checks, via the pure-Python bitset path (not the numpy path used by the
    solver's constraint construction):
      1. weight equals the claimed value;
      2. it lies in ker(Hc)  -- i.e. commutes with the opposite-type checks;
      3. it is NOT in rowspace(Hs) -- i.e. not a stabilizer;
      4. as a full symplectic vector it is a nontrivial logical of the code.
    """
    w = np.zeros(n, dtype=np.uint8)
    w[support] = 1
    in_kernel = not matmul(Hc, w[:, None]).any()
    is_stab = row_space_contains(rows_to_bitsets(Hs), n, rows_to_bitsets(w[None, :])[0])
    v = np.zeros(2 * n, dtype=np.uint8)
    if kind == "X":
        v[:n] = w
    else:
        v[n:] = w
    return {
        "weight": int(w.sum()),
        "commutes_with_opposite_checks": bool(in_kernel),
        "is_stabilizer": bool(is_stab),
        "is_nontrivial_logical": bool(code.is_logical(v)),
        "symplectic_weight": symplectic_weight(v),
        "support": sorted(int(j) for j in support),
    }


def run(name: str, search_cap: int | None, time_limit: float, workers: int) -> dict:
    spec = BRAVYI_BB[name]
    HX, HZ = build_bb(spec)
    code = bb_stabilizer(spec)
    n = HX.shape[1]
    t0 = time.time()
    res = exact_distance_css(HX, HZ, time_limit_s=time_limit, workers=workers,
                             upper_bound=search_cap)
    checks: dict = {}
    ok = True
    for tag, Hc, Hs in (("X", HZ, HX), ("Z", HX, HZ)):
        sol = res.get(f"d_{tag}_witness")
        if sol is None:
            checks[tag] = {"witness": None,
                           "note": "no witness produced -- lower bound only"}
            ok = False
            continue
        v = verify_witness(sol, Hc, Hs, n, tag, code)
        v["matches_reported_value"] = (v["weight"] == res[f"d_{tag}"])
        checks[tag] = v
        ok = ok and v["commutes_with_opposite_checks"] and not v["is_stabilizer"] \
            and v["is_nontrivial_logical"] and v["matches_reported_value"]

    certified = bool(res["d_exact"] and ok)
    return {
        "code": name, "ell": spec.ell, "m": spec.m, "A": spec.A, "B": spec.B,
        "n": spec.n, "k": code.k,
        "d_X": res["d_X"], "d_Z": res["d_Z"], "d": res["d"],
        "d_lower_bound": res["d_lower_bound"],
        "solver_reports_exact": res["d_exact"],
        "witnesses_independently_verified": ok,
        "CERTIFIED_EXACT": certified,
        "witness_X": checks.get("X"), "witness_Z": checks.get("Z"),
        "all_sectors_decided": {"X": res["d_X_all_sectors_decided"],
                                "Z": res["d_Z_all_sectors_decided"]},
        "search_cap": search_cap,
        "wall_time_s": round(time.time() - t0, 2),
        "method": ("CP-SAT exact IP per logical sector; lower bound by proven "
                   "INFEASIBLE below the value, upper bound by an explicit "
                   "witness verified through an independent GF(2) bitset path"),
        "solver": "ortools CP-SAT",
        "env": {"python": platform.python_version(), "platform": platform.platform()},
    }


if __name__ == "__main__":
    name = sys.argv[1]
    cap = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != "none" else None
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 3600.0
    wk = int(sys.argv[4]) if len(sys.argv) > 4 else 14
    rec = run(name, cap, tl, wk)
    print(json.dumps({k: v for k, v in rec.items()
                      if k not in ("witness_X", "witness_Z")}, indent=2))
    for tag in ("X", "Z"):
        w = rec[f"witness_{tag}"]
        if w and "support" in w:
            print(f"  witness_{tag}: weight {w['weight']}, nontrivial="
                  f"{w['is_nontrivial_logical']}, support={w['support'][:12]}"
                  f"{'...' if len(w['support']) > 12 else ''}")
    if not rec["CERTIFIED_EXACT"]:
        print("  *** NOT CERTIFIED EXACT -- this is a bound, not a distance ***")
    safe = name.replace("[", "").replace("]", "").replace(",", "_")
    (OUT / f"bb_distance_{safe}.json").write_text(json.dumps(rec, indent=2))
