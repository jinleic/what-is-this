"""EXP-011: can any PBB code escape its CSS shadow?

Theorem 2 gives, for every PBB code Q, a CSS code Q' on the same n qubits with
the same k and d_Z(Q') >= d(Q).  Q' dominates Q in [[n,k,d]] as soon as

        d_X(Q') >= d(Q).

The only way a PBB code can be genuinely better than its own shadow is
d_X(Q') < d(Q).  This script decides that inequality exactly, for the small
catalogue codes where CP-SAT terminates.

Reported quantities (all with explicit exactness flags):
  d_Z_shadow  = min weight of a nontrivial Z-logical of Q'
              = min weight of a nontrivial pure-Z logical of Q      (Theorem 2)
  d_X_shadow  = min weight of a nontrivial X-logical of Q'
  d_shadow    = min(d_X_shadow, d_Z_shadow)
  escapes     = d_shadow < d(Q)         -> Q beats its own shadow
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import PBBSpec, build_pbb  # noqa: E402
from qec_research.codes.pbb_theory import analyse_pbb, css_shadow  # noqa: E402
from qec_research.distance.exact import exact_distance_css  # noqa: E402
from qec_research.distance.sectors import min_pure_z_logical  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def analyse(r: dict, tl: float, workers: int) -> dict:
    t0 = time.time()
    spec = PBBSpec(r["ell"], r["m"],
                   [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                   [tuple(t) for t in (r["C_terms"] or [])],
                   [tuple(t) for t in (r["D_terms"] or [])])
    st = analyse_pbb(spec)
    HX, SZ = css_shadow(spec)
    res = exact_distance_css(HX, SZ, time_limit_s=tl, workers=workers)
    pbb = build_pbb(spec)
    pz = min_pure_z_logical(pbb, time_limit_s=tl, workers=workers)

    dZ, dX = res["d_Z"], res["d_X"]
    d_shadow = res["d"]
    D = int(r["d"])
    # Theorem 2 consistency: the shadow's Z-distance must equal the PBB's
    # minimum pure-Z logical weight.
    thm2_ok = (dZ is not None and pz.weight is not None and dZ == pz.weight)
    return {
        "code_id": r.get("code_id"), "n": r["n"], "k": r["k"], "d_reported": D,
        "trust_level": r.get("trust_level"), "delta": st.delta,
        "d_Z_shadow": dZ, "d_Z_shadow_exact": res["d_Z_exact"],
        "d_X_shadow": dX, "d_X_shadow_exact": res["d_X_exact"],
        "d_shadow": d_shadow, "d_shadow_exact": res["d_exact"],
        "pbb_min_pure_z": pz.weight, "pbb_min_pure_z_exact": pz.exact,
        "theorem2_dZ_identity_holds": bool(thm2_ok),
        "shadow_dominates": bool(d_shadow is not None and d_shadow >= D),
        "escapes_shadow": bool(d_shadow is not None and res["d_exact"] and d_shadow < D),
        "wall_time_s": round(time.time() - t0, 1),
    }


def _job(a):
    r, tl, w = a
    try:
        return analyse(r, tl, w)
    except Exception as e:
        return {"code_id": r.get("code_id"), "n": r.get("n"), "error": repr(e)}


def main() -> None:
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 72
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    per = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    only_delta_pos = (sys.argv[5].lower() != "all") if len(sys.argv) > 5 else True

    rows = [json.loads(l) for l in open(CATALOG) if json.loads(l)["n"] <= max_n]
    if only_delta_pos:
        keep = []
        for r in rows:
            spec = PBBSpec(r["ell"], r["m"],
                           [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                           [tuple(t) for t in (r["C_terms"] or [])],
                           [tuple(t) for t in (r["D_terms"] or [])])
            if analyse_pbb(spec).delta > 0:
                keep.append(r)
        rows = keep
    print(f"analysing {len(rows)} codes (n<={max_n}, delta>0 only={only_delta_pos})",
          flush=True)

    out = []
    with ProcessPoolExecutor(max_workers=procs) as ex:
        futs = [ex.submit(_job, (r, tl, per)) for r in rows]
        for i, f in enumerate(as_completed(futs)):
            rec = f.result()
            out.append(rec)
            if "error" in rec:
                print(f"[{i+1}/{len(rows)}] {rec['code_id']} ERROR {rec['error'][:90]}", flush=True)
            else:
                tag = "ESCAPES" if rec["escapes_shadow"] else (
                    "dominated" if rec["shadow_dominates"] else "undecided")
                print(f"[{i+1}/{len(rows)}] {rec['code_id']} n={rec['n']} k={rec['k']} "
                      f"d={rec['d_reported']} delta={rec['delta']} | shadow "
                      f"dZ={rec['d_Z_shadow']} dX={rec['d_X_shadow']} d={rec['d_shadow']} "
                      f"exact={rec['d_shadow_exact']} thm2={rec['theorem2_dZ_identity_holds']} "
                      f"-> {tag}", flush=True)

    (OUT / f"exp011_shadow_escape_n{max_n}.json").write_text(json.dumps(out, indent=2))
    good = [x for x in out if "error" not in x]
    print(f"\nSUMMARY over {len(good)} codes")
    print(f"  Theorem 2 d_Z identity holds : "
          f"{sum(1 for x in good if x['theorem2_dZ_identity_holds'])}/{len(good)}")
    print(f"  shadow dominates             : {sum(1 for x in good if x['shadow_dominates'])}")
    print(f"  ESCAPES the shadow           : {sum(1 for x in good if x['escapes_shadow'])}")
    print(f"  undecided                    : "
          f"{sum(1 for x in good if not x['shadow_dominates'] and not x['escapes_shadow'])}")
    esc = [x for x in good if x["escapes_shadow"]]
    for x in esc:
        print(f"    ESCAPER {x['code_id']} d_pbb={x['d_reported']} d_shadow={x['d_shadow']} "
              f"(dX={x['d_X_shadow']}, dZ={x['d_Z_shadow']})")


if __name__ == "__main__":
    main()
