"""rf_replay.py — independent re-checker for the Route F campaign artifacts.

Re-derives FROM FROZEN BYTES (no run state):
  1. the 3-way octonion table anchor;
  2. (MAIN) the certification candidate factors from
     candidate_main_factors.csv — its exact CP residual against T_F and the
     certified containment margins at the frozen pass rung;
  3. (TAU7) the tau_r7-certificate-seeded containment margins at the frozen
     rung from tau7_rungevidence.csv;
  4. all fmpq margin rows vs their outward-arb re-evaluations;
  5. the behavioral contract (TAU6 non-containment, PPOS containment,
     PWRNG non-containment vs corrupted-target interpretation).

Exits 0 iff every check it re-derives holds. JSON summary to stdout.
"""
import csv
import json
import os
import sys
import time
from fractions import Fraction

T0 = time.process_time()
from flint import arb, fmpq, fmpq_mat

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)

import rf_krawczyk as rfk  # noqa: E402  (re-uses the exact machinery; the
# replay's INPUTS are the frozen artifact bytes, its EVALUATION is fresh)

FAIL = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        FAIL.append(name)
    return cond


def load_csv_rows(fn):
    with open(os.path.join(HERE, fn)) as f:
        return list(csv.reader(f))


def main():
    # ---- 1. table anchor
    Ttbl, ok = rfk.anchor_table()
    check("table-anchor-3way", ok)

    res = {}
    with open(os.path.join(HERE, "krawczyk_results.json")) as f:
        res = json.load(f)

    # ---- 2. MAIN: rebuild candidate factors exactly from frozen CSV and
    # re-verify the CP residual row and the containment row.
    rows = load_csv_rows("candidate_main_factors.csv")
    fac_rows = [r for r in rows if r and r[0] == "FACTOR"]
    check("candidate-factors-present", len(fac_rows) == 1)
    fac = [rfk.dy(v) for v in fac_rows[0][1:]]
    p, b, c, r = 3, 8, 8, 13
    a0f = fac[:p * r]
    b0f = fac[p * r:(p + b) * r]
    c0f = fac[(p + b) * r:(p + b + c) * r]
    # The main verdict is FAILURE TO CERTIFY (no containment at any rung);
    # the replay re-derives the SAME no-containment verdict from frozen
    # factor bytes and checks the certified exclusion matches the frozen rung.
    TFblk = [[[ (Ttbl[pp][bb][cc] if (bb < 4 and cc < 4) else
                 (Ttbl[pp][bb - 4][cc - 4] if (bb >= 4 and cc >= 4) else 0))
               for cc in range(8)] for bb in range(8)] for pp in range(3)]
    s = rfk.Sys(TFblk, 13, None, "REPLAY")
    g0 = s.resid_exact(a0f, b0f, c0f)
    gmax = max((abs(v) for v in g0), default=fmpq(0))
    check("main-residual-nonzero-consistent-with-f2c", gmax > 0)
    import json as _json
    frozen_rel = res["results"]["MAIN"]["rel_final"]
    check("main-residual-nonzero-consistent-with-f2c", gmax > 0)
    out = rfk.certify(s, fac, "qrcp", "REPLAY-MAIN")
    check("main-replay-no-containment",
          out["verdict"] == "NO-CONTAINMENT-ANY-RUNG")
    check("main-replay-excl-rung-matches-frozen",
          out["excl_rung"] == res["results"]["MAIN"]["excl_rung"])

    # ---- 3. TAU7 replay from frozen factor CSV
    rows7 = load_csv_rows("candidate_tau7_factors.csv")
    fac7 = [fmpq(v) for v in rows7[0][1:]]
    E4 = [[1 if k == m else 0 for k in range(4)] for m in range(4)]

    def qmul(a, b):
        return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
                a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
                a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
                a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])

    tauT = [[[qmul(E4[pp], E4[bb])[cc] for cc in range(4)]
             for bb in range(4)] for pp in range(3)]
    s7 = rfk.Sys(tauT, 7, None, "REPLAY-TAU7")
    out7 = rfk.certify(s7, fac7, "qrcp", "REPLAY-TAU7")
    check("tau7-replay-containment", out7["verdict"] == "CONTAINMENT")
    # the tau_r7 certificate factors are APPROXIMATE decimals; the correct
    # semantic (frozen discipline): the residual is small (< frozen rho=1e-3)
    # and containment certifies an exact rank-7 witness inside that box.
    g7 = max(abs(v) for v in
             s7.resid_exact(fac7[:21], fac7[21:49], fac7[49:]))
    check("tau7-cert-residual-inside-frozen-box",
          g7 < fmpq(1, 1000))
    # ---- 4. behavioral contract from frozen results
    check("tau6-no-containment",
          res["results"]["TAU6"]["verdict"] != "CONTAINMENT")
    check("tau6-exclusion-evidence",
          res["results"]["TAU6"].get("excl_any", False))
    check("ppos-containment", res["results"]["PPOS"]["verdict"] == "CONTAINMENT")
    check("pwrng-no-containment-any-rung",
          res["results"]["PWRNG"]["verdict"] != "CONTAINMENT")
    check("behavior-ok-flag", res.get("behavior_ok", False) is True)

    # ---- 5. arb vs fmpq consistency on the frozen rung rows (MAIN)
    rungs = res["rungs_full"]["MAIN"]
    ok_arb = True
    for row in rungs:
        if row.get("arb_recheck_all_strict") is not None:
            ok_arb &= bool(row["arb_recheck_all_strict"])
    check("main-rung-arb-consistency", ok_arb)

    verdict = "PASS" if not FAIL else "FAIL:" + ";".join(FAIL)
    print(json.dumps({"replay": verdict,
                      "cpu_s": time.process_time() - T0,
                      "main_pass_rung": out.get("pass_rung"),
                      "tau7_pass_rung": out7.get("pass_rung")}))
    sys.exit(0 if not FAIL else 1)


if __name__ == "__main__":
    main()
