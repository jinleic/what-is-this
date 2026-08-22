"""Standalone replay of proofs/shell_stress.md: uniform shell-defect stress test.

Recomputes the PM-S comparison from first-party stored data only:
results/bounds/paired_momentum.json (certified intervals) and
results/shell_stress.json (this front's per-mode table).  Asserts: (1) the four
recorded K*=1/4 min shell defects are reproduced exactly from the source artifact;
(2) all are negative and far below eps*(1/4,2); (3) the structural aggravation
eps*(1/4,2) > 1 - 2K*/I3 holds on certified intervals; (4) the incumbent-endpoint
largest-box shell minimum is below eps*(63/250,2).
Run: PYTHONPATH=src python tests/test_shell_stress.py
"""
import json
import os
import sys
from fractions import Fraction

sys.set_int_max_str_digits(300000)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    src = json.load(open(os.path.join(ROOT, "results", "bounds",
                                      "paired_momentum.json")))
    own = json.load(open(os.path.join(ROOT, "results", "shell_stress.json")))
    ok = []

    def check(name, cond):
        ok.append((name, bool(cond)))
        print(("  ok: " if cond else "FAIL: ") + name)

    eps = {k: [Fraction(v) for v in vv["eps_star_interval"]]
           for k, vv in src["data"]["shell_thresholds_eps_star"].items()}

    # 1-2. recompute K*=1/4 min shell defects per box, exact equality with our table
    recomputed = {}
    for L in src["data"]["lattices"]:
        ev = L["evaluations"]["K=1/4"]
        Khi = Fraction(ev["K_interval"][1])
        best = None
        for md in ev["modes"]:
            lam = Fraction(md["lambda_T"])
            if lam < 2:
                continue
            defect_lo = 1 - 2 * Khi * lam * Fraction(md["ghat"][1])
            best = defect_lo if best is None or defect_lo < best else best
        recomputed[L["name"]] = best
    recorded = {k: Fraction(v) for k, v in
                own["data"]["stress_summary"]["K=1/4,Lambda=2"]
                ["min_shell_defect_by_box"].items()}
    check("four boxes recomputed == recorded", recomputed == recorded)
    check("all four negative", all(v < 0 for v in recomputed.values()))
    check("largest-box minimum is -0.020384-class",
          -Fraction(3, 100) < recomputed["3x3x2"] < -Fraction(1, 100))
    e142 = eps["K=1/4,Lambda=2"][0]
    check("every box violates eps*(1/4,2)", all(v < e142 for v in recomputed.values()))

    # 3. structural aggravation: eps*(1/4,2) exceeds the Gaussian saturant defect
    I3lo, I3hi = [Fraction(v) for v in src["data"]["certified_inputs"]["I3_interval"]]
    sat_defect_lo = 1 - 2 * Fraction(1, 4) / I3lo
    sat_defect_hi = 1 - 2 * Fraction(1, 4) / I3hi
    check("Gaussian saturant defect < eps*(1/4,2) certified",
          sat_defect_hi < eps["K=1/4,Lambda=2"][0] and sat_defect_lo > 0)
    print(f"  info: saturant defect in [{float(sat_defect_lo):.6f},"
          f" {float(sat_defect_hi):.6f}] vs eps* {float(e142):.6f}")

    # 4. 0.252-target: incumbent-endpoint 3x3x2 shell minimum < eps*(63/250,2)
    ev = next(L for L in src["data"]["lattices"]
              if L["name"] == "3x3x2")["evaluations"]["incumbent"]
    Khi = Fraction(ev["K_interval"][1])
    best = min(1 - 2 * Khi * Fraction(md["lambda_T"]) * Fraction(md["ghat"][1])
               for md in ev["modes"] if Fraction(md["lambda_T"]) >= 2)
    check("incumbent 3x3x2 shell min defect < eps*(63/250,2)",
          best < eps["K=63/250,Lambda=2"][0])
    check("the recorded note of the missing 63/250 per-lattice tables is honest",
          "K=63/250" not in
          next(L for L in src["data"]["lattices"] if L["name"] == "3x3x2")
          ["evaluations"])

    # 5. verdict recorded
    check("verdict FALSIFIED-for-purpose recorded",
          "FALSIFIED" in own["data"]["verdict"])

    # 6. wave-17 Lambda=3,4 extension: re-derive eps* from certified inputs and
    #    confirm the box shells fail it at every evaluation point.
    lam34 = json.load(open(os.path.join(ROOT, "results", "shell_stress_lambda34.json")))
    Jlo = {int(k): Fraction(v) for k, v in
           src["data"]["shell_weights"]["J_lo"].items()}
    dts = src["data"]["delta_star_thresholds"]
    Kstr = {"incumbent": "63/250", "K=1/4": "1/4",
            "K=49/200": "49/200", "K=6/25": "6/25"}
    for rec in lam34["eps_star_rows"]:
        K = Fraction(rec["K"]); Lam = rec["Lambda"]
        lo_d, hi_d = (Fraction(x) for x in dts[f"K={rec['K']}"][:2]) \
            if f"K={rec['K']}" in dts else (Fraction(rec["K"]), Fraction(rec["K"]))
        e_lo = lo_d * (2 * K) / Jlo[Lam]
        check(f"re-derived eps*(K={rec['K']},L={Lam}) matches artifact",
              e_lo == Fraction(rec["eps_star_interval"][0]))
    box_rows = lam34["box_shell_rows"]
    check("all 48 (box, eval, Lambda) combinations recorded", len(box_rows) == 48)
    check("all 48 premise statuses FALSIFIED, 0 ALIVE",
          all(r["premise_status"].startswith("FALSIFIED") for r in box_rows))
    check("tightest combination is 3x3x2 incumbent (shared across shells)",
          max(Fraction(r["weakest_mode_defect_upper_end"]) for r in box_rows)
          == Fraction("+0.005196")
          and all(r["evaluation"] == "incumbent" and r["lattice"] == "3x3x2"
                  for r in box_rows
                  if Fraction(r["weakest_mode_defect_upper_end"])
                  == Fraction("+0.005196")))
    check("concentration factor Lambda=3 ~4x, Lambda=4 ~11.5x",
          abs(2 * Fraction(1, 4) / Jlo[3] - 4.02) < 0.02
          and abs(2 * Fraction(1, 4) / Jlo[4] - 11.4) < 0.05)
    check("lambda34 verdict closes the PM-S family",
          "closed" in lam34["verdict"])

    n_ok = sum(1 for _, c in ok if c)
    print(f"\n{n_ok}/{len(ok)} checks passed")
    assert n_ok == len(ok), [n for n, c in ok if not c]


if __name__ == "__main__":
    main()
    print("PASS")
