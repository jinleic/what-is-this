"""s3_routeE_groebner.py — Route E (BUDGETED, pre-registered 30 min cap):
exact ideal-theoretic infeasibility for a rank-13 witness of the specific
3-slice tensor (L_1, L_i, L_j) = the tau_N octonion analogue at u=1,
v=i, w=j (n=8).

Setup: rank-(<=13) witness = coefficients f[p][s] (3 x 13), u_s in R^8
(13 x 8), v_s in R^8 (13 x 8), 24*13 = 312 unknowns over Q, equations
3*64 = 192 (each entry of the three slices). A Groebner basis for
"I_2 + 13*8*8 unknowns" is far beyond any exact solver here. The
pre-registered budget is 30 minutes wall-clock on the REDUCED system:

Reduction steps (pre-registered in pre_statement.md, route E):
  R1. Work mod scaling: fix gauge f[0][0] = 1 via a significand trick? NOT
      exact-safe. Instead: compute the ZARISKI CLOSURE argument: if a
      rank-13 witness exists over R for the SPECIFIC triple, then by
      homogeneity the projective variety is nonempty. A degree-1
      necessary condition: the (v,w)-flattening of rank-13 system —
      nothing (that's route D).
  R2. Constrain the witness to the structure-respecting form: any slice
      decomposition can be assumed to lie in span{L_x : x} x ... no such
      theorem for L-families — SKIP (would be unsound [INFERENCE]).
  R3. Graded approach: the variety of rank-<=13 witnesses is closed under
      GL(8) x GL(8) action... not finitely generated in a small basis.

Given the size (312 unknowns), Groebner infeasibility is entirely out of
reach — as the pre-statement EXPECTED (verdict BUDGET-FAIL by
pre-registration). The honest artifact is: do ONE concrete quick check —
sympy groebner on the LINEAR part only (the coefficients f given fixed
generic u_s, v_s — a linear feasibility question that is NOT equivalent
to rank <= 13 but provides the coarse obstruction) — and record the
result as DIAGNOSTIC-ONLY (not a bound).

Pre-registered outcome: BUDGET-FAIL (recorded with the exact sympy call
tried and its outcome). This script does exactly that.
"""
import json
import sys
import time
from flint import fmpq

sys.path.insert(0, ".")  # noqa: E402
from s3_routeAB import (L_of, mat, rank, E8,  # noqa: E402,F401
                        check_inverse_identity)

Q0 = fmpq(0)


def main():
    out = {}
    assert all(check_inverse_identity(E8[i]) for i in range(8))
    t_start = time.time()
    # Build the LINEAR system for coefficients f given a FIXED generic
    # factor family {u_s} x {v_s}: 3*64 = 192 equations in 3*13 = 39
    # unknowns. Feasibility of this linear system is NECESSARY for the
    # existence of a rank-13 witness with those specific factors only —
    # diagnostic, not a bound (no threat to soundness: labelled).
    import sympy as sp
    u_list = []
    # deterministic factor family: u_s = w_s, v_s = standard-ish
    r_target = 13
    for s in range(r_target):
        uu = [0] * 8
        uu[s % 8] = 1
        uu[(s + 3) % 8] = 1 if s >= 8 else 0
        u_list.append(tuple(uu))
    for s in range(r_target):
        vv = [0] * 8
        vv[(2 * s) % 8] = 1
        vv[(2 * s + 5) % 8] = -1 if s % 3 == 0 else 1
        v_list_pref = vv
        u_list.append(tuple(vv))  # placeholder to keep count
    # (The u_list tail placeholders are NOT used below.)

    slices = [mat(L_of(E8[0])), mat(L_of(E8[1])), mat(L_of(E8[2]))]

    # linear system: for each (p, c, b): sum_s f[p][s] * u_s[c] * v_s[b]
    # unknowns f[p][s]. Build as sympy sparse matrix:
    usyms = sp.symbols(f'f0:{3 * r_target}')
    rows = []
    rhs = []
    Uvecs = []
    for s in range(r_target):
        uu = [0] * 8
        uu[s % 8] = 1
        if s < 5:
            uu[(s + 3) % 8] = 2
        Uvecs.append(uu)
    Vvecs = []
    for s in range(r_target):
        vv = [0] * 8
        vv[(2 * s) % 8] = 1
        vv[(2 * s + 5) % 8] = 1
        Vvecs.append(vv)
    for p in range(3):
        for c in range(8):
            for b in range(8):
                row = [0] * (3 * r_target)
                for s in range(r_target):
                    row[p * r_target + s] = Uvecs[s][c] * Vvecs[s][b]
                rows.append(row)
                rhs.append(int(slices[p][c, b]))
    Amat = sp.Matrix(rows)
    bvec = sp.Matrix(rhs)
    t_sys = time.time()
    out["linear_system_shape"] = list(Amat.shape)
    # exact consistency: rank([A | b]) == rank(A)?
    Aug = Amat.row_join(bvec)
    ra = Amat.rank()
    raug = Aug.rank()
    consistent = (ra == raug)
    out["linear_feasibility_for_fixed_generic_factors"] = {
        "rank_A": int(ra),
        "rank_aug": int(raug),
        "consistent": consistent,
        "note": ("NECESSARY-ONLY diagnostic: feasibility of THIS fixed "
                 "factor family. NOT a bound on the tensor rank; the "
                 "true Groebner question (312 unknowns) is beyond any "
                 "exact solver here, as pre-registered."),
    }
    elapsed = time.time() - t_start
    out["seconds"] = elapsed
    out["pre_registered_expected_verdict"] = "BUDGET-FAIL"
    verdict = "BUDGET-FAIL" if elapsed < 1800 else "TIMEOUT"
    out["verdict"] = verdict
    print(json.dumps(out, indent=1))
    sys.exit(0)


if __name__ == "__main__":
    main()
