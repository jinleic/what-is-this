"""gate_c_lower_rank4.py — exact-rational lower-bound scaffolding for the
quaternion vector multiplication tensor tau_H (4x4x4, |entries| <= 1).

Machine-checked facts (all exact, python-flint fmpq):

  F1  Each 4x4 slice L_1, L_i of tau_H has matrix rank exactly 4, and the
      slices are pairwise non-proportional over Q.  Hence no rank argument
      via proportional-slice reduction applies, and any 3-term claim would
      need a genuine CP lower bound (NOT established here — see
      pre_statement.md C1(a): the lower bound 7 is the paper's peel+pencil
      argument, HUMAN-AUDITED territory, not a flattening computation).

  F2  For every homogeneous sign pattern s[p,b,c] = sigma[p] * t[c]
      (sigma in {+-1}^3, t in {+-1}^4 — the only sign flips that are
      automorphism-free and slice-homogeneous), the 12x4 flattening of the
      sign-flipped tensor still has rank exactly 4: there is NO sign-flip
      of tau_H whose flattening rank drops to 3 or below.  Flattening rank
      <= CP rank, so no sign-conjugate of tau_H has CP rank < 4 by this
      route either.

SCOPE HONESTY: flattening rank bounds are necessary conditions only.
None of F1/F2 is a CP lower bound.  This file existso that the claim
"R_R(tau) = 7" is re-verified only through (a) Krawczyk upper (7 <=) and
(b) the paper's peel+pencil lower argument whose Lean artifact replay
constitutes the machine check of (>= 7)  [CITED-DEPENDENCY].
"""
import itertools
import json

from flint import fmpq, fmpq_mat


def qmul(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])


E4 = [[1 if k == m else 0 for k in range(4)] for m in range(4)]
w = (E4[0], E4[1], E4[2])   # slices: left multiplication by 1, i, j
tau = [[[qmul(w[p], E4[b])[c] for c in range(4)] for b in range(4)]
       for p in range(3)]


def flat_rank(rows):
    m = fmpq_mat(len(rows), len(rows[0]),
                 [fmpq(int(v)) for row in rows for v in row])
    return int(m.rref()[1])


ranks = {
    "slice_1": flat_rank(tau[1]),
    "slice_i": flat_rank(tau[2]),
}

T1 = fmpq_mat(4, 4, [fmpq(int(v)) for row in tau[1] for v in row])
T2 = fmpq_mat(4, 4, [fmpq(int(v)) for row in tau[2] for v in row])
P = T2 * T1.inv()
ranks["slices_pairwise_nonprop"] = not all(
    P[i, j] == P[0, 0] for i in range(4) for j in range(4))

# F2: all 2^3 * 2^4 = 128 homogeneous sign patterns
low = []
total = 0
for t in itertools.product([1, -1], repeat=4):
    for bits in range(8):
        sigma = [1 if (bits >> p) & 1 == 0 else -1 for p in range(3)]
        s = [[[sigma[p] * t[c] for c in range(4)] for b in range(4)]
             for p in range(3)]
        flat = [[fmpq(int(tau[p][b][c])) * s[p][b][c] for c in range(4)]
                for p in range(3) for b in range(4)]
        rk = flat_rank(flat)
        total += 1
        if rk < 4:
            low.append({"sigma": sigma, "t": list(t), "rank": rk})

verdict = {
    "facts": {
        "F1": {"slice_ranks": ranks,
               "note": "flattening/slice ranks exactly 4; not a CP bound"},
        "F2": {"homogeneous_sign_patterns_total": total,
               "flattening_below_4": len(low),
               "examples": low[:8]},
    },
    "scope": "scaffolding; lower bound R_R(tau)>=7 is NOT claimed from these",
}
ok = ranks["slice_1"] == 4 and ranks["slice_i"] == 4 \
    and ranks["slices_pairwise_nonprop"] and len(low) == 0
verdict["all_exact_checks_pass"] = bool(ok)
print(json.dumps(verdict, indent=1))
