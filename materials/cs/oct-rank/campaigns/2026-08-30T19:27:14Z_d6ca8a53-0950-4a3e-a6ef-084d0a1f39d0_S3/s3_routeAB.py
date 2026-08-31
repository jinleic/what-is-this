"""s3_routeAB.py — Routes A and B: Strassen commutator + Blaeser/Lickteig
improvements, EXACTLY, on 3-slice L-families (u, v, w independent in R^8).

FIX of s3_routeA_strassen.py (rule 5 retraction, recorded in the VERDICT):
that script normalized A^{-1}B as N(u)^{-1} L_{ubar} L_v and then computed
rank([X, Y]) on the HALF-scaled pair directly; the derived quantity fed to
the Strassen bound must be the commutator of the NORMALIZED matrices
[X, Y] where X = A^{-1}B, Y = A^{-1}C — that part was right — but the
formula printed as 8 + rank/2 used the WRONG rank: with X, Y 8x8 and
[X, Y] skew of rank 8 everywhere on the swept families, the bound reads
8 + 4 = 12, NOT 16.  The bound that uses rank commutator divided by 2 is
Strassen's R >= n + rank([X,Y]) for TRIPLE products A,B,C; the /2 form is
the PAIR (pencil) case.  This script computes BOTH:
    s_pair := rank([X, Y])           (used in the /2 pencil variant)
    gives  n + s_pair/2  and  n + s_pair  candidates,
and reports each honestly. The commutator rank itself is the only exact
input; both formula variants are reported (no dressing).

Families swept (pre-registered in pre_statement.md S3 routes):
  F1: all (e_a, e_b, e_c), a<b<c — 56 triples.
  F3: (e_a+e_b, e_a, e_c) and (e_a, e_a+e_b, e_c) — 112 triples.
  F4: 4 fixed mixed-coefficient triples (list in code, deterministic).
Control at n=4: quaternions (1,i,j) and (i,j,k), plus the worst-case
across all basis triples for H.

All ranks exact (fmpq_mat.rref-based rank over Q).
"""
import json
import sys
from itertools import combinations
from flint import fmpq, fmpq_mat

Q0 = fmpq(0)


# ---------------- CD octonions, verified identical to upstream/gate-A table
def cd_conj(x):
    n = len(x)
    if n == 1:
        return x
    m = n // 2
    return cd_conj(x[:m]) + tuple(-v for v in x[m:])


def cd_mul(x, y):
    n = len(x)
    if n == 1:
        return (x[0] * y[0],)
    m = n // 2
    a, b = x[:m], x[m:]
    c, d = y[:m], y[m:]
    left = tuple(p - r for p, r in zip(cd_mul(a, c),
                                       cd_mul(cd_conj(d), b)))
    right = tuple(p + r for p, r in zip(cd_mul(d, a),
                                        cd_mul(b, cd_conj(c))))
    return left + right


E8 = [tuple(1 if k == i else 0 for k in range(8)) for i in range(8)]
T_tab = [[cd_mul(E8[i], E8[j]) for j in range(8)] for i in range(8)]


def L_of(u):
    n = 8
    return [[sum(u[p] * T_tab[p][b][c] for p in range(n)) for b in range(n)]
            for c in range(n)]


def mat(rows):
    return fmpq_mat([[fmpq(v) for v in row] for row in rows])


def rank(M):
    return M.rref()[1]


def conj(u):
    return (u[0],) + tuple(-x for x in u[1:])


def comm_rank(u, v, w):
    """rank([L_ubar L_v, L_ubar L_w]) computed exactly."""
    Lub = mat(L_of(conj(u)))
    X = Lub * mat(L_of(v))
    Y = Lub * mat(L_of(w))
    return rank(X * Y - Y * X)


def check_inverse_identity(u):
    n = 8
    Nu = sum(x * x for x in u)
    prod = mat(L_of(conj(u))) * mat(L_of(u))
    for i in range(n):
        for j in range(n):
            want = fmpq(Nu) if i == j else Q0
            if prod[i, j] != want:
                return False
    return True


def main():
    out = {}
    assert all(check_inverse_identity(E8[i]) for i in range(8)), \
        "L_ubar L_u = N(u) I failed — table defect"
    # validate convention: L_of(e_1) column 1 == e_1 products checked by
    # norm multiplicativity (already asserted above via identity check).

    # n=4 control (quaternions, assoc => [L_v, L_w] has rank = |v x w|-ish)
    def qmul(a, b):
        return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
                a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
                a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
                a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])
    E4 = [tuple(1 if k == i else 0 for k in range(4)) for i in range(4)]
    qtab = [[qmul(E4[i], E4[j]) for j in range(4)] for i in range(4)]

    def L4(u):
        return [[sum(u[p] * qtab[p][b][c] for p in range(4))
                 for b in range(4)] for c in range(4)]

    def cr4(u, v, w):
        ub = conj(u)
        A = mat(L4(ub)) * mat(L4(v))
        B = mat(L4(ub)) * mat(L4(w))
        return rank(A * B - B * A)

    h_ctrl = [cr4(E4[a], E4[b], E4[c])
              for a, b, c in combinations(range(4), 3)]
    out["control_H_commutator_ranks_all_basis_triples"] = h_ctrl
    tau_c = cr4(E4[0], E4[1], E4[2])
    out["control_tau_1ij"] = tau_c
    # On H: (L_i, L_j) normalized pencil = (I, L_{ij}) on imaginary part;
    # [L_{ubar}L_v, L_{ubar}L_w] = [L_{ubar v}, L_{ubar w}] = [L_k_assoc...]:
    # with u=1: [L_i, L_j] = L_i L_j - L_j L_i = L_k - L_-k = (matrix) 2 L_k
    # rank 4. Strassen twin forms: 4+4=8 and 4+2=6 (true tau rank is 7:
    # BOTH miss; boundaries: twin-6 < 7 < 8).

    # ---- n=8: F1 all basis triples
    f1 = [comm_rank(E8[a], E8[b], E8[c])
          for a, b, c in combinations(range(8), 3)]
    out["F1_comm_rank_min"] = min(f1)
    out["F1_comm_rank_max"] = max(f1)
    dist = {}
    for r in f1:
        dist[r] = dist.get(r, 0) + 1
    out["F1_distribution"] = dist

    # ---- n=8: F3 mixed
    vals = []
    for a, b, c in combinations(range(8), 3):
        u1 = tuple(E8[a][i] + E8[b][i] for i in range(8))
        vals.append(comm_rank(u1, E8[a], E8[c]))
        vals.append(comm_rank(E8[a], u1, E8[c]))
    out["F3_min"] = min(vals)
    out["F3_max"] = max(vals)

    # ---- n=8: F4 fixed
    fixed = [
        ((1, 1, 0, 0, 0, 0, 0, 0), (0, 1, 1, 0, 0, 0, 0, 0),
         (0, 0, 1, 1, 0, 0, 0, 0)),
        ((1, 1, 1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 1, 0, 0, 0),
         (0, 0, 1, 0, 0, 1, 0, 0)),
        ((2, 1, 0, 0, 0, 0, 0, 0), (1, -1, 3, 0, 0, 1, 0, 0),
         (0, 4, -2, 1, 1, 0, 1, 0)),
        ((1, 0, 0, 0, 1, 0, 0, 0), (0, 1, 0, 0, 0, 1, 0, 0),
         (0, 0, 1, 0, 0, 0, 1, 2)),
        ((1, 2, 3, 4, 5, 6, 7, 8), (8, 7, 6, 5, 4, 3, 2, 1),
         (1, -1, 2, -2, 3, -3, 4, -4)),
    ]
    vals4 = [comm_rank(u, v, w) for (u, v, w) in fixed]
    out["F4_values"] = vals4
    out["F4_min"] = min(vals4)

    lo = min(min(f1), min(vals), min(vals4))
    out["overall_comm_rank_min"] = lo
    out["bounds_reported"] = {
        "pencil_form_8_plus_half": 8 + lo // 2,
        "strassen_8_plus_rank": 8 + lo,
    }
    print(json.dumps(out, indent=1))
    sys.exit(0)


if __name__ == "__main__":
    main()
