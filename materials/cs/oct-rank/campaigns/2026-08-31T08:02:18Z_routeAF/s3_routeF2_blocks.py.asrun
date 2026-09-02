"""s3_routeF2_blocks.py — Route F completion: the exact 4x4 block algebra
of the (L_1, L_i, L_j) family in the H' (+) H'^perp splitting (H' =
span{e0..e3}, complement = span{e4..e7} — a CD doubling, exactly).

Established exactly by s3_routeF_subalgebra.py:
  span{e0..e3} is a quaternion subalgebra, and for x in H' (e0,e1,e2),
  L_x maps H' -> H' and H'^perp -> H'^perp (columns 0..3 in H', 4..7 in
  the complement) — block DIAGONAL; for x in the complement, L_x maps
  H' -> comp and comp -> H' — block OFF-diagonal.

Therefore the three slices decompose in the 2 x 2 block form over the
splitting R^8 = H' (+) C:
   L_1 = (I4,   0  )     L_i = (B_i, 0)      L_j = (B_j, 0)
         (0,    I4' )           (0, B_i')           (0, B_j')
with B_x the 4x4 quaternion left-multiplication on H' and B'_x its
copy on the complement (BUT with the opposite orientation coming from
the CD conjugation — need EXACT block extraction; the complement action
could be R rather than L, giving the transpose).

RANK CONSEQUENCE (the interesting question): with the family block
diagonal as (A1 (+) A1', B1 (+) B1', C1 (+) C1') the rank of the 3-slice
tensor is... NOT simply the sum of the block ranks (the terms in a CP
decomposition can MIX blocks). What IS immediate:
  lower: rank >= max over the two blocks of the corresponding bound:
    >= max( rank((I4, B_i, B_j)), rank((I4', B_i', B_j')) ) = max(7, 7)
    = 7 if the blocks are ± the quaternion tau (rank exactly 7).
  upper: rank <= 7 + 7 = 14 IF the two blocks can be decomposed
    independently (always: concatenate witnesses) — so rank <= 14 and
    possibly < 14 via cross terms.

Compute EXACTLY here:
  1. B_i, B_j on H' in the (e0,e1,e2,e3) basis: the quaternion table.
  2. B'_i, B'_j on the complement in (e4..e7): verify whether B'_x =
     B_x, -B_x, B_x^T or other — EXACT entrywise.
  3. rank((I4, B_i, B_j)) upper=7: the frozen tau witness gives 7; here:
     only re-verify the block IS the tau tensor (then rank = 7 known).
  4. Hence block-diagonal 3-family (I_8, L_i, L_j) has rank EXACTLY
     <= 14 and >= 7. For the S3 question the WORST CASE (independent
     triple) needs more; this specific triple — u=e0-ish — is fixed.
"""
import json
import sys
from flint import fmpq, fmpq_mat

sys.path.insert(0, ".")  # noqa: E402
from s3_routeAB import (L_of, mat, rank, conj, E8, cd_mul,  # noqa: E402,F401
                        check_inverse_identity)
sys.path.insert(0, "../..")  # noqa: E402


def qmul(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])


def main():
    out = {}
    assert all(check_inverse_identity(E8[i]) for i in range(8))
    Hbasis = [E8[i] for i in range(4)]
    Cbasis = [E8[i] for i in range(4, 8)]

    def block_of(x, basis):
        cols = []
        for h in basis:
            prod = cd_mul(x, h)
            vec = list(prod)
            coords = [sum(a * b for a, b in zip(vec, bh)) for bh in basis]
            cols.append(coords)
        return [[cols[b][c] for b in range(4)] for c in range(4)]

    B1 = block_of(E8[1], Hbasis)
    B2 = block_of(E8[2], Hbasis)
    Bc1 = block_of(E8[1], Cbasis)
    Bc2 = block_of(E8[2], Cbasis)
    B1m, B2m, Bc1m, Bc2m = (fmpq_mat([[fmpq(v) for v in row] for row in M])
                            for M in (B1, B2, Bc1, Bc2))
    # quaternion table check on H':
    # cleaner: direct quaternion left-mult tables:
    def Lq(qv):
        E4 = [tuple(1 if k == m else 0 for k in range(4)) for m in range(4)]
        return [[qmul(tuple(qv), E4[b])[c] for b in range(4)]
                for c in range(4)]
    Lqi = Lq((0, 1, 0, 0))
    Lqj = Lq((0, 0, 1, 0))
    out["Bi_equals_quaternion_Li"] = (B1 == Lqi)
    out["Bj_equals_quaternion_Lj"] = (B2 == Lqj)
    out["Bi_complement_equals_Li"] = (Bc1 == Lqi)
    out["Bi_complement_equals_negLi"] = (Bc1 == [[-v for v in row]
                                                 for row in Lqi])
    out["Bi_complement_equals_Li_transpose"] = (
        Bc1 == [[Lqi[c][b] for c in range(4)] for b in range(4)])
    out["Bi_complement_equals_negLi_transpose"] = (
        Bc1 == [[-Lqi[c][b] for c in range(4)] for b in range(4)])
    out["Bi_comp_eq_Bj_comp_type"] = (Bc1 == Bc2, Bc1, Bc2)

    # exact block ranks:
    out["rank_Bi"] = rank(B1m)
    out["rank_Bj"] = rank(B2m)
    out["rank_Bci"] = rank(Bc1m)
    out["comm_block_rank"] = rank(B1m * B2m - B2m * B1m)
    out["comm_block_comp_rank"] = rank(Bc1m * Bc2m - Bc2m * Bc1m)
    print(json.dumps(out, indent=1))
    sys.exit(0)


if __name__ == "__main__":
    main()
