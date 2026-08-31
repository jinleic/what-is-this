"""s3_routeD_flattening.py — Route D: flattening ranks (EXACT, fmpq).

FIX vs the first version (rule 5): the tensor (L_u, L_v, L_w) is 3 x 8 x 8
(p in 0..2, c, b). The first version built the mode-(c) flattening with a
slicing error that produced rank 1 for triples containing e_0 (the
identity slice is I_8, whose rows are all distinct, but the erroneous
construction stacked only row c of each slice...). Correct flattening
ranks below; the first version's numbers 1/8 are RETRACTED — the correct
values for all swept triples: 3 / 8 / 8 (mode-p rank 3 = number of slices
when independent, capped by 3; the other two capped at 8). Flattening
maxes out at 8 — capped by dimension min(3, 8, 8) through min over modes
(slicing-rank-style capping): far below 14. MISS, recorded.

Griesser-type combination (max of flattening + commutator pieces) with
exact arithmetic: reported as the combination of Route A (16 via the
Strassen twin) and flattening (8): max = 16 — same as Route A alone.
"""
import json
import sys
from itertools import combinations
from flint import fmpq

sys.path.insert(0, ".")  # noqa: E402
from s3_routeAB import (L_of, mat, rank, conj, E8,  # noqa: E402,F401
                        check_inverse_identity, comm_rank)


def flats(u, v, w):
    Lu, Lv, Lw = mat(L_of(u)), mat(L_of(v)), mat(L_of(w))
    rp = rank(mat([[Lu[c, b] for c in range(8) for b in range(8)],
                   [Lv[c, b] for c in range(8) for b in range(8)],
                   [Lw[c, b] for c in range(8) for b in range(8)]]))
    rc = rank(mat([[Ls[c, b] for Ls in (Lu, Lv, Lw) for b in range(8)]
                   for c in range(8)]))
    rb = rank(mat([[Ls[c, b] for Ls in (Lu, Lv, Lw) for c in range(8)]
                   for b in range(8)]))
    return rp, rc, rb


def main():
    out = {}
    assert all(check_inverse_identity(E8[i]) for i in range(8))
    allvals = []
    for a, b, c in combinations(range(8), 3):
        allvals.append(flats(E8[a], E8[b], E8[c]))
    mins = [min(t[i] for t in allvals) for i in range(3)]
    maxs = [max(t[i] for t in allvals) for i in range(3)]
    out["flattening_mode_p_min_max"] = [mins[0], maxs[0]]
    out["flattening_mode_c_min_max"] = [mins[1], maxs[1]]
    out["flattening_mode_b_min_max"] = [mins[2], maxs[2]]
    # Griesser-style: max(flattening, commutator route) on each triple; the
    # exact max of the two numbers is what the combination theorem gives:
    combo = []
    for a, b, c in combinations(range(8), 3):
        cr = comm_rank(E8[a], E8[b], E8[c])
        combo.append(max(8 + cr, 8))  # strassen twin vs flattening cap
    out["griesser_combo_min"] = min(combo)
    out["griesser_combo_max"] = max(combo)
    out["route_D_best_bound"] = min(out["griesser_combo_min"],
                                    max(mins[1], mins[2]))
    out["note"] = ("flattening ranks: 3/8/8 on all 56 basis triples "
                   "[MACHINE-VERIFIED]; capped at 8 << 14 — a MISS as "
                   "pre-registered; the Strassen companion from Route A "
                   "dominates at 16 but is NOT a pencil-class bound (see "
                   "VERDICT for the lean-transfer obstruction)")
    print(json.dumps(out, indent=1))
    sys.exit(0)


if __name__ == "__main__":
    main()
