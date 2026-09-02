"""Plants stage — M12ALPHA-DIVONLY (standalone, writes plants_result.json).

CF-1/CF-2/CF-3 exactly as amended (AMENDMENTS 2 and 3, committed before
any plant run).
"""
import json
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mceliece/src")

from instance import Instance  # noqa: E402
from gfield import GF  # noqa: E402
from fastfield import EField  # noqa: E402
import census as _c  # noqa: E402

OUT = ("/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/"
       "2026-09-01T03-28-10Z_M12ALPHA_DIVONLY")
rec = {"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

# ---------------- CF-1: duplicate support (registered dup a12 := a11)
cf1 = {"plant": "duplicate_support", "base": "(11,2048,48,6211)"}
m, n, t, seed = 11, 2048, 48, 6211
gf = GF(m)
import random as _r
rng = _r.Random(seed)
tries = 0
while True:
    tries += 1
    G = [rng.randrange(gf.q) for _ in range(t)] + [1]
    if gf.pis_irreducible(G):
        break
idx = list(range(gf.q))
rng.shuffle(idx)
support = idx[:n]
cf1["pre_dup_positions_10_11"] = [support[10], support[11]]
support[11] = support[10]          # registered duplication: a_12 := a_11
cf1["g_tries_rng_consumed"] = tries
Pi = [1]
for a in support:
    Pi = gf.ptrim(gf.pmul(Pi, [a, 1]))
PiD = gf.pderiv(Pi)
PiDa = [gf.peval(PiD, a) for a in support]
dup_zero = [i for i, v in enumerate(PiDa) if v == 0]
cf1["PiD_zero_at"] = dup_zero[:8]
cf1["PiD_zero_count"] = len(dup_zero)
cf1["rejected"] = bool(dup_zero)
cf1["rejected_by"] = ("Pi'(a)=0 at duplicated position(s); the build's "
                      "distinct-support assert (instance.py, PiDa check) "
                      "fires" if dup_zero else "NOT REJECTED")
rec["CF1_duplicate_support"] = cf1
print("CF1:", cf1, flush=True)

# ---------------- shared m=11 anchor instance (for CF-2, CF-3)
inst11 = Instance(11, 2048, 48, 6211)
gf11, ef11 = inst11.gf, inst11.ef
D11, k11, n11 = inst11.D, inst11.k, inst11.n
Dp1 = D11 + 1
Cmat = np.zeros((k11, Dp1), dtype=np.uint16)
for j, f in enumerate(inst11.F):
    Cmat[j, :len(f)] = f
Pi11 = np.asarray(inst11.Pi, dtype=np.uint16)
PiD11 = np.asarray(inst11.PiD, dtype=np.uint16)
G11 = np.asarray(inst11.G, dtype=np.uint16)
nz11 = np.nonzero(Cmat.any(axis=0))[0]


def alpha_vals_matrix(C, i):
    """alpha LHS/RHS value vectors at support point i for coefficient
    matrix C (k x Dp1, uint16). Derivative: slot s <- c_{s+1}, kept iff
    s even (s+1 odd)."""
    a = inst11.support[i]
    w0 = _c.lucas_w(ef11, Dp1 - 1, 0, a)
    Fa = np.bitwise_xor.reduce(ef11.MUL[C[:, nz11], w0[nz11][None, :]], axis=1)
    Dm2 = np.zeros_like(C)
    Dm2[:, 0: Dp1 - 1] = C[:, 1:Dp1]
    Dm2[:, 1::2] = 0
    Fpa = np.bitwise_xor.reduce(ef11.MUL[Dm2[:, nz11], w0[nz11][None, :]], axis=1)
    return _finish(Fa, Fpa, a, ef11, Dp1 - 1, k11)


def _finish(Fa, Fpa, a, ef, Ddeg, kk):
    wp = _c.lucas_w(ef, len(Pi11) - 1, 0, a)
    pia = int(np.bitwise_xor.reduce(ef.MUL[Pi11, wp]))
    w1 = _c.lucas_w(ef, len(PiD11) - 1, 0, a)
    pida = int(np.bitwise_xor.reduce(ef.MUL[PiD11, w1]))
    Ga = int(np.bitwise_xor.reduce(ef.MUL[G11, _c.lucas_w(ef, len(G11) - 1, 0, a)]))
    lhs = ef.MUL[pia, Fpa] ^ ef.MUL[pida, Fa]
    rhs = ef.MUL[ef.MUL[Ga, Ga], ef.MUL[Fa, Fa]]
    return lhs, rhs


def alpha_vals_lists(flist, i):
    """alpha at support point i from coefficient LISTS (any length;
    exact, no truncation). flist[j] = f_j coefficient list."""
    a = inst11.support[i]
    Fa = [gf11.peval(f, a) for f in flist]
    dflist = [gf11.pderiv(f) for f in flist]
    Fpa = [gf11.peval(df, a) for df in dflist]
    return _finish(np.array(Fa, dtype=np.uint16),
                   np.array(Fpa, dtype=np.uint16), a, ef11, D11, k11)


# sanity 1: uncorrupted matrix path, alpha OK at point 0
lhs0, rhs0 = alpha_vals_matrix(Cmat, 0)
rec["sanity_uncorrupted_matrix_alpha_point0"] = bool(np.array_equal(lhs0, rhs0))
# sanity 2: list path agrees with matrix path on rows at point 0
lhs0l, rhs0l = alpha_vals_lists(inst11.F, 0)
rec["sanity_list_equals_matrix_point0"] = bool(
    np.array_equal(lhs0, lhs0l) and np.array_equal(rhs0, rhs0l))
# sanity 3: list-path derivatives via pderiv equal shift logic at point 5
lhs5l, rhs5l = alpha_vals_lists(inst11.F, 5)
lhs5m, rhs5m = alpha_vals_matrix(Cmat, 5)
rec["sanity_point5_agreement"] = bool(
    np.array_equal(lhs5l, lhs5m) and np.array_equal(rhs5l, rhs5m))
print("sanity:", {k: v for k, v in rec.items() if k.startswith("sanity")}, flush=True)

# ---------------- CF-2: f0 -> f0 + Pi (degree plant) — LIST path (deg 2048
# exceeds the Dp1-truncated matrix on purpose; that is the plant)
cf2 = {"plant": "f0_plus_Pi", "base": "(11,2048,48,6211)"}
f0p = gf11.ptrim(gf11.padd(inst11.F[0], inst11.Pi))
cf2["deg_f0"] = gf11.pdeg(inst11.F[0])
cf2["deg_f0_plus_Pi"] = gf11.pdeg(f0p)
cf2["D"] = D11
cf2["deg_ok_expected_false"] = cf2["deg_f0_plus_Pi"] > D11
Fcorrupt = [f0p if j == 0 else inst11.F[j] for j in range(k11)]
vfail = 0
vfail_first = []
pts = list(range(0, n11, 97))
for i in pts:
    lhs, rhs = alpha_vals_lists(Fcorrupt, i)
    if not np.array_equal(lhs, rhs):
        vfail += 1
        vfail_first.append(i)
cf2["value_check_points"] = len(pts)
cf2["value_check_failures"] = vfail
cf2["value_check_first_fail_at"] = vfail_first[:3]
cf2["matches_amendment2_invisibility_prediction"] = (vfail == 0)
cf2["ABORT_semantics"] = ("degree condition max_j deg f_j <= D FAILED "
                          f"({cf2['deg_f0_plus_Pi']} > {D11}): per ADDENDUM "
                          "2(i) NO delta verdict is produced")
cf2["verdict"] = "REJECTED via degree gate (ABORT, no delta verdict)"
rec["CF2_f0_plus_Pi"] = cf2
print("CF2:", cf2, flush=True)

# ---------------- CF-3: row0 -> 3*f0 (alpha plant, AMENDMENT 3) — matrix path
cf3 = {"plant": "row0_scalar_x_3", "base": "(11,2048,48,6211)"}
c3 = 3
cf3["c"] = c3
cf3["c_squared"] = int(gf11.mul(c3, c3))
cf3["c_sq_ne_c"] = cf3["c_squared"] != c3
C3 = Cmat.copy()
C3[0] = ef11.MUL[c3, Cmat[0]]
fails = []
checked = 0
for i in range(n11):
    lhs, rhs = alpha_vals_matrix(C3, i)
    badj = np.nonzero(lhs != rhs)[0]
    checked += 1
    if badj.size:
        fails.append({"i": i, "rows_failing_first5": badj[:5].tolist()})
        if len(fails) >= 3:
            break
cf3["points_checked"] = checked
cf3["first_failing_points"] = fails
cf3["rejected_by_alpha_value_check"] = bool(fails)
cf3["verdict"] = ("REJECTED via alpha value check" if fails
                  else "PASSED (would be a finding)")
rec["CF3_row_scalar"] = cf3
print("CF3:", {k: v for k, v in cf3.items() if k != "first_failing_points"}, flush=True)

rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
json.dump(rec, open(f"{OUT}/plants_result.json", "w"), indent=1, default=str)
print("plants_result.json written", flush=True)
