"""Plants stage — T96 cell (standalone, writes plants_result.json).

MECHANICAL RE-INSTATIATION of M12ALPHA_DIVONLY/plants_stage.py with the
t=64 campaign's amendment history inherited verbatim (AMENDMENT 2:
CF-2 value-level invisibility prediction; AMENDMENT 3: CF-3 as
scalar-multiply — row-XOR proven F2-inert a priori, never
re-registered), PLUS the pre-statement section 4 P-4 add-on: a bounded
128-point t=96-native direction check using the SAME verifier ops —
the built (12,3488,96,16384) instance must be alpha-clean on the probe
points, and its row-0 x3 corruption must FAIL at least one of them.

All plants run through the SAME inherited verifier ops (gfield/fastfield
scalar path + the exact alpha value check). No outside ad-hoc checker.
"""
import json
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/20260901T100025Z_bd55fdfd_8fbba9be5ab3/code")

from instance import Instance  # noqa: E402
from gfield import GF  # noqa: E402
from fastfield import EField  # noqa: E402
import census as _c  # noqa: E402

RUN = "20260901T100025Z_bd55fdfd_8fbba9be5ab3"
OUT = f"/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/{RUN}"
rec = {"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

# ---------------- P-1: duplicate support (registered dup a12 := a11), base (11,2048,48,6211)
p1 = {"plant": "duplicate_support", "base": "(11,2048,48,6211)"}
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
p1["pre_dup_positions_10_11"] = [support[10], support[11]]
support[11] = support[10]          # registered duplication: a_12 := a_11
p1["g_tries_rng_consumed"] = tries
Pi = [1]
for a in support:
    Pi = gf.ptrim(gf.pmul(Pi, [a, 1]))
PiD = gf.pderiv(Pi)
PiDa = [gf.peval(PiD, a) for a in support]
dup_zero = [i for i, v in enumerate(PiDa) if v == 0]
p1["PiD_zero_at"] = dup_zero[:8]
p1["PiD_zero_count"] = len(dup_zero)
p1["rejected"] = bool(dup_zero)
p1["rejected_by"] = ("Pi'(a)=0 at duplicated position(s); the build's "
                     "distinct-support assert (instance.py, PiDa check) "
                     "fires" if dup_zero else "NOT REJECTED")
rec["P1_duplicate_support"] = p1
print("P1:", p1, flush=True)

# ---------------- shared m=11 anchor instance (for P-2, P-3)
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
    return _finish(Fa, Fpa, a, ef11, k11)


def _finish(Fa, Fpa, a, ef, kk):
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
                   np.array(Fpa, dtype=np.uint16), a, ef11, k11)


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
print("sanity:", {k: v for k, v in rec.items() if k.startswith("sanity")},
      flush=True)

# ---------------- P-2: f0 -> f0 + Pi (degree plant) — LIST path
p2 = {"plant": "f0_plus_Pi", "base": "(11,2048,48,6211)"}
f0p = gf11.ptrim(gf11.padd(inst11.F[0], inst11.Pi))
p2["deg_f0"] = gf11.pdeg(inst11.F[0])
p2["deg_f0_plus_Pi"] = gf11.pdeg(f0p)
p2["D"] = D11
p2["deg_ok_expected_false"] = p2["deg_f0_plus_Pi"] > D11
Fcorrupt = [f0p if j == 0 else inst11.F[j] for j in range(k11)]
vfail = 0
vfail_first = []
pts = list(range(0, n11, 97))
for i in pts:
    lhs, rhs = alpha_vals_lists(Fcorrupt, i)
    if not np.array_equal(lhs, rhs):
        vfail += 1
        vfail_first.append(i)
p2["value_check_points"] = len(pts)
p2["value_check_failures"] = vfail
p2["value_check_first_fail_at"] = vfail_first[:3]
p2["matches_amendment2_invisibility_prediction"] = (vfail == 0)
p2["ABORT_semantics"] = ("degree condition max_j deg f_j <= D FAILED "
                         f"({p2['deg_f0_plus_Pi']} > {D11}): per ADDENDUM "
                         "2(i) NO delta verdict is produced")
p2["verdict"] = "REJECTED via degree gate (ABORT, no delta verdict)"
rec["P2_f0_plus_Pi"] = p2
print("P2:", p2, flush=True)

# ---------------- P-3: row0 -> 3*f0 (alpha plant, AMENDMENT 3 form) — matrix path
p3 = {"plant": "row0_scalar_x_3", "base": "(11,2048,48,6211)"}
c3 = 3
p3["c"] = c3
p3["c_squared"] = int(gf11.mul(c3, c3))
p3["c_sq_ne_c"] = p3["c_squared"] != c3
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
p3["points_checked"] = checked
p3["first_failing_points"] = fails
p3["rejected_by_alpha_value_check"] = bool(fails)
p3["verdict"] = ("REJECTED via alpha value check" if fails
                 else "PASSED (would be a finding)")
rec["P3_row_scalar"] = p3
print("P3:", {k: v for k, v in p3.items() if k != "first_failing_points"},
      flush=True)

# ---------------- P-4: t=96-native direction check (bounded probe grid)
# SAME verifier ops at the registered cell: the UNCORRUPTED instance must
# be alpha-clean on the first 128 support points; its row-0 x3 corruption
# must fail >= 1 of them. Instrument-direction check at m=12 scale.
p4 = {"plant": "t96_row0_scalar_x3_probe128",
      "base": "(12,3488,96,16384)"}
inst96 = Instance(12, 3488, 96, 16384)
gf96, ef96 = inst96.gf, inst96.ef
D96, k96, n96 = inst96.D, inst96.k, inst96.n
Dp196 = D96 + 1
Cm96 = np.zeros((k96, Dp196), dtype=np.uint16)
for j, f in enumerate(inst96.F):
    Cm96[j, :len(f)] = f
Pia96 = np.asarray(inst96.Pi, dtype=np.uint16)
PiDa96 = np.asarray(inst96.PiD, dtype=np.uint16)
Ga96 = np.asarray(inst96.G, dtype=np.uint16)
nz96 = np.nonzero(Cm96.any(axis=0))[0]


def a96(C, i):
    a = inst96.support[i]
    w0 = _c.lucas_w(ef96, Dp196 - 1, 0, a)
    Fa = np.bitwise_xor.reduce(ef96.MUL[C[:, nz96], w0[nz96][None, :]], axis=1)
    Dm = np.zeros_like(C)
    Dm[:, 0: Dp196 - 1] = C[:, 1:Dp196]
    Dm[:, 1::2] = 0
    Fpa = np.bitwise_xor.reduce(ef96.MUL[Dm[:, nz96], w0[nz96][None, :]], axis=1)
    wp = _c.lucas_w(ef96, len(Pia96) - 1, 0, a)
    pia = int(np.bitwise_xor.reduce(ef96.MUL[Pia96, wp]))
    w1 = _c.lucas_w(ef96, len(PiDa96) - 1, 0, a)
    pida = int(np.bitwise_xor.reduce(ef96.MUL[PiDa96, w1]))
    Ga = int(np.bitwise_xor.reduce(ef96.MUL[Ga96, _c.lucas_w(ef96, len(Ga96) - 1, 0, a)]))
    lhs = ef96.MUL[pia, Fpa] ^ ef96.MUL[pida, Fa]
    rhs = ef96.MUL[ef96.MUL[Ga, Ga], ef96.MUL[Fa, Fa]]
    return lhs, rhs


PROBE = 128
clean_fail = []
cor_fail = []
CmC = Cm96.copy()
CmC[0] = ef96.MUL[3, Cm96[0]]  # same x3 corruption mechanism as P-3
for i in range(PROBE):
    l0, r0 = a96(Cm96, i)
    if not np.array_equal(l0, r0):
        clean_fail.append(i)
    l1, r1 = a96(CmC, i)
    if not np.array_equal(l1, r1):
        cor_fail.append(i)
p4["probe_points"] = PROBE
p4["uncorrupted_clean_on_all_probes"] = (len(clean_fail) == 0)
p4["uncorrupted_failures"] = clean_fail[:5]
p4["corrupted_failing_points"] = cor_fail[:5]
p4["corrupted_fail_count"] = len(cor_fail)
p4["c"] = 3
p4["c_squared"] = int(gf96.mul(3, 3))
p4["rejected"] = bool(cor_fail) and (len(clean_fail) == 0)
p4["verdict"] = ("REJECTED via alpha value check at t=96 scale "
                 "(same ops, same verifier)" if cor_fail
                 and not clean_fail else
                 "PASSED (would be a finding)")
rec["P4_t96_probe"] = p4
print("P4:", {k: v for k, v in p4.items()}, flush=True)

rec["all_plants_rejected"] = bool(
    p1["rejected"] and p2["deg_ok_expected_false"] and p2["matches_amendment2_invisibility_prediction"]
    and p3["rejected_by_alpha_value_check"] and p4["rejected"])
rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
json.dump(rec, open(f"{OUT}/plants_result.json", "w"), indent=1, default=str)
print("ALL PLANTS REJECTED =", rec["all_plants_rejected"], flush=True)
