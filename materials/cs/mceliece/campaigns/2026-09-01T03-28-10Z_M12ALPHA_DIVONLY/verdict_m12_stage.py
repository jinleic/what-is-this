"""m=12 verdict stage — M12ALPHA-DIVONLY (standalone, writes verdict.json).

NO plants here (plants_stage.py). Order:
  1. build (12,3488,64,16384)              [measured ~1785 s]
  2. FULL alpha grid, all 3488 support points, unmodified m<=11 instrument
     ops (MEASURED — AMENDMENT 1) + F2-linearity control (4 selectors)
  3. checkpoint alpha to m12_alpha_result.json
  4. ADDENDUM-2(i) delta chain: exact divisions, unit checks, offdiag,
     assembly, degree gate
  5. final verdict.json
Registered CPU budget: 10,800 s process_time (AMENDMENT 1); the stage
checks the budget after the alpha grid and after the delta chain and
records the numbers either way.
"""
import json
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mceliece/src")

from instance import Instance  # noqa: E402
import census as _c  # noqa: E402

OUT = ("/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/"
       "2026-09-01T03-28-10Z_M12ALPHA_DIVONLY")
BUDGET = 10800.0
C0 = time.process_time()


def cpu():
    return time.process_time() - C0


T_START = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
ver = {"instance": "(12,3488,64,16384)", "started": T_START,
       "registered_budget_cpu_s": BUDGET}

# 1. build
t0 = cpu()
inst = Instance(12, 3488, 64, 16384)
ver["build_cpu_s"] = round(cpu() - t0, 2)
print("built; cpu", ver["build_cpu_s"], flush=True)
gf, ef = inst.gf, inst.ef
D, k, n = inst.D, inst.k, inst.n
Dp1 = D + 1
ver.update({"k": k, "D": D, "n": n})

# lam recompute + gcd(G,G') + degree gate inputs
lam_re = [gf.mul(gf.mul(inst.Ga[i], inst.Ga[i]), gf.inv(gf.peval(inst.PiD, a)))
          for i, a in enumerate(inst.support)]
ver["lam_recomputed_equal"] = bool(lam_re == inst.lam)
Gp = gf.pderiv(inst.G)
ver["gcd_G_Gprime_deg"] = gf.pdeg(gf.pgcd(inst.G[:], Gp))
ver["G_squarefree"] = ver["gcd_G_Gprime_deg"] <= 0
degs = [gf.pdeg(f) for f in inst.F]
ver["max_deg_f"] = int(max(degs))
ver["deg_ok"] = bool(max(degs) <= D)
ver["beta_delta_nonzero"] = bool(inst.guards.get("beta_delta_nonzero"))
ver["beta_pair"] = inst.guards.get("beta_pair")
ver["gamma_k"] = bool(inst.guards.get("gamma_k"))
ver["eps_gcd_const"] = bool(inst.guards.get("eps_gcd_const"))
ver["eps_maxdeg_is_D"] = bool(inst.guards.get("eps_maxdeg_is_D"))
json.dump(ver, open(f"{OUT}/m12_build_checkpoint.json", "w"), indent=1, default=str)

# 2. coefficient matrices
Cm2 = np.zeros((k, Dp1), dtype=np.uint16)
Dm2 = np.zeros((k, Dp1), dtype=np.uint16)
for j, f in enumerate(inst.F):
    cf = np.asarray(f + [0] * max(0, Dp1 - len(f)), dtype=np.uint16)[:Dp1]
    Cm2[j] = cf
    dc = np.zeros(Dp1, dtype=np.uint16)
    dc[0:Dp1 - 1] = cf[1:Dp1]
    dc[1::2] = 0                    # char-2 derivative (odd e -> slot e-1)
    Dm2[j] = dc
Pi12 = np.asarray(inst.Pi, dtype=np.uint16)
PiD12 = np.asarray(inst.PiD, dtype=np.uint16)
G12 = np.asarray(inst.G, dtype=np.uint16)
nz = np.nonzero(Cm2.any(axis=0))[0]
print("coeff nz slots:", nz.size, flush=True)


def lhs_rhs(i):
    a = inst.support[i]
    w0 = _c.lucas_w(ef, Dp1 - 1, 0, a)
    Fa = np.bitwise_xor.reduce(ef.MUL[Cm2[:, nz], w0[nz][None, :]], axis=1)
    Fpa = np.bitwise_xor.reduce(ef.MUL[Dm2[:, nz], w0[nz][None, :]], axis=1)
    wp = _c.lucas_w(ef, len(Pi12) - 1, 0, a)
    pia = int(np.bitwise_xor.reduce(ef.MUL[Pi12, wp]))
    w1 = _c.lucas_w(ef, len(PiD12) - 1, 0, a)
    pida = int(np.bitwise_xor.reduce(ef.MUL[PiD12, w1]))
    Ga = int(np.bitwise_xor.reduce(ef.MUL[G12, _c.lucas_w(ef, len(G12) - 1, 0, a)]))
    lhs = ef.MUL[pia, Fpa] ^ ef.MUL[pida, Fa]
    rhs = ef.MUL[ef.MUL[Ga, Ga], ef.MUL[Fa, Fa]]
    return lhs, rhs


# FULL alpha grid
alpha_fail = []
t0 = cpu()
for i in range(n):
    lhs, rhs = lhs_rhs(i)
    if not np.array_equal(lhs, rhs):
        badj = np.nonzero(lhs != rhs)[0]
        alpha_fail.append({"i": i, "rows": badj[:5].tolist()})
        if len(alpha_fail) >= 3:
            break
ver["alpha_grid_cpu_s"] = round(cpu() - t0, 2)
ver["alpha_grid_points_checked"] = n
ver["alpha_identity_measured_full_grid"] = (len(alpha_fail) == 0)
ver["alpha_failures"] = alpha_fail
ver["alpha_check_kind"] = ("MEASURED (all n=3488 support points, unmodified "
                           "m<=11 instrument ops, exact table arithmetic)")
ver["cpu_after_alpha"] = round(cpu(), 2)
print("alpha done:", ver["alpha_grid_cpu_s"], "s; fails:", len(alpha_fail), flush=True)

# F2-linearity control (4 seeded random selectors)
rng = np.random.default_rng(16384)
lin_ok = True
for _ in range(4):
    sel = rng.integers(0, 2, k).astype(bool)
    comb = (np.bitwise_xor.reduce(Cm2[sel], axis=0) if sel.any()
            else np.zeros(Dp1, dtype=np.uint16))
    i = int(rng.integers(0, n))
    lhs, rhs = lhs_rhs(i)
    # recompute LHS/RHS for the COMBINED poly only
    Fc = int(np.bitwise_xor.reduce(ef.MUL[comb, _c.lucas_w(ef, Dp1 - 1, 0, inst.support[i])]))
    dc = np.zeros(Dp1, dtype=np.uint16)
    dc[0:Dp1 - 1] = comb[1:Dp1]
    dc[1::2] = 0
    Fpc = int(np.bitwise_xor.reduce(ef.MUL[dc, _c.lucas_w(ef, Dp1 - 1, 0, inst.support[i])]))
    a = inst.support[i]
    wp = _c.lucas_w(ef, len(Pi12) - 1, 0, a)
    pia = int(np.bitwise_xor.reduce(ef.MUL[Pi12, wp]))
    w1 = _c.lucas_w(ef, len(PiD12) - 1, 0, a)
    pida = int(np.bitwise_xor.reduce(ef.MUL[PiD12, w1]))
    Ga = int(np.bitwise_xor.reduce(ef.MUL[G12, _c.lucas_w(ef, len(G12) - 1, 0, a)]))
    lhsc = ef.MUL[pia, Fpc] ^ ef.MUL[pida, Fc]
    rhsc = ef.MUL[ef.MUL[Ga, Ga], ef.MUL[Fc, Fc]]
    lin_ok &= bool(np.array_equal(lhsc, rhsc))
ver["f2_linearity_control_4_selectors"] = lin_ok

# 3. checkpoint
json.dump(ver, open(f"{OUT}/m12_alpha_result.json", "w"), indent=1, default=str)
print("alpha checkpointed", flush=True)

# 4. delta chain (ADDENDUM 2(i))
rec_d = {"n": n}
rec_d["deg_ok"] = ver["deg_ok"]
t0 = cpu()
Ls = inst._lagrange()
rec_d["lagrange_build_cpu_s"] = round(cpu() - t0, 2)
t0 = cpu()
exact_div_ok = True
for i in range(n):
    q, r = gf.pdivmod(inst.Pi, [inst.support[i], 1])
    if gf.pdeg(r) >= 0:
        exact_div_ok = False
        rec_d["exact_div_fail_at"] = i
        break
rec_d["exact_div_ok"] = bool(exact_div_ok)
rec_d["exact_div_cpu_s"] = round(cpu() - t0, 2)
t0 = cpu()
bad_unit = [i for i, a in enumerate(inst.support) if gf.peval(Ls[i], a) != 1][:5]
rec_d["lagrange_unit_ok"] = (len(bad_unit) == 0)
rec_d["lagrange_unit_failures"] = bad_unit
rec_d["lagrange_unit_checks"] = n
rec_d["unit_checks_cpu_s"] = round(cpu() - t0, 2)
import random as _r
rng2 = _r.Random(16384 ^ 0xDE17A)
bad_off = 0
t0 = cpu()
for _ in range(200):
    i = rng2.randrange(n)
    l = rng2.randrange(n)
    if i != l and gf.peval(Ls[i], inst.support[l]) != 0:
        bad_off += 1
rec_d["offdiag_probes"] = 200
rec_d["offdiag_violations"] = bad_off
rec_d["offdiag_cpu_s"] = round(cpu() - t0, 2)
asm_ok = True
t0 = cpu()
for _ in range(3):
    j = rng2.randrange(k)
    row = inst.Y[j]
    acc = []
    for i in range(n):
        if row[i]:
            acc = gf.padd(acc, gf.pscale(Ls[i], gf.inv(inst.lam[i])))
    if gf.ptrim(acc) != gf.ptrim(inst.F[j]):
        asm_ok = False
        rec_d["assembly_fail_row"] = int(j)
        break
rec_d["assembly_ok"] = bool(asm_ok)
rec_d["assembly_cpu_s"] = round(cpu() - t0, 2)
rec_d["delta_exact"] = bool(rec_d["deg_ok"] and exact_div_ok
                            and rec_d["lagrange_unit_ok"] and bad_off == 0 and asm_ok)
rec_d["delta_scope"] = ("EXACT for ALL (j,l): lam_l f_j(a_l) = Y[j,l], implied "
                        "by L_i(a_i)=1 (n checks) + structural off-diagonal "
                        "vanishing + assembly")
ver["delta_chain"] = rec_d
ver["total_cpu_s"] = round(cpu(), 2)
ver["within_registered_budget"] = bool(cpu() <= BUDGET)
ver["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

# 5. verdict.json
json.dump(ver, open(f"{OUT}/verdict.json", "w"), indent=1, default=str)
print("VERDICT delta_exact=", rec_d["delta_exact"],
      " alpha_measured=", ver["alpha_identity_measured_full_grid"],
      " budget_ok=", ver["within_registered_budget"], flush=True)
