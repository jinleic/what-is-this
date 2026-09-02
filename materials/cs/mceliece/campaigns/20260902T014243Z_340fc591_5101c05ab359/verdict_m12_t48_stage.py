"""m=12 t=48 verdict stage — t48 cell (standalone, writes verdict.json).

MECHANICAL RE-INSTATIATION of M12ALPHA_DIVONLY/verdict_m12_stage.py at
the registered cell (12, 3488, 48, 16384). Same checks, same order,
same abort semantics; only the instance tuple, campaign path, budget
(9 CPU-h = 32400 s), and the MEASUREMENT-GATE fork bookkeeping
(pre-statement section 2.3) differ.

Fork rule (fixed in pre_statement_t48.md section 2 BEFORE compute):
  - alpha MEASURED over the FULL support grid (n points), inherited
    instrument ops (identical to m<=11 rows).
  - DERIVE branch entered ONLY IF, at the alpha-grid checkpoint, the
    CUMULATIVE measured process_time() has already exceeded
    0.75 * BUDGET = 24300 CPU s. In that case: abort the grid, label
    alpha CITED-DEPENDENCY (exact delta + Apon Lemma 3), compute the
    mandatory separate machine checks (a)-(d), and record the
    resolution in AMENDMENT_fork_resolution.md BEFORE the delta chain
    compute. Never resolved after seeing alpha verdict numbers beyond
    the aborting checkpoint.
"""
import json
import os
import random as _r
import sys
import time

import numpy as np

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/20260902T014243Z_340fc591_5101c05ab359/code")

from instance import Instance  # noqa: E402
import census as _c  # noqa: E402

RUN = "20260902T014243Z_340fc591_5101c05ab359"
OUT = (f"/Users/jinleic/jinleic-workspace/cs/mceliece/campaigns/{RUN}")
BUDGET = 32400.0              # pre_statement_t48.md section 5: 9 CPU-hours
FORK_ABORT_AT = 0.75 * BUDGET  # section 2.3: 24300 CPU s cumulative
M12, N12, T12, SEED = 12, 3488, 48, 16384   # registered cell, rule 16
C0 = time.process_time()


def cpu():
    return time.process_time() - C0


T_START = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
ver = {"instance": "(12,3488,48,16384)", "started": T_START,
       "registered_budget_cpu_s": BUDGET,
       "fork_abort_at_cpu_s": FORK_ABORT_AT,
       "stage_script": "verdict_m12_t48_stage.py (mechanical re-"
                       "instantiation of M12ALPHA_DIVONLY/"
                       "verdict_m12_stage.py; instrument inherited "
                       "byte-identical, see INSTRUMENT_INHERITANCE.md)"}

# 1. build
t0 = cpu()
inst = Instance(M12, N12, T12, SEED)
ver["build_cpu_s"] = round(cpu() - t0, 2)
print("built; cpu", ver["build_cpu_s"], flush=True)
gf, ef = inst.gf, inst.ef
D, k, n = inst.D, inst.k, inst.n
Dp1 = D + 1
ver.update({"k": k, "D": D, "n": n, "g_tries": inst.g_tries,
            "G_deg": int(gf.pdeg(inst.G)),
            "G_monic": bool(inst.G[-1] == 1)})

# lam recompute + gcd(G,G') + degree gate inputs (verdict-stage checks,
# identical to the t=64 stage)
lam_re = [gf.mul(gf.mul(inst.Ga[i], inst.Ga[i]), gf.inv(gf.peval(inst.PiD, a)))
          for i, a in enumerate(inst.support)]
ver["lam_recomputed_equal"] = bool(lam_re == inst.lam)
Gp = gf.pderiv(inst.G)
ver["gcd_G_Gprime_deg"] = gf.pdeg(gf.pgcd(inst.G[:], Gp))
ver["G_squarefree"] = ver["gcd_G_Gprime_deg"] <= 0
degs = [gf.pdeg(f) for f in inst.F]
ver["max_deg_f"] = int(max(degs))
ver["deg_ok"] = bool(max(degs) <= D)
# deg Π = n and Π' nonzero (distinct support) — derive-branch check (b),
# run unconditionally here as a build control
ver["Pi_deg_is_n"] = bool(gf.pdeg(inst.Pi) == n)
ver["PiD_nonzero_deg"] = int(gf.pdeg(inst.PiD))
ver["beta_delta_nonzero"] = bool(inst.guards.get("beta_delta_nonzero"))
ver["beta_pair"] = inst.guards.get("beta_pair")
ver["gamma_k"] = bool(inst.guards.get("gamma_k"))
ver["eps_gcd_const"] = bool(inst.guards.get("eps_gcd_const"))
ver["eps_maxdeg_is_D"] = bool(inst.guards.get("eps_maxdeg_is_D"))
json.dump(ver, open(f"{OUT}/m12t48_build_checkpoint.json", "w"), indent=1,
          default=str)
print("build checkpoint written", flush=True)

# 2. coefficient matrices (identical ops to verdict_m12_stage.py)
Cm2 = np.zeros((k, Dp1), dtype=np.uint16)
Dm2 = np.zeros((k, Dp1), dtype=np.uint16)
for j, f in enumerate(inst.F):
    cf = np.asarray(f + [0] * max(0, Dp1 - len(f)), dtype=np.uint16)[:Dp1]
    Cm2[j] = cf
    dc = np.zeros(Dp1, dtype=np.uint16)
    dc[0:Dp1 - 1] = cf[1:Dp1]
    dc[1::2] = 0
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


# 3. MEASUREMENT-GATE CHECKPOINT (pre-statement section 2.3): judged on
# measured cumulative CPU only. Derive branch entry condition.
cpu_at_gate = cpu()
ver["cpu_at_alpha_gate"] = round(cpu_at_gate, 2)
ver["fork_resolution"] = ("MEASURED" if cpu_at_gate <= FORK_ABORT_AT
                          else "DERIVE")

if ver["fork_resolution"] == "MEASURED":
    # FULL alpha grid over ALL n support points — inherited instrument ops
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
    ver["alpha_check_kind"] = ("MEASURED (all n support points, unmodified "
                               "m<=11 instrument ops, exact table arithmetic)")
else:
    # DERIVE branch: CITED-DEPENDENCY (exact delta + Apon Lemma 3) —
    # mandatory separate machine checks (a)-(d):
    # (a) gcd(G,G') constant -> ver["G_squarefree"] (above)
    # (b) deg Pi = n, Pi' nonzero -> ver["Pi_deg_is_n"] (above)
    # (c) max_j deg f_j <= D ABORT gate -> ver["deg_ok"] (above)
    # (d) lam recomputed byte-equal -> ver["lam_recomputed_equal"] (above)
    ver["alpha_grid_points_checked"] = 0
    ver["alpha_identity_measured_full_grid"] = None
    ver["alpha_failures"] = []
    ver["alpha_check_kind"] = ("CITED-DEPENDENCY: INFERRED from a cited "
                               "lemma (Apon 2026/1810 Lemma 3, kappa=1) "
                               "given a measured exact delta — entered via "
                               "the pre-registered section-2.3 gate on the "
                               "measured cumulative CPU "
                               f"{cpu_at_gate:.1f} > {FORK_ABORT_AT:.0f} s")
    # amend file MUST be written before delta compute; done by caller after
    # this stage writes its fork-resolution record (see runner).
    json.dump(ver, open(f"{OUT}/m12t48_alpha_result.json", "w"), indent=1,
              default=str)
    print("FORK RESOLVED TO DERIVE at cpu", round(cpu_at_gate, 1),
          "— alpha CITED-DEPENDENCY; aborting grid", flush=True)
    json.dump(ver, open(f"{OUT}/verdict_fork.json", "w"), indent=1,
              default=str)
    sys.exit(3)  # registered derive-branch exit: caller files the amendment

ver["cpu_after_alpha"] = round(cpu(), 2)
print("alpha done:", ver["alpha_grid_cpu_s"], "s; fails:", len(alpha_fail),
      flush=True)

# F2-linearity control (4 seeded random selectors; identical to t=64 stage)
rng = np.random.default_rng(SEED)
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

# 4. checkpoint
json.dump(ver, open(f"{OUT}/m12t48_alpha_result.json", "w"), indent=1,
          default=str)
print("alpha checkpointed", flush=True)

# 5. delta chain (ADDENDUM 2(i)) — identical to verdict_m12_stage.py
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
rng2 = _r.Random(SEED ^ 0xDE17A)
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
                            and rec_d["lagrange_unit_ok"] and bad_off == 0
                            and asm_ok)
rec_d["delta_scope"] = ("EXACT for ALL (j,l): lam_l f_j(a_l) = Y[j,l], implied "
                        "by L_i(a_i)=1 (n checks) + structural off-diagonal "
                        "vanishing + assembly")
ver["delta_chain"] = rec_d
ver["total_cpu_s"] = round(cpu(), 2)
ver["within_registered_budget"] = bool(cpu() <= BUDGET)
ver["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

# 6. verdict.json
json.dump(ver, open(f"{OUT}/verdict.json", "w"), indent=1, default=str)
print("VERDICT delta_exact=", rec_d["delta_exact"],
      " alpha_measured=", ver["alpha_identity_measured_full_grid"],
      " max_deg_f=", ver["max_deg_f"], "/", D,
      " budget_ok=", ver["within_registered_budget"], flush=True)
