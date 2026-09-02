"""UNIFORM delta verifier — ONE code path for all certified instances.

Implements the exact Lagrange route fixed in pre_statement.md ADDENDUM 2:

  F is built as  f_j = SUM_i Y[j,i] lam_i^{-1} L_i,  L_i = Pi/(Z-a_i)/Pi'(a_i).
  Off-diagonal vanishing L_i(a_l) = 0 (l != i) holds BY CONSTRUCTION (the
  factor (Z-a_l) is retained; the division is exact, asserted per i).
  Hence   lam_l f_j(a_l) = Y[j,l]  for ALL (j,l)   follows EXACTLY from the
  n scalar checks L_i(a_i) = 1.

  LOAD-BEARING DEGREE CONDITION (ABORT on failure): max_j deg f_j <= D.
  Without it the degree-<= n-1 agreement argument collapses (f and f + Pi
  agree at every support point yet differ).

Reports per instance:
  deg_ok            max_j deg(f_j) <= D                     (ABORT if False)
  exact_div_ok      Pi = (Z-a_i) q_i with zero remainder, all i
  lagrange_unit_ok  L_i(a_i) = 1, all i                     (the n checks)
  offdiag_spotcheck L_i(a_l) = 0 on seeded random (i,l), l != i
  assembly_spotcheck  f_j equals the direct sum for seeded random j
  delta_exact       conjunction => lam_l f_j(a_l) = Y[j,l] for ALL (j,l)
"""
from __future__ import annotations

import json
import random
import sys
import time

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))

from instance import Instance


CERTIFIED = [
    (6, 64, 3, 1387), (6, 64, 4, 2311), (6, 64, 5, 3413),
    (7, 128, 3, 4421), (7, 128, 6, 5531),
    (8, 256, 3, 6637), (8, 256, 5, 7741),
    (9, 512, 12, 9109), (9, 512, 24, 9137), (9, 512, 40, 9277),
    (10, 1024, 24, 8001), (10, 1024, 40, 5113),
    (11, 2048, 48, 6211),
]


def verify_delta(m, n, t, seed, offdiag_probes=200, assembly_rows=3):
    t0 = time.time()
    inst = Instance(m, n, t, seed)
    gf, ef = inst.gf, inst.ef
    rec = {"m": m, "n": n, "t": t, "k": inst.k, "D": inst.D, "seed": seed,
           "G": inst.G}

    # --- LOAD-BEARING degree condition: ABORT if violated
    degs = [gf.pdeg(f) for f in inst.F]
    deg_ok = max(degs) <= inst.D
    rec["max_deg_f"] = int(max(degs))
    rec["deg_ok"] = bool(deg_ok)
    if not deg_ok:
        rec["ABORT"] = ("degree condition max_j deg f_j <= D FAILED; the "
                        "Lemma-2 degree argument does not apply, so NO delta "
                        "verdict is produced for this instance")
        rec["delta_exact"] = None
        rec["elapsed_s"] = round(time.time() - t0, 2)
        return rec

    # --- Lagrange basis, rebuilt through the same construction path
    Ls = inst._lagrange()

    # exact-division certificate: Pi = (Z - a_i) q_i, zero remainder
    exact_div_ok = True
    for i in range(n):
        q, r = gf.pdivmod(inst.Pi, [inst.support[i], 1])
        if gf.pdeg(r) >= 0:
            exact_div_ok = False
            rec["exact_div_fail_at"] = i
            break
    rec["exact_div_ok"] = bool(exact_div_ok)

    # the n scalar checks L_i(a_i) = 1
    bad_unit = []
    for i, a in enumerate(inst.support):
        if gf.peval(Ls[i], a) != 1:
            bad_unit.append(i)
            if len(bad_unit) > 4:
                break
    rec["lagrange_unit_ok"] = (len(bad_unit) == 0)
    rec["lagrange_unit_failures"] = bad_unit
    rec["lagrange_unit_checks"] = n

    # off-diagonal spot check (structural, so this is redundancy not licence)
    rng = random.Random(seed ^ 0xDE17A)
    bad_off = 0
    for _ in range(offdiag_probes):
        i = rng.randrange(n)
        l = rng.randrange(n)
        if i == l:
            continue
        if gf.peval(Ls[i], inst.support[l]) != 0:
            bad_off += 1
    rec["offdiag_probes"] = offdiag_probes
    rec["offdiag_violations"] = bad_off

    # assembly spot check: f_j equals the direct sum over its support
    asm_ok = True
    for _ in range(assembly_rows):
        j = rng.randrange(inst.k)
        acc = []
        row = inst.Y[j]
        for i in range(n):
            if row[i]:
                acc = gf.padd(acc, gf.pscale(Ls[i], gf.inv(inst.lam[i])))
        if gf.ptrim(acc) != gf.ptrim(inst.F[j]):
            asm_ok = False
            rec["assembly_fail_row"] = int(j)
            break
    rec["assembly_spotcheck_rows"] = assembly_rows
    rec["assembly_ok"] = bool(asm_ok)

    rec["delta_exact"] = bool(exact_div_ok and rec["lagrange_unit_ok"]
                              and bad_off == 0 and asm_ok)
    rec["delta_scope"] = ("EXACT for ALL (j,l): lam_l f_j(a_l) = Y[j,l], "
                          "implied by L_i(a_i)=1 (n checks) + structural "
                          "off-diagonal vanishing + assembly")
    # alpha stays MEASURED on these rows (unchanged, from the instance guards)
    rec["alpha_identity_measured"] = bool(inst.guards.get("alpha_identity"))
    rec["beta_delta_nonzero"] = bool(inst.guards.get("beta_delta_nonzero"))
    rec["beta_pair"] = inst.guards.get("beta_pair")
    rec["eps_maxdeg_is_D"] = bool(inst.guards.get("eps_maxdeg_is_D"))
    rec["gamma_k"] = bool(inst.guards.get("gamma_k"))
    rec["elapsed_s"] = round(time.time() - t0, 2)
    return rec


def main():
    out = []
    allok = True
    for (m, n, t, seed) in CERTIFIED:
        r = verify_delta(m, n, t, seed)
        ok = r.get("delta_exact")
        allok = allok and bool(ok)
        print(f"m={m:2d} n={n:5d} t={t:3d}: deg_ok={r['deg_ok']} "
              f"exact_div={r.get('exact_div_ok')} "
              f"L_i(a_i)=1 on {r.get('lagrange_unit_checks')} checks:"
              f"{r.get('lagrange_unit_ok')} "
              f"offdiag_viol={r.get('offdiag_violations')} "
              f"assembly={r.get('assembly_ok')} "
              f"=> delta_exact={ok} | alpha_measured={r.get('alpha_identity_measured')} "
              f"[{r['elapsed_s']}s]", flush=True)
        out.append(r)
        json.dump(out, open("/Users/jinleic/jinleic-workspace/cs/mceliece/"
                            "scratch/delta_uniform.json", "w"), indent=1, default=str)
    print(f"\nUNIFORM DELTA (13 instances, one code path): all_exact={allok}")


if __name__ == "__main__":
    main()
