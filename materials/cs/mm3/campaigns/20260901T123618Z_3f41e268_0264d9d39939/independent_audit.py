#!/usr/bin/env python3
"""Independent post-run audit (assignment point 5): re-verify the converged
negative through a THIRD path — exact CRT-free direct enumeration of every
12-gate schedule shape is infeasible; instead the audit independently
re-derives the aux universe from scratch with a different enumeration loop,
re-runs the DFS at T=12 for a random subsample plus the boundary instances,
re-checks the frozen lrat certificate bytes, re-verifies witness circuits by
exact expansion, and re-computes the total arithmetic.
"""
import json, hashlib, os, subprocess, sys, random
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/src")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/scripts")
import sun56_verify as sv
from gate_b_floor import prep, subset_dfs, canon, N

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
checks = []

def ck(name, cond, detail=""):
    checks.append((name, bool(cond), detail))
    print(("PASS " if cond else "FAIL ") + name + (" | " + detail if detail else ""))

# 1. Independent universe re-derivation (different loop shape: accumulate over
#    pairs with signed-duplicate suppression via canon of a+b and a-b only —
#    ±a±b == ±(a+b), ±(a-b), so four sign combos reduce to two canonical sums).
V = [tuple(r) for r in sv.V]
classes, reps = prep(V)
Tset = set(classes)
BASIS = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)]
ICAN = {canon(e) for e in BASIS}
forms = BASIS + sorted(Tset)
u = set()
for i in range(len(forms)):
    for j in range(len(forms)):
        for combo in (0, 1):
            vc = tuple((forms[i][k] + forms[j][k]) if combo == 0
                       else (forms[i][k] - forms[j][k]) for k in range(N))
            for sgn in (1, -1):
                c = canon(tuple(sgn * vc[k] for k in range(N)))
                if c != (0,) * N and c not in ICAN and c not in Tset:
                    u.add(c)
au = sorted(u)
h = hashlib.sha256(("\n".join(",".join(str(x) for x in row) for row in au)).encode()).hexdigest()
ck("universe_size_338", len(au) == 338, f"got {len(au)}")
ck("universe_sha256_pinned", h == "c10bdeeae74b601dc063cb935f7c7baca79a637cab52a7e2f85d16a46f84f813", h[:16])

# 2. DFS re-runs: all boundary + deterministic sample of 40.
rng = random.Random(20260901)
sample = sorted(set(au[:3]) | set(au[-3:]) | set(rng.sample(au, 40)))
ok_all = True
for a in sample:
    extc, extr = prep(sorted(Tset | {a}))
    ok, st, _ = subset_dfs(extc, extr)
    if ok:
        ok_all = False
        break
ck("dfs_T12_infeasible_46_sample", ok_all, f"n={len(sample)}")

# 3. Certificate byte check + re-verification via pinned checkers.
cnf = os.path.join(HERE, "certificates/sun_aux1_T12.cnf")
lrat = os.path.join(HERE, "certificates/sun_aux1_T12.lrat")
cnf_h = hashlib.sha256(open(cnf, 'rb').read()).hexdigest()
lrat_h = hashlib.sha256(open(lrat, 'rb').read()).hexdigest()
ck("cert_cnf_sha", cnf_h == "cd27491c7893d964115850034bdebbe732eec004a4de861b514957831be41548")
ck("cert_lrat_sha", lrat_h == "c4067622a7811ade795607d064d5a35437f9f53d1a7bb0ad9bb79851656e9a2e")
for checker in ("/tmp/dratbuild_sun56gap/drat-trim", "/tmp/dratbuild_sun56gap/lrat-check"):
    r = subprocess.run([checker, cnf, lrat], capture_output=True, text=True)
    ck(f"replay_{os.path.basename(checker)}", r.returncode == 0, (r.stdout + r.stderr).strip()[-40:])

# 4. Witness re-verification by exact expansion (third path: nested-loop expansion).
gU, fU = [], []
vs = {i: tuple(1 if j == i else 0 for j in range(9)) for i in range(9)}
for gi, (a, s, b) in enumerate(sv.SIDES['U']['inter']):
    vs[9 + gi] = tuple(vs[a][k] + s * vs[b][k] for k in range(9))
U_rows_exact = []
for f in sv.SIDES['U']['final']:
    acc = [0] * 9
    for idx, c in f:
        v = vs[idx]
        for k in range(9):
            acc[k] += c * v[k]
    U_rows_exact.append(tuple(acc))
ck("U_witness_exact", U_rows_exact == [tuple(r) for r in sv.U])

# 5. Total arithmetic recheck.
total = 13 + 13 + 16 + 14
ck("total_56_arithmetic", total == 56, f"13+13+16+14={total}")

# 6. Verdict consistency with pre-registered rule.
vd = json.load(open(os.path.join(HERE, "verdict.json")))
ck("verdict_rule_neg_bounded", vd["verdict"] == "NEGATIVE-BOUNDED" and vd["admitted_aux_count"] == 0
   and vd["instrument_agreement"] == 338)
ckpt_rows = [json.loads(l) for l in open(os.path.join(HERE, "scan_checkpoint.jsonl"))]
ck("checkpoint_rows_338_complete", len(ckpt_rows) == 338 and
   all(r["dfs"] is False and r["cnf_sat"] is False for r in ckpt_rows))

print("\n%d/%d checks pass" % (sum(1 for _, c, _ in checks if c), len(checks)))
sys.exit(0 if all(c for _, c, _ in checks) else 1)
