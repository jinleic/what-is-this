
#!/usr/bin/env python3
"""Gate C continuation: DRAT/LRAT certificates for all 18 floor-impossibility
instances of the named set S survivors ({paper55, sun56} x sigma-orbit x 3 sides).
Encoding: src/gate_b_encodings.py:build_floor_cnf (the exact encoding
checker-certified in campaign 2026-08-30T031544Z). Solver: kissat 4.0.4
(binary DRAT). Conversion: drat-trim -L -> LRAT. Checkers: drat-trim and
lrat-check built from tools_snapshot 2e3b2dc."""
import sys, json, subprocess, time, hashlib, os
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/src")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/scratch")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/scripts")
from gatec_decomps import load_paper55, load_sun56
import gatec_sweep as gs
from gate_b_floor import prep, subset_dfs
from gate_b_encodings import build_floor_cnf

CAMP = "/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T174243Z_edbdd408-840f-42b4-9be0-192a95783ab9"
CNFDIR = os.path.join(CAMP, "cnf"); PR = os.path.join(CAMP, "proofs"); CHK = os.path.join(CAMP, "checker_output")
for d in (CNFDIR, PR, CHK): os.makedirs(d, exist_ok=True)
KISSAT = "/opt/homebrew/bin/kissat"
DRATTRIM = "/tmp/dratbuild/drat-trim"
LRATCHECK = "/tmp/dratbuild/lrat-check"
TMPD = "/tmp/dratbuild"

def wr(path, txt):
    with open(path, "w") as f: f.write(txt)

results = {}
instances = []
instances = []
for D, loader in (("paper55", load_paper55), ("sun56", load_sun56)):
    U, V, W = loader()
    orbit = gs.sigma_orbit(U, V, W)
    names = ["L", "R", "Ofac"]
    for sp, triple in enumerate(orbit):
        for si, targets in enumerate(triple):
            instances.append((D, sp, names[si], targets))
for (D, sp, side, targets) in instances:
    classes, reps = prep(list(targets))
    d = len(classes)
    ok, stats, order = subset_dfs(classes, reps)
    if ok:
        results[f"{D}|s{sp}|{side}"] = {"d": d, "floor_achievable": True, "note": "no UNSAT instance; witness exists"}
        print(D, sp, side, "d=", d, "FLOOR ACHIEVABLE (witness) — no certificate needed", flush=True)
        continue
    clauses, varmap, nv = build_floor_cnf(classes, reps, d)
    dimacs = [f"p cnf {nv} {len(clauses)}"] + [" ".join(map(str, c)) + " 0" for c in clauses]
    base = f"{D}_s{sp}_{side}_d{d}"
    cnfp = os.path.join(CNFDIR, base + ".cnf")
    wr(cnfp, "\n".join(dimacs) + "\n")
    prop = os.path.join(PR, base + ".drat")
    t0 = time.time()
    r1 = subprocess.run([KISSAT, "-q", cnfp, prop], capture_output=True, text=True, timeout=7200)
    dt = time.time() - t0
    unsat = "s UNSATISFIABLE" in r1.stdout or "UNSATISFIABLE" in r1.stdout
    lr = os.path.join(PR, base + ".lrat")
    r2 = subprocess.run([DRATTRIM, cnfp, prop, "-L", lr], capture_output=True, text=True)
    r3 = subprocess.run([DRATTRIM, cnfp, lr], capture_output=True, text=True)
    r4 = subprocess.run([LRATCHECK, cnfp, lr], capture_output=True, text=True)
    key = f"{D}|s{sp}|{side}"
    results[key] = {
        "d": d, "floor_achievable": False, "vars": nv, "clauses": len(clauses),
        "kissat": r1.stdout.strip()[-80:], "unsat": unsat, "kissat_s": round(dt,2),
        "drattrim_lrat_conv": (r2.stdout + r2.stderr).strip()[-120:],
        "drattrim_verify": (r3.stdout + r3.stderr).strip()[-80:], "drattrim_rc": r3.returncode,
        "lratcheck_verify": (r4.stdout + r4.stderr).strip()[-80:], "lratcheck_rc": r4.returncode,
        "sha256_cnf": hashlib.sha256(open(cnfp,"rb").read()).hexdigest(),
        "sha256_lrat": hashlib.sha256(open(lr,"rb").read()).hexdigest(),
    }
    print(key, "d=", d, "vars", nv, "UNSAT", unsat, "| drat:", r3.stdout.strip()[-20:], "rc", r3.returncode,
          "| lrat:", r4.stdout.strip()[-20:], "rc", r4.returncode, f"({dt:.1f}s)", flush=True)

with open(os.path.join(CAMP, "certificates_result.json"), "w") as f:
    json.dump(results, f, indent=1)
print("ALL DONE")
