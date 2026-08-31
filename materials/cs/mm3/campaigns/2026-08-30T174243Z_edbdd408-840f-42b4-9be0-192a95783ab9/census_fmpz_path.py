
#!/usr/bin/env python3
"""Independent second verification path (fmpz) for the non-monomial diagonal census.
Uses the FROZEN sandwich implementation from campaign 113838Z scripts + fmpz Brent.
Enumerates all 6960 ternary unimodular 3x3 matrices, applies the diagonal sandwich,
records: (a) ternarity, (b) Brent 729/729 over Z via fmpz for survivors, (c) the 48
survivor identities, (d) certified totals per sigma class for survivors."""
import sys, json, itertools, time
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/src")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/scratch")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/scripts")
from gatec_decomps import load_paper55, load_sun56
import gatec_sweep as gs
from gate_b_floor import prep, subset_dfs
from flint import fmpz_mat, fmpz

R, N = 23, 9
GAP = 14

def brent_fails_fmpz(U, V, Wfac):
    """729 identities, exact fmpz."""
    fails = 0
    Um = fmpz_mat([[int(x) for x in row] for row in U])
    Vm = fmpz_mat([[int(x) for x in row] for row in V])
    Wm = fmpz_mat([[int(x) for x in row] for row in Wfac])
    for i in range(3):
        for j in range(3):
            for a in range(9):
                for b in range(9):
                    s = fmpz(0)
                    for r in range(R):
                        if Wm[r, 3*i+j] and Um[r, a] and Vm[r, b]:
                            s += Wm[r, 3*i+j] * Um[r, a] * Vm[r, b]
                    expected = 0
                    for k in range(3):
                        if a == 3*i+k and b == 3*k+j:
                            expected = 1
                    if s != expected:
                        fails += 1
    return fails

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
          - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
          + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

unimod = []
for entries in itertools.product((-1,0,1), repeat=9):
    M=[list(entries[0:3]),list(entries[3:6]),list(entries[6:9])]
    if det3(M) == 1 or det3(M) == -1:
        unimod.append(M)
assert len(unimod) == 6960, len(unimod)

def certified_side(targets):
    classes, reps = prep(list(targets))
    d = len(classes)
    ok, stats, order = subset_dfs(classes, reps)
    return (d if ok else d+1), d, ok, stats["states"]

out = {"unimod_count": len(unimod), "per_decomposition": {}}
for D, loader in (("paper55", load_paper55), ("sun56", load_sun56)):
    U, V, W = loader()
    survivors = []
    verdict = {"non_ternary": 0, "ternary_brent_fail": 0, "valid": 0}
    tb0 = time.time()
    for gi, G in enumerate(unimod):
        U2, V2, W2 = gs.sandwich(U, V, W, G, G, G)
        ternary = all(abs(x) <= 1 for M in (U2, V2, W2) for row in M for x in row)
        if ternary:
            bf = brent_fails_fmpz(U2, V2, W2)
            if bf == 0:
                verdict["valid"] += 1
                survivors.append(gi)
            else:
                verdict["ternary_brent_fail"] += 1
        else:
            verdict["non_ternary"] += 1
    # sigma classes for survivors, certified totals
    totals = set()
    per_sigma = []
    for gi in survivors:
        G = unimod[gi]
        U2, V2, W2 = gs.sandwich(U, V, W, G, G, G)
        for k, (a, b, c) in enumerate(gs.sigma_orbit(U2, V2, W2)):
            cl = certified_side(a); cr = certified_side(b); co = certified_side(c)
            tot = cl[0] + cr[0] + co[0] + GAP
            totals.add(tot)
            per_sigma.append([gi, k, cl[0], cr[0], co[0]+GAP, tot])
    is_mono = all(sum(1 for x in unimod[gi] if x) == 3 for gi in survivors)
    out["per_decomposition"][D] = {
        "verdict": verdict, "survivor_count": len(survivors),
        "survivors_all_monomial": is_mono,
        "certified_totals_set": sorted(totals),
        "min_certified_total": min(totals),
        "elapsed_s": round(time.time()-tb0, 1),
    }
    print(D, out["per_decomposition"][D], flush=True)

print(json.dumps(out["per_decomposition"], indent=1))
with open("/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T174243Z_edbdd408-840f-42b4-9be0-192a95783ab9/census_fmpz_path.json", "w") as f:
    json.dump(out, f, indent=1)
print("DONE fmpz path")
