"""parseval_direct — Main's decider: ||D3H||^2_{L^2(T)} DIRECTLY as a Parseval sum of H's own
coefficients b_m on the circle, CERTIFIED interval arithmetic. Zero overlap with the contested
Phi3/chain machinery (no Phi3, no c-table, no kernels, no theta-quadrature).

Parseval identity used: paper line 884 (eq 26), verbatim: "||D^3H||^2_{L^2(T)} = sum_{m>=1,
m odd} m^6 |b_m|^2" (from Parseval applied to D^3 = m^3 multiplier on odd modes; D = t d/dt).
Wait: D^3 H has coefficients (im)^3 b_m, so |D^3H|^2 = sum (m^3)^2 |b_m|^2 = sum m^6 |b_m|^2,
as the paper writes. So the paper's identity with the m^6 weight appears at line 884 (26).

b_m sources:
  * m in [1, 251]: EXACT from the certified A-grid (gate_A_final.compute_b_m machinery),
    with interval widths in Arb at outward rounding.
  * m > 251: bound via Lemma 6.2-type reasoning? NO — we need the TRUE m^6-weighted tail,
    not a tail on |b_m|. Instead bound tail directly:
      sum_{m>251, odd} m^6 |b_m|^2: from the paper's own Lemma 6.2 the |b_m| tails are
      controlled by ||D^3H|| itself — circular. INSTEAD: bound each b_m for m > 251 via the
      coefficient-extraction bound: b_m = (pi/2) sum_{(a,b)} (-1)^b A^2 [t^{m-b}] rho^a only
      for (a,b) in the grid — but the grid has a+b <= 251, while [t^{m-b}]rho^a requires
      m-b <= 5a i.e. a >= (m-b)/5: for m > 251 the coefficient [t^{m-b}]rho^a needs a > 50ish.
      The GRID A_{a,b} only covers a+b <= 251: for b <= 251-a: still exact for a+b <= 251!
      So the grid gives b_m EXACTLY on the window where a+b <= 251 AND a <= (m-b)/1...
      For m > 251: rows with a+b <= 251 but m-b <= 5a+b are needed: a >= (m - b - ...).
      For m = 253: a+b <= 251 with m-b in rho^a support: a + b <= 251 AND m - b <= 5a.
      For a=51, b=0: m - 0 = 253 <= 5*51 = 255 ✓. So rows DO contribute at m = 253.
      The missing mass comes from (a,b) with a+b > 251 — bounded by notes section 9 machinery
      (per-slice masses mu_a -> 0, falling-factorial bounds). That is a big lift; BUT for the
      DIRECT decision we only need ||D3H||^2 to within a factor of 3: 208.58 vs 624.61. The
      GRID-ONLY partial sum (m <= 251 with rows a+b <= 251) is a LOWER BOUND on ||D3H||^2
      (all terms m^6 |b_m|^2 >= 0). If the lower bound alone EXCEEDS 208.583976, the paper's
      conclusion ||D3H|| < 14.4424 is contradicted, regardless of the tail. If it does not,
      we must bound the tail to decide.

CERTIFIED ARITHMETIC: every b_m in [1, 251] from the Arb ball grid, outward-rounded; the sum
of m^6 |b_m|^2 accumulated with interval arithmetic; the result is a certified LOWER bound.
"""
import math, sys, time, json
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
from flint import arb, fmpq
import gate_A_final as G
import core

PREC = 256

import os
CACHE = '/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/parseval_bm_cache.json'
def load_grid_balls():
    A_raw = G.load_grid()
    A = {}
    for k, v in A_raw.items():
        mid, half = G.A_to_arb(v)
        A[k] = (mid, half)
    return A



if __name__ == "__main__":
    print("loading grid...")
    t0 = time.time()
    A = load_grid_balls()
    print("computing b_m EXACTLY for m = 1..251 (Arb, outward)...")
    t0 = time.time()
    B = G.compute_b_m(A, m_max=251)
    print(f"b_m done ({time.time()-t0:.1f}s)")
    # Parseval: sum m^6 |b_m|^2 with |b_m|^2 as interval square:
    # lower bound of |b_m|^2: (|mid| - half)^2 if positive else 0; upper: (|mid|+half)^2
    lo_sum = arb(0)
    hi_sum = arb(0)
    with core.Prec(PREC):
        for m in range(1, 252, 2):
            mid, half = B[m]
            mabs = abs(mid)
            up2 = (mabs + half) ** 2
            dn = (mabs - half)
            dn2 = (dn * dn) if not (dn < 0) else arb(0)
            w = arb(m) ** 6
            lo_sum += w * dn2
            hi_sum += w * up2
    print("Parseval grid-part (m <= 251, odd):")
    print("  LOWER bound = ", lo_sum.str(40))
    print("  UPPER bound = ", hi_sum.str(40))
    print("  sqrt(LOWER) = ", lo_sum.sqrt().str(30))
    print("  sqrt(UPPER) = ", hi_sum.sqrt().str(30))
    paper_budget = core.Prec(PREC) and None
    with core.Prec(PREC):
        budget = arb(208) + arb("0.583976")
    print(f"paper's conclusion budget: ||D3H||^2 <= {budget.str(20)}")
    print(f"LOWER > budget?  {(lo_sum > budget)}")
    print(f"our chain implies 624.611141: grid-part share = ", end="")
    with core.Prec(PREC):
        six_hundred = arb("624.611141")
    print((lo_sum / six_hundred).str(12))
    json.dump({
        "lower": lo_sum.str(60),
        "upper": hi_sum.str(60),
        "paper_budget": "208.583976",
        "lower_exceeds_paper_budget": bool(lo_sum > budget),
    }, open('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/parseval_gridpart.json', 'w'), indent=1)
    print("saved scratch/logs/parseval_gridpart.json")
