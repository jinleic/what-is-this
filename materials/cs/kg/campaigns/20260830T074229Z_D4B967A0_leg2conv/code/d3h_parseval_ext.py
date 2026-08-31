"""parseval_ext — extend the certified Parseval sum to the full grid-supported window
m in (251, 1255] (rows a+b <= 251 with m-b <= 5a still contribute EXACTLY from the grid).
Identity: paper line 884 (26): ||D3H||^2 = sum_{m>=1 odd} m^6 |b_m|^2 (D^3 multiplies the
m-th Fourier coefficient of H by (im)^3... wait: D^3(t^m) = m^3 t^m (Section 6 line 879),
so D^3H = sum m^3 b_m t^m over odd m and Parseval gives sum m^6 |b_m|^2 — paper line 884.
Tail beyond m = 1255: no grid row can contribute (window argument, notes §4), so the m>1255
part needs nongrid (a,b) rows — BOUNDED, not computed, in this script (report as (T2) bound
candidate, not part of the certified sum)."""
import math, sys, time, json
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
from flint import arb
import gate_A_final as G
import core

PREC = 256
M_EXT = 1255

if __name__ == "__main__":
    print("loading grid (balls)...")
    A_raw = G.load_grid()
    A = {}
    for k, v in A_raw.items():
        mid, half = G.A_to_arb(v)
        A[k] = (mid, half)
    print(f"rows = {len(A)}")
    print(f"computing b_m EXACTLY for m = 1..{M_EXT} odd (Arb, outward)... this is the expensive part")
    t0 = time.time()
    B = G.compute_b_m(A, m_max=M_EXT)
    print(f"b_m done ({time.time()-t0:.1f}s)")
    lo_sum = arb(0)
    hi_sum = arb(0)
    with core.Prec(PREC):
        for m in range(1, M_EXT + 1, 2):
            mid, half = B[m]
            mabs = abs(mid)
            up2 = (mabs + half) ** 2
            dn = (mabs - half)
            dn2 = (dn * dn) if not (dn < 0) else arb(0)
            w = arb(m) ** 6
            lo_sum += w * dn2
            hi_sum += w * up2
    print("Parseval grid-supported window (m <= 1255, odd), CERTIFIED:")
    print("  LOWER = ", lo_sum.str(50))
    print("  UPPER = ", hi_sum.str(50))
    print("  sqrt(LOWER) = ", lo_sum.sqrt().str(35))
    print("  sqrt(UPPER) = ", hi_sum.sqrt().str(35))
    with core.Prec(PREC):
        budget = arb("208.583976")
        ours = arb("624.611141")
    print(f"paper conclusion budget 208.583976: LOWER > budget?  {(lo_sum > budget)}")
    print(f"our-chain 624.611141: LOWER > ours? {(lo_sum > ours)}")
    json.dump({
        "window": "m <= 1255 odd, rows a+b <= 251",
        "lower": lo_sum.str(80),
        "upper": hi_sum.str(80),
        "sqrt_lower": lo_sum.sqrt().str(60),
        "sqrt_upper": hi_sum.sqrt().str(60),
        "lower_exceeds_paper_budget": bool(lo_sum > budget),
        "lower_exceeds_our_chain": bool(lo_sum > ours),
    }, open('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/parseval_ext1255.json', 'w'), indent=1)
    print("saved scratch/logs/parseval_ext1255.json")
