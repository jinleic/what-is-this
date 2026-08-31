import sys, math, time
sys.path.insert(0, ".")
from d3h_inner import S_float
for t in [0.0, 0.25, 0.5, 0.75, 0.9]:
    t0 = time.time()
    a, b, s = S_float(t)
    print(f"S({t}) = {a:.6f} + 4·{b/4:.6f} = {s:.6f}   ({time.time()-t0:.0f}s)", flush=True)
print("budget (π²/6)·126.804 =", math.pi**2/6*126.80385221)
