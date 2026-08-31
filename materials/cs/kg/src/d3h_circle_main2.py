# Slower-t-focused continuation: θ in {1.0, 1.2, 1.4, 1.5, 1.57} — where S decays; plus endpoint probes
# via |t|→1 circle points near ±i where |1−ρ²| gets small.
import sys, math, time
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/kg/src")
from d3h_circle import S_c
print("fine θ sweep:", flush=True)
for th in [0.5, 1.0, 1.3, 1.5, 1.57, 1.8, 2.2]:
    t = complex(math.cos(th), math.sin(th))
    t0 = time.time()
    a, b, s = S_c(t)
    r = (t - (0.34101124**2)*t**3 + (0.05276111**2)*t**5)/V if False else None
    print(f"θ={th:.3f}: S = {a.real:.4f} + 4·{b.real/4:.4f} → {s.real:.4f} [{time.time()-t0:.0f}s]", flush=True)
