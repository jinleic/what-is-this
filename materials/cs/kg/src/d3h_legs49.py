"""legs49 — recompute leg1 (‖Φ3‖²) and leg2 (‖∂x∂yΦ3‖²) on the 49-point θ grid; save npy."""
import sys, math, time
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/kg/src")
import numpy as np
from numpy.polynomial.hermite import hermgauss
from d3h_phi3 import G_pa_terms, _IO_LIST
from d3h_kernels import mult_terms
from d3h_dxdy import weight_eval
from flint import acb

V = 1+0.34101124**2+0.05276111**2
VT = 0.136419125/math.sqrt(V)
def wmid(w, xf): return float(weight_eval(w, acb(xf)).mid().real)
def rho_c(t): return (t-(0.34101124**2)*t**3+(0.05276111**2)*t**5)/V
def cfun(pa, t):
    s3 = 0.34101124; s5 = 0.05276111
    r1 = (1-3*s3*s3*t*t+5*s5*s5*t**4)/V
    r2 = (6*s3*s3*t*t+20*s5*s5*t**4)/V
    r3 = (12*s3*s3*t*t+80*s5*s5*t**4)/V
    return {(1,0): r3, (2,0):3*r1*r2, (3,0):r1**3, (0,1):-t, (0,2):3*t*t, (0,3):-(t**3),
            (1,1):-3*t*(r1+r2), (2,1):-3*t*r1*r1, (1,2):3*t*t*r1}[pa]

NN = 140
yg, wg = hermgauss(NN)
xs = math.sqrt(2)*yg; w = wg/math.sqrt(math.pi)
msk = np.abs(xs) <= 14
xs, w = xs[msk], w[msk]
ux = VT*(xs**3-3*xs)/math.sqrt(6)
Wij = np.outer(w, w)
_wv = {}
def wvec(wt):
    if wt not in _wv: _wv[wt] = np.array([wmid(wt, x) for x in xs])
    return _wv[wt]

th49 = np.load("/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/d3h_circle_thfine49.npy")
L1 = np.zeros(len(th49)); L2 = np.zeros(len(th49))
t0 = time.time()
for idx, th in enumerate(th49):
    t = complex(math.cos(th), math.sin(th))
    r = rho_c(t); s = 1-r*r
    Km = np.exp(-(ux[:,None]**2-2*r*ux[:,None]*ux[None,:]+ux[None,:]**2)/(2*s))/(2*math.pi*np.sqrt(s))
    Ksq = np.sqrt(Km)
    Ph1 = np.zeros((len(xs), len(xs)), dtype=complex)
    Ph2 = np.zeros((len(xs), len(xs)), dtype=complex)
    for (p, a) in _IO_LIST:
        cval = cfun((p, a), t)
        for (aeff, store) in [(a, Ph1), (a+1, Ph2)]:
            for (wx, wy, kind, tp) in G_pa_terms(p, aeff):
                M = np.zeros((len(xs), len(xs)), dtype=complex)
                for (au, av, ar, pw), cc in mult_terms(kind).items():
                    rp = (r**ar) if ar else 1.0
                    sp = (s**(-pw)) if pw else 1.0
                    up = (ux**au) if au else np.ones_like(ux)
                    vp = (ux**av) if av else np.ones_like(ux)
                    M += cc*rp*sp*np.outer(up, vp)
                store += cval*(2*math.pi)*(VT**tp)*np.outer(wvec(wx), wvec(wy))*M*Ksq
    L1[idx] = float(np.real(np.sum(Wij*np.abs(Ph1)**2)))
    L2[idx] = float(np.real(np.sum(Wij*np.abs(Ph2)**2)))
    if idx % 8 == 0:
        print(f"{idx}/49 θ={th:.3f} l1={L1[idx]:.3f} l2={L2[idx]:.3f} ({time.time()-t0:.0f}s)", flush=True)
np.save("/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/d3h_leg1_49.npy", L1)
np.save("/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/d3h_leg2_49.npy", L2)
h = th49[1]-th49[0]
a1 = h*(L1[0]/2 + L1[-1]/2 + L1[1:-1].sum())/math.pi
a2 = h*(L2[0]/2 + L2[-1]/2 + L2[1:-1].sum())/math.pi
print(f"DONE {(time.time()-t0):.0f}s  avg(leg1) = {a1:.4f}  avg(leg2) = {a2:.4f}  avg(S) = {a1+4*a2:.4f}")
print(f"leg1 share of S: {a1/(a1+4*a2):.6f}")
print(f"B3 from leg1 only: {math.sqrt((math.pi**2/6)*a1):.6f}")
