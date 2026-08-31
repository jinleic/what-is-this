"""Validate Φ1 = ρ3·G10 + 3ρ1ρ2·G20 + ρ1³·G30 − t·G01 against the recursion
Φ1 = D_t[Φ0] − t·∂x∂yΦ0 with Φ0 = C_{ρ(t)}, by complex-FD of Φ0 in t (∂r-direction via ρ) and
dxdy-machinery for the ∂x∂y part. Cross-checks the paper's c-table END-TO-END."""
import sys, math, cmath
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/kg/src")
from d3h_phi3 import G_pa_eval, rho_j, C_TABLE, _IO_LIST
from d3h_kernels import vartheta
from d3h_dxdy import dxdy_eval
from flint import acb
from core import Prec

# C_ρ(x,y) at (x,y,t): use dxdy_eval? Need C ITSELF (a=0 base): use the series C_r = (π/2)Σ r^b q_bq_b
# — for |r|<1 converges; at t=0.3: ρ(0.3)=0.27 fine. Use float series... do arb-lite:
import math
def psi_n(n, xf):
    hp, hc = 1.0, xf
    if n == 0: return 1.0
    if n == 1: return xf
    for k in range(1, n):
        hp, hc = hc, xf*hc - k*hp
    return hc/math.sqrt(math.factorial(n))
Vv = 1+0.34101124**2+0.05276111**2
vt0 = 0.136419125/math.sqrt(Vv)
def q_ser(b, u):
    if b == 0: return math.erf(u/math.sqrt(2))
    return 2*math.exp(-u*u/2)/math.sqrt(2*math.pi)*psi_n(b-1, -u)/math.sqrt(b)
def C_series(xf, yf, tf, NB=60):
    r = (tf - (0.34101124**2)*tf**3 + (0.05276111**2)*tf**5)/Vv
    ux = vt0*(xf**3-3*xf)/math.sqrt(6); uy = vt0*(yf**3-3*yf)/math.sqrt(6)
    s = 0.0
    for b in range(NB):
        s += (r**b)*q_ser(b, ux)*q_ser(b, uy)
    return (math.pi/2)*s

def dCdr(xf, yf, tf, h=1e-5):
    # ∂r C at FIXED (u,v)... but FD in t moves r through ρ: ∂rC = ∂tC/ρ1(t):
    return (C_series(xf, yf, tf+h) - C_series(xf, yf, tf-h))/(2*h)/((1 - 3*(0.34101124**2)*tf*tf + 5*(0.05276111**2)*tf**4)/Vv)

def dxdyC_series(xf, yf, tf, h=1e-3):
    f = lambda xx, yy: C_series(xx, yy, tf)
    return (f(xf+h, yf+h) - f(xf+h, yf-h) - f(xf-h, yf+h) + f(xf-h, yf-h))/(4*h*h)

xf, yf, tf = 0.4, -0.7, 0.3
# LHS recursion: Φ1 = DΦ0 − t∂xyΦ0 with DΦ0 = ρ1·∂rC (chain: D[C_ρ(t)] = ρ1·∂rC):
rho1 = (1 - 3*(0.34101124**2)*tf*tf + 5*(0.05276111**2)*tf**4)/Vv
Phi1_rec = rho1*dCdr(xf, yf, tf) - tf*dxdyC_series(xf, yf, tf)

# RHS c-table: Φ1 = ρ3·G10 + 3ρ1ρ2·G20 + ρ1³·G30 − t·G01:
rho2 = (6*(0.34101124**2)*tf*tf + 20*(0.05276111**2)*tf**4)/Vv
rho3 = (12*(0.34101124**2)*tf*tf + 80*(0.05276111**2)*tf**4)/Vv
with Prec(150):
    vt = vartheta()
    x, y, t = acb(xf), acb(yf), acb(tf)
    rhs = (rho3*float(G_pa_eval(1,0,x,y,t,vt,150).mid().real)
           + 3*rho1*rho2*float(G_pa_eval(2,0,x,y,t,vt,150).mid().real)
           + rho1**3*float(G_pa_eval(3,0,x,y,t,vt,150).mid().real)
           - tf*float(G_pa_eval(0,1,x,y,t,vt,150).mid().real))
print(f"Φ1 recursion (FD): {Phi1_rec:.8g}")
print(f"Φ1 c-table ( closed): {rhs:.8g}")
print(f"ratio: {rhs/Phi1_rec:.6f}")
