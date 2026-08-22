
import numpy as np, math
from ising.lattices import chain, square, cubic, hyperrect
from ising.exact_enumeration import dos_bonds, even_subgraph_polynomial, joint_dos

# --- 1D chain, periodic: Z = 2^L (cosh^L K + sinh^L K) exactly ---
for L in [3,4,5,6,7,8]:
    lat = chain(L, True)
    g = dos_bonds(lat)
    P = even_subgraph_polynomial(lat)
    # even subgraphs of a cycle: empty set and full cycle -> P = 1 + v^L
    expect = [0]*(L+1); expect[0]=1; expect[L]=1
    assert P == expect, (L, P)
print("1D periodic chain even-subgraph polynomial 1+v^L : OK for L=3..8")

# --- 1D open chain: only the empty subgraph -> P = 1
for L in [3,5,8]:
    P = even_subgraph_polynomial(chain(L, False))
    assert P == [1]+[0]*(L-1), (L,P)
print("1D open chain P(v)=1 : OK")

# --- 2D 4x4 torus ---
lat = square(4,4,True)
P = even_subgraph_polynomial(lat)
print("2x2 square torus bonds:", square(2,2,True).n_bonds)
print("4x4 torus N,B =", lat.n_sites, lat.n_bonds)
print("P(v) for 4x4 torus:", P)
print("sum coefficients P(1) =", sum(P), " (= number of even subgraphs = 2^{B-N+1} = %d)" % (2**(lat.n_bonds-lat.n_sites+1)))

# --- 2x2x2 cube torus ---
c = cubic(2,2,2,True)
print("2x2x2 torus:", c.describe(), "degrees", set(c.degree_sequence()))
g = dos_bonds(c)
print("DOS:", {b:int(v) for b,v in enumerate(g) if v})
Pc = even_subgraph_polynomial(c)
print("P(v) 2x2x2:", Pc, "P(1)=",sum(Pc), "expect 2^(B-N+1)=",2**(c.n_bonds-c.n_sites+1))
