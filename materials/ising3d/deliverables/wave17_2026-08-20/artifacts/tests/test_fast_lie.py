
import time, json
import numpy as np
from ising.clifford import pauli_x, zz
from ising.clifford.fast_lie import GeneratedAlgebra, P1, P2

def grid(a,b,per_a=False,per_b=False):
    n=a*b; idx=lambda x,y: x*b+y; bonds=[]
    for x in range(a):
        for y in range(b):
            if x+1<a: bonds.append((idx(x,y),idx(x+1,y)))
            elif per_a and a>2: bonds.append((idx(x,y),idx(0,y)))
            if y+1<b: bonds.append((idx(x,y),idx(x,y+1)))
            elif per_b and b>2: bonds.append((idx(x,y),idx(x,0)))
    return list(range(n)), bonds, n

def onsager_dim(a,b,per_a=False,per_b=False,p=P1,max_dim=None,report=None):
    sites,bonds,n = grid(a,b,per_a,per_b)
    xs=[pauli_x(n,i) for i in sites]; zs=[zz(n,i,j) for i,j in bonds]
    singles = xs+zs
    A={g:1 for g in xs}; B={g:1 for g in zs}
    alg=GeneratedAlgebra([A,B],singles,n)
    alg.p=p
    gA=list(range(len(xs))); gB=list(range(len(xs),len(singles)))
    d=alg.closure_dim([A,B],[gA,gB],max_dim=max_dim,report=report)
    return d, alg.m

# regression against the slow implementation
for L,exp in [(2,4),(3,9),(4,16),(5,25),(6,36),(7,49),(8,64)]:
    d,m = onsager_dim(1,L)
    assert d==exp, (L,d,exp)
print("chain n^2 law reproduced by fast engine for n=2..8  OK")
for L,exp in [(3,8),(4,11),(5,14),(6,17),(7,20)]:
    d,m = onsager_dim(1,L,per_b=True); assert d==exp,(L,d,exp)
print("ring 3n-1 law reproduced  OK")
d,m = onsager_dim(2,3); print("2x3 grid:", d, "(slow engine gave 263)  support size", m)
assert d==263
print("REGRESSION OK")
