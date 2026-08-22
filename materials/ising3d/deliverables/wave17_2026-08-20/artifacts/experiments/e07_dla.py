
import time, json, numpy as np
from ising.clifford import (tfim_generators, dla_pauli_closure, frustration_graph,
                            find_induced_claw, symplectic_form)

def grid(a,b,per_a=False,per_b=False):
    n=a*b; sites=list(range(n)); bonds=[]
    idx=lambda x,y: x*b+y
    for x in range(a):
        for y in range(b):
            if x+1<a: bonds.append((idx(x,y),idx(x+1,y)))
            elif per_a and a>2: bonds.append((idx(x,y),idx(0,y)))
            if y+1<b: bonds.append((idx(x,y),idx(x,y+1)))
            elif per_b and b>2: bonds.append((idx(x,y),idx(x,0)))
    return sites,bonds,n

rows=[]
print("=== 1D chains (2D classical Ising layer operator) ===")
for L in range(2,13):
    sites,bonds,n = grid(1,L)
    g,lab,n = tfim_generators(sites,bonds,n)
    t0=time.time(); S,sat = dla_pauli_closure(g,n); dt=time.time()-t0
    adj = frustration_graph(g,n); claw = find_induced_claw(adj)
    print(f"  open chain n={L:2d}: |gens|={len(g):3d}  dim DLA={len(S):6d}   n(2n-1)={n*(2*n-1):6d}  claw={'YES' if claw else 'no'}  ({dt:.2f}s)")
    rows.append(dict(kind="chain_open",n=L,gens=len(g),dla=len(S),so2n=n*(2*n-1),claw=bool(claw)))

print("=== 1D rings (periodic chain) ===")
for L in range(3,11):
    sites,bonds,n = grid(1,L,per_b=True)
    g,lab,n = tfim_generators(sites,bonds,n)
    S,sat = dla_pauli_closure(g,n)
    adj = frustration_graph(g,n); claw = find_induced_claw(adj)
    print(f"  ring n={L:2d}: dim DLA={len(S):6d}  2*n(2n-1)={2*n*(2*n-1):6d}  claw={'YES' if claw else 'no'}")
    rows.append(dict(kind="chain_periodic",n=L,gens=len(g),dla=len(S),so2n=n*(2*n-1),claw=bool(claw)))

print("=== 2D layers (3D classical Ising layer operator) ===")
for (a,b) in [(2,2),(2,3),(2,4),(3,3)]:
    sites,bonds,n = grid(a,b)
    g,lab,n = tfim_generators(sites,bonds,n)
    t0=time.time()
    try:
        S,sat = dla_pauli_closure(g,n,max_size=1<<23)
        dt=time.time()-t0
        adj = frustration_graph(g,n); claw = find_induced_claw(adj)
        maxposs = 2**(2*n-1)-1
        print(f"  grid {a}x{b} (n={n:2d}): dim DLA={len(S):8d}  max possible 2^(2n-1)-1={maxposs:9d}  so(2n) dim={n*(2*n-1):5d}  claw={'YES' if claw else 'no'} ({dt:.1f}s)")
        rows.append(dict(kind="grid_open",a=a,b=b,n=n,gens=len(g),dla=len(S),so2n=n*(2*n-1),span_max=maxposs,claw=bool(claw)))
    except Exception as e:
        print("  grid",a,b,"FAILED",e)
json.dump(rows, open("results/dla_table.json","w"), indent=1)
