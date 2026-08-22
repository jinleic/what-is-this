from __future__ import annotations
from fractions import Fraction
import math,time
from collections import defaultdict
Coord=tuple[int,int,int]

def sites_l1(r): return {(x,y,z) for x in range(-r,r+1) for y in range(-r,r+1) for z in range(-r,r+1) if abs(x)+abs(y)+abs(z)<=r}
def sites_linf(r): return {(x,y,z) for x in range(-r,r+1) for y in range(-r,r+1) for z in range(-r,r+1)}
def sites_l2(r): return {(x,y,z) for x in range(-r,r+1) for y in range(-r,r+1) for z in range(-r,r+1) if x*x+y*y+z*z<=r*r}
D=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
def graph(sites):
 vs=list(sites); idx={p:i for i,p in enumerate(vs)}; es=[]
 for i,p in enumerate(vs):
  for d in D[::2]:
   q=tuple(p[j]+d[j] for j in range(3))
   if q in idx: es.append((i,idx[q]))
 return vs,es

def minfill_order(n,edges,keep=None):
 adj=[set() for _ in range(n)]
 for a,b in edges: adj[a].add(b);adj[b].add(a)
 alive=set(range(n)); order=[]; maxw=0; keep=set() if keep is None else set(keep)
 while alive-keep:
  best=min(alive-keep,key=lambda v:(sum(1 for a in adj[v] for b in adj[v] if a<b and b not in adj[a]),len(adj[v]),v))
  nei=adj[best]&alive; maxw=max(maxw,len(nei))
  for a in nei:
   adj[a].update(nei-{a})
   adj[a].discard(best)
  alive.remove(best);order.append(best)
 return order,maxw

def eliminate(sites,v):
 vs,edges=graph(sites); n=len(vs); order,w=minfill_order(n,edges)
 factors=[]
 # tuple(vars), values index little endian
 for a,b in edges: factors.append(((a,b),(1.0,v,v,1.0)))
 q=[6-sum(tuple(p[j]+d[j] for j in range(3)) in sites for d in D) for p in vs]
 oi=vs.index((0,0,0))
 # expectation separately for every boundary vertex (slow), use a marked numerator via dual-number factor values (Z, sum q sigma_o sigma_x)
 # represent values pairs, initial edge factors scalar; insert global observable is sum of unary/pair terms, handle by derivative marker accumulated via semiring product
 # factors values (z,m), product (z1z2,z1m2+m1z2); add component observables as factors? exp(epsilon O)=1+eps O, unary pair between o,x.
 # One first-order marker factor: (1 + eps * sum_x q_x sigma_o sigma_x).
 # This factor is dense on all boundary variables, so instead compute each
 # two-point numerator separately below; the prototype is only for small widths.
 def partition(extra):
  ff=[]
  for var,val in factors: ff.append((var,tuple((x,0.0) for x in val)))
  ff.extend(extra)
  return ff
 for x in order:
  hit=[f for f in ff if x in f[0]]; ff=[f for f in ff if x not in f[0]]
  union=sorted(set().union(*(set(f[0]) for f in hit))); rest=[u for u in union if u!=x]; pos={u:i for i,u in enumerate(union)}
  out=[]
  for mask in range(1<<len(rest)):
   vals=[]
   for xb in (0,1):
    assign={u:(mask>>j)&1 for j,u in enumerate(rest)};assign[x]=xb; z,m=1.,0.
    for vars_,arr in hit:
     j=sum(assign[u]<<k for k,u in enumerate(vars_)); zz,mm=arr[j]; z,m=z*zz,z*mm+m*zz
    vals.append((z,m))
   out.append((vals[0][0]+vals[1][0],vals[0][1]+vals[1][1]))
  ff.append((tuple(rest),tuple(out)))
 # remaining constants multiply
 z,m=1.,0.
 for _,arr in ff: zz,mm=arr[0];z,m=z*zz,z*mm+m*zz
 return v*m/z,w,n,len(edges)

v=math.tanh(0.21221190116616784)
for kind,fn in [('l1',sites_l1),('l2',sites_l2),('box',sites_linf)]:
 for r in range(1,7):
  s=fn(r); t=time.time()
  if len(s)>1000: break
  try: val,w,n,e=eliminate(s,v); print(kind,r,n,e,w,val,time.time()-t,flush=True)
  except MemoryError: print('MEM',kind,r,flush=True);break
