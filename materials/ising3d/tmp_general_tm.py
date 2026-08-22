from __future__ import annotations
import numpy as np, math,time
D=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
def build(sites:set[tuple[int,int,int]],axis=2):
 zs=sorted({p[axis] for p in sites}); trans_axes=[i for i in range(3) if i!=axis]
 trans=lambda p:(p[trans_axes[0]],p[trans_axes[1]])
 layers=[]
 for z in zs:
  ps=sorted((p for p in sites if p[axis]==z),key=trans); idx={p:i for i,p in enumerate(ps)}; masks=[]
  for i,p in enumerate(ps):
   for d in ((1,0,0),(0,1,0),(0,0,1)):
    q=tuple(p[j]+d[j] for j in range(3))
    if q in idx: masks.append((1<<i)|(1<<idx[q]))
  b=[]
  for p in ps:
   b.append(6-sum(tuple(p[j]+d[j] for j in range(3)) in sites for d in D))
  layers.append((ps,idx,masks,b))
 return layers,trans

def eval_tm(sites,v,axis=2):
 layers,trans=build(sites,axis); prev=None; plain=fixed=marked=None
 for li,(ps,idx,masks,boundary) in enumerate(layers):
  n=len(ps); N=1<<n; ids=np.arange(N,dtype=np.uint32)
  if li==0:
   plain=np.zeros(N);fixed=np.zeros(N);marked=np.zeros(N);plain[0]=fixed[0]=1
  else:
   # previous arrays indexed by prev layer; since nearest-neighbor slice sets may differ, restrict states to common coordinates and remap
   pidx={trans(p):i for i,p in enumerate(prev)}; cidx={trans(p):i for i,p in enumerate(ps)}; common=sorted(set(pidx)&set(cidx)); outside=set(pidx)-set(common)
   incoming=np.zeros(N);infix=np.zeros(N);inmark=np.zeros(N)
   valid=np.ones(len(plain),dtype=bool)
   for t in outside: valid &= ((np.arange(len(plain),dtype=np.uint32)>>pidx[t])&1)==0
   old=np.flatnonzero(valid).astype(np.uint32); new=np.zeros(len(old),dtype=np.uint32)
   for t in common: new |= ((old>>pidx[t])&1).astype(np.uint32)<<cidx[t]
   incoming[new]=plain[old];infix[new]=fixed[old];inmark[new]=marked[old]
   plain,fixed,marked=incoming,infix,inmark
  for mask in masks:
   sh=ids^mask;plain=plain+v*plain[sh];fixed=fixed+v*fixed[sh];marked=marked+v*marked[sh]
  target=1<<idx[(0,0,0)] if (0,0,0) in idx else 0
  if li==len(layers)-1:
   return v*(marked[target]+sum(m*fixed[target^(1<<j)] for j,m in enumerate(boundary)))/plain[0],max(len(x[0]) for x in layers)
  nxt=layers[li+1][0]; ncoords={trans(p) for p in nxt}; keep=[j for j,p in enumerate(ps) if trans(p) in ncoords]
  # outgoing edge mask uses current ordering, and only common sites can carry edges
  valid=np.ones(N,dtype=bool)
  for j in set(range(n))-set(keep): valid &= ((ids>>j)&1)==0
  weights=np.zeros(N);weights[valid]=np.power(v,np.bitwise_count(ids[valid]))
  sh=ids^target; nm=marked[sh].copy()
  for j,m in enumerate(boundary): nm+=m*fixed[sh^(1<<j)]
  marked=nm*weights;plain=plain*weights;fixed=fixed[sh]*weights;prev=ps
 raise AssertionError

def geom(kind,r):
 return {(x,y,z) for x in range(-r,r+1) for y in range(-r,r+1) for z in range(-r,r+1) if (kind=='l1' and abs(x)+abs(y)+abs(z)<=r) or (kind=='l2' and x*x+y*y+z*z<=r*r) or (kind=='linf')}

if __name__=='__main__':
 v=math.tanh(.21221190116616784)
 for k,r in [('l1',1),('l1',2),('l1',3),('l2',1),('l2',2),('l2',3)]:
  s=geom(k,r);t=time.time();print(k,r,len(s),eval_tm(s,v),time.time()-t,flush=True)
