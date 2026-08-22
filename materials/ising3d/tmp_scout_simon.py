from __future__ import annotations
import math,time
import numpy as np

Coord=tuple[int,int,int]

def shape_sites(kind:str,r:int):
    if kind=='l1':
        return {(x,y,z) for x in range(-r,r+1) for y in range(-r,r+1) for z in range(-r,r+1) if abs(x)+abs(y)+abs(z)<=r}
    if kind=='l2':
        return {(x,y,z) for x in range(-r,r+1) for y in range(-r,r+1) for z in range(-r,r+1) if x*x+y*y+z*z<=r*r}
    if kind=='linf':
        return {(x,y,z) for x in range(-r,r+1) for y in range(-r,r+1) for z in range(-r,r+1)}
    raise ValueError

def tm(sites:set[Coord],v:float):
    # choose propagation direction giving smallest max layer width
    choices=[]
    for a in range(3):
      layers={}
      for p in sites: layers.setdefault(p[a],[]).append(p)
      choices.append((max(map(len,layers.values())),a,layers))
    width,a,layers=min(choices)
    zs=sorted(layers)
    # normalize transverse coordinates as tuples of remaining axes, global union ordering per adjacent layer
    trans=lambda p: tuple(p[i] for i in range(3) if i!=a)
    data=[]
    dirs=[(1,0,0),(0,1,0),(0,0,1)]
    for z in zs:
      layer=sorted(layers[z],key=trans); idx={p:i for i,p in enumerate(layer)}
      masks=[]
      for p in layer:
       for d in dirs:
        q=tuple(p[i]+d[i] for i in range(3))
        if q in idx: masks.append((1<<idx[p])|(1<<idx[q]))
      boundary=[]
      for p in layer:
       deg=sum(tuple(p[i]+sgn*d[i] for i in range(3)) in sites for d in dirs for sgn in (-1,1))
       boundary.append(6-deg)
      data.append((layer,idx,masks,boundary))
    prev_layer=[]; plain=np.array([1.]); fixed=np.array([1.]); marked=np.array([0.])
    origin=(0,0,0)
    for li,(layer,idx,masks,boundary) in enumerate(data):
      n=len(layer); N=1<<n; incoming=np.zeros(N); infix=np.zeros(N); inmark=np.zeros(N)
      if li==0: incoming[0]=1; infix[0]=1
      else:
       # map prior outgoing bits by same transverse coordinate into current incoming mask
       pidx={trans(p):i for i,p in enumerate(prev_layer)}; cidx={trans(p):i for i,p in enumerate(layer)}
       common=set(pidx)&set(cidx)
       # every previous state having bits outside common invalid; remap bits
       for s in range(len(plain)):
        if any((s>>pidx[t])&1 for t in set(pidx)-common): continue
        cs=sum(((s>>pidx[t])&1)<<cidx[t] for t in common)
        incoming[cs]=plain[s]; infix[cs]=fixed[s]; inmark[cs]=marked[s]
      ids=np.arange(N,dtype=np.uint32); plain=incoming; fixed=infix; marked=inmark
      for mask in masks:
       sh=ids^mask; plain=plain+v*plain[sh]; fixed=fixed+v*fixed[sh]; marked=marked+v*marked[sh]
      target=(1<<idx[origin]) if origin in idx else 0
      if li==len(data)-1:
       P=plain[0]; Q=marked[target]+sum(m*fixed[target^(1<<j)] for j,m in enumerate(boundary)); return v*Q/P,width
      next_layer=data[li+1][0]; nidx={trans(p):i for i,p in enumerate(next_layer)}
      common=set(map(trans,layer))&set(map(trans,next_layer)); outgoing_positions=[idx[p] for p in layer if trans(p) in common]
      # states with outgoing on current common vertices, remap kept bits into CURRENT state convention; next loop remaps
      weights=np.zeros(N)
      for s in range(N):
       if any((s>>j)&1 for j in range(n) if j not in outgoing_positions): continue
       weights[s]=v**s.bit_count()
      sh=ids^target; nm=marked[sh].copy()
      for j,m in enumerate(boundary): nm+=m*fixed[sh^(1<<j)]
      marked=nm*weights; plain=plain*weights; fixed=fixed[sh]*weights
      prev_layer=layer
    raise AssertionError

v=math.tanh(0.21221190116616784)
for kind in ('l1','l2'):
 for r in range(1,6):
  s=shape_sites(kind,r); t=time.time()
  try: val,w=tm(s,v); print(kind,r,len(s),w,val,time.time()-t,flush=True)
  except MemoryError: print('MEM',kind,r,flush=True); break
