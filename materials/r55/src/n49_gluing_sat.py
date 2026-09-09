
import os
import sys, json, time, os
from pysat.solvers import Cadical195
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

X1, X2 = [4094, 10363617, 13533905, 12388753, 10007597, 10676499, 14238983, 15877135, 14783081, 4667655, 3954897, 6669865, 10938128, 7950898, 12795104, 5708084, 1764168, 2740902, 10001934, 2565214, 4564186, 4865448, 3206084, 283134], [4094, 13314137, 11847737, 13161223, 11613319, 15369413, 16057123, 14426929, 12207305, 5183689, 3578161, 517639, 8906264, 7133330, 7551244, 5848108, 3725906, 13131554, 11550404, 2335658, 4572628, 877940, 1237738, 397822]
def bits(x):
    while x:
        b = x & -x; yield b.bit_length()-1; x ^= b
def compl(n, adj):
    full=(1<<n)-1; return [full & ~a & ~(1<<v) for v,a in enumerate(adj)]
def common(adj,u,v): return bin(adj[u]&adj[v]).count("1")
def k4s(n, adj):
    out=[]
    for u in range(n):
        for du in bits(adj[u]>>(u+1)):
            v=u+1+du; B=adj[u]&adj[v]&~((1<<(v+1))-1)
            for w in bits(B):
                C=B&adj[w]&~((1<<(w+1))-1)
                for x in bits(C): out.append((u,v,w,x))
    return out

def solve(Xa, tag):
    n=24; Db=compl(n,Xa)
    pool=IDPool(); M=[[pool.id(("m",w,x)) for x in range(n)] for w in range(n)]
    cnf=[]
    for w in range(n):
        cnf+=CardEnc.equals([M[w][x] for x in range(n)],12,vpool=pool,encoding=EncType.seqcounter).clauses
    for x in range(n):
        cnf+=CardEnc.equals([M[w][x] for w in range(n)],12,vpool=pool,encoding=EncType.seqcounter).clauses
    def inner(vecs,adjg,k_adj,k_non,tagg):
        cls=[]
        for u in range(n):
            for v in range(u+1,n):
                t=(k_adj if (adjg[u]>>v)&1 else k_non)-common(adjg,u,v)
                ys=[]
                for x in range(n):
                    y=pool.id(("y",tagg,u,v,x)); a,b=vecs[u][x],vecs[v][x]
                    cls+=[[-y,a],[-y,b],[y,-a,-b]]; ys.append(y)
                cls+=CardEnc.equals(ys,t,vpool=pool,encoding=EncType.seqcounter).clauses
        return cls
    cnf+=inner(M,Xa,10,11,"r")
    cnf+=inner([[M[w][x] for w in range(n)] for x in range(n)],Db,11,12,"c")
    # row picks induce exactly 22+l_X(w) edges in D; column picks 22+l_D(x)... derive:
    ledge=lambda adj,w:[e for e in ((a,b) for a in bits(adj[w]) for b in bits(adj[w]) if a<b) if (adj[e[0]]>>e[1])&1]
    Dedges=[(x,y) for x in range(n) for y in range(x+1,n) if (Db[x]>>y)&1]
    Xedges=[(u,v) for u in range(n) for v in range(u+1,n) if (Xa[u]>>v)&1]
    for w in range(n):
        lw=len(ledge(Xa,w)); ys=[]
        for (x,y) in Dedges:
            z=pool.id(("re",w,x,y)); a,b=M[w][x],M[w][y]
            cnf+=[[-z,a],[-z,b],[z,-a,-b]]; ys.append(z)
        cnf+=CardEnc.equals(ys,22+lw,vpool=pool,encoding=EncType.seqcounter).clauses
    for x in range(n):
        lx=len(ledge(Db,x)); ys=[]
        for (u,v) in Xedges:
            z=pool.id(("ce",x,u,v)); a,b=M[u][x],M[v][x]
            cnf+=[[-z,a],[-z,b],[z,-a,-b]]; ys.append(z)
        # e(X[col_x]) = 132 - e_local_D(x) - sum_{y~x}(11-c_D) ; compute exact:
        s=sum(11-common(Db,x,y) for y in bits(Db[x]))
        cnf+=CardEnc.equals(ys,132-lx-s,vpool=pool,encoding=EncType.seqcounter).clauses
    # K5/I5 exclusions across the cut
    cX=compl(n,Xa); cD=compl(n,Db)
    triX=[(u,v,w) for u in range(n) for v in range(u+1,n) for w in range(v+1,n)
          if (Xa[u]>>v)&1 and (Xa[u]>>w)&1 and (Xa[v]>>w)&1]
    triD=[(x,y,z) for x in range(n) for y in range(x+1,n) for z in range(y+1,n)
          if (Db[x]>>y)&1 and (Db[x]>>z)&1 and (Db[y]>>z)&1]
    for (u,v,w) in triX:
        for (x,y) in Dedges:
            cnf.append([-M[u][x],-M[u][y],-M[v][x],-M[v][y],-M[w][x],-M[w][y]])
    for (u,v) in Xedges:
        for (x,y,z) in triD:
            cnf.append([-M[u][x],-M[u][y],-M[u][z],-M[v][x],-M[v][y],-M[v][z]])
    for u in range(n):
        for q in k4s(n,Db): cnf.append([-M[u][q[0]],-M[u][q[1]],-M[u][q[2]],-M[u][q[3]]])
    for q in k4s(n,cX):
        for x in range(n): cnf.append([M[q[0]][x],M[q[1]][x],M[q[2]][x],M[q[3]][x]])
    i3X=[(u,v,w) for u in range(n) for v in range(u+1,n) for w in range(v+1,n)
         if (cX[u]>>v)&1 and (cX[u]>>w)&1 and (cX[v]>>w)&1]
    i3D=[(x,y,z) for x in range(n) for y in range(x+1,n) for z in range(y+1,n)
         if (cD[x]>>y)&1 and (cD[x]>>z)&1 and (cD[y]>>z)&1]
    nD=[(x,y) for x in range(n) for y in range(x+1,n) if (cD[x]>>y)&1]
    nX=[(u,v) for u in range(n) for v in range(u+1,n) if (cX[u]>>v)&1]
    for (u,v,w) in i3X:
        for (x,y) in nD:
            cnf.append([M[u][x],M[u][y],M[v][x],M[v][y],M[w][x],M[w][y]])
    for (u,v) in nX:
        for (x,y,z) in i3D:
            cnf.append([M[u][x],M[u][y],M[u][z],M[v][x],M[v][y],M[v][z]])
    t0=time.time()
    with Cadical195(bootstrap_with=cnf) as s:
        res=s.solve()
        model = s.get_model() if res else None
    out={"combo":tag,"clauses":len(cnf),"result":"SAT" if res else "UNSAT","secs":round(time.time()-t0,1)}
    print(json.dumps(out), flush=True)
    return out

results=[solve(X1,"X1/complX1"), solve(X2,"X2/complX2")]
json.dump(results, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'n49_gluing_sat.json'),"w"), indent=1)
print("DONE", flush=True)
