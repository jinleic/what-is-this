import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import sys,time
import numpy as np
from ldpc import BpOsdDecoder
sys.path.insert(0,'experiments')
sys.path.insert(0,'third_party/qcode-discovery')
import exp029_circuit_distance as e
from qec_research.codes.bicycle import PBBSpec,build_pbb
from qec_research.symplectic.core import symplectic_weight,symplectic_product_matrix
from qec_research.gf2.linalg import nullspace_np, rref_np

def _compute_achievable_basis(checks, logicals):
    kernel = nullspace_np(checks)
    image = (logicals @ kernel.T) & 1
    rows, _ = rref_np(image.T)
    rows = rows[rows.any(axis=1)]
    return None if len(rows) >= len(logicals) else rows

def _safe_for_bposd(matrix):
    return bool(matrix.size and matrix.any(axis=0).all())
row=e.load_catalogue_row();code=build_pbb(PBBSpec(row['ell'],row['m'],row['A_terms'],row['B_terms'],row['C_terms'],row['D_terms']))
n=code.n; S=code.H; L=code.logical_basis();Sx,Sz=S[:,:n],S[:,n:];Lx,Lz=L[:,:n],L[:,n:]
rng=np.random.default_rng(20260812); best=None; bw=999;t0=time.time()
channels=[]
channels.append(('X',Sz,Lz,np.zeros(n,dtype=np.uint8)))
channels.append(('Z',Sx,Lx,np.ones(n,dtype=np.uint8)))
channels.append(('Y',(Sx^Sz),(Lx^Lz),'Y'))
for i in range(1000):
 m=rng.integers(0,2,n,dtype=np.uint8)
 channels.append((f'M{i}',(Sz*(1-m)^Sx*m).astype(np.uint8),(Lz*(1-m)^Lx*m).astype(np.uint8),m))
for ci,(name,C,A,mask) in enumerate(channels):
 active=A[A.any(axis=1)]; checks=C[C.any(axis=1)];
 if not len(active):continue
 E=np.vstack([checks,active]).astype(np.uint8)
 if not _safe_for_bposd(E):continue
 basis=_compute_achievable_basis(checks,active)
 dec=BpOsdDecoder(E,error_rate=.05,bp_method='product_sum',osd_method='osd_0',osd_order=0,max_iter=100)
 synd=np.zeros(len(checks)+len(active),dtype=np.uint8)
 for trial in range(300):
  if basis is None:
   bits=np.zeros(len(active),dtype=np.uint8)
   while not bits.any():bits=rng.integers(0,2,len(active),dtype=np.uint8)
  else:
   co=np.zeros(len(basis),dtype=np.uint8)
   while not co.any():co=rng.integers(0,2,len(basis),dtype=np.uint8)
   bits=np.bitwise_xor.reduce(basis[co.astype(bool)],axis=0)
  synd[:len(checks)]=0;synd[len(checks):]=bits
  r=np.asarray(dec.decode(synd),dtype=np.uint8)
  if np.array_equal((E@r)%2,synd) and r.any():
   v=np.zeros(2*n,dtype=np.uint8)
   if isinstance(mask,str):v[:n]=r;v[n:]=r
   else:v[:n]=r*(1-mask);v[n:]=r*mask
   w=symplectic_weight(v)
   if w<bw:
    bw=w;best=v.copy();print('channel',name,'trial',trial,'best',bw,'wall',time.time()-t0,flush=True)
    if bw<=12:break
 if bw<=12:break
 if ci%25==0:print('progress',ci,'best',bw,'wall',time.time()-t0,flush=True)
print('FINAL',bw,'wall',time.time()-t0)
if best is not None:
 print('vector',''.join(map(str,best.tolist())))
 print('x',np.flatnonzero(best[:n]).tolist());print('z',np.flatnonzero(best[n:]).tolist())
 print('commutes',not symplectic_product_matrix(best[None,:],code.H).any(),'logical',code.is_logical(best))
