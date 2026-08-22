import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import sys
import numpy as np
from ldpc import BpOsdDecoder
sys.path.insert(0, 'experiments')
import exp029_circuit_distance as e
from qec_research.codes.bicycle import PBBSpec, build_pbb
from qec_research.symplectic.core import lambda_swap, symplectic_weight, symplectic_product_matrix

row=e.load_catalogue_row(); code=build_pbb(PBBSpec(row['ell'],row['m'],row['A_terms'],row['B_terms'],row['C_terms'],row['D_terms']))
# Syndrome equations for error e: H Lambda e^T; logical action similarly.
checks = lambda_swap(code.H)
logicals = lambda_swap(code.logical_basis())
effective = np.vstack([checks, logicals]).astype(np.uint8)
ns = checks.shape[0]; nl = logicals.shape[0]
print('matrix',effective.shape,'rank?',flush=True)
dec=BpOsdDecoder(effective,error_rate=.05,bp_method='product_sum',osd_method='osd_0',osd_order=0,max_iter=100)
rng=np.random.default_rng(20260812)
best=None; bw=999
for trial in range(200000):
    bits=np.zeros(nl,dtype=np.uint8)
    while not bits.any(): bits=rng.integers(0,2,nl,dtype=np.uint8)
    syndrome=np.concatenate([np.zeros(ns,dtype=np.uint8),bits])
    v=np.asarray(dec.decode(syndrome),dtype=np.uint8)
    if np.array_equal((effective@v)%2,syndrome):
        w=symplectic_weight(v)
        if w < bw:
            bw=w;best=v.copy();print('trial',trial,'best',bw,flush=True)
            if bw <= 12: break
if best is None: raise SystemExit('no witness')
print('weight',bw)
print('vector',''.join(map(str,best.tolist())))
print('x',np.flatnonzero(best[:code.n]).tolist())
print('z',np.flatnonzero(best[code.n:]).tolist())
print('valid syndrome',not symplectic_product_matrix(best[None,:],code.H).any())
print('logical',code.is_logical(best))
