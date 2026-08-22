from __future__ import annotations
import numpy as np, math,time

def initial_tensors(nx,ny,nz,v):
 # Tensor network high-temp: each edge bit weight sqrt(v) at endpoints; vertex parity delta.
 # Contract along z exactly into a 2D PEPS layer? Boundary-MPS contraction across x/y/z with SVD truncation is complex.
 pass
