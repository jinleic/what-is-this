import sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/src')
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/scripts')
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/scratch')
from gatec_decomps import load_stapleton60
from gate_b_floor import prep
U, V, W = load_stapleton60()
c, _ = prep(U)
assert len(c) == 15, 'planted wrong d'
