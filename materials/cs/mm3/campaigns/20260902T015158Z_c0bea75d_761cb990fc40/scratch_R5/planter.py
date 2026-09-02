import sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/src')
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/scripts')
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/scratch')
from gatec_decomps import load_paper55
from gate_b_floor import prep
U, V, W = load_paper55()
c, _ = prep(W)
assert len(c) == 14, 'planted wrong d'
