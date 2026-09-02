import sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/src')
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56/scripts')
import sun56_verify as sv
from gate_b_floor import prep, subset_dfs
V = [tuple(r) for r in sv.V]
c, _ = prep(V)
assert len(c) == 12, 'planted wrong d'
