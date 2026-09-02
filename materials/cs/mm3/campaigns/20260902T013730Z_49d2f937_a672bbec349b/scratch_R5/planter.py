import sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/src')
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/mm3/scratch')
from pathlib import Path
lines = Path('/Users/jinleic/jinleic-workspace/cs/mm3/scratch/mws59_layout.txt').read_text().splitlines()
data = lines[236:265]
blocks, cur = [], []
for L in data:
    s = L.strip()
    if s.startswith('#'):
        blocks.append(cur); cur = []
    elif s and not s.startswith(('A','T')):
        try: cur.append([int(x) for x in s.split()])
        except ValueError: pass
if cur: blocks.append(cur)
b1,b2,b3 = blocks
V = [tuple(b2[k][r] for k in range(9)) for r in range(23)]
from gate_b_floor import prep
c, _ = prep(V)
assert len(c) == 13, 'planted wrong d'
