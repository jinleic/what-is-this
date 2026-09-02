"""K1 middle-term root-cause scratch: what was wrong in the assert chain.

K1 (corrected): the machine content is the pair
    rank(stacked C_i ⊗ A ⊗ A generator rows) == dim V == N - prod(s_i - t_i),
verified directly. The FIVE-way chain in the first draft included a middle
term r_q := prod(s_i - t_i) compared against dim V == 37 — wrong:
prod(s_i-t_i) = 27 is dim Q_0⊗Q_1⊗Q_2 (the QUOTIENT image), not dim V.
rank(pi0 ⊗ pi1 ⊗ pi2) = 27 (pure-tensor span of the images) — never 37, and
never asserted as == dim V in any correct statement of L2.
This scratch documents the root cause for the report; no re-derivation of the
theorem is implied (see theorem_corrected.md L2 for the corrected line).
"""
import sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/rs-pe3d')
from src.rs import Inst
from src.field import rank_mod

I = Inst(17, (4, 4, 4), (1, 1, 1))
rows = []
for i in range(3):
    rows.extend(I.lift_basis(i))
print('rank_gen =', rank_mod(rows, 17))
print('dim ker pred = N - prod(s-t) =', 64 - 27)
print('dim quotient image = prod(s-t) =', 27)
