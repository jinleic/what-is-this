#!/usr/bin/env python3
"""STANDALONE machine-checkable certificate for the n=18 Q=3 Theorem-M legs
(pure stdlib; exact integer arithmetic).

Solver functions copied VERBATIM from the audited
scratch/kobon/n14/escape/appendix_independent_check.py (mult_solutions /
feasible). System (scratch/kobon/n18/obstruction_n18.md): under
near-saturation the multiplicities m_i satisfy
    m_{i-1} + m_i = 2 s_i,     0 <= m_i <= s_i * s_{i+1},    m_i integer,
indices mod R on the direction cycle with class-size list (s_0..s_{R-1}).

Legs (clearing obstruction_n18.md §7 pending queue):
  (A') n=18 paired: three size-2 + twelve singles, R=15 (odd cycle):
       all C(15,3)=455 placements infeasible; gap classification subtotals
       140/90/225 (all-even / two-odd e=0 / two-odd e>0); 19 dihedral orbit
       representatives; reflection-only quotient count 231.
  (B') n=18 triad: one size-3 + fifteen singles, R=16 (even cycle):
       alternating-sum inconsistency at every t (|a|=4); all 16 placements
       infeasible.
  (C') controls that MUST stay feasible: all-single odd/even cycles and the
       hexagon control [2,2,2] (six lines in three pairs).
"""
from itertools import combinations


def mult_solutions(sizes):
    """Exact solution SET of m_{i-1}+m_i = 2 s_i on the R-cycle as
    m_i(t) = a_i * t + b_i for integer t; returns (a, b, fixed_t or None).
    For even R the closure forces either 'inconsistent' (fixed_t=() ) or a
    free parameter t."""
    R = len(sizes)
    d = [2 * s for s in sizes]
    a, b = [1], [0]
    for i in range(1, R):
        a.append(-a[-1])
        b.append(d[i] - b[-1])
    if R % 2 == 1:                       # a_{R-1} = +1: closure: 2 t = d_0 - b
        num, den = d[0] - b[-1], 2
        if num % den:
            return a, b, ()              # no integer t
        t0 = num // den
        return a, b, t0
    else:                                # a_{R-1} = -1: closure: b_{R-1} == d_0
        if b[-1] != d[0]:
            return a, b, ()
        return a, b, None                # free integer parameter t


def feasible(sizes, verbose=False):
    R = len(sizes)
    U = [sizes[i] * sizes[(i + 1) % R] for i in range(R)]
    a, b, tfix = mult_solutions(sizes)
    if tfix == ():
        return False, 'no rational solution (even-cycle alternating sum != 0)'
    if tfix is not None:
        cand = [a[i] * tfix + b[i] for i in range(R)]
        return (all(0 <= cand[i] <= U[i] for i in range(R)),
                f'unique candidate {cand}, bounds {U}')
    # free integer t: intersect per-edge intervals
    lo, hi = -10**9, 10**9
    for i in range(R):
        if a[i] == 1:
            lo, hi = max(lo, -b[i]), min(hi, U[i] - b[i])
        else:
            lo, hi = max(lo, b[i] - U[i]), min(hi, b[i])
    for t in range(lo, hi + 1):
        cand = [a[i] * t + b[i] for i in range(R)]
        if all(0 <= cand[i] <= U[i] for i in range(R)):
            return True, f'free-parameter solution at t={t}: {cand}'
    return False, f'no integer t in [{lo},{hi}]'


def gap_triple(S, R):
    """Sorted cyclic gap composition: gaps count INTERMEDIATE single classes
    between consecutive doubled positions, so the three gaps sum to R-3
    (12 for R=15), matching obstruction_n18.md's convention r1+r2+r3=12."""
    p = sorted(S)
    gaps = [((p[(i + 1) % 3] - p[i]) % R) - 1 for i in range(3)]
    return tuple(sorted(gaps))


print('== (Cprime) controls: must stay feasible ==')
for sizes, tag in [([1] * 9, 'Blanc-21 at n=9'), ([1] * 15, 'all-single R=15'),
                   ([1] * 17, 'all-single R=17'), ([1] * 5, 'pentagon'),
                   ([2, 2, 2], 'hexagon (n=6, Q=3, three 2-classes)'),
                   ([1, 1, 1, 1], 'even R=4 all-single (m=1s)')]:
    ok, why = feasible(sizes)
    assert ok, f'control {tag} wrongly infeasible: {why}'
    print(f'  {tag}: feasible [{why[:60]}]')

print('== (Aprime) n=18 paired: three size-2 + twelve singles, R=15 ==')
cnt = 0
buck = {'all_even': 0, 'two_odd_e0': 0, 'two_odd_ep': 0}
orbits = set()
for S in combinations(range(15), 3):
    sizes = [2 if i in S else 1 for i in range(15)]
    ok, why = feasible(sizes)
    assert not ok, f'paired placement {S} FEASIBLE: {why}'
    cnt += 1
    g = gap_triple(S, 15)
    orbits.add(g)
    par = [x % 2 for x in g]
    if par == [0, 0, 0]:
        buck['all_even'] += 1
    elif sum(par) == 2:
        buck['two_odd_e0' if 0 in g else 'two_odd_ep'] += 1
    else:
        raise AssertionError(f'unclassified gap parity {g}')
print(f'  {cnt}/455 paired placements infeasible [odd cycle]')
assert buck == {'all_even': 140, 'two_odd_e0': 90, 'two_odd_ep': 225}, buck
print(f'  gap subtotals 140/90/225 confirmed: {buck}')
# exact 19 dihedral-orbit representatives with raw multiplicities, per the
# obstruction table: all-distinct entries -> 30; exactly two equal -> 15;
# (4,4,4) -> 5.
expected_orbits = {
    (0, 0, 12): 15, (0, 1, 11): 30, (0, 2, 10): 30, (0, 3, 9): 30,
    (0, 4, 8): 30, (0, 5, 7): 30, (0, 6, 6): 15,
    (1, 1, 10): 15, (1, 2, 9): 30, (1, 3, 8): 30, (1, 4, 7): 30,
    (1, 5, 6): 30, (2, 2, 8): 15, (2, 3, 7): 30, (2, 4, 6): 30,
    (2, 5, 5): 15, (3, 3, 6): 15, (3, 4, 5): 30, (4, 4, 4): 5,
}
from collections import Counter
cnt_orbits = Counter()
for S in combinations(range(15), 3):
    cnt_orbits[gap_triple(S, 15)] += 1
assert dict(cnt_orbits) == expected_orbits, (
    {k: (cnt_orbits.get(k, 0), v) for k, v in expected_orbits.items()
     if cnt_orbits.get(k, 0) != v})
assert sum(expected_orbits.values()) == 455 and len(expected_orbits) == 19
for rep, m in expected_orbits.items():
    neq = len(set(rep))
    want = 30 if neq == 3 else (15 if neq == 2 else 5)
    assert m == want, (rep, m, want)
print(f'  19 dihedral orbit representatives confirmed with exact multiplicities')
# reflection-only quotient: (455 + 7) / 2 = 231 (7 reflection-fixed placements:
# one fixed vertex + one transposed pair among 3-subsets of a 15-cycle)
refl_fixed = 0
pivot = 0  # reflection fixing vertex 0, pairing i <-> -i mod 15
for S in combinations(range(15), 3):
    if tuple(sorted(((-x) % 15) for x in S)) == tuple(sorted(S)):
        refl_fixed += 1
assert refl_fixed == 7, f'reflection-fixed placements {refl_fixed} != 7'
assert (455 + 7) // 2 == 231
print(f'  reflection-fixed placements = 7; reflection-only quotient = 231 confirmed')

print('== (Bprime) n=18 triad: one size-3 + fifteen singles, R=16 ==')
signs = {}
for t in range(16):
    sizes = [3 if i == t else 1 for i in range(16)]
    ok, why = feasible(sizes)
    assert not ok, f'triad placement t={t} solution: {why}'
    a = sum((-1) ** i * (2 * s) for i, s in enumerate(sizes))
    assert abs(a) == 4, f'alt sum {a} at t={t}, expected +/-4'
    signs[t] = a
# document the sign pattern reported by the checker
pat = ''.join('+' if signs[t] > 0 else '-' for t in range(16))
print(f'  16/16 triad placements infeasible; alternating sums |a|=4, pattern {pat}')
# reflection-only quotient for triad placements on R=16: reflection t -> -t
# fixes t=0 and t=8, so (16 + 2)/2 = 9 orbits.
triad_refl_fixed = [t for t in range(16) if (-t) % 16 == t]
assert triad_refl_fixed == [0, 8], triad_refl_fixed
assert (16 + 2) // 2 == 9
print('  triad reflection-only quotient = 9 confirmed ((16+2)/2)')

print('CERTIFICATE OK: n=18 Q=3 legs (Aprime)+(Bprime) hold.')
