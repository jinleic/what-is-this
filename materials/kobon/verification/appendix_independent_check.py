#!/usr/bin/env python3
"""STANDALONE machine-checkable certificate for Theorem M
(pure stdlib; parametrized cycle solver with exact integer arithmetic).

Target system (derived in escape_obstruction.md): the multiplicities
m_i = #(double-extreme vertices at direction classes i, i+1 on the direction
cycle of a saturated (n, Q, simple) arrangement) satisfy

    m_{i-1} + m_i = 2 s_i,     0 <= m_i <= s_i * s_{i+1},    m_i integer.

Certificate legs:
  (A) n=14 paired layout: three size-2 classes + eight size-1, R=11 (odd
      cycle, unique integer solution): all C(11,3)=165 placements infeasible.
  (B) n=14 TRIAD layout: one size-3 class + ELEVEN size-1, R=12 (EVEN cycle):
      consistency condition  sum_i (-1)^i d_i == -4 != 0  for every placement
      of the unique 3-class => NO solution, not even rational.
      (An earlier revision of this file mistakenly used R=11 with one size-3
      => 13 lines; thanks to Main for the catch. The 14-line triad is now the
      certified object.)
  (C) controls which MUST stay feasible (guard against an over-strong
      checker): known equality arrangements + even-cycle positive examples.
"""
from fractions import Fraction as Fr
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


print('== (C) controls: must stay feasible ==')
for sizes, tag in [([1] * 9, 'Blanc-21 at n=9'), ([1] * 15, 'rot-sym-65 at n=15'),
                   ([1] * 5, 'pentagon'), ([1] * 3, 'triangle'),
                   ([2, 2, 2], 'hexagon (n=6, Q=3)'),
                   ([1, 1, 1, 1], 'even R=4 all-single (m=1s)')]:
    ok, why = feasible(sizes)
    assert ok, f'control {tag} wrongly infeasible: {why}'
    print(f'  {tag}: feasible [{why[:60]}]')

print('== (A) n=14 paired: three size-2 + eight size-1, R=11 ==')
cnt = 0
for S in combinations(range(11), 3):
    sizes = [2 if i in S else 1 for i in range(11)]
    ok, why = feasible(sizes)
    assert not ok, f'paired pattern {S} FEASIBLE: {why}'
    cnt += 1
print(f'  {cnt}/165 placements infeasible')

print('== (B) n=14 triad: one size-3 + eleven size-1, R=12 (14 lines) ==')
for t in range(12):
    sizes = [3 if i == t else 1 for i in range(12)]
    ok, why = feasible(sizes)
    assert not ok, f'triad placement {t} FEASIBLE: {why}'
print(f'  12/12 placements infeasible [even-cycle alternating sum = +-4]')
print('  note: this is the 14-LINE triad (one 3-class + 11 singles); an earlier')
print('        13-line leg (R=11, one 3-class) has been retired.')
print('CERTIFICATE OK: Theorem M holds for the full Q=3 branch at n=14.')
