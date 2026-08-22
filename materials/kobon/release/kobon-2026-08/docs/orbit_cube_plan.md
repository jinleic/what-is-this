# n=14 concurrency-orbit CNF plan (next 12; no builds performed)

This is the low-cost planted-triple queue for the exact parallel-pattern tree.
The three already-running cubes are recorded first so that they are not queued a
second time.  A row's **signature** is the exact admissible lattice row used to
choose its capacity budget.  The positive `C` units plant the displayed
witnesses; all other `C` variables are intentionally left free.  Therefore a
SAT model can have a richer actual signature (the coverage/nesting argument in
`concurrency_cubes.md` requires this monotonicity).

## 1. Existing running cubes (do not rebuild)

`concurrency_cubes.md` records these as built and solving:

| cube | exact P units | C+ unit | orbit interpretation |
|---|---|---|---|
| `n14c-q0tp1` | `P+ = {Pij : 0 <= i < j < 14}` | `C(0,1,2)` | Q=0, orbit `()` |
| `n14c-q1tp1d` | `P01=0`; every other `Pij=1` | `C(2,3,4)` | Q=1, pattern `(2)`, orbit `()` |
| `n14c-q1tp1s` | `P01=0`; every other `Pij=1` | `C(0,2,3)` | Q=1, pattern `(2)`, orbit `(2)` |

Here `Pij=1` means the pair crosses (not parallel), and `Pij=0` means the
pair is parallel.  Thus Q=0 has one pattern and one triple orbit, while Q=1
has one pattern and two triple orbits.  The three running cubes exhaust all
Q<=1 planted-triple orbits.

## 2. Orbit count and exact representatives

The orbit key is the sorted multiset of parallel-class sizes used by the
planted triple.  A triple uses at most one line from each parallel class.
Labels below are the canonical labels used by `orbit_units`: parallel classes
come first, then singleton lines.

* Q=2: `2 = C(2,2)+C(2,2) = 1+1`, so the only pattern is `(2,2)`.
  Its available singleton count is `14-2-2=10`; the possible chosen-class
  counts are r=0,1,2.  Hence the three orbit keys are
  `()`, `(2)`, `(2,2)`.
* Q=3: either `3 = 1+1+1`, giving `(2,2,2)`, or `3=C(3,2)`, giving
  `(3)`.  Pattern `(2,2,2)` has singleton count `14-2-2-2=8` and r=0,1,2,3,
  hence four orbit keys `()`, `(2)`, `(2,2)`, `(2,2,2)`.  Pattern `(3)` has
  singleton count `14-3=11` and r=0,1, hence `()`, `(3)`.

Thus the Q=2/Q=3 queue has `3+4+2=9` distinct planted-triple orbit cubes.
The final three rows below are N3=2 variants of shared-line orbits, giving
`9+3=12` queued cubes.

For compact exact unit notation, let

```
ALLP = {Pij : 0 <= i < j < 14}          (|ALLP| = C(14,2) = 91)
P-(A) = ALLP \ A                        (negative P units)
```

For each row, `P+ = A` means `Pij=0` for every pair in A and `P-=P-(A)`
means `Pij=1` for every pair outside A.  `C+` lists all positive C units;
`C- = empty` means that no negative C units are asserted.  This is the
complete P/C unit list without expanding the 89 or 88 complement literals in
every row.

## 3. The next 12 orbit cubes

The command column is intentionally an execution line, not an executed build.
Before launch, Main must register the row's `(classes, C+)` recipe in
`concurrency_cases.py:CUBES` using the same pattern as the three running
entries.  Then `concurrency_cases.py build NAME` followed by
`solve_sat_capture.sh NAME` is the existing build/DRAT/model pipeline.

The size estimates use the low-cost Q2/Q3 per-line budget family: the checked
Q2 placement reference is `332,861 vars / 5,790,173 clauses`, and the checked
Q3-paired reference is `324,175 vars / 5,772,809 clauses`.  Consequently every
row is expected to be approximately `340k vars / 5.8M clauses`; exact headers
for the planted-C variant are **PENDING-VERIFICATION**.

| # / orbit index | cube name | exact signature `(Q,N3,N4,N5,N6,N7)` | parallel classes; orbit key; representative triple(s) | exact units (`P+`; `C+`) | expected vars / clauses | runnable build + solve line |
|---:|---|---|---|---|---|---|
| Q2-o00 | `n14c-q2p22tp1d` | `(2,1,0,0,0,0)` | `[[0,1],[2,3]]`; `()`; `C(4,5,6)` | `P+={01,23}`; `C+={C(4,5,6)}` | `~333k / ~5.79M` (reference `332861/5790173`) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q2p22tp1d && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q2p22tp1d` |
| Q2-o01 | `n14c-q2p22tp1s1` | `(2,1,0,0,0,0)` | `[[0,1],[2,3]]`; `(2)`; `C(0,4,5)` | `P+={01,23}`; `C+={C(0,4,5)}` | `~333k / ~5.79M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q2p22tp1s1 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q2p22tp1s1` |
| Q2-o02 | `n14c-q2p22tp1s2` | `(2,1,0,0,0,0)` | `[[0,1],[2,3]]`; `(2,2)`; `C(0,2,4)` | `P+={01,23}`; `C+={C(0,2,4)}` | `~333k / ~5.79M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q2p22tp1s2 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q2p22tp1s2` |
| Q3p-o00 | `n14c-q3p222tp1d` | `(3,1,0,0,0,0)` | `[[0,1],[2,3],[4,5]]`; `()`; `C(6,7,8)` | `P+={01,23,45}`; `C+={C(6,7,8)}` | `~324k / ~5.77M` (reference `324175/5772809`) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p222tp1d && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p222tp1d` |
| Q3p-o01 | `n14c-q3p222tp1s1` | `(3,1,0,0,0,0)` | `[[0,1],[2,3],[4,5]]`; `(2)`; `C(0,6,7)` | `P+={01,23,45}`; `C+={C(0,6,7)}` | `~324k / ~5.77M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p222tp1s1 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p222tp1s1` |
| Q3p-o02 | `n14c-q3p222tp1s2` | `(3,1,0,0,0,0)` | `[[0,1],[2,3],[4,5]]`; `(2,2)`; `C(0,2,6)` | `P+={01,23,45}`; `C+={C(0,2,6)}` | `~324k / ~5.77M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p222tp1s2 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p222tp1s2` |
| Q3p-o03 | `n14c-q3p222tp1s3` | `(3,1,0,0,0,0)` | `[[0,1],[2,3],[4,5]]`; `(2,2,2)`; `C(0,2,4)` | `P+={01,23,45}`; `C+={C(0,2,4)}` | `~324k / ~5.77M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p222tp1s3 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p222tp1s3` |
| Q3t-o00 | `n14c-q3p3tp1d` | `(3,1,0,0,0,0)` | `[[0,1,2]]`; `()`; `C(3,4,5)` | `P+={01,02,12}`; `C+={C(3,4,5)}` | `~324k / ~5.77M` (PENDING exact triad header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p3tp1d && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p3tp1d` |
| Q3t-o01 | `n14c-q3p3tp1s3` | `(3,1,0,0,0,0)` | `[[0,1,2]]`; `(3)`; `C(0,3,4)` | `P+={01,02,12}`; `C+={C(0,3,4)}` | `~324k / ~5.77M` (PENDING exact triad header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p3tp1s3 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p3tp1s3` |
| Q2-o01-N3=2 | `n14c-q2p22tp2s1` | `(2,2,0,0,0,0)` | `[[0,1],[2,3]]`; shared-line orbit `(2)`; `C(0,4,5)`, `C(0,6,7)` | `P+={01,23}`; `C+={C(0,4,5),C(0,6,7)}` | `~333k / ~5.79M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q2p22tp2s1 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q2p22tp2s1` |
| Q3p-o03-N3=2 | `n14c-q3p222tp2s3` | `(3,2,0,0,0,0)` | `[[0,1],[2,3],[4,5]]`; shared-line orbit `(2,2,2)`; `C(0,2,4)`, `C(1,3,5)` | `P+={01,23,45}`; `C+={C(0,2,4),C(1,3,5)}` | `~324k / ~5.77M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p222tp2s3 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p222tp2s3` |
| Q3t-o01-N3=2 | `n14c-q3p3tp2s3` | `(3,2,0,0,0,0)` | `[[0,1,2]]`; shared-line orbit `(3)`; `C(0,3,4)`, `C(1,5,6)` | `P+={01,02,12}`; `C+={C(0,3,4),C(1,5,6)}` | `~324k / ~5.77M` (PENDING exact planted-C header) | `VPY=/Users/jinleic/jinleic-workspace/scratch/kobon-audit/venv/bin/python; "$VPY" /Users/jinleic/jinleic-workspace/scratch/kobon/n14/concurrency_cases.py build n14c-q3p3tp2s3 && /Users/jinleic/jinleic-workspace/scratch/kobon/n14/solve_sat_capture.sh n14c-q3p3tp2s3` |

In every row, the omitted P units are exactly `Pij=1` for `Pij in P-(A)`;
there are no omitted parallel pairs.  The C units are all sorted triples.
The two C triples in each N3=2 row intersect in at most one line and share no
pair, so they do not force an A3 closure into N4:

```
{0,4,5} ∩ {0,6,7} = {0};
{0,2,4} ∩ {1,3,5} = empty;
{0,3,4} ∩ {1,5,6} = empty.
```

## 4. Exact budget checks for the tagged signatures

For a signature with no N4..N7 bulk points, the three relevant inequalities
reduce to

```
B1: C + 2Q - 3*N3 <= 6,
B2: Q + N3 <= 24,
B3: 3*N3 <= 91-Q.
```

Base Q=2 rows `(Q,N3)=(2,1)`:

```
B1: C + 2*2 - 3*1 <= 6  => C + 1 <= 6 => C <= 5;
B2: 2 + 1 = 3 <= 24;
B3: 3*1 = 3 <= 91-2 = 89.
```

Base Q=3 rows `(Q,N3)=(3,1)`:

```
B1: C + 2*3 - 3*1 <= 6  => C + 3 <= 6 => C <= 3;
B2: 3 + 1 = 4 <= 24;
B3: 3*1 = 3 <= 91-3 = 88.
```

N3=2 Q=2 row `(2,2,0,0,0,0)`:

```
B1: C + 2*2 - 3*2 <= 6  => C - 2 <= 6 => C <= 8;
B2: 2 + 2 = 4 <= 24;
B3: 3*2 = 6 <= 91-2 = 89.
```

N3=2 Q=3 rows `(3,2,0,0,0,0)`:

```
B1: C + 2*3 - 3*2 <= 6  => C + 0 <= 6 => C <= 6;
B2: 3 + 2 = 5 <= 24;
B3: 3*2 = 6 <= 91-3 = 88.
```

The exact orbit enumeration and these integer checks are hand-derived from
the formulas in `concurrency_cubes.md`; no build, solver, or local subprocess
was run for this plan.

## 5. Dead/duplicate pattern ledger

* **Duplicate, do not enqueue:** Q=0 `()` and Q=1 `(2)` orbits are exactly the
  three `n14c-q*tp1*` cubes already built and solving.  Rebuilding any of those
  would add no coverage.
* **Not evidence of death:** `q2pos_index.json` and `q3pos_index.json` describe
  no-concurrency cubes (`C(t)=0` for every triple).  Their prior Q=2/Q=3
  placement work does not refute any row here, because every row above asserts
  at least one positive C unit.  In particular, the Q3 certificate-equivalent
  placement skipped by `q3_positions.py` is only a no-concurrency skip, not a
  dead planted-triple orbit.
* **No selected Q2/Q3 planted orbit is dead by current evidence.**  All nine
  base orbits are retained for coverage; the three N3=2 variants are the
  closest shared-line analogues to the running `q1tp1d/q1tp1s` cubes.

## 6. Verification queue (when CPU slots free)

1. Register the twelve `(classes, C+)` recipes in `concurrency_cases.py:CUBES`
   under the names above.  Preserve the Q2/Q3 low-cost per-line budget mode;
   the old generic all-concurrency recipe has a larger header and is not the
   `~340k/~5.8M` estimate used here.
2. For each row, run the exact command in the table.  Record the generated
   DIMACS header (`cnf.nv`, `len(cnf.clauses)`) and replace the corresponding
   PENDING-VERIFICATION marker only after the build.  No such build was run in
   preparing this file.
3. Run `solve_sat_capture.sh NAME`; on SAT, extract/straighten and run
   `engine.verify_selection(14, ms, bs, sel, minimum=54)`, then compute the
   model's actual `(Q,N3..N7)` signature.  On UNSAT, record discovery UNSAT;
   do not promote it to a theorem until the capacity-theorem audit and full
   orbit coverage are complete.
