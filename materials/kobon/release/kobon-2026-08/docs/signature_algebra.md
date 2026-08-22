# n=20, T=117 — degeneracy signature algebra (hand derivation)

Status: **hand-derived; no local enumeration or solver execution**. The
verification queue at the end names the one computation still required.

## 1. Budgets and exact slack

Let `Q` be the number of parallel pairs, let `N_k` be the number of
multiplicity-`k` finite points, and let `C` be the selected-triangle /
outside-line open-interior incidence count. For `n=20`,

\[
\binom{20}{2}=190,\qquad 20(20-2)=360,
\qquad 3T=3\cdot117=351,
\]
so the segment slack is

\[
\Delta=360-351=9.
\]

For the requested `k=3,\dots,7` coordinates, the three integer budgets are

\[
\begin{aligned}
\tag{B1}
 C+2Q-3N_3+0N_4+5N_5+12N_6+21N_7&\le 9,\\
\tag{B2}
 Q+N_3+3N_4+6N_5+10N_6+15N_7&\le 54,\\
\tag{B3}
 3N_3+6N_4+10N_5+15N_6+21N_7&\le 190-Q.
\end{aligned}
\]

Here `54` is the face slack
\[
\binom{19}{2}-117=171-117=54.
\]
The corresponding bounded-face expression is
\[
F_b=171-Q-\bigl(N_3+3N_4+6N_5+10N_6+15N_7+\cdots\bigr).
\]

The corrected distinct-finite-point count (the point-count correction from
Sections 1–3 of the n=14 note) is
\[
\begin{aligned}
\mathrm{pts}
 &=190-Q-\sum_{k\ge3}\left(\binom{k}{2}-1\right)N_k\\
 &=190-Q-(2N_3+5N_4+9N_5+14N_6+20N_7
             +27N_8+35N_9+44N_{10}+54N_{11}).
\end{aligned}
\]
The `B2` coefficients are instead `\binom{k-1}{2}`; confusing these two
penalties undercounts the point reduction by `k-2` per multipoint.

For completeness, the general coefficients needed if a point has `k≥8`
are
\[
\begin{array}{c|rrrrrrrrr}
k&3&4&5&6&7&8&9&10&11\\ \hline
k(k-4)&-3&0&5&12&21&32&45&60&77\\
\binom{k-1}{2}&1&3&6&10&15&21&28&36&45\\
\binom{k}{2}&3&6&10&15&21&28&36&45&55\\
\binom{k}{2}-1&2&5&9&14&20&27&35&44&54
\end{array}
\]
Thus a formula written only through `N_7` is the requested core lattice; the
`N_8,N_9,N_{10}` terms above give the complete n=20 extension rather than
silently assuming the n=14 cutoff `k≤7`.

## 2. Complete symbolic parameterization of the signature lattice

Write the bulk vector as
\[
b=(N_4,N_5,N_6,N_7)
\tag{bulk}
\]
and put `m=N_3`. Define its three bulk costs
\[
 s(b)=5N_5+12N_6+21N_7,
\quad f(b)=3N_4+6N_5+10N_6+15N_7,
\quad p(b)=6N_4+10N_5+15N_6+21N_7.
\]
For every nonnegative integer `m`, the exact permitted `Q` interval is
\[
\boxed{
0\le Q\le q_b(m):=
\min\left\{
\left\lfloor\frac{9+3m-s(b)}2\right\rfloor,
54-f(b)-m,
190-p(b)-3m
\right\}.
}
\tag{2.1}
\]
It is nonempty exactly for
\[
L_b:=\max\left(0,\left\lceil\frac{s(b)-9}{3}\right\rceil\right)
\le m\le
U_b:=\min\left(54-f(b),\left\lfloor\frac{190-p(b)}3\right\rfloor\right).
\tag{2.2}
\]
Equations (2.1)–(2.2), with `N_4,N_5,N_6,N_7,m` integral, are a complete
hand-checkable enumeration of all requested `(Q,N_3,\dots,N_7)` signatures:
for each allowed `(b,m)`, take every integer `Q` from zero through `q_b(m)`.
The corresponding available crossing slack is, exactly,
\[
0\le C\le 9-2Q+3m-s(b),
\tag{2.3}
\]
so every signature admitted by (2.1) has at least the algebraic witness
`C=0`.

For a bulk vector including `N_8,N_9,N_{10}` (and any hypothetical `N_{11}`),
use the same equations after replacing
\[
\begin{aligned}
 s(b)&\mapsto s(b)+32N_8+45N_9+60N_{10}+77N_{11},\\
 f(b)&\mapsto f(b)+21N_8+28N_9+36N_{10}+45N_{11},\\
 p(b)&\mapsto p(b)+28N_8+36N_9+45N_{10}+55N_{11}.
\end{aligned}
\tag{2.4}
\]
This is the symbolic completion beyond the displayed `N_3..N_7`
coordinates.

### Maximal `Q`

From `B1` alone (set `C=0`),
\[
2Q\le 9+3m-s(b),
\qquad
Q\le\left\lfloor\frac{9+3m-s(b)}2\right\rfloor.
\tag{2.5}
\]
In the bulk-free slice this is the exact B1-only maximum
\[
Q\le\left\lfloor\frac{9+3m}{2}\right\rfloor
 =4+\left\lceil\frac{3m}{2}\right\rceil
 =4+m+\left\lfloor\frac{m+1}{2}\right\rfloor.
\]
Combining all three budgets in that slice gives
\[
Q\le\min\left\{
4+\left\lceil\frac{3m}{2}\right\rceil,
54-m,
190-3m\right\}.
\tag{2.6}
\]
The maximum over all bulk-free signatures is therefore `Q=34`, attained at
`m=N_3=20` (both the B1 and B2 bounds equal 34 there). Positive bulk costs
only decrease the bound according to (2.1); in particular, `N_4` lowers the
face bound but not B1.

### Largest possible multiplicity `k`

The face budget for one `k`-fold point gives
\[
\binom{k-1}{2}\le54.
\]
Since `\binom{10}{2}=45` and `\binom{11}{2}=55`, this face-only test says
`k≤11`.
The segment budget refines it. With one `k`-fold point, `m=N_3` triple
points, and no other positive-cost bulk, `Q=C=0` is the least-cost case and
\[
k(k-4)\le9+3m,
\qquad
m\ge\left\lceil\frac{k(k-4)-9}{3}\right\rceil.
\tag{2.7}
\]
In particular, with `m=0`, `k(k-4)≤9`; the positive root is
`2+\sqrt{13}<6`, hence the exact integer cutoff is `k≤5`.
Checking the face upper bound `m≤54-\binom{k-1}{2}` against (2.7) gives

\[
\begin{array}{c|c|c|c}
k&k(k-4)&m\text{ required by B1}&m\text{ allowed by B2}\\ \hline
8&32&m\ge8&m\le33\\
9&45&m\ge12&m\le26\\
10&60&m\ge17&m\le18\\
11&77&m\ge23&m\le9\quad(\text{impossible})
\end{array}
\]
Thus the coupled n=20 algebra has **maximum multiplicity `k=10`**, not 11:
`k=10` survives at `(N_{10},N_3)=(1,17)` or `(1,18)` with `Q=0` (and the
first row has `C=0` exactly). Pure high-k checks, useful for auditing the
extension in (2.4), are
\[
\begin{array}{c|c|c|c}
\text{bulk}&N_3\text{ range}&Q\text{ maximum}&\#(Q,N_3)\\ \hline
(N_8=1)&8..33&15&216\\
(N_9=1)&12..26&8&72\\
(N_{10}=1)&17..18&0&2.
\end{array}
\]
These counts are direct sums of the intervals in (2.1), not machine output.
`N_{11}=1` is impossible by the displayed `m≥23` versus `m≤9` clash.
Mixed high-k vectors remain covered exactly by (2.4).

## 3. Bulk-class volume table (hand arithmetic)

Here `volume` means the number of `(Q,N_3)` pairs in a fixed bulk class,
computed as `\sum_{m=L_b}^{U_b}(q_b(m)+1)`. The table gives the ten largest
classes in the requested `k≤7` lattice. The `p` column is included to show
that `B3` was checked; for every row below the active minimum is the B1/B2
minimum, and the pair bound is slack throughout. The arithmetic in the last
column is the actual sum, so no total-lattice claim is being made.

In the `active q_b` column, `A=\lfloor(9+3m-s)/2\rfloor` and
`F=54-f-m`; an interval written `A[a..b]; F[c..d]` means those are the
minimizing bounds on the indicated integer `m` ranges.

| bulk vector `b` | `(s,f,p)` | `N_3` range | active `q_b(m)` | `Q` range | volume (exact sum) |
|---|---:|---:|---|---:|---:|
| `()` | `(0,0,0)` | `0..54` | `A[0..19]; F[20..54]` | `0..34` | `390+630=1020` |
| `(4,1)` | `(0,3,6)` | `0..51` | `A[0..18]; F[19..51]` | `0..32` | `356+561=917` |
| `(4,2)` | `(0,6,12)` | `0..48` | `A[0..17]; F[18..48]` | `0..30` | `324+496=820` |
| `(5,1)` | `(5,6,10)` | `0..48` | `A[0..18]; F[19..48]` | `0..29` | `309+465=774` |
| `(4,3)` | `(0,9,18)` | `0..45` | `A[0..16]; F[17..45]` | `0..28` | `293+435=728` |
| `(4,1)(5,1)` | `(5,9,16)` | `0..45` | `A[0..17]; F[18..45]` | `0..27` | `279+406=685` |
| `(4,4)` | `(0,12,24)` | `0..42` | `A[0..15]; F[16..42]` | `0..27` | `276+378=654` |
| `(4,2)(5,1)` | `(5,12,22)` | `0..42` | `A[0..16]; F[17..42]` | `0..26` | `251+351=602` |
| `(6,1)` | `(12,10,15)` | `1..44` | `A[1..18]; F[19..44]` | `0..25` | `243+351=594` |
| `(4,5)` | `(0,15,30)` | `0..39` | `A[0..13]; F[14..39]` | `0..25` | `210+351=561` |

Near-tie omitted from the ten-row cutoff: `(5,2)` has
`(s,f,p)=(10,12,20)`, `N_3=1..42`, `q_b=A[1..17];F[18..42]`, and
`234+325=559` signatures. This is recorded to make the ranking auditable,
not as a claim about the unlisted classes.

## 4. Eight cheapest corner cubes / readiness plan

The existing n20 builder already has seven exact parallel-pattern/no-
concurrency corners. They are the cheapest structural corners because they
use only `P` units (and the all-negative concurrency units already emitted by
the builder). In the table, `P(i,j)-` means the pair is parallel and every
unlisted pair receives `P+`; `C(t)-` means the concurrent-triple literal is
forbidden. `C` in the budget is the incidence count, not the `C(t)` literal.

| priority | cube corner | exact units / signature corner | readiness when a solver slot frees |
|---:|---|---|---|
| 1 | `q0` | all `P+`; all `C(t)-`; `(Q,N_3)=(0,0)` | CNF already present as `n20_cube_q0_t117.cnf`; launch first if the no-parallel corner is free. |
| 2 | `q1` | `P(0,1)-`, all other `P+`; all `C(t)-`; `(1,0)` | Existing `n20_cube_q1_t117.cnf`; one parallel pair, `C≤7`. |
| 3 | `q2` | `P(0,1)-,P(2,3)-`, all other `P+`; all `C(t)-`; `(2,0)` | Existing `n20_cube_q2_t117.cnf`; two disjoint pairs, `C≤5`. |
| 4 | `q3paired` | `P(0,1)-,P(2,3)-,P(4,5)-`, all other `P+`; all `C(t)-`; `(3,0)` | Existing `n20_cube_q3paired_t117.cnf`; three pair classes, `C≤3`. |
| 5 | `q3triad` | all three `P` pairs inside `{0,1,2}` negative, all other `P+`; all `C(t)-`; `(3,0)` | Existing `n20_cube_q3triad_t117.cnf`; one parallel 3-class, `C≤3`. |
| 6 | `q4four` | four negative pairs `(0,1),(2,3),(4,5),(6,7)`, all other `P+`; all `C(t)-`; `(4,0)` | Existing `n20_cube_q4four_t117.cnf`; tight no-concurrency corner, `C≤1`. |
| 7 | `q4pairtriad` | `P(0,1)-` plus all three internal pairs of `{2,3,4}` negative, all other `P+`; all `C(t)-`; `(4,0)` | Existing `n20_cube_q4pairtriad_t117.cnf`; same `C≤1`, different parallel-class shape. |
| 8 | `q0-tp1` | all `P+`; `C(0,1,2)+` and all other `C(t)-`; exact `(Q,N_3)=(0,1)` | Not built by the current no-concurrency recipe. Add the one planted `C` unit to the q0 `P` pattern and queue first concurrency run; algebraic `C≤12`. The N3=2 sibling is the next cheap upgrade: add also `C(3,4,5)+` (all other `C(t)-`), giving `C≤15`. |

The first seven entries are the ordered next runs requested for the already
built CNFs. Entry 8 is the cheapest planted-concurrency signature; if two CNFs
are available, run its disjoint N3=2 sibling immediately after it. Shared-line
planted triples are intentionally deferred because closure can create an
N4-or-higher point and no longer gives a clean N3=2 corner.

## 5. Verification queue (one item; intentionally not run here)

- **PENDING-VERIFICATION:** run the exact `run-concurrency_cases.py-20`
  enumeration over the integer inequalities (B1)–(B3), including the
  `k=8,9,10` extension, and compare its per-bulk `(Q,N_3)` ranges and volumes
  against (2.1)–(2.4) and the hand rows above. This is the sole queued
  computation; no execution was performed for this artifact.
