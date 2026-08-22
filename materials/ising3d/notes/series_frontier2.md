# Second-order finite-prefix series-structure frontier

## Exact input and fixed search protocol

[COMPUTATION] The input artifacts are the exact simple-cubic free-energy prefixes `results/series/extended2_sc_ht_free_energy.json` through `v^22` (23 coefficients, including exact zeros) and `results/series/extended2_sc_lt_free_energy.json` through `x^32` (33 coefficients). Their complete coefficient arrays and SHA-256 hashes are copied into `results/series/frontier2.json`. The searches use the original variables `v` and `x`; no critical-point value or other external numerical input is used for fitting, budget selection, or validation.

[LEMMA] Let an ansatz have `U` homogeneous rational coefficients and let its known coefficient equations form a matrix over `Q`. A nonzero `U x U` minor proves exact full column rank and therefore excludes every nonzero relation in that finite ansatz space. The experiment computes rational ranks and maximal-minor determinants by rowwise denominator clearing and fraction-free integer Bareiss elimination.

[COMPUTATION] Every frontier row uses the value-independent rule inherited from `e43_series_structure.py`: train on coefficient equations `0,...,U+1` and reserve every later derivative-safe equation as strict holdout. Thus every training system is overdetermined by two equations in its nominal budget. A training kernel is never reported as a discovery: later coefficients either eliminate it exactly, or the remaining directions receive an explicit truncation-unobservability certificate. Degree grids, the `U+2` split, and holdout orders depend only on the ansatz dimension and available series length.

## Second-order differential-algebraic frontier

The second-order ansatz is

\[
P(t,F,F',F'')=\sum_{a,b,c,d}p_{abcd}t^aF^b(F')^c(F'')^d=0.
\]

[LEMMA] If `N` coefficients of `F` through `t^(N-1)` are known, the occurrence of `F''` leaves exactly `N-2` derivative-safe coefficient equations, through order `t^(N-3)`. The next equation would require an unknown coefficient of `F`.

[COMPUTATION] Two nested grids were exhausted subject to `U+2 <= N-2`:

1. total degree `a+b+c+d <= D`;
2. jet bidegree `0 <= a <= d_t` and `b+c+d <= D_jet`.

Every space contains an `F''` monomial. The complete exact frontier is:

| series | safe equations | rows | largest `U` | `NO_RELATION_AT_BUDGET` | `CANDIDATE_REFUTED_BY_HOLDOUT` | `TRUNCATION_UNOBSERVABLE` | unexplained |
|---|---:|---:|---:|---:|---:|---:|---:|
| HT, `t=v` | 21 | 7 | 16 | 6 | 1 | 0 | 0 |
| LT, `t=x` | 31 | 12 | 28 | 8 | 4 | 0 | 0 |

[COMPUTATION] The HT `(d_t,D_jet)=(0,2)` training kernel is eliminated first at withheld order `v^12`; this holdout refutation is independently reconstructed by `tests/test_series_frontier2.py`. On the LT side, four training kernels are eliminated at withheld orders `x^6`, `x^10`, `x^14`, and `x^26`. All remaining second-order rows already have exact full-rank training-prefix minors.

[COMPUTATION] Consequently no observable second-order candidate survives in either searched frontier. This bounds only the displayed finite total-degree and jet-bidegree spaces. It does **not** prove that either infinite free energy is differentially transcendental or even non-differentially-algebraic.

## Linear ODE frontier through the true data budget

The homogeneous Euler-form ansatz is

\[
\sum_{j=0}^{r}Q_j(t)\theta^jF(t)=0,
\qquad \theta=t\frac{d}{dt},\quad r\le 2,
\]

with a common polynomial bound `deg Q_j <= d`. All degrees satisfying `(r+1)(d+1)+2 <= N` were tested. This extends the earlier HT reduced-series scan beyond its parameter budget 11 and also supplies the full LT table.

[COMPUTATION] The maximal tested polynomial degrees at orders `r=0,1,2` are respectively `(20,9,6)` for HT and `(30,14,9)` for LT. The complete table is:

| series | equations | rows | largest `U` | `NO_RELATION_AT_BUDGET` | `CANDIDATE_REFUTED_BY_HOLDOUT` | `TRUNCATION_UNOBSERVABLE` | unexplained |
|---|---:|---:|---:|---:|---:|---:|---:|
| HT, `t=v` | 23 | 38 | 21 | 29 | 8 | 1 | 0 |
| LT, `t=x` | 33 | 56 | 31 | 0 | 48 | 8 | 0 |

[COMPUTATION] The apparent top-budget kernels are interpolation forced by sparse series support, not ODE discoveries. For example, the HT series is supported on even powers. In the order-two, degree-six row, the even and odd polynomial-coefficient sectors decouple; the surviving exact direction vanishes through `v^22`, while its first unfixed supported coefficient is `v^24`. The JSON certificate partitions columns by their exact support congruence class, stores exact component ranks and kernel bases, and identifies the first unavailable supported order. The LT prefixes have analogous parity-supported directions; four order-zero top-degree rows are also literal zero-column cases because their first possible terms lie beyond `x^32`.

[COMPUTATION] Therefore every linear-ODE kernel is either killed by an exact holdout or certified as truncation-unobservable. This finite statement does **not** prove non-D-finiteness of either true free energy.

## Mahler-type frontier

The exotic functional ansatz is

\[
M(t,F(t),F(t^2))=
\sum_{a,b,c}m_{abc}t^aF(t)^bF(t^2)^c=0,
\qquad b+c\le D,
\]

with dependent total degree `D=1` or `2`. For each `D`, every polynomial degree in `t` allowed by `U+2 <= N` was tested.

[COMPUTATION] The complete low-degree Mahler-type frontier is:

| series | equations | rows | largest `U` | `NO_RELATION_AT_BUDGET` | `CANDIDATE_REFUTED_BY_HOLDOUT` | `TRUNCATION_UNOBSERVABLE` | unexplained |
|---|---:|---:|---:|---:|---:|---:|---:|
| HT, `t=v` | 23 | 10 | 21 | 8 | 2 | 0 | 0 |
| LT, `t=x` | 33 | 15 | 30 | 6 | 8 | 1 | 0 |

[COMPUTATION] The sole full-system LT kernel, at `(degree_t,D)=(4,2)`, is exactly truncation-unobservable under its certified parity support: the known prefix fixes every available supported direction, and the first unfixed supported order is beyond `x^32`. No Mahler candidate remains for parent review. This is only a finite low-dependent-degree exclusion, not a functional-transcendence theorem.

## Certificates, holdouts, and scope

[COMPUTATION] Across all six HT/LT family tables, 138 exact budget rows were completed: 57 have full-rank training-prefix minors, 71 training candidates are refuted by withheld exact coefficients, and 10 kernels are certified truncation-unobservable. There are zero unexplained or observable survivors. Every row points to one certificate under `/data/certificates/Cxxxx` in `results/series/frontier2.json`; negative certificates store nonzero rational maximal minors, while support-unobservable certificates store exact component ranks, kernel bases, support progressions, and first unfixed orders.

[COMPUTATION] `.venv/bin/python tests/test_series_frontier2.py` independently reconstructs four representative `NO_RELATION_AT_BUDGET` matrices and their stored rational determinants, reruns the HT second-order holdout refutation from scratch with a separate exact `Fraction` elimination, checks a parity-support truncation certificate, validates source hashes and the result envelope, and confirms zero unexplained survivors.

[UNRESOLVED] These negative frontiers bound **structure searches**, not the mathematical structure of the unknown infinite functions. Nothing here proves non-D-finiteness, non-differential-algebraicity, differential transcendence, absence of higher-degree Mahler relations, or an exact solution of the three-dimensional Ising model.
