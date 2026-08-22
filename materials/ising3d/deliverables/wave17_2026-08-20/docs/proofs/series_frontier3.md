# Third-frontier holdout evaluation of the frozen e56 Euler-ODE survivor

## Exact inputs and frozen candidate

[COMPUTATION] The canonical input is the exact simple-cubic high-temperature
reduced free energy `phi = log 2 + sum_n f_n v^n` through `v^28`
(`results/series/ht_v28.json`, 29 exact rational coefficients). The stored
prefix agrees coefficient-for-coefficient with the frozen e56 training prefix
`results/series/extended2_sc_ht_free_energy.json` through `v^22`, and the new
holdout coefficients agree with their own production artifacts:
`f_24 = 2135670379057/8`
(`results/series/ht_v24.json`) and
`f_26 = 115377512914251/26`
(`results/series/ht_v26.json`). SHA-256 hashes of all six source artifacts
are recorded in `results/series/frontier3.json`.

[COMPUTATION] The evaluated candidate is the sole HT row of the frozen e56
frontier not closed at the `v^22` budget: the homogeneous Euler-form linear
ODE `sum_j Q_j(v) theta^j F = 0` with `theta = v d/dv`, order `r = 2` and
`deg Q_j <= 6` (21 monomial columns). Its e56 verdict was
TRUNCATION_UNOBSERVABLE with first unfixed supported order `v^24`. The
producer verified it is the unique HT truncation-unobservable row across all
three e56 HT families (1 row in total).

## Predeclared holdout protocol

[LEMMA] The orders `23..28` were unavailable at the e56 freeze, so under the
inherited U+2 discipline they are strict holdouts: the training prefix of the
candidate remains exactly the equations at orders `0..22`, and no coefficient
beyond `v^22` may enter any refit, reselection, or kernel reconstruction.

[COMPUTATION] The decision tree was fixed before any holdout evaluation:
(i) a nonzero residual at order 24 classifies
`CANDIDATE_REFUTED_BY_HOLDOUT` with `v^26`/`v^28` recorded as confirmation
only; (ii) a zero residual at 24 triggers sequential evaluation at 26 then
28; (iii) if all three vanish, a predeclared expanded prefix-rank grid
(cells `(r,d)` with `(r+1)(d+1)+2 <= 29`, training orders `0..U+1`, holdouts
`U+2..28` preserved) is recomputed with no survival discovery and no
D-finiteness claims. The branch taken was
`REFUTED_AT_FIRST_UNOBSERVABLE_ORDER`; the expanded grid was not
executed.

## Reconstruction of the exact relation

[COMPUTATION] Rebuilding the 23x21 training matrix from the raw `v^0..v^22`
coefficients with the e56 entry rule
`[v^n](v^p theta^q F) = (n-p)^q f_(n-p)` gives exact rank
20 and a one-dimensional kernel
spanned by the primitive integer vector stored in the frozen artifact
(span agreement certified; the vector annihilates all 23 training equations
exactly). In polynomial form the candidate relation is
`Q0(v) F + Q1(v) theta F + Q2(v) theta^2 F = 0` with

* `Q0 = [415411587595614524886324598429103082240, -17542892800829250901062241025978616364800, 155942196950027554724099723299859143578240, -562403771025920547345517448419249201947520]` on `v^0, v^2, v^4, v^6`,
* `Q1 = [-271312938457991527181927220060150534042, 10416309270194065879763218917315987064674, -100711595067991261198116533001037767870252, 461536851620893363874832770964225826846920]`,
* `Q2 = [31803572330092132369382460422799496461, -721817054292921798409008630259147497042, 3765219061709894039282288376258559012611, -53418410543932338351276479457442506002918]`.

## Decisive holdout evaluation

[LEMMA] (Parity.) Every column with a nonzero kernel weight has even `p`,
and every odd-index coefficient of the HT series is exactly zero, so the
residual of the candidate at any odd order is identically zero before any
arithmetic: orders 23, 25, 27 carry no information.

[COMPUTATION] The first informative holdout order is 24. The exact residual
of the frozen kernel there is

    [v^24] (Q0 + Q1 theta + Q2 theta^2) F = 499752451221372610239349461368647432702945848808/11

which is a nonzero rational number: the candidate is refuted by the first
previously unobservable coefficient. Writing the residual as `S f_24 + T`,
where `S = sum_q c_(0,q) 24^q = 12222758726736886117284368520518000226768` is the exact weight the
relation puts on `f_24` and `T = -35892202921745856477746014532342520613720483516398134/11` is fixed by the training
prefix alone, the relation would force
`f_24 = 17946101460872928238873007266171260306860241758199067/67225172997052873645064026862849001247224`, whereas the exact coefficient is
`f_24 = 2135670379057/8`.

[COMPUTATION] The surviving-kernel dimension is 1 through order 23 and 0
from order 24 onward, so `first_refuting_holdout_order = 24`. The full
29-equation system has exact rank 21 (full column rank); the stored
nonzero 21x21 maximal minor on row orders
[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24] has exact determinant
`-31265665004583691032828566614210634515125637906910453679736021091221406669327547/284710580000000`. Orders 26 and 28 give further
nonzero exact residuals (500071378398815869213964078020051854907899278078520/143 and
22667554452529238570074088272879080635241164323868544/143) but are flagged
confirmation-only: they played no role in the classification.

[COMPUTATION] An independent evaluator (the coordinating parent agent)
recomputed the three residuals from the artifacts before this script ran;
all three match the values above exactly (recorded in
`results/series/frontier3.json` under `independent_crosscheck`).

## Controls on the strict holdout machinery

[COMPUTATION] All 8 strict
`CANDIDATE_REFUTED_BY_HOLDOUT` rows of the frozen e43 first-order
differential-algebraic LT frontier were re-audited from the stored reduced
coefficients: the e43 U+2 budget ladder reproduces every stored first
refuting holdout order and every stored nonzero maximal-minor determinant
exactly. One frozen e43 HT `NO_RELATION_AT_BUDGET` row was likewise
reproduced with full training-prefix rank and its stored minor. These
controls demonstrate the exact three-way classification semantics on known
data.

## Verdict semantics and scope

* `NO_RELATION_AT_BUDGET`: a nonzero exact maximal minor of the training
  prefix proves no nonzero relation exists in the ansatz space at this
  budget.
* `CANDIDATE_REFUTED_BY_HOLDOUT`: a training kernel exists but later exact
  coefficients kill it; here first killed at withheld order `v^24`.
* `CANDIDATE_SURVIVES`: an observable exact kernel remains after every
  known coefficient; no row of this experiment carries this verdict.

[UNRESOLVED] This is a finite-prefix statement about one frozen ansatz row.
It does not prove non-D-finiteness or non-differential-algebraicity of the
true 3D Ising free energy, does not exclude relations outside the searched
spaces, and is not an exact solution of the model. No survival discovery
was performed and none is claimed.
