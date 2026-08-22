# Finite-prefix algebraic and differential-algebraic structure certificates

## Series and exact design rule

[COMPUTATION] The inputs are the extended exact reduced free-energy series in `results/series/extended2_sc_ht_free_energy.json` and `results/series/extended2_sc_lt_free_energy.json`.  Odd powers vanish, so the searches use

\[
f_{\rm HT}(z)=\sum_{n=0}^{11}[v^{2n}]\phi_{\rm HT}\,z^n,\qquad z=v^2,
\]

with 12 known coefficients and valuation one, and

\[
f_{\rm LT}(u)=\sum_{n=0}^{16}[x^{2n}](\phi_{\rm LT}-3K)\,u^n,\qquad u=x^2,
\]

with 17 known coefficients and valuation three.  The source-file SHA-256 hashes and every coefficient are copied into `results/series/structure_certificates.json`.

[LEMMA] For an ansatz with `U` unknown rational coefficients, fewer than `U` homogeneous coefficient equations force a nonzero kernel by dimension alone.  Such a kernel is interpolation, not evidence for a relation.

[COMPUTATION] This sweep uses the stronger fixed rule: train on the first `U+2` exact coefficient equations, so every training matrix is overdetermined by at least two equations, and use every later known coefficient as a strict holdout.  The prefix length depends only on `U`; no holdout residual, candidate, critical benchmark, or external constant enters model selection.  Every rectangular degree tuple satisfying `U <= E-2`, where `E` is the number of derivative-safe equations, is included.  Thus the largest allowed budgets have no holdout left, but remain overdetermined by two equations; any kernel there is explicitly a `CANDIDATE_SURVIVES(!)` rather than a conclusion.

[LEMMA] If the `E x U` exact rational coefficient matrix has rank `U`, then no nonzero polynomial in that displayed ansatz annihilates the known prefix.  Each negative row in the result stores an explicit nonzero `U x U` rational minor.  If the training kernel is nonzero, the result also stores a primitive integer basis and its exact residual on every holdout equation.

## Algebraicity sweep

The algebraic ansatz is

\[
P(t,f)=\sum_{a=0}^{d_x}\sum_{b=0}^{d_f}p_{ab}t^a f^b,
\qquad U=(d_x+1)(d_f+1).
\]

[COMPUTATION] The full exact frontier is:

| series | equations | budget rows | `NO_RELATION_AT_BUDGET` | `CANDIDATE_REFUTED_BY_HOLDOUT` | `CANDIDATE_SURVIVES(!)` |
|---|---:|---:|---:|---:|---:|
| HT, `z=v^2` | 12 | 27 | 27 | 0 | 0 |
| LT, `u=x^2` | 17 | 45 | 28 | 7 | 10 |

[COMPUTATION] Therefore all 27 tested HT bidegrees have exact full-column-rank certificates.  On the LT side, 28 training matrices already have full column rank and seven prefix kernels are eliminated by later exact coefficients.

[COMPUTATION] **The ten LT survivors are not algebraic discoveries.**  In every survivor row, the full kernel is exactly the span of matrix columns that are identically zero at the available truncation.  Since `f_LT(u)=u^3+O(u^5)`, for example `f_LT(u)^6=O(u^18)`, while the data stop at `u^16`.  Thus the apparent polynomial `P=f^6` vanishes for the tautological reason that its first possible term has not yet been reached.  The JSON flags every such row with `CANDIDATE_SURVIVES(!)`, lists its zero columns and full kernel basis, and requests parent audit against new coefficients.

[COMPUTATION] Absence at one listed bidegree means only: **no polynomial of that bidegree explains these known coefficients**.  It is not a proof that either infinite series is transcendental, and it says nothing about a larger bidegree.

## First-order differential-algebraic probe

The differential-algebraic ansatz is

\[
Q(t,f,f')=\sum_{a=0}^{d_x}\sum_{b=0}^{d_f}
\sum_{c=0}^{d_{f'}}q_{abc}t^a f^b(f')^c,
\qquad d_{f'}\geq1,
\]

with rectangular parameter budget
`U=(d_x+1)(d_f+1)(d_fprime+1)`.  Powers of `f'` are allowed, but no derivative above first order occurs.

[LEMMA] If `N` coefficients of `f` are known, only `N-1` general coefficient equations involving `f'` are safe: the next equation would depend on the unknown next coefficient of `f`.  Hence the HT and LT probes have 11 and 16 equations respectively.

[COMPUTATION] The full exact first-order frontier is:

| series | derivative-safe equations | budget rows | `NO_RELATION_AT_BUDGET` | `CANDIDATE_REFUTED_BY_HOLDOUT` | `CANDIDATE_SURVIVES(!)` |
|---|---:|---:|---:|---:|---:|
| HT, `z=v^2` | 11 | 21 | 21 | 0 | 0 |
| LT, `u=x^2` | 16 | 45 | 29 | 8 | 8 |

[COMPUTATION] Every tested HT first-order differential-algebraic budget has an exact full-rank certificate.  On the LT side, 29 training matrices have full rank and eight prefix candidates fail strict holdouts.

[COMPUTATION] **The eight surviving LT rows are again truncation artifacts, not differential-algebraic discoveries.**  Here `f'_LT(u)=3u^2+O(u^4)`.  Monomials such as `(f')^8=O(u^16)` are invisible to derivative-safe equations ending at `u^15`, and mixed high powers such as `f^5 f'=O(u^17)` are likewise invisible.  Exact rank confirms that the survivor space is precisely the span of these zero columns.  All 72 stored survivor-basis vectors across the algebraic and differential-algebraic tables were independently re-evaluated on every available equation by `tests/test_series_structure.py`.

[COMPUTATION] These finite negative certificates extend the earlier linear D-finite probes only to the displayed first-derivative polynomial spaces.  They do not decide differential algebraicity, do not test higher derivatives, and do not constrain larger polynomial budgets.

## Interlayer extension through one more order

[COMPUTATION] An independent two-dimensional finite-lattice calculation was implemented specifically for the overlap susceptibility.  For each open rectangle it enumerates edge subsets by exact boundary-parity polynomial, forms every exact finite-box correlation ratio, sums the ordered pairs `sum_(r,s) G(r,s)^2`, and performs rectangular finite-lattice Möbius inversion.

[LEMMA] After the two partition-function denominators are removed by finite-lattice inversion, a contributing doubled-correlation cluster is a connected even multigraph: the two color copies have the same odd source pair, so their union has even degree at every vertex.  Any connected even multigraph spanning an `a x b` bounding rectangle contains a closed traversal that crosses each coordinate range in both directions and therefore uses at least `2((a-1)+(b-1))` edges, counted with multiplicity.  Hence every contribution through `v^8` has `(a-1)+(b-1) <= 4`.  All such rectangles are included, and adding one extra span unit changes no coefficient.

[COMPUTATION] This independent exact calculation gives

\[
\sum_rG(r)^2
 =1+4v^2+36v^4+236v^6+1556v^8+O(v^{10}).
\]

[COMPUTATION] Consequently,

\[
c_2^{\rm tot}(v)
 =\frac12+2v^2+18v^4+118v^6+778v^8+O(v^{10}),
\]

and

\[
c_2^{\rm res}(v)
 =2v^2+18v^4+118v^6+778v^8+O(v^{10}).
\]

[COMPUTATION] The independent anisotropic three-dimensional finite-lattice calculation `anisotropic_flm_c2_series(8)` gives the identical exact `Fraction` array.  Adding one span unit also leaves that array unchanged.  Thus the previous agreement through `v^6` is extended by one nonzero order, with exact new coefficient `[v^8]c_2=778`.

[THEOREM] In an open stack, an even number of selected interlayer bonds crosses every layer cut.  Equivalently, for an even-layer periodic stack, flipping every other layer maps `w=tanh K_z` to `-w` while leaving the in-plane measure invariant.  Therefore the thermodynamic formal free energy is even in `w`, and

\[
c_3(v)=[w^3](\phi_{3D}-\phi_{2D})=0
\]

identically as a formal function of `v`.  A fixed odd periodic layer count can have a wrapping term at order equal to that count; it is excluded before taking the layer-direction thermodynamic limit.

[COMPUTATION] As a finite exact cross-check, every one of the 15 anisotropic boxes needed through `v^8` has identically zero `w^1` and `w^3` polynomial columns.  Hence the first two requested coefficients, `[v^0]c_3` and `[v^2]c_3`, are both zero, as are all coefficients through `v^8` in the finite check.

[COMPUTATION] No `c_4` claim is made.  Unlike `c_2`, it involves connected four-spin layer correlations and logarithmic subtractions, so it is not supplied by the squared two-point-correlation machinery used here.

## Reproduction and scope

[COMPUTATION] `experiments/e43_series_structure.py` uses only exact Python `Fraction` and SymPy `Rational` arithmetic, enforces a hard 1800-second alarm on each exact linear solve, and completed all 138 budget rows without a timeout.  It writes the full frontier tables, exact minors, candidate bases, holdout residuals, interlayer arrays, and computed checks to `results/series/structure_certificates.json`.

[COMPUTATION] `.venv/bin/python tests/test_series_structure.py` independently reconstructs representative HT and LT algebraic and differential-algebraic matrices using a separate `Fraction` elimination routine, rechecks every surviving basis vector, re-derives the two-dimensional `v^8` coefficient, and cross-checks the anisotropic result.

[COMPUTATION] Nothing in this finite-prefix sweep proves algebraic transcendence, differential transcendence, non-D-finiteness at untested budgets, an analytic continuation in interlayer coupling, or an exact solution of the isotropic three-dimensional Ising model.
