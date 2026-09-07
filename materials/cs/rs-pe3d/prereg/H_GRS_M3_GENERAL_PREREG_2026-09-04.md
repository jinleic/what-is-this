# PREREG-RSPE3D-18-GRS-M3-GENERAL — gate **H-GRS-M3-GENERAL**

Written 2026-09-04 UTC by agent `RsPe3dConverse`, immediately after
H-GRS-CONVERSE was frozen, closed FROZEN-CERTIFIED (run
`20260904T020255Z_45d37a74_920d291c4cbe`), and committed. This file is
path-scoped committed before any new parameter-dependent computation. Binding
order: preregistration -> commit -> `campaign.py init` -> byte-identical
in-run copy with source commit and SHA-256 -> controls before censuses ->
analytic record -> `freeze` -> exactly one `close --verdict`. Write only under
`cs/rs-pe3d/`; do not edit any frozen run or root ledger.

**Provenance of the candidate.** Main (owner) IRC-signalled the m>=3
candidate formula (K = 2+sum(r_i-1); capacity-partition separators; pairwise
swap cross-ratios; simultaneous PGL graphs; count sum binom(k_M, K)) as a
candidate to independently re-derive, not authority. This preregistration is
this agent's independent derivation. **One correction is registered up
front:** the forward (graph => dependent) degree bound sum_i(r_i-1) = K-2 is
obtained with the **individual reduced denominators** D_j of each Moebius map
M_j = N_j/D_j (deg <= 1 each). Main's wording "common-denominator" would put
degree up to m-1 in the denominator power and break the count; the corrected
computation is Theorem M-F step 2 below. Everything else re-derives as
signalled.

## 1. Setting

m >= 3 GRS/Vandermonde factors over F_p, sizes r_i x n_i, r_i >= 2, distinct
evaluation points inside each factor, nonzero multipliers suppressed:
a_{z}=(1,z,...,z^{r_i-1})^T. Product column at a point z=(z_1,...,z_m):
a_{z_1}(x) a_{z_2} (x) ... (x) a_{z_m}, the evaluation vector of the
monomials prod_i z_i^{e_i}, 0 <= e_i < r_i, in row-major multi-index order,
dimension prod_i r_i. V = (x)_{i=1}^m F[z_i]_{<r_i}; the pairing
<Gamma, (x)_i a_{z_i}> = F_Gamma(z_1,...,z_m) transports dependence of
columns to syzygies sum_t c_t F(s_t) = 0 for all F in V (gate-15 transport
generalized; obligation (O7) below verifies the m-factor column form,
associativity, and the vec identity on structured test systems).

**All-distinct** = coordinate-injective: within every coordinate i the
z_i-values of the support are pairwise distinct. **K := 2 + sum_{i=1}^m
(r_i - 1)**. Tied strata are out of scope. Claims are field-general; in-run
verification is over GF(7), GF(11), GF(13) as registered.

## 2. Theorem M-E — independence below K (candidate proof)

**M-E.** Let S be all-distinct with k = |S| < K. Then S is independent.

**Proof.** Fix t; W = S\{s_t}, |W| = k-1 <= sum_i(r_i-1) = K-2. Assign the
victims greedily to coordinates: disjoint V_i subseteq W with sum |V_i| =
k-1 and |V_i| <= r_i-1. Then g_t = prod_i prod_{u in V_i}(z_i - z_i(u)) has
multidegree |V_i| <= r_i-1 in each coordinate, so g_t in V; it kills every
s_u, u != t (u is killed in the coordinate whose V_i contains it) and
g_t(s_t) != 0 by all-distinctness. The {g_t} form a Lagrange system:
evaluation V -> F^k is surjective, rank = k, columns independent. □

Corollary: every proper subset of a size-K all-distinct set is independent;
any syzygy on such a set has full support; a dependent size-K all-distinct
set is a circuit.

## 3. Theorem M-C — converse at K (candidate proof)

**M-C.** Let S be all-distinct, |S| = K >= 5 (m >= 3, r_i >= 2), dependent.
Then there is a **unique** tuple (M_2,...,M_m) of nondegenerate Moebius maps
with z_j(s) = M_j(z_1(s)) for every s in S, j = 2..m.

**Proof.**
(i) Full support (corollary of M-E).
(ii) **Split relations.** For distinct t,t' let O = S\{s_t,s_t'}, |O| = K-2 =
sum_i(r_i-1). For every ordered partition O = (x)_i U_i with |U_i| = r_i-1:
F = prod_i prod_{u in U_i}(z_i - z_i(u)) in V kills O and is nonzero at
s_t, s_t'; the syzygy gives c_t A(t) + c_t' A(t') = 0 where
A(w) = prod_i prod_{u in U_i}(z_i(w) - z_i(u)).
(iii) **Two-coordinate swap.** Fix coordinates i,j, distinct u_1,u_2 in O.
The two partitions identical except U_i^(1) = W_i (+){u_1}, U_j^(1) = W_j
(+){u_2} versus U_i^(2) = W_i (+){u_2}, U_j^(2) = W_j (+){u_1}, with
|W_i| = r_i-2, |W_j| = r_j-2 and all other U_l fixed (l != i,j), exist
exactly: |O\{u_1,u_2}| = K-4 = (r_i-2)+(r_j-2)+sum_{l!=i,j}(r_l-1).
Dividing the two syzygy equations cancels every common factor and leaves
cr(z_i(t),z_i(t');z_i(u_1),z_i(u_2)) = cr(z_j(t),z_j(t');z_j(u_1),z_j(u_2))
for every ordered distinct quadruple and every coordinate pair (i,j).
(iv) **Reconstruction.** For each j >= 2, coordinates 1 and j satisfy the
cross-ratio equality on every quadruple; the 3-seed reconstruction of
H-GRS-CONVERSE Theorem C step (iv) (M_j fixing the three seed images;
fourth-slot injectivity on the z_j-line) gives z_j = M_j(z_1) on S, M_j
unique and nondegenerate (three distinct images by all-distinctness), no
pole on the z_1-values. □

**Reducible / repeated-component clause.** As in the two-factor gate: no
factorization or component hypothesis is used anywhere; degenerate Delta=0
relations force constant coordinates (excluded by all-distinctness), and
multiple PGL representatives of one graph are a counting convention collapsed
by per-coordinate canonical normalization and in-run set semantics.

## 4. Theorem M-F — forward at K (candidate proof, corrected denominator
device)

**M-F.** Let S be all-distinct, |S| = K, with z_j(s) = M_j(z_1(s)) for
nondegenerate M_j = N_j/D_j (coprime, deg <= 1, D_j(z) != 0 on the
z_1-values of S). Then S is dependent — hence a circuit by M-E.

**Proof.** 1. For F = sum_e f_e prod_i z_i^{e_i} in V (e_i < r_i):
F(z, M_2(z),...,M_m(z)) = P_F(z)/W(z) with W(z) = prod_{j>=2} D_j(z)^{r_j-1}
and
P_F(z) = sum_e f_e z^{e_1} prod_{j>=2} N_j(z)^{e_j} D_j(z)^{r_j-1-e_j}.
2. deg P_F <= e_1 + sum_{j>=2}[e_j deg N_j + (r_j-1-e_j) deg D_j]
<= e_1 + sum_{j>=2}(r_j-1) <= (r_1-1) + sum_{j>=2}(r_j-1) = K-2,
using deg N_j = deg D_j = 1 bound (PGL2 = degree-<=1 rational maps) and
**individual** denominators D_j. So the numerator space {P_F} lies in
F[z]_{<=K-2}, dimension <= K-1 < K.
3. The K evaluation functionals of that space at the K distinct z_1-values
are dependent: there is c != 0 with sum_t c_t P(z_t) = 0 for all P in the
space. Then sum_t c'_t F(s_t) = 0 for all F in V, where
c'_t = c_t W(z_t) != 0: a full-support syzygy. □

## 5. Count (candidate formula)

With M-F and M-C, the all-distinct size-K circuits of A_1 (x) ... (x) A_m are
EXACTLY the simultaneous-graph supports, and
$$N = sum_{(M_2,...,M_m) in PGL(2,p)^{m-1}} binom(k_M, K),$$
k_M = #{z in X_1 : M_j(z) in X_j for all j >= 2}. The parameterization is
injective for K >= 3 (three points determine each M_j; here K >= 5), with
per-coordinate canonical normalization (first nonzero of the 4-tuple to 1).

## 6. In-run obligations

On **every** enumerated dependent all-distinct K-set: **(O1)** rank K-1 by
two independent rank algorithms, nullspace dimension exactly 1, full support;
**(O2)** split relations for all pairs and ALL partitions (multinomial count)
with |U_i| = r_i-1; **(O3)** cross-ratio equality for all coordinate pairs
(i,j) and all ordered distinct quadruples; **(O4)** reconstruction of
(M_2,...,M_m) from the seed triple (indices 0,1,2), each M_j nondegenerate,
pole-free, fitting all points, and the independent seed (1,2,3) giving the
same canonical tuple; **(O7)** m-factor transport: column form (column at z
equals the multi-index evaluation vector, Kronecker associativity), vec
identity on structured test systems, multipliers rank-free. **(O5)**
toolkit identities (Moebius invariance of cr, fourth-slot injectivity) per
prime. **(O6)** m-factor separators for sampled sets at sizes 2 <= k <= K-1
only — the construction exists exactly there (k-1 <= sum(r_i-1)); size-K
sets are NOT separator-verified (dependent sets exist at K; the all-but-two
split relations of (O2) are the K-size object, verified on every dependent
set).

## 7. Registered anchors (exhaustive; owner-independent prediction checks)

Grids X_i = {1..n} mod p; exhaustive all-distinct K-set sweeps; direct count
must equal the PGL-tuple sum (eq-b) AND the PGL-graph support set (eq-a) at
every registered cell; every dependent set passes (O1)-(O4).

1. (2,2,2), K=5, n=5: all three primes 13, 11, 7; 14400 sets per prime;
   sizes 2..4 all-distinct dependent count = 0 (exhaustive) at all primes.
2. (2,2,2), K=5, n=6: 3,110,400 sets; primes 13 and 11.
3. (2,2,3), K=6, n=6: 518,400 sets; primes 13 and 11; size 5 dependent
   count = 0 (exhaustive) at both primes.
4. Regression anchor: the two-factor (2,4) n=6 GF(13) size-6 count 2 is
   reproduced by the shared library as a cross-check of the m=2 special
   case.

## 8. Sampled verification (no count pins; constructive + random)

- (3,3,3), K=8, n=8, primes 11, 13: 200 random all-distinct 8-sets per prime
  (rank 8 asserted = independent); PLUS constructive dependent sets from
  enumerated simultaneous PGL triples (at least 8 per prime): full
  (O1)-(O4).
- (2,2,2,2), K=6, n=6, prime 11: 200 random all-distinct 6-sets (rank 6
  asserted); constructive simultaneous-graph sets: full (O1)-(O4).
- Separator samples (O6) at (2,2,2) GF(7) n=5, (2,2,3) GF(11) n=6,
  (3,3,3) GF(11) n=8, (2,2,2,2) GF(11) n=6 for sizes 2..K-1.

## 9. Control battery (before any census; every ACCEPT paired with a REJECT)

1. **m-factor ACCEPT**: a planted simultaneous-graph 5-set at (2,2,2) GF(13)
   (e.g. z_2 = 2 z_1^{-1}, z_3 = 3 z_1 + 1 on a common 5-subset of
   domains): circuit; (O1)-(O4) pass end to end.
2. **Non-graph REJECT**: an all-distinct 5-set on {1..5}^3 off every
   (M_2,M_3) tuple: rank 5 exactly; no nondegenerate pair of relations fits
   (kernel of the 2-relation system contains no nondegenerate pair).
3. **Sub-K plant REJECT**: ACCEPT 5-set minus one point: rank 4; corrupted
   variant: rank 4; plus the exhaustive sub-K zeros of anchor 1.
4. **Concat-vs-Kronecker guard (m=3)**: true column dimension 8 and
   row-major multi-index coordinates at (2,2,2); block-concatenation
   impostor rejected on dimension and coordinates.
5. **Ambient/dimension plant**: full 3x3x3 evaluation space at (2,2,2) GF(13)
   has rank exactly 8; a planted 9-cell set has rank exactly 8.
6. **Wrong-seed REJECT**: reconstructions from a corrupted seed triple fail
   on the true point set (or are degenerate/absent), for (M_2,M_3)
   simultaneously.
7. **Uniqueness control**: two independent seed triples give the same
   canonical (M_2,M_3) on every ACCEPT set (also enforced in (O4)).
8. **m=2 regression**: (2,4) n=6 GF(13) size-6 anchor count 2 (anchor 4).

## 10. Stage plan for `run_grs_m3_general.py`

- **A** m-factor library: multi-index Kronecker columns, dual ranks,
  nullspace/syzygies, per-coordinate canonical PGL, cr utilities, tuple
  enumeration, transport obligations (O7), toolkit identities (O5).
- **F** control battery (before censuses).
- **B** PGL layer: |PGL(2,p)| = p(p^2-1); tuple domains k_M.
- **C** anchors (section 7), primes in order 13, 11, 7.
- **D** sampled verification (section 8).
- **E** separators (O6).
- **G** wrap up: budgets, results, assertion count.

## 11. Arithmetic, runtime discipline, budgets, defects

Exact GF(p) integer arithmetic only; no floats, no NumPy/sympy;
sys.dont_write_bytecode=True first; launched with PYTHONDONTWRITEBYTECODE=1
and all five thread caps = 1; nice -n 10 asserted; RLIMIT_CPU soft = hard =
5400 before the first long stage. **Hard cap 5400 CPU s / 6000 wall s**,
checked per family; no adaptive extension. Any defect: disclosed in
defect_log.json, broken output preserved under a distinct name, affected
family rerun from its beginning (complete-rerun precedent). Prereg
amendments, if any, follow the A1/A2 transparency rule: original hash
preserved, both hashes in provenance, never claimed unamended.

## 12. Promotion rule and verdict

**FROZEN-CERTIFIED** iff: M-E, M-C, M-F are registered as proved with the
in-run obligations verified as specified; every anchor of section 7
reproduced exactly (eq-a set equality and eq-b counts); section 8 sampled
verification clean (random sets independent; constructive graph sets
dependent with full pipeline); the full control battery passes with every
planted REJECT firing. Otherwise **FROZEN-INCONCLUSIVE**, naming the first
failed proof step or control while certifying inside the record every part
that did pass. Nothing incomplete is promoted; any proof gap is named at the
exact failing step.
