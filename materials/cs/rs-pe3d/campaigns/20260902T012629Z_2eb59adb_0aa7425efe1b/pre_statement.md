# PREREG-RSPE3D-8-ALLDISTINCT-PGLCOUNT — gate **H-ALLDISTINCT-PGLCOUNT**

Written 2026-09-01 (local date), agent `RsPe3dH2`, after predecessor
H-ALLDISTINCT-CROSSRATIO was frozen, closed FROZEN-CERTIFIED, independently
verified by Main, and committed. This file is path-scoped committed before any
parameter-dependent compute for this run. Binding order: preregistration →
`git commit --only` → `campaign.py init --gate H-ALLDISTINCT-PGLCOUNT` →
byte-identical prereg copy in the minted run directory with source commit and
SHA-256 → exact controls and censuses → proof record → `freeze` →
`close --verdict` with exactly one verdict. No frozen campaign directory is
modified.

The principle, specializations, and numeric anchors below were supplied by
Main as derived and externally verified, not as authority. This run must
re-derive every counting step and independently reproduce every anchor.

## 1. Fixed theorem and notation

Let $X,Y\subseteq F_p$ be distinct evaluation sets for two 2-row GRS factors,
with columns $(1,x)^{\mathsf T}$, $x\in X$, and $(1,y)^{\mathsf T}$,
$y\in Y$. The predecessor theorem proves that an all-distinct size-four
support is a circuit iff its four paired points are the graph of a unique
Möbius transformation $M\in\mathrm{PGL}(2,p)$, equivalently iff the pairing
preserves cross-ratio.

For $M\in\mathrm{PGL}(2,p)$ define

$$D_M(X,Y):=X\cap M^{-1}(Y),$$

where $Y$ contains only finite points, so a finite pole of $M$ is automatically
excluded. The target theorem is

$$\boxed{\quad N_{\rm all}(X,Y)=
\sum_{M\in\mathrm{PGL}(2,p)}\binom{|D_M(X,Y)|}{4}.\quad}$$

The support/pair bijection is

$$S\longleftrightarrow(M,U),\qquad M\in\mathrm{PGL}(2,p),\quad
U\in\binom{D_M(X,Y)}4,$$

with $S=\{(x,M(x)):x\in U\}$. Uniqueness of $M$ is load-bearing and must be
proved, not inferred from the count.

Two closed specializations are fixed:

1. **Full field $X=Y=F_p$.**
   $$N_{\rm all}=p(p-1)\binom p4+p^2(p-1)\binom{p-1}4.$$
2. **Multiplicative group $X=Y=F_p^*$, $q=p-1$.**
   $$N_{\rm all}=2q\binom q4+4q^2\binom{q-1}4+
   q^2(q-1)\binom{q-2}4.$$

## 2. Proof obligations

- **P1 — circuit supports biject with $(M,U)$.** From the predecessor
  cross-ratio theorem, prove every all-distinct circuit support has a Möbius
  graph. Prove the map is unique because two Möbius transformations agreeing
  on three distinct projective points are equal. Conversely, for every
  $M$ and four-subset $U\subseteq D_M(X,Y)$, prove $M$ has no pole on $U$,
  is injective there, gives an all-distinct support, preserves cross-ratio,
  and therefore gives a circuit. Prove two pairs cannot yield the same support.
- **P2 — the general principle.** Count the four-subsets independently for
  each unique $M$, obtaining the displayed sum. Explain precisely why its
  field dependence is only through $p$, $X$, $Y$, and the induced PGL action.
- **P3 — PGL size and pole/zero fibers.** Prove
  $|\mathrm{PGL}(2,p)|=p(p^2-1)$. Every $M$ has distinct projective preimages
  `pole` $=M^{-1}(\infty)$ and `zero` $=M^{-1}(0)$. Prove every ordered pair
  of distinct projective points carries exactly $p-1$ maps (fix a third domain
  point and choose its image in $F_p^*$, or an equivalent exact argument).
- **P4 — full-field specialization.** A map with pole $\infty$ keeps all $p$
  finite points; there are $p(p-1)$ such affine maps. A map with finite pole
  keeps exactly $p-1$ finite points; there are $p^2(p-1)$ such maps. Derive
  formula 1 without overcounting.
- **P5 — multiplicative specialization.** Put $S=F_p^*$, $q=p-1$. Since
  $\mathbb P^1(F_p)\setminus S=\{0,\infty\}$, a point of $S$ fails to map into
  $S$ iff it is the pole or the zero. If exactly $k\in\{0,1,2\}$ of the ordered
  pair `(pole,zero)` lie in $S$, prove $|D_M(S,S)|=q-k$. Derive the map-class
  sizes
  $$2q\quad(k=0),\qquad4q^2\quad(k=1),\qquad q^2(q-1)\quad(k=2),$$
  their identity
  $$2q+4q^2+q^2(q-1)=p(p^2-1),$$
  and formula 2. Explicitly reject the known slip $q^2(q-2)$ for the $k=2$
  class.

**Promotion rule.** FROZEN-CERTIFIED only if P1–P5 are complete with no gap,
every listed numeric anchor is reproduced by three independent routes, all
feasible support sets agree exactly, the class-size identity passes, and every
planted REJECT fires. Otherwise close FROZEN-INCONCLUSIVE naming the first
proof/control gap.

## 3. Exact control matrix

All mathematical decisions use Python integers modulo $p$, exact
Gauss–Jordan rank/determinant, exact projective PGL enumeration, stdlib only.
PGL matrices are enumerated modulo scalar by normalizing the first nonzero
entry to one; enumeration size and uniqueness are asserted.

### T1 — PGL structural census

For $p\in\{5,7,11,13\}$ enumerate every canonical nondegenerate $2\times2$
projective matrix. Assert exactly $p(p^2-1)$ maps. Compute each map's unique
pole and zero by its action on $\mathbb P^1(F_p)$; assert they are distinct,
every ordered distinct `(pole,zero)` pair has exactly $p-1$ maps, and the
multiplicative-group class sizes are exactly
$2q,4q^2,q^2(q-1)$. Assert their sum equals the enumerated PGL size.

### T2 — full-field formulas by three routes

For $X=Y=F_p$ and $p\in\{5,7\}$:

1. **Product census:** build the actual pure-tensor columns
   $(1,y,x,xy)$ and scan every $\binom{p^2}{4}$ four-subset, retaining genuine
   profile-$(4,4)$ circuits by exact rank and all proper-subset independence.
2. **Cross-ratio route:** enumerate all
   $\binom p4^2 4!$ all-distinct pairings and retain cross-ratio equality.
3. **Direct PGL route:** for every canonical $M$, compute $D_M(F_p,F_p)$,
   enumerate all four-subsets, and materialize their graph supports; assert no
   support is emitted twice.

Assert support-set equality among all three routes and formula 1. Pin counts:
$p=5\mapsto200$, $p=7\mapsto5880$.

### T3 — multiplicative-group formulas by three routes

For $X=Y=F_p^*$ and $p\in\{5,7,11\}$:

1. exhaust every all-distinct support and retain actual product-column
   circuits by exact rank/proper-subset checks;
2. retain cross-ratio equality independently;
3. enumerate every PGL map (120/336/1320 maps), materialize every valid graph
   support, and assert uniqueness.

Assert all three support sets equal formula 2. Pin counts:
$p=5\mapsto8$, $p=7\mapsto1080$, $p=11\mapsto117600$.

### T4 — controls in both directions

1. **Wrong-class-coefficient REJECT.** At $p=7$, replace the correct $k=2$
   class $q^2(q-1)$ by $q^2(q-2)$. Assert the planted class sum is 300 rather
   than $|\mathrm{PGL}(2,7)|=336$, and the planted formula gives 1044 rather
   than 1080.
2. **Selected-pole REJECT.** Use $M(x)=1/x$ over GF(7), whose finite pole is
   $0$. Plant $U=\{0,1,2,3\}\subset F_7$ and assert it is rejected because
   $M(0)=\infty\notin F_7$; the direct enumerator must never emit that pair.
3. **Pure-tensor ACCEPT / non-tensor REJECT.** The identity graph on four
   finite points is a genuine pure-tensor circuit. Replace actual $h(0,0)$ by
   $h(0,0)+h(0,1)$ without changing labels and assert actual rank rises to 4
   while the label/PGL predicate still accepts, rejecting the corrupted
   instrument.
4. **Duplicated-column REJECT.** Clone one factor column; assert spark drops
   from 3 to 2 and a dependent pair appears.
5. **Concat-vs-Kronecker guard.** On a 2-row by 3-row factor pair, assert true
   row-major Kronecker coordinates and dimension 6; planted block
   concatenation has dimension 5 and different coordinates.

## 4. Arithmetic, budgets, defects, and kill criteria

- Exact mathematical arithmetic only: Python integers modulo $p$; inverses by
  `pow(x,-1,p)`; no floating-point mathematical decisions, NumPy, or external
  algebra package. Timing does not enter a theorem decision.
- `sys.dont_write_bytecode=True`; no frozen import; execute with
  `PYTHONDONTWRITEBYTECODE=1`.
- One process at `nice -n 15`, with
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
- Hard cap 600 CPU seconds and 900 wall seconds, checked during every exhaustive
  family. No adaptive parameter extension.
- Any defect is disclosed, broken output preserved, and every affected family
  rerun from its beginning. No assertion is weakened or silently repinned.
- Kill on any P1–P5 gap; incomplete PGL/support enumeration; support-set
  mismatch; formula/anchor mismatch; wrong class identity accepted; selected
  pole accepted; tensor/duplicate/concat guard failure; runtime or budget
  violation; or unresolved defect. Close exactly one verdict.

## 5. Fixed scope and nonclaims

Certified if promoted: the general finite-set PGL sum principle for two 2-row
GRS factors and the two stated full-field/multiplicative-group closed forms.
Not certified: a closed simplification for arbitrary subsets $X,Y$ beyond the
PGL sum; higher-row factors beyond the predecessor's bi-Vandermonde rank
statement; spark-2 branches; three or more factors; asymptotics; extension
fields not explicitly represented by the prime-field controls (the group-
action proof is stated for $F_p$ as preregistered).
