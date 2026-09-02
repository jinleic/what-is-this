# Closed PGL counts for two all-distinct GRS families
## gate H-ALLDISTINCT-PGLCOUNT, run `20260902T012629Z_2eb59adb_0aa7425efe1b`

Agent `RsPe3dH2`, 2026-09-01/02 UTC. Prereg:
`prereg/H_ALLDISTINCT_PGLCOUNT_PREREG_2026-09-01.md`, path-scoped commit
`bddaa259bc2e7f175453c42e978d91159fc8ac73`, sha256
`ffb5360100e19c6d2ca28d336f91d4116a223411f80935483b3c30338f6daca7`,
byte-identical run copy `pre_statement.md` (cmp-verified). The predecessor
H-ALLDISTINCT-CROSSRATIO run was frozen and closed before this campaign was
initialized; no frozen directory was modified and no parameter-dependent
compute preceded the prereg commit. The supplied principle and formulas were
candidate inputs, not authority; §§1–5 derive them. Exact machine record:
`controls_results.json` (151 aggregate asserts, 0 failures; CPU 15.578/600 s,
wall 16.886/900 s; `nice` 15; exact GF($p$) arithmetic; thread caps 1;
bytecode disabled).

## 0. Setting and theorem

Let $X,Y\subseteq F_p$ be distinct evaluation sets for two 2-row GRS factors,
with columns $(1,x)^{\mathsf T}$, $x\in X$, and $(1,y)^{\mathsf T}$,
$y\in Y$. The frozen predecessor theorem says an all-distinct size-four
support is a circuit exactly when its paired points are the graph of a
nondegenerate Möbius map, equivalently when the pairing preserves cross-ratio.
For $M\in\mathrm{PGL}(2,p)$ put

$$D_M(X,Y)=X\cap M^{-1}(Y).$$

Here $X,Y$ contain finite points; if $M$ has a finite pole in $X$, that pole is
not in $D_M$ because $M(\mathrm{pole})=\infty\notin Y$.

**Theorem PGL-SUM.**

$$\boxed{N_{\rm all}(X,Y)=
\sum_{M\in\mathrm{PGL}(2,p)}\binom{|D_M(X,Y)|}{4}.}$$

Two closed corollaries are:

1. for $X=Y=F_p$,
   $$N_{\rm all}=p(p-1)\binom p4+p^2(p-1)\binom{p-1}4;$$
2. for $X=Y=F_p^*$ and $q=p-1$,
   $$N_{\rm all}=2q\binom q4+4q^2\binom{q-1}4+
   q^2(q-1)\binom{q-2}4.$$

Thus for arbitrary finite point sets the PGL sum is an exact evaluation
principle; the count's field dependence is exhausted by $p$, $X$, $Y$, and
the induced PGL action.

## 1. P1 — bijection between circuit supports and pairs $(M,U)$

Let $S$ be an all-distinct circuit support. Its row projection is a unique
four-subset $U\subseteq X$, and the predecessor theorem supplies a Möbius map
$M$ with

$$S=\{(x,M(x)):x\in U\},\qquad M(U)\subseteq Y.$$

The map is unique: if $M$ and $M'$ give the same support, they agree on the
four points of $U$, hence on any three distinct points. A nonidentity Möbius
map cannot fix three distinct projective points (equivalently, three distinct
source/image pairs determine one projective linear transformation), so
$M=M'$. In particular, this is not a many-to-one orbit count.

Conversely let $M\in\mathrm{PGL}(2,p)$ and
$U\in\binom{D_M(X,Y)}4$. By definition, $M(x)$ is finite and in $Y$ for every
$x\in U$, so $U$ contains no pole. A projective linear transformation is a
bijection of $\mathbb P^1(F_p)$; hence distinct $x\in U$ have distinct images.
Therefore

$$S(M,U)=\{(x,M(x)):x\in U\}$$

is all-distinct. Möbius maps preserve cross-ratio, so the frozen cross-ratio
theorem makes $S(M,U)$ a genuine circuit. If
$S(M,U)=S(M',U')$, their row projections give $U=U'$, and then $M,M'$ agree on
four (in particular three) points, so $M=M'$. Thus

$$S\longleftrightarrow(M,U),\qquad
M\in\mathrm{PGL}(2,p),\quad U\in\binom{D_M(X,Y)}4$$

is a bijection in both directions.

## 2. P2 — the general PGL sum

For fixed $M$, the bijection gives exactly
$\binom{|D_M(X,Y)|}{4}$ supports: choose any four permitted domain points.
The uniqueness proved in §1 makes the classes for different $M$ disjoint.
Summing over all maps proves Theorem PGL-SUM.

This also closes the predecessor's “simpler evaluation” loop at the correct
level of generality. No coefficient, determinant, or additional quotient is
missing: the terms depend only on the finite field and the two point sets
through the intersection cardinalities $|X\cap M^{-1}(Y)|$. Consequently all
field dependence is precisely the dependence of those cardinalities under the
finite PGL action. A further closed simplification for arbitrary subsets would
require additional structure on $X,Y$; none is claimed here.

## 3. P3 — size of PGL and the pole/zero fibers

An invertible $2\times2$ matrix is chosen by a nonzero first column
($p^2-1$ choices) and a second column outside its span ($p^2-p$ choices). Thus

$$|\mathrm{GL}(2,p)|=(p^2-1)(p^2-p).$$

Two matrices define the same projective map exactly when they differ by a
nonzero scalar; the scalar center has size $p-1$. Therefore

$$|\mathrm{PGL}(2,p)|=
\frac{(p^2-1)(p^2-p)}{p-1}=p(p^2-1).$$

Every $M$ is a bijection of $\mathbb P^1(F_p)$, so it has unique, distinct
preimages

$$r=M^{-1}(\infty)\quad\text{(pole)},\qquad
z=M^{-1}(0)\quad\text{(zero)}.$$

They are distinct because $0\ne\infty$ and $M$ is injective. Conversely fix an
ordered pair $(r,z)$ of distinct projective points and a third point
$t\notin\{r,z\}$. A projective transformation is uniquely determined by the
images of three distinct points. Requiring

$$M(r)=\infty,\qquad M(z)=0,
$$

leaves $M(t)$ free to be any projective point other than $0,\infty$, i.e. any
of the $p-1$ values in $F_p^*$. Each choice gives one map, so **every ordered
distinct pole/zero pair carries exactly $p-1$ maps**. There are
$(p+1)p$ ordered distinct pairs, recovering
$(p+1)p(p-1)=p(p^2-1)$ as a cross-check, not a separate assumption.

## 4. P4 — full-field closed form

Take $X=Y=F_p$. If the pole is $\infty$, no finite domain point maps to
$\infty$, so $D_M=F_p$ and $|D_M|=p$. Such maps are precisely affine maps.
Counting by §3: the pole is fixed at $\infty$, the zero may be any of the $p$
finite projective points, and each ordered pair carries $p-1$ maps. Hence
there are $p(p-1)$ affine maps.

If the pole is finite, exactly that one point of $F_p$ maps outside $Y$;
all other finite points have finite images. Thus $|D_M|=p-1$. There are $p$
choices of finite pole, $p$ choices of a distinct zero in the remaining
projective line (including possibly $\infty$), and $p-1$ maps per pair:
$p^2(p-1)$ maps. Substitution into PGL-SUM gives

$$N_{\rm all}(F_p,F_p)=
 p(p-1)\binom p4+p^2(p-1)\binom{p-1}4.$$

The two map classes total $p(p-1)+p^2(p-1)=p(p^2-1)$, so no maps are omitted
or counted twice.

## 5. P5 — multiplicative-group closed form

Let $S=F_p^*$ and $q=|S|=p-1$. Its complement in the projective line is
exactly $\{0,\infty\}$. Since $M$ is a bijection, a point $x\in S$ fails to
map into $S$ exactly when $M(x)$ is one of these two missing values. Their
unique preimages are the zero $z=M^{-1}(0)$ and pole $r=M^{-1}(\infty)$.
Therefore, if exactly $k$ of the distinct ordered pair $(r,z)$ lie in $S$,

$$|D_M(S,S)|=q-k,\qquad k\in\{0,1,2\}.$$

Count the ordered pole/zero pairs and multiply by the $q=p-1$ maps per pair:

- $k=0$: both entries are the two distinct points in $\{0,\infty\}$. There are
  2 ordered pairs, hence $2q$ maps.
- $k=1$: choose which of pole/zero lies in $S$ (2), its value in $S$ ($q$),
  and the other value in $\{0,\infty\}$ (2): $4q$ ordered pairs, hence
  $4q^2$ maps.
- $k=2$: choose an ordered distinct pair in $S$: $q(q-1)$ pairs, hence
  $q^2(q-1)$ maps.

The class identity is

$$2q+4q^2+q^2(q-1)=q^3+3q^2+2q
=q(q+1)(q+2)=p(p^2-1).$$

This identity is load-bearing: the tempting but wrong coefficient
$q^2(q-2)$ omits $q^2$ maps. Finally substitute the valid-domain sizes
$q,q-1,q-2$ into PGL-SUM:

$$N_{\rm all}(F_p^*,F_p^*)=
2q\binom q4+4q^2\binom{q-1}4+q^2(q-1)\binom{q-2}4.$$

This derives, rather than fits, the second closed form.

## 6. Exact in-run evidence

The instrument independently enumerated projective matrix classes by
normalizing the first nonzero matrix entry to one. For each map it computed
the action on all $p+1$ projective points, located the unique pole and zero,
and materialized every graph support. Duplicate support emission from two
$(M,U)$ pairs was a hard failure; none occurred.

### T1 — PGL structure

| $p$ | enumerated $|\mathrm{PGL}|$ | ordered pole/zero pairs | maps/pair | $k=0$ | $k=1$ | $k=2$ | class sum |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 120 | 30 | 4 | 8 | 64 | 48 | 120 |
| 7 | 336 | 56 | 6 | 12 | 144 | 180 | 336 |
| 11 | 1320 | 132 | 10 | 20 | 400 | 900 | 1320 |
| 13 | 2184 | 182 | 12 | 24 | 576 | 1584 | 2184 |

Every row exactly matches $p(p^2-1)$ and
$(2q,4q^2,q^2(q-1))$; every one of the $p(p+1)$ ordered distinct pole/zero
pairs occurred exactly $p-1$ times.

### T2 — full field, three independent routes

For $X=Y=F_p$, the product route scanned **every** product four-subset, not
only matchings, and retained genuine profile-$(4,4)$ circuits by actual
pure-tensor rank and all proper-subset checks. The cross-ratio route swept all
matchings independently. The PGL route materialized graph supports directly.

| $p$ | all product 4-subsets | all-distinct pairings | product rank | cross-ratio | PGL graphs | formula |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 12,650 | 600 | 200 | 200 | 200 | 200 |
| 7 | 211,876 | 29,400 | 5,880 | 5,880 | 5,880 | 5,880 |

All three **support sets**, not just counts, were equal. The PGL valid-domain
histograms independently expose the proof's two classes: at $p=5$, 20 maps
keep 5 points and 100 keep 4; at $p=7$, 42 keep 7 and 294 keep 6. The full
product censuses also recorded the separate crossing populations (900 and
11,025), never mixed into the all-distinct count.

### T3 — multiplicative group, three independent routes

The actual product-column route and cross-ratio route each swept every
all-distinct pairing; the PGL route enumerated all 120/336/1320 maps and every
valid four-subset.

| $p$ | all pairings | product rank | cross-ratio | PGL graphs | formula | valid-domain histogram |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 24 | 8 | 8 | 8 | 8 | $|D|=4:8,\ 3:64,\ 2:48$ |
| 7 | 5,400 | 1,080 | 1,080 | 1,080 | 1,080 | $|D|=6:12,\ 5:144,\ 4:180$ |
| 11 | 1,058,400 | 117,600 | 117,600 | 117,600 | 117,600 | $|D|=10:20,\ 9:400,\ 8:900$ |

Again all three support sets were exactly equal. This reproduces every supplied
anchor and directly verifies the map-class/domain-size correspondence used in
§5.

### T4 — fail-loud controls

- **Wrong $k=2$ coefficient REJECT (GF(7)).** Planting $q^2(q-2)=144$ gives
  class sum 300, rejected against enumerated $|\mathrm{PGL}(2,7)|=336$; the
  resulting count 1044 is rejected against 1080. This catches exactly the
  supplied hand-derivation slip.
- **Selected pole REJECT.** $M(x)=1/x$ has pole 0 and zero $\infty$. The planted
  $U=\{0,1,2,3\}$ is not contained in its valid domain
  $\{1,2,3,4,5,6\}$ and is never emitted.
- **Pure-tensor ACCEPT / non-tensor REJECT.** The identity graph has actual
  rank 3. Replacing $h(0,0)$ by $h(0,0)+h(0,1)$ while retaining labels raises
  actual rank to 4 although the labels still preserve cross-ratio; the
  mismatch rejects the corrupted instrument.
- **Duplicate-column REJECT.** A cloned factor column drops spark from 3 to 2
  and exposes dependent pair $(0,7)$.
- **Concat guard.** True 2-row by 3-row Kronecker columns have dimension 6 and
  independently checked multiplicative coordinates; block concatenation has
  dimension 5 and different coordinates.

Runtime controls asserted bytecode disabled, all thread caps 1, and nice value
15. The exact run used CPU 15.578 s and wall 16.886 s, below the fixed
600/900 s caps.

## 7. Scope, defects, and verdict

**Certified:** the general PGL-SUM principle for arbitrary finite
$X,Y\subseteq F_p$, the full-field formula, and the multiplicative-group
formula for two 2-row GRS factors. This completes the registered request for a
simpler exact evaluation in the two structured point-set families and closes
the loop on field dependence: only $p$, the point sets, and their PGL
incidence profile enter.

**Not certified:** a further symbolic simplification of the PGL sum for
arbitrary unstructured subsets; higher-row factors beyond the frozen
bi-Vandermonde rank statement; spark-2 converse branches; three or more
factors; extension-field analogues or asymptotics.

No run defect, failed assertion, parameter amendment, or prereg amendment
occurred. The wrong $q^2(q-2)$ derivation was a supplied historical slip and was
used only as a pre-registered REJECT plant; the run never adopted it. P1–P5 are
complete, all three routes agree as support sets in all five formula rows,
all structural identities and anchors pass, and every plant fires. Under the
promotion rule this supports exactly one verdict: **FROZEN-CERTIFIED**.
