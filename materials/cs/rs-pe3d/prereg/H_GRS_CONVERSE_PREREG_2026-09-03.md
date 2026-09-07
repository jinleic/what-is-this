# PREREG-RSPE3D-17-GRS-CONVERSE — gate **H-GRS-CONVERSE**

Written 2026-09-03 UTC by agent `RsPe3dConverse`, after H-44ROW-LAYER was
frozen, closed FROZEN-CERTIFIED (run `20260902T135324Z_516630ab_ed03d008f0f6`),
and after H-GRS-GENERAL-SIZE (run `20260902T113614Z_a999d1b0_5daa598b6a04`,
FROZEN-CERTIFIED) recorded C1's converse and general all-distinct emptiness
below $r_A+r_B$ as the two surviving scoped open leads of the general-size
theorem. This file is path-scoped committed before any new parameter-dependent
computation. Binding order: preregistration -> commit -> `campaign.py init` ->
byte-identical in-run copy with source commit and SHA-256 -> controls before
censuses -> analytic record -> `freeze` -> exactly one `close --verdict`. Write
only under `cs/rs-pe3d/`; do not edit any frozen run or root ledger.

The two candidate claims below (Theorems E and C) were derived analytically
before registration, in full, by a syzygy/separator route that avoids the
slice-device boundary recorded in PREREG-RSPE3D-15 §3; they are not authority.
The run must re-verify every proof obligation line by line and confirm every
positive claim in exact arithmetic, each with a known-true ACCEPT and a planted
REJECT. Counts outside proved constructions remain census data.

## 1. Setting (identical to PREREG-RSPE3D-15 §1)

Two GRS/Vandermonde factors over $F_p$, sizes $r_A\times n_A$, $r_B\times n_B$,
distinct evaluation points, nonzero multipliers suppressed:
$a_x=(1,x,\ldots,x^{r_A-1})^{\mathsf T}$, $b_y=(1,y,\ldots,y^{r_B-1})^{\mathsf T}$.
Product column at $(x,y)$: $h(x,y)=a_x\otimes b_y=(x^iy^j)_{0\le i<r_A,0\le
j<r_B}$, row-major Kronecker, dimension $r_Ar_B$. For a support $S$ of grid
cells, $E_S$ is the $r_Ar_B\times|S|$ evaluation matrix with columns
$h(x,y)$. $V=F[x]_{<r_A}\otimes F[y]_{<r_B}$ = polynomials of bidegree
$\deg_x<r_A$, $\deg_y<r_B$; the pairing $\langle\Gamma,h(x,y)\rangle=
F_\Gamma(x,y)$ identifies column dependence of $A\otimes B$ on $S$ with
syzygies $\sum_{t}c_tF(s_t)=0\ \forall F\in V$. **All-distinct** =
coordinate-injective ($x$'s pairwise distinct and $y$'s pairwise distinct).
Circuit predicate (C) as in PREREG-15. $r_A,r_B\ge2$ throughout; the tied
strata are out of scope exactly as before.

## 2. Theorem E — emptiness below $r_A+r_B$ (candidate proof, to be verified)

**Theorem E.** Over any field, let $S$ be all-distinct with $k=|S|<
r_A+r_B$. Then the $k$ columns $h(s_t)$, $s_t\in S$, are linearly independent;
equivalently $S$ imposes $k$ independent conditions on $V$; equivalently no
all-distinct circuit has size $<r_A+r_B$.

**Proof (registered route: explicit separators).** For each $t$ let
$W_t=S\setminus\{s_t\}$, $|W_t|=k-1\le(r_A-1)+(r_B-1)$. Choose
$U_t\subseteq W_t$ with $|U_t|=\min(k-1,r_A-1)$ and $V_t=W_t\setminus U_t$;
then $|V_t|=k-1-|U_t|\le r_B-1$ (if $k-1\le r_A-1$ then $V_t=\emptyset$;
else $|V_t|=k-r_A\le r_B-1$). Define
$$g_t(x,y)=\prod_{u\in U_t}(x-x_u)\prod_{w\in V_t}(y-y_w).$$
Bidegree $(|U_t|,|V_t|)\le(r_A-1,r_B-1)$, so $g_t\in V$. For $u\ne t$: if
$u\in U_t$ the factor $(x_u-x_u)$ kills $g_t(s_u)$; otherwise $u\in V_t$ and
the factor $(y_u-y_u)$ kills it. And $g_t(s_t)=\prod_{u\in U_t}(x_t-x_u)
\prod_{w\in V_t}(y_t-y_w)\ne0$ by all-distinctness. Hence the evaluation map
$V\to F^k$, $F\mapsto(F(s_t))_t$, is surjective, so $\operatorname{rank}E_S=k$
and the columns are independent. $\square$

The proof is constructive and field-general; it uses no sweep, no primality,
no degree hypothesis beyond $k<r_A+r_B$. In-run obligations (O6) below verify
the construction line by line; exhaustive sub-$k$ emptiness censuses at the
anchors corroborate the statement as a whole.

## 3. Theorem C — C1 converse in general (candidate proof, to be verified)

**Theorem C.** Over any field, let $S=\{s_t=(x_t,y_t)\}_{t=1}^k$ be
all-distinct with $k=r_A+r_B\ge4$ and dependent (i.e.
$\operatorname{rank}E_S\le k-1$). Then there is a **unique** nondegenerate
Möbius map $M$ with $y_t=M(x_t)$ for every $t$; i.e. $S$ lies on the graph of
a single Möbius map. Nondegenerate means: with the gate-15 convention
$c_3xy+c_1y+c_2x+c_0=0$, $\Delta=c_0c_3-c_1c_2\ne0$ and $M(x)=-(c_2x+c_0)/
(c_3x+c_1)$ is finite at every $x_t$.

**Proof (registered route: two-point separators force all cross-ratios).**

(i) *Full support.* By Theorem E every proper subset of $S$ (size $<k$) is
independent. A nonzero syzygy $c$ (exists by dependence) therefore has
$c_t\ne0$ for every $t$: otherwise its support would be a dependent proper
subset.

(ii) *Split relations.* Fix distinct $t,t'$ and let $O=S\setminus\{s_t,s_{t'}\}$,
$|O|=k-2=r_A+r_B-2$. For every $U\subseteq O$ with $|U|=r_A-1$ (exists since
$|O|\ge r_A-1$), put $V=O\setminus U$ ($|V|=r_B-1$) and
$F_{t,t',U}(x,y)=\prod_{u\in U}(x-x_u)\prod_{w\in V}(y-y_w)\in V$ (bidegree
$(r_A-1,r_B-1)$). Then $F(s_u)=0$ for $u\in O$ while $F(s_t),F(s_{t'})\ne0$,
so the syzygy gives, with $A_U=\prod_{u\in U}(x_t-x_u)$,
$A'_U=\prod_{u\in U}(x_{t'}-x_u)$, $B_V,B'_V$ the $y$-analogues:
$$c_tA_UB_V=-c_{t'}A'_UB'_V.\tag{$\ast$}$$

(iii) *Cross-ratio equality.* Let $u_1\ne u_2\in O$ and
$W\subseteq O\setminus\{u_1,u_2\}$, $|W|=r_A-2$ (exists since
$|O|-2=r_A+r_B-4\ge r_A-2$). Apply $(\ast)$ to $U_2=W\cup\{u_1\}$ and
$U_3=W\cup\{u_2\}$ and divide (all factors nonzero by (i) and
all-distinctness). Cancellation of the common $W$/$V_1$ factors leaves
$$\frac{(x_t-x_{u_1})(x_{t'}-x_{u_2})}{(x_{t'}-x_{u_1})(x_t-x_{u_2})}
=\frac{(y_t-y_{u_1})(y_{t'}-y_{u_2})}{(y_{t'}-y_{u_1})(y_t-y_{u_2})},$$
i.e. $\operatorname{cr}(x_t,x_{t'};x_{u_1},x_{u_2})=
\operatorname{cr}(y_t,y_{t'};y_{u_1},y_{u_2})$ for **every ordered tuple of
distinct indices** $(t,t',u_1,u_2)$, where
$\operatorname{cr}(a,b;c,d)=\frac{(a-c)(b-d)}{(b-c)(a-d)}$.

(iv) *Reconstruction lemma.* Pick distinct $p,q,r$; let $M$ be the unique
Möbius map with $M(x_p)=y_p$, $M(x_q)=y_q$, $M(x_r)=y_r$ (exists and is
nondegenerate: a degenerate $\Delta=0$ relation is a constant map and cannot
take three distinct values). For any $s\notin\{p,q,r\}$, (iii) and Möbius
invariance of $\operatorname{cr}$ give
$\operatorname{cr}(x_p,x_q;x_r,x_s)=\operatorname{cr}(y_p,y_q;y_r,y_s)
=\operatorname{cr}(x_p,x_q;x_r,M(x_s))$; since $z\mapsto\operatorname{cr}
(x_p,x_q;x_r,z)$ is a bijection of $\mathbb P^1$ (fractional linear in $z$,
$x_p\ne x_q$), $x_s=M(x_s)$, i.e. $y_s=M(x_s)$ (finite: $y_s\in F$). The
three seed indices are covered by construction. Uniqueness: two Möbius maps
agreeing on three distinct points coincide. $\square$

**Reducible / repeated-component clause.** The proof never factors a vanishing
polynomial: no irreducibility, squarefreeness, or component hypothesis is used
anywhere, so there are no exceptional reducible or repeated-component cases to
handle — the statement is unconditional on the all-distinct stratum. The
alternative slice route's graph relation $q_1(x)-q(x)y$ and its possible
redundant factorisations correspond to (a) the same graph counted by several
PGL representatives (a counting convention, collapsed by the canonical
normalization "first nonzero of $(c_3,c_1,c_2,c_0)$ to 1" — the in-run set
semantics asserts this collapse), and (b) degenerate $\Delta=0$ relations,
which live in the tied strata (they force equal $y$-coordinates) and are
excluded here by all-distinctness; the reconstructed $M$ is automatically
nondegenerate (three distinct seed images). The slice-device boundary recorded
in PREREG-15 §3 (composition degree count, closing only
$k>\deg R\deg R'+1$) is **not used** and is superseded by this route.

**Forward direction.** Not re-proved here: Möbius-graph sets at size
$r_A+r_B$ are dependent of rank $k-1$ per gate 15 (C1 forward, certified);
the run reproduces it at the anchors (every enumerated dependent set has rank
exactly $k-1$ with all deletions independent) and re-uses the certified
per-prime PGL support-set equality (eq-a) and PGL sum (eq-b) as the bridge to
the frozen record.

## 4. In-run proof obligations (line-by-line verification)

On **every** all-distinct dependent $k$-set found in any sweep (all of them,
not a sample), and before any count is recorded from it:

- **(O1)** rank $E_S=k-1$ by two independent rank algorithms; nullspace
  dimension exactly $1$; the unique syzygy has full support (O1 = proof step
  (i) + Theorem E's corollary that nullity is 1).
- **(O2)** the split relations $(\ast)$ for **all** pairs $\{t,t'\}$ and
  **all** $U\subseteq O$, $|U|=r_A-1$ (proof step (ii)).
- **(O3)** the cross-ratio equalities for **all** ordered distinct quadruples
  (proof step (iii)).
- **(O4)** reconstruction: solve the seed system at indices $0,1,2$; assert a
  nondegenerate solution exists ($\Delta\ne0$), assert every $x_t$ is off its
  pole, assert $y_t=M(x_t)$ for all $t$; assert a second, independent seed
  triple yields the same canonical map (single-map clause, step (iv)).
- **(O5)** toolkit identities per prime: Möbius invariance of $\operatorname{cr}$
  and bijectivity of the fourth slot, on deterministic test tuples.
- **(O6)** Theorem E separators: for sampled all-distinct sets at every
  registered configuration and size (including $k=r_A+r_B$ itself), build the
  $g_t$ of §2, assert the bidegree bounds $|U_t|\le r_A-1$, $|V_t|\le r_B-1$,
  the kill pattern $g_t(s_u)=0$ ($u\ne t$), the nonzero diagonal
  $g_t(s_t)\ne0$, and the resulting exact rank $k$.

A dependent set failing any obligation freezes the affected family, preserves
broken output, and reruns that family from the beginning (§9).

## 5. Certified anchors (exact reproduction required; gate-15 §6 and gate-12)

Points $\{1..n\}\bmod p$; exhaustive sweeps over all-distinct sets; primes in
order 13, 11, 7. "Circuits" = dependent with all $k$ deletions of rank $k-1$.

1. $(r_A,r_B)=(2,4)$, $k=6$, $n=6$: dependent $=$ circuits $=$ PGL sum $=$
   **2 / 4 / 12** at GF(13)/GF(11)/GF(7); eq-a support-set equality; sub-$k$
   sizes 2–5 all-distinct dependent count $=0$ at $n=6$, all three primes.
2. $(2,4)$, $k=6$, $n=7$: circuits $=$ **42 / 72 / 588**; eq-a; size-5
   dependent count $=0$ at $n=7$, all three primes.
3. $(3,3)$, $k=6$, $n=6$: circuits $=$ **2 / 4 / 12**; eq-a; sizes 2–5 empty
   at $n=6$, all three primes.
4. $(3,3)$, $k=6$, $n=7$: circuits $=$ **42 / 72 / 588**; eq-a; size-5 empty
   at $n=7$, all three primes.
5. $(3,4)$, $k=7$, $n=7$: circuits $=$ **2 / 2 / 42**; eq-a; size-6
   (sub-$k$) dependent count $=0$, all three primes.
6. $(4,4)$, $k=8$, $n=8$ (full matchings): circuits $=$ **4 / 8** at
   GF(13)/GF(11), reproducing gate-12's certified all-distinct Möbius-8
   bijection counts; eq-a; eq-b against $\sum_M\binom{k_M}{8}$.

Obligations (O1)–(O4) run on **every** dependent set enumerated in 1–6.

## 6. New general-coverage instances (census; internal two-route checks only)

The proof is general, so these instances are corroboration, not new pinned
claims; their counts are cross-checked two ways (direct exhaustive rank sweep
vs PGL-family enumeration $\sum_M\binom{k_M}{k}$ and eq-a set equality) and
every dependent set passes (O1)–(O4). $n=k$ full matchings
(permutations of $n$), points $\{1..n\}\bmod p$:

- $(2,5)$, $k=7$, $n=7$, primes 7, 11, 13; plus exhaustive sub-$k$ size-6
  emptiness ($=0$ expected) at all three primes; sizes 2–5 sampled.
- $(3,5)$, $k=8$, $n=8$, primes 11, 13; sizes 2–7 sampled (300 per size,
  seeded).
- $(2,6)$, $k=8$, $n=8$, primes 11, 13; sizes 2–7 sampled.
- $(4,5)$, $k=9$, $n=9$, primes 11, 13; sizes 2–8 sampled.

These are exactly the unbalanced regimes where the superseded slice-device
boundary failed to close; a converse failure there would refute Theorem C, so
any nonzero population of non-Möbius dependent sets, any (O1)–(O4) failure, or
any eq-a/eq-b mismatch aborts the gate (§9).

## 7. Control battery (before any census; every ACCEPT paired with a REJECT)

1. **Converse ACCEPT on a planted known-true graph.** A size-6 all-distinct
   set on one PGL graph at (2,4) GF(13) (e.g. $M(x)=2x^{-1}$ on a 6-subset of
   its domain): dependent, rank 5, all deletions rank 5, then (O1)–(O4) pass
   end to end.
2. **Non-Möbius REJECT (the converse's contrapositive, planted).** An
   all-distinct $k$-set chosen off every PGL graph at (2,4) GF(13) $n=6$ and
   at one more configuration: rank $=k$ exactly (square case: $\det\ne0$);
   no nondegenerate Möbius relation exists (solution space of
   $c_3xy+c_1y+c_2x+c_0=0$ contains no $\Delta\ne0$ vector).
3. **Sub-$k$ plant REJECT ("a dependent set below $r_A+r_B$ must not
   exist").** Delete one point from the ACCEPT circuit (size $k-1$):
   independent, rank $k-1$; plus corrupted variants (one $y$ moved to an
   unused value): independent. Exhaustive sub-$k$ emptiness at the anchors
   (§5) is the global form of this control.
4. **Concat-vs-Kronecker guard.** True row-major Kronecker dimension and
   coordinates at (2,4); the 5-dimensional block-concatenation impostor has
   wrong dimension/coordinates.
5. **Ambient/dimension plant.** Full evaluation space of a factor pair has
   dimension exactly $r_Ar_B$ (8 at (2,4)): any 9 distinct cells are
   dependent; a planted 9-cell set has rank exactly 8 in the 8-dimensional
   ambient, and a 9-cell all-distinct sub-$k$ set at (2,4) (where $9>k=6$)
   confirms rank cannot exceed the ambient while staying below $k+1$
   exactly as measured.
6. **Single-map / wrong-seed REJECT.** On the ACCEPT circuit: corrupting one
   coordinate of the seed triple makes the reconstruction fail on the
   remaining points (no nondegenerate map fits), demonstrating the map is
   determined by the data, not by the instrument.
7. **(2,2) corner.** On $\{1,2,3,4\}^2$ full matchings over GF(7), GF(11),
   GF(13): dependent 4-sets $\iff$ $\det[1,x,y,xy]=0$ $\iff$ cross-ratio
   equality $\iff$ reconstruction succeeds — the $k=4$ corner of Theorem C,
   consistent with gate H-ALLDISTINCT-CROSSRATIO.

## 8. Stage plan for `run_grs_converse.py`

- **A — library + toolkit identities (O5).** Row-major Kronecker columns;
  dual ranks (batch RREF + incremental basis) on every computation;
  nullspace; canonical PGL enumeration; cr utilities; obligation (b)
  identity of gate 15 re-verified.
- **F — control battery (§7) before any census.**
- **B — PGL layer:** $|\mathrm{PGL}(2,p)|=p(p^2-1)$ = 2184 / 1320 / 336;
  per-map domains $D_M$; eq-b sums.
- **C — certified anchors (§5)** in primes order 13, 11, 7; eq-a, eq-b,
  sub-$k$ emptiness, (O1)–(O4) on every dependent set.
- **D — new coverage (§6)** with the same pipeline; budget-checked per
  family.
- **E — (O6) separator verifications** across registered configurations and
  sampled sets.
- **G — wrap up:** budgets, results file, assertion count.

Per-stage checkpoints to `checkpoint.json`; final results to
`controls_results.json`.

## 9. Arithmetic, runtime discipline, budgets, defects

Exact GF(p) integer arithmetic only: Python ints mod $p$, inverses by
`pow(x,-1,p)`; no floats, no NumPy/sympy, no external algebra;
`sys.dont_write_bytecode=True` before any import; launched with
`PYTHONDONTWRITEBYTECODE=1` and all five thread caps $=1$; `nice -n 10`
asserted; `RLIMIT_CPU` soft $=$ hard $=5400$ set before the first long stage.
**Hard cap 5400 CPU seconds / 6000 wall seconds**, checked by every
exhaustive family via the shared `check_budget`; no adaptive extension. The
launcher polls artifact files rather than blocking indefinitely and verifies
the process is not OS-napped (R state, ~0 CPU growth) before relaunching.

Any defect is disclosed in `defect_log.json`, broken output preserved under a
distinct name, and every affected family rerun from its beginning (gate-15
precedent: five disclosed aborts, all remediated pre-freeze). No assertion
weakened, no count silently repinned, prereg not amended. If the prereg must
be amended after init, the provenance record carries both hashes.

## 10. Promotion rule and verdict

**FROZEN-CERTIFIED** iff: Theorem E and Theorem C are registered as proved
with their in-run obligations (O1)–(O6) verified as specified (every
enumerated dependent set passing (O1)–(O4); separator checks (O6) passing);
all six certified anchors of §5 reproduced exactly with eq-a support-set
equality and eq-b counts; the §6 new-coverage sweeps are internally
consistent (direct $=$ PGL-family, all dependent sets Möbius); and the full
control battery of §7 passes with every planted REJECT firing. Otherwise
**FROZEN-INCONCLUSIVE**, naming the first failed obligation or control while
certifying inside the record every part that did pass.

Both statements are proved analytically in §§2–3, so no partial-outcome
branch is anticipated; if an in-run failure exposes a proof gap, the gap is
named exactly (step (i)–(iv) or the failing obligation), nothing incomplete
is promoted, and the affected statement is recorded as open with that exact
failing step.
