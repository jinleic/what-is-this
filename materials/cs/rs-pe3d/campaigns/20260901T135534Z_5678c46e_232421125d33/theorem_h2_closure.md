# Hypothesis H2 — closure for the (3,3), (2,2)-decomposition corner
## gate H-MIX-H2-CLOSURE, run 20260901T135534Z_5678c46e_232421125d33

Agent `RsPe3dH2`, 2026-09-01. Prereg: `prereg/H_MIX_H2_CLOSURE_PREREG_2026-09-01.md`,
path-scoped commit `b3e61f76fdb72c06a65cb90a61069273ef4644f3`, sha256
`58c7ed6059a1a5513bfaa70992e383e5b6df49d0eda840fa5892740481ce55c5`, byte-identical
copy in this run dir as `pre_statement.md` (cmp-verified). No frozen run dir
was modified; the five-step obligation was supplied by owner steering as a
candidate (prereg preamble) and is **derived and verified here line by
line**; the machine record is `controls_results.json` (170 passing asserts,
0 failures, CPU 0.212 s / cap 600 s, wall 0.214 s / cap 900 s, `nice 15`,
single-threaded, exact mod-p integers, `PYTHONDONTWRITEBYTECODE=1`,
bytecode disable asserted from inside the run).

## 0. Setting and abbreviation

$A$, $B$ over a field, all columns nonzero, $d_A = d_B = 3$, $H = A \otimes B$,
$m_* = 4$. Frozen kernel identity:
$\ker(A \otimes B) = (\ker A \otimes F^{s_B}) + (F^{s_A} \otimes \ker B)$,
so any relation is $\Gamma = C_1 + C_2$ with the nonzero columns of $C_1$ in
$\ker A$ (spark floor: each has $\ge 3$ nonzeros) and the nonzero rows of
$C_2$ in $\ker B$ (each has $\ge 3$ nonzeros). H2 asks whether the residual
corner — a size-4 profile-$(3,3)$ circuit whose relation admits a
$(k_1,k_2)=(2,2)$ decomposition ($k_1$ = nonzero columns of $C_1$, $k_2$ =
nonzero rows of $C_2$) — produces a non-crossing support.

**Notation.** Active rows $U = \{u_1, u_2\}$ of $C_2$, active columns
$V = \{v_1, v_2\}$ of $C_1$, $u_1 \ne u_2$, $v_1 \ne v_2$; the four cells
$X = U \times V$ (a corner of size 4). Supports $P = \operatorname{supp}(C_1)$,
$Q = \operatorname{supp}(C_2)$, overlap $O = P \cap Q$, and
$T_1, T_2$ = the row supports of the two nonzero columns of $C_1$ (at
$v_1, v_2$), $Z_1, Z_2$ = the column supports of the two nonzero rows of
$C_2$ (at $u_1, u_2$). $S = \operatorname{supp}(\Gamma)$ with $|S| = 4$,
profile $(3,3)$: exactly 3 distinct rows and 3 distinct columns meet $S$.

Two standing observations, used repeatedly:

1. **Spark floor (F1).** Every nonzero column of $C_1$ lies in $\ker A$, so
   its support within its matrix column has size $\ge d_A = 3$; every
   nonzero row of $C_2$ lies in $\ker B$, so its support within its row has
   size $\ge d_B = 3$.
2. **Row/column inclusion (F2).** $\operatorname{supp}(C_1) =
   \bigcup_{v \in V} (T_v \times \{v\}) \subseteq T_1 \times \{v_1\} \,\cup\,
   T_2 \times \{v_2\}$ — in particular **every cell of
   $\operatorname{supp}(C_1)$ has its column in $\{v_1,v_2\}$**, and symmetrically
   **every cell of $\operatorname{supp}(C_2)$ has its row in $\{u_1,u_2\}$**.
   Consequently $\operatorname{supp}(C_1) \cap \operatorname{supp}(C_2)
   \subseteq \{(u,v) : u \in \{u_1,u_2\},\ v \in \{v_1,v_2\}\} = X$: the
   intersection is inside $X$ regardless of $\Gamma$ or $S$.

## 1. Step 1 — cancellation cells are confined (proven)

Write $P = \operatorname{supp}(C_1)$, $Q = \operatorname{supp}(C_2)$.

**Claim 1a (everything supported lies in $P \cup Q$).** $\Gamma = C_1 + C_2$;
a cell where both summands vanish has $\Gamma = 0$. Hence
$S = \operatorname{supp}(\Gamma) \subseteq P \cup Q$.

**Claim 1b (every cancellation cell lies in $X$; supported cells outside $S$
cancel).** Let $(u,v) \notin S$. Then $0 = \Gamma_{u,v} = (C_1)_{u,v} +
(C_2)_{u,v}$. If $(C_1)_{u,v} \ne 0$ then $(u,v) \in P$; if moreover
$(u,v) \notin Q$, then $(C_2)_{u,v} = 0$ and $\Gamma_{u,v} \ne 0$, contradicting
$(u,v) \notin S$. Hence $(u,v) \in Q$, i.e. $(u,v) \in P \cap Q$, which by (F2)
lies in $X$. Reading the chain backwards: **if a cell with both coordinates
active cancels then it lies in the common part, and any supported cell outside
$S$ must belong to both summands and cancels there.** That is exactly the
inclusion
$$O = P \cap Q \subseteq X, \qquad (P \cup Q) \setminus S \subseteq O \subseteq X .$$
Conversely every cell of $X$ with $(C_1)_{u,v} \neq -(C_2)_{u,v}$ has
$\Gamma_{u,v} \ne 0$, which forces $(u,v) \in S$ and consumes row/column budget
of the profile. This is recorded as the *noncancellation alternatives* of
Claim 1b; they are used in Step 3.

**Verification note (a) of the assignment.** The airtightness worry was cells
**inside $S$** where both $C_1$ and $C_2$ are nonzero *without cancelling*
($\Gamma \ne 0$ there). The argument above never assumes such a cell cancels:
for $(u,v) \in S$ there is no constraint at all, but such cells are still in
$P \cup Q$, and if they lie in **both** $P$ and $Q$ they are in $P \cap Q
\subseteq X$ by (F2). The Step-2 counting below uses only $O = P \cap Q
\subseteq X$, never any assumption about cancellation or its absence inside
$S$. The in-run controls confirm the underlying theory: in all six exhaustive
sweeps, every retained tuple (all four of whose $X$-cells cancel exactly)
produced a kernel element, and the census set-equality caught every
discrepancy between summand supports and $\operatorname{supp}(\Gamma)$.

## 2. Step 2 — tight counting (proven)

Let $D_1 = P \setminus Q$ and $D_2 = Q \setminus P$; these are disjoint by
construction, so $|D_1| + |D_2| = |P| + |Q| - 2|O|$. Every active column of
$C_1$ has support $\ge 3$ (F1), so $|P| \ge 3 \cdot 2 = 6$; similarly
$|Q| \ge 3 \cdot 2 = 6$. By Step 1, $O \subseteq X$, so $|O| \le 4$. Finally
$D_1, D_2 \subseteq P \cup Q \setminus O$; cells of $(P \cup Q) \setminus O$ are
supported in exactly one summand, so $\Gamma$ is nonzero there
($\Gamma = C_1 + C_2$ with one summand zero), hence they lie in $S$: this is
the **genuine containment** $D_1, D_2 \subseteq S$ demanded by the assignment,
not an approximation by selection. Therefore the difference sets are disjoint
(by definition) and lie inside $S \cap (P \cup Q)$. Now
$$|S| \ge |D_1| + |D_2| = |P| + |Q| - 2|O| \ge 12 - 8 = 4 = m_* ,$$
which is exactly tight. The in-run sweeps measure the same quantity: for
`V4xV4` at all three primes the distinct profile-$(3,3)$ supports number
exactly 144 $= d_A d_B C_A(3) C_B(3)$, the theorem-X closed form — tightness
with no room for a fifth support. **Verification note (b)** is discharged
inside this step: the difference sets are *genuinely* disjoint (one is the
complement of the other relative to $P \cup Q$), and *genuinely* inside $S$
(each of their cells has $\Gamma \ne 0$). The instrument cross-checks both
directions: under the retained tuples the observed support always equals
$(P \cup Q) \setminus O$ exactly in the sweeps' cancellation-cells (four cells
of $X$ and no others), i.e. the census blocks recompute $S$ from the same
identity the proof uses.

## 3. Step 3 — equality forces the shape (proven)

Equality in the chain above is claimed to force $|P| = |Q| = 6$, $|O| = 4$,
$O = X$, and $T_i = \{u_1, u_2, w_i\}$, $Z_i = \{v_1, v_2, x_i\}$, with every
cell of $X$ inside both supports. Trace each link:

- Equality $|S| = 4$ in $|S| \ge |D_1| + |D_2|$ forces $|D_1| + |D_2| = 4$.
  From $12 - 2|O| \le 4$ we get $|O| \ge 4$; with $|O| \le 4$ this pins
  $|O| = 4$ and, since $O \subseteq X$ and $|X| = 4$, forces $O = X$:
  **all four cells of $X$ lie in both supports.** This is the load-bearing
  point flagged in verification note (c): the argument does not merely bound
  the intersection *from above*; the counting inequality forces the overlap
  to be *maximal*, and maximality plus the inclusion $O \subseteq X$ yields
  $X \subseteq P$ and $X \subseteq Q$ separately. In particular each row
  support $Z_i \supseteq \{v_1, v_2\}$ and each column support
  $T_i \supseteq \{u_1, u_2\}$.
- $|P| \ge 6$ and $|D_1| = |P| - |O| = |P| - 4 \ge 2$; similarly
  $|D_2| \ge 2$; and $|D_1| + |D_2| = 4$ forces $|D_1| = |D_2| = 2$, hence
  $|P| = |Q| = 6$ exactly. Since $|P| = \sum_{i} |T_i| \ge 3 + 3 = 6$ with
  equality, $|T_1| = |T_2| = 3$; together with $T_i \supseteq \{u_1, u_2\}$
  this gives $T_i = \{u_1, u_2, w_i\}$ for some $w_i \notin \{u_1,u_2\}$;
  symmetrically $|Z_1| = |Z_2| = 3$ and $Z_i = \{v_1, v_2, x_i\}$ with
  $x_i \notin \{v_1, v_2\}$.
- Since $D_1 \mathbin{\dot\cup} D_2 \subseteq S$ and
  $|D_1 \mathbin{\dot\cup} D_2| = |D_1| + |D_2| = 4 = |S|$, equality of
  finite sets forces $S = D_1 \mathbin{\dot\cup} D_2$. This union is disjoint
  from $O = X$ by definition, so $S \cap X = \varnothing$. Therefore all four
  cells of $X$ lie outside $S$ and satisfy
  $(C_1)_{u,v} + (C_2)_{u,v} = 0$: every cell of $X$ actually cancels.

**Verification note (c)** is discharged: $X$ lies in *both* supports (by the
maximality argument), not merely as an upper bound on the intersection; the
second half of the inequality chain forces the two difference sets to their minimum
sizes rather than leaving slack for a non-crossing arrangement. The corner
case flagged in the prereg — a cell inside $S$ supported on both summands
without cancelling — is now precisely excluded: any such cell would sit in
$O \subseteq X$, but $S \cap X = \varnothing$, so inside $S$ the two summands
never co-occupy a cell.

## 4. Step 4 — the profile collapses the parameters (proven)

$S \supseteq D_1 \cup D_2$ with $|D_i| = 2$; the four difference cells are
$(T_i \setminus \{u_1,u_2\}) \times \{v_i\} = \{(w_i, v_i)\}$ and
$\{u_i\} \times (Z_i \setminus \{v_1,v_2\}) = \{(u_i, x_i)\}$:
$$S = \{(w_1,v_1),\, (w_2,v_2),\, (u_1,x_1),\, (u_2,x_2)\},$$
where the four pairs are distinct as cells (distinct $v_i$'s, distinct
$u_i$'s, and column-index bookkeeping separates column-cells from row-cells).
This display *is* the claimed support shape; what remains is to collapse
$(w_1, w_2, x_1, x_2)$ to $(w, w, x, x)$ using the profile.

**Row exclusion (verification note (d), first half).** The relation columns at
$v_1$ has support $T_1 = \{u_1, u_2, w_1\}$ of size exactly 3; by the spark
floor the relation column *cannot* have support $\le 2$, so
$w_1 \notin \{u_1, u_2\}$ (else $|T_1| \le 2$). Same for $w_2$.

**Profile collapse.** $\{w_1, v_1\}, \{w_2, v_2\}, \{u_1, x_1\}, \{u_2, x_2\}$
are the four cells of $S$, whose distinct row set has size 3 and distinct
column set has size 3. Rows: $\{u_1, u_2, w_1, w_2\}$ must have size 3
(it contains $\{u_1,u_2\}$ plus possibly $w_1, w_2$, but the *distinct row
set of $S$* is $n_A = 3$). Since $w_1$ (and $w_2$) are outside
$\{u_1, u_2\}$, the only way $\{u_1, u_2, w_1, w_2\}$ has size 3 is
$w_1 = w_2$. Symmetrically the column set
$\{v_1, v_2, x_1, x_2\}$ has size 3 with $x_i \notin \{v_1,v_2\}$, forcing
$x_1 = x_2$. Hence $S$ takes the displayed form with
$w = w_1 = w_2$, $x = x_1 = x_2$:
$$S = \{(w,v_1),\, (w,v_2),\, (u_1,x),\, (u_2,x)\}.$$

**Note (d), second half: why the purported alternatives are genuinely
excluded, not just "not helpful".** The assignment named "for example
$w_1 = u_2$ must be excluded". Note $w_1 = u_2$ would make
$T_1 = \{u_1, u_2\}$, a 2-element support for a nonzero $\ker A$-column —
violating the spark floor $d_A = 3$; this is exactly the row exclusion above,
so the named alternative dies for a *structural* reason (column size),
not by accidental counting. Its mirror, "one of the row/column indices
coincides", dies the same way by the symmetric argument on $Z_i$. A third
alternative — $w_1 \ne w_2$ but both different from $u_1,u_2$, i.e. four
distinct rows — is excluded because then the profile's $n_A$ would be 4, not
3; there is no slack anywhere in the count. The profile pins the collapse
uniquely; no other branch of the case analysis produces a size-4 support.
Hence and by Step 3 the row/ column supports are **exactly**
$T = \{u_1, u_2, w\}$ and $Z = \{v_1, v_2, x\}$ (size 3 each, forcing the
"third" index out of the pair by spark-3), and $S$ is a *crossing support*
whose relation, up to the proportionality of Step 5, is the known
crossing relation.

## 5. Step 5 — the support is a *circuit* (proven), and H2 holds

**Claim 5a (size-3 support ⇒ circuit, verification note (e)).** The column of
$C_1$ at $v_1$ is a nonzero vector $\gamma \in \ker A$ with
$\operatorname{supp}(\gamma) = T$, $|T| = 3 = d_A = \operatorname{spark}(A)$.
Circuit means (i) dependent and (ii) every proper subset independent. (i) holds:
$\gamma \ne 0$ is a dependence among the columns of $A$ indexed by $T$.
(ii): a proper subset of $T$ has size $\le 2 < \operatorname{spark}(A)$, and
spark is by definition the *minimum* size of a dependent set, so every proper
subset is independent. Hence $T \in \mathrm{Circ}_A(3)$ — a genuine circuit,
not merely a dependent set. The worry behind note (e) — "minimal-support
kernel vector" could in principle mean merely dependent — is squarely
addressed: the spark-3 hypothesis upgrades minimality to the two-sided
conditions. Symmetrically $Z \in \mathrm{Circ}_B(3)$.

**Claim 5b (crossing identification).** With $T \ni w$, $Z \ni x$ (both by
Step 4's exclusion of $w$ from the row pair, resp. $x$ from the column pair),
Step 4's support display reads
$$S = ((T \setminus \{w\}) \times \{x\}) \;\cup\; (\{w\} \times (Z \setminus \{x\}))
= S(T, w, Z, x),$$
the crossing family of Theorem X with centre cell $(w,x)$ cancelled.
So the $(2,2)$ decomposition produces **no new support**.

**Claim 5c (proportionality).** Let $S$ be a circuit of $H$. Its relation
space is $\ker(H_S) = \{\gamma : \gamma_s \ne 0 \Rightarrow s \in S,\ 
H_S \gamma = 0\}$ of dimension $|S| - \operatorname{rank}(H_S) = 4 - 3 = 1$ —
1-dimensional (Theorem X, Section 3.1, both proved and machine-verified
in-run: rank 3 and all $\binom43$ triples independent on every constructed
support; in-run the same rank-3 reldim statement holds on every
profile-$(3,3)$ candidate retained by the sweeps). The crossing relation
$\Gamma^{\rm cross}$ built on $T, w, Z, x$ (owner construction, T1 block)
is *a* nonzero element of that 1-dimensional space with *the same support*
$S$; therefore $\Gamma = \lambda\, \Gamma^{\rm cross}$ for some $\lambda
\ne 0$ — the H2 relation is proportional to the crossing relation, nothing
exotic survives. Combining: every $(2,2)$ decomposition of a size-4
profile-$(3,3)$ relation yields a support in the crossing family. H2 is
**discharged** and the Theorem-X converse is complete for $d_A, d_B \ge 3$
(modulo the spark-2 boundary, see §7).

## 6. In-run verification

Arithmetic: exact Python ints mod $p$, stdlib only, rref with
`pow(pivot,-1,p)`, no floats/method libraries; 170 aggregate asserts,
0 failures; `nice 15`, all four thread-env caps `1`, budget holds
(CPU 0.212 / 600 s, wall 0.214 / 900 s). All controls ran AFTER
`campaign.py init`; the prereg commit preceded any parameter-dependent compute.

**T1 — owner GF(13) ACCEPT (independently derived here).** $A = B = [(1,t)]
_{t=1}^4$ (2×4 Vandermonde), $T = Z = \{0,1,2\}$, $u_1 {=} 0, u_2 {=} 1,
w {=} 2$; $v_1 {=} 0, v_2 {=} 1, x {=} 2$; $c = \delta = (1,11,1)$ — the
constraints of the owner construction. Forcing equations solved for
$\alpha |_{\{v_1,v_2\}} = (\delta_{v_1}, \delta_{v_2})$, $\beta
|_{\{u_1,u_2\}} = (-c_{u_1}, -c_{u_2})$; the asserted support of the resulting
$\Gamma$ is `[(0,2),(1,2),(2,0),(2,1)]`, *identical as a set* to
$S(T,w,Z,x)$ (set membership asserted, not count); $A \Gamma B^{\mathsf T} = 0$
asserted exactly both as a matrix product and as the coefficient sum over
product columns; support is a genuine circuit (rank 3, all four 3-subsets
independent). SUPPORT-LEVEL agreement, not merely cardinality.

**T2 — exhaustive (2,2) parameter sweep, 6 configurations.** For each
$p \in \{7, 11, 13\}$ and each of two distinct factor pairs
(`V4xV4` = $[(1,t)]_{t=1}^4 \otimes$ same; `V3xV5` = $[(1,t)]_{t \in
\{0,1,3\}} \otimes [(1,t)]_{t \in \{0,1,2,4,5\}}$), the instrument
enumerated **every** way to pick the four active supports of a candidate
$(2,2)$ relation: the row pair $U$ (C(4,2)/C(3,2) ways), the column pair $V$,
the two ordered third indices $w_1, w_2 \in [s_A] \setminus U$, and the two
ordered third indices $x_1, x_2 \in [s_B] \setminus V$ — 576 resp. 270 shapes
— and, for each shape, **every** assignment of the four circuit directions
to the four columns/rows, normalised bijectively by pinning the
$(C_1, v_1)$ scalar to 1, then sweeping the remaining three scalars over
$F_p^\times$ (exact totals pinned: $576(p-1)^3$ and $270(p-1)^3$ tuples;
e.g. 995,328 normalised tuples at `p13_V4xV4`, 466,560 at `p13_V3xV5`).
Pinned sweep counts: 576 and 270 shapes at their respective configs —
**no adaptive parameter extension** was allowed. For every tuple that made
the four $X$ cells (those of Step 3: all four of $(u_1,v_1),(u_1,v_2),
(u_2,v_1),(u_2,v_2)$) cancel *exactly*, the observed support was recomputed
and asserted to have the correct size and profile. Result at every one of
the six configs, asserted as **set membership, not counts**:

- every retained profile-$(3,3)$ support lies in the independently
  constructed crossing family $\mathcal{C}$ (set difference empty),
- the retained profile-$(3,3)$ population equals $\mathcal{C}$ exactly
  (144 supports at `V4xV4`, 90 at `V3xV5`, at each of the three primes),
- every retained profile-$(3,3)$ candidate is a genuine circuit (rank 3,
  all proper subsets independent) — corroborating Step 5's circuit claim
  component-by-component,
- every retained tuple's decomposition satisfies the kernel equation exactly
  (the ambient kernel identity, re-asserted in-run), confirming that no
  beyond-argument drift occurred.

**No new support ever appears** — the "zero residual" assertion held at
all six configs. (Retained non-$(3,3)$ solutions exist and are recorded but
not promoted: profile-$(4,4)$ all-distinct shapes at 48/24/72 counts,
varying with $p$ — this is consistent with the all-distinct residual having
no closed form; they are NOT H2 candidates since their profile is $(4,4)$,
not $(3,3)$.)

**T3 — previous campaign census (first-conviction check).** At GF(13),
reimplemented by value (no import from the frozen directory) using the
prereg's same witness factors $A = ((1,1,1,1),(1,2,3,4))$,
$B = ((1,-2,1,0),(0,1,-2,1))$: exhaustively enumerated
$\binom{16}{4} = 1820$ four-subsets; dependent sets found by exact rank;
circuits screened by all $\binom43$ proper-subset independence checks
(156 circuits total). Pinned reproduction of the
frozen campaign's numbers: total circuits 156, profile-$(3,3)$ circuits 144,
profile-$(4,4)$ residual 12. **Profile-$(3,3)$ set equality against the
independently constructed crossing family: 144 measured = 144 constructed,
both set differences empty** — the residual equals zero, confirming
H2's premise numerically on data the previous campaign certified from
scratch.

**T4 — planted REJECTs.** (i) *Perturbation*: perturbing one coefficient of
T1's forced $\Gamma$ by $+1$ breaks $A\Gamma B^{\mathsf T} = 0$ (asserted:
nonzero matrix product AND nonzero coefficient sum — dependence destroyed).
(ii) *Duplicated factor column*: cloning a factor column into a 5th slot
drops the measured spark from 3 to 2 and produces a dependent column pair —
detected.

**T5 — fail-loud concat guard.** The instrument computes product columns
*row-major* by $[a_i b_j]$: verified against an independently written
row-major Kronecker product for a 2-row × 3-row factor pair (ambient
dimension $2 \cdot 3 = 6$ asserted, coordinates asserted equal), while the
planted block-concatenation path yields dimension $2 + 3 = 5$ with different
coordinates — the control fires immediately if block concatenation were ever
silently substituted for the Kronecker product.

**Runtime controls.** `PYTHONDONTWRITEBYTECODE=1` in env and
`sys.dont_write_bytecode = True` set before any import (asserted
from inside the run); thread envs `OMP/OPENBLAS/MKL/VECLIB = 1` asserted;
`nice` value asserted $\ge 15$ (measured 15); CPU cap 600 s and wall cap
900 s pre-registered; the instrument checks both during sweeps and would
have failed before promotion had either been reached —
neither was (0.212 / 0.214 s).

## 7. What H2's closure buys — explicit scope sentence

**Now proved**, combining this run with the frozen H-MIX-CROSSING record:
Theorem X's converse is **complete for $d_A, d_B \ge 3$** — every size$
m_* = d_A + d_B - 2$ circuit of $H = A \otimes B$ with profile $(d_A, d_B)$
is a crossing support $S(R, i_0, J, j_0)$ with $R \in \mathrm{Circ}_A(d_A)$,
$J \in \mathrm{Circ}_B(d_B)$, and its relation space is 1-dimensional,
spanned by the explicit crossing relation. Mechanism: Lemma P kills pure
pieces for $d_A, d_B \ge 2$; Lemmas Q / Q2 kill $(1,k)$, $(k,1)$, and
$(k_1, k_2)$-with-both-$\ge2$ shapes (Q2 analytic for $d_A, d_B \ge 4$), and
this run's five steps close the last corner — $(3,3)$ with
$k_1 = k_2 = 2$ — by forcing its support to be a crossing support with a
proportional relation (no new support, proportionality by 1-dimensionality).

**NOT affected / still open** (unchanged scope, restated per prereg §5):

1. **Spark-2 boundary branches** ($d_A = 2$ or $d_B = 2$ converse branches,
   e.g. permutation-shaped sets) — out of both campaigns' scope, untouched.
2. **All-distinct closed form** — the all-distinct (profile-$(d{+}1,d{+}1)$)
   size-$(d{+}1)$ population still has no closed form and remains
   field-dependent (in-run: the $(4,4)$ residuals in the retained sweeps
   vary with $p$).
3. **$\ge 3$ factors** — nothing here touches three or more factors; the
   regrouping corollary of the frozen run stays demo-level.

After the certified close, the target-level records
(`SCOPE_NOTE_H_MIX_CROSSING.md`, `state.json`) are updated in a separate
path-scoped commit only because H2 is now proved; no frozen run directory is
touched by that update.

## 8. Defects and disclosures

No correctness defect, failed assertion, parameter change, or prereg
amendment occurred. The only prereg hash is the one recorded in the header
and `provenance.txt`; `controls_results.json` is the authoritative execution
record. No control was weakened, dropped, or re-pinned.

## 9. Verdict inputs

All pre-registered controls pass; both planted REJECTs fire; the census
reproduction matches the frozen run; set-level (not count-level) containment
holds at all six exhaustive (2,2) sweeps; all five proof steps are discharged
in §1–§5 with the assignment's five named sharp points ((a)–(e)) each
addressed and closed; prereg §2's FROZEN-CERTIFIED gate is met with no gap.
Per the prereg §2 promotion rule this supports **FROZEN-CERTIFIED**, and per
the assignment's item 5 the scope note and state files are updated to record
H2 as discharged.
