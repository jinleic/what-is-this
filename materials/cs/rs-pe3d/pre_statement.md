# PRE-STATEMENT — rs-pe3d — ECCC TR26-150 Conjecture 4.2, k=3

Committed: 2026-08-29 (UTC), by agent `RsPe3d`, BEFORE any computation run.
Status of this document: anti-hindsight contract. Thresholds, precision,
certificate formats, and labelling rules below are FIXED prior to observing
any result. If a gate must be re-scoped, the re-scoping is added as a dated
addendum and the original text is left visible. No numeric thresholds in this
document may be edited after runs begin.

Primary source read first-hand (full text): ECCC TR26-150, Mitali Bafna &
Nikhil Vyas, *"Private PCPs from Product Expansion"* (2026-08-19),
https://eccc.weizmann.ac.il/report/2026/150/ . All verbatim quotes below were
taken from the PDF text extraction of the full report (stored re-paging copy
in-session, ~1600 lines); every load-bearing quote was re-checked against the
extraction this session. Provenance: [REPRODUCED] unless noted.

---

## 1. The claim under study, verbatim

### 1.1 Definition 4.1 (product expansion) [REPRODUCED]

> **Definition 4.1 (Product expansion).** A tuple $(C_i)_{i\in I}$, where
> $I \subseteq [k]$, is $\rho$-product-expanding if every
> $M \in \sum_{i\in I} L_i(C_i)$ has a decomposition $M = \sum_{i\in I} M_i$,
> with $M_i \in L_i(C_i)$, such that
> $\mathrm{wt}(M) \geqslant \rho \sum_{i\in I} s_i \ell_i(M_i).$

Here $L_i(C_i) = \mathbb{F}_q^{S_1} \otimes \cdots \otimes C_i \otimes \cdots
\otimes \mathbb{F}_q^{S_k}$ (the tensor code lifted in axis $i$, identity in
all other axes), $s_i = |S_i|$, $\ell_i(M_i)$ = number of nonzero $i$-axis
lines of $M_i$, and $\mathrm{wt}$ is Hamming weight on
$\Omega = S_1 \times \cdots \times S_k$.

### 1.2 Conjecture 4.2 [REPRODUCED — the exact text is load-bearing]

> **Conjecture 4.2.** For every $\eta \in (0,1)$, balance parameter
> $\beta \geqslant 1$, and integer $r \geqslant 2$, there is
> $\rho = \rho(r,\eta,\beta) \in (0,1]$ such that the following holds. Let $q$
> be prime, and let $S_1,\dots,S_r \leqslant \mathbb{F}_q^\times$ be
> multiplicative subgroups of pairwise coprime orders $s_i = |S_i|$
> satisfying $\beta^{-1} \leqslant s_i/s_j \leqslant \beta$ $(i,j \in [r])$.
> For integers $t_i$ with $1 \leqslant t_i \leqslant (1-\eta)s_i$, and
> invertible diagonal matrices $\Lambda_i$ indexed by $S_i$, the tuple
> $(\Lambda_i \cdot \mathrm{RS}(S_i,t_i))_{i=1}^r$ is $\rho$-product-expanding.
> Note that, without loss of generality, we may assume that $\rho(r,\eta,\beta)$
> is nonincreasing in $r$.

(Extraction notes: "⩾" inside the balance hypothesis was rendered
"$\leqslant$" in raw extraction both times as $s_i/s_j$ comparisons of the
form $\beta^{-1}\le s_i/s_j\le\beta$ — the two-sided $\beta^{-1} \cdots \beta$
form is the only sensible reading and is consistent with the README's
paraphrase and the paper's own Theorem 4.11 hypothesis
"$K^{-1} \leqslant n_1/n_2 \leqslant K$". Also the conjecture's WLOG remark
on nonincreasingness in $r$ is included.)

### 1.3 RS code definition (Section 3) [REPRODUCED]

> define $\mathrm{RS}(S,a) = \{(f(\alpha))_{\alpha \in S}: f \in
> \mathbb{F}_q[X], \deg f < a\}$.

So the parameter $t_i$ equals the **dimension**; evaluation map is injective
when $t_i \le s_i$.

### 1.4 Theorem 2.1 (proved, 2-dimensional) [REPRODUCED]

> **Theorem 2.1.** Fix $\varepsilon > 0$. There exists $\rho = \rho(\varepsilon) > 0$
> [such that] Let $S_i \leqslant \mathbb{F}_q^\times$ have coprime orders
> $n_i = |S_i|$ with $|S_1| \approx |S_2|$ and let $C_i = \mathrm{RS}(S_i,d_i)$,
> where $d_i \leqslant (1-\varepsilon)n_i$. Every matrix $M$ that is a sum of
> $C_1$-columns and $C_2$-rows admits a decomposition $M = A + B$ such that
> $\mathrm{wt}(M) \geqslant \rho \cdot n_1 \cdot \#\{\text{nonzero columns of } A\}
> + n_2 \cdot \#\{\text{nonzero rows of } B\}$.

### 1.5 Theorem 4.11 (proved, author-verified) [REPRODUCED]

> **Theorem 4.11.** Fix $\varepsilon \in (0,1)$ and $K \geqslant 1$. There is a
> constant $\rho = \Omega(\varepsilon^6/K^3) \in (0,1]$ [such that] Let $q$ be
> prime, and let $S_1, S_2 \leqslant \mathbb{F}_q^\times$ be multiplicative
> subgroups of coprime orders $n_i = |S_i|$ satisfying
> $K^{-1} \leqslant n_1/n_2 \leqslant K$. For integers $k_i$ satisfying
> $1 \leqslant k_i \leqslant (1-\varepsilon)n_i$, let
> $C_i = \mathrm{RS}(S_i,k_i)$. Then $(C_1,C_2)$ is $\rho$-product-expanding.
> The same conclusion holds after independently rescaling the coordinates of
> $C_1$ and $C_2$ by nonzero field elements.

The last sentence is exactly the $\Lambda_i$-rescaling robustness, for $k=2$.

### 1.6 The AI candidate claim (intro) [REPRODUCED]

> We have a candidate AI-generated proof of this conjecture which gives
> $\rho = \varepsilon^{\exp(k)}$, but it has not been verified by the authors.

Context: Conjecture 1.4 (informal higher-dim version); $\rho$ a function of
$\varepsilon = 1 - \max_i t_i/s_i$-style rate slack and $k$ only.

### 1.7 AI Disclosure [REPRODUCED]

> ChatGPT was used throughout the research and writing process. The proof of
> Theorem 4.11 and the main ideas in Sections 4.3 and 6.2 were generated by
> ChatGPT 5.6 Sol. The authors independently verified all AI-generated
> arguments and take full responsibility for all the results in this paper.

**Consequence (evidence-handling rule):** Theorem 4.11 is author-verified;
Conjecture 4.2's higher-dimensional candidate proof ("$\rho =
\varepsilon^{\exp(k)}$") is NOT verified. Escalation to Main on
refutation-appearance applies to the former (hypotheses met + certificate in
hand); for the latter, disconfirming data is a normal finding (Gate C) but is
still reported to Main promptly. [Restated from the campaign order.]

### 1.8 Related statements to keep in scope [REPRODUCED; one reading DERIVED]

- **Corollary 4.3:** with $\rho=\rho(k,\eta,\beta)$, every
  $M \in K^\perp \setminus C^\perp$ satisfies
  $\mathrm{wt}(M) \geqslant \eta\,\rho^{\,k-1} N$ where $N = \prod s_i$.
  (Exponent "$k-1$" read from Lemma 4.6's proof structure, which is fully
  reproduced in Section 4.2 and yields $\delta\rho^{k-r}\prod_{i\ge t} s_i$ at
  stage $t$; the corollary display itself is partially garbled in extraction —
  reading tagged **DERIVED**, not REPRODUCED.) Canonical paper instantiation:
  $\rho(k) = \rho(k, 1/128, 4)$ (Corollary 4.5 usage).
- **Lemma 4.6** (Relative dual distance): define
  $\delta_i = \min\{\mathrm{wt}(w): w \in W_i \setminus U_i\}/s_i$,
  $\delta = \min_i \delta_i$; if every suffix tuple
  $(W_t, U_{t+1},\dots,U_k)$ is $\rho$-product-expanding then
  $\mathrm{wt}(M) \geqslant \delta\rho^{k-1}N$ for $M \in K^\perp\setminus C^\perp$,
  where $K_t = K_t \otimes C_{t+1} + V_t \otimes K_{t+1}$, $C_t = \bigotimes V_i$.
- **Lemma 4.8:** $P(\rho,r) = \dfrac{\rho}{r(2+1/\rho^{\,r})}$ [full display:
  the second "$r$" exponent placement is garbled in the extraction; reading
  used in proof: $|B| \le \frac{2+1/\rho^r}{\rho} |A|$ giving
  $\sum_i s_i \ell_{i,\cup} \le \frac{2+1/\rho^r}{\rho}\mathrm{wt}_\cup$; tagged
  **DERIVED** reading, load-bearing only for Gate A commentary, not used in any
  certificate.]
- **Lemma 4.9** uses $s_r\lambda_i \leqslant \frac{P}{\rho}(t_i+t_r)$;
  **Corollary 4.4/4.10:** $C(\rho,k) = \exp(O(k^3+k^2\log(1/\rho)))$.
- **Theorem 4.12** (specialized Corvaja–Zannier): hypotheses include
  $N < q$, $x^{n_1}, y^{n_2}$ multiplicatively independent modulo constants;
  conclusion $\Gamma_F \leqslant C N^{1/3} (ab)^{2/3} + ab$.
- **Lemma 4.13**: constant $c_0 > 0$; $0 \leqslant k_i < n_i$; notation as in
  paper; $t \leqslant c_0 h^3 / N$ regime.
- **Facts 3.2/3.3:** dual of RS; Combinatorial Nullstellensatz.

### 1.9 Falsifiability stance [DERIVED — restated as binding]

For every finite instance $(q,S_i,t_i,\Lambda_i)$, the per-instance optimum
$\rho_{\mathrm{inst}}$ (defined below) satisfies $0 < \rho_{\mathrm{inst}} \le 1$:
proof, $\mathrm{wt}(M) \le \sum_i s_i \ell_i(M_i)$ for every decomposition, so
no finite instance can refute the ∃ρ-infinite-family Conjecture 4.2 itself.
Finite campaigns can produce: (a) first exact values (Gate A); (b) evidence
about balance sensitivity / a $\beta^\ast$ collapse (Gate B); (c) audit of the
AI-claimed *form* $\rho=\varepsilon^{\exp(k)}$ (Gate C). A Gate-B "collapse"
is evidence that the conjecture **statement** needs a balance hypothesis
(publishable correction of the statement, per target README), never a
refutation of its truth. The ONLY escalation trigger is an instance whose
certificate contradicts the author-verified Theorem 4.11 or 2.1 **as stated
with hypotheses met** (2-dimensional, balanced, coprime orders, $k_i \le
(1-\varepsilon)n_i$, diagonal rescalings) — that triggers STOP + hub to Main
with certificate before any record.

---

## 2. The exact quantity computed (Gate A/B/C input)

Fix an instance $q$, pairwise-coprime $s_1,s_2,s_3$ (with $s_i \mid q-1$),
$t = (t_1,t_2,t_3)$ with $1 \le t_i \le (1-\eta)s_i$, $\Lambda_i$. Define

$$
\rho_{\mathrm{inst}} \;=\; \min_{M \in V \setminus \{0\}} \frac{\mathrm{wt}(M)}{\delta(M)},
\qquad
\delta(M) = \min\Big\{ \sum_{i=1}^{3} s_i\,\ell_i(M_i) \;:\; M_i \in L_i(\Lambda_i\mathrm{RS}(S_i,t_i)),\; \textstyle\sum_i M_i = M \Big\},
$$

where $V = L_1 + L_2 + L_3 \subseteq \mathbb{F}_q^{N}$, $N = s_1 s_2 s_3$.

**Exactness chain:**

- $(\Lambda_i\mathrm{RS}(S_i,t_i))$ = coordinatewise-rescaled RS; membership in
  $L_i$ is determined per $i$-axis line: an $i$-line with coordinate
  $x_{-i} \in \prod_{j \ne i} S_j$ carries the restriction of some polynomial
  $f \in \mathbb{F}_q[X]$, $\deg f < t_i$, evaluated at $S_i$, rescaled pointwise
  by $\Lambda_i$'s diagonal entries indexed by $S_i$.
- **Rounding/precision:** all computation in exact $\mathbb{F}_q$ arithmetic
  ($q$ prime $\le$ ~700 in-scope; elements as Python ints mod $q$; no
  floating-point anywhere in the verifier; rationals output as exact
  fractions via `fractions.Fraction`/flint `fmpq`, printed `"p/q"`).
- Pass/fail per instance vs a conjectured form value $\rho_{\mathrm{claim}}$:
  PASS iff $\rho_{\mathrm{claim}} \le \rho_{\mathrm{inst}}$ EXACTLY (integer
  comparison $p_1 q_2 \le p_2 q_1$), i.e. the conjectured $\rho$ is compatible
  with the instance; FAIL iff strict opposition; INFEASIBLE-INSTANCE if no
  coprime triple satisfies the constraint.

**Fixed thresholds (pre-registered):**

- $0 < \rho_{\mathrm{inst}} \le 1$: sanity assertion. Any computed
  $\rho_{\mathrm{inst}} > 1$ is a **code bug**: $\mathrm{wt}(M) \le \sum_i
  s_i\ell_i(M_i)$ for every decomposition, so the min ratio cannot exceed 1.
  If the pipeline ever outputs $>1$: STOP, debug, do not record.
- $\rho_{\mathrm{inst}} = 0$: impossible ($\mathrm{wt} \ge 1$ for $M \ne 0$, and
  $\delta$ finite ⇒ ratio $> 0$). If observed: bug, same STOP.
- Refutation thresholds for Gate C: fixed AFTER the mapping
  $\varepsilon \mapsto$ numeric value is decided (Section 5), before any run.

---

## 3. Instance enumeration (Gate A)

- **Primes (only primes — Conjecture 4.2 requires $q$ prime).** The target
  README's illustrative sweep includes $q = 25$, which is **not prime**;
  recorded here as a README discrepancy, not silently followed. Replaced by
  the next applicable primes.
- **Coprimality:** orders $s_i \mid q-1$, pairwise coprime, all $s_i \ge 2$
  (a subgroup of order 1 makes line weights degenerate; also excluded by the
  conjecture's implicit use — recorded as interpretation). Vacuous primes:
  $q=13$ ($12 = 4\cdot3$: no triple with all $s_i\ge2$ pairwise coprime...
  $\{4,3\}$ only 2 factors), $q=37$ ($36=4\cdot9$: same). First active:
  $q=61$ ($60 = 4\cdot3\cdot5$).
- **Sweep (in-scope, subject to compute budget):**
  $q \in \{31, 61, 109, 181, 241, 421, 601\}$ — full coprime-triple census
  per $q$ (divisor triples of $q-1$, pairwise coprime, all $\ge 2$), with
  $q=31$ ($s=(2,3,5)$) as the **pilot**.
- **$\eta \in \{1/128, 1/64, 1/32\}$** intersected with feasibility
  $t_i \le \lfloor(1-\eta)s_i\rfloor$ (for $s_i \le 127$, $\eta=1/128$ is
  equivalent to $\eta=1/64$ constraint-wise at these sizes — dedupe rows,
  recorded).
---

## 9. Re-scoping addendum (2026-08-29, before first computation)

The original Anchor 3 above stated that the $q=31,s=(2,3,5),t=(1,1,1)$
space could be fully enumerated because the displayed **upper bound** on
$\dim V$ is $31$. This was a mistaken use of an upper bound: the actual
dimension is at most $N=30$ and, for $t=(1,1,1)$, the exact tensor-dual
dimension calculation gives
$\dim V = N-(s_1-t_1)(s_2-t_2)(s_3-t_3)=30-1\cdot2\cdot4=22$
(this identity is also checked directly by exact row reduction in `src/`).
Enumerating $q^{22}$ vectors is not a feasible finite campaign. The original
claim is retained above as historical text and is **withdrawn**.

Replacement scope, fixed before any computation:

1. For every reported witness $M$, $\delta(M)$ is either closed by an
   exhaustive exact line-support enumeration below the returned cost, or is
   marked solver-supported COMPUTATIONAL-EVIDENCE. The independent verifier
   always checks all $\mathbb F_q$ equations, line counts, weight, and the
   exact fraction.
2. A global $\rho_{\rm inst}$ is called MACHINE-VERIFIED only if the outer
   minimization is closed by a complete finite certificate (for example an
   exhaustive support/value census, or an independently checked exact
   universal certificate). Otherwise the reported quantity is explicitly
   $\rho_{\rm inst}^{\rm upper}$, the best verified ratio in the stated
   candidate class, and the unsearched complement is OPEN.
3. The first executable campaign uses $q=31,s=(2,3,5)$, identity $\Lambda$,
   all feasible $t$ profiles at $\eta=1/32$, and a pre-registered candidate
   class consisting of exact kernel representatives for all supports of
   weight at most $W_{\rm cap}=3$, plus the line-support MILP witnesses
   listed in each manifest. The weight-$\le3$ census is complete as a
   support census; it is not a claim about heavier $M$.
4. The earlier sentence in Section 4 saying that the full $V$ enumeration
   has $15+10+6=31$ elements is corrected in interpretation: those are
   dimensions of generating lift spaces, not a number of field elements.
   No computation will use that sentence as a completeness argument.

This addendum is a correction of feasibility, not an after-the-fact numeric
threshold change. The numerical threshold and fraction comparison rules
remain unchanged.
- **$(s_1,s_2,s_3)$ profiles:** all coprime triples per $q$; **balanced** and
  the **most unbalanced** valid triples per $q$ are always included in the
  table (highest-value probe per README).
- **$\Lambda$ policy (FIXED):** headline Gate-A values at
  $\Lambda_1=\Lambda_2=\Lambda_3=\mathrm{Id}$; a separate **$\Lambda$-probe**
  run (random invertible diagonals, fixed seed documented in campaign) on a
  representative subset recorded as COMPUTATIONAL-EVIDENCE for the
  $\forall \Lambda_i$ quantifier. $\Lambda$ is NOT trivially removable by
  support invariance: the equation $\sum_i M_i = M$ mixes per-axis rescalings,
  so line supports change.
- **$t$-sweep policy (FIXED):** pilot $q=31$: full $t$-grid
  $t_i \in \{1,\dots,\lfloor(1-\eta)s_i\rfloor\}$ for each $\eta$. For larger
  $q$: pre-registered corner stratification $t_i \in \{1,\ \lfloor s_i/2
  \rfloor,\ \lfloor(1-\eta)s_i\rfloor\}$ per axis (27-profiles/triple/η),
  full grid only where compute allows. Conjecture 4.2 quantifies
  $\forall t_i$; any sampled row is labeled evidence-only (see §6 labels).

---

## 4. Solver formulation (FIXED before first run)

**Primary formulation: exact-$\delta$ enumeration + min-is-attained witness.**

For a fixed instance, the exact algorithm is:

1. Enumerate line-support candidate sets $B_i \subseteq \{x_{-i} :
   x_{-i} \in \prod_{j\ne i}S_j\}$ (there are $N/s_i$ possible $i$-lines).
   For small instances ($N/s_1 \le ~20$ per axis group, pilot scale),
   enumerate all triples $(B_1,B_2,B_3)$ with $\sum_i s_i|B_i| \le D_{\max}$
   for a fixed outer-search budget $D_{\max}$
   (pilot: $D_{\max} = \max(s_i)\cdot\sum_i |B_i|$-bounded, exact).
2. For each $(B_1,B_2,B_3)$: the space
   $W(B) = \{A_1+A_2+A_3: A_i \in L_i(C_i),\ A_i$ zero off $B_i\}$ has dim
   $\le \sum_i t_i |B_i|$; test $M \in W(B)$ by exact Gaussian elimination
   in $\mathbb{F}_q$; $\delta(M) = \min\{\sum_i s_i|B_i| : M \in W(B)\}$
   (minimize over all $B$ with $M \in W(B)$; ties broken by enumerating all
   $B$ at the minimum size — exhaustive, hence EXACT per $M$).
3. $\rho_{\mathrm{inst}}$: minimize $\mathrm{wt}(M)/\delta(M)$ over
   $M \in V/\{0\}$ reached by the campaign's search (enumeration of
   line-support patterns of $M$ + lowest-weight console via the
   $W(B)$-coset structure). Full exhaustive minimization over all of
   $V\setminus\{0\}$ is only claimed when the search provably dominated
   every $(B_1,B_2,B_3)$ pattern; otherwise the instance row states exactly
   what was searched (see evidence labels).
4. **Rational LP/ILP cross-check (HiGHS via `highspy`, version recorded**
   from `highspy.__version__` / HiGHS banner): rewrite
   $\rho_{\mathrm{inst}}$ as min-ratio = LP for a FIXED witness $M$:
   $\min \sum_i s_i \ell_i$ s.t. $\sum_i M_i = M$,
   $\ell_i$ line-indicator sums, nonzero-integer field coefficients
   mod-$q$-generalized (ILP over $\mathbb{Z}_{\ge 0}$ with modular equality
   rows). Dual returned by HiGHS is **re-solved as an exact-rational LP**
   with `fmpq`-backed simplex (flint rational arithmetic) to machine-verify
   the per-instance optimum WHEN the LP relaxation is exact; if the returned
   solution is integral VERIFIED by exact re-decoration and dual closure,
   row label = **MACHINE-VERIFIED**; if solver-only without exact dual
   closure, row label = **COMPUTATIONAL-EVIDENCE** (per repo evidence rules;
   README explicitly requires exact-rational dual verification for
   "optimum" claims).

**Anti-overclaim rule (FIXED):** if for a row the global min over
$M \in V\setminus\{0\}$ is only reached by heuristic search (not exhaustive
proof of optimality), the row reports $\rho_{\mathrm{inst}}^{\mathrm{upper}}$
= the best ratio found (an upper bound on the true min), the label
COMPUTATIONAL-EVIDENCE, and the search protocol. A claimed exact
$\rho_{\mathrm{inst}}$ listing REQUIRES the exhaustive pattern argument
(Section 4 step 3) plus the LP dual closure at the found witness. No
middle-ground claim.

**Certificate format (FIXED, per instance):**
JSON in `campaigns/<UTC>_<uuid>/` with fields: $q$; $(s_1,s_2,s_3)$;
$(t_1,t_2,t_3)$; $\eta$; $\beta$ as used; $\Lambda$ id/seed;
$\rho_{\mathrm{inst}}$ exact string `"p/q"`; witness $M \in \mathbb{F}_q^N$
as integer list; decomposition $(M_1,M_2,M_3)$ certificates (integer lists
$\bmod q$); line counts $\ell_i$; $\delta$ value; solver name+version;
exact-dual-closure yes/no; code hash12; seed; UTC timestamp; evidence label.

**Independent verifier (FIXED):** a separate script re-reads the
certificate and re-checks, in exact $\mathbb{F}_q$: (a) each
$M_i \in L_i(\Lambda_i\mathrm{RS}(S_i,t_i))$ per $i$-line existence of a
polynomial of degree $< t_i$ matching (Gaussian elimination on the
Vandermonde subsystem, scaled by $\Lambda_i$ at $S_i$); (b) $\sum M_i = M$;
(c) recomputed $\ell_i$, recomputed wt; (d) the inequality
$\mathrm{wt}(M) \ge \rho_{\mathrm{inst}} \sum_i s_i\ell_i(M_i)$ with exact
fraction arithmetic. Verifier runs are part of the campaign; campaign
directory is frozen (append-only history via new files) after the verifier
passes.

**Unit anchors (must pass before any real instance):**

- Anchor 1 (full-space degenerate): $C_i = \mathbb{F}_q^{S_i}$ (i.e.
  $t_i = s_i$, outside conjecture range but a valid sub-routine test):
  then $L_i = \mathbb{F}_q^N$, $V = \mathbb{F}_q^N$, and for any $M$, the
  decomposition $M_i = M$ (projection onto one line each) must give
  $\delta(M) = \min(s_1 \ell_1(M), s_2\ell_2(M), s_3\ell_3(M))$... HAND-CHECK
  on a toy $q=7$, $s=(2,3,4)$ toy grid (coprimality not needed for the
  units).
- Anchor 2: $t=(1,1,1)$; each $L_i$ = line-spanned grid; cross-check exact
  $\delta$ and $\rho$ against brute-force enumeration on a tiny grid
  ($q=5$, $s=(2,2,2)$ forced non-coprime for the unit) — cross-tool
  agreement only, labeled unit-test, never a campaign row.
- Anchor 3: on the pilot $q=31$, $s=(2,3,5)$: for $t=(1,1,1)$, verify by
  direct enumeration over the FULL space $V$ (dim $\le 1\cdot(30/2) +
  1\cdot(30/3) + 1\cdot(30/5) = 15+10+6 = 31$ elements → enumerate all
  nonzero $M \in V$ exactly) that the reported $\rho_{\mathrm{inst}}$ is
  achieved and minimal. This is a complete enumeration ⇒ a MACHINE-VERIFIED
  exact value for that row.

---

## 5. Gate C: the AI-claim mapping (FIXED before any run)

The AI candidate proof's claim, verbatim: "$\rho = \varepsilon^{\exp(k)}$".
To instantiate at measured instances, fix the mapping:

- $\varepsilon := \eta$ (the Conjecture 4.2 rate-slack parameter:
  $t_i \le (1-\eta)s_i$ = $d_i \le (1-\varepsilon)n_i$ mapping; this is an
  **[INFERENCE]** — the AI claim does not spell out its $\varepsilon$; the
  paper's Theorems 2.1/4.11 use exactly the $(1-\varepsilon)n$ form, so
  $\varepsilon = \eta$ is the canonical reading; alternatives explicitly
  listed here for sensitivity: $\varepsilon := 1 - t_i/s_i$ per axis with
  the min over $i$, which equals $\eta$ at the boundary $t_i =
  (1-\eta)s_i$ and is larger otherwise — recorded as a secondary mapping
  $\varepsilon' = \min_i(1-t_i/s_i) \ge \eta$).
- $k := r = 3$; $\exp(k)$ read as $e^k$ (natural exponential, standard
  reading; recorded).
- $\rho_{\mathrm{claim}}(\eta) = \eta^{e^3} = \eta^{20.0855\ldots}$.
  For $\eta \in \{1/32,1/64,1/128\}$:
  $\rho_{\mathrm{claim}}(1/32) = 32^{-20.0855\ldots} \approx 10^{-30.2}$;
  $\rho_{\mathrm{claim}}(1/64) = 64^{-20.0855\ldots} \approx 10^{-36.2}$;
  $\rho_{\mathrm{claim}}(1/128) = 128^{-20.0855\ldots} \approx 10^{-42.3}$.
  (All < any plausible $\rho_{\mathrm{inst}}$, so the AI form is trivially
  satisfiable at these $\eta$ — recorded as the expected outcome; the
  interesting numeric is the margin
  $\rho_{\mathrm{inst}}/\rho_{\mathrm{claim}}$ per instance.)
- **Gate C verdict rules (FIXED):** consistent-with-form iff
  $\rho_{\mathrm{claim}}(\eta) \le \rho_{\mathrm{inst}}$ (exact rational
  comparison to the rational rounding of $\eta^{e^3}$ computed at 50-digit
  precision, then re-exactified as the fraction closest strictly below —
  documented per row); "inconsistent" iff
  $\rho_{\mathrm{claim}}(\eta) > \rho_{\mathrm{inst}}$ — i.e. the instance
  witnesses that the AI-claimed $\rho$ (if the proof were valid) cannot
  hold. Since the AI claim is unverified, ANY Gate-C outcome is
  COMPUTATIONAL-EVIDENCE about the claim, never about the conjecture.

---

## 6. Evidence labels (FIXED usage)

- **MACHINE-VERIFIED**: the row's exactness chain completes (exhaustive
  pattern minimization OR complete-enumeration instance; exact $\mathbb{F}_q$
  verification of the witness; exact-rational dual closure where an LP is
  used). Includes the 2-anchor-complete rows.
- **COMPUTATIONAL-EVIDENCE**: solver-located optimum without exhaustive
  global-min proof; includes all $\Lambda$-probe rows and all sampled-$t$
  rows at large $q$.
- **HUMAN-AUDITED**: none at campaign time (this doc + code review by Main
  may add it later; not asserted by this agent).
- **CITED-DEPENDENCY**: TR26-150 (primary), flint/Highs/sympy (tools),
  Polishchuk–Spielman context.
- **OPEN**: rows not covered by the budget.
- **FAILED**: pipeline crashes/abortions with cause.
- Provenance: [REPRODUCED] / [DERIVED] / [REPORTED] / [INFERENCE] tagged
  per quote/claim above. **[REPORTED]** items: none in this pre-statement
  (no second-hand load-bearing facts used).

---

## 7. Environment and resource policy (FIXED)

- Interpreter: `/Users/jinleic/jinleic-workspace/cs/.venv/bin/python`
  (3.14.3). Tools: python-flint 0.9.0 (`fmpq`, `fmpz`, `fq_default` if
  needed), sympy 1.14.0 (`primitive_root`), python-sat 1.9.dev15, highspy
  (HiGHS), z3-solver 5.1.0, mpmath 1.3.0, numpy 2.5.2. No new installs
  anticipated; any install recorded in the campaign log.
- Every run: `nice -n 10`, and
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1`. One heavy job at a time;
  bounded searches with explicit stop conditions; no daemonized loops.
- Campaign directories: `campaigns/` frozen after creation
  (`campaigns/<ISO-UTC>Z_<short-uuid>/` containing `manifest.json`,
  per-instance `*.json`, `run.log`, `codehash.txt`); rehearsals in
  `scratch/` or `campaigns-smoke/` only.

---

## 8. Re-scoping log (append-only)

- 2026-08-29: original pre-statement committed. No re-scoping yet.
 
## 10. Placement note for the re-scoping addendum

The correction in Section 9 is intentionally retained at its original
insertion point so that it is adjacent to the affected enumeration text.
The later original sections remain historical fixed policy text; where they
say “exact” without the finite-window qualification, Section 9 controls.
 
## 11. Pilot execution budget (fixed before the first run)

The campaign will run the complete weight-$\le3$ support census for every
feasible $(t_1,t_2,t_3)$ profile at $(q,s,\eta)=(31,(2,3,5),1/32)$.
Because the higher-$t$ kernel intersections contain thousands of support
representatives, the fixed first-run optimization budget is 120 CPU seconds:
the MILP plus exact cheaper-support closure is run for the full candidate
class only for the baseline profile $t=(1,1,1)$. For every other profile the
support census is retained, but no ratio is reported unless a witness is
explicitly completed within this budget. This is a pre-registered resource
cut, not a post-result threshold change; omitted profile/value searches are
OPEN.
 
## 12. Unit-anchor feasibility correction

My pre-statement's illustrative Anchor 1 uses $q=7,s=(2,3,4)$, but $4\nmid
(7-1)$, so a multiplicative subgroup of order $4$ does not exist in
$\mathbb F_7^\star$. This is a unit-grid typo, not a campaign parameter.
Before relying on the anchor, the executable unit check substitutes
$q=13,s=(2,3,4)$ (all three orders divide $12$), keeps $t_i=s_i$, and
records the substitution in its smoke artifact. Coprimality remains disabled
for this unit under the pre-statement's unit convention. The pilot campaign
parameters are unchanged. The one-direction equality tested in this anchor
was my derived scratch premise; its failure is not attributed to the README
or to TR26-150.
 
## 13. Gate-B finite sweep scope (fixed before execution)

To obtain a fixed-$q$ balance probe without claiming an asymptotic
conclusion, use $q=61$, $\eta=1/32$, $t=(1,1,1)$, identity $\Lambda$, and
the two valid coprime triples
$s=(2,3,5)$ and $s=(3,4,5)$. Their balance parameters under
$\beta=\max_i s_i/\min_i s_i$ are respectively $\beta=5/2$ and
$\beta=5/3$. For each triple, exhaustively census all point supports of
weight at most $W_{\rm cap}=3$ and optimize the normalized reduced-basis
candidate class; each accepted delta has the same exact Fq plus cheaper
line-support closure as Gate A. This is a finite-window Gate-B comparison at
fixed $q,\eta,t,\Lambda$; all heavier supports, non-basis values, other
$t,\Lambda$, and any claim about monotone $\rho(\beta)$ remain OPEN.
 
## 14. Status addendum (2026-08-30) — machine-verifiable finite rows

Section 13's Gate-B scope was executed and extended. The frozen campaign is
`campaigns/2026-08-30T01-04-49Z_exacttable/`. For each of its nine instances
the best-in-window ratio was equipped with a complete per-witness proof:
every line-support tuple of cost $<\delta$ was exactly checked absent and
$\ge1$ tuple of cost $=\delta$ was checked present, closing $\delta(M^*)$ for
the witness. Ratios so obtained are
$(13,(2,2,4))$: $1/2,\ 3/8,\ 1/4$ at
$t=(1,1,1),(1,1,2),(1,1,4)$;
$(31,(2,3,5)),(31,(2,3,10)),(31,(3,10,2)),(61,(2,3,5))$: $3/5$;
$(31,(2,6,5)),(61,(3,4,5))$: $1$.
These are exact, machine-verified per-witness quantities inside the stated
weight-$\le3$ window; they are NOT global $\rho_{\rm inst}$ values, and the
balance conclusion is that no collapse at $\beta\in\{5/3,5/2\}$ was found
(heavier-support global behavior remains OPEN).
