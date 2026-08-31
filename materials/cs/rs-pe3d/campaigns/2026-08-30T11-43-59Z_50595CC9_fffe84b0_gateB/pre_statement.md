# PRE-STATEMENT — Gate B — rs-pe3d — fixed-(q, eta) balance sweep

Committed: 2026-08-30T11:43:59Z (UTC), by agent `RsPe3dGateB`, BEFORE any gate-B
computation run. This file, `../pre_statement.md` (parent, Sections 1–14),
`scratch/gateB_premise_worksheet.json`, and `scratch/gateB_anchor_plan.json`
are the anti-hindsight contract. Code hash at commitment time
(`src.campaign._code_hash()`):
`fffe84b00b9ca48ba5b3a6970b5523c9063bb61a559a88f81328a501e7b8bd56`.

When this directory was created it contained ONLY this file. No sweep numbers
existed anywhere in my session at commitment time. Any later edit to this file
must be a dated addendum; the original text stays visible.

---

## GB1. What balance means in TR26-150 — owner-read verbatim quotes

Source owner-read this session, first-hand, full text, PDF taken from
`https://eccc.weizmann.ac.il/report/2026/150/download` (title page:
TR26-150 | 19th August 2026, Bafna & Vyas, *Private PCPs from Product
Expansion*). Quotes re-checked against my own extraction (`/tmp/tr26-150.txt`,
pdftotext -layout) line by line. Provenance [REPRODUCED] unless tagged.

### GB1.1 Definition 4.1 (Product expansion) [REPRODUCED]

> **Definition 4.1 (Product expansion).** A tuple $(C_i)_{i\in I}$, where
> $I \subseteq [k]$, is $\rho$-product-expanding if every
> $M \in \sum_{i\in I} L_i(C_i)$ has a decomposition $M = \sum_{i\in I} M_i$,
> with $M_i \in L_i(C_i)$, such that
> $\mathrm{wt}(M) \geqslant \rho \sum_{i\in I} s_i \ell_i(M_i).$

with the paper's setup: $s_i = |S_i|$, $\Omega = \prod_i S_i$, an $i$-axis line
fixes every coordinate except the $i$-th, $\ell_i(M_i)$ counts nonzero $i$-axis
lines of $M_i$, $\mathrm{wt}$ is Hamming weight on $\Omega$. This is exactly
the quantity computed on disk since Gate A ($\rho_{\rm inst}$ in the parent
pre-statement §2).

### GB1.2 Conjecture 4.2 [REPRODUCED — quantifiers are load-bearing]

> **Conjecture 4.2 (Product Expansion of Punctured Reed–Solomon Codes).** For
> every $\eta \in (0, 1)$, balance parameter $\beta \geqslant 1$, and integer
> $r \geqslant 2$, there is $\rho = \rho(r, \eta, \beta) \in (0, 1]$ such that
> the following holds. Let $q$ be prime, and let
> $S_1, \dots, S_r \leqslant \mathbb{F}_q^\times$ be multiplicative subgroups
> of pairwise coprime orders $s_i = |S_i|$ satisfying
> $$ \beta^{-1} \leqslant \frac{s_i}{s_j} \leqslant \beta \quad (i, j \in [r]). $$
> For integers $t_i$ with $1 \leqslant t_i \leqslant (1 - \eta)s_i$, and
> invertible diagonal matrices $\Lambda_i$ indexed by $S_i$, the tuple
> $(\Lambda_i \mathrm{RS}(S_i, t_i))_{i=1}^r$ is $\rho$-product-expanding.
>
> Note that, without loss of generality, we may assume that $\rho(r, \eta,
> \beta)$ is nonincreasing in $r$.

Quantifier structure, made explicit [INFERENCE from the displayed text]:

- $\forall \eta \in (0,1)\ \forall \beta \ge 1\ \forall r \ge 2\
  \exists \rho = \rho(r,\eta,\beta) \in (0,1]$;
- then $\forall$ prime $q$, $\forall$ subgroups with pairwise coprime orders
  satisfying the two-sided ratio bound, $\forall t_i \in [1, (1-\eta)s_i]$,
  $\forall$ invertible diagonal $\Lambda_i$: the tuple is
  $\rho$-product-expanding **per Definition 4.1**.
- **Balance is a HYPOTHESIS inside the conjecture**, not a conclusion: the
  conjecture's conclusion holds only on the balanced family
  $\beta^{-1} \le s_i/s_j \le \beta$. The balance parameter is an input.
  The conjecture is SILENT (omits any claim) for triples outside the bound.
- The conjecture as stated is $\exists\rho$ per fixed $(r,\eta,\beta)$ over an
  infinite family of finite instances — no finite instance can refute it
  (parent pre-statement §1.9 falsifiability stance, retained).

### GB1.3 Theorem 2.1 (proved 2-D, balance hypothesis explicit) [REPRODUCED]

> **Theorem 2.1 (2d Product expansion).** Fix $\varepsilon > 0$. There exists
> $\rho = \rho(\varepsilon) > 0$ such that the following holds. Let
> $S_i \leqslant \mathbb{F}_q^\times$ have coprime orders $n_i = |S_i|$ with
> $|S_1| \approx |S_2|$ and let $C_i = \mathrm{RS}(S_i, d_i)$, where
> $d_i \leqslant (1-\varepsilon)n_i$. Every matrix $M$ that is a sum of
> $C_1$-columns and $C_2$-rows admits a decomposition $M = A + B$ such that
> $$ \mathrm{wt}(M) \geqslant \rho \left( n_1 \cdot \#\{\text{nonzero columns of } A\} + n_2 \cdot \#\{\text{nonzero rows of } B\} \right). $$

Here "$|S_1| \approx |S_2|$" is INFORMAL. The paper's own formalization of the
same hypothesis appears in:

### GB1.4 Theorem 4.11 (proved, author-verified, formal balance) [REPRODUCED]

> **Theorem 4.11 (Two-dimensional product expansion of subgroup Reed–Solomon
> codes).** Fix $\varepsilon \in (0, 1)$ and $K \geqslant 1$. There is a
> constant $\rho = \Omega(\varepsilon^6/K^3) \in (0, 1]$ with the following
> property. Let $q$ be prime, and let $S_1, S_2 \leqslant
> \mathbb{F}_q^\times$ be multiplicative subgroups of coprime orders
> $n_i = |S_i|$ satisfying $K^{-1} \leqslant n_1/n_2 \leqslant K$. For
> integers $k_i$ satisfying $1 \leqslant k_i \leqslant (1-\varepsilon)n_i$,
> let $C_i = \mathrm{RS}(S_i, k_i)$. Then $(C_1, C_2)$ is
> $\rho$-product-expanding. The same conclusion holds after independently
> rescaling the coordinates of $C_1$ and $C_2$ by nonzero field elements.

### GB1.5 Balance in the paper's own parameter choices [REPRODUCED]

- Lemma 5.4, item 2 (Domain sizes): "$\frac{1}{4} < \frac{s_i}{s_j} < 4$
  $(i, j \in [k])$" — i.e. the paper instantiates the conjecture at
  **$\beta = 4$** ("4-balanced", Corollary 4.5 usage
  $\rho(k) = \rho(k, 1/128, 4)$).
- Section 6 discussion: "Privacy comes from the degree slack left after fixing
  a polynomial on $D$, together with our product expansion conjecture, whose
  hypotheses require the side lengths of $\Omega$ to be pairwise coprime and
  **balanced**."
- The proved 2-D dependence is $\rho = \Omega(\varepsilon^6/K^3)$ — balance
  enters the 2-D constant as $K^{-3}$; the hypothesis itself must hold.

### GB1.6 The AI candidate claim [REPRODUCED]

Paper (intro): "We have a candidate AI-generated proof of this conjecture
which gives $\rho = \varepsilon^{\exp(k)}$, but it has not been verified by
the authors." The linked AI PDF (`https://mitalibafna.github.io/higher-dimension-
product-expansion.pdf`, owner-read this session, header: "Disclaimer: This is
an AI-generated proof that has not been verified.") states its own Theorem with
a **$K$-comparability hypothesis** "$K^{-1} \le n_i/n_j \le K$ for all
$i, j \in [d]$" and proves
$\rho = \dfrac{1}{2dK^d}\left(\dfrac{\varepsilon}{(C_0 d K)^{10d}}\right)^{3d\,6^{d-2}}$
— **balance-dependence is carried inside the AI claim too**, through the
$K^{-\exp}$ factors. Gate C (done, on disk) audited the $\varepsilon$-form at
$\eta$-grid points; gate B probes the $\beta$-dependence that BOTH the 2-D
theorem and the AI proof carry but that the *informal* higher-D statement
(Conjecture 1.4: "$\rho$ ... only a function of $\varepsilon$ and $k$") drops.
This distinction is the entire point of gate B.

---

## GB2. The imbalance parameter — fixed BEFORE any computation

The conjecture's own scalar is a two-sided bound:
$\beta_{\rm conj}(s) = \max_{i,j} s_i/s_j = (\max_i s_i)/(\min_i s_i)$ over the
triple. This is a HYPOTHESIS threshold: $\beta_{\rm conj} \le \beta$ means the
conjecture's conclusion applies; $\beta_{\rm conj} > \beta$ means the
conjecture is silent. Gate B therefore sweeps the hypothesis parameter itself.

**Fixed definition [COMMITTED BEFORE RUNNING]: for a subgroup-order triple
$s = (s_1, s_2, s_3)$ define**

$$ \beta(s) \;=\; \frac{\max_i s_i}{\min_i s_i}. $$

- $\beta(s) = 1$ iff all orders equal (the perfectly balanced case; note
  pairwise coprime + equal forces $s_i = 1$, which the construction excludes —
  so perfect balance is not realizable in this family; the minimum achievable
  $\beta$ grows with the number of pairwise-coprime orders $> 1$).
- Specializes correctly to the 2-D theorem's hypothesis: for $r = 2$,
  $\beta(s) = \max(n_1, n_2)/\min(n_1, n_2)$, and
  $K^{-1} \le n_1/n_2 \le K \iff \beta(s) \le K$. So Theorem 4.11's $K$ is
  exactly this $\beta$, reading the 2-D hypothesis as
  "$\beta(s_1, s_2) \le K$".
- Matches Conjecture 4.2's quantifier: the conjecture's hypothesis
  $\beta^{-1} \le s_i/s_j \le \beta\ \forall i,j$ is equivalent to
  $\beta(s) \le \beta$.
- Anti-hindsight note: this is the conjecture's OWN parameter read off its
  displayed hypothesis, not a quantity invented after seeing data. Alternatives
  (e.g. $\log$-spread $\beta_{\log} = \max \log s_i / \min \log s_i$, or
  $\max_i s_i / \mathrm{median}$) are computable post hoc from the same triples
  if adjudication needs them; they are secondary and no threshold may be
  defined through them first.

---

## GB3. Fixed-(q, eta) anchor-sweep design (what "hold (q, eta) fixed" means here)

The contract says: hold $q, \eta$ fixed, vary the subgroup triple across the
pre-registered imbalance range. The tension: at a single prime $q = 61$ there
are exactly TWO coprime triples ($\beta = 5/2$ and $\beta = 5/3$ — a
two-point $\beta$-grid). A wider $\beta$ range at fixed $q$ does not exist
within a single prime. Pre-registered resolution (fixed now, before running):

- **Primary anchor sweep (fixed $q = 61$, $\eta = 1/32$, $t = (1,1,1)$,
  $\Lambda = \mathrm{Id}$):** triples $(2,3,5)$ [$\beta = 5/2$] and
  $(3,4,5)$ [$\beta = 5/3$]. This is the contract's literal fixed-$(q,\eta)$
  sweep and it reproduces both frozen gate-A rows as its anchors.
- **Grid-extension tier (each point still a full fixed-$(q,\eta,t,\Lambda)$
  instance; the tier as a whole is labeled as a cross-prime comparison, not a
  single-prime sweep):** add pre-registered instances at in-scope primes
  $q \in \{13, 31, 241\}$ frozen in the gate-A-style instance list below.
  Rationale committed NOW: the dichotomy (graceful vs collapse) is about the
  sign/structure of $\rho(\beta)$; the fixed-$q$ tier gives exact two-point
  within-prime ordering; the cross-prime tier extends the $\beta$ range to
  $[5/3, 8/3] \cup \{3, 4, 8/3, 16/3, 6, ...\}$ as the instance list fixes.
  Every cross-prime row will be labeled separately and the fixed-$q$ rows will
  be the primary evidence.

### GB3.1 The full pre-registered instance grid (COMMITTED)

All instances: $\Lambda_i = \mathrm{Id}$, weight-$\le 3$ support census,
normalized reduced-basis candidate class, exhaustive cheaper-line-tuple closure
per accepted witness (the frozen `optcheck._proof_delta` procedure), exact
rational ratios.

Group 1 — fixed-(q, eta) anchor sweep, q = 61, eta = 1/32 (contract's core):

| id | q | s | t | beta(s) | role |
|----|---|---|---|---------|------|
| B01 | 61 | (2,3,5) | (1,1,1) | 5/2 | gate-A row reproduce (anchor 1 of 2) |
| B02 | 61 | (3,4,5) | (1,1,1) | 5/3 | gate-A row reproduce (anchor 2 of 2); more balanced |

Group 2 — cross-prime beta extension, same eta = 1/32, t = (1,1,1),
Lambda = Id (pre-registered for range, labeled cross-prime):

| id | q | s | t | beta(s) | source of q |
|----|---|---|---|---------|-------------|
| B03 | 31 | (2,3,5) | (1,1,1) | 5/2 | gate-A frozen |
| B04 | 31 | (2,3,10) | (1,1,1) | 5 | gate-A frozen |
| B05 | 31 | (2,6,5) | (1,1,1) | 3 | gate-A frozen |
| B06 | 31 | (3,10,2) | (1,1,1) | 5 | gate-A frozen |
| B07 | 13 | (2,2,4) | (1,1,1) | 2 | gate-A frozen |
| B08 | 13 | (2,2,4) | (1,1,2) | 2 | gate-A frozen |
| B09 | 13 | (2,2,4) | (1,1,4) | 2 | gate-A frozen |
| B10 | 241 | (3,5,8) | (1,1,1) | 8/3 | new q, wide-beta end |
| B11 | 241 | (3,5,16) | (1,1,1) | 16/3 | new q, max beta in scope |
| B12 | 241 | (2,3,5) | (1,1,1) | 5/2 | new q, balanced-side control |
| B13 | 241 | (3,4,5) | (1,1,1) | 5/3 | new q, most balanced in scope |

(In-scope primes with $\ge 2$ coprime triples: 61 has $\{(2,3,5),(3,4,5)\}$;
241 has $\{(2,3,5),(3,4,5),(3,5,8),(3,5,16)\}$; 421 has 19 triples but
$N = s_1 s_2 s_3$ up to $210$ and window cost makes several of them
infeasible — B10–B13 use 241's full four-triple census, which is complete,
so 241 is itself a full fixed-$q$ mini-sweep. 421's 19-triple census is NOT
pre-registered (budget, GB5); recorded OPEN.)

Wait — correction to B13 row: (3,4,5) has $\gcd(4, 5)=1$, $\gcd(3, 4)=1$,
$\gcd(3,5)=1$ — valid, divides $240 = 3 \cdot 8 \cdot 5$: yes
($3 \mid 240$, $4 \mid 240$, $5 \mid 240$). OK. And (2,3,5) likewise. Good.

### GB3.2 Why these exact instances

- B01–B02: the contract's literal anchor sweep at fixed $q = 61$; also the
  pre-run anchors (must reproduce 3/5 and 1/1 before extending).
- B03–B09: already frozen gate-A values; recomputed with THIS campaign's code
  as a second-route check (second-route rule: two structurally independent
  histories for the same numbers — the frozen gate-A campaign and this one).
- B10–B13: 241 gives the widest $\beta$ range among in-scope primes with a
  complete multi-triple census; $16/3$ is the pre-registered extreme.

### GB3.3 The t-profile policy (FIXED)

$t = (1,1,1)$ everywhere in gate B. Reason fixed now: $t$ interacts with
$\delta(M)$'s composition space and multiplies cost; conjecture quantifies
$\forall t_i$ but gate B's question is the $s$-dependence at fixed slack, and
$\eta = 1/32$ with $t_i = 1$ is feasible at every listed instance. Any $t$
variation recorded: none (OPEN).

### GB3.4 The Lambda policy (FIXED)

$\Lambda = \mathrm{Id}$ on all axes (parent pre-statement headline policy).
The $\forall \Lambda_i$ quantifier is NOT swept in gate B; all rows are
$\Lambda = \mathrm{Id}$ and this is part of the swept-set statement.

---

## GB4. The exact quantity per instance (unchanged from parent pre-statement §2)

$\rho_{\rm inst}^{\rm upper}$ = min over the campaign's candidate class
(weight-$\le 3$ supports, normalized reduced-basis vectors) of
$\mathrm{wt}(M)/\delta(M)$, with $\delta(M)$ per-witness-closed exactly by the
exhaustive cheaper-tuple sweep. This is what twin gate-A campaigns did; gate B
inherits the identical window and closure and adds only new instances. No
global-$\rho$ claim is made anywhere (Rule 7, GB7).

Sanity thresholds (inherited, parent §2): any $\rho_{\rm inst} > 1$ or
$= 0$ is a bug — STOP, do not record.

## GB5. Budget (FIXED before running)

- Per-candidate MILP limit: 120 s (as gate-A table campaign).
- Instance-level outer deadline: 4 h wall; the sweep checkpoints to
  `scratch/gateB_checkpoint.json` after EVERY instance (write-first protocol;
  two predecessors lost work to mid-run failures with nothing on disk).
- Phase order (COMMITTED):
  1. **Anchor phase first**: unit anchor $(5,(2,2,2),(1,1,1))$ exhaustive
     $\rho = 1/3$ check, then B01, B02 — all three must reproduce frozen
     values (1/3, 3/5, 1/1) with THIS campaign's code before any B03+ row is
     computed. If an anchor fails: STOP, debug, escalate to Main, no B03+ run
     at all.
  2. Cross-prime tier B03–B09 (cheap, frozen-value reproduction; a mismatch
     with frozen gate A stops everything and triggers escalation per the
     second-route rule — disagreement between two routes on the same quantity
     is an anchor failure).
  3. New-q tier B10–B13.
- Skip conditions: none pre-registered (all 13 rows are committed; if the
  budget runs out mid-tier, the checkpoint records exactly which rows are
  done and the rest stay OPEN).
- Concurrency: single-threaded-ish per parent §7 (`OMP_NUM_THREADS=1` etc.),
  `nice -n 10`.

---

## GB6. Adjudication protocol (FIXED before seeing any B03+ value)

Data: pairs $(\beta(s), \rho)$ across the grid. Verdict rules:

- **Graceful degradation:** across the swept grid, $\rho$ stays bounded away
  from the universal floor (i.e. $\rho \ge c > 0$ with $c$ not shrinking as
  $\beta$ grows within the grid) — conjecture plausible without balance
  hypothesis, at grid resolution.
- **Collapse at $\beta^*$:** there is a last grid point $\beta^-$ with
  $\rho$ at its balanced-typical value and a next grid point $\beta^+$
  (the smallest swept $\beta$ above $\beta^-$) where $\rho$ drops by a
  pre-registered break factor — defined now as: $\rho(\beta^+) < \frac{1}{3}
  \cdot \min$-over-balanced-grid ($\beta \le 5/2$) — then the conjecture NEEDS
  a balance hypothesis. $\beta^*$ is then located only to grid resolution:
  reported as $\beta^* \in (\beta^-, \beta^+]$, uncertainty = grid spacing.
- **Neither:** if $\rho(\beta)$ is non-monotone in a way that fits neither
  shape (e.g. rises then falls), or falls but stays far above the floor, or
  the window effect (wt/delta granularity) dominates the trend — report as a
  third outcome with the full data; do NOT force the dichotomy.
- Any verdict that accuses the conjecture statement (i.e. "collapse"): freeze
  campaign, hub `Main` with exact discrepancy + campaign path, and STOP that
  thread before writing any claim into any README/index.

---

## GB7. Rule-7 statement (fixed now; restated verbatim in the final report)

The swept set, exactly:
$$ \{ (q, s, t, \eta, \Lambda) : \eta = 1/32,\ t = (1,1,1),\ \Lambda = \mathrm{Id},\ (q,s) \in \text{GB3.1 table} \} $$
with witnesses restricted to weight-$\le 3$ supports and the normalized
reduced-basis candidate class of $V \cap \mathbb{F}^S$, $\delta$ exactly
closed per accepted witness.

NOT swept (stated plainly in the report): heavier supports (weight $\ge 4$),
non-basis field values inside a support (dim-$>1$ intersections contribute
only their reduced basis), $t \ne (1,1,1)$, $\eta \ne 1/32$, non-identity
$\Lambda$, primes other than $\{13, 31, 61, 241\}$, the 421 nineteen-triple
census, continuous $\beta$ (only the discrete grid $\{2, 5/3, 5/2, 3, 8/3, 5,
16/3\}$ intersected with realizable triples), and every asymptotic statement
in $\beta$. All gate-B conclusions are statements about THIS swept set;
$\beta^*$, if seen, carries grid resolution as its uncertainty.

---

## GB8. Escalation and labels

- Labels per parent §6: rows fully closed in-window = MACHINE-VERIFIED
  (finite-window); anything else COMPUTATIONAL-EVIDENCE; float rows
  COMPUTATIONAL-EVIDENCE (none planned — everything is exact rationals);
  budget-cut rows FAILED/OPEN with cause.
- Escalation trigger: per assignment — a collapse verdict, or ANY numeric
  disagreement with a frozen gate-A row or with TR26-150's statements, goes to
  Main FIRST, before any README claim.
- The unit-anchor premise lesson is inherited from
  `scratch/anchor1_counterexample.json`: no derived premise goes untested on a
  feasible substitution. Gate B's own candidate premises (window stays
  exhaustively closable at extreme imbalance) are themselves tested on the
  B11 instance $(3,5,16)$ before any interpolation claim.

---

## GB9. ADDENDUM — pre-registered prediction for the in-flight row B11

Committed 2026-08-30T13:39Z (UTC), by agent `RsPe3dGateB`, **while the B11
census was running at 45.6% (1,050,000 of 2,304,200 supports) and BEFORE any
B11 ratio existed anywhere** — no B11 candidate had been optimized, the
`scratch/gateB_b11_result.json` file did not exist. The original GB1–GB8 text
above is unmodified. This section adds a falsifiable prediction so that the
pattern below is a TEST rather than a post-hoc story.

### GB9.1 The observed pattern on the twelve completed rows

All twelve completed rows (exact rationals, per-witness $\delta$ closed):

| id | q | s | t | beta(s) | rho (exact) | contains both 2 and 3? |
|----|---|---|---|---------|-------------|------------------------|
| B01 | 61 | (2,3,5) | (1,1,1) | 5/2 | 3/5 | yes |
| B02 | 61 | (3,4,5) | (1,1,1) | 5/3 | 1 | no |
| B03 | 31 | (2,3,5) | (1,1,1) | 5/2 | 3/5 | yes |
| B04 | 31 | (2,3,10) | (1,1,1) | 5 | 3/5 | yes |
| B05 | 31 | (2,6,5) | (1,1,1) | 3 | 1 | no |
| B06 | 31 | (3,10,2) | (1,1,1) | 5 | 3/5 | yes |
| B07 | 13 | (2,2,4) | (1,1,1) | 2 | 1/2 | no |
| B08 | 13 | (2,2,4) | (1,1,2) | 2 | 3/8 | no |
| B09 | 13 | (2,2,4) | (1,1,4) | 2 | 1/4 | no |
| B10 | 241 | (3,5,8) | (1,1,1) | 8/3 | 1 | no |
| B12 | 241 | (2,3,5) | (1,1,1) | 5/2 | 3/5 | yes |
| B13 | 241 | (3,4,5) | (1,1,1) | 5/3 | 1 | no |

**Pattern P (hypothesis, NOT a claim):** at $t=(1,1,1)$, $\Lambda=\mathrm{Id}$,
$\{2,3\} \subseteq \{s_1,s_2,s_3\} \iff \rho = 3/5$ in this window. It holds on
all twelve rows with no exception. Provenance [DERIVED from my own rows];
suggested by Main; explicitly to be **killed, not confirmed** — Main records
that structural hypotheses offered this session are 0-for-8.

**Exclusion E (a conclusion, not a hypothesis):** no monotone function of
$\beta$ can fit the completed rows. At $q = 31$ the ratios are INTERLEAVED in
$\beta$: $\beta = 5/2 \to 3/5$, $\beta = 3 \to 1$, $\beta = 5 \to 3/5$ — a
value of $1$ sits strictly between two values of $3/5$ as $\beta$ increases.
Same shape at $q = 241$: $5/3 \to 1$, $5/2 \to 3/5$, $8/3 \to 1$. Since these
are exact rationals with per-witness closure, monotone-in-$\beta$ is EXCLUDED
for $\rho^{\rm window}$ on the swept set — this holds regardless of what B11
does.

### GB9.2 The two opposing predictions for B11 (committed)

B11 $= (q,s,t) = (241,(3,5,16),(1,1,1))$, $\beta = 16/3$ — the largest $\beta$
in the grid; its orders contain $3$ but NOT $2$.

- **Under pattern P:** $\rho(\mathrm{B11}) = 1$ (no $2$ in the orders).
- **Under any $\beta$-driven degradation story:** B11 must be the LOWEST ratio
  in the whole grid, i.e. $\rho(\mathrm{B11}) < 1/4$ (strictly below B09's
  $1/4$, the current grid minimum), since $\beta = 16/3$ is the most
  imbalanced point swept.

These are incompatible on the same row, so B11 discriminates. Pre-committed
readings of the three possible outcomes:

1. $\rho = 1$ → P survives one real test; $\beta$-degradation is dead as a
   driver at the top of the swept range (together with E).
2. $\rho = 3/5$ → P is FALSIFIED (a row without $2$ giving $3/5$); the
   $\beta$ story survives only in the weak sense that the top-$\beta$ point is
   not maximal-$\rho$; still not monotone, by E.
3. anything else ($4/5$, $2/3$, $1/2$, $< 1/4$, …) → BOTH P and the simple
   $\beta$ story die; report the number and treat the shape question as OPEN.

Reporting rule (committed): **the number first, the interpretation second**;
the verdict sentence must name which of the three branches occurred.

### GB9.3 Verdict-language guard

Nothing in GB9 changes GB6. In particular a "collapse at $\beta^*$" verdict
still requires the GB6 break-factor condition, still gets escalated to Main
before any README wording, and P is a statement about MY window, never about
TR26-150.
