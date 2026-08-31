# PRE-REGISTERED DESIGN — Gate B: failure rate of Apon's `Δ_{p,q} ≠ 0` genericity hypothesis

**Committed 2026-08-30T17:35Z (UTC), BEFORE the first computation of this campaign.**
Agent `MceliececelGateB`. Inherits the Gate-A conventions of
`campaigns/2026-08-30T14-09Z_57200ADD/pre_statement_with_addendum2.md` and the
addendum-2 label policy. This document is committed to disk (git add) before
the first search execution; nothing below is chosen after seeing data (repo
rule 16, domain-shopping prohibited).

## 0. Object measured.

Apon's Lemma 9 (ePrint 2026/1810, p. 11) needs a pair p<q with
`Δ_{p,q} = f_p f_q' + f_q f_p' != 0` to prove `A ↦ T(A) mod EF` injective; his
§3.6 "Scope of the result" states verbatim: *"If every Δ_{p,q} is zero, Lemma 9
does not apply."* **That is the hypothesis whose failure rate we measure.**
Gate A observed the guard passing on 13/13 census instances — in-sample, chosen
for the census, not evidence of a rate. Vedenev (1747) imposes no such
condition, but his v2/§9 costs (q^{m(c-2)}) do not cover the degenerate branch
either; see Reports/interpretation. no disputant publishes a rate.

Instance object (inherited verbatim from Gate A, code `instance.py` of
campaign `57200ADD`): given (m, n, t, seed) build the binary Goppa instance on
E=F_{2^m}: G a seeded random monic irreducible of degree t (Rabin-tested);
support = prefix of a seeded shuffle of the full enumeration of E;
Π = ∏(Z−a_i); λ_i = G(a_i)²/Π′(a_i); Y = first k rows of the RREF of the
binary parity-check nullspace (deterministic function of column order);
F = (f_1..f_k), deg ≤ D = n−2t−1, with λ_i f_j(a_i) = Y[j,i].

**`Δ_{p,q} = f_p f_q' + f_q f_p'` (char-2 Wronskian; coefficients in E;
"≠0" = the polynomial is not identically zero, decided EXACTLY by
`gf.ptrim(...) != []`).**

## 1. Degeneracy of the instance = (definition, fixed in advance)

An instance (m,n,t,seed) is **DEGENERATE** iff

> for ALL pairs p<q in {1..k}:  `Δ_{p,q} = 0` EXACTLY (zero polynomial in E[Z]).

One failing pair suffices to be nondegenerate (Apon needs only one).

## 2. Search space (pre-registered)

Two layers, reported SEPARATELY; never merged into one number.

- **(A) Sampled layer** — models "random Goppa public key, then apply Gate-A
  instance construction". Seeds only; G, support order, Y, F all determined by
  the seed stream. **Boundary t = 1 pre-registered here:** at m=6, t=1 the
  census build *aborts* (G's unique root must miss the n<n_full support —
  at n = 2^m, i.e. full support, no degree-1 irreducible G can avoid the
  support: root ∈ E = support, so G(a_i)=0 for the root a_i... BUT the build
  asserts G(a_i)≠0 on the support and the assertion fires. Therefore **t=1
  instances do not exist** at n = 2^m and t=1 is a nonce-boundary, not
  sampled). Consequently sampled t ranges:

  m=6:  t ∈ {2,3}  (n=64, k=60 or 58, D = 59 or 57)
  m=7:  t ∈ {2,3}  (n=128, k=122 or 120, D = 123 or 121)
  m=8:  t ∈ {2,3}  (n=256, k=250 or 248, D = 251 or 249)
  m=9:  t ∈ {2,3}  (n=512, k=506 or 504, D = 507 or 505)
  m=10: t = 3      (n=1024, k=994, D = 1017)
  m=11: t = 3      (n=2048, k=2024, D = 2041)

  N = 25 seeds per (m,t) cell (fixed stream: seeds 1000..1024 per cell,
  advance-on-guard-failure substitution logged, same rule as Gate A §6).
  Cell sizes chosen for ≥25 dominance-decaying CI half-width; the floor
  budget (rule 2 of the main text) is the believer's assurance that a
  grossly common failure (say, 10%) is caught with overwhelming probability.

- **(B) Adversarial layer** — hand-built configurations chosen to stress the
  hypothesis structurally, FIXED before any of their runs:

  * (B1) *t=2 at every m ∈ {6..11}* (smallest nontrivial t). Any t=2 Goppa
    poly factors over E as (Z−β)(Z−β^{2^{m-1}}) + ... — i.e. degree-2
    irreducible = quadratic — with the root {β, β^{2^{m/2}}} structure when m
    even. We test whether the *pair* (f_p, f_q) can be proportional because
    g has G-fold structure. 25 seeds per cell.
    — Wait, this is also covered under sampled layer A above; (B1) is the
    **force-degenerate attempt** at those cells: we *search* seeds until we
    find a degenerate instance or exhaust the stream (cap 200 seeds/cell,
    recorded as "not found within cap").
    — Re-stated cleanly in §2-fix below.
  * (B2) *k=1 and k=2* instances: with k=1 there are no pairs at all, so the
    hypothesis is VACUOUSLY... ruled by the strongest possible vacuity note,
    not a rate event. k=2: exactly one Δ; if λ-interpolation keeps f_1, f_2
    equal (up to scalar) the instance degenerates. Pre-registered built
    directly at (m, n=2^m, t=n−k over m — feasibility k=n−mt>0 forces
    t=(n−k)/m ≤ (2^m−k)/m, so k=2 ⇒ t ≈ (2^m−2)/m, large t). **Budget:**
    m=6 (t=10), m=7 (t=18), m=8 (t=31), m=9 (t=63), m=10 (t=102), m=11
    (t=186). 25 seeds each.
  * (B3) *repeated-order / structured support orderings*: the SAME support set
    with 25 different seeded orderings (orderings ∈ {seeded shuffle,
    reversed shuffle, "sort-by-G-value, seeded tiebreak"}), at
    (m,n,t)=(6,64,3) and (7,128,6). 25 shuffles + 25 reversed + 25 sorted = 75
    orderings per (m,t), all with G fixed from the first 25-seed stream.
  * (B4) *small-field / small-t boundary sweep*: (m, n, t) with the *smallest
    growing* families: (m,n,t) ∈ {(4,16,t): t∈{1,2,3}}, exhaustively census
    all irreducible G and all support orders — no seed at all, FULL
    enumeration where the space is small enough. For m=4: |E|=16. t=1:
    G linear, root must avoid support; support n<16 (t=1 ⇒ k=n−4·1>0 ⇒
    n ≥ 5). Witness: (4, 5..15, 1): enumerate G ∈ {Z+β: β∈E}, all supports
    n ∈ {5..15} (subsets of size n — sampled, not exhaustive), all
    orderings. For n=15 one β fails (the excluded element IS the root —
    when β ∈ E∖support, G has a root off-support → valid instance), 1 of
    16 β; for n=15 exactly the excluded element works. Rate over the
    enumerated family is REPORTED but this is a boundary, not a rate.
  * (B5) *k=2 forced-degenerate attempts via Y-manipulation*: at (6,64,58)
    (m=6, t=3), take the Gate-A construction, then REPLACE the RREF-selected
    Y with {row, row ⊕ row'} (a binary combination that keeps Y a valid
    generator of the same code) such that the two interpolants are designed
    proportional. If ANY such Y exists making Δ=0, that is a degeneracy
    EVENT (the hypothesis fails on an actual code of the family). Pre-state
    the search: enumerate {Y = (u, u+v) : u, v ∈ C'} where C' = the k-dim
    code, cap 200 seeds. 25 seeds.

  etc. — the point of layer B is *hunting*, not sampling: some of these will
  find nothing; reporting the miss rate of the hunts is the deliverable.
  All handful (B1..B5) are pre-registered WITH their caps; a hunt that hits
  its cap without finding degeneracy is a REPORTED MISS, not suppression.

## 2-fix. Cleaned adversarial operator list (authoritative; §2 (B1) wording above is superseded)

*(Written after the first commit was discovered ambiguous — ambiguity itself
is a risk; the superseded wording remains visible above per rule 5.)*

The pre-registered adversarial probes are:

| probe | cell | what varies | cap |
|---|---|---|---|
| ADV-1 | any (m,t) | seed stream, brute force for one degenerate instance | 200 seeds/cell |
| ADV-1b | (m, n=2^m, t∈{2,3}) m∈{6..11} | t small, seeds | 200/cell |
| ADV-2 | (m, n, t=((n−k)/m)) for k=2 | k=2 (one-Wronskian) instances | 25 seeds |
| ADV-3 | (6,64,3),(7,128,6) | support orderings: 25 shuffle + 25 reversed + 25 sort-by-G | 75 orderings |
| ADV-4 | m=4, t∈{1,2,3}, n as constrained | full enumeration of irreducible G at complete-support | exhaustive |
| ADV-5 | (6,64,58)+(7,128,120) | Y-manipulation: Y=(u,u+v) proportional pair | 200 seeds |

Nothing else is searched; anything else discovered opportunistically DURING
these probes is reported but does not change the rate denominators. A
degeneracy event in layer A or B is an event FOR THE RATE (layer A) or a
CERTIFIED EXHIBIT (layer B).

## 3. Budget (fixed)

- 1 nice-10 process, BLAS/OMP = 1 thread (repo resource policy).
- Sampled layer: 6 m-values × ≤2 t-cells × 25 seeds = 275 instances.
- Adversarial layer: ≤ 200 seeds/cell × ≤ 8 cells + 75 orderings + full m=4
  enumeration + 25 k=2 cells ≤ 3200 instance-builds worst case.
- Wall-clock cap for the whole campaign: **2 hours**. If the cap is hit
  before layer A finishes, report the completed cells only, with the exact
  cut stated (rule-7 scope sentence).
- Stopping: layer A fully swept (275) → then adversarial probes until done or
  cap. No data-dependent stopping anywhere in layer A.

## 4. Adjudication rule (fixed BEFORE any run)

> Let N_A = number of layer-A instances built (all events counted, including
> abort-substituted seeds — each build that survives guards is one
> denominator unit). Let D_A = number of layer-A instances with ALL
> Δ_{p,q} = 0. The reported rate is D_A / N_A with an exact Clopper–Pearson
> 95% interval reported alongside (computed by mpmath, not normal approx).

Decision thresholds (fixed now, before the number exists):

- **(i) "Generic in practice."** D_A = 0 in layer A **and** every adversarial
  probe comes back empty. Report: `0/N_A (95% CI [0, 3/N_A])`, i.e. the
  measurement bounds the failure rate above by ≈3/N_A at 95% — three-event
  rule of thumb, stated with the exact CP bound, labelled
  COMPUTATIONAL-EVIDENCE over the enumerated domain, and the hole in Apon
  §3.6 remains OPEN at the measured scale.
- **(ii) "A real hole."** Any layer-A event, or an adversarial exhibit of a
  degenerate instance that is a *legitimate output of the Gate-A
  construction* (guards all pass). Report exact rate + exhibit. **Escalate to
  Main BEFORE writing the exhibit anywhere as established**
  (hard constraint, repo rule; the exhibit touches a published theorem's
  stated hole).
- **(iii) "Boundary only."** Degeneracy arises only at degenerate-parameter
  cells with k ≤ 2 (vacuous or single-Wronskian). Reported as a boundary
  phenomenon; the rate statement keeps k ≥ 3 separately.
- **(iv) "Non-interpretable."** Guard-chain failures at rates that swamp the
  measurement (>20% build-abort in a cell) — report the abort rate itself as
  the finding; the rate over the surviving population is then NOT an
  estimate of the rate over the construction's image, and the report says
  so.

Anaointment of the label: every Δ evaluation is exact F_{2^m} polynomial
arithmetic (`gf.pmul/padd/pderiv/ptrim`) — the exact same code path as the
frozen gate-A guard (MACHINE-VERIFIED there on 13/13). No floats, no
ball arithmetic, no tolerances. The Clopper-Pearson CI is computed with
`mpmath` (binary search on the incomplete beta function) or by explicit
combinatorial summation — the interval endpoints are exact rational sums of
binomial terms, no float in the *interval* computation.

## 5. What would make this measurement wrong (pre-mortem, committed)

- **The RREF Y selection is deterministic given (m,n,t,seed)** — a degenerate
  instance found at some (m,n,t,seed) is bit-exactly reproducible from the
  frozen instance records, so an exhibit is a certificate, not an anecdote.
- If the adversarial probes systematically confound "degenerate" with
  "guard-failing" (e.g. all degenerate candidates die at the `eps` gcd
  guard), that is a FINDING about the construction (Apon's guards themselves
  exclude the degenerate branch), reported as (iv)-adjacent interpretation,
  not as a suppressible miss.
- Rule-14 self-audit, committed in advance: what could layer A NOT see?
  (a) degenerate instances that the build REJECTS (guards) — counted under
  (iv); (b) degeneracy that depends on held-count c or flag depth R —
  Δ_{p,q} is a fixed polynomial pair of the instance, c- and R-independent,
  so this family of failures is structurally invisible to the Δ-guard —
  that is exactly why the audit is out of scope for the Δ-rate and stated
  here; (c) non-binary-Goppa codes — outside the target claim; (d) t=1
  full-support — build-impossible, boundary case; (e) k<2 — vacuous, not a
  Δ-failure. (b) is the honest limit of what a Δ-rate can say about the
  waterfall: it prices one input to Lemmas 9/7, not the whole Step-3 rank
  behavior.

## 6. Commitment

This file is committed (git) BEFORE the first instance is built.�inely
Numbers recorded later go in the campaign dir (producer-generated name), not
here. If the design must change, the change is a visible ADDENDUM below,
never an edit above the original text.

— MceliececelGateB, 2026-08-30T17:35Z.

---

# ADDENDUM 1 — registered BEFORE the first computed instance (17:52Z, post-commit of the base text, pre-computation)

The base text above is committed and stands (rule 5: superseded wording stays
visible). This addendum fixes, BEFORE any instance is built, exactly three
things and changes nothing else.

## A1.1 — Algebraic reduction of the search space (proven below, then machine-verified)

At full support $n = 2^m$ (every cell registered above is full-support):

(a) $\Pi'(a_i) = 1$ for all $i$, since $\Pi = Z^q + Z$ over $\mathbb F_q$ and
$\Pi' = qZ^{q-1} + 1 = 1$. Hence $\lambda_i = G(a_i)^2$.

(b) The code equals the value-image of the admissible polynomials:
$C = \{(\lambda_i f(a_i) \bmod 2)_i : f \in E[Z]_{\le D}\}$, and
$V := \mathrm{span}_E\{f_y : y \in C\}$ (equivalently the E-span of the
interpolants of ANY $F_2$-basis of $C$) equals
$V_0 = \{f \in E[Z]_{\le D} : G(a_i)^2 f(a_i) \in F_2 \ \forall i\}$.
Derivation: $y \in C$ iff $\exists f$ of degree $\le D$ with
$\lambda_i f(a_i) = y_i$ for all $i$; invertibility of evaluation on degree
$\le D < n$ makes $f$ unique per $y$; and the parity checks transport exactly:
$\sum_i y_i a_i^j/G(a_i) = \sum_i G(a_i) f(a_i) a_i^j / \Pi'(a_i)$, which is
the degree-$(n{-}2{-}j)$ Lagrange sum, zero for $j \le t-1$ because
$\deg(Gf) = n-1$ and $t + D + j \le n-2 \iff j \le t-1$.

(c) RREF (row reduction) preserves the row space, so the $k$ selected rows
span $C$ for EVERY support ordering; interpolation is $F_2$-linear, so the
E-span of the interpolants is $V_0$ for every ordering and every basis; and
$W$ is $E$-bilinear, so the all-pairs Wronskian verdict is invariant.

**Consequence 1 (search-space collapse, to be machine-verified
in-campaign):** at full support the degeneracy verdict is a pure function of
$(m, t, G)$. Support orderings and Y-basis choices CANNOT change it. The
"support-ordering axis" of the registered search space is provably
degenerate; the ordering/basis probes (ADV-3, ADV-5) are re-scoped from
"hunting" to VERIFICATION of this invariance (they must return equal
verdicts; a violation of invariance would falsify (b)/(c) and is itself a
finding to escalate). The report states the rate over distinct G, with the
per-seed denominator also reported.

**Consequence 2 (dedup):** different seeds can draw the same G (the seed
stream rejects non-irreducible draws). The pre-registered per-seed
denominator $N_A$ is reported as registered; the semantically correct rate
denominator is the number of DISTINCT G values drawn, $N_G$, also reported,
with the dedup count. If they disagree materially the discrepancy is
reported explicitly, not resolved silently.

## A1.2 — Exact verdict algorithm (instrument upgrade, replacing the O(k^2)-gcd cascade as the front end)

The cascade (gcd of each pair, coprime parts must be squares) remains the
DEFINITIVE exact test but costs $O(k^2 \cdot D^2)$ field ops — unaffordable
at $k \approx 2000$. The engine therefore runs, in order:

1. **All-squares sufficient certificate** (exact, $O(kD)$): if every
   coordinate $f_j$ has all odd-degree coefficients zero, then every
   $f_j' \equiv 0$ (char-2 derivative kills squares), so every
   $\Delta_{p,q} \equiv 0$: **DEGENERATE**, certified directly by the
   recorded zero-odd-coefficient matrices (no cascade needed).
2. **Pointwise necessary filter** (exact, $O(kh)$ per point, $h \approx D/2$):
   $\Delta_{p,q}(a_i) = v_i[p] w_i[q] + w_i[p] v_i[q]$ with $v_i[j] =
   f_j(a_i) = Y[j,i]\lambda_i^{-1}$ (known exactly from Y, lambda) and
   $w_i[j] = f_j'(a_i)$ (one MUL-gather pass per point). If at ANY point the
   $k \times 2$ matrix $[v_i | w_i]$ has rank 2, some pair has
   $\Delta(a_i) \ne 0$ so $\Delta_{p,q} \not\equiv 0$: **NONDEGENERATE**,
   with the witness $(p, q, i)$ recorded. (Deg $W \le 2D-1$ makes pointwise
   vanishing only NECESSARY, not sufficient — that is why the certifying
   direction here is the RANK-2 one-sided witness; one nonzero point value
   certifies $\Delta \not\equiv 0$ exactly.)
3. **Full exact cascade** only for instances that pass all points' rank-1
   test without being all-squares (expected rare; affordable at small m;
   at large m an explicit budget line is recorded when it triggers).

Both certificates are exact F_{2^m} polynomial facts. No float enters any
verdict.

## A1.3 — Registered scope extensions and re-scopes (computed only after this commit)

* **E1 — m=6, t=2 EXHAUSTIVE:** all monic irreducible degree-2 polynomials
  over F_64 — exactly $(64^2-64)/2 = 2016$ of them — full census of the
  population (forced-G builds). Expected cost ~10 min. This yields an exact
  census rate at $(m,t) = (6,2)$ complementing ADV-4's m=4 census.
* **E2 — special-G families (adversarial constructions):** at registered
  cells, forced-G families with special shape:
  (i) $G = Z^3 + c$: reducible iff $c$ is a cube (gcd(3, q-1) = 3 at even
  m); enumerate ALL $c \in F_{64}$, keep the irreducibles (~42 expected).
  (ii) $G = Z^3 + Z + c$, all $c \in F_{64}$, keep irreducibles.
  (iii) $G = Z^2 + Z + c$ at m=7 (Tr(c) = 1: irreducible quadratics of the
  monic depressed family; 64 candidates).
  These are G-space probes the random stream cannot hit (measure-zero
  shapes); they are labelled ADVERSARIAL, and their rates are reported per
  family, never merged into the sampled rate.
* **ADV-3, ADV-5 re-scoped** to invariance VERIFICATION (ordering-invariance
  at (6,64,3) with 3 orderings incl. identity; basis-invariance by swapping
  two RREF rows at the same cell). Their registered caps stand; their
  purpose is stated as verification.
* **k=1 proven unreachable:** (2^m - 1)/m is not an integer for any
  2 <= m <= 11 (2^m = 1 mod m forces m | 1), so the k=1 vacuous cell never
  existed in this regime; recorded as a settled boundary, not a search.
* **k=2 exact cells recomputed:** k=2 requires $mt = 2^m - 2$; exact at
  m=7 (t=18) and m=11 (t=186); other registered k=2-intending cells
  actually give k = 4, 6, 8 (listed in the report per cell as built).

— MceliececelGateB, 2026-08-30T17:52Z, before the first computed instance.

---

# ADDENDUM 2 — instrument self-test caught and fixed an invalid iff (2026-08-30T18:2xZ, still BEFORE the first campaign instance)

The FIRST version of the O(k^2) gcd-cascade instrument (gcd the pair, strip
the common factor, require both coprime parts to be squares) was falsified
by the property-based self-test (2909/3000 disagreements): its converse is
FALSE. Exact statement of the correct char-2 algebra, now machine-checked
(selftest.py, ALL_OK, 3000/3000 agreement, 5/5 counterfactual plants):

* Write f = f_e(Z^2) + Z·f_o(Z^2) (even/odd parts). Then char-2 algebra gives
  Delta_{f,g} = f_e·g_o + f_o·g_e as a polynomial in x = Z^2 — the
  Z^{odd+odd} and Z^{even+even} terms cancel in W = f g' + f' g.
* So all-pairs coprime-parts-squares is a SUFFICIENT certificate for
  degeneracy (route all_squares-style), but the general degenerate condition
  is (f_e, f_o) proportional to (g_e, g_o) over E(x) — NOT both-squares.
* The engine's rank-scan route never used the broken iff and was correct as
  built; its completeness lemma (deg_x Delta <= D-1 < n at full support,
  distinct squares x = a_i^2) is exact and was re-verified in the same run.
* Instrument fix: gcd-cascade replaced by exact brute Wronskian polynomials
  as the small-k ground truth (three instruments: engine rank-scan,
  brute Wronskian, and known-answer/plant tests — all agreeing on 3000
  property cases and all plants).
* The synthetic-F harness of selftest v2 had its OWN artifact (a Y-synthesis
  inconsistent with the F passed in — rule 14 applied to the tester); v3
  synthesizes every case through ONE canonical lift so engine, brute force,
  and cascade act on the identical object. The v2 mismatches were an
  artifact of that harness, not of the engine; v2 was never used for any
  campaign number, and no number from any selftest enters the rate.

— MceliececelGateB.
