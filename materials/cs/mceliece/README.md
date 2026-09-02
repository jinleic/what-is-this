# `mceliece/` — the Classic McEliece hold-out "waterfall" dispute

**Status: BENCHMARK (no campaign run yet).** A four-paper dispute is live on
IACR ePrint (all four abstracts owner-read 2026-08-29, **V-owner**):

1. ePrint **2026/1747**, Kirill Vedenev, *"Extending Distinguishing to Key
   Recovery for Subfield Subcodes of GRS codes"* — turns the "hold-out
   distinguisher" of Ghoshal–Ishai–Jain–Sun (ePrint **2026/1630**, not yet
   owner-read) into a proposed key-recovery attack for GRS subcodes including
   Goppa, validated experimentally on Goppa codes over $\mathbb F_4$, and
   *conjectured* to extend to binary Goppa.
2. ePrint **2026/1810**, Daniel Apon, *"An Algebraic-Geometry Lower Bound
   against the ePrint:2026/1747 McEliece Key-Recovery Attack"* — argues the
   attack's linear equations are **not** independent across held positions: for
   a binary Goppa polynomial of degree $t$, an explicit linear map constructs a
   **$(2t+3)$-dimensional family** modulo the true solution, and at each held
   support point the entire derivative-flag block restricts on this family to
   **at most one** ordinary evaluation condition. Under a concrete nondegeneracy
   condition on the hidden vector polynomial $\mathbf F$ and $\mathbf F'$:
   $$c_{\rm need} > 2t+3,$$
   where $c_{\rm need}$ is the number of sampled held positions required at the
   critical step — and the attack cost depends on $c_{\rm need}$ **in the
   exponent**. For `mceliece8192128` ($t=128$, NIST Category 5) this yields
   $c_{\rm need}>259$ and *">2^1500 bit operations"*.
3. ePrint **2026/1786**, Markku-Juhani O. Saarinen, *"Bit Operation Cost of
   'Holdout' Key-Recovery Attacks Against Classic McEliece"* — a *living
   costing paper*: relation-generation reductions (different vanishing orders
   per column; discarding redundant binary derivative rows; injective
   zero-padding to square), a width-64 relation-generation floor of
   $2^{112.35}$ bit ops, conditional subtotal $2^{132.22}$ for
   `mceliece348864`, and — fatally for a clean cost estimate — *"End-to-end
   attack complexity for this singleton tangent/anchor route is not yet
   established."*
4. Vedenev **v2** with a rejection criterion (V-scout, not yet owner-fetched).

The empirical crux both sides appeal to: the **"waterfall"** — in experiments
on proper binary Goppa instances, just before the required independent
equation count, additional held positions *"sharply drop in value, providing
just a single, independent equation"* rather than the $\approx{k\choose 2}$
from early positions. Apon explains it via the $(2t+3)$-family; Vedenev v2
responds with a rejection criterion. **Nobody has published a controlled
exact-arithmetic measurement of the collapse on synthetic binary Goppa
instances.** That is this target's gate — one hour of exact $\mathbb F_{2^m}$
linear algebra decides it, with a certificate either way.

## Why a certificate decides a security dispute

Every quantity in the dispute is finite $\mathbb F_{2^m}$-linear algebra for a
*fixed, explicit* instance: Hasse-derivative flag spaces at held positions,
their kernels, and the rank/nullity bookkeeping across positions. An echelon
form over $\mathbb F_{2^m}$ computed by `python-flint` is an exact,
reproducible artifact — no floats, no heuristics, no decoder choice. The gate's
outcome either corroborates the $c_{\rm need}>2t+3$ bound's mechanism on real
instances or exhibits a concrete counterexample instance, which is precisely
what a v3 of this dispute would need in either direction.

## Gates

1. **Gate A (waterfall census — refute-or-confirm, ~1 h).**
   Build Hasse-derivative flag systems on **nine synthetic binary-Goppa
   instances** and measure the collapse directly. Instance ladder — every pair
   satisfies the binary-Goppa feasibility constraint $k=n-mt>0$ with
   $n\le 2^m$:

   | $m$ | $n$ | $t$ | $mt$ | $k=n-mt$ | $2t+3$ |
   |---:|---:|---:|---:|---:|---:|
   | 10 | 1024 | 24 | 240 | 784 | 51 |
   | 10 | 1024 | 32 | 320 | 704 | 67 |
   | 10 | 1024 | 40 | 400 | 624 | 83 |
   | 11 | 2048 | 32 | 352 | 1696 | 67 |
   | 11 | 2048 | 48 | 528 | 1520 | 99 |
   | 11 | 2048 | 64 | 704 | 1344 | 131 |
   | 12 | 3488 | 48 | 576 | 2912 | 99 |
   | 12 | 3488 | 64 | 768 | 2720 | 131 |
   | 12 | 3488 | 96 | 1152 | 2336 | 195 |

   The last three rows sit at the real `mceliece348864` field size ($m=12$,
   $n=3488$), whose shipped parameter is $t=64$; $t=96$ probes above it.
   **Correction to the scan input:** the scout proposed
   $\lambda\in\{9,10,11\}$ with $t\in\{72,80,88\}$ (thresholds
   $2t+3=147,163,179$). Those are infeasible — at $m=9,t=72$ we get
   $mt=648>512\ge n$, so $k<0$ and no such Goppa code exists. The ladder above
   preserves the intended threshold range while staying constructible.

   For each instance and each held-count $c\in\{2t+1,\dots,2t+4\}$, assemble
   the derivative-flag relation system in exact $\mathbb F_{2^m}$ arithmetic
   and record (i) the nullity of the system restricted to the
   $(2t+3)$-dimensional family and (ii) the number of *new* independent
   equations contributed as the $(c+1)$-st held position is added.
   **Pass (Apon):** in all nine instances each additional held position in the
   critical window contributes **exactly one** new independent equation, and
   the count stabilizes exactly at $2t+3$.
   **Refutation (Vedenev):** any instance where a held position contributes
   $\ge2$ new independent equations below the $2t+3$ threshold — exhibit the
   instance, its Goppa polynomial, its support ordering, and the echelon
   certificate.
   Runtime basis: 9 instances × 4 held-counts, each a
   $\approx(2t+3)\times(\text{flag width})$ echelon over $\mathbb F_{2^m}$ via
   `python-flint`; the largest is $195\times O(k)$ at $m=12$. Exact
   finite-field elimination at this size is seconds per system, so the census
   is well under an hour on one low-priority core. Seeds recorded in the
   campaign inventory.
2. **Gate B (nondegeneracy stress — try to break the bound, days).**
   Apon's $c_{\rm need}>2t+3$ holds *under a concrete nondegeneracy condition*
   on $(\mathbf F,\mathbf F')$. Search for binary Goppa polynomials (and held
   sets) where the condition degenerates: grid/ILP-guided search over
   square-free irreducible $g$ and support orderings, measuring whether the
   flag restriction ever yields **more** than one condition per position.
   **Refutation:** a certified instance with $c_{\rm need}\le 2t+3$ — this
   would meaningfully re-price the mooted attack and is exactly the "v3"
   input either side would want. **Pass:** the collapse is uniform across the
   searched degeneracy profiles.
3. **Gate C (symbolic audit of Apon, days).**
   Independently re-derive the $(2t+3)$-dimensional family from the stated
   linear-map construction in exact rational/symbolic arithmetic (sympy +
   flint), verify: (i) dimension count $2t+3$ modulo the true solution,
   (ii) the "at most one ordinary evaluation condition per held point"
   restriction, (iii) the equivalence to the stated nondegeneracy condition on
   $(\mathbf F,\mathbf F')$. Small-$t$ exhaustive checks ($t\le 12$) + symbolic
   derivation. **Pass:** Apon's algebra confirmed end-to-end at all small $t$;
   **Refutation:** a dimension other than $2t+3$, or a position yielding more
   than one condition without nondegeneracy failing.

## Pre-campaign requirements (before any campaign directory exists)

- Owner first-hand read of the **full texts** (abstracts are read): Vedenev
  2026/1747 (algorithm, and v2's rejection criterion), Apon 2026/1810 (the
  explicit linear map and the nondegeneracy condition statement), Saarinen
  2026/1786 (relation-model reductions — needed for gate A's system assembly
  to match the dispute's convention), and the upstream distinguisher paper
  ePrint **2026/1630** (Ghoshal–Ishai–Jain–Sun). Extract: the exact definition
  of the Hasse-derivative flag used by Vedenev, the exact sense in which Apon
  restricts "the entire derivative-flag block" to the family, and the exact
  rejected-vs-accepted criterion in v2.
- One-page pre-statement for gate A committed here before the first run: the
  precise flag construction, the held-set convention, the nullity metric, the
  nine instances with their polynomials (seeded, recorded in the campaign
  inventory), and the numeric pass/fail criterion as above.
- Cross-check gate A's system assembly against Saarinen's reductions so a
  "refutation" cannot be dismissed as a convention mismatch.

## Disjointness

No overlap with `../math/` (classical error-correcting *code design* —
"Automated co-design of quantum LDPC codes" — is quantum-side) or `../physics/`
(decoders for quantum codes; this is the classical McEliece public-key
cryptosystem, complexity and finite-field algebra only). No attack is executed
on real keys; all instances are synthetic and published in the campaign
inventory.

## Layout

- `src/` — implementation (not yet created).
- `campaigns/` — immutable run snapshots (not yet created).
- `scratch/` — non-authoritative exploration.

## Running

```sh
cd cs
./.venv/bin/python -c "import flint; print(flint.__version__)"   # toolchain smoke
```

## Current state (2026-08-30, agent `Mceliece`)

**Gate A: PASS — Apon's waterfall mechanism confirmed by measurement.**
Evidence label: MACHINE-VERIFIED (exact $\mathbb F_{2^m}$ arithmetic; no
floating point, no ball arithmetic anywhere in the census).

**Claim of record — the enumerated instance list.** No universal over $m$
or $t$ is claimed; "$m\le N$ complete" is not asserted for any $N$.

| $m$ | $n$ | $t$ | $k$ | $2t{+}3$ | seed | verdict |
|---:|---:|---:|---:|---:|---:|:--|
| 6 | 64 | 3 | 46 | 9 | 1387 | PASS |
| 6 | 64 | 4 | 40 | 11 | 2311 | PASS |
| 6 | 64 | 5 | 34 | 13 | 3413 | PASS |
| 7 | 128 | 3 | 107 | 9 | 4421 | PASS |
| 7 | 128 | 6 | 86 | 15 | 5531 | PASS |
| 8 | 256 | 3 | 232 | 9 | 6637 | PASS |
| 8 | 256 | 5 | 216 | 13 | 7741 | PASS |
| 9 | 512 | 12 | 404 | 27 | 9109 | PASS |
| 9 | 512 | 24 | 296 | 51 | 9137 | PASS |
| 9 | 512 | 40 | 152 | 83 | 9277 | PASS |
| 10 | 1024 | 24 | 784 | 51 | 8001 | PASS |
| 10 | 1024 | 40 | 624 | 83 | 5113 | PASS |
| 11 | 2048 | 48 | 1520 | 99 | 6211 | PASS |

On every instance and every held count $c$ in the swept range:
$N_{\rm fam}(c) = 2t+3-c$ **exactly**, reaching $0$ precisely at $c = 2t+3$;
every per-(point, order) family-restricted block has rank **exactly 1**
(rank 1 at even $j$, 0 at odd $j$ — the $(j{+}1)$ parity in
$(j{+}1)\partial^{[j+1]}\mathbf F$); **zero refutation events**; the
$\Delta_{p,q}\neq0$ nondegeneracy guard passes; the lower-bound certificate
($T(A_S)$ jets inside the flags at all $c = 2t+2$ held points) passes.

**Not run** (the boundary is part of the result): $m=12$ — ALL $t$; and at
each swept $m$, every $t$ outside the table above. Support is always the
full field $n = 2^m$; no other support size was tested at any $m$. No
instance with $t\ge49$ was run at any $m$, so `mceliece8192128` ($t=128$)
is untouched and $c_{\rm need}>2t+3$ at cryptographic scale is **not**
established here. Depth-uniformity of Lemma 7 beyond $R=4$ remains
CITED-DEPENDENCY (Apon Lemma 7).

**Guard chain** (uniform across all 13 rows, one code path):
$\delta$ MACHINE-VERIFIED exact via the Lagrange route with an
abort-on-failure degree condition $\max_j\deg f_j\le D$;
$\alpha$ (the differential identity) **MEASURED per instance** — it is the
independent check on Apon Lemma 3 itself and is deliberately *not* replaced
by a citation to that lemma; $\beta,\gamma,\varepsilon$ measured per
instance.

**Campaigns.** `2026-08-30T01-42Z_9D0FC9E7` (small-$t$ census),
`2026-08-30T04-47Z_98C586C7` (ladder addendum), `2026-08-30T12-03Z_DCAB8F20`
($m=9$, closes the ladder hole), `2026-08-30T13-57Z_41B2CA1B`
(consolidated + boundary tables), `2026-08-30T14-09Z_57200ADD`
(uniform $\delta$ across all 13, one code path).

**Source tracking.** Vedenev 2026/1747 was revised in place on 2026-08-30,
detected by content hash (245,667 B → 264,802 B; md5
`cd044d4c…` → `f27d10be…`; 1633 → 1835 text lines). Its new §9 changes the
*cost* framing of Frobenius alignment; §8.1, $c^\star = 2\delta+3$, Apon's
Theorem 1, criterion (24) and Conjecture 5 are unchanged, so the object this
census measures is unchanged and none of v3 is folded into the verdict. The
superseded copy is retained beside the revision. The hold-out dispute
remains four-paper (1630, 1747, 1786, 1810) as of 2026-08-30; ePrint
2026/1778 (HOVER, Saarinen) is a fifth McEliece paper in the same wave but
on the higher-order-vanishing route, not the hold-out route.

## Current state (2026-08-30, agent `MceliececelGateB`) — Gate B: the `Δ_{p,q}` genericity rate

**Measurable object.** Apon's Lemma 9 (ePrint 2026/1810) needs some pair
$p<q$ with $\Delta_{p,q}=f_pf_q'+f_qf_p'\neq0$ (char-2 Wronskian) to make
$A\mapsto T(A)\bmod EF$ injective; his §3.6 states verbatim that if every
$\Delta_{p,q}$ is zero the proof does not apply. Gate A had observed the
guard on 13/13 census instances — in-sample, hence no evidence about a
RATE (the same in-sample/out-of-sample defect that governs `rs-pe3d`'s
pattern P). This campaign measured the out-of-sample failure rate,
adversarially and by sampling, with the adjudication rule frozen in
advance (pre_registration committed before the first computation, git
`02166a6` + addenda `03ecafd`/`b8e8273`).

**Result — three separate findings, plus a negative.** No degeneracy
event was found anywhere in the searched space, and NO configuration
forcing all-$\Delta$-zero was found: the §3.6 hole remains OPEN — this
campaign did NOT close it, it measured how hard it is to hit (Main
confirmed the framing mid-run).

1. **EXHAUSTIVE THEOREM at $(m,t)=(6,2)$** — MACHINE-VERIFIED. All 2016
   monic irreducible degree-2 polynomials over $\mathbb F_{64}$
   (population exhausted; count $(q^2-q)/2$ exact) built at complete
   support $n=64$: **0/2016 degenerate**. For this cell, EVERY binary
   Goppa polynomial satisfies Apon's genericity condition. An exhaustive
   census carries no confidence interval — attaching one would be a
   category error.
2. **SAMPLED layer: 0/250** — COMPUTATIONAL-EVIDENCE, and a rate over the
   sampling distribution only, NOT over the space (that qualification is
   part of the claim). Ten $(m,t)$ cells × 25 registered seeds
   (m ∈ {6..11}, smallest admissible $t$ per m; $t=1$ at full support is
   construct-impossible and $k=1$ unreachable for $m\le11$): 0
   degenerate, 0 guard failures. Certified one-sided 95% Clopper–Pearson
   upper bound **0.0119** (Arb-ball, radius $1.2\times10^{-16}$;
   tail$(\text{upper})=0.05$ verified; mpmath cross-check agrees to
   $10^{-16}$).
3. **ADVERSARIAL special-G families: 0/127** — COMPUTATIONAL-EVIDENCE
   over these named families: $Z^3{+}c$ (42/64 irreducible),
   $Z^3{+}Z{+}c$ (21/64), all depressed quadratics $Z^2{+}Z{+}c$ at
   $m=7$ (Tr$(c)=1$; 64/128). CP95 upper bound **0.0233**.
4. **Ordering/basis invariance (new exact algebra, proven AND
   machine-checked).** At full support $\Pi'\equiv1$ and
   $\lambda_i=G(a_i)^2$; the code is the value-image of $V_0=\{f:
   G(a_i)^2f(a_i)\in\mathbb F_2\}$ independent of column order, RREF
   row-selection, and any $E$-basis change; $W$ is $E$-bilinear; hence
   the all-pairs verdict is a pure function of $(m,t,G)$. The
   support-ordering axis of the originally-registered space is therefore
   provably degenerate and was re-scoped (pre-registered ADDENDUM 1) to
   verification probes, which passed. Degenerate cells that do not
   exist: $t=1$ (root in support), $k=1$ ($2^m\not\equiv1\bmod m$).

**Instrument (rule 14 discharged with counterfactual plants).** Verdicts
by a NEW exact engine: pointwise rank-scan of $[f(a_i)\,|\,f'(a_i)]$ with
the completeness lemma $\deg_x\widetilde\Delta\le D-1<n$ at full support
($x=Z^2$, even/odd split $\Delta=f_eg_o+f_og_e$, Frobenius injective).
Three exact instruments agree on 3000 property cases; FIVE counterfactual
plants (all-squares, duplicate rows, shared-non-square, zero-row,
witness-constructive) prove the engine can return DEGENERATE, not merely
fail to return it. The self-test caught a real iff-bug in my first gcd
cascade (coprime-squares only SUFFICIENT for degeneracy; the coprime
degenerate condition is $(f_e,f_o)\parallel(g_e,g_o)$ over
$E(x)$) BEFORE any campaign number existed. The frozen Gate-A ladder was
rebuilt byte-exactly from seeds before use (13/13); the $\Psi=C\cdot J$
factorization is NOT on the Gate-B verdict path (waterfall census only);
the interpolation identity $\lambda_iF(a_i)=Y_{j,i}$ was re-derived and
re-verified. Inherited-quantity disclosure given to Main mid-run and
confirmed.

**Scope of what was swept (rule 7).** Gate B, exact-$\mathbb F_{2^m}$
$\Delta_{p,q}$ verdicts on the Gate-A construction at complete support:
the full population of irreducible $G$ at $(m,t,n)=(6,2,64)$ — 2016
polynomials, zero parameters unswept inside that cell; sampled seeds
1000–1024 at the ten cells $(m,t)\in\{(6,2),(6,3),(7,2),(7,3),(8,2),
(8,3),(9,2),(9,3),(10,3),(11,3)\}$ with support $n=2^m$ and seed-stream
regeneration rule as registered; the three structured-G families named in
(3); five invariance probes; NOT searched: any $m>11$, any $t$ beyond the
cells above at any $m$ (in particular $t\ge12$ at $m\le11$ and
cryptographic scale $t=69/128$), any support size $n<2^m$, non-full-rank
or non-Apon-style instance builds, held-position or flag-dependence of
degeneracy ($\Delta_{p,q}$ is $c$- and $R$-independent by construction),
and non-binary-Goppa codes. Differing constructions excluded by design:
$\Psi=C\cdot J$-based waterfall census quantities were not merged into
this rate, and no quantity measured here is asserted about any
construction other than the Gate-A instance build.

**Not decided by this campaign:** whether any degenerate binary-Goppa
instance exists at ANY size (open; bounds now much tighter at the swept
cells); whether Apon's $\Delta$-condition can fail on a sub-support
instance; the $m=12$ and cryptographic-scale waterfall (explicit
non-goal, inherited from the predecessor's stop decision).

**Campaign.** `campaigns/2026-08-30T21-01Z_B6EB2171/` (frozen,
content-hashed; `manifest.md` with the incident log incl. the reverted
frozen-file touch and the 9-instance stale-module rebuild; evidence
labels per row inside).

**Amendment to the section above (2026-08-30, post-acceptance, Main-directed
sharpening — appended beside, not replacing, the committed text):**

1. *k=1 cell:* the congruence $2^m\not\equiv1\pmod m$ has **no solution for
   any $m>1$ whatsoever** (for any prime $p\mid m$, $\mathrm{ord}_p(2)$
   divides both $m$ and $p-1$, contradiction; direct check at $m=2$). The
   $k=1$ boundary is therefore unreachable **unconditionally**, not merely
   in the swept range.
2. *Sampled-layer semantics (now index-canonical):* since the Gate-B
   invariance theorem makes the verdict a pure function of $(m,t,G)$ at
   full support, the registered seeds vary $G$ and **nothing else**. The
   0/250 layer is therefore a sample over **irreducible polynomials $G$**
   within each $(m,t)$ cell, and the CP bound 0.0119 bounds the degeneracy
   rate **over that within-$G$-distribution** — a consequence of the
   campaign's own invariance theorem, not an assumption. The $t=1$
   full-support impossibility is likewise unconditional ($\beta\in E$ =
   support); only the sub-support case is excluded-from-scope rather than
   impossible.

**Amendment 2 (2026-08-30, closes Main's last review item — appended
beside the committed text; no number changes):**

*K=1 nonexistence, airtight form (owner-supplied tightening of my parity
sketch; adopted verbatim as the canonical proof):* let $p$ be the
**smallest** prime dividing $m>1$. If $p\mid 2^m-1$ then
$\mathrm{ord}_p(2)$ divides both $m$ and $p-1$ (Fermat). But
$\mathrm{ord}_p(2)\le p-1<p$, while every prime factor of $m$ is
$\ge p$ — so $\mathrm{ord}_p(2)$ shares no prime factor with $m$,
forcing $\mathrm{ord}_p(2)=1$, whence $p\mid 2-1=1$: contradiction.
Main numerically verified no solutions to $m\mid 2^m-1$ at every
$m\le 20{,}000$. Hence the $k=1$ cell is unreachable for **every**
$m>1$, unconditionally.

*Sampled-layer semantics (owner-supplied sharpening of my own invariance
theorem — a consequence, not a new measurement):* the verdict is a pure
function of $(m,t,G)$ at full support; the registered seeds therefore vary
$G$ and nothing else. The sampled layer is a sample over irreducible
Goppa polynomials $G$ within each $(m,t)$ cell, and the CP bound 0.0119
bounds the degeneracy rate over that **within-cell $G$-distribution**.

— `MceliececelGateB`, final item closed; standing down.

## Current state (2026-08-31, agent `MceliecelGeneric`) — source-freshness correction and route cutover

**Rule-17a correction (MACHINE-VERIFIED).** The preceding source-tracking
sentence naming `pe1747_v3_2026-08-30.pdf` as the superseded copy is
**FAILED as a path-role statement**: on disk that path contained the current
264,802-byte PDF while canonical `pe1747.pdf` contained the previous
245,667-byte PDF. No content was deleted. Canonical `prim/pe1747.pdf` now
contains current IACR version `20260830:010812`, SHA-256
`f3e83fe71596367d6544048fa5ac0c2c4400309a9b7227a389cf4c565d26321a`;
the previous IACR version is preserved byte-exactly as
`prim/pe1747_superseded_20260828T221256.pdf`, SHA-256
`e9fc4fff89a27a579840b2c2a437714bc9596b3162fbaf8ad66e96b3de544525`.

**ePrint 2026/1786 changed materially, not cosmetically.** The locally
cached 556,874-byte PDF is exact IACR version `20260828:000118` and is now
preserved as `prim/pe1786_superseded_20260828T000118.pdf`, SHA-256
`ee0c447b251d5c18735f4687ac8b476e9aeb5dcd57c5660a4fc3fb57e4da9d21`.
Canonical `prim/pe1786.pdf` is current version `20260830:172206`,
340,507 bytes, SHA-256
`12e244c760a068d74bc2784c79cfe0b0a6eef698d7bf10994fadb5c7bc50731d`.
Current Appendix B says Sections 2–5 replace the singleton tangent/anchor
route with a one-chart direct locator-recovery route, remove superseded
singleton tables and artifact functions, and revise Section 7 and Assumption
1. The full old/current text diff and all hashes/roles are frozen in
`prim/source_diff_addendum_2026-08-31.md` and
`prim/provenance_manifest_2026-08-31.json`.

**Semantic map (HUMAN-AUDITED against current 1630/1786 and frozen Gate-B
code): disjoint.** Gate B decides Apon's full-support polynomial predicate
`all Delta_{p,q}=0` for interpolation polynomials `F`. Current 1786
Assumption 1 instead quantifies a shortened multivariate homogeneous-jet
matrix `E_l`, its kernel factorization `C*=U*V*`, accepted projective labels,
and full/block-projected kernel ranks. Gate B constructs none of those
objects; current 1786 contains no `Delta`/Wronskian premise. Therefore the
Gate-B census tests **none of current Assumption 1(1)–(4)**. Its frozen
theorem and sampled/adversarial bounds remain valid in their stated Apon
domains, but are not evidence for the live direct route. The planned
genericity extension was stopped before computation.

**No hole was closed.** Apon's Section 3.6 all-`Delta`-zero branch remains
OPEN; the current 1786 revision changes routes rather than proving it
impossible. The current direct route is separately conditional on its four
finite-instance premises, including four-holdout singleton rigidity and
full/projected rank admission. No refutation of either paper is claimed.

**Next campaign.** Pre-register and exactly test the smallest nontrivial
binary-Goppa direct-route cell
`(m,n,t,k,D,d,s,h)=(3,6,1,3,3,5,3,4)`, where the binary-Goppa no-slack
identity is `15=dD+1-t`, `s-1=2` is a power of two, the public jet matrix is
only `21 x 60`, and every locator block has size 9. This replaces the old
`Delta`-extension campaign; no matrix was constructed before the new
pre-statement.

## Current state (2026-08-31, Main) — direct-route toy census exhausted; registered genericity fails

The scalar parent campaign
`campaigns/2026-08-31T08-44Z_DA168C79/` passed all mandatory controls and
durably completed the contiguous `beta=8,...,21` prefix. Its registered
7,200-second wall cap then interrupted `beta=22` inside the public-phase exact
RREF. The 14 complete cells all pass P1_U3 and P1_U4 and all fail P2, P3, and
P4. This prefix is preserved as **COMPUTATIONAL-EVIDENCE / NOT EXHAUSTIVE** in
`partial_result.json`; `checksums_partial.sha256` verifies all 22 result,
control, log, and plant artifacts.

The separately preregistered exact continuation
`campaigns/2026-08-31T16-12-46Z_vectorized-rref-replay/` replaced only scalar
finite-field elimination with table-indexed `uint16` Gauss–Jordan operations
and cached each cell's unchanged `E` and `ker(E)`. All 26 controls passed.
In particular, the complete RREF tuple, nullspace basis order, and rank agree
on 1,004 deterministic systems; a wrong multiplication table is detected;
and cached and uncached `beta=8` payloads agree with the hard-hashed parent.
An additive strict audit also compares the full completed parent prefix
`beta=8,...,21`, removing only `header` and `t_utc`: all 14 semantic payloads,
including `complete=true`, are exactly equal.

The continuation exhausts all 24 registered monic degree-one polynomials over
the fixed GF(32) support. P1_U3 and P1_U4 pass in all 24 cells and all 120
branches of each size. P2 fails in all 24: every cell accepts 33 projective
labels rather than the expected five. P3 fails in all 24 with 28 unmapped
labels per cell. P4 fails in all 24 with
`rank(E)=36`, `nullity=244`, and both the full and block predicates false.
Therefore the conjunction of P1–P4, and hence registered `genericity_pass`,
is **false in every cell of this toy population**. Runner exit code 4 is the
predefined mathematical-failure exit, not an instrument or integrity failure.
Exact result: `report.json`; state manifest SHA-256
`ac21cdf3d11cb2844138e0bf459f6b53df3b457b78b81ab27285a40173ee0a79`;
all static, state, and final checksum ledgers pass.

**Protocol disclosure.** Frozen control C26 removed `complete` together with
`header` and `t_utc`, while the pre-statement permitted removal of only the
last two fields. The frozen control is not relabeled. Final assembly already
requires `complete=true` in every successor cell, and the additive post-run
strict audit checks `complete` across all 14 parent-prefix cells and passes.
This repairs the omitted comparison after the run; it is not represented as a
retroactive pre-production control.

**Rule-7 boundary.** The exhausted domain is only the 24 polynomials
`G_beta=Z+beta`, `beta=8,...,31`, on ordered support `(0,...,7)` at
`(m,n,t,k,ell,n_ell,k_ell,D_ell,d,s,h)=(5,8,1,3,0,8,3,5,7,5,4)`.
It is outside current ePrint 1786 Table 1 and Assumption 1's five Classic
McEliece cells. The result validates and falsifies predicates only on this toy
population; it does **not** refute the paper's conditional assumption, verify
any NIST cell, or close Apon's unrelated all-`Delta`-zero hole.

**Evidence grade.** The successor's exact finite-field census, predicate
verdicts, controls, and equality audits are **MACHINE-VERIFIED**. Parent and
successor wall times are **COMPUTATIONAL-EVIDENCE** only.

## Current state (2026-09-01, agent `MceliecelM12`) — m=12 counter-swept at the NIST cell; alpha MEASURED, not derived

**Verdict.** The pre-registered bounded m=12 campaign
`campaigns/2026-09-01T03-28-10Z_M12ALPHA_DIVONLY/` (pre-statement committed
`7fb44ae` BEFORE any m=12 compute) ran and PASSED at the single registered
instance $(m,n,t,\text{seed}) = (12, 3488, 64, 16384)$ ($k=n-mt=2720$,
$D=n-2t-1=3359$, $2t+3=131$). Exact result artifacts: `verdict.json`
(fields `alpha_identity_measured_full_grid=true`,
`alpha_grid_points_checked=3488`, `alpha_failures=[]`,
`delta_chain.delta_exact=true`), incremental checkpoints
`m12_build_checkpoint.json` / `m12_alpha_result.json`, all inside
`campaigns/2026-09-01T03-28-10Z_M12ALPHA_DIVONLY/` with `checksums.sha256`
(17/17 verified).

**ADDENDUM 2 relationship — implemented as written, fork resolved by
measurement, no amendment to the design itself.** Clause (iii)'s first
path was AFFORDABLE at the NIST cell: in-process calibration
(`calib.json`, CPU by `time.process_time()` only, per Main's steering
that `ps`-based rates are broken on this box) measured the unmodified
alpha instrument at 0.1504 CPU s per support point at m=12 scale,
projecting 524.6 CPU s for the full 3488-point grid — inside the
pre-registered 10,800 CPU-s budget with 20x headroom, so the DERIVE
fork (alpha as CITED-DEPENDENCY) was NEVER reached and NEVER used. The
fork resolution (`AMENDMENT_1_fork_resolution.md`, committed `30b6c8f`)
records the measured numbers BEFORE the verdict compute. Two further
amendments, both committed before the corresponding compute: AMENDMENT_2
corrected my own CF-2 plant expectation pre-run (F+Pi is INVISIBLE to
the value-level alpha check at support points — both cross terms cancel
in char-2 and Pi vanishes on the support — so its rejection is the
DEGREE-GATE ABORT, exactly ADDENDUM 2 clause (i)'s load-bearing
warning), and AMENDMENT_3 proved the originally registered CF-3 plant
(row-XOR) mathematically inert a priori (both sides of the identity are
F_2-additive in each coordinate, so no row-XOR can fail alpha) and
substituted a scalar-multiply plant (row 0 -> 3*f_0; 3^2=5!=3 verified
in-run). All three amendments are quoted-and-reasoned files in the
campaign dir, dated, checksummed, committed before their targets ran.

**alpha at m=12 is MEASURED (MACHINE-VERIFIED arithmetic) — the headline
defect does not occur here.** The m<=11 alpha instrument is
`cs/mceliece/src/instance.py::_guards`, step (alpha): for each support
point, the k-vector LHS $\Pi(a)F'(a)+\Pi'(a)F(a)$ and RHS
$G(a)^2F(a)^2$ are evaluated by exact uint16 log-table gathers
(`fastfield.EField.MUL`), vectorized over all k coordinates; char-2
derivative = shift-left with odd-coefficient parity. This campaign ran
that IDENTICAL operation sequence over ALL n=3488 support points —
byte-identical code path to m=6..11: the code snapshot in the campaign
dir's `code/` is byte-diffed against the frozen 14-09Z snapshot
(`campaigns/2026-08-30T14-09Z_57200ADD/code/`) and differs only by
`instance.py`'s inert `forced_G` ctor parameter (line-range 84-91), an
unused-input artifact of the Gate-B tooling; `fastfield.py`, `gfield.py`,
`census.py`, `verify_delta_uniform.py` are byte-identical. Result:
identity holds at every support point, 0 failures
(`verdict.json:alpha_failures=[]`), 700.53 CPU s measured
(`verdict.json:alpha_grid_cpu_s`). The label is MEASURED, not
CITED-DEPENDENCY and not derived-imported-relabelled; the one-sentence
addendum sentence is therefore NOT needed for m=12, and no m<=11 row's
label moved.

**delta at m=12: MACHINE-VERIFIED exact.** The full ADDENDUM-2(i)
5-part chain (exact division Pi=(Z-a_i)q_i zero-remainder on all 3488,
n unit checks L_i(a_i)=1 3488/3488 pass, 200 off-diagonal probes 0
violations, 3/3 assembly spot-checks, degree gate max_j deg f_j = 3359
= 3359 exactly) — `verdict.json:delta_chain.*`, each part with its own
measured CPU cost. Hypothesis hygiene also machine-checked: gcd(G,G')
constant (square-free), lam recomputation byte-equal, F_2-linearity
control on 4 seeded selectors.

**Anchors (hard gate, pre-statement section 2): 3/3 byte-equal.** Rows
(11,2048,48,6211), (10,1024,40,5113), (6,64,3,1387) re-run through THIS
campaign's code path and compared field-by-field against the frozen
14-09Z records
(`campaigns/2026-08-30T14-09Z_57200ADD/delta_uniform_13instances.json`):
every field byte-equal except `elapsed_s` (declared timing exclusion);
one additive field `guards_build` disclosed and not part of the frozen
set. Artifacts: `anchor_11_2048_48_6211.json`,
`anchor_10_1024_40_5113.json`, `anchor_6_64_3_1387.json` (field
`byte_equal_vs_frozen=true` in each; independently re-derived 26/26 in
the owner-side cross-check).

**Counterfactual plants (rule 14): 3/3 REJECTED**
(`plants_result.json`). CF-1 duplicate-support: Pi'(a)=0 at positions
[10,11] after the registered duplication a_12:=a_11 — the build's
distinct-support guard fires. CF-2 f_0 -> f_0 + Pi: deg 2048 > 1951 = D
— degree-gate ABORT, NO delta verdict produced, and the value-level
invisibility prediction CONFIRMED (0 failures / 22 probe points), a
live demonstration that the ADDENDUM-2(i) degree condition is
load-bearing. CF-3 row_0 -> 3*f_0: alpha value check FAILS on row 0
(first failing points 1, 2, 5), c^2=5!=3 in-run.

**Budget (registered pre-compute, measured in-process).** Build 2175.43
CPU s (`verdict.json:build_cpu_s`), alpha grid 700.53
(`verdict.json:alpha_grid_cpu_s`), delta chain 15.0
(`verdict.json:delta_chain.*` CPU fields), total 2893.37
(`verdict.json:total_cpu_s`) of the registered 10,800 CPU-s budget
(`verdict.json:within_registered_budget=true`); wall ~51 min. One
disclosed calibration miss: the AMENDMENT-1 projection under-estimated
the alpha grid (524.6 vs actual 700.53, factor 1.33); absorbed by
budget headroom, recorded for future calibration honesty.

**Rule-7 scope.** Swept: the single instance (12,3488,64) at seed 16384
with alpha over the FULL support (3488 points, identical instrument as
the certified ladder) and the exact 5-part delta chain; the 3 anchor
rows in byte-equality mode; 3 counterfactual plants. NOT swept: ANY
other t at m=12 (t=48 and t=96 remain UNREACHED), any n != 3488, any
support ordering other than the single seed-16384 shuffle, any m > 12,
any census row (N_fam / per-point rank / held-count c-ladder) at m=12,
any sub-support construct, any non-irreducible G, no G population
exhausted (one seeded irreducible, scalar-Rabin-confirmed). This is a
BOUNDED FINITE statement about ONE instance: NOT "m=12 complete", NOT a
verification of any NIST attack-cost figure, and it does NOT close
Apon's section-3.6 all-Delta-zero hole (still BOUNDED by the Gate-B
rates, not closed).

**Named next action.** Second m=12 row at t=96 (the README ladder's
probe-above-ship row), fresh pre-statement, budgeted from this
campaign's measured rates (build ~2175 s, alpha grid ~700 s, delta
chain ~15 s => ~1 h CPU per row at t=96 scale).

— `MceliecelM12`, 2026-09-01T06:09Z; campaign frozen at commit `240820d` + ledger fix `6027427`.
