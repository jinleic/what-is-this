# Theoretical CS results — authoritative index

**Current through 2026-09-01.** Repository opened 2026-08-29; a second wave of gate-C
campaigns and a repo-wide source audit landed 2026-08-30, followed by the dated
frontier and mutable-source audit on 2026-08-31, and on 2026-09-01 the
`rs-pe3d/` H-MINLINE falsification with the certified T-DGE replacement and the
`omega/` branch-pinning closure. Evidence labels
are literal and defined in [`README.md`](README.md#evidence-labels):
**MACHINE-VERIFIED** is a named executable fact; **HUMAN-AUDITED** is ordinary
checked mathematics; **COMPUTATIONAL-EVIDENCE** is non-proof numerical
evidence; **CITED-DEPENDENCY** is imported; **OPEN** is unresolved; **FAILED**
means the wording or inference is not established.

## Headline result

**A certified upper bound on the Grothendieck constant that is stronger than
the one its source paper states.**

$$K_G \;\le\; 1.7818413238804372880\ldots$$

improving arXiv:2608.11158's stated $1.7818666069360661$ by
$2.528\times10^{-5}$, i.e. an improvement over Krivine's constant of
$3.72654\times10^{-4}$ where the paper claims $\ge3.4737\times10^{-4}$.

**Position against the full known field, 2026-08-30 (owner, Arb at 400 bits).** The
competitive set was incomplete until today: `2608.14817`'s bibliography exposed four 2026
preprints on $K_G$ that this repository did not know (recorded in
[the scan's post-scan audit](docs/CS_FRONTIER_SCAN_2026-08-29.md), finding 9). Every
published upper bound now known, weakest first:

| bound | value | source |
|---|---|---|
| Krivine 1977 | $1.782213978191369112$ | `Kri77` |
| $\text{Kri}-10^{-500}$ | $1.782213978191369112$ | BMMN13 (numerically invisible — qualitative) |
| $\text{Kri}-10^{-217}$ | $1.782213978191369112$ | Heilman `2606.00247` (degree-3 Hermite, analytic) |
| $\text{Kri}-1.013\times10^{-5}$ | $1.782203848191369112$ | Heilman `2606.00247` **Thm 1.9**, exact-$\eta$ interval certificate ($\eta=0.04249900400783211$, $d=3$, Sage, 54 h Rouché) |
| $\text{Kri}-10^{-5}$ | $1.782203978191369112$ | Heilman's own abstract (rounded); Li et al. `2606.03991v3` |
| $\text{Kri}-10^{-4}$ | $1.782113978191369112$ | `2608.11158` v2 abstract |
| certificate | $1.781866606936066100$ | `2608.11158` shipped `d3h_certificate.json` |
| **ours** | $\mathbf{1.781841323880437288}$ | this repository |

So the headline is **the strongest upper bound on $K_G$ presently known**, $3.62524311\times
10^{-4}$ below the best previously published *rigorous computer-assisted* bound — Heilman's
Theorem 1.9, whose exact-$\eta$ interval certificate gives $\text{Kri}-1.013\times10^{-5}$,
$1.3\times10^{-7}$ **stronger than his own abstract's rounded $10^{-5}$**. Citing his theorem
rather than his abstract is the same discipline this index applies to its own source: quote the
strongest number a paper actually establishes, in both directions. The $10^{-217}$ and $10^{-500}$ figures are
qualitative existence statements, not competitive numbers, and are listed only so the field is
complete. On the lower side $6\pi/11=1.71359599286715995$ still dominates both 2026
improvements over $K_{\rm DR}$ — Heilman `2603.22616` at $K_{\rm DR}+10^{-26}$ and
Jones–Malavolta `2603.30039` at $K_{\rm DR}+10^{-12}$ — by $0.0366393187$.

**Standing caveat, honestly stated.** Heilman `2606.00247` reaches its bound by *interval
arithmetic on Hermite thresholding at degrees three and five*, which is the same instrument and
the same degrees as `kg/` gate C. That gate is **held** pending a written duplication finding
from its agent. This does not affect the headline above, which rests on slack in
`2608.11158`'s own certificate, but it does mean the *method* is not ours alone and the repo
never knew it until today.

**Source-version audit, 2026-08-30 — and its gap found a live revision.** The owner's sweep covered the five **arXiv** sources by API but explicitly excluded `mceliece/`'s IACR **ePrint** primaries, since ePrint exposes no equivalent dated field. A real in-place revision then turned up exactly there: **ePrint 2026/1747 (Vedenev) was revised on 2026-08-30, after our copy**, detected by **byte comparison rather than a listing date** — 245,667 bytes md5 `cd044d4c…` became 264,802 bytes md5 `f27d10be…`, text growing 1633 to 1835 lines. The addition is a new Section 9 on extended multiplicity shortening and **Frobenius alignment**, which he states removes the $m^{c-1}$ alignment factor from Step 3's cost, leaving only the $q^{m(c-2)}$ support enumeration. **It does not touch what `mceliece/` measures:** Section 8.1 is intact, $c^*=2\delta+3$ for binary Goppa is still asserted (line 1404), Apon's Theorem 1 is still the quoted lower bound, and the v2 rejection criterion (24) and Conjecture 5 are unchanged — so the 13-instance ladder stands as measured and nothing was folded in. Separately, **ePrint 2026/1778 (Saarinen, HOVER)** is a fifth McEliece cryptanalysis paper in the same wave but on the **higher-order-vanishing** route, not the hold-out route, so the hold-out dispute proper remains four-paper (1630, 1747, 1786, 1810); it reports end-to-end key recovery on five TII challenges including TII-252 while stating it does not threaten Classic McEliece parameters in its present form. **Method now binding as README rule 17a: sources without a dated API are tracked by CONTENT HASH, with the superseded copy kept alongside rather than overwritten.**

All five arXiv source papers re-checked
for post-scan revisions: `2608.16884` (omega), `2607.28676` (mm3), `2608.16649` (oct-rank)
and `2305.07156` (delcap) are all still at **v1**. `2608.11158` (this headline's source) is
at **v2, updated 2026-08-12** — before the owner read of 2026-08-29, so what we read *was*
v2 and nothing here is stale. No source has been superseded since the scan.

**CORRECTION — 2026-08-31, owner:** the final sentence above is valid for the five
arXiv sources it names, but false as a statement about all seven targets. A live
byte/version-history audit found that ePrint `2026/1786` changed **four times** after
the cached `20260828:000118` PDF. Current `20260830:172206` bytes (340,507 bytes,
sha256 `12e244c…`) replace the singleton tangent/anchor route with a one-chart direct
locator-recovery route and report complete **conditional** cost tallies for all five
Classic McEliece parameter sets. The released `mccost.py` totals replay exactly, and
its solver-failure extension degrees were independently certified by exact rational
comparison plus 300-bit Arb; that verifies arithmetic, **not** the four finite-instance
premises in current Assumption 1. Existing `mceliece/` gate B tests a different
Wronskian premise and cannot be cited for the new route. The full provenance matrix,
hashes, conditional scope, and exact threshold table are in the newest
[`PROGRESS.md`](PROGRESS.md) entry. Remote ePrint `1810`, `1778`, and `1630`, ECCC
`TR26-150`, and all 14 queried arXiv primary/adjacent records were unchanged; `1747`
also exposed a local filename-role inversion, now being repaired without deleting
either byte version.

**Citation hazard, locked deliberately.** The v2 *abstract* states the upper bound only in
the conservative closed form $\pi/\bigl(2\log(1+\sqrt2)\bigr)-10^{-4}$, which evaluates to
$1.7821139781913691118$ — **weaker than the paper's own certificate value**
$1.7818666069360661$. The full ordering, Arb-verified at 300 bits:

$$\underbrace{1.78184132}_{\text{ours}}\;<\;\underbrace{1.78186661}_{\text{paper certificate}}\;<\;\underbrace{1.78211398}_{\text{v2 abstract}}\;<\;\underbrace{1.78221398}_{\text{Krivine}}$$

So our margin over "the paper" is $2.5283056\times10^{-5}$ against its certificate but
$2.7265431\times10^{-4}$ against its abstract — a **$10.8\times$ inflation available for
free** by citing a different line of the same paper. This index cites the certificate, the
strongest number the paper actually establishes. That choice is deliberate and must not be
"improved": comparing against a source's most conservative published statement instead of
its best one is how a real result gets turned into an overclaim.

The mechanism is that the paper states a conservative $\gamma$. Its sufficiency
condition is $\gamma+\Delta_H<b_1$, so the largest admissible value is
$\gamma^\*=b_1-\Delta_H$; taking it gives the bound above. Certified in
`python-flint`/Arb at 256 bits, outward-rounded, by agent `Kg`:

```
b1_lower       = 0.88157382204959954858188766422240828   (half-width 2.46e-36)
gamma_paper    = 0.881545409                              (exact)
Delta_H_upper  = 1.59045454374832e-5  = head 1.13288599276897e-5 + tail 4.57569e-6
margin_low     = b1_lower - gamma_paper - Delta_H_upper = 1.2508504162e-5  > 0
gamma*         = b1_lower - Delta_H_upper = 0.881557917504162
K_G <= pi/(2 gamma*) = 1.781841323871362
```

**Owner-verified independently** (2026-08-30, this session, Arb at 300 bits):
$\pi/(2\gamma_{\rm paper})=1.7818666069360660912$, which reproduces the
paper's own `certificate.json` value $1.7818666069360661$ **exactly** — so
$\gamma_{\rm paper}$ is correctly identified and $\gamma^\*$ is real unclaimed
slack, not an arithmetic artifact. The re-derived $\gamma^\*$ bound is
$1.7818413238804372880$, and the resulting interval
$[6\pi/11,\ 1.78184132\ldots]$ is non-empty.

**Three caveats, load-bearing.** (i) This is **"PASS with tail input CITED"**:
the tail rests on $\lVert D^3H\rVert\le14.44243664663976457$ imported
`CITED-DEPENDENCY` from the authors' frozen `d3h_certificate.json`, as is their
certified A-grid `grid_M251.json`. The improvement is therefore
arithmetic-level — same scheme, same coefficients, tighter extraction of
$\gamma$ — not an independent certification of the scheme. (ii) The dominant
enclosure width is that imported tail, $4.58\times10^{-6}$, only $2.7\times$
below the margin rather than the $10\times$ the owner required; recorded as
`width_qualification` in the artifact rather than hidden, and re-deriving
$\lVert D^3H\rVert$ is the first task of gate C. (iii) Everything else — $b_m$,
head, tail arithmetic, margin, $\gamma^\*$, the $K_G$ bound — is `DERIVED` here.

Also machine-verified, and unchanged: the **tenths digit of $K_G$ is 7**, since
both $6\pi/11=1.713595992867159948\ldots$ and the upper bound lie strictly
inside $(1.7,1.8)$.

**All seven targets now carry committed evidence, and as of 2026-08-30 a second
wave of gate-C campaigns has landed.** `mm3/` and `oct-rank/` are complete
through gate C — `oct-rank/` additionally closed its **gate C4** $n=8$
peeling-constant probe, which had stood as the repo's last `HUMAN-AUDIT-PENDING`
item with *nothing filed*, as a certified negative. `omega/` has gate A's
verdict and gate C **PARTIAL**, but both historical gate-B rung certificates
are **SUSPENDED** after an outward-upper defect. Its separately repaired
five-point fixed-m ray replay is **MACHINE-VERIFIED**: every move worsens and
the verdict remains **OPEN**. It is not a feasible construction or exponent
result.
`mceliece/` has gate A on **13** instances with $m$ now **contiguous** $6$–$11$,
after the owner found and the agent closed a hole at $m=9$. `kg/` has gate A plus
an independently corroborated tail, and its gate C was held on a methodological
collision with newly-found literature, then cleared. `rs-pe3d/` has the first
exact 3-D expansion table; `delcap/` has a certified Blahut–Arimoto pipeline plus
two documented negative findings.

**Three currently accepted results are, to repo knowledge, firsts:** the
$K_G$ bound above; the first machine-checked per-orientation optimality map for
rank-23 $3\times3$ matrix multiplication (`mm3/`); and — added 2026-08-30 —
`mm3/`'s **monomial-transfer theorem**, that $C(F)=C(gF)$ for every
signed permutation $g$ of the 9 input coordinates, which makes the entire
$48^3$ sandwich group cost-invariant. Structurally, sign classes map
bijectively, input directions are fixed, and synthesis problems are isomorphic;
machine-confirmed by identical $d$, identical floor verdicts and identical DFS
state counts across 1,658,880 swept orientations. It discharges the monomial half
of gate C's contract **by proof** and localizes the remaining open direction to
the non-monomial ternary group, exactly $6960-48=6912$ generators per factor and
$\approx10^{12}$ orientations — a wall, not a sweep.

**Two second-order findings about the literature, both about checkability
rather than correctness:** the current $\omega$ world record has no released
parameters and so is not independently checkable from its paper as published
(re-audited 2026-08-30, still negative); and, by contrast, the $K_G$ paper ships
frozen Arb certificates, per-coefficient enclosures and an explicit margin. Same
field, same month, both AI-assisted, opposite ends of the artifact-quality
spectrum.

**A third literature finding, from the 2026-08-30 source audit
([details](docs/CS_FRONTIER_SCAN_2026-08-29.md#post-scan-audit--2026-08-30-owner)):
the repo's own frontier scan had two coverage defects, both from stating coverage
in a coordinate the source does not order by.** ECCC report numbers are **not
chronological** — 11 of 19 adjacent pairs in TR26-142..161 are date-inverted — so
the scan's "TR26-001..154" ceiling silently excluded seven reports dated 12–28
August; and the arXiv side missed two papers at the month boundary. Reading one of
those missed papers rather than filing it exposed **five** unknown 2026 works on
$K_G$, including a rigorous interval-arithmetic bound by Hermite thresholding at
exactly the degrees `kg/` gate C was about to sweep. **The miss was one paper; the
gap was five.** Rule now binding on every scan here: **state source coverage by
DATE against a dated listing, never by report or identifier number.**

**A field-level finding about printed tables, from `delcap/`'s three primaries
(2026-08-30).** Printed numeric tables in finite-length deletion-channel work
differ from certified recomputation for **four distinct reasons**, and **none
touches a theorem**: Morozov–Duman's Table III is a conservative round-up (a
single offset $\le7.1\times10^{-7}$ plus 5-dp ceiling reproduces all 36 rows);
Pinto–Ribeiro state theirs outright as "BA rate + tolerance"; Tavakoli–Nguyen–Bose's
Table I formula columns are *inconsistently* correctly-rounded (no single offset
is feasible — the spread $1.024\times10^{-3}$ exceeds one 3-dp ulp, where MD's
$9.689\times10^{-6}$ fell inside a 5-dp ulp, and the same criterion therefore
separates the two papers); and TNB's $C_{q,n}$ column is an unconverged
Blahut–Arimoto estimate. **So certified recomputation is not redundant with
reading the paper.** This replaces an owner hypothesis that conservative round-up
was a field convention, which the third paper falsified.

**Three apparent problems in the literature were raised and all three
dissolved** on investigation — see [Resolved threads](#resolved-threads-2026-08-30).
No refutation of any published claim is asserted anywhere in this repository.

### The price of rigour on the $\omega$ ladder (2026-08-30)

A second-order but reusable finding, and the first time anyone has measured it.
Both released rungs of the combination-loss line were enclosed here under
identical semantics, and both cost about the same to make rigorous:

| rung | published | certified enclosure at the released point | cost of rigour |
|---|---|---|---|
| Alman et al. 2025 | $2.371339$ | $\omega\le2.3713400836689$ | $+1.084\times10^{-6}$ |
| Vassilevska Williams–Xu–Xu–Zhou 2024 | $2.37155181$ | $\omega\le2.3715538358544617$ | $+2.026\times10^{-6}$ |

In both cases the gap is dominated by the **ln-space dual residual of the
shipped Lagrange multipliers** in the Lemma-1 max-entropy certificate — for
VXXZ24, worst $3.0767\times10^{-6}$ in glob region 2, $2.1502\times10^{-6}$ at
$r_0$, $1.1433\times10^{-6}$ at $r_1$ — measured against interval widths of
$1.08\times10^{-19}$. **The arithmetic is five orders of magnitude tighter than
the gap, so it is nowhere near the bottleneck.**

**Gate C Part 1 settles this by measurement, 2026-08-30 (agent `OmegaGateC`, owner-verified
arithmetic).** The above attribution was previously an inference from residual magnitudes. It is
now a falsifiable test that came back the hard way. Running the **same interval code path** at
MID precisions $300/200/128/96/64$ bits leaves the certified endpoint **bit-identical at all
five**, with the with-arm gap fixed at $+2.0258544615181506\times10^{-6}$ and the without-arm
arm $3.794237812826395\times10^{-9}$ *below* published at every precision; $R_{\rm sum}$ widths
move from $1.57\times10^{-89}$ to $1.73\times10^{-18}$ and the total rounding contribution is
bounded by $8.3\times10^{-19}$ of $\omega$ even at 64 bits. **A gap that does not move across
four octaves of precision is not an arithmetic gap.** Owner-recomputed at 300 bits: the frozen
enclosure sits $2.02585446\times10^{-6}$ above published, matching the agent's figure exactly.

**And the release is not defective.** Under our stricter certified standard the shipped
multipliers leave an ln-space Lemma-1 constraint-equality residual of $3.60\times10^{-11}$ —
comfortably *inside* VXXZ24's own stated tolerance of $1.1\times10^{-9}$. The release therefore
meets its own specification; the gap is the price of holding it to a stricter one. That
distinction — **stricter standard, not defect** — is why no escalation was warranted and none
was made.

**Gate C Part 2 — a retraction the agent made on itself, and a no-go that replaced it.** The
pre-registered box was $B=\{p^*+\delta:\lvert\delta_i\rvert\le10^{-7}\}$ over all $6759$
coordinates, stage (a) float then stage (b) interval on the same path; the $B_0$ anchor
reproduced at $2.3715538358544612$. Fourteen signed structured boxes were screened and **no box
certified below the published $2.37155181$, so there is no record claim.** Three boxes did appear
to certify below $B_0$ — `region_prop+` at $-8.15\times10^{-7}$, `glob_dist_2+` at
$-4.19\times10^{-7}$, `glob_dist_0+` at $-1.25\times10^{-7}$, the best closing $40.23\%$ of the
$2.026\times10^{-6}$ gap (owner-recomputed at 300 bits). **All three are RETRACTED: a
feasibility audit in exact dyadic arithmetic found every one of them INFEASIBLE.**
`region_prop+` breaks the $\sum r_p=1$ row by $+3.000000000086\times10^{-7}$;
`glob_dist_0+/1+/2+` each break their region's dist-simplex sum by $+7.0\times10^{-7}$ (seven
coordinates at $10^{-7}$) and also the dist$\leftrightarrow$dist$_{\max}$ marginal equality rows.
They are certified objective evaluations at exact-infeasible points and **bound nothing about
$\omega$**. The certified rung-2 value therefore remains the frozen
$\omega\le2.3715538358544617$. Retractions are inline in the campaign's
`feasibility_audit_addendum.md`, `stage_b_boxes_record.md` and `part1_gap_decomposition.md`, with
the README snapshot and checksums sealed.

**What replaced the retracted number is stronger than it was: an exact no-go on the `rp`
direction.** The $r_p$ block is pinned by **three** rows — $\sum r_p=1$ together with the
asymmetry rows $Y-Z=0$ and $X-Y=0$ — which force a uniform exact triple. Any asymmetric offset
violates an asymmetry row; any uniform offset violates the simplex row. The exactly-feasible
$\lvert\Delta r\rvert$ is $\approx5.5\times10^{-17}$, **ten orders below** the
$1.94\times10^{-11}$ defect absorption. So **no exactly-feasible `rp` move exists at radius
$10^{-7}$**, and the retracted $-8.15\times10^{-7}$ is unrealizable there in principle rather
than merely unrealized. Simplex-repaired `glob_dist` probes are the one remaining live path; they
were **not** run, being outside the pre-registered family, and are approved only with a written
pre-statement addendum fixing the repair rule before any evaluation.

**A latent defect in a validation accessor, found by this audit — provenance changes, numbers do
not.** `ParamManager._lin_beq` is **replaced**, not extended, on every `add_lincon_eq` call, so
`pm.linear_violations()` compares all $3174$ rows against only the **last** recorded $b$. The
stage-a screen's "lin $=3.0\times10^{-7}$" was an artifact of that ledger bug. This matters
beyond gate C: **both frozen rungs quoted "linear violations $0.0$" through the same accessor.**
Those claims are **RETAINED** — an independent exact-dyadic per-row replay at $p^*$ also returns
$0.0$ — but they now rest on the replay rather than on the accessor, which is a change in
provenance even though no number moves. The fix is applied in `src/` only; frozen campaign
scripts are immutable and stay untouched.

**Third instance this session of one failure shape.** A validation shares a component with the
thing it validates, and therefore cannot see a defect in that component: `kg/`'s
finite-difference audit checked a wrong function against the analytic derivative of *that same
wrong function*; `mm3/`'s three-way floor agreement is three solvers over **one** modelling
premise; and here an accessor compared every row to a single $b$. Recorded as a standing rule —
**independent instrument, or it is not an audit.**

**Gate C VERDICT — a certified negative, which is the outcome its contract named.** Under the
approved Decision-1 addendum a pre-registered simplex repair (exact rule: $+7\times10^{-8}$ on
the smallest coordinate, $-1\times10^{-8}$ on each of the seven largest, sum **exactly** zero,
owner-verified) was run on six boxes. **All six were rejected at feasibility before any objective
evaluation.** The repair reproduces $p^*$'s simplex residual *bitwise* — dyadic cancellation
works, residuals identical to the last digit — and $\mathrm{ceq}_{\max}$ stays exactly at the
$p^*$ level $3.5953524879506205\times10^{-11}$; but the dist$\leftrightarrow$dist$_{\max}$
marginal rows move by $\pm7.0\times10^{-8}$, worst $7.000097472856237\times10^{-8}$, which
exceeds the $1.1\times10^{-9}$ gate by a factor of **63.64** (owner-recomputed). The structural
reason is now proven rather than observed: the marginals are **linear** in dist while
dist$_{\max}$ was held fixed, so *any* strict-inside dist move at $10^{-7}$ scale breaks the 27
share-marginals rows.

**RETRACTED 2026-08-30, same day, by the owner — the local-optimality conclusion was an
overclaim and the framing error was mine, not the agent's.** Recorded per rule 5, not deleted.

**What was written and is now withdrawn:** "at radius $\pm10^{-7}$ every single-block direction
is infeasible", "glob-dist is jointly frozen by simplex-versus-marginals", "the released point is
the unique feasible point of the swept single-block family inside the box", and "local optimality
by constraint rigidity … none exists". The agent's original wording was the weaker
"certified-locally-optimal … within every feasible single-block family **tried**"; **I sharpened
it into the universal form**, which is the stronger and unsupported claim. Owner structural
judgement now stands at 0-for-13. The thirteenth was **H-COP**, reframed from the owner's own coprimality tabulation of the frozen `rs-pe3d` rows and killed by the P2 battery. **One counterweight, recorded at its earned strength:** the owner's `min(s)\le3` window-emptiness gate (H-GATE) **survived its first prospective test** — but the agent volunteered, unprompted, that the crisp iff-form crystallized only after four rows had landed, so its prospective content is **exactly one** committed read (PN10), not twenty. So: 13 falsified, 1 surviving on a single prospective test. The twelfth was a mechanism the owner falsified himself before the agent could act on it (the `mm3` per-factor alphabet collapse: 23.55% of non-monomial ternary unimodular pairs have a ternary product, so independent factors cancel generically and the diagonal $6960\to48$ result does not extend). **Three of the last four misses were attribution or mechanism, not arithmetic** — the owner's numbers survive checking, the owner's *reasons* do not. The eleventh, `mm3/` 2026-08-30, was again an **attribution** error rather than a wrong guess: the owner's commissioning brief asserted that the monomial-transfer theorem delivers the record attack's orientation collapse. It does not — the $145\times$ collapse is the ternary alphabet, an instance-exact census on the diagonal slice, and the theorem covers only the 48 monomials. The agent declined to inherit the framing. **Two of the last three owner misses were attribution rather than hypothesis, which is specific enough to act on: when the owner supplies a *reason why* something is tractable, that clause is the least reliable part of the instruction.** The tenth was a **misattribution** rather than a wrong hypothesis: a competent static audit of a **vendored external** Lean checkout, followed by adding it to a sibling repository's Targets table — which would have claimed an outside unconditional formalization of $H_1\le246$ as that repository's own result. Caught and reverted within the same edit sequence; the directory had **already** been documented as third-party by the owner hours earlier. Recorded as rule 17c.

**Why it fails.** Six pre-registered probe directions were tested — one $\pm$ repair pattern per
region. **Exact rank computation (agent, SVD, superseding the owner's estimate):** each
`glob_dist` block carries **45 coordinates**; the 27 margin rows are three 9-row blocks, each
full-rank 9, but their spans **overlap** so the stacked margin rank is **24, not 27**; and the
simplex row of ones lies **inside the margin span**, adding zero rank. Total **24** independent
equality rows, so the nullspace is **exactly 21** dimensions. (The owner's $45-28=17$ assumed the
28 rows independent — the very assumption the owner had instructed the agent to check and had not
checked itself.)

**A symptom of rank deficiency was read as a symptom of rigidity.** The repaired probes cleared
the simplex row *bitwise*, which both owner and agent took as the repair working correctly. It was
in fact a **trace of the dependency**: the simplex row is in the margin span, so any move
satisfying the margins satisfies it automatically. The kernel is simplex-orthogonal to machine
precision ($\max\lvert v\cdot\mathbf 1\rvert=5.55\times10^{-16}$), so simplex preservation on
the kernel needs no repair at all. **Finitely many rejected
directions cannot certify a positive-dimensional space.** The mechanism argument that was offered
in support — "the marginals are linear in dist while dist$_{\max}$ was held fixed, so any dist
move breaks the marginal rows" — proves only that *a dist move which changes the marginals breaks
them*, which is a tautology. The step actually required is that **every** dist move changes the
marginals, and that is false precisely when the marginal map has nontrivial kernel on the block —
which the dimension count says it does. Directions inside that kernel preserve all 27 marginal
rows and, if also summing to zero, the simplex row too; none was searched.

**What genuinely survives.** (i) **The `rp` no-go stands, and the owner re-derived it
independently rather than accepting it.** The region_prop block is the coordinate triple
$\{12,13,14\}$ carrying rows $\sum=1$, $x_{13}-x_{14}=0$, $x_{12}-x_{13}=0$, i.e.

$$A=\begin{pmatrix}1&1&1\\0&1&-1\\1&-1&0\end{pmatrix},\qquad \det A=-3\neq0 .$$

Rank $3$ on $3$ coordinates $\Rightarrow$ **nullspace dimension exactly $0$**, so the feasible
set is a single point, and solving gives it explicitly as the uniform triple
$(\tfrac13,\tfrac13,\tfrac13)$ — matching the agent's "uniform rational triple" exactly.
Exactly-feasible $\lvert\Delta r\rvert\approx5.5\times10^{-17}$, ten orders below the
$1.94\times10^{-11}$ absorption. **This is what a rigidity proof looks like: the rank closes the
space.** It is precisely the ingredient `glob_dist` lacks, where $45$ coordinates against $28$
rows leave $\ge17$ dimensions unsearched. (ii) **The six probes are
infeasible, as measured**, with the marginal rows moving $\pm7.0\times10^{-8}$ against a
$1.1\times10^{-9}$ gate (worst $7.000097472856237\times10^{-8}$, exceeding it $63.64\times$;
owner-recomputed). (iii) **No box certified below the published $2.37155181$**, so there was never
a record claim. (iv) Part 1's precision-invariance result is untouched — it depends on none of
this.

**Correct scope of gate C, FINAL:** *no feasible improvement was found among the six
pre-registered single-block probe directions at radius $\pm10^{-7}$, and the region_prop block
admits no feasible motion at all.* That is a real negative and it is all that was established.

**A full certification over the true domain was then attempted and FAILED TO CERTIFY — reported
as such, not patched and not shrunk to close.** Parameterizing by the coordinate perturbation over
$D=\{\delta:A\delta=0,\ \lVert\delta\rVert_\infty\le10^{-7}\}$, the naive two-sided
certified enclosure gives $\omega\in[2.3715431556360005,\ 2.3715624359361276]$ — width
$1.9280300127\times10^{-5}$ — against a first-order LP signal of
$\min G\cdot\delta=-1.5816497000997742\times10^{-7}$: the enclosure is **$121.9\times$ looser
than the signal** (owner-verified). Subdivision cannot rescue it, since widths are linear in radius
and closing the gap needs $(10^3)^{45}=10^{135}$ boxes. The missing machinery is named precisely: a
**certified gradient/slope enclosure** through the transcribed tree plus a second-order remainder.

**What the attempt nevertheless established.** An exact rational kernel basis with
$A\cdot V_{\rm exact}=0$ verified identically over $\mathbb Q$; the exact rank $24$ / nullspace
$21$ proof; the three named structural dependencies (one grand-total relation plus a 2-dimensional
**moment** family, the second total relation being their sum — so the laser-method margins are
*moment-coupled*, not merely total-coupled); a rigorous if weak enclosure valid over **all** of
$D$; and a reusable bound instrument requiring **no** dual-feasibility check, since $A\delta=0$
gives $G\cdot\delta=(G-A^{\!\top}y)\cdot\delta\ge-R\lVert G-A^{\!\top}y\rVert_1$ for
*arbitrary* $y$.

**And it reframed the open question.** The best conceivable gain over the true domain is
$1.58\times10^{-7}$, only $\mathbf{12.8\times}$ short of the gap to published — against
$810\times$ for the SVD-basis inscribed cube and $10\,677\times$ for the rref cube. Whether
`glob_dist` admits a feasible improving direction remains **OPEN**, but it is now one order short
rather than four, with the exact basis, the rank proof and the instrument all in hand.

**Caveat on the dimension itself, added 2026-08-30:** the count $21$ rests on a **float SVD**
rank decision ("three zero singular values"), which does **not** prove exact rank $24$. The risk is
asymmetric: a true rank *above* 24 only means the basis spans some infeasible directions, which
strengthens a cap; but a true rank *below* 24 means the kernel is **larger than 21** and the basis
**misses a direction**, so a cap over its span would be silent about that direction. The owner's
$17$ assumed row independence it had not checked; the agent's $24$ rests on a float rank test —
**same class of error, and the exact rank is still owed.** An exact rational nullspace of the
$27\times45$ coefficient matrix, plus an exact verification that the simplex row lies in the margin
span, has been ordered before any certification proceeds.

**The OPEN question is genuinely open, not merely unproven.** Agent-measured structure on the
21-dimensional kernel: every basis vector is two-sided at $p^*$ (support minimum
$1.28\times10^{-7}$ exceeds the margins), so $p^*\pm\epsilon v$ is feasible — box, open simplex
and all 27 marginals exact — for small $\epsilon$; the feasible set inside the box around dist is
therefore **at least 21-dimensional**. And the certified objective is **not constant** on that
kernel: the marginals are invariant, but $H(\mathrm{dist})$, the $p_{\rm comp}$ complete-split
terms and the `part_frac` flow into part-level blocks all vary, with float slack sensitivity along
kernel directions at $\epsilon=10^{-7}$ spanning $-2.5\times10^{-9}$ to $+6.8\times10^{-9}$ —
**both signs**. So descent directions exist at float level and a feasible certified improvement is
not excluded by invariance.

**But its size is bounded in advance, and this scoping is load-bearing.** At that sensitivity a
*fully successful* certification improves our own enclosure by $\approx10^{-9}$ at most. The gap
from our $2.3715538358544617$ to the published $2.37155181$ is $2.02585446\times10^{-6}$, a
**shortfall factor of $810.34$**; reaching the published value would need
$\epsilon\approx8.103\times10^{-5}$, three orders **outside** the pre-registered $10^{-7}$ box
(owner-computed). **A kernel search therefore cannot reach the published value, let alone the
unpublished record — its value is the certification, not the number.**

**`omega/sms` checked rather than assumed, and it lands on the second branch too.** No linear
equality row touches $\omega$ or `single_mat_size` at all ($K$ is pinned by bounds
$\mathrm{lb}=\mathrm{ub}=1$, outside box scope). The only acting constraint is the Schönhage line
$\mathrm{sms}\cdot\omega=T-R$ — **one** equation on **two** box coordinates, codimension 1, so
they are **not** frozen. But the total feasible $\omega$ excursion inside the box along that line
is only $\approx4.5\times10^{-8}$, under half the box radius, and any real motion must also keep
the $R$-side rows satisfied — which is precisely the compensating-$R$ multi-block class already
marked OPEN. **No free second no-go, and no live single-coordinate descent either.**

**Exactly what is NOT established, so "locally optimal" cannot be misread.** Coordinated
multi-block moves are **entirely unexamined** and are the live class: a dist move can only be
feasible if dist$_{\max}$ moves in the exact opposite linear pattern, and omega/sms moves need a
compensating $R$. Neither was pre-registered or run. The statement covers the **VXXZ24 rung at
$2.37155181$**, not the record $\omega<2.371177$, whose parameters remain unpublished as of
2026-08-30 and were therefore never searched. Construction-changing work — `2605.21738`
(asymptotic-rank speedups), `2608.27434` (centroids) — lies outside every statement here.

**Two accessor bugs, the second found *by* the mandated re-validation.** Fixing bug 1 (the
trailing-$b$ comparison) made the accessor immediately report $1.478849\times10^{-4}$ at $p^*$,
exposing **bug 2**: `_lc_rows` read `coeff[r,c]` with $r$ as the output row, silently
**transposing every matrix-coefficient row set** — precisely the j2m share-marginals equalities,
the only rows carrying real coefficient matrices. The MATLAB source was read first-hand
(`ParamManager.m:142-152`: `lin_Aeq(rows, group) = coeff'`, transposed) and the fix matches it
verbatim. **Bug 1 masked bug 2** — the trailing-$b$ comparison made the corrupt transposed rows
report clean, so the accessor was unreliable in a way that concealed a second unreliability.
Closed with **three** instruments agreeing on one number: corrected accessor, independent per-row
replay, and a direct structural recomputation of $\max\lvert Ad-Ad_{\max}\rvert$ over all 63
Parts $\times$ 3 margins, all at $3.5258684860650646\times10^{-11}$ — inside the release
tolerance. Rung 1 is untouched: `alman25_float` has **no** linear ledger at all (no
`AddLinearConstraintEq` in its transcription), so its frozen numbers came through the
$c/\mathrm{ceq}$ path and never met either bug. Stage-(a) headline numbers unchanged; the stage-b
anchor re-ran at $2.3715538358544617$ unchanged.

**This is a statement about released certificate data, not about the validity of
either bound, and not about any error in either paper.** The decisive evidence:
with tight multipliers the published VXXZ24 value **is** recoverable — the
no-arm diagnostic gives $2.3715518062057623$, below $2.37155181$ with margin
$3.7942377\times10^{-9}$. So the released artifacts are *incomplete* rather than
*insufficient*: the constants are float-optimal points whose accompanying duals
are loose at the $\sim10^{-6}$ level, and rigorously enclosing the released data
costs one to two parts in $10^6$.

**The positive result leads:** the **two-rung ladder order survives rigorous
arithmetic.** Rung 1 certifies to $2.3713400836689$ and rung 2 to
$2.3715538358544617$, a gap of $2.1375219\times10^{-4}$, so Alman25's rung
genuinely improves on VXXZ24's even when both are enclosed under identical
semantics with explicit slack — the first two-rung rigorous comparison in this
line, holding under both treatments tested, and both endpoints also sit below
DWZ23's $2.371866$. All figures owner-verified at 200 bits. Machinery validated
before use: rung-2 stage (a) clean at 6759/6759 parameters
($c_{\max}=4.63\times10^{-12}$, $c_{\rm eq,max}=3.60\times10^{-11}$, Schönhage
line exact), containment 9/9 num_block + 3/3 penalty + 1512/1512 Lagrange
form-consistency, precision sweep live
($1.08\times10^{-19}\to4.77\times10^{-18}$).

**Remaining certification gap — CLOSED 2026-08-30.** All three of `mm3/`'s
floor-impossibility results now carry proof logs accepted by a pinned external
checker: **kissat 4.0.4** DRAT proofs, verified by **drat-trim** (Heule/Wetzler
build 2024-04-21, commit `2e3b2dc`) **and** independently by `lrat-check`, in
about 0.02 s each. Frozen at
`mm3/campaigns/2026-08-30T031544Z_e0f3f117_c9df97a8bf3a/` with the CNFs, the
binary and LRAT proofs, the checker output verbatim, checksums, tool sources and
a manifest. So $C(U)=13$, $C(V)=14$, $C(W\text{-factor})=14$ are now
**checker-certified**, not merely agreed by three procedures.

**Encoding-adequacy control passed, and its shape matters.** A proof that a CNF
is UNSAT only certifies the floor claim if the CNF encodes the intended
question, so a SAT control was required. The raw floor encoder is UNSAT at both
$T=d$ and $T=d+1$ — which is correct rather than a bug, because the paper's
witness circuits use one auxiliary intermediate — and the aux-1 encoding is
**SAT** at $d+1$, exhibiting the witness class. That is consistent with agent
`Mm3`'s independent aux-1 finding, and it is the control in the only form that
is meaningful here.

The earlier pseudo-Boolean blocker was routed around rather than fought: VeriPB
v0.1.0's `recordclass>=0.15` API break made the PB route the fallback, and the
DRAT route (already proven in this workspace by `../math/`'s Kobon result)
carried it instead.

### Proof-log certification (2026-08-30)

Per-instance detail for the closure recorded above. Solver
`OMP_NUM_THREADS=1 nice -n 10 kissat <cnf> <proof.bin>`; checkers `drat-trim`
and `lrat-check` from the same pinned tree, sources frozen in the campaign's
`tools_snapshot/`.

| instance | vars | clauses | proof (bin) | proof (LRAT) | core clauses | verdict |
|---|---:|---:|---:|---:|---:|---|
| paper-$U$, $d=12$ | 396 | 1440 | 719 B | 6111 B | 13 | `s VERIFIED` / `c VERIFIED` |
| paper-$V$, $d=13$ | 468 | 1807 | 148 B | 7946 B | 14 | `s VERIFIED` / `c VERIFIED` |
| paper-$W$-factor, $d=13$ | 481 | 1820 | 668 B | 8011 B | 14 | `s VERIFIED` / `c VERIFIED` |

One core lemma each; ~0.02 s per check. Cross-checked with a second solver:
**cadical 3.0.1** emitted ASCII DRAT for paper-$U$, also verified. **Encoding
soundness stated explicitly:** UNSAT of the over-approximating floor encoding
implies no circuit exists at $T=d(F)$; the encoding mirrors the paper's floor
lemma (every gate value lies in $\pm$ the needed target set), with clause groups
for at-most-one-per-slot, coverage, $y\to\bigvee p$, and $p\to$ operand
availability. **SAT control:** the aux-1 extended encoder is satisfiable at
$d+1$ — paper-$U$ at $T=13$ (533 vars / 1975 clauses) and paper-$V$ at $T=14$
(616 / 2379) — admitting the paper's own witnesses, whose auxiliary
intermediates are exactly $u_{12}=A_1-A_4$ and $v_9=v_5+B_7$.

## Results table

| target | problem | status | headline result | evidence |
|---|---|---|---|---|
| [`mceliece/`](mceliece/README.md) | Classic McEliece hold-out "waterfall": is the derivative-flag collapse at $c=2t+3$ real on binary Goppa? | **gate A PASS on 13 instances, $m$ CONTIGUOUS $6$–$11$, $\delta$ guard UNIFORM and exact; **gate B FROZEN 2026-08-30 — an exhaustive THEOREM at $(m,t)=(6,2)$ plus a bounded genericity rate; Apon's §3.6 hole bounded, NOT closed** | **Apon's mechanism is confirmed by measurement, exactly, with zero refutation events.** Small-$t$ artifact: seven instances $(m,n,t)$ = $(6,64,3)$, $(6,64,4)$, $(6,64,5)$, $(7,128,3)$, $(7,128,6)$, $(8,256,3)$, $(8,256,5)$; ladder addendum: $(m{=}10,t{\in}\{24,40\})$ and $(m{=}11,t{=}48)$ on full-field supports with the same engine, seeds and rules. On **every** instance and **every** held count $c$: $N_{\rm fam}(c)=2t+3-c$ **exactly**, reaching $0$ precisely at $c=2t+3$; every per-(point, order $j<4$) family-restricted block has rank **exactly 1**. Two quantities are measured and compared, not one: $N_{\rm fam}(c)$ from all explicitly-built family-restricted flag rows, against $N_{\rm eval}(c)$ from the evaluation/Vandermonde matrix on the family, which *is* the Lemma-7 prediction. **Why this has force:** the family-restricted rows are built explicitly from Apon's Definition-4 map $T(A)=Q_A+\rho_A Z F$ (eqs. 17–20), so Lemma 7's claim — that flag conditions on the family reduce to $A(a)=0$, hence at most one condition per held point — is **derived by measurement, not assumed**. Five machine-checked instance guards per row: Apon's differential identity $\Pi F'+\Pi'F=\kappa G^2F^{(2)}$ coefficientwise; $\Delta_{p,q}\ne0$ with the exhibited pair recorded; $k=n-mt$; $\lambda_i F(a_i)$ exactly matching the public $Y$ columns in $\mathbb F_2^k$; constant coordinate gcd with max coordinate degree exactly $D$. Convention risk closed in advance by adopting Apon's own §2 perfect-oracle protocol, with Saarinen's eqs. (29)/(30) as cross-check anchor and his relation-generation reductions correctly excluded from Step-3; ambient rank is never measured, so the spurious-refutation mode is designed out. | MACHINE-VERIFIED (exact $\mathbb F_{2^m}$ arithmetic, uint16 bitmask field with multiplication/inverse tables; **zero ball-arithmetic paths**), cross-validated against `python-flint` `fq_default` (500/500 scalar products, 100/100 polynomial products, 8/8 irreducibility verdicts on $\mathbb F_{64}$ and degree-8 over $\mathbb F_{4096}$). Two frozen campaigns: `2026-08-30T01-42Z_9D0FC9E7` (small $t$) and `2026-08-30T04-47Z_98C586C7` (ladder addendum), each carrying the committed pre-statement, polynomials, supports, seeds and guard outcomes. **Scope, per rule 9:** a finite theorem over 10 instances. The $m=12$ rows ($n=3488$, $t\in\{48,64,96\}$) are compute-bound at $\approx1$ h/row in the current engine and continue in a further artifact — a **performance** limit, not an algebraic one. **$c_{\rm need}>2t+3$ at `mceliece8192128` scale remains NOT established**, and the artifact says so. **Self-corrections preserved, not deleted** (`PROBE_RESULTS.md`): a Rabin irreducibility test using Frobenius by $2$ rather than $q=2^m$, which silently passed split polynomials over extension fields and would have voided the squarefreeness guard; and the wrong GRS multiplier $\lambda$ instead of $\lambda'=G(a_i)^2/\Pi'(a_i)$, which interpolates a different ambient code. The conclusion drawn from the second bug ("$F'\equiv0$ identically") is explicitly **VOIDED** there. The corrected analysis shows $F$ is unique (evaluation injective since $D=n-2t-1<n$), so $F'\equiv0$ is impossible for a nonzero row — which makes Apon's $\Delta_{p,q}\ne0$ a **genuine** genericity condition and gate B's actual subject: a measured degeneracy rate over seeds and support orderings, a quantity neither disputant has. |
| [`kg/`](kg/README.md) | Grothendieck constant $K_G$ | **gate A PASS, tail input INDEPENDENTLY CORROBORATED; **gate B PARTIAL and two prior band PASSes RETRACTED 2026-08-30** — predecessor envelope was unsound in the permissive direction; gate C COMPLETE — certified cap on the septic family + a certified ceiling on the degree-5 program** (agents `Kg`, `KgLeg2Conv`, `KgGateC-2`) — gate C closed 2026-08-30: the septic family does **not** beat $3.47\times10^{-4}$ and the paper's own $(s_3,s_5)$ at $s_7=0$ is the best certified point in the sweep; separately, a hyperplane-attained $b_5$-strip is **provably impossible** for $\lambda>\lambda^*=0.9362333111198615$ | **See the headline.** $K_G\le1.7818413238804372880\ldots$, stronger than arXiv:2608.11158's stated $1.7818666069360661$ by $2.528\times10^{-5}$; improvement over Krivine $3.72654\times10^{-4}$ against the paper's claimed $\ge3.4737\times10^{-4}$. Certified at 256 bits outward-rounded: $b_1^{\rm low}=0.88157382204959954858188766422240828$ (half-width $2.46\times10^{-36}$), $\Delta_H^{\rm upper}=1.59045454374832\times10^{-5}$ = head $1.13288599276897\times10^{-5}$ (`[DERIVED]`, half-width $1.77\times10^{-36}$) + tail $4.57569\times10^{-6}$ (Lemma 6.2 at $N=251$ with the **cited** $\lVert D^3H\rVert$), giving $\mathrm{margin}_{\rm low}=b_1^{\rm low}-\gamma_{\rm paper}-\Delta_H^{\rm upper}=1.2508504162\times10^{-5}>0$ and $\gamma^\*=0.881557917504162$. **CORRECTION to this row's earlier text, from the agent's final report:** it previously said an independent single-cell computation gave $b_1=0.881598$, making the paper's $0.881573822049$ "a rigorous lower bound $2.4\times10^{-5}$ below the true value". **That was wrong** — the $0.881598$ figure was a *partial* evaluation missing the $A_{0,1}$ term. The correct closed-form identity $b_1=\frac{\pi}{2}\!\left(A_{1,0}^2/V-A_{0,1}^2\right)$ reproduces $b_1^{\rm low}$ to $2.34\times10^{-36}$ from the imported grid bytes and matches the paper's printed value to $4.86\times10^{-17}$ **absolute**. So the paper's $b_1$ is essentially exact, not conservative; the unclaimed slack is entirely in their choice of $\gamma$ ($0.881545409$ against an admissible $0.881557917504162$), which is where the whole $2.528\times10^{-5}$ improvement comes from. The headline is unaffected; its mechanism is now stated correctly. **Gate B — PARTIAL, and explicitly not a refutation.** Interval **cover** (no sampling) of the $(p,c)$ threshold family per the declared design (even-coefficient box subdivision, Christoffel envelopes dominating odd parts, Lemma D.4 monotonicity). Certified bands: $c\in[1.30,1.45]$ (856 boxes, $J(c,p)\le d(c)$ with positive certified margin on every box), $[3.50,4.083]$, $[4.083,6.0]$. Open bands where the agent's envelope did not close within budget: $[1.0,1.3]$, $[1.45,1.75]$, $[6.0,12.0]$. A spot-check confirms the paper's inequality itself holds there — corner cubic $a_0=-0.98$, $a_2=-0.05$, odd norm $0.105$ gives $J(6)=0.0605\pm3.1\times10^{-14}\le d(6)=0.0991$ — so **the loose link is the box certificate, not the paper's Lemma D.2/D.4**. Reported OPEN rather than faked-pass. | MACHINE-VERIFIED for the derived chain. Gate A frozen at `kg/campaigns/20260830T010702Z_2a373482_385dbaab2fd3/`; gate B at `kg/campaigns/20260830T013720Z_e24aa159_afa515bbc532/`. **CITED-DEPENDENCY inputs** (SHA-256 in `provenance.json`): `grid_M251.json` and `d3h_certificate.json` ($\lVert D^3H\rVert\le14.44243664663976457$) from `github.com/trishullab/grothendieck-bounds` — so the improvement is arithmetic-level on the authors' own scheme, **not** an independent certification of it. **Interval hygiene verified** with the corrected diagnostic: precision sweep at 512/256/128/96 bits gives widths $6.10\times10^{-121}$, $4.85\times10^{-78}$, $1.65\times10^{-39}$, $7.08\times10^{-30}$ — monotone growth with margin positive throughout, so no float leak; taint test confirms A-grid half-width inflation propagates to the $b_1$ output width; a two-argument `arb(lo,hi)` audit found **zero** misconstructions in executed paths. **Stated width qualification:** dominant width is the imported tail $4.58\times10^{-6}$, $2.7\times$ below the margin rather than the required $10\times$; kept as `width_qualification`. **Cost correction, attributed to the agent, replacing the owner's estimate:** $b_m$ loop measured 544.8 s (Apple M3 Ultra, `nice -n 10`, python-flint 0.9.0, 256 bits); fresh-user end-to-end $\approx$10 min reusing the authors' certified A-grid; deriving the A-grid from scratch is a separate larger task, not completed. $K_{\rm DR}$ `[REPORTED]` (`bound2.dvi` dead, no Wayback 200; Davie 1984 unobtainable) — no gate depends on it. Gate C not started; handed to a fresh agent. |
| [`omega/`](omega/README.md) | matrix multiplication exponent $\omega$: is the record rigorously established as stated? | **gate A verdict + gate B rungs 1 & 2 CERTIFIED + gate C PARTIAL** (agents `Omega`, `OmegaRung2`, `OmegaGateC`) — gate C 2026-08-30: no feasible improvement among **six** pre-registered single-block probes at $\pm10^{-7}$; exact `rp` no-go (3 coords, 3 rows, zero-dimensional); three infeasible probes retracted; two accessor bugs fixed with provenance recorded. **A wider local-optimality claim was written and RETRACTED the same day — owner framing error**; `glob_dist` has a $\ge17$-dimensional unsearched nullspace, so its freeze is **OPEN** | **Gate A — the current world record is not independently checkable from the published paper, dated 2026-08-29.** arXiv:2608.16884v1 gives the reformulated program (Eq. 11, $q=5$, $\ell^\*=4$, ~7M parameters, Lemma-1 $H^{\max}$ scheme) but **no parameter table and no ancillary files**, stating only that the authors *"are preparing a repository in which we will release the verification code and our discovered solution."* Six-path auditable negative: e-print tarball carries no ancillary data (`assets/` = 3 logo PNGs), TeX has no parameter table, all ten co-authors' GitHub accounts enumerated, `google-deepmind/alphaevolve_results` holds only the May-2025 notebook, GitHub full-text search for `2.371177` = 0, no OSF project, web text search 0. **A dated statement about artifact availability, not a claim the bound is wrong.** **Gate B rung 1 — the first rigorous interval enclosure of any rung of the combination-loss ladder.** At the exact released Alman25 parameter point (24855 float64 values, SHA `f2336913…b46f`): $$\omega\;\le\;2.3713400836689$$ = raw interval endpoint $2.371340083602922$ plus $6.5978\times10^{-11}$ of absorbed feasibility slack. $\varepsilon$ decomposition, **corrected** after an owner-flagged double-count (see evidence): the Lemma-1 max-entropy residual of the shipped Lagrange multipliers contributes $1.0836\times10^{-6}$ and is carried **inside** the raw endpoint, not added again; the three stage-(a) feasibility defects — inequality violation $1.137\times10^{-10}$, equality violation $2.390\times10^{-11}$, Schönhage half-ulp $-6.2\times10^{-15}$ — absorb to $6.5978\times10^{-11}$ total. Region-wise Lemma-1 residuals: $1.33,1.05,0.74,0.63,2.06,1.64\times10^{-6}$, max $2.062\times10^{-6}$. **Consequence, stated plainly:** the certified endpoint sits $2.1173\times10^{-4}$ **below** VXXZ24's published $2.37155181$, so Alman25's rung genuinely improves on its predecessor even under rigorous arithmetic with explicit slack — independently of the unreleased $2.371177$ artifact. It correctly does **not** certify $2.371177$, and it remains $1.0837\times10^{-6}$ **looser** than the published display $2.371339$, which is what an honest enclosure of a float-optimal point costs. | Gate A HUMAN-AUDITED, frozen at `omega/campaigns/2026-08-30T00:26:00Z__4528cddb80c2/`. Gate B rung 1 MACHINE-VERIFIED, frozen at `omega/campaigns/2026-08-30T03:20:00Z_c9d4e2a8/` (republishable JSON, with a dated correction note). **Owner-flagged and owner-verified correction:** the first reported endpoint was $2.3713411672715115$, which **conservatively double-counted** the Lemma-1 gap — `raw − published` and the applied slack were near-identical ($1.0836029$ vs $1.0836686\times10^{-6}$), which the owner flagged as too coincidental to be independent. The agent's decisive test toggled the Lemma-1 charge and showed the raw endpoint moves by exactly that amount, proving the gap was already inside it. Double-counting slack can only loosen an upper bound, never invalidate it, so nothing was ever wrong — but the corrected slack term is **16,424.7× smaller** and the endpoint improved by $1.0836026\times10^{-6}$, which also widened the margin below VXXZ24 from $2.1064\times10^{-4}$ to $2.1173\times10^{-4}$. All figures re-derived independently by the owner at 200 bits. **Interval hygiene, now repo knowledge:** an earlier apparent zero-width collapse was `float()` JSON printing, not a real collapse (entropies carry genuine $1.8\times10^{-20}$ widths at Arb precision); the original phantom's root cause was python-flint's `arb(lo, hi)` constructing (midpoint, **radius**) rather than an interval — broadcast to all agents, and `Delcap` then found the complementary trap that a tight radius reads as $0.0$ if width is computed by subtracting float64 endpoints. Per-block containment 18/18 + 6/6 + 5/5 + 1/1; taint test (`dist[r0]` → `num_block` width $3.09\times10^{-8}$) confirmed every value-chain edge live. DWZ23's rung remains MATLAB MCOS-opaque and unevaluated. |
| [`oct-rank/`](oct-rank/README.md) | real tensor rank of octonion multiplication, window $18\le \mathrm R_{\mathbb R}(T_{\mathbb O})\le25$ | **gate A PASS; gate B clean negative; gate C PASS; gate C4 FILED as a certified negative; **S3 CLOSED as a filed negative 2026-08-30 — best certified LB 13, gap to 14 exactly 1, and the substitution route PROVEN unable to exceed 13**** (agents `OctRank`, `OctRankGateC`) — the $n=8$ peeling-constant probe, previously `HUMAN-AUDIT-PENDING` with **nothing filed**, is closed 2026-08-30 | **Gate A — arXiv:2608.16649 independently re-verified end to end.** Exact-rational re-check of the rank-25 Krawczyk certificate: $K=0.319556737<1$, $\mathrm{off}_{\max}=1.314\mathrm e{-14}$, $(1-K)\rho-\mathrm{off}_{\max}=+6.804\mathrm e{-7}>0$ at $\rho=10^{-6}$, **all 512 per-equation margins strictly positive**; own Cayley–Dickson tensor matches the upstream table entrywise; $\lVert YJ_S-I\rVert_\infty=2.978\mathrm e{-12}$; outward-rounded Arb re-evaluation strict (worst-margin ball $[6.8044\mathrm e{-7}\pm9.89\mathrm e{-22}]$). Radius re-derived `[DERIVED]`: hypotheses hold for all $\rho<\rho^\*\approx3.1293\mathrm e{-6}$ by 60-iteration bisection, contraction-limited rather than off-vector-limited, so the paper's $\rho=10^{-6}$ carries $\approx3.1\times$ slack. Lean 4 replay: `leanprover/lean4:v4.29.1`, mathlib `5e932f9`, `lake build` 2398 jobs clean; axiom audit of every key theorem = exactly `[propext, Classical.choice, Quot.sound]`, no `sorry`; `oct_rank_ge_eighteen : IsSlicesOfRank r T_O → 18 ≤ r` matches the prose. **Gate B — no rank-24 candidate found, and the negative is well-powered.** 273 probes: 256 random seeds (best rel $2.631\mathrm e{-3}$, max factor-column norm $30.4$, all 300 LM rounds consumed), compress-mode merges (best $1.76\mathrm e{-2}$), small-alphabet snap (rel $\approx3$), and a deflation-continuation curve from the frozen rank-25 point. Control for contrast: $r=25$ reaches rel $1.85\mathrm e{-14}$ with norms $\le9.5$ in 287 rounds. **Border-rank discrimination was decisive and is why this negative is trustworthy:** on the continuation, residual stays $\approx1.3\mathrm e{-15}$ while norms grow $13.7\to1135.7$ as $s\to0.01$ — but those rows carry **25** columns, so that is inter-term re-weighting compensating a shrinking term, *not* border-rank drift toward a 24-term limit; norms snap back at $s=0$, and the true 24-column endpoint runs away ($0.563\to0.519$) with bounded norms. An earlier $s=0.9$ run that reached rel $1.747\mathrm e{-15}$ was **self-escalated and correctly rejected** by the agent as re-convergence to the rank-25 certificate, since scaling a column by $0.9$ does not remove the term. **Gate C — PASS.** $\mathrm R_{\mathbb R}(\tau)=7$ independently re-verified in both directions (exact, worst margin $9.946\mathrm e{-6}$, $K=0.005352$), and the general $\frac52n-2$ lower bound confirmed **sharp** at both small even dimensions by exact witnesses: $n=2$ rank 3, $n=4$ rank 8. | Gate A MACHINE-VERIFIED (re-verification) `[REPRODUCED]` / CITED-DEPENDENCY for the mathematics, frozen at `oct-rank/campaigns/2026-08-30T00:04:39ZZ_gateA_A6A4B0CB/`. Gate B **COMPUTATIONAL-EVIDENCE about the search under a declared budget only** — explicitly **no** claim that $\mathrm R_{\mathbb R}(T_{\mathbb O})>24$; the window $18..25$ remains open at both ends. Gate C MACHINE-VERIFIED. The $n=8$ peeling-constant probe is HUMAN-AUDIT-PENDING and stays in `scratch/`. |
| [`rs-pe3d/`](rs-pe3d/README.md) | high-dimensional product expansion for Reed–Solomon tensor codes (TR26-150 Conjecture 4.2) | **gate A first exact 3D table frozen; gate B COMPLETE — 13/13 exact rows, verdict NEITHER; **P2 battery COMPLETE 2026-08-30 — pattern P survives 9/9 on the strong arm, H-COP falsified, and the window's domain gate identified**** (agents `RsPe3d`, `RsPe3dGateB`) — gate B closed 2026-08-30: $\beta$ is **not** the driver, falsified at the opposite extreme from its prediction | **The first exact 3-dimensional product-expansion ratios ever computed**, 9 instances, each $\delta$ exact over $\mathbb F_q$ with an exhaustive cheaper-support check. Window: all point supports of weight $\le3$, normalized reduced-basis candidates. Ratios — $q=13$, $s=(2,2,4)$: $1/2$, $3/8$, $1/4$ at $t=(1,1,1),(1,1,2),(1,1,4)$; $q=31$: $s=(2,3,5)\to3/5$, $(2,3,10)\to3/5$, $(2,6,5)\to1$, $(3,10,2)\to3/5$; $q=61$: $(2,3,5)\to3/5$, $(3,4,5)\to1$. **Honest labelling:** every row attains its best $\delta$ via a weight-$\le3$ candidate but none is support-proven globally minimal over the $M$-side, so each value is $\rho_{\rm inst}^{\rm upper}$ and the global $\rho$ stays OPEN per row. **Balance wedge, first data:** at $q=61$, $\beta=1.67\to\rho=1$ versus $\beta=2.5\to\rho=3/5$ — visible degradation as imbalance grows, explicitly **not** a collapse claim. **RETRACTED 2026-08-30 — this framing was the repository's error, not the paper's.** It said the proved 2-D Theorem 2.1/4.11 carries a $\lvert S_1\rvert\approx\lvert S_2\rvert$ hypothesis "that Conjecture 4.2 drops for higher $r$". **Conjecture 4.2 does not drop it.** First-hand read by agent `RsPe3dGateB`, quoted with line numbers in `campaigns/…_gateB/finding2_balance_hypothesis_provenance.md`: Conjecture 4.2 (lines 740–748) carries *"balance parameter $\beta\ge1$ … $\beta^{-1}\le s_i/s_j\le\beta$"*; Theorem 4.11 (1280–1285) carries $K^{-1}\le n_1/n_2\le K$; Lemma 5.4 item 2 (1632) carries $1/4\le s_i/s_j\le4$; Section 6 (1604) states the hypotheses require side lengths *"pairwise coprime and balanced"*; and **even the unverified AI proof carries it** (its Theorem line 8 $K^{-1}\le n_i/n_j\le K$, Prop 5 line 252 "$K$-comparable", constant $\rho=(2dK^d)^{-1}\bigl(\varepsilon/(C_0dK)^{10d}\bigr)^{d/\gamma_d}$ with $d/\gamma_d=3d\cdot6^{d-2}$, which at $d=3$ degrades as $\mathbf{K^{-1623}}$ — `[DERIVED]` from the printed display, unverified proof, no claim about its truth. **An earlier $K^{-4863}$ here is RETRACTED, and the owner's "confirmation" of it was worse than the error itself.** The exponent is pinned by the proof's own definitions: line 245 sets $\gamma_1=1$ and $\gamma_e=1/(3\cdot6^{e-2})$, line 489 gives the recurrence $\gamma_e=\gamma_{e-1}/6$, and line 933 states the exponent *is* $d/\gamma_d$. Exact arithmetic then yields $6,54,432,3240$ at $d=2,3,4,5$, matching the **linear** $3d\cdot6^{d-2}$ at every $d$ and excluding a tower $3^d\cdot6^{d-2}$ ($9,162,2916,52488$). The owner instead inferred the reading **from the agent's stated K-power** — a quantity derived out of the very formula in question — so the check **matched the error instead of catching it**: two apparently independent confirmations of one wrong number, because one was derived from the other. Caught by the agent re-deriving from the definitions. **This is discipline rule 14's shared-component failure, committed by the owner while enforcing it.**). Balance is **absent only from Conjecture 1.4's informal phrasing** (lines 224–226, *"$\rho$ being only a function of $\varepsilon$ and $k$"*). So the gap is **informal restatement versus every formal statement**, the wording correction belongs to Conjecture 1.4, and Conjecture 4.2 needs no correction on this axis. No claim is made against any published statement. | **Gate B, 2026-08-30: the balance question is answered in the negative, and it is the contract's third outcome.** 13 of 13 pre-registered rows exact, $\beta$-degradation **falsified at the opposite extreme** — it predicted $\rho(\mathrm{B11})<1/4$, below the grid minimum, and the observed value is $1$, the grid **maximum**. At $q=241$ in increasing $\beta$: $5/3\to1$, $5/2\to3/5$, $16/3\to1$, so $\beta$ grows $3.2\times$ with $\rho$ returning to its maximum while the *intermediate* $\beta$ sits lowest — no monotone function of $\beta$ can do that. **Conclusion: on the swept set $\rho^{\rm window}$ is not a function of $\beta$; there is no $\beta^*$ to locate.** A candidate pattern, owner-tabulated and pre-registered with a cryptographic order proof before the deciding row landed, is that $\{2,3\}\subseteq\operatorname{orders}(s)\iff\rho=3/5$ — it fits all 13 rows and survived **one** out-of-sample test. Recorded as evidence, **not** established: 12 of 13 rows are in-sample, the $1/\binom{13}{5}=7.8\times10^{-4}$ separation figure is therefore an in-sample statistic, and the three $q=13$ rows show $\rho$ varying with $t$ at fixed orders ($1/2$, $3/8$, $1/4$), so the pattern cannot be the whole story. Falsifier named. Frozen at `rs-pe3d/campaigns/2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB/` with the pre-statement committed before compute and witness certificates regenerated inside the freeze. Gate A frozen at `rs-pe3d/campaigns/2026-08-30T01-04-49Z_exacttable`. Exact $\mathbb F_q$ arithmetic; per-row optimality attestation being strengthened by exhaustive line-support minimization. Two premise failures documented as dead ends rather than deleted: a $q^{22}$ enumeration that was infeasible, and a false "Anchor 1" equality disproved by an exact counterexample certificate (the agent's own derived unit premise, not a claim of this repo or of TR26-150 — attribution corrected in its `pre_statement.md`). |
| [`mm3/`](mm3/README.md) | additive complexity of rank-23 $3\times3$ matmul; is the per-orientation optimality claim correct, and can 55 be beaten? | **gates A+B+C complete; UNSATs checker-certified; gate C extended with a THEOREM; OFF-DIAGONAL PROGRAMME COMPLETE 2026-08-30 — 3,387,432,960 valid orientations all certified $\ge55$, zero $\le54$, and the published orientation is the UNIQUE minimizer** | **Gate A — arXiv:2607.28676 reproduced exactly** `[REPRODUCED]`: all 729 Brent identities over $\mathbb Z$ (27 unit, 702 zero, **0 failures**), counts $13/14/28=55$, three independent presentations agreeing, Perminov's `cr58_cn122` provenance closed to 0 mismatches. Full predecessor ladder also re-verified in one fixed convention: Stapleton-60 (729/729 + 2000-trial randomized end-to-end, printed-structure recount exactly 60), MWS-59 (729/729 on both their Table-3 artifact and a reconstruction of printed Table 2), Perminov-58, Sun-56. **Gate B — the paper's "provably optimal for this fixed orientation" claim is machine-decided, CONFIRMED, and now CERTIFIED.** $C(U)=13$, $C(V)=14$, $C(W\text{-factor})=14$, output $=28$, established by **three independent decision procedures** (exhaustive floor DFS, HiGHS ILP infeasible, CaDiCaL 1.5.3 UNSAT) plus the paper's own circuits as $d(F)+1$ witnesses, and since 2026-08-30 by **kissat 4.0.4 DRAT proofs accepted by drat-trim and `lrat-check`** (see [Proof-log certification](#proof-log-certification-2026-08-30)). Certificate structure after the agent's own retraction of an over-reading: **[no circuit at the floor $d(F)$] + [witness at $d(F)+1$]** $\Rightarrow C(F)=d(F)+1$; aux-1 at $d+1$ *is* possible, as the witnesses require and as the SAT control confirms. **Gate C — no $\le54$ total exists in the valid orientation orbit, and the landscape is mapped.** Swept set stated exactly: the $\sigma$-orbit $\{\mathrm{Id},\sigma,\sigma^2\}$, $\sigma:(U,V,W)\mapsto(V,W^\top,U^\top)$ of order 3; pure $(U,V)$ swaps **excluded by proof, not budget** ($B\cdot A\ne A\cdot B$). Certified totals: $\mathrm{Id}\to13/14/28=\mathbf{55}$, $\sigma\to14/14/29=57$, $\sigma^2\to14/13/30=57$ — the paper's orientation is the **certified optimum of its own orbit**, and $\{55,57,57\}$ is the **first machine-checked per-orientation optimality map for rank-23 $3\times3$**. Sun's decomposition separately certified on *his* factors at $13/13/30=56$. **Two cross-paper alarms raised and both dissolved:** Sun-56 is a *different* rank-23 decomposition, not a reorientation (zero common product triples, and multiset equality fails up to signs, under the $T$ involution, and transposed); and the 59-vs-57 flag was the agent's **own tally slip**, machine-recounted to exactly $15+15+29=59$, retracted inline, with record history $60\to59\to58\to56\to55$ unchanged. Bonus: the same procedure independently confirms Sun's own lower bounds, including his reported `V_aux1_at_12_possible: False` over all 338 auxiliary directions. | Gates A/B/C MACHINE-VERIFIED across **six** frozen campaigns, the last being `mm3/campaigns/2026-08-30T031544Z_e0f3f117_c9df97a8bf3a/` (CNFs, binary + LRAT proofs, verbatim checker output, checksums, tool sources, manifest) and `2026-08-30T012035Z_a7ea8e8e_cbdf8e63fa94/` (the 59-vs-57 resolution). **Scope, per rule 7:** gate C is *complete* over the valid $\sigma$-orbit but *partial* over abstract orientation space, so "no 54" is not a universal claim. **Certification gap CLOSED:** the UNSATs no longer rest on procedure agreement alone. VeriPB v0.1.0's `recordclass` API break made the pseudo-Boolean route a fallback; the DRAT route carried it instead, with cadical 3.0.1 as a second solver on paper-$U$. |
| [`delcap/`](delcap/README.md) | binary/$q$-ary deletion channel: is the published constant chain rigorously enclosable, and can the finite-length bounds be tightened? | **gate A closed at max certifiable $n$; gate B enclosure table delivered; **Open Item 1 RESOLVED 2026-08-30 — 20/20 $q{=}3$ rows strictly improve BOTH ends of the published TNB sandwich**; gate C TARGET ACHIEVED — certified improvement on the published $d=1/2$ sandwich, both ends, all 13 rows** (agents `Delcap`, `DelcapGateC`) — **row added 2026-08-30: this target had NO row in this table while the summary above asserted its evidence, the same SSOT drift found in `../physics/` the same day** | **GATE C TARGET HIT: a certified finite-blocklength improvement on Tavakoli–Nguyen–Bose's published $d=1/2$ sandwich, at BOTH ends, on all 13 rows** ($q\in\{2,3,4\}$, $n$ up to 8). Truncation order **NONE** — the full exact $q^n\times\sum_k q^k$ channel enumerated, no window, no alphabet cut — and tail bound **exactly $0$**, Arb 400 bits. Their sandwich is $\mathrm{LB}_1=(1-d)\log_2q-h_2(d)$, $\mathrm{LB}^+$ their Cor. 1 tightening, $\mathrm{UB}=(1-d)\log_2q$. Example rows: $q{=}2,n{=}2$ their $[0.375,0.5]$ against certified $[0.415241012,0.415241013]$ (width $6.7\times10^{-10}$); $q{=}2,n{=}8$ their $[0.18324295,0.5]$ against $[0.265375771,0.265382499]$; $q{=}4,n{=}4$ their $[0.63826289,1.0]$ against $[0.707114444,0.707114895]$. Verdict on every row `CERT_LOWER_BEATS_LBplus` **and** `CERT_UPPER_BEATS_UB`, decided by Arb-endpoint inequalities with margins $0.04$–$0.08$ and $0.08$–$0.29$ bits/symbol — four to six orders above the widths. `MACHINE-VERIFIED`. **Owner-verified at 300 bits:** $h_2(1/2)=1$ exactly, their closed forms reproduce their printed $\mathrm{UB}$/$\mathrm{LB}_1$ columns for $q=2,3,4$, and three spot-checked rows match the reported gains to the printed digit with the certified interval strictly inside $[\mathrm{LB}^+,\mathrm{UB}]$. Their own points: 15 further rows, every certified interval strictly inside their sandwich, **zero exclusions**. Extension rows with no printed analogue (MD show these $\delta$ only as Fig. 2 curves): **30** LO-CVB rows at $\delta\in\{1/20,1/2,4/5\}$, ball radius $\le2.84\times10^{-120}$, plus **30** Pinto–Ribeiro $C_{n,k}$ rows at $n=6..9$ — first certified values, explicitly *not* comparisons. **Rule 7 on the extension rows:** only the top 24 $\Lambda$ subsets per row are Arb-certified, so each is a *valid* converse bound but **not** a proven global optimum over $\Lambda$; $E(m,w)$ is `CITED-DEPENDENCY` except $w\in\{0,1,2,m\}$, re-derived and matching. **Also first rigorous interval enclosure of the published deletion-channel constants.** **And tighter certified lower bounds than the published BA column:** TNB's $C_{q,n}$ column sits below our certified enclosure in **9 of 15** rows by more than one 3-dp ulp ($+1.3\times10^{-3}$ to $+1.1\times10^{-2}$), from the primal certificate alone. **No theorem is contradicted** — their ordering $\mathrm{LB}_1\le\mathrm{LB}^+\le C_{q,n}\le\mathrm{UB}$ holds and our intervals lie strictly inside their own sandwich in all 15 rows, confirming it. The excess is monotone increasing in **both** $n$ and $d$ (owner-verified), the signature of an unconverged Blahut–Arimoto whose primal climbs to capacity **from below** — so their $C$ column is an unconverged *estimate*, not a claim to printed precision (`INFERENCE`). Certified-Arb Blahut–Arimoto: mpmath 150-dps locating step $\to$ exact rational snap on the simplex $\to$ outward-rounded Arb evaluation of **both** the primal $I(p^*)$ and the Csiszár–Tusnády dual $\max_x \mathrm{KL}(W(\cdot\mid x)\Vert D_{p^*})$, so the interval $[\text{primal},\text{dual}]$ contains the true capacity **regardless of BA convergence** — convergence controls only tightness. 12 Fertonani–Duman Table II entries enclosed exactly, plus a new certified $U(28,0.68)=0.13168$. Gate C: **all 36 Morozov–Duman Table III LO-CVB rows independently recomputed and certified**, zero-width intervals at 400 bits, agreeing with the published values to within $7.1\times10^{-7}$ absolute under the paper's evident **round-up** printing convention — the conservative and correct convention for an upper bound. `MACHINE-VERIFIED` for our values; the convention identification is `INFERENCE`. **No published claim is contradicted.** Two negative findings recorded: the proposed type-channel symmetry reduction was **falsified by the target itself** ($S_n$-equivariance of the channel law fails), and the Fertonani–Duman apparent discrepancy dissolved into their stated round-up rule on an unconverged dual. Solver audit clean — zero LP/QP/MILP in the certified path — and a **live invalid-certificate branch was caught before it bit**: at $q=2,n=10,d=1/20$ the empty-output mass $\sim10^{-13}$ sits *below* the snap resolution $2^{-40}=9.09\times10^{-13}$, so $D'(y)$ would have snapped to zero where $W(y\mid x)>0$, making $\mathrm{KL}=+\infty$; silently skipping that term would have produced a too-small "upper bound". `cert_dual` returns $+\infty$ (no claim) instead, with full support guaranteed two ways. |


**Current-status override for the table (2026-08-31).** The historical
`omega/` gate-B rung labels in the row above are **SUSPENDED** after the shared
outward-upper defect; only the separately replayed five-point candidate ray is
restored, with all five moves strictly worsening and an **OPEN** verdict. The
two historical rungs and box enclosures await corrected replay. For `kg/`, the
direct-D.4 restart has now re-certified `[1,1.3]`, `[1.30,1.45]`, and
`[1.45,1.75]`; `[1.75,3.5]` is partial with 25/36 tiles,
`[3.5,4.083]` and the refined `[4.083,6]` remain open at their panel caps, and
the widened interior refinement is running. Later correction sections below
are the controlling evidence.

## `oct-rank/` Routes A/F closed to their honest boundary, 2026-08-31

**Route A is dead as a route to 14.** For every linearly independent real
octonion triple $(u,v,w)$,

$$[L_{\bar u}L_v,L_{\bar u}L_w]^2
  =-4N(u)\det\operatorname{Gram}(u,v,w)\,I.$$

The Gram determinant is positive, so the commutator has rank eight. The proof
uses the coefficientwise polar identity
$L_{\bar x}L_y+L_{\bar y}L_x=2\langle x,y\rangle I$, orthogonally decomposes
$v,w$ from $u$, and squares the normalized commutator. Exact rational
coefficient controls cover all basis pairs, 56 independent triples, 28
dependent plants, and the quaternion calibration; the universal deduction is
**HUMAN-AUDITED**, not overlabelled as a machine proof.

Owner first-hand reading of Landsberg §6.1, Theorem 6.1.1 confirms the actual
Strassen-form inequality

$$\operatorname{rank}[T_{\alpha,\alpha_1},T_{\alpha,\alpha_2}]
  \le 2(r-b).$$

Thus $b=8$ and commutator rank eight give only $r\ge12$, **not 16**. The
quaternion control gives $r\ge6$, consistent with true rank seven; the naked
$b+\operatorname{rank}$ reading would falsely give eight and is retracted.
Landsberg and Koiran were read first-hand; attribution to Strassen's 1983
primary remains `CITED-DEPENDENCY / PRIMARY UNESTABLISHED`.

**Route F remains exactly $\{13,14\}$.** The conjugated $(1,i,j)$ tensor is
the shared-first-factor block duplication
$\tau\boxtimes(e_0\otimes I_2)$: two identical $4\times4$ quaternion slice
blocks. The frozen rank-seven $\tau$ witness supplies upper 14, while the
existing substitution/pencil chain supplies lower 13. Five numerical CP
starts and three commuting-extension starts found no witness; their best
relative residuals, $1.33\times10^{-3}$ and $1.00\times10^{-4}$, are
`COMPUTATIONAL-EVIDENCE`, never infeasibility. Ordinary direct-sum additivity
does not apply after identifying the first factors, and the available
matrix-pencil multiplicativity theorem has first dimension two over
$\mathbb C$, not three over $\mathbb R$.

Frozen campaign:
`oct-rank/campaigns/2026-08-31T08:02:18Z_routeAF/`; owner checksums passed all
23 listed payloads. Two numerical searches violated `kg/`'s exclusive CPU
slot (105.30 s and 32.67 s wall), so timing/performance interpretation is
forbidden. The agent also deleted two untracked drafts without required user
confirmation; both were restored byte-exactly from the edit transcript and
are explicitly `UNEXECUTED / NOT EVIDENCE`. No lasting file loss was found.
The full octonion multiplication window remains
$18\le R_{\mathbb R}(T_{\mathbb O})\le25$.

## `omega/` Gate-C slope replay repaired — V1 passes; full-box sign remains OPEN, 2026-08-31

**The old $+0.9008715$ residual is retracted as nonreproducible.** The
advertised frozen source, SHA-256
`4504365a7917f80362424a7252203cd770b7c51392061e517f8c10e62fc946fa`,
raises `NameError: s_const` before its first aggregate checkpoint; no frozen
producer for raw $3.272425321778391$ exists. The corrected replay matches all
seven frozen $R$ branches, $R=2.8170035674609757$,
$M=2.0942543887102634$, and raw
$\Omega=2.3715538358350807$, residual
$4.440892098500626\times10^{-16}$: **V1 PASS**.

The inherited penalty hypothesis was false. Zeroing only the stated Lemma-1
contribution shifts raw $\Omega$ by
$2.029648699330977\times10^{-6}$, not $0.9008715$; there is no $3\times10^5$
amplification. The old number is `REPORTED` history only, and its implied
$R=0.9303495043602614$ is `INFERENCE`.

Four certificate-path defects were fixed before decision: an extra $/r$
shrunk the effective LP radius from $10^{-7}$ to $10^{-14}$; 21
variable-bound marginals were incorrectly treated as a 45-coordinate dual;
the dense SVD basis was loaded instead of the exact integer JSON kernel basis;
and the arbitrary-$y$ full-box instrument was absent. Exact Arb endpoints now
replace inward binary64 subtraction, and every interval-overlapping minimum
branch is hulled.

**Honest result:** the exact-basis midpoint LP gives a numerical candidate
$-1.577332901516852\times10^{-7}$ (`COMPUTATIONAL-EVIDENCE`). The rigorous
arbitrary-$y$ enclosure over the complete named 45-coordinate box is
$[-5.841014555619999,\,5.841014555619999]\times10^{-6}$, width
$1.168202911124\times10^{-5}$ (`MACHINE-VERIFIED`), so it straddles zero.
Verdict: **FAILURE TO CERTIFY / OPEN**—neither an improving witness nor local
optimality.

Campaign:
`omega/campaigns/2026-08-31T07:57:00Z_OmegaGateCDiag_4859327/`; owner
`sha256sum -c` accepted all 17 payloads. Scope is only the region-0
`glob dist[0]` nullspace at radius $10^{-7}$ in the transcribed
$q=5$, max-level-3 program. It bears on neither the published
$2.37155181$ endpoint nor the $2.371177$ record. A held successor will test
the single numerical direction with exact dyadic feasibility and direct Arb
point comparisons.

### Correction — the registered ray preserves a non-unit center mass

The held successor described above did run its five fixed points and remained
**OPEN**, but “exact dyadic feasibility” was too strong.  Its exact guard
proves `A*d=0`, `sum(d)=0`, positivity, and preservation of the released
center's mass.  It does not prove that mass is one.  Exact summation of the
registered coordinates in
`omega/campaigns/2026-08-31T09:18:39Z_OmegaGateCCandidateWitness_e9c6414/protocol_static_corrected.json`
gives

$$
\sum_i p_i
=\frac{37778931862949057407573}{37778931862957161709568}
=1-\frac{8104301995}{37778931862957161709568}
<1.
$$

All five old interval differences still contain zero, so their **OPEN**
verdict is unchanged and no false witness was promoted.  The points form a
fixed-mass affine ray at the released binary64 center, not a feasible ray for
the registered sum-one constraint.  Any later strict sign on those same
points can establish only descent of the explicitly repaired
generalized-entropy/absorbed-defect objective; it is not a feasible
construction, a constrained-local-optimality counterexample, or an exponent
improvement.

### Correction — shared Omega interval enclosures suspended

Frozen-source audit found that `omega/src/interval_core.py::from_two` used
`b.lower()` for its upper endpoint instead of `b.upper()`.  A nonzero-radius
Arb upper expression can therefore be rounded inward before later buffering.
Pending a corrected replay, all Omega results whose rigorous label depends on
that core are **SUSPENDED**: the two-rung interval endpoints, certified-box
endpoints, full-box sign enclosure, and five candidate-ray interval
differences.  Stored numerical outputs remain evidence of what the old code
computed, not accepted formal enclosures.  Float reproductions and exact
integer/rational rank, mass, and feasibility statements are unaffected.

The first generalized-entropy repair run (v9) aborted at the base point on a
transposed feature matrix before producing any certificate.  Its static v10
successor fixed that interface but was rejected before production: its
binary64-only exact-point converter cannot represent the registered moved
dyadics, and it retained the defective interval dependency.  Neither campaign
established a sign.  Restoration requires an outward-upper fix, exact Arb
dyadic extraction without binary64 rounding, a fully frozen code/input
perimeter, and replay.

### Partial restoration — candidate ray replayed under repaired core

Campaign
`omega/campaigns/2026-08-31T14:29:40Z_OmegaGateCEntropySolver_v12_d2d3d4d5/`
replays the base and all five registered divisors with the outward-upper fix,
exact finite-Arb dyadic extraction, hard-hashed local inputs, and an
independently checked KKT gate. The intervening v11 production is preserved as
**ABORTED** at a stale old-objective raw anchor; its diagnostics motivated a
least-squares residual solve after the L-BFGS-B warm start.

All 1,152 entropy certificates passed. The maximum outward KKT upper bound is
`1.4456822252392066e-10 < 1/2^32`. The repaired full-objective differences
(moved minus base), in divisor order 16, 8, 4, 2, 1, are strictly positive:
`[6.475289389514407e-9,6.475289390521749e-9]`,
`[1.2954036629144681e-8,1.2954036630152025e-8]`,
`[2.5911541940662355e-8,2.5911541941669702e-8]`,
`[5.182659614230212e-8,5.182659614330947e-8]`, and
`[1.0365687869954472e-7,1.0365687870055208e-7]`
(`MACHINE-VERIFIED`). Thus none of these five moves is a descent direction for
the repaired fixed-m objective.

Verdict remains **OPEN**: the fixed-mass points are not sum-one feasible, and
one worsening ray proves neither local optimality nor an exponent bound.
Result SHA-256 is
`0403e652099a217d421c9bb1aa7fd56f07fb6daff72ddf14273d6a6dff063122`.
This restores only the candidate-ray statements. The historical two-rung and
box enclosures remain **SUSPENDED** pending their separate corrected replay.

### Corrected-core two-rung replay measured, but formal endpoints remain suspended

The separate two-rung campaign ran all three frozen child computations. Its
aggregate launcher then aborted with `KeyError: 'omega_cert_upper'` because the
rung-1 schema contains distinct raw and absorbed-slack fields, not that key. A
hash-gated additive recovery compared both without silently selecting a
replacement: rung 1 raw/absorbed are `2.371340083602922` /
`2.3713411672715115`, while rung 2 with Lemma is
`2.3715538358544617`; both orderings hold and all are below `2.371866`.
Recovered artifact SHA-256:
`8160e846c4b64a167624dad8460e83bc29c4e69ac16115a7015da208a838adda`.

This is **COMPUTATIONAL-EVIDENCE**, not restoration. The original launcher
abort is preserved, its missing field is unresolved, and the formal two-rung
endpoints stay **SUSPENDED** pending an end-to-end outward aggregate free of
decision-relevant float extraction and a proof for the feasibility-absorption
term.

## `delcap/` orbit-invariance correction — 19 prior $q=3$ rows suspended, 2026-08-31

**Repository-certificate defect found:** the frozen largest-remainder
per-output-word snap can split an output orbit, while the old dual certificate
then maximized only one input representative per orbit. At $(q,n,d)=(3,6,1/5)$,
flat output orbit 51 (size 12) receives exact integer counts $578{,}028$ and
$578{,}029$ at denominator $2^{30}$. The reference law $D$ is therefore not
output-orbit invariant, so the representative-only maximum is inadmissible. A separate
$(3,3,1/2)$ plant proves the failure can be permissive: the compressed value is
strictly smaller than the true full-alphabet dual.

Mandatory audit stopped at the first failure: $(3,6,1/2)$ passed exact orbit
equality, $3{,}279$ generator checks, all $729$ inputs, and full-versus-
representative agreement; $(3,6,1/5)$ failed; the other 18 rows did not run.
Consequently only the first row currently retains certification. The other
19 intervals in the table and summary are preserved but **RETRACTED /
SUSPENDED** until replay from exact orbit-coordinate masses. No external
published claim is affected.

The clean fix is not a filter: snap total dual-reference mass in
**output-orbit** coordinates, set $D_y=D(O)/|O|$ exactly, prove positivity,
normalization, and invariance, then compare the representative maximum with
every full-alphabet input. If the primal is symmetry-compressed, its input
distribution and full-alphabet value need separate equality checks. The replacement
campaign must re-certify all prior 20 rows before extending to $q=4$.

### Correction replay complete — all 20 $q=3$ rows restored, 2026-08-31

**The suspension above is resolved by a replacement certificate, not by
reinstating the defective endpoints.** Campaign
`delcap/campaigns/2026-08-31T09:13:25Z_7dc5babe-5e1b-44fa-a9e5-03d2f10183b6_q3-invariance-correction/`
snaps exact total mass in input- and output-orbit coordinates, expands each
mass uniformly within its orbit, and directly replays the primal and dual on
every positive-mass input word. All 20 replacement rows at
$n\in\{6,7,8,9,10\}$ and
$d\in\{1/2,1/5,1/10,1/20\}$ certify both
`CERT_LOWER_BEATS_LBplus` and `CERT_UPPER_BEATS_UB`.

The old defect was widespread: 16/20 accepted `ba_word` vectors split an
output orbit; only four were invariant. Every corrected row instead selects
`ba_total_orbit_mass`. The replay performed 1,058,508 direct full-input KL
evaluations, 5,821,704 exact generator checks, and 237,744,816
positive-$W$/positive-$D$ term checks, with zero positive-mass words skipped.
Owner verification accepted all 17 checksummed payloads. Conservative
400-bit-Arb interval widths range from
$1.2226475645090283\times10^{-8}$ to
$2.220013228610702\times10^{-5}$; binary64 fields remain report-only.

**Scope:** this restores exactly the 20 finite-$n$, $q=3$ statements above.
It does not restore the old endpoints, prove an asymptotic capacity claim, or
certify any $q=4$ row. The separately frozen $q=4$, $n=5,\ldots,10$ extension
remains **NOT RUN / FAILURE TO CERTIFY** pending its released compute slot.

### $q=4$ startup anchor released and passed

The frozen extension has now completed its mandatory
`q=4,n=5,d=1/2` anchor. Its 400-bit outward per-symbol interval is
`[0.6664806108007938,0.6664806978587714]`, width upper
`8.705797740570955e-8`; all 3,072 input and 4,095-per-candidate output
generator checks passed, and both dual candidates were evaluated on all 1,024
inputs. Verdict:
`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB` (**MACHINE-VERIFIED**).
Anchor SHA-256:
`c47aeb35250eade511e3cfe8e65701e9302545b18d46d995952cdd5f073b794b`.
The anchor writes no row artifact, so the fixed 24-row `n=5..10` grid remains
**COMPUTE PENDING**; no other $q=4$ row or asymptotic claim is certified.

## `rs-pe3d/` H-GATE — support-window mechanism proved, 2026-08-31 (owner-derived, agent-audited)

**Outcome:** for the implemented lifted-line sum of diagonal generalized
Reed–Solomon codes,

$$d_{\rm Ham}(V)=\min_i(s_i-t_i+1).$$

Hypotheses: distinct evaluation points, nonzero diagonal entries, and
$1\le t_i\le s_i$; a zero component code is handled separately. Precisely:
every support $S$ smaller than this distance has
$V\cap\mathbb F_q^S=\{0\}$, and a minimum GRS word lifted on one line supplies
**some** support at equality. Hence the **union** of all supports of weight at
most $w$ is empty iff $w<d_{\rm Ham}(V)$. No individual-support “iff” is
claimed above the threshold.

At $t=(1,1,1)$ this proves the observed H-GATE exactly: the weight-$\le3$
window is empty iff $\min_i s_i\ge4$. It also replaces the sampled predictor
outside that slice: the general predicate is $\min_i(s_i-t_i+1)$, not
$\min(s)$. Pairwise coprimality, repeated/equal orders, $s_i=1$, and any
nonzero $\Lambda$ do not change the distance theorem.

**Proof:** put $Q_i=A_i/C_i$. Right exactness identifies
$V=\ker(\bigotimes_i\pi_i)$. MDS distance makes every set of fewer than
$d(C_i)$ coordinate cosets independent. Quotient-dual functionals, assigned
coordinatewise, isolate any chosen point of a smaller support; their tensor
pullback annihilates $V$, forcing every coefficient to zero. A
degree-$(t_i-1)$ polynomial with $t_i-1$ evaluation roots gives the sharp
line word. Evidence strength: **HUMAN-AUDITED, machine-anchored**, not a formal
proof.

Frozen campaign:
`rs-pe3d/campaigns/2026-08-31T08-26-54Z_hgateMechanism/`; four checksums
owner-accepted. Exact controls report the dimensions $13/22/37$, the quotient
kernel dimension $37=64-27$, complete below-distance emptiness on
$16/30/1{,}830$ supports, zero hits on all $82{,}160$ weight-3 supports at
$(41,(4,4,5))$, and both $s_i=1$ corners nonempty on every tested support.

**Citation lock:** the frozen prose contains seven errors caught in owner
review—most importantly a false per-support biconditional, two incorrect GRS
sentences, a wrong kernel-dimension formula, and a false preregistration
claim. It is authoritative only through the appended
[`rs-pe3d/README.md`](rs-pe3d/README.md) owner corrections. The campaign had no
separate pre-compute commit; its timing is session-chronology/`REPORTED`.
One 52.9-second control also overlapped `kg/`'s heavy slot, so wall time is not
performance evidence. None of this alters the corrected linear-algebra proof.

**Untouched:** Conjecture 4.2, every global $\rho_{\rm inst}$, non-GRS
constructions, and optimization in the first nonempty weight-$\ge4$ window.
For the old empty window, $\rho^{\rm window}$ is **undefined**, not zero.

## `mm3/` gate C record attack — certified negative, 2026-08-30 (agent `Mm3RecordAttack`, owner-verified)

**Outcome: no 54-addition scheme was found, and the certified minimum is exactly 55 over the 288 monomial orientations actually DECIDED — checker-verified, and valid under any tensor automorphism.** The wider claim over all of $S$ **does not stand**: the exclusion of the other 6912 was computed under a map that is an automorphism only on the orthogonal subgroup, and under the corrected action `sun56` admits **576** valid diagonal orientations, **528 of them non-monomial and never decided**. See the re-opened audit below. **No published claim is affected** — the record 55 was never beaten here, and the gap is in this repository's coverage, not in arXiv:2607.28676. The published record 55
(arXiv:2607.28676) is **not beaten**, nothing was escalated, and no published claim is
contradicted. `pre_statement.md` was committed (`fb731d4`) before any computation, with $S$, the
budget, the wall-clock cap and the adjudication rule fixed in advance.

**The named set $S$, exactly:** $\{\texttt{paper55},\texttt{sun56}\}\times\{6960$ ternary
$3\times3$ matrices with $|\det|=1\}\times$ diagonal sandwich $(X,Y,Z)=(G,G,G)\times
\sigma\text{-orbit }3$ = **41,760 triples, every one decided exactly, no sampling.**

**Owner-verified by independent brute force over all $3^9=19683$ ternary matrices** — every count
reproduced: 11808 invertible, **6960** with $|\det|=1$, **48** monomial, **6912** non-monomial;
$41760=2\cdot6960\cdot3$; 288 valid orientations $=2\cdot48\cdot3$; alphabet factor
$6960/48=\mathbf{145}$; unswept off-diagonal $6960^3\cdot3=\mathbf{1{,}011{,}460{,}608{,}000}$.

**⚠⚠ AUDIT RE-OPENED AND THE CLAIM RE-NARROWED, 2026-08-30 — THE OWNER'S EARLIER CLOSURE WAS
PREMATURE AND IS RETRACTED HERE (rule 5).** The closure below rested on a measurement comparing the
frozen map against the agent's **first** derivation of the honest action. That derivation was
subsequently **disproved by the agent's own control** (failed Brent $729$ on $10/10$ random
non-monomial triples). **Two invalid maps agreeing is not corroboration** — precisely the failure
the owner had flagged in writing hours earlier about the campaign's "two independent paths", and
then walked into.

**The corrected automorphism family, owner-verified symbolically and numerically.** The matmul
trilinear form is $F(U,V,W)=\sum_{i,k,j}U_{ik}V_{kj}W_{ij}=\langle UV,W\rangle=\mathrm{tr}((UV)^{\mathsf T}W)$
— **not** $\mathrm{tr}(UVW)$, which is what the first derivation assumed. The correct family is
$$u'=P^{-\mathsf T}UQ^{-1},\qquad v'=QVR^{-1},\qquad w'=PWR^{\mathsf T}$$
invariant for **all** invertible $P,Q,R$. **Owner checks: symbolic identity confirmed; 200/200
random integer $(U,V,W,P,Q,R)$ invariant for the corrected family and 0/200 for the naive one.**

**Consequence — a concrete counterexample population, not a hypothetical.** At the diagonal the
corrected action is $u'=G^{-\mathsf T}UG^{-1}$, $v'=GVG^{-1}$, $w'=GWG^{\mathsf T}$; the $v$-side
differs from the frozen $GVG^{\mathsf T}$ off the orthogonal subgroup. Re-censused:

| decomposition | frozen map | corrected map |
|---|---|---|
| `paper55` | 6912 / 0 / **48** (all monomial) | 6912 / 0 / **48** — agree |
| `sun56` | 6912 / 0 / **48** | 6384 / 0 / **576**, of which **528 are NON-monomial** |

So the record attack's *"the 48 survivors are exactly the monomials"* is **map-relative**: true for
the frozen map and for `paper55`, **false for `sun56` under the true automorphism**.

**What stands and what does not.**
- **The 288 monomial decisions STAND unconditionally.** Monomial orientations are valid under *any*
  automorphism, and those were Brent-verified and floor-certified.
- **The exclusion of the 6912 does NOT stand.** There is now a concrete address for the gap: **528
  `sun56` non-monomial diagonal orientations that are ternary AND Brent-valid** — genuine
  decompositions of the same tensor, inside the ternary alphabet, **never decided**. Their certified
  lower bounds are **UNKNOWN**, and whether any has total $\le54$ is exactly the open question.
- **Second-order finding:** the record-attack manifest's structural note that *"each $G$ is a tensor
  automorphism with integral inverse"* is **wrong off the orthogonals** — 8/8 random non-monomial
  diagonals fail all 729 Brent identities under the frozen map. They were never Brent-tested there
  because the ternarity gate excluded them first.

**◆ DISPUTE RESOLVED — AND THE ANSWER IS THAT "TERNARY" IS FAMILY-RELATIVE, 2026-08-30.** Neither
side had a bug. Owner replayed the agent's witness (`sun56`, $G=[[-1,-1,-1],[0,1,1],[0,0,1]]$,
$\det=-1$, non-monomial), **anchoring the Brent implementation on the UNMAPPED decomposition first
at 0 failures of 729**, and got **Brent 729/729 PASS but max $|x|=4$** (10/23, 10/23, 11/23 rows out
of range). The agent got max $|x|=1$ on all three. Both passed Brent.

**The decisive test.** The owner enumerated all $4^6=4096$ letter assignments over
$\{G,G^{\mathsf T},G^{-1},G^{-\mathsf T}\}$ and asked which give **all-ternary AND Brent 729/729** on the
owner's arrays. **Exactly one survives:**
$$U\mapsto G^{-1}UG^{-\mathsf T},\qquad V\mapsto G^{\mathsf T}VG^{-\mathsf T},\qquad W\mapsto G^{\mathsf T}WG$$
which is precisely the **transpose-dual** of the family the owner applied. Since
$(LMR)^{\mathsf T}=R^{\mathsf T}M^{\mathsf T}L^{\mathsf T}$, a transpose-dual map is the same map applied to
**per-block transposed** arrays. **Diagnosis: the two sides store the $3\times3$ blocks with opposite
row/column orientation.** Each applies a genuine automorphism and each gets Brent 729/729 within its
own convention.

**The load-bearing consequence.** Brent-validity is convention-invariant; **ternarity is not.** Of
4096 variants only one is all-ternary, so a block transpose moves a survivor from inside the
alphabet to outside it. Therefore *"the set of ternary orientations"* is well-defined **only relative
to a named automorphism family and a named block convention** — and any count entering this index
must name both or it is not a well-posed claim. The owner's position: the **frozen `gatec_sweep`
reshape is the reference**, because the ternary alphabet, the $d$-counters, the floor CNFs and the
entire published-comparison ladder are all expressed in it; a census in the other orientation is
internally consistent but **not comparable** to the frozen landscape or to the paper's 55.

**The agent's DFS result, accepted as a conditional and reassuring in every branch:** certified
totals for all **576** valid `sun56` diagonal orientations, **minimum 55, ZERO $\le54$**, all decided
in 109 s, worst row `gi=1068` at $13{+}12{+}16{+}14=55$, coinciding with the known `sun56`
achievable-floor structure. So the population the divergence exposed **does not contain a 54**. It is
recorded as certified $\ge55$ over the 576 **under the named family and convention**, pending the
reference decision — not as an unconditional statement.

**⚠ SUPERSEDED DISPUTE ENTRY, retained per rule 5:** Owner replayed the agent's supplied witness independently: `sun56`,
$G=[[-1,-1,-1],[0,1,1],[0,0,1]]$ ($\det=-1$, non-monomial, $G^{-1}$ also ternary), corrected map
applied blockwise. **Brent implementation anchored first on the UNMAPPED decomposition: 0 failures
of 729.** Then on the mapped witness:

- **Brent 729/729 PASS, zero failures — CONFIRMED.** The witness is a genuine rank-23 decomposition,
  and since a wrong map would fail Brent, the pass corroborates that the owner's replay implements
  the agent's map faithfully.
- **"Ternary throughout" NOT REPRODUCED, and not marginally.** $U'$: 10 of 23 rows carry an entry
  with $|x|>1$, max $|x|=4$. $V'$: 10 of 23, max 4. $W'$: 11 of 23, max 4. Example
  $V'$ row 0 $=[-2,-4,0,1,2,0,0,0,0]$.

**The two outcomes are opposite and the difference is one predicate.** If the 528 are Brent-valid
but **non-ternary**, they lie outside the ternary-alphabet model that scopes this target, **the
exclusion of the 6912 stands after all**, and the record attack's coverage is fine. If they **are**
ternary and the owner's replay is wrong, the re-narrowing above stands. Resolution requested as a
single number (max $|x|$ per factor on this witness); the floor-DFS over the 576 is held until the
selection predicate is settled, to avoid certifying a mis-selected set.

**The index deliberately does NOT swing a third time on unresolved evidence.** This claim has moved
twice today — narrowed, restored, re-narrowed — and a third move before the predicate is settled
would make the record less trustworthy than the stated uncertainty. What is **not** in dispute: the
**288 monomial decisions stand unconditionally**, and **no published claim is affected in any
branch** — the record 55 was never beaten here.

**Assigned and in progress:** certified totals for all 576 valid `sun56` diagonal orientations
(minutes), then the off-diagonal census under the corrected map. Any total $\le54$ escalates.

**✗ SUPERSEDED CLOSURE, retained per rule 5 — the owner's premature "restored" entry:** The exclusion was challenged, the owner refused an argument and required a measurement, and
the measurement came back confirming it. **Result (agent `Mm3OffDiag`, pass A):** the diagonal
census re-run over all 6960 under **both** maps in the same enumeration order gives
**frozen $6912/0/48$ and honest $6912/0/48$, survivor sets identical bit-for-bit (48 members, all
monomial), Brent-fail cell $0$ in all four runs, for both decompositions.** So the certified
minimum **55 over the full named set $S$ of 41,760 triples stands**, and the exclusion of the 6912
is sound. It is now verified under two distinct maps, one of which is provably the correct action —
evidence the original campaign did not have. The audit trail below is retained per rule 5.

**⚠ WHAT WAS AUDITED, AND WHY IT WAS NOT A FALSE ALARM (2026-08-30).**
The frozen `gatec_sweep.inv_sp()` implements the sandwich inverse as the **TRANSPOSE**. For an
integer matrix $G^{\mathsf T}=G^{-1}$ holds exactly when $G$ is orthogonal, i.e. exactly for the
signed permutations — **owner-verified: the two sets coincide, 48 of 6960, bit-for-bit**. So:

- **The certified $\ge55$ STANDS, but only over the 288 decided orientations.** Those are the
  monomial ones, where transpose *is* the inverse exactly, so they were decided under the correct
  action. Nothing about them moves.
- **What is under audit is the EXCLUSION of the 6912**, which was computed under
  $U\mapsto GUG^{\mathsf T}$ rather than $U\mapsto GUG^{-1}$ — a genuinely different map for every
  one of them. **Owner measurement: 4608 non-monomial $G$ have a ternary $G^{-1}$**, so an
  inverse-ternarity argument does **not** exclude them; 4656 of 6960 have both $G$ and $G^{-1}$
  ternary, of which only 48 are monomial. Example: $G=[[-1,-1,-1],[-1,-1,0],[-1,0,-1]]$ has
  $G^{-1}=[[1,-1,-1],[-1,0,1],[-1,1,0]]$, both ternary, $G$ non-monomial.
- **ROOT CAUSE SHARPENED, 2026-08-30 (agent `Mm3OffDiag`, owner-re-derived independently).** The
  owner's first diagnosis — *"transpose used where inverse was required"* — was the **symptom**. The
  cause: reading `gatec_sweep.py:97-110`, the implemented factor map is
  $u'=X_i^{\mathsf T}UY^{\mathsf T}$. Matching it against the tensor-preserving family
  $u'=P^{\mathsf T}UQ^{-\mathsf T},\ v'=Q^{\mathsf T}VR^{-\mathsf T},\ w'=R^{\mathsf T}WP^{-\mathsf T}$
  — which fixes $\mathrm{tr}(ABC)$ for **any** invertible $P,Q,R$ — the $U$-side forces $Q=G^{-1}$
  while the $V$-side forces $Q=G^{\mathsf T}$. **Both hold only if $G^{\mathsf T}=G^{-1}$, i.e. $G$
  orthogonal, i.e. monomial.** So the frozen formula **is a genuine tensor automorphism exactly on
  the 48 and on nothing else** — not a miscomputed inverse, but a map valid only on the orthogonal
  subgroup, applied outside it with the domain assumption invisible at the call site.
- **The frozen count is probably sound-by-accident, and is being measured rather than argued.** The
  *intended* map $u'=G^{-\mathsf T}UG^{\mathsf T}$ equals the honest action at
  $(P,Q,R)=(G^{-1},G^{-1},G^{-1})$ — a relabeling by the bijection $G\mapsto G^{-1}$ — and inversion
  preserves monomiality, so the survivor set should be the 48 either way and $6912/0/48$ should
  stand as a count. **That is currently an argument, so per rule 14 the index stays narrowed until
  it is a measurement:** the assigned control reports honest-action versus frozen-formula survivor
  vectors elementwise over all 6960, with the falsifiable prediction that the frozen formula
  **fails** Brent on random non-monomial $(X,Y,Z)$ while the honest action passes 729/729.
- **No frozen code is reusable off-diagonal.** The frozen $W$-side pattern
  $(X^{-\mathsf T},Z^{\mathsf T})$ has $X$ and $Z$ roles swapped against the honest cyclic structure
  $R^{\mathsf T}WP^{-\mathsf T}$, so for $X\ne Z$ it is not a conjugation at all.
- **Owner-verified census design (every figure re-derived).** Right cosets of the 6 permutation
  matrices over the 6960: **exactly 1160, every orbit size exactly 6, action free**
  ($6960=1160\times6$). Permutations factor out of $u'=P_p^{\mathsf T}(a_p^{\mathsf T}Ua_q^{-\mathsf T})P_q$,
  so ternarity, $d$, floors and Brent depend on **data-pairs only** — reducing the sweep to a
  triangle count over 1160 nodes with three $1160^2=1{,}345{,}600$-entry pair tables, times exactly
  $648=216\times3$ schemes per surviving triple. $1160^3\times648=1{,}011{,}460{,}608{,}000$,
  **identical** to $6960^3\times3$, so the factorization loses nothing. **Built-in control:** the
  monomial floor is $2^3=8$ data nodes, $8^3=512$ triples, $512\times648=\mathbf{331{,}776}$ — which
  equals $48^3\times3$ exactly, forcing the machinery to reproduce a figure derived a completely
  different way before it may report anything new.
- **Owner pinned the call chain, and the precondition is stated in the code it was violated
  against.** `gatec_sweep.py:115-122` — `inv_sp(M)`, docstring *"Inverse of a 3x3 signed
  permutation matrix"*, body `out[j][i] = M[i][j]` with the comment *"transpose of monomial gives
  inverse"*. The function is **correct and its domain is documented.**
  `gatec_sweep.py:93` — `sandwich()` hard-codes `Xi, Yi = inv_sp(X), inv_sp(Y)`, making the
  precondition **invisible to every caller**. Inside gate C's own sweep the domain is respected:
  `gatec_sweep.py:80-81` builds only `SP` with `assert len(SP) == 48`. The violation is in the
  record attack: `census_fmpz_path.py:52` asserts `len(unimod) == 6960` and lines 66-67 / 82-83
  call `gs.sandwich(U, V, W, G, G, G)` across all of them. **A precondition violation at a call
  site, not a wrong formula** — which is exactly why owner count-checking had no power against it.
- **A second consequence, and the more damaging one: the "two independent paths" claim is in
  question.** `census_fmpz_path.py:4` states it *"Uses the FROZEN sandwich implementation"*, so it
  inherits `inv_sp` verbatim; the kernel path is *"my own re-implementation of sandwich"*. If that
  re-implementation also used the transpose — which exact agreement on 6912/0/48 suggests — the two
  paths differ in **arithmetic backend** (pure int vs `fmpz`) while sharing **one wrong map**. Two
  implementations agreeing on the same defective semantics is not corroboration. Assigned for
  explicit determination, with the possibility recorded that the transpose substitution may not
  change the alphabet verdicts at all, which would largely rescue the census conclusion.
- **Also flagged for the record:** that campaign's `pre_statement.md:43-55` justifies the diagonal
  restriction by asserting the diagonal *"contains every non-monomial sandwich that keeps the
  factor blocks ternary with the SAME rank-23 structure"* — an unproved claim about the very
  off-diagonal region now being swept. Not inherited; falsifiable by the running census.
- **The narrowing this triggered is now LIFTED** (measurement above). At the time it was correct
  to impose: had the two maps disagreed, $S$ would have been smaller than claimed and $\le54$
  could have hidden in the wrongly-excluded region. If any of the 4608 preserves the alphabet on
  the actual decomposition, the swept set was smaller than claimed and $\le54$ could hide in the
  region wrongly excluded. The re-run is assigned; verdict pending. **No published claim is
  affected either way** — the record 55 was never beaten here.

**Per-decomposition census (MACHINE-VERIFIED under the frozen transpose map, two independent
arithmetic backends agreeing — pure-int kernel and frozen-code `fmpz`; the map itself is the
subject of the audit above):** exactly **6912 of 6960** non-monomial $G$ **exit the
$\{-1,0,1\}$ alphabet** under the diagonal action (entries $\pm2$ appear); **0** kept ternarity
but failed Brent; **48 survive, and they are exactly the monomials — under the frozen map.** **Qualification added 2026-08-30:** this survivor identification is **map-relative**. It holds for the frozen map and for `paper55` under the corrected automorphism, but is **false for `sun56` under the corrected family**, where 576 survive and 528 are non-monomial. The *count* $6912/0/48$ is confirmed correct under all three maps tested; it is the *“exactly the monomials”* reading that does not generalise. See the off-diagonal programme above.

**The two reductions are split honestly, and this correction matters.** The $145\times$ collapse
is the **ternary alphabet**, not the monomial-transfer theorem: it is a diagonal-slice,
instance-exact census and **not a theorem**. The theorem covers **exactly the 48 monomial
orientations**; the 6912 are excluded by the alphabet. The owner's framing when commissioning this
work — that the theorem delivers the reduction — was **wrong**, and the agent corrected the
attribution rather than accepting it. The non-monomial space remains **a wall, not a sweep**, now
with the diagonal slice exactly measured.

**Certificate standard held:** 15 non-achievable side-instances carry kissat 4.0.4 DRAT→LRAT
proofs accepted by **both** `drat-trim` and `lrat-check` (frozen build `2e3b2dc`), plus one
**CaDiCaL 3.0.1** cross-verification; 3 achievable Sun floors carry explicit witnesses. Anchors
were reproduced **first** — 55/58/56/59/60, 729/729 Brent identities each. All validity checks
over $\mathbb Z$ (`fmpz`); no $\mathbb F_2$ encoding entered the set, closing the
arXiv:2607.29291 under-constraining-ring hazard by construction.

**Agent's own rule-5 correction, owner-confirmed:** *"ternary unimodular group"* is a **misnomer**
— the 6960-element set does **not** close under multiplication. The owner found an explicit
witness product with an entry of **3**, i.e. outside $\{-1,0,1\}$ entirely. The 48 monomials
**do** close, as the signed permutation group $B_3$ of order $2^3\cdot3!=48$. The per-element
orientation action remains sound; only the word "group" was wrong.

**◆ OFF-DIAGONAL CENSUS COMPLETE (`paper55`, $\sigma=\mathrm{Id}$) — "the non-monomial space is a
wall" is FALSE, and the wall was the owner's sentence too.** Agent `Mm3OffDiag`, owner-verified.

**Correction to an owner figure first:** the data-node count is $6960/48=\mathbf{145}$, not the
owner's 1160 — the owner quotiented by the 6 permutations only, while the full monomial group is
$B_3$ with $|B_3|=2^3\cdot3!=48$, so the 8 sign matrices fold into the node as well. $145^3 =
3{,}048{,}625$ data-triples, $145^2=21{,}025$ pair-table cells.

**Method, and the factorization is PROVEN not sampled.** Outer-monomial invisibility is structural
via $(a_1p_1)^{-\mathsf T}=p_1^{-\mathsf T}a_1^{-\mathsf T}$, so ternarity, $d$, floors and Brent are
functions of data-pairs only. Controls still run: **two complete $145\times145$ rows enumerated
directly through the real map — 42,050 triples, ZERO mismatches** against table prediction; plus
200/200 random-triple concordance and 25/25 independent factorization tests. $\sigma$-invariance of
the all-ternary predicate proven structurally, 40/40 sampled.

**Counts (all owner-re-verified exactly):** survivor data-triples $=\mathbf{5{,}700}$ of
$3{,}048{,}625$ ($0.187\%$). Classification by monomial-node count: **1** all-monomial, **49** with
two, **825** with one, **4,825** with none — summing to 5,700, so **5,699 survivors involve at least
one non-monomial factor.** Touchstone: the single all-monomial triple $\times48^3=110{,}592$,
matching the owner's independently derived figure exactly. Scheme instances: $5{,}700\times110{,}592
= 630{,}374{,}400$ per decomposition at $\sigma=\mathrm{Id}$; $\times3\ \sigma = 1{,}891{,}123{,}200$;
$\times2$ decompositions $=\mathbf{3{,}782{,}246{,}400}$ over the named set.

**The monomial orientations are $1/5700 = 0.0175\%$ of the valid population.** The owner's falsified
per-factor prediction ($331{,}776$ per decomposition) is off by a factor of **exactly 5,700** — not a
coincidence, but the survivor-data-triple count itself, since that prediction was precisely the one
all-monomial triple. **The "wall" framing was inherited verbatim by three successive campaigns and
by the owner's own commissioning brief; the census replaced it with a number.**

**THE CORRECTED ISOTROPY FAMILY — the single most reusable artifact of this thread.** For the
$3\times3$ matmul tensor, whose trilinear form is
$F(U,V,W)=\sum_{i,k,j}U_{ik}V_{kj}W_{ij}=\langle UV,W\rangle$ — **not** $\mathrm{tr}(UVW)$ — the
tensor-preserving action, in **standard matrix-product convention**, blockwise on the $23$ $3\times3$
blocks, is
$$U'=P^{-1}UQ^{-\mathsf T},\qquad V'=Q^{\mathsf T}VR^{-\mathsf T},\qquad W'=P^{\mathsf T}WR$$
for **any** invertible $P,Q,R$. **Independently confirmed unique** by two separate $4^6=4096$-sheet
enumerations — the agent's and the owner's — as the **only** all-ternary + Brent-passing three-letter
assignment on the shared `load_sun56()` arrays. Owner also verified invariance symbolically with
fully symbolic matrices and numerically at **200/200** random integer sextuples, against **0/200**
for the naive cyclic-conjugation family. **Record the convention with the formula:** `gatec_sweep`
uses transposed-$L$ semantics ($\sum L[i][i_2]M[i][k]R[k_2][k]=L^{\mathsf T}MR^{\mathsf T}$), and a full
day was spent on a phantom disagreement that was purely this transcription difference — both sides
had the same map and the same arrays throughout.

**◆◆ OFF-DIAGONAL PROGRAMME COMPLETE, BOTH DECOMPOSITIONS — the published orientation is the
UNIQUE minimizer of the entire landscape, 2026-08-30.** Owner-verified arithmetic throughout.

| | `paper55` | `sun56` |
|---|---|---|
| survivor data-triples | 5,700 | 4,510 |
| decided | all, 929 s | all, 809 s |
| **minimum certified total** | **55** | **55** |
| triples attaining 55 | **exactly 1** (all-monomial) | **exactly 1** (all-monomial) |
| minimum among non-monomial-involving | **56** (5,699 triples) | **56** (4,509 triples) |
| median | 66 | 68 |
| spread | 55-74 | 55-75 |
| scheme-instances | 1,891,123,200 | 1,496,309,760 |

**Histograms sum exactly** (5,700 and 4,510, owner-checked). **ZERO data-triples at $\le54$ in
either.** Named-set total **3,387,432,960** scheme-instances — a **5,105$\times$** extension of the
certified no-go surface over the inherited monomial-only $48^3\times3\times2=663{,}552$
(owner-verified exact: $663{,}552\times5105=3{,}387{,}432{,}960$).

**The structural result, and it is the interesting one.** In **both** decompositions the value 55 is
attained by **exactly one** data-triple — the all-monomial one — and **every** orientation involving
any non-monomial factor certifies $\ge56$. So **the published record orientation is the unique
minimizer of the whole off-diagonal ternary landscape in the named set**, with the medians sitting
11 and 13 additions above it. The region that three campaigns and the owner had called *"a wall"*
turns out to be wide open, exactly countable, and uniformly worse than the known scheme.

**Verification battery before freeze, against the FROZEN reference functions:** frozen `d_count`
cross-check 12/12, floor-verdict re-run 6/6, survivor full-map (ternary + Brent) recheck 8/8 — all
green. Anchor **asserted at startup** rather than checked afterwards, per the instruction issued
after the `paper55` bookkeeping alarm.

**Undecided remainder: empty** for these two decompositions — both decision layers ran to completion
inside caps. **Not searched, and no universal no-54 claim is made:** other decompositions,
non-ternary alphabets, and $\mathrm{GL}(3,\mathbb Q)$ sandwiches. **No escalation condition fired
anywhere in the programme.**

**◆ `paper55` decision layer detail, 2026-08-30.** All 5,700 decided in 929 s, well inside the 24 h cap.

**Histogram of certified totals** (owner-verified to sum to exactly 5,700): 55:1, 56:2, 57:3, 58:13,
59:26, 60:54, 61:129, 62:248, 63:388, 64:670, 65:824, 66:895, 67:829, 68:689, 69:470, 70:268,
71:126, 72:46, 73:15, 74:4. **Minimum $=55$, attained by exactly ONE data-triple — the all-monomial
one, i.e. the known 55-addition scheme. ZERO data-triples at $\le54$; 5,700 of 5,700 are $\ge55$.**
**Median $=66$**, so the 55-point is the extreme low tail attained precisely by the published scheme,
and the surrounding off-diagonal landscape is markedly worse.

**So the certified no-go extends from a single monomial point to the entire 5,700-triple off-diagonal
survivor set for `paper55`** — $630{,}374{,}400$ scheme-instances at $\sigma=\mathrm{Id}$,
$1{,}891{,}123{,}200$ across the $\sigma$-classes — with certified minimum exactly the published
record value.

**A false alarm the agent raised against itself and retracted inline (rule 5), caught by the
mandatory anchor.** The decision loop wrote the per-side sum into the `total` field **without the
$+14$ transposition gap**, so the live $\le54$ trigger fired repeatedly on side-sums (minimum 41).
What exposed it was the anchor requirement: the all-monomial data-triple **must** certify to exactly
55, and the code printed 41 for it. **The escalation protocol never had to fire** — no true total was
at or below 54, and the DFS decisions ($d$-values and floor verdicts) were unaffected; only the
summary field was wrong. This is the second time today an anchor-on-a-known-value caught a defect
that the surrounding numbers looked consistent with.

**Decision layer green-lit** under pre-registered caps ($\approx$2 h estimated, 24 h DFS cap), with
`sun56` to follow on identical machinery. Any total $\le54$ escalates before anything is written.

**Named next opening:** the off-diagonal sandwich product — $6960^3\times3=1{,}011{,}460{,}608{,}000$ **per decomposition**, hence $\mathbf{2{,}022{,}921{,}216{,}000}$ across both `paper55` and `sun56` of the named set (**scope corrected 2026-08-30**: the owner carried the per-decomposition figure into a sentence describing the two-decomposition set, understating the unswept region by exactly $2\times$; caught by agent `Mm3OffDiag` before its census was designed) — sandwich
product is **the only remaining place inside the ternary alphabet where $\le54$ could hide**, and
it was pre-registered as budget-out rather than quietly omitted. Frozen at
`mm3/campaigns/2026-08-30T174243Z_edbdd408-840f-42b4-9be0-192a95783ab9/`.

### 2026-08-31 extension — `mws59` and `stapleton60` fully decided; no new LB-55 row

The exact standard-action census now extends the off-diagonal programme to the
two remaining published decompositions selected for this gate:

| decomposition | survivor data triples | valid orientations | LB range | upper median | minimizers | nonmonomial min | at most 54 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `mws59` | 5,796 | 1,922,973,696 | 56–74 | 67 | 2 | 56 | 0 |
| `stapleton60` | 4,800 | 1,592,524,800 | 58–76 | 67 | 9 | 58 | 0 |
| aggregate | 10,596 | 3,515,498,496 | 56–76 | 67 | 2 | 56 | 0 |

No new lower-bound-55 data triple enters from `mws59` or `stapleton60`.
`paper55` remains the only **verified 55-addition circuit** in the expanded
swept set, while `sun56` remains LB 55 / UB 56. The emphasized distinction is
load-bearing: `paper55` is not the unique LB-55 data-triple minimizer because
`sun56` also has an LB-55 all-monomial triple.

All values are exact integer enumeration plus complete finite subset DFS.
They are lower bounds, not synthesized circuits. The run reproduces all five
published totals 55/58/56/59/60 and 729/729 Brent identities before the new
counts. An independent 476-assertion audit regenerates pair masks, survivor
sets, side bounds, histograms, medians, minimizers, crosschecks, and aggregate
counts. Its first owner-authored attempt made a postcompute schema assumption,
aborted without a verdict, and left an empty write-once file; the corrected
retry and failed attempt are both preserved. Final checksum replay passes
31/31.

Scope is exactly `mws59` and `stapleton60`, every ternary
determinant-`+/-1` `(P,Q,R) in T^3`, `|T|=6960`, and all three sigma powers
under
$U'=P^{-1}UQ^{-\mathsf T}$,
$V'=Q^{\mathsf T}VR^{-\mathsf T}$,
$W'=P^{\mathsf T}WR$.
Other decompositions, alphabets, general rational/integer sandwiches, other
actions, anti-cyclic swaps, and upper-bound synthesis remain outside scope.
No universal no-54 statement follows. Evidence:
`mm3/campaigns/2026-08-31T083323Z_55447c11-5887-4fdd-aa2f-30e31f2332e8_b3046e6a6c50/`;
final-manifest SHA-256
`f33c416199bf9e9f661d5733cd0c23400932943b9f7e7ad5456dfb852d276bcc`.

## `kg/` gate B — two prior band PASSes RETRACTED, 2026-08-30 (agent `KgGateBClose`, owner-verified)

**Outcome: the three OPEN bands were NOT closed, and two bands previously recorded as PASS are
withdrawn.** The cause is a predecessor defect **unsound in the permissive direction**, found by
the agent, confirmed exactly by the owner. **No paper claim is affected; nothing approaches a
refutation of Lemma D.2/D.4; zero escalations.**

**Defect 1, the load-bearing one.** `j_integrand_bound` capped $|e|\le1$. But at the paper's own
D.4 reproducing kernel $p_0=(3-s^2)/\sqrt6$ (`paper_full.txt` line 2036),
$|e(0)|=|a_0-a_2/\sqrt2|=\sqrt6/2=\sqrt{3/2}\approx1.2247$. **Owner-verified exactly in sympy:**
$a_0=\sqrt6/3$, $a_2=-\sqrt3/3$, and $a_0-a_2/\sqrt2=\sqrt6/2$ identically. A cap that
**understates** $|e|$ understates $J$ and therefore makes $J\le d$ appear to hold — it can
manufacture false PASSes. **Defect 2:** 3-D circumradius inflation used in place of the paper's
2-D even-box envelope (lines 1972-1981).

**Band verdicts after correction:**

| band | before | after | evidence |
|---|---|---|---|
| $[1.30,1.45]$ | PASS | **PASS re-confirmed** | mandatory anchor re-run under the corrected envelope: certified margin $+0.1930$, cell $(0.8,-0.6)\pm10^{-4}$, 512 panels, Arb 256-bit outward |
| $[3.50,4.083]$ | PASS | **RETRACTED → PARTIAL** | 8 pre-declared geometric $c$-tiles, 2257 certified evaluations; kernel-adjacent tiles certify to $+7.5\times10^{-5}$, but corner cell $[0.833,1.0]\times[-0.333,-0.167]$ is **STUCK at certified margin $-2.0$ to $-4.7\times10^{-5}$** under full branch-and-bound |
| $[4.083,6.0]$ | PASS | **RETRACTED → UNVERIFIED** | not re-run under the corrected envelope, so its PASS rests on the unsound cap |
| $[6.0,12.0]$ | OPEN | OPEN, one cell certified | razor cell at $c$-pair $(6.0,6.05)$: certified margin $+3.675\times10^{-5}$, $n=8192$ — the tightest mid-band certification in this repo, and consistent in magnitude with the paper's band-8 Table-1 margin $+2.539\times10^{-5}$ `[REPORTED]` |
| $[1.0,1.3]$, $[1.45,1.75]$, $[1.75,3.5]$ | OPEN | OPEN | not run, budget exhausted |

**CORRECTION — 2026-08-31, owner:** the first two “after” rows and the
$mq=t^2$ diagnosis below do not survive the frozen-source audit. Campaign
`kg/campaigns/20260831T080542Z_kg_band_close2/` was pre-registered in a
one-file commit, then stopped at its mandatory anchor: the exact as-run source gives

$$d(1.45)_{\rm lo}-U(1.30,B)_{\rm hi}
=-0.149911438966541\ldots<0,$$

not $+0.1930$. Owner checksum verification passed for all ten frozen artifacts,
and an independent replay failed at the same Arb interval. Stronger mechanism:
the box center is itself the exact unit cubic
$p=\frac45\psi_0-\frac35\psi_2=A-Bs^2$. Splitting its Gaussian integral at the
two explicit roots gives, at 300-bit Arb,
$J(1.30,p)=0.3748305545887635700\ldots$ while
$d(1.45)=0.3726816997945499091\ldots$. Hence
$d(1.45)-J(1.30,p)=-0.0021488547942136609\ldots$: **no sound upper envelope
on a box containing its center can possibly have positive margin for that coarse
endpoint pair.** A separate 80-digit direct quadrature agrees with the closed
form. This does not challenge $J(c,p)\le d(c)$ at equal $c$:
$d(1.30)-J(1.30,p)=+0.0322365204873065\ldots$.

The predecessor's frozen `b5_result.json` also has **zero** certified
$[3.50,4.083]$ tiles and unresolved margins about $-0.0283$ through $-0.0194$,
not seven of eight tiles with $10^{-5}$-scale misses. Thus `[1.30,1.45]` is
**NOT CERTIFIED**, `[3.50,4.083]` is **UNVERIFIED/OPEN**, and no full band was
run on 31 August. Exact CDF panel masses were already present; $\phi(p_b)$ is a
lower, not upper, density rectangle on positive panels; and Lemma D.4 directly
permits $f(E_P,W_P;t)$, so the unsupported $mq=t^2$ branch is unnecessary.
Nothing here is a paper refutation. A new immutable campaign is restarting from
the analytic equal-$c$ control and ratio-1.02 tiles; no result from it is imported
before its mandatory anchors pass.

**DIRECT-D.4 RESTART RESULT — 2026-08-31, owner-inspected.** Immutable campaign
`kg/campaigns/20260831T082425Z_kg_direct_d4_restart/` replaces every superseded
band statement above. The direct 256-bit, 132-root-box driver gives:

| band | restart result | worst/closest evidence |
|---|---|---|
| `[1,1.3]` | **PASS 14/14** | worst certified lower margin `+9.495739844715277e-7` |
| `[1.30,1.45]` | **PASS 6/6** | `+1.95705024245827645e-6` |
| `[1.45,1.75]` | **PASS 10/10** | `+3.10559377267758493e-6` |
| `[1.75,3.5]` | **PARTIAL 25/36** | certified prefix through `2.81476518658164445646615405049` plus terminal tile from `3.49980671715929653996950394581`; closest open `-1.6374391446870114e-6` |
| `[3.5,4.083]` | **OPEN 0/8** | all eight hit panel cap; closest open `-1.5986209269591551e-5` |
| `[4.083,6]` | **OPEN 0/20 after cached refinement** | all reached depth 22/panel cap; closest open `-1.9726188208686084e-7` |

Negative numbers in the OPEN rows are lower bounds on an upper envelope's
margin, not lower bounds on `J-d`; they certify neither failure of Lemma D.4
nor a counterexample. No refutation trigger fired. The `[3.5,4.083]` result
SHA-256 is
`35030d91f4a220f348991cd6731951118f8b7887d5bd47991e219b510980968e`.
The `[4.083,6]` cached result SHA-256 is
`7e144344c11127ef3a069ca28053719b50ba1b9d730de616c22981c349f61b3c`.
The registered widened-interior refinement is now the sole
low-priority/thread-one process. Gate B is still partial.

Process disclosure: the cached-highband staging agent deleted three
agent-generated `__pycache__` files without required confirmation. They were
not restored; the campaign correction records the deletion and makes clear
that those compiled caches and staging probes are not evidence.

**MIDBAND REFINEMENT RESULT — 2026-08-31.** The widened unresolved-interior
campaign certified 6/12 tiles. Five adjacent tiles extend the continuous
`[1.75,3.5]` certified prefix through
`3.107728208020454953573248`; one terminal tile
`[3.49980671715929099891398176911,3.49980671715930]` also certifies. The
remaining contiguous gap is
`[3.107728208020454953573248,3.49980671715929653996950394581]`.
Smallest certified lower margin:
`+3.51851846436958939716186293650e-9`; closest open margin:
`-4.52960421419799285256261953577e-8`. Result SHA-256:
`f97604585aba167e583e8d3c8859e25442db5e09c7a482516715af02d540592d`.
The band remains **PARTIAL**, not refuted.

Mechanism correction: the named “add a panel-wise `K_o`-refined odd radius”
successor is **RETRACTED as new work**. Both baseline and cached runners
already apply that factor panel by panel. A successor must strengthen the
enclosure itself, for example by certified hinge splitting or coupled
even/odd maximization; repeating the existing radius is not a new mechanism.

**The stuck cells are our looseness, not a violation, and this is measured rather than asserted.**
Independent float probes at the stuck corner give $\max_p J(3.75,p)=0.14467$ against
$d(3.75)=0.15684$ — true spare $+0.0122$, roughly **600$\times$** the paper's own typical margin.
So the envelope excess is $\approx0.012$ while the truth has ample room; **no cell anywhere is near
refuting the paper.** Dominant remaining looseness is the fold-inflation $q$-hinge firing
whole-panel when $mq$ crosses $t^2$ **inside** a panel; the named fix is panel subdivision at the
$mq=t^2$ crossing, **pre-registered and deliberately NOT adopted mid-campaign** per rule 16.

**Disclosure the agent volunteered:** Addendum 2 (ordered-pair $q$-scheme) postdates the `b5 v2`
run it governs, and the addendum says so. A pre-registration timing gap declared by the agent
against itself.

## `oct-rank/` S3 — filed negative with a quantified obstruction, 2026-08-30 (agent `OctRankS3`, owner-verified)

**Outcome: 14 not reached; S3 stays OPEN but is now a *filed* negative with a named obstruction
profile rather than an untouched item.** Best certified lower bound **13** (1 substitution peel +
pencil floor 12 — exactly the published chain's arithmetic), **gap to 14 = 1**. Reaching 14 would
have lifted the published $18\le \mathrm R_{\mathbb R}(T_{\mathbb O})\le25$ to 19; the window is
**untouched** and no escalation was triggered, correctly, since a negative is not an improvement.

**Anchors reproduced first (MACHINE-VERIFIED):** $n{=}2\Rightarrow3$, $n{=}4\Rightarrow8$,
$\mathrm R_{\mathbb R}(\tau)=7$. **Convention stated and owner-re-verified:** Cayley-Dickson doubling
$\mathbb H\oplus\mathbb H\ell$, basis $(1,i,j,k,\ell,i\ell,j\ell,k\ell)$, byte-identical to the
gate-A table; owner independently confirmed $L_{\bar u}L_u = N(u)I$ exactly.

| route | computed | reaches 14? |
|---|---|---|
| **A** Strassen commutator | $\mathrm{rank}[L_{\bar u}L_v,L_{\bar u}L_w]=8$ on **173/173** triples; owner reproduced 8 on 5 independent basis triples | **No certificate.** Twin readings $8{+}8=16$ or $8{+}4=12$; **12 is below the chain's own 13**, so the safe reading yields nothing, and the naked twin that would give 16 **overshoots the $\tau$ control at $n{=}4$ (predicts $4{+}4=8$ against the true 7)**. **Agent's precision, adopted:** the $\tau$ control refutes the $n{=}4$ twin reading; the $16$-at-$n{=}8$ claim fails on **provenance** (unread Strassen/Bläser form) and on the **universal quantifier**, i.e. a *missing certificate* rather than a measured contradiction at $n{=}8$. A successor chasing $\ge16$ needs **both** a first-hand reading of the exact 3-slice hypotheses **and** a rank argument (not a sweep) for the commutator bound, with the $\tau$ control as the mandatory cross-check. Reported as the **only live candidate** for a future $\ge16$ route |
| **B** AFT substitution pivots | pivot-residual independence exact, 3 pivots $\times$ 56 triples; double-peel 10 | No — ceiling 13 |
| **C** residual-pencil floor | residual pencils satisfy the irreducible quadratic $C^2-2aC+bI=0$, $a^2<b$, entrywise in `fmpq` | No — and **this is the campaign's strongest structural result: the substitution route CANNOT exceed 13**, by an exact 181-checks-per-triple proof. A no-go on the *route*, which tells future work to abandon the avenue rather than retry it |
| **D** flattening ranks | $3/8/8$ exact on all 56 basis triples | No — **structurally** capped: for an $8\times8\times3$ tensor the caps are $\min(8,24)$, $\min(8,24)$, $\min(3,64)$ = 8, 8, 3, which the owner confirmed matches exactly |
| **E** exact ideal-theoretic infeasibility | BUDGET-FAIL as pre-registered (312-unknown Gröbner) | Not attempted to completion |
| **F** (post-hoc) quaternion-subalgebra blocks | global $T$-conjugation block-diagonalizes; the specific triple $(1,i,j)$ has rank $\le14$ exactly ($7{+}7$ concatenation of two $\tau$ copies) and $\ge13$ | No — an upper bound on one triple does not lift a universal statement; **$(1,i,j)$ is 13-or-14, open** |

**Labelling discipline worth noting:** the universal statement *"commutator rank $\ge8$ for ALL
independent triples"* is labelled **`COMPUTATIONAL-EVIDENCE`**, not machine-verified — 173/173 is
equivalence-class evidence, **not a universal quantifier**. The agent drew that line itself.
Strassen/Bläser/Lickteig twin forms are labelled `CITED-DEPENDENCY`, **unverified**, because they
were not read first-hand — and that honesty is exactly what kept the $16$ reading out of the record.

## `delcap/` Open Item 1 RESOLVED — 20 rows, both ends improved, 2026-08-30 (agent `DelcapNextGate`, owner-verified from the frozen artifact)

**Outcome: all 20 pre-registered rows certified, and every one strictly improves BOTH ends of the
published TNB sandwich.** `MACHINE-VERIFIED`: exact rational snap plus Arb 400-bit outward-rounded
certificates, 9 sha256 checksums verified. **No published value is excluded by any interval, so no
escalation — the intervals sit strictly inside the published sandwich, which is what improvement on
both ends means.**

**Owner verification, read from `orbit_rows.jsonl` rather than the report:** for all 20 rows,
`cert_lo` $>\max(\texttt{lb1},\texttt{lbplus})$ and `cert_hi` $<\texttt{ub}$ — **20/20 on each
side**, matching the artifact's own `verdict` field
`CERT_UPPER_BEATS_UB+CERT_LOWER_BEATS_LBplus`. Certified width per symbol ranges
$9.862\times10^{-8}$ to $1.320\times10^{-4}$; all widths strictly positive.

| $(q,n,d)$ | certified interval (per symbol) | width | published $\texttt{lbplus}$ / $\texttt{ub}$ |
|---|---|---|---|
| $(3,10,1/2)$ — the wall row | $[0.42369976,\ 0.42371229]$ | $1.253\times10^{-5}$ | $0.32826401$ / $0.79248125$ |
| $(3,6,1/20)$ — tightest | $[1.41580723,\ 1.41580733]$ | $9.862\times10^{-8}$ | $1.41473940$ / $1.50571438$ |
| $(3,8,1/2)$ — widest | $[0.44916027,\ 0.44929226]$ | $1.320\times10^{-4}$ | $0.35735985$ / $0.79248125$ |

**Method:** sparse $G=S_q\times C_2$ orbit certificates on the **symbol-value** group — explicitly
**not** the falsified position/type reduction. Equivariance was anchored at runtime and
**independently confirmed by the owner** in exact rationals at $(n,q,d)=(3,2,1/3),(4,2,1/5),(3,3,2/7)$,
with a **position-permutation control that BREAKS** in the same test that confirms the value action —
so the live reduction and the dead one are separated by measurement, not assertion.

**Three bugs found and fixed mid-run, recorded inline — the second is the important one.**
(i) A snap remainder loop spun $\sim10^{11}$ iterations on partial-alphabet marginals; replaced by a
cycle-decomposed equivalent verified identical on 200 random + 60 adversarial trials.
(ii) **A rehearsal produced a NEGATIVE-width certificate** (dual below primal) because per-word $D'$
was snapped over only representative-emittable words, making it sub-stochastic — *"the exact
'invalid certificate that looks clean' shape,"* caught by a **width-sign check** and fixed by
extending the marginal to the full output alphabet via the orbit identity.
(iii) The first $+\infty$ re-verification used orbit-**aggregated** masses, which never vanish,
whereas the frozen failure mode lives at **per-word** granularity — i.e. the check was being run at
a granularity that could not see the thing it was testing, and the agent caught that itself.

**Rule-7 scope:** $q=3$ only; $(3,10)$ and $(3,n)$ for $n\in\{6,7,8,9\}$ at $d\in\{1/2,1/5,1/10,1/20\}$;
plus 2 dense-vs-orbit cross-validation rows. NOT swept: $q\ge4$ at $n>5$, any $n>10$, MD $m>23$, the
$\Lambda$-optimality gap at $m\in\{22,23\}$, Pinto-Ribeiro $n\in\{29,31\}$. **Every row is a finite-$n$
theorem; the asymptotic side remains open.**

## `mceliece/` gate B — Apon's $\Delta_{p,q}\ne0$ genericity, measured, 2026-08-30 (agent `MceliecelGateB`, owner-verified)

**Outcome: zero degeneracy events anywhere, across ~570 exact $\Delta$-verdicts. Apon's §3.6 hole is
BOUNDED, not closed** — no forcing configuration was found, and that is stated as plainly as the
positive findings. Neither disputant in the four-paper dispute had this rate.

**Three layers, deliberately NOT merged — and the first is a theorem, not a rate.**

| layer | result | evidence |
|---|---|---|
| **E1, adversarial, EXHAUSTIVE** | **0/2016** irreducible degree-2 Goppa polynomials over $\mathbb F_{64}$ at $(m,t)=(6,2)$ | `MACHINE-VERIFIED`, **no confidence interval** — owner verified the population is exactly $(q^2-q)/2=(4096-64)/2=2016$, so this is the **complete** cell. For $(m,t)=(6,2)$ the genericity condition holds for **every** irreducible Goppa polynomial: a finite exhaustive **theorem**, and a CI on it would be a category error |
| **sampled** | 0/250, exact one-sided 95% Clopper-Pearson $\le\mathbf{0.0119}$ | owner reproduced: $1-0.05^{1/250}=0.011911$. Bounds the rate over the **sampling distribution**, not over the space |
| **adversarial special-$G$ families** | 0/127, CP95 $\le\mathbf{0.0233}$ | owner reproduced: $0.023312$ |

Invariance V1 (ordering/basis) 5/5. Script corpus and all checkpoints content-hashed.

**Two owner sharpenings of the agent's own results.** (i) The $k=1$ boundary cell was reported
unreachable *"for $m\le11$, since $2^m\not\equiv1\bmod m$."* That congruence has **no solution for
any $m>1$** — no integer $m>1$ divides $2^m-1$ — so the cell is unreachable **unconditionally**, not
merely in the swept range. (ii) The agent **proved** that at full support $\Pi'\equiv1$,
$\lambda_i=G(a_i)^2$, the code is the value-image of $V_0=\{f: G(a_i)^2f(a_i)\in\mathbb F_2\}$, and
the Wronskian is $E$-bilinear — hence **the verdict is a pure function of $(m,t,G)$.** That
sharpens its own sampled layer: the seeds vary **neither support, basis, nor column order**, only
$G$, so the CP bound is a bound on the degeneracy rate over the **distribution of irreducible
Goppa polynomials** within each cell — a far more specific statement than "over the sampling
distribution", and a consequence of the agent's theorem rather than an assumption.

**Sampled layer composition:** 10 cells $(m,t)\in\{(6,2),(6,3),(7,2),(7,3),(8,2),(8,3),(9,2),(9,3),(10,3),(11,3)\}$
$\times$ seeds 1000-1024 $=250$ (owner-checked). **Special-$G$ families:** $Z^3+c$ (42/64 irreducible),
$Z^3+Z+c$ (21/64), $Z^2+Z+c$ with $\mathrm{Tr}(c)=1$ at $m=7$ (64/128) $=127$ (owner-checked).
**Boundary cells settled:** $t=1$ at full support is construct-impossible (the root of a degree-1 $G$
lies in $E$).

**The instrument was built new and tested against planted falsehoods — the control most campaigns
skip.** Rather than reuse the inherited route, the agent built a fresh verdict instrument (rank-scan
plus the completeness lemma $\deg_x\Delta\le D-1<n$) and self-tested it on **3000 property cases and
5 counterfactual PLANTS** against brute-Wronskian ground truth. **It caught a genuine iff-bug in its
own first cascade draft before any campaign number existed.** A test that only ever sees true
instances cannot distinguish a correct predicate from one that always returns `True`; planting
falsehoods is what gives a zero its meaning. It also rebuilt all 13 frozen instances **byte-exactly
from seeds** rather than trusting them, proved *and* machine-checked ordering/basis invariance, and
correctly identified that $\Psi=C\!\cdot\!J$ belongs to the waterfall census and is **not** on the
gate-B verdict path (re-deriving $\lambda_iF(a_i)=Y$ instead).

**Incident log kept, including against itself:** a reverted frozen-file touch (byte-verified) and 9
$m{=}11$ stale-module rebuilds (all nondegenerate), plus the seed-1018 dedup that moved the sampled
denominator $228\to250$.

### 2026-08-31 direct-route correction — exhaustive toy population fails P2–P4

The section above remains a theorem about Apon's older
$\Delta_{p,q}$ route. It is **not evidence** for current ePrint 1786
Assumption 1. A source-faithful direct-route campaign now exhausts a separate
toy population: all 24 monic `G_beta=Z+beta`, `beta=8,...,31`, over GF(32)
on ordered support `(0,...,7)` at
`(m,n,t,k,ell,n_ell,k_ell,D_ell,d,s,h)=(5,8,1,3,0,8,3,5,7,5,4)`.

All 24 cells and all 240 P1 branches pass P1_U3/P1_U4. Every cell fails each
remaining premise: P2 accepts all 33 projective labels rather than the five
expected; P3 leaves 28 labels unmapped; and P4 has
`rank(E)=36`, `nullity=244`, with both full and block admission false.
Thus `genericity_pass=false` in **24/24** cells. This is a complete finite
theorem for the registered toy population and an implementation stress test,
not a counterexample to the paper: the toy parameters are outside current
Table 1 and Assumption 1's five conditional Classic McEliece cells.

The scalar parent first passed its controls and recorded `beta=8,...,21`,
then stopped exactly at its 7,200-second cap. The exact table-indexed successor
completed all 24 in 1,622.94 seconds after 26/26 controls passed. Its RREF
tuple, nullspace order, and rank match the scalar reference on 1,004 systems;
cached/uncached `beta=8` matches the hard-hashed parent; and an additive strict
audit matches every semantic field across all 14 completed parent-prefix
cells. The strict audit is also the disclosed correction for a frozen C26
perimeter deviation: C26 omitted `complete` along with `header,t_utc`, whereas
the pre-statement allowed only the latter two. Both versions have
`complete=true`; the post-run audit checks it but is not relabeled as a
pre-production control.

Evidence:
`mceliece/campaigns/2026-08-31T16-12-46Z_vectorized-rref-replay/report.json`;
state manifest SHA-256
`ac21cdf3d11cb2844138e0bf459f6b53df3b457b78b81ab27285a40173ee0a79`.
All static, state, parent-partial, and final checksum ledgers pass.
Evidence grade: exact finite-field census, predicate, control, and equality
statements are **MACHINE-VERIFIED**; wall times are
**COMPUTATIONAL-EVIDENCE** only.

## `rs-pe3d/` P2 battery — 19/19 adjudicated, 2026-08-30 (agent `RsPe3dPatternP`, owner-verified)

`pre_statement.md` committed as the campaign's **first** file, with all 19 P-predictions and H-COP
fixed in advance. Every number is $\rho^{\rm window}$, **never** $\rho_{\rm inst}$.

**1. Pattern P SURVIVES — and is still not established.** Zero misses across 19. The arms are
reported separately because they are not equally informative:
- **PP arm (strong — predicts *exactly* $3/5$): 9/9.** New primes $q\in\{13,43,97,109,67,157,601,193\}$,
  $s=(2,3,k)$ for $k$ up to 32, $\beta$ up to 16; every one $\rho^{\rm window}=3/5$ exactly with
  $\mathrm{wt}=3$, $\delta=5$, per-witness exhaustive cheaper-tuple closure proven. Illustrative
  labelled null $(1/5)^9=5.12\times10^{-7}$ — **explicitly not a p-value**.
- **PN arm (weak — predicts merely "not $3/5$"): 7/7 landed.** PN8 $(11,(2,2,5))$ and PN9
  $(17,(2,2,8))$ both gave $1/2$, the battery's only new sub-1 value.
- **Out-of-sample record: B11 + 16 landing rows, 0 failures. The PN direction remains thin**, and for
  a biconditional the weak direction is where it usually dies.

**2. H-COP FALSIFIED — the owner's own coprimality observation, reframed and killed.** H-COP.a holds
(landed coprime rows lie in $\{1,3/5\}$); **H-COP.b fails** — landed non-coprime rows take
$\{3/5,1,1/2\}$, i.e. **two** distinct sub-1 values against a predicted $\ge3$. Falsification scoped
to the landed set. **Owner structural record: 0-for-13.**

**3. H-GATE PROSPECTIVELY CONFIRMED — and it is the battery's most durable output.** The
weight-$\le3$ window is non-empty **iff $\min(s)\le3$**, holding across all rows ever run with zero
exceptions (owner-verified: non-empty rows have $\min(s)\in\{2,3\}$, empty rows $\min(s)\in\{4,5\}$).
**PN10 was a committed prospective test** — addendum and order proof written at 49% census with zero hits so far — and came back **EMPTY as predicted**. **Agent-supplied precision that narrows this, recorded because it weakens the claim:** the crisp iff-form crystallized only *after* PN4/PN5/PN6/PN7 had already landed, so **H-GATE's prospective content is exactly and only the PN10 read** — one committed data point, not twenty. PN4/PN5/PN6/PN7 are labelled **retrospective** in the target README and the boundary is kept explicit there. The 20-row consistency is real but retrofitted; the single prospective test is what the claim rests on. Three complete emptiness censuses, owner-verified
to reproduce exactly as $\binom{L}{1}+\binom{L}{2}+\binom{L}{3}$: PN4 $457{,}450$ ($L=140$), PN5
$288{,}100$ ($L=120$), PN10 $3{,}658{,}900$ ($L=280$), **zero** nonzero supports in each. The
product-magnitude confound is **broken by the data**: $(3,5,16)$ has product 240 and a non-empty
window while $(4,5,7)$ has product 140 and is empty, so the larger instance is the non-empty one.

**Rule-15 significance, stated in the target README:** if H-GATE holds, this target's entire reported
domain is **$\{$instances with an order $\le3\}$** — every ratio ever reported here, including the
frozen 13, is conditional on that predicate. Mechanism from the dual product codes is **OPEN**;
corner cases $s_i=1$ untested. This is the same evidential status pattern P had after B11: one
committed prospective test, promising, **not** established.

**4. All three owner checks on the frozen table CONFIRMED** from the frozen witnesses themselves
rather than payload strings: $\rho=\mathrm{wt}/\delta$ in every row; P holds within each coprimality
class; and the OC4 description of the frozen 13. Clarification recorded: B04-B06 are **pairwise**
non-coprime (mutual gcd 1, a shared pair factor), and the classification is pairwise.

**5. Two of the agent's own reporting strings were stale post-freeze** (PN count across landing
waves; the illustrative null figure) — corrected as **rule-5 inline addenda** with frozen originals
untouched, and the vectorized empty-window kernel was cross-checked exactly against the frozen census
code on non-empty instances before being relied on.

**Two named next campaigns, agent-priced.** (a) The weight-$\ge4$ support window at a $\min(s)\le3$ instance — **the only way to grow the domain while H-GATE stands** — needing the vectorized kernel plus a checkpointing budget ($\approx11{,}700$ s census-alone at PN4's $N=140$, more at larger $N$). (b) The **dual product-code mechanism** for $\min(s)\ge4$ emptiness: a proof of the gate at $(4,5,7)$-class instances **would replace H-GATE's 20-row evidence with a theorem**.

**Not swept:** weight $\ge4$ supports — **everywhere open, and where the $\min(s)\ge4$ instances
actually live**; non-basis values in dim $>1$ intersections; other $t/\eta/\Lambda$; $\rho_{\rm inst}$
everywhere. TR26-150 untouched; no escalation ever warranted.

## Campaign ladder

All seven targets committed a `pre_statement.md` before their first computation.
Frozen campaign artifacts to date:

| target | frozen | decisive artifacts | verdict |
|---|---|---|---|
| `mm3/` | 7 | `2026-08-30T031544Z_e0f3f117_c9df97a8bf3a/` — kissat 4.0.4 DRAT proofs for all three floor CNFs, accepted by drat-trim **and** `lrat-check`, with CNFs, binary + LRAT proofs, verbatim checker output, checksums, tool sources, manifest; `2026-08-30T012035Z_a7ea8e8e_cbdf8e63fa94/` — the 59-vs-57 tally resolution; `2026-08-30T010534Z_610e7116_bf045517802b/` — three agreeing decision procedures at $T=d(F)$ | gates A+B+C complete; UNSATs **checker-certified** |
| `rs-pe3d/` | 6 | `2026-08-30T01-04-49Z_exacttable/` — 9-instance exact 3-D $\rho$ table with exhaustive weight-$\le3$ cheaper-support check; two documented premise failures preserved as dead ends | first exact 3-D table frozen; refinement running |
| `delcap/` | 5 + `enclosure_table.csv` | `2026-08-30T01:43:09Z_eabfc720f2ef/` — 12 Fertonani–Duman Table II entries enclosed in exact interval arithmetic; `2026-08-30T02:47:50Z_1811e9be088b/` — new certified $U(28,0.68)=0.13168$, plus the falsified type-channel reduction | gate A closed honestly at max certifiable $n$; enclosure table delivered |
| `omega/` | 4 | `2026-08-30T03:20:00Z_c9d4e2a8/` — rung 1 certified endpoint $2.3713400836689$ with the corrected $\varepsilon$ decomposition; `2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/` — rung 2 at $2.3715538358544617$, both arms' JSON byte-identical on re-run | gate A verdict; gate B rungs 1+2 CERTIFIED |
| `oct-rank/` | 4 | `2026-08-30T00:04:39ZZ_gateA_A6A4B0CB/` — 512 exact-rational Krawczyk margins, radius bisection to $\rho^\*\approx3.129\mathrm e{-6}$, Lean 4 v4.29.1 replay + axiom audit; `2026-08-30T01:55:00ZZ_gateB_ladder/` — 273-probe ladder with border-rank instrumentation | gate A PASS; gate B clean negative; gate C PASS |
| `kg/` | 3 | `20260830T010702Z_2a373482_385dbaab2fd3/` — gate A, the repo headline certificate; `20260830T013720Z_e24aa159_afa515bbc532/` — gate B interval cover with open bands stated | gate A PASS (tail input CITED); gate B PARTIAL |
| `mceliece/` | 2 | `2026-08-30T01-42Z_9D0FC9E7/` — 7 small-$t$ instances; `2026-08-30T04-47Z_98C586C7/` — the $m\in\{10,11\}$ ladder addendum on full-field supports | gate A PASS on 10 instances |

**31 frozen campaigns across seven targets.** One further `kg/` directory,
`20260830T030633Z_405bd000_d3hR3start/`, is a live provisional checkpoint for the
$\lVert D^3H\rVert$ re-derivation, not a verdict: its contents are explicitly
labelled COMPUTATIONAL-EVIDENCE, non-certified.

**Proof-log certification leg — CLOSED 2026-08-30**, superseding this section's
earlier "no target has yet produced a proof-log certificate". See
[Proof-log certification](#proof-log-certification-2026-08-30): `mm3/`'s three
floor-impossibility results carry kissat 4.0.4 DRAT proofs accepted by drat-trim
and independently by `lrat-check`. VeriPB v0.1.0's `recordclass` API break made
the pseudo-Boolean route a fallback; the DRAT route carried it instead.

**`kg/`'s imported tail is now INDEPENDENTLY CORROBORATED — resolved 2026-08-30.** The
headline's $\lVert D^3H\rVert\le14.44243664663976457$ entered as a CITED-DEPENDENCY, so
the bound was arithmetic-level on the source paper's own scheme rather than independently
reconstructed. Owner-computed robustness, at 300 bits: the tail budget for the headline to
survive is $b_1^{\rm low}-\gamma_{\rm paper}-\mathrm{head}=1.708418967\times10^{-5}$
against a cited tail of $4.57569\times10^{-6}$ — slack $3.7336860\times$, i.e. **the
headline holds for any $\lVert D^3H\rVert\le53.923523$**, and still beats Krivine up to
$596.25224$. Since $\gamma_{\rm paper}=\pi/(2\times1.7818666069360661)$ to
$4.3\times10^{-18}$, "margin $>0$" and "beats the source paper's stated bound" are the same
condition; there is no separate cliff between them.

**Our own corrected route now lands below the paper's boundary.** After the Euler-operator
defect was found and fixed (see [the resolution](#resolution-the-3times-was-ours-and-the-paper-stands)),
our pipeline gives avg$(S)=123.377609\le126.80385221$ and $B_3=14.2459830<14.4424366$,
raising break-even slack from the cited $3.7336860\times$ to $3.7851739\times$. The cited
dependency is therefore corroborated from inside rather than merely consumed — at
`COMPUTATIONAL-EVIDENCE` grade only, the $\theta$-quadrature being uncertified float.

**Superseded, retained per Rule 5, uncitable**: the earlier reading of this item as an open
exposure with `KgD3H`'s float sweep at $2.8583\times$ above $126.80385221$ implying
$\lVert D^3H\rVert\approx24.4$, together with the fallback statement that "at $B_3=24.43$
the tail is $7.73998\times10^{-6}$ and $K_G\le1.78184771968824648$". Those figures are
properties of the plain-derivative pipeline. The fallback reasoning was nonetheless sound
and its conclusion held on every branch, which is why the headline never moved.

**Owner hypothesis FALSIFIED, recorded per rule 5.** I predicted the excess would
converge to exactly $3.000$ as a multiplicity/measure bookkeeping factor, on the
grounds that dividing the then-current $3.0488$ by $3$ or by $\pi$ landed within
$1.5\%$ of the cited constant. Refinement went $3.0488$ (9-point) $\to2.8604$
(13-point) $\to2.8583$ (49-point), passing **through** $3$ downward and settling at
$2.86$ — neither $3$ nor $\pi$. The agent independently killed the same story from
the other side: the per-$(p,a)$ decomposition gives diagonal norms $48.95$ plus
PD-verified cross-terms $+79.3$, which are legitimate $\langle G_i,G_j\rangle$
overlaps rather than a counting fault, and the $(48)$-bridge series identity is
exact including the $b=0$ term (Mehler-confirmed; the earlier apparent violation
was the agent's own omission in its analytic comparison, documented in its notes).

**Quadrature objection raised and ANSWERED; the discrepancy is now localized to
one additive term.** I objected that the $\theta$-sweep could not be trusted: the
integrand carries $494\times$ dynamic range with layers at $\theta=0,\pi$, and in
the 9-point rule the endpoint panels supplied $31.76\%$ of the mass from 2 of 9
nodes, so an unresolved boundary layer — or an integrable singularity, which the
agent's own T-trace route exhibits at $\lvert\lambda\rvert\to1$, the same endpoint
— would make uniform-grid "stability" an artifact rather than convergence. **The
agent answered it and I withdrew the objection.** $S$ *saturates*: $S(0)=334.347$,
$S(\pi)=1424.281$, flat to 6 digits from $\theta=10^{-3}$ inward, so there is no
singularity in the $\mu\times\mu$-norm integrand; the T-route divergence belongs to
a diagonal-Hermite **line** functional, a different object, so the apparent
contradiction was in my reading and not in the mathematics. Endpoint panel shares
fall linearly ($15.4/40.7\to8.1/28.7\to1.9/8.1\%$ at $h=\pi/8,\pi/12,\pi/48$),
first-order as a trapezoid on a smooth integrand must, and the three integrals
agree to $1\%$. Converged: avg$(S)=362.5\pm1$, excess $2.85875$,
$\lVert D^3H\rVert\approx24.419$.

With the quadrature cleared, the disagreement lives in the *object*, and the agent
has now measured the leg split over the full 49-point grid:
avg$(\lVert\Phi_3\rVert^2)=18.2592$, avg$(\lVert\partial_x\partial_y\Phi_3\rVert^2)=86.0457$,
summing to $18.2592+4(86.0457)=362.4420$ against the reported avg$(S)=362.4418$.
The first leg is only $5.038\%$ of the mass and gives $B_3=5.48044$ alone, so
**the entire excess is one scalar on one object**: the paper's certificate tolerates
avg$(\lVert\partial_x\partial_y\Phi_3\rVert^2)\le(126.80385221-18.2592)/4=27.1362$,
against the agent's $86.0457$ — a ratio of $3.17089$.

**Two further owner hypotheses RETRACTED, recorded per rule 5.** (i) I argued the
multiplicity $4$ in $S=\lVert\Phi_3\rVert^2+4\lVert\partial_x\partial_y\Phi_3\rVert^2$
was the unexplained constant, since third-derivative multinomial coefficients in two
variables are $1,3,3,1$. **Wrong — a category error.** The $4$ is the Cauchy–Schwarz
constant of the source paper's own chain: from
$\lvert T[\Psi]\rvert\le\lVert\Psi\rVert+(\pi/\sqrt6)\lVert\partial_x\partial_y\Psi\rVert$,
writing $x+\pi y/\sqrt6=\langle(1,\pi/(2\sqrt6)),(x,2y)\rangle$ gives
$(x+\pi y/\sqrt6)^2\le(1+\pi^2/24)(x^2+4y^2)$, and $(\pi/(2\sqrt6))^2=\pi^2/24$
**exactly** (verified at 200 bits), with $1+\pi^2/24=1.41123351671\le\pi^2/6=1.64493406685$.
The $4$ pairs with the $2y$ in the split; the multinomial coefficients govern the
$D^3$ split, not this pairing. The agent had already concluded the 4-factor structure
was right, and I overrode a correct finding of its own. (ii) The **double-count**
hypothesis is dead: it required the mixed term to be redundant with
$\lVert\Phi_3\rVert^2$, but line 1854 defines the two legs independently. (iii) My
$\theta=0$ extrapolation — projecting avg(leg1)$\,=139$ from a single node's $38.3\%$
share — was refuted by measurement at $5.038\%$, since leg1 varies $128\to1.6\to14.7$
across $\theta$. Projecting a $494\times$-dynamic-range integrand from one node was
the same error I had just warned the agent against.

**Three audits were run and all passed — and (1) and (2) are now VOID as evidence,
because audit passes do not transfer across a defect in the layer they exercise.**
Both exercised `rho_eval`/`rho_j`, the evaluator layer where defect (b) below lives, and a
finite-difference check of a wrong function against the analytic derivative of *that
same wrong function* passes by construction: it tests internal consistency, not
correctness. They may not be cited as evidence for anything; any future claim resting
on them must re-run them on corrected code. Recorded rather than deleted, per rule 5.
(1) The
$a\to a{+}1$ shift: coefficients $c_{p,a}$ read unshifted against $G_{p,a+1}$, shift
applied exactly once, index set extended rather than reused or truncated, no $c$-row
re-read at the shifted index — PASS. (2) The $\partial_x\partial_y$ operator: $a{=}4$
validated against finite-difference-of-$a{=}3$ to $2.7\times10^{-8}$, and
$\mathbb E[\psi_3'^2]=3$ exactly — PASS. (3) The agent's own flagged item $(\gamma)$,
leg2's inner 2-D quadrature convergence, which the owner reopened after it was
skipped in favour of object-level forensics: avg(leg2) $=88.244,86.429,86.046,86.002$
at $N=70,100,140,200$, **identical to 4 decimals at both $\lvert x\rvert\le13$ and
$\lvert x\rvert\le16$** (Hermite weights beyond 13 are below $10^{-37}$), so
truncation error is zero and the $N$-drift is a monotone $2.541\%$. Owner-verified by
geometric extrapolation: last-two decrement ratio $0.1149$, tail $-0.005711$, limit
$85.996289$, consistent with the agent's $86.00\pm0.06$. **The excess is not a
quadrature artifact *of that implementation* — but see the defect finding below,
which removes the implementation itself from evidence.** Also resolved: the $86.0457$ figure came from $N=140$,
$\lvert x\rvert\le14$, and the 13-vs-16 identity makes the earlier record ambiguity
harmless. **Scope of (3), stated precisely:** the 8-cell table covers the
$\theta$-*averaged* leg2, whereas item $(\gamma)$ as written in the agent's notes
(lines 302–304) named $\theta=0$. That point is separately settled and
cross-consistent — its partial trajectory $55.98\to52.68\to51.52$ has flattening
decrements ($-5.90\%$, then $-2.20\%$) and its endpoint agrees to $0.033\%$ with the
$51.5368$ implied independently by the reported $S(0)$ and leg1$(0)$. What has never
been checked at any node count is $\theta=\pi$, where leg2 $\approx352.395$ — roughly
$7\times$ the $\theta=0$ value, on the largest-magnitude part of the integrand. A
converged average largely implies it, since the average is dominated by that region,
but that is inference rather than measurement, and `KgLeg2Conv` is now required to
report per-$\theta$ node-convergence at both endpoints as part of its independent
route.

**No claim that the paper's tail constant is understated is available or asserted.**
The paper states only the aggregate at line 1864 with no leg-level or intermediate
values anywhere in lines 1827–1900, so the narrowest defensible finding — if it ever
lands — is "the certificate is not reconstructible at leg level from lines 1827–1854
plus the printed $c$-table", **not** that the source object differs. What is still
missing is the thing that made every other result in this repository hold up: a
second structurally independent route. leg1 has one ($\lVert G_{0,1}\rVert^2$,
`inner_c` versus direct, reconciled to $3.3\%$); leg2 = $86.00$ rests on a single
code path, and self-convergence of one implementation is not validation of it.

**The independent route was built, and it found two defects in our own machinery —
so the $3.17\times$ excess is withdrawn as an established quantity.** Agent
`KgLeg2Conv` implemented a structurally disjoint second route (Mehler $b$-series with
an exact $q$-ladder $\mathrm dq/\mathrm du=\sqrt{b+1}\,q_{b+1}$ plus Faà di Bruno
jets), self-validated harder than anything else on this target: the (47)/(48) bridges
to $6\times10^{-16}$, $\lVert G_{10}\rVert^2$ against the closed form
$4\pi^2\mathbb E[K^2]$ to 14 digits on three quadrature systems, per-$b$
Cauchy-integral jets exact to $10^{-15}$, F-matrix vs Gram vs direct-2D at
$10^{-16}$, flat window and node sweeps. Its values, COMPUTATIONAL-EVIDENCE:
leg1$(0)=1.807176$, leg2$(0)=34.960447$, leg1$(\pi)=51.137718$,
leg2$(\pi)=415.069782$. **Three assemblies mutually disagree** at $\theta=0$ — leg2
$=34.960$ / $42.442$ (`inner_c`) / $51.525$ (`legs49`, the frozen path behind the
$86.00$), a spread of $1.474\times$; and leg1 $=1.807$ / $127.97$ / $128.25$, a
spread of $71\times$. Two concrete defects in the predecessor pipeline are pinned:
**(a)** `d3h_legs49.cfun` row $(1,0)$ computes $\rho_3=(+12s_3^2+80s_5^2)/V=+1.44599$
where paper eq (8) gives $(-12s_3^2+80s_5^2)/V=-1.04798$ — a sign flip whose
amplitude ratio $1.3797878$ squares to $1.9038143$, matching the agent's stated
$1.904$ block inflation (legs are squared norms); **(b)** `d3h_kernels.rho_eval` and
`d3h_phi3.rho_j` implement $\rho(t)=t+C_3t^3+C_5t^5$ with $\rho(1)=0.89857$ instead
of the paper's $(t-s_3^2t^3+s_5^2t^5)/V$ with $\rho(1)=0.79217$, so the *pointwise*
evaluators sit on a different correlation curve from the *grid* evaluators, and the
(48)-bridge closes only at $0.79217$. Both are a sign on the same $s_3^2$ term at two
call sites; one upstream convention error is the parsimonious reading.

**DECIDER LANDED: one defect, two symptoms — confirmed.** `KgLeg2Conv` fixed both call
sites in a clean copy, reproducing the as-frozen anchor leg1$(0)=128.2480$,
leg2$(0)=51.5248$ exactly first as a control, and the corrected `legs49` assembly then
lands on the independently-validated route at relative $1.5\times10^{-7}$
(leg1$(0)$), $2.5\times10^{-7}$ (leg2$(0)$), $1.6\times10^{-10}$ (leg1$(\pi)$) and
$1.5\times10^{-12}$ (leg2$(\pi)$), with pointwise assembly-vs-jet sup-difference
$3.7\times10^{-14}$. **The circularity objection to testing against the second
route's own values is answered by overdetermination, not by trust:** the correction
has zero degrees of freedom — one sign and one missing $1/V$ — so it cannot be dialled
to hit four independent quantities at $10^{-10}$. **Defect (b) closes from a source
predating the investigation:** notes line 44 records, from `KgD3H`'s own first-hand
paper read, $\rho$-eff coefficients $(0.8936,-0.1039,+0.002487)$; owner-recomputed at
120 bits, $1/V-s_3^2/V+s_5^2/V=0.792187$ matches the paper value and
$1-s_3^2/V+s_5^2/V=0.898587$ matches the code value, the gap being exactly
$1-1/V=0.1064$. The defect is therefore precisely *the $1/V$ normalisation applied to
the cubic and quintic terms but not to the linear one*. Corrected endpoints
(COMPUTATIONAL-EVIDENCE): leg2$(0)=34.9604$, leg2$(\pi)=415.0698$. Three agent
self-retractions accepted: the $316$ residual was its own probe contaminated by
`rho_j`, the $218/409/7048$ intermediates were comparator-script incoherence, and the
$0.59\%$ herm gap is quadrature error on a squared non-polynomial integrand. So
$86.0457$, $18.2592$ and $362.44$ are properties of a defected pipeline, not of the
paper's $\Phi_3$. **The one live exposure is a read, not a compute:** both routes are
structurally disjoint in machinery but share the eq (8) transcription, corroborated by
notes line 44 yet still one document read twice. A verbatim eq (8) quote with line
number, resolving whether $1/V$ attaches to all three terms, is required before any
corrected sweep value may be stated against the paper.

**Transcription exposure CLOSED; corrected sweep ESCALATED, not recorded.** The eq (8)
check passed: `paper_full.txt` lines 687, 243 and 1721 all print
$\rho(t)=(t-s_3^2t^3+s_5^2t^5)/V$ with $V:=1+s_3^2+s_5^2$ at line 682, the LaTeX
`\frac` spanning the whole numerator, so $1/V$ attaches to **all three** terms — not
ambiguous, and matching notes line 44 from an earlier independent read. The corrected
49-point sweep then gave avg$(S)=379.7180$ against the paper's $126.80385221$, a ratio
of $2.9945305$; owner-verified at 150 bits, as are $B_3^2=624.61107$,
$B_3=24.992220$, and $2.15761\times$ slack to break-even. **[ALL FOUR OF THESE NUMBERS
ARE SUBSEQUENTLY RETRACTED — see "Resolution: the $3\times$ was ours" below. They are
properties of a pipeline using plain $\partial_t^j\rho$ where the paper's line 1826
requires the Euler $\rho_j=D_t^j\rho$. Corrected: avg$(S)=123.377609$,
$B_3=14.2459830$, slack $3.7851739\times$. Retained here per Rule 5, uncitable.]**
**The agent escalated without
concluding, per instruction, and nothing paper-side is asserted.** Three defects in the
*evidence* were then flagged by the owner: (i) the $S$-sweep **lost route independence
in its interior** — the endpoints carry two disjoint routes agreeing to
$10^{-10}\!-\!10^{-7}$, but the interior carries only the closed-form assembly (the
$b$-series diverges at $\lvert r\rvert=1$, $\theta=\pi/2$) with one spot-check at
$0.1\%$, so the $3\times$ discrepancy sits exactly in a single-implementation region —
the same failure mode that produced the retracted $86.0457$; (ii) the reported
$379.72/379.72$ at five significant figures cannot distinguish a converged rule from a
stalled one, so unrounded values and the $\theta$-rule family are required; (iii) if
the *object* rather than merely its series is singular at $\theta=\pi/2$, the
documented uniform-rule-near-integrable-singularity trap applies and the convergence
table is not yet evidence. **The pending decider is free:** the paper's $126.804$ lies
*below* our own $S(0)=141.65$, and for nonnegative $S$ no probability measure over
$\theta$ can average below the pointwise minimum — so if $\min_\theta S$ over the 49
nodes exceeds $126.804$, then no weighting, Jacobian or refinement can reconcile the
two, $\theta$-quadrature is eliminated, and the discrepancy is pointwise or
definitional in the $\Phi_3$-assembly layer (c-table, $D_t$ bookkeeping, the (52)
recursion) that the kernel-level $10^{-14}$ validation never touched. The ratio
$2.9945$ is **banned as a lead** and permitted only to order tests: it is $0.18\%$ from
$3$, which is not $3$, and owner guesses of exactly this shape are nought for two
today.

**Min-test run: the owner's elimination FAILED, and the same data delivered a stronger
one in the opposite direction.** $\min S=2.5945$ at $\theta=1.7671$ with 21 of 49 nodes
*below* $126.804$, so the pointwise-minimum argument does not fire and $\theta$-shape
remains admissible. Two of the three owner flags are then **withdrawn on the agent's own
evidence**: (a) the trapezoid decrement ratio is $71{,}794$ where $O(h^2)$ predicts $4$
and $O(h^4)$ predicts $16$ — that is *spectral* accuracy, exactly what a trapezoid gives
for a smooth function of $\cos\theta$ (even about both $\theta=0$ and $\theta=\pi$, so
$[0,\pi]$ with endpoints is half a periodic-analytic trapezoid on $[0,2\pi]$), so the
rule is in the right accuracy class and the singularity trap is closed; (b) the
single-route interior is **exonerated by magnitude** — band $[1.5,2.1]$ is $19.1\%$ of
the interval but carries $S\approx3$, contributing $0.573$ to a mean of $379.718041$
when reconciliation requires removing $252.914$, so it supplies at most $0.227\%$ of
what is needed and a $100\times$ error there moves the mean $15\%$. Conversely
$[1.7671,\pi]$ is $43.75\%$ of the interval and a linear proxy alone accounts for
$374.95$ of the mean: **the average lives almost entirely in the best-validated region**
($\theta=\pi$ dual-route at $6.3\times10^{-10}$, $3\pi/4$ at $0.25\%$). Independently,
our $\min\lvert1-\rho^2\rvert=0.3725$ reproduces the **paper-stated** $0.3724$,
corroborating the $\rho$/kernel layer against the source rather than against ourselves.
**One candidate survives on our side:** the $\Phi_3$-assembly layer (c-table, $D_t$
bookkeeping, the (52) recursion), consumed by *both* routes and never validated against
anything. Ordered next: first a verbatim read of the paper's own definition of the
quantity it certifies at $\le126.80385221$ — uniform $\theta$-average or some other
measure — since $S$ spans $660\times$ and a convention difference buys a factor of $3$
for free; then an independent $\Phi_3$ derivation from eq (52) with its own c-table,
compared pointwise at $\theta=0,3\pi/4,\pi$, a constant ratio indicting a multiplicity
or normalisation constant and a $\theta$-varying one indicting the recursion.
**Owner calibration, recorded:** process demands this session paid (the second-route
requirement found a real defect; escalation, audit-non-transfer and the numerology ban
all held), while owner structural hypotheses are **nought for three**. Effort
reallocated to verifying agent arithmetic and enforcing discipline over proposing
explanations.

**Paper's definition read and matched; then the $\Phi_3$ derivation found two printed
rows it disagrees with.** The certified object is confirmed as the *uniform*
$\theta$-average of $S$ over $[0,\pi]$ with no Jacobian and no weight —
`paper_full.txt` line 1854 defines
$S(\theta):=\lVert\Phi_3\rVert^2_{L^2(\mu)}+4\lVert\partial_x\partial_y\Phi_3\rVert^2_{L^2(\mu)}$,
line 1812 fixes $\mu$ as the product standard Gaussian, and lines 1861/1864 give
$(1/\pi)\int_0^\pi S\le126.80385221$. So our averaging convention is theirs and the
$660\times$ range of $S$ is intrinsic to their own object; the cheapest escape is
closed. The independent derivation from eq (52), with its own c-table read from lines
1835–1852, then reproduced **seven** printed rows exactly — $(1,0),(2,0),(3,0),(2,1),(1,2),(0,3),(0,2)$ —
and **disagreed on two**: $(0,1)$, printed $-t$ against a recursion that generates $0$;
and $(1,1)$, printed $-3t(\rho_1+\rho_2)$ against a derived $-3\rho_1-3t\rho_2$.
**Owner finding — the disagreement is undetectable exactly where it was tested.**
Assumption-free, printed minus derived on $(1,1)$ equals $3\rho_1(1-t)$: identically
zero at $t=+1$ *by construction*, so $\theta=0$ never had power to discriminate the two
forms, and equal to $6\rho_1$ at $t=-1$. Since $\theta$ near $\pi$ supplies roughly
$375$ of the $379.718041$ mean with the $4\times$ weight on leg2$(\pi)=415.069782$, a
row correct at $t=+1$ and wrong at $t=-1$ is invisible to the strongest existing check
and decisive for the headline number. **The agent's own conclusion that the discrepancy
is "not in the c-table print" is therefore withdrawn as self-contradictory**, and its
supporting evidence — that the paper-verbatim table reproduces the validated endpoints —
is circular, since those endpoints were themselves computed from that table. Pending
test, cheap and ordered first: recompute the full sweep with the *derived* rows and
report the direction of movement before any interpretation. The $\mu$-normalisation
search is **deferred and must be pre-registered** when it runs — enumerate candidates
from the paper's text with line numbers, commit the list, then report every value
including misses, because a search halting when it hits $126.804$ is numerology at the
level of conventions. The hypothesis that the paper's constant errs **may not be held**
while our own derivation disagrees with the paper's printed table.

**Derived-rows sweep run: direction is UP, so the c-table question is second-order.**
avg$(S)$ moves $379.718041\to397.667200$, ratio $2.9945308\to3.1360814$ — a movement of
only $4.73\%$, so the $3\times$ is **robust to the table** and lives in a layer neither
variant touches. The owner's row algebra was confirmed exactly: $6\rho_1(-1)=3.56573539$
against the agent's reported $3.5657354$, with $\rho_1(-1)=+0.5942892309$ and
$\rho_2(-1)=-0.5737406266$. The agent also retracted its own "identical at $t=1$ by
accident" reading; the agreement there is structural.

**Owner scoping error, corrected.** The paper makes *two separable* claims and this
thread has been conflating them. The **intermediate** (line 1864) is
avg$(S)\le126.80385221$. The **conclusion** (line 1861) is
$\lVert D^3H\rVert^2\le(\pi^2/6)\,$avg$(S)$, i.e. $\le208.583976$ and
$\lVert D^3H\rVert\le14.442437$ — which reproduces the paper's printed $14.4424$
exactly, so its own chain arithmetic is internally sound. **Line 1861 is an
inequality**, so a false intermediate does *not* imply a false conclusion: the bound may
simply be loose. Only the conclusion carries mathematical weight, and every hour spent
here has tested the intermediate instead. Accordingly the deferred **direct Parseval
route is now the instrument**: compute $\lVert D^3H\rVert^2_{L^2(\mathbb T)}$ as a
Parseval sum over $H$'s own coefficients $b_m$, sharing none of the contested machinery —
no $\Phi_3$, c-table, kernel layer, $\theta$-quadrature or $\mu$-normalisation —
discriminating $208.583976$ (paper) against the $624.611141$ our chain implies. Required:
the exact identity with its paper line number, certified widths, and an explicit
truncation order and tail bound as first-class numbers. The $\mu$-normalisation search
remains deferred and pre-registered; the "certified constant errs" hypothesis remains
unadvanced and, after this correction, was never a single hypothesis at all.

**Direct Parseval route run, and its refutation of the paper is VOID — owner-diagnosed.**
The agent reported $\sum_{m\le1255,\,m\ \rm odd} m^6|b_m|^2 = 290.3819745$ as a certified
*lower* bound, hence $\lVert D^3H\rVert\ge17.0406 > 14.4424$, contradicting the paper's
conclusion. **It is not a lower bound and there is no contradiction.** Per the standing
rule, a claimed refutation of a published claim was re-verified independently before it
could move, and the defect was found. The imported grid is not $a{+}b\le251$; it is the
odd-parity half of the **square** $a\le251$, $b\le251$ — $31752 = 252^2/2$ rows,
$\max(a{+}b)=501$. The recursion is $b_m=(\pi/2)\sum_{a+b\le m,\ a+b\ \rm odd}
(-1)^b A_{a,b}^2\,[t^{m-b}]\rho^a$, so the true $b_m$ requires **every** row with
$a{+}b\le m$. For $m\le251$ all such rows satisfy $a,b\le251$ and are present, so $b_m$
is exact; for **every** $m>251$ rows with $a>251$ or $b>251$ are required and missing,
and they carry sign $(-1)^b$, i.e. mixed sign. A partial-grid $b_m$ is therefore neither
an upper nor a lower bound on $|b_m|$, and $m^6|b_m^{\rm partial}|^2$ bounds nothing.
The agent's window argument ran the wrong direction: the question was never which $m$ the
grid reaches, but whether rows beyond the truncation feed the $m$ already summed.

Owner audit, `kg/scratch/owner_parseval_audit.py`, reproduces the artifact one truncation
level down where the correct answer is known. Validated first: the reimplemented recursion
returns $b_1 = 0.8815738220495995485818877$, matching the paper/A-grid value to
$2.65\times10^{-36}$ relative. Self-check: for $m\le125$ the beyond-truncation bucket
vanishes identically, $0$ rows nonzero, so the exactness edge is real and located as
claimed.

|Parseval window|grid $a,b\le125$|grid $a,b\le251$|status|
|---|---|---|---|
|$m\le125$|$0.777362899$|$0.777362899$|both exact, agree|
|$m\le251$|$8755.035274$|$0.777362899$|left column **is the agent's procedure**|

Truncating at $125$ and summing to $251$ — precisely what was done at $251$ summing to
$1255$ — inflates $0.777362899$ into $8755.04$, a factor of $1.1\times10^4$. Mechanism:
the true $b_m$ collapse by near-total inter-row cancellation ($b_{145}=1.55\times10^{-18}$,
$b_{251}=1.11\times10^{-25}$, exact), truncation destroys the cancellation and leaves an
uncancelled residue twelve orders too large ($b_{145}^{\rm trunc}=5.27\times10^{-6}$),
and $m^6$ then amplifies it — that single $m$ contributes $257.9$ against the agent's
entire reported total of $290.38$. **Retracted, recorded not deleted**: $290.3819745$,
$\lVert D^3H\rVert\ge17.0406$, and the implied contradiction of line 1861.

**What survives, and it is a genuine certified result**:
$\lVert D^3H\rVert^2_{L^2(\mathbb T)} \ge 0.777362899$ (`COMPUTATIONAL-EVIDENCE` pending
frozen certified widths), from the exact window $m\le251$ with all omitted terms
nonnegative. This is the agent's own `parseval_gridpart.json` value $0.77737$, computed
correctly and then discarded by extending past the edge. It is consistent with the paper's
$\le208.583976$.

**Second owner error, recorded.** Calling the Parseval route "the decider" was wrong.
Line 1861 is an **upper** bound, so a small true $\lVert D^3H\rVert^2$ is consistent with
*any* value of avg$(S)$ whatever: the route has exactly **zero** power over
avg$(S)\le126.80385221$. It tests the conclusion only, and the conclusion **survives**,
merely loose. The intermediate question is untouched and fully open. Owner structural
hypotheses now stand at 0-for-3 with two instrument-specification errors alongside.

**Next instrument, ordered as a test and not as an explanation**: unit-test the
$L^2(\mu)$ norm layer — the one layer both routes share beyond $\rho$ — against exact
Gaussian moments, feeding monomials $x^iy^j$ and comparing to
$\mathbb E[x^{2i}]\mathbb E[y^{2j}]=(2i-1)!!(2j-1)!!$, reporting every value including
every match. Ordering rationale, no more: $379.718041/3 = 126.5726803$ sits $0.183\%$
**below** the paper's $126.80385221$, exactly the slack expected of a certified upper
bound, so a clean factor of $3$ there would reconcile everything at once; and a $3$ has an
obvious manufacturing route, since $\mathbb E[x^4]=3$ for monomials and
monomial-versus-normalised-Hermite bookkeeping gives ratios $1$, $1.5$, $2.5$ at degrees
$1$, $2$, $3$. Motive and opportunity justify a measurement, not a belief.

**$L^2(\mu)$ norm layer tests CLEAN — owner hypothesis falsified by measurement.** All
$14$ monomial cases $\lVert x^iy^j\rVert^2$ match $(2i-1)!!(2j-1)!!$ at quadrature
resolution ($\lVert1\rVert^2=1$ to $4\times10^{-9}$, $\lVert x^2\rVert^2=3$ to
$9.5\times10^{-7}$, $\lVert x^3\rVert^2=15$ to $7.3\times10^{-6}$,
$\lVert x^3y^3\rVert^2=225$ to $1.5\times10^{-5}$, $\lVert x^4\rVert^2=105$ to
$4.0\times10^{-5}$); the kernel-weighted path matches an mpmath 30-dps ground truth at
$10^{-15}$; the moment ladder $\int x^k\,\mathrm d\mu$, $k=0..12$, is clean. An apparent
$97.5\%$ miss mid-test was the agent's own missing $\varphi(u)^2=e^{-u^2}/(2\pi)$
reference factor, caught and corrected inside the run and recorded rather than reported as
a finding. **Owner structural hypotheses now 0-for-4.**

**Discrepancy PROVABLY localized to leg2, by nonnegativity alone — no hypothesis.** From
the frozen 49-point sweep (`leg2asm_*`, which reproduces avg$(S)=379.718041$ exactly):

|quantity|value|share of avg$(S)$|
|---|---|---|
|avg(leg1) $=$ avg$\lVert\Phi_3\rVert^2$|$6.404064$|$1.69\%$|
|avg(leg2) $=$ avg$\lVert\partial_x\partial_y\Phi_3\rVert^2$|$93.328494$|—|
|$4\,$avg(leg2)|$373.313976$|$98.31\%$|

Both legs are squared norms, hence nonnegative. Setting leg1 **identically zero** — its
minimum possible value — still leaves avg$(S)=373.313976 = 2.9440\times$ the paper's
$126.80385221$; reconciling through leg1 would require avg(leg1) $=-246.51$. Therefore
**no correction to leg1 of any kind can close this gap**, and the entire discrepancy sits
in $\lVert\partial_x\partial_y\Phi_3\rVert^2$, which must shrink by at least
$3.100620\times$. This retires leg1 and every leg1-side c-table question permanently.

Note the corollary: the **only** validation the $\partial_x\partial_y$ layer ever had is
the finite-difference check already marked **VOID as evidence** (it exercised
`rho_eval`/`rho_j`, the defected layer, checking a wrong function against the analytic
derivative of that same wrong function). The layer where the discrepancy is now provably
localized therefore carries **zero surviving validation**.

**Test ordered, genuinely independent**: $\Phi_3$ itself is dual-route validated pointwise
(assembly vs jet, sup-diff $3.7\times10^{-14}$), so finite differences taken on $\Phi_3$
directly yield an independent $\partial_x\partial_y\Phi_3$ sharing no code with the leg2
path. (a) Pointwise ratio of FD-mixed-difference against the leg2 path's
$\partial_x\partial_y\Phi_3$ at several $(x,y)$ and $\theta\in\{0,\text{mid},\pi\}$ —
$\sqrt3$ in amplitude identifies the mechanism, $1.0$ exonerates the derivative
implementation and moves the fault to what gets normed. (b) $\lVert\partial_{xx}\Phi_3
\rVert^2$, $\lVert\partial_{xy}\Phi_3\rVert^2$, $\lVert\partial_{yy}\Phi_3\rVert^2$ as
three distinct numbers beside exactly what leg2 returns. Ordering rationale only: $3$ is
the number of distinct second-derivative components in two variables and the Hessian
Frobenius norm is $\lVert\partial_{xx}\rVert^2+2\lVert\partial_{xy}\rVert^2+
\lVert\partial_{yy}\rVert^2$, so if leg2 is silently such an aggregate that is the whole
factor; if leg2 is $\lVert\partial_{xy}\rVert^2$ alone, both candidates die and the next
seat is the derivative counting inside the assembly.

**Dead artifacts, do not cite**: `kg/scratch/logs/leg2fix_leg1_49.npy` and
`leg2fix_leg2_49.npy` are stale defected-pipeline output (avg(leg1) $=131481.8$,
avg(leg2) $=5.58\times10^7$). The live sweep is `leg2asm_*`.

**Both derivative tests came back clean, and then the ROOT CAUSE was found in the paper —
it is our defect.** (a) Pointwise: an independent mpmath-30-dps route taking FD mixed
differences of the closed-form $\Phi_3$ against the leg2 path's
$\partial_x\partial_y\Phi_3$ gives ratios $0.99999$–$1.00002$ across 12 points at
$\theta\in\{0,\pi/2,1.7671,\pi\}$ — not $\sqrt3$, not anything but $1.0$. The derivative
implementation is correct. (b) The three second-derivative norms, computed on independent
tensor grids with $\lVert\partial_{xx}\rVert^2=\lVert\partial_{yy}\rVert^2$ symmetry
passing exactly, give leg2$/\lVert\partial_{xy}\rVert^2\to1.001$, while the Hessian
Frobenius aggregate ratio is $5.33$ at $\theta=0$ and $3.72$ at $\pi$ — not constant, not
$3$. So leg2 **is** $\lVert\partial_x\partial_y\Phi_3\rVert^2$ alone: both owner candidates
dead, and the fault is above the derivative and norm layers, inside the $\Phi_3$ object.

**The paper's recursion is not underdetermined; it was being read in the wrong place.**
Line **1816**: $\Phi_0(t;x,y):=C_{\rho(t)}(x,y)$, $\Phi_{k+1}:=D_t\Phi_k-t\,\partial_x
\partial_y\Phi_k$. Line **857**: $D:=t\,\dfrac{d}{dt}$, the **Euler operator**, which
"multiplies the coefficient of $t^m$ by $m$". Line **1826**, immediately before the
c-table: "Write $\rho_j=D_t^j\rho$." Since $D_t^jt^m=m^jt^m$ the closed form is exact:

$$\rho_j(t)=\bigl(t-3^j s_3^2 t^3+5^j s_5^2 t^5\bigr)/V.$$

Anchor for this formula, independent of the dispute: it yields $\rho_0(i)=i$ exactly with
$\lvert\rho_0(i)\rvert=1$, reproducing the known $\lvert r\rvert=1$ divergence at
$\theta=\pi/2$ that the $b$-series hits.

**The code uses plain $\partial_t^j\rho$, evaluated at the wrong point.** The agent's
reported pair at $\theta=\pi$ ($t=-1$) was $\rho_1=+0.5942892309$,
$\rho_2=-0.5737406266$; **both** match plain derivatives at $t=+1$ to eleven digits
($\rho'(1)=0.59428923089$, $\rho''(1)=-0.57374062658$). Two compounding errors: wrong
operator and wrong evaluation point.

|$j$|Euler $\rho_j$ at $t=\pm1$|code, plain $\partial_t^j\rho$|discrepancy|
|---|---|---|---|
|$1$|$\pm0.59428923089$|$+0.59428923089$|coincides at $t{=}+1$, **sign-flipped** at $t{=}-1$|
|$2$|$\pm0.020548604311$|$\mp0.57374062658$|code **$27.92\times$ too large**|
|$3$|$\mp1.6011717695$|$-0.47423912061$|code $3.376\times$ too small, sign-flipped at $t{=}-1$|

**This accounts for every observation.** It sits above the derivative and norm layers,
which test clean; it lives inside $\Phi_3$ exactly where nonnegativity localized it; it
feeds the $\partial_x\partial_y$-activated rows, hence leg2 and not leg1; and it is shared
by **both** c-table variants, which is precisely why swapping the entire table moved
avg$(S)$ by only $4.73\%$. It also explains the $12\times$ asymmetry between $\theta=0$ and
$\theta=\pi$ at equal $\lvert\rho\rvert$: the code's $\rho_j$ do not sign-flip as Euler
derivatives must.

**The paper's printed c-table is vindicated and the row dispute is CLOSED in the paper's
favour.** $D_t$ carries the $t$, so the printed $-3t(\rho_1+\rho_2)$ is right; the agent's
derived $-3\rho_1-3t\rho_2$ is missing a $t$ on the $\rho_1$ term — the exact signature of
applying $\partial_t$ where $D_t=t\,\partial_t$ was required. The derived-rows variant is
explained and dead; the algebra was sound, the operator was not. This is the **second**
defect located in the `rho_j`/`rho_eval` layer, which is exactly why the voided
finite-difference audit could never have caught it.

Unaffected: the Parseval survivor consumes $[t^{m-b}]\rho^a$, i.e. $\rho$ itself and never
$\rho_j$, so $\lVert D^3H\rVert^2\ge0.77736289857822393$ stands. **On this evidence the
$2.94\times$ is ours and no paper-side claim will be made.** Pending measurement: re-run
the 49-point sweep with the Euler $\rho_j$ and the paper's printed table, reporting
avg(leg1), avg(leg2), avg$(S)$, the ratio, and $S(0),S(\pi/2),S(\pi)$ — a miss to be
reported as readily as a hit, and the $(47)/(48)$ bridges and assembly-vs-jet anchors
re-established on the corrected $\rho_j$.

### Resolution: the $3\times$ was ours, and the paper stands

**Euler run, definitive.** Paper-printed c-table verbatim plus Euler
$\rho_j=(t-3^js_3^2t^3+5^js_5^2t^5)/V$, 49-point uniform trapezoid:

|quantity|defected|corrected (Euler)|paper|
|---|---|---|---|
|avg(leg1)|$6.404064$|$1.939156$|—|
|avg(leg2)|$93.328494$|$30.359613$|—|
|avg$(S)$|$379.718041$|$\mathbf{123.377609}$|$\le126.80385221$|
|ratio to paper|$2.9945308$|$\mathbf{0.97297998}$|—|
|$S(0)$, $S(\pi)$|$141.649$, $1711.417$|$171.6916$, $171.6916$|—|

Our pipeline now reproduces the paper's certified intermediate **from inside**, with
$2.7\%$ slack — what a certified upper bound should look like. **Hypothesis (iii) is DEAD,
not merely unadvanced**: both the paper's intermediate and its conclusion stand.

**Owner arithmetic correction.** The agent reported $B_3=14.2202$; that is wrong.
$(\pi^2/6)\times123.377609 = 202.948032$, so $B_3 = 14.2459830$ (the agent's value squares
to $202.214088$). The conclusion is unaffected — $14.2459830 < 14.4424366$, the paper's
certified boundary — but the digits were about to enter the index wrong.

**Headline strengthens.** Break-even is $\lVert D^3H\rVert\le53.9235230$. The old
$B_3=24.9922200$ was the *defected* value and its $2.1576124\times$ slack is **RETRACTED**.
Corrected: $B_3 = 14.2459830$, slack $\mathbf{3.7851739\times}$. Our own route now lands
below the paper's boundary rather than merely consuming its cited tail, so the
`CITED-DEPENDENCY` tail $4.57569\times10^{-6}$ is now **independently corroborated** by our
pipeline — at `COMPUTATIONAL-EVIDENCE` grade only, since the $\theta$-quadrature is
uncertified float. It may not be labelled stronger.

**Independent structural confirmation, from the paper's own averaging step.** Every
exponent in $\rho_j=(t-3^js_3^2t^3+5^js_5^2t^5)/V$ is odd, so $\rho_j(-t)=-\rho_j(t)$ for
**every** $j$. Since $C_{-r}(x,y)=C_r(x,-y)$ and product standard Gaussian $\mu$ is
symmetric under $y\mapsto-y$, flipping every $\rho_j$ leaves both $L^2(\mu)$ norms
invariant, hence $S(\theta+\pi)=S(\theta)$: **$S$ is $\pi$-periodic**. Line 1861 replaces
$(1/2\pi)\int_0^{2\pi}$ by $(1/\pi)\int_0^{\pi}$ — a step **valid only if $S$ is
$\pi$-periodic**. So the paper's own averaging is a structural proof that $\rho_j$ must be
Euler derivatives, and plain derivatives cannot satisfy it ($\rho'$ even, $\rho''$ odd,
$\rho'''$ even — mixed parity). Measurement matches exactly: Euler gives
$S(0)=S(\pi)=171.6916$; the defected run gave $141.649$ vs $1711.417$, violating the
symmetry by $12\times$. This is a second wholly independent line of evidence, and it
reclassifies the $12\times$ asymmetry from curiosity to symmetry violation.

**All three closures LANDED, and the numbers are owner-verified from the frozen arrays
rather than from the agent's report.** (1) **Reproducibility** — the Euler c-table is now
committed inside `d3h_sweep49_asm.py` (the runtime patch removed; the defected
`ctab_PLAIN_DERIVS` retained behind a DO-NOT-USE banner), and a clean-process rerun gives
avg$(S)=123.475199$, $123.377610$, $123.377609$ at $13$, $25$, $49$ nodes — the $25\to49$
increment is $10^{-6}$, the same spectral class established earlier. (2) **The flagged
anchor** — a fresh independent mpmath-40-dps evaluator (`euler_fresh_anchor.py`) puts
assembly-vs-jet at worst relative $6.6\times10^{-15}$, the class of the pre-Euler
$3.7\times10^{-14}$; the agent's earlier $1.2\times10^{-3}$ was its own probe truncation and
is retracted. (3) **$\pi$-periodicity is now a standing assertion** of the committed sweep,
holding to $3.96\times10^{-15}$. The $(47)/(48)$ bridges re-verified on Euler bits
($1.7\times10^{-16}$) and the norm layer is unaffected ($224.9967$ vs $225$).

**Owner independent recomputation from `leg2euler_*.npy`** — every figure reproduced, and
one check the agent did not run:

|quantity|owner recompute|agent reported|
|---|---|---|
|avg(leg1)|$1.939156$|$1.939156$|
|avg(leg2)|$30.359613$|$30.359613$|
|avg$(S)$|$123.377609$|$123.377609$|
|$B_3$|$14.2459830$|$14.2459830$|
|slack to break-even|$3.7851739\times$|$3.7851739\times$|
|ratio to $126.80385221$|$0.97297997$|$0.972980$|

**Stronger symmetry check, owner-derived.** $\pi$-periodicity alone constrains only
$S(0)=S(\pi)$ — a single equality. But $S$ is also *even*: $\theta\mapsto-\theta$
conjugates $t=e^{i\theta}$, the $\rho_j$ have real coefficients, and the $L^2(\mu)$ norms
are conjugation-invariant. Even **and** $\pi$-periodic forces the full reflection
$S(\theta)=S(\pi-\theta)$, i.e. symmetry about $\theta=\pi/2$ — $24$ independent equalities
across the $49$-node grid rather than one. Measured on the frozen array: maximum relative
deviation $7.936\times10^{-15}$, with $S(0)=S(\pi)$ to $9.93\times10^{-16}$. Consistently,
$\min S = 1.266545$ falls **exactly** at $\theta=\pi/2$, the symmetry point, as a smooth
symmetric function's critical point must; $\max S = 310.806960$ is attained at an interior
pair, which is why avg$(S)=123.38$ sits well above the endpoint/centre midpoint of $86.5$.

**Label.** The row stands at `COMPUTATIONAL-EVIDENCE`: the $\theta$-quadrature is
uncertified float, so despite spectral self-convergence and $10^{-15}$ structural
agreement it may **not** be labelled `MACHINE-VERIFIED`. Promotion requires certified
interval widths end to end through the quadrature.

**Retracted as properties of a defected pipeline**: $86.0457$, $18.2592$, $362.44$,
avg$(S)=379.718041$ and $397.667200$, $B_3=24.992220$, slack $2.1576\times$, and the
entire $(0,1)$/$(1,1)$ c-table dispute. **Retained**:
$\lVert D^3H\rVert^2\ge0.77736289857822393$ over the exact $m\le251$ window, which never
touched $\rho_j$.

**Standing scope on this whole thread.** These are findings about in-repo code, **not
about the source paper**, and nothing about the paper's certificate is asserted.
Note what this vindicates: the reason a second
independent route was demanded is exactly that a converged single implementation can
be a converged *wrong* implementation, and that is what happened. For the record,
$3.17$ sits within $1\%$ of $\pi$ ($0.93\%$), $\sqrt{10}$ ($0.27\%$) and $19/6$
($0.13\%$) simultaneously, so it carries **zero** evidential weight; two prior owner
guesses of exactly this shape were falsified. Remaining per-target gaps are stated in
each row above and in each target README.

## Resolved threads (2026-08-30)

Two apparent problems in the literature were raised by agents, escalated before
publication, and **both dissolved on investigation**. Recording them, because a
correctly-killed false positive is calibration evidence for everything else here.

1. **`mm3/` — Sun-56 versus the 55-paper: decomposition mismatch, not
   contradiction.** Sun's 13-addition left circuit is real (verified twice: his
   shipped `verify.py` SIDES and his printed §4; 729 identities, 0 failures) but
   computes *his own* $U$-map. His decomposition is not a reorientation of
   `cr58_cn122`-as-printed — factor-column multiset equality fails, and still
   fails up to per-column signs, under the $T$ involution, and entry-transposed,
   with **zero** common product triples against $\sigma$ applied to the
   55-paper's factors. Both decompose $M\langle3,3,3\rangle$ from different
   points in scheme space.
2. **`mm3/` — the 59-versus-57 flag: the paper is right, the agent's recount was
   wrong.** A hand tally of Mårtensson–Stankovski Wagner–Stapleton's printed
   Table 2 came to 57 (and once 58) against their claimed 59. Convention was
   checked first and matched exactly — their footnote 1 charges additions *and*
   subtractions, and the table captions state negation is avoided by rearranging
   terms, which is the model the recount used. A dict-driven machine recount
   then gave left $15$ + right $15$ + output $29 = \mathbf{59}$, an exact match
   with the paper; the earlier figures were the agent's own tally slips
   (undercounted multi-term operands `M15R`/`M8R`/`M22R`, and a C-assembly row
   count of 19 against an actual 20). **Retracted inline** per rule 5.
   Corroborated two further ways: a certified lower bound on their Table-3
   artifact ($d(U)=13$, $d(V)=12$, $d(W)=14$, all floor-impossible over
   $132/81/672$ states, giving a certified sum LB of $56\le59$), and a
   constructive transposition check (reverse-accumulation transposed circuit is
   exactly 15 additions, $15+14=29$, matching the output recount). **Record
   history unchanged: $60\to59\to58\to56\to55$, no underclaim by any author.**
   Frozen at `mm3/campaigns/2026-08-30T012035Z_a7ea8e8e_cbdf8e63fa94/`.

3. **`delcap/` — Fertonani–Duman Table II: valid upper bounds, not an error.**
   An apparent discrepancy ($f(3,2)=1.48$ published against a certified
   $1.46978\ldots$) was measured across 12 tabulated pairs: the offset is
   *irregular*, $+0.0006$ to $+0.0127$, never uniform, and all 12 published
   2-dp values sit above the certified intervals. Mechanism reproducing all 12:
   their stated round-up rule applied to a not-fully-converged Blahut–Arimoto
   **dual**, which approaches capacity from above. No literature correction
   warranted.

**Supersedes** the earlier note in `mm3/`'s results row that the 60- and
59-rungs were "pending artifact availability": the full predecessor ladder is
now machine-checked over $\mathbb Z$ in one fixed convention —
Stapleton-60 (729/729, ternary, 2000-trial randomized end-to-end PASS,
printed-structure recount exactly 60 = 12 pre-gates + 20 inline + 9 v-gates +
19 C-assembly), MWS-59 (729/729 on both their Table 3 artifact and a
reconstruction of printed Table 2), Perminov-58, Sun-56, and the 55-paper.

## Adoption rules

- A row's status moves only with committed evidence (campaign artifacts or a
  proof certificate accepted by a pinned checker) **and** a matching
  `PROGRESS.md` entry.
- Retractions are recorded inline, never silently removed.
- `BENCHMARK` = harness/protocol phase; `ACTIVE` = campaigns running;
  `BLOCKED`/`PAUSED` carry the reason in the target README.
- A finite-$n$ verification is reported as a finite theorem with its $N$; it
  never occupies the headline row as an asymptotic claim.
- **Refutations of published claims never enter this table from a subagent
  report alone.** Owning agents are instructed to escalate to the owner
  session, which re-checks independently before any row moves.

## `delcap/` q=4 correction — 24/24 orbit-total certificates complete, 2026-08-31

This append-only result supersedes the earlier **COMPUTE PENDING** q=4
status. The frozen `c9c3a842…` runner completed the preregistered
`q=4`, `n=5..10`,
`d in {1/2,1/5,1/10,1/20}` grid. All 24 rows are `CERTIFIED`; all 24
select the exact `ba_total_orbit_mass` dual; and all 24 conservative
intervals strictly improve both Tavakoli-Nguyen-Bose finite-$n$ endpoints
(`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`,
**MACHINE-VERIFIED**).

Outward width upper bounds range from
`6.165071504436472e-09` at `(5,1/10)` to
`2.3125135592425783e-05` at `(8,1/2)` bits/symbol. Even the tightest
published-margin cell, `(5,1/20)`, retains positive serialized Arb
lower endpoints of approximately `0.0007402311927961205` on the lower
side and `0.08766315500354434` on the upper side. The exact Arb balls
and exact outward binary rationals—not these display decimals—are the
certificates.

The 24 rows record 11,182,080 direct full-input dual evaluations,
3,585,643,216 positive-$W$/positive-$D$ dual term checks,
1,792,821,608 direct primal conditional entries, 16,773,120 expanded
input generator checks, and 44,728,272 expanded output generator
checks, with zero full-word/representative overlap failures. Post-run
resume validation passed with 24 rows; owner checks verified all 24
newline-inclusive row hashes and the exact CSV projection
(**MACHINE-VERIFIED**). Runtime was `4651.466 s` wall / `4649.094 s`
CPU, maximum row-recorded RSS `1,441,169,408` bytes, with no resource
stop (**COMPUTATIONAL-EVIDENCE**).

Canonical report:
[`delcap/.../report.md`](delcap/campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/report.md).
Final checksum-ledger SHA-256:
`d2660432da9a4d79fb3d9c0e134a249cef6c546fe199f2a3bce672dce905eb30`.
The old pre-production manifest/ledger and post-hold amendments remain
preserved; additive final artifacts close the campaign.

**Scope:** this is a complete finite theorem only for the registered
24-cell q=4 box. It makes no q=3, q>=5, n>=11, other-$d$, or
asymptotic-capacity claim. Named next campaign:
**q=3,n=11 total-output-orbit-mass extension**, requiring its own
prospective freeze.

## Final current-status override for the seven targets, 2026-08-31

This is the controlling status summary for the historical table at lines
515–523 and its earlier override. Detailed claims remain in the target
sections and frozen campaign artifacts.

- **`omega/`: PARTIAL / FAILURE TO CERTIFY (V1 closed).** The slope replay now
  reproduces the frozen raw anchor under the corrected v11 interval core:
  replay `2.3715538358350807` vs frozen `2.3715538358350803`, residual
  `4.440892098500626e-16` — **exactly one ULP**, inside the `5e-13` V1 gate,
  with all seven R branches, `R_sum=2.8170035674609757` and
  `M_low=2.0942543887102634` float-equal to frozen (**MACHINE-VERIFIED**;
  owner-re-derived, and the replay path is `gate_c_slope_pass.run_slope_pass`,
  distinct from the `importlib`-loaded frozen control, so the agreement is
  informative rather than tautological). The requested Lemma-1 penalty
  diagnostic is answered: the lemma1-off shift is `2.029648699330977e-6`, so no
  `~2.9e5` amplification exists and the penalty path never carried a 0.9-scale
  residual. The LP/ANY-y stage ran and returned **FAILURE TO CERTIFY**: the
  rigorous enclosure `[-5.841014555620028e-6, +5.841014555620028e-6]` straddles
  zero at `73.85977508485684x` the `1.5816497000997742e-7` signal (enclosure
  **MACHINE-VERIFIED**, the ratio **INFERENCE**). Gate C stays **PARTIAL**; the
  21-dimensional kernel question is **UNCHANGED and OPEN**; and the guardrail is
  re-verified arithmetically — the gap to the published `2.37155181` is
  `2.0258350805768544e-6`, i.e. `12.808x` the best available signal, so this
  route can never produce a record. The earlier `+0.9008714859433109` residual
  is **retracted, non-reproducible** (`NameError: name 's_const'` before its
  first checkpoint) and survives as **REPORTED** history only.
  A follow-on pre-registered feasibility gate then closed the named next action as
  a **quantified obstruction** (campaign
  `2026-09-01T03:15:00Z_OmegaGateCBnBGate_9d2e4a7c`): the enclosure's L1 width is
  diffuse, not concentrated — `k* = 36` of 45 coordinates are needed to reach 90%
  of the width, the top six carry only `21.8493%`, and 40 coordinates each carry
  at least `1%` (owner-re-derived from `stage1_decomposition.json`; the widths sum
  to `113.66562530936517`, matching `GATE_Z.live_l1` exactly). At `~6.2`
  bisections per binding coordinate that implies `2^223.2 ~ 1.548e67` subboxes, so
  **branch-and-bound over D at this enclosure technology is infeasible at this
  cost** — a first-class negative, with the pre-registered threshold honored and
  zero subbox compute spent. Every R branch's minimizer identity is undecided, the
  widest candidate span being `R_glob[0]` at `3.2123059392596964e-05` = `203.1x`
  the signal. The obstruction binds this technology over this D only; a different
  certified-gradient technology and analytic branch-pinning remain open.
- **`kg/`: PARTIAL / FAILURE TO CERTIFY the full band; the retraction stands.**
  Exact direct-D.4 intervals certify the continuous prefix through
  `3.107728208020454953573248` and a terminal overlap through `3.5`; the
  contiguous interior gap remains open (**MACHINE-VERIFIED** finite-domain
  statements). The proposed panel-wise `K_o` mechanism was already present and is
  retracted as new. **Band `[4.083, 6.0]` re-run under the corrected 2-D even-box
  envelope, 2026-09-01** (campaign
  `20260831T231500Z_KgEvenBoxRerun_campA_corr_env`, commit `8d2cb26`): tiles 1-2
  PASS with certified margins `+1.04499016763434518e-6` and
  `+4.71904853878616890e-6` (radii `~1e-36`, Arb 256), tile 3
  `[4.1485893120, 4.1817780265]` OPEN at open-leaf margin
  `-1.76291787028287693e-5` with `95/251` leaves certified at depth 14 and the
  `131072` panel cap reached, so by the pre-registered stop rule the band is
  **FAILURE TO CERTIFY** and tiles 4-49 were not run. The corrected instrument
  does not reproduce the retracted `kgcloseB3` PASSes, so that retraction
  **STANDS** and no new PASS is earned. The failing cell is a disk-ring where
  `r_o -> 0` and the binding cost is panel-sup `|e|`-hinge slack: it is a failure
  of our envelope at our caps, **not** a lower bound on `J - d` and **not** a
  paper refutation. Startup battery 16/16 (**MACHINE-VERIFIED**) including
  `e(0) = sqrt(3/2)` admitted and tight (slack `1.87e-96`), the cap-at-1
  counterfactual rejected with clipping quantified at `0.2247448713915890`, the
  3-D-inflation counterfactual confirming the old form was looser, the
  byte-verified `d3h` certificate, and a three-way validation of the margin
  predicate itself.
- **`rs-pe3d/`: corrected finite-slice theorem.** On `t=(1,1,1)`, the union
  support window is empty exactly below `min_i d(C_i)` and first becomes
  nonempty at that weight. This is not an individual-support biconditional.
  The kernel formula, coset argument, minimum-word wording, and
  non-uniqueness of triple-sum representations are corrected in the target
  README. Algebraic deduction is a human proof anchored by exact machine
  checks; separate-commit preregistration is only **REPORTED** chronology.
  **Corrected-form theorem re-derived with a commit-anchored pre-statement,
  2026-09-01** (campaign `2026-09-01T02-33-35Z_mech2_mechanismCorrected`,
  pre-statement committed pre-compute at `29b0dfe`, closing the predecessor's
  **REPORTED**-grade chronology gap). Theorem T, corrected form: (1) for EVERY
  support with `|S| < d`, `V cap F^S = {0}`; (2) there EXISTS a support of size
  `d` carrying a codeword; (3) NO per-support claim at `|S| >= d`. Census
  semantics: the union window is empty iff `w < d`, first nonempty exactly at
  `d`. At `t=(1,1,1)`, `d = min_i s_i`, so H-GATE's min-form is a **theorem** on
  that slice and the `(4,5,7)` class sits inside it. Controls 16/16 exact `F_q`.
  **Owner independent verification**, by a different construction
  (`V = ker(H_0 (x) H_1 (x) H_2)` from GRS parity checks, emptiness via
  `dim(V cap F^S) = |S| - rank(pi[:,S])`, versus the agent's rref of stacked line
  classes): the L5 separator argument audited line by line and valid; PN4
  `dim V = 68`, the `w <= 3` census `457,450` supports with `0` nonzero, and
  `35/35` direction-0 minimum lines all reproduced exactly. The owner also
  extended past the campaign's controls, whose emptiness blocks were all
  `t=(1,1,1)`: parts 1 and 2 both hold at `(13,(5,4,4),(2,2,2))`,
  `(13,(4,4,4),(2,1,3))` and `(31,(5,5,6),(3,2,2))` with random `Lambda`, so the
  machine anchoring is no longer `t=1`-only. Independently, two kernels
  (V-basis and parity) agree `15/15` on the 13 frozen gate-B tallies plus
  `PN6 3486` and `PP2 917`.
  **New owner finding, sharper than the theorem's clause (3).** Exhaustive exact
  censuses at weight exactly `d` return candidate sets that are SET-EQUAL to the
  minimum lines in argmin-`d` directions, with zero non-line dependencies: PN4
  `35 = 35` over all `C(140,4) = 15,329,615` supports, `(17,(4,4,4))` `48 = 48`
  (all three directions minimum), `(41,(4,4,5))` `40 = 40` (directions 0 and 1
  only, direction 2 correctly absent at `s_2 = 5 > d`). Status: **MACHINE-VERIFIED**
  as a finite statement about those three instances; **CONJECTURE** ("H-MINLINE")
  as a general claim — three instances are three instances. If proved it would
  upgrade clause (2) from existence to exact identification. Campaign W (PN4
  weight-exactly-4 window) is frozen mid-flight with its pre-statement at
  `8114a8f`, validation and a demonstrated kill-restart resume complete, and a
  finishing agent holds the resume protocol.
- **`oct-rank/`: Route A CLOSED-NEGATIVE, provenance-clean; Route F stays OPEN at
  `{13,14}`.** The Strassen 1983 primary (LAA 52/53:645-685, DOI
  `10.1016/0024-3795(83)80041-X`) was obtained and read FIRST-HAND, closing the
  predecessor's `UNESTABLISHED` provenance flag; owner-verified from the retained
  text layer and matching publisher metadata. Theorem 4.1 with the proof's `(4.2)`
  step carries a **one-half** factor, so the family yields
  `rank >= n + (1/2) rank(comm)`. Since `comm` is `8x8`, `rank(comm) <= 8`
  trivially, capping the entire theorem family at `8 + 4 = 12 < 13`, the standing
  certified floor — the closure needs no further hypothesis. The naked
  `n + rank(comm) = 16` reading is dead twice over: it has **no primary source**
  (Bläser is uncited in a 1983 paper, "commutator" is never named there), and it
  is **arithmetically false** — on the `tau` control it would assert `4+4=8`
  against a true rank of `7`, and a lower bound cannot exceed the true rank.
  Corollary 4.2 is inapplicable to `O` because its associativity clause fails.
  Controls: `tau` gives `4 + 4/2 = 6 <= 7` (**MACHINE-VERIFIED**), `tau-box-s`
  gives `12 <= 13`. An owner objection that the universal `rank(M) = 8` wording
  outran its certificate (C3 verifies `M^2 = -4 detGram/N(u)^3 I` at six named
  triples) is **RETRACTED**: the universal statement is a valid HUMAN PROOF in
  `campaigns/2026-08-31T08:02:18Z_routeAF/commutator_derivation.txt`, audited
  line by line by the owner — the anchor identities are bilinear, so basis
  checks plus bilinearity give them for all arguments, and the orthogonal
  decomposition yields `D^2 = -4 U detGram I`, nonzero exactly on independent
  real triples. The six C3 points are cross-check anchors, not the quantifier's
  basis. The verdict is independent of this either way, resting on the trivial
  `rank <= 8` cap. Numerical
  rank-13 misses remain **COMPUTATIONAL-EVIDENCE**, never nonexistence. Two
  untracked drafts were deleted without approval and then restored byte-exactly;
  neither is evidence.
  **Route F attempted and closed as FAILURE TO CERTIFY, 2026-09-01** (campaign
  `2026-09-01T04:30:00Z_routeF_kraw`, pre-statement `739e63f` before compute).
  The agent first caught a rule-14 defect in the inherited plan — a
  gradient-system Krawczyk certifies critical points of the residual, and no
  interval method can prove a residual is exactly zero — and replaced it with an
  existence-carrying square-slice instrument, counted explicitly: `247`
  parameters (`13 x 19`) minus `55` frozen gauge coordinates equals `192`
  equations (`8 x 8 x 3`), owner-checked. Outcome: the pre-registered
  `5 x 5000`-nfev polish bottoms at relative `8.816619e-06`, containment `0/192`
  at ALL 15 fixed rungs, and certified no-root exclusion is COMPLETE (`192/192`)
  for every `rho <= 1e-6` around both polished candidates — so no exact on-slice
  rank-13 decomposition exists within `1e-6` of either seed
  (**MACHINE-VERIFIED**, exact `fmpq` with outward-`arb` re-verification). This is
  explicitly **not** a rank-13 impossibility: off-slice roots are not excluded and
  nothing global is claimed. The instrument discriminates in both directions,
  owner-replayed from frozen bytes (15/15 PASS, 42.6 s): `tau`-7 containment
  climbs `1/48 -> 45/48 -> 48/48` as `rho` tightens (**accepts** a known-true
  case), while the globally infeasible `tau`-6 plant, a genuine near-miss
  corrupted-target candidate, and both main candidates are all **rejected** with
  certified exclusion. `rank((L_1, L_i, L_j))` stays **OPEN in {13,14}**: absence
  of a 13-witness is not evidence for 14, and absence of an impossibility argument
  is not evidence for 13. Four defects disclosed, all caught by pre-registered
  anchors or startup halts BEFORE certification data existed — including a
  target-convention error (raw slices versus blockdiag conjugate) on run 1.
- **`mceliece/`: exact toy direct-route census complete, paper-scale route
  OPEN.** All 24 registered GF(32) cells pass P1 and fail P2–P4, so the toy
  `genericity_pass` is false (**MACHINE-VERIFIED**). The population lies
  outside current ePrint 1786 Table 1 and Assumption 1 and proves nothing
  about a NIST cell. Parent/successor wall times are
  **COMPUTATIONAL-EVIDENCE**.
  **`m=12` REACHED at a real NIST cell, 2026-09-01** (campaign
  `2026-09-01T03-28-10Z_M12ALPHA_DIVONLY`, pre-statement `7fb44ae` before
  compute, checksum ledger 17/17): instance `(12, 3488, 64, seed 16384)` — the
  `mceliece348864` parameter set, where `m=12` had previously been UNREACHED and
  frozen as such. `alpha` is **MEASURED**, not derived: the unmodified `m<=11`
  instrument evaluated all `n = 3488` support points with zero failures, and the
  ADDENDUM-2 derive fork (**CITED-DEPENDENCY**) was never reached because the
  measured projection came in far under the registered `10800` CPU-s budget.
  Identity of the instrument is licensed by 3/3 byte-equal anchors against the
  frozen gate-A records at `(11,2048,48,6211)`, `(10,1024,40,5113)` and
  `(6,64,3,1387)`. The `delta` chain is exact (**MACHINE-VERIFIED**): exact
  division `3488/3488`, Lagrange unit checks `3488/3488`, structural off-diagonal
  vanishing (definitional for a Lagrange basis, with 200 redundant probes), and
  the degree gate closing exactly at `max deg f = D = 3359 = n - 2t - 1`, with
  `k = 2720` owner-cross-checked. Three counterfactual plants were REJECTED by
  three distinct mechanisms. **Scope, not inflated:** a BOUNDED FINITE statement
  about ONE instance — `t = 48` and `t = 96` at `m = 12` remain UNREACHED, no `G`
  population is exhausted, this is not an `m=12`-complete result and not a NIST
  attack-cost verification, and **Apon's section-3.6 hole remains BOUNDED, not
  closed**. Six defects disclosed, three of them amendments committed before the
  verdict compute.
- **`mm3/`: expanded named landscape complete.** Exact censuses put the
  `mws59` and `stapleton60` lower-bound minima at 56 and 58, with zero
  admitted triples at most 54 (**MACHINE-VERIFIED**). `paper55` is the only
  verified 55-addition circuit in the expanded swept set, but is not the
  unique LB-55 data-triple minimizer; `sun56` remains LB 55 / UB 56. No
  universal no-54 claim follows.
- **`delcap/`: q=3 correction and q=4 extension complete.** The corrected
  total-output-orbit-mass construction certifies all 20 registered q=3 rows
  and all 24 registered q=4 rows, each strictly inside both published
  finite-$n$ endpoints (**MACHINE-VERIFIED**). These are finite theorems, not
  asymptotic-capacity claims.

The forward arXiv/ECCC delta found no collision with these targets. The
matching-vector candidate was the already-audited TR26-156 under a second
identifier and had no load-bearing target-scale finite gate; **no decorative
eighth target was adopted**. The Monday 2026-08-31 window (owner-re-derived:
`cs.CC` 16 items, max published `2026-08-31T12:04:08Z`; `cs.IT` 37 items, max
`2026-08-31T11:03:56Z`; ECCC zero new reports, ceiling TR26-162 unchanged) is
likewise **no collision, no gate affected, no target opened**. That sweep also
exposed a third coverage defect in this repository's own instrument — one
`max published` printed for a two-category sweep hid three `cs.IT` items, all
triaging clean — recorded as rule 17d and in `PROGRESS.md`; no result changes.

## `delcap/` q=3,n=11 extension — 4/4 orbit-total certificates complete, 2026-08-31

This append-only result closes the cell that the q=4 campaign excluded and
named as its successor, and supersedes the immediately preceding delcap
bullet only by adding `q=3,n=11`; no earlier delcap claim is withdrawn. The
frozen `08ecbc49…` runner completed the preregistered `q=3`,
`n=11`, `d in {1/2,1/5,1/10,1/20}` box. All four rows are
`CERTIFIED`; all four select the exact `ba_total_orbit_mass` dual; and all
four conservative intervals strictly improve both Tavakoli-Nguyen-Bose
finite-$n$ endpoints
(`CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`,
**MACHINE-VERIFIED**).

| $d$ | certified interval (outward) | width upper |
|:---:|:---|---:|
| 1/2 | `[0.4136415126778377, 0.4136475665528263]` | `6.053874988461827e-06` |
| 1/5 | `[0.9165946848818924, 0.9166017532024916]` | `7.068320598920295e-06` |
| 1/10 | `[1.2046904010391062, 1.2046925848595933]` | `2.1838204867062363e-06` |
| 1/20 | `[1.3817159992481873, 1.3817177042129642]` | `1.7049647766059455e-06` |

The tightest published margin is at $d=1/20$, where the archived
Arb lower-margin ball has lower endpoint approximately
`0.00129169637501813941901314966788` and the upper-margin ball
approximately `0.12399667147213431625505939351440`. The exact Arb balls and
exact outward binary rationals—not these display decimals—are the
certificates.

All four rows reproduce the exact census tuple
`(14884, 22450, 6148309, 3573542, 9721851)`, orbit-size sums
$177147=3^{11}$ and $265720=\sum_{k=0}^{11}3^k$, and the independent
Burnside counts. The rows record 1,417,176 direct full-input dual
evaluations, 634,415,384 positive-$W$/positive-$D$ dual term checks,
316,716,419 direct primal conditional entries, 2,125,764 expanded input
generator checks, and 9,565,920 expanded output generator checks, with zero
full-word/representative overlap failures. Post-run resume validation passed
with four rows; owner checks re-derived all four newline-inclusive row
hashes, re-verified each outward decimal against its archived exact binary
rational, re-checked `width <= 1/500`, and confirmed the exact CSV
projection (**MACHINE-VERIFIED**). Runtime was `1151.072 s` wall /
`1150.401 s` CPU, maximum row-recorded RSS `929,677,312` bytes, with no
resource stop (**COMPUTATIONAL-EVIDENCE**).

Canonical report:
[`delcap/.../report.md`](delcap/campaigns/2026-08-31T20-34-52Z_04ed25bd-455f-4814-8288-5ee0106e0db8_q3-n11-total-output-orbit-mass/report.md).
Final checksum-ledger SHA-256:
`0ba5629a13854ba4eb1cdff32aad83fa0a385ec9a3da6b33bef5f94ac28b7761`.
The pre-statement, static-guard evidence, and static ledgers remain
byte-preserved; additive final artifacts close the campaign.

**Scope:** this is a complete finite theorem only for the registered
four-cell `q=3,n=11` box. It makes no q>=4, `q=3` with `n!=11`,
other-$d$, or asymptotic-capacity claim. Named next campaign:
**q=3,n=12 total-output-orbit-mass extension** (44,530 input orbits over
531,441 words; 66,980 output orbits over 797,161 words), requiring its own
prospective freeze.

## `rs-pe3d/` H-MINLINE falsified, and theorem T-DGE certified with a sharp hypothesis, 2026-09-01

**Two separate verdicts, deliberately not merged.** The conjecture H-MINLINE is
**FALSE**; the corrected statement it was reaching for is now a proved theorem
with a hypothesis shown to be sharp.

**H-MINLINE is FALSIFIED (MACHINE-VERIFIED).** At exactly the reported slice
$q=13$, $s=(2,2,4)$, $t=(1,1,1)$ ($N=16$, $d=2$, $\dim V=13$,
$\dim V^{\perp}=3$) the weight-$d$ carriers number **24 = 16 minimum lines + 8
off-line diagonal pairs**. The off-line pair is exhibited, not merely counted:
the support $\{(0,0,0),(1,1,0)\}$ carries an exact word built as a sum of two
axis-constant words that cancels exactly on one flat, certified in $V$ by two
space-exact dual routes and independently by $\operatorname{rank}(G+[w]) =
\operatorname{rank}(G)$. Controls pass in both directions: unique-minimum
$s=(2,3,4)$ gives 12/12 lines with 0 off-line, tied $s=(3,3,4)$ ($d=3$) gives
24/24 with 0 off-line, a known line is accepted at dimension 1, and both
planted noncarriers are rejected at dimension 0. One root-cause instrument
defect was found by the agent's own cross-check and disclosed: the
preregistered route-A space was **not** $V$'s annihilator, coincidentally
agreeing on every $\le 4$-weight support of all three instances; the
falsification stands on rebuilt space-exact duals.

**Theorem T-DGE (HUMAN-AUDITED proof; counts MACHINE-VERIFIED).** Let
$H=H_1\otimes\cdots\otimes H_n$ over any field, every $H_i$ with nonzero
columns and spark $d_i$, and $d=\min_i d_i\ge 2$. Then
$\operatorname{spark}(H)=d$; and **if at most one $d_i$ equals 2**, a set of
exactly $d$ columns of $H$ is dependent iff it is an axis fiber — all index
tuples agree outside one coordinate $i$ with $d_i=d$, and the varying indices
form a size-$d$ circuit of $H_i$. If two or more $d_i$ equal 2 then $d=2$ and
the characterization fails, the dependent pairs being exactly the distinct
tuples that agree componentwise up to parallel class.

**The hypothesis is sharp in both directions.** "$d\ge3$" is strictly too
weak: the theorem also covers $d=2$ with a single spark-2 factor, which is
precisely why $s=(2,3,4)$ and $s=(2,4,4)$ satisfy the characterization while
$s=(2,2,4)$ does not. Sharpness at the other end is machine-checked: every
configuration with two spark-2 factors exhibits at least one off-fiber
diagonal.

**Counts closed-form and exhaustively verified.** At $d=2$ the dependent-pair
count is exactly $(\prod_i A_i - N)/2$ with $A_i=\sum_k m_{i,k}^2$ over the
parallel-class multiplicities of $H_i$, and the fiber subcount is
$\sum_i (\prod_{j\ne i} n_j)\sum_k \binom{m_{i,k}}{2}$; in the fiber regime
the carrier count is $\sum_{i:\,d_i=d} (\prod_{j\ne i} n_j)\cdot
\#\{\text{size-}d\text{ circuits of }H_i\}$. For the diagonal GRS instances
$d_i=s_i-t_i+1$, classes are trivial unless $d_i\le2$, and $d_i=1$ means that
factor is the zero map so every product column is a carrier.

**The "minimum line" phrasing is a $t_i=1$ artifact.** For $t_i>1$ the
minimum-weight supports of $C_i$ are proper $d$-subsets of a line, so the line
wording fails even where the theorem holds: at $s=(4,4,4)$, $t=(2,1,1)$ all
**64** carriers are axis-0 fibers on 3-subsets of a 4-point line, against only
16 full lines. The LINE corollary is therefore stated only at $t_i=1$.

In-run controls for the proof campaign: M1 six structured $\mathrm{GF}(13)$
fiber-regime configurations exact (100/50/192/64/180/90, zero off-fiber), M2
seven diagonal-regime rows against the pair formula with measured
multiplicities, M3 21/21 randomized non-MDS $\mathrm{GF}(7)$ configurations
with zero violations and $\operatorname{spark}(H)=\min_i d_i$ throughout, M4
both plants rejected fail-loud with a pristine accept, M5 15/15 census
reconstruction rows, M6 boundary sharpness 6/6. Four instrument defects were
caught by first execution, disclosed pre-freeze, and every affected family
re-run from scratch.

**Owner-independent verification.** All headline counts were reproduced in a
from-scratch $\mathrm{GF}(p)$ Kronecker instrument written independently of the
campaign code: the six structured rows, 21 randomized non-MDS configurations
across 2–4 factors, the pair formula on 33/33 randomized multi-class rows, the
15/15 census closed forms, and the two control plants. The frozen census file
was cross-read and its $\dim V + \dim V^{\perp} = N$, 24 = 16 + 8, 12/12/0 and
24/24/0 rows agree exactly.

| run | gate | verdict | checksum-ledger SHA-256 |
|---|---|---|---|
| `20260901T103743Z_d86a1fa2_f72db722778d` | H-MINLINE | FROZEN-NEGATIVE | `3e636c27cc98b6bd1ad15928c339fb73910861a60696f06484b0bf97b061a778` |
| `20260901T104923Z_fba446bf_5aee86fef866` | H-MINLINE-DGE3 | FROZEN-INCONCLUSIVE | `ea89fe1bc7df19892bf1ac2fac4e8b8671dbb5ef53837256a150dd940461bfdf` |
| `20260901T111136Z_99745efb_610bb61dc5a8` | H-MINLINE-DGE3-PROOF | FROZEN-CERTIFIED | `1fd92c879a416aedbfc80c8b2ce6f75f06d02c8eb6ff0429477c965f40112cb1` |

Canonical proof artifact:
[`rs-pe3d/.../theorem_t_dge.md`](rs-pe3d/campaigns/20260901T111136Z_99745efb_610bb61dc5a8/theorem_t_dge.md).
The middle run is the honest record of the intermediate state: the corrected
theorem was registered with its proof obligations and **not** promoted until
they were discharged in the third run.

**Scope.** T-DGE is a theorem about size-exactly-$d$ dependent sets. It says
nothing at $|S|>d$, and the naive extension there is false: with two factors of
spark 3 and four columns each ($d=3$, no spark-2 factor, so inside the fiber
regime) there are **156** size-4 circuits and **not one of them is a fiber**
(144 with three distinct indices per side, 12 pure diagonals); with five
columns each, 984 size-4 circuits, again none a fiber, while all size-3
circuits remain fibers (32 and 100 respectively). An explicit witness over
$\mathrm{GF}(13)$: $A$ with columns $(1,k)$, $k=1..4$, and $B$ with columns
$(1,0),(-2,1),(1,-2),(0,1)$ give $\sum_k a_k\otimes b_k=0$ on the diagonal.
A necessary condition for such an all-distinct size-$(d+1)$ circuit is
$\rho_A+\rho_B\le d+1$, hence $(d_A-1)+(d_B-1)\le d+1$, so a tied minimum
forces $d\le3$ — verified with zero mismatches on six configurations. The
enabling structural identity is
$\ker(A\otimes B) = (\ker A\otimes F^{n_B}) + (F^{n_A}\otimes \ker B)$, an
equality by dimension count, which places the $|S|>d$ question in the same
"sum of axis-constant tensor spaces" frame as $V$ itself. Named next campaign:
the size-$(d+1)$ classification under that identity.

## `omega/` Gate-C analytic branch-pinning falsified; $\varepsilon_0$ certified bit-constant and dual-coupled, 2026-09-01

**Outcome: FROZEN-NEGATIVE, with a new certified structural fact.** The
registered analytic branch-pinning route is **not** viable: over the $10^{-7}$
`dist[0]` box (a superset of the exact-kernel domain $D$) both pre-registered
order inequalities straddle zero,
$G_1 \in [-2.9602186479446094\times10^{-5},\,
+2.6689872684820378\times10^{-5}]$ and $G_2$ within $4\times10^{-15}$ of it,
a width of $5.63\times10^{-5}$ = **355.9×** the $1.5816497000997742\times10^{-7}$
signal. The $p_{\rm comp}$-versus-penalty minimizer identity is therefore not
pinnable by freezing $\varepsilon$.

**What is now certified (MACHINE-VERIFIED).** $\varepsilon_0$ is
**bit-constant** over the box — $2.1502086942121845\times10^{-6}$ at both
$r=0$ and $r=10^{-7}$ — and fully dual-coupled: the C2 drift boundary moves by
exactly $4\,\text{drift}$ (coupling gap $4.41\times10^{-16}$) and
$R_{\rm glob}[0]$ responds by exactly $-2\,\text{prop}\cdot4\,\text{drift}$
(gap $2.94\times10^{-16}$). So the parent obstruction's $R_{\rm glob}[0]$
straddle is a $p_{\rm comp}$-vs-penalty box-motion race, **not** a Lemma-1
artifact, and $\varepsilon$ leaves the movers list for every successor.

Controls all pass in the pre-registered directions: GATE-Z reproduces both
frozen records (R-branches to $3.33\times10^{-16}$, $\Omega$ to
$8.88\times10^{-16}$), C1 accepts the identity at gap
$5.79\times10^{-13}$, C2 rejects on dual drift, C3 rejects a wrong-index dual
lookup. Attempt 1 aborted at C2 and its traceback is preserved per-attempt;
the authoritative in-directory replay is deterministic with byte-identical
decision fields.

Run `20260901T102048Z_b555441e_d3303018f163`, gate `OmegaGateCBranchPin`,
checksum-ledger SHA-256
`b42449bfe144206d12ebd3f3d592c37a1c0f10cc9a54ecd6b25bbb560d97d1d4`
(9 files). Preregistration commit `754b66d` preceded all compute and
Amendment 1 (`68a65d9`) preceded the first runner execution; the in-directory
prereg copy is the amended file and the manifest records the init-time hash,
both disclosed.

**Scope.** The 21-dimensional kernel question is **UNCHANGED and OPEN** and
gate C stays PARTIAL. This is a route closure plus a structural fact, not a
bound: the 12.808× guardrail still means this line can never produce a record.
Named next gate, frozen in the manifest: a successor must pin a
$p_{\rm comp}$-vs-penalty ordering, or pin a different branch.

## `mm3/` sun56 gap closed upward — LB = UB = 56 for the fixed orientation, 2026-09-01

**A registered gap decided.** `sun56` stood at LB 55 / UB 56, the
all-monomial data triple being the unique LB-55 triple of the frozen sun56
census. It is now **exactly 56** for its fixed orientation, so Sun's published
56-addition scheme is **optimal for its own orientation**. This closes our own
registered gap; it is not a new record against the literature.

**Exact reduction (HUMAN-AUDITED, and the load-bearing step is a proof, not an
assumption).** With $C(U)=13$ and $C(W_{\rm fac})=16$ already exact in the
frozen landscape and output accounted by the transposition principle,
$$\text{total} \;=\; C(U)+C(V)+C(W_{\rm fac})+14 \;=\; 13 + C(V) + 30 ,$$
so a 55-total exists iff $C(V)\le 12$. The closure argument then makes the
question finite: 12 gates covering the 11 required $\tau$ classes leave
**exactly one** non-$\tau$ gate, and that gate's operands must already be
available, hence its class is a signed sum of two elements of
$\text{inputs}\cup\tau(V)$. So $C(V)\le12$ iff some **single** auxiliary
direction $a\notin\tau(V)\cup\text{inputs}$ makes $\tau(V)\cup\{a\}$
schedulable in 12 slots — a universe of exactly **338** candidates, enumerated
by two independent loops and hash-pinned in the pre-statement *before* any
compute (CSV SHA-256 `c10bdeea…`).

**Decision (MACHINE-VERIFIED).** All 338 instances are INFEASIBLE, decided
twice per instance — a complete memoized subset-DFS at $T=12$ and a
freshly written slot-availability CNF solved by kissat 4.0.4 — with
**agreement on every instance** and abort-on-disagreement armed. Owner
re-verification of the frozen checkpoint: 338/338 rows, indices complete
$0..337$, `dfs` and `cnf_sat` both false in every row, zero disagreements, 338
distinct auxiliary vectors with no zero vector and coordinates confined to
$\{-2,-1,0,1,2\}$ exactly as signed pair-sums require. Hence $C(V)\ge13$;
Sun's own 13-gate $V$-circuit (expanded exactly, reproducing all 23 $V$ rows,
13 distinct classes = 11 $\tau$ + 2 auxiliary) gives $C(V)\le13$. Therefore
$C(V)=13$ and the minimum is $13+13+16+14=\mathbf{56}$.

Controls, all in-run and in both directions: ACCEPT — A1 the known-true
paper55 circuit passes 729/729 Brent over $\mathbb Z$ in two independent exact
paths with split recount 13/13/30 = 56; A2 witness circuits reproduce all 23
$U$ and 23 $V$ rows by exact expansion; A3 $\tau$ together with Sun's two
auxiliaries is schedulable at $T=13$ in both instruments; A4 a frozen
landscape row is re-derived from raw factors through the same code path.
REJECT — R1 a one-addition-deleted plant and R2 a one-operand-perturbed plant
are both rejected by exact expansion; R3 floor@11 is infeasible in both
instruments with a fresh DRAT→LRAT certificate accepted by two pinned
checkers (`drat-trim` `s VERIFIED`, `lrat-check` rc 0); R4 the $d$-count
planter fires its assertion. An independent third-path audit re-derived the
338-element universe with a different enumeration loop and re-verified the
certificate with both checkers, 11/11 pass.

**Disclosed defect.** A first CNF encoder omitted negative unit clauses for
non-input classes and admitted 39/338 spurious SATs against the DFS; it was
found in pre-scouting by replaying an extracted SAT model, fixed *before*
preregistration, and never used for any verdict.

Run `20260901T123618Z_3f41e268_0264d9d39939`, gate `sun56-gap-total55`,
**FROZEN-NEGATIVE** for the 55-move, checksum-ledger SHA-256
`f434d19794eea4b8cc5d3e0d00cf28d9970223562237c2edc3d090d360598dbc`
(18 files, 18/18 verified). Measured cost 73.5 s wall of a pre-registered
3.0 CPU-h budget (phases: anchor 0.192 s, controls 1.415 s, scan 71.479 s,
certificates 0.403 s). Canonical report:
[`mm3/.../report.md`](mm3/campaigns/20260901T123618Z_3f41e268_0264d9d39939/report.md).

**Scope.** The statement is: for `sun56` in the $\sigma^0$ monomial
orientation, three-stage linear SLP model, additions including subtractions,
free inputs and free sign changes, output via transposition accounting, **no
circuit with $\le55$ additions exists**. Other orientations of the same
tensor, other decompositions and other models are untouched, and no unbounded
"no such circuit exists" claim is made anywhere. Combined with the frozen
census (campaign `192400Z`, all non-monomial sun56 triples $\ge56$), the
frozen swept set now contains exactly one 55-addition scheme — the published
`paper55` one.

## `rs-pe3d/` size-$(d+1)$ circuits — fiber extension falsified, rank criterion proved, 2026-09-01

**The natural extension of T-DGE to $|S|=d+1$ is FALSE, and the sharp
replacement is proved.** T-DGE settles size-exactly-$d$ dependent sets; this
campaign settles what happens one step up, and the answer is that the fiber
picture does not merely weaken, it collapses.

**Falsification (MACHINE-VERIFIED, explicit witness).** Over
$\mathrm{GF}(13)$ let $A$ have columns $(1,k)$, $k=1,\dots,4$ and $B$ have
columns $(1,0),(-2,1),(1,-2),(0,1)$. Both have spark 3, so $d=3$ and **no
factor has spark 2** — squarely inside T-DGE's fiber regime. Then
$\sum_{k=1}^{4} a_k\otimes b_k = 0$ exactly, the set $\{(k,k)\}$ has rank 3,
every one of its four 3-subsets is independent, and all indices are distinct on
both sides. It is therefore a genuine size-4 circuit that is **not** a fiber.
The construction is not accidental: $B$'s columns are a basis of $\ker R$
where $R$ has rows $(c_k\,a_k[i])$.

**Exhaustive census.** In that product, all $\binom{16}{4}=1820$ subsets were
swept: **156** size-4 circuits, **zero** of them fibers, splitting as 144 with
three distinct indices per side and 12 with all four distinct. With five
columns per factor: 984 = 900 + 84, again zero fibers. Meanwhile every size-3
circuit is a fiber (32 and 100 respectively), exactly as T-DGE requires.

**Thm-CRIT (proved).** Any all-distinct size-$(d+1)$ circuit satisfies
$\rho_U+\rho_V\le d+1$, where $\rho$ are the ranks of the selected columns on
each side; since spark forces $\rho_X\ge\min(d+1,d_X-1)$, the condition
$$d_A+d_B\;\le\;d+3$$
is necessary. In particular for a **tied** minimum $d_A=d_B=d$ this reads
$d\le3$, so the all-distinct channel is **closed for every tied minimum with
$d\ge4$**. Verified with zero mismatches: $(3,3)\to156$ and $(3,3)$ at five
columns $\to984$ are feasible; $(4,3)$, $(3,4)$, $(3,5)$ and $(4,4)$ at $d=4$
all give zero all-distinct circuits, the last inside $\binom{25}{5}$.

**Thm-FIB (proved).** A size-$(d+1)$ axis fiber is a circuit **iff** the
varying factor columns form a $(d+1)$-circuit of that factor, giving
$\sum_i (\prod_{j\ne i} s_j)\,C_i(d+1)$; for an MDS factor $C_i(d+1)>0$ iff
$d_i=d+1$ exactly, and then it is $\binom{s_i}{d+1}$. Checked on non-MDS
factors and on three-factor RS anchors against both census and formula:
$(3,3,4)\to24$, $(3,4,4)\to16$, $(3,5,5)$-type $(3,3,3)\to27$,
$(4,4,4)\to48$.

**The count is field-dependent — registered OPEN, not papered over.** The
all-distinct channel is determinantal, so its cardinality varies with the
field while the mixed channel does not. Owner-independent reproduction of the
same witness product across primes: $p=7\to152$, $p=11\to148$,
$p=13\to156$, $p=17\to148$, i.e. the $(3,3)$ channel is **constant at 144**
and only the $(4,4)$ channel moves (8, 4, 12, 4); the agent's own rows
($\mathrm{GF}(7)\to152$, $\mathrm{GF}(31)\to148$) agree exactly where they
overlap. An exact all-distinct count function is therefore left OPEN, together
with the repeated-mixed census and the $\ge3$-factor case.

**Enabling identity (MACHINE-VERIFIED on every census, zero failures).**
$\ker(A\otimes B) = (\ker A\otimes F^{s_B}) + (F^{s_A}\otimes\ker B)$, an
equality by dimension count since $\operatorname{rank}(A\otimes B)=\rho_A\rho_B$.
This puts the $|S|>d$ question in the same "sum of axis-constant tensor
spaces" frame as $V$ itself, and it is the mechanism behind the non-fiber
circuits: cancellation between the two summands, the same mechanism as the
H-MINLINE diagonal witness one size down.

**Boundary sharpness.** T-DGE's "at most one spark-2 factor" hypothesis stays
necessary here: two spark-2 factors at $d=2$ give 40 size-3 circuits of which
only 16 are fibers, and the rank criterion does not bar them ($2\le3$).

Controls in both directions: pure-tensor sweeps accept with full primal
expectations; a planted non-tensor column produces 5 non-fiber dependent
triples against 0 in the pristine matrix; a planted duplicated factor column
drops the measured spark $3\to2$ with a merged parallel class and 4 dependent
pairs below $d$. Four instrument defects were disclosed, including a
three-factor column builder that used block concatenation instead of a true
Kronecker product — caught because the pinned anchor $(3,3,4)\to24$ returned
0 — and the superseded assertion record is preserved in-run.

Run `20260901T124213Z_b77f0809_378bd1faf353`, gate `H-DP1-CIRCUIT`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`b9f1c02e14d11aa6d59469de2ef1973118b422578f91d371c533ddcb64f5654c`
(7 files, 7/7 verified), 34.7 s of a 45-minute cap. Canonical artifact:
[`rs-pe3d/.../theorem_hd1_circuits.md`](rs-pe3d/campaigns/20260901T124213Z_b77f0809_378bd1faf353/theorem_hd1_circuits.md).

**Scope.** Thm-CRIT and Thm-FIB are two-factor statements with measured
sparks at subset size exactly $d+1$; the necessary condition is necessary, not
sufficient, and no exact all-distinct count is claimed. Named next campaign:
the repeated-index (mixed-profile) census and the $\ge3$-factor extension.

## `kg/` band [4.083,6.0] tile 3 certifies at panel cap $2^{19}$ — the registered open frontier closes, 2026-09-01

**The cell that failed is now certified, by reverting a disclosed weakening
rather than by shopping the domain.** Campaign A originally pre-registered
panel cap $2^{19}$, then weakened it to $2^{17}$ in its own disclosed
Addendum 2 purely to cut cost; at $2^{17}$ tile 3 of band $[4.083,6.0]$ stayed
OPEN. Campaign B then measured that on the failing ring cell the weakened cap
loses strictly while the original wins strictly. Campaign C restored exactly
that one knob.

**Phase-1 PASS (MACHINE-VERIFIED).** Tile 3 (c-pair
$[4.1485893120,\,4.181778026496]$) certifies **completely**: 484/484
disk-intersecting leaves certified, **zero open**, DFS stack empty, 1120 boxes
and 4140 envelope evaluations. Worst certified margin
$+5.06389732008514710540159866096\times10^{-7}\pm2.76\times10^{-37}$ —
strictly positive lower end, the pre-registered threshold — at the leaf
$a_0\in[-0.8203125,-0.8177083\ldots]$, $a_2\in[+0.5703125,+0.5729166\ldots]$
with 16384 panels; maximum panels used across leaves 524288. The refutation
gate never fired: no certified $J>d$ interval anywhere.

**Gate rigour held in the right order.** The startup battery is 16/16 with the
artifact byte-identical to A and B, and the reconciliation gate reproduced A's
cap-$2^{17}$ frontier margin **bit-exactly before any $2^{19}$ number was
computed** (radius 0, difference $-4.19\times10^{-35}$, bound $<2^{-40}$).
The knob ladder was asserted as $[4096,16384,65536,262144,524288]$ before any
cell evaluation.

**Owner-independent recompute.** Both numbers were recomputed from the
byte-frozen Campaign A instrument in a separate process, with Campaign B's
frozen frontier cell: at cap $2^{17}$,
$-1.7629178702828769324466747654642\times10^{-5}$, and at cap $2^{19}$,
$+4.4630909176417023964274202614533\times10^{-6}$ — positive, hence
certifying — agreeing with the campaign's strings to $4.19\times10^{-35}$ and
$3.26\times10^{-36}$ respectively (the residual is display truncation of the
30-significant-figure claim).

**Phase 2 deferred as a quantified negative, not a failure.** Tiles 4–49 at
cap $2^{19}$ project to **78.8 CPU-h** on the optimistic model (157.6
like-tile-3, 265 on the measured multiplier) against the pre-registered
**40 CPU-h** budget — at least $1.97\times$ over. Zero Phase-2 CPU was spent
and the deferral is recorded in `phase2_deferral.json`. Named next campaign:
tiles 4–49 under a $\ge160$ CPU-h prereg with the now-proven tile-granular
checkpointing, or a jointly derived (even,odd) coefficient-space tightening
first.

Campaign `20260901T055808Z_KgCampC_pan2p19_campC` (frozen at commit
`95922c8`; predates the lifecycle tool and was finished as-is), 16-entry
checksum ledger replays green (owner-verified 16/16).

**Scope.** Swept: tile 3 only, over the full even disk (132-root cover, all
compatible odd parts) under A's byte-shared 2-D even-box envelope with closed
tail at 256-bit outward Arb, plus the $2^{17}$ reconciliation cell and the
$2^{19}$ frontier datum. **No band PASS is claimed for $[4.083,6.0]$ beyond
tile 3**; A's tile-1/2 PASSes stand as A's, A's tile-3 OPEN stands as a
statement about A's own $2^{17}$ instrument, the `kgcloseB3` retraction
stands, and nothing here bears on the source paper's correctness. The headline
$K_G$ bound is untouched by this campaign.

## `omega/` $p_{\rm comp}$-vs-penalty ordering not pinnable on the kernel ray — the straddle is not a box artifact, 2026-09-01

**Second route closed the same day, and it closes a class rather than a
guess.** After branch-pinning by freezing $\varepsilon$ was falsified, the
registered successor was to pin the $p_{\rm comp}$-versus-penalty ordering.
Restricting the domain from the whole $10^{-7}$ box to the **certified
exact-kernel ray** $\{p^*+\rho s:\rho\in[0,10^{-7}]\}$ — with $s$ the frozen
candidate direction, $A s=0$ over $\mathbb Z$, $\sum s=0$, kernel membership
re-proved in-run on all 27 rows with $A$'s rank re-established as exactly 24
over $\mathbb Q$ — does **not** pin it. At the certified endpoint both order
inequalities straddle zero:
$G_1\in[-1.8897456060889468\times10^{-5},\,+1.6396874550250076\times10^{-5}]$,
width $3.5294330611\times10^{-5}$ = **223.15×** the
$1.5816497000997742\times10^{-7}$ signal, and $G_2$ within $4\times10^{-15}$
of it.

**The decisive quantity is the survival fraction.** The kernel ray retains
**62.7%** of the parent's whole-box $G_1$ width
($3.5294\times10^{-5}$ of $5.6292\times10^{-5}$) and the candidate hull span
drops only $1.52\times$ ($3.2123\times10^{-5}\to2.1109\times10^{-5}$, still
133.46× signal). So the parent obstruction was **not** an artifact of taking
the whole box, and domain restriction is removed as a mover. All four ratios
were recomputed independently by the owner and match to the digits quoted.

$\varepsilon_0$ is again bit-constant along the ray
($2.1502086942121845\times10^{-6}$ at both $\rho=0$ and $\rho=10^{-7}$), a
third independent confirmation, so $\varepsilon$ is exonerated for every
future analysis and the straddle is carried by $p_{\rm comp}$-vs-penalty
motion itself. Controls pass in both directions: GATE-Z reproduces both frozen
records, C1 accepts the identity at gap $5.79\times10^{-13}$, C2 rejects a
wrong-index dual lookup, C3 rejects dual drift with exact $4\,$drift coupling.
Five runner defects were disclosed with per-attempt logs, and Amendment 1
(commit `2bc7628`) — which corrected a transcription error in the expected
positivity wall — was committed **before** any verdict compute; the exact wall
was then re-derived from first principles as
$\mathrm{Fraction}(4771454744492531,\,295147905179352825856)$, bit-identical
to the frozen record.

Run `20260901T132102Z_f17bb4cb_117823c32d87`, gate `OmegaPcVersusPen`,
**FROZEN-NEGATIVE**, checksum-ledger SHA-256
`025b41b85135c0bbee32ab7fc4fdbae3c8cd9c1dca3f6b523ae4b5c32f2d6d8b`
(18 files, owner-verified 18/18).

**Provenance correction (owner-found).** That run's
`PREREG_PROVENANCE.json` records the init-time prereg hash `443fc302…` with
`byte_identical: true`, but the in-directory `pre_statement.md` hashes
`03fae37c…` because Amendment 1 edited it after the copy. The frozen
`sha256s.txt` pins the amended bytes, so no evidence is lost and the chain is
reconstructible, but the provenance file reads misleadingly on its face. Rule
going forward: **a prereg amended after `campaign.py init` must have both
hashes recorded in the provenance.**

**Scope.** The 21-dimensional kernel question is **UNCHANGED and OPEN**, gate
C stays PARTIAL, and the 12.808× guardrail means this line still cannot
produce a record. What is now certified is narrower and more useful than a
guess: two movers ($\varepsilon$-freezing, domain restriction) are eliminated
with measured factors. Named next route: a certified derivative-sign attack,
or the other six branches, whose whole-box spans are $\le2.06\times10^{-6}$ —
roughly 13× signal, about 15× narrower than $R_{\rm glob}[0]$.

## `rs-pe3d/` the crossing channel — closed-form count, field-independence, and coherence with Thm-CRIT, 2026-09-01

**The third and structurally decisive channel at $|S|>d$ is now classified.**
T-DGE handles $|S|=d$; run 4 falsified the fiber picture at $|S|=d+1$ and
proved the rank criterion; this run explains *where* the non-fiber circuits
come from and counts them exactly.

**Construction.** Using the kernel identity
$\ker(A\otimes B)=(\ker A\otimes F^{s_B})+(F^{s_A}\otimes\ker B)$, take one
column supported on a $d_A$-circuit $R$ of $A$ and one row supported on a
$d_B$-circuit $J$ of $B$; if the crossing row $i_0\in R$ and crossing column
$j_0\in J$, the scaling $\sigma=-c_{i_0}/\delta_{j_0}$ (available over any
field, both entries being nonzero) cancels the centre cell exactly.

**Theorem X (proved directions).** The resulting support has
$$|S| \;=\; (d_A-1)+(d_B-1) \;=\; d_A+d_B-2 \;=:\; m_*$$
with profile exactly $(d_A,d_B)$, is a genuine circuit with a
1-dimensional relation space, distinct parameter tuples give distinct
supports, and the count is
$$N_{\rm mix} \;=\; d_A\,d_B\,C_A(d_A)\,C_B(d_B).$$
**Field-independence** holds for fixed factor circuit spectra, because the
construction uses only the existence and full support of factor circuit
relations, a scaling possible in any field, and a field-free injectivity
argument.

**Coherence — the two channels meet exactly at the criterion's boundary.**
The crossing channel emits at its own size $m_*$, so it lands in the
$|S|=d+1$ window **iff** $d_A+d_B-2=d+1$, i.e. **iff $d_A+d_B=d+3$** — which
is precisely the tightness case of run 4's necessary condition
$d_A+d_B\le d+3$. For tied factors this is $2d-2=d+1$, i.e. $d=3$, which is
why the earlier tied-$d{=}4$ configuration had *no* size-5 circuits at all:
the all-distinct channel is barred by Thm-CRIT and the crossing channel sits
at size 6 instead.

**Set equality, not merely matching counts.** At six configurations the
census population with profile $(d_A,d_B)$ was proved **equal as a set of
supports** to the constructed family: the witness 144, a $2\times3$
Vandermonde row 36, $2\times5$ 900, $2\times6$ 3600, the $3\times4$ MDS pair
16, and a composite regrouping 36; the relation-space dimension $m_*-1$ was
verified on every member. Owner-independent reproduction: $9\binom{n}{3}^2$
gives 144, 900 and 3600 at $n=4,5,6$, and the $3\times4$ pair at $d=4$ has
**zero** size-5 circuits with exactly **16** $=4\cdot4\cdot1\cdot1$ size-6
circuits, all of profile $(4,4)$.

**Cross-prime contrast, independently reproduced.** At
$p\in\{7,11,13,17,31\}$ the factor spectra stay $\{3\!:\!4\}$ and the mixed
population is **144 at every prime**, while the all-distinct residual varies
$8/4/12/4/4$ with totals $152/148/156/148/148$. The owner's own sweep gives
$p=7\to152\,(8)$, $p=11\to148\,(4)$, $p=13\to156\,(12)$, $p=17\to148\,(4)$ —
agreeing exactly. This is the sharp separation between a field-free
combinatorial channel and a determinantal one.

Controls, 66/66 pinned asserts with zero failures over $\ge163{,}000$
enumerated subsets: **T3a** perturbs each of the 144 crossing relations in
every single coefficient (6336 perturbations) and *all* destroy dependence,
proving the cancellation is essential rather than cosmetic; **T3c** rejects a
duplicated factor column (spark $3\to2$, four size-2 dependent pairs below
$d$); **T3d** rejects a non-tensor column plant (5 non-fiber dependent
triples against 0 pristine); **T8** is a dedicated guard against the
block-concatenation-versus-Kronecker defect that bit the previous run, and it
separates them cleanly (true 3-factor anchors 24/16/27/48 versus the defect
path's divergent 0/0/1/1732). Four instrument/pin defects were disclosed,
including three mis-stated sweep constants in the committed prereg.

Run `20260901T132343Z_4f169cbb_2f81dd8f2d8d`, gate `H-MIX-CROSSING`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`7afdc91a1701e86a65a054965f5aab77321f32361201b72daa418cf61f957c15`
(8 files, owner-verified 8/8), 50.3 s of a 3600 s cap.

**Scope, stated at the boundary the proof actually reaches.** The forward
construction, the closed form, the injectivity and the field-independence
claim are general. The **converse** — that every size-$m_*$ circuit of
profile $(d_A,d_B)$ is of this form — is proved for $d_A,d_B\ge3$ by Lemmas
P, Q, Q2 and R **except one corner**: $(d_A,d_B)=(3,3)$ with a
$(2,2)$-decomposition, where Lemma Q3 does not discharge the case
analytically. That corner is closed **by exhaustive census at every swept
$(3,3)$ configuration** (set equality, not counts) and is registered as
**Hypothesis H2**; the general $(3,3)$ converse is therefore **OPEN**. The
spark-2 boundary branches, all-distinct closed forms, $\ge3$-factor general
theorems, characteristic 2, and $|S|$ beyond the swept rows are all out of
scope.

## `oct-rank/` both registered routes to $\{13,14\}$ closed with quantified evidence, 2026-09-01/02

**The frontier did not move, and that is now a measured fact rather than a
lack of trying.** $\operatorname{rank}((L_1,L_i,L_j))$ stays **OPEN in
$\{13,14\}$**: floor 13 and upper 14 both certified, Route A's $12<13$
closure intact, the published window $18\le R_{\mathbb R}(T_O)\le25$
untouched. What changed is that the two routes the owner re-scoped — exact
CP completion downward and Gröbner infeasibility upward — are each frozen with
the exact cost at which they fail.

**N1, exact CP completion (FROZEN-INCONCLUSIVE, MACHINE-VERIFIED
negative-of-search).** Does an exact rank-13 CP decomposition of
$T_F=\operatorname{blockdiag}(\tau,\tau)$ — rank-equivalent to
$(L_1,L_i,L_j)$ by two frozen entrywise diagonal-similitude identities — exist,
findable from 40 declared seeds (32 random, 7 $\tau_{r7}$-merge, 1 frozen
restart) via $3\times5000$-nfev TRF followed by exact $\mathbb Q$ Newton
admission? **No candidate was admitted**: 40/40 seeds completed in 3769.7 s
CPU, exact residual range $1.4\times10^{-5}$ to $7.8\times10^{-4}$, best
relative residual $1.414\times10^{-6}$ — four orders above the pre-registered
exact gate $\|g_0\|_{\max}\le10^{-7}$ — and the E4 kill criterion fired
exactly as registered. Owner re-read of the frozen `n1_results.json`: 40
seeds, minimum `rel_final` $1.413946\times10^{-6}$, zero at or below the gate.
All four controls in their expected directions through the byte-identical
frozen instrument: $\tau_7$ CONTAINMENT at $\rho=10^{-3}$ (accepts a known
decomposition), $\tau_6$ NO-CONTAINMENT at every rung (rejects a globally
infeasible one), synthetic P-POS CONTAINMENT at $10^{-4}$, and P-WRNG
NO-CONTAINMENT with complete certified exclusion for a one-coordinate
perturbation. A pre-outcome instrument defect ($\tau^{\mathsf T}$ builder
emitting 4-vectors) was fixed and disclosed with the broken output preserved.
Run `20260901T103338Z_06de66a3_d0b7bc2dc6a5`, checksum-ledger SHA-256
`a137e5eae5ea5d15b37c7a6c9f7170d8e72a3d92578da6a48f8bdb543203a510` (11 files,
owner-verified).

**N2, Gröbner infeasibility (FROZEN-NEGATIVE as a quantified resource
negative).** The pre-registered protocol required calibrating msolve on
known-status systems before touching the 273-variable gauge-extended rank-13
systems $E$ (218 equations) and $E_1$ (220). The **smallest** calibration
system, the known-infeasible $\tau$ at $r=6$ with 66 variables and 48
equations — four times smaller than the main systems — blew its 30-minute S1
cap while still inside F4 degree 9: matrix $8{,}080{,}476\times74{,}499{,}764$,
22.65 GB RSS at kill (35.93 GB in an earlier attempt), degree 8 alone costing
3802 s wall / 3593 s CPU, i.e. more than twice the whole cap in a single
round. Three attempts, all 0-byte output; K-1 ($\tau$ at $r=7$) had already
exceeded its 120 s probe. Per the registered S1 rule, S2 on $E$ and $E_1$ at
primes 32003, 65537 and 2147483647 was **not run**; primes, monomial order,
systems and caps are unchanged from the prereg. Frozen for reproduction:
msolve 0.8.0 pinned to the Sage 10.7 tarball (SHA-256 `319ba0de…badfd`,
binary `1e8a6725…bc21`), build logs, all serialized integer systems with
hashes (owner-verified 8/8), the deterministic generator with generation-time
anchors, exact commands, full F4 traces, an OS memory snapshot, and the
0-byte outputs as witnesses of non-termination (owner-verified: 5 present).
Disclosed: $E_1$ has 273 variables, the prereg's "275" being a miscount; and
the verdict label FROZEN-NEGATIVE was set by owner instruction where the
prereg had registered FROZEN-INCONCLUSIVE for the same BUDGET-FAIL semantics
— recorded in `VERDICT.md` rather than silently substituted. Run
`20260901T115801Z_8d6683fb_39c752e7eadd`, checksum-ledger SHA-256
`c8d027ac1b8728faff08ae9a91be0204c0fbf953f2e355c76242328c1328c891` (47 files,
owner-verified).

**Scope, verbatim from the frozen verdict.** "Absence of a 13-witness is NOT
evidence for 14; absence of an impossibility argument is NOT evidence for
13." A resource negative is not a rank claim. What the two runs establish is
narrower and useful: local TRF-plus-exact-Newton from 40 seeds does not reach
the exact gate, and F4 Gröbner elimination is infeasible by orders of
magnitude even on a system four times smaller than the target. Any successor
must change technology — structured or symmetry-reduced elimination, or a
different certificate of infeasibility — rather than rerun either route.

## `mm3/` mws59 lower bound 56 → 58 — and a correction to frozen session-3 evidence, 2026-09-02

**Frontier move.** `mws59`, the widest remaining gap in the frozen landscape
at $[56,59]$, is now **$\{58,59\}$** for its fixed orientation. Nothing here
improves on the published 59-addition scheme; the lower bound moved up by two
with certificates.

**Exact stagewise decision (MACHINE-VERIFIED).** With
$\text{total}=C(U)+C(V)+C(W_{\rm fac})+14$ and output bounded by the
transposition principle (preconditions audited in-run: 23/23 nonzero rows per
side, $\operatorname{rank}(W_{\rm fac})=9$ exact over $\mathbb Q$):
- $C(U)=14$ **exactly**: floor@13 re-certified impossible in-run, and the
  single admitted auxiliary among 389 hash-pinned candidates yields a 14-gate
  witness verified by exact expansion (owner re-read of the frozen scan: 389
  rows, exactly 1 feasible).
- $C(W_{\rm fac})=15$ **exactly**: floor@14 impossible, and a 15-gate witness
  synthesized in-run from the Table-3 $W$ and verified by exact expansion.
- $C(V)\ge15$: single-auxiliary scan 0/366 at $T=13$ and a **complete**
  two-auxiliary pair scan 0/79,728 at $T=14$ (owner re-read: 366 outer rows,
  79,728 pairs done, 0 hits).
Hence LB $=14+15+15+14=58$ and UB $=59$ (published 15/15/29 recounted in-run).

Every one of the 80,483 instances was decided twice — the frozen memoized
subset-DFS and a fresh slot-availability CNF under kissat 4.0.4 — with
**zero disagreements**; five DRAT→LRAT certificates for the load-bearing
UNSATs were accepted by both pinned checkers; controls A1–A4 accept and R1–R4
reject all pass, including a cross-solver CaDiCaL confirmation of the SAT
control; an independent post-run audit passed 7/7 (third-path universe
re-enumeration, independent pair count 79,728, checkpoint completeness,
boundary and random DFS re-runs, verdict arithmetic, witness re-verification,
certificate replay). Cost 1.33 of a 24.0 CPU-h budget.

Run `20260901T133118Z_8c14c3dd_7970c5318d5d`, gate `mws59-gap-total59`,
**FROZEN-NEGATIVE** for the 59-total question as registered, checksum-ledger
SHA-256 `d9075921dc4fbc23c029c4dd6d50a0d313a67fecf1a263770c039e2c1457e7a4`
(37 files, owner-verified 37/37).

**CORRECTION to frozen evidence (owner-adjudicated).** The session-3
campaign `2026-08-30T012035Z_a7ea8e8e_cbdf8e63fa94` supplied
$C(W_{\rm fac}^{\rm mws59})\le15$ through a "constructive transposition
check", `transpose_check.py`, whose `mws59_resolution.json` and manifest
record "15 additions, `computes_Wfac_exactly: true`". That script is
**self-circular**: its `chain` dictionary is created and never populated
(owner-verified: five references, no item assignment), so it registers no
output consumers, and its own stdout reports "transposed gates: 0 / Wfac
circuit additions = 0 / computes exactly the W factor map: False". The frozen
verification was therefore **void**. The number 15 happens to be correct and
is now established independently by this run's exact-expansion-verified
witness; the frozen directory is untouched, and this paragraph is the
authoritative adjudication of that artifact. Lesson recorded: a verifier's
own stdout must be compared against the JSON it emits before either is
frozen.

**Scope.** Fixed orientation, three-stage linear SLP model, additions
including subtractions, free inputs and sign changes, output via transposition
accounting. The residual 58-versus-59 decision needs a $U$ two-auxiliary or
$V$ three-auxiliary scan at $T=15$ — projected $10^6$–$10^7$ instances at the
measured 21 ms each, explicitly unregistered, and possibly infeasible; if so
it is to be frozen as a measured negative. `laderman23` is not decidable by
this machinery at all (no verified factors). Next on the ladder:
`stapleton60` at $[58,60]$.

## `rs-pe3d/` Theorem X's converse completed (H2), and the all-distinct channel is cross-ratio preservation, 2026-09-02

**Two certified closures finish the size-$(d+1)$ theory for two-factor
products with sparks $\ge3$.**

**H2 discharged (HUMAN-AUDITED proof, MACHINE-VERIFIED sweeps).** The one
corner the crossing theorem's converse had left to census — $(d_A,d_B)=(3,3)$
with a $(2,2)$-decomposition, two $\ker A$ columns at $v_1\ne v_2$ and two
$\ker B$ rows at $u_1\ne u_2$ — is now proved to produce **no new support**.
The argument is a tight count: every cancellation cell lies in
$X=\{u_1,u_2\}\times\{v_1,v_2\}$ (a cell outside $S$ where $C_1\ne0$ forces
$C_2\ne0$ there); $|\operatorname{supp}C_1|,|\operatorname{supp}C_2|\ge6$
with intersection inside $X$, so $|S|\ge(6-4)+(6-4)=4=m_*$ is tight;
equality forces both column supports to be $\{u_1,u_2,w_i\}$ and both row
supports $\{v_1,v_2,x_i\}$ with all four cells of $X$ cancelling; profile
$(3,3)$ then forces $w_1=w_2=w$ and $x_1=x_2=x$; and
$S=\{(w,v_1),(w,v_2),(u_1,x),(u_2,x)\}=S(T,w,Z,x)$ with $T=\{u_1,u_2,w\}$ a
3-circuit of $A$ and $Z=\{v_1,v_2,x\}$ a 3-circuit of $B$ — exactly a
crossing support, with the relation proportional to the crossing relation
because a circuit's relation space is one-dimensional. Owner derivation
confirmed numerically before dispatch (forced $\Gamma$ has support identical
to $S(T,w,Z,x)$, $A\Gamma B^{\mathsf T}=0$, genuine circuit); the campaign
discharged all five steps including the inside-$S$ overlap and the
spark-to-circuit step, with 170 exact asserts, six exhaustive $(2,2)$ sweeps
over $p\in\{7,11,13\}$ and two factor pairs giving profile-$(3,3)$ set
equality (144 and 90) with zero new supports, both plants rejected and the
concat guard passing. **Theorem X's converse is therefore complete for all
$d_A,d_B\ge3$.** Run `20260901T135534Z_5678c46e_232421125d33`, gate
`H-MIX-H2-CLOSURE`, **FROZEN-CERTIFIED**, checksum-ledger SHA-256
`389a74dc3b38bdccf0bc7680c300e6d240472c648860854e7df5a981fbf52e1a` (6 files,
owner-verified); `SCOPE_NOTE_H_MIX_CROSSING.md` and `state.json` updated at
target level.

**The all-distinct channel is cross-ratio preservation (HUMAN-AUDITED proof,
MACHINE-VERIFIED set equalities).** For $2\times n$ GRS factors with columns
$(1,x_u)$ and $(1,y_v)$ and distinct evaluation points, $r_Ar_B=4=d+1$, so
an all-distinct size-4 support $\{(u_k,v_k)\}$ carries a circuit iff a single
determinant vanishes, $\det[1,\;y_k,\;x_k,\;x_ky_k]=0$, and the chain of
equivalences is now proved:
$$\det=0 \iff \exists\,F=c_0+c_1y+c_2x+c_3xy\ne0 \text{ vanishing at the four points}
\iff \text{the pairs lie on one nondegenerate Möbius graph } y=-\tfrac{c_2x+c_0}{c_3x+c_1}
\iff \operatorname{CR}(x_{u_1},\dots,x_{u_4})=\operatorname{CR}(y_{v_1},\dots,y_{v_4}),$$
with the three degeneracies handled explicitly: $D=c_3x+c_1\not\equiv0$, the
factor-degenerate case $c_1c_2=c_0c_3$ (a vertical-plus-horizontal line pair
that cannot hold four points with distinct coordinates), and the absence of a
selected pole; conversely every nondegenerate Möbius map $y=(ax+b)/(cx+d)$
with $ad-bc\ne0$ yields such an $F$. Triple-independence of the three
proper subsets gives circuitness, so the predicate is exactly "the pairing is
a Möbius restriction". Hence
$$N_{\rm all\text{-}distinct}=\#\{(U,\sigma):\operatorname{CR}(U)=\operatorname{CR}(\sigma(U))\},$$
and the field dependence of this channel is **exactly** the field dependence
of cross-ratio coincidences — a field equation, which is why the crossing
channel is field-free and this one is not. Owner-independent verification
preceded dispatch: the determinantal count reproduces the frozen $(4,4)$
census exactly at $p=7,11,13,17$ (8/4/12/4) and predicts 4 at $p=31$ where
the frozen total is $148=144+4$; cross-ratio equality agreed with
$\det=0$ on all 120 bijections across those five primes and on all 600
$(U,\sigma)$ pairs of an asymmetric $5\times5$ configuration at $p=13$ and
$17$, zero disagreements. In-run: 176 exact asserts; determinant, cross-ratio
and census populations equal as **sets** on five symmetric primes (8/4/12/4/4)
and two asymmetric $5\times5$ runs (84 and 40); on every accepted support a
nonzero bilinear left-null vector was solved and its $\Delta\ne0$, pole-free
Möbius values and line-or-hyperbola normal form asserted; every planted
cross-ratio mismatch, non-tensor column, duplicated column and concat guard
fired. A launch-wrapper defect (shell `exec` unavailable, exit 127, before
any Python ran) was disclosed and the unchanged run replayed. Run
`20260902T011204Z_7ece99fe_a108b246f87d`, gate `H-ALLDISTINCT-CROSSRATIO`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`53eea991a0a4634dcaf5eb73fb3819fe28ca928aa0fa8e50ed1c0e1885c37e01` (7
files, owner-verified).

**The two-factor picture, now complete for sparks $\ge3$.** At $|S|=d$:
fibers over factor circuits (T-DGE). At $|S|=d+1$: fibers over
$(d+1)$-circuits (Thm-FIB); crossing supports of size $d_A+d_B-2$, counted by
$d_Ad_BC_AC_B$, field-free, present iff $d_A+d_B=d+3$ (Theorem X with H2);
and, for $2\times n$ GRS factors, all-distinct supports exactly when the
pairing preserves cross-ratio, field-dependent (this run). The rank criterion
$d_A+d_B\le d+3$ (Thm-CRIT) is the boundary all three respect.

**Scope.** Registered open, in `state.json`: the spark-2 boundary branches
of Theorem X; a simpler closed evaluation of the cross-ratio coincidence count
for arbitrary point sets and any analogue for higher-row bi-Vandermonde
factors, where $r_Ar_B>d+1$ makes the predicate a rank condition on
$\operatorname{span}\{x^iy^j:i<r_A,\,j<r_B\}$ rather than one determinant;
and three or more tensor factors.

## `mm3/` stapleton60 gap closed — LB = UB = 60 for the fixed orientation, 2026-09-02

**Second published scheme certified optimal for its own orientation.**
`stapleton60` in the $\sigma^0$ orientation goes from $[58,60]$ to
**exactly 60**. The frozen lower bound 58 was valid but not tight: both
side-stage $d{+}1$ values (15, 15) are unattainable.

**Reduction and decisions (MACHINE-VERIFIED).** $\text{total}=C(U)+C(V)+
C(\text{output})$ with $C(\text{output})=28$ exact (frozen floor@13 on
$W_{\rm fac}$ gives output $\ge C(W_{\rm fac})+14\ge28$; the printed output
stage re-expands in-run to exactly 28, reproducing all 9 $C$-entries over the
23 products). That reduces the whole gap to two single-auxiliary questions at
$T=d+1=15$, with the closure argument pre-registered: a 15-gate circuit for
14 $\tau$ classes has at most one non-$\tau$ gate value, duplicate-class gates
are deletable, and zero auxiliaries would strip to the refuted 14-gate floor.
Results: $C(U)=16$ exact (0 of 428 admitted at $T=15$; the printed 16-gate
left stage re-expanded exactly), $C(V)=16$ exact (0 of 426), and
$C(W_{\rm fac})=14$ exact (8 of 398 admitted at $T=14$; the first admission
replayed as an explicit 14-gate circuit and verified by exact expansion,
`witness_Wfac_14gates.json`). Owner re-read of the frozen checkpoints:
428/0, 426/0, 398/8, total 1,252 instances, **zero disagreements** between
the frozen subset-DFS and the fresh kissat 4.0.4 CNF.

Universes were hash-pinned before the prereg commit by enumeration only, with
**zero** feasibility instances decided pre-commit. Controls A1–A4 accept and
R1–R4 reject all pass on the first pass, including a tri-instrument SAT
control (DFS / kissat / CaDiCaL 3.0.1), exact transposition preconditions
($\operatorname{rank}(W_{\rm fac})=9$ over $\mathbb Q$), and 729/729 Brent in
both `fmpz` and pure integers; **37 DRAT→LRAT certificates** each accepted by
both pinned checkers and replayed post-run (owner count: 74 certificate
files); independent post-run audit 7/7. Cost 110.7 CPU-s of a 6.0 CPU-h cap,
so the negative is complete, not truncated. No transposition-derived witness
is used anywhere, honouring the session-3 defect lesson structurally.

Run `20260902T012056Z_25e1634d_f146ed17f287`, gate
`stapleton60-gap-total60`, **FROZEN-NEGATIVE** for the sub-60 question,
checksum-ledger SHA-256
`129ee31432ab78f94b7d5b8b37f28dc8a26d4a67e16520ce034b92c615b5a858`
(129 files, owner-verified 129/129).

**Fixed-orientation ladder after today.** `paper55` $55=55$ (the record,
untouched), `sun56` $56=56$, `mws59` $\in\{58,59\}$ (LB 58 certified),
`stapleton60` $60=60$. Two published schemes are now certified optimal for
their own orientations. **Scope** as before: fixed orientation, three-stage
linear SLP model; other orientations, decompositions and models are
untouched. Remaining: the mws59 58-versus-59 bit, `paper55` $\sigma^1/\sigma^2$
at $[55,57]$, and `laderman23` (no verified factors).

## `rs-pe3d/` the all-distinct count in closed form — a $\mathrm{PGL}(2,p)$ sum, 2026-09-02

**The registered open item "a closed evaluation of the cross-ratio
coincidence count" is resolved for full-support GRS factors.** Since each
all-distinct size-4 circuit support is the graph of a **unique** Möbius map
restricted to its four $x$-points (existence by the cross-ratio theorem,
uniqueness because three points determine a Möbius map), and conversely each
pair $(M,U)$ with $U\subseteq X$, $|U|=4$, $M(U)\subseteq Y$ and no pole in
$U$ yields a distinct support,
$$N_{\rm all\text{-}distinct}(X,Y)\;=\;\sum_{M\in\mathrm{PGL}(2,p)}\binom{|X\cap M^{-1}(Y)|}{4},$$
with $M^{-1}(Y)$ excluding the pole, $|\mathrm{PGL}(2,p)|=p(p^2-1)$, and
exactly $p-1$ maps for each ordered pair of distinct (pole, zero) in
$\mathbb P^1$. Classifying maps by where their pole and zero fall relative to
the point set gives two closed forms (HUMAN-AUDITED, MACHINE-VERIFIED):
- full field $X=Y=\mathbb F_p$: the $p(p-1)$ affine maps keep all $p$
  points and the other $p^2(p-1)$ lose their finite pole, so
  $$N=p(p-1)\binom{p}{4}+p^2(p-1)\binom{p-1}{4};$$
- multiplicative group $X=Y=\mathbb F_p^*$, $q=p-1$: classes of size $2q$,
  $4q^2$, $q^2(q-1)$ by the number $k\in\{0,1,2\}$ of $\{\text{pole},\text{zero}\}$
  inside $\mathbb F_p^*$, summing to $p(p^2-1)$, each keeping $q-k$ points, so
  $$N=2q\binom{q}{4}+4q^2\binom{q-1}{4}+q^2(q-1)\binom{q-2}{4}.$$

Owner-independent verification preceded dispatch: brute-force product census
with full $2\times p$ Vandermonde factors gives 200 at $p=5$ and 5880 at
$p=7$, equal to the first form; on $\mathbb F_p^*$ the census gives 8, 1080
and 117600 at $p=5,7,11$, equal to the second form **and** to a direct
enumeration of every Möbius map (336 at $p=7$, 1320 at $p=11$) under the
principle. An owner hand-derivation slip — $q^2(q-2)$ for the $k=2$ class,
giving 1044 instead of 1080 at $p=7$ — was caught by that brute force and
turned into a control: in-run the wrong coefficient is rejected by the
class-size identity ($300\ne336$) and by the count ($1044\ne1080$). In-run:
151 exact asserts, zero failures; PGL sizes and class sizes at
$p=5,7,11,13$; three-route support-set equality (census, cross-ratio
predicate, PGL enumeration) at 200/5880 and 8/1080/117600, including
1,058,400 pairings at $p=11$; selected-pole, non-tensor, duplicated-column
and concat guards all fired.

Run `20260902T012629Z_2eb59adb_0aa7425efe1b`, gate
`H-ALLDISTINCT-PGLCOUNT`, **FROZEN-CERTIFIED**, checksum-ledger SHA-256
`1ff319febd6d3da614513981b9bf265b88012057b3f0c091f3bc5c8b26df6096` (6 files,
owner-verified).

**Scope.** Two $2$-row GRS factors with evaluation sets inside $\mathbb F_p$;
the sum is general in $X,Y$, the closed forms are for the two full supports.
With this, the field dependence of the all-distinct channel is fully
explained: it enters only through $p$ and the point sets, via
$\mathrm{PGL}(2,p)$. Still open: the spark-2 boundary branches of Theorem X,
higher-row bi-Vandermonde factors, and three or more tensor factors.

## `mceliece/` second $m=12$ row reached: $t=96$ with $\alpha$ measured, not derived, 2026-09-02

**A bounded finite statement about one more instance, exactly as registered.**
At $(m,n,t,\text{seed})=(12,3488,96,16384)$ — $k=n-mt=2336$,
$D=n-2t-1=3295$, $2t+3=195$ (owner-recomputed) — the binary-Goppa Guard-A
construction passes with the same standard as the $t=64$ NIST cell:
- $G$ monic irreducible of degree 96 (312 seeded tries, scalar-Rabin
  confirmed), square-free by $\gcd(G,G')=0$;
- **$\alpha$ MEASURED**: the fork resolved at the pre-registered gate on
  measured cumulative CPU $3353.46\ll24300$, so the unmodified $m\le11$
  instrument ran over **all 3488 support points with 0 failures** (648.43
  CPU-s); the DERIVE branch (Apon Lemma 3 as CITED-DEPENDENCY) was never
  entered;
- $\delta$ chain **exact** 5/5: zero-remainder division 3488/3488, unit
  checks $L_i(a_i)=1$ 3488/3488, 0/200 off-diagonal violations, 3/3 assembly,
  degree gate $\max_j\deg f_j=3295=D$ exactly;
- controls: $\lambda$ recompute byte-equal, four-selector $\mathbb F_2$
  linearity, $\beta$ pair $(0,1)$, $\varepsilon$ max-degree; anchors 3/3
  byte-equal against the frozen $m\le11$ records at $(11,2048,48)$,
  $(10,1024,40)$, $(6,64,3)$; plants 4/4 rejected by the same verifier —
  duplicated support ($\Pi'=0$), $f_0+\Pi$ caught by the degree-gate abort
  ($2048>1951$) with value-level invisibility confirmed 0/22, row-0 $\times3$
  failing $\alpha$ ($c^2=5$), and a $t{=}96$-native $\times3$ probe failing
  60/128 corrupted against 0 clean.
Cost 4017.55 of a 32400 CPU-s budget. Run
`20260901T100025Z_bd55fdfd_8fbba9be5ab3`, **FROZEN-CERTIFIED**,
checksum-ledger SHA-256
`f8e7d407ef9d121a4a37b91fe34f04b35c5a1b3569144c5ea8684d0c64c936d0`
(20 files, owner-verified 20/20; `verdict.json` fields
`alpha_identity_measured_full_grid: true`, `alpha_failures: []`).

**Governance.** An earlier same-cell run, `20260901T094531Z_ae227cfd_fe99571a47fa`,
was closed **REHEARSAL** and all its outputs discarded because tool-level
$t{=}96$ probes preceded the prereg commit; its `REHEARSAL_PROVENANCE.md`
records the timeline. The certified run above started its first $t=96$
computation only after commit, init and claim.

**Scope.** ONE instance, one seeded $G$, one seeded support ordering.
**Not** "$m=12$ complete"; $t=48$ remains unreached and is the named next
row; no census rows at $m=12$; the Apon section-3.6 hole remains **bounded,
not closed**; no NIST-cost claim.

## `mm3/` the mws59 bit is 58 — an explicit 58-addition circuit, one below the published 59, 2026-09-02

**The ladder is now exact at every rung.** `mws59`'s fixed orientation has
minimum **exactly 58** additions: $C(U)+C(V)+C(\text{output})=14+15+29$. The
published 59-addition scheme is therefore **one addition above the optimum of
its own orientation** — its left stage spends 15 where 14 suffice. No
published claim is contradicted (59 additions do suffice) and 58 is **not a
record**: `paper55`'s 55 stands untouched.

**How the two sides met.** Lower bounds are entirely the frozen previous run:
$C(U)=14$ exact, $C(V)\ge15$ (single-auxiliary 0/366 and complete
two-auxiliary 0/79,728), $C(W_{\rm fac})=15$ exact, output $\ge29$. The upper
bound became a certified object by assembling session-10's frozen 14-gate
left circuit with the paper's printed 15-addition right stage and 29-addition
output stage (Table 2, transcribed verbatim from the pinned layout, recounted
in the fixed model as $15+15+29=59$, and matched to the Table-3 blocks by an
exact product correspondence: identity coordinate maps, a bijection on the 23
products, per-product signs, nine negated products compensated in the output
column, all 621 entry equalities exact). The decisive artifact
`assembled_58_endtoend.json` expands the assembled circuit **symbolically over
all 81 monomials** $A_iB_j$ and finds all nine outputs equal to
$C[i][j]=\sum_l A[3i+l]B[3l+j]$ exactly, with the total recounted from
explicit gate lists as 58 — independent of the DFS/CNF searches, of the Brent
check on the factor blocks, and of any transposition argument.

**Owner-independent re-expansion.** With a fresh expander over exact
integers: the published Table 2 gives 9/9 outputs exact at $15+15+29=59$; the
frozen 14-gate witness (`witness_U_14gates.json`, each recorded gate value
re-derived from its operands) realises all 23 left operands of the printed
products, 19 exactly and 4 (products 15, 16, 18, 21) up to sign; absorbing
those signs into the output stage, which is free in the model, the assembled
algorithm gives **9/9 outputs exact at $14+15+29=58$**.

**Route S never ran, by its own pre-registered numbers.** The registered
fallback — a three-auxiliary closure enumeration — projected $6.30\times10^7$
ordered chains at the measured 59–120 ms per instance, i.e. 173–2,100 CPU-h
against a 4.0 CPU-h cap; Route P settled the bit first, so no prefix was run,
nothing was sampled, and no Route-S claim exists. Measured cost 2.2 CPU-s.

**Amendment 1, disclosed and committed before the amended compute.** The
first correspondence matcher forced $\varepsilon\varepsilon'=+1$, forbade a
compensating output-column sign and assumed one shared $3\times3$ flattening,
and found no correspondence although it recounted 15/15/29 correctly. It was
widened to $\sigma_{A,B,C}\in\{\mathrm{id},T\}$ with independent per-product
signs, accepting only on exact equality of every entry, with the
wrong-correspondence plant strengthened to the widened space; the failed
attempt is preserved verbatim (`verdict_preamendment.json`). Only the
$\varepsilon\varepsilon'=-1$ case was needed. `PREREG_PROVENANCE.json`
records both prereg hashes (owner-verified).

Controls A1–A4 accept (729/729 Brent in `fmpz` and integers, frozen row, both
session-10 witnesses re-verified by exact expansion, tri-instrument SAT
control, transposition preconditions) and R1–R5 reject (gate deleted,
operand perturbed, wrong correspondence, frozen floor@12 re-certified with a
fresh dual-checker DRAT→LRAT, $d$-count planter); independent post-run audit
6/6. Run `20260902T013730Z_49d2f937_a672bbec349b`, gate
`mws59-bit-58vs59`, **FROZEN-CERTIFIED**, checksum-ledger SHA-256
`59718ab491e66b196ddd8f79f1eb33880f04fe859fe9290b8d8f03b0ddd8d661` (22
files, owner-verified 22/22).

**Exact fixed-orientation ladder, published → certified minimum.**
`paper55` $55\to55$ (the record), `sun56` $56\to56$ (published optimal),
`mws59` $59\to\mathbf{58}$ (published $+1$), `stapleton60` $60\to60$
(published optimal). **Scope** unchanged: fixed orientations, three-stage
linear SLP model; other orientations, decompositions and the global sub-55
question are untouched. Remaining in `mm3/`: `paper55` $\sigma^1/\sigma^2$
at $[55,57]$ and `laderman23` (no verified factors).

## `rs-pe3d/` size-$(d+1)$ circuits of $n$-fold products vary in at most two coordinates, 2026-09-02

**The capstone of today's Kronecker-circuit theory, and it closed a gap I
had flagged as the honest risk.** For $H=H_1\otimes\cdots\otimes H_n$ with
every spark $d_i\ge3$ and $d=\min_i d_i$: **every size-$(d+1)$ circuit
varies in at most two coordinates.** It is either a fiber over a
$(d+1)$-circuit of a single factor, or a regrouped two-factor circuit —
crossing or all-distinct — with every other coordinate fixed. No genuinely
three-dimensional $(d+1)$-circuit exists. Corollary: for **$d\ge4$** every
size-$(d+1)$ circuit is a fiber.

**Completeness of the two-factor classification, proved in general
(HUMAN-AUDITED).** The step the earlier campaigns had certified only by
census at $d=3$ is now analytic: a non-fiber two-factor circuit must use at
least $d_A$ distinct $A$-indices and $d_B$ distinct $B$-indices (dual
isolation, as in T-DGE); writing the dependency as
$A_E\,\operatorname{diag}(\gamma)\,B_E^{\mathsf T}=0$ gives
$\rho_A+\rho_B\le d+1$, hence $d_A+d_B\le d+3$, which with $d_A,d_B\ge3$
forces $d=d_A=d_B=3$. The residual profiles $(3,4)$ and $(4,3)$ are then
excluded using the two singleton rows and the one-dimensional factor-circuit
kernel, leaving exactly the crossing profile $(3,3)$ and the all-distinct
profile $(4,4)$. The $n$-factor statement follows by induction through the
composite $B=\bigotimes_{j\ge2}H_j$, with the branches $d_B=d$, $d_B=d+1$ and
$d_B>d+1$ handled separately, using T-DGE on $B$ to force a size-$d_B$
circuit of $B$ to be an axis fiber and the rank inequality to force the
all-distinct $B$-columns onto one axis line.

**Exact count, with a real double-counting trap avoided.** The number of
size-$(d+1)$ circuits is the sum over factors of fiber counts plus the sum
over unordered pairs $\{i,j\}$ of $(\prod_{k\ne i,j}s_k)\cdot
N^{\rm nonfiber}_2(H_i,H_j)$; using all-pair counts instead of non-fiber pair
counts double-counts fibers, and the planted "all-pair" route is rejected
in-run at raw 656 versus the true 624 with the 32-element duplicate excess
identified.

**Owner-independent evidence preceded and followed dispatch.** Three
$2\times3$ Vandermonde factors over $\mathrm{GF}(7)$: 27 size-3 circuits all
fibers, and **81** size-4 circuits of profiles $(1,3,3)/(3,1,3)/(3,3,1)$, 27
each, $=3$ pairs $\times$ 3 fixed indices $\times$ 9 crossings, zero
three-dimensional. Three $2\times4$ factors: **1824** size-4 circuits,
$576\times3$ crossings and $32\times3$ all-distinct, $=3\times4\times152$
exactly, zero three-dimensional. Post-freeze, a $3\times4$ MDS factor with
two $2\times4$ factors: **624** $=576$ crossings $+32$ all-distinct $+16$
fibers over the $3\times4$ factor's single 4-circuit, zero
three-dimensional. In-run: 78 exact asserts after amendment, support-set
equality at 81 (+27 fibers), 1824, 624 and 156 ($3\times4$ with
$2\times3$ and $2\times4$), every plant fired.

**Amendment, scope-only, both preregs bound.** The original prereg's P1
sentence over-scoped the frozen GRS Möbius normal form to arbitrary spark-3
factors; the general theorem needs only the abstract all-distinct two-factor
channel. The sentence was narrowed without changing theorem, pins,
instrument or promotion rule; the original prereg (commit `1b90473`, hash
`209c0f6a…`, bound in the manifest) and the amended one (commit `0cc460f`,
hash `a34cf7a7…`) are both preserved in the run directory with the
pre-amendment controls file, and `provenance.txt` records both — exactly the
rule set after the omega provenance defect.

Run `20260902T014115Z_25757c71_9c6065a4173c`, gate `H-DP1-THREEFACTOR`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`acd1866f611cb5612f46885d9ec310bf6be8fb7ff0e21770c54c6dc7f6e28303`
(9 files, owner-verified).

**State of the theory after one day.** Size $d$: T-DGE (fibers iff at most
one spark-2 factor). Size $d+1$, all sparks $\ge3$: at most two varying
coordinates; the two-factor channels are fibers (Thm-FIB), crossings (Theorem
X with H2, count $d_Ad_BC_AC_B$, field-free, present iff $d_A+d_B=d+3$) and,
at $d=3$ only, all-distinct pairs (Thm-CRIT; for $2\times n$ GRS exactly
cross-ratio preservation, counted by a $\mathrm{PGL}(2,p)$ sum). Registered
open: spark-2 boundary branches, higher-row bi-Vandermonde all-distinct
counts, and sizes $\ge d+2$.

## `mm3/` paper55's $\sigma^1$ and $\sigma^2$ orientations both exactly 55 — the ladder is complete, 2026-09-02

**Not a new record; the record value at two further orientations of the
same tensor.** The frozen $[55,57]$ intervals for `paper55` in the
$\sigma^1=(V,T(W),T(U))$ and $\sigma^2=(T(W),U,T(V))$ orientations both
collapse to $\{55\}$: $\sigma^1=14+14+27$, $\sigma^2=14+13+28$. The "best
known 57" rows are superseded by 2; $\sigma^0$'s 55 stands and nothing
sub-55 is claimed.

**Why only the output stages were open.** Gate B's frozen exact side costs
($C(U)=13$, $C(V)=14$, $C(W_{\rm fac})=14$) transport to the $\sigma$
orientations by the free-relabelling symmetry $C(T(X))=C(X)$ — a fixed
input-coordinate permutation renames free wires at zero cost, re-checked
numerically in-run for all three blocks — so each orientation's only open
quantity was whether its output stage attains its bound (27 and 28). Both do.

**Constructed and verified (MACHINE-VERIFIED).** A 14-gate $W_{\rm fac}$
circuit was **synthesized in-run** by the aux-1 dual instrument over a
hash-pinned 391-element universe (two independent enumeration loops
agreeing), admitted at instance 307 and verified by exact expansion — so no
transposition-derived witness is used anywhere, honouring the session-3
lesson; owner re-read: 391 distinct instances, three passes each, exactly
instance 307 feasible in every pass, zero instrument disagreements, and all
14 recorded gate values re-derived exactly from their operands. The printed
13-gate left and 14-gate right SLPs are reused under free renaming and
re-verified; each output stage is built by reverse-mode transposition as
explicit adjoint accumulations with additions counted exactly and checked
entry-by-entry ($23\times9$ equalities) against its orientation's
$W$-combination map; and each assembled scheme is expanded **end-to-end over
all 81 monomials** with 9/9 outputs exact, both orientations.

Controls A1–A5 accept (729/729 Brent in `fmpz` and integers for all three
orientations, printed-SLP recount $13/14/28=55$, frozen landscape rows,
relabelling symmetry, transposition preconditions, tri-instrument SAT
control) and R1–R5 reject (gate deleted, operand perturbed, frozen
$W_{\rm fac}$ floor@13 re-certified infeasible with a fresh dual-checker
DRAT→LRAT, $d$-count planter); independent audit 6/6 including a
from-scratch re-derivation of both assemblies. Two instrument defects were
caught by pre-registered checks before any claim — a relabelling that
rewrote stored values (refused by the exact-expansion verifier) and
transposed outputs emerging in wire order rather than coordinate order
(caught by the output-map equality check) — with the intermediate
INCONCLUSIVE verdict preserved verbatim. Cost about 81 CPU-s of a 2.0 CPU-h
cap. Run `20260902T015158Z_c0bea75d_761cb990fc40`, gate
`paper55-sigma12-total55`, **FROZEN-CERTIFIED**, checksum-ledger SHA-256
`70818358cf61d1a454e7f6710ccbe30157c249197c29346f276b7bad7362c5c2`
(23 files, owner-verified 23/23).

**The fixed-orientation ladder, complete.** `paper55`
$\sigma^0/\sigma^1/\sigma^2=55$ (the record, three orientations),
`sun56`$=56$ (published optimal), `mws59`$=58$ (published 59 is $+1$),
`stapleton60`$=60$ (published optimal). Every named orientation with a frozen
gap now has an exact per-orientation minimum. **Scope** unchanged; the
global sub-55 question is untouched by any campaign here, and `laderman23`
remains undecidable without verified factors.

## `rs-pe3d/` size $d+2$ classified completely for $2\times n$ GRS pairs — circuits are degree-$\le2$ bipartite graphs with minor conditions, 2026-09-02

**A complete characterization one size further, and it changes the
language.** For two $2\times n$ GRS factors ($d=3$), a 5-set of product
columns is a circuit **iff** (i) its bipartite support graph has maximum
row and column degree $\le2$, (ii) it contains no 4-edge crossing, and
(iii) every all-distinct 4-minor is nonzero. Exhaustive templates by profile:
$(3,3)$ a 4-cycle plus an edge; $(3,4)$ a 5-path plus an edge; $(3,5)$ two
same-side 3-paths plus an edge, with mirrors; $(4,4)$ a 4-path plus two
edges carrying one minor; $(4,5)/(5,4)$ a 3-path plus three edges carrying
two minors; $(5,5)$ a perfect matching carrying five minors.

**Closed forms (HUMAN-AUDITED, MACHINE-VERIFIED).** The field-free
$(1,2)/(2,1)$-decomposition channels — one $\ker A$ column on a 3-circuit
$R$ and two $\ker B$ rows at $u_1,u_2\in R$ on 3-circuits $Z_1,Z_2\ni j_0$
with both crossings cancelled — give, with injectivity proved,
$$N_{33}=9\binom{n}{3}^2,\qquad N_{34}=N_{43}=72\binom{n}{3}\binom{n}{4},\qquad
N_{35}=N_{53}=18n\binom{n}{3}\binom{n-1}{4},$$
according to $Z_1=Z_2$, $|Z_1\cap Z_2|=2$, $|Z_1\cap Z_2|=1$, and the
coincidence $N_{33}$ = the size-4 crossing count is given a structural
reason. The determinantal channels are exact sums: $N_{44}=12(M_4-N_{\rm all4})$
in terms of the 4-matchings and the all-distinct size-4 count, with explicit
$P_{45}$ and $P_{55}$ minor sums, and canonical decompositions
$(T_1,T_2,r,c)$ for every survivor — $(6,6,4,3)$ for $(4,4)$, $(6,7,4,4)$
for $(4,5)$, $(9,8,6,6)$ for $(5,5)$. The general size law is
$|S|=T_1+T_2-r-c$ with $r$ the support overlap and $c$ the cancelled
overlap; the owner's proposed $-2c$ law is only the tight full-cancellation
corollary, corrected in-run and turned into a "wrong-law" plant.

**A sharpness result inside the result.** The factor-rank criterion
$\rho_A+\rho_B\le5$ is **insufficient** for the all-distinct $(5,5)$
channel: at $n=5$, $p=13$, all 120 matchings satisfy it while only 60 are
circuits — the five 4-minors are load-bearing.

**Owner-independent evidence.** Pre-dispatch censuses: $2\times4$ over
$\mathrm{GF}(13)$ 864 with $(3,3)$:144, $(3,4)$:288, $(4,3)$:288, $(4,4)$:144;
over $\mathrm{GF}(7)$ 912 with $(4,4)$:192; $2\times5$ over $\mathrm{GF}(13)$
17880 with $(3,3)$:900, $(3,4)$:3600, $(3,5)$:900, $(4,4)$:6192,
$(4,5)$:864, $(5,5)$:60 — and the three closed forms above derived and
matched before dispatch. Post-freeze: $2\times4$ over $\mathrm{GF}(11)$ gives
**960** with $(4,4)$:240, matching the run's new anchor; every one of the 864
circuits at $p=13$ has maximum degree $\le2$ and contains no dependent
4-subset. In-run: 170 exact asserts; direct, structural and constructed
**set** equality at $p\in\{7,11,13\}$ for $n=4,5$; new totals $n=5$:
17124/17560/17880 at $p=7/11/13$ with only the determinantal profiles
varying; every plant fired including an inserted size-4 crossing and the
wrong-law route.

Run `20260902T021721Z_7629001e_215016bf40c3`, gate `H-DP2-SIZE5`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`e5d03e06999d18b67de989074955ae1e39c7b111fb594b7e271168f58567c82c` (6
files, owner-verified).

**Scope.** Two $2\times n$ GRS factors at $d=3$, size exactly 5. The
degree-$\le2$ condition is the general face of "no factor circuit inside"
($\deg\le d-1$), which suggests the general-size statement; that is the
named next gate, not a claim.

## `delcap/` q=3,n=12 box complete — 4/4 orbit-total certificates, every row beats both published endpoints, 2026-09-02

**A parameter frontier moved with replayable certificates.** The registered
$q=3$, $n=12$ extension is complete: all four cells are CERTIFIED with
conservative outward intervals that **strictly beat both** published
Tavakoli–Nguyen–Bose endpoints, exactly as the $n=11$ box did.

| $d$ | certified interval (outward, bits/symbol) | width upper |
|:---:|:---|---:|
| 1/2 | `[0.40486950906723745, 0.4048828292730583]` | `1.332e-05` |
| 1/5 | `[0.9060415197452321, 0.9060644026068896]` | `2.288e-05` |
| 1/10 | `[1.1968923784214915, 1.1968988416160002]` | `6.463e-06` |
| 1/20 | `[1.376928469567818, 1.3769335302038865]` | `5.061e-06` |

All widths inside the $1/500$ target; every dual is `ba_total_orbit_mass`;
LB$^+$ margins $\ge9.77\times10^{-2}$ / $2.10\times10^{-2}$ /
$5.29\times10^{-3}$ / $1.31\times10^{-3}$ and UB margins
$\ge3.88\times10^{-1}$ / $3.62\times10^{-1}$ / $2.30\times10^{-1}$ /
$1.29\times10^{-1}$, all as Arb lower bounds $>0$ (**MACHINE-VERIFIED**; the
exact Arb balls and outward binary rationals in `orbit_rows.jsonl` are the
certificates, not these display decimals).

**Pre-registered counts asserted, not adjusted.** 44,530 input orbits over
$3^{12}=531{,}441$ words and 66,980 output orbits over 797,161; Burnside
$B(0..12)=1,1,2,4,10,25,70,196,574,1681,5002,14884,44530$ reproduces both
(owner check: $\sum_{k\le12}B(k)=66{,}980$ and the $k\le11$ prefix
$22{,}450$ equals the frozen $n=11$ anchor); census tuple
$(44530,66980,30666848,18311678,48978526)$ identical across both static
guards and all four rows. Exact orbit-mass identities verified per row
(input numerators sum to $2^{30}$ exactly; the full-support bump fires as
$2^{30}+66{,}980$ where expected; orbit-uniform dual sums to 797,161).

**Controls in both directions, in each campaign.** ACCEPT: the frozen
$n=11$ row $0{:}3{:}11{:}1/2$ reproduced **bit-exactly** on all 35
certificate fields — census tuple, all 14,884 input and 22,450 output
numerators, every Arb ball, the whole conservative interval as exact binary
rationals, the verdict, even the BA float locator. REJECT: R1a mass total
$2^{30}+1$, R1b/R1c negative numerators, R2 negative width, R3 a planted
split output law rejected before compression while orbit-uniform is
admitted — each at its pre-registered assertion point. A 180-assertion
independent re-derivation from the frozen bytes passed; both `--validate-resume`
checks PASS.

**Two campaigns, one of them CRASHED and closed honestly.** The first run
`20260901T121346Z_d1db3155_4354a957b8ba` certified and byte-pinned rows 0–1
($d=1/2,1/5$), then died mid-row-2 with `BrokenPipeError` because its
launcher piped stdout through `tail`, deviating from the frozen command;
`FAILURE.json` is preserved, its rule forbade a rerun under the same
prereg, and it was frozen and closed **CRASHED** (13 files,
checksum-ledger SHA-256 `5c4094b5c6ad580be6a6eb81162a62df3040556637b3c80a51a589930cb1f460`).
The continuation `20260902T010546Z_adc6518e_e85083dd9888`, gate
`delcap-q3-n12-continuation-d10-d20-v1`, certified rows 2–3 under its own
prereg and closed **FROZEN-CERTIFIED** (17 files, checksum-ledger SHA-256
`1ff2b5edecf9c7bc0b4b4918e4291d98ee3131a2dc6354c98dfc904a5f918867`;
owner-verified both ledgers clean, both preregs bound). Row CPU 7664.5 s =
2.13 CPU-h of a 6 CPU-h budget; max RSS 4.40 GB of a 96 GiB cap; no
resource stop (**COMPUTATIONAL-EVIDENCE**). Two close-out checking errors by
the agent (width versus display-double difference; post-bump denominator)
are disclosed in its report.

**Scope.** A finite-$n$ theorem for exactly these four $q=3,n=12$ cells.
**No asymptotic-capacity claim.** Named next campaign: $q=3,n=13$
($B(13)=133{,}225$ input orbits over $3^{13}$; 200,205 output orbits over
2,391,484), estimated 6–12 CPU-h and 12–14 GiB, requiring its own prereg
with a larger budget and never a resume under changed caps.

## `rs-pe3d/` the two-row theory is complete; the three-row frontier is opened and honestly left open, 2026-09-02

**Certified subresult (HUMAN-AUDITED, MACHINE-VERIFIED).** Product columns
of $A\otimes B$ live in $F^{r_Ar_B}$, so every set of $r_Ar_B+1$ columns is
dependent and **no circuit has size $>r_Ar_B+1$**, with equality
attainable. For two $2\times n$ GRS factors this caps circuits at size 5,
and the classification is therefore exhaustive: size 3 fibers (T-DGE), size
4 crossing and all-distinct (Theorem X with H2, Thm-CRIT, cross-ratio,
PGL count), size 5 the six templates (H-DP2-SIZE5) — and nothing else.
Owner check: all $\binom{16}{5}=4368$ 5-sets of the $2\times4$ pair are
dependent, so size-6 circuits are empty. The degree bound is scoped
precisely: for the tied two-row layer at sizes $>3$ every circuit has
maximum degree $\le d-1$; in the unequal $(d_A,d_B)=(3,4)$ case the
non-fiber bounds are $A$-vertex degree $\le d_B-1=3$ and $B$-vertex degree
$\le d_A-1=2$, the size-4 $B$-fiber being the expected exception.

**Three-row frontier, opened.** For $(2\times n)\otimes(3\times n)$ GRS
pairs — $d_A=3$, $d_B=4$, $d=3$, ambient dimension 6, circuit sizes 3..7 —
exact censuses at $n=4,5$ over $\mathrm{GF}(7)$ and $\mathrm{GF}(13)$
passed 121 asserts, and the field-free channels were proved: the size-5
crossing at profile $(3,4)$ with count $3\cdot4\cdot C_A(3)\,C_B(4)$ (owner
check: exactly 48 size-5 circuits for the $2\times4\otimes3\times4$ pair,
all of profile $(3,4)$, $=3\cdot4\cdot4\cdot1$); three field-free size-6
families with counts $18\binom{n}{3}\binom{n}{4}$, $144\binom{n}{4}^2$ and
$180\binom{n}{5}\binom{n}{4}$; and the size-7 $K_{2,3}+K_2$ family with
$12\binom{n}{3}\binom{n}{4}$.

**Left OPEN, by verdict.** The whole gate closed **FROZEN-INCONCLUSIVE**
because the determinantal channels — size-5 $(5,5)$, size-6 $(4,5)$ and
$(5,5)$, and the residual size-7 channels — were not characterised by
explicit minor conditions in-window. Nothing incomplete was promoted; the
two-row capstone inside the run is fully proved.

Run `20260902T023120Z_8ac826db_96f766bfe6cc`, gate
`H-2ROW-COMPLETE-3ROW-OPEN`, checksum-ledger SHA-256
`d2a2177a0f8f6cf4656834918ba32dd78482481b0cd6b84435525af1824624b5` (7
files, owner-verified); target scope note and `state.json` updated.

**Where the theory stands after nine campaigns in one day.** Kronecker
circuits are now understood completely for two-row GRS pairs (sizes
$3,4,5$; field-free versus determinantal channels separated; the
determinantal ones reduced to cross-ratio and $\mathrm{PGL}(2,p)$ counting),
for $n$-fold products at size $d+1$ under sparks $\ge3$ (at most two varying
coordinates), and partially for mixed two/three-row pairs. Named next gate:
minor characterisations of the three-row determinantal channels, in the
style of the size-5 $(4,4)$ and $(5,5)$ treatment.

## `mm3/` laderman23 source-locked from the primary and its ladder closed at exactly 62, 2026-09-02

**The item blocked since session 8 is resolved, in two lifecycle
campaigns.**

**Phase 1, source lock (CITED-DEPENDENCY for the text, MACHINE-VERIFIED
for validity).** The factors come from the **primary**: Laderman, *Bull.
Amer. Math. Soc.* 82(1), 1976, 126–128, open-access publisher PDF, 236,490
bytes, SHA-256 `a1deb200f2e270b13b870bbcf8ce8c669d7a09a4acdeddc012df5d6c90042a5b`
(owner-verified), pinned with 300/400 dpi rasterisations. All 23 products
and 9 output equations were transcribed verbatim and the paper's own
criterion — "729 nonlinear algebraic equations involving 621 unknowns" —
machine-verified: 729/729 Brent over $\mathbb Z$ in pure integers and
`fmpz`, agreeing on every identity; owner re-check on
`laderman23_factors.json`: 23 products, **0/729 failures**. Page fidelity
verified twice, and an independent pinned reproduction (arXiv:1108.2830
§2.4, SHA-256 `d768596b…`) also verifies 729/729 with 9/9 identical output
supports, differing only by free per-product signs — strong corroboration
that the object is Laderman's. Defect D1: `pdftotext` truncated the
$m_3/m_{11}$ tails at the page edge and one entry was typed as $b_{33}$
where the scan prints $b_{32}$; caught by control A1 at 727/729 with no
verdict, and Amendment 1 (committed before the amended run) made the
rasterised scan authoritative. The prereg's FROZEN-NEGATIVE branch was
deliberately not invoked: the typing failed, not the source. Run
`20260902T022208Z_8c279ad1_d0276dfee63d`, gate `laderman23-source-lock`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`081df488e998a2d06974a951a1087365805ba8d76dfb04407ec3f3a349dee0a6`
(33 files, owner-verified; provenance carries both prereg hashes).

**Phase 2, the ladder (MACHINE-VERIFIED).** Laderman's printed "basic form"
recounts to exactly **98** additions $=28+28+42$ (recounted twice
independently). The certified minimum for the fixed orientation is exactly
**62** $=16+16+30$, LB $=$ UB. So the certified minimum beats the printed
count by exactly 36 — the honest framing being that it **quantifies the
slack the paper itself advertises** ("the number of additions could be
greatly reduced, but it is being given in its more basic form"): the
reducible amount is exactly 36 and 62 is the floor; no claim that Laderman
erred. Per side $C(U)=C(V)=C(W_{\rm fac})=16$ exactly: floors at 14 all
impossible, then the complete hash-pinned single-auxiliary census at $T=15$
admitted **zero of $3\times408=1{,}224$** instances, memoized subset-DFS
and a fresh kissat 4.0.4 CNF agreeing on every one, universes pinned before
any feasibility compute and re-derived by two enumeration loops, and all
1,224 infeasibilities carrying DRAT→LRAT certificates accepted by both
pinned checkers (1,224/1,224 at rc 0/0; retention policy disclosed). Upper
bounds by randomized greedy common-subexpression elimination (240 restarts
per side, 240/240 proposals verified, 0 discarded) at 16 gates per side,
each verified by exact expansion in code sharing nothing with the search;
the output stage constructed by reverse-mode transposition at exactly 30
additions and checked against all 207 $W_{\rm fac}$ map entries; the
assembled scheme expanded over all 81 monomials with 9/9 outputs exact. No
transposition-derived witness anywhere. Controls A1–A6 accept and R1–R5
reject; CaDiCaL cross-solver agreement; independent audit 8/8 with its own
re-transposition and end-to-end expansion; about 7.5 CPU-min of a 2.0 CPU-h
cap. Defect D1: control R5 was mis-specified by the agent (UNSAT exactly
because the floor is), the first run aborted with no verdict, and
Amendment 1 replaced it with a construction-known-SAT chain required to
pass in both instruments; both aborted artifacts preserved. The unrun
two-auxiliary route is quantified (~90,000 pairs per side, ~2.0 CPU-h per
side) and nothing is claimed from it. Run
`20260902T023856Z_28ff8823_43330ec64df9`, gate `laderman23-ladder-total`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`5d91d160b1bd963e6b559712fab0e9a9ce578bfd53f8ed497829b4aa44a9a28d`
(68 files, owner-verified).

**The exact fixed-orientation ladder, every named decomposition.**
`paper55` $\sigma^{0/1/2}=55$ (record), `sun56`$=56$, `mws59`$=58$,
`stapleton60`$=60$, `laderman23`$=62$ (printed 98). Laderman's 62 is above
the record and not a competitor for it. **Scope** unchanged; what remains in
`mm3/` is the global sub-55 question, which needs a genuinely different
search space and its own prereg, and the standing owner adjudication of the
void session-3 `transpose_check.py` recorded above.

## `rs-pe3d/` the three-row determinantal channels characterised — Möbius again, one row up, 2026-09-02

**The residual left open one campaign earlier is closed.** For
$(2\times n)\otimes(3\times n)$ GRS pairs the determinantal channels now
have explicit predicates and counts (HUMAN-AUDITED, MACHINE-VERIFIED as
sets):
- **size 5, profile $(5,5)$**: five all-distinct $B$-coordinates reduce,
  via Vandermonde-kernel barycentric weights, to a **unique Möbius graph**,
  so a matching is a circuit iff its pairing is a Möbius restriction and
  $N_{55}=\sum_{M\in\mathrm{PGL}(2,p)}\binom{k_M}{5}$ with
  $k_M=|X\cap M^{-1}(Y)|$ — the same law as the two-row case, one size up.
  Owner check: the $(2\times5)\otimes(3\times5)$ pair over $\mathrm{GF}(13)$
  has exactly 6 such circuits and the PGL sum gives 6;
- **size 6, one repeated $B$-coordinate**: the singular determinant
  factors as $(-1)^{b-a-1}(x_b-x_a)\prod(y_i-y_0)\,D_4(R)$ on the four
  singleton edges, giving $N_{45}=6(|Y|-4)\sum_M\binom{k_M}{4}$ and
  $N_{55}=4\sum_M\binom{k_M}{4}\big((|X|-4)(|Y|-4)-(k_M-4)\big)$ with a
  central-deletion exclusion;
- **size 7**: fully counted at $q=4$; for $q=5,6,7$ exact necessary-and-
  sufficient deletion-minor lists.
Direct and predicted populations agree as **sets** — overall, by profile,
and across all 17 size-7 signatures — for $n=4,5$ and $p\in\{7,11,13\}$;
261 exact asserts; 83 CPU-s. Two instrument defects were preserved,
disclosed and fixed with all families rerun.

Run `20260902T025427Z_1b1f202c_6247a0994bce`, gate `H-3ROW-DETERMINANTAL`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`a564b92bfa62440465cfc3c9f212a88c98fc88df42e69833792d4f113f871843`
(10 files, owner-verified).

**Scope and next gate.** Two-row and mixed two/three-row GRS pairs are now
classified with the field-free/determinantal separation explicit at every
size, and the recurring principle is that all-distinct channels are Möbius
restrictions counted through $\mathrm{PGL}(2,p)$. Named next gate, not a
claim: $(3\times n)\otimes(3\times n)$ pairs ($d=4$, ambient dimension 9,
circuit sizes 4..10), where a second-order analogue of the Möbius law is
expected but unproved.

## `rs-pe3d/` tied three-row pairs: sizes 4–7 settled, the Möbius law holds at size 6, size 8 is the first open layer, 2026-09-02

**Layer-by-layer, promoted only where proved.** For $(3\times n)\otimes(3\times n)$
GRS pairs ($d_A=d_B=d=4$, ambient dimension 9, circuit sizes 4..10):
- **size 4**: fibers only, $N_4=2n\binom{n}{4}$ (T-DGE);
- **size 5 $=d+1$**: **empty** — exactly Thm-CRIT's prediction for a tied
  minimum with $d\ge4$, now confirmed in the first configuration where the
  crossing channel sits strictly above $d+1$;
- **size 6**: **complete** — crossings $16\,C_A(4)\,C_B(4)$ plus the
  all-distinct channel, which is again a **unique Möbius graph** with count
  $\sum_{M\in\mathrm{PGL}(2,p)}\binom{k_M}{6}$, every other profile excluded.
  The analytic route: for six distinct $y$, the Vandermonde-kernel form gives
  $q_0,q_1,q_2$ of degree $\le2$ with $q_1^2=q_0q_2$, forcing $x=q_1/q_0$ to
  be linear-fractional; the converse restricts the bidegree-$(2,2)$ space to
  a 4-dimensional one, so any five of the six points are independent;
- **size 7**: **empty**;
- **size 8**: an exact injective reduced rational-degree-2 subfamily is
  proved, but a **residual** remains — dedicated full-matching sweeps at
  $n=8$ leave 560 and 416 unexplained circuits at $\mathrm{GF}(11)$ and
  $\mathrm{GF}(13)$ — so this is the first unresolved structural layer;
- **sizes 9 and 10**: exact curve-plus-deletion and deletion-determinant
  predicates.
Owner check on the $(3\times4)\otimes(3\times4)$ pair over $\mathrm{GF}(13)$:
size-5 circuits $=0$, size-6 circuits $=16$, all of profile $(4,4)$,
$=16\binom{4}{4}^2$ crossings with the all-distinct channel structurally
empty at $n=4$ (six distinct indices cannot fit in four).

In-run: the final corrected full matrix passed 143 assertions in 618 CPU-s;
a first pass and a cancelled partial are preserved and disclosed with the
complete matrix rerun; independent math and compliance audits found no
missing obligations. Run `20260902T031527Z_2d19c5c4_718416983678`, gate
`H-33ROW-OPEN`, **FROZEN-INCONCLUSIVE** for the whole gate with sizes 4–7
certified inside, checksum-ledger SHA-256
`43707d2c02a25f24aa46ff7e2440c7962a7a25f88f52144bd011a30bece7b53c`
(11 files, owner-verified).

**Recurring principle, now seen four times.** For GRS factor pairs, the
first all-distinct channel at each row-configuration is a Möbius
restriction counted through $\mathrm{PGL}(2,p)$: two-row at size 4, mixed
two/three-row at size 5, tied three-row at size 6. Named next gate: the
size-8 residual.

## `rs-pe3d/` the size-8 residual is a complete intersection — the second-order Möbius law, 2026-09-02

**The last open layer of the tied three-row theory is explained, and the
explanation is geometry.** For $(3\times n)\otimes(3\times n)$ GRS pairs,
an **all-distinct size-8 circuit** is exactly one of two things: an
exact-degree-2 graph of class $(1,2)$ or $(2,1)$, or a **reduced complete
intersection of two bidegree-$(2,2)$ curves** — eight points whose pencil of
$(2,2)$-curves is two-dimensional with **no fixed component**, and whose two
normalised resultants equal the 8-node projection polynomials. Bézout on the
bidegree gives $(2,2)\cdot(2,2)=8$, so there is no ninth point. On the
$n=8$ point sets both graph orientations are empty, so the residual —
**560 at $\mathrm{GF}(11)$ and 416 at $\mathrm{GF}(13)$** — equals the
complete-intersection family as exact support sets.

**Why rank alone is not the predicate.** "Evaluation matrix of rank 7"
over-accepts: 2192 and 1344 pencils, of which exactly 1632 and 928 are the
common-$(1,1)$-component sets — pencils whose base locus contains a
six-point Möbius subcircuit. The size-6 Möbius law and the size-8
complete-intersection law are thus the first- and second-order rungs of one
ladder: a $(1,1)$-curve through six points, a pencil of $(2,2)$-curves
through eight.

In-run: 72 controls in 18.9 CPU-s, including a wrong-degree predicate, an
inserted lower-size circuit, the ambient plant and the concat guard; a
pre-execution syntax error was preserved and disclosed with the complete
family rerun; independent proof and compliance audits found no substantive
issue. Run `20260902T041030Z_64e0c4ef_b5de16e49f8a`, gate
`H-33ROW-SIZE8-RESIDUAL`, **FROZEN-CERTIFIED**, checksum-ledger SHA-256
`bbd14c3bea49a1a1850eaf0a0dae4889049d72f6e7fc5416aefb3ffaabd836cd`
(10 files, owner-verified; anchors 560/416/2192/1344/1632/928 present).

**Scope.** The all-distinct size-8 gate of the tied three-row pair;
repeated-index size-8 supports are outside it. **The target rests** for
this session after twelve campaigns: eleven certified statements and one
honestly inconclusive gate, from the falsification of H-MINLINE at
$q=13,s=(2,2,4)$ this morning to complete intersections of $(2,2)$-curves
tonight.

## `mceliece/` third $m=12$ row: $t=48$ — the ladder's last unreached cell, $\alpha$ measured, 2026-09-02

**Three bounded instances now sit at $m=12$: $t\in\{48,64,96\}$, one seed
each.** At $(m,n,t,\text{seed})=(12,3488,48,16384)$ — $k=2912$, $D=3391$,
$2t+3=99$ (owner-recomputed) — the Guard-A construction passes with the same
standard as the two earlier rows: $G$ monic irreducible of degree 48 (17
seeded tries, scalar-Rabin), square-free; **$\alpha$ MEASURED** on all 3488
support points with 0 failures (798.70 CPU-s), the fork resolved at the
registered gate on measured cumulative CPU $1966.59\ll24300$ so the DERIVE
branch was never entered; $\delta$ chain exact 5/5 with the degree gate
closing at $\max_j\deg f_j=3391=D$ exactly; controls, anchors 3/3
byte-equal against the frozen $m\le11$ records, and plants 4/4 rejected
including a $t{=}48$-native probe failing 66/128 corrupted against 0 clean.
Cost 2781.02 of a 32400 CPU-s budget. Run
`20260902T014243Z_340fc591_5101c05ab359`, **FROZEN-CERTIFIED**,
checksum-ledger SHA-256
`a096403d55136cad7e4a28fa8c11cfeee555e6ccabeca7dd2a63aa228978a78d`
(24 files, owner-verified 24/24; `verdict.json`
`alpha_identity_measured_full_grid: true`, `alpha_failures: []`,
`alpha_grid_points_checked: 3488`).

**Process note, disclosed in `RLIMIT_DIAGNOSIS.md`.** Two silent deaths of
the verdict stage at exactly 7200 s CPU were SIGXCPU from a soft
`RLIMIT_CPU` in tool-spawned contexts, not memory pressure as first
suspected; the fix was an exec wrapper raising the soft limit to the hard
limit, with the stage code unchanged and no amendment needed.

**Scope.** Still **not** "$m=12$ complete"; no $m=12$ waterfall census rows
exist; the Apon section-3.6 hole remains **bounded, not closed**; no
NIST-cost claim. Owner decision: the target stands down for this session with
the three bounded instances; the named next campaign, if any, is the $m=12$
waterfall census under its own prereg and a budget sized from the measured
rates recorded in these three runs.

## `rs-pe3d/` the general-size theorem — every Möbius law is one instance, 2026-09-02

**Unification, certified.** For GRS factors with distinct evaluation points,
let $V=F[x]_{<r_A}\otimes F[y]_{<r_B}$, of dimension $r_Ar_B$. Because the
product column of $(u,v)$ is the evaluation vector $(x_u^iy_v^j)_{i,j}$, an
all-distinct $k$-set $S$ of grid points is a circuit of $A\otimes B$ **iff**
$S$ imposes exactly $k-1$ independent conditions on $V$ and every proper
subset imposes independent conditions — the circuit definition transported
through the evaluation map (HUMAN-AUDITED). The corollaries carry the
content:
- **C1, the unified Möbius law.** Points on the graph of a Möbius map lie on
  the $(1,1)$-curve $c_3xy+c_2x+c_1y+c_0=0$, whose multiples inside $V$ form
  a space of dimension $(r_A-1)(r_B-1)$, so any $r_A+r_B$ such points impose
  at most $r_A+r_B-1$ conditions and are dependent. The first all-distinct
  circuits therefore appear at size $r_A+r_B$ and, on every swept cell, are
  exactly Möbius restrictions counted by
  $\sum_{M\in\mathrm{PGL}(2,p)}\binom{k_M}{r_A+r_B}$ — sizes 4 ($2{\times}2$
  rows), 5 ($2{\times}3$), 6 ($3{\times}3$ and $2{\times}4$), 7 ($3{\times}4$) are
  one statement;
- **C2, top sizes.** At $k=r_Ar_B$ the points lie on a unique curve of $V$
  with no special proper subset; at $k=r_Ar_B+1$ every set is dependent and
  circuits are the sets all of whose $r_Ar_B$-subsets are independent —
  the two-row size-5 and tied three-row size-9/10 predicates as instances;
- **C3, pencils.** At $k=r_Ar_B-1$ the vanishing space is a pencil; the
  tied three-row size-8 complete-intersection law is the instance $r_A=r_B=3$.

**Evidence (MACHINE-VERIFIED, set equality at all 15 anchor cells).**
Owner-derived anchors before dispatch, over $\mathrm{GF}(13)$ with points
$1..n$: $(2\times n)\otimes(4\times n)$ has no all-distinct size-5 dependency
at $n=6$, and its size-6 circuits number 2 at $n=6$ and 42 at $n=7$, all
Möbius, equal to the PGL sums; $(3\times7)\otimes(4\times7)$ has no size-6
all-distinct dependency and exactly 2 size-7 circuits, both Möbius, PGL sum 2.
In-run: these plus $p=11$ and $p=7$ (588 at $n=7$, $k=6$, where the residues
of $1..7$ fold to $\{0,\dots,6\}$ — an expectation first pinned from a
wrong-residue probe was corrected independently through the PGL sum and
disclosed), the tied three-row re-expressions $12/4/2$ and $588/72/42$, C2 on
$4\times4$ at sizes 9/10 $=0/16$ over three primes and on $5\times5$ over
$\mathrm{GF}(7)$ at $36378/181960$, and the size-8 bootstrap
$560/416$ complete intersections against $2192/1344$ rank-7 pencils with
$1632/928$ common-$(1,1)$-component sets. Owner post-freeze check at a new
prime: $(2\times7)\otimes(4\times7)$ over $\mathrm{GF}(11)$ has 72
all-distinct size-6 circuits, all Möbius, PGL sum 72. Controls: ambient
plant discriminating rank $=r_Ar_B$ exactly, a planted non-Möbius
$(r_A{+}r_B)$-set that must be independent, wrong-coefficient and
selected-pole replays, corrupted column, duplicated-column spark drop, concat
guard — every accept/reject pair fired. 103 asserts, 1627.5 CPU-s of a 5400 s
cap under `RLIMIT_CPU` raised to the cap. Five pre-freeze launch aborts and
one post-run control-audit gap were disclosed in `defect_log.json` and
remediated by a full second complete run; the prereg was never amended.

Run `20260902T113614Z_a999d1b0_5daa598b6a04`, gate `H-GRS-GENERAL-SIZE`,
**FROZEN-CERTIFIED**, checksum-ledger SHA-256
`dc9d41031de17e82951413821eb69fa46794b4b0851e59ae7b6fe0086c304f77`
(13 files, owner-verified).

**Scope, with the honest leads recorded in `state.json`.** The forward
direction of C1 (Möbius $\Rightarrow$ dependent) and all set equalities are
proved and verified; the **converse** of C1 (dependent at size $r_A+r_B$
$\Rightarrow$ Möbius) is proved only for the swept row-configurations, the
emptiness of all-distinct circuits below $r_A+r_B$ is a slice boundary rather
than a general theorem, and C3's no-fixed-component characterisation is
stated for $r_A=r_B=3$ only. Named next layer: $(4\times n)\otimes(4\times n)$.

## `delcap/` q=3,n=13 box complete — 4/4 rows beat both published endpoints, one clean run, 2026-09-02

**The frontier moves one more step, and this time without a crash.** All
four $q=3,n=13$ cells are CERTIFIED in a single run with the budget not
binding, every conservative outward interval strictly beating **both**
Tavakoli–Nguyen–Bose endpoints.

| $d$ | certified interval (outward, bits/symbol) | width upper |
|:---:|:---|---:|
| 1/2 | `[0.39713666568100003, 0.397161636764393]` | `2.497e-05` |
| 1/5 | `[0.8966732662961678, 0.8967669676702715]` | `9.370e-05` |
| 1/10 | `[1.1898717660767768, 1.1898901620041096]` | `1.840e-05` |
| 1/20 | `[1.3725688358567392, 1.3725840117520043]` | `1.518e-05` |

All widths inside $1/500$, all duals `ba_total_orbit_mass`; LB$^+$ margins
$\ge9.86\times10^{-2}$ / $2.13\times10^{-2}$ / $5.37\times10^{-3}$ /
$1.33\times10^{-3}$ and UB margins $\ge3.95\times10^{-1}$ /
$3.71\times10^{-1}$ / $2.37\times10^{-1}$ / $1.33\times10^{-1}$
(**MACHINE-VERIFIED**). Owner re-check from the frozen `orbit_rows.jsonl`:
for every row the display decimals are outward of the exact binary
rationals and the exact width is $\le1/500$.

**Counts asserted, not adjusted.** $B(13)=133{,}225$ input orbits over
$3^{13}=1{,}594{,}323$ and 200,205 output orbits over 2,391,484 (owner:
Burnside table extended by one term reproduces both); the $k\le11$ and
$k\le12$ prefixes reproduce the frozen 22,450 and 66,980; dense census
$(133225,200205,152963378,93991184,246954562)$ byte-identical across the
static guard, both production builds and all rows. Per-row orbit-mass
identities hold (input $2^{30}$; BA post-bump $2^{30}+200{,}205$; uniform
2,391,484). Controls: ACCEPT the frozen $n=12$ row $0{:}3{:}12{:}1/10$
bit-exactly on all 41 certificate fields; REJECT R1a/R1b/R1c/R2/R3 all fired.
An AST diff proves the 22 certificate functions byte-identical to the $n=12$
instrument; the 180-assertion re-derivation passed.

Cost 38,668 s CPU = 10.74 of 24 CPU-h, wall 11.28 of 26 h, maximum row RSS
16.67 GB of 32 GiB (**COMPUTATIONAL-EVIDENCE**). Run
`20260902T024656Z_90c757de_c4e9905d05a1`, gate
`delcap-q3-n13-total-output-orbit-mass-v1`, **FROZEN-CERTIFIED**,
checksum-ledger SHA-256
`b362203bca3acf9fdf01aecbd5c75c7f5a33a5709f76ec95177bef1c2d0be6d5`
(14 files, owner-verified 14/14).

**Scope, and why the family stops here.** A finite-$n$ theorem for exactly
these four cells; **no asymptotic-capacity claim.** $n=14$ was evaluated
and **not opened**: measured per-row cost scales $287.6\to1916.1\to9573.8$ s
across $n=11,12,13$, projecting ~13.9 h per row, ~55.5 CPU-h per box and a
35–40 GB live peak — above the 32 GiB envelope this instrument was
registered under. The targets are recorded ($B(14)=399{,}310$ input orbits,
599,515 output) for a future fresh prereg with ~64 GiB / ~64 CPU-h caps or a
separately frozen instrument change. The family rests at a complete
certified box.

## `rs-pe3d/` the tied four-row layer confirms every prediction of the general theorem, 2026-09-02

**A prediction test, and the theory passed it.** For $(4\times n)\otimes(4\times n)$
GRS pairs ($d=5$, ambient dimension 16, circuit sizes 5..17), each layer in
scope was predicted from the certified theorems before compute and then
certified (HUMAN-AUDITED conditions, MACHINE-VERIFIED sets):
- **size 5**: circuits are exactly the fibers, $2n\binom{n}{5}$ — 10 at
  $n=5$, 72 at $n=6$ — as support-set equality over both primes (owner
  census: exactly 10 size-5 circuits, all fibers);
- **size 6 $=d+1$**: **empty** at $n=5,6$ — Thm-CRIT's tied channel needs
  $d_A+d_B-2=8\ne6$ and fiber lifts are circuits only at $k=\text{spark}=5$;
  the only dependent 6-sets are fiber-plus-lone (200 and 2160), each
  rejected by a deletion witness;
- **size 7**: **empty** at $n=5,6$, by a full $\binom{36}{7}=8.35$M sweep in
  which all 31,320 dependent sets contain a fiber;
- **size 8**: at $n=5$ the exhaustive 1.08M census finds exactly **25**
  circuits $=25\binom{5}{5}^2$ crossings and nothing else; at $n=6$ all 900
  constructed crossings verify as circuits; and the first all-distinct
  circuits appear at $r_A+r_B=8$ as **Möbius restrictions**, equal as sets
  to $\sum_M\binom{k_M}{8}$ — 8 over $\mathrm{GF}(11)$ and 4 over
  $\mathrm{GF}(13)$ (owner PGL sums: 8 and 4).
Controls: a $4\times4$ sub-grid tensor basis with rank exactly 16, a 17-set
dependent at rank 16, a 15-plus-duplicate 16-set of corank 1 with singular
deletion rejected, corrupted column, concat guard (16 versus 7), crossing
circuit and Möbius-8 restriction-minimality — every pair fired. 49 asserts,
2463.8 CPU-s of a 5400 s cap. Defects disclosed: one control-plant slip
fixed pre-launch, and one SIGXCPU kill when the first launch enumerated
$\binom{64}{8}\approx4.4\times10^9$ tuples in the $n=8$ sweep, replaced by
the $8!=40320$ bijection enumeration and rerun from the beginning; the
prereg was never amended.

Run `20260902T135324Z_516630ab_ed03d008f0f6`, gate `H-44ROW-LAYER`,
**FROZEN-CERTIFIED** within its registered scope, checksum-ledger SHA-256
`1dc80d247d61c5e8d1da69abb37ee615faca0fca0a293c326f07f96601e79445`
(10 files, owner-verified).

**Scope.** Tied size-8 exhaustiveness at $n\ge6$ (30.3M subsets) and sizes
9–17 are explicit non-claims. The `rs-pe3d/` target rests for this session
after **fourteen** campaigns.
