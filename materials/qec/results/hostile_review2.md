# Hostile review — round 2 (post-repair shipped candidate)

**Date:** 2026-08-18 · **Target:** `math/qec/reports/paper_pbb_nogo.md` (584 ln, edited 2026-08-18 15:43) and the shipped `paper_pbb_nogo.tex`/`.pdf` (15:47) · **Reviewer:** HostileReview2
**Stance:** adversarial; every check below was re-derived from the artifacts or recomputed in GF(2) numpy (≤20 min single-threaded total). No SAT/stim-scale reruns.

## Verdict

**NO_CRITICAL_FINDING.** No theorem, headline number, or table cell is wrong: the six 2026-08-18 edits are directionally sound and every checked numerical claim matches its artifact (see §Rerun table). However, the shipped candidate is **not submittable as-is**:

- **2 HIGH** — a false justification sentence inside the proof of Theorem G(ii) (the repaired proof replaced one false sentence with another; the identity itself is fine), and a visible text corruption in the shipped PDF (`extbfThe X sector.`).
- **3 MEDIUM** — scope/coverage misdescriptions in the two new §6 paragraphs and a Conjecture-B′ naming collision.
- **5 LOW** — wording/staleness nits.

---

## Findings

### HIGH-1 — Theorem G(ii)'s replacement proof contains a false statement ("row space extends")

**Quote** (md 131–133; tex 173–176):
> "For the dimension we count generators directly: the row space of $H_Q$ extends that of $H_P$ exactly by the $[C\;D]$-parts of first-block combinations taken modulo $S_Z$, so $\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\bar\Delta$ …"

**Why wrong.** `rowspace(H_Q) ⊇ rowspace(H_P)` is **false in general**, and it fails on the paper's own showcase rows. Recomputed exactly (GF(2), this session):

| row | rank H_Q | rank H_P | rank[H_Q;H_P] | rowspace(H_P) ⊆ rowspace(H_Q)? | parent row (e₀[A B]\|0) ∈ rowspace(H_Q)? |
|---|---|---|---|---|---|
| `phase2_58` | 68 | 64 | **92** | **no** | **no** |
| `12_6_0193` | 132 | 132 | **178** | **no** | **no** |

The parent's first-block stabilizer rows are *destroyed*, not extended — this is exactly the "demoted parent stabilizers" phenomenon the paper itself quotes in §6 (J.2/J.4; the X-sector note records 255 catalogue rows with nonzero demotion space). The sentence is the same defect class as the pre-repair "X-sector of the quotient is unchanged" sentence it replaced.

**What survives (verified):** the *conclusion* is correct — `rank H_Q = rank H_P + dim bar Δ` holds on both rows (= `rank[A;B] + dim(S_Z+Δ)`: 68 = 32+36, 132 = 132), `k_Q = k_P − dim bar Δ` holds (phase2_58: k_Q=4=k_P−4, k_P=8; 12_6_0193: k_Q=12=k_P−0 — both match catalogue k), and the sector-identity derivation in the same paragraph ("Equivalently, via the sector identity …") is complete and correct on its own, as is the 368/368 machine check. So Theorem G(ii) stands; only the "count generators directly" sentence is wrong.

(The parent's rowspace is *not* contained in $Q$'s — parent X-stabilizers are demoted — so any phrasing implying containment must go.)

### HIGH-2 — Shipped PDF renders corrupted text `extbfThe X sector.`

**Quote** (`reports/paper_pbb_nogo.tex:525`, the only tab in the file):
```
	extbf{The $X$ sector.} The full analysis is in \texttt{notes/theorem\_j\_xsector.md};
```
(`\textbf` lost its `\t` to a literal TAB.)

**Why wrong.** Confirmed in the shipped PDF itself (decompressed content stream): the paragraph heading typesets as literal body text **"extbfThe X sector."** — visible garbage at the start of the new §6 paragraph in every reader's first pass. The `.md` heading (`**The X sector.**`, md 414) is fine; only tex/pdf are affected.

**Exact fix.** `reports/paper_pbb_nogo.tex:525` → `\textbf{The $X$ sector.}` and rebuild the PDF. (Grep confirms this is the only tab-corrupted macro in the tex.)

### MEDIUM-1 — Residual paragraph misdescribes EXP-042's enumeration scope and quotient

**Quote** (md 426–428; tex 536–539):
> "among all $4{,}621$ distinct connected trinomial $k_P\ge2$ parents with $\ell\cdot m\le 24$ (exhaustive enumeration, quotiented by the $\mathbb Z_\ell\times\mathbb Z_m$ translation symmetry and the Hadamard $\ell\leftrightarrow m$ swap), **none** has $T<k_P/2$."

**Why wrong (two mismatches vs `results/processed/exp042_residual_probe.json`).**
1. **"trinomial" is the wrong word.** The artifact's `support_families` are unordered support-size pairs **(2,2), (2,3), (3,3)** — i.e. $A,B$ each of weight 2–3. Per-family counts recomputed from `lattices[].families`: 2×2 → 855, 2×3 → 1,012, 3×3 → **2,754**. Strictly trinomial (3×3) parents number 2,754, not 4,621; the quoted 4,621 is the union including 1,867 binomial parents. (Codebase usage agrees: `exp040_saturation_probe.py` calls (3,3) "trinomial" and (2,2) the "binomial fallback".)
2. **The stated quotient is not the one used.** `protocol.parent_equivalence` = "supports identified under **independent** $\mathbb Z_\ell\times\mathbb Z_m$ translations of the A and B blocks; $\{A,B\}$ **unordered**" — not a single diagonal translation, and not an $\ell\leftrightarrow m$ Hadamard identification (all **37 ordered** $(\ell,m)$ shapes were scanned separately; `shapes: 37`). The count 4,621 is the fingerprint count under the *artifact's* equivalence, which is coarser than the paper's parenthetical implies.

**What survives (verified):** `verdict = SMALL_CLASS_EMPTY`, `distinct_parent_fingerprints = 4621`, `bucket_T_lt_half_k = 0`, `residual_parents = []`, `unresolved_parents = 0` on all 37 shapes, `upper_law_breaks = 0`; buckets T=k/2: 3,908, T=k: 713 (3,908+713=4,621). Since the scanned set ⊇ any reading of "trinomial", the substantive claim (no parent with T<k_P/2 at ℓm≤24) is true and, correctly, the paragraph does not claim anything beyond ℓm≤24. This is a scope/ provenance misdescription, not an overclaim of strength.

**Exact fix** (md 426–428; tex 536–539): "among all $4{,}621$ distinct connected $k_P\ge2$ parents with $A,B$ of weight $\le 3$ each (support families $2{\times}2,2{\times}3,3{\times}3$; exhaustive over all 37 ordered lattices $\ell,m\ge2$, $\ell\cdot m\le24$, quotiented by independent block translations with $\{A,B\}$ unordered) …".

### MEDIUM-2 — X-sector paragraph overstates catalogue coverage of the $d_X$ check

**Quote** (md 419–421; tex 530–532):
> "zero violations in $3.16$M valid perturbations across $12$ lattices and on all $368$ catalogue rows"

**Why wrong.** The catalogue comparison $d_X(Q)\ge d_X(P)$ is only fully determined on the rows where **both** quantities exist. Recomputed from `results/partial_runs/xsector/catalogue_scan.json` (368 records): $d_X(Q)$ exact on **276** (gated by $\dim\mathrm{Xcen}(Q)\le 22$), certified $d_Z(P)=d_X(P)$ on **279**, both on **221** (17 at equality, **0 violations**). On the remaining 147 rows no violation could have been checked; "on all 368 catalogue rows" asserts a coverage the anchor does not have. The companion note states this precisely ("221 rows carry both quantities"); the paper compressed it incorrectly. (The 12-lattice/3.16M half is correct: hunt specs union to exactly 12 distinct $(\ell,m)$ pairs, and 3,228+1,029,827+2,131,004+368 = 3,164,427 ≈ 3.16M.)

**Exact fix** (md 419–421; tex 530–532): "… and on the $221$ catalogue rows where both $d_X(Q)$ (exact) and $d_X(P)=d_Z(P)$ (certified) are known ($17$ at equality, $0$ violations); $d_X(Q)$ is exact on $276$ of $368$ rows in total."

### MEDIUM-3 — Two different statements named "Conjecture B′"; §6's parenthetical is self-contradictory

**Quotes.** §6 (md 430):
> "Whether Conjecture B′ (all residual parents satisfy $T\ge k_P/2$) holds beyond $\ell\cdot m\le 24$, or whether a second exception exists at larger lattices, is under computation."

§7 (md 482–483):
> "**Conjecture B′ (residual saturation).** For parents with $T<k_P/2$ (one known: `9a7638586033` …), distance increase still forces $\dim\bar\Delta=T$."

**Why wrong.** (a) "Residual parents" are by definition the $T<k_P/2$ parents, so "all residual parents satisfy $T\ge k_P/2$" is vacuous — no residual parent can satisfy it. The intended §6 statement is evidently "the residual class is empty beyond the known exception" (i.e. no second exception), which is a *different* conjecture from §7's B′ (which presupposes residual parents exist and asserts saturation on them). (b) Same label, two statements: EXP-042's `b_prime_distance_increase_off_t` counterexample key tracks §7's meaning. (c) Additionally, the shipped **tex omits the §6 sentence entirely** (the tex paragraph ends at "far outside this scope.") — an md↔tex divergence introduced by this edit round, so md and PDF currently make different conjecture claims.

**Exact fix.** In md §6 rewrite to: "Whether the residual class (parents with $T<k_P/2$) is empty beyond `9a7638586033` at larger lattices is under computation; Conjecture B′ (§7) governs what happens on any residual parent that does exist." Then re-sync tex/md so both carry the corrected sentence.

### LOW-1 — "both counts drift monotonically with the running sweep" is false for the escapee count

**Quote** (md 438): "… both counts drift monotonically with the running sweep."
**Why wrong.** Escapee history: 21 → 35 → 30 (claims_matrix.md R8 records 84 pending/21 escapees pre-revision and 60/35 at its check; the live artifact now has 60/30). Pending is monotone (84→60→60); escapees are not — an escapee becomes dominated whenever the CSS pool gains a certificate.
**Exact fix.** "…; the pending count only decreases, and the escapee count shrinks as the CSS pool grows — both are provisional."

### LOW-2 — §5 summary table's "30" lacks the "of 279 classified" qualifier

**Quote** (md 306): "| catalogue rows Theorem H permits (candidates for a genuine increase) | $30$ |"
**Why wrong.** Across all 368 rows, 89 sit on uncertified parents ($T$ a lower bound) and are *unclassified*, not "permitted"; 249+30=279≠368, while the neighbouring row says "249 of 368". §5.1 has the correct statement ("exactly 30 of 279").
**Exact fix.** "| catalogue rows Theorem H permits, of the 279 classified | $30$ |".

### LOW-3 — [5,12] sentence: two wording compressions (substance verified)

**Quote** (md 441–443): "weight-$4$ meet-in-the-middle exclusion completed on both circuits (no weight-$\le4$ single-mechanism logical failure; structural-probability DEM)".
**Why imprecise.** (a) The weight-4 exclusion is *two* exhaustive methods — a pair-pair XOR MITM over the perfect-matching class and a per-centre exact-cover over the star $K_{1,3}$ class (`exp040_circuit_distance.json`, `w4_partition_argument`); only the first is meet-in-the-middle. (b) "no weight-≤4 single-mechanism logical failure" reads as if only single mechanisms were excluded; the exclusion covers **all** combinations of ≤4 DEM mechanisms (w1/w2 singles, w3 pairs, w4 triples/quads). The word "certified" is **not** an overstatement: both circuits carry `lb{value:5, certified:true}`, w3/w4-matching/w4-star all `completed`+`certified` with no timeouts, and `ub{value:12}` from EXP-029's certified weight-12 data-logical witnesses on both — "tied at the certified interval [5,12] on each" is exactly what the artifact supports (top-level verdict additionally certifies `d_DEM_mech ≥ 5` on both).
**Exact fix.** "weight-$4$ exclusion completed on both circuits (perfect-matching class by meet-in-the-middle, star class by exact cover; no logical failure from any set of $\le4$ DEM mechanisms; structural-probability DEM)".

### LOW-4 — Environment line's suite count is stale by one

**Quote** (md 512; tex 627): "Full suite: $897$ passing tests, $1$ skipped."
**Why wrong.** `pytest --collect-only` over `tests/` now collects **899** tests; 897+1 accounts for 898. All 748 tests in the three non-SAT-heavy files (`test_paper_claims`, `test_sector_identity`, `test_dimension_identity_368`) pass this session. Either the count predates a test added on 2026-08-18 (`test_sector_identity.py`, 01:50) or one test currently fails — a fresh full-suite run is needed before quoting a number. (All package versions in the line verified live: Python 3.13.9, numpy 2.4.6, scipy 1.18.0, stim 1.16.0, pymatching 2.4.0, sinter 1.16.0, ldpc 2.4.1, galois 0.4.11, ortools 9.15.6755, python-sat 1.9.dev13 — all exact matches.)
**Exact fix.** Re-run the suite, update to the observed pass/skip counts.

### LOW-5 — "(ii) The $n=360$ parent is not yet closed" understates the gap

**Quote** (md 438). §5's own table says 39 parents at n=360 with only 1/39 certified and 0 closed; 38 are not even certified. The singular suggests one outstanding parent.
**Exact fix.** "(ii) The $n=360$ parents remain open (38 of 39 not yet certified)."

---

## Rerun / verification table (everything recomputed this session)

| # | Item | Source / method | Result vs paper |
|---|---|---|---|
| 1 | rowspace(H_Q) vs rowspace(H_P); rank identity | GF(2) numpy rebuild of `phase2_58`, `12_6_0193` | rank identity ✓ (68=64+4; 132=132+0; = rank[A B]+dim(S_Z+Δ)); **inclusion false** (HIGH-1) |
| 2 | 134 certs: exact/closed/T≥k_P/2, exception, per-n | all `exp039/parent_*.json` | 134/134 exact; 61 closed; 133 ≥ k_P/2; exception `9a7638586033` (n=144, k=12, T=4, d_Z=6) ✓; by n {36:3,72:8,108:68,144:38,180:16,360:1} ✓; closed by n {3,5,28,21,4,0} ✓; slowest 16,492.95 s ✓ |
| 3 | 11 catalogue [[144,12,12]] parents all T=12 closed | certs filtered n=144,k=12,d_Z=12 | 11, all T=12, all family_closed ✓; `4c7eb964` 220.43 s, 3 witnesses ✓; Gross `aebb6495` T=12 exact closed ✓ |
| 4 | §5.1 cross-tab + spread + capped-by-n | `exp039_nogo_module.json` rows | cells {(1,1):109,(1,½):67,(1,⅓):1,(⅚,1):11,(⅚,½):10,(¾,1):25,(¾,½):13,(½,1):13,(½,½):30} ✓; spread −24…0 ✓; capped by n {47,32,80,65,24,1} ✓; counts 202/134/61/249/30 ✓ |
| 5 | Escapees/pending | `exp037_envelope_classification.json` | {DOMINATED:278, ESCAPEE:30, PENDING:60}, sum 368 ✓; 7 reversals all DOMINATED with U ∈ {10,10,6,6,6,6,6} ✓ |
| 6 | §6 domination table | `exp036_envelope_check.json` | 7/7 dominated; 162 candidates ✓; kd²/n CSS {12.00, 7.41, 6.00} vs PBB {2.78, 1.85, 2.00, 0.67} ✓; factors 4.32/4.00/3.70/11.11/3.00 ✓; strongest-by-kd² at 72/108/144 = [[72,12,6]]/[[108,8,10]]/[[144,12,12]] ✓; 5/7 parent-profitable, max gain 1.125≤1.13 ✓ |
| 7 | §5.2 table / §7 item 2 | `exp040_saturation_probe.json` | 64 parents (16×4 lattices); 26,898 δ>0; <T/=T/>T = 23,634/3,264/0 ✓; 496 increases all at =T (counterexamples empty) ✓; all 64 T∈{k/2,k} ✓; containment verified (760 contained); 7/7 dim bar Δ=T, hypothesis certified, T=k_P/2 ✓ |
| 8 | Residual paragraph | `exp042_residual_probe.json` | verdict SMALL_CLASS_EMPTY; 4,621; bucket_T_lt_half_k=0; residual_parents=[]; upper_law_breaks=0; 37 shapes; unresolved 0; exception ℓm=72 ✓ — but scope words wrong (MEDIUM-1) |
| 9 | [5,12] sentence | `exp040_circuit_distance.json` | both circuits lb=5 certified (w3 MITM, w4 matching MITM + star exact-cover, all completed, no timeouts, no witnesses); ub=12 via EXP-029 certified witnesses; stagnation all true ✓ (wording LOW-3) |
| 10 | Environment | live `.venv` | Python 3.13.9, numpy 2.4.6 + 8 packages all exact ✓; **899 collected vs 897+1 quoted** (LOW-4); 748/748 pass on the 3 cheap files |
| 11 | §8 check counts | `pytest --collect-only` | test_pbb_survival 15 ✓, test_pbb_nogo 17 ✓ |
| 12 | §4 `12_6_0193` certificate | `results/certificates/pbb_12_6_0193_distance.json` + sector files | MITM ≤5 ✓; 72/72 weight-6 = stored pure-Z rows ✓; 4 UNSAT sectors at 501.98/673.75/421.97/575.90 s (cadical195) → "502/674/422/576" ✓; weight-12 witness two-path ✓; CERTIFIED_EXACT d=12 ✓ |
| 13 | 3.16M / 12 lattices | hunt scripts + JSONs | union of hunt lattice specs = exactly 12 distinct (ℓ,m); 3,228+1,029,827+2,131,004+368 = 3,164,427 ≈ 3.16M ✓ |
| 14 | d_X catalogue coverage | `xsector/catalogue_scan.json` | 368 rows; d_X(Q) exact 276; d_Z(P) cert 279; both 221; 17 equal; 0 violations — paper says "all 368" (MEDIUM-2) |
| 15 | PDF text | decompressed content streams | "extbfThe X sector." rendered ✓ (HIGH-2); tex contains exactly one tab |
| 16 | §6/§7 arithmetic | by hand | kd²/n cells, +2 parity gap (10−8, 6−4), 1.78→2.00, 16,493 s, "2.5× loose" (40/16) all ✓ |
| 17 | ρ_X statement vs note | text compare | paper §6 ρ_X=dim((R_CD+S_Z)/S_Z)≥dim bar Δ ≡ note J.1 (Δ=L·[C;D]⊆R_CD) ✓; J.0/J.4 citations ✓; G/H/I proofs use pure-Z survivors only — nothing earlier relies on d_X(Q)≥d_X(P) ✓ |
| 18 | Theorem-block re-read | hand proof-check | G(i),(ii)-statement,(iii), Lemma 1/2/3, Cor. 2, Thm H, Thm I proofs all check out; Lemma 1's ideal argument correct (L ideal; Δ=L·[C;D]; S_Z rowspace closed under monomial left-multiplication) ✓ |

**Not re-verified within budget (no finding implied):** full 897-test suite execution (SAT/stim scale; only LOW-4's off-by-one observed at collection), EXP-042 re-execution (would overwrite the artifact; 61.6 s wall per its own record), EXP-029 latency re-measurement (ratios traced exactly to `exp013_isolated_latency.json`: p50 10.32×, p95 16.60×, p99 6.31× — the p99<p50 ordering is legitimate, each ratio is same-percentile across the two circuits), §3.1 solver-internals (6.2×/77.8→12.5 s), and the demotion-hunt re-runs.

## Bottom line

Round 2 finds **no critical hole**: the rank identity behind Theorem G(ii), the new §6 paragraphs' substance, the [5,12] intervals, the environment versions, and the 30/60 sweep numerology are all artifact-true. But the G(ii) proof sentence "the row space of $H_Q$ extends that of $H_P$" is demonstrably false on the paper's own examples (fix: two-line projection argument, HIGH-1), and the shipped PDF contains visible corrupted text (HIGH-2). Fix both, correct the three MEDIUM scope/wording items, and re-sync md↔tex; nothing else blocks submission.
