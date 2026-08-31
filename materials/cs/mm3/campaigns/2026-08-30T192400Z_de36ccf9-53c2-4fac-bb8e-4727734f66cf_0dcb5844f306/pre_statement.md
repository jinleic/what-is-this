# Pre-statement — off-diagonal sandwich census + certified decision (agent `Mm3OffDiag`)

**Committed 2026-08-30T19:24Z, BEFORE any computation of this campaign.**
Folder: `cs/mm3/`. Campaign: this directory.

## 0. Prior state this campaign inherits, and the two corrections it must honour

Frozen prior artifacts: gate-C sweep
`2026-08-30T113838Z_c9532c07-61e4-40e1-8b8f-b599d1fad6dc_7e0e84591c56`
(monochrome/monomial sandwiches; 55/58/56/59/60 anchors), record attack
`2026-08-30T174243Z_edbdd408-840f-42b4-9be0-192a95783ab9` (diagonal non-monomial
census). Two material corrections occurred BETWEEN those freezes and this
commit, both confirmed jointly by Main and this agent:

1. **The frozen `gatec_sweep.inv_sp()` is the TRANSPOSE, exact only on the 48
   monomial matrices.** The record attack fed all 6960 into it (precondition
   violation at the call site, extensively confirmed by Main).
2. **Derivation (this agent, pass 0, pre-compute):** the frozen map inside
   `sandwich()` is `u' = X⁻ᵀ·U·Yᵀ`, `v' = Y⁻ᵀ·V·Zᵀ`, `w' = X⁻ᵀ·W·Zᵀ` — not
   `X·U·Y⁻¹`. The honest tensor-preserving conjugation family is
   `tr((P A Q⁻¹)(Q B R⁻¹)(R C P⁻¹)) = tr(ABC)` for any invertible `P,Q,R`, whose
   factor-block map is
   `u' = Pᵀ·U·Q⁻ᵀ`, `v' = Qᵀ·V·R⁻ᵀ`, `w' = Rᵀ·W·P⁻ᵀ` (cyclic structure:
   W-side is `Rᵀ·W·P⁻ᵀ`, SYMMETRIC in the roles of P,Q — unlike the frozen
   formula, whose W-side `(X⁻ᵀ, Zᵀ)` transposes the roles and is provably NOT a
   tensor conjugation for `X ≠ Z`).
   On the DIAGONAL `(X,Y,Z) = (G,G,G)` the frozen map equals the honest map at
   `(P,Q,R) = (G⁻¹,G⁻¹,G⁻¹)`; `G ↦ G⁻¹` is a bijection of the 6960 preserving
   the monomial subset, so the frozen diagonal census 6912/0/48 and the
   certified ≥55 over the 288 monomial orientations **survive as counts and as
   valid claims about the honest action**. This survives-by-relabeling claim is
   itself re-verified below (anchor A2), not inherited.

Per repo rule 5 both corrections are recorded in the README session entry this
campaign appends; the frozen campaign directories themselves are not edited.

**Owner-hypothesis provenance (transparency, rule 17):** the owner's per-factor
collapse hypothesis was FALSIFIED BY THE OWNER before this commit (his own
300k-sample proxy: ~23.55% of non-monomial×non-monomial products stay ternary,
exact witness pair whose product is monomial; owner-supplied figures labelled
REPORTED, never an input). This pre-statement is designed WITHOUT the collapse
assumption; the falsification is recorded rather than repaired. The owner's
`4608` = # of non-monomial G with `G⁻¹` also ternary is likewise a
necessary-condition SCREEN only — a candidate upper bound, never assumed.

## 1. Falsifiable outcomes, fixed before compute (rule 14)

Exactly one of:

- **(COUNT)** An exact integer count of the ternary-alphabet-valid orientations
  over the off-diagonal `(P,Q,R)` sandwich product, reported per (D, σ) class,
  with explicit factorization verdict (subcube test). If pairwise factorization
  fails, the count is reported as an UPPER BOUND, labelled as such, never as a
  count.
- **(REC)** A ≤54-addition scheme found among decided survivors — all 729
  Brent identities EXACT over ℤ in the frozen convention, two independent
  arithmetic paths, ESCALATE TO MAIN BEFORE WRITING ANYWHERE.
- **(NOGO)** Certified `total_lb ≥ 55` (d-counter + floor-DFS + transposition
  model, the frozen convention) over an explicitly named, fully-enumerated
  decided subset of survivors, with every non-achievable floor instance
  carrying dual-checker (drat-trim + lrat-check) DRAT/LRAT certificates and
  verbatim frozen checker output. Any UNDECIDED remainder at campaign wrap named
  EXACTLY (per (D,σ): predefined enumeration order, prefix decided, remainder
  count) — a partial decision certifies only its prefix.

A failed search over an unnamed set certifies anything and is prohibited.

## 2. Search space, fixed NOW (rule 16, no movement after results)

- Decompositions: `paper55`, `sun56` (LOADERS in the frozen
  `gatec_decomps.py`), identical to the record attack's named set. No other
  decomposition enters this campaign, even if budgets allow.
- Sandwich: `(P,Q,R)` each independently over the 6960-element ternary
  unimodular set T; the FULL off-diagonal product, ordered 6960³ = 337,153,536
  triples per (D, σ-power), ×3 σ-powers ×2 decompositions = 2,022,921,216
  scheme-instances total per decision point. Both per-(D,σ) and total figures
  will be reported (the two readings that conflated the last sessions' totals).
- Action (honest, derived above): `u' = Pᵀ·U·Q⁻ᵀ`, `v' = Qᵀ·V·R⁻ᵀ`,
  `w' = Rᵀ·W·P⁻ᵀ`, rows as 3×3 blocks reshaped to the 23×9 convention; then the
  frozen σ-orbit {Id, σ, σ²}, σ = (U,V,W) → (V, Wᵀ, Uᵀ), MAXIMAL cycle per the
  frozen note (non-maximal pair-swaps compute B·A; excluded as before).
- Validity predicate (fixed): a scheme-instance is VALID iff (a) all 23 rows
  of all three blocks are ternary {−1,0,1}, and (b) all 729 Brent identities
  hold EXACTLY over ℤ. Both checked by census machinery; no other admission
  criterion. We do NOT require P⁻¹, Q⁻¹, R⁻¹ ternary (the sandwich's products
  are what must be ternary; that P/Q/R inverses happen to be ternary is
  neither assumed nor required — observation recorded in artifacts only).

## 3. Census design (factorization to be PROVEN then machine-tested)

**Algebra to be established before counting (pass A, before census table
launch):**

- T decomposes as T = A·𝔐 where A = 1160 canonical ternary-unimodular
  DATA-right matrices with non-decreasing rows (representatives of
  T/𝔐-mod-right, 𝔐 = the 48 monomials), and 𝔐 the monomial matrices, since
  6·1160 = 6960. Verify |orbits| exactly 6, 1160 total, or FALL BACK to full
  6960³ (abandon factorization, enumerate all 48.4B data-triples NO — fall
  back = report infeasibility honestly; the count design dies with the
  factorization).
- Under the honest map, u' = Pᵀ·U·Q⁻ᵀ has row-spaces/permutation content
  governed by P's LEFT part only: writing P = a·p, Q = b·q (a canonical,
  p monomial), u' = pᵀ·(aᵀ·U·b⁻ᵀ)·q. Ternarity, d, floor-DFS results, and
  Brent validity are invariant under right-monimal factors: **the alphabet
  predicate and the bound layer must depend only on the DATA-TRIPLE
  (a₁,a₂,a₃)**. This is the factorization claim, stated as a theorem-shaped
  claim with a completeness argument (permutation content pulled out
  exactly), and machine-checked by the subcube test below.
- Census outputs (per decomposition):
  - `alpha_U[a1][a2]` ∈ {0,1}: is `a1ᵀ·U·a2⁻ᵀ` ternary rowwise,
  - `alpha_V[a2][a3]`, `alpha_W[a3][a1]` similarly on the *data-pair* sides,
  - Hamming-weight / d-count summaries cached per pair,
  - exact survivor count = #{(a1,a2,a3) ∈ A³ : alpha_U[a1,a2] ∧ alpha_V[a2,a3]
    ∧ alpha_W[a3,a1]} × 648 (the 2·2 monomial-quotient factor 2² × 3 σ-powers
    per decomposition... to be reconciled exactly: 48 monomial triples per data
    triple ×3 σ = 432; 216 per σ... — the exact per-permutation instance
    arithmetic is fixed in code comments before table launch and cross-checked
    against the frozen 48³×3 = 331,776 monomial count and against an
    independent recount).
- Instrument control (rule 14, mandated by Main): DIRECT enumeration of the
  left-anchored subcube {(1,a2,a3)}: 1·1160·1160 ≈ 1.35M sandwich evaluations,
  role-set exact trifecta test — run DIRECT-3 and DIRECT-3', plus sample-level
  comparison of direct vs table-derived alternation tests. If ANY discrepancy:
  the factorization claim is FALSIFIED, the count is reported as an upper
  bound, and options (b)/(c) are NOT pursued this campaign (rule 16: no
  redesign after seeing results; a redesign is a separate pre-statement).
- RING GUARD: all Brent and ternarity checks in EXACT ℤ arithmetic
  (Python ints for census; fmpz for the anchor reproduction and for any
  witness re-verification). No F₂ encoding anywhere; F₂ validity does not
  imply ℤ validity.

## 4. Anchors mandatory BEFORE any new number (order, not suggestion)

- **A1:** Reproduce the frozen diagonal census claim 6960 → (6912, 0, 48) per
  decomposition under BOTH maps (frozen transpose-as-inverse AND honest), in
  the SAME script, same G enumeration order (combinatorial
  `itertools.product` order as in `census_fmpz_path.py`), and report the
  verdict-vector diff as a hard number. Expectation from the relabeling
  argument: both give 6912/0/48 with the 48-element survivor sets differing by
  inversion (same set, as monomials are inversion-stable — verify exactly,
  don't assume: check set equality AND sample the non-monomial per-G verdicts
  elementwise).
- **A2:** On ≥10 pseudorandom non-monomial (P,Q,R) triples with seed fixed
  (seed=20260830): frozen-style map (X⁻ᵀ, Yᵀ formula) vs honest map —
  Brent-failure counts for both. PREDICTION (pre-registered): frozen map FAILS
  Brent for all sampled triples with X≠Z (its W-side is not a conjugation);
  honest map PASSES on all. If the honest map fails Brent on ANY triple, the
  derivation is wrong: STOP, re-derive, do not run the census.
- **A3:** 55/58/56/59/60 anchors from the frozen harness
  (`verify_anchors.py` convention) re-run under my own code path (fresh Brent,
  fresh loaders untouched) — 729/729 each.
- Only AFTER A1–A3 pass does the census table launch.

## 5. Decision layer over survivors (pre-registered caps, all-or-nothing)

- Per VALID data-triple (a1,a2,a3) in a random-but-FIXED order (rng seed
  20260830), for each sigma-power: compute d(L), d(R), d(O_fac) via the frozen
  `prep` + floor-DFS convention (all exact ℤ); total_lb =
  Σ(sides) + 14 (transposition gap; transposition-activity precondition
  checked per instance; if INACTIVE the gap-JUSTIFY step fails and the
  instance is reported individually rather than summed — see §7 adjudication).
- Enumeration ORDER (rule 16, fixed now): lexicographic in (a1,a2,a3)
  canonical indices; survivors processed in lex order of canonical index
  triples; per (D,σ).
- Decide: dfs-floor-decide each of L/R/O_fac per data-triple, caching per-side
  (d, floor_ok, DFS states); total_lb computed identically to the frozen
  `certified_landscape.json` row model; escalation on total_lb ≤ 54 iff the
  rows are check-satisfied: ALL d/floor-checked + witness-verified. A
  total_lb ≤ 54 row does NOT create a ≤54 scheme claim by itself: it flags an
  orientation for LAYER-S3 exact scheduling (gates considered, aux-1, per the
  frozen `build_floor_cnf` + aux-extension conventions, kissat 4.0.4,
  dual-checker per the frozen `cert_gen.py` pipeline).
- **Wall-clock caps (fixed before compute; hard):**
  - Anchor + derivation + subcube controls (pass A): 2 h.
  - Census table construction (pass B): 8 h wall; expected ≪1 h at the
    measured rate.
  - Survivor decision DFS (pass C): 24 h wall. If unfinished at cap,
    remaining survivors named exactly (count + interval of lex-indices),
    reported as UNDECIDED; never claimed as excluded.
  - LAYER-S3 (only if a ≤54 flag arises): 24 h solver wall total.
  - Total campaign hard cap: 48 h. Partial = prefix-named remainder, per §1.
- **Primary deliverable priority order (fixed now):** (1) exact survivor
  count per (D,σ) with factorization verdict; (2) diagonal dual-map diff
  (already part of A1); (3) survivor decision prefix + certified totals; (4)
  remainder naming if cap hit. A partial pass-C with exact remainder naming
  is a REPORTABLE DELIVERABLE (Main-directed); it is not a failure.

## 6. Certificates and toolchain (frozen, dual-checker mandate)

- Solver (only if LAYER-S3 triggers): kissat 4.0.4 (/opt/homebrew/bin/kissat,
  patched env identical to the frozen pipeline); CaDiCaL 3.0.1 cross-check on
  the audited instance; binary DRAT proofs converted with `drat-trim -L`;
  verification by BOTH checkers; binaries pinned to the frozen
  snapshot-tree 2e3b2dc SHA-256 drat-trim
  `111b0405566d55629f5d391b80b3220cc7a3ebc3cd406a35895ff65cc9acd5e4`,
  lrat-check `b4bdebfcc40da664be2fd416451b43f9e764e5e92383ce2b4e9d1f05148e9b07`
  (locally rebuilt from `tools_snapshot/2026-08-30T031544Z_*/tools_snapshot/`
  and byte-verified against the pinned SHAs before use on 2026-08-30T19:20Z).
- Encoding: ONLY the frozen `build_floor_cnf` (+ aux-1 extension pattern of
  campaign 031544Z controls); encoding audit inherited
  (`gate_b_encoding_audit.md`, no positive-unit under-constraint, both-direction
  model ⇔ circuit). Any NEW encoding (should a genuinely new floor question
  arise) would need its own audit BEFORE use — and per rule 16 I do not intend
  one; the floor question set is fixed.
- python 3.14.3, `.venv`, python-flint 0.9.0 (fmpz for anchor Brent and all
  witness re-verification); numpy 2.5.2 census tables (uint8 + uint16 counts;
  popcount via `np.bitwise_count`; collision-freedom verified by an exact
  Python-int audit script). Boundary: census counting in uint64 accumulators
  with modular verification — numpy sums EXCEED 2⁶⁴? None: max triangle count
  ≤ 1160³·6⁴ ≈ 3.2·10¹² < 2⁶⁴ ≈ 1.8·10¹⁹; overflow impossible by argument
  recorded here; the exact Python re-audit recomputes the same count by
  different accumulation order as an independent path (this census has two
  count paths, in the spirit of the record attack's two-path discipline, but
  BOTH now implement the same honest map — map-independence comes from the
  adversarial test battery instead).

## 7. Adjudication (fixed)

- If survivor count = 48³·2-basic×...== exactly 331,776 per (D,σ)-class
  (= 48³ per D-σ... exact reconciliation to be locked in code comments BEFORE
  table launch, target: 48³ monomial-only per (D,σ) = 110,592 per D,
  331,776 with σ): hold the diagonal verdict, update README with the exact
  boundary statement.
- If survivor count > 331,776·(D,σ equivalents): report exact count, declare
  the owner's hypothesis (already falsified) WITH NUMBERS from the actual
  action, proceed to decision layer on the fixed caps.
- If factorization fails subcube test: count = upper bound (labelled), report,
  do NOT proceed to decision layer this campaign.
- If total_lb ≤ 54 flags arise: LAYER-S3; found scheme ⇒ ESCALATE (rule as in
  §1 REC); all misses reported.
- Transposition-inactive orientation: excluded from +14 totals (labelled
  individually in artifacts; not silently folded into totals). Precondition:
  activity audit from `gatec_audits.py` method.

## 8. Non-goals (unchanged from the assignment)

No re-run of gates A/B/C as such; no monomial-transfer re-derivation; no other
target; no edits to `cs/RESULTS.md`, `cs/PROGRESS.md`, or any frozen campaign
directory; README changes are APPEND-only below the contract; escalation to
Main before any write bearing on a published result (the 55 record or the
published optimality claims).

## 9. Escalation policy (fixed)

A ≤54 scheme or ANY disagreement with a published claim (e.g. a defect in
arXiv:2607.28676's own published bounds under the honest maps) → report to
Main via hub BEFORE writing anything. Certified no-gos and census counts over
the named set enter the README session entry after the run; those do not
require pre-clearance.
