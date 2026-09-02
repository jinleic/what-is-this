# VERDICT.md — Route A provenance campaign: Strassen 1983 3-slice commutator
# theorem — primary import + clause checks

Campaign: `campaigns/2026-09-01T03:21:00Z_routeA_prov/`
Pre-statement committed BEFORE compute: commit `dd40d2b`
("oct-rank: preregister route A provenance campaign (OctRankRouteA) …",
master), directory named with the actual UTC start `2026-09-01T03:21:00Z`.
All scripts frozen here with SHA256 in `checksums.sha256`; the verification
output of the exact check is `clause_check.out` (exit 0).

## Headline

**ROUTE A: CLOSED-NEGATIVE at the provenance layer, and now provenance-clean.
The first-hand Strassen 1983 text yields, for every linearly independent
(u, v, w) in R^8 exactly the bound**

    rank_R(L_u, L_v, L_w) ⩾ 12      (12 = 8 + 8/2)

**— strictly weaker than the standing certified floor 13. The "naked"
n + rank(commutator) reading (16) does not exist in the primary text.**
No 3-slice theorem in Strassen 1983 exceeds n + ½·rank(commutator); the
Blaeser "twin" has no primary support (Blaeser is not even cited in the
paper; Lickteig is only mentioned as an independent announcement). The
≥ 16 lower-bound route is dead as a theorem-family route. S3 stays OPEN at
13 with gap exactly 1 to 14, and the published window
18 ≤ R_R(T_O) ≤ 25 is untouched (nothing here bears on it; no refutation
of published work is claimed).

## What changed relative to Route AF (delta, no re-litigation)

Route AF closed Route A with Landsberg Thm 6.1.1 + Koiran (secondary
readings) and flagged "Strassen primary provenance UNESTABLISHED". THIS
campaign closes exactly that gap: the primary (Strassen 1983 LAA 52/53
645-685) was obtained, read, transcribed, hashed, and the theorem imported
from it. The conclusion (12, not 16) is CONFIRMED — now from the primary.
No number of any prior campaign changed.

## Provenance record (rule 1/17a discipline)

Primary text: `sources/Strassen_1983_Rank_and_optimal_computation_of_generic_tensors_LAA52-53.pdf`
SHA256 `f8481709de3a93e5747ce45e5c19aae81225139f39db347ef373a82df5bd0376`
(1,936,012 bytes). Acquisition chain (full detail: `sources/ACQUISITION_LOG.md`):
ScienceDirect 403 (challenge page even via Wayback 2024 capture); EuDML
zero-hit (content-layer negative archived); OpenAlex/SemanticScholar
bronze-OA pointer 403s; **scholar.archive.org → CORE 82110360 → Wayback
2019-04-16 capture (mimetype application/pdf) → PDF obtained**.
Identity verified from the title page: "Rank and Optimal Computation of
Generic Tensors / V. Strassen / … Universität Zürich / Dedicated to
Alexander M. Ostrowski … / LINEAR ALGEBRA AND ITS APPLICATIONS 52/53:645-685
(1983)". Text layer (OCR-degraded): `sources/Strassen_1983_text.txt`
SHA256 `b138a8b896c1b9a8056510adda3777aba70c1808c7ce7add2bba795b3e106e49`.
VERBATIM theorem transcriptions + glyph-level verification method:
`theorem_import.md`. OCR landmine: the text layer renders inequalities and
fractions as ASCII garbage ("R(t)~n++rank(...)"); every load-bearing
formula was verified against 400 dpi page rasters, glyph by glyph
(`verification_images/`), including a pixel-ink test distinguishing the
theorem's BORDER-RANK symbol R̲ (underbar present: ink rows 1011-1014 below
the R baseline) from the proof's plain-R display (4.2) (0 below-baseline
ink pixels) on the same page.

## The exact imported statements (see theorem_import.md for full quotes)

- **Theorem 4.1** (p. 673, Part 4 "A DETERMINANT FOR 3-SLICE TENSORS",
  dim U = dim V = n, dim W = 3, slices A, B, C):
  "Let A be invertible. Then **R̲(t) ⩾ n + ½ rank(BA⁻¹C − CA⁻¹B)**"
  (R̲ = border rank; algebraically closed base field k).
- **The proof's first display (4.2) is the plain-rank version**:
  "We first show **R(t) ⩾ n + ½ rank(BA⁻¹C − CA⁻¹B)**", proven by pure
  coefficient comparison on an existing rank-r decomposition
  (rank(BC − CB) ≤ 2(r − n), his (4.4), plus induction over vanishing
  diagonal entries) — NO algebraic closure, NO Zariski closure in this
  step, so it is valid over ℝ; the printed proof's final step passes to
  the border-rank form by closure/irreducibility only because the paper's
  Theorem 4.6 needs a border-rank hypersurface statement.
- **Corollary 4.2** (p. 676): for a finite-dimensional ASSOCIATIVE algebra
  Λ with unit, a module N, dim N = n, b, c ∈ Λ:
  **R̲(N) ⩾ n + ½ rank_N(bc − cb)** (underlined R; subscript N on rank —
  verified twice at glyph level).
- Negative first-hand results: no factor-free n + rank(commutator)
  bound exists anywhere in the paper; "commutator" never appears as a
  name for BC − CB; Blaeser uncited; Lickteig only mentioned.

## Clause-by-clause hypothesis check (rule: a hypothesis you cannot verify FAILS)

For our 3-slice tensor t = (L_u, L_v, L_w) (dim U = dim V = 8, dim W = 3):
[verified in `clause_check.md`, arithmetic in `clause_check.out`, all exact fmpq]

| clause | hypothesis | status for D2 (independent triples) | type |
|---|---|---|---|
| 1 | 3 slices, dim U = dim V = n, dim W = 3 | SATISFIED (8, 8, 3) | HUMAN-AUDITED (definitional + exact shape check) |
| 2 | A = slice 1 invertible | SATISFIED: N(u) = ‖u‖² > 0 over R; L_{ū}L_u = N(u)I verified ENTRYWISE EXACTLY on all 8 basis slices and a generic rational triple | MACHINE-VERIFIED (exact) |
| 3 | conclusion uses M = BA⁻¹C − CA⁻¹B over the base field | COMPUTED: M = (L_vL_u⁻¹)(L_wL_u⁻¹) − (L_wL_u⁻¹)(L_vL_u⁻¹); M² = (−4·det Gram/N(u)³)·I verified ENTRYWISE EXACTLY on 5 basis triples + 1 generic triple; rank(M) = 8 in every case | MACHINE-VERIFIED (exact) |
| 4 | (Corollary 4.2 only): Λ ASSOCIATIVE | NOT SATISFIABLE with Λ = O (octonions are non-associative): the corollary does NOT apply to the octonion family; the reach to our tensor is exclusively via Theorem 4.1's 3-slice tensor statement, which needs no algebra on the slices | HUMAN-AUDITED (definitional) |
| 5 | real-rank transport: Strassen's k algebraically closed | The (4.2)-step is field-independent (printed proof uses no closure); complexification alternative R_R ≥ R̲_k also consistent; the campaign's bound rests on the ℝ-valid (4.2) argument, NOT on the closure step | HUMAN-AUDITED (line-by-line reading of the printed proof) |

## The derived bound, the control, and the adjudication

- Universal bound (all independent u, v, w ∈ R^8):
  rank ⩾ 8 + ½·8 = **12** [HUMAN-AUDITED deduction; inputs: CITED-DEPENDENCY
  theorem (primary read) + MACHINE-VERIFIED rank(M) = 8 arithmetic].
- **tau CONTROL (n = 4) — passed**: M_tau = L_iL_j − L_jL_i = 2L_k exactly,
  rank 4, bound 4 + 4/2 = 6 ≤ 7 = true rank. Consistent, non-tight. The
  naked reading 4 + 4 = 8 WOULD have refuted the control — no such form
  exists in the primary text (confirming Route AF's retraction).
- **tau ⊠ s CONTROL (8 × 8 × 3)**: M = blockdiag(2L_k, 2L_k) exactly,
  rank 8, bound 12 ≤ 13 ≤ true rank ∈ [13, 14]. Consistent.
- Wy the theorem cannot give 13 or 14 here: rank(M) = 8 is proven
  (Route AF's exact identity D² = −4N(u)·detGram·I with detGram > 0 —
  MACHINE-VERIFIED anchors re-run in this campaign), so the commutator
  term is capped at 8 and the bound at 12 on the whole declared class.
  No parameter of the theorem changes this.

## Evidence labels (strict)

- Theorem text/transcription: CITED-DEPENDENCY with PRIMARY READ
  (provenance now ESTABLISHED first-hand; upgrade of Route AF's
  UNVERIFIED flag).
- Clause 2/3 arithmetic (inverses, M², ranks on 6 triples; tau and tau⊠s
  controls): MACHINE-VERIFIED (exact fmpq, exit 0, deterministic).
- The deduction "rank(M) = 8 on ALL independent triples ⇒ universal
  floor 12": HUMAN-AUDITED (deduction rests on Route AF's frozen
  commutator_derivation.txt (4) whose coefficient anchors are
  MACHINE-VERIFIED, plus the CITED theorem).
- The claim "(4.2) is valid over ℝ": HUMAN-AUDITED reading of the printed
  proof; NOT machine-verified; independently the same numeric bound is
  what our own exact rank argument would give via contraposition, so no
  load-bearing step rests on the human audit alone.
- No COMPUTATIONAL-EVIDENCE is load-bearing anywhere in this verdict.

## Retractions / corrections (rule 5)

- None of this campaign's numbers retracts any prior verdict. The prior
  retraction (Route AF: "naked twin 16") is CONFIRMED at the primary
  layer: the factor-free reading has no source.
- This campaign's Route-A pre-statement draft contained two draft-thought
  sentences ("Wait— stricken sentence…", "附加") removed before freeze;
  the frozen theorem_import.md carries the cleaned bookkeeping. No
  artifact output was ever based on those sentences.

## Resource note

Compute: 1 exact script, 0.09 s wall (in-process accounting per Main's
ruling: `time.process_time()` equivalents in `time` output; no ps-derived
numbers used or needed). No heavy slots touched; no breaches.

## Rule-7 scope sentence

Covered: the Route-A provenance question over the declared domains D1-D4 —
first-hand acquisition/quantification of Strassen 1983's 3-slice
statements (Theorem 4.1 + the (4.2) proof-step + Corollary 4.2, with
negative results on any stronger form), clause-by-clause hypothesis
verification for (L_u, L_v, L_w) with u,v,w independent in R^8 (the full
S3-declared class; 8 basis-slice invertibility checks + 1 generic-triple,
6 exact M²/rank checks, tau control, tau⊠s control), and the adjudicated
consequence 12 < 13. NOT covered, hence not excluded by anything here: a
rank ≥ 14 argument for 3-slice L-families from any OTHER theorem family
(none found in this primary; other primaries unswept by this campaign),
rank-13 witness/infeasibility for the specific (1, i, j) triple (that is
the Route-F question), any statement about R_R(T_O) itself beyond the
untouched window 18 ≤ R_R(T_O) ≤ 25, border-rank work over ℂ beyond the
imported statements, dependent triples beyond the frozen square-zero
observation, and any n outside {4, 8}.

## Named next action

Route F job (second assignment item): decide rank((L_1, L_i, L_j)) ∈
{13, 14} — exact 13-witness or impossibility, per the af-triple13-existence
plan already named in the Route AF verdict. Route A offers no input to it
(12 < 13): the campaign proceeds to Route F only if time remains in this
session; otherwise the named next action stands for the next run.
