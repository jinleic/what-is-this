# Pre-statement — Gate C orientation sweep (record attack, ≤54)

**Committed 2026-08-30T11:38Z, BEFORE any sweep computation.**
Agent: Mm3GateC. Folder: `cs/mm3/`. Campaign: this directory.

## 0. Assignment and provenance of this statement

Main directed gate C: sweep ternary orientations of public rank-23 decompositions
over {−1,0,1} against the total-addition objective. All five predecessor anchor
values are frozen in the repo (README session updates + campaigns 1–7); this sweep
must reproduce them with the sweep counter BEFORE trusting any new number.

Two papers were scanned on Main's direction before this statement was committed
(records in §8); neither adds a member to the decomposition set.

## 1. Decompositions entering the sweep (exact, with sources)

| ID | Source file (reproducible from) | Ring | Anchor total (frozen) |
|---|---|---|---|
| `paper55` | `cs/mm3/src/tensor_data.py` (printed §5.2 blocks; gate-A verified 729/729 over Z) | Z, ternary | 55 = 13/14/28 |
| `perminov58` | `cs/mm3/scratch/cr58_cn122_ZT_reduced.json` (commit 98ba522), decoded with his loader conventions + CPERM [0,3,6,1,4,7,2,5,8] | Z, ternary | 58 (CSE convention; his JSON) |
| `sun56` | `cs/mm3/scratch/sun56_verify.py` SIDES (also `scratch/sun56_blocks.json`), factors already verified 729/729 over Z | Z, ternary | 56 = 13/13/30 |
| `mws59` | `cs/mm3/scratch/mws59_layout.txt` lines 232–260 = MWS paper (arXiv:2601.05272) Table 3, three `#`-separated 23×9 ternary blocks; 729/729 verified in campaign `2026-08-30T011537Z_*` | Z, ternary | 59 = 15/15/29 (machine recount, campaign `2026-08-30T012035Z_*`) |
| `stapleton60` | arXiv:2508.03857 Appendix A scheme — **not yet on disk**; to be fetched BEFORE sweep start; ENTRY CONDITIONS: (i) fetched full text; (ii) my transcription expands to 729/729 Brent over Z; (iii) recount = 60. If any condition fails, the decomposition is recorded `FAILED-ENTRY` with the reason and is NOT swept | Z, ternary | 60 (frozen from campaign `2026-08-30T011537Z_*` manifest) |

Every member is a rank-23 decomposition of M⟨3,3,3⟩ verified over **Z** (729/729
Brent identities, exact integer arithmetic). No F2-only scheme enters the set: a
scheme valid over F2 is not automatically valid over Z, and none of the four
on-disk sources is F2-restricted. Smirnov's cr58_cn122 and Perminov's cr58_cn122
are the same object at the same provenance commit (gate A: identical tensor after
CPERM, 0 mismatches) — they are ONE member (`perminov58`), recorded under both
names for posterity.

**Set-extension rule (pre-registered here):** the set is EXTENDED only if a new
public Z-ternary rank-23 decomposition appears that satisfies the same entry
conditions, and the extension is APPENDED here with date, source, and reason
BEFORE its sweep runs. The set never grows because a target was hit — only because
a new public decomposition appeared.

## 2. The group action (exact)

The matmul tensor's isotropy: for any invertible X, Y, Z ∈ GL(3,R), the action

  (U, V, W) ↦ (X·U·Y⁻¹, Y·V·Z⁻¹, Z·W·X⁻¹)

preserves the tensor (hence Brent 729/729) but NOT the addition count. Pre-registered
sub-group of orientations, per decomposition:

- X, Y, Z range independently over the 3×3 **signed permutation matrices**
  (48 each; triple product 48³ = 110,592), composed with the sigma-orbit
  sigma: (U,V,W) → (V, Wᵀ, Uᵀ), {Id, sigma, sigma²} (order 3, cyclic trace
  symmetry — precedent: arXiv:2607.29291 over F2 uses GL(3,2)³ × cyclic trace;
  here over Z the finite analogue is signed permutations).
- **Total swept orientations per decomposition: 110,592 × 3 = 331,776.**
- Equivalence: two orientations related by summand relabeling (column permutation
  of the 23 products, consistently across U, V, W) have EQUAL counts; the counter
  is column-order invariant by construction (§3), so no matching step is needed
  beyond canonical column ordering — verified against the anchors (§4): each
  anchor reproduces from its frozen factor block in its printed column order.
- **(U,V)-swap elements are EXCLUDED with stated reason:** swapping U↔V (left↔right
  roles without transposition) computes B·A, not A·B — the earlier gate-C session
  verified that such variants FAIL the Brent identities as decompositions of A·B.
  The sigma action (V, Wᵀ, Uᵀ) is valid (verified 57 in session 2) and is included.

**What is NOT in the group (stated plainly, rule 7):** general ternary
{-1,0,1} invertible matrices (GL₃ of the 3×3 ternary unimodulars contains
6,960−48 = 6,912 non-monomial elements per side — counted exactly on 2026-08-30,
identically sized group gives 6960³ × 3 ≈ 1.01×10¹¹ orientations — NOT swept,
budgeted out), and the continuous GL(3,Z)/GL(3,Q) part (infinite; not swept).

## 3. Addition-count model (verbatim from the target convention)

Identical to the gate-A/gate-B counter (campaigns 1–2, REPRODUCED anchor status):

- Inputs A_0..A_8 / B_0..B_8 (and the 23 products for the output map) are FREE.
- Each gate computes x + y or x − y of previous quantities; cost 1 per gate.
  Additions include subtractions. Negations are avoided by rearranging terms —
  sign changes and copies are FREE.
- Total = left additions + right additions + output additions.

**Per-orientation counting procedure** (no solver, pure deterministic bound):
for each of L/R/O factor maps F: 9 → 23 (target columns of the oriented factors):
  d(F) := |{ sign-classes of the 23 target vectors }  reduced  { ±e_i input directions }|.
  d(F) is a rigorous LOWER BOUND: every gate creates ≤ 1 new output direction.
  (gate B's certified theorems: floor-DFS impossibility at d(F) for the three
  55-paper blocks + the paper's own d(F)+1-witness ⇒ C = d+1 there; for other
  factor maps the DFS is re-run per map, with the SAME completeness argument.)
  The counter reports d(F) per side; when a witness or certified count exists
  (from a gate-A-verified SLP or transposition of one), the certified count is
  used instead and LABELLED as such. No float anywhere; all arithmetic over
  Python ints / fmpz; no interval needed — everything exact.

**Column-order invariance guarantee (goal of anchor test):** d(F) depends only on
the SIGN-CLASS SET of the target columns. Column permutations of the 23 products
are absorbed silently. Verified by the anchor test (§4): each frozen decomposition
must reproduce its frozen count from its printed column order AND from a
random 23-column relabeling (test: exact match of d on ≥1 random permutation).

## 4. Anchor reproduction (MUST pass before any new orientation is swept)

| Anchor | expected | source |
|---|---|---|
| paper55 | 13 / 14 / 28 (d-counter: 12/13/13 + witnesses) | gate A count; gate B floors |
| sun56 | 13 / 13 / 30 (d: 12/12/…) | Sun's own verify.py; session-2 recheck |
| perminov58 | 58 (CSE accounting: nonzero-coeff − 2·23 − 81 from HIS JSON) | gate A `check_perminov` |
| mws59 | 15 / 15 / 29 | machine recount, campaign `2026-08-30T012035Z_*` |
| stapleton60 | 60 total recount from printed Appendix A | campaign `2026-08-30T011537Z_*` |

Any anchor miss ⇒ STOP, fix the counter, re-anchor. The sweep counter is the SAME
code path as gate B's d(F) counter (module `gate_b_floor.py`, `factor_targets`,
`canon`, `prep`) extended to accept arbitrary factor blocks; the anchor test
proves this extension reproduces all frozen values before new data flows.

## 5. Validity check (mandatory, per orientation family)

For any orientation family that achieves d(F) total ≤ 54 on ANY triple:
- re-verify Brent 729/729 over Z (exact fmpz) on the oriented factors, and
- independently verify that the certified/witness circuit (if claimed) computes
  exactly those oriented factors, and
- the valid-orientation checker (U∘Y⁻¹ vs V∘Z⁻¹ etc.) re-derives the factors.

A count WITHOUT a passing Brent check on the oriented factors is DISCARDED —
count never outranks validity.

## 6. Budget

- 5 decompositions × 331,776 orientations = **1,658,880 (U,V,W) triples** swept.
- Per-orientation cost is dominated by 3 sign-class-count reductions over 23
  ternary vectors in Z⁹ — O(23) each; target wall-clock ≤ 2 h single-thread
  (M3 Ultra); checkpoint every ~10⁴ triples to `sweep_checkpoint.jsonl` (crash-safe).
- If budget expires, the sweep is reported as PARTIAL with the exact swept prefix
  enumerated; no universal claim ever issues from a partial sweep (rule 7).

## 7. Escalation policy (fixed before launch)

Any total ≤ 54 ⇒ freeze, re-verify Brent independently, then `hub` message Main
with the exact numbers; **write no claim anywhere** (README, RESULTS, index)
before Main's independent re-verification.

## 8. Notes on Main-directed paper scans (2026-08-30, recorded here per the
set-extension rule because they happened BEFORE commit; neither added a member)

- **arXiv:2608.27434v1** (Kassabov–Landsberg–Souza–Speegle, "Disjoint and nearly
  disjoint sums of matrix multiplication tensors and their centroids", 2026-08-27):
  full HTML scanned (section titles 1–10 enumerated; text searched for
  "3×3", "rank 23", "Laderman", "Smirnov", "M⟨3,3,3⟩", "addition", "additive":
  zero hits on all seven). The paper's explicit decompositions are border-rank
  decompositions of Strassen's T_Str,q, Schönhage's T_Sch,u,v and
  big-centroid tensors (T_cen_2,3,4, T_cen_3,3,3), at exponent window
  2.46016–2.548 — no rank-23 additive decomposition of the 3×3 tensor is
  produced. **Negative recorded; set extended by NO member.**
- **arXiv:2607.29291** (Palladinos, "SAT Certificates for the Matrix-Multiplication
  Challenges over F2…", 2026-07-31): abstract read first-hand. Sanity-check on
  the gate-C design (GL-isotropy + cyclic trace + slot-matching is published
  live technique — over F2, GL(3,2)³; here over Z, the signed-permutation
  analogue). Also a warning: encoding-completeness defects can flip benchmark
  expectations wholesale (their ten "expected-UNSAT" were all SAT; cause =
  positive unit clauses that REQUIRED selected incidences without FORBIDDING
  extra ones). One bounded audit of gate B's encoding is scheduled in this
  campaign's notes (§9) — gate B's three-way agreement checks the same modelling
  premise three times, not the premise itself.

## 9. In-campaign obligations (beyond the sweep itself)

- (a) written gate-B encoding audit: which clauses FORBID vs REQUIRE (the
  encoding must not merely require selected additions without forbidding extras);
  (b) one-line two-direction soundness argument for the floor CNFs
  (assignment ⇔ circuit-of-that-count); (c) any under-constraint → ESCALATE,
  no quiet patch; (d) if sound, one-line statement that gate B is strengthened.
- Full landscape per decomposition (min + distribution), not just the minimum.
- Frozen campaign at close: scripts, raw sweep output, checksums, manifest.
