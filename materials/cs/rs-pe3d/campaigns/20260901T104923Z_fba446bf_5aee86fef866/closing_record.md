# Run 2 closing record — gate H-MINLINE-DGE3 — verdict FROZEN-INCONCLUSIVE

Campaign: 20260901T104923Z_fba446bf_5aee86fef866. Prereg committed pre-compute
at git `27001a3` (sha256 965000475aed0a9aea9fbbd54325f16611db93e4eef06196dee64c6feaebd
...c1, hash-bound in manifest.json). Per Main's verdict-semantics instruction
(a positive universal theorem is never buried under the negative gate's single
status) and the prereg's fixed verdict rule ("promote only if P1-P5 are all
discharged with no gap AND all controls pass exactly; any gap: theorem stays
OPEN, campaign closes FROZEN-INCONCLUSIVE"), Run 2 closes FROZEN-INCONCLUSIVE.

## What was completed here

- Preregistration committed (theorem T-DGE statement with the parallel-class
  hypothesis, proof obligations P1-P5, controls M1-M6, verdict rule).
- Matrix-model instrument `matrix_model_controls.py`: **M1 structured
  fiber-regime configs 6/6 EXACT PASS** over GF(13) — [(3,5),(3,5)] -> 100,
  [(3,5),(4,5)] -> 50, [(3,4)]^3 -> 192, [(3,4),(4,4),(4,4)] -> 64,
  [(4,6)]^2 -> 180, ([(4,6),(5,6)]) -> 90, every one with off = 0 and
  spark(product) = min spark verified by exhaustive below-d search plus a
  witnessed dependent d-set. These corroborate T-DGE regime 1 and Main's
  independent numbers (dep=100/50/192/64/180/90, off=0).

## What is incomplete (exact gaps; theorem NOT promoted)

1. **Proof artifact P1-P5 not written/discharged.** The 2-factor independence
   lemma (P1), the size-d r=d_A circuit step (P2, including the zero-
   coefficient handling of the isolation functional flag by Main's fable
   review), the no-parallel-columns product lemma (P3), the regime-2
   classification + count (P4), and the RS/GRS specialization incl.
   Lambda-invariance of proportionality classes (P5) are REGISTERED as
   obligations in the prereg but no complete proof document exists yet. Per
   Main's instruction ("encode in theorem artifact only if complete") no
   theorem artifact was written. Theorem T-DGE status: **OPEN (proof
   outstanding, evidence-only)**.
2. **M2 diagonal-regime structured configs ambiguous**, run halted on
   [(2,5),(2,5)] (instrument got 300 = all C(25,2) pairs with distinct-column
   spark-2 factors; Main's reported 120 for "[2,2]" decomposes as 48 fibers +
   72 diagonals, which matches spark-2 factors with 4 columns each —
   N = 16, product rank 1, all C(16,2) = 120 pairs dependent, fibers
   2 x C(4,2)*4 = 48). Config semantics (column counts and parallel-class
   structure of Main's instrument) must be pinned before the assertion is
   meaningful; instrument assertion needs the regime-2 count including
   MULTI-class factors A_i = Sum_k m_{i,k}^2 formula:
   unordered distinct dependent pairs = (Prod_i A_i - N)/2, fibers vs off
   split per the prereg's D-subset sum.
3. **M3 randomized (21 x GF(7)), M4 RS reconciliations (13 rows incl.
   t != (1,1,1) circuit-form rows), M5/M6 plants and Run-1 cross-derivation:
   not run.**

## Next gate (ready to run, price minutes-hours)

Re-open a fresh campaign (same prereg content + pinned M2 configs with the
A_i = Sum m^2 pair-count formula and Main's unambiguous [2,2](4,4) row),
write theorem_t_dge.md discharging P1-P5 (P1/P2 per Main's steering route
with the audited corrections; P4 count formula validated against the Run-1
census 24 = 16 + 8 and Main's 15 census instances), run M1-M6, then freeze +
close. The Run-1 falsification (FROZEN-NEGATIVE) is INDEPENDENT of this run
and unaffected.
