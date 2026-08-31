# Convention note — unique-family identification and semantics reconciliation

**Timestamp 2026-08-30T21:0xZ. Machine-verified on the shared arrays; resolves the
Main–agent ternarity dispute in `provenance_resolution.md` (which located the W-swap
but mis-named the mechanism).**

## Facts (all replayed in this campaign's artifacts)

1. **Arrays:** both parties load IDENTICAL sun56 rows (verbatim row comparison: all
   four cited rows equal). No storage divergence of the factor data.
2. **Semantics of my map code:** `map_block(M, L, R)` computes
   `out[i2][k2] = Σ_{i,k} L[i][i2]·R[k2][k]·M[i][k]`, i.e. in standard matrix product
   form it is `Lᵀ·M_blk·Rᵀ`. My textual transcription "`u2 = P⁻ᵀ·U·Q⁻¹`" described the
   DELTA-COLLAPSE derivation (symbols-as-matrices); the implemented census family in
   standard semantics is
     **U-side: `P⁻¹·U_blk·Q⁻ᵀ`  V-side: `Qᵀ·V_blk·R⁻ᵀ`  W-side: `Pᵀ·W_blk·R`**,
   which is the same automorphism family as the derived text with generators
   `(Pᵀ,Qᵀ,Rᵀ)` — one family, two transcriptions.
3. **Unique-family proof (independent, both parties):** replaying Main's 4⁶ = 4096
   letter enumeration on the shared arrays yields EXACTLY ONE all-ternary AND
   Brent-729/729 assignment: `(G⁻ᵀ, G⁻¹) | (G, G⁻¹) | (G, Gᵀ)` — my family. Main's
   independent enumeration agreed on uniqueness; the letter-list difference was pure
   L-vs-Lᵀ transcription, verified bit-for-bit (my winning U-block equals Main's
   `(G⁻¹, G⁻ᵀ)` under standard matmul; same for V and W).
4. **Main's original non-ternary replay** used standard `G·W·Gᵀ` for the W-side — the
   roles-swapped variant, which my `w_swap` probe reproduced byte-for-byte including
   all cited rows (max|entry| = 4). The swapped variant is NOT the unique family and is
   not Brent-relevant (my swapped-variant Brent test on paper55 fails; Main's ternarity
   rows came from this variant).
5. **No second family exists:** the transpose-dual concern (Main) is absorbed because
   the dual of my family is my family at transposed generators — the SAME unique
   survivor. There is exactly one ternary+Brent 3-letter diagonal family on these
   arrays, found independently by both parties.

## Reference convention (fixing item 2)

The **frozen `gatec_sweep` reshape** (row-major 3×3 blocks inside the 23×9 rows;
Brent index `U[r][3i+k]`) is the repo's semantics. My implementation is verified
IDENTICAL to the frozen sandwich at monomials bit-for-bit on all three blocks — only
possible under the same reshape. All census outputs (576/55 diagonal; the off-diagonal
tables) are therefore **comparable to the frozen landscape and to the published 55**
without conversion.

## Status of the numbers

- 576 valid sun56 diagonal orientations (48 monomial + 528 non-monomial), all
  ternary AND Brent-valid under the unique family in the frozen convention:
  **unconditional**. Floor-decided: min certified total = 55, zero ≤ 54
  (`sun56_diag_honest2_floor.json`). Per rule 5: my earlier
  `provenance_resolution.md` named the W-swap as the replay defect; the precise
  mechanism is the L/Lᵀ semantics documented here — the byte-reproduction stands,
  the naming is corrected.
- The record-attack's exclusion of the 6912 remains **overturned as stated**: under
  the correct automorphism family, non-monomial diagonal survivors exist (528 on
  sun56) and are now CERTIFIED ≥ 55 — the gap it exposed is real but safe on the
  diagonal slice. The off-diagonal census proceeds under the same family.
