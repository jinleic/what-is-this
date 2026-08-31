# Provenance-resolution record — the ternarity dispute on witness gi=134 (settled)

**Timestamp 2026-08-30T20:4xZ. All claims below machine-checked in this campaign's
artifacts; scripts byte-frozen at wrap.**

## Dispute
Main's independent replay of the corrected honest map on witness
G = [[-1,-1,-1],[0,1,1],[0,0,1]] (sun56, diagonal) reproduced Brent 729/729 but
NOT ternarity: Main's U2/V2/W2 max|entry| = 4, mine = 1. Main suspended the 528
figure pending resolution (correctly, per rule 14: a measurement whose reference
is unvalidated is not evidence — the same lesson as miss #13).

## Resolution
Main's replay applied the W-block as `w2 = Rᵀ·W·P` (Gᵀ on the LEFT). The corrected
family's W-side is `w2 = P·W·Rᵀ` — G untransposed on the left, Gᵀ on the right.
Computed in this campaign (`crosscheck_sun56_G134.json`):

- `w_swap` variant (= Main's replay): W2 row0 = (2,−2,−2,−1,1,1,−1,1,1), max|entry| = 4
  — **byte-identical to Main's replayed rows, including all cited rows**.
- Corrected family: W2 row0 = (0,0,−1,0,0,−1,0,0,0), max|entry| = 1 (U2, V2, W2 all 1).

The witness's raw sun56 factor rows (for cross-loading) are frozen in
`crosscheck_sun56_G134.json`: U[0] = (0,0,0,0,1,0,0,1,1), U[1] = (0,0,0,1,−1,−1,1,−1,0),
V[0] = (−1,1,1,−1,1,1,0,0,0), W[0] = (0,0,1,0,0,0,0,0,1).

## Status after resolution
- The **528 figure STANDS**: 528 non-monomial sun56 diagonal orientations under the
  corrected honest map are ternary AND Brent-valid (729/729 over ℤ). Plus 48 monomial
  = 576 valid diagonal sun56 orientations under the corrected action (vs 48 under the
  frozen non-automorphism).
- **Floor decisions over all 576 completed** (109 s wall): **minimum certified total
  = 55, ZERO orientations at ≤ 54** (checkpoint `sun56_diag_honest2_floor.json`).
  Worst row gi=1068: sides (d=12, floor-impossible), (d=11, floor-impossible),
  (d=16, floor-achievable) → 13+12+16+14 = 55, matching the known sun56 structure.
- open question moved OFF-diagonal: the same corrected action on the full 6960³
  product is the census now launching.
- NOTE for the record: Brent 729/729 passing under Main's W-swapped variant is
  NOT expected (the swapped map is not an automorphism); if Main's replay truly
  passed Brent with the swapped W2, the replay's Brent and ternarity consumed
  different arrays — flagged to Main for a side-by-side re-run; does not affect
  the resolution above, which is exact-reproduction based.
