# Report — Gate-B certified close-out restart: three full bands PASS, three FAILURE TO CERTIFY

Campaign: `campaigns/20260831T082425Z_kg_direct_d4_restart/`.
Pre-statement first-file commit: `87a43d2b9a9cd30a58b8ecdad74375719d0635aa`.
The stopped predecessor campaign `20260831T080542Z_kg_band_close2` was left
immutable.

## Headline

All analytic and direct-envelope startup controls passed, and the preregistered
six-band sweep then ran to its declared protocol.  The band table is now:

| band | verdict | worst certified leaf margin (Arb 256) | tiles certified/open | envelope evals |
|---|---|---|---|---|
| `[1.0,1.3]` | **PASS** | `+9.49573984471527704909052347604e-7` | 14/0 | 13,718 |
| `[1.30,1.45]` | **PASS** | `+1.95705024245827645177570572783e-6` | 6/0 | 6,340 |
| `[1.45,1.75]` | **PASS** | `+3.10559377267758493229905832146e-6` | 10/0 | 11,650 |
| `[1.75,3.5]` | **PARTIAL — FAILURE TO CERTIFY** | `+2.30038358179570808985477029928e-7` | 25/11 | 57,551 |
| `[3.50,4.083]` | **FAILURE TO CERTIFY** | — | 0/8 | 3,900 |
| `[4.083,6.0]` | **FAILURE TO CERTIFY** | — | 0/20 | 6,016 |

All margins are **MACHINE-VERIFIED** Arb `256`-bit values with zero added
ball radius; the run counters are `COMPUTATIONAL-EVIDENCE`.  PASS means every
ratio-`1.02` c-tile of the band is certified by separated endpoints over the
whole preregistered parameter cover; it is not a claim inherited from any prior
campaign and it exactly contradicts nothing in the paper.  The three non-PASS
bands are upper-envelope failures (unresolved leaves), not refutations; the
paper's refutation trigger never fired.

## Mandatory startup controls, all passed

Recording facts first, then their strengths:

1. Analytic unit-cubic control.  The center
   `p_*(s)=(4/5)psi_0-(3/5)psi_2=A-Bs^2`, with
   `A=4/5+3/(5sqrt(2))`, `B=3/(5sqrt(2))`, satisfies
   `(4/5)^2+(3/5)^2=1`; on `[0,infinity)` the active set splits at the roots
   `r1=(sqrt(c^2+4AB)-c)/(2B)`, `r2=(sqrt(c^2+4AB)+c)/(2B)`, and the closed
   form uses only `P(r1)`, `Q(r2)`, densities, and
   `r phi(r)` corrections summed over both active components.
   At 300 Arb bits: `J(1.30,p_*)=[0.374830554588763570024533304529760488130
   +/- 4.49e-62]`; `d(1.45)-J=[-0.00214885479421366092840253289474137
   +/- 4.59e-63]`; `d(1.30)-J=[0.0322365204873065325951173054080068
   +/- 2.01e-62]`.  **MACHINE-VERIFIED.**  An independent 80-digit mpmath
   quadrature gave 0.3748305545887635700245333045297604881300360210880756,
   matching; it is `COMPUTATIONAL-EVIDENCE` only, not campaign evidence.
2. Root semantics: exact sign signs at `0,r1,r2` upon the two quadratic
   factors and ordering `0<r1<sqrt(A/B)<r2`, with `A-B=4/5` exactly.
   **MACHINE-VERIFIED.**
3. Rule-17b count: the exact `12x12` rational cover has exactly `132`
   disk-intersecting boxes, asserted every tile run.  **MACHINE-VERIFIED
   against the independently pre-anchored count.**
4. Weights/plants: `phi(1)<Phi(1)-Phi(0)<phi(0)`; the pseudo-upper
   `J-1e-6` is rejected; margin plants around the actual `d(1)` all behave.
   **MACHINE-VERIFIED.**
5. Direct-envelope anchors on the `0.0001`-radius cell at `(0.8,-0.6)`,
   `512` panels: equality c-pair margin
   `+0.02515143610672980214527072481596129` and the ratio-`1.02` pair
   `(1.30,1.326)` margin `+0.0187889631784633691674141804635343`, both zero
   additional radius; the analytic ball is contained (`<= +2.81e-51` width)
   and rejection holds.  **MACHINE-VERIFIED.**
6. The preregistered coarse pair is declared impossible as a PASS anchor,
   since even the center cubic has margin `-0.00214885...` there.  It was
   therefore never used as a positive control.

## What each band run establishes

For a PASS band: for every adjacent ratio-`1.02` tile `[c_L,c_U]`, every
closed rectangle in the 132-root even-disk cover (refined to side
`>= 1/1024` where needed), and every unit vector compatible with the odd
bound `sqrt(1-dist(0,B)^2)`, the direct panel envelope

`(E+W-c_L a)_+ + (|E-W|-c_L a)_+`, exact weight `Phi(b)-Phi(a)`,

after the `[8,infinity)` tail majorant, gives
`d(c_U).lower() - U.upper() > 0`.  By exact c-monotonicity of
`J(c,p)` and `d(c)`, this certifies the full named band over the full cover.
A single tile or cell never covers a whole band, and no c-monotonicity shortcut
across tiles was used.

Open tiles report their unresolved leaf, the floor-`1/1024` depth reached,
and the exact unresolved interval lower margin at 32,768 panels.  The band
`[1.75,3.5]` certifies 25/36 tiles; the defective ring-region unresolved
cell floor is around `(-0.87,-0.5 of a vector)`, matching prior campaigns'
uneven-corner observation but now with exact rectangles and no inherited
claim.  In `[4.083,6.0]` the unresolved margin nearest zero remains
`-2.63399664328978880111736054560e-7`; in `[3.50,4.083]` it is
`-1.59862092695915505593006178343e-5`.  No unresolved leaf was ever proven to
exceed `d` (no refutation trigger); all are envelope slack.

## Timing confounds, honestly registered

Recorded wall times per band (COMPUTATIONAL-EVIDENCE, never performance
evidence): `[1.0,1.3]` 2583.990 s (CSI, ~4-minute McEliece overlap);
`[1.30,1.45]` 1193.090 s; `[1.45,1.75]` 2034.288 s; `[1.75,3.5]` 9640.884 s
(CSI, minor pre-base Omega overlap 0.98 s); `[3.50,4.083]` 660.247 s;
`[4.083,6.0]` 1482.250 s (CSI, 52.9-second RS census plus Oct Route-F
numsearch overlap).  All heavy runs were `nice -n 10` with
OMP/OpenBLAS/MKL/vecLib/NumExpr pinned to 1 and single-slot ownership
confirmed by the parent before launch.  No correctness claim here depends on
any timing.

## Fixing the historical `+0.1930` defect

The inherited claim of a positive margin on the coarse pair `(1.30,1.45)` is
impossible, because the cell's own cubic has negative margin there.  This does
not contradict the paper: at equal `c=1.30`, `d(1.30)-J(1.30,p_*)` is
`+0.0322365...`, and sub-band equal-c and 1.02-ratio anchors both certify
positive margins.  No prior artifact had established the historical number,
and this campaign replaces it with the equal-c and ratio-1.02 checks instead
of pretending to reconstruct it.

## Yes, the paper is fully consistent

Nothing here refutes paper Lemma D.2/D.4 or Table 1; every `MACHINE-VERIFIED`
comparison in this campaign that the paper also asserted (tile-adjacent
endpoints plus the exact `J<=d` implied by (56) on the three closed bands) is
strictly consistent.  The several PASS minimums in the 1.45-to-3.5 region are
about three-number-magnitude smaller than Table-1's +2.26e-5 REPORTED figure
and should not be compared as paper-contradicting; they are envelope margins,
not analytic J values.

## Rule 7 scope sentence

This campaign swept exactly all unit cubics in the orthonormal-probabilist
Hermite span `psi_0,...,psi_3` through the full even-coefficient disk via its
exact 132-root cover and certified odd-radius bound, over the six full c-bands
`[1.0,1.3]`, `[1.30,1.45]`, `[1.45,1.75]`, `[1.75,3.5]`, `[3.50,4.083]`, and
`[4.083,6.0]`, using adjacent ratio-`1.02` c-tiles, floor side `1/1024`,
panel ladder to 32768, 256-bit Arb panels on `[0,8]`, and the closed
`[8,infinity)` tail majorant; it also evaluated only the explicit cubic
`p_*=(4/5)psi_0-(3/5)psi_2` at c `1.30` and the d values at `1.30`/`1.45`,
the radius-`0.0001` anchor cell at c-pairs `(1.30,1.30)` and `(1.30,1.326)`,
and startup control plants; it did **not** search c in `[0.993405,1.0)`,
c above `6.0`, coefficient cells outside the disk-cover floors, non-unit
cubics, degrees above three, the paper's C1/C3/splice certificates,
Gate-C septic/Krivine scheme families, or any construction outside the
paper's unit-cubic threshold family.

## Frozen evidence inventory

- `pre_statement.md` — preregistration, closes the stopped predecessor's stop
  at the impossible historical anchor and names every control above.
- `code/startup_controls.py`, `code/gate_direct.py`, `code/core.py` and their
  `asrun/` byte-copies — `frozen` on completion (checksums listed in
  `checksums.sha256`), including the direct Lemma-D.4 envelope source, CLI,
  and no post-run modification.
- Per-band artifacts (six `logs/band_*_result.json` files and progress logs)
  plus `artifacts/band_summary.json` with an independently cross-checked
  `99,175` envelope-evaluation count and band-by-band verdicts.
- Every band run used the same preregistered CLI, ratio, budget, and
  Arb precision.  No dry-run result was silently rolled backward.

## Named next campaign (priced)

`Gate-B near-vector refinement toward the 3.5-6 ring` — target closing the
`[3.5,4.083]`/`[4.083,6.0]` unresolved ring rectangles by adding the
per-panel odd-radius refinement `K_o` (only `phi(1.4-1.6)` region above)
at the same floor and direct-envelope arithmetic.  Based on this campaign's
floor depths/panel histogram, estimated cost is **2-4 CPU-hours**
(COMPUTATIONAL-EVIDENCE projection; not measured here).
