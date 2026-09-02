# Amendment 1 — the authoritative primary read is the RASTERIZED scan, not the text layer

**Created 2026-09-02, AFTER the first (aborted) run and BEFORE any amended
compute.** Committed path-scoped before the amended runner executes.
Campaign: `20260902T022208Z_8c279ad1_d0276dfee63d`, gate
`laderman23-source-lock`, prereg sha
`154f807854c7335587aea99a9427b128a801dfc6cbe2626297ee9597c2210c79`.

## What happened (defect in MY instrument, not in the source)

The first run aborted at control **A1** exactly as pre-registered: the
transcription satisfied only **727 of 729** Brent identities, both instruments
agreeing on the count (int = fmpz = 2 failures). **No verdict was produced and
no claim was made.** The aborted output is preserved verbatim as
`verdict_preamendment_aborted.json`, and the pre-amendment runner as
`laderman_lock_run.py.preamendment`.

Root cause: pre-statement §0 recorded that the primary is a **page scan** whose
text layer was extracted with `pdftotext -layout`. That text layer **truncates
the final term of two long products at the right page edge**, rendering both
`m_3` and `m_11` as `... - b_31 + b_3` with the second subscript digit lost:

    m   =                           +
     ll     «32(- *U    *13 + Z?21 - ^22 - ^23 ~ *31 + *3

I typed `b_33` for both. The primary scan actually prints **`b_33` for `m_3`**
and **`b_32` for `m_11`**. That single wrong coefficient is the entire 2-identity
failure.

This is a **transcription defect in my own instrument**, caught by the
pre-registered ACCEPT control before any claim — it is NOT a fact about
Laderman's paper, and pre-statement §6's FROZEN-NEGATIVE branch ("no
trustworthy source is reachable, or the primary transcription fails Brent") must
not be read as licensing the claim "Laderman's published scheme fails Brent".
Freezing that would be exactly the overclaim this project forbids: my typing
failed, the source did not. Hence this amendment rather than a negative verdict.

## What changes (method, fixed before the amended run)

1. **Authoritative read.** The primary source of truth for the transcription is
   the **rasterized scan** of the pinned PDF (`pdftoppm -r 300/400 -png`), read
   line by line. The `pdftotext -layout` text layer is DEMOTED to a
   cross-checking aid; it is known-lossy at the page edge.
2. **Every line re-read.** All 23 products and all 9 output equations were
   re-read from the rasterized primary before this amendment was written. The
   readings are recorded verbatim in §"Recorded primary read" below, and the
   rasterizations themselves are pinned in `sources/raster/` with sha256
   prefixes.
3. **New control A5 (added, gating).** The runner must assert that its
   transcription table equals, term for term and sign for sign, the recorded
   primary read below. Any drift aborts.
4. **The single coefficient corrected**: `m_11`'s last term `b_33 -> b_32`.
   Nothing else in the transcription changes.

Everything else in the pre-statement stands unchanged: convention, both
instruments, controls A1-A4 and R1-R4, the reproduction cross-check being a
non-gating control, the 0.5 CPU-h budget, the verdict rule, and the explicit
list of what this phase may never claim.

## Recorded primary read (from the rasterized scan, pages 127-128)

    m1  = (a11 + a12 + a13 - a21 - a22 - a32 - a33) b22
    m2  = (a11 - a21)(-b12 + b22)
    m3  = a22(-b11 + b12 + b21 - b22 - b23 - b31 + b33)
    m4  = (-a11 + a21 + a22)(b11 - b12 + b22)
    m5  = (a21 + a22)(-b11 + b12)
    m6  = a11 b11
    m7  = (-a11 + a31 + a32)(b11 - b13 + b23)
    m8  = (-a11 + a31)(b13 - b23)
    m9  = (a31 + a32)(-b11 + b13)
    m10 = (a11 + a12 + a13 - a22 - a23 - a31 - a32) b23
    m11 = a32(-b11 + b13 + b21 - b22 - b23 - b31 + b32)      <-- b32, not b33
    m12 = (-a13 + a32 + a33)(b22 + b31 - b32)
    m13 = (a13 - a33)(b22 - b32)
    m14 = a13 b31
    m15 = (a32 + a33)(-b31 + b32)
    m16 = (-a13 + a22 + a23)(b23 + b31 - b33)
    m17 = (a13 - a23)(b23 - b33)
    m18 = (a22 + a23)(-b31 + b33)
    m19 = a12 b21
    m20 = a23 b32
    m21 = a21 b13
    m22 = a31 b12
    m23 = a33 b33

    C11 = m6 + m14 + m19
    C12 = m1 + m4 + m5 + m6 + m12 + m14 + m15
    C13 = m6 + m7 + m9 + m10 + m14 + m16 + m18
    C21 = m2 + m3 + m4 + m6 + m14 + m16 + m17
    C22 = m2 + m4 + m5 + m6 + m20
    C23 = m14 + m16 + m17 + m18 + m21
    C31 = m6 + m7 + m8 + m11 + m12 + m13 + m14
    C32 = m12 + m13 + m14 + m15 + m22
    C33 = m6 + m7 + m8 + m9 + m23

All output equations carry `+` on every term, exactly as printed.

## Independent corroboration of the corrected glyph (recorded, not gating)

- The pinned independent reproduction (arXiv:1108.2830 §2.4) also has `b_32`
  as `P11`'s last term and `b_33`-family signs for `P03`, i.e. it agrees with
  the rasterized primary and not with the truncated text layer.
- The coefficient is additionally **forced algebraically**: with `m6, m7, m8,
  m12, m13, m14` fixed as printed, `C31 = a31 b11 + a32 b21 + a33 b31` forces
  `m11 = a32(-b11 + b13 + b21 - b22 - b23 - b31 + b32)` uniquely. The same
  argument independently forces `m3`'s last term to be `+b33`, confirming that
  the two truncated tails genuinely differ.

Neither corroboration is used as evidence for the verdict: the verdict rests on
the in-run 729/729 dual-arithmetic check of the transcription recorded above.
