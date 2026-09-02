# Campaign report — gate `laderman23-source-lock` (Phase 1) (run 20260902T022208Z_8c279ad1_d0276dfee63d)

## Question (prereg sha `154f8078…`, commit `425c8a5`; amendment 1 commit `cb319be`)

Is a trustworthy, verbatim-transcribed, machine-verified rank-23 decomposition
of the 3x3 matmul tensor obtainable from **Laderman's own published paper**, and
does it satisfy the exact 729-identity Brent system that the paper itself states
as its correctness criterion (p.128)?

## Answer — verdict FROZEN-CERTIFIED (`CERTIFIED-SOURCE-LOCK`)

**Yes.** The primary was obtained, hash-pinned, transcribed verbatim, and
**729/729 Brent identities hold exactly** in two independent arithmetics that
agree on every identity.

## The locked source

| role | document | bytes | sha256 |
|---|---|---|---|
| **PRIMARY** | J. D. Laderman, *"A noncommutative algorithm for multiplying 3x3 matrices using 23 multiplications"*, **Bull. Amer. Math. Soc. 82(1), Jan 1976, 126-128** (open-access publisher PDF) | 236,490 | `a1deb200f2e270b13b870bbcf8ce8c669d7a09a4acdeddc012df5d6c90042a5b` |
| reproduction (control) | Courtois–Bard–Hulme, arXiv:1108.2830 §2.4 "The Laderman Solution From 1975" (Maple-checkable form) | 131,107 | `d768596b226198996170cc645814ebef93eac0181291c54db34ae21b95b78707` |

Both are stored in the run dir under `sources/`, together with the pinned
300/400 dpi rasterizations (`sources/raster/`) that are the **authoritative
read** per amendment 1, and the lossy `pdftotext -layout` extraction kept only
as a cross-check.

## Provenance labels (both required, not interchangeable)

- *"These formulas are Laderman's published scheme"* — **CITED-DEPENDENCY**
  on the hash-pinned open-access publisher PDF. This campaign does not and
  cannot machine-verify authorship or publication.
- *"The transcribed triple is a valid rank-23 decomposition of the 3x3 matmul
  tensor over Z"* — **MACHINE-VERIFIED** in-run: 729/729, dual arithmetic.
- *"The transcription matches the printed page"* — **MACHINE-VERIFIED**: the 7
  single-product terms and all 9 output supports equal the printed equations
  (control A3), and every product and output equation equals the committed
  amendment's recorded primary read term-for-term and sign-for-sign (A5).

## The locked object

Convention (identical to this target's frozen convention, reused verbatim):
`U[r][3i+k]` = coeff of `a_{i+1,k+1}`; `V[r][3k+j]` = coeff of `b_{k+1,j+1}`;
`W[r][3i+j]` = coeff of product r in `C_{i+1,j+1}`; Brent
`sum_r U[r][3i+k] V[r][3k'+j] W[r][3i'+j'] = [i=i'][j=j'][k=k']`.

Triple sha256 **`522ba07f4f1784ac…`**, alphabet exactly {-1, 0, 1}, 23 products,
no zero rows on any side. Stored as `laderman23_factors.json` together with the
printed products and output equations.

## Controls — all passed, both directions

| control | result |
|---|---|
| **A1** primary transcription, Brent | **729/729**, int and fmpz agreeing on every identity |
| **A2** frozen-row re-verification (frozen `paper55` triple through the SAME path) | 729/729 — a known-true object is accepted |
| **A3** printed-page fingerprint | 7 single products + all 9 output supports match the page |
| **A4** alphabet / rank audit | alphabet {-1,0,1}, 23 products, 0 zero rows |
| **A5** (amendment 1) equality with the committed recorded primary read | 23/23 products, 9/9 outputs, zero drift |
| **R1** one U sign flipped | REJECTED (1 identity broken) |
| **R2** one product deleted | REJECTED (7 identities broken) |
| **R3** two V rows swapped | REJECTED (35 identities broken) |
| **R4** planted wrong count | assertion fired |

## Reproduction cross-check (non-gating, and a genuinely interesting finding)

The independent reproduction transcribes to a **different but equally valid**
point: it also satisfies **729/729**, and its nine output supports are
**identical 9/9** to the primary's printed equations, while **per-product signs
differ** (e.g. it prints `P06 = a11*(-b11)` and compensates with `C11 = -P06 +
P14 + P19`). Negating a product's operand together with its output coefficients
is a free equivalence, so the two are the same algorithm presented with
different sign conventions. The support agreement across two independently
typeset sources is strong corroboration that the transcription is Laderman's
scheme and not some other rank-23 decomposition.

## Instrument defects (disclosed; both caught by pre-registered controls, no claim affected)

- **D1 — transcription defect, amended.** The first run aborted at control A1
  with **727/729** (both instruments agreeing on 2 failures) and produced **no
  verdict**. Cause: the `pdftotext` text layer truncates the final term of the
  two long products `m_3` and `m_11` at the right page edge (`... - b_31 + b_3`)
  and I typed `b_33` for both; the scan prints **`b_33` for `m_3`** and
  **`b_32` for `m_11`**. Amendment 1 was written and committed BEFORE the
  amended run: it demotes the text layer, makes the rasterized scan
  authoritative, records the full line-by-line read of all 23 products and 9
  output equations, corrects that one coefficient, and adds gating control A5.
  The corrected glyph is independently corroborated (the reproduction also has
  `b_32`, and `C31 = a31b11 + a32b21 + a33b31` forces it uniquely) — but the
  verdict rests only on the in-run 729/729 check. The aborted verdict and the
  pre-amendment runner are preserved verbatim
  (`verdict_preamendment_aborted.json`, `laderman_lock_run.py.preamendment`).
  **The prereg's FROZEN-NEGATIVE branch was deliberately NOT invoked**: my
  typing failed, the source did not, and freezing "Laderman's scheme fails
  Brent" would have been exactly the kind of overclaim this project forbids.
- **D2 — control A5's own parser, fixed in place.** On its first execution A5
  tokenized a trailing prose annotation (`<-- b32, not b33`) as two extra
  formula terms and aborted. The parser now truncates annotations before
  tokenizing. Neither the recorded read nor the transcription changed; A5 found
  its own bug and refused to pass, which is the correct behaviour.

## Cost

Measured **0.15 s** wall / well under 1 CPU-second, of a pre-registered 0.5
CPU-hour cap. Fully deterministic; no randomness, no sampling, no seeds.

## What this phase explicitly does NOT claim

No addition count, no stage cost, no floor, no optimality or minimality, no
comparison to any published operation count, and no authorship verification.
Those are Phase 2 (which requires its own fresh pre-registration and commit) or
out of scope entirely.

## exactly-one verdict

**FROZEN-CERTIFIED** — the source lock holds: Laderman's own paper is pinned,
the transcription is verbatim and matches the printed page, and the resulting
rank-23 triple is machine-verified 729/729 over Z in two independent
arithmetics. **Phase 2 is unlocked** and requires a fresh pre-registration.
