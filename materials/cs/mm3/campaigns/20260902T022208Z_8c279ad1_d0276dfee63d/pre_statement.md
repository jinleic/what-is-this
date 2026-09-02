# Pre-statement — `laderman23_source_lock` (Phase 1: source lock only, no search)

**Created 2026-09-02, BEFORE any target compute of this campaign.**
Target: `cs/mm3/`. Gate: `laderman23-source-lock`. Agent: `Mm3NextGap`
(subagent of Main). Committed path-scoped BEFORE `campaign.py init`; the run
dir receives a byte-identical copy with source commit + sha256 recorded.

## 0. Source reads performed BEFORE this statement (disclosed, per rule 5)

Two documents were fetched first-hand and are hash-pinned NOW, before any
compute of this campaign:

1. **PRIMARY.** Julian D. Laderman, *"A noncommutative algorithm for
   multiplying 3x3 matrices using 23 multiplications"*, Bulletin of the
   American Mathematical Society **82**(1), January 1976, 126-128.
   Publisher PDF (open access), fetched 2026-09-02 from
   `https://www.ams.org/journals/bull/1976-82-01/S0002-9904-1976-13988-2/S0002-9904-1976-13988-2.pdf`
   — 236,490 bytes, **sha256
   `a1deb200f2e270b13b870bbcf8ce8c669d7a09a4acdeddc012df5d6c90042a5b`**.
   This is a page scan; the text layer was extracted with `pdftotext -layout`
   (poppler) and pages 2-3 were also rasterized at 400 dpi for human reading.
2. **INDEPENDENT REPRODUCTION (control only).** N. T. Courtois, G. V. Bard,
   D. Hulme, *"A New General-Purpose Method to Multiply 3x3 Matrices Using Only
   23 Multiplications"*, arXiv:1108.2830, Section 2.4 "The Laderman Solution
   From 1975", which prints the scheme in directly machine-checkable Maple
   form. Fetched 2026-09-02 from `https://arxiv.org/pdf/1108.2830` —
   131,107 bytes, **sha256
   `d768596b226198996170cc645814ebef93eac0181291c54db34ae21b95b78707`**.

**Scouting disclosure (pre-registration honesty).** While reading the primary
scan I hand-expanded **three** of its nine output equations (C13, C23, C32)
against the printed products to confirm the text layer was legible and my
reading of the ambiguous glyphs was right; all three closed exactly. I also
observed that the reproduction's product-index supports appear to match the
primary's printed output equations. **None of that is evidence**: every claim
this campaign makes is re-derived in-run by the pinned code path, and the
hand-checks are superseded by the full 729-identity verification below.

No factor matrices for Laderman-23 exist anywhere on disk (session 10 recorded
"`laderman23_source_lock` (no verified factors on disk)"); nothing is inherited.

## 1. Fixed question and falsifiable outcome

**Is a trustworthy, verbatim-transcribed, machine-verified rank-23
decomposition of the 3x3 matrix-multiplication tensor obtainable from
Laderman's own published paper — and does it satisfy the exact Brent identity
that the paper itself states as its correctness criterion?**

The paper states its own criterion on page 128: an integer solution of "the
following system of 729 nonlinear algebraic equations involving 621 unknowns",
`sum_t x_{ght} y_{rst} z_{ikt} = delta_{ig} delta_{hr} delta_{sk}`. That is
exactly the Brent identity this target already uses. So the paper's own
acceptance test is the one applied here.

Falsifiable: if the transcription fails even one of the 729 identities in either
arithmetic, the source lock FAILS and Phase 2 does not start.

## 2. Model, conventions, and the transcription contract (fixed)

- **Factor convention (identical to the frozen target convention, reused
  verbatim from `campaigns/2026-08-30T113838Z_.../scripts/verify_anchors.py`):**
  `U[r][3i+k]` = coefficient of `a_{i+1,k+1}` in product r's left operand;
  `V[r][3k+j]` = coefficient of `b_{k+1,j+1}` in product r's right operand;
  `W[r][3i+j]` = coefficient of product r in output `C_{i+1,j+1}`.
  Brent: `sum_r U[r][3i+k] V[r][3k'+j] W[r][3i'+j'] = [i=i'][j=j'][k=k']`,
  729 identities, all indices 1..3, row-major.
- **Verbatim transcription.** The 23 products and the 9 output combinations are
  typed from the primary's own printed formulas, with the printed sign of every
  term. No re-derivation, no simplification, no reordering, no sign
  normalization. The transcription is stored as an explicit literal table in
  the runner and hashed.
- **Provenance labels (fixed now).** The *source text* — that these are
  Laderman's published formulas — is **CITED-DEPENDENCY** (owner-readable
  hash-pinned PDF; this campaign does not and cannot machine-verify
  authorship). That the transcribed triple *is* a valid rank-23 decomposition
  of the 3x3 matmul tensor over Z is **MACHINE-VERIFIED** (729/729, two
  arithmetics). These labels are not interchangeable and the report must carry
  both.
- **No search, no synthesis, no circuit model in this phase.** Addition counts,
  stage costs, floors, auxiliary universes, SAT and DFS instruments are all
  OUT OF SCOPE for Phase 1. Phase 1 produces exactly one artifact: a
  hash-pinned, Brent-verified factor triple plus its provenance record.

## 3. Instruments (pinned before use)

- **Instrument 1 — pure-integer Brent evaluator.** Python `int` only, exact,
  no floats, no solver: the 6-fold loop above evaluated directly.
- **Instrument 2 — `fmpz` Brent evaluator.** The same loop over
  `python-flint` `fmpz` integers (flint 0.9.0), an independently implemented
  arithmetic. **Both must agree on every one of the 729 identities**; any
  disagreement aborts with no verdict (FROZEN-INCONCLUSIVE).
- No SAT solver, no DFS, no certificates are needed or used in Phase 1.

## 4. Controls — both directions, in-run, AFTER init, BEFORE any claim

**ACCEPT (must pass):**
- **A1.** The primary transcription satisfies **729/729** Brent identities in
  BOTH arithmetics.
- **A2.** Frozen-row re-verification through the SAME code path: the frozen
  `paper55` factor triple (loaded from the frozen
  `scripts/gatec_decomps.py::load_paper55`) must also give 729/729 in both
  arithmetics. A known-true object must be accepted by this verifier.
- **A3.** Structural fingerprint against the primary's printed page: the 7
  single-product terms must be exactly
  `m6 = a11*b11`, `m14 = a13*b31`, `m19 = a12*b21`, `m20 = a23*b32`,
  `m21 = a21*b13`, `m22 = a31*b12`, `m23 = a33*b33`, and the nine output
  supports must be exactly the printed index sets
  C11{6,14,19}, C12{1,4,5,6,12,14,15}, C13{6,7,9,10,14,16,18},
  C21{2,3,4,6,14,16,17}, C22{2,4,5,6,20}, C23{14,16,17,18,21},
  C31{6,7,8,11,12,13,14}, C32{12,13,14,15,22}, C33{6,7,8,9,23}.
  A mismatch means my transcription drifted from the printed page and aborts.
- **A4.** Alphabet and rank audit: every U, V, W entry lies in {-1,0,1};
  exactly 23 products; no all-zero product row on any side.

**REJECT (must fail, and the failure count must be reported):**
- **R1.** One sign flipped in one U entry — must fail Brent with a nonzero
  count in both arithmetics.
- **R2.** One product deleted (its W row zeroed, leaving 22 effective
  products) — must fail.
- **R3.** Two products' V rows swapped — must fail (guards against a verifier
  that is insensitive to product identity).
- **R4.** A planted assertion error (deliberately wrong expected 729 count)
  must fire, proving the harness's own asserts are live.

**Transcription cross-check (reported, not gating).** The independent
reproduction (source 2) is transcribed separately and (a) Brent-checked through
the same instruments, and (b) compared to the primary on the nine output
supports. Its per-product SIGNS are expected to differ (negating a product's
operand and its output coefficients is a free equivalence); a support
disagreement or a Brent failure on the reproduction is reported as a finding
about the *reproduction*, never as a defect of the primary, and never blocks a
verdict that rests on the primary alone.

## 5. Budget, seeds, determinism

- Hard cap **0.5 CPU-hours** total. Expected actual: seconds (729 identities
  x 2 arithmetics x a handful of objects). `nice -n 15`,
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1`,
  `PYTHONDONTWRITEBYTECODE=1`, and `sys.dont_write_bytecode = True` before
  importing anything from a frozen campaign dir.
- Fully deterministic: no randomness, no sampling, no timing-dependent
  behaviour. No seeds are needed; none are used.
- No checkpointing needed at this scale; the runner writes `verdict.json` and
  `runner_output.json` once, at the end.

## 6. Verdict rule (fixed before compute) — exactly one terminal verdict

- **FROZEN-CERTIFIED** iff A1-A4 all pass and R1-R4 all fire correctly. Claim
  then permitted, and no more: *"the factor triple transcribed verbatim from
  the hash-pinned primary is a rank-23 decomposition of the 3x3
  matrix-multiplication tensor over Z, machine-verified 729/729 in two
  independent arithmetics; its identification as Laderman's published scheme is
  a CITED-DEPENDENCY on the pinned PDF."* This unlocks Phase 2, which requires
  its own fresh pre-registration and commit.
- **FROZEN-NEGATIVE** iff no trustworthy source is reachable, or the primary
  transcription fails Brent. The campaign then freezes the failure with the
  exact reachable-source list and the exact failing identities, and Phase 2
  does NOT start.
- **FROZEN-INCONCLUSIVE** on any control failure, instrument disagreement, or
  budget exhaustion.
- What this phase may NEVER claim: any addition count, any stage cost, any
  optimality, any statement about Laderman's minimum, any comparison to the
  published operation count, or authorship verification. Those are Phase 2 or
  out of scope entirely.

## 7. Reproduction command (fixed)

    cd cs/mm3/<run_dir>
    OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 15 \
      ../../../cs/.venv/bin/python laderman_lock_run.py

Phases: A anchors (transcribe + Brent both arithmetics) -> B controls
(A2-A4, R1-R4) -> C reproduction cross-check -> D verdict. Aborts (exit 1) on
any control failure or instrument disagreement.
