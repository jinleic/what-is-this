# Route A pre-statement — Strassen 3-slice commutator theorem: primary import
# and clause-by-clause hypothesis check for T_O-relevant 3-slice L-families

Campaign: `campaigns/2026-09-01T03:21:00Z_routeA_prov/` (directory named with
the ACTUAL UTC start time of the campaign, 2026-09-01T03:21Z).
Owner agent: `OctRankRouteA`.

**Committed BEFORE any load-bearing computation of this campaign.** Rule 16:
domains fixed here; no domain-shopping. Any re-scoping goes into an appended
"Re-scoping log" with the original text retained.

## 0. What this campaign is and is not

Object: the 3-slice lower bound available for the L-family slices of T_O via
the commutator machinery of Strassen's 1983 paper, with the theorem imported
FIRST-HAND from the primary text (not from Landsberg's survey, not from
Koiran's restatement — those were the provenance chain of the Route AF
campaign, which found Primary-provenance NOT ESTABLISHED; this campaign
closes exactly that gap or freezes it).

Questions to settle, in order:

A1 (provenance). Read Strassen 1983 first-hand; quote the theorem(s)
   verbatim, with page, and record how the text was obtained (source, hash).
A2 (hypotheses). For EACH hypothesis of each quoted theorem, check it
   clause by clause against our tensors ((L_u, L_v, L_w) 3-slice families
   from T_O in the fixed Cayley-Dickson basis; the specific tensor
   tau ⊠ s = block-diagonal doubled tau, and tau itself as control),
   with the verification named (exact arithmetic vs. structural argument).
A3 (consequence). Derive the exact implied number (if any) for
   rank(L_u, L_v, L_w) from the theorem AS STATED, and re-audit the tau
   control at n = 4 with the same text. The outcome must be consistent with
   the known true value R_R(tau) = 7: a reading that gives 8 at tau is
   DEAD and is reported as refuted, not used.
A4 (upper discipline). NOTHING in this campaign touches the published
   window 18 <= R_R(T_O) <= 25. If any reading produced here appeared to
   bear on that window, escalate to Main BEFORE recording any claim; no
   refutation of published work is claimed by this campaign.

Pre-registered non-goals: no new numerical searches; no Route-F 13-vs-14
compute (that is the second campaign, only if Route A resolves or
hard-blocks); no claims about R_R(T_O) itself; no claims about dependent
triples beyond what the algebra forces.

## 1. Fixed domains (rule 16)

D1. Algebra: Cayley-Dickson octonions on H(+)H, basis
    (1, i, j, k, l, i.l, j.l, k.l), product
    (a,b)(c,d) = (ac - conj(d)b, da + b conj(c)) — byte-compatible with
    gate_a_verify.py A1 and the S3/RouteAF campaigns.
D2. Tensor class: T_{u,v,w} = three-slice real tensor with slices
    (L_u, L_v, L_w), L_x[c,b] = (x e_b)_c, for (u,v,w) ranging over
    linearly independent triples in R^8 (the S3-declared class).
D3. Control class: tau = (L_1, L_i, L_j) on the quaternion subalgebra
    H' = span(e0..e3) (4-dimensional slices), true rank exactly 7
    [MACHINE-VERIFIED in gate C; replayed by anchors].
D4. The specific conjugated target tau ⊠ s (slices blockdiag copies of
    tau slices) is admitted as a SECONDARY control for the
    block-duplication structure; its rank is OPEN in [13,14] — it is NOT a
    known value and MUST NOT be used as a calibration fit. It is used only
    as a shape on which the theorem's hypotheses must be checked.

No other domains. If a route needs a new domain, it is a Re-scoping log
entry, reported alongside, not a replacement.

## 2. Arithmetic discipline

- All matrix ranks and identities in fmpq/sympy exact arithmetic.
- Floats allowed only for orientation; no float may support any claim.
- Hypothesis checks of structural type (e.g. "dim W = 3") are written
  arguments pinned to exact dimension computations.

## 3. Verbatim quotes required in the theorem import (A1)

The pre-statement is written BEFORE the theorem text was fully digested
(the acquisition is ongoing; the PDF has been obtained and rasterized).
The VERBATIM theorem statements, with PDF page and transcription method
(text layer / pixel-verified image read), go into `theorem_import.md`
in the NEXT write, after the reading, and the campaign compute follows
only after that file exists. This staging is deliberate: the hunt for the
text (Wayback/CORE route) is provenance work and preceded this file only
as acquisition, not analysis.

## 4. Adjudication (fixed now)

- A2 find = every clause of THE theorem (the exact statement, not a
  paraphrase) verified for (L_u, L_v, L_w): the campaign derives the
  certified lower bound N_theorem for the class, and reports it with its
  label (the import is CITED-DEPENDENCY with primary read = rule-1 clean;
  the clause checks here are MACHINE-VERIFIED / HUMAN-AUDITED as labeled).
- tau misreading = refuted reading; reported as FAILED (rule 5 retraction
  mechanism; the naked 16 was already retracted by Route AF).
- If the primary text is found to state the inequality for border rank
  only (R with underbar), the campaign reports the border-rank bound and
  the real-rank consequence separately, naming the complexification
  argument explicitly with its exact compatibility (real tensor rank >=
  complex border rank after base change) — and labels that step
  CITED-DEPENDENCY [INFERENCE] unless proven in-text.
- Rule 7 scope sentence drafted at close; every route's numeric outcome,
  including misses, reported.

## 5. Budget (fixed)

- Reading + transcription: no cap (provenance gate).
- A2/A3 compute: <= 60 min CPU (in-process time.process_time() accounting
  per Main's ruling; no ps-derived numbers anywhere).
- No heavy unbounded processes; single-threaded BLAS.

## 6. Rule-7 scope sentence (skeleton, numbers at close)

"Covered: [clauses verified on D2/D3/D4 with exact numbers]. NOT swept:
[...]."
