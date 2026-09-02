# Theorem import — Strassen 1983, first-hand transcription

Campaign: `2026-09-01T03:21:00Z_routeA_prov`. This file is written AFTER the
pre-statement commit (dd40d2b) and BEFORE any hypothesis-checking compute.

## Source provenance

- Primary text obtained 2026-09-01 (UTC), PDF, 41 pages:
  `sources/Strassen_1983_Rank_and_optimal_computation_of_generic_tensors_LAA52-53.pdf`
  SHA256 `f8481709de3a93e5747ce45e5c19aae81225139f39db347ef373a82df5bd0376`
  (1,936,012 bytes, PDF 1.2, application/pdf).
- Acquisition chain (all steps logged in sources/ACQUISITION_LOG.md and in
  the session transcript):
  1. ScienceDirect direct: HTTP 403 (both `/pdf` and `/pdfft` endpoints;
     bronze-OA per OpenAlex/S2 but Cloudflare-challenged).
  2. EuDML: zero hits (content-layer negative, evidence file retained).
  3. Wayback CDX: no `application/pdf` capture of any ScienceDirect PDF
     URL for this PII (only challenge-page HTML).
  4. scholar.archive.org (fatcat) listed a work whose archived location is
     CORE download/pdf/82110360; Wayback CDX for that URL had a 2019-04-16
     capture with mimetype `application/pdf`; fetched
     `https://web.archive.org/web/20190416143337if_/https://core.ac.uk/download/pdf/82110360.pdf`.
  5. Identity check of the fetched PDF: title page reads
     "Rank and Optimal Computation of Generic Tensors / V. Strassen /
     Institut für Angewandte Mathematik / Universität Zürich ... Dedicated
     to Alexander M. Ostrowski on the occasion of his ninetieth birthday /
     Submitted by Walter Gautschi", running head "LINEAR ALGEBRA AND ITS
     APPLICATIONS 52/53:645-685 (1983)"; Elsevier Science Publishing Co.
     This matches the bibliographic record (LAA 52/53 (1983) 645-685,
     DOI 10.1016/0024-3795(83)80041-X, Crossref metadata verified).
     Text layer: `sources/Strassen_1983_text.txt`
     SHA256 `b138a8b896c1b9a8056510adda3777aba70c1808c7ce7add2bba795b3e106e49`.
- IMPORTANT METHODOLOGICAL NOTE (glyph fidelity): the 1983 PDF's text layer
  is OCR-degraded (e.g. "inuertible" for "invertible", statements rendered
  with `~` for the inequality signs and `+` for the fraction bar). ALL
  displayed formulas quoted below were verified against rasterized page
  images, glyph by glyph (coordinates and comparison crops frozen in
  `verification_images/`), NOT taken from the text layer. Where the text
  layer and the image disagree, the image wins; the pixel evidence is
  cited per formula.

## Field convention of the paper (verbatim, page 648)

"Throughout this paper k denotes an algebraically closed field."
[text layer, clean ASCII, high confidence; consistent with (1.5) defining
border rank via Zariski closure over k.]

So Strassen 1983 works over an ALGEBRAICALLY CLOSED field k, and rank /
border rank are complex (k-) ranks. There is no real-rank statement in the
paper. (Section headers: "2. FIBRES OF THE COMPUTATION MAP", "3. PERFECT
SHAPE", "4. A DETERMINANT FOR 3-SLICE TENSORS".)

## THEOREM 4.1 (page 673; PDF page 29 of 41)

Setup (same page, text layer cleaned of OCR slips; setup confirmed in
image):

"4. A DETERMINANT FOR 3-SLICE TENSORS

Here we assume dim U = dim V = n, dim W = 3. We choose bases
(e_1, ..., e_n) for U, (f_1, ..., f_n) for V, (g_1, g_2, g_3) for W. If
t ∈ U ⊗ V ⊗ W we write [t = Σ_{i,j,l} α_ijl e_i ⊗ f_j ⊗ g_l  (4.1)]
and define the n × n matrices
A = (α_ij1),  B = (β_ij) ,  C = (γ_ij),
A, B, C are the three 'slices' of t."

[The α, β, γ superscript conventions in the printed definition are garbled
in the OCR text layer; from the display image A = (α_ijl) with distinct
third indices per slice; the PROOF uses slices A = l_1 = I_n, B = l_2,
C = l_3 w.r.t. the g-basis, consistent with standard usage. The exact
reading of the definition line is NOT load-bearing for us: only the
theorem's hypothesis (A invertible) and conclusion enter.]

Statement (image-verified; `verification_images/thm_display_wide.png`,
`lineA_thm_disp.png`, pixel checks in R_compare.png):

"THEOREM 4.1. Let A be invertible. Then

    R̲(t) ⩾ n + (1/2) rank(B A^{-1} C - C A^{-1} B)."

CRITICAL GLYPH FACTS (verified at pixel level, and cross-checked by two
independent image reads plus a layout-model read):

- The inequality is "greater or equal", pointing from the rank symbol
  toward n + ....
- The coefficient of the rank term is the literal fraction 1/2.
- The rank argument is the matrix B A^{-1} C - C A^{-1} B (plain products
  and minus; NO commutator brackets [ , ] printed).
- The rank symbol in the THEOREM DISPLAY is R with an UNDERBAR — Strassen's
  border-rank notation R̲. Evidence: the underbar segment is present below
  the R baseline in the theorem display (ink rows 1011-1014 below glyph
  rows 960-1000, columns 821-926 ≈ R-bowl span) and ABSENT in the proof's
  display (4.2) for the plain-rank R (0 below-baseline ink pixels in
  rows 1267-1292, columns 825-863). Side-by-side: `R_compare.png`,
  `micro_R.png`. Both the layout-image read and the adversarial second
  read concurred: theorem display = underline; (4.2) = plain R.
- The PROOF's first display, tagged (4.2), is the PLAIN-rank version:

      R(t) ⩾ n + (1/2) rank(B A^{-1} C - C A^{-1} B).   (4.2)

  ("We first show R(t) ⩾ ... (4.2)" — text-layer line 16 of p.29 with the
  sign read from the image; plain R, underline ABSENT — verified.)

- Why the theorem displays R̲ while the proof shows (4.2) with plain R —
  the printed proof structure, read line by line from the image:
  1. "We first show R(t) ⩾ n + ½ rank(BA⁻¹C − CA⁻¹B)  (4.2)": a direct,
     base-field-independent argument: W.1.o.g. A = I_n ("The tensor with
     slices I_n, A⁻¹B, A⁻¹C is isomorphic to t"); from a rank-r
     decomposition t = Σ u_ρ ⊗ v_ρ ⊗ w_ρ compare coefficients to get
     H'Γ_1Z = I_n, H'Γ_2Z = B, H'Γ_3Z = C (his (4.5)); augment Z to an
     r × r matrix; derive BC + BΘ₁₃ − CB − CΘ₁₃ = 0, hence
     rank(BC − CB) = rank(CΘ̄₁₃ − BΘ̄₁₃) ≤ 2(r − n) — his displayed (4.4).
     The vanishing-diagonal case is handled by induction on p (# of
     e_ρ1 = 0), each vanishing entry adding ≤ 2 via a rank-1 matrix:
     rank(BC − CB) ≤ rank(B'C' − C'B') + 2 ≤ 2(r−n).
     This step uses NO algebraic closure, NO Zariski closure: it is pure
     linear algebra on an existing rank-r decomposition, valid over any
     base field. Contrapositive: rank(BC − CB) = 2m ⇒ no decomposition
     with r < n + m, i.e. R(t) ≥ n + ½ rank(BC − CB).
  2. The final step passes to BORDER rank by Zariski closure/irreducibility
     ("Since the right side of this inclusion is closed, the left side may
     be replaced by its closure {t : R̲(t) ≤ r}"), giving the theorem's
     displayed R̲(t) form — needed by the paper because Theorem 4.6 is a
     statement about the border-rank hypersurface X_{(3n−1)/2}.
- Real-rank bookkeeping (rule 14, pre-registered at A3): R̲ ≤ R (paper's
  p. 648 line "Since R̲ ≤ R"), so the PLAIN-rank inequality (4.2) is the
  STRONGER statement at each individual tensor, and it is the one whose
  printed proof is field-independent. For OUR use on real tensors the
  load-bearing content is exactly (4.2) with the plain R; the underlined
  theorem-display form is weaker once moved to a real tensor but is the
  proper statement over the paper's algebraically closed k.

[Of these two glyph facts: the underlined R in the theorem display and the
plain R in (4.2) were each confirmed twice by independent image reads and
once by pixel-ink analysis. The (4.2)-hypothesis "no e_ij1 vanishes" case
split is in the printed proof.]

## COROLLARY 4.2 (page 676; PDF page 32)

Setup (text layer, clean): "Let A be a finite dimensional associative
algebra (with unity) and N an n-dimensional (unitary, left) A-module...
For a ∈ A let rank_A(a) denote the rank of the linear map l_a : N → N,
x ↦ a x."

[The map definition: l_a : N → N, x ↦ ax — text-layer shows "l_a : N → N,"
with the RHS in the display; image-verified.]

Statement (image-verified twice, `verification_images/cor42_v3.png` crop
and p32 full page):

"COROLLARY 4.2. Let N be a Λ-module, dim N = n, b, c ∈ Λ. Then

    R̲(N) ⩾ n + (1/2) rank_N(b c − c b)."

with the glyph facts:
- leading rank symbol R̲ = underlined (border rank) — underline under the
  R of R(N) only, confirmed by two independent reads (HIGH confidence in
  both), matching the theorem-display convention;
- the printed algebra letter is capital Lambda Λ ("two-stroke apex glyph
  with no crossbar"; also used in "1, b, c linearly independent over k"
  and in U = N*, V = N, W = (k·1 + kb + kc)*) — i.e. the corollary is
  stated for a finite-dimensional associative k-algebra Λ with k
  algebraically closed, acting on a left module N;
- the rank in the display carries a SUBSCRIPT: rank_N(b c − c b) — the
  subscript N (the module) is genuine ink below the baseline; two reads
  agree; the text-layer "rank,(a)" is the same subscript N garbled by OCR;
- terminal period after (bc − cb).

Proof chain inside the corollary (for the record): W.1.o.g. 1, b, c
linearly independent over k; sets U = N*, V = N, W = (k·1 + kb + kc)*;
dual map α : A* → W; the structural tensor t of N maps to t' whose slices
(over basis of V, dual basis for U, dual of (1, b, c) for W) are the
matrices of l_1, l_b, l_c (in particular A = 1_n); then
R̲(N) = R̲(t) ⩾ R̲(t') ⩾ n + ½ rank(BC − CB) = n + ½ rank_N(bc − cb)
(image-verified chain on p. 677 top: `p33_top.png`; all three R's
underlined; inequalities direction re-verified in second read).

## THEOREM 4.6 / Proposition 4.7 (recorded for completeness; NOT used)

Theorem 4.6 (p. 681): for n ≥ 3 odd, dim U = dim V = n, dim W = 3, the
polynomial F generates the ideal of the hypersurface X_{(3n−1)/2}; r =
(3n−1)/2. Proposition 4.7: n even ⇒ R̲(n,n,3) = 3n/2. These are the paper's
own 3-slice border-rank results. NOT load-bearing for this campaign; quoted
only for context of what the theorem machinery was built for.

## What is NOT in Strassen 1983 (negative search for other variants)

- Searched full text for every Theorem/Lemma/Corollary/Proposition header
  (24 headers enumerated): Section 4 contains ONLY Theorem 4.1, Corollary
  4.2, Corollary 4.3, Lemma 4.4, Lemma 4.5, Theorem 4.6, Proposition 4.7.
- NO third theorem gives a 3-slice bound of the form
  rank ≥ n + rank(commutator) (no factor 1/2). The only 3-slice rank
  lower bounds in the paper are Theorem 4.1 / Corollary 4.2 (both with
  the ½), plus the shape bounds (1.3). The "n + rank" reading previously
  entertained ("Strassen twin", retracted by Route AF) has NO primary-text
  support in Strassen 1983 — first-hand-confirmed negative here.
- The word "commutator" as a NAME for BC − CB does not appear in the
  paper (the expression is always written as a difference of products;
  "commutative diagram" is its only use of the root "commut...").
  The identification [X, Y] = BC − CB (when A⁻¹... ) is OURS to make and
  is trivial algebra, but the paper does not name it.
- Lickteig is mentioned (independent Jacobian-criterion application,
  announcement) — no Lickteig theorem in this paper.
- Blaeser is not cited anywhere in the paper (bibliography checked by
  scan: no 2000s entries; the paper is 1983).

## Real-rank standing of the import (pre-registered reading, applied at A3)

Strassen's k is algebraically closed (k = ℂ use). For a REAL tensor t_R:
every real rank-r decomposition complexifies to a k-rank-r decomposition,
so complex rank R_k(t) ≤ R_R(t); also border rank R̲_k(t) ≤ R_k(t).
Therefore R̲-form: R̲_k(t) ≥ n + ½ rank(...) ⇒ R_R(t) ≥ n + ½ rank(...).
The REAL-rank conclusion 18-window-relevant is the same number with the
plain-R proof (4.2) if its hypotheses extend to ℝ — the paper's proof of
(4.2) uses Zariski-closure/irreducibility ONLY in the final underlined
step, not in (4.2) itself; (4.2)'s coefficient-comparison argument is
purely linear-algebraic over the base field, so it is valid over ℝ as
well. THIS claim (validity of (4.2) over arbitrary base field incl. ℝ) is
a HUMAN-AUDITED reading of the printed proof: the proof of (4.2) uses
"R(t) ≤ r" decompositions, matrix rank inequalities, and induction on p =
#{vanishing diagonal entries}; none of it needs k algebraically closed or
any closure operation. Bookkeeping label: the (4.2)-over-ℝ step is
HUMAN-AUDITED (line-by-line reading of the printed proof); it is NOT
machine-verified (no formalization), and the campaign does not
over-reach: the lemma-level content that executes over ℝ is re-proven by
OUR OWN exact-arithmetic rank argument in the clause check (A2), so the
bound we report does not REST on the human audit alone — see clause check
C2/C3 below.
