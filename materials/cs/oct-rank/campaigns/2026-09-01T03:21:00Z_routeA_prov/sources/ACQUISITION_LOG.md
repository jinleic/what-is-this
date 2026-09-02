# Sources directory — contents and acquisition log

Campaign: Route A provenance (pre-statement committed before this log is
complete; the hunt below starts AFTER frozen scope).

Goal: FIRST-HAND text of the exact hypotheses/conclusions of the 3-slice
(Strassen-type) lower-bound theorems relevant to

        rank(T) >= n + (1/2) * rank([A1^-1 A2, A1^-1 A3])

with full citation (author, title, venue, year, theorem number), plus the
n + rank-flavored and Blaeser variants if they exist as separate theorems.

## Known bibliographic facts (from Route AF campaign source_audit.txt)

- V. Strassen, "Rank and optimal computation of generic tensors",
  Linear Algebra and its Applications 52/53 (1983), 645-685,
  DOI 10.1016/0024-3795(83)80041-X.
- J. M. Landsberg, "Geometry and the Complexity of Matrix Multiplication"
  (survey), Theorem 6.1.1 attributes a 3-slice commutator inequality to
  [52] = Strassen 1983; exact transcription already frozen in
  campaigns/2026-08-31T08:02:18Z_routeAF/source_audit.txt.
- P. Koiran, "On tensor rank and commuting matrices", arXiv:2006.02374v2,
  Theorem 1 (page 1) and Lemma 6 (page 9) — secondary restatement with the
  1/2 factor, read first-hand by Route AF.
- M. Blaeser, "On the complexity of the multiplication of matrices of
  small formats", Theor. Comput. Sci. 306 (2003) — status per prior
  campaign: no primary found supporting a separate 3-slice formula;
  must be established or left UNVERIFIED.

## Access attempts (log, retractions inline)

1. ScienceDirect direct PDF (`.../pdfft?isDTMRedir=true&download=true`):
   HTTP 403 (2026-09-01T03:2xZ). Same as Route AF's attempts.
2. OpenAlex API: marks the article BRONZE open-access with pdf_url the
   ScienceDirect `/pdf` endpoint; that endpoint 403s (Cloudflare-class
   challenge; the body fetched was an HTML challenge document).
3. Semantic Scholar API: openAccessPdf = same ScienceDirect URL (bronze).
   No repository copy.
4. EuDML: exact-phrase title search returns ZERO documents (content-layer
   read). Evidence: evidence_eudml_search_negative.{html,md}.
5. Wayback Machine CDX for the article URL: only landing pages (200 in
   2022-2024, 403 from 2025 on) and ONE /pdf capture 20240423085041 which
   is itself a bot-challenge HTML page, not a PDF. Domain-wide CDX sweep
   for PII-bearing captures running at time of writing.
6. zenodo / university archives / institutional copies: pending.

Verdict so far: primary text NOT yet obtained. Continue until every
reasonable venue is exhausted, then freeze the provenance blocker per the
assignment ("If you cannot obtain the primary text, STOP and freeze that
as the blocker").
