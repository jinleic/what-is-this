# EuDML provenance probe — content-layer negative

Provenance checkpoint (owner note, 2026-09-01T03:21Z session):

- Query sent, exactly: `"Rank and optimal computation of generic tensors"`
  (percent-encoded) to BOTH:
  - `https://eudml.org/search/page?q=...` (direct)
  - `https://eudml.org/searchForm?search-box-1=...&SCHEME=general&expression-1=EQUALS_all`
    (followed its 302 to the canonical `search/page` URL
    `q=sc.general*op*l_0*c_0all_0eq%3A1.%22Rank+and+optimal+computation+of+generic+tensors%22&qt=SEARCH`)
- Result layer read: DOCUMENT CONTENT. The `#primary-content` block of the
  returned page says verbatim:
  `<p>Your search did not match any documents.<BR>Suggestions:...`
- Interpretation: EuDML does not host the digitized Strassen 1983 LAA 52/53
  article under this title. This rules out EuDML as a full-text source.
- The raw HTML of the deciding page is retained as
  `evidence_eudml_search_negative.html` (annotated in place with the note).
