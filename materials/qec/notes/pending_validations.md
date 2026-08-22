# Pending validations — arXiv kit (created 2026-08-18T22:03:39Z)

Created by the arXiv-kit follow-up (file did not previously exist). Items the kit
leaves open, each mapped to its satisfying test/artifact. Paths relative to `math/qec/`.

| # | Item | Status | Blocking on | Satisfying test/artifact |
|---|---|---|---|---|
| 1 | Author/affiliation/submitter block | OPEN (hard blocker) | human | arXiv named-human submitter policy; `reports/paper_pbb_nogo.tex` line 45 `\author{}` must be filled |
| 2 | License choice (recommend CC-BY 4.0) | OPEN | human decision | `reports/arxiv_package.md` §7, `reports/arxiv_metadata.json`.license |
| 3 | Quickstart step 3 runtime (`pytest tests/test_pbb_theorems.py -q`, exp 0, 4 passed) | UNVERIFIED this session | session CPU policy forbids SAT launches | `reports/claims_matrix.md` C51/C52 reproducers |
| 4 | Quickstart step 4 runtime (`pytest tests/test_pbb_nogo.py tests/test_pbb_survival.py -q`, exp 0, 17+15) | UNVERIFIED this session | same policy | `reports/claims_matrix.md` C01-C08 rows |
| 5 | Hostile-review disposition recorded from session context (NO_CRITICAL_FINDING, LOW-1 fixed) | [INFERENCE] — review artifact not vendored under reports/ | vendor the review report if it must be citable | session context; `reports/claims_matrix.md` ledger |
| 6 | Ancillary exp039 artifact path | RESOLVED as discrepancy | none — documented | template said `results/processed/exp039_nogo_module.json` (absent); real path `results/partial_runs/exp039_nogo_module.json` (66,430 B, paper §8) |

Verified-closed this session (for contrast): abstract byte-verbatim substring check
PASS; PDF counts measured from build log (10 pages, 7 tables, 0 figures); gate re-run
exit 0 in 1.5 s. See `notes/arxiv_submission_snapshot.md`.
