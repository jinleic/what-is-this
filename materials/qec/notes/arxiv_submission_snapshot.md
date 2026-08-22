# arXiv submission snapshot — paper_pbb_nogo (pre-submission capture)

Captured 2026-08-18T22:03:39Z by the arXiv-kit agent, BEFORE any arXiv submission or
on-arXiv recompilation. Staging only — **do not post**. Full kit: `reports/arxiv_package.md`.
Machine-readable twin of this capture: `reports/arxiv_metadata.json`.

## Source fingerprint (what this snapshot describes)

| File | Bytes | MD5 |
|---|---|---|
| `reports/paper_pbb_nogo.md` | 33,262 | `4340160bb3e821359c0c8743a1ed6735` |
| `reports/paper_pbb_nogo.tex` | 37,483 | `ab90deef04a9522236184a0c0de278c5` |
| `reports/paper_pbb_nogo.bib` | 2,772 | `91ecec291e989e46ecf55f6acfd22f34` |
| `reports/paper_pbb_nogo.pdf` | 622,680 | `579d2d193bd05f90d68a7b646b33ab2d` |

Build log: `reports/paper_pbb_nogo.log (pdflatex, 2026-08-17 23:07)` → "Output written on paper_pbb_nogo.pdf
(10 pages, 619338 bytes)". Paper md mtime 2026-08-17 23:07.

## Metadata fields as staged

- **Title:** Perturbing a bivariate-bicycle code cannot buy distance for free
- **Subtitle (masthead, not an arXiv field):** An exact, perturbation-independent rate–distance obstruction for perturbed bivariate-bicycle (PBB) codes, with machine-checked certificates
- **Primary category:** quant-ph — cross-lists: cs.IT, math.CO
- **Comments:** 10 pages, 7 tables, 0 figures; two boxed main results; machine-checked with exact SAT certificates (CaDiCaL 1.9.5) and GF(2) replay; full provenance ledger and per-claim reproducers in the ancillary claims matrix
- **Measured counts:** 10 pages / 7 tables / 0 figures (from build log + tex grep, not the md)
- **License:** CC-BY 4.0 (recommended) vs arXiv default non-exclusive license
- **Authors/affiliations/contact/submitter:** [TO FILL BY HUMAN] — arXiv requires a named
  human submitter; `paper_pbb_nogo.tex` line 45 currently ships `\author{}` empty.

## Abstract verification (acceptance check)

Abstract = `## Abstract` body of `reports/paper_pbb_nogo.md`, lines 13-63
(3,158 chars, 51 lines, TWO boxed results). Checks run
2026-08-18T22:03:39Z:

- exact substring of the paper md: **PASS** (occurs exactly once)
- embedded copy in `reports/arxiv_package.md` §1 equals source bytes: **PASS**
- abstract string inside `reports/arxiv_metadata.json` equals source bytes: **PASS**

## Quickstart evidence (this session)

- `exp039_nogo_module.py gate`: **exit 0 in 1.5 s**, 7/7 certified reversals permitted,
  slack 0 on all, `theorem_refuted: false` (run 2026-08-18, single thread, GF(2) replay only).
- Timing ledger from stored certs: n≤72 → 3 s (11 certs); n≤144 → 37,242 s (~10.3 h);
  all 134 → 104,392 s (~29 h). The ≤30-min quickstart therefore uses `run --ns 72`, not
  the paper §8 `run --ns 144`.
- pytest steps listed in the kit are UNVERIFIED this session (SAT launches forbidden by
  session CPU policy); expected exit 0 per claims-matrix reproducers.

## Open blockers

1. Author block [TO FILL BY HUMAN] (hard blocker for submission).
2. License pick + submission are human actions.
