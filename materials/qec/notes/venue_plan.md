# Venue plan — PBB no-go paper (2026-08-19)

**SUBMISSION HELD per user decision (2026-08-19).** Do NOT post to arXiv; author block
remains empty; no external publication action is authorized. Internal artifacts
(package, metadata, snapshot, this plan) stay current as the paper evolves.

Status of the artifact: `reports/paper_pbb_nogo.md` + `.tex` + `.pdf` (11 pp, 898/899
suite green, md↔tex synced, hostile-review rounds 1–2 closed with NO_CRITICAL).
Author block is user-HOLD (see `reports/arxiv_metadata.json`, `author_placeholders`) —
resume the route below only on explicit user instruction.

## Route (recommended)

1. **arXiv preprint (immediate)**. Categories: **quant-ph** (primary), cross-list
   **cs.IT** (secondary). License CC-BY 4.0 per `arxiv_metadata.json`. Ancillary
   ancillary files listed there; `paper_pbb_nogo.bib` ships in-tree.
   Rationale for posting before the matched-MC (W6) verdict lands: the paper's
   theorem layer is complete and self-contained; the circuit section is explicitly
   scoped "not settled", so the preprint claims nothing the running experiment could
   contradict. If W6 lands later it becomes v2 material.
2. **QIP 2027 poster/talk submission** (deadline historically **mid-October**, four-page
   extended abstract; non-archival, compatible with journal submission). The exact
   obstruction theorem + machine-checked certificates fit QIP's profile; public
   deadline check needed when the CFP posts.
3. **Journal submission within ~4–6 weeks**, ideally including the W6 matched-MC
   verdict if the cells complete:
   - **Quantum** (first choice): open-access, explicit support for negative/results
     with reproducible-machine evidence; theory paper profile matches (proofs +
     certificates + explicit scope). No page fee issue beyond standard OJ fees.
   - Alternatives: **PRX Quantum** (if the circuit narrative lands strong),
     **TQC** proceedings talk (theory-fit, deadline-dependent), **IEEE Transactions
     on Quantum Engineering** (engineering-leaning, weaker fit), **npj Quantum
     Information** (general audience, weaker fit for a proof-centric negative).

## Non-blocking open items that would strengthen v2 / journal version

- W5/W6 of EXP-040: production-decoder CP intervals + final separation verdict
  (matched MC running, `exp040-mc`; stage-persistent).
- n=180/360 parent certification sweep (60 rows pending, background shards).
- Conjecture B′ residual-class extension beyond ℓ·m ≤ 40 (budget wall at 42+).
- J.5/J.7 (unconditional d_X(Q) ≥ d_X(P)) — currently OPEN-typed in the paper;
  EXP-044 evidence layer is in (0 violations, 6.47M perturbations).

## What must NOT ship uncorrected

- The `[5,12]` circuit-distance claim is **structural-probability DEM only**; any
  v2 claiming more needs the probability-informed certificate (out of scope today).
- The figure PDFs (`fig_rate_distance.pdf`, `fig_trade_law.pdf`) are generated and
  referenced in `notes/figure_provenance.md`; if included in the submission they
  need \includegraphics wiring (currently text-only paper).
