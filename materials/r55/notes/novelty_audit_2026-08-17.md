# Novelty Audit 2026-08-17 (primary-source verified)

Method: every claim verified by direct fetch of the primary source on
2026-08-17 (arXiv API/abs pages, journal DL pages, author PDFs, arXiv atom
feed, CrossRef). Full machine-readable table: `.superpowers/sdd/novelty-audit`
agent record; citations below are verbatim with quote and URL.

## Verified status of the field

- **R(5,5) ∈ {43, 46}** — DS1 revision #18 (2026-04-24), Table Ia/B: lower 43
  (Exoo, JGT 13, 1989, 97-98, key [Ex5]); upper 46 (Angeltveit–McKay [AnM3],
  now also journal-published: J. Graph Theory 11x (2026) 1-11, arXiv:2409.15709v2).
- **AM46 method**: m=2-excess LP ("the m = 2 case of [8, Theorem 2.2]") +
  gluing; NO m≥3 rows, NO proof certificates; trust = dual independent
  implementation; cost ~30 + ~50 CPU-years.
- **MR97** (r55.pdf): Thm 2.2 general display verified; (I2)/(I3) verified;
  Thm 3.2's n=49 kill via m=4 verified verbatim (right side 132 vs
  s(K4, H̄₁)=144, s(K4, H̄₂)=138 — the figure our n=49 replay reproduces as
  1584=12·132 with i4∈{138,144}); §5 LP(s,t,n) (pp.11-14) with degree and
  per-vertex (order,edges) bins, exact rationalization protocol, proves
  R(4,6) ≤ 41.
- **Engström arXiv:1002.4304v3** (terminal): Thm 2.18 degree-1 normal form;
  proves MR97 Conjecture 1 for K=4; his Conjectures 1.1/1.2 remain open.
- **Lidicky–Pfender 2021** (SIAM JDM 35:2328-2344) and **SIAM Review 68(2):
  385-403 (2026)**: no R(5,5) movement; flag/SDP lane only.
- **arXiv sweep**: ti:"Ramsey", 2025-01-01..2026-08-18, 263 entries read in
  full. The ONLY abstract printing "R(5,5)" is Tamburini 2508.16699v2 —
  heuristic quantum diagnostics, not a proof. No exact-certificate or
  neighborhood-state LP cone work. Certificate-lane precedent elsewhere:
  Gauthier–Brown (HOL4), MathCheck/DRAT for R(3,8)/R(3,9) (IJCAI 2025).

## What remains unclaimed (campaign delta, stands after audit)

1. Exact motif-refined per-vertex LP cones (q-refined states) with exact
   rational affine certificates — closest precedent MR97 §5.
2. Fully disjoint independent checker reconstructing the entire evidence
   chain from hash-pinned catalogs (AM trust model is dual-implementation,
   not certificate-verification).
3. Any identity-row (m=2/3/4) certified restriction on hypothetical
   Ramsey(5,5,45) graphs — zero rival publications attempt such restrictions.

## Scoop watch

- Braginsky–Yellenki–Zhu et al. SAT-attacking srg(45,22,10,11) (existence
  lane; would prove R(5,5)=46 the "hard" way) — orthogonal.
- Angeltveit arXiv:2602.11459 (R(K5,K5-e)=30) — machinery off-diagonal, not
  (5,5).
- Tamburini heuristic — preempt "already shown?" pushback with one footnote.

Scope honesty: arXiv-title sweep + DS1 backstop; a fresh non-arXiv claim
between 2026-04-24 and 2026-08-17 is possible but unlikely; SIAM Review full
text unfetched (403), assessed via CrossRef + 69 references.
