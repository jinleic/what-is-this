# Open-status record — freshness check 2026-08-17

**Scope.** Novelty/open status of the Theorem G/H headline claims and two collateral
claims, refreshed 2026-08-17. "Open status supported" below means the most recent
external literature state we could verify does NOT contain the claim; it is not a
claim of exhaustive world-literature coverage. Supersedes nothing in
`notes/open_status.md`; this file is the dated novelty annex.

**Method audit trail (all executed this session, 2026-08-17).**

Web searches (Google-style):
1. `"perturbed bivariate bicycle" quantum code` -> 2 hits, both the catalogue's own
   qcode-discovery repo. NEGATIVE.
2. `bivariate bicycle code perturbation distance increase quantum` -> BB-generic
   results (Postema-Kokkelmans, ECZoo, Bravyi GitHub); no PBB theory. NEGATIVE.
3. `"Gross code" 144 12 12 distance quantum error correction` -> ECZoo gross page,
   architectures, Bravyi et al.; no perturbation distance-increase theory. NEGATIVE.
4. `arxiv 2606.02418 perturbed bivariate bicycle catalogue` -> self, repo, author's
   page only. NEGATIVE.
5. `citing "2308.07915" bivariate bicycle 2026 perturbation non-CSS` -> catalogue +
   surveys only. NEGATIVE.
6. `"non-CSS" "144,12,12" quantum code` -> mirror codes 2603.05496, catalogue,
   unrelated IEEE. NEGATIVE for PBB theory.
7. `"2606.02418" citing OR citation perturbation theorem` -> self/repo/aggregators;
   Google Scholar shows "Cited by 3". NEGATIVE.
8. `"perturbed bivariate bicycle" theorem dimension k truncation no-go` -> nothing.
9. `bivariate bicycle "dimension formula" OR "logical dimension" perturbation kernel
   symplectic 2026` -> Gröbner k for CSS BB (2502.17052) only. NEGATIVE.
10. `"PBB code" quantum` -> nothing relevant. NEGATIVE.
11. `bivariate bicycle "envelope" distance rate upper bound quantum LDPC` -> only the
    catalogue's own rate-distance-envelope remark. NEGATIVE.
12. `arxiv 2026 non-CSS quantum code bicycle "logical dimension" theorem stabilizer
    rank` -> self-dual BB 2510.05211, SBB 2605.04151, qudit twisted-torus 2602.04443,
    mirror 2603.05496; none carries a PBB dimension identity. NEGATIVE.
13. `site:arxiv.org bivariate bicycle perturbation dimension` -> catalogue self +
    2607.05724 (QCNN) only. NEGATIVE.
14. `"perturbed bicycle" OR "dressing space" quantum code dimension` -> zero hits.
15. `no-go theorem quantum LDPC perturbation "logical qubits" parent dimension bound`
    -> only unrelated transversal-gadget no-go (2602.13395). NEGATIVE.

arXiv API (export.arxiv.org, 2026-08-17):
16. `all:"perturbed bivariate bicycle"` -> exactly 1 record: arXiv:2606.02418 itself.
    NEGATIVE.
17. `all:"bivariate bicycle"` sortBy=submittedDate desc, top 40 -> 10 post-June-2026
    papers; every abstract read: 2608.15754 (RL decoder co-design), 2608.11516
    (trapping sets), 2608.09913 (Koszul/SET, CSS-only), 2608.09115 (univariate
    bicycle search), 2608.02773 (Cornucopia codes), 2607.19563 (syndrome
    post-selection), 2607.06177 (Tanner routing), 2607.05724 (QCNN), 2607.04462
    (encoder synthesis), 2606.27119 (neural decoders). NONE on PBB theory. NEGATIVE.

Citation listings:
18. INSPIRE `refersto:recid:3164530` -> 1 citing record: Heitritter, Brown, Hardikar,
    "Evolving Quantum Error-Correcting Encodings for Molecular Simulation"
    (recid 3172672, arXiv-created 2026-06-25); methodological citation of the
    catalogue's discovery workflow; no PBB theory. NEGATIVE.
19. OpenAlex W7163132656 -> cited_by_count 0; preprint-only record. NEGATIVE.
20. Semantic Scholar API -> UNREACHABLE via available tooling (endpoint parse-blocked);
    compensated by (17)+(18)+(19)+Google Scholar "Cited by 3".

Revision metadata and full-text checks:
21. arXiv:2606.02418: submission history shows v1 (2026-06-01) ONLY; no v2, no
    errata, no journal version (OpenAlex: not published, not retracted).
22. arXiv:2502.17052: v1 2025-02-24, v2 2025-02-25, v3 2025-05-01, v4 2026-04-30
    (latest). Full-text of v4 contains ZERO occurrences of "perturb" -> says nothing
    about perturbation k-truncation; its dimension results are CSS-BB Gröbner theory.
23. Full-text greps: 2603.05496v1 (mirror codes) -> no "perturb", no "no-go", no
    closed-form k beyond the definition k := log dim C(S); 2510.05211v2 -> no
    "perturb"; 2605.04151v1 -> no "perturb".

## Status table

| Claim | Exact statement | Most recent source supporting open status | Checked | Confidence | Residual risk |
|---|---|---|---|---|---|
| (a) Perturbation-independent no-go / closure theorem | For PBB child Q over CSS BB parent P: d_Q > d_Z(P) => k_Q <= k_P - T(P), T(P) an exact parent-only invariant; family-closed when T(P) = k_P (holds for 61/202 catalogue parents, incl. Gross T = 12 = k_P and all 11 catalogue [[144,12,12]] parents) | Catalogue arXiv:2606.02418v1 (2026-06-01, only version) contains only an A=B BB distance-trap theorem and empirical PBB distances; entire post-June-2026 BB arXiv batch (through 2608.15754, 2026-08-16) has no such theorem; searches 1,8,13-17 all negative | 2026-08-17 | HIGH — NOVEL_LIKELY | Differently-phrased equivalent (e.g. in mirror-code/2BGA group-algebra language) outside arXiv phrase coverage; Semantic Scholar citation graph unreachable this session; journal-latency papers not yet indexed |
| (b) Exact dimension identity | k_Q = k_P - dim(Delta_bar), where Delta = {lambda[C;D] : lambda[A;B] = 0} is the dressing R-submodule, R = F2[x,y]/(x^l-1, y^m-1) | Nearest external dimension statements: Postema-Kokkelmans arXiv:2502.17052v4 (2026-04-30; Gröbner k for CSS BB, zero perturbation content in v4) and mirror codes arXiv:2603.05496v1 (2026-03; definition-level k, no closed form, no perturbation framing); searches 2,9,12,14,16,22,23 negative | 2026-08-17 | HIGH — NOVEL_LIKELY | A 2BGA/group-algebra paper could state the same rank computation without calling it an identity; non-arXiv classical-coding venues under-covered; S2 unreachable |
| (c) Saturation observation | All 7 certified distance-increasing PBB perturbations satisfy k_Q = k_P - T(P) with slack exactly 0; each CSS-dominated at equal n | Internal: results/processed/exp036_envelope_check.json (all_reversals_css_dominated_at_equal_n: true). External: catalogue v1 (2026-06-01) reports reversals as isolated distance values with no parent-level invariant; nothing published since (searches 5-8,16-19 negative) | 2026-08-17 | HIGH — NOVEL_LIKELY as an external observation (its truth rests on our exact-T certificates, tests/test_pbb_nogo.py) | Same channels as (a); observation could be independently re-derived quickly once (a) is public |
| (d) CSS envelope stands for delta>0 PBB catalogue | Every delta>0 PBB catalogue row is dominated at equal n by a certified-exact CSS code (249/368 rows parent-capped without search; all 7 reversals CSS-dominated; no PBB child beats the certified CSS envelope) | Internal: results/processed/exp036_envelope_check.json + EXP-037 certified CSS pool. External: catalogue arXiv:2606.02418v1 (2026-06-01) is the most recent systematic PBB distance source and reports no envelope-beating child; post-June BB batch clean; searches 3,6,11 negative | 2026-08-17 | HIGH (as a computation); MEDIUM-HIGH (that no published non-CSS beater exists) | Catalogue's own upper-bound rows are untrusted distances — a hidden beater could sit in unclosed n in {180,360} MILP rows; non-arXiv code tables |
| (e) Circuit-level end-to-end PBB [[144,12,12]]: open | No published end-to-end circuit-level fault-tolerant memory demonstration establishing any advantage of a PBB [[144,12,12]] over the CSS Gross; internal: EXP-024 NEGATIVE_STRUCTURAL (multi-ancilla gadget route), EXP-029 interval 4 <= d_DEM_mech <= 12, EXP-030 Gross-vs-surface mapping | Internal ledger rows checked 2026-08-13 (EXP-024/029/030 artifacts); external freshest: mirror codes arXiv:2603.05496v1 (2026-03) runs circuit-level experiments on its own non-CSS codes (pseudothreshold ~0.2%) but not on catalogue PBB; Tour de Gross arXiv:2506.03094v1 (2025-06-03) for Gross FT context | 2026-08-17 | HIGH that it remains open | Mirror-code group could land a [[144,12,12]]-class non-CSS circuit-level demo; internal EXP-025/EXP-026 still in flight and could change our own end-to-end verdict |

## Watch list (papers that come close)

1. **arXiv:2603.05496 v1 — mirror codes (2026-03).** Nearest neighbor overall.
   Non-CSS stabilizer codes from (G,A,B) with mixed Z(A g) X(B g^{-1}) checks;
   contain BB codes up to permutation/local Clifford. Has a well-definedness kernel
   characterization (Prop. 3.2) but no dimension identity, no distance no-go, no
   parent-perturbation framing. Risk channel: a v2 extending the group-algebra
   formalism to k formulas. Monitored: v1 is current as of 2026-08-17.
2. **arXiv:2608.09913 — Koszul-complex SET rigidity (Song, 2026-08-10).**
   Homological machinery computing BB k via Ext/Tor; CSS translation-invariant only.
   One conceptual step (a short exact sequence for the perturbed check complex) from
   a (b)-type identity; does not take that step.
3. **arXiv:2606.02418 + qcode-discovery repo.** Any v2 of the catalogue could add
   theory; v1 (2026-06-01) is the only version as of 2026-08-17. Repo README already
   states the A=B distance trap; no closure theory.
4. **arXiv:2502.17052 v4 (2026-04-30), Postema-Kokkelmans.** Definitive CSS-BB
   dimension/existence theory; v4 still contains no perturbation content (grep-verified).
   A v5 could in principle extend to PBB.
5. **arXiv:2608.09115 — divisor-driven univariate bicycle search (2026-08-10).**
   Exact distance certification via the F4-additive correspondence on the cyclic
   (univariate) family; signals fast.algebraization of bicycle-type families.
6. **arXiv:2510.05211 v2 — self-dual BB with transversal Cliffords (npj QI 2026).**
   BB k theory on twisted tori; no perturbation content (grep-verified).
7. **Heitritter-Brown-Hardikar GSE (arXiv-created 2026-06-25; INSPIRE 3172672).**
   Sole independently tracked citation of the catalogue; methodological only. Monitor
   for a follow-up in which a discovery team bolts perturbation theory onto the catalogue.

## Verdicts

- (a) Perturbation-independent no-go/closure theorem: **NOVEL_LIKELY**.
- (b) Exact dimension identity k_Q = k_P - dim(Delta_bar): **NOVEL_LIKELY**.
- (c) Saturation k_Q = k_P - T(P) on certified reversals: **NOVEL_LIKELY**.

All three: at least four independent search channels negative each (phrase-level arXiv
API, Google-style web, citation listings INSPIRE/OpenAlex, post-June arXiv BB batch).
Named residual risk for all three: Semantic Scholar citation graph unreachable this
session; non-arXiv venues and differently-phrased formulations under-covered.
