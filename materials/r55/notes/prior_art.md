# Prior Art: Computer-Assisted / Certified Ramsey-Number Computation

Recon notes for the R(5,5) <= 45 campaign. Compiled 2026-08-13.

Evidence discipline used in this file:

- Every claim carries a source URL.
- **[VERBATIM-VERIFIED]** = quote checked by me directly against raw page/PDF text
  (curl + local extraction), not just against a search snippet.
- **[FETCHED-ABSTRACT]** = quote returned by an automated fetch of the primary page
  (arXiv abstract page). Low risk, but re-check before citing in a paper.
- **[REPORTED]** = secondary source or search-result claim not verified against a
  primary source; treat as a lead, not a fact.
- I distinguish what a source PROVES (theorem + artifact) from what it merely reports.

---

## 0. Current status of R(5,5) (baseline for the campaign)

**Best published bounds: 43 <= R(5,5) <= 46.**

- Lower bound 43: Exoo (1989). Reported in Angeltveit-McKay 2024: "The lower bound of
  43, which was established by Exoo [5] in 1989, is still the best." [VERBATIM-VERIFIED
  from https://arxiv.org/html/2409.15709v2]
- Upper bound 46: Angeltveit & McKay, arXiv:2409.15709 (v1 2024-09-24, v2 2025-09-01).
  https://arxiv.org/abs/2409.15709
- Widely-believed value: 43. Gasarch (blog, 2024-09-29): "It is widely believed that
  R(5)=43." [REPORTED — commentary]
  https://blog.computationalcomplexity.org/2024/09/progress-on-r5-will-there-be-more.html
- DeepMind `formal-conjectures` issue #2364 (opened 2026-02-19) requests a Lean
  statement of R(5,5) bounds; it still quotes the stale bounds "43 ≤ R(5,5) ≤ 48",
  confirming no formal statement/verification of the 46 bound exists there yet.
  [VERBATIM-VERIFIED via GitHub API]
  https://github.com/google-deepmind/formal-conjectures/issues/2364

**No published or announced proof of R(5,5) <= 45 was found** (searches across arXiv,
Google, blogs, sat4math.com index, 2024-2026). Closest live activity is the CMU
srg(45,22,10,11) SAT search (Section 3.1), which attacks the n=45 case from the
*existence* side (would prove R(5,5)=46, i.e., would make <= 45 false).

### 0.1 Angeltveit-McKay R(5,5) <= 46: proof structure and compute (the thing we would reproduce/beat)

Source: https://arxiv.org/abs/2409.15709 and https://arxiv.org/html/2409.15709v2
(all quotes below [VERBATIM-VERIFIED] from the v2 HTML).

- Abstract (v2 body): "All of the computational parts of the proof were independently
  implemented by both authors, with consistent results." (Abstract page wording:
  "All of the computations were independently implemented by both authors, with
  consistent results." [FETCHED-ABSTRACT])
- Skeleton: any F in R(5,5,46) has all degrees in {21,22,23,24}; each vertex
  neighbourhood is a (4,5)-Ramsey graph and each non-neighbourhood a (5,4)-Ramsey
  graph; LP ("excess" counting) reduces which parts of the (4,5) catalogs are needed;
  the rest is pointed-graph gluing.
- Catalog facts: "|R(4,5,24)| = 352,366" (census completed in the 2017 <=48 paper,
  started by McKay-Radziszowski); estimates "|R(4,5,21)| ≈ 5.5×10^17",
  "|R(4,5,22)| ≈ 1.9×10^15", "|R(4,5,23)| ≈ 10^11" — "Hence any approach that relies
  on the complete catalogue of these Ramsey graphs is impractical. Instead we use
  linear programming to reduce the number of graphs we need to consider."
- R(4,5,23,e>=119) census alone: "Finding all of these graphs was the second most time
  consuming part of the project, taking approximately 5 years of CPU time."
- Scale: "in total we computed approximately 2 trillion gluing operations."
- **SAT is already inside this proof.** First author: "Each gluing problem can be
  encoded as a SAT problem whose clauses forbid cliques or independent sets of size 5.
  We also added symmetry breaking clauses to distinguish the vertices of L1." and
  "Our SAT solver was a very simple special purpose solver without advanced features
  like clause learning and restarts." Second author: "the reduced SAT program was
  handed to the SAT solver Glucose [3] with a time limit of 1 minute. Glucose decided
  satisfiability in 1 second on average." plus "As a sanity check, a number of cases
  were run as well using Kissat [4] instead of Glucose."
  (Gasarch's blog also says they "used Glucose, a SAT Solver" [REPORTED]; the paper
  confirms Glucose only for the second author's gluing pipeline.)
- Compute: "We estimate that completing the census of R(4,5,n,e≥e0) took approximately
  15 years of CPU time while gluing the necessary graph took another 15 years of CPU
  time. Hence the whole project took about 30 years of CPU time for the first author to
  complete." … "This replication took about 50 years of additional CPU time."
- Trust model: dual independent implementation, **no proof certificates**: "the two
  approaches were deliberately intended to be disjoint to avoid the well-known axiom of
  so[ftware]…" (sentence continues in paper); "Given the theory and the two independent
  computations, Theorem 1.1 has now been firmly established."
- Their outlook on 45: "We believe that future improvements to the upper bound on
  R(5,5) will need new theoretical insights, as an excessive amount of computer time
  would be required to apply the same method."
- Gasarch's cost guess: "Getting R(5)≤45 may take 3000 years of CPU time. So it may not
  be for a while or ever." [REPORTED — blog commentary via automated fetch]
  https://blog.computationalcomplexity.org/2024/09/progress-on-r5-will-there-be-more.html
- Data artifact: Angeltveit & McKay, "Edge-extremal Ramsey (4,5)-graphs",
  https://users.cecs.anu.edu.au/~bdm/data/ramsey.html (cited as [2] in the paper) —
  relevant to Gate 1 (independent validation of the R(4,5,24) catalog).
- Predecessor: Angeltveit & McKay, "R(5,5) <= 48", arXiv:1703.08768 (2017), J. Graph
  Theory 89(1):5-13, 2018. Abstract: "We improve the upper bound on the Ramsey number
  R(5,5) from R(5,5) ≤ 49 to R(5,5) ≤ 48. We also complete the catalogue of extremal
  graphs for R(4,5)." [FETCHED-ABSTRACT] https://arxiv.org/abs/1703.08768

**Campaign takeaway:** the 46 proof is LP + catalog + ~10^12 SAT-shaped gluing subproblems,
verified only by dual implementation. A certified SAT reproduction (or a 45 push) is
exactly the gap the certificate tooling in Section 4 was built to fill — but the
LP/counting layer and the catalog-completeness layer are *outside* plain DRAT and need
either VeriPB-style PB reasoning or theorem-prover glue.

---

## 1. SAT-based Ramsey results

### 1.1 Codish-Frank-Itzhakov-Miller: R(4,3,3) = 30

- arXiv:1510.08266, "Computing the Ramsey Number R(4,3,3) using Abstraction and
  Symmetry breaking" (2015); journal version Constraints 21(3), 2016,
  DOI 10.1007/s10601-016-9240-3.
  https://arxiv.org/abs/1510.08266 / https://link.springer.com/article/10.1007/s10601-016-9240-3
- Abstract claims: "compute the value R(4,3,3)=30"; intermediate step "the previously
  unknown set R(3,3,3;13) consisting of 78,892 Ramsey colorings." [FETCHED-ABSTRACT]
- Companion: arXiv:1409.5189, "Breaking Symmetries in Graph Coloring Problems with
  Degree Matrices: the Ramsey Number R(4,3,3)=30". https://arxiv.org/abs/1409.5189v2
- PROVES (computationally): R(4,3,3)=30 via SAT + degree-matrix abstraction + symmetry
  breaking. REPORTS: no independently checkable proof certificate (pre-dates practical
  certified symmetry breaking) — trust rests on the toolchain. [REPORTED — I did not
  find a certificate artifact for it; verify if it matters.]
- Relevance: the classic template for "abstraction (degree sequences) + symmetry
  breaking + SAT" on multicolor Ramsey; the abstraction layer is the same kind of
  uncertified gap our campaign must close with VeriPB/Lean.

### 1.2 Heule methodology: cube-and-conquer + DRAT at extreme scale

- Cube-and-conquer origin: Heule, Kullmann, Wieringa, Biere, "Cube and Conquer:
  Guiding CDCL SAT Solvers by Lookaheads", HVC 2011, DOI 10.1007/978-3-642-34188-5_8.
  [VERBATIM-VERIFIED citation from SMS survey reference list]
- Boolean Pythagorean Triples: Heule, Kullmann, Marek, "Solving and Verifying the
  boolean Pythagorean Triples problem via Cube-and-Conquer", arXiv:1605.00723 (SAT
  2016). https://arxiv.org/abs/1605.00723 — produced "a proof in the DRAT format,
  which is almost 200 terabytes in size", compressed certificate 68 GB; solved on 800
  cores in ~2 days. [REPORTED — figures from search summary of the paper; primary
  paper is the citation]
- Schur Number Five: Heule, AAAI 2018, arXiv:1711.08076. https://arxiv.org/pdf/1711.08076
  — the "largest proof ever" style methodology with certified result. [REPORTED]
- Relevance: this is the canonical pipeline for our conquer phase: cube-and-conquer +
  per-cube DRAT/LRAT + a *cube-cover completeness* certificate (see LRAT-Catcher,
  Section 4.5, which packages exactly this for Lean).

### 1.3 Neiman-Mackey-Heule: directed Ramsey / tournaments

- Neiman, Mackey, Heule, "Tighter Bounds on Directed Ramsey Number R(7)",
  arXiv:2011.00683; Graphs and Combinatorics 38, #81 (2022),
  DOI 10.1007/s00373-022-02560-5. https://arxiv.org/abs/2011.00683
- Result: 34 <= R(7) <= 47 (transitive subtournaments); includes computer-assisted
  proof of a Sanchez-Flores conjecture (all TT6-free tournaments on 24-25 vertices are
  subtournaments of ST27). [REPORTED — from search summary of abstract]
- Relevance: Heule-style SAT applied to a Ramsey-type extremal problem with
  catalog+extension ("glue-like") structure; closest Heule-group analog of the
  Angeltveit-McKay workflow.

### 1.4 MathCheck (SAT+CAS): verified certificates for R(3,8) and R(3,9)

- Li, Duggan, Bright, Ganesh, "Verified Certificates via SAT and Computer Algebra
  Systems for the Ramsey R(3,8) and R(3,9) Problems", arXiv:2502.06055; IJCAI 2025
  (https://www.ijcai.org/proceedings/2025/0292.pdf).
  https://arxiv.org/abs/2502.06055
- Abstract [FETCHED-ABSTRACT]: "we use the software MathCheck to generate certificates
  for Ramsey problems R(3,8) and R(3,9) … Our SAT+CAS approach significantly
  outperforms traditional SAT-only methods … solves R(3,8) … sequentially in 59 hours
  … while a SAT-only approach using state-of-the-art CaDiCaL solver times out after 7
  days. … Our results provide the first independently verifiable certificates for
  these Ramsey numbers, ensuring both correctness and completeness of the exhaustive
  search process."
- PROVES: R(3,8)=28 / R(3,9)=36 recomputations *with* independently checkable
  certificates (CAS isomorph-free generation logged so the completeness of the split
  is auditable). Uses parallel cube-and-conquer for R(3,9).
- Relevance: the closest existing "certified Ramsey number" pipeline; the paper to
  benchmark our certificate design against. Note it is at R(3,k) scale — orders of
  magnitude below R(5,5) gluing.

### 1.5 SAT Modulo Symmetries (Kirchweger-Szeider) and Ramsey applications

- Original: Kirchweger & Szeider, "SAT Modulo Symmetries for Graph Generation",
  CP 2021; extended: ACM TOCL 2024, DOI 10.1145/3670405.
  https://dl.acm.org/doi/10.1145/3670405
- Co-Certificate Learning: Kirchweger, Peitl, Szeider, IJCAI 2023,
  DOI 10.24963/IJCAI.2023/216. https://www.ijcai.org/proceedings/2023/0216.pdf
- Survey: Szeider, "SAT Modulo Symmetries: A Survey", CEUR-WS Vol-4116.
  https://ceur-ws.org/Vol-4116/invited1.pdf
  Quotes [VERBATIM-VERIFIED from the PDF]:
  - Ramsey enumeration: "Using a direct encoding of these properties, SMS enumerated
    the complete sets for R(3,5,n) and R(4,4,n). The framework's ability to produce
    verifiable proofs of completeness is particularly valuable in this context. The
    DRAT proof from the SAT solver certifies the search, while nc-certificates for
    each learned symmetry clause can be independently verified [4]."
  - Proof architecture: "For the symmetry-breaking component, SMS outputs an
    nc-certificate (a permutation π) for each learned symmetry-breaking clause. An
    independent checker can verify in polynomial time that this permutation indeed
    witnesses the non-canonicity of the pruned partial graph."
  - Two-pass DRAT trick: "After an initial run that generates all necessary
    symmetry-breaking clauses, these clauses are added to the problem encoding. A
    second run with a proof-producing SAT solver outputs a DRAT proof, which can be
    checked by existing proof verification tools."
  - Tooling: pysms `--ramsey 3 5` built-in; "SMS supports cube-and-conquer
    parallelization [23] via –simple-assignment-cutoff … and –cube-file".
- Smart cubing for SMS: "Smart Cubing for Graph Search: A Comparative Study",
  arXiv:2501.17201. https://arxiv.org/pdf/2501.17201 [REPORTED — not read]
- Relevance: SMS is the strongest *dynamic* isomorph-free search with a designed
  certificate story (nc-certificates + DRAT re-run); directly applicable to
  neighborhood-catalog generation in an R(5,5) attack. See also Section 4.6 for the
  Lean-verified SMS pipeline (IJCAR 2026).

### 1.6 Itzhakov-Codish: complete symmetry breaks for small Ramsey instances

- Itzhakov & Codish, "Breaking symmetries in graph search with canonizing sets",
  Constraints 21(3), 2016, DOI 10.1007/s10601-016-9244-z.
  https://link.springer.com/article/10.1007/s10601-016-9244-z
- Itzhakov & Codish, "Incremental Symmetry Breaking Constraints for Graph Search
  Problems", AAAI 2020. https://ojs.aaai.org/index.php/AAAI/article/view/5513
  [REPORTED from search summaries]: complete+compact symmetry break for order-10
  graphs with 7,853 lex constraints; first complete break for order 11; incremental
  extension n -> n+1.
- Relevance: complete symmetry breaks make small Ramsey UNSAT instances (e.g.
  R(4,4)<=18, R(3,5)) solvable directly and — crucially for us — a *statically added*
  symmetry break keeps the whole proof inside DRAT (unlike dynamic SMS pruning), at
  the price of scaling limits (order <= ~11 for complete breaks). For 45-vertex
  instances only partial/orbit-based breaks are feasible, which is why VeriPB
  (Section 4.3) matters.

### 1.7 Other recent SAT-for-Ramsey work (2024-2026)

- Wesley, "Lower Bounds for Book Ramsey Numbers", arXiv:2410.03625 (rev. 2025-09-05).
  https://arxiv.org/pdf/2410.03625 — SAT for Ramsey-type lower bounds. [REPORTED]
- "New bounds for some small multicolor Ramsey numbers", arXiv:2509.03784 (2025).
  https://arxiv.org/html/2509.03784 [REPORTED — not read]
- Przybocki, Mackey, Heule, Subercaseaux, "Doubly Saturated Ramsey Graphs: A Case
  Study in Computer-Assisted Mathematical Discovery", arXiv:2604.21187 (2026-04-23).
  https://arxiv.org/abs/2604.21187 — abstract [FETCHED-ABSTRACT]: "We present a method
  combining SAT solving with bespoke LLM-generated code to discover infinite families
  of such graphs, answering a question of Grinstead and Roberts from 1982. In
  addition, we use LLMs to generate and formalize correctness proofs in Lean."
  Relevance: the Heule/Mackey group's current SAT+LLM+Lean workflow on R(5,5)-adjacent
  objects (Ramsey-good graphs); shows who else is circling this territory.
- Kalfus & Lidický, "An automated proof that R(B8,B10)=37", arXiv:2606.05629
  (2026-06-04). https://arxiv.org/abs/2606.05629 — abstract [FETCHED-ABSTRACT]: "The
  problem as well as the proof were found with AutoMath, an AI-assisted mathematical
  discovery workflow … A Lean formalization of the upper-bound argument is available
  in the accompanying repository." Relevance: 37-vertex book-Ramsey upper bound with
  Lean artifact — evidence that ~40-vertex Ramsey-type upper bounds with formal
  artifacts are now being done.

---

## 2. Formal verification of Ramsey numbers

### 2.1 Gauthier & Brown: R(4,5)=25 in HOL4 (the direct precedent)

- arXiv:2404.01761; ITP 2024, DOI 10.4230/LIPIcs.ITP.2024.16.
  https://arxiv.org/abs/2404.01761 /
  https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ITP.2024.16
- Abstract [FETCHED-ABSTRACT]: "we prove this theorem in the interactive theorem
  prover HOL4 limiting the uncertainty to the small HOL4 kernel. Instead of verifying
  their algorithms directly, we rely on the HOL4 interface to MiniSat SAT to prove
  gluing lemmas. To reduce the number of such lemmas and thus make the computational
  part of the proof feasible, we implement a generalization algorithm. We verify that
  its output covers all the possible cases by implementing a custom SAT-solver
  extended with a graph isomorphism checker."
- PROVES: R(4,5)=25 end-to-end inside HOL4 (McKay-Radziszowski 1995 made formal).
- Relevance: the *only* fully formal exact classical Ramsey number at this difficulty;
  its architecture (SAT-proved gluing lemmas + verified case-cover via
  generalization + isomorphism checking) is the blueprint for formalizing an
  Angeltveit-McKay-style argument one level up.

### 2.2 R(3,3,3)=17 and small Ramsey numbers in Lean 4

- Narváez, Song, Zhang, "Formalizing Finite Ramsey Theory in Lean 4", CICM 2024,
  DOI 10.1007/978-3-031-66997-2_6 [VERBATIM-VERIFIED via dblp lookup].
- Repo: https://github.com/cruisesong7/formal_ramsey — README [VERBATIM-VERIFIED]:
  "The more general theory is `Ramsey.lean` which proves statements about multicolor
  Ramsey numbers. It also contains the proof of `R(3,3,3)=17`." Also: verified CNF
  encoders — "To obtain verifiable correct encodings of the Ramsey property as a
  Boolean formula in CNF run `lake exe RamseyEncoder <N> <s> <t>` … This is proven in
  the trestle framework." (trestle: https://github.com/FormalSAT/trestle)
- Relevance: Lean 4 + verified Ramsey CNF encoder already exists — a candidate base
  for our encoding-correctness layer.

### 2.3 Isabelle: diagonal Ramsey asymptotics (not exact values)

- Paulson, "Formalising New Mathematics in Isabelle: Diagonal Ramsey",
  arXiv:2501.10852. https://arxiv.org/abs/2501.10852 — formalizes the
  Campos-Griffiths-Morris-Sahasrabudhe exponential improvement on diagonal Ramsey
  upper bounds. [REPORTED — from abstract-level search] Relevance: shows ITP capacity
  for hard Ramsey *theory*; not directly a tool for exact small values.

### 2.4 Empty hexagon: the model for "verified encoding + checked SAT proof"

- Subercaseaux, Nawrocki, Gallicchio, Codel, Carneiro, Heule, "Formal Verification of
  the Empty Hexagon Number", ITP 2024, DOI 10.4230/LIPIcs.ITP.2024.35;
  arXiv:2403.17370. https://arxiv.org/pdf/2403.17370
- [REPORTED from search summary]: formalizes in Lean the reduction "if the CNF is
  unsatisfiable then every set of 30 points … contains an empty convex hexagon"; the
  17,300 CPU-hour SAT UNSAT result is checked by certificate, with the encoding
  verified in Lean. Relevance: the reference architecture for our Gate structure —
  Lean-verified reduction + industrial SAT certificate checking, at 10^4 CPU-hour scale.

---

## 3. Known attacks on R(5,5) <= 45 / reproductions of Angeltveit-McKay

### 3.1 CMU: srg(45,22,10,11) SAT search (live, 2026)

- Braginsky, Yellenki, Zhu (CMU; joint with John Mackey and Zachary Battleman),
  "Exploring srg(45,22,10,11) and implications for the Ramsey number R(5,5)", talk
  abstract, 57th Southeastern Intl. Conf. on Combinatorics, Graph Theory & Computing
  (March 2026). https://www.math.fau.edu/combinatorics/abstracts/braginsky-57.pdf
- Abstract [VERBATIM-VERIFIED from the PDF]: "We investigate the existence of an
  R(5,5)-good strongly regular graph with parameters (v,k,λ,µ) = (45,22,10,11)
  (equivalently, a K5-free simple graph whose complement is also K5-free); such a
  graph would imply R(5,5) = 46. We encode the search as a Boolean satisfiability
  (SAT) problem and, using symmetry reduction, decompose it into 313 SAT encodings,
  141 of which we certify as unsatisfiable. Additionally, we introduce a new encoding
  for enforcing K5-freeness that runs empirically faster than the standard approach."
- Status: PARTIAL (141/313 sub-encodings certified UNSAT as of the abstract). This is
  the existence side: finding such a graph would *disprove* our target R(5,5) <= 45.
  Conversely, ruling out this srg is a natural early sub-goal of any <= 45 proof
  (a (5,5,45)-graph need not be this srg, so completion of their search does not by
  itself give <= 45). Mackey is a frequent Heule collaborator — expect follow-up.
- Relevance: direct competitor/complement on the exact n=45 battlefield; their
  K5-freeness encoding trick and their symmetry decomposition are immediately
  reusable; also a scoop risk in the =46 direction.

### 3.2 No SAT reproduction of Angeltveit-McKay found

- Searches (Aug 2026) for a SAT-based or certified reproduction of R(5,5) <= 46 or
  <= 48 found nothing published or announced. The sat4math.com community index
  ("SAT for Mathematics is a community-maintained research index", 126 papers,
  2003-2026, https://sat4math.com/) lists no R(5,5) project. [REPORTED — absence of
  evidence from index + searches, not proof of absence]
- The A-M 46 paper itself contains SAT gluing but no certificates (Section 0.1),
  and states the authors believe the same method cannot reach 45 without "new
  theoretical insights" [VERBATIM-VERIFIED].

### 3.3 Commentary / expectations

- Gasarch (2024): "Getting R(5)≤45 may take 3000 years of CPU time. So it may not be
  for a while or ever." [REPORTED — blog]
  https://blog.computationalcomplexity.org/2024/09/progress-on-r5-will-there-be-more.html
  Note for calibration: A-M spent ~80 CPU-years total (30 + 50 replication)
  [VERBATIM-VERIFIED], so "3000 CPU-years" is a guess ~40x their spend — large but not
  out of reach for a cluster campaign if the method is SAT-parallelizable.
- Radziszowski's dynamic survey "Small Ramsey Numbers" (ElJC DS1) is the canonical
  bounds ledger; the 2024 edition is cited by A-M as [9] ("the 2024 edition of
  Radziszowki's survey"). https://www.cs.rit.edu/~spr/ElJC/ejcram18.pdf (site was
  unreachable at recon time; entry content [REPORTED]).

---

## 4. Proof-certificate tooling relevant at this scale

### 4.1 DRAT / LRAT trimming and checking

- DRAT-trim: Wetzler, Heule, Hunt, "DRAT-trim: Efficient Checking and Trimming Using
  Expressive Clausal Proofs" (SAT 2014); tool https://github.com/marijnheule/drat-trim.
  [REPORTED — standard reference]
- LRAT (annotated, machine-checkable-fast format) is the target format emitted by
  DRAT-trim and checked by verified checkers; modern solvers (CaDiCaL/Kissat) emit
  LRAT natively. [REPORTED — standard knowledge; verify exact solver flags at build time]
- Relevance: at A-M scale (10^12 subproblems), per-subproblem DRAT->LRAT trimming and
  aggregate bookkeeping is the throughput bottleneck to engineer around.

### 4.2 cake_lpr: verified checker with binary-level guarantees

- Tan, Heule, Myreen, "cake_lpr: Verified Propagation Redundancy Checking in CakeML",
  TACAS 2021, DOI 10.1007/978-3-030-72013-1_12; tool
  https://github.com/tanyongkiam/cake_lpr.
  https://link.springer.com/content/pdf/10.1007/978-3-030-72013-1_12.pdf
- [REPORTED from search summary + SAT competition checker doc
  https://satcompetition.github.io/2025/downloads/checkers/cakelpr.pdf]: checks LRAT
  and LPR; "verified using CakeML's binary code extraction toolchain, which yields
  correctness guarantees for its machine code (binary) implementation"; DRAT/DPR
  supported via DRAT-trim/DPR-trim preprocessing.
- Relevance: the default trusted checker for the conquer phase; used by SAT Competition.

### 4.3 VeriPB / pseudo-Boolean proof logging — the symmetry-breaking fix

- Bogaerts, Gocht, McCreesh, Nordström, "Certified Symmetry and Dominance Breaking
  for Combinatorial Optimisation", AAAI 2022, arXiv:2203.12275; extended JAIR 77,
  pp. 1539-1589 (2023). https://arxiv.org/pdf/2203.12275
- Tool: https://github.com/StephanGocht/VeriPB ; publications list
  https://veripb.org/publications.html
- Why it matters [REPORTED from search summaries of the above]: DRAT cannot
  efficiently express the *derivation* of symmetry-breaking constraints (adding them
  is sound only up to satisfiability-equivalence arguments outside DRAT); VeriPB's
  cutting-planes + redundance/dominance rules certify "automated static symmetry
  breaking in SAT, manual static symmetry breaking in CP, and automated dynamic
  dominance handling in maximum clique solving." This is exactly the gap in a
  symmetry-reduced R(5,5) encoding.
- Verified backend: CakePB — "CakePB is a kernel proof checker that has been formally
  verified in the HOL4 theorem prover using the CakeML suite"; workflow: "use VeriPB
  as a preprocessor to compile augmented proofs into kernel format for formal
  verification by CakePB." [REPORTED from SAT Competition docs
  https://satcompetition.github.io/2023/downloads/proposals/veripb.pdf and
  https://satcompetition.github.io/2025/downloads/checkers/veripb.pdf]
- PB<->DRAT bridge: Bryant/Heule et al., "Translating Pseudo-Boolean Proofs" (PBIP),
  FMCAD 2024, https://www.cs.cmu.edu/~mheule/publications/PBIP.pdf [REPORTED — not read].

### 4.4 PBLean: VeriPB certificates into Lean 4 (2026)

- Szeider, "PBLean: Pseudo-Boolean Proof Certificates for Lean 4", arXiv:2602.08692
  (v1 2026-02-09, v2 2026-04-02). https://arxiv.org/abs/2602.08692
- Abstract [FETCHED-ABSTRACT]: "a method for importing VeriPB pseudo-Boolean (PB)
  proof certificates into Lean 4. Key to our approach is reflection … Our checker
  supports all VeriPB kernel rules, including cutting-plane derivations,
  proof-by-contradiction subproofs, and redundance-based reasoning for symmetry
  breaking. … To derive theorems about the original combinatorial problems rather
  than about PB constraints alone, we support verified encodings."
- Relevance: closes solver->theorem gap for *symmetry-broken* encodings; strongest
  single tool for making an R(5,5) sub-result a composable Lean lemma.

### 4.5 LRAT-Catcher: LRAT + cube-and-conquer into Lean 4 (2026)

- Szeider, "LRAT-Catcher: Importing SAT Solver Certificates into Lean4 by
  Reflection", arXiv:2607.00815 (2026-07-01). https://arxiv.org/abs/2607.00815
- Abstract [FETCHED-ABSTRACT]: "a standalone, general-purpose tool that imports a
  DIMACS formula together with an LRAT certificate into Lean 4 as a theorem" using
  "the formally verified LRAT checker from Lean core as compiled native code via
  reflection"; supports "cube-and-conquer solving runs entirely inside Lean,"
  combining per-cube refutations "with a cover-completeness certificate into a single
  unsatisfiability proof"; evaluated on "establishing the Schur number S(4) = 44 and
  the Ramsey number R(4,4) = 18 as Lean theorems."
- Relevance: exactly the cube-cover + per-cube-certificate aggregation an R(5,5)
  campaign needs, already demonstrated on a (small) Ramsey number.

### 4.6 Verified SMS pipeline (2026)

- Kirchweger, Manrique, Szeider, "Formally Verified Graph Generation with SAT Modulo
  Symmetries and Lean", IJCAR 2026, DOI 10.1007/978-3-032-32589-1_8
  [VERBATIM-VERIFIED via dblp]. https://doi.org/10.1007/978-3-032-32589-1_8
- Relevance: makes SMS's nc-certificate + DRAT story land in Lean — candidate for
  certifying isomorph-free neighborhood catalogs (the analog of the R(4,5,n) censuses).

### 4.7 Scale reference points

- Pythagorean triples: ~200 TB DRAT, 68 GB compressed certificate (Section 1.2).
- Empty hexagon: 17,300 CPU-hours, certificate-checked, encoding verified in Lean
  (Section 2.4).
- A-M R(5,5)<=46: ~80 CPU-years, ~2×10^12 gluing ops, no certificates (Section 0.1).
- Implication: a certified <= 45 attempt is plausibly O(10^2-10^3) CPU-years of
  solving plus a comparable certificate-checking budget; certificate *volume*
  management (trim-on-the-fly, per-cube checking, no global proof file) is a
  first-class design constraint.

---

## 5. erdosproblems.com AI-collaboration workflow (process prior art)

Source: "Problem 728 and the use of AI on Erdős problems", Kevin Barreto, 26 Jan 2026,
https://www.erdosproblems.com/forum/thread/blog:2 (fetched raw; quotes
[VERBATIM-VERIFIED]).

Workflow described (GPT-5.2 Pro proof generation + Harmonic's Aristotle
autoformalization to Lean):

1. "Prompt the model with the problem and see if it finds any relevant literature or
   makes progress on it."
2. Re-prompt disguised as a competition problem: "This is a complex competition-style
   math problem. Solve the problem and give a rigorous proof or disproof. Do not
   search the internet." — because "if it discovers that a problem is an open research
   problem online, it will refuse to make a good attempt".
3. "Ask it to write its solution, formatted as a publishable maths paper, in a LaTeX
   code block, and then pass that TeX file to Aristotle to attempt to autoformalise."
4. Iterate Aristotle runs; "Once you've reached a full Lean file and Aristotle
   comments that it has formalised the solution to the problem at the top, then one
   checks the final main statement for accuracy to ensure it proved what was intended."

Gate rule: "I would not post an AI-generated proof unless I was very certain it was
correct and had a Lean formalisation of the proof".

**Verification pitfalls documented in the thread (each one maps onto our campaign):**

- Literature-search failure / novelty error, twice: problem [481] "it was later noted
  that it had been previously solved, which I was oblivious to"; problem [333] "I had
  a slip of judgement in announcing this prematurely on Twitter … I quickly retracted
  the claim when KoishiChan discovered it had been reported in the literature
  previously." -> our analog: verify live bounds before/at announcement (scoop check).
- Problem-statement misreading (quantifier direction): "it was noted that C was meant
  to be taken arbitrarily large, and I myself misread the problem, interpreting both
  C and ε as being sufficiently small." -> our analog: the formal statement of
  "R(5,5) <= 45" and of every gate lemma must be reviewed independently of the proof;
  a Lean proof of the wrong statement is the failure mode.
- Formalization checks the proof, not the theorem choice: the workflow's final check
  is explicitly on "the final main statement" — the residual trust gap.
- Suspicious lemma internals can survive informal review: "there were some things like
  a lemma where the model's proof seemingly had k being fixed in a proof where it
  should have been growing logarithmically that seemed a bit suspicious to me" —
  informal plausibility review is not verification.
- Context poisoning in autoformalizers: "passing the Lean file for [728] as extra
  context, thinking it could just copy relevant lemmata from there, actually just
  confuses Aristotle, rather than helping it."
- Site disclaimer: "All comments are the responsibility of the user. Comments
  appearing on this page are not verified for correctness."
- Community calibration (Nat Sothanaphan's reply): "at some level all of mathematics
  is such combination; the distinction is the difficulty level of the combination" —
  and novelty assessments shifted post hoc (Pomerance overlap).

---

## 6. Synthesis: strongest available tooling for this campaign

1. **Search/scale engine:** cube-and-conquer (Heule et al.) with CaDiCaL/Kissat
   emitting LRAT; SMS (or static Itzhakov-Codish-style breaks where complete breaks
   are feasible) for isomorph-free catalog phases; MathCheck-style SAT+CAS if we need
   nauty-grade canonization inside the loop.
2. **Certification spine:** VeriPB for anything symmetry-broken or counting/LP-flavored
   (cutting planes covers the A-M "excess" LP layer natively — DRAT does not);
   CakePB/cake_lpr as verified checkers; DRAT-trim/LRAT for plain UNSAT cubes.
3. **Formal capstone:** Lean 4 via PBLean (VeriPB->Lean) + LRAT-Catcher
   (LRAT+cube-cover->Lean) + trestle/formal_ramsey verified Ramsey encoders;
   HOL4 Gauthier-Brown architecture as the fallback pattern for gluing-lemma
   formalization.
4. **Known unknowns to resolve early:** (a) whether the A-M LP layer can be re-derived
   as PB cutting-planes reasoning at acceptable size; (b) certificate volume at
   10^11-10^12 subproblems; (c) whether the R(4,5,24) catalog (352,366 graphs,
   https://users.cecs.anu.edu.au/~bdm/data/ramsey.html) can be independently
   re-generated with an SMS/SAT+CAS certified pipeline (Gate 1); (d) CMU
   srg(45,22,10,11) progress — coordinate or race.
