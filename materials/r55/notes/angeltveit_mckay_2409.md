# Angeltveit–McKay, "R(5,5) <= 46" (arXiv:2409.15709) — proof structure recon

Recon notes for the R(5,5) campaign (first target: prove R(5,5) <= 45).
Date: 2026-08-13. Author of notes: Claude (subagent recon).

**Evidence discipline.** Unless tagged [REPORTED], every quote below was verified
directly against a primary source downloaded on 2026-08-13:

- Abstract page: <https://arxiv.org/abs/2409.15709> (v1 2024-09-24, v2 2025-09-01; no journal ref, no ancillary files listed).
- Full text (HTML, v2): <https://arxiv.org/html/2409.15709v2>.
- LaTeX source (ground truth for all quotes): <https://arxiv.org/e-print/2409.15709v2>
  (file `R55le46v3.tex`; v1 source `R55le46v1.tex` from <https://arxiv.org/e-print/2409.15709v1>).
- 2018 predecessor paper source: <https://arxiv.org/e-print/1703.08768> (file `Ramsey_paper2.0.tex`),
  abstract page <https://arxiv.org/abs/1703.08768>.
- Data page: <https://users.cecs.anu.edu.au/~bdm/data/ramsey.html> and
  <https://users.cecs.anu.edu.au/~bdm/data/r45extreme.tar.gz> (downloaded; graph counts re-counted locally).

Quotes keep the authors' LaTeX math notation. `\cR` = mathcal R (set of Ramsey graphs).

---

## 1. Main theorem and framing

> "Theorem 1.1. The Ramsey number R(5,5) is less than or equal to 46."
> — arXiv:2409.15709v2, Sec. 1. <https://arxiv.org/html/2409.15709v2>

Abstract (verbatim, v2):

> "We prove that the Ramsey number R(5,5) is less than or equal to 46. The proof uses a
> combination of linear programming and checking a large number of cases by computer.
> All of the computational parts of the proof were independently implemented by both
> authors, with consistent results."
> — <https://arxiv.org/abs/2409.15709>

Context (Sec. 1):

> "The lower bound of 43, which was established by Exoo [5] in 1989, is still the best.
> The upper bound of 49 was proved by the second author and Radziszowski [8], and this
> was improved by the authors [1] to 48."

Notation: `\cR(s,t,n)` = isomorphism classes of graphs on n vertices with no s-clique and no
independent t-set; `\cR(s,t,n,e>=e0)` etc. restrict edge counts; `e(s,t,n)`/`E(s,t,n)` = min/max
edge counts. `F_v^+` = neighbourhood (induced subgraph on neighbours of v), `F_v^-` = "dual
neighbourhood" (induced on vertices neither v nor adjacent to v).

Degree constraint driving everything (Sec. 1, using R(4,5)=25 [McKay–Radziszowski 1995]):

> "It follows that any vertex v of a graph in \cR(5,5,m) must have degree m-25 <= d(v) <= 24."

So for a hypothetical `F ∈ \cR(5,5,46)`: **every vertex has degree in {21,22,23,24}** (Sec. 2:
"For a hypothetical graph F ∈ \cR(5,5,46) every vertex must have degree d(v) ∈ {21,22,23,24}").

The subgraph-counting identity (from the m=2 case of McKay–Radziszowski 1997, Thm 2.2), Eq. (1.2):

> "excess(F) = 0, where excess(F) = \sum_{v∈V(F)} ( e(F_v^-) - e(F_v^+) - (1/2) d(v)(n-2d(v)) )."

## 2. Method: what gets glued to what (Sec. 2 + Figure 1)

Figure 1 structure (verified by rendering `r55pic46b.svg`/`r55pic46b.pdf` from the source):
two adjacent vertices a,b; K = common neighbourhood of a and b; H = N(a) = {b} ∪ K ∪ B;
G = N(b) = {a} ∪ K ∪ A; L = vertices adjacent to neither a nor b, with distinguished subsets
L1 (and L2 ⊆ L used in the second computation). Heavy dashed lines (the unknown edges the
gluing must fill in) run: A–B, K–L1, A–L1, B–L1, and inside L1.

> "Consider two adjacent vertices a,b. The neighbourhood of a is H = {b} ∪ K ∪ B and that of b
> is G = {a} ∪ K ∪ A. Both of these neighbourhoods are in \cR(4,5). L is the part of the graph
> adjacent to neither a nor b. The subgraph consisting of G and H, overlapping in K, plus the
> edge {a,b} but not the edges between A and B, will be denoted by G ∪_K H."

> "Using a mixture of theory and computation we compile a collection of pairs {(G,a),(H,b)}
> such that G_a^+ is isomorphic to H_b^+ and every graph in \cR(5,5,n) necessarily contains a
> pair in our collection overlapped as in the figure with a adjacent to b. We call a pair like
> (G,a) a pointed graph."

> "Next, for each pair of pointed graphs in our collection, we determine all the ways to fill
> in edges in the places indicated by heavy dashed lines in the figure without creating cliques
> or independent sets of size 5. We call this gluing along the edge ab."

Why L1 is nonempty here (unlike the 2018 <=48 proof):

> "In our previous work [1], our strategy was to start with L1 = ∅, so that only edges between
> A and B were sought. ... With n=46, choosing L1 = ∅ is impractical as the number of valid
> ways to add edges between A and B is extremely large. Instead, we chose a larger L1 and found
> that the number of solutions became manageable."

Summary of the 2018 <=48 proof as restated in this paper (Sec. 2):

> "First we determined the complete catalogue of \cR(4,5,24), which contains a total of
> 352,366 graphs. Given a hypothetical graph F ∈ \cR(5,5,48), either F or its complement must
> have a pair of adjacent vertices a,b of degree 24 whose neighbourhoods intersect in some
> subgraph K ∈ \cR(3,5,d) for d <= 11. ... Hence it suffices to consider all ways of gluing
> G ∪_K H for K ∈ \cR(3,5,d) where d <= 11 and G,H ∈ \cR(4,5,24). There is one gluing
> operation required for each pair of pointed graphs of type K and for each automorphism of K,
> and in total we computed approximately 2 trillion gluing operations."

Why the full-catalog approach cannot work at n=46 (Sec. 2, citing estimates from [7], updated in
their Appendix):

> "in [7] the authors estimated that |\cR(4,5,21)| ≈ 5.5×10^17, |\cR(4,5,22)| ≈ 1.9×10^15 and
> |\cR(4,5,23)| ≈ 10^11. ... Hence any approach that relies on the complete catalogue of these
> Ramsey graphs is impractical. Instead we use linear programming to reduce the number of
> graphs we need to consider. The basic idea is to determine \cR(4,5,n,e>=e0) for n=21,22,23
> and suitable e0 (which depends on n) and exclude these graphs by gluing along an edge as
> described above. Then we can use linear programming to finish the proof of Theorem 1.1."

> "In more detail, we determine the sets \cR(4,5,23,e>=119), \cR(4,5,22,e>=113) and
> \cR(4,5,21,e>=107). We also use \cR(4,5,24,e>=127), determined in [1]."

Two technical obstacles and their fixes (Sec. 2, enumerated):

> "(1) The method used in [7,1] to determine \cR(4,5,24) is too slow to determine \cR(4,5,n)
> for n=21,22,23 ... We explain our approach in Section 3.
> (2) The method used in [1] to glue two graphs along an edge produces far too many output
> graphs ... We get around that by adding extra vertices straight away when gluing along an
> edge as in [1]. In addition, we show in Section 5 that we do not have to perform all possible
> gluing operations."

### Role of "linear programming"

The LP is the excess-identity bookkeeping of Sec. 4 ("Some linear programming"): rewrite
excess(F)=0 grouped by degree class,

> "excess(F) = \sum_{d(v)=24}( (e(F_v^-)-104) + (127-e(F_v^+)) + 1 )
>            + \sum_{d(v)=23}( (e(F_v^-)-119) + (118-e(F_v^+)) + 1 )
>            + \sum_{d(v)=22}( (e(F_v^-)-135) + (112-e(F_v^+)) + 1 )
>            + \sum_{d(v)=21}( (e(F_v^-)-149) + (106-e(F_v^+)) + 1 )."

> "Given such an F ∈ \cR(5,5,46), suppose each F_v^+ ∉ E and each F_v^- ∉ \bar{E}. Then each
> vertex v ∈ V(F) contributes at least 1 to excess(F), so excess(F) >= 46 and hence we cannot
> have excess(F) = 0. This suggests a proof strategy: Deal with the relatively small number of
> graphs in E (and \bar{E}) separately, and then use the above equation for excess(F) to finish."

The weighted counting in Prop. 5.3's proof (α = 5m1+2m2+m3, β dual; α+β >= 46) is the actual
"LP"-style inequality juggling. No industrial LP solver is mentioned anywhere in the paper.

## 3. The sets E and the strata (Sec. 4–5)

Definitions (verbatim from LaTeX source, Sec. 4):

```
A   = \cR(4,5,24, e >= 127)
B_1 = \cR(4,5,23, e >= 121)   B_2 = \cR(4,5,23, e = 120)   B_3 = \cR(4,5,23, e = 119)
C_2 = \cR(4,5,22, e = 114)    C_3 = \cR(4,5,21, e = 113)   [sic — see typo note below]
D_3 = \cR(4,5,21, e = 107)
E   = A ∪ B_1 ∪ B_2 ∪ B_3 ∪ C_2 ∪ C_3 ∪ D_3
E_1 = A ∪ B_1,   E_2 = B_2 ∪ C_2,   E_3 = B_3 ∪ C_3 ∪ D_3
```

**Typo note (my analysis, not the paper's):** the source literally reads
`C_3 &= \cR(4,5,21, e = 113)` in both v1 and v2, but E(4,5,21)=107 (their Sec. 3.4 and Table 1),
so no 21-vertex Ramsey(4,5) graph has 113 edges. From Sec. 3.3 ("for the proof of Theorem 1.1
it suffices to consider \cR(4,5,22,e>=113)"), the excess rewrite (`112 - e(F_v^+)` for d(v)=22),
and |\cR(4,5,22,e=113)| = 30,976 published in the data files, C_3 must mean
**\cR(4,5,22, e = 113)**. Tagged: paper typo, meaning unambiguous from context.

`\bar A` etc. are complement sets, e.g. "\bar{A} = \cR(5,4,24, e <= 149)". Pointed graphs:
P(G) = the (up to) n pointed graphs of G; P_k(G) = P(G) after discarding the k-1 heuristically
most difficult pointed vertices:

> "Let P_k(G) denote the set of pointed graphs obtained as follows: First, find the n pointed
> graphs for G (without removing isomorphic graphs). Then sort them according to a heuristic
> measure of difficulty, and throw away the k-1 most difficult pointed graphs."

Supporting lemmas (all "explicit calculation"):

> "Lemma 5.1. Any graph G ∈ \cR(5,4,17) has a vertex of degree at least 8."
> Proof quote: "we can check this by completing a census of \cR(5,4,17,e<=59). There are 7147
> such graphs, and all of them have at least one vertex of degree greater than or equal to 8."

> "Lemma 5.2. Suppose G ∈ \cR(5,5,21) has two non-adjacent vertices of degree at most 4. Then
> either G contains a vertex of degree at least 8 or G contains a 4-clique {w1,w2,w3,w4} with
> deg(w1)+deg(w2)+deg(w3)+deg(w4) <= 24."
> Proof quotes: "there are 2029 graphs in \cR(5,4,16) with maximum degree at most 7 and a
> vertex of degree at most 4"; "there are 148 graphs in \cR(5,5,21) with maximum degree at
> most 7 and two non-adjacent vertices of degree at most 4."

Key structural reduction:

> "Proposition 5.3. Given F ∈ \cR(5,5,46), either F[E] or the complement \bar{F}[E] must
> contain one of the following:
> (1) Some vertex in E_1 adjacent to at least 1 other vertex in E.
> (2) Some vertex in E_2 adjacent to at least 1 other vertex in E_2.
> (3) Some vertex in E_2 adjacent to at least 5 other vertices in E.
> (4) Some vertex in E_3 \ D_3 adjacent to at least 8 other vertices in E.
> (5) Some vertex in D_3 adjacent to at least 8 other vertices in E \ D_3."

The proof of Prop. 5.3 is the pen-and-paper/LP heart: weights α = 5m1+2m2+m3 and β (dual),
excess(F)=0 forces α+β >= 46 (plus a refinement via n21 = number of degree-21 vertices in F[E]:
"α+β >= 46 + max(n21-\bar m1, 0) + max(\bar n21 - m1, 0). (This was the reason for considering
\cR(4,5,24,e>=127) rather than \cR(4,5,24,e>=128).)"), then case analysis with Ramsey-type and
inclusion-exclusion arguments.

**The five gluing strata** (this is the whole computational case split):

> "Theorem 5.4. To prove Theorem 1.1 it suffices to do the following gluings:
> (1) Glue P(E_1) to P(E);
> (2) Glue P(E_2) to P(E_2);
> (3) Glue P_5(E_2) to the rest of P(E);
> (4) Glue P_8(E_3 \ D_3) to P(E).
> (5) Glue P_8(D_3) to P(E \ D_3)."

Sizes of the underlying graph sets (Sec. 3; all verified against the published data files, see §6):

| set | definition | count |
|---|---|---|
| A | R(4,5,24,e>=127) | 2 + 3 + 32 + 147 + 843 + 3,401 = 4,428 (per-e counts below) |
| B1 | R(4,5,23,e>=121) | 2 + 119 = 121 |
| B2 | R(4,5,23,e=120) | 7,800 |
| B3 | R(4,5,23,e=119) | 332,778 |
| C2 | R(4,5,22,e=114) | 133 |
| C3 | R(4,5,22,e=113) [paper writes n=21; typo] | 30,976 |
| D3 | R(4,5,21,e=107) | 31 |
| E (total) | union | 376,267 graphs |

Per-e counts, Sec. 3.1–3.4 (verbatim numbers):
`|R(4,5,24,e=132)|=2, |e=131|=3, |e=130|=32, |e=129|=147, |e=128|=843, |e=127|=3,401` (from [1]);
`E(4,5,23)=122; |e=122|=2, |e=121|=119, |e=120|=7,800, |e=119|=332,778`;
`E(4,5,22)=114; |e=114|=133, |e=113|=30,976`;
`E(4,5,21)=107; |e=107|=31`.

Census method for these edge-extremal sets (Sec. 3.2):

> "As in [7], the idea is to glue \cR(3,5,p) to \cR(4,4,q) for p+q+1=23. ... consider all
> tuples (d_0,...,d_{q-1}) so that if we glue G to H with v_i adjacent to d_i vertices in G,
> the results lie in \cR(4,5,23,e>=119). Now we can order the vertices of H in a clever way
> ... and organise the set of such pairs (H,(d_0,...,d_{q-1})) into a tree."

> "Finding all of these graphs [R(4,5,23,e>=119)] was the second most time consuming part of
> the project, taking approximately 5 years of CPU time."

Flexibility note (Sec. 3): "we could consider \cR(4,5,22,e=114) only at the price of having to
consider \cR(4,5,23,e>=118) instead of \cR(4,5,23,e>=119). It also seems like it would suffice
to consider \cR(4,5,24,e>=128) rather than \cR(4,5,24,e>=127). But if we did that we would have
to weaken Proposition 5.3 and perform more of the more difficult gluing operations."

## 4. The two independent gluing computations

### First computation (Angeltveit, Sec. 6)

- Encoding: "Each gluing problem can be encoded as a SAT problem whose clauses forbid cliques
  or independent sets of size 5. We also added symmetry breaking clauses to distinguish the
  vertices of L1."
- Target size 37: "If |G| + |H| - |K| >= 37 we chose L1 to be empty, and if |G| + |H| - |K| < 37
  we chose L1 so that the output graphs (if any) have exactly 37 vertices."
- Heuristic: pick v ∈ K of maximal degree d in G ∪_K H, make v have degree 21 so its
  neighbourhood must lie in \cR(4,5,21) ("this provided an additional bottleneck").
- Solver: "Our SAT solver was a very simple special purpose solver without advanced features
  like clause learning and restarts."
- **Result:** "These gluing operations produced a total of 8,485,247 graphs in \cR(5,5,37).
  None of those extended to \cR(5,5,38)."

### Second computation (McKay, Sec. 7)

- Assumptions on L1 by size (|L1| = 1..5), using v ∈ K of min degree k_min and w ∈ K of max
  degree < 21; e.g. "|L1|=3 and |L|>=9: Since R(3,4)=9, L contains a triangle which we can take
  as L1"; "|L1|=4 and |L|>=13: By direct computation we find that every graph in \cR(5,4,13)
  contains an induced subgraph consisting of two triangles sharing an edge"; "|L1|=5 and
  |L|>=14: By direct computation, every graph in \cR(5,4,14) contains either B1 or B2 as an
  induced subgraph" (B1 = bowtie: two triangles sharing a vertex; B2 = bowtie + one extra edge;
  shown in the `bowtie.pdf` figure — note these B1/B2 are 5-vertex graphs, a notation collision
  with the sets B_1,B_2 of Sec. 4). Also: "none of the 413 graphs in \cR(5,4,14) without B2 as
  an induced subgraph extend to \cR(5,4,24)."
- Phased search with |L2| = 10 (auto-degrading to 9): (1) propagation + limited backtracking;
  (2) residual SAT instance to Glucose with 1-minute limit ("Glucose decided satisfiability in
  1 second on average"); (3) rare leftovers by Glucose + backtracking. "As a sanity check, a
  number of cases were run as well using Kissat [4] instead of Glucose."
- **Result:** "In the first phase, the total number of (G,H)-overlaps was about 12 million, of
  which 5.6 million were processed by Glucose. Glucose found 15,248 satisfiable cases, and the
  remainder were unsatisfiable including 77 which timed out and needed step (3). The 15,248
  satisfiable cases went to the second phase, where only 17 cases were found to be satisfiable
  and were sent to the third phase where they were unsatisfiable. This proves that none of the
  initial configurations can be extended to \cR(5,5,46)."

### Independent-verification setup (Sec. 2 + Conclusions)

> "In the interests of confidence, all the computations were repeated by the second author
> using independent programs and usually with different methods. This replication took about
> 50 years of additional CPU time. For some tasks, such as creation of the catalogues described
> in the following section, the outputs of the two implementations could be directly compared.
> For the remaining much more complicated stages, the two approaches were deliberately intended
> to be disjoint to avoid the well-known axiom of software engineering that two programmers
> implementing the same algorithm tend to make errors in the same places. Moreover, if a graph
> \cR(5,5,46) exists, each implementation should have found many of its subgraphs by their
> completely different approaches."

(v1 additionally claimed the replication "produced identical results"; that phrase was dropped
in v2 in favour of the fuller explanation above — verified by diffing the two arXiv sources.)

> "Given the theory and the two independent computations, Theorem 1.1 has now been firmly
> established. We believe that future improvements to the upper bound on R(5,5) will need new
> theoretical insights, as an excessive amount of computer time would be required to apply the
> same method." — Sec. 8.

## 5. CPU-cost attribution: 2 trillion gluings (2018, <=48) vs 50 CPU-years (2024, <=46)

**Resolved with primary quotes from both papers:**

- **~2 trillion gluing operations = the 2018 proof of R(5,5) <= 48** (arXiv:1703.08768,
  J. Graph Theory 89(1):5–13, 2018). The 2018 paper itself:
  > "The proof of Theorem 1 is via computer verification, checking approximately two trillion
  > separate cases." (`Ramsey_paper2.0.tex`, Introduction)
  > "we needed to solve approximately 2 trillion gluing problems. While that is certainly a
  > lot, we were able to perform hundreds of thousands of such gluings per second per core.
  > The whole calculation took approximately six core-months for one implementation and two
  > core-months for the other." (Sec. on gluing)
  Also from 2018: completing the R(4,5,24) census "took about 1.5 core-years of computer time
  and discovered 1462 new graphs in \cR(4,5,24)" plus "another 6 core-months to sanity-checking
  of the completed catalogue." So the entire 2018 computation was of order ~2–3 core-years.
  The 2024 paper repeats the attribution: "in total we computed approximately 2 trillion gluing
  operations" (describing the [1]=2018 proof).

- **The CPU-decades belong to the 2024 proof of R(5,5) <= 46** (arXiv:2409.15709):
  > "We estimate that completing the census of \cR(4,5,n,e>=e0) took approximately 15 years of
  > CPU time while gluing the necessary graph took another 15 years of CPU time. Hence the
  > whole project took about 30 years of CPU time for the first author to complete."
  > "This replication took about 50 years of additional CPU time."
  So: first author's primary computation ~30 CPU-years (~15 census + ~15 gluing); second
  author's independent replication ~50 CPU-years; total ~80 CPU-years for the <=46 proof.
  Within the census, R(4,5,23,e>=119) alone took "approximately 5 years of CPU time" and was
  "the second most time consuming part" (the gluing phase being the largest).

- The 2024 paper does NOT publish a total gluing-operation count for the <=46 computation; the
  only published volume figures for its gluing are the second computation's "about 12 million"
  first-phase (G,H)-overlaps and the first computation's 8,485,247 output graphs in R(5,5,37).

## 6. Code and data availability

- **Paper's own pointer** (Sec. 3.4): "the graphs themselves, including additional classes, are
  available on the internet [2]" where [2] = "Vigleik Angeltveit and Brendan D. McKay.
  Edge-extremal Ramsey (4,5)-graphs. <https://users.cecs.anu.edu.au/~bdm/data/ramsey.html>."
- **Data page contents** (fetched 2026-08-13, verbatim):
  - "In 1995, McKay and Radziszowski proved that there are no Ramsey(4,5)-graphs with more than
    24 vertices and found 350904 of them with 24 vertices. The remainder were found in 2016 by
    McKay and Angeltveit. There are 352366 altogether, see r45_24.g6." (350,904 + 1,462 = 352,366
    matches the 2018 paper's "discovered 1462 new graphs".)
  - "We also provide, for each number of vertices 4-23, the complete sets of Ramsey(4,5)-graphs
    with the smallest few edge counts and the largest few edge counts. See r45extreme.tar.gz."
  - Complete catalogs of R(3,5,n) (n=1..13) and R(4,4,n) (n=1..17) — the gluing inputs for the
    census method — plus R(3,4), R(3,6)..R(3,9), and "r55_42some.g6 contains 328 of these
    [Ramsey(5,5,42)] graphs; the other 328 are their complements."
- **I downloaded `r45extreme.tar.gz` (~90.6 MB) and re-counted the key files** (graph6 format,
  one graph per line): `r4521.107.g6` = 31, `r4522.113.g6` = 30,976, `r4522.114.g6` = 133,
  `r4523.119.g6` = 332,778, `r4523.120.g6` = 7,800, `r4523.121.g6` = 119, `r4523.122.g6` = 2.
  **All match the paper's published counts exactly.** (`r4521.106.g6` also present, matching
  Table 1's N(emax-1)=10,188 for n=21.)
- **No code published.** The arXiv listing has no ancillary files (checked abs page and both
  e-print tarballs: only .tex + figure PDFs). The data page offers no programs ("gluing"
  programs, the special-purpose SAT solver, and the census code are not distributed there).
  McKay's site distributes nauty/Traces separately at <https://users.cecs.anu.edu.au/~bdm/nauty/>
  [REPORTED — well-known location; not re-verified in this session]. Third-party SAT solvers
  used: Glucose (Audemard–Simon) and Kissat (Biere et al.), both public.
- **Journal status:** no journal reference or DOI beyond the arXiv DOI is listed on the abs
  page as of 2026-08-13; a web search found no published journal version. [REPORTED — absence
  of evidence from arXiv page + search, not proof.]

## 7. Appendix A: census of R(4,5) (Table 1 highlights)

> "Floating-point numbers in the final column are estimates obtained by a statistical method
> described by the first author [6]. ... In total, |\cR(4,5)| ≈ 2.93×10^19."

Selected exact rows (n; e_min; e_max; N(e_min); N(e_min+1); N(e_max−1); N(e_max); |R(4,5,n)|):

```
21; 77; 107;  83; 5,940; 10,188;  31; ≈5.6×10^17
22; 88; 114;   3;    94; 30,976; 133; ≈1.8×10^15
23; 101; 122;  1;    76;    119;   2; ≈9×10^10
24; 116; 132;  9;    90;      3;   2; 352,366 (exact)
```

(Updated estimates vs [7]: R(4,5,21) ≈ 5.6×10^17, R(4,5,22) ≈ 1.8×10^15, R(4,5,23) ≈ 9×10^10.)

## 8. Reproduction targets (my assessment)

Ranked smallest-first, each with a published ground-truth count:

1. **Lemma 5.1 census: R(5,4,17, e<=59), published count 7,147 graphs, all with max degree >= 8.**
   Self-contained, no gluing infrastructure needed beyond isomorph-free generation; a good
   warm-up but not a "stratum" of the main computation.
2. **D3 = R(4,5,21, e=107) census: published count 31** (paper Sec. 3.4: "|R(4,5,21,e=107)|=31";
   "this calculation was faster still"). This is the smallest genuine stratum of the proof's
   catalog stage; inputs (complete R(3,5,p) and R(4,4,q) catalogs, p+q+1=21) are on the data
   page, and the answer can be diffed graph-by-graph against `r45extreme/r4521.107.g6`
   (31 graphs, verified locally). **Recommended first reproduction.**
3. **C2 = R(4,5,22, e=114): 133 graphs** (`r4522.113.g6`/`r4522.114.g6` for the e>=113 set),
   then B2/B1 at n=23, then the big one B3 (332,778; ~5 CPU-years by their method).
4. Smallest gluing stratum: **Theorem 5.4(2), glue P(E2) to P(E2)** with |E2| = 7,800 + 133
   = 7,933 graphs. Caveat: the paper publishes no per-stratum output counts for the gluing
   stage, so the only checkable aggregate numbers are the first computation's total
   8,485,247 graphs in R(5,5,37) (across all strata) and the second computation's phase counts
   (~12M overlaps / 15,248 / 17 / 0).
5. 2018-style sanity anchor: the R(4,5,24) catalog (352,366 graphs, `r45_24.g6`) and
   E(4,5,24)=132 with |R(4,5,24,e=132)|=2 — relevant to campaign task "Gate 1".

## 9. Implications for an R(5,5) <= 45 attack (brief)

- At n=45 every vertex of F ∈ R(5,5,45) has degree in {20,...,24}; the excess identity now has
  five degree classes and the e0 thresholds must be re-derived (the d(v)=20 term is
  e(F_v^-) - e(F_v^+) - (1/2)·20·(45-40) = ... - 50, and F_v^+ ∈ R(4,5,20), F_v^- ∈ R(5,4,24)).
  R(4,5,20) is the worst regime (≈8.57×10^18 graphs, e_max=100 with N(100)=1 — Table 1), so the
  edge-extremal slices needed at n=20 may be the new bottleneck.
- The authors' own verdict (Sec. 8): future improvements "will need new theoretical insights,
  as an excessive amount of computer time would be required to apply the same method."

## 10. Local artifacts from this recon

Scratchpad (session-local, may be garbage-collected):
`/private/tmp/claude-501/-Users-jinleic-jinleic-workspace/affa935d-1114-4410-b38c-41022bc8bb31/scratchpad/r55/`
contains `paper.html`/`paper.txt` (arXiv HTML + extracted text), `srcv1/`, `srcv2/` (LaTeX
sources of 2409.15709 v1/v2), `src48/` (LaTeX source of 1703.08768), `r45extreme.tar.gz` +
extracted counts, `fig1.svg(.png)`, `bowtie.pdf(.png)`, `ramsey_page.html`.
