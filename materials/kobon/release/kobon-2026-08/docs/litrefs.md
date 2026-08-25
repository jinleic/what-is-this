# Literature pins for the Kobon capacity paper

Source-audited 2026-08-20 (scout `LiteraturePin`). VERBATIM-SOURCE standard:
quote only text returned from the directly linked page/PDF; preserve source
typos and label corrections separately.

## 1. K(14) >= 53 provenance and attribution

**Verdict:** the lower-bound construction is attributed to **Johannes
Bader**, not to Gilles Clément. The correct author names in the upper-bound
paper are **Gilles Clément and Johannes Bader**. No accessible source in
this audit uses "F.-Z. Bader" or "R. Clément"; treat those initials as
transcription/OCR errors and do not cite them.

Primary Bader page:
https://www.sop.tik.ee.ethz.ch/people/baderj/other.html. Verbatim page
text: "The following table lists other configurations I found (click on the
images to enlarge them):" followed by "n: 10 n: 14", "triangles: 25
triangles: 53", "best known: 25 best known: ?", and "upper bound: 26 upper
bound: 56". The same page links `kobonTriangles/14_KobonTriangle.pdf` and
`14_KobonTriangle_large.png`.

Bader cached PDF: https://oeis.org/A006066/a006066_1.pdf, page 2. Extracted
table text says: "Other Configurations ... n: 10 n: 14 triangles: 25
triangles: 53 best known: 25 best known: ? upper bound: 26 upper bound:
56". It is an image/vector-style PDF; its text extraction contains no
Cartesian line-equation list.

Independent vector/image corroboration:
https://zegalur.github.io/line-order/gallery/kobon.html has the exact
heading "14-Line Solution (53 Triangles) by Johannes Bader" and publishes a
14-row line-order table. The associated straight-line SVG is
https://raw.githubusercontent.com/zegalur/line-order/master/gallery/imgs/
kobon_14_53tri_lines.svg; its header literally contains "(14 lines, 53
triangles, 54 upper bound)" and "by Johannes Bader", followed by 14 SVG
`<line>` elements with 500x500 rendering coordinates and triangle paths.
This is an accessible vector rendering/coordinate-like artifact generated
by LineOrder, not a canonical exact/rational coordinate list.

OEIS plain-text source: https://oeis.org/search?q=id:A006066&fmt=text. Its
current table header is "The known values a = a(n) and upper bounds U
(usually A032765(n)) with name of discoverer of the arrangement when known
are as follows:" and "n a U [Found by]"; current row is `14 >= 53 54
[Bader]`. The current `%E` history note says "Updated with results from
Johannes Bader (johannes.bader(AT)tik.ee.ethz.ch), Dec 06 2007, who says
\"Acknowledgments and dedication to Corinne Thomet\"." This is the
strongest OEIS provenance trail for the 53 arrangement.

Wikipedia https://en.wikipedia.org/wiki/Kobon_triangle_problem currently
does not list n=14/53 in its prose. It cites C-B for strictness in n = 0
or 2 mod 6, BBL for the stronger even bound, OEIS A006066, and Savchuk;
Bader appears only in the archived external link. Therefore Wikipedia is
corroboration/context, not the originating source for 53.

**Coordinate-source status:** no accessible PDF with a textual Cartesian
or exact-rational coordinate list for Bader's 53-triangle arrangement was
found. Mark the requested coordinate item `<no accessible source>`. The
Bader image/PDF and LineOrder SVG are valid visual/vector evidence; do not
call their decimal drawing coordinates exact or canonical.

**Current-vs-proposed caveat, corrected 2026-08-20:** published OEIS
A006066 remains approved and says `>=53 [Bader]`. OEIS history records
proposed Aug-2026 edits for a(14)=54 by Andrea Maiorana, with exact
rational data at https://github.com/rufio72/kobon_triangles_k14
(commit e47c7cfd); that README says "Status: NOT yet independently
confirmed" — **however, this campaign independently verified the claim on
2026-08-20**: SHA-pinned sol1 (83fc2666...) plus full sol2..sol15 sweep,
all exact-rational counts = 54, including the audited campaign engine
check (`scratch/kobon/n14/maiorana_verify.log`, exit-0 transcript;
certificate `scratch/kobon/n14/Kgen14_ge54_certificate.md`). Statement for
our paper: `K_gen(14) >= 54` is verified local evidence; the source's own
README status and the unapproved OEIS history remain exactly that — do not
report the OEIS row as changed; report the arrangement as independently
verified by us.

## 2. Upper-bound sources and exact statements

### Blanc (correct identifier is arXiv:0801.2845, not 1012.1931)

Official publication list:
https://algebra.dmi.unibas.ch/blanc/publications.html lists "The best
polynomial bounds for the number of triangles in a simple arrangement of n
pseudo-lines", *Geombinatorics* 21 (2011), no. 1, 5-14. Full paper:
https://algebra.dmi.unibas.ch/blanc/articles/bestbound.pdf. ArXiv metadata:
https://arxiv.org/abs/0801.2845. The mistaken identifier
https://arxiv.org/abs/1012.1931 is unrelated: "Charged particle elliptic
flow in p+p collisions at LHC energies in a transport model PACIAE"
(Dai-Mei Zhou et al., nucl-th).

Verbatim Theorem 1 (mathematical display transcribed faithfully): "Theorem
1 (polynomial bounds for affine arrangements). Let A be an affine
arrangement of n pseudo-lines. Then"

    a_3(A) <= n(n-5/2)/3      if n == 0,4 (mod 6),
              (n(n-2)-2)/3    if n == 1   (mod 6),
              (n(n-5/2)-2)/3  if n == 2   (mod 6),
              n(n-2)/3        if n == 3,5 (mod 6).

"Furthermore, for any integer 0 <= k <= 5, there exist infinitely many
integers n such that the bound above is reached, for n congruent to k (mod
6)."

Verbatim Theorem 3 content: "The bound of Theorem 1 is reached for any
integer n <= 30, n != 11,12; and we have a_3^s(11) = 32 and a_3^s(12) =
37." The paper's a_3 counts bounded triangular regions in **simple affine
arrangements of pseudo-lines** (and distinguishes affine/projective and
line/pseudoline cases). This is not a theorem for the broad
crossed/degenerate K_gen convention. In particular, do not cite Blanc as
proving broad K(14)=53.

### Rasukaru 2005 web analysis (newly located 2026-08-22)

Primary archived Japanese pages, dated September 22--October 8, 2005:

- https://web.archive.org/web/20061111125706/http://www10.plala.or.jp/rascalhp/nlines.htm
- https://web.archive.org/web/20221209001801/http://www10.plala.or.jp/rascalhp/nlines2.htm
- https://web.archive.org/web/20221208231339/http://www10.plala.or.jp/rascalhp/nlines3.htm

Honma's page, which Clément--Bader cite as reference [6], explicitly credits
the analysis to forum participant “らすかる” (Rasukaru) and dates the
conclusions:
https://web.archive.org/web/20171111045123/https://www004.upp.so-net.ne.jp/s_honma/triangle/triangle2.htm.

Part 1 uses a sequential/per-line multipoint charge: adding a third line
through an existing crossing destroys at least three segments and creates at
most two new shared-side credits; later lines lose 5, 7, ... segments. It
then gives the endpoint-loop parity obstruction for even n. Part 2 claims
K(11)=32 and the mod-6 -1 bound. Part 3 gives a detailed n=10 case analysis
claiming K(10)=25 and states the even formula

    floor((n(n-2) - floor((n+2)/4)) / 3),

which gives upper bound 54 at n=14. The final summary explicitly says the
conditions/reductions in 5-1 and 5-3 are “probably” ignorable, so the general
formula is not presented with a complete proof. The pages are public prior
art, not refereed papers or proof certificates.

Priority consequence: do not call the campaign's n=6,8,10 results first
mathematical determinations. They are the first machine-checkable,
all-degeneracy certificates located by this audit. The web argument motivates
the face-specific square penalty

    3F <= n(n-2) - 2Q - sum_p (k_p-2)^2,

but the penalty is false even for triangular faces. A generic exact-rational
chart of the classical simplicial arrangement A(12,1) has Q=0, N3=15, N6=1,
and 30 triangular bounded faces, so its two sides are 90 and 89. The exact
replay is `scratch/kobon/square_penalty_counterexample.py`; the earlier
4,000-arrangement audit missed this structured simplicial family. More
generally A(2m,1), defined as the m side lines and m mirror axes of a regular
m-gon, violates the penalty for every m >= 6. The canonical source read
directly is Branko Grünbaum, *A catalogue of simplicial arrangements in the
real projective plane*, Ars Math. Contemp. 2 (2009), 1-25,
doi:10.26493/1855-3974.88.e12. Pages 2 and 4 define the simplicial family
R(1)=A(2m,1) and give, for m>3, t2=m, t3=m(m-1)/2, and t_m=1. For even m,
half the axes join opposite vertices and half join opposite edge midpoints;
for odd m, every axis joins a vertex to the opposite edge midpoint. The proof
uses only m>=6, and the exact m=6 replay constructs three axes of each even
type.

Identifier correction: arXiv:0904.1244 is an astrophysics paper, not
Grünbaum's catalogue; arXiv:1011.1862 is a lattice-QCD paper, not Cuntz.
Cuntz's *Simplicial arrangements with up to 27 lines* is arXiv:1108.3000.
The paper cites the verified published Grünbaum source, not any of these
misassigned identifiers.

### Clement-Bader 2007 ETH draft

PDF: https://www.sop.tik.ee.ethz.ch/publicationListFiles/cb2007a.pdf; OEIS
cache: https://oeis.org/A006066/a006066.pdf. The first page is headed
"Tighter Upper Bound for the Number of Kobon Triangles", "Draft Version,
Timestamp: December 21, 2007 19:03", "Gilles Clement" and "Johannes Bader".

Verbatim abstract claim: "In this paper we present the concept of perfect
configuration of lines which are then used to proof that the bound known
of Tamura can not be reached for all n with n congruent to 0 mod 6 and n
congruent to 2 mod 6." It then displays the new bound as n(n-2)/3 for n mod
6 in {3,5}, (n+1)(n-3)/3 for n mod 6 in {0,2}, and (n-1)^2/3 for n mod 6
in {1,4}.

Verbatim Theorem 1 prose/display: "The maximal number of Kobon triangles
K(n) for a given number of n straight lines in a plane is upper bounded by"

    K(n) <= floor(n(n-2)/3) - I_{n: n mod 6 in {0,2}}(n),

"where I_A(x) denotes the indicator function. In other words the upper
bound known by Saburo Tamura cannot be reached for all n with n congruent
to 0 mod 6 and n congruent to 2 mod 6."

The draft's final n congruent to 1,4 display is a typo:
https://mathworld.wolfram.com/KobonTriangle.html explicitly states, "where
the above expression (n^2-2n-2)/3 was incorrectly written as (n-1)^2/3 in
the paper." Use the corrected value (n^2-2n-2)/3 only as a corrected
transcription, never as a verbatim C-B quote. For n=14, C-B gives at most
55. Its Lemma 1 proof says one multipoint belongs to at most two shared-side
pairs; the campaign's exact cevian-flower construction gives 2k such pairs at
a k-fold point using genuine triangular faces. This invalidates the stated
proof step, not the numerical bound. Rasukaru's earlier sequential/per-line
argument is different and is not refuted by the flower alone; the
$A(2m,1)$ family above does refute the coefficient-1 square inequality that
the sequential credit suggests. Neither counterexample by itself refutes the
stated mod-6 numerical bound. C-B is an unpublished draft; do not treat its
theorem as a certified general-arrangement upper bound.
The exact counterexample changes the proposed proof route, not this source-status conclusion.

### BBL context for the current even upper bounds

https://arxiv.org/pdf/0706.0723, Theorem 1.1 verbatim: "If n is an even
integer, then a_3^s(n) <= floor(n(n - 7/3)/3). The bound of Theorem 1.1 is
reached for 4, 6, 10, 16 pseudo-lines but not for 8, 12, 14 pseudo-lines."
This is a simple affine pseudoline triangular-face quantity, not broad
K_gen. For n=14 it yields 54; for n=18 it yields 94; for n=20 it yields
117.

### Parpalak-Utkin 2026 (two preprints)

arXiv:2607.29236 (Parpalak, Utkin; 2026-07-31) "Enumeration and
Classification of Triangle-Maximal Pseudoline Arrangements": exhaustive
enumeration of simple pseudoline arrangements (n odd) attaining
floor(n(n-2)/3), w0-wiring-diagram DFS over even-indexed generators.
Verbatim: "The zero entry at n = 11 shows the known non-existence of
perfect arrangements in this case." Table 1 exact counts: n=13 -> 47;
n=17 -> 85; n=19 -> 107; n=27: 85,562,064 wiring diagrams in 56,646
projective classes. Corollary 7.1: a bound-attaining simple pseudoline
arrangement exists for every odd n <= 89 except n=11. Convention note
(verbatim): "Stretchability is not considered: the results are
combinatorial, leaving the straight-line realization as a separate
question." Scope: bound-attaining family only — its n=11 #D=0 rules out
the upper target 33 at pseudoline level but does NOT by itself state
K_pseudo(11)=32 (needs the cross-source lower bound 32, e.g. Savchuk's
constructed 32).

arXiv:2604.22035 (same authors; 2026-04-23, NOT July) "The 18*2^t+1
Triangle-Maximal Series of Straight Lines": Theorem 4.1: explicit 19-line
simple straight arrangement with exactly 107 bounded triangular faces,
certified by computer-assisted interval/combinatorial checks; Theorem
5.1: infinite series n=18*2^t+1 with a3 = 108*4^t - 1 (n=19: 107; 37:
431; 73: 1727). One-step examples: n=41 -> 533, n=45 -> 645, both
matching the upper bound; n=21/23/27 base-config searches returned zero
(strong negative evidence). Scope: simple straight lines face count.

The 2607 paper explicitly excludes the present branch: "We restrict to odd
\(n\); the even case (pairs of parallel lines in the Arnold problem,
non-simple arrangements in the Kobon problem) is beyond the scope of this
paper." Section 3.1 further says even optimal arrangements involve parallel
pairs or triple points and that simple even arrangements have many unavoidable
defects. Its defect and wall-debt methods are therefore research prompts, not
universal clauses for \((12,39)\).

### Rote NumPSLA 2025 and colored-pseudoline signotopes 2026

Rote, arXiv:2503.02336, "NumPSLA — An experimental research tool for
pseudoline arrangements and order types", enumerates small simple
x-monotone pseudoline arrangements by incremental cutpaths. The abstract says
12-point abstract order types and 11-pseudoline arrangements are practical;
Section 1.1 explicitly assumes general position. The code's exclude database
is a set difference against an order-type database, not a nonrealizability
proof certificate. It was used for representation comparison only.

Radtke--Keszegh--Lauff, arXiv:2601.20574, "On Triangles in Colored
Pseudoline Arrangements", works with simple pairwise-crossing pseudolines.
Definition 1 records the rank-three signotope four-set sign-change axiom and
the text identifies triangular-face flips with changing one sign. This
motivated a zero-aware projective chirotope audit, but none of its colored
triangle theorems supplies an even non-simple Kobon upper bound.

### Savchuk 2025

Paper: https://arxiv.org/abs/2507.07951; full PDF:
https://arxiv.org/pdf/2507.07951. Verbatim abstract method/result claims:
"We present new methods and results for constructing optimal Kobon triangle
arrangements. First, we introduce a compact table notation for describing
arrangements of pseudolines, enabling the representation and analysis of
complex cases, including symmetrical arrangements, arrangements with
parallel lines, and arrangements with multiple-line intersection points.
Building on this, we provide a simple heuristic method and tools for
recovering straight-line arrangements from a given table, with the ability
to enforce additional properties such as symmetries. ... we develop a tool
that transforms the search for optimal Kobon arrangement tables into a SAT
problem, allowing us to leverage modern SAT solvers (specifically Kissat)
... Using these techniques, we find new optimal Kobon arrangements for 23
and 27 lines, along with several other new results."

Verbatim Conclusion: "The approach described in this paper has proven
productive, enabling the discovery of new maximal Kobon arrangements (for
n in {23,27}). For n = 11 we confirm that arrangement with 33 triangles
cannot be built even with pseudolines. For n in {3,5,9,15,17} we enumerate
all possible Kobon arrangements."

Methods stated in the paper: tables encode each line's intersection order
and support parallel/multiple-line cases; CNF variables include G
(immediately-after), X (somewhere-after), A (column), and optional M
(missing triangle); Kissat solves generated DIMACS instances; straightening
uses line equations x cos(a_i)+y sin(a_i)+C_i=0 and constrained
SciPy/NumPy minimization; LineOrder generates the displayed straight-line
images. The paper makes no claim settling n=14.

## 3. Current OEIS A006066 statuses (published text)

Use plain text URL https://oeis.org/search?q=id:A006066&fmt=text or HTML
https://oeis.org/A006066. Current exact rows under `n a U [Found by]` are:

- `14 >= 53 54 [Bader]` - lower bound 53, upper bound 54; open in the
  published record.
- `15 65 65 [Suzuki]` - exact because a=U.
- `18 >= 93 94 [Bader]` - lower bound 93, upper bound 94; open.
- `20 >= 116 117 [Wood]` - lower bound 116, upper bound 117; open.

The `>=` is part of the a-column and means a construction/lower bound
only; U is the separate upper-bound column. A bare equal a and U records
an exact value; `?` records unknown. `[Found by]` is attached to the
arrangement in the row, not proof of the upper bound. The proposed 2026 54
record is **independently verified by this workspace (2026-08-20)** (SHA-pinned
sol1 + replayable 15/15 sweep; `scratch/kobon/n14/Kgen14_ge54_certificate.md`),
but remains **unapproved OEIS history**: report the approved row text and our
local verification separately, exactly as done in §1's corrected caveat.

## 4. Recommended citation entries

```bibtex
@misc{BaderKobonTriangles,
  author       = {Johannes Bader},
  title        = {Kobon Triangles},
  howpublished = {ETH Zuerich Systems Optimization archived web page},
  url          = {https://www.sop.tik.ee.ethz.ch/people/baderj/other.html},
  note         = {14-line, 53-triangle configuration; page also links
                  14_KobonTriangle.pdf and image}
}

@unpublished{ClementBader2007,
  author = {Gilles Cl{\'e}ment and Johannes Bader},
  title  = {Tighter Upper Bound for the Number of Kobon Triangles},
  note   = {Draft Version, December 21, 2007},
  url    = {https://www.sop.tik.ee.ethz.ch/publicationListFiles/cb2007a.pdf}
}

@article{Blanc2011,
  author  = {J{\'e}r{\'e}my Blanc},
  title   = {The best polynomial bounds for the number of triangles in a
             simple arrangement of n pseudo-lines},
  journal = {Geombinatorics},
  volume  = {21}, number = {1}, pages = {5--14}, year = {2011},
  url     = {https://algebra.dmi.unibas.ch/blanc/articles/bestbound.pdf},
  note    = {Preprint arXiv:0801.2845}
}

@incollection{BartholdiBlancLoisel2008,
  author    = {Nicolas Bartholdi and J{\'e}r{\'e}my Blanc and
               S{\'e}bastien Loisel},
  title     = {On simple arrangements of lines and pseudo-lines in P^2 and
               R^2 with the maximum number of triangles},
  booktitle = {Surveys on Discrete and Computational Geometry},
  series    = {Contemporary Mathematics}, volume = {453},
  pages     = {105--116},
  publisher = {American Mathematical Society}, year = {2008},
  doi       = {10.1090/conm/453/08797},
  eprint    = {0706.0723}
}

@misc{Savchuk2025,
  author        = {Pavlo Savchuk},
  title         = {Constructing Optimal Kobon Triangle Arrangements via
                   Table Encoding, SAT Solving, and Heuristic
                   Straightening},
  year          = {2025}, eprint = {2507.07951}, archivePrefix = {arXiv},
  url           = {https://arxiv.org/abs/2507.07951}
}

@misc{OEISA006066,
  author       = {{OEIS Foundation Inc.}},
  title        = {A006066: Kobon triangles},
  howpublished = {The On-Line Encyclopedia of Integer Sequences},
  url          = {https://oeis.org/A006066},
  note         = {Published text accessed 2026-08-20; rows retain >=
                  lower-bound notation}
}
```

Optional vector corroboration:
`@misc{SavchukLineOrderKobon14, author={Pavlo Savchuk}, title={LineOrder
Gallery #1 -- Kobon triangle problem},
url={https://zegalur.github.io/line-order/gallery/kobon.html}, note={14-line
solution heading credits Johannes Bader; associated SVG has decimal
rendering coordinates}}`

## 5. Verification standard

- **VERBATIM-SOURCE:** quote only text returned from the directly linked
  page/PDF; preserve source typos and label corrections separately.
- **ATTRIBUTION:** call 53 a Bader construction only with Bader's
  page/PDF/vector artifact plus current OEIS `[Bader]`; cite C-B only for
  the upper-bound refinement, never as the 53 constructor.
- **THEOREM-SCOPE:** label Blanc/BBL quantities as simple affine/projective
  pseudoline triangular-face counts; do not promote them to the broad
  crossed/parallel/multiple-point K_gen theorem without a separate proof.
- **STATUS:** interpret OEIS `a`, `U`, `>=`, and `?` exactly as the table
  convention; published approved text and proposed history must be reported
  separately.
- **COORDINATES:** exact Cartesian/rational coordinates require a source
  that actually publishes them. The Bader PDF/image and LineOrder SVG are
  visual/vector evidence only; record `<no accessible source>` for a
  textual exact-coordinate list.
- **VERIFIED-LOCAL (2026-08-20):** this workspace independently verified
  Maiorana's 54 candidates (SHA-pinned `sol1` hash `83fc2666…`, replayable
  15/15 sweep with exact-rational counts; transcript
  `scratch/kobon/n14/sweep_all_verify.log`; sol1 deep certificate with the
  audited `engine.verify_selection`, `scratch/kobon/n14/maiorana_verify.log`).
  This is *local* verification: it backs the claim `$K_{\rm gen}(14)\ge54$`
  in our documents but does NOT change the approved OEIS text. Only the OEIS
  row itself and any Parpalak/Utkin independent 54 claim remain
  PENDING-VERIFICATION-in-literature.
