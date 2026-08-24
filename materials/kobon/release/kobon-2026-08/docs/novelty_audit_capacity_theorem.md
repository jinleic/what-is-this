# Novelty audit: the crossing-refined capacity theorem

Date: 2026-08-17. Scope: is the inequality

> **(Cap)** For an arrangement of $n$ distinct straight lines in the real affine plane —
> parallels and multiple points allowed, selected triangles required only to have
> pairwise disjoint **open** interiors, other arrangement lines **may cross** a selected
> triangle —
> $$3|S| + C \;\le\; n(n-2) - 2Q - \textstyle\sum_p k_p(k_p-4),$$
> $C$ = number of (selected triangle, outside arrangement line) incidences where the
> line meets the triangle's open interior; $Q$ = number of parallel pairs; $p$ ranges
> over multiple intersection points, $k_p$ = line multiplicity.

proved anywhere in the published literature? Also requested: weaker equivalents in the
face-only convention of type $3T \le n(n-2) - 2Q - \sum_p (k_p-2)^2$.

Internal provenance of (Cap): `math/kobon/report.md` line 309
("With multipoints: $3T \le n(n-2) - 2Q - \sum_p k_p(k_p-4)$"), `math/kobon/engine.py`
lines 624–625 (linewise form $\sigma_r + \gamma_r \le n-2-q_r-\sum_{p\text{ on }r} k_p(k_p-4)$),
and the no-concurrency crossing budget `3·target + crossings ≤ Σ_r (n−2−q_r)`
(engine.py 690–708). **Status flag (per campaign constraint): (Cap) is still under
independent internal audit and is used as a search heuristic, NOT an established theorem.**
The naive per-line charging steps fail in the broad convention (report.md §6 records three
exact-rational counterexamples; audit_shard.py records 256/500 random configs with
concurrency-free opposite-side segment sharing). This audit concerns **novelty only**,
NOT correctness.

## Addendum — 2026-08-22: non-published 2005 prior art located

The original audit was correct only for its stated **published literature**
scope; it missed important public web prior art. Three Japanese essays by
Rasukaru, dated 22 September--8 October 2005 and preserved on the Honma page
later cited by Clément--Bader, use a sequential/per-line multipoint charge,
claim $K(6)=7$, $K(8)=15$, $K(10)=25$, $K(11)=32$, and argue for the
face-specific square penalty
$$3T\le n(n-2)-2Q-\sum_p(k_p-2)^2.$$
Part 3 states the resulting even bound (54 at $n=14$), but explicitly calls
two perturbation reductions “probably” ignorable. It is an informal web
argument, not a refereed theorem or checkable certificate. Primary archives:

- https://web.archive.org/web/20061111125706/http://www10.plala.or.jp/rascalhp/nlines.htm
- https://web.archive.org/web/20221209001801/http://www10.plala.or.jp/rascalhp/nlines2.htm
- https://web.archive.org/web/20221208231339/http://www10.plala.or.jp/rascalhp/nlines3.htm

This changes priority language: the campaign's $n=6,8,10$ results are the
first machine-checkable all-degeneracy certificates known to us, not the first
claims or geometric arguments for the values. It does **not** subsume (Cap):
Rasukaru counts triangular faces, has no crossing-incidence term $C$, gives no
proof in the broad crossed-triangle convention, and does not contain the
$k_p(k_p-4)$ charge. Statements below saying “no source” for the square form
must be read as “no located **published** source”; the web precursor is now
part of the provenance record.

## Addendum — 2026-08-22: the face square penalty is refuted

The square form is not merely unproved. The exact-rational replay
`scratch/kobon/square_penalty_counterexample.py` constructs a generic affine
chart of the classical simplicial arrangement $A(12,1)$ and verifies, with two
independent face predicates,
$$Q=0,\quad N_3=15,\quad N_6=1,\quad F=30,$$
so $3F=90$ while the proposed right-hand side is
$120-(15+16)=89$. Its complete edge census is
$(M,E_2,E_0)=(51,39,0)$, making the equivalent failure
$E_2-E_0=39>2\sum_p(k_p-2)=38$ explicit.

This is an infinite obstruction, not an isolated coordinate accident.
Grünbaum's canonical catalogue, Ars Math. Contemp. 2 (2009), pp. 2 and 4,
defines the simplicial family $R(1)=A(2m,1)$ from the $m$ side lines and $m$
mirror axes of a regular $m$-gon and gives, for $m>3$,
$t_2=m$, $t_3=m(m-1)/2$, and $t_m=1$. For even $m$, the axes split between
opposite vertices and opposite edge midpoints; for odd $m$, each joins a
vertex to the opposite edge midpoint. The proof uses only $m\ge6$. A generic
affine chart violates the square form by at least $(m^2-7m+8)/2$ for every
$m\ge6$, so no fixed coefficient greater than $2/3$ can replace the
coefficient $1$ universally. This refutes the proposed proof path, not the
prior mod-$6$ numerical bound.

## Citation corrections in the assignment (verified)

| task citation | actually is | correct source |
|---|---|---|
| "Blanc 2011 … arXiv:1012.1931" | arXiv:1012.1931 = *Charged particle elliptic flow in p+p collisions* (Zhou, Yan, Dong, Li, Wang, Cai, Sa) — particle physics | Blanc's paper is *The best polynomial bounds for the number of triangles in a simple arrangement of n pseudo-lines*, **arXiv:0801.2845** (pub. Monatsh. Math. 2011); the title "*The maximum number of triangles in arrangements of pseudolines*" belongs to **Roudneff 1996**, JCTB 66:44–74 |
| "BBL 2017 arXiv:1705.07183" | arXiv:1705.07183 = MIMO power normalization (Sadeghi et al.) — EE | BBL = Bartholdi–Blanc–Loisel, **arXiv:0706.0723** (2007; Contemp. Math. 453, 2008) |
| "Bader-Clement 0706.0726" | arXiv:0706.0726 = gravitational-wave radiation reaction (Kidder, Blanchet, Iyer) — gr-qc | Clément–Bader draft is not on arXiv: [cb2007a.pdf](https://www.sop.tik.ee.ethz.ch/publicationListFiles/cb2007a.pdf) (ETH, Dec 21 2007); cached at oeis.org/A006066/a006066.pdf |
| "OEIS A008765" | A008765 = a partition generating-function sequence, unrelated to Kobon | companion bound sequence is **A032765** = ⌊n(n+2)/3⌋ |
| "Grünbaum catalogue arXiv:0904.1244" | arXiv:0904.1244 = *A Spectroscopic Study of Young Stellar Objects in the Serpens Cloud Core and NGC 1333* — astrophysics | Grünbaum, *A catalogue of simplicial arrangements in the real projective plane*, Ars Math. Contemp. 2 (2009), 1–25, [doi:10.26493/1855-3974.88.e12](https://doi.org/10.26493/1855-3974.88.e12); no arXiv identifier used |
| "Cuntz arXiv:1011.1862" | arXiv:1011.1862 = *Quark masses with Nf=2 twisted mass lattice QCD* | Cuntz, *Simplicial arrangements with up to 27 lines*, **arXiv:1108.3000** |

## Source-by-source table

| source | exact statement | hypotheses: pseudoline/straight, simple/general, crossed allowed? | formula | URL |
|---|---|---|---|---|
| **Tamura** (via Clément–Bader draft, and Wikipedia's gloss) | "Saburo Tamura proved that ⌊n(n−2)/3⌋ provides an upper bound on the maximal number [K(n)]" | straight lines; *stated* for general arrangements; **no Q, no multipoint penalty**; proof route is segment charging ("Proof … directly follows from the proof of Lemma 1"), whose local steps are the ones our report §6 refutes under crossed+degenerate conventions | $K(n)\le\lfloor n(n-2)/3\rfloor$ (pointwise integer bound, not a degeneracy-refined inequality) | https://mathworld.wolfram.com/KobonTriangle.html ; https://en.wikipedia.org/wiki/Kobon_triangle_problem |
| **Clément–Bader 2007 draft** (*Tighter Upper Bound for the Number of Kobon Triangles*, ETH, **unpublished draft**) | Thm 1: $K(n)\le \lfloor n(n-2)/3 \rfloor - \mathbf 1_{n\bmod 6\in\{0,2\}}$ (piecewise $B(n)$: $n(n-2)/3$ for $n\equiv3,5$; $n(n-2)/3-1$ for $n\equiv0,2$; $(n-1)^2/3$ for $n\equiv1,4$) | straight lines ("configuration of n straight lines"); degeneracies touched ("if a line intersects an existing point then the number of points decreases by 2 and the number of segments by 3") but **only via perfect-configuration charging**: Lemma 1 step 2 "If a line segment is the side of two triangles then the corresponding line intersects at one of the two endpoints an existing point which belongs to both triangles" and step 3 "every intersection point with more than two corresponding lines is part of at most two pairs of triangles that share a common side" — **exactly the steps report.md §6 refutes with exact-rational witnesses** (n=4 double-covered segment with simple endpoints; n=6 with six shared-side pairs at one triple point). **No Q-term, no Σk_p(k_p−4) term anywhere.** | pointwise integer bound only; **no symbolic degeneracy-refined inequality** | https://www.sop.tik.ee.ethz.ch/publicationListFiles/cb2007a.pdf |
| **BBL** — Bartholdi, Blanc, Loisel 2007/2008 | bound for triangular **faces**: $a_3\le \lfloor n(n-7/3)/3\rfloor$ even affine; projective $p_3$ bounds | **simple arrangements** by definition ("two … intersect transversally into exactly one point, and this one does not belong to any other curve"); pseudoline and straight; triangles are arrangement faces ⇒ no crossing, no parallels, no multipoints by fiat | $a_3(A)\le n(n-7/3)/3$ (even affine) | https://arxiv.org/abs/0706.0723 |
| **Blanc 2011** (*best polynomial bounds*, arXiv:0801.2845) | intro: "Since a bounded segment may not delimit two different triangles, and the number of bounded segments is n(n−2), we have a₃(A) ≤ n(n−2)/3"; Cor 2.4: even affine $a_3\le\lfloor n(n−5/2)/3\rfloor$; Prop 2.3 unused-segment parity argument | **simple** pseudoline arrangements; face count $a_3(A)$ ("polygons delimited by three bounded segments"); crossing/parallels/concurrency excluded by definition. The proof's hinge — "a bounded segment may not delimit two different triangles" — is true for faces of simple arrangements and **false** for disjoint-open-interior selections in the broad convention (report.md 230–239) | $a_3(A)\le n(n-2)/3$; $a_3\le\lfloor n(n-5/2)/3\rfloor$ even | https://arxiv.org/abs/0801.2845 ; https://algebra.dmi.unibas.ch/blanc/articles/bestbound.pdf |
| **Roudneff 1996** (*The maximum number of triangles in arrangements of pseudolines*, JCTB 66:44–74) | projective triangular faces, simple; $p_3^s(n)\le n(n-1)/3$-type bounds | simple projective pseudolines; faces; no degeneracies | $p_3^s$ mod-6 bounds | JCTB 66 (1996) 44–74 (cited in Blanc §1; handbook chap5) |
| **Forge–Ramírez-Alfonsín 1998** (*Straight Line Arrangements in the Real Projective Plane*, DCG 20:155–161) | construction method: from a projective arrangement with the maximum $\tfrac13 k(k-1)$ triangles to one with $2k-2$ lines (doubling); gives $K=2…$ constructions | straight lines, projective; constructions in the simple/perfect regime; cited by Blanc only for the doubling method; **no capacity inequality of any kind** | none (construction) | DCG 20 (1998) 155–161; see https://en.wikipedia.org/wiki/Kobon_triangle_problem |
| **Füredi–Palásti 1984** (*Arrangements of lines with a large number of triangles*, PAMS 92) | lower-bound constructions with ∼n²/3 triangles | straight lines; constructions; no upper-bound inequality | none (lower bound) | PAMS 92 (1984) 561–566 |
| **Savchuk 2025** (arXiv:2507.07951) | SAT/table-encoding constructions; a(11)=32 "confirmed that no optimal solution exists in the 11-line case" by SAT case exhaustion; finds new optima n=23 (161), n=27 (225). The ONLY segment-level statement is inherited: "when n mod 6 ∈ {3,5}, arrangements meeting the upper bound N(n) = (1/3)n(n−2) will have every finite non-overlapping segment as a side of a non-overlapping triangle [5]" | pseudoline tables + straightening; **represents parallels and multiple-line intersections in the encoding** (broad-convention exposure!) but the only upper bounds used are Tamura's ⌊⌋ and C-B; **no analytic degeneracy-refined inequality** | none; upper bounds quoted from [5] | https://arxiv.org/abs/2507.07951 |
| **Parpalak–Utkin 2026** (arXiv:2604.22035; arXiv:2607.29236) | "Given n lines **in general position** in the plane, how many **bounded triangular faces** can the arrangement have?" — constructs the 18·2^t+1 straight series; enumerates triangle-maximal **simple** pseudoline arrangements | **general position / simple**: faces only; no parallels, no multipoints, no crossing | inherits BBL/Blanc bounds | https://arxiv.org/abs/2604.22035 ; https://arxiv.org/abs/2607.29236 |
| **OEIS A006066** (Kobon sequence) | table of constructions and bounds; cites the same Tamura / C-B / BBL bounds; "Savchuk 2025 proved maximal" for a(11)=32; companion bound sequence is A032765 = ⌊n(n+2)/3⌋ | broad problem stated; **no formula with Q or multipoint penalties appears anywhere on the page** | pointwise table only | https://oeis.org/A006066 |
| **Alkauskas 2025** (arXiv:2510.22584, *Triangle unions with maximal number of sides*) | extremal problem for **unions** of triangles (edge/side counts of the union), not arrangement capacity | different problem entirely | none applicable | https://arxiv.org/abs/2510.22584 |
| **Felsner handbook ch.5** (Arrangements of pseudolines) | bibliography lists exactly the simple-convention items above as the known triangle-maximum literature | corroborates: the entire published "maximum triangles" theory lives in the simple/face convention | — | http://www.csun.edu/~ctoth/Handbook/chap5.pdf |
| **Grünbaum 2009 catalogue** | Defines the regular simplicial family $R(1)=A(2m,1)$ and gives its $t$-vector for $m>3$ | straight projective line arrangements; all chambers triangular | $t_2=m$, $t_3=m(m-1)/2$, $t_m=1$ | https://doi.org/10.26493/1855-3974.88.e12 |

Also checked: "BBL prop 2.1" in the assignment ≈ Blanc's Prop 2.3/Cor 2.4 and the
intro segment observation quoted above — all under the simple hypothesis. No
published source contains the surviving broad candidate
$\sum_p k_p(k_p-4)$. The stronger square extension
$\sum_p (k_p-2)^2$ is now refuted both for broad selections and, by the
$A(12,1)$ replay above, for triangular faces. No published source contains
the crossing-incidence term $C$ at all: in the simple face convention
$C\equiv0$ by definition, so the term is structurally invisible to all prior
work.
Grünbaum's published catalogue supplies the canonical $A(2m,1)$ family,
simpliciality, parity-sensitive axis description, and multiplicity census,
but states neither (Cap) nor the square inequality tested here. Cuntz
(arXiv:1108.3000) independently treats simplicial-arrangement classification.

## Verdict

**(Cap) is new mathematical content** — as a *statement*. Specifically:

1. **Not subsumed.** No located published source proves any upper bound for
   pairwise-open-interior-disjoint selections in the broad convention (crossed triangles
   + parallels + concurrencies) that is *refined by arrangement degeneracies*
   ($Q$, $k_p$) *and* *charged for crossed triangles* ($C$). Every published upper-bound
   proof is either (a) in the simple face convention where $Q=0$, all $k_p=2$, $C=0$, so
   the inequality degenerates to the classical $3T\le n(n-2)$ — which is a strict special
   case, not a subsummation — or (b) the Clément–Bader **draft**, which states only
   pointwise mod-6 integer bounds via charging steps that are **demonstrably false**
   in exactly the convention (Cap) targets (quote-level evidence above; refutation
   documented in `math/kobon/report.md` §6 with exact-rational witnesses).

2. **Not contradicted.** Nothing located asserts anything incompatible with (Cap);
   the classical bounds agree with it in their common domain (simple ⇒
   $\Sigma k_p(k_p-4)=0$, $C=0$ ⇒ $3|S|\le n(n-2)$).

3. **Novelty is concentrated exactly where the generality is.** The $.C$ term and the
   multipoint penalty are responses to failure modes (opposite-side segment sharing with
   no multipoint; per-line over-count under concurrency — report.md §6 items 1–3) that
   **cannot occur** in the published simple/face regime. Hence no prior framework even
   encounters, let alone prices, them.

4. **Required caveats before attaching names.**
   (a) (Cap) is **not yet proven**: it is under independent internal audit; the
   campaign constraint instructs treating it as a search heuristic. The hard part of the
   proof is precisely the globally-true-but-locally-false segment accounting that no
   published argument handles.
   (b) The **simple special case** ($3|S|\le n(n-2)$ for selections of interior-disjoint
   triangles in simple arrangements, with sides injecting into used bounded segments) is
   classical/easily-derivable — a referee will likely call that slice folklore. The
   claimable novelty is the **broad-convention inequality with the $C$ and $\sum k_p(k_p-4)$
   refinements and its proof**.
   (c) Tamura's ⌊n(n−2)/3⌋ is *stated* for the broad problem (Wikipedia/MathWorld) but
   the only locateable proofs are segment-charging arguments whose local steps fail in
   the broad convention — i.e., the published record for the broad Tamura bound itself
   has the same gap (already recorded in report.md §2 as the "scope gap"). This
   strengthens, not weakens, the novelty of a correct proof of (Cap).

**Bottom line: new — not subsumed, not contradicted — but publishable only once the
internal audit completes; the claim must be scoped to the degeneracy- and
crossing-refined inequality, whose special cases are classical.**
