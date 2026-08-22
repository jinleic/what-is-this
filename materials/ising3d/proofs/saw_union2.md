# Exact enlarged finite SAW union at \(n=36\)

Artifacts: `experiments/e120_saw_union2.py`,
`tests/test_saw_union2.py`, and `results/bounds/saw_union2.json`.

## [LEMMA] Scope and notation

Let \(c_r\) be the number of rooted \(r\)-step nearest-neighbour
self-avoiding walks (SAWs) on \(\mathbb Z^3\), and let \(a_r\) count the
rooted \(r\)-step SAWs that avoid \(e_1=(1,0,0)\).  Put

\[
h(x,y,z)=-x+y+z,\qquad
\mathcal M=\{-e_1,+e_2,+e_3\}.
\tag{1}
\]

The three letters in \(\mathcal M\) have \(\Delta h=+1\), while
\(+e_1,-e_2,-e_3\) have \(\Delta h=-1\).  Every lower-family statement
below concerns only the finite length \(r=35\), unless it is explicitly
labelled otherwise.

## [EXTERNAL] The exact input count

The certificate uses

\[
c_{36}=2941370856334701726560670,
\tag{2}
\]

from Schram--Barkema--Bisseling, Table I.  The producer and standalone
verifier hash the cached PDF before consuming (2):

```text
898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12
```

No numerical estimate of \(K_c\) is used in choosing, fitting, or checking
any family below.

---

## [LEMMA] One positive-height DFS produces all required blocks

For integers \(n,q>0\), let \(B(n,q)\) count the \(n\)-step SAWs that
start at height \(0\), end at height \(q\), and have every interior height
strictly in \((0,q)\).  Let \(T(n)\) count the \(n\)-step SAWs whose
non-origin vertices have positive height.

A positive-height walk ending at height \(q\) is in \(B(n,q)\) if and only
if \(q\) is its strict record height.  Indeed, positivity gives the lower
interior inequality, and strict-record status gives the upper one; the
converse is immediate.  Therefore a visited-set enumeration of all
positive-height SAWs through depth 13 simultaneously computes every
\(B(n,q)\) and \(T(n)\) in the stated range.

## [COMPUTATION] The exact profile

The direct packed-coordinate DFS gives

\[
(T(0),\ldots,T(13))=
(1,3,9,45,171,837,3411,16425,69525,331947,1436643,
6827355,29971017,142016661).
\tag{3}
\]

The block entries used at length 11 are

\[
\begin{array}{c|rrrrr}
q&3&5&7&9&11\\ \hline
B(11,q)&5832&411156&729000&551124&177147.
\end{array}
\tag{4}
\]

The complete exact profile through depth 13 is stored in
`data.positive_height_profile` of the result JSON.  Its standalone verifier
uses an iterative tuple-coordinate DFS rather than the producer's recursive
packed-coordinate DFS.

---

## [COMPUTATION] Full bounded exact split search

The result JSON literally stores every row, rather than only an optimizer:

* `data.split_search.two_block_table` has all **392** schedules
  \((\ell_1,q_1,\ell_2,q_2,t)\) with
  \(9\le\ell_1,\ell_2,t\le13\),
  \(\ell_1+\ell_2+t=35\), and every nonzero profile block value;
* `data.split_search.three_block_table` has all **1651** schedules with
  \(5\le\ell_1,\ell_2,\ell_3\le9\), \(8\le t\le13\),
  \(\ell_1+\ell_2+\ell_3+t=35\), and every nonzero profile block value.

For a two-block row, the exact product is

\[
B(\ell_1,q_1)B(\ell_2,q_2)T(t).
\tag{5}
\]

The unique maximum in the declared two-block table is the pre-existing
schedule

\[
(\ell_1,q_1,\ell_2,q_2,t)=(11,7,11,7,13),
\qquad
729000^2\,142016661=75473476338501000000.
\tag{6}
\]

The three-block table includes asymmetric length and rise variants.  Its
maximum product is

\[
17606452560245513280,
\tag{7}
\]

attained by the six permutations of the displayed block multiset in the
artifact, for example \((5,5,8,6,9,7,13)\).  Thus no selected
three-block row exceeds (6) in this exact bounded search.

## [UNRESOLVED] Search scope

Equations (5)--(7) are maxima only over the two explicitly stored finite
profile-feasible domains.  They are not claimed to optimize arbitrary
35-step block/tail schedules, whose complete exact enumeration was not
attempted.

---

## [THEOREM] A 25-member pairwise-disjoint height-slab grid

For \(q_1,q_2\in\{3,5,7,9,11\}\), define \(H_{q_1,q_2}\) by concatenating

1. a block in \(B(11,q_1)\),
2. a spatial translate of a block in \(B(11,q_2)\), and
3. a spatial translate of a tail in \(T(13)\).

The three respective height bands are disjoint except at their splice
vertices.  Hence the concatenation is self-avoiding, and fixed split points
11 and 22 recover its three components.  Since every non-origin vertex has
positive height and \(h(e_1)=-1\), it avoids \(e_1\).  Consequently

\[
|H_{q_1,q_2}|=B(11,q_1)B(11,q_2)T(13)
\quad\hbox{and}\quad H_{q_1,q_2}\subseteq a_{35}.
\tag{8}
\]

## [LEMMA] The 25 slab schedules are disjoint

A member of \(H_{q_1,q_2}\) has

\[
h_{11}=q_1,\qquad h_{22}=q_1+q_2.
\tag{9}
\]

If two index pairs differ in \(q_1\), their first values in (9) differ; if
they have the same \(q_1\) and differ in \(q_2\), their second values
differ.  Thus the 25 sets are pairwise disjoint.  The old wave-11 height
family is exactly \(H_{7,7}\), so this grid contains it.

## [THEOREM] Exact \(C\)-overlaps of the grid

Let \(C\) be the coordinate-monotone family from `proofs/saw_union.md`.
For a signed-coordinate orthant \(Q\) (at most one sign on each axis), let
\(B_i(Q)\) be the restricted block word count for the \(i\)-th block and
let \(T(Q)\) be the restricted tail word count.  A word over \(Q\) is
coordinate-monotone and therefore self-avoiding.  If \(e(Q)\) is the
number of unused coordinate slots, Möbius inversion gives

\[
|C\cap H_{q_1,q_2}|
 =\sum_Q(-1)^{e(Q)}B_1(Q)B_2(Q)T(Q).
\tag{10}
\]

The positive-height condition excludes an initial \(+e_1\), so (10) counts
the required \(e_1\)-avoiding coordinate-monotone words without any
additional subtraction.  The 26 nonempty exact DP rows for every schedule
are stored in `data.height_grid.cm_mobius_tables`.

## [THEOREM] Exact \(L\)- and triple-overlaps of the grid

Let \(L\) be the seven-family one-defect construction from
`proofs/saw_union.md`, with

\[
|L|=644233352324156721,
\qquad
|C\cap L|=50032713330104219.
\tag{11}
\]

For each of the seven disjoint regular languages defining \(L\), a finite
height DP records the step position, current height, and (for the triple)
the three coordinate-sign slots.  Since every \(L\)-word is already a
self-avoiding \(e_1\)-avoiding walk, satisfying the slab height constraints
is equivalent to membership in the corresponding slab set.  This gives all
pairwise and triple terms exactly.  The only nonzero \(L\)-overlap rows are

\[
\begin{array}{c|r|r}
(q_1,q_2)&|L\cap H_{q_1,q_2}|&|C\cap L\cap H_{q_1,q_2}|\\ \hline
(9,11)&121063985671653612&240518168576\\
(11,9)&121063985671653612&240518168576\\
(11,11)&248922378702430641&50031923056121755.
\end{array}
\tag{12}
\]

All other 22 rows vanish.  The exact total terms are

\[
\begin{aligned}
\sum_s|C\cap H_s|&=51027099125602155,\\
\sum_s|L\cap H_s|&=491050350045737865,\\
\sum_s|C\cap L\cap H_s|&=50032404092458907.
\end{aligned}
\tag{13}
\]

Here and below \(s\) ranges over the 25 ordered pairs.  Therefore the
exact net contribution of the entire disjoint grid, relative to \(L\cup C\),
is

\[
\begin{aligned}
G
&=\sum_s\bigl(|H_s|-|C\cap H_s|-|L\cap H_s|
                  +|C\cap L\cap H_s|\bigr)\\
&=498390727964942607428.
\end{aligned}
\tag{14}
\]

---

## [THEOREM] Broad one-defect coordinate-monotone candidates are redundant

Define \(R_y\) to consist of words with exactly one \(-e_2\), every other
letter in \(\{-e_1,+e_3\}\), and define \(R_z\) symmetrically.  These
are self-avoiding and avoid \(e_1\), because each uses one sign per
coordinate and has nonpositive physical \(x\)-coordinate.  Their exact
sizes are

\[
|R_y|=|R_z|=35\,2^{34}=601295421440.
\tag{15}
\]

Both are already subsets of \(C\), so they cannot enlarge a union that
contains \(C\).  Their exact intersections are

\[
|R_y\cap C|=|R_z\cap C|=601295421440,
\quad
|R_y\cap L|=|R_z\cap L|=34\,2^{33}=292057776128,
\tag{16}
\]

and each has exact total grid intersection

\[
\sum_s|R_y\cap H_s|=\sum_s|R_z\cap H_s|=429496729600.
\tag{17}
\]

Those grid intersections are also their \(L\)-triples: positivity rules out
a defect at the first step, and every remaining word lies in the applicable
corner family of \(L\).  The per-schedule values are stored in
`data.candidate_audit.one_defect_coordinate_monotone`.

## [THEOREM] The tested single-fold candidate is already \(D_1\subset L\)

The candidate language

\[
u(-e_1,+e_2,+e_1)w,
\tag{18}
\]

with \(u\in\mathcal M^*\) and \(w\) either empty or beginning in
\(\{+e_2,+e_3\}\), is exactly the existing \(D_1\) language: the
post-fold restriction prevents its unique immediate return.  It is therefore
contained in \(L\), has no coordinate-monotone members, and adds no new
union mass.  Its exact total grid intersection is

\[
\sum_s|D_1\cap H_s|=31501343210481297.
\tag{19}
\]

## [UNRESOLVED] Other single-fold designs

Equation (18) audits the concrete safe one-fold candidate used here; it does
not classify all conceivable one-fold languages outside the fixed
\(\mathcal M\) background.

---

## [THEOREM] Two disjoint prescribed two-defect languages

Let \(N=35-6=29\).  Define a fold gadget and a corner gadget by

\[
F=(-e_1,+e_2,+e_1),
\qquad
Q=(+e_2,+e_3,-e_2).
\tag{20}
\]

Define \(\mathcal F\) as all words

\[
uFvFw,
\tag{21}
\]

where \(u\in\mathcal M^*\), \(v\in\mathcal M^+\) begins in
\(\{+e_2,+e_3\}\), and \(w\) is empty or begins in that same two-letter
set.  Define \(\mathcal Q\) as all words \(uQvQw\), with the analogous
restriction that \(v\) and any nonempty \(w\) begin in
\(\{-e_1,+e_3\}\).

### [LEMMA] Self-avoidance and \(e_1\)-avoidance

Use coordinates \((a,b,c)=(-x,y,z)\).  Contracting each \(F\) in (21)
to one \(+b\) step leaves an \(\mathcal M\)-monotone skeleton.  The two
interior vertices of a gadget are the side vertices of that skeleton edge.
A later skeleton vertex can hit the second side vertex only by taking an
immediate \(+a\) step after the gadget; the required first post-gadget
advance rules this out, and monotonicity then rules it out forever.  Distinct
gadgets have distinct side vertices because the required intervening advance
separates their contracted skeleton edges.  Earlier skeleton vertices cannot
meet a side vertex by coordinate monotonicity.  Thus \(\mathcal F\) is
self-avoiding.  Its \(a\)-coordinate never becomes negative, so its physical
\(x\)-coordinate never becomes positive and it avoids \(e_1\).

The identical argument for \(Q\) contracts it to a \(+c\) skeleton edge;
the forbidden immediate return is now \(+b\).  Hence \(\mathcal Q\) is
self-avoiding.  It has no \(+e_1\) letter and therefore avoids \(e_1\).

### [THEOREM] Exact counts and mutual disjointness

For each class, the free lengths in (21) sum to \(N\).  If the final free
segment is empty, there are \(2N3^{N-1}\) words.  If it is nonempty, there
are \(2N(N-1)3^{N-2}\) words.  Therefore

\[
|\mathcal F|=|\mathcal Q|
=2N3^{N-1}+2N(N-1)3^{N-2}
=13710824278006626.
\tag{22}
\]

Each \(\mathcal F\)-word has two \(+e_1\) defect letters, so it is
disjoint from \(L\), whose nonmonotone families have one defect.  It has
both \(x\)-signs, so it is disjoint from \(C\).  Each \(\mathcal Q\)-word
has two \(-e_2\) defect letters and both \(y\)-signs, giving the analogous
disjointness from \(L\) and \(C\).  Finally \(\mathcal F\cap\mathcal Q\)
is empty because every \(\mathcal F\)-word has \(+e_1\), whereas no
\(\mathcal Q\)-word does.

### [COMPUTATION] Exact grid intersections

An exact finite word DP over the three free-segment lengths gives

\[
\sum_s|\mathcal F\cap H_s|
=\sum_s|\mathcal Q\cap H_s|
=7732355849776818.
\tag{23}
\]

The only nonzero schedule rows for either class are

\[
\begin{array}{c|r}
(q_1,q_2)&|\mathcal F\cap H_{q_1,q_2}|=|\mathcal Q\cap H_{q_1,q_2}|\\ \hline
(7,11)&183014339639688\\
(9,9)&1494617107057452\\
(9,11)&2455442390165814\\
(11,7)&183014339639688\\
(11,9)&2455442390165814\\
(11,11)&960825283108362.
\end{array}
\tag{24}
\]

All mixed triples containing one of these two-defect classes vanish by the
pairwise disjointness just proved.  Thus each class contributes exactly

\[
13710824278006626-7732355849776818=5978468428229808
\tag{25}
\]

new walks beyond \(L\cup C\cup\bigcup_sH_s\).

---

## [THEOREM] Exact inclusion--exclusion certificate

The established base union has

\[
|L\cup C|=|L|+|C|-|L\cap C|=927743929390000207.
\tag{26}
\]

The grid is internally disjoint, all its pairwise/triple intersections with
\(L,C\) are (10)--(13), and the two new languages have only the grid
overlaps (23).  Therefore exact inclusion--exclusion gives

\[
\begin{aligned}
U
={}&|L|+|C|-|L\cap C|\\
&+\sum_s\bigl(|H_s|-|C\cap H_s|-|L\cap H_s|
                 +|C\cap L\cap H_s|\bigr)\\
&+|\mathcal F|-\sum_s|\mathcal F\cap H_s|
 +|\mathcal Q|-\sum_s|\mathcal Q\cap H_s|\\
={}&499330428831189067251.
\end{aligned}
\tag{27}
\]

Every constituent is a subset of \(a_{35}\), hence

\[
\boxed{a_{35}\ge499330428831189067251.}
\tag{28}
\]

The standalone verifier contains traps showing that omitting any of the
\(C\cap H\), \(L\cap H\), \(C\cap L\cap H\), or gadget-grid overlap
classes changes (27).  In particular, omitting the triple term changes the
answer by \(50032404092458907\), and omitting the two gadget-grid terms
changes it by \(2\cdot7732355849776818\).

---

## [THEOREM] Chaining to a connective-constant bound

The first-backtrack injection proved in `proofs/saw_union.md` gives, for all
\(m,n\ge1\),

\[
c_{m+n}\le c_m\bigl(c_n-a_{n-1}\bigr).
\tag{29}
\]

Iterating (29) at fixed \(n\) and applying Fekete's lemma along multiples
of \(n\) yields

\[
\mu^n\le c_n-a_{n-1}.
\tag{30}
\]

The producer independently checks (29) by exact SAW backtracking for every
split with \(m+n\le8\); that finite check is an implementation audit, not a
substitute for the all-size injection proof.

Substituting (2) and (28) into (30) at \(n=36\) gives

\[
\begin{aligned}
M&=c_{36}-499330428831189067251\\
 &=2940871525905870537493419,\\
\mu^{36}&\le M.
\end{aligned}
\tag{31}
\]

This is strictly stronger than the wave-11 base
\(2941294455272074779839351\).

## [THEOREM] The resulting rigorous Ising lower endpoint

The SAW correlation-domination proof in `proofs/kc_bounds.md`, §A gives
\(v_c\ge1/\mu\).  By monotonicity of \(\operatorname{atanh}\), (31) yields

\[
K_c\ge\operatorname{atanh}(M^{-1/36}).
\tag{32}
\]

## [COMPUTATION] Directed rounding

At 120 decimal digits of directed interval arithmetic, the producer obtains
an interval for (32) contained in the 40-place bin beginning at

```text
0.2122129322754723621039731646196408804906
```

and therefore certifies the downward floor

\[
\boxed{K_c\ge0.2122129322754723621039731646196408804906.}
\tag{33}
\]

The standalone verifier recomputes the enclosure at 150 decimal digits and
checks

\[
\operatorname{floor}_{40}\le\text{interval lower}
\le\operatorname{atanh}(M^{-1/36})
\le\text{interval upper}<\operatorname{floor}_{40}+10^{-40}.
\tag{34}
\]

## [UNRESOLVED] Remaining scope

The finite union in (28) does not construct an all-\(r\) family strengthening
\(L_r\), does not settle a Kesten-ratio bound, and does not constitute an
exact solution of the three-dimensional Ising model.
