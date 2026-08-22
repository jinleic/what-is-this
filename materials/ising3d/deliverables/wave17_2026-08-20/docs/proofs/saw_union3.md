# Variable-cut and three-block exact SAW union at \(n=36\)

Artifacts: `experiments/e123_saw_union3.py`,
`tests/test_saw_union3.py`, and `results/bounds/saw_union3.json`.

## [LEMMA] Scope and notation

Let \(c_r\) be the number of rooted \(r\)-step nearest-neighbour
self-avoiding walks (SAWs) on \(\mathbb Z^3\), and let \(a_r\) count the
rooted walks avoiding \(e_1=(1,0,0)\).  Put

\[
h(x,y,z)=-x+y+z,\qquad
\mathcal M=\{-e_1,+e_2,+e_3\}.
\tag{1}
\]

[LEMMA] Each letter of \(\mathcal M\) raises \(h\) by one; each of the
other three directed steps lowers it by one.  Every finite family below has
length \(35\) and is proved only to be a subset of \(a_{35}\).

## [EXTERNAL] Exact length-36 input

[EXTERNAL] The sole external count is

\[
c_{36}=2941370856334701726560670,
\tag{2}
\]

from Schram--Barkema--Bisseling, Table I.  Both producer and verifier hash the
cached source before using (2):

```text
898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12
```

[EXTERNAL] No numerical value or estimate of \(K_c\) selects, tunes, or
validates any family in this certificate.

---

## [LEMMA] Strict-record blocks and positive tails

For \(n,q>0\), let \(B(n,q)\) count \(n\)-step SAWs from height zero to
height \(q\), with every interior height in \((0,q)\), and let \(T(n)\)
count \(n\)-step SAWs whose non-origin vertices have positive height.

[LEMMA] A positive-height terminal at \(q\) lies in \(B(n,q)\) exactly when
\(q\) is a strict record height.  Thus one raw positive-height visited-set DFS
through depth 13 gives all required block and tail counts.

## [COMPUTATION] Exact raw profile used here

[COMPUTATION] The producer's packed-coordinate DFS and the verifier's
independent iterative tuple-coordinate DFS give the following selected values:

\[
\begin{array}{c|rrr}
n&T(n)& &\\ \hline
11&6\,827\,355&&\\
12&29\,971\,017&&\\
13&142\,016\,661&&
\end{array}
\qquad
\begin{array}{c|rrrrr}
n\backslash q&3&5&7&9&11\\ \hline
7&432&2916&2187&&\\
9&1566&34668&43740&19683&\\
11&5832&411156&729000&551124&177147\\
13&20952&4815450&11544444&11573604&6377292
\end{array}
\tag{3}
\]

[COMPUTATION] The complete strict-record profile through depth 13 is stored in
`data.positive_height_profile` and is recomputed from raw lattice steps by the
standalone verifier.

---

## [THEOREM] Height-slab concatenation

For a two-block schedule \(s=(\ell_1,q_1,\ell_2,q_2,t)\), concatenate a
block in \(B(\ell_1,q_1)\), its spatial translate in
\(B(\ell_2,q_2)\), and a translated tail in \(T(t)\).  For a three-block
schedule \(s=(\ell_1,q_1,\ell_2,q_2,\ell_3,q_3,t)\), concatenate the three
blocks and tail analogously.

[THEOREM] If the displayed lengths sum to 35, then each such concatenation is
self-avoiding and avoids \(e_1\), and its cardinality is respectively

\[
\prod_{i=1}^{2} B(\ell_i,q_i)T(t),
\qquad
\prod_{i=1}^{3} B(\ell_i,q_i)T(t).
\tag{4}
\]

[LEMMA] The proof is the strict-band argument: the first block occupies
heights strictly between its two splice heights; every next block occupies the
next open height band; the tail lies strictly above the final splice.  These
bands meet only at prescribed splice vertices, so components cannot otherwise
intersect.  The fixed cut positions recover all components, proving
injectivity.  Every non-origin vertex has positive height while
\(h(e_1)=-1\), proving \(e_1\)-avoidance.

---

## [THEOREM] The selected variable-cut and three-block grids

Let the old 25-member grid be
\(H_{q_1,q_2}=(11,q_1,11,q_2,13)\), where
\(q_1,q_2\in\{3,5,7,9,11\}\).  Its exact union with the pre-existing
\(L,C,F_2,Q_2\) families is used only as a certified base below.

[THEOREM] The new variable-cut grid is

\[
V_q=(13,q,11,3,11),\qquad q\in\{5,7,9,11\}.
\tag{5}
\]

[THEOREM] The new three-block grid is

\[
J_{q_1,q_3}=(7,q_1,7,3,9,q_3,12),
\quad q_1\in\{3,5,7\},\quad q_3\in\{3,5,7,9\}.
\tag{6}
\]

[THEOREM] All rises in (5)--(6) are drawn from the declared admissible odd
set \(\{3,5,7,9,11\}\).  The new cuts are genuinely different from the old
cuts \(11,22\): (5) cuts at 13 and 24, while (6) cuts at 7, 14, and 23.

## [COMPUTATION] Complete pairwise grid-disjointness certificate

For a schedule \(s\), let \(A_s(j)\) be its exact set of allowed heights at
step \(j\), restricted to \(\{1,\ldots,35\}\); this restriction loses no
walk because a 35-step nearest-neighbour walk cannot have larger height.

[LEMMA] If \(A_s(j)\cap A_{s'}(j)=\varnothing\) for one \(j\), then the
underlying schedule families are disjoint.

[COMPUTATION] There are 41 schedules in the combined old/new grid, hence
\(\binom{41}{2}=820\) pair terms.  The artifact literally stores all 820 rows
`(left,right,step,left_allowed,right_allowed)` in
`data.schedule_disjointness.pairwise_separator_rows`; every row has disjoint
finite allowed-height sets.

[THEOREM] Therefore every pairwise intersection of distinct height-slab
families is empty.  Every triple or higher height-slab intersection is also
empty because it is contained in one of the listed empty pair intersections.
No hidden schedule-schedule inclusion--exclusion term remains.

---

## [THEOREM] Exact overlaps of every new height-slab row

Let \(C\) be the globally coordinate-monotone family and \(L\) the seven
old one-defect regular languages.  For a signed-coordinate orthant \(Q\), let
\(B_{\ell,q}(Q)\) and \(T_t(Q)\) be exact finite word-DP counts over that
orthant.  Then the signed-orthant Möbius formula is

\[
|C\cap H_s|=
 \sum_Q(-1)^{e(Q)}\left(\prod_{(\ell,q)\in s}B_{\ell,q}(Q)\right)T_t(Q),
\tag{7}
\]

where \(e(Q)\) is the number of unused coordinate slots and the product ranges
over the block components of \(s\).

[THEOREM] Equation (7) is exact because inclusion--exclusion over the 27
signed-coordinate slot patterns counts every globally coordinate-monotone
word once.  The positive first height automatically excludes an initial
\(+e_1\) word.

[COMPUTATION] The exact finite DPs give the complete variable-cut table.  In
every listed row, \(|L\cap V_q|=|C\cap L\cap V_q|=|F_2\cap V_q|=|Q_2\cap
V_q|=0\); the net column is therefore \(H-C\cap H\).

\[
\begin{array}{c|r|r|r}
q&|V_q|&|C\cap V_q|&\text{net}\\ \hline
5&191737419653862000&676984409640&191736742669452360\\
7&459666677859495840&1624493301552&459665053366194288\\
9&460827745497433440&1673651036592&460826071846396848\\
11&253925492416953120&886041378936&253924606375574184\\ \hline
\text{sum}&1366157335427744400&4861170126720&1366152474257617680
\end{array}
\tag{8}
\]

[COMPUTATION] The exact finite DPs give the complete three-block table.  In
every listed row, \(|L\cap J|=|C\cap L\cap J|=|F_2\cap J|=|Q_2\cap J|=0\).

\[
\begin{array}{cc|r|r|r}
q_1&q_3&|J_{q_1,q_3}|&|C\cap J_{q_1,q_3}|&\text{net}\\ \hline
3&3&8759125145968128&36392828076&8759088753140052\\
3&5&193908908403846144&582248679264&193908326155166880\\
3&7&244651426490833920&727799427900&244650698691406020\\
3&9&110093141920875264&291118630428&110092850802244836\\
5&3&59124094735284864&218343254724&59123876392030140\\
5&5&1308885131725961472&3493437253920&1308881638288707552\\
5&7&1651397128813128960&4366779456420&1651392762033672540\\
5&9&743128707965908032&1746710075628&743126961255832404\\
7&3&44343071051463648&145559885580&44342925491578068\\
7&5&981663848794471104&2328949043424&981661519845427680\\
7&7&1238547846609846720&2911183459380&1238544935426387340\\
7&9&557346530974431024&1164473100648&557345366501330376\\ \hline
\text{sum}&&7141848962632019280&18012995095392&7141830949636923888
\end{array}
\tag{9}
\]

[THEOREM] From (8)--(9), the exact contribution of all 16 new grid rows,
relative to the complete old union, is

\[
G_{\rm new}=8507983423894541568.
\tag{10}
\]

[COMPUTATION] The artifact records every row of (8)--(9), every zero pair or
triple term, and every signed-orthant summand of (7), rather than only the
aggregates.

---

## [THEOREM] Two new three-defect languages

Let

\[
F_y=(-e_1,+e_2,+e_1),\qquad
Q_y=(+e_2,+e_3,-e_2).
\tag{11}
\]

Define \(F_3\) to be all words

\[
uF_yvF_ywF_yx,
\tag{12}
\]

where \(u\in\mathcal M^*\), each internal segment \(v,w\in\mathcal M^+\)
begins in \(\{+e_2,+e_3\}\), and \(x\) is empty or begins in that same
set.  Define \(Q_3\) by replacing \(F_y\) by \(Q_y\) and the required first
letters of the internal/nonempty final segments by \(\{-e_1,+e_3\}\).

[LEMMA] Contracting each \(F_y\) to a \(+e_2\) skeleton edge leaves an
\(\mathcal M\)-monotone skeleton.  The required first post-gadget advance
excludes the unique immediate return to the gadget side vertex; coordinate
monotonicity excludes every later collision and separates side vertices of
distinct gadgets.  Thus every \(F_3\) word is self-avoiding.

[LEMMA] In coordinates \((a,b,c)=(-x,y,z)\), an \(F_y\) gadget never makes
\(a\) negative, so \(F_3\) avoids \(e_1\).  The identical contracted-skeleton
argument proves \(Q_3\) self-avoiding, and \(Q_3\) has no \(+e_1\) letter, so
it also avoids \(e_1\).

[LEMMA] For \(k\) prescribed gadgets, put \(N=35-3k\).  The exact segment sum is

\[
\sum_{\substack{n_0\ge0,\ n_1,\ldots,n_{k-1}\ge1,\ n_k\ge0\\
 n_0+\cdots+n_k=N}}
3^{n_0}\prod_{i=1}^{k-1}\bigl(2\,3^{n_i-1}\bigr)
\begin{cases}
1,&n_k=0,\\
2\,3^{n_k-1},&n_k>0.
\end{cases}
\tag{13}
\]

[COMPUTATION] Equation (13) at \(k=3\) gives

\[
|F_3|=|Q_3|=2325336517026900.
\tag{14}
\]

[COMPUTATION] The verifier directly generates every short word for both
families through the nontrivial small-length range, checks self-avoidance and
\(e_1\)-avoidance vertex by vertex, and independently checks (13).

[THEOREM] \(F_3\) is disjoint from \(L,C,F_2,Q_2,Q_3\): it has exactly three
\(+e_1\) defect letters, both \(x\)-signs, while the named competing families
have respectively at most one, at most one sign, two, zero, and zero
\(+e_1\) letters.  Symmetrically, \(Q_3\) is disjoint from
\(L,C,Q_2,F_2,F_3\) by its exactly three \(-e_2\) defect letters, both
\(y\)-signs, and absence of \(+e_1\).

[COMPUTATION] The exact language DP against all 41 pairwise-disjoint height
schedules has precisely the following nonzero rows, identically for \(F_3\)
and \(Q_3\):

\[
\begin{array}{c|r}
H_{q_1,q_2}&|F_3\cap H_{q_1,q_2}|=|Q_3\cap H_{q_1,q_2}|\\ \hline
H_{7,9}&31632108085872\\
H_{7,11}&51967034712504\\
H_{9,7}&31632108085872\\
H_{9,9}&424397450152116\\
H_{9,11}&166068567450828\\
H_{11,7}&51967034712504\\
H_{11,9}&166068567450828\\
H_{11,11}&9790890598008
\end{array}
\tag{15}
\]

[COMPUTATION] All 16 intersections of each of \(F_3,Q_3\) with the new
variable-cut and three-block rows are exactly zero and are retained in the
artifact.  The total of (15) is

\[
933523761248532,
\qquad
|F_3\setminus\textstyle\bigcup H_s|
=|Q_3\setminus\textstyle\bigcup H_s|
=1391812755778368.
\tag{16}
\]

---

## [THEOREM] Exact inclusion--exclusion certificate

[COMPUTATION] The independently recomputed inherited base is

\[
U_{\rm old}=499330428831189067251.
\tag{17}
\]

[THEOREM] The old base already contains \(L,C\), its 25 disjoint height
slabs, and the disjoint two-gadget languages \(F_2,Q_2\).  Equations
(8)--(10), (15)--(16), the disjointness theorem, and the stated family
separations give the complete exact expression

\[
\begin{aligned}
U={}&U_{\rm old}
 +\sum_{s\in\{V,J\}}\bigl(|H_s|-|C\cap H_s|-|L\cap H_s|
       +|C\cap L\cap H_s|-|F_2\cap H_s|-|Q_2\cap H_s|\bigr)\\
&+|F_3|-\sum_{s\in\{H,V,J\}}|F_3\cap H_s|
 +|Q_3|-\sum_{s\in\{H,V,J\}}|Q_3\cap H_s|\\
={}&507841195880595165555.
\end{aligned}
\tag{18}
\]

[THEOREM] Every constituent of (18) is a subset of \(a_{35}\), so

\[
\boxed{a_{35}\ge507841195880595165555.}
\tag{19}
\]

[COMPUTATION] The standalone verifier has omission traps: dropping the new
\(C\cap H\) terms changes (18) by
\(4861170126720+18012995095392\), and dropping the two three-gadget grid
overlap totals changes it by \(2\cdot933523761248532\).

---

## [THEOREM] Chaining and the Ising lower endpoint

[THEOREM] The first-backtrack injection proved in `proofs/saw_union.md` and
used in `proofs/saw_union2.md` gives, for every \(m,n\ge1\),

\[
c_{m+n}\le c_m\bigl(c_n-a_{n-1}\bigr),
\qquad
\mu^n\le c_n-a_{n-1}.
\tag{20}
\]

[THEOREM] Substitution of (2) and (19) at \(n=36\) gives

\[
M=c_{36}-a_{35}^{\rm lower}
=2940863015138821131395115,
\qquad \mu^{36}\le M.
\tag{21}
\]

[THEOREM] The SAW correlation-domination proof in `proofs/kc_bounds.md`,
section A, gives \(v_c\ge1/\mu\).  Monotonicity of \(\operatorname{atanh}\)
therefore yields

\[
K_c\ge\operatorname{atanh}\bigl(M^{-1/36}\bigr).
\tag{22}
\]

[COMPUTATION] At 120 decimal digits of directed interval arithmetic, the
producer enclosure is

```text
[0.21221294985163282536308213395043835554375109820294088436712143642598023712200841018763813390035329507018572587440834334308152,
 0.21221294985163282536308213395043835554375109820294088436712143642598023712200841018763813390035329507018572587440834334342037]
```

[COMPUTATION] The independent verifier repeats the calculation at 150 digits,
checks directed-rounding containment, and certifies

\[
\boxed{K_c\ge0.2122129498516328253630821339504383555437.}
\tag{23}
\]

## [COMPUTATION] Resource record

[COMPUTATION] The producer observed 43.78457750000234 seconds wall time and a
29,720,576-byte peak RSS before artifact serialization.
[UNRESOLVED] The artifact separately records the unlaunched preflight wall of
300 seconds and RSS wall of 8 GiB for the depth-13 profile, 41 grid rows, 820
separators, and 82 three-gadget schedule intersections.

## [UNRESOLVED] Scope not established

[UNRESOLVED] The selected grids are explicit exact disjoint subfamilies, not a
global optimization over all 35-step block/tail schedules.

[UNRESOLVED] This finite lower bound does not provide an all-size recursion
stronger than the existing \(L_r\) family, does not prove a Kesten-ratio bound,
and does not solve the three-dimensional Ising model.
