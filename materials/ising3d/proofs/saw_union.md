# Exact finite SAW-family union at \(n=36\)

Artifacts: `experiments/e116_saw_union.py`,
`tests/test_saw_union.py`, and `results/bounds/saw_union.json`.

## Scope and notation

Let \(c_r\) be the number of rooted \(r\)-step nearest-neighbour
self-avoiding walks (SAWs) on \(\mathbb Z^3\), and let \(a_r\) be the
number of rooted \(r\)-step SAWs avoiding \(e_1=(1,0,0)\).  Let
\(\mu=\lim_{r\to\infty}c_r^{1/r}\) be the cubic-lattice connective
constant.  All finite family counts below are subsets of \(a_{35}\);
none is claimed to be an all-\(r\) family unless explicitly stated.

Put

\[
 h(x,y,z)=-x+y+z,
 \qquad
 \mathcal M=\{-e_1,e_2,e_3\}.
 \tag{1}
\]

The three letters of \(\mathcal M\) have \(\Delta h=+1\).  The three
defect letters used below, \(e_1,-e_2,-e_3\), have \(\Delta h=-1\).

The external exact input is

\[
 c_{36}=2941370856334701726560670,
 \tag{2}
\]

from Schram--Barkema--Bisseling, Table I.  The producer and clean-room
verifier hash the cached PDF; its SHA-256 is

```text
898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12
```

No numerical estimate of \(K_c\) is used to select, fit, or validate any
certificate below.

---

## 1. The height-slab family \(H\)

Let \(B\) be the set of 11-step SAWs which start at height \(0\), end at
height \(7\), and have every interior height strictly between \(0\) and
\(7\).  Let \(T\) be the set of 13-step SAWs from height \(0\) whose
non-origin vertices have positive height.

**[COMPUTATION]** Exact visited-set enumeration gives

\[
 |B|=729000,
 \qquad |T|=142016661.
 \tag{3}
\]

Define \(H\) by concatenating one block from \(B\), a translated second
block from \(B\), and a translated tail from \(T\).  Their height bands
are respectively

\[
 \{0\}\cup\{1,\ldots,6\}\cup\{7\},
 \quad
 \{7\}\cup\{8,\ldots,13\}\cup\{14\},
 \quad
 \{14\}\cup\{15,16,\ldots\}.
 \tag{4}
\]

**[THEOREM]** This concatenation is an injection into the walks counted by
\(a_{35}\), and

\[
 |H|=|B|^2|T|=75473476338501000000.
 \tag{5}
\]

*Proof.*  The height bands in (4) are disjoint except for their prescribed
splice vertices.  Hence the three pieces cannot intersect elsewhere; the
concatenation is self-avoiding.  Fixed split positions 11 and 22 recover
the three components, proving injectivity.  Every non-origin vertex has
positive height, while \(h(e_1)=-1\), so the path avoids \(e_1\). \(\square\)

---

## 2. The existing defect family \(L\)

The finite defect construction from `proofs/saw_refine.md` consists of:

* all \(\mathcal M\)-words;
* four one-defect corner families, with one \(-e_2\) or \(-e_3\) step
  and the prescribed preceding/following local directions; and
* two fold families with the local pattern
  \((-e_1,e_2,e_1)\) or \((-e_1,e_3,e_1)\).

With

\[
 s(0)=1,\qquad s(j)=2\,3^{j-1}\quad(j\geq1),
\]

write

\[
 G_r=\sum_{k=0}^{r-2}3^k s(r-2-k),
 \qquad
 D_r=\sum_{k=0}^{r-3}3^k s(r-3-k),
 \qquad
 L_r=3^r+4G_r+2D_r.
 \tag{6}
\]

**[THEOREM]** \(a_r\ge L_r\) for every \(r\ge0\).

*Proof.*  A monotone \(\mathcal M\)-word is self-avoiding and avoids
\(e_1\).  For a corner defect, the restricted first suffix letter prevents
the only possible immediate revisit; monotonicity before and after the one
defect prevents any other revisit.  For a fold defect, the displayed
three-step local pattern has the same property.  Each defect has a unique
type and prescribed preceding local pattern, so the seven families are
disjoint.  All paths retain nonpositive first coordinate: in a fold, the
sole \(+e_1\) is preceded by \(-e_1\), so it cannot reach \(e_1\).  The
position sums are exactly (6). \(\square\)

**[COMPUTATION]** Formula (6), independently classified against exhaustive
SAWs through \(r=10\), gives

\[
 L_{35}=644233352324156721,
 \qquad
 L_{29}=741377533262625.
 \tag{7}
\]

---

## 3. The coordinate-monotone family \(C\)

A word is *globally coordinate-monotone* if, for each coordinate, it uses
at most one of the two signs.  Such a word is automatically self-avoiding:
if it revisited a vertex, every monotone coordinate would have to remain
constant on the intervening nonempty segment.

**[LEMMA]** A globally coordinate-monotone word visits \(e_1\) if and only
if its first step is \(+e_1\).

*Proof.*  To reach \((1,0,0)\), the \(x\)-coordinate must have chosen the
positive sign and have taken exactly one \(+e_1\) step.  Since the other
two coordinates are monotone from zero, they can still be zero only if no
\(y\)- or \(z\)-step has occurred.  Thus the prefix is exactly the first
step \(+e_1\).  The converse is immediate. \(\square\)

Define \(C\) as the globally coordinate-monotone 35-step words whose
first step is not \(+e_1\).  By the lemma, \(C\subseteq a_{35}\).

**[COMPUTATION]** A finite automaton with state

\[
 (s_x,s_y,s_z)\in\{\text{unused},-,+\}^3
 \tag{8}
\]

counts \(C\).  Appending a letter is permitted exactly when its coordinate
slot is unused or has the same sign; its slot is then fixed.  The one
initial transition \(+e_1\) is forbidden.  This 27-state automaton gives

\[
 |C|=333543290395947705.
 \tag{9}
\]

The standalone verifier independently obtains (9) by summing, over each
used signed-coordinate pattern, surjective word counts with
inclusion--exclusion and subtracting patterns whose first letter is
\(+e_1\).  It also compares both descriptions with direct SAW enumeration
through length 7.

---

## 4. Exact overlaps

### 4.1 \(H\cap L=\varnothing\)

**[LEMMA] Height partition.**

\[
 H\cap L=\varnothing.
 \tag{10}
\]

*Proof.*  Every \(H\)-walk has \(h=7\) exactly after its first 11 steps.
A monotone \(L\)-word has only \(\Delta h=+1\) letters.  Each of the six
nonmonotone defect families has exactly one defect letter, and that letter
has \(\Delta h=-1\); every other letter is in \(\mathcal M\).  Therefore
for a word in \(L\), if \(d\) is the number of its (at most one) defects
among the first 11 steps,

\[
 h_{11}=11-2d\in\{11,9\}.
 \tag{11}
\]

This cannot equal 7. \(\square\)

**[COMPUTATION] Audited disjoint-union constant.**
\[
|H\cup L|=|H|+|L|
=75473476338501000000+644233352324156721
=76117709690825156721.
\]

This is stricter than the weaker observation that both families have
positive height.  The forced slab endpoint at step 11 is essential.  The
producer and verifier additionally enumerate three direct mini slab
geometries and find the corresponding intersections with the defect family
empty.

### 4.2 \(C\cap L\)

**[THEOREM]**

\[
 |C\cap L|
 =3^{35}+4\,(34\,2^{33})
 =50032713330104219.
 \tag{12}
\]

*Proof.*  All \(3^{35}\) monotone \(\mathcal M\)-words lie in \(C\).
Consider, for example, a \(-e_2\) corner defect.  Coordinate monotonicity
forbids every \(+e_2\) occurrence, so all nondefect letters must be in
\(\{-e_1,+e_3\}\).  The defect is at one of the 34 noninitial positions,
and its prescribed predecessor is fixed.  The remaining 33 positions are
arbitrary over two letters, giving \(34\,2^{33}\).  This applies to both
\(-e_2\) corner types and, by exchanging \(e_2,e_3\), to both \(-e_3\)
corner types.

Each fold word contains both \(-e_1\) and \(+e_1\), so it cannot be
coordinate-monotone.  The original defect families are disjoint, proving
(12). \(\square\)

The producer derives the four terms with a defect-position/previous-letter
automaton.  The verifier independently sums the 34 positions and directly
classifies all SAWs through length 9.

### 4.3 \(C\cap H\)

For a signed-coordinate letter set \(Q\) containing at most one sign for
each coordinate, let

\[
 \begin{aligned}
 B_Q&=\#\{\text{11-step words over }Q:
                 h:0\to7,\ 0<h_j<7\ (1\le j<11)\},\\
 T_Q&=\#\{\text{13-step words over }Q:
                 h_j>0\ (1\le j\le13)\}.
 \end{aligned}
 \tag{13}
\]

A word over \(Q\) is coordinate-monotone and therefore self-avoiding, so
these are exact finite height DPs, not relaxations.  The transition is
simply

\[
 D_{j+1}(q)=\sum_{i\in Q}D_j(q-\Delta h_i),
 \tag{14}
\]

with the stated strict height filters.  The second slab block is a height
translate of the first; hence

\[
 M_Q=B_Q^2T_Q
 \tag{15}
\]

counts members of \(H\) whose steps all lie in \(Q\).

Let \(e(Q)\) be the number of coordinate slots unused by \(Q\).  There
are \(3^3=27\) choices of \(Q\), including the empty set.  Exact
Möbius inversion gives

\[
 |C\cap H|=\sum_Q(-1)^{e(Q)}M_Q.
 \tag{16}
\]

Indeed, a coordinate-monotone walk whose actual used letter set is \(P\)
is counted in every \(Q\supseteq P\).  For each coordinate unused by
\(P\), the contributions from choosing it unused, positive, or negative
are \(-1+1+1=1\); thus its total coefficient is exactly one.

**[COMPUTATION]** The exact DP (13)--(16) gives

\[
 |C\cap H|=157640944278888.
 \tag{17}
\]

For every nonempty \(Q\), the producer independently re-enumerates
\(B_Q\) and \(T_Q\) with a restricted visited-set DFS.  The standalone
verifier uses only that DFS and a per-used-pattern form of Möbius inversion.
Both also compare to direct mini-slab enumerations:

\[
\begin{array}{c|r|r}
(\ell,q,t) & |H_{\ell,q,t}| & |C\cap H_{\ell,q,t}|\\ \hline
(5,3,4)&1{,}994{,}544&24{,}618\\
(5,3,5)&9{,}762{,}768&73{,}872\\
(5,3,6)&39{,}785{,}904&209{,}184.
\end{array}
\tag{18}
\]

In all three cases the direct and Möbius counts agree, and the direct
intersection with the corresponding defect family is empty.

---

## 5. The exact union and a conservative fallback

Since (10) makes the triple overlap empty, inclusion--exclusion using
(5), (7), (9), (12), and (17) gives

\[
\begin{aligned}
 |H\cup L\cup C|
 &=|H|+|L|+|C|-|H\cap L|-|C\cap L|-|C\cap H|\\
 &=76401062626946721319.
\end{aligned}
\tag{19}
\]

**[THEOREM]**

\[
 a_{35}\ge76401062626946721319.
 \tag{20}
\]

*Proof.*  All three families are subsets of \(a_{35}\), and (19) counts
their union exactly. \(\square\)

For auditability, the producer also derives the conservative subfamily
specified in the initial task.  Let \(R_y\) contain words with exactly one
\(-e_2\) and every other letter in \(\{-e_1,+e_3\}\), and define \(R_z\)
symmetrically.  Each has \(35\,2^{34}\) words.  Define

\[
 E=C\setminus(\mathcal M^{35}\cup R_y\cup R_z).
 \tag{21}
\]

**[THEOREM]** \(E\cap L=\varnothing\), and

\[
\begin{aligned}
 |E|&=333543290395947705-3^{35}-2(35\,2^{34})\\
     &=283510542706105118,\\
 |E\cap H|&=|C\cap H|,\\
 |E\setminus H|&=283352901761826230.
\end{aligned}
\tag{22}
\]

*Proof.*  Every coordinate-monotone corner-defect word is in \(R_y\) or
\(R_z\), while a coordinate-monotone fold is impossible; this is the
argument of §4.2 without the predecessor restriction.  Thus (21) is
disjoint from \(L\).  A word in the removed monotone class has no
\(\Delta h=-1\) step, and a word in either removal class has one; its
step-11 height is therefore 11 or 9, never the slab value 7.  Hence the
removed classes do not meet \(H\), proving \(|E\cap H|=|C\cap H|\).
\(\square\)

This reproduces the deliberately weaker fallback

\[
 |H\cup L\cup E|=76401062592586982951,
 \tag{23}
\]

and its corresponding base

\[
 M_E=2941294455272109139577719.
 \tag{24}
\]

It is not the primary result: the exact union (19) is larger by
\(4\,2^{33}=34359738368\).

---

## 6. First-backtrack chaining and the critical-coupling bound

**[THEOREM] First-backtrack inequality.**  For all \(m,n\ge1\),

\[
 c_{m+n}\le c_m\bigl(c_n-a_{n-1}\bigr).
 \tag{25}
\]

*Proof.*  Split an \((m+n)\)-SAW after step \(m\), translating the suffix
to start at the origin.  This embeds all \((m+n)\)-SAWs into valid pairs
\((\omega,\eta)\) of an \(m\)-SAW and an \(n\)-SAW whose translated
concatenation is self-avoiding.

Fix \(\omega\) and let \(\sigma=\omega_m-\omega_{m-1}\).  Choose any
\((n-1)\)-SAW \(\zeta\) avoiding the neighbour \(\sigma\); cubic
symmetry gives exactly \(a_{n-1}\) choices.  Make an \(n\)-walk
\(\eta\) by taking first step \(-\sigma\), then following the translate
\(-\sigma+\zeta\).  Since \(\zeta\) avoids \(\sigma\), this is itself a
self-avoiding \(n\)-walk.  Its translated concatenation with \(\omega\)
hits \(\omega_{m-1}\), so it is invalid.  Dropping the forced first step
recovers \(\zeta\), proving injection.  Therefore at least \(a_{n-1}\)
of the \(c_n\) suffixes are invalid for every \(\omega\), which proves
(25). \(\square\)

**[THEOREM] Chained consequence.**  For every \(n\ge1\),

\[
 \mu^n\le c_n-a_{n-1}.
 \tag{26}
\]

*Proof.*  Put \(X_n=c_n-a_{n-1}\).  Applying (25) with
\(m=(k-1)n\) inductively gives

\[
 c_{kn}\le c_nX_n^{k-1}.
 \tag{27}
\]

Ordinary SAW concatenation gives submultiplicativity, so Fekete's lemma
identifies the limit of \(c_r^{1/r}\), including the subsequence \(kn\),
with \(\mu\).  Taking \(kn\)-th roots of (27) and letting \(k\to\infty\)
yields (26). \(\square\)

**[COMPUTATION]** Both programs exhaustively verify (25) for every split
with \(m+n\le9\) or \(10\), respectively, and verify small instances of
(27).  These finite tests are implementation checks, not replacements for
the all-size injection proof.

Using (20), (26) at \(n=36\), and (2), define

\[
\begin{aligned}
 M
 &=c_{36}-76401062626946721319\\
 &=2941294455272074779839351.
\end{aligned}
\tag{28}
\]

**[THEOREM]**

\[
 \mu^{36}\le M.
 \tag{29}
\]

The inequality is strict relative to the wave-10 base

\[
 c_{36}-|H|=2941295382858363225560670>M.
\tag{30}
\]

The high-temperature SAW-domination argument proved in
`proofs/kc_bounds.md`, §A, gives \(v_c\ge1/\mu\).  Since
\(K_c=\operatorname{atanh}(v_c)\), (29) implies

\[
 \boxed{K_c\ge\operatorname{atanh}(M^{-1/36}).}
 \tag{31}
\]

**[COMPUTATION]** Directed interval arithmetic at 120 decimal digits gives

\[
\operatorname{atanh}(M^{-1/36})
=0.212212058921446585933461933042941687454109711371115920\ldots,
\]

and downward rounding to 40 decimal places gives the certified endpoint

```text
0.2122120589214465859334619330429416874541
```

The standalone verifier independently repeats the enclosure at 150 digits
and verifies

\[
\text{floor}_{40}\le \text{interval lower endpoint}
\le \text{true value}
\le \text{interval upper endpoint}<\text{floor}_{40}+10^{-40}.
\]

---

## 7. Finite selection record and open scope

**[COMPUTATION]** The producer compares the recorded finite lower-family
bases for the external counts \(c_{30},\ldots,c_{36}\) by exact integer
powers only.  The \(n=36\) union certificate is the strongest recorded
instance.  This is not a fit to, or comparison against, a numerical
critical-coupling benchmark.

**[UNRESOLVED]** The stronger union is a finite \(a_{35}\) certificate; no
all-\(r\) strengthening of the defect-family theorem has been proved.

**[UNRESOLVED]** The ratio route \(\mu^2\le c_{32}/c_{30}\) remains open.
A sufficient lower bound would be

\[
 a_{29}\ge c_{30}-\left\lfloor\frac{c_{32}}{c_{30}}\right\rfloor
 =270569905525454674592,
\]

whereas the explicit defect family supplies only
\(L_{29}=741377533262625\).

## Reproduction

From the repository root:

```sh
PYTHONPATH=src .venv/bin/python experiments/e116_saw_union.py
PYTHONPATH=src .venv/bin/python tests/test_saw_union.py
```
