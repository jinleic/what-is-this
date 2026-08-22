# Ladder two-generator growth: exact L=8 patch and boundary-transfer obstruction

## 1. Outcome and scope

For the open two-leg ladder with \(L\) rungs, write

\[
 A_L=\sum_v X_v,
 \qquad
 B_L=\sum_{(u,v)\in E_L} Z_uZ_v,
\]

and let \(D_k\) be the dimension of the generator-left-nested-word filtration through
bracket depth \(k\).  All brackets below are divided by two in the Pauli
\(Q_{(a\mid b)}=X^aZ^b\) convention; this changes nonzero scalars but not any rank.

**[COMPUTATION]** There is an explicit family of 256 generator-left-nested words of
maximum depth 16 on the eight-rung ladder whose last-three-rung Pauli projection has
rank 256 modulo both

\[
 p_1=2147483647,\qquad p_2=2147483629.
\]

Therefore the 256 full operators are linearly independent over \(\mathbb Q\), and

\[
 \boxed{D_{16}\ge256=2^8.}
\]

This replaces the dead 256-word AB/BB certificate at \(L=8\), whose exact rank is
248, with a finite full-rank certificate.  It is not an all-length recursion.

**[LEMMA/COMPUTATION]** A complete exact search also rules out the requested
last-three-rung triangular transfer mechanism for the following natural stationary
class: four depth-\(\le4\) seed words and two distinct depth-two extension words.  All
144 families in that class which are full rank through \(L=7\) have exact
last-three-rung projection rank 248 at \(L=8\).  Thus none can have a
256-by-256 triangular boundary transfer with nonzero diagonal at that step.

**[UNRESOLVED]** Neither result proves \(D_{2L}\ge2^L\) for arbitrary \(L\).  In
particular, the boundary obstruction does not upper-bound the full Pauli rank of every
stationary candidate, and the patched \(L=8\) family is not claimed to recur.

## 2. Exact finite \(L=8\) replacement family

Let \(W_8(s)\), \(s\in\{0,1\}^8\), be the prior AB/BB family:

\[
 W_2(00)=A,\quad W_2(01)=BA,\quad W_2(10)=ABA,\quad W_2(11)=BABA,
\]

\[
 W_{L+1}(s0)=ABW_L(s),\qquad W_{L+1}(s1)=BBW_L(s).
\]

Words are written outermost generator first.  Sparse integer echelonization of its
last-three-rung projection retains 248 rows and rejects precisely the labels

\[
 R=\{10111011,10111101,10111110,11111010,
       11111100,11111101,11111110,11111111\}.
\]

For a six-bit string \(t=t_1\cdots t_6\), put \(m_0=AB\), \(m_1=BB\), and define

\[
 V(t)=m_{t_6}\cdots m_{t_1}\,BAB.
\]

Use the eight extension labels

\[
 T=\{000000,100000,010000,110000,001000,000100,100100,010100\}.
\]

Their literal words are

\[
\begin{array}{c|l}
 t&V(t)\\ \hline
000000&\mathtt{ABABABABABABBAB}\\
100000&\mathtt{ABABABABABBBBAB}\\
010000&\mathtt{ABABABABBBABBAB}\\
110000&\mathtt{ABABABABBBBBBAB}\\
001000&\mathtt{ABABABBBABABBAB}\\
000100&\mathtt{ABABBBABABABBAB}\\
100100&\mathtt{ABABBBABABBBBAB}\\
010100&\mathtt{ABABBBABBBABBAB}
\end{array}
\]

and every one has depth 15.  The certificate family is

\[
 \mathcal P_8=\{W_8(s):s\notin R\}\cup\{V(t):t\in T\}.
\]

It has 248 retained AB/BB rows and eight new rows.  Let \(\pi_{8,3}\) keep only
Pauli coordinates that are identity on the first five rungs, i.e. coordinates supported
on the final three full rungs.  Exact integer expansion and sparse modular echelon give

\[
 \operatorname{rank}_{\mathbb F_{p_1}}\pi_{8,3}(\mathcal P_8)
 =\operatorname{rank}_{\mathbb F_{p_2}}\pi_{8,3}(\mathcal P_8)=256.
\]

**[LEMMA — modular lower bound].** A rank-\(r\) integer matrix modulo a prime has
rational rank at least \(r\).

**Proof.** Its nonzero modular \(r\)-minor is an integer minor not divisible by that
prime, hence nonzero over \(\mathbb Q\). \(\square\)

Apply the lemma to the projected family.  Projection cannot increase rank, so the full
256 rows have rational rank at least 256; since there are only 256 rows, their rational
rank is exactly 256.  This proves the displayed finite lower bound for \(D_{16}\).
The artifact records the ordered word digest

```text
05e423043c744258d9b1f211b23ee98f4e7a39d806c1b179e6b81eb66921a8c8
```

and all 256 projected-row SHA-256 digests.

## 3. The eight old relations are not a removable reflection or parity artifact

Let \(\tau\) exchange ladder rows and \(\rho\) reverse rungs.  Both graph
automorphisms fix \(A_L\) and \(B_L\), hence fix every literal bracket word:

\[
 \tau W=W,\qquad \rho W=W.
\]

Thus row exchange and rung reversal act trivially on every member of the old AB/BB
family; they do not split its rows into \(+\) and \(-\) sectors that could separate the
relations.

There are, however, two useful character gradings.  Conjugation by
\(U_Z=\prod_vZ_v\) maps \(A\mapsto-A\), \(B\mapsto B\).  On the bipartite open
ladder, conjugation by \(U_E=\prod_{v\in E}X_v\), for either color class \(E\), maps
\(A\mapsto A\), \(B\mapsto-B\).  Therefore a literal word with \(n_A\) occurrences
of \(A\) and \(n_B\) occurrences of \(B\) has character

\[
 (n_A\bmod2,n_B\bmod2).
\]

**[COMPUTATION]** The eight primitive integer residual-zero relations from the prior
AB/BB certificate are homogeneous in these characters.  Their distribution, and the
resulting exact rank decomposition of the old 256 rows, is

\[
\begin{array}{c|c|c|c}
(n_A,n_B)\bmod2&\text{family rows}&\text{independent exact relations}&
 \operatorname{rank}_{\mathbb Q}\\ \hline
(0,0)&64&3&61\\
(0,1)&64&0&64\\
(1,0)&64&3&61\\
(1,1)&64&2&62
\end{array}
\]

For each sector, the relation-coefficient rows have the displayed exact rational rank,
and direct substitution into full integer Pauli expansions leaves an empty residual
dictionary.  The two-prime boundary ranks give the matching lower bounds.  Hence the
sector ranks are exact, not an inference from coincident modular ranks.

Consequently, projecting to a fixed character sector preserves each relation that lies
there; it cannot repair the defect.  A fixed sector also contains only 64 original rows.
This answers the symmetry question negatively: the eight relations are not caused by
mixing rung-reversal, row-swap, or character-parity sectors.

## 4. Exact stationary two-letter boundary-transfer obstruction

The searched class is deliberately finite and depth-compatible.

* Seeds are every unordered four-word subset of the 30 literal \(A/B\) words of
  depths 1 through 4 whose \(L=2\) vectors have exact rational rank four.
* Moves are every unordered distinct pair from
  \(\{AA,AB,BA,BB\}\).  A repeated move creates duplicate children and is already
  incapable of doubling.
* For each parent word \(w\), both children are \(m_0w\) and \(m_1w\).  A seed of
  depth at most four then remains at depth at most \(2L\).

The exact \(L=2\) seed-rank histogram over all \({30\choose4}=27405\) subsets is

\[
\begin{array}{c|rrrrr}
\text{rank}&0&1&2&3&4\\ \hline
\text{number of bases}&1001&8269&12708&4991&436.
\end{array}
\]

For every branch with a modular full minor, that minor is already a rational certificate.
For every modular-deficient branch, the experiment instead runs exact `Fraction`
echelonization, so the following is an exact exhaustive filter rather than a
prime-dependent screen:

\[
\begin{array}{c|c|l|c}
L&\text{input candidates}&\text{exact rank histogram}&\text{full-rank survivors}\\ \hline
3&2616&3:8,\ 4:88,\ 5:520,\ 6:920,\ 7:672,\ 8:408&408\\
4&408&14:48,\ 15:32,\ 16:328&328\\
5&328&30:64,\ 31:32,\ 32:232&232\\
6&232&61:40,\ 63:48,\ 64:144&144\\
7&144&128:144&144
\end{array}
\]

At \(L=5\), every survivor uses exactly the unordered move pair \(\{AB,BB\}\).
The 144 surviving base sets have SHA-256 digest

```text
c6a687a640feb0ed34462551bf06fbe0c59ecb4342e879d2668afd3181b7ad2a
```

For each of these candidates, \(\pi_{L,3}\) has full rank \(2^L\) at each checked
length \(2\le L\le7\).  At the next step, exact rational echelonization of every
256-row boundary matrix gives the single histogram

\[
 \operatorname{rank}_{\mathbb Q}\pi_{8,3}=248\qquad\text{for all 144 candidates}.
\]

A representative has rank 248 modulo both recorded primes as an arithmetic cross-check;
the upper bound 248 here is from exact rational echelon, not modular arithmetic.

**[LEMMA/COMPUTATION — restricted boundary-transfer obstruction].** No member of this
stationary two-letter class can supply a full-rank triangular last-three-rung transfer
at \(L=8\).

**Proof.** A triangular 256-by-256 boundary transfer with nonzero diagonal has boundary
rank 256.  The exhaustive exact calculation above gives boundary rank 248 for every
candidate that remains eligible through \(L=7\). \(\square\)

This lemma is intentionally narrow: it rules out the proposed boundary-transfer proof
mechanism in the stated class, not all two-generator word families and not the full
operator rank of every candidate.

## 5. Reproduction and limitations

From the repository root, run

```sh
timeout 3600 .venv/bin/python experiments/e82_ladder_all_l.py
timeout 3600 .venv/bin/python tests/test_ladder_all_l.py
```

The producer measured 1273.865094 seconds and the independent standalone verifier
measured 1084.06 seconds.  Both use only exact integer Pauli arithmetic, exact
`Fraction` echelon where an upper bound is required, and the two recorded primes only
for lower-bound certificates.

The standalone verifier does not import the producer.  It regenerates the complete
natural-class enumeration, all 144 exact \(L=8\) boundary ranks, the eight full-Pauli
relation residuals, the symmetry checks, and the 248-plus-8 patched family.
