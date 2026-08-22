# Two-leg ladder growth: finite certificates and a failed pivot rule

## 1. Result and scope

Let \(\Lambda_L\) be the open \(2\times L\) ladder, with sites numbered row-major, and put

\[
 A_L=\sum_{v\in\Lambda_L}X_v,
 \qquad
 B_L=\sum_{(u,v)\in E(\Lambda_L)}Z_uZ_v.
\]

For the depth filtration

\[
 F_1=\operatorname{span}\{A_L,B_L\},\qquad
 F_{k+1}=F_k+[A_L,F_k]+[B_L,F_k],\qquad D_k=\dim F_k,
\]

the requested all-length assertion is

\[
 D_{2L}(\Lambda_L)\ge 2^L\quad(L\ge2). \tag{1}
\]

**[COMPUTATION] Outcome (case (b)).** The fixed, translation-compatible boundary-lexicographic
leading-term rule defined in Section 4 is refuted at the finite length \(L=4\). It supplies enough
distinct raw leading Pauli words for \(L=2,3\), but only 15 distinct leaders among *all* nonzero
bracket words of depth at most 8 on \(\Lambda_4\), whereas a leading-term injection would require
16. An explicit collision inside the tested recursive binary family is

\[
 s=0101,\quad w(s)=\mathtt{BBABBA},
 \qquad
 t=0110,\quad w(t)=\mathtt{ABBBBA},
\]

and both words have leader

\[
 \operatorname{lead}(w(s))=\operatorname{lead}(w(t))
   =\mathtt{II\ II\ ZI\ XY}. \tag{2}
\]

The two bracket vectors in (2) are not proportional: their exact reductions have pair rank 2
modulo \(2^{31}-1\). Thus (2) is specifically a collision of the proposed raw pivot rule, not a
linear dependence of those two words.

**[CONJECTURE] The all-\(L\) inequality (1) remains open.** This work neither proves nor refutes
(1). In particular, the independent modular calculations still establish (1) for every
\(2\le L\le6\). No finite counterexample to the inequality itself was found.

No fitted dimension law or finite-size growth base is used anywhere in this conclusion.

## 2. Exact independent bracket engine

**[COMPUTATION]** `experiments/e28_ladder_proof.py` is independent of
`experiments/e20_algebra_growth.py`: it does not import that module, its orbit machinery, or its
basis builder. It uses the full Pauli basis.

Write

\[
 Q_{(a\mid b)}=X^aZ^b,
\]

encoded by the integer `a | (b << n)`. Every commutator is divided by 2. For Pauli words \(Q_g,Q_v\),

\[
 \frac12[Q_g,Q_v]=
 \begin{cases}
  0,&\langle g,v\rangle=0,\\
  \varepsilon(g,v)Q_{g\oplus v},&\langle g,v\rangle=1,
 \end{cases}
\]

where \(\varepsilon(g,v)\in\{+1,-1\}\) is evaluated exactly from the binary symplectic phases.
Starting from the exact integer coefficient dictionaries for \(A_L\) and \(B_L\), the code expands
every actual generator-left-nested word up to length \(2L\). A word string is outermost-first:
`ABBA` denotes \([A,[B,[B,A]]]\).

**[LEMMA]** These words span \(F_{2L}\).

**Proof.** The filtration recurrence applies either \(\operatorname{ad}_{A_L}\) or
\(\operatorname{ad}_{B_L}\) to every vector in the previous filtration. Induction on depth shows
that its spanning set is exactly the set of generator-left-nested words of length at most that
depth. Conversely every such word is produced by the recurrence. \(\square\)

The code eliminates the exact-integer coefficient rows modulo each of

\[
 p_1=2147483647,
 \qquad p_2=2147483629.
\]

**[LEMMA]** A rank \(r\) found modulo either prime is a rigorous lower bound on the rank over
\(\mathbb Q\).

**Proof.** Modular rank \(r\) exhibits an \(r\times r\) minor whose determinant is nonzero modulo
\(p\). Its integer determinant is therefore nonzero, hence the same rows have rank at least \(r\)
over \(\mathbb Q\). \(\square\)

Two-prime agreement is recorded only as an independent arithmetic cross-check; it is not described
as a rational upper-bound certificate.

## 3. Re-derived finite certificate and actual supports

**[COMPUTATION]** Full-Pauli elimination reproduces, independently and at both primes,

| \(L\) | depth | \(D_{2L}\) mod \(p_1\) | \(D_{2L}\) mod \(p_2\) | \(2^L\) |
|---:|---:|---:|---:|---:|
| 2 | 4 | 6 | 6 | 4 |
| 3 | 6 | 16 | 16 | 8 |
| 4 | 8 | 38 | 38 | 16 |
| 5 | 10 | 85 | 85 | 32 |
| 6 | 12 | 197 | 197 | 64 |

Hence \(D_{2L}\ge2^L\) over \(\mathbb Q\) for exactly this verified range, \(2\le L\le6\).

The JSON artifact stores an actual maximal independent list at each length: 6, 16, 38, 85, and
197 word records respectively. Each record contains

1. the concrete left-nested bracket word and depth;
2. every Pauli support code occurring with nonzero exact coefficient;
3. an SHA-256 digest of the exact `(Pauli code, integer coefficient)` support vector;
4. its raw leader under the order of Section 4;
5. its modular-echelon pivot under descending integer Pauli-code order; and
6. a histogram of touched-rung masks.

**[COMPUTATION] Structural observation.** The actual supports do not grow by simply appending one
Pauli letter at the right boundary. Already for \(L=4\), the first independent rows include

| depth | word | raw boundary leader | modular pivot | support size |
|---:|---|---|---|---:|
| 1 | `A` | `II II II XI` | `II II II IX` | 8 |
| 1 | `B` | `ZZ II II II` | `II II IZ IZ` | 10 |
| 2 | `BA` | `II II ZI YI` | `II II IZ IY` | 20 |
| 3 | `ABA` | `II II II YY` | `II II IY IY` | 20 |

The rung-mask histograms in the artifact show translated copies, disconnected touched-rung masks,
and overlapping boundary patterns. At larger depth, elimination pivots can differ from every raw
leader chosen before elimination. This is the concrete obstruction to reading an induction directly
from the finite ranks.

## 4. The tested translation-compatible leading rule

Fix a monomial order on Pauli strings uniformly for every \(L\):

1. compare rungs from right to left;
2. within a rung compare the top site, then the bottom site; and
3. at each site use the local order \(X>Y>I>Z\).

For a nonzero bracket vector \(W=\sum_q c_qQ_q\), define

\[
 \operatorname{lead}(W)=\max\{Q_q:c_q\ne0\}
\]

in this order. The rule is translation-compatible in the relevant finite-ladder sense: adding an
identity rung on the left leaves all existing comparisons unchanged, while the same local ordering
is used at every rung.

The binary recursive family tested exhaustively is

\[
\begin{array}{c|c}
 s & w_2(s)\\ \hline
00&A\\
01&[B,A]\\
10&[A,[B,A]]\\
11&[B,[A,[B,A]]]
\end{array}
\]

and for each added bit

\[
 w_{L+1}(s0)=[A,[B,w_L(s)]],
 \qquad
 w_{L+1}(s1)=[B,[B,w_L(s)]]. \tag{3}
\]

The maximum base depth is 4, and every recursion step adds 2, so every word in (3) has depth at
most \(2L\), exactly as required.

**[COMPUTATION]** Despite raw-leader collisions, the \(2^L\) rows in (3) have modular rank
\(2^L\) at both primes for every \(L=2,3,4,5,6\). This is useful finite evidence for (1), but it is
not an induction: the echelon pivots depend on eliminations against earlier full rows and no stable
all-\(L\) pivot formula was proved.

## 5. Exact failure certificate

There are two nested finite failures.

### 5.1 Capacity failure among all words

**[COMPUTATION]** Exhaustive enumeration of every nonzero left-nested word of depth at most \(2L\)
gives the following number of distinct raw leaders under the fixed order:

| \(L\) | nonzero words enumerated | distinct leaders | required \(2^L\) |
|---:|---:|---:|---:|
| 2 | 16 | 5 | 4 |
| 3 | 64 | 9 | 8 |
| 4 | 256 | 15 | 16 |
| 5 | 1024 | 27 | 32 |
| 6 | 4096 | 40 | 64 |

At \(L=4\), only 15 leaders exist among all eligible words. Therefore no injection whatsoever from
16 binary labels to eligible bracket words can have distinct leaders under this order. This is a
finite counterexample to the proposed pivot rule, not merely a poor choice of the family (3). The
largest verified length for the rule is \(L=3\).

### 5.2 Explicit collision in the recursive family

For \(L=4\), applying (3) gives

\[
 w_4(0101)=\mathtt{BBABBA},
 \qquad
 w_4(0110)=\mathtt{ABBBBA}.
\]

Exact integer expansion gives the shared maximum (2), whose encoded Pauli integer is 33928. The
artifact stores both full support dictionaries and their hashes. The common maximum remains the
same before any modular reduction, so the collision is characteristic-independent. Reduction of
the pair modulo \(p_1\) has rank 2; thus cancellation internal to either bracket has already been
accounted for, and the failure is precisely non-distinct leading monomials.

## 6. Why no leading-term induction is written

For ordinary polynomial products, one often uses
\(\operatorname{lead}(fg)=\operatorname{lead}(f)\operatorname{lead}(g)\). That identity cannot be
assumed for Pauli commutators of sums. If

\[
 U=\sum_u a_uQ_u,\qquad V=\sum_v b_vQ_v,
\]

then a target Pauli word \(Q_t\) in \([U,V]\) receives contributions from every anticommuting pair
with \(u\oplus v=t\). Distinct source pairs can land on the same \(t\) with opposite phase signs,
and the individually largest source terms may commute and contribute zero. A valid induction must
prove both uniqueness of the maximal target and non-cancellation of its coefficient for the
specific recursive family.

**[COMPUTATION]** Equation (2), and more strongly the 15-versus-16 capacity failure, shows that the
fixed boundary-lex rule does not have the required uniqueness property. Writing an induction for it
would therefore be false. The finite modular independence of (3) does not repair that missing
all-length argument.

## 7. Controls and reproducibility

**[COMPUTATION]** The independent engine also reproduces the required controls at both primes:

- open chains for \(2\le n\le6\) saturate at \(n^2\);
- rings for \(3\le n\le6\) saturate at \(3n-1\); and
- the total two-generator algebra for the \(2\times2\) ladder saturates at 11.

These are finite control computations, not new all-\(n\) equality theorems.

Run from the repository root:

```sh
.venv/bin/python experiments/e28_ladder_proof.py
.venv/bin/python tests/test_ladder_proof.py
```

Both programs print a final `PASS` on success. The machine-readable source is
`results/algebra_growth/ladder_proof.json`; its ranks, actual supports, collision, controls,
provenance, and checks are regenerated by the experiment.
