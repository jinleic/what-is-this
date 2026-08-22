# Two-leg ladder Lie growth: hard-core-boson structure, state-space certificates, and the remaining sector-saturation lemma

## 0. Scope, notation, and status

Let \(\Lambda_L\) be the open two-leg ladder with \(L\) rungs, \(n=2L\) sites,
and \(3L-2\) edges.  With site \(2i\) on the top and \(2i+1\) on the bottom
of rung \(i\), set

\[
 A_L=\sum_{v\in\Lambda_L}X_v,\qquad
 B_L=\sum_{(u,v)\in E(\Lambda_L)}Z_uZ_v,\qquad
 \mathfrak g_L=\langle A_L,B_L\rangle_{\rm Lie}.
\]

All abstract statements below are over a field of characteristic zero unless a
finite field is named explicitly.  The finite certificates are computed over
\(\mathbb F_p\), \(p=2147483647\).  A modular lower rank is a valid
characteristic-zero lower rank: reduction modulo \(p\) cannot increase the
rank of a finite family of integral Lie words.

**[UNRESOLVED]** This note does not prove \(\dim\mathfrak g_L\ge 2^L\) for
all \(L\).  It proves a stronger all-\(L\) implication conditional on one
precisely named sector-saturation lemma, and gives an independent finite
certificate through \(L=8\).

The machine-readable evidence is
`results/algebra_growth/ladder_alll_proof.json`.  It uses the repository
`provenance` / `data` / `checks` envelope, records all launched phase walls
and process peak RSS values, and distinguishes the unlaunched longer-word
preflight from measured computations.

---

## 1. [THEOREM] Hard-core-boson decomposition and grading

Work in the simultaneous \(X\)-eigenbasis \(\{|x\rangle:x\in\{0,1\}^n\}\),
where

\[
 X_v|x\rangle=(-1)^{x_v}|x\rangle,
 \qquad
 N=\sum_v\frac{1-X_v}{2},
 \qquad
 \psi=|0\rangle=|+\rangle^{\otimes 2L}.
\]

Thus an occupied site has \(x_v=1\).  On every ladder edge define the
hard-core operations

\[
\begin{aligned}
 D&=\sum_{(uv)\in E}b_u^\dagger b_v^\dagger,\\
 F&=\sum_{(uv)\in E}\bigl(b_u^\dagger b_v+b_v^\dagger b_u\bigr),\\
 D^\dagger&=\sum_{(uv)\in E}b_ub_v.
\end{aligned}
\]

Here \(D\) fills an empty edge, \(F\) hops across a mixed edge, and
\(D^\dagger\) empties a full edge.

**Theorem 1.** For every \(L\ge2\),

\[
 A_L=2L-2N,\qquad B_L=D+F+D^\dagger,
\]

and \(D,F,D^\dagger\) are the \(\operatorname{ad}_{A_L}\)-eigencomponents
of \(B_L\), with eigenvalues \(-4,0,+4\), respectively.  In particular,

\[
\boxed{
 D=\frac1{32}[A,[A,B]]-\frac18[A,B],\quad
 D^\dagger=\frac1{32}[A,[A,B]]+\frac18[A,B],\quad
 F=B-\frac1{16}[A,[A,B]].
}
\]

Consequently \(D,F,D^\dagger\in F_3(\mathfrak g_L)\), and
\(\mathfrak g_L\) is graded by particle-number change:

\[
 \mathfrak g_L=\bigoplus_{r\in2\mathbb Z}\mathfrak g_L^{[r]},
 \qquad
 \mathfrak g_L^{[r]}\mathcal H_k\subseteq\mathcal H_{k+r}.
\]

*Proof.* \(Z_uZ_v\) flips both endpoint occupations.  Partitioning the four
endpoint patterns into empty-empty, mixed, and full-full gives
\(B=D+F+D^\dagger\).  Since
\([N,D]=2D\), \([N,F]=0\), and \([N,D^\dagger]=-2D^\dagger\), the displayed
\(\operatorname{ad}_A\) eigenvalues follow from \(A=2L-2N\).  Lagrange
interpolation in the three eigenvalues gives the three formulas.  The bracket
of homogeneous elements adds particle-number changes. \(\square\)

**[COMPUTATION]** The standalone verifier reconstructs the three pieces both
as exact rational vectors in the integer Pauli basis
\(Q_{(a\mid b)}=X^aZ^b\) and from their configuration action.  It checks
\(B=D+F+D^\dagger\), all three displayed formulas, and all three
\(\operatorname{ad}_A\) eigenvalue equations at \(L=2,3,4\).

---

## 2. [THEOREM] Evaluation and stabiliser-module bounds

**Theorem 2 (evaluation).** For every state \(\phi\),

\[
 \dim\mathfrak g_L\ge\dim(\mathfrak g_L\phi).
\]

More generally,

\[
 \dim\mathfrak g_L\ge
 \operatorname{rank}\bigl[T\mapsto(T\phi_1,\ldots,T\phi_s)\bigr]
\]

for any finite family of states.  This is the rank-nullity bound for the
indicated linear map.

**Theorem 3 (stabiliser module).** Let

\[
 \mathfrak s_\psi=\{S\in\mathfrak g_L:S\psi=0\}.
\]

Then \(\mathfrak s_\psi\) is a Lie subalgebra, \(\mathfrak g_L\psi\) is
invariant under every element of \(\mathfrak s_\psi\), and hence under the
unital associative algebra it generates.  In particular, for
\(v\in\mathfrak g_L\psi\),

\[
 \dim\mathfrak g_L\ge\dim\bigl(\mathcal S_\psi v\bigr).
\]

*Proof.* If \(S\psi=S'\psi=0\), then \([S,S']\psi=0\).  For
\(S\in\mathfrak s_\psi\) and \(T\in\mathfrak g_L\),

\[
 S(T\psi)=[S,T]\psi+T(S\psi)=[S,T]\psi\in\mathfrak g_L\psi.
\]

The last bound follows from Theorem 2. \(\square\)

At the vacuum, \(F\psi=D^\dagger\psi=0\), \(B\psi=D\psi\), every
negative-grade element annihilates \(\psi\), and every grade-zero element
acts by a scalar on the one-dimensional zero-particle sector.  The finite
certificate also constructs selected grade-\(+2\) linear combinations that
annihilate \(\psi\).

**[COMPUTATION]** Those selected annihilation relations are certified in
\(\mathbb F_p\), not asserted as rational relations.  This is sufficient:
the entire module calculation takes place in \(\mathfrak g_{L,\mathbb F_p}\psi\),
so its rank lower-bounds \(\dim\mathfrak g_{L,\mathbb F_p}\), which in turn
lower-bounds the characteristic-zero Lie dimension.  The standalone verifier
reconstructs every recorded generator, checks its vacuum annihilation, checks
the grade-zero scalar condition, and includes the negative control
\(D\psi\ne0\).

---

## 3. [THEOREM] Exact symmetry/parity ceiling

Let \(\tau\) exchange the two legs and \(\rho\) reverse the rungs.  Both
fix \(A_L\), \(B_L\), and \(\psi\), so every vector in
\(\mathfrak g_L\psi\) is \(\langle\tau,\rho\rangle\)-invariant.  The
operators \(A_L\) and \(B_L\) preserve particle-number parity.  Therefore

\[
 \mathfrak g_L\psi\subseteq
 \mathcal K_L:=\{\text{even-particle, \(\langle\tau,\rho\rangle\)-invariant vectors}\}.
\]

**Theorem 4.**

\[
 \dim\mathcal K_L
 =\frac{2^{2L-1}+3\,2^L}{4}
 =2^{2L-3}+3\,2^{L-2}
 \ge2^{2L-3}\ge2^L\qquad(L\ge3).
\]

*Proof.* Burnside's lemma applies to the four-element group generated by
\(\tau\) and \(\rho\), acting on even configurations.  The identity fixes
\(2^{2L-1}\) configurations.  Each of \(\tau\), \(\rho\), and
\(\tau\rho\) fixes \(2^L\) even configurations: for odd \(L\), the free
middle rung contributes two even states after the parity condition.  Averaging
gives the formula. \(\square\)

**[COMPUTATION]** Direct disjoint orbit-sum construction verifies the formula
for \(L=2,\ldots,9\):

| \(L\) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| \(\dim\mathcal K_L\) | 5 | 14 | 44 | 152 | 560 | 2144 | 8384 | 33152 |

The artifact records the direct orbit dimensions in every even sector, and
the verifier independently reconstructs those orbit-sum bases.

**[COMPUTATION]** Complete Pauli closure over \(\mathbb F_p\) gives

\[
\begin{array}{c|cc}
L&\dim\mathfrak g_{L,\mathbb F_p}&\dim(\mathfrak g_{L,\mathbb F_p}\psi)\\
\hline
2&11&3\\
3&263&14.
\end{array}
\]

Thus \(\dim(\mathfrak g_3\psi)=14=\dim\mathcal K_3\).  The verifier
independently recomputes this saturation and the \(263\) closure; it also
cross-checks the established stored \(L=3\) dimension in
`results/algebra_growth/exact_dimensions.json`.

---

## 4. [UNRESOLVED] The named missing sector-saturation/cyclicity lemma

For a fixed \(L\), write \(\mathcal K_L^{(k)}\) for the invariant
\(k\)-particle sector.

**Missing Lemma — SectorSaturationCyclicity(\(L\)) [UNRESOLVED].** For every even
\(k\in\{0,2,\ldots,2L\}\),

\[
 \mathcal K_L^{(k)}\subseteq\mathfrak g_L\psi.
\]

Equivalently,

\[
 \mathfrak g_L\psi=\mathcal K_L.
\]

This is the precise missing lemma.  It is **not** replaced here by an
irreducibility assertion about \(\mathfrak g_L^{[0]}\): such an assertion
would require a separate proof that the stabiliser-generated associative
algebra has the needed sector action and that suitable nonzero seeds occur in
every sector.

### Exact finite Krylov data relevant to the lemma

A finite exact check for one \(L\) can be made in the direct orbit-sum basis
of every \(\mathcal K_L^{(k)}\): choose actual seeds
\(v_j=X_j\psi\in\mathfrak g_L\psi\) in that sector and grade-zero
stabilisers \(S_i\in\mathfrak s_\psi\).  If the matrix with columns

\[
 S_{i_1}\cdots S_{i_r}v_j
\]

has exact rank \(\dim\mathcal K_L^{(k)}\) in every sector, then
SectorSaturationCyclicity(\(L\)) follows.  This is the exact state-space
Krylov data needed to settle a fixed \(L\); an all-\(L\) theorem would still
need a uniform proof of those ranks.

The artifact also records a deliberately narrower diagnostic: the
\(F\)-Krylov orbit of the selected seed
\(\operatorname{ad}_D^m(F)\psi\).  The exact modular ranks below are compared
with \(\dim\mathcal K_L^{(2m)}\):

| \(L\) | sector | seed | \(F\)-Krylov rank / \(\dim\mathcal K_L^{(2m)}\) |
|---:|---:|:---|---:|
| 3 | 2 | \(DF\psi\) | \(4/6\) |
| 3 | 4 | \(DDF\psi\) | \(4/6\) |
| 3 | 6 | \(DDDF\psi\) | \(0/1\) (zero seed) |
| 4 | 2 | \(DF\psi\) | \(8/10\) |
| 4 | 4 | \(DDF\psi\) | \(16/22\) |
| 4 | 6,8 | selected seed | \(0/10,0/1\) |
| 5 | 2 | \(DF\psi\) | \(12/15\) |
| 5 | 4 | \(DDF\psi\) | \(48/60\) |
| 5 | 6,8,10 | selected seed | \(0/60,0/15,0/1\) |

**[COMPUTATION]** The standalone verifier independently reproduces every
entry in this table.  These shortfalls rule out only this *selected-seed,
single-\(F\)* construction.  They neither disprove
SectorSaturationCyclicity(\(L\)) nor rule out other grade-zero stabilisers or
other seeds.

---

## 5. [THEOREM] Conditional consequence of sector saturation

**Theorem 5 (conditional on SectorSaturationCyclicity).** If
SectorSaturationCyclicity(\(L\)) holds for every \(L\ge3\), then, for every
\(L\ge3\),

\[
 \dim\mathfrak g_L
 \ge\dim\mathcal K_L
 \ge2^{2L-3}
 \ge2^L.
\]

*Proof.* The named hypothesis gives
\(\mathfrak g_L\psi=\mathcal K_L\).  Apply Theorem 2 and then Theorem 4.
\(\square\)

This is quadratically stronger than the original \(2^L\) lower-bound target,
but the hypothesis remains unresolved beyond the verified \(L=3\) saturation.

---

## 6. [COMPUTATION] Finite state-space certificate

For each row the producer uses all left-nested words of length at most four in
\(D,F,D^\dagger\) as seeds, selected negative-, zero-, and positive-grade
stabiliser combinations, and a deterministic signed fold of state coordinates.
The fold can only decrease rank.  Every used stabiliser is checked to annihilate
\(\psi\) in \(\mathbb F_p\) before use.

| \(L\) | target | observed modular lower rank | status |
|---:|---:|---:|:---|
| 2 | 4 | 3 | not enough; complete closure gives \(\dim\mathfrak g_2=11\) |
| 3 | 8 | 8 | certified |
| 4 | 16 | 16 | certified |
| 5 | 32 | 32 | certified |
| 6 | 64 | 64 | certified |
| 7 | 128 | 128 | certified |
| 8 | 256 | 256 | certified |
| 9 | 512 | 420 | non-certificate |

For \(L=3,\ldots,8\), the search deliberately stops when the target is
reached.  Thus the entries equal to \(2^L\) are certified lower bounds, not
claims that \(\dim(\mathfrak g_L\psi)=2^L\).

**[COMPUTATION]** The \(L=9\) run uses the same small recipe, has an explicit
335-second per-row wall budget with a one-second guard, and records the
observed wall/RSS in the artifact.  It finds the valid lower rank
\(420<512\), but the budget expires before a target certificate.  This is
not an upper bound on \(\dim(\mathfrak g_9\psi)\), not a failure of
SectorSaturationCyclicity(9), and not evidence against it.  A separate,
previous literal-Pauli certificate in `ladder_nonstat.json` remains the
repository's \(L=9\) \(512\)-dimensional certificate.

The finite recipe closure, with early stop disabled and no fold loss at
\(L\le6\), has modular ranks

\[
3,\ 13,\ 31,\ 71,\ 142\qquad(L=2,3,4,5,6),
\]

inside \(\mathcal K_L\) of dimensions \(5,14,44,152,560\).  These are
facts about the explicitly truncated stabiliser recipe only; they are not
claims about the full \(\mathfrak g_L\psi\).

---

## 7. [COMPUTATION] Cheap L=9 improvement attempts

### Three-state multi-vector evaluation

The direct map

\[
 T\longmapsto\bigl(T\psi_{+},T\psi_{-},T\psi_{\rm rung}\bigr)
\]

was evaluated over \(\mathbb F_p\), where \(\psi_+\) is all-plus,
\(\psi_-\) is all-minus, and \(\psi_{\rm rung}\) alternates empty and
fully occupied rungs in the \(X\) basis.  On the cheap length-\(\le4\)
pool of 120 left-nested \(D,F,D^\dagger\) words, the rank is **26** under
three deterministic signed folds of width 1024 each.

This is a valid multi-vector lower bound on \(\dim\mathfrak g_9\), but it
does not close the \(512\) gap.  Independently of the measured rank, its
120-column pool cannot reach 512.  The first unrestricted all-word pool with
at least 512 candidates is length \(\le6\), with 1092 words.

**[UNRESOLVED]** That longer multi-vector preflight was not launched: the
recursive state actions were not cheap within the explicit wall/RSS limits.
No conclusion is drawn from its absence.

### Rung-lex triangular screen

For the same bounded L=9 module search, state coordinates were ordered by

\[
(\text{particle number},\ \text{rung-adapted lexicographic order}).
\]

Vectors with distinct minimum-support coordinates in this order form a
triangular independent family.  Among the 420 projected-independent candidate
vectors examined by the bounded recipe, only **7** distinct such leaders occur.

**[COMPUTATION]** This exact leader screen therefore does not close the L=9
gap.  It is only a lower-bound screen for the tested recipe; it is not a
no-go theorem for other state-space triangular constructions.

---

## 8. [UNRESOLVED] Earlier obstructions remain scoped as documented

The earlier ladder results remain in force on the classes they address:

* `proofs/ladder_growth_family.md` documents the single-Pauli-word and
  related leading-word injection obstruction;
* `proofs/ladder_representation.md` records the AB/BB family rank
  \(248<256\) at \(L=8\);
* `proofs/ladder_all_l.md` records the stationary two-letter
  boundary-transfer obstruction;
* `proofs/ladder_nonstat.md` records its erase-a-rung obstruction and its
  separate literal-Pauli L=9 certificate.

Those results do not become false here.  They constrain their respective
operator-word constructions.  The present certificates instead rank state
vectors generated through a stabiliser module, so none of those scoped
obstructions alone proves or disproves SectorSaturationCyclicity(\(L\)).

---

## 9. Reproduction and claim ledger

From the repository root, run exactly:

```sh
PYTHONPATH=src .venv/bin/python experiments/e129_ladder_alll_proof.py
PYTHONPATH=src .venv/bin/python tests/test_ladder_alll_proof.py
```

The test does not import the producer.  It independently rebuilds the integer
Pauli identities, direct invariant orbit-sum bases, small full Lie closures,
all recorded modular stabiliser recipes and ranks, and the selected-seed
Krylov table.

| statement | status |
|---|---|
| hard-core decomposition, depth-3 formulas, and particle grading | **[THEOREM]**; exact checks at \(L=2,3,4\) |
| evaluation and stabiliser-module bounds | **[THEOREM]** |
| symmetry/parity ceiling and its Burnside formula | **[THEOREM]**; direct values through \(L=9\) |
| \(\dim(\mathfrak g_3\psi)=14=\dim\mathcal K_3\), with \(\dim\mathfrak g_3=263\) cross-check | **[COMPUTATION]** |
| \(\dim\mathfrak g_L\ge2^L\) for \(2\le L\le8\) | **[COMPUTATION]** |
| 420/512 at \(L=9\) for the small state-space recipe | **[UNRESOLVED]** non-certificate, not a negative result |
| SectorSaturationCyclicity(\(L\)) for all \(L\) | **[UNRESOLVED]** |
| \(\dim\mathfrak g_L\ge2^{2L-3}\ge2^L\) for all \(L\ge3\) | **[THEOREM]** (conditional on SectorSaturationCyclicity) |
