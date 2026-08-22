# Onsager quotients: the closed-ideal classification and a conditional 2x3 no-go

Claim labels in this note are **[EXTERNAL]** (proved in a cited full source),
**[EXTERNAL-ABSTRACT]** (only an abstract was available and no full-text theorem is used),
**[LEMMA]** or **[THEOREM]** (proved here), **[COMPUTATION]** (an exact finite
calculation with an artifact), and **[UNRESOLVED]**.  The application to the 2x3
layer is deliberately labeled **[THEOREM, CONDITIONAL]**.

## 1. Scope and the retracted shortcut

The Onsager algebra below is the complex Lie algebra `OA`.  For a Lie algebra
initially defined over `Q` or `R`, “characteristic-zero Levi quotient” means the
semisimple quotient after extension of scalars to `C`.

**[UNRESOLVED]** This note does **not** determine the characteristic-zero Levi
structure of the 2x3 layer algebra.  That computation belongs to the separate
`CharZeroLevi` task.  The result proved here is the exact external implication
which that computation may, or may not, trigger.

The earlier Davies-based shortcut is retracted for a precise reason:
classifying finite-dimensional **irreducible representations** does not classify
arbitrary quotient Lie algebras.  An arbitrary quotient can have a radical and
its kernel need not be a “closed ideal” in the terminology of Date--Roan.  The
argument below never makes that inference.  Instead, it applies the closed-ideal
classification only to the canonical **semisimple quotient**, whose kernel is
closed for an elementary reason (Section 5).

## 2. Source audit

### Full sources used for mathematical statements

**[EXTERNAL]** E. Date and S.-S. Roan, *The structure of quotients of the
Onsager algebra by closed ideals*, J. Phys. A: Math. Gen. **33** (2000),
3275--3296, DOI
[`10.1088/0305-4470/33/16/316`](https://doi.org/10.1088/0305-4470/33/16/316),
arXiv [`math/9911018v3`](https://arxiv.org/abs/math/9911018).
The cached full PDF is
`sources/fulltext/date_roan2000_onsager.pdf`, SHA-256

```
c1334b58ed30c77b1592e628d8bf26c40d52ef63b9226ecad1ddec0e96771b06
```

The exact locations used are: the definition following Theorem 1 in Section 2
(loop fixed points); Lemma 1 and Theorem 2 in Section 3 (Chinese remainder and
closed ideals); Proposition 1 and Theorem 3 in Section 4 (non-fixed local
factors); and Proposition 4 plus the paragraph following it in Section 5
(fixed-point factors).  Section 6, Lemma 14 and Theorem 6, is consistent with
but is not needed for the quotient theorem proved here.

**[EXTERNAL]** S.-S. Roan, *Onsager's Algebra, Loop Algebra and Chiral Potts
Model*, MPIM preprint MPI/91-70 (1991), official
[MPIM scan](https://archive.mpim-bonn.mpg.de/3714/1/preprint_1991_70.pdf).
The cached full PDF is `sources/fulltext/roan1991_onsager.pdf`, SHA-256

```
82988421997ae8fb388c79c4f46381def00ee15ecaaed04922bd7ab94c352313
```

Proposition 1 and the following paragraph (scan PDF pages 11--12) give the loop
realization and fixed-point description.  Section 2, especially Lemma 7,
describes simultaneous evaluation at distinct non-fixed reciprocal orbits.
Date--Roan 2000 is the source used for the complete closed-ideal statement.

Both hashes and both manifest records are recomputed by
`tests/test_oa_quotient.py`; they also appear in
`results/algebra_structure/oa_quotient.json`.  The records follow the schema in
`sources/manifest.yaml` under keys `date_roan2000_onsager` and
`roan1991_onsager`.

### Abstract-only historical sources

**[EXTERNAL-ABSTRACT]** B. Davies, *Onsager's algebra and superintegrability*,
J. Phys. A **23** (1990), 2245--2261, DOI
[`10.1088/0305-4470/23/12/010`](https://doi.org/10.1088/0305-4470/23/12/010).
The IOP abstract was obtained; the full article was paywalled.  The abstract
says that irreducible representations are considered.  No theorem from its
unobtained full text is used here, and in particular that representation result
is not promoted to a statement about arbitrary quotients.

**[EXTERNAL-ABSTRACT]** B. Davies, *Onsager's algebra and the Dolan--Grady
condition in the non-self-dual case*, J. Math. Phys. **32** (1991), 2945--2950,
DOI [`10.1063/1.529036`](https://doi.org/10.1063/1.529036).  The Crossref
abstract was obtained, but not the full article.  Its abstract-level statement
is historical support only; the universal presentation used below is stated in
the full Date--Roan source.

These access limitations are recorded without hashes as `obtained:
abstract_only` under manifest keys `davies1990_onsager` and
`davies1991_onsager`.

## 3. Exact loop realization and terminology

Let `sl2` have basis `e,f,h` and brackets

\[
[e,f]=h,\qquad [h,e]=2e,\qquad [h,f]=-2f.
\]

Define

\[
\theta(e)=f,\qquad \theta(f)=e,\qquad \theta(h)=-h
\]

and the Chevalley-type loop involution

\[
\widehat\theta\bigl(p(t)\otimes x\bigr)
   =p(t^{-1})\otimes\theta(x)
\]

on `L(sl2)=C[t,t^{-1}] tensor sl2`.

**[EXTERNAL] (Date--Roan, Section 2; Roan, Proposition 1.)** The Onsager
algebra is the fixed-point algebra

\[
 OA=L(\mathfrak{sl}_2)^{\widehat\theta}.
\]

In the normalization of those sources,

\[
 A_m=2t^m e+2t^{-m}f,\qquad
 G_m=(t^m-t^{-m})h,
\]

and `OA` is universal on `A_0,A_1` subject to

\[
[A_1,[A_1,[A_1,A_0]]]=16[A_1,A_0],\qquad
[A_0,[A_0,[A_0,A_1]]]=16[A_0,A_1].                 \tag{DG}
\]

An element of `OA` is uniquely of the form

\[
p(t)e+p(t^{-1})f+q(t)h,\qquad q(t)=-q(t^{-1}).
\]

For a reciprocal polynomial `P`, Date--Roan write

\[
 I_P=\{X\in OA:P\text{ divides all three Laurent coefficients of }X\}.
\]

Their word **closed is algebraic, not topological**: for an ideal `I`, put

\[
 Z(I)=\{x\in OA:[x,OA]\subseteq I\}.
\]

Then `Z(I)/I` is the centre of `OA/I`, and `I` is called closed exactly when
`Z(I)=I`, equivalently when `OA/I` has zero centre.

## 4. What the closed-ideal classification actually says

For `a in C*`, set

\[
 U_a(t)=
 \begin{cases}
 (t-a)(t-a^{-1})=t^2-(a+a^{-1})t+1,&a\ne\pm1,\\
 t-a,&a=\pm1.
 \end{cases}
\]

Thus `U_a=U_(a^{-1})`.

**[EXTERNAL] (Date--Roan, Theorem 2.)** A proper ideal `I` of `OA` is closed
if and only if

\[
 I=I_P
\]

for a reciprocal polynomial `P` whose zero multiplicities at `t=1` and
`t=-1` are even.  Consequently a finite-codimensional closed ideal has, up to a
nonzero scalar, a unique factorization

\[
P(t)=(t-1)^{2r_+}(t+1)^{2r_-}
       \prod_{j=1}^{n}U_{a_j}(t)^{L_j},             \tag{1}
\]

where `r_+,r_- >= 0`, `L_j>=1`, `a_j != +/-1`, and the `a_j` represent
distinct inversion orbits.

**[EXTERNAL] (Date--Roan, Lemma 1.)** The Chinese remainder map for pairwise
coprime reciprocal factors is a Lie-algebra isomorphism.  Applied to (1), it
gives

\[
\begin{split}
OA/I_P\cong{}& OA/I_{(t-1)^{2r_+}}
 \oplus OA/I_{(t+1)^{2r_-}}\\
&\oplus\bigoplus_{j=1}^{n}OA/I_{U_{a_j}^{L_j}},     \tag{2}
\end{split}
\]

with absent (`r=0`) factors omitted.  Date--Roan prove surjectivity directly by
the polynomial Chinese remainder theorem while preserving the fixed-point
condition `q(t)=-q(t^{-1})`.

### 4.1 Points away from `t=+/-1`

Expand at `t=a_j` with `u=t-a_j`.

**[EXTERNAL] (Date--Roan, Proposition 1 and Theorem 3.)** If `a_j != +/-1`,
Taylor expansion induces

\[
 OA/I_{U_{a_j}^{L_j}}
   \cong \bigl(\mathbb C[u]/(u^{L_j})\bigr)
                 \otimes\mathfrak{sl}_2.            \tag{3}
\]

For `L_j=1`, this is the ordinary evaluation quotient `sl2`.

**[LEMMA 1 — radical of a truncated current factor.]** In the algebra in (3),

\[
N_j=u\bigl(\mathbb C[u]/(u^{L_j})\bigr)\otimes\mathfrak{sl}_2
\]

is exactly the solvable radical, and the quotient by it is `sl2`.

*Proof.* It is an ideal and

\[
[u^r x,u^s y]=u^{r+s}[x,y].
\]

Therefore its lower central series satisfies
`gamma_k(N_j) subseteq u^k C[u]/(u^{L_j}) tensor sl2` and terminates at
`k=L_j`; hence `N_j` is nilpotent.  Evaluation at `u=0` has kernel `N_j` and
image `sl2`.  Since the quotient is semisimple, any solvable ideal maps to zero
there, so no radical larger than `N_j` is possible. `[]`

The exact-integer bracket enumeration in
`experiments/e47_oa_quotient_sanity.py` checks this filtration for
`1 <= L <= 6`; it is only a finite sanity check, while the displayed degree
argument proves the lemma for every `L`.

### 4.2 The fixed points `t=+/-1`

Write `OA_(+1,L)=OA/I_(t-1)^L`; the involution `t -> -t` identifies the
`+1` and `-1` cases.

**[EXTERNAL] (Date--Roan, Lemma 4 and Proposition 4.)**

* `OA_(+1,1)` is the one-dimensional abelian algebra `C(e+f)`.
* For every `L>=1`, `OA_(+1,L)` is solvable of dimension
  `L+floor(L/2)`.  Date--Roan give a basis `X_k` (`0<=k<L`) together with
  `Y_j` for odd `j` (`0<=j<L`).
* For the even exponents that occur in a **closed** ideal, `L=2r>=2`, the
  quotient by the derived algebra is one-dimensional (spanned by `X_0`).
  Hence the factor has dimension `3r`, is solvable and non-abelian, and has
  zero Levi part.

Thus the fixed points do not hide an additional simple factor: their entire
algebra is radical.  The only semisimple contribution of (2) is one `sl2` from
each non-fixed reciprocal orbit, independent of the truncation exponent.

## 5. The critical corollary for arbitrary finite-dimensional quotients

**[THEOREM 2 — Levi corollary.]** Let `q` be any finite-dimensional complex
Lie algebra which is a quotient of `OA`.  Then

\[
 q/\operatorname{Rad}(q)\cong
   \mathfrak{sl}_2^{\oplus n}                       \tag{4}
\]

for some `n>=0`.  Equivalently, every Levi subalgebra of `q` is a direct sum of
copies of `sl2`; it has no simple factor of rank at least two and no simple
factor of any other rank-one type.

*Proof.* Let `pi:OA -> q` be surjective and let
`S=q/Rad(q)`.  If `S=0`, (4) holds with `n=0`.  Otherwise compose the maps:

\[
 OA\longrightarrow q\longrightarrow S.
\]

Its kernel `J` is **closed**, because `OA/J isomorphic S` and a semisimple Lie
algebra has zero centre.  This is the only point at which closedness is needed;
we make no assertion that `ker(pi)` itself is closed.

Apply Date--Roan Theorem 2 and the decomposition (2)--(3) to `OA/J`.  The
radical of an off-fixed factor is the positive-`u` ideal of Lemma 1, and a
fixed-point factor is entirely solvable.  Since `OA/J isomorphic S` is itself
semisimple, no fixed-point factor can occur and every off-fixed exponent must
be `L_j=1`.  Hence `S` is a direct sum of the resulting evaluation copies of
`sl2`.  Finally, the Levi--Malcev theorem identifies any Levi subalgebra of
`q` with `q/Rad(q)`. `[]`

This proof is the repair of the converse trap.  The classification of closed
ideals is not silently extended to arbitrary kernels; instead, the radical is
removed first, making the new kernel closed by definition.

## 6. Conditional theorem for the 2x3 layer algebra

Let `g_Q` be the characteristic-zero Lie algebra generated by the two 2x3
layer generators `A,B`, and put `g_C=g_Q tensor_Q C`.  (Equivalently, begin
with the complexification of the real dynamical Lie algebra.)

**[THEOREM, CONDITIONAL — Onsager-quotient no-go.]** If

\[
 g_C/\operatorname{Rad}(g_C)
\]

contains a simple ideal not isomorphic to `sl2`, then `g_C` is not isomorphic
to **any** quotient of the Onsager algebra.  In particular there is no pair of
elements of `g_C` which both generates all of `g_C` and satisfies the
normalized Dolan--Grady relations (DG).  Thus the generating pair `(A,B)` has
no Onsager-algebra realization after conjugation, change of basis,
automorphism, or replacement by a different generating pair: the obstruction
is to the abstract Lie algebra, not merely to the displayed generators.

*Proof.* Suppose `g_C isomorphic OA/I` for some ideal `I`, with no closedness
assumption.  Theorem 2 applied to `q=g_C` says that its semisimple quotient is
`sl2^n`, whose only simple ideals are copies of `sl2`.  This contradicts the
hypothesis.

For the final assertion, any generating pair `(X,Y)` satisfying (DG) induces,
by the universal presentation of `OA`, a Lie homomorphism `OA -> g_C` sending
`A_0` to `X` and `A_1` to `Y`.  Because `(X,Y)` generates `g_C`, that map is
surjective, again contradicting the first assertion. `[]`

The phrase “in any disguise” has exactly this abstract meaning.  It does not
claim that every conceivable integrability algebra must be Onsager, nor does it
address non-Lie or infinite-dimensional extensions.

## 7. Why the two modular `sp(14)` witnesses do not fire the theorem

**[COMPUTATION]** Reading the already-certified artifact
`results/algebra_structure/structure.json` (without recomputing the heavy
closure), `experiments/e47_oa_quotient_sanity.py` verifies all of the following
at each of the two recorded primes:

| prime | sector | module dimension | image dimension | exact finite-field type |
|---:|:---:|---:|---:|:---|
| 2147483647 | `000` | 14 | 105 | `C7 = sp(14)` |
| 2147483629 | `000` | 14 | 105 | `C7 = sp(14)` |

At both primes the stored invariant-form computation has a one-dimensional
space of alternating forms, form dimension 14, and kernel dimension zero.
These are exact statements over the two named finite fields.  The script also
checks that the source artifact retains `levi_factors: null`,
`levi_decomposition: UNDETERMINED`, and the explicit statement that no
characteristic-zero lift is claimed.

**[UNRESOLVED]** The modular witnesses do not prove the antecedent of the
conditional theorem.  A nonzero modular minor gives a lower bound on a
characteristic-zero rank; agreement at two primes is not a deterministic
rational upper certificate, and a simple modular image is not by itself a
certified characteristic-zero simple ideal.  In particular, it would be
incorrect to replace “if the characteristic-zero Levi quotient contains ...”
by “because a `C7` image appears modulo two primes.”

A sufficient next artifact would be an exact characteristic-zero radical/Levi
certificate, or an exact characteristic-zero surjection onto a certified
simple non-`sl2` algebra.  That is the sole missing antecedent and is being
computed independently; it is not assumed here.

## 8. Reproduction and artifact scope

Run only the dedicated standalone test:

```bash
.venv/bin/python tests/test_oa_quotient.py
```

It reruns `experiments/e47_oa_quotient_sanity.py`, recomputes both primary PDF
hashes, validates the result envelope and all computed checks, independently
reads both modular `C7` witnesses, enforces the characteristic-zero caveat, and
re-derives the small truncated-current lower-central dimensions using exact
integer brackets.  The generated certificate is
`results/algebra_structure/oa_quotient.json`.

**Exact scope.** Theorem 2 is an all-finite-dimensional-quotients statement
over `C`, conditional only on the cited Date--Roan classification.  The 2x3
conclusion is conditional on a separate characteristic-zero hypothesis.  No
claim about the 3D thermodynamic free energy, critical coupling, or an exact
solution follows from this note.

## 9. Synthesis: the antecedent is certified and the theorem is unconditional

*Added by the lead agent after both wave-6 computations landed; this section is
the only part of this note that depends on `proofs/char0_structure.md`.*

The sole missing antecedent named in Section 7 has now been produced.
`proofs/char0_structure.md` (experiments/e45, `tests/test_char0_levi.py`, all
checks passing) certifies, with exact rational arithmetic throughout:

* `dim_Q g = 263` by two independent exact closures (full-orbit Pauli sums and
  a faithful matrix closure), with the orbit-equivariance lemma proved and the
  earlier unsound one-representative shortcut replaced and cross-checked;
* `rad(g) = Z(g)` of dimension 1 (Killing rank 262, Cartan's criterion plus
  direct centrality);
* the semisimple quotient is `C7 + C3 + C3 + A8 + A5`
  (`sp(14) + sp(6) + sp(6) + sl(9) + sl(6)`, dimensions `105+21+21+80+35=262`),
  each factor certified isomorphic over `Q` to the named **split** classical
  algebra; the `C7` factor is distinguished from `B7` by a faithful
  14-dimensional rational module preserving a nondegenerate alternating form
  of determinant `16384`, and `B7 = so(15)` has no nontrivial module of
  dimension `< 15`.

**Base change.** Because every factor is certified isomorphic over `Q` to a
split classical algebra, complexification is immediate:
`g_C/Rad(g_C) ~ sp(14,C) + sp(6,C) + sp(6,C) + sl(9,C) + sl(6,C)`.  (For the
`C7` leg alone one can argue directly: the kernel of a `Q`-defined
representation base-changes as `ker(rho tensor C) = (ker rho) tensor C = 0`,
so the complexified image is a 105-dimensional subalgebra of `sp(14,C)`,
hence equals it.)  In particular `g_C/Rad(g_C)` contains the simple ideal
`sp(14,C)`, which is not isomorphic to `sl2`.

**[THEOREM — Onsager-quotient no-go, unconditional.]** The complexified 2x3
open-layer transfer Lie algebra `g_C` is not isomorphic to any quotient of the
Onsager algebra.  No pair of elements of `g_C` both generates `g_C` and
satisfies the normalized Dolan--Grady relations: the pair `(A,B)` admits no
Onsager-algebra realization after any change of basis, automorphism, or
replacement of generators.  *Proof.* The antecedent of the Theorem of
Section 6 holds by the certificate above; apply that theorem. `[]`

**Scope.** This is a statement about the single finite `2x3` layer algebra
(and, a fortiori, a statement that the 3D layer construction does not carry a
hidden Onsager structure at this size).  It is not an all-sizes theorem,
although the DG-defect theorem makes the same conclusion expected at every
size with a vertex of degree `>= 3`.  It does not restrict non-Lie,
infinite-dimensional, or non-Onsager integrability mechanisms.
