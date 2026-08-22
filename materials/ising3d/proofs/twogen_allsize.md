# Two-generator branching layers: exact counterexamples and the grade-separation wall

**[THEOREM — achieved negative statement].** Let

\[
 \alpha=iA=i\sum_{v\in V(\Gamma)}X_v,
 \qquad
 \beta=iB=i\sum_{\{u,v\}\in E(\Gamma)}Z_uZ_v,
 \qquad
 \mathfrak h_\Gamma=\operatorname{Lie}_{\mathbb Q}\langle\alpha,\beta\rangle
\]

act on the full physical `2^n`-dimensional spin space.  Branching, bipartiteness, and a Hamiltonian path do **not** by themselves force
`dim h_Gamma > n(2n-1)`: the Hamiltonian bipartite branching graphs `K_{2,3}` and `K_{3,3}` have exact dimensions

\[
 \dim\mathfrak h_{K_{2,3}}=44<45=5(2\cdot5-1),
 \qquad
 \dim\mathfrak h_{K_{3,3}}=63<66=6(2\cdot6-1).
\]

**[UNRESOLVED — exact `n_0`].** No universal eventual threshold is proved or disproved.  If an `n_0` with the requested quantifier exists, the two counterexamples force `n_0>=7`.  A complete finite census proves the inequality for every such graph at `n=7`, but one size is not an induction and is not reported as `n_0=7`.

**[THEOREM — contrast with local terms].** This negative result concerns the two **sums** `A,B`.  It does not weaken the exact local-term path/cycle/branching trichotomy in `proofs/clifford_grade_classification.md` (`H480`--`H487`), where every `X_v` and every `Z_uZ_v` is an independent generator and every branching graph has the stated exponentially large local-term algebra.  The distinction between those two generating sets is load-bearing.

**[THEOREM][EXTERNAL — proved grid scope].** For the physical open-grid family, the dimension obstruction is known for every open `2xL` grid with `L>=3` by `proofs/ladder_w8.md` (`H460`--`H462`).  Beyond the ladder, it is certified for `3x2`, `3x3`, and `3x4`; no all-`3xL` or all-rectangle induction is proved here.

**[COMPUTATION].** The machine-readable artifact is `results/algebra_growth/twogen_allsize.json`; its producers are `experiments/e182_adA_spectrum.py`, `experiments/e183_twogen_closure.py`, and `experiments/e184_twogen_allsize.py`.  The standalone verifier `tests/test_twogen_allsize.py` imports none of those producers and rebuilds all decisive quantities from raw generators.

## 1. Exact `ad_(iA)` spectrum on every Clifford grade

**[LEMMA — grade preservation].** Fix a Hamiltonian-path Jordan--Wigner order and write `N=2n`.  The conventions and phases are those already proved in `proofs/clifford_grade_classification.md`:

\[
 iX_r=-\gamma_{2r}\gamma_{2r+1},
 \qquad
 Z_pZ_q=i^{q-p}\gamma_{2p+1}\cdots\gamma_{2q}.
\]

Put `D=ad_alpha`.  Since `alpha` has Clifford grade two, its commutator with one Majorana is a linear combination of one Majorana.  The derivation rule then makes `D` act on `Lambda^k(Q^N)` without changing `k`; hence every Clifford grade is `D`-invariant.

**[LEMMA — exact spectrum and multiplicities].** Over `C`, define

\[
 \eta_r^+=\gamma_{2r}+i\gamma_{2r+1},
 \qquad
 \eta_r^-=\gamma_{2r}-i\gamma_{2r+1}.
\]

Direct Clifford multiplication gives

\[
 D\eta_r^+=-2i\eta_r^+,
 \qquad
 D\eta_r^-=+2i\eta_r^-.
\]

Choose `a` of the `n` minus-weight vectors and `b` of the `n` plus-weight vectors in an exterior monomial, with `a+b=k`.  The resulting eigenvalue and multiplicity are

\[
 \lambda_{a,b}=2i(b-a),
 \qquad
 m_{a,b}=\binom na\binom nb.
\]

Thus the grade-`k` eigenvalues are

\[
 2i\ell,
 \qquad
 \ell=-r,-r+2,\ldots,r,
 \qquad
 r=\min(k,2n-k),
\]

and the multiplicities sum exactly to `C(2n,k)` by Vandermonde's binomial identity.

**[LEMMA — coincidences across grades].** For even `k`, the eigenvalues are integer multiples of `4i`.  Every noncentral grade `2<=k<=2n-2` contains the three eigenvalues `-4i,0,+4i`; higher grades also share every eigenvalue allowed by the smaller radius.  Therefore `D`-eigenvalue separation is not, by itself, grade separation.

## 2. What Vandermonde separation does and does not recover

**[LEMMA — Vandermonde separation].** Let a diagonalizable operator `T` act on a finite-dimensional characteristic-zero complex vector space, and write

\[
 v=\sum_{j=1}^s v_{\lambda_j},
 \qquad
 Tv_{\lambda_j}=\lambda_jv_{\lambda_j},
\]

for the distinct eigenvalues occurring in `v`.  The coefficient matrix expressing `v,Tv,...,T^(s-1)v` in the vectors `v_{lambda_j}` is the Vandermonde matrix `(lambda_j^r)`.  Its determinant is `prod_{j<k}(lambda_k-lambda_j)`, so every `v_{lambda_j}` is a polynomial in `T` applied to `v`.  Applied to `T=D` and `v=beta`, every distinct `D`-eigencomponent of `beta` lies in the complexification of `h_Gamma`.

**[LEMMA — exact coincidence rule].** If grade components `v_{k,lambda}` and `v_{k',lambda}` share the same eigenvalue, Vandermonde inversion returns only their sum

\[
 v_\lambda=\sum_k v_{k,\lambda}.
\]

No polynomial in `D` distinguishes the summands at that eigenvalue.  This is the exact limitation relevant below.

## 3. Every Ising edge occupies the same three eigenvalues

**[LEMMA — endpoint occupancy].** For a path edge or chord at positions `p<q`, the Clifford monomial

\[
 \gamma_{2p+1}\gamma_{2p+2}\cdots\gamma_{2q}
\]

contains one Majorana from site pair `p`, both Majoranas from every interior site pair, and one Majorana from site pair `q`.  Each interior product `gamma_(2r)gamma_(2r+1)` commutes with `alpha`.  Only the two endpoints rotate under `D`.  Consequently **every** `Z_pZ_q`, independent of its Clifford grade `2(q-p)`, has nonzero `D`-eigenvalue support exactly

\[
 \{-4i,0,+4i\}.
\]

**[LEMMA — exact rational components].** In Pauli notation, factoring out the common leading `i` from the anti-Hermitian generator, the derivation satisfies

\[
 D(Z_uZ_v)=2(Y_uZ_v+Z_uY_v),
 \qquad
 D^2(Z_uZ_v)=8(Y_uY_v-Z_uZ_v),
\]

and `D(D^2+16)(Z_uZ_v)=0`.  Therefore the three real rational components of the full bond sum are

\[
\begin{aligned}
 B_0&=B+\frac1{16}D^2B
     =\frac12\sum_{uv\in E}(Z_uZ_v+Y_uY_v),\\
 B_c&=-\frac1{16}D^2B
     =\frac12\sum_{uv\in E}(Z_uZ_v-Y_uY_v),\\
 B_s&=\frac14DB
     =\frac12\sum_{uv\in E}(Y_uZ_v+Z_uY_v).
\end{aligned}
\]

They obey `DB_0=0`, `DB_c=4B_s`, `DB_s=-4B_c`, and `B=B_0+B_c`.  Over `C`, `B_c +/- iB_s` are the `-/+4i` components.  The producer and verifier check these brackets term by term for `2x3`, `2x4`, `3x3`, `K_{3,3}`, `K_{3,4}`, and the `2x2` control.

**[THEOREM — exact refutation of the proposed method].** Write `B=B_path+B_high` by Clifford grade.  Every polynomial `p(D)` that kills the three independent eigencomponents of all path bonds has

\[
 p(-4i)=p(0)=p(4i)=0.
\]

The same polynomial therefore kills every chord.  Conversely, every Vandermonde projector above retains the path and chord contributions together.  Thus the proposed inference `B_high in h_Gamma` does not follow; neither does subtraction to obtain the grade-two path-bond sum.

**[UNRESOLVED — failed induction step].** Without `B_high` or the separate grade-two path part, the suggested brackets cannot isolate the coordinate Clifford monomials needed for the requested explicitly counted all-size Pauli-string family.  This is the exact step that fails.  No finite computation is promoted to an induction in its place.

## 4. Complete-bipartite symmetry and exact counterexamples

**[LEMMA — symmetry container].** Every graph automorphism `g` fixes both sums `A` and `B`, so it fixes every Lie word in them.  If `L_Gamma` denotes the local-term algebra generated by the individual fields and bonds, then

\[
 \mathfrak h_\Gamma\subseteq(\mathfrak L_\Gamma)^{\operatorname{Aut}(\Gamma)}.
\]

The fixed space has one basis vector for each automorphism orbit of Pauli strings in the local-term support.  This supplies an exact finite upper container and explains why complete bipartite graphs, whose automorphism groups are large, can have much smaller two-sum algebras than their local-term algebras.

**[COMPUTATION — exact rational closures].** The producer and clean-room verifier enumerate those orbit coordinates, apply both adjoints with exact integer structure constants, and close over `Fraction` for the following small graphs:

| tag | graph | `n` | `|Aut|` | fixed local-term upper | `dim_Q h_Gamma` | `n(2n-1)` |
|---|---|---:|---:|---:|---:|---:|
| **[COMPUTATION]** | `K_{2,3}` | 5 | 12 | 49 | **44** | 45 |
| **[COMPUTATION]** | `K_{3,3}` | 6 | 72 | 63 | **63** | 66 |
| **[COMPUTATION]** | `K_{3,4}` | 7 | 144 | 174 | **167** | 91 |

**[THEOREM — symmetry explanation for `K_{3,3}`].** The `K_{3,3}` local-term algebra has 63 Pauli-string orbits under `(S_3 x S_3) semidirect C_2`, so symmetry gives `dim h<=63`.  The exact rational closure reaches all 63 orbit coordinates, proving equality.  This is a proved mechanism, not merely a fitted correlation between symmetry and a small dimension.

**[THEOREM — the two load-bearing negatives].** Both `K_{2,3}` and `K_{3,3}` are simple, bipartite, carry the alternating Hamiltonian path `0-1-...-(n-1)`, and have maximum degree at least three.  Their exact dimensions lie below `dim so(2n)`, refuting the claim that the graph hypotheses alone force the desired dimension obstruction at those sizes.

**[COMPUTATION — complete small-size census].** Relabel any chosen Hamiltonian path as `0-1-...-(n-1)`.  Bipartiteness permits exactly the odd-distance chords.  Exhausting every chord subset and retaining `Delta>=3` gives:

| tag | `n` | branching chord graphs | certified above ceiling | exceptions |
|---|---:|---:|---:|---|
| **[COMPUTATION]** | 4 | 0 | 0 | none exist |
| **[COMPUTATION]** | 5 | 3 | 2 | `K_{2,3}`, exact `44<45` |
| **[COMPUTATION]** | 6 | 14 | 13 | `K_{3,3}`, exact `63<66` |
| **[COMPUTATION]** | 7 | 63 | 63 | none at this size |

**[LEMMA — modular direction].** Each capped positive row contains `n(2n-1)+1` actual integer Lie words independent modulo the good prime `2147483647`; the corresponding rational words are independent.  The exceptional modular closures saturate at 44 and 63 and are independently closed over `Q`, so both inequalities are exact.

**[UNRESOLVED].** The `n=7` census supplies no monotonicity under adding vertices or changing the graph, and complete-bipartite symmetry persists at all balanced sizes.  It therefore does not settle any `n>=8` universal quantifier.

## 5. What is proved for open grids

**[COMPUTATION][EXTERNAL — exact finite anchors].** Fresh raw-generator good-prime closures reproduce the four requested two-generator dimensions.  Existing characteristic-zero sector certificates supply the opposite inequalities for the three larger rows:

| tag | open grid | `n` | `dim_Q h` | `n(2n-1)` | conclusion |
|---|---:|---:|---:|---:|---|
| **[COMPUTATION]** | `2x2` | 4 | 11 | 28 | dimension method inconclusive |
| **[COMPUTATION][EXTERNAL]** | `2x3` | 6 | 263 | 66 | fixed-basis quadratic pair excluded |
| **[COMPUTATION][EXTERNAL]** | `2x4` | 8 | 2952 | 120 | fixed-basis quadratic pair excluded |
| **[COMPUTATION][EXTERNAL]** | `3x3` | 9 | 8034 | 153 | fixed-basis quadratic pair excluded |

**[EXTERNAL].** The exact characteristic-zero upper/lower certificates are `proofs/char0_structure.md`, `proofs/char0_complete_2x4.md`, and `proofs/char0_3x3.md`.  This front independently rebuilds their raw modular closure ranks; it does not duplicate their sector and factor classifications.

**[THEOREM][EXTERNAL — infinite ladder family].** `proofs/ladder_w8.md` proves

\[
 \dim_{\mathbb Q}\mathfrak h_{2\times L}>2L(4L-1)
 \qquad\text{for every }L\ge3.
\]

Its `L=2` row is the genuine control `11<=28`.  This note cites the ladder theorem and does not claim its word/transport construction as part of the present front.

**[COMPUTATION][EXTERNAL — nonladder finite range].** `proofs/dla_3xl.md` and the later exact `3x3` certificate give `3x2=263`, `3x3=8034`, and `3x4>=1242>276`.  These prove the desired inequality for `3xL` only at `L=2,3,4`.

**[CONJECTURE — open-grid target].** Every open `a x b` grid with `a,b>=2` and `ab>=6` has

\[
 \dim\mathfrak h_{a\times b}>ab(2ab-1).
\]

**[UNRESOLVED — exact grid gap].** The cited ladder theorem settles all `2xL`, but `3xL` for `L>=5` and every family with `min(a,b)>=4` remain open.  The finite `3xL` binary word ranks in `proofs/dla_3xl.md` do not have a proved all-size block recurrence and are not used as an induction.

## 6. Fixed-basis transfer consequence and its limits

**[LEMMA — quadratic dimension ceiling].** Bilinears in the physical `2n` Majoranas form `so(2n)`, of dimension `n(2n-1)`.  If one fixed change of basis made both `alpha=iA` and `beta=iB` Majorana-bilinear, every iterated commutator would remain in that same `so(2n)`.  Hence

\[
 \dim\mathfrak h_\Gamma>n(2n-1)
 \quad\Longrightarrow\quad
 (A,B)\text{ are not simultaneously quadratic in a fixed physical basis}.
\]

**[THEOREM — transfer-factor corollary].** For every grid row certified above its ceiling, the displayed factorization

\[
 V=\exp(K^*A)\exp(KB)
\]

is not a standard free-fermion transfer factorization in which both exponent generators are quadratic in one fixed basis of the physical `2^n`-dimensional space.

**[UNRESOLVED — below-ceiling cases].** The inequalities `44<45`, `63<66`, and `11<28` do **not** prove that `K_{2,3}`, `K_{3,3}`, or `2x2` is free-fermionic.  They say only that this dimension obstruction is inconclusive.

**[THEOREM — precise non-spectral scope].** The Lie-dimension argument constrains the pair `(A,B)`, not the spectrum of one product.  It does not exclude `V` at a single coupling, or its eigenvalue multiset, from being Gaussian by an unrelated or coupling-dependent coincidence.  The pair-product spectral certificates in `proofs/allsize_gaussian.md` and the parity-half certificate in `proofs/gaussian_2x6.md` are logically independent results.

**[THEOREM — global scope].** Nothing in this note computes a thermodynamic free energy or critical coupling, and nothing solves the three-dimensional Ising model.

## 7. Verification

**[COMPUTATION].** From the repository root:

```sh
.venv/bin/python experiments/e182_adA_spectrum.py
.venv/bin/python experiments/e183_twogen_closure.py
.venv/bin/python experiments/e184_twogen_allsize.py
.venv/bin/python tests/test_twogen_allsize.py
```

**[COMPUTATION].** The verifier uses a separately written orbit-coordinate closure.  It reconstructs `A` and `B` from raw `X_v` and `Z_uZ_v` strings; reproduces `11,263,2952,8034`; closes `K_{2,3}`, `K_{3,3}`, and `K_{3,4}` over exact fractions; repeats all 80 branching chord cases at `n=5,6,7`; re-derives every displayed endpoint bracket for `n=6,8,9` and two non-grid controls; recomputes the exact factor-dimension arithmetic; checks artifact hashes; and prints final `PASS` under the declared process-time and RSS walls.

## Final status ledger

| tag | requested obligation | achieved status |
|---|---|---|
| **[LEMMA]** | grade preservation and exact spectrum | proved, including multiplicities and cross-grade coincidences |
| **[LEMMA]** | Vandermonde separation | proved only by distinct eigenvalue, exactly as valid |
| **[THEOREM — refuted method]** | recover `B_high` | impossible by `ad_A` alone: all edges occupy the same three eigenvalues |
| **[THEOREM][COMPUTATION]** | universal branching lower bound | refuted at `n=5,6` by exact complete-bipartite closures |
| **[UNRESOLVED]** | exact universal `n_0` | no value proved; if it exists, `n_0>=7` |
| **[THEOREM][EXTERNAL]** | positive physical family | all open `2xL`, `L>=3`; nonladder grids only `3x2,3x3,3x4` |
| **[THEOREM]** | transfer consequence | fixed simultaneous quadratic basis excluded exactly when the certified dimension clears the ceiling |
| **[THEOREM]** | spectral caveat | no conclusion about a single product spectrum from Lie dimension alone |
