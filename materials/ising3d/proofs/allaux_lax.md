# All-auxiliary Lax audit: exact method limitation and the valid Schur theorem

**[THEOREM — headline].** Exponential dimension and irreducibility of a physical local-term algebra do **not**, by themselves, imply the absence of local Lax, RLL, or RLLL intertwiners at arbitrary auxiliary dimension.  There are two exact all-dimension counterfamilies below: a strict nonidentity permutation RLLL solution for every common leg dimension, and the rational Yang RLL solution for every auxiliary/quantum dimension.  In the latter family the entries of the actual local Lax operator generate the full one-site matrix algebra and their chain envelope has dimension exponential in the number of sites.

**[THEOREM — strongest sound replacement].** Schur's lemma gives a zero intertwiner space only after one proves that the candidate intertwines **every element** of one named algebra between two named nonisomorphic simple representations.  A single RLLL product equation is only a restricted one-operator intertwiner equation.  Neither the dimension nor the irreducibility of a different physical local-term algebra supplies the missing all-elements equations.

**[UNRESOLVED — Ising scope].** No unrestricted higher-auxiliary Ising local factor is classified here.  The arbitrary-
\(8\times8\) theorem in `proofs/tetra_ungraded.md` and the fixed-order one-line-cabled \(16\times16\) theorem in `proofs/tetra16_locus.md` retain exactly their stated scopes; neither is broadened.

**[COMPUTATION].** The exact integrated artifact is `results/integrability/allaux_lax.json`.  Its producers are `experiments/e212_allaux_intertwiner.py`, `experiments/e213_allaux_controls.py`, and `experiments/e214_allaux_certificate.py`.  The clean-room verifier is `tests/test_allaux_lax.py` and imports none of the producers.

## 1. What the RLLL equation actually says

**[LEMMA — exact spaces].** Let \(k\) be a characteristic-zero field, put

\[
 A=A_1\otimes A_2\otimes A_3,
 \qquad
 Q=Q_4\otimes Q_5\otimes Q_6,
\]

**[LEMMA].** Let the three embedded local factors have the tetrahedral placements \(145,246,356\), and define

\[
 F=L_{145}L_{246}L_{356},
 \qquad
 G=L_{356}L_{246}L_{145}.
\]

**[LEMMA].** The vertex-type identical-\(L\) equation is

\[
 (R\otimes I_Q)F=G(R\otimes I_Q),
 \qquad R\in\operatorname{End}_k(A). \tag{1}
\]

**[LEMMA — coefficient map].** Equation (1) is the kernel of the linear map

\[
 \Phi_L:\operatorname{End}_k(A)\longrightarrow
 \operatorname{End}_k(A\otimes Q),
 \qquad
 R\longmapsto (R\otimes I_Q)F-G(R\otimes I_Q). \tag{2}
\]

**[LEMMA].** Thus the exact fixed-ansatz question is whether \(\ker\Phi_L\) contains a nonzero or invertible matrix.  This is precisely the question solved for the fixed binary-leg Ising tensor by `proofs/tetra_ungraded.md`.

**[LEMMA — one-operator interpretation].** If the tensor-factor restriction is temporarily dropped and \(T\) is allowed to range over all of \(\operatorname{End}_k(A\otimes Q)\), then \(TF=GT\) says exactly that \(T\) is a \(k[t]\)-module homomorphism when \(t\) acts by \(F\) on the source and by \(G\) on the target.  Restoring the RLLL ansatz gives

\[
 \ker\Phi_L
 =\operatorname{Hom}_{k[t]}((A\otimes Q,F),(A\otimes Q,G))
   \cap (\operatorname{End}_k(A)\otimes I_Q). \tag{3}
\]

**[THEOREM — missing bridge].** Equation (3) is not an intertwiner equation for the wave-18 physical local-term Lie algebra

\[
 \mathfrak g_\Gamma
 =\operatorname{Lie}_k\{iX_v,\ iZ_uZ_v\}.
\]

**[THEOREM].** To apply Schur to \(\mathfrak g_\Gamma\), one would first have to derive named representations \(\rho,\sigma\) and all equations

\[
 (R\otimes I)\rho(x)=\sigma(x)(R\otimes I)
 \quad\text{for every }x\in\mathfrak g_\Gamma. \tag{4}
\]

**[THEOREM].** No such derivation follows from the one equality (1).  Replacing (1) by (4) without proof is the exact logical gap in the conjectured inference.

**[LEMMA — auxiliary dimension is part of the ansatz].** The dimensions of \(A_1,A_2,A_3\) occur in the domains and codomains of the local factors.  Changing an auxiliary dimension therefore requires a newly specified family of local operators \(L^{(d)}\).  There is no canonical arbitrary-\(d\) coefficient matrix obtained by merely enlarging the 64 unknown entries of the fixed binary-leg \(8\times8\) problem.

## 2. The correct Schur and tensor-commutant statements

**[THEOREM — Schur zero criterion].** Let \(B\) be a unital associative \(k\)-algebra and let \(P,Q\) be finite-dimensional simple \(B\)-modules with representations \(\rho,\sigma\).  If

\[
 T\rho(b)=\sigma(b)T\qquad(b\in B), \tag{5}
\]

**[LEMMA].** Every nonzero \(T\) is therefore an isomorphism, and consequently

\[
 P\not\cong Q\quad\Longrightarrow\quad
 \operatorname{Hom}_B(P,Q)=0. \tag{6}
\]

**[LEMMA — proof].** The kernel of \(T\) is a \(B\)-submodule of \(P\), and its image is a \(B\)-submodule of \(Q\).  Simplicity makes a nonzero \(T\) injective and surjective, hence an isomorphism.  This proves (6) without any dimension estimate.

**[THEOREM — equivalent simple modules].** If, in addition, \(P\) is absolutely simple in the needed sense \(\operatorname{End}_B(P)=k\), and if \(S:P\to Q\) is one module isomorphism, then

\[
 \operatorname{Hom}_B(P,Q)=kS. \tag{7}
\]

**[THEOREM].** Every nonzero element of (7) is invertible.  Thus irreducibility can produce a one-dimensional space of invertible intertwiners; it does not inherently produce a no-go.

**[LEMMA — proof].** For any \(T\in\operatorname{Hom}_B(P,Q)\), the map \(S^{-1}T\) lies in \(\operatorname{End}_B(P)=k\), so \(T\) is a scalar multiple of \(S\).

**[THEOREM — tensor multiplicity].** Let \(U,V\) be finite-dimensional vector spaces on which \(B\) acts trivially.  Then

\[
 \operatorname{Hom}_B(U\otimes P,V\otimes Q)
 \cong
 \operatorname{Hom}_k(U,V)\otimes\operatorname{Hom}_B(P,Q). \tag{8}
\]

**[THEOREM].** In particular, if \(P=Q\) and \(\operatorname{End}_B(P)=k\),

\[
 (I_U\otimes\rho(B))'
 =\operatorname{End}_k(U)\otimes I_P,
 \qquad
 \dim (I_U\otimes\rho(B))'=(\dim U)^2. \tag{9}
\]

**[LEMMA — block proof].** Choose bases of \(U,V\) and write an operator \(T:U\otimes P\to V\otimes Q\) as a matrix of blocks \(T_{ij}:P\to Q\).  Equation (5) holds for \(T\) exactly when it holds for every block.  This gives (8), and (9) follows from (7).

**[LEMMA — double-commutant warning].** Over an algebraically closed field, Burnside's theorem can turn absolute irreducibility of an associative image into \(\rho(B)=\operatorname{End}(P)\), whose physical commutant is scalar.  Over \(\mathbb Q\), irreducibility alone does not even guarantee the scalar-commutant hypothesis; it must be stated separately.  After an auxiliary multiplicity space is included, the commutant is the full auxiliary algebra in (9), not zero.

**[THEOREM — valid application gate].** A Schur-based RLLL no-go is sound under the named assumptions: (i) the component equations imply (5) for every element of a specified algebra \(B\); (ii) the source and target are simple \(B\)-modules; and (iii) they are nonisomorphic.  The exponential wave-18 dimension formula supplies none of these three statements for \(F,G\) in (1).

## 3. Strict arbitrary-dimension RLLL counterfamily

**[THEOREM — permutation RLLL solution].** Let all six legs be copies of \(V=k^d\), for any \(d\ge2\).  On a three-leg space define

\[
 L=P_{12}\otimes I_3,
 \qquad
 R=P_{12}\otimes I_3, \tag{10}
\]

**[THEOREM].** Here \(P\) swaps the indicated tensor factors.  In the tetrahedral embedding this means

\[
 R_{123}=P_{12},\quad
 L_{145}=P_{14},\quad
 L_{246}=P_{24},\quad
 L_{356}=P_{35}. \tag{11}
\]

**[THEOREM].** These nonidentity permutation matrices satisfy

\[
 R_{123}L_{145}L_{246}L_{356}
 =L_{356}L_{246}L_{145}R_{123} \tag{12}
\]

**[THEOREM].** Equation (12) holds for every \(d\).

**[LEMMA — exact proof].** In the symmetric group on the six factor labels,

\[
 P_{12}P_{14}P_{24}=P_{24}P_{14}P_{12}, \tag{13}
\]

**[LEMMA].** The factor \(P_{35}\) commutes with all three factors in (13).  Multiplying (13) by that commuting factor proves (12).  Equality in the permutation group gives equality in every tensor-permutation representation, so this is an all-\(d\) theorem rather than finite extrapolation.

**[THEOREM — nonzero and invertible].** The auxiliary \(R=P_{12}\otimes I_3\) is nonscalar for \(d\ge2\), satisfies \(R^{-1}=R\), and has determinant

\[
 \det R=(-1)^{d^2(d-1)/2}\in\{+1,-1\}. \tag{14}
\]

**[THEOREM].** Thus (12) is a strict invertible RLLL counterexample to any bare assertion that auxiliary dimension alone forces a zero kernel.

**[LEMMA — local-entry size].** Relative to the first leg of \(L=P_{12}\otimes I_3\), its \((r,s)\) matrix entry is

\[
 L_{rs}=E_{sr}\otimes I_3. \tag{15}
\]

**[LEMMA].** The entries therefore contain the full active-leg algebra \(\operatorname{End}(V)\).  On \(n\) independent active physical legs, their associative envelope contains every tensor product of matrix units and is \(\operatorname{End}(V^{\otimes n})\), of exact dimension \(d^{2n}\).

**[THEOREM — interpretation].** The strict family proves that neither a large local-entry envelope nor arbitrary auxiliary dimension is an abstract obstruction to RLLL.  A theorem for the Ising tensor must use the specific way its physical local terms are encoded in its specific \(L^{(d)}\), not only the abstract dimension of a separate local-term algebra.

## 4. Known rational integrable counterfamily

**[EXTERNAL — identification].** The following permutation solution is the standard rational \(GL_d\)/XXX Yang \(R\)-matrix.  The external name is not used in the proof; the artifact verifies the polynomial identity exactly in \(\mathbb Z[u,v][S_3]\).

**[THEOREM — rational Yang identity].** On auxiliary copies \(V_a,V_b\) and a quantum copy \(V_q\), set

\[
 L_{aq}(u)=uI+P_{aq},\qquad
 L_{bq}(v)=vI+P_{bq},\qquad
 R_{ab}(u-v)=(u-v)I+P_{ab}. \tag{16}
\]

**[THEOREM].** For every \(d\ge2\),

\[
 R_{ab}(u-v)L_{aq}(u)L_{bq}(v)
 =L_{bq}(v)L_{aq}(u)R_{ab}(u-v). \tag{17}
\]

**[LEMMA — exact group-algebra proof].** Put \(A=P_{ab}\), \(B=P_{aq}\), and \(C=P_{bq}\).  The transpositions obey

\[
 AB=CA=BC,
 \qquad
 BA=AC=CB,
 \qquad
 ABC=CBA. \tag{18}
\]

**[LEMMA].** Expanding both sides of (17), grouping by the six permutations in \(S_3\), and using \(u-v-u+v=0\) makes every coefficient agree.  The producer stores the complete six-term normal form and a zero residual in \(\mathbb Z[u,v][S_3]\).

**[THEOREM — invertible exact sample].** At \(u=2,v=-1\), the auxiliary operator is \(R=3I+P\).  The swap has eigenvalue \(+1\) on the symmetric subspace of dimension \(d(d+1)/2\) and \(-1\) on the alternating subspace of dimension \(d(d-1)/2\).  Hence

\[
 \det(3I+P)
 =4^{d(d+1)/2}2^{d(d-1)/2}\ne0. \tag{19}
\]

**[THEOREM].** The solution is nonscalar and invertible for every \(d\ge2\).

**[LEMMA — actual Lax entries].** Relative to the auxiliary basis,

\[
 L(u)=\sum_{r,s}E_{rs}^{(a)}\otimes
       \bigl(u\delta_{rs}I+E_{sr}^{(q)}\bigr),
 \qquad
 L_{rs}(u)=u\delta_{rs}I+E_{sr}. \tag{20}
\]

**[LEMMA].** For \(r\ne s\), the entries give every off-diagonal matrix unit.  For any \(r\) and any \(s\ne r\), a product \(E_{rs}E_{sr}=E_{rr}\) gives the diagonal unit.  Therefore the associative algebra generated by the entries is exactly \(\operatorname{End}(V)\), and it acts absolutely irreducibly.

**[THEOREM — exponential chain envelope].** Embed these local entries independently at \(n\) quantum sites.  Products of the local matrix units give every

\[
 E_{r_1s_1}\otimes\cdots\otimes E_{r_ns_n}.
\]

**[THEOREM].** They form a basis of \(\operatorname{End}(V^{\otimes n})\), so the exact associative dimension is

\[
 d^{2n}. \tag{21}
\]

**[THEOREM].** This exponentially large irreducible envelope coexists with the exact local intertwiner (17) for every \(d,n\).

**[THEOREM — entry-derived full Lie control].** Let \(H=E_{11}-E_{22}\), take all onsite \(\mathfrak{sl}_d\) generators supplied by linear combinations of entries in (20), and add the nearest-neighbour products \(H^{(j)}H^{(j+1)}\), which lie in the range-two entry-product envelope.  For every \(d\ge2,n\ge1\), these local operators generate

\[
 \mathfrak{sl}_{d^n},
 \qquad
 \dim\mathfrak{sl}_{d^n}=d^{2n}-1. \tag{22}
\]

**[LEMMA — induction proof].** The base case is onsite \(\mathfrak{sl}_d\).  Suppose the first \(m\) sites generate \(\mathfrak{sl}_N\) with \(N=d^m\).  The edge coupling is \(X\otimes H\), where \(0\ne X=I^{\otimes(m-1)}\otimes H\in\mathfrak{sl}_N\).  Because \(\mathfrak{sl}_N\) and \(\mathfrak{sl}_d\) are simple in characteristic zero, iterated adjoints of their elements on \(X\) and \(H\) span the two full algebras.  Independent adjoints therefore generate \(\mathfrak{sl}_N\otimes\mathfrak{sl}_d\).  Together with the onsite summands,

\[
 \mathfrak{sl}_{Nd}
 = (\mathfrak{sl}_N\otimes I)
   \oplus(I\otimes\mathfrak{sl}_d)
   \oplus(\mathfrak{sl}_N\otimes\mathfrak{sl}_d), \tag{23}
\]

**[LEMMA].** This closes the induction.

**[LEMMA — physical caveat].** Equation (22) is a method-limitation theorem about what local-algebra size can imply when no compatibility theorem is supplied.  It does not claim that adding arbitrary controls to the XXX Hamiltonian preserves the XXX commuting transfer family.

## 5. Known-model decomposition controls

**[EXTERNAL].** For spin \(1/2\), the rational permutation density satisfies

\[
 P_{j,j+1}=\tfrac12(I+X_jX_{j+1}+Y_jY_{j+1}+Z_jZ_{j+1}). \tag{24}
\]

**[EXTERNAL].** This is the local XXX density associated with (16).

**[COMPUTATION — finite exact audit].** If the three Pauli summands \(XX,YY,ZZ\) on every open-chain bond are treated as separate local generators, exact characteristic-zero Pauli closure gives:

| tag | sites \(n\) | exact Lie dimension | exact reached labels |
|---|---:|---:|---|
| **[COMPUTATION]** | 3 | 15 | all noncentral labels with even \(X\)- and \(Z\)-parity |
| **[COMPUTATION]** | 4 | 60 | all noncentral labels with even \(X\)- and \(Z\)-parity |
| **[COMPUTATION]** | 5 | 255 | all noncentral labels with even \(X\)- and \(Z\)-parity |
| **[COMPUTATION]** | 6 | 1020 | all noncentral labels with even \(X\)- and \(Z\)-parity |

**[LEMMA — exact finite characterization].** For each displayed \(n\), the verifier independently constructs the even-even binary symplectic subspace.  Its radical contains only zero when \(n\) is odd and the four global \(I,X,Z,Y\) labels when \(n\) is even.  The raw commutator closure equals that explicit noncentral set coefficient-for-coefficient, not merely in dimension.

**[UNRESOLVED — no finite promotion].** The four rows are exact finite controls and are not promoted to an all-\(n\) Lie-closure theorem.  The all-size conclusions used in this note are instead the proved associative-envelope theorem (21) and the separately proved full-control theorem (22).

**[THEOREM — decomposition sensitivity].** A local-term algebra depends on which summands are declared independent generators, whereas RLL and RLLL are identities of the specified \(L\) and \(R\).  Therefore a dimension-only argument must first prove that its chosen generator decomposition is functorially encoded by the relevant intertwiner equations.

## 6. Consequences for the Ising claims

**[THEOREM — what wave 18 proves].** `proofs/clifford_grade_classification.md` proves an exact path/cycle/branching trichotomy for the physical Lie algebra generated by every \(iX_v\) and every \(iZ_uZ_v\) separately.  On the Hamiltonian branching branch its dimension is exponential.  That result is retained without change.

**[THEOREM — what wave 18 does not prove].** The trichotomy neither defines higher-dimensional local Ising tensors \(L^{(d)}\) nor turns the single product equation (1) into the all-elements equation (4).  Consequently it cannot be combined with Schur or a double-commutant slogan to conclude \(\ker\Phi_{L^{(d)}}=0\) for arbitrary \(d\).

**[THEOREM — fixed \(8\times8\) boundary].** `proofs/tetra_ungraded.md` proves that for the displayed identical binary-leg Ising tensor and every characteristic-zero \(q\notin\{0,1,-1\}\), the exact \(4096\times64\) coefficient matrix has rank 64 and its arbitrary \(8\times8\) auxiliary-\(R\) kernel is zero.  The present theorem does not weaken, reprove, or enlarge that result.

**[THEOREM — cabled \(16\times16\) boundary].** `proofs/tetra16_locus.md` proves the complete specialization locus only for its fixed-order one-line-cabled \(16\times16\) ansatz.  It explicitly does not classify a free \(\mathbb C^4\) local leg.  The present theorem does not replace that named cabling by an unrestricted auxiliary theorem.

**[UNRESOLVED — exact remaining problem].** For any proposed higher auxiliary dimension, one must first specify a genuine Ising local factor \(L^{(d)}\), its physical reconstruction assumptions, and the allowed equivalences/projections.  One must then either compute \(\ker\Phi_{L^{(d)}}\) directly or prove the all-elements, simple, nonisomorphic representation hypotheses required by (6).  No such unrestricted family is presently classified.

## 7. Exact controls and reproduction

**[COMPUTATION — universal certificates].** The artifact stores a zero six-permutation normal-form residual for (17) in \(\mathbb Z[u,v][S_3]\), and the exact \(S_6\) factor-permutation equality for (12).  The latter common factor permutation is `[3,1,4,0,2,5]` on zero-based labels.

**[COMPUTATION — finite RLL/RLLL controls].** At \(d=2,3,4\), the exact rational Yang three-factor matrices have sizes \(8,27,64\), zero residuals, and auxiliary determinants \(128,32768,67108864\).  The strict permutation RLLL six-factor matrices have sizes \(64,729,4096\), zero residuals, and auxiliary determinants in \(\{+1,-1\}\).

**[COMPUTATION — algebra controls].** Exact local-entry ranks are \(4,9,16\) for \(d=2,3,4\).  Tensor products give full chain ranks 64 at \((d,n)=(2,3)\) and 81 at \((3,2)\).  Small full-control qubit closures have dimensions \(3,15,63,255\) at \(n=1,2,3,4\).

**[COMPUTATION — Schur controls].** Exact rational row reduction gives commutant nullity one for full \(M_2\) and \(M_3\), intertwiner nullity one with an explicit determinant-one witness for equivalent simple \(M_2\)-modules, nullity zero for two nonisomorphic one-dimensional modules, and nullities \(4,9,4\) for tensor multiplicities \((\dim U,\dim P)=(2,2),(3,2),(2,3)\).  The one-operator identity equation on a two-dimensional space has nullity four while the full \(M_2\) intertwiner system has nullity one.

**[COMPUTATION — resource discipline].** The integrated producer records `time.process_time()` usage and peak RSS, uses exact integer/rational/symplectic arithmetic, and stays below the 2 GiB cap.  No modular rank is promoted to a rational equality and the benchmark value \(K_c\) is unused.

**[COMPUTATION — producer commands].** Regenerate the three exact producer layers with:

```bash
.venv/bin/python experiments/e212_allaux_intertwiner.py
.venv/bin/python experiments/e213_allaux_controls.py
.venv/bin/python experiments/e214_allaux_certificate.py
```

**[COMPUTATION — independent verification command].** Run the clean-room verifier with:

```bash
.venv/bin/python tests/test_allaux_lax.py
```

**[LEMMA — verifier independence].** The verifier imports no producer.  It separately rebuilds both permutation identities, all six finite RLL/RLLL controls, the local-entry and tensor-product bases, exact rational intertwiner ranks, tensor commutants, both packed-Pauli control families, source hashes, and the prior-theorem scope guards.

**[THEOREM — final scope].** The achieved theorem is an exact refutation of the algebra-size/irreducibility method at every auxiliary dimension, plus the precise Schur theorem that would make such a no-go valid under stronger named assumptions.  It is not an unrestricted higher-auxiliary Ising RLLL classification and does not solve the three-dimensional Ising model.
