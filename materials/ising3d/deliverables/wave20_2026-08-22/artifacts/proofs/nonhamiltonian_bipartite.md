# Non-Hamiltonian connected bipartite local-term algebras

**[THEOREM — headline].** Let `Gamma` be a finite connected simple bipartite graph on `n>=2` vertices, with colour-class sizes `r,s`, and let

```text
g_Gamma = Lie_Q { i X_v, i Z_u Z_v : v in V(Gamma), uv in E(Gamma) }.
```

**[THEOREM — new branch].** If `Delta(Gamma)>=3`, then `g_Gamma` has the Pauli-label basis

```text
{(a,b) in F_2^n x F_2^n : |b| is even and Q_c(a,b)=1}
```

**[THEOREM — qualification].** The global label `P_X=(1,...,1;0)` is omitted from that displayed set when `n` is odd.

**[THEOREM — quadratic form].** Here `c` is either bipartition indicator and

```text
Q_c(a,b) = a.b + |a| + c.b  (mod 2).
```

**[THEOREM — exact dimension].** Consequently,

```text
dim g_Gamma = 2^(2n-2)-1                              if n is odd,
              2^(2n-2)-(-1)^r 2^(n-1)                 if n is even.
```

**[LEMMA — parity is well-defined].** When `n` is even, `r` and `s` have the same parity, so either colour-class size gives the same sign.

**[THEOREM — requested non-Hamiltonian classification].** Every connected bipartite graph without a Hamiltonian path has `Delta>=3`, so every such graph is in the preceding quadratic-root branch.

**[THEOREM — full connected-bipartite trichotomy].** Combining the new branch with the path/cycle results in `proofs/clifford_grade_classification.md` gives

| tag | connected bipartite graph | exact dimension |
|---|---|---:|
| **[THEOREM]** | path `P_n` | `n(2n-1)` |
| **[THEOREM]** | even cycle `C_n` | `2n(2n-1)` |
| **[THEOREM]** | `Delta>=3` | the parity formula above |

**[COMPUTATION — artifact].** The integrated exact artifact is `results/algebra_growth/nonhamiltonian_bipartite.json`; its producers are `experiments/e224_nonhamiltonian_census.py`, `experiments/e225_matching_grade.py`, and `experiments/e226_nonhamiltonian_theorem.py`.

**[COMPUTATION — independent verifier].** The clean-room verifier is `tests/test_nonhamiltonian_bipartite.py`; it imports no producer and rebuilds graph isomorphism reduction, all-pairs Pauli closure, matching, path cover, the quadratic root set, and Jordan--Wigner grades.

## 1. [LEMMA] Exact Pauli closure and the quadratic container

**[LEMMA — packed labels].** Encode a Pauli monomial modulo its nonzero scalar by `(a,b) in F_2^n x F_2^n`, where `a` records `X` support and `b` records `Z` support.

**[LEMMA — bracket rule].** For labels `p=(a,b)` and `q=(a',b')`, their Pauli commutator is nonzero exactly when

```text
B(p,q) = a.b' + b.a' = 1  (mod 2),
```

and its label is `p+q`.

**[LEMMA — label space].** Connectedness makes the bond vectors `e_u+e_v` span the even-weight subspace of `F_2^n`; together with all field vectors, the raw generators span

```text
H = {(a,b): |b| is even},     dim_F2 H = 2n-1.
```

**[LEMMA — colour-swap independence].** Replacing `c` by the complementary colouring adds `|b| mod 2` to `Q_c`, which vanishes on `H`; the quadratic form therefore does not depend on which colour class is labelled one.

**[LEMMA — polar identity].** Direct expansion gives

```text
Q_c(p+q) = Q_c(p)+Q_c(q)+B(p,q).
```

**[LEMMA — raw roots].** Every field has `Q_c(e_v,0)=1`, and every bond has

```text
Q_c(0,e_u+e_v)=c_u+c_v=1
```

because every graph edge crosses the bipartition.

**[LEMMA — containment].** If two `Q_c=1` labels have a nonzero bracket, their sum also has `Q_c=1`; hence every generated label lies in the `Q_c=1` level.

**[LEMMA — radical].** The symplectic radical of `H` is exactly

```text
rad(H) = span{P_X},     P_X=(1,...,1;0).
```

**[LEMMA — global parity cannot be generated].** If `p+q=P_X`, then `B(p,q)=B(p,p+P_X)=0`, so no nonzero commutator produces `P_X`; it is not a raw generator.

**[LEMMA — value on the radical].** The formula gives `Q_c(P_X)=n mod 2`, explaining why only odd `n` requires removing a `Q_c=1` radical label.

## 2. [LEMMA] A branching spanning tree supplies an `E6` root basis

**[LEMMA — branch-preserving spanning tree].** If `Delta(Gamma)>=3`, choose three edges incident to one degree-at-least-three vertex; those three edges form a forest and can be extended to a spanning tree `T` of `Gamma`.

**[LEMMA — raw tree basis].** The `n` fields and `n-1` bonds of `T` are linearly independent in `H`: the fields occupy independent `a` coordinates, while the tree incidence vectors are an independent basis of the even-weight `b` subspace.

**[LEMMA — commutation graph].** The commutation graph of that raw basis is the subdivision graph of `T`: a field `X_v` is adjacent precisely to the bond generators incident to `v`, while two fields or two bonds commute.

**[LEMMA — induced `E6`].** Let `v` be the retained branch vertex with neighbours `u_1,u_2,u_3`; the six basis elements

```text
X_v,
Z_v Z_(u_1), Z_v Z_(u_2), Z_v Z_(u_3),
X_(u_1), X_(u_2)
```

have commutation edges

```text
01, 02, 03, 14, 25,
```

which are the Dynkin tree `E6` as an induced subgraph.

**[COMPUTATION — finite certificate audit].** The artifact constructs such a spanning tree and checks GF(2) basis rank `2n-1`, all six pairings, and `Q_c=1` for every tree generator in each of the `20` branching isomorphism classes through `n=6`.

## 3. [EXTERNAL] The transvection-orbit input and its exact application

**[EXTERNAL — cited theorem].** Theorem 3.2 of Ahmet Seven, *Orbits of groups generated by transvections over F_2*, arXiv:math/0303098v1 (2003), states that when a basis is equivalent to a tree containing induced `E6`, the generated symplectic-transvection group has exactly two orbits outside the radical, cut out by the associated quadratic values `0` and `1`: <https://arxiv.org/abs/math/0303098>.

**[LEMMA — transvections are commutator steps].** For a raw basis label `g`, the transvection is

```text
tau_g(p)=p+B(p,g)g.
```

**[LEMMA — orbit equals generator-adjoint closure].** If `B(p,g)=1`, the move `p -> p+g` is exactly the label of the nonzero commutator with the raw generator `g`; if `B(p,g)=0`, the transvection fixes `p`.

**[LEMMA — the basis lies in one orbit].** The raw commutation graph is connected, and adjacent roots `x,y` satisfy `tau_y tau_x(y)=x`, so all raw basis labels lie in one transvection orbit.

**[LEMMA — exact tree closure].** The cited theorem, the induced-`E6` certificate, and `Q_c(g)=1` therefore identify the raw-tree orbit and exact generator-adjoint closure with every nonradical `Q_c=1` label.

**[LEMMA — extra edges do not change the answer].** Generator monotonicity gives `g_T subseteq g_Gamma`, while the quadratic containment gives `g_Gamma subseteq {Q_c=1}` and the radical argument removes `P_X`; hence both inclusions are equalities.

**[THEOREM — all-size status].** The all-size classification follows from the cited orbit theorem plus the explicit application above; it does not follow from the finite census.

## 4. [LEMMA] Exact counting without an Arf-invariant black box

**[LEMMA — count setup].** For fixed even-weight `b`, the dependence on `a` is

```text
Q_c(a,b) = a.(b+1) + c.b.
```

**[LEMMA — nonspecial fibres].** If `b` is not the all-ones vector, then `b+1` is a nonzero linear functional of `a`, so exactly `2^(n-1)` of the `2^n` choices of `a` have `Q_c=1`.

**[LEMMA — odd n].** If `n` is odd, the all-ones `b` has odd weight and is absent from `H`; all `2^(n-1)` allowed `b` fibres are balanced, giving `2^(2n-2)` roots before removing the single radical root `P_X`.

**[LEMMA — even n special fibre].** If `n` is even, the all-ones `b` is allowed and its constant quadratic value is

```text
c.(1,...,1) = r mod 2.
```

**[LEMMA — even n count].** The other `2^(n-1)-1` fibres contribute `2^(2n-2)-2^(n-1)` roots, while the special fibre contributes `0` when `r` is even and `2^n` when `r` is odd; this is exactly

```text
2^(2n-2)-(-1)^r 2^(n-1).
```

## 5. [THEOREM] Consequences for the requested graph families

**[THEOREM — non-Hamiltonian implication].** A connected graph of maximum degree at most two is a path or a cycle, and both have a Hamiltonian path; therefore every connected non-Hamiltonian graph has `Delta>=3`.

**[THEOREM — complete bipartite family].** A Hamiltonian path in `K_(r,s)` must alternate colours, so it exists exactly when `|r-s|<=1`; for every `|r-s|>=2`, `K_(r,s)` is an infinite non-Hamiltonian subclass governed by the quadratic-root formula.

**[THEOREM — trees].** Every non-path tree has a branch vertex and therefore has the same exact quadratic-root classification; its dimension depends on `n` and, for even `n`, the parity of a colour-class size, not on its detailed branching shape.

**[THEOREM — one-defect path covers].** A connected bipartite graph with minimum vertex-disjoint path-cover number two is classified by the same formula whenever it is non-Hamiltonian; path-cover number introduces no extra Lie-algebra branch.

## 6. [COMPUTATION] Deliberate invariant falsification

**[COMPUTATION — two six-vertex trees].** Define

```text
T_(2,4): edges 04,14,05,25,35,
T_(3,3): edges 04,14,25,35,45.
```

**[COMPUTATION — shared candidate invariants].** Both trees have `n=6`, no Hamiltonian path, maximum matching size `2`, matching deficiency `2`, and minimum path-cover number `2`.

**[COMPUTATION — decisive mismatch].** Their colour classes and exact dimensions are

```text
T_(2,4): parts (2,4), dim = 992,
T_(3,3): parts (3,3), dim = 1056.
```

**[THEOREM — narrowed invariant].** Matching deficiency together with path-cover number is therefore not a complete invariant; at even `n`, the dimension formula retains colour-class parity.

**[COMPUTATION — imbalance magnitude is unnecessary].** The star `K_(1,5)` has imbalance `4` and the balanced non-Hamiltonian tree `T_(3,3)` has imbalance `0`, but both have odd smaller-colour parity and exact dimension `1056`.

**[THEOREM — narrowed colour law].** Colour-imbalance magnitude is finer than the dimension requires; for fixed even `n`, only the parity of either colour-class size enters.

**[COMPUTATION — order dependence of grades].** For `T_(3,3)`, the identity vertex order gives the exact grade histogram

```text
{2:34, 4:256, 6:476, 8:256, 10:34},
```

while order `(0,2,1,3,5,4)` gives

```text
{2:66, 6:924, 10:66}.
```

**[THEOREM — grade limitation].** Reachable coordinate Clifford grades are therefore presentation data and cannot be an intrinsic graph invariant without specifying the Jordan--Wigner order.

**[COMPUTATION — grade-set counterexample].** In the identity order, both `T_(2,4)` and `T_(3,3)` occupy grades `{2,4,6,8,10}`, but their exact dimensions are `992` and `1056`; even the occupied-grade set at fixed `n` does not determine the algebra dimension.

## 7. [COMPUTATION] Exhaustive bounded census and exact limits

**[COMPUTATION — graph scope].** The producer enumerates every labelled simple graph through `n=6`, filters connected bipartite graphs, and reduces by the minimum edge bit-code over all vertex permutations, obtaining `27` isomorphism classes.

**[COMPUTATION — non-Hamiltonian scope].** Exactly `12` of those `27` representatives have no Hamiltonian path, and every one closes to its complete nonradical `Q_c=1` level.

**[COMPUTATION — exact closure].** Every stored dimension is the number of distinct raw Pauli labels reached by generator-adjoint breadth-first closure, and every reached set is independently checked under all nonzero pairwise brackets.

**[COMPUTATION — grade scope].** For each of the `12` non-Hamiltonian representatives, the producer audits all `n!` vertex orders and stores every distinct exact grade histogram with its order multiplicity.

**[COMPUTATION — verifier method].** The standalone verifier instead enumerates isomorphism signatures as adjacency-bit tuples, saturates new Pauli labels against all previously reached labels, derives Jordan--Wigner grades by a backward suffix recurrence, and compares all dimensions, digests, graph invariants, and order-profile multiplicities.

**[COMPUTATION — reproduction].** From the repository root, the producer and verifier commands are

```text
.venv/bin/python experiments/e226_nonhamiltonian_theorem.py
.venv/bin/python tests/test_nonhamiltonian_bipartite.py
```

**[UNRESOLVED — finite evidence].** The exhaustive census ends at `n=6`; no finite observation is promoted to an all-size theorem, and the all-size step explicitly uses the cited transvection-orbit theorem.

**[UNRESOLVED — excluded algebras].** This note classifies the algebra generated by every individual local field and bond; it does not classify the two-generator algebra generated by `sum_v X_v` and `sum_e Z_uZ_v`.

**[UNRESOLVED — Ising scope].** The classification gives no three-dimensional partition function, free energy, critical point, transfer spectrum, thermodynamic limit, or exact solution of the three-dimensional Ising model.
