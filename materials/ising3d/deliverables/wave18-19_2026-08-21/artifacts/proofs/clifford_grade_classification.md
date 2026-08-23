**[THEOREM — exact refutation].** The proposed universal classification is false: the bipartite Hamiltonian graph `C_6` has one non-path edge, but its local-term algebra has grade set `{2,10}` on `N=12` Majoranas and exact dimension `132=2*C(12,2)`, not the proposed `1056`.

# Clifford grades and the exact path/cycle/branching trichotomy

**[COMPUTATION].** The integrated artifact is `results/algebra_growth/clifford_grade.json`; its producers are `experiments/e178_clifford_setup.py` through `experiments/e181_clifford_grade.py`, and its independent verifier is `tests/test_clifford_grade.py`.

## 0. Corrected classification and scope

**[THEOREM — exact trichotomy].** Let `Gamma` be a finite simple bipartite graph on `n` vertices that has a Hamiltonian path, let `N=2n`, and let

```text
g_Gamma = Lie_R { i X_v, i Z_u Z_v : v in V(Gamma), uv in E(Gamma) }
```

act on the full `2^n`-dimensional spin Hilbert space.  Relative to the Hamiltonian-path Jordan--Wigner order, exactly one of the following holds.

| tag | graph branch | generated Clifford grades | exact dimension |
|---|---|---|---:|
| **[THEOREM]** | `Gamma` is the path | `{2}` | `C(N,2)=n(2n-1)` |
| **[THEOREM]** | `Gamma` is the even cycle | `{2,N-2}` | `2C(N,2)=N(N-1)=2n(2n-1)` |
| **[THEOREM]** | `Delta(Gamma)>=3` | every `k=2 mod 4`, with `k=N` omitted | `2^(2n-2)-(-1)^(n/2)2^(n-1)` for even `n`; `2^(2n-2)-1` for odd `n` |

**[THEOREM — relation to the requested statement].** The requested hypotheses include at least one non-path edge and therefore exclude the path branch, but they do not exclude the cycle branch.  Every even cycle `C_n` with even `n>=6` is a counterexample.  The `2x2` grid is `C_4`; it is in the cycle branch, but its two grades `{2,6}` accidentally are the whole `2 mod 4` class for `N=8`, so its dimension `56` did not expose the missing hypothesis.

**[LEMMA — branching equivalence].** For a graph carrying a Hamiltonian path, the following are equivalent: some non-path edge is not the pair of path endpoints; some chord has grade `k<=N-4`; and `Delta(Gamma)>=3`.

**[LEMMA — proof].** Every interior path vertex already has path degree two.  A chord other than the endpoint pair meets at least one interior vertex and raises its degree to at least three.  Conversely, if `Delta<=2`, no chord can meet an interior path vertex, so the graph is either the path or the path plus its endpoint edge, namely a cycle.  The endpoint edge has path distance `n-1` and grade `2(n-1)=N-2`; every other chord has distance at most `n-2` and hence grade at most `N-4`. `[]`

## 1. Obligation 1: exact Jordan--Wigner setup

**[LEMMA — status: proved].** Number path positions from `r=0` to `n-1` and use Hermitian Pauli matrices with `XY=iZ`, `YZ=iX`, and `ZX=iY`.  Define

```text
gamma_(2r)   = X_0 X_1 ... X_(r-1) Z_r,
gamma_(2r+1) = X_0 X_1 ... X_(r-1) Y_r.
```

Then `gamma_a gamma_b + gamma_b gamma_a=2 delta_ab`, with no suppressed phase.

**[LEMMA — field and path-bond phases].** Direct local multiplication gives

```text
X_r       = i gamma_(2r) gamma_(2r+1),
Z_r Z_(r+1) = i gamma_(2r+1) gamma_(2r+2).
```

Thus the anti-Hermitian generators are

```text
iX_r          = -gamma_(2r) gamma_(2r+1),
iZ_rZ_(r+1)   = -gamma_(2r+1) gamma_(2r+2).
```

**[LEMMA — arbitrary edge phase].** For path positions `p<q`, put `d=q-p`.  Grouping consecutive Majoranas into path-bond pairs gives

```text
gamma_(2p+1) ... gamma_(2q)
 = product_(r=p)^(q-1) [gamma_(2r+1) gamma_(2r+2)]
 = (-i)^d product_(r=p)^(q-1) Z_r Z_(r+1)
 = (-i)^d Z_p Z_q.
```

Consequently

```text
Z_p Z_q = i^d gamma_(2p+1) ... gamma_(2q),
```

which is one Clifford monomial of exact grade `2d=2(q-p)`.  In particular, for odd `d`,

```text
i Z_p Z_q = i^(d+1) gamma_(2p+1) ... gamma_(2q)
           = (-1)^((d+1)/2) gamma_(2p+1) ... gamma_(2q).
```

**[LEMMA — path generates grade two].** The fields and path bonds give every consecutive bilinear `gamma_a gamma_(a+1)`, for `0<=a<N-1`.  For distinct `a,b,c`,

```text
[gamma_a gamma_b, gamma_b gamma_c] = 2 gamma_a gamma_c.
```

Induction on index distance therefore gives every `gamma_a gamma_b`, so the path terms generate the full grade-two component `Lambda^2(Q^N)`, isomorphic to `so(N)` and of dimension `C(N,2)=n(2n-1)`. `[]`

**[COMPUTATION — phase audit].** `e178` checks all Clifford relations, every field identity, and every edge identity for all `n<=9` using exact powers of `i`; all 197 setup/formula checks in the integrated artifact pass.

## 2. Obligation 2: the bipartite grade lemma and a non-bipartite test

**[LEMMA — status: proved].** Along a Hamiltonian path the vertex colours alternate.  In a bipartite graph every edge joins opposite colours, so its path positions differ by an odd number.  Every non-path edge therefore has grade

```text
2(q-p) = 2 mod 4.
```

A non-path edge has distance at least three, because distance one is already a path edge, so every chord grade is at least six.

**[LEMMA — non-bipartite converse].** The parity of path position is itself a valid two-colouring of the path.  If the full graph is non-bipartite, some additional edge violates this colouring; it joins positions of the same parity, has even path distance, and therefore has grade `0 mod 4`.

**[THEOREM — exact interpretation].** One `0 mod 4` generator already proves that the algebra is not contained in the `2 mod 4` span.  It does not by itself force the full even Clifford commutator algebra: for example, the bare odd cycle `C_5` has only grades `{2,8}` and dimension `90`.  The stronger full-even outcome below uses both the grade-six chord and the grade-eight odd-cycle edge.

**[COMPUTATION — required falsification].** Take the Hamiltonian path `0-1-2-3-4` and the graph `C_5` plus chord `(0,3)`.  Its two non-path edges have distances three and four, hence grades six and eight.  Exact raw-generator closure gives

```text
grade 2:  C(10,2) = 45
grade 4:  C(10,4) = 210
grade 6:  C(10,6) = 210
grade 8:  C(10,8) = 45
total:                 510.
```

The `2 mod 4` prediction is `45+210=255`; all even monomials including the identity number `512`; omitting only the identity gives `511`; and the full noncentral even commutator algebra, which omits both identity and volume, has dimension `512-2=510`.  The measured closure is exactly the last alternative.  The standalone verifier recomputes the 510 strings, their full histogram, pairwise bracket saturation, and the artifact digest.

## 3. Obligation 3: grade generation, exact failure, and corrected lemma

**[LEMMA — monomial bracket rule].** For ordered supports `I,J` with sizes `a,b` and overlap `j=|I intersection J|`, Clifford multiplication gives a scalar sign times `gamma_(I symmetric-difference J)`.  The two multiplication orders differ by the sign

```text
(-1)^(ab-j).
```

For even `a,b`, the commutator is nonzero exactly when `j` is odd, and its output grade is

```text
a+b-2j.
```

Supports with the prescribed sizes and overlap exist exactly when

```text
0 <= j <= min(a,b),       a+b-j <= N.
```

**[THEOREM — status: the requested lemma is refuted].** Grade two plus one grade `k` with merely `6<=k<=N` does not always generate the full class.  If `N=0 mod 4` and `k=N-2`, then the only generated grades are `{2,N-2}`.

**[THEOREM — endpoint proof].** Two `(N-2)`-subsets are complements of two-element sets and therefore overlap in at least `N-4` indices.  The possible overlaps are `N-4`, `N-3`, and `N-2`; because `N` is even, only `N-3` is odd.  Its bracket grade is

```text
2(N-2)-2(N-3)=2.
```

A grade-two monomial bracketed with a grade-`N-2` monomial either commutes or stays at grade `N-2`.  Hence `{2,N-2}` is bracket-closed. `[]`

**[COMPUTATION — load-bearing counterexample].** For `C_6`, `N=12` and the endpoint chord has `k=10=N-2`.  Exact raw Pauli closure and independent pairwise saturation give the histogram `{2:66,10:66}` and dimension `132`; the proposed full class would add `C(12,6)=924` and have dimension `1056`.

**[THEOREM — corrected generation lemma].** Let `N` be even and let `k=2 mod 4` satisfy

```text
6 <= k <= N-4.
```

Grade two plus one coordinate grade-`k` monomial generates every noncentral grade `2 mod 4`.

**[LEMMA — proof, descent to grade six].** By the grade-two orbit lemma in the next section, one grade-`k` monomial supplies the whole grade-`k` coordinate space.  Choose two `k`-subsets overlapping in exactly `k-3` indices.  Their union has size `k+3<=N-1`, so such supports exist; `k-3` is odd; and their bracket has grade

```text
2k-2(k-3)=6.
```

**[LEMMA — proof, raising from grade six].** If grade `r=2 mod 4` is present, choose a six-subset and an `r`-subset overlapping in one index.  Their bracket has grade `r+4`, and the required union has size `r+5`.  For `N=0 mod 4`, this reaches the last class grade `N-2` from `r=N-6`, with union size `N-1`.  For `N=2 mod 4`, it reaches the last noncentral class grade `N-4` from `r=N-8`, with union size `N-3`; the central top grade `N` is neither needed nor reachable.  This proves the lemma. `[]`

**[LEMMA — small cases].** At `N=8`, the exceptional set `{2,N-2}={2,6}` already is the whole class, explaining `C_4`.  At `N=10`, the corrected seed is `k=6=N-4` and the noncentral target is only `{2,6}`.  At `N=12`, seed `k=6` generates `{2,6,10}`, whereas the endpoint seed `k=10` generates only `{2,10}`.

**[COMPUTATION].** `e179` exhausts every corrected seed through `N=80`, separately exhausts the endpoint failures through `N=64`, and directly checks the support formula through `N=10`; all 361 checks pass.

## 4. Obligation 4: orbit coverage, including middle grade

**[LEMMA — elementary orbit proof over `Q`].** Let `I` be a `k`-subset with `0<k<N`.  If `a` is in `I` and `b` is not in `I`, then

```text
[gamma_a gamma_b, gamma_I] = +/-2 gamma_((I-{a}) union {b}).
```

Thus grade-two brackets perform the edges of the Johnson graph `J(N,k)`.  That graph is connected, so the `so(N)`-cyclic span of one coordinate monomial contains every coordinate grade-`k` monomial.  This proof works over `Q`, includes middle grade without an exception, and also supplies explicit bracket words.

**[EXTERNAL — standard representation description].** Over an algebraic closure of characteristic zero, `Lambda^k(V)` for the standard `so(2n)` module is irreducible when `k!=n`.  At `k=n` it splits into two Hodge-star eigensummands, over a quadratic extension if the star eigenvalues are not already in the ground field.

**[LEMMA — requested middle-grade argument].** The middle grade belongs to the `2 mod 4` class exactly when `n=2 mod 4`.  For a coordinate decomposable monomial `gamma_I`, Hodge star sends

```text
*gamma_I = c_I gamma_(I^c),       c_I != 0.
```

The supports `I` and `I^c` are disjoint, so the two vectors are linearly independent and both Hodge projections of `gamma_I` are nonzero.  The two Hodge summands are non-isomorphic simple `so(2n)` modules; hence a cyclic submodule with nonzero projection to each is their full direct sum.  After scalar extension the cyclic space therefore has dimension `C(2n,n)`.  Scalar extension preserves dimension, so the original `Q`-cyclic space has the same dimension and is all of `Lambda^n(Q^(2n))`. `[]`

**[COMPUTATION].** The artifact records six exact grade-two swaps taking `{0,1,2,3,4,5}` to its disjoint Hodge complement `{6,7,8,9,10,11}` at `N=12`, plus non-middle orbit samples at `N=18` and `N=24`.

## 5. Obligation 5: centre and the odd-`n` correction

**[THEOREM — status: the requested scalar claim is refuted].** The top monomial

```text
omega = gamma_0 gamma_1 ... gamma_(N-1)
```

is central in the even Clifford algebra, but it is not scalar on the full `2^n`-dimensional irreducible Clifford representation.  Under the conventions of Section 1,

```text
P_X := product_(r=0)^(n-1) X_r = i^n omega.
```

Thus `omega=i^(-n)P_X`; global `X` parity is a nonidentity, traceless Pauli string on the full physical space.  It becomes a scalar only after restriction to either parity/half-spin sector, each of dimension `2^(n-1)`.

**[LEMMA — correct reason for subtraction].** When `n` is odd, `N=2 mod 4`, so `omega` belongs formally to the `2 mod 4` class.  It is absent from every generator: bipartiteness forbids the endpoint chord because its distance `n-1` is even, so every chord grade is at most `N-4`.  It also cannot be produced by a bracket.  If `I symmetric-difference J` were all `N` indices, then `J=I^c`, so `I intersection J` would be empty; even monomials on disjoint complementary supports commute.  Hence the generated Lie algebra omits `omega`. `[]`

**[THEOREM — dimension convention].** This repository counts linearly independent anti-Hermitian Pauli strings in the generated Lie algebra on the full physical space.  It excludes the identity and excludes global parity because parity is not generated, not because parity is a scalar or because it is traceful.  This convention gives `65535`, not `65536`, for the `3x3` layer.

**[LEMMA — binomial sum].** Let

```text
S_N = sum_(k=2 mod 4) C(N,k).
```

A fourth-root-of-unity filter gives

```text
S_N = (2^N-(1+i)^N-(1-i)^N)/4
    = 2^(N-2)-2^(N/2-1) cos(N*pi/4).
```

With `N=2n`, this is

```text
2^(2n-2)-(-1)^(n/2)2^(n-1),   n even,
2^(2n-2),                       n odd.
```

Subtracting the nongenerated volume in the odd case gives `2^(2n-2)-1`. `[]`

**[EXTERNAL — structural identification, not needed for generation].** Over an algebraic closure, Clifford reversal acts on an even grade-`k` monomial by `(-1)^(k(k-1)/2)`, so its `-1` eigenspace is exactly the `k=2 mod 4` span.  For even `n`, its action on the two half-spin spaces yields the familiar `so+so` case when `n=0 mod 4` and `sp+sp` case when `n=2 mod 4`; for odd `n`, removing the central volume leaves the corresponding `sl(2^(n-1))` algebra.  The elementary support proof above establishes the exact vector space without relying on this representation-theoretic identification.

## 6. Obligation 6: containment, independently of generation

**[LEMMA — status: proved].** Let `A_N` be the span of all coordinate Clifford monomials of grade `2 mod 4`, with `omega` omitted when `N=2 mod 4`.  Then `A_N` is a Lie subalgebra.

**[LEMMA — proof].** If two basis monomials in `A_N` have grades `a=b=2 mod 4`, a nonzero bracket requires odd overlap `j` and has grade

```text
a+b-2j = 2 mod 4.
```

The bracket cannot equal `omega` by the complementary-support argument in Section 5.  Thus the bracket remains in `A_N`.  Section 1 and bipartiteness put every local generator in `A_N`, so

```text
g_Gamma subseteq A_N.
```

This containment uses no orbit-saturation or raising argument. `[]`

**[THEOREM — reverse inclusion only on the branching branch].** If `Delta(Gamma)>=3`, Section 0 supplies a chord grade `6<=k<=N-4`; Section 4 supplies its whole coordinate orbit; and Section 3 supplies every grade of `A_N`.  Therefore `A_N subseteq g_Gamma`, independently completing equality.  For a cycle, the only chord grade is the exceptional `N-2`, and the reverse inclusion stops exactly at `{2,N-2}`.

## 7. Obligation 7: exact verification

**[COMPUTATION — exact BFS method].** For monomial Pauli generators, every Lie word is a linear combination of left-normed words.  The exact finite closure can therefore be found by starting from the raw packed generators and repeatedly applying `v -> v XOR g` only when the binary symplectic product with a raw generator `g` is one.  Distinct packed Pauli strings are linearly independent, so the number reached is the exact Lie dimension.  The producer and verifier contain separate implementations and the verifier imports none of `e178`--`e181` or the repository Clifford package.

**[COMPUTATION — known values, all raw closures].** The independent verifier obtains:

| tag | graph | `n` | exact grade histogram | exact dimension |
|---|---|---:|---|---:|
| **[COMPUTATION]** | `2x2=C4` | 4 | `2:28, 6:28` | 56 |
| **[COMPUTATION]** | `2x3` | 6 | `2:66, 6:924, 10:66` | 1056 |
| **[COMPUTATION]** | `2x4` | 8 | `2:120, 6:8008, 10:8008, 14:120` | 16256 |
| **[COMPUTATION]** | `2x5` | 10 | `2:190, 6:38760, 10:184756, 14:38760, 18:190` | 262656 |
| **[COMPUTATION]** | `3x3` | 9 | `2:153, 6:18564, 10:43758, 14:3060` | 65535 |

**[COMPUTATION — trichotomy controls].** The same raw engine gives `P_6:66` with grade set `{2}`, `C_6:132` with histogram `{2:66,10:66}`, and `C_8:240` with histogram `{2:120,14:120}`.  The verifier additionally checks every pairwise bracket for the two cycle closures.

**[COMPUTATION — three requested new graphs].** The new cases are:

| tag | graph | method | chord grades | result |
|---|---|---|---|---:|
| **[COMPUTATION]** | `3x3` minus edge `(0,3)` | independent raw Pauli closure | contains grade 6 | 65535 |
| **[THEOREM]** | open `3x4` grid | non-materialized grade characterization plus orbit sample | `14,10,6,6,10,14` | 4192256 |
| **[THEOREM]** | `C_10` plus bipartite chord `(0,3)` | non-materialized grade characterization plus orbit sample | `18,6` | 262656 |

**[THEOREM — large-case certification scope].** The `3x4` row does not materialize 4,192,256 strings.  Its snake path is checked edge by edge; all chords have grade `2 mod 4`; a displayed grade-six chord invokes the proved corrected generation lemma; an exact Johnson-walk sample verifies coordinate-orbit action; and the dimension is independently summed from the allowed binomial grades.  The `C_10`-plus-chord row uses the same method, with the grade-six chord providing saturation while the endpoint grade-18 chord is redundant.

**[COMPUTATION — non-bipartite check].** The required `C_5`-plus-chord computation gives dimension `510`, versus `255` for the bipartite grade-class prediction and `510` for the full noncentral even-grade alternative, as detailed in Section 2.

## 8. Obligation 8: consequence, restricted to what is proved

**[LEMMA — Gaussian dimension bound].** Majorana bilinears on `m` fermionic modes span `so(2m)`, of dimension

```text
m(2m-1).
```

Therefore any realization of a generated algebra of dimension `D` by `m`-mode bilinears must satisfy

```text
m >= m_min(D) := ceil((1+sqrt(1+8D))/4).
```

**[THEOREM — branching consequence].** On the `Delta>=3` branch,

```text
D = Theta(4^n),       m_min(D)=Theta(2^n).
```

Thus the local terms are not simultaneously Majorana-quadratic after any unitary change of basis on the physical `n`-mode Fock space, and any bilinear realization allowed to enlarge the mode space needs exponentially many modes.

**[THEOREM — universal exponential claim refuted].** The same conclusion is false for the cycle counterfamily.  An even cycle has

```text
D=2n(2n-1)=Theta(n^2),
```

so the dimension argument gives only `m_min=Theta(n)`.  This is the ring/free-fermion-solvable branch, naturally resolved into parity sectors or represented after only a linear-size enlargement.  The open path is exactly `so(2n)` and is quadratic in `n` modes.

**[THEOREM — full-physical-space nuance].** Although a ring is free-fermion solvable sector by sector, its full-space algebra has dimension `2n(2n-1)>n(2n-1)`.  Hence no single unitary basis change on the full `n`-mode physical Fock space puts all ring local terms into one `so(2n)` bilinear algebra.  This statement is compatible with, and more precise than, sectorwise solvability.

**[THEOREM — comparison with the claw theorem].** The existing frustration-graph induced-claw theorem already gives a term-wise yes/no bilinearity obstruction whenever `Delta(Gamma)>=3`.  The Clifford-grade route adds exact Lie dimensions, an exponential quantitative mode lower bound on that branch, and the precise location of the dimension jump.  It also computes the `Delta<=2` controls, including the `2x2=C4` algebra of dimension `56`, where the claw argument does not apply.  For cycles it gives an exact quadratic-size algebra and only a linear mode count; it does not falsely promote that branch to the exponential no-go.

## Final obligation ledger

| tag | obligation | status |
|---|---:|---|
| **[LEMMA]** | 1 | proved, with exact phases |
| **[LEMMA]+[COMPUTATION]** | 2 | proved and tested non-bipartitely |
| **[THEOREM]** | 3 | refuted as stated; corrected `k<=N-4` lemma proved |
| **[LEMMA]** | 4 | proved over `Q`, including middle grade |
| **[THEOREM]** | 5 | requested scalar rationale refuted; correct central-generation argument proved |
| **[LEMMA]** | 6 | containment proved independently |
| **[COMPUTATION]** | 7 | all required and additional cases completed |
| **[THEOREM]** | 8 | exponential only for branching; cycle counterclaim exactly refuted |

**[UNRESOLVED].** No obligation in the corrected trichotomy remains open.  The two assertions left false are precisely the original endpoint-grade generation lemma and the resulting universal exponential-mode conclusion; both are replaced above by exact counterexamples and branchwise theorems.
