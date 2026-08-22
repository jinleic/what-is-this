# Beyond the Hamiltonian trichotomy: exact boundary, non-bipartite branch, and box corollaries

**[THEOREM — headline].** The Hamiltonian-path hypothesis in `proofs/clifford_grade_classification.md` cannot be deleted for connected bipartite graphs: the claw `K_{1,3}` has `Delta=3` but its individual-local-term Lie algebra has exact dimension `72`, not the Hamiltonian branching value `56`.  On the other hand, the non-bipartite branch admits a Hamiltonian-free classification: a connected simple non-bipartite graph is either an odd cycle, with dimension `2n(2n-1)`, or has `Delta>=3`, in which case its local-term algebra is the full noncentral even Pauli/Clifford algebra and has dimension `2^(2n-1)-2`.

**[COMPUTATION].** The integrated exact artifact is `results/algebra_growth/trichotomy_extend.json`.  Its producers are `experiments/e185_order_counterexample.py`, `experiments/e186_nonbipartite.py`, and `experiments/e187_trichotomy_extend.py`.  The clean-room verifier `tests/test_trichotomy_extend.py` imports none of those producers and starts every finite closure from raw `X_v` and `Z_uZ_v` labels.

**[THEOREM — scope].** Throughout, for a finite simple graph `Gamma` on `n` vertices,

```text
g_Gamma = Lie_Q { i X_v, i Z_u Z_v : v in V(Gamma), uv in E(Gamma) }.
```

The same dimension statements hold over `R`.  This is the many-generator **local-term** algebra.  It is not the two-generator algebra `Lie<A=sum_v X_v,B=sum_e Z_uZ_v>`.

## 1. What an arbitrary Jordan--Wigner order does

**[LEMMA — arbitrary-order edge grade].** Let `sigma(0),...,sigma(n-1)` be any linear order of the vertices and define

```text
gamma_(2r)   = product_(s<r) X_(sigma(s)) Z_(sigma(r)),
gamma_(2r+1) = product_(s<r) X_(sigma(s)) Y_(sigma(r)).
```

If `u=sigma(p)`, `v=sigma(q)`, and `p<q`, then

```text
Z_u Z_v = i^(q-p) gamma_(2p+1) gamma_(2p+2) ... gamma_(2q).
```

Thus the edge is one coordinate Clifford monomial of exact grade `2(q-p)`.

**[LEMMA — proof].** Relabeling the phase calculation in `proofs/clifford_grade_classification.md` Section 1 gives

```text
gamma_(2r+1) gamma_(2r+2) = -i Z_(sigma(r)) Z_(sigma(r+1)).
```

Multiplying from `r=p` through `q-1` cancels every intermediate `Z` and gives `(-i)^(q-p) Z_uZ_v`; rearranging proves the displayed identity.  No adjacency assumption on consecutive ordered vertices is used. `[]`

**[THEOREM — exact order criterion].** For a connected graph, every edge has grade `2 mod 4` in the order `sigma` if and only if

```text
c_sigma(v) = position_sigma(v) mod 2
```

is a proper two-colouring.  Consequently, for a connected bipartite graph with colour-class sizes `r` and `s`, such an order exists if and only if `|r-s|<=1`.

**[LEMMA — proof and count].** An edge `uv` has grade `2 mod 4` exactly when its endpoint positions differ oddly, which is exactly `c_sigma(u) != c_sigma(v)`.  A connected bipartite graph has a unique bipartition up to swapping colours.  The even and odd positions have sizes `ceil(n/2)` and `floor(n/2)`, so they can carry the two colour classes exactly when `|r-s|<=1`.  The number of compatible orders is

```text
2(r!)(s!)  if r=s,
(r!)(s!)   if |r-s|=1,
0          otherwise.
```

In the unequal case the larger class is forced onto the larger position parity. `[]`

**[LEMMA — two notions of a working order].** A *grade-compatible* order only makes every edge grade `2 mod 4`.  The wave-18 saturation proof needs more: every consecutive pair in the order must itself be an edge, so that the fields and those `n-1` bonds supply all consecutive Majorana bilinears.  Such an order is precisely an oriented Hamiltonian-path order.  A grade-compatible order need not be Hamiltonian and does not by itself generate the full grade-two component.

**[THEOREM — order independence of the algebra].** The matrix Lie algebra `g_Gamma` and its dimension are independent of the Jordan--Wigner order.  Changing `sigma` changes only the Clifford coordinates used to express the fixed physical Pauli matrices `X_v` and `Z_uZ_v`; it neither adds nor removes a physical generator.  Individual displayed Clifford grades are therefore presentation data, whereas the raw Pauli closure and its dimension are intrinsic.

**[COMPUTATION — exhaustive finite order audit].** The producer checks every connected bipartite labelled graph through `n=5`, all `218` graphs and every vertex order, and obtains exactly the colour-class count above in all cases.  The standalone verifier independently repeats the exhaustive audit through `n=4` (`23` graphs).

## 2. The Hamiltonian-path hypothesis is genuinely necessary

**[THEOREM — exact counterexample].** Let `Gamma=K_{1,3}`.  It is connected, simple, bipartite, and has `Delta=3`, but it has no Hamiltonian path.  Its bipartition sizes are `1` and `3`, so Section 1 also proves that no linear order can put every bond in grade `2 mod 4`.  Exact local-term closure gives

```text
dim g_(K1,3) = 72,
```

whereas blindly applying the Hamiltonian branching formula at `n=4` gives `56`.  Therefore the wave-18 bipartite trichotomy does not extend to all connected bipartite graphs.

**[COMPUTATION — exact closure certificate].** Encode a Pauli monomial by `(a,b) in F_2^n x F_2^n`, with raw generators

```text
X_v       = (e_v,0),
Z_u Z_v   = (0,e_u+e_v).
```

For labels `p,q`, a nonzero commutator occurs exactly when their binary symplectic pairing is one, and then its label is `p+q`.  Generator-adjoint breadth-first closure reaches `72` distinct labels for `K_{1,3}`.  An independent all-pairs saturation check finds no missing nonzero bracket.  In the representative order `(0,1,2,3)`, the three raw bond grades are `2,4,6` and the reached grade histogram is

```text
grade 2: 16,   grade 4: 40,   grade 6: 16.
```

The producer and verifier independently obtain closure digest `1a170f73a6401c7318b7255289a83370f26d486af521fd6c239833b607288d2b`.

**[COMPUTATION — second control].** The six-vertex branch tree with edges

```text
01, 12, 13, 34, 45
```

is connected and bipartite with colour-class sizes `2` and `4`, has no Hamiltonian path and no grade-compatible order, and closes exactly at dimension `992`, not the misapplied Hamiltonian value `1056`.  Thus the claw mismatch is not the only measured non-Hamiltonian failure.

**[UNRESOLVED — exact remaining bipartite scope].** This front does not classify every connected bipartite graph without a Hamiltonian path.  The two exact counterexamples settle the requested removal attempt negatively; they do not license a replacement formula based only on `n`, `Delta`, or colour-class balance.

## 3. The full non-bipartite branch

### 3.1 The order-free upper algebra

**[LEMMA — parity container].** Put `P_X=product_v X_v`.  Every raw local term commutes with `P_X`.  In binary Pauli coordinates, its centraliser is the span of labels `(a,b)` with `|b|` even and therefore contains `2^(2n-1)` coordinate Pauli strings.

**[LEMMA — the two central labels are not generated].** Neither the identity label `0` nor the global parity label `(1,...,1;0)` can arise.  Equal Pauli labels commute, so a nonzero commutator cannot have label zero.  Global `X` parity is the radical of the symplectic form on the even-`b` label space: if `p+q=P_X`, then

```text
<p,q> = <p,p+P_X> = 0,
```

so that proposed bracket also vanishes.  Since neither label is raw, Pauli-basis independence excludes both from `g_Gamma`.

**[THEOREM — full noncentral even algebra].** Define

```text
E_n = span_Q { Pauli(a,b) : |b| even } minus span_Q{I,P_X}.
```

Then `E_n` is a Lie algebra of dimension

```text
dim E_n = 2^(2n-1)-2.
```

Under any Jordan--Wigner realization, it is exactly the span of all coordinate even Clifford grades `2,4,...,2n-2`; equivalently it is the derived, noncentral part of the even Clifford algebra.  Section 3.1 proves `g_Gamma subseteq E_n` for every connected graph with `n>=2`.

### 3.2 Two support-generation lemmas

**[LEMMA — an interior `0 mod 4` seed].** Let `N` be even.  Grade two together with one coordinate grade `h` satisfying

```text
h = 0 mod 4,    4 <= h <= N-4
```

generates every noncentral even coordinate grade.

**[LEMMA — proof].** The grade-two orbit/Johnson-graph lemma from `proofs/clifford_grade_classification.md` supplies the entire grade-`h` coordinate space.  Two `h`-subsets with overlap `h-3` have odd overlap, fit because their union has size `h+3<=N-1`, and bracket to grade six.  Grade six raises every present grade `r=2 mod 4` to `r+4` by an overlap-one bracket; this recovers every noncentral `2 mod 4` grade.

If `h=4`, grade four is already present.  If `h>=8`, the preceding step supplies grade `h-2`; an `h`-subset and an `(h-2)`-subset overlapping in `h-3` indices fit in `h+1<=N-3` indices and bracket to grade four.  Finally, overlap-one brackets of grade six with grades `4,8,12,...` raise by four until every noncentral `0 mod 4` grade is present.  Hence grades `2,4,...,N-2` all occur. `[]`

**[LEMMA — endpoint `0 mod 4` plus an interior `2 mod 4` seed].** Suppose `N=2 mod 4`.  Grade two, one coordinate grade `k` with

```text
k = 2 mod 4,    6 <= k <= N-4,
```

and one coordinate grade `N-2` generate every noncentral even grade.

**[LEMMA — proof].** The corrected wave-18 grade lemma gives every noncentral grade `r=2 mod 4`, and the grade-two orbit gives all coordinate `(N-2)`-monomials.  Choose an `(N-2)`-subset and an `r`-subset with overlap `r-1`.  The overlap is odd, the union has size `N-1`, and their bracket has grade

```text
(N-2)+r-2(r-1) = N-r.
```

As `r` runs through `2,6,...,N-4`, the complements `N-r` run through all grades `N-2,N-6,...,4`, which are exactly the noncentral `0 mod 4` grades. `[]`

**[COMPUTATION — support audit].** Independent finite support closures exhaust the first lemma in `361` `(N,h)` cases and the second in `171` `(N,k)` cases through `N=80`.  These finite checks audit the overlap arithmetic; the preceding support constructions, not the finite range, prove the lemmas for all `N`.

### 3.3 First the Hamiltonian case

**[THEOREM — Hamiltonian non-bipartite dichotomy].** Let `Gamma` be connected, simple, non-bipartite, and Hamiltonian on `n` vertices, and put `N=2n`.  Exactly one of the following holds:

| tag | graph branch | exact Clifford grades in a Hamiltonian-path order | exact dimension |
|---|---|---|---:|
| **[THEOREM]** | odd cycle `C_n` | `{2,N-2}` | `2*C(N,2)=2n(2n-1)` |
| **[THEOREM]** | `Delta(Gamma)>=3` | every even grade `2,4,...,N-2` | `2^(2n-1)-2` |

**[LEMMA — proof].** Path-position parity is the path's two-colouring.  Non-bipartiteness supplies a same-parity chord, hence a chord of even distance and grade `0 mod 4`.  If such a chord is not the endpoint pair, its grade lies between `4` and `N-4`, and the first support lemma gives `E_n`.

If the endpoint pair is the only same-parity chord, then that edge exists, its distance `n-1` is even, and therefore `n` is odd and `N=2 mod 4`.  When `Delta>=3`, the branching equivalence in the wave-18 note supplies another nonendpoint chord.  By assumption it joins opposite path parities, so its grade is `2 mod 4` and lies between `6` and `N-4`.  The second support lemma again gives `E_n`.

If `Delta<=2`, connectedness and non-bipartiteness force an odd cycle.  Its only chord relative to the opened Hamiltonian path is the endpoint edge of grade `N-2`.  The grade set `{2,N-2}` is bracket-closed: grade two preserves the `(N-2)` orbit, and two `(N-2)` supports can bracket nontrivially only back to grade two.  Both complete coordinate orbits occur, giving `2*C(N,2)`. `[]`

**[THEOREM — necessary odd-cycle exception].** The sentence “one `0 mod 4` edge forces the full even algebra” is false without the branching qualification.  The bare `C_5` has exact dimension `90` and grades `{2,8}`, whereas `dim E_5=510`.  The triangle `C_3` is an accidental small coincidence: `{2,4}` already comprises every noncentral even grade at `N=6`.

### 3.4 Removing Hamiltonicity in the non-bipartite branch

**[LEMMA — lollipop seed].** An odd cycle with one new leaf attached to a cycle vertex has full local-term algebra `E_m` on its `m` vertices.

**[LEMMA — proof].** If the odd cycle has length `ell`, order the vertices as the leaf, its cycle neighbour, and then the remaining cycle vertices around the cycle.  This is a Hamiltonian path.  The omitted closing-cycle edge has order distance `ell-1`, hence grade `2(ell-1)=N-4`, which is `0 mod 4`.  The first support lemma and the parity upper bound give equality with `E_m`. `[]`

**[LEMMA — leaf extension].** If a graph `H` on `m>=3` vertices has local-term algebra `E_m`, then adjoining a new leaf by one edge produces a graph whose local-term algebra is `E_(m+1)`.

**[LEMMA — proof].** The full algebra `E_m` contains every old field `X_v` and every old pair bond `Z_uZ_v`.  After adding the leaf field and its attachment bond, the generated algebra therefore contains the local-term algebra of `K_m` with a pendant leaf.  That comparison graph has a Hamiltonian path beginning at the leaf, is non-bipartite, and has maximum degree at least three, so Section 3.3 gives its full algebra `E_(m+1)`.  Generator monotonicity gives the lower inclusion, while Section 3.1 gives the matching upper inclusion. `[]`

**[LEMMA — spanning construction].** Every connected non-bipartite graph that is not itself an odd cycle has a spanning subgraph obtainable from an odd-cycle-plus-leaf seed by repeated leaf additions.

**[LEMMA — proof].** Choose an odd cycle `C`.  If a vertex lies outside `C`, connectedness supplies a path leaving `C`; retain its first edge to obtain the seed.  If `C` is spanning but the graph is not exactly `C`, choose an extra chord.  The chord and one of the two cycle arcs form a proper odd cycle, and the first edge of the omitted arc attaches a vertex outside that smaller odd cycle, again giving the seed.  Contract the seed and choose a spanning tree of the remaining connected graph.  Add its vertices in parent-before-child order, retaining only each parent edge; every step is a leaf addition. `[]`

**[THEOREM — general connected non-bipartite classification].** Let `Gamma` be any connected simple non-bipartite graph on `n>=3` vertices.  Then

```text
Gamma is an odd cycle
    => g_Gamma has grades {2,2n-2} and dim = 2n(2n-1);
Gamma is not an odd cycle (equivalently Delta(Gamma)>=3)
    => g_Gamma = E_n and dim = 2^(2n-1)-2.
```

**[LEMMA — proof].** The odd-cycle row is Section 3.3.  Otherwise the spanning construction, the lollipop seed, repeated leaf extension, and generator monotonicity give `E_n subseteq g_Gamma`; the parity container gives the reverse inclusion.  A connected simple graph of maximum degree at most two is a path or a cycle, and the only non-bipartite option is an odd cycle, proving the stated equivalence. `[]`

**[COMPUTATION — independent graph checks].** Exact raw-generator closures verify four different full-even controls and the odd-cycle exception:

| tag | graph | Hamiltonian? | decisive chord/structure | exact grade histogram | dimension |
|---|---|---|---|---|---:|
| **[COMPUTATION]** | bare `C5` | yes | endpoint grade `8` only | `2:45, 8:45` | `90` |
| **[COMPUTATION]** | paw = `C3` plus leaf | yes | interior grade `4` | `2:28, 4:70, 6:28` | `126` |
| **[COMPUTATION]** | `C5` plus chord `(0,3)` | yes | grades `8` and `6` | `2:45, 4:210, 6:210, 8:45` | `510` |
| **[COMPUTATION]** | triangle plus two leaves at one vertex | no | one leaf-extension step | `2:45, 4:210, 6:210, 8:45` | `510` |
| **[COMPUTATION]** | `C6` plus chord `(0,2)` | yes | interior grade `4` | `2:66, 4:495, 6:924, 8:495, 10:66` | `2046` |

**[COMPUTATION].** The producer pairwise-saturates all five closures.  The standalone verifier independently rebuilds all five from raw fields and bonds, re-inverts its own Jordan--Wigner basis, and reproduces every dimension and histogram.

## 4. Rectangular-grid corollaries and repository-status audit

**[THEOREM — all open rectangular grids].** Let `G_(r,s)=P_r square P_s` be the open `r x s` grid and put `n=rs`, `N=2n`.  A row snake is a Hamiltonian path and checkerboard parity is a bipartition.  Therefore:

```text
min(r,s)=1:
    grades {2}, dim = n(2n-1);
(r,s)=(2,2):
    grades {2,N-2}, dim = 2n(2n-1)=56;
r,s>=2 and (r,s)!=(2,2):
    every noncentral grade k=2 mod 4,
    dim = 2^(2n-2)-(-1)^(n/2)2^(n-1)  if n is even,
          2^(2n-2)-1                   if n is odd.
```

**[LEMMA — proof].** The snake visits adjacent horizontal vertices within a row and crosses vertically at each row end.  The grid is a path when one side is one, is the cycle `C4` at `2x2`, and otherwise has a degree-three or degree-four vertex.  The wave-18 exact bipartite Hamiltonian trichotomy applies in the three respective branches. `[]`

**[EXTERNAL — repository audit rule].** Before assigning a novelty status, this front searched the repository for each named geometry and exact dimension, then checked `proofs/alll_saturation.md`, `proofs/clifford_grade_classification.md`, `results/algebra_growth/clifford_grade.json`, and `proofs/char0_4x4.md`.  “Previously implicit” below means the wave-18 theorem already implied the value, not that this front proved a mathematically new classification.

**[THEOREM — concrete corollaries with prior status].** The exact local-term values and their pre-front repository status are:

| tag | open layer | `n` | exact local-term grade space | exact dimension | previously open / previously known |
|---|---:|---:|---|---:|---|
| **[THEOREM]** | `2x6` | 12 | noncentral `k=2 mod 4` | `4,192,256` | **previously known exactly**: full `L=6` closure in `proofs/alll_saturation.md`, `H431`--`H433` |
| **[THEOREM]** | `3x4` | 12 | noncentral `k=2 mod 4` | `4,192,256` | **previously known exactly**: explicit wave-18 row in `proofs/clifford_grade_classification.md`, `H486` |
| **[THEOREM]** | `3x5` | 15 | noncentral `k=2 mod 4`, top volume omitted | `268,435,455` | **previously implicit, not graph-specifically tabled**: wave-18 formula plus stored `n=15` binomial check |
| **[THEOREM]** | `4x4` | 16 | noncentral `k=2 mod 4` | `1,073,709,056` | **previously implicit for local terms**: wave-18 formula plus stored `n=16` check; the explicitly open `>=1,794` object was the different two-generator algebra |

**[THEOREM — the `4x4` scope correction].** `proofs/char0_4x4.md` defines

```text
A=sum_v X_v,   B=sum_e Z_uZ_v,   g=Lie_Q<A,B>
```

and proves only `dim g>=1794`.  The `1,073,709,056` value above instead belongs to the algebra generated by every individual `iX_v` and `iZ_uZ_v`.  This front neither improves the two-generator lower bound nor proves its equality, type, centre, or radical.

**[COMPUTATION — raw corollary certificates].** For every concrete row, the artifact stores the raw field/bond generator digest, the complete snake path, all chord grades, the degree, and both the direct binomial sum and closed-form dimension.  The standalone verifier reconstructs all of them.  In particular, it derives from raw graph generators

```text
3x5 chord grades = (18,14,10,6,6,10,14,18),
4x4 chord grades = (14,10,6,6,10,14,14,10,6),
```

finds a valid interior grade-six seed in each, and independently obtains `268435455` and `1073709056`.  It rejects any artifact corollary row lacking a recognized raw graph/path/grade certificate or a sourced previous-status field.

## 5. Simple-cubic finite boxes

**[LEMMA — explicit box Hamiltonian path].** Let `B_(a,b,c)=P_a square P_b square P_c`.  Snake through one `a x b` layer.  In successive `c`-layers alternate that two-dimensional snake with its reversal.  The end of one layer and the start of the next then have identical in-layer coordinates and differ by exactly one vertical edge.  This constructs a Hamiltonian path for every positive `a,b,c`.

**[THEOREM — exact finite-box local-term algebra].** Put `n=abc` and `N=2n`.  If `a,b,c>=2`, the box is bipartite, Hamiltonian, and has `Delta>=3`; hence

```text
g_(a,b,c)
  = span_Q of all coordinate Clifford monomials of grade k=2 mod 4,
    with the central top grade N omitted when n is odd,

dim g_(a,b,c)
  = 2^(2n-2)-(-1)^(n/2)2^(n-1)  if n is even,
    2^(2n-2)-1                   if n is odd.
```

**[THEOREM — degenerate boxes].** The same explicit path gives the complete boundary cases.  If at most one side exceeds one, the box graph is a path and has dimension `n(2n-1)`.  If the multiset of sides exceeding one is `{2,2}`, it is `C4` and has dimension `56`.  Every other open box is in the branching formula above.

**[COMPUTATION — finite geometry audit].** The producer and verifier independently build every box with sides at most four and volume at most `24`: `53` boxes total, including `40` branching cases.  Every displayed path is edge-valid, every graph is connected and bipartite, and every branching chord census lies in grade `2 mod 4`.  This finite audit checks the construction; the explicit alternating-layer proof establishes it for all side lengths.

**[THEOREM — what the box corollary does not mean].** This is an exact statement about the Lie algebra generated by the individual local terms of one finite open box.  It is **not** a statement about the thermodynamic limit, the transfer spectrum, the two-generator algebra, or solvability of the three-dimensional Ising model.  In particular, it does not provide a partition function, critical point, spectrum, or 3D free-energy formula.

## 6. Verification ledger and exact limits

**[COMPUTATION — reproduction].** From the repository root, the producer and standalone verifier are

```text
.venv/bin/python experiments/e187_trichotomy_extend.py
.venv/bin/python tests/test_trichotomy_extend.py
```

The verifier prints `PASS`, independently re-derives all five non-bipartite closures, the two bipartite counterexamples, every concrete grid corollary, and the general finite box path audit.  It imports no producer conclusion.

**[COMPUTATION — resources].** The integrated producer records a peak RSS below `27 MB`; the standalone verifier records a peak below `25 MB`.  Both enforce a `2,000,000,000`-byte cap and budget loops with `time.process_time()`.  These measurements are host-specific resource observations, not mathematical theorems.

**[THEOREM — exact achieved scope].** The achieved all-size statements are the arbitrary-order criterion, the impossibility of deleting Hamiltonicity from the bipartite wave-18 formula, the Hamiltonian-free connected non-bipartite cycle/full-even dichotomy, the full open-grid corollary, and the finite simple-cubic-box local-term corollary.

**[UNRESOLVED — excluded problems].** The local-term classification does not settle the remaining connected bipartite non-Hamiltonian graphs, any two-generator all-size problem, any transfer-spectrum problem, any thermodynamic limit, or the 3D Ising model.  No finite closure or finite-box Lie dimension is promoted here to a solution of those problems.
