# Boundary matchgate defect for the open `2x3x3` box

Artifacts:

- producer: `experiments/e247_boundary_matchgate_defect.py`;
- clean-room verifier: `tests/test_boundary_matchgate_defect.py`;
- producer output: `results/integrability/boundary_matchgate_defect.json`;
- inherited topology source: `experiments/e241_box_genus_bound.py`.

## 1. Exact status

`[THEOREM]` below means a self-contained algebraic conclusion. `[COMPUTATION]`
means a bounded integer, rational, or polynomial enumeration whose complete
coverage is recorded in the artifact and independently reconstructed by the
verifier.

Let

\[
G=P_2\mathbin\square P_3\mathbin\square P_3
\]

have open boundary conditions, let every vertex be colored by
\(x+y+z\pmod 2\), and set \(v=\tanh K\).

**[THEOREM — fixed-boundary matchgate defect].** For each of the two choices
of checkerboard color to eliminate, retain the five sites of the opposite
color on the unique `3x3` \(x\)-face containing five such sites and average
the other four retained-color spins. In the Walsh basis and after normalizing
by the empty Walsh entry, the resulting five-leg signature is not the tensor
of even principal Pfaffians of the pair matrix for **any** of the \(5!=120\)
boundary-leg permutations, for every

\[
0<v<1.
\]

This is an all-\(v\) statement, not merely a check at selected rational
points. It concerns exactly this open-box contraction, this Walsh basis, and
the complete ordering class of its five displayed legs.

The theorem is deliberately narrower than an arbitrary-basis or global
statement. It is **not** a no-go for arbitrary hidden auxiliaries, is **not** a
global Pfaffian-count lower bound, does **not** cover all box sizes or all
lengths \(L\), and is **not** a solution of the three-dimensional Ising model.

## 2. Explicit geometry and the two boundary choices

Use coordinates

\[
(x,y,z)\in\{0,1\}\times\{0,1,2\}\times\{0,1,2\}.
\]

Two vertices are adjacent exactly when one coordinate changes by one. Thus
there are 18 vertices and

\[
(2-1)3\cdot3+2(3-1)3+2\cdot3(3-1)=9+12+12=33
\]

edges. Every edge flips \(x+y+z\pmod2\). Direct coordinate enumeration gives
nine even and nine odd vertices, so this is an explicitly checked \(9+9\)
bipartition.

The two retained boundary orders are as follows.

1. Eliminate the even color and retain the odd color. The odd-color counts on
   \(x=0\) and \(x=1\) are respectively 4 and 5. Select \(x=1\), with leg
   order

   \[
   (1,0,0),(1,0,2),(1,1,1),(1,2,0),(1,2,2).
   \]

   The averaged odd sites are

   \[
   (0,0,1),(0,1,0),(0,1,2),(0,2,1).
   \]

2. Eliminate the odd color and retain the even color. The even-color counts on
   \(x=0\) and \(x=1\) are respectively 5 and 4. Select \(x=0\), with leg
   order

   \[
   (0,0,0),(0,0,2),(0,1,1),(0,2,0),(0,2,2).
   \]

   The averaged even sites are

   \[
   (1,0,1),(1,1,0),(1,1,2),(1,2,1).
   \]

In both cases the nine eliminated stars have degree histogram

\[
3^4,\quad 4^4,\quad 5^1.
\]

The producer records every center coordinate and its neighbors in increasing
coordinate order. Reflection in the middle of the \(x\)-direction identifies
the two resulting polynomial signatures, but the program constructs and
gates them separately before checking that equality.

## 3. Checkerboard elimination and the boundary Walsh polynomial

For signs \(s_i\in\{-1,1\}\), eliminating one star center gives the normalized
all-degree factor

\[
\begin{aligned}
f_d(s_1,\ldots,s_d)
 &=\frac{\prod_i(1+v s_i)+\prod_i(1-v s_i)}{2}\\
 &=\sum_{\substack{A\subseteq[d]\\ |A|\text{ even}}}
       v^{|A|}\prod_{i\in A}s_i. \tag{3.1}
\end{aligned}
\]

Let \(U\) be the nine eliminated sites, \(R\) the nine retained sites,
\(B\subset R\) the five boundary legs, and \(I=R\setminus B\) the four
averaged sites. For \(A\subseteq B\), define the unnormalized boundary Walsh
entry

\[
B_A(v)=2^{-9}\sum_{\sigma_R\in\{-1,1\}^R}
 \left(\prod_{b\in A}\sigma_b\right)
 \prod_{u\in U}f_{\deg(u)}(\sigma_{N(u)}). \tag{3.2}
\]

The factor \(2^{-9}=2^{-4}2^{-5}\) combines the average over \(I\) and the
Walsh transform over \(B\).

Expanding (3.1) gives an equivalent integer edge-boundary formula:

\[
B_A(v)=
\sum_{\substack{F\subseteq E:\\
 \deg_F(u)\equiv0\ (u\in U),\\
 \deg_F(i)\equiv0\ (i\in I),\\
 \{b\in B:\deg_F(b)\equiv1\}=A}}
 v^{|F|}. \tag{3.3}
\]

Indeed, the local expansion at each \(u\) chooses an even subset of its
incident edges. Averaging a retained spin kills precisely the terms with odd
incidence there, while the boundary Walsh character selects precisely the
boundary syndrome \(A\). This also shows that every odd-cardinality boundary
entry is zero.

The producer evaluates (3.3) by a parity-state dynamic program over the nine
stars. As an independent semantic control, on the open `2x2x2` cube it compares
all \(2^4=16\) retained-spin assignments in (3.2) with all \(2^{12}=4096\)
edge subsets in (3.3). Every boundary polynomial agrees exactly. The
clean-room verifier rebuilds the target `2x3x3` signatures by the direct spin
sum (3.2), not by importing or replaying the producer's star dynamic program.

The shared empty entry of the two target signatures is

\[
\begin{aligned}
D(v)=B_\varnothing(v)={}&1+20v^4+78v^6+402v^8+1640v^{10}
 +4909v^{12}+11636v^{14}\\
&+17685v^{16}+16840v^{18}+9374v^{20}+2542v^{22}
 +376v^{24}+32v^{26}+v^{28}. \tag{3.4}
\end{aligned}
\]

Every coefficient is nonnegative and the constant term is one, hence
\(D(v)>0\) for \(0<v<1\). The normalized signature is

\[
g_A(v)=\frac{B_A(v)}{D(v)},\qquad g_\varnothing=1. \tag{3.5}
\]

All 32 raw entries, including the 16 exact zeros in odd parity, are stored as
coefficient lists in the artifact.

## 4. Why the five-leg Pfaffian test is complete

Fix a boundary permutation

\[
\pi=(\pi_0,\pi_1,\pi_2,\pi_3,\pi_4)\in S_5.
\]

Here `pi` is an ordinary unsigned relabeling of the set-indexed Walsh tensor.
It does not attach an active fermionic-SWAP or Koszul phase to a subset.
Phase-decorated fermionic rewiring is a different operation and lies outside
this fixed-Walsh/order theorem.

If the normalized even signature is a principal-sub-Pfaffian tensor in this
order, its skew pair matrix is forced by its two-leg entries:

\[
M^{(\pi)}_{ij}=g_{\{\pi_i,\pi_j\}}\quad(i<j),\qquad
M^{(\pi)}_{ji}=-M^{(\pi)}_{ij}. \tag{4.1}
\]

There is no freedom left in \(M\): every \(2\times2\) principal Pfaffian is
its corresponding upper-triangular entry. For positions \(a<b<c<d\), the
remaining condition is

\[
g_{\{\pi_a,\pi_b,\pi_c,\pi_d\}}=
 M_{ab}M_{cd}-M_{ac}M_{bd}+M_{ad}M_{bc}. \tag{4.2}
\]

A five-leg even signature has only subsets of sizes 0, 2, and 4. Therefore its
one empty, ten pair, and five four-leg equations are all of its even principal
Pfaffian equations. Checking all 16 for every \(\pi\in S_5\) is complete.

Clear the positive denominator in (4.2):

\[
\begin{aligned}
R_{\pi;abcd}(v)={}&D(v)B_{\{\pi_a,\pi_b,\pi_c,\pi_d\}}(v)\\
&-\bigl(
 B_{\{\pi_a,\pi_b\}}B_{\{\pi_c,\pi_d\}}
 -B_{\{\pi_a,\pi_c\}}B_{\{\pi_b,\pi_d\}}
 +B_{\{\pi_a,\pi_d\}}B_{\{\pi_b,\pi_c\}}
 \bigr). \tag{4.3}
\end{aligned}
\]

The normalized residual is \(R_{\pi;abcd}/D^2\), so it has the same sign as
\(R_{\pi;abcd}\) on the stated interval.

**[COMPUTATION — complete coverage].** For each color choice the producer
checks

\[
120\times(1+10+5)=1920
\]

even principal Pfaffians symbolically. The 1320 empty-and-pair residuals are
identically zero by construction. Every one of the 600 four-leg residuals is
nonzero. Thus every ordering has exactly five symbolic failures. The verifier
repeats the 120-order enumeration by factoradic unranking and uses the explicit
four-term formula (4.3), rather than the producer's recursive polynomial
Pfaffian routine.

## 5. Exact factorization and strict sign on `0<v<1`

Order the five four-position subsets lexicographically, so the first is
\((0,1,2,3)\). It is already nonzero for every one of the 120 permutations.
Across all first residuals, and also across all 600 residual occurrences, only
four polynomials occur. With each \(Q_j\) irreducible over \(\mathbb Z[v]\),
the complete positive-interval factorizations are

\[
\begin{aligned}
R_{4,36}(v)&=4v^4(1-v)^6(1+v)^6(1+v^2)^2Q_{36}(v),\\
R_{4,38}(v)&=2v^4(1-v)^6(1+v)^6(1+v^2)Q_{38}(v),\\
R_{6,28}(v)&=4v^6(1-v)^8(1+v)^8(1+v^2)^2Q_{28}(v),\\
R_{10,20}(v)&=8v^{10}(1-v)^{11}(1+v)^{11}(1+v^2)Q_{20}(v). \tag{5.1}
\end{aligned}
\]

The remainder polynomials are

\[
\begin{aligned}
Q_{36}(v)={}&1+9v^2+86v^4+508v^6+2389v^8+8704v^{10}
 +25555v^{12}+60683v^{14}\\
&+115407v^{16}+178677v^{18}+224201v^{20}+214498v^{22}
 +137351v^{24}\\
&+59430v^{26}+17373v^{28}+3287v^{30}+388v^{32}+28v^{34}+v^{36},\\[2mm]
Q_{38}(v)={}&1+17v^2+190v^4+1484v^6+7961v^8+34263v^{10}
 +115785v^{12}+315551v^{14}\\
&+692059v^{16}+1211057v^{18}+1691429v^{20}+1802207v^{22}
 +1394971v^{24}\\
&+751317v^{26}+283139v^{28}+72885v^{30}+12784v^{32}
 +1426v^{34}+81v^{36}+v^{38},\\[2mm]
Q_{28}(v)={}&1+30v^2+199v^4+990v^6+3395v^8+8616v^{10}
 +17510v^{12}+27180v^{14}\\
&+30691v^{16}+24110v^{18}+13079v^{20}+4406v^{22}
 +785v^{24}+76v^{26}+4v^{28},\\[2mm]
Q_{20}(v)={}&25+203v^2+874v^4+2266v^6+3718v^8+4018v^{10}
 +2734v^{12}\\
&+1178v^{14}+293v^{16}+47v^{18}+4v^{20}. \tag{5.2}
\end{aligned}
\]

The ordering inventory is:

| residual class | one representative permutation | first residuals | all four-subset occurrences |
|---|---:|---:|---:|
| \(R_{4,36}\) | `(0,1,2,4,3)` | 64 | 320 |
| \(R_{4,38}\) | `(0,1,3,4,2)` | 16 | 80 |
| \(R_{6,28}\) | `(0,1,2,3,4)` | 32 | 160 |
| \(R_{10,20}\) | `(0,1,4,3,2)` | 8 | 40 |
| total | — | 120 | 600 |

The producer factors each residual over \(\mathbb Z[v]\), reconstructs the
unfactored coefficient list from the factor list, and checks the displayed
multiplicities. The verifier independently strips the elementary factors,
reconstructs the positive remainder, and uses exact \(\mathbb Z[v]\)
factorization to check its irreducibility.

Most importantly, the nonvanishing proof needs no numerical root finder and no
floating-point sign inference. Every scalar in (5.1) is positive; on
\(0<v<1\), each of \(v\), \(1-v\), \(1+v\), and \(1+v^2\) is positive; and
every coefficient in every polynomial (5.2) is strictly positive. Exact
coefficient positivity therefore gives

\[
R_{\pi;abcd}(v)>0
\quad\text{for every tested }(\pi;abcd)\text{ and every }0<v<1. \tag{5.3}
\]

Together with \(D(v)>0\), equations (4.2)--(4.3) prove the theorem. There are
no exceptional roots inside the interval to report. The endpoint zeros caused
by powers of \(v\) and \(1-v\) are outside the open domain.

## 6. Exact `v=1/3` and `v=1/2` controls

The rational points are controls, not the basis of the all-\(v\) proof. The
empty entries are

\[
D(1/3)=\frac{33287590051840}{22876792454961},\qquad
D(1/2)=\frac{2387190145}{268435456}. \tag{6.1}
\]

For each color and each point, all 120 pair matrices are reconstructed from the
normalized two-leg entries and all 1920 even principal Pfaffians are evaluated
with `Fraction`. Exactly 1320 entries match (the empty and pair entries) and
all 600 four-leg entries fail, with five failures in every ordering.

For additional exact scale information, the four normalized positive
residuals \(R/D^2\) are:

| class | at \(v=1/3\) | at \(v=1/2\) |
|---|---:|---:|
| \(R_{4,36}\) | \(635276288146432/10319646927149569\) | \(9428963830002468/227947071535404841\) |
| \(R_{4,38}\) | \(3089366431234048/51598234635747845\) | \(68017283373030714/1139735357677024205\) |
| \(R_{6,28}\) | \(113895401877504/10319646927149569\) | \(1820091846874176/227947071535404841\) |
| \(R_{10,20}\) | \(62354588663808/51598234635747845\) | \(1099489220642304/1139735357677024205\) |

Every number in the table is positive, consistently with (5.3).

## 7. Planar alternating-cycle positive control

A failure result is meaningful only if the sign and ordering convention also
accepts a known planar composition. Take an alternating 10-cycle

\[
b_0,u_0,b_1,u_1,b_2,u_2,b_3,u_3,b_4,u_4,b_0,
\]

where each eliminated \(u_i\) is adjacent to boundary spins \(b_i\) and
\(b_{i+1\bmod5}\). Eliminating the five degree-two stars gives

\[
\prod_{i=0}^4(1+v^2 b_i b_{i+1}). \tag{7.1}
\]

Equivalently, select any subset of the five effective cycle edges; the
boundary mask is its mod-two boundary and its weight is \(v^{2|F|}\). The
empty Walsh entry is

\[
D_{C_{10}}(v)=1+v^{10}>0. \tag{7.2}
\]

The producer feeds this signature to the **same** permutation, pair-matrix,
and principal-Pfaffian implementation as the target. Exactly the ten rotations
and reversals of `(0,1,2,3,4)` pass all 16 equations symbolically. The other
110 orders demonstrate the expected dependence on boundary order. At both
\(v=1/3\) and \(v=1/2\), the complete 1920-entry audit has 1520 matches and
400 failures, with the same ten fully passing orders. The clean-room verifier
reconstructs (7.1) by a direct 32-state spin sum.

Thus the target failure is not an implementation that rejects every input or
a convention with the Pfaffian signs reversed.

## 8. Inherited nonplanarity provenance and terminology

The graph count and nonplanarity statement are inherited from the current
source of `e241`, whose SHA-256 hash is stored and gated. For
\((a,b,c)=(2,3,3)\), the `e241` lower-bound numerator is

\[
abc-ab-bc-ca+4=18-6-9-6+4=1,
\]

so its orientable-genus lower bound is

\[
\left\lceil\frac14\right\rceil=1.
\]

This proves nonplanarity. It does **not** determine the exact orientable genus,
and no exact-genus claim is made here. The boundary calculation above is an
independent finite contraction on that graph; it is not presented as a proof
that nonplanarity is the unique cause of the residual.

The box is open. No periodic seams, seam signs, or surface-sector sum is used.
In particular, flat graph-\(H^1\) cochains are not identified with surface
spin structures here. The proof uses only the explicit edge set, the
checkerboard stars, and the five displayed boundary legs.

Recent work on fermionic tensor-network contraction makes global ordering and
swap structure explicit even on arbitrary geometries
([arXiv:2410.02215v3](https://arxiv.org/abs/2410.02215v3)). The
Gaussian-plus-interaction viewpoint in
[arXiv:2412.04216v2](https://arxiv.org/abs/2412.04216v2) and the
Gaussian-plus-universal-quartic decomposition of parity-preserving circuits in
[arXiv:2504.19317v1](https://arxiv.org/abs/2504.19317v1) motivated looking for
the first explicit residual after locally Gaussian stars are composed.
[arXiv:2505.07804v2](https://arxiv.org/abs/2505.07804v2) provides a separate
topological-complexity perspective on matchgate simulation. These papers
motivate the question only: equations (3.1)--(5.3) are a self-contained exact
proof and import no external numerical or theorem claim.

## 9. Provenance, resource gates, and exact arithmetic

The producer records current SHA-256 hashes of itself, the verifier, this proof,
and `experiments/e241_box_genus_bound.py`. The verifier requires exact equality
with all four current hashes and checks a canonical hash of the certificate.
It imports neither producer.

All mathematical objects are Python integers, `Fraction`s, or exact
\(\mathbb Z[v]\) polynomials. SymPy is used only for exact polynomial
factorization over \(\mathbb Z\); no floating-point root or sign calculation
enters any claim. No benchmark coupling or fitted critical value is used.

The bounded inventories are:

- 18 target vertices and 33 target edges;
- nine eliminated stars and at most 512 retained-parity states per color;
- 120 boundary permutations, 16 even principal subsets per permutation, and
  600 four-leg residuals per signature;
- 4096 edge subsets and 16 retained-spin assignments in the `2x2x2` semantic
  control;
- 120 permutations and 1920 principal-Pfaffian checks for the planar control.

Both programs impose a process-CPU budget of 120 seconds and an RSS cap of
2 GiB. On Darwin they measure
`mach_task_basic_info.resident_size_max` through `task_info` flavor 20; on
other supported systems they convert `getrusage(RUSAGE_SELF).ru_maxrss` to
bytes. The producer raises before writing the JSON if any semantic,
provenance, symbolic, rational-control, or resource gate fails.

## 10. Exact scope boundary

What has been proved is:

1. the two explicit five-leg boundary Walsh signatures for the open `2x3x3`
   box;
2. the strict positive four-leg Pfaffian obstruction for all 120 orders and all
   \(0<v<1\);
3. the exact rational controls at \(v=1/3\) and \(v=1/2\);
4. acceptance of the planar alternating-cycle signature in exactly its ten
   dihedral orders.

What has not been proved is:

1. an obstruction under arbitrary local basis transformations;
2. an obstruction for arbitrary hidden auxiliary spins, gadgets, or
   unrestricted tensor-network rewrites;
3. a lower bound on the number of Pfaffians in an arbitrary global formula;
4. the same statement for other boundaries, other boxes, or all \(L\);
5. global Pfaffian solvability or nonsolvability of the three-dimensional Ising
   partition function;
6. a three-dimensional Ising solution, thermodynamic free energy, critical
   coupling, or critical exponents.

The finite boundary defect is therefore a rigorous obstruction in its stated
basis and ordering class, not a promotion to any of the excluded global or
auxiliary claims.
