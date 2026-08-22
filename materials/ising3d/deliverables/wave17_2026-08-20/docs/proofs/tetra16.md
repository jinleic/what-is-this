# One-line-cabled $16\times16$ auxiliary-$R$ RLLL certificate

[COMPUTATION] This note studies one precisely defined higher-auxiliary ansatz.  It does **not** solve the three-dimensional Ising model, and it does not classify a free local $\mathbb C^4$ leg.

## Corrected cabled problem

[DEFINITION] Let the binary Ising local operator be

$$
L(q)_{o,i}=\begin{cases}
q^{(|o|+|i|)/2},&|o|+|i|\equiv0\pmod2,\\
0,&|o|+|i|\equiv1\pmod2,
\end{cases}
$$

for $o,i\in\{0,1\}^3$.  The convention is
$L[\mathrm{output},\mathrm{input}]=T[\mathrm{input},\mathrm{output}]/2$.

[DEFINITION] Cable only auxiliary line $1$:

$$
\mathbb C^4_1=\mathbb C^2_{i_1}\otimes\mathbb C^2_{i_2},\qquad
L_4(q)_{145}:=L(q)_{i_1,4,5}\,L(q)_{i_2,4,5}. \tag{1}
$$

The displayed factor order in (1) is part of the definition.  It is material:
at $q=1/4$ the two embedded factors do not commute exactly.

[DEFINITION] In seven zero-based binary positions

$$
(i_1,i_2,2,3,4,5,6)=(0,1,2,3,4,5,6),
$$

let an unrestricted $R$ act on the first four bits and define

$$
R_{i_1 i_2 2 3}L_{4,145}L_{246}L_{356}
=L_{356}L_{246}L_{4,145}R_{i_1 i_2 2 3}. \tag{2}
$$

Thus $R$ is an arbitrary $16\times16$ matrix with 256 independent entries.
The global space is $16\cdot 2^3=128$ dimensional, so expanding (2) gives

$$
B(q)\operatorname{vec}(R)=0,\qquad B(q)\in\operatorname{Mat}_{16384\times256}(\mathbb Z[q]). \tag{3}
$$

[COMPUTATION] The previously quoted $65536\times256$ count is not the size of this canonical cabled formulation; (3) is the corrected count.

## Exact parity-sector reduction

[LEMMA] The 256 columns of $B(q)$ split into 128 parity-preserving and 128 parity-reversing $R$ coordinates, with disjoint row support.

Each binary $L$ preserves parity on its three legs.  Therefore the cabled product (1), and hence both ordered products in (2), preserve total seven-bit parity.  A matrix unit
$R_{o,i}$ contributes only to component rows with parity change
$|o|-|i|\pmod2$.  Ordering all parity-preserving entries first gives two independent $8192\times128$ row systems.

[COMPUTATION] At the exact control $q=1/4$, streamed construction finds all 16,384 component rows nonzero and exactly zero off-sector coefficients.  No $16384\times256$ dense matrix is materialized.  The selected $128\times128$ minor supports are connected bipartite graphs with respectively 1,955 and 1,980 nonzero entries.  The structured forward and reverse $128\times128$ products each have exact rank 16 at this control.

## Generic characteristic-zero rank

[THEOREM] Over the rational-function field $\mathbb Q(q)$,

$$
\operatorname{rank} B(q)=256,\qquad \ker B(q)=\{0\}. \tag{4}
$$

[THEOREM] Proof of (4).  For each parity sector, fix the 128 component rows stored in
`results/integrability/tetra16.json`.  They form a literal $128\times128$
minor $M_\epsilon(q)\in\mathbb Z[q]^{128\times128}$.  At $q=2$ modulo
$p=101$, exact finite-field elimination gives rank 128 for both sectors.  The
same independent result holds modulo $p=103$.  Hence each
$\det M_\epsilon(q)$ is not the zero polynomial: if it were zero in
$\mathbb Z[q]$, every integral specialization and every reduction modulo a
prime would vanish.  The two off-sector blocks vanish identically, so the
block diagonal product is a nonzero $256\times256$ minor of $B(q)$.  Since
$B$ has 256 columns, (4) follows.  The standalone verifier reconstructs the
operators and both modular eliminations without importing the producer.

[UNRESOLVED] Equation (4) is a generic-rank theorem, not an explicit
all-specialization theorem.  The zero locus of the two polynomial minors has
not been factored, so this work does **not** name the finite algebraic
exception set or assert rank 256 at every nonexceptional-looking $q$.

## Exact rational specializations

[COMPUTATION] For $q=a/d$, the calculation uses the integer local matrix
$\widehat L=d^3L(a/d)$.  Both ordered products and every coefficient row are
then scaled by $d^{12}$; a $128\times128$ sector determinant is scaled by
$d^{1536}$.  The following are complete prime factorizations of the nonzero
integer determinants $\det\widehat M_0$ and $\det\widehat M_1$ of the
stored sector minors:

| $q$ | parity-preserving $\det\widehat M_0$ | parity-reversing $\det\widehat M_1$ |
|---:|:---|:---|
| $1/4$ | $2^{1346}3^{290}5^{274}17^2\,1801$ | $2^{1332}3^{290}5^{272}\,11\,17^2\,463\,4373$ |
| $1/3$ | $2^{854}3^{668}5^2\,37$ | $2^{859}3^{665}5^2$ |
| $2/5$ | $2^{306}3^{290}5^{668}7^{274}\,13\,29^2\,367$ | $2^{292}3^{290}5^{664}7^{272}\,29^2\,31\,53\,18911\,20611$ |

Every listed factor is independently checked prime by the standalone test,
and the product is checked against the reconstructed fraction-free
$128\times128$ determinant.  Both sector determinants are nonzero at each
listed $q$, so the full rank is exactly 256 and the nullspace is zero:

$$
\begin{array}{c|c|c}
q&\operatorname{rank}B(q)&\dim\ker B(q)\\ \hline
1/4&256&0\\
1/3&256&0\\
2/5&256&0
\end{array}
$$

These are finite exact certificates, not an interpolation claim.

## Symbolic-minor wall

[UNRESOLVED] The same two $128\times128$ sparse-derived sector minors selected
at $q=1/4$ were rebuilt over $\mathbb Z[q]$ and sent to SymPy's
fraction-free domain determinant routine.  Each sector was given a hard
120-second wall.  Both timed out after 120.005033 and 120.005257 seconds,
for 240.158974 seconds total.  Process peak resident memory recorded by
Darwin `getrusage` was 157,990,912 bytes (83,968,000 bytes before the probe).
No symbolic determinant, factorization, complementary minor, or Bezout
identity was obtained.  This is the precise current obstruction to turning
(4) into an explicit all-specialization exception theorem.

## Reproduction and scope

[CHECK] Rebuild the artifact and its bounded symbolic search:

```bash
timeout 900 .venv/bin/python experiments/e83_tetra16.py
```

[CHECK] Run only the standalone verifier:

```bash
timeout 900 .venv/bin/python tests/test_tetra16.py
```

[SCOPE] The result concerns the fixed-order one-line cabling (1), three
otherwise identical binary Ising $L(q)$ factors, and an arbitrary ungraded
$16\times16$ $R$.  Cabling all three auxiliary lines is not attempted here.
An unconstrained $L$ on a genuine $\mathbb C^4$ leg, a free $\mathbb C^4$
local factor, nonidentical/spectral factors, IRF/dynamical relations, and
singular dimension-changing constructions are different problems and are
outside this certificate.
