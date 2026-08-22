# Support-size-stratified extensive charges beyond the wave-12 box frontier

## 1. Scope and finite class

Put, on \(\mathbb Z^d\),

\[
H=\sum_xX_x+\sum_{x\in\mathbb Z^d}\sum_{i=1}^d Z_xZ_{x+e_i},
\qquad D=\tfrac12[H,\,\cdot\,].
\]

The ordered-Pauli convention is \(X^aZ^b\).  As in
`proofs/extensive_charges.md`, \(\pi\) identifies every finite Pauli word
with all of its lattice translates.  Thus \(\pi Dq=0\) is exactly the
finite-density divergence condition for the formal translation sum of \(q\).
No infinite-volume operator sum is used as an element of the compact local
algebra.

For

\[
B_R=\{0,\ldots,R\}^d,\qquad
C_R=\{-1,\ldots,R+1\}^d,
\]

and an integer \(s\ge0\), define the **support-size class**

\[
\mathcal C_{d,R,s}=
 \operatorname{span}_{\mathbb Q}\{X^aZ^b:
 \operatorname{supp}(X^aZ^b)\subseteq B_R,
 \ |\operatorname{supp}(X^aZ^b)|\le s\}.
\tag{1}
\]

Here support size counts sites with a nonidentity local Pauli, so an \(XZ\)
site counts once.  This is **not** a point-group-invariant subclass: it
contains arbitrary asymmetric linear combinations of every word meeting the
bound.  It is only restricted by word support size.

Let

\[
S_{d,R,s}=\pi|_{\mathcal C_{d,R,s}},\qquad
M_{d,R,s}=\pi D|_{\mathcal C_{d,R,s}}.
\tag{2}
\]

All claims below are at the fixed ratio \(a=b=1\).

---

## 2. Exact structural reduction

**[LEMMA — support-size-stratified finite presentation].** Every output word
of \(D\) applied to a word in \(\mathcal C_{d,R,s}\) lies in \(C_R\) and
has support size at most \(s+1\).  Hence \(M_{d,R,s}\) is presented exactly
by the finite set of translation orbits of such output words that are actually
hit by the class columns.  Omitting all unhit orbits changes neither its kernel
nor its rank.

*Proof.*  In the repository convention,

\[
\tfrac12[X_u,X^aZ^b]=[b_u=1]X^{a+e_u}Z^b,
\tag{3}
\]

so a field term leaves the word support unchanged: the \(Z\)-bit at \(u\)
remains set.  For a bond \(uv\),

\[
\tfrac12[Z_uZ_v,X^aZ^b]
=-[a_u\mathbin{\mathrm{xor}}a_v=1]X^aZ^{b+e_u+e_v}.
\tag{4}
\]

The endpoint carrying the unique \(X\)-bit remains in the support because the
\(X\)-mask is unchanged.  At the other endpoint the new \(Z\)-bit either
adds one previously trivial site or removes a \(Z\)-only site.  Thus a bond
term has support size at most \(s+1\).  Every neighbor of \(B_R\) lies in
\(C_R\), proving the spatial claim.  Since each class column has finitely
many terms and the class has finitely many columns, only finitely many orbit
rows are hit; every other coordinate of \(\pi Dq\) is identically zero.
\(\square\)

**[LEMMA — exact class quotient].** The finite class quotient by trivial
translation shifts and the identity is

\[
\dim\frac{\ker M_{d,R,s}}
 {\ker S_{d,R,s}+\mathbb Q I}
 =\operatorname{rank}S_{d,R,s}-\operatorname{rank}M_{d,R,s}-1.
\tag{5}
\]

*Proof.*  The divergence criterion in `proofs/extensive_charges.md` gives
\(\ker S_{d,R,s}=\mathcal C_{d,R,s}\cap\ker\pi\): these are exactly the
class densities whose formal translation sum vanishes.  The identity is in
\(\mathcal C_{d,R,s}\), has \(DI=0\), and is not in \(\ker S_{d,R,s}\)
because its empty-word orbit coefficient is nonzero.  Rank-nullity applied to
(2) therefore gives (5).  This argument quotients by all finite divergence
densities that themselves lie in the stated class; it does not impose a
separate truncation on the divergence criterion. \(\square\)

**[LEMMA — analytic rank upper bound].** For \(R\ge1\) and \(s\ge2\),

\[
\operatorname{rank}_{\mathbb Q}M_{d,R,s}
\le\operatorname{rank}S_{d,R,s}-2.
\tag{6}
\]

*Proof.*  Besides \(I\), the density

\[
h=X_0+\sum_{i=1}^dZ_0Z_{e_i}
\tag{7}
\]

lies in the class because its individual words have sizes one and two.  Its
formal translate sum is \(H\), hence \(\pi Dh=0\); the scripts also
substitute (3)--(4) directly and obtain the empty orbit image.  The images of
\(I\) and \(h\) under \(S_{d,R,s}\) are independent: \(I\) occupies the
empty orbit while \(h\) has nonempty one-site and bond-orbit coordinates.
Thus \(\dim\ker M_{d,R,s}\ge2\), which is (6). \(\square\)

**[LEMMA — modular sandwich].** For an integer matrix \(A\),
\(\operatorname{rank}_{\mathbb F_p}A\le\operatorname{rank}_{\mathbb Q}A\).
Consequently, a modular rank reaching the right side of (6) certifies the
rational rank, rather than merely suggesting it.

*Proof.*  A nonzero minor modulo \(p\) is a nonzero integer minor, hence is
nonzero over \(\mathbb Q\). \(\square\)

---

## 3. Certified finite class census

**[COMPUTATION — exact arithmetic and certificates].**
`experiments/e125_extensive_3d_r2.py` enumerates the raw class columns in
blocks, builds the exact translation-orbit rows prescribed by the first lemma,
and stores each column image in a streamed CSR representation.  The smallest
new class is eliminated over \(\mathbb Q\) with `Fraction`; all larger classes
are eliminated at both primes

\[
p_1=2\,147\,483\,647,\qquad p_2=2\,147\,483\,629.
\]

A finite-field rank is used only as the lower inequality in the preceding
lemma.  The JSON includes each rank, deterministic pivot-trace SHA-256, and
full compact traces when at most 300,000 pivots occur.  The independent
standalone verifier recomputes all listed ranks at the third prime
\(1\,000\,000\,007\), whose primality it establishes by deterministic trial
division; it uses a separate bytes-based translation-orbit encoding and
recomputes all images from raw inputs.

| lattice / class | density columns | \(\operatorname{rank}S\) | active output orbit rows | exact or both modular ranks of \(M\) | upper bound from (6) | certified quotient |
|---|---:|---:|---:|---:|---:|---:|
| \(\mathbb Z^3,\ B_2,\ s\le2\) | 3,241 | 562 | 4,837 | 560 over \(\mathbb Q\) | 560 | 1 |
| \(\mathbb Z^3,\ B_2,\ s\le3\) | 82,216 | 29,749 | 343,783 | 29,747 at both \(p_1,p_2\) | 29,747 | 1 |
| **\(\mathbb Z^3,\ B_2,\ s\le4\)** | **1,503,766** | **822,334** | **11,093,578** | **822,332 at both \(p_1,p_2\)** | **822,332** | **1** |
| \(\mathbb Z^2,\ B_3,\ s\le4\) | 163,669 | 83,164 | 671,033 | 83,162 at both \(p_1,p_2\) | 83,162 | 1 |
| \(\mathbb Z^2,\ B_3,\ s\le5\) | 1,225,093 | 790,294 | 6,463,183 | 790,292 at both \(p_1,p_2\) | 790,292 | 1 |

The exact rational rank in the first row is 560.  Every other displayed
rational rank follows from the modular sandwich, not from an unsupported
identification of a modular and rational rank.

**[THEOREM — finite support-size class certificate].** For each row in the
table, every \(q\in\mathcal C_{d,R,s}\) satisfying \(\pi Dq=0\) has the
form

\[
q=d+\alpha I+\beta h,
\qquad d\in\mathcal C_{d,R,s}\cap\ker\pi,
\quad \alpha,\beta\in\mathbb Q.
\tag{8}
\]

In particular, within the clearly stated class \(\mathcal C_{3,2,4}\), there
is no non-Hamiltonian translation-covariant extensive density: modulo finite
translation shifts and the identity, the only class is \(h\).

*Proof.*  The table and the sandwich certify that the quotient in (5) has
dimension one.  The image of \(h\) is nonzero in that quotient because it has
nonempty Pauli-orbit coordinates whereas \(I\) does not.  It therefore spans
the one-dimensional quotient, which is equivalent to (8).  This is a
finite-box, finite-support-size theorem only. \(\square\)

**[COMPUTATION — regression anchors and defect traps].** The same producer
pipeline reproduces the frozen wave-12 full-box values

\[
(\mathbb Z^3,R=1):\quad
4^8=65\,536,\quad \operatorname{rank}S=64\,813,\quad
\operatorname{rank}M=64\,811,
\]

and

\[
(\mathbb Z^2,R=2):\quad
4^9=262\,144,\quad \operatorname{rank}S=254\,209,\quad
\operatorname{rank}M=254\,207.
\]

The standalone verifier independently repeats both values.  It also checks,
for every class column and every generated output term, the support-size and
current-box assertions in the first lemma; verifies \(\pi Dh=0\); checks
\(D(Z_0)=+X_0Z_0\) and that \(D(X_0)\) has six coefficient-\(-1\) bond
terms in three dimensions; and verifies that the one-dimensional oriented
chain current is conserved only in one dimension and has nonzero projected
commutator when transverse bonds are included.  These are deliberate traps
for a field/bond sign error, an omitted direction, a support off-by-one, or a
lost translation orbit.

---

## 4. Resource observations and the symmetry route

**[COMPUTATION — completed resource record].** The headline
\(\mathbb Z^3,B_2,s\le4\) enumeration took 147.605438 s; its two modular
eliminations took 13.799390 s and 12.929278 s.  Its recorded conservative
process-lifetime peak was 3,302.469 MiB.  The finished
\(\mathbb Z^2,B_3,s\le5\) enumeration took 84.568285 s and its two
eliminations took 10.550260 s and 10.388093 s, under the same conservative
peak.  The final process-lifetime maximum including the deliberately halted
probe was 5,133.094 MiB; that probe triggered when the observed value crossed
the 5 GiB RSS wall.

**[COMPUTATION — observed \(\mathbb Z^3,B_2,s\le5\) wall].** The next class
has exactly

\[
\sum_{k=0}^5\binom{27}{k}3^k=21\,121\,156
\tag{9}
\]

columns.  It was actually launched, not merely projected.  The enumeration
was stopped by the 5 GiB RSS wall after 171.824079 s and 1,654,784 columns,
at a measured 9,630.687 columns/s and 5,133.094 MiB process-lifetime RSS.
A linear extrapolation of the enumeration alone at that measured rate is
2,193.1 s; it is **not** a lower bound and excludes any subsequent rank work.
No rank or quotient is claimed for this class.

**[LEMMA — full-box anchored-orbit count].** For a full \(B_R\) box with
\(R\ge1\), the number of translation orbits met by its ordered-Pauli basis is

\[
N_{d,R}=1+\sum_{k=0}^d(-1)^k\binom dk
 4^{R^k(R+1)^{d-k}}.
\tag{10}
\]

*Proof.*  Apart from the identity, choose the unique translate whose support
has minimum coordinate zero in every direction.  Inclusion-exclusion over the
\(d\) coordinate-zero faces counts patterns touching every such face.  If
\(k\) chosen faces are avoided, support is restricted to
\(R^k(R+1)^{d-k}\) sites.  The initial one restores the identity orbit.
\(\square\)

**[COMPUTATION — symmetry reduction is quantitatively insufficient here].**
Equation (10) gives

\[
4^{27}=18\,014\,398\,509\,481\,984,
\qquad N_{3,2}=18\,014\,192\,401\,317\,889,
\tag{11}
\]

for the full \(3\times3\times3\) box, and

\[
4^{16}=4\,294\,967\,296,
\qquad N_{2,3}=4\,261\,675\,009
\tag{12}
\]

for the full \(4\times4\) square box.  These are input-size preflights; the
full systems were not launched.

**[LEMMA — cubic point-group lower bound].** The bond-preserving geometric
point group of the \(3\times3\times3\) cube has \(3!\,2^3=48\) signed
coordinate permutations.  Even a perfect reduction to its invariant
subspace leaves at least

\[
\left\lceil\frac{4^{27}}{48}\right\rceil
=375\,299\,968\,947\,542
\tag{13}
\]

orbit-summed density coordinates.

*Proof.*  The signed coordinate permutations preserve the cube and its
nearest-neighbor bonds.  Burnside's count for their action on the \(4^{27}\)
Pauli words is the average of nonnegative fixed-point counts and is therefore
at least \(4^{27}/48\); integrality gives the ceiling. \(\square\)

Thus symmetry reduction alone cannot make the full \(\mathbb Z^3,R=2\) box
finite work under the stated resource wall.  The support-size reduction in
(1), unlike symmetry restriction, retains all point-group behavior of the
bounded-size words.

---

## 5. Boundaries of the result

**[UNRESOLVED]** The headline theorem does **not** decide the full
\(\mathbb Z^3,R=2\) box: it leaves open any density that requires at least
one word with five or more nonidentity sites.  The observed \(s\le5\) wall
is a measured resource result, not evidence that no such density exists.

**[UNRESOLVED]** The full \(\mathbb Z^2,R=3\) box, all larger support-size
classes, all larger boxes, quasilocal densities, and non-translation-covariant
charges remain undecided.

**[UNRESOLVED]** The ranks certify only \(a=b=1\).  They do not identify or
exclude exceptional nonzero coupling ratios.

**[UNRESOLVED]** Nothing in this finite census proves non-integrability, a
solution, or an all-size theorem for the three-dimensional Ising model.
