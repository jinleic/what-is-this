# Translation-covariant extensive charges of the transverse-field Ising generator

## 1. Scope and status

This note addresses the extensive class deliberately outside the compact-local
scope of `proofs/local_commutant.md`.  Put, on \(\mathbb Z^d\),

\[
 H_{a,b}=a\sum_xX_x+b\sum_{x\in\mathbb Z^d}\sum_{i=1}^d Z_xZ_{x+e_i},
 \qquad ab\ne0,
\]

and consider formal translation sums

\[
 Q(q)=\sum_{x\in\mathbb Z^d}\tau_x(q)
\]

with a finite-support density \(q\).  Neither bare sum is asserted to be an
operator in the compact local algebra; only the finite density relation below
is used.

**[THEOREM — finite-density divergence criterion].** For a finite density
\(q\), the formal derivation \([Q(q),H_{a,b}]\) vanishes if and only if its
finite local density is a lattice divergence:

\[
 [q,H_{a,b}]_{\rm local}
  =\sum_{i=1}^d\bigl(s_i-\tau_{e_i}s_i\bigr)                 \tag{1}
\]

for finite-support densities \(s_i\).  This is a statement about formal
translation sums of finite Pauli densities, not an all-size solution of the
3D Ising model.

**[LEMMA — proof of the criterion].** Let \(\mathcal D\) be the rational
vector space with basis all finite ordered-Pauli words on \(\mathbb Z^d\), and
let \(\pi\) identify all translates of a word.  In each translation orbit a
finite density is a finitely supported function \(f:\mathbb Z^d\to\mathbb Q\).
If \(f=\sum_i(g_i-\tau_{e_i}g_i)\), its coefficient sum is zero by telescoping.
Conversely, if \(\sum_xf(x)=0\), choose one root site and route the coefficient
at every other occupied site to that root along a coordinate path.  The finite
edge flow on these paths gives finite \(g_i\) with the stated divergence.  Do
this independently in every Pauli translation orbit.  Thus

\[
 \ker\pi=\sum_i(1-\tau_{e_i})\mathcal D.                   \tag{2}
\]

Since the local commutator is finite and translations commute with the
derivation, \(\pi([q,H]_{\rm local})=0\) is equivalent to (1).  The script
uses \(D=\tfrac12[H,\cdot]\), whereas the displayed convention has
\([\cdot,H]=-2D\); multiplication by \(-2\) changes neither condition.
\(\square\)

## 2. Exact finite presentation and the trivial quotient

For an explicit radius convention, the calculation uses

\[
 B_R=\{0,\ldots,R\}^d,\qquad
 C_R=\{-1,\ldots,R+1\}^d.                                  \tag{3}
\]

Let \(V_R\) be the rational span of every \(X^aZ^b\) supported in \(B_R\).
It has \(4^{(R+1)^d}\) ordered-Pauli density coefficients.  The local
commutator of a word in \(V_R\) has support in \(C_R\).  In the repository
convention,

\[
 \frac12[X_s,X^aZ^b]=[b_s=1]X^{a+e_s}Z^b,
\quad
 \frac12[Z_uZ_v,X^aZ^b]=-[a_u\mathbin{\mathrm{xor}}a_v=1]X^aZ^{b+e_u+e_v}.
                                                                    \tag{4}
\]

**[LEMMA — finite current matrix].** For every active output Pauli translation
orbit, the experiment joins all output translates by a rooted coordinatewise
tree.  It introduces exactly one coefficient of a current \(s_i\) per
oriented tree edge.  The incidence columns are precisely
\(P-\tau_{e_i}P\), lie in \(C_R\), and span the zero-total-coefficient
subspace of that orbit.  Therefore the sparse finite equation system

\[
 Dq-\sum_i(1-\tau_{e_i})s_i=0                              \tag{5}
\]

is exact for the divergence question, even though it uses a tree basis rather
than every redundant current word in \(C_R\).

Let \(S_R=\pi|_{V_R}\) and \(M_R=\pi D|_{V_R}\).  Equation (2) says that
\(\ker S_R\) is exactly the pure-divergence density shift space.  The identity
density is separately quotiented, exactly as required in the task.  Hence the
reported quotient is

\[
 \dim\frac{\ker M_R}{\ker S_R+\mathbb QI}
 =\operatorname{rank}S_R-\operatorname{rank}M_R-1.          \tag{6}
\]

This proves both that the quotient is neither too small (it removes every
finite zero translation sum) nor too large (a surviving class has a nonzero
formal translation sum).

**[COMPUTATION — audit data].** The JSON records for every completed case the
number of density unknowns, tree-current unknowns by direction, active
ordered-Pauli equation rows, ranks, nullities, timings, and the conservative
process-lifetime maximum RSS at measurement end.  The corresponding
full-system rank is `current_unknown_count + rank(M_R)`, because every
rooted-tree incidence map is injective and its image is disjoint from the
chosen orbit-sum quotient.

The standalone verifier `tests/test_extensive_charges.py` rebuilds the
per-orbit rooted trees independently: it confirms that every incidence column
is a pure difference \(t-\tau_{e_i}t\), that each tree has incidence rank
\(|V|-1\) so that its column space is exactly the zero-orbit-sum subspace,
that a nonzero orbit sum is unsolvable, and that the per-direction current
counts and equation counts match the stored presentation for \(d=1,2\) at
\(R=1\), with both transverse directions exercised in \(d=2\).

## 3. Mandatory one-dimensional positive control

**[EXTERNAL]** The transverse-field Ising chain is free-fermionic and has a
local charge hierarchy.  At the first nontrivial box one familiar member is
the oriented energy-current density

\[
 j=X_0Z_0Z_1-X_1Z_0Z_1
   =-i\bigl(Y_0Z_1-Z_0Y_1\bigr),                             \tag{7}
\]

where the last equality only translates the repository's \(X^aZ^b\)
convention to Hermitian \(Y\)'s.

**[COMPUTATION]** The same parameterized ordered-Pauli code gives
\(M_Rj=0\) at \(d=1,R=1\), not a separate one-dimensional implementation.
Exact `Fraction` elimination gives:

| \(d=1\) radius \(R\) | \(4^{R+1}\) density coefficients | \(\operatorname{rank}S_R\) | \(\operatorname{rank}M_R\) | quotient dimension |
|---:|---:|---:|---:|---:|
| 0 | 4 | 4 | 3 | 0 |
| 1 | 16 | 13 | 10 | 2 |
| 2 | 64 | 49 | 44 | 4 |
| 3 | 256 | 193 | 186 | 6 |
| 4 | 1,024 | 769 | 760 | 8 |
| 5 | 4,096 | 3,073 | 3,062 | 10 |
| 6 | 16,384 | 12,289 | 12,276 | 12 |
| 7 | 65,536 | 49,153 | 49,138 | 14 |
| 8 | 262,144 | 196,609 | 196,592 | 16 |
| 9 | 1,048,576 | 786,433 | 786,414 | 18 |

Thus the control detects both the Hamiltonian density and an independent
current already at \(R=1\), and the quotient continues to grow through the
largest completed exact chain box, \(R=9\).  At \(R=10\), the 4,194,304-column
ordered-Pauli preparation reached its explicit 120-second wall after
1,960,704 columns at 2,258.953 MiB process peak RSS, so every chain radius
\(R\ge10\) remains unresolved here.  The control would fail if the code
omitted the telescoping quotient, used the opposite Pauli sign incorrectly,
or omitted a bond direction.

## 4. Two- and three-dimensional finite results

The universal density

\[
 h=X_0+\sum_{i=1}^dZ_0Z_{e_i}                                \tag{8}
\]

has translation sum \(H_{1,1}\).  It and the identity are independent
classes in \(V_R/\ker S_R\) for \(R\ge1\), and direct exact substitution
shows \(M_Rh=0\).  Consequently

\[
 \operatorname{rank}_{\mathbb Q}M_R\le \operatorname{rank}S_R-2. \tag{9}
\]

**[COMPUTATION — exact-Q small boxes].**

| lattice | \(R\) | density coefficients | \(\operatorname{rank}S_R\) | \(\operatorname{rank}_{\mathbb Q}M_R\) | quotient | active equations / tree currents |
|---|---:|---:|---:|---:|---:|---:|
| \(\mathbb Z^2\) | 0 | 4 | 4 | 3 | 0 | 10 / 0 |
| \(\mathbb Z^2\) | 1 | 256 | 229 | 227 | 1 | 1,282 / 163 |
| \(\mathbb Z^3\) | 0 | 4 | 4 | 3 | 0 | 14 / 0 |

These three ranks were obtained by exact rational Gaussian elimination.
At \(d=2,R=1\), the one surviving class is (8).  In particular the embedded
one-dimensional density (7) has nonzero projected commutator once the
transverse bonds are included.

**[COMPUTATION — finite-field lower bounds plus an independent exact upper
bound].** The larger cases were reduced modulo each of
\(2\,147\,483\,647\) and \(2\,147\,483\,629\).  A modular rank is used
only as \(\operatorname{rank}_{\mathbb F_p}M_R\le
\operatorname{rank}_{\mathbb Q}M_R\).  In each row below it reaches the
separate upper bound (9); this sandwich, not an unsupported modular equality,
certifies the displayed rational rank.

| lattice | \(R\) | density coefficients | \(\operatorname{rank}S_R\) | ranks at both primes | upper bound (9) | certified quotient | active equations / tree currents |
|---|---:|---:|---:|---:|---:|---:|---:|
| \(\mathbb Z^2\) | 2 | 262,144 | 254,209 | 254,207 | 254,207 | 1 | 1,835,646 / 66,431 |
| \(\mathbb Z^3\) | 1 | 65,536 | 64,813 | 64,811 | 64,811 | 1 | 852,106 / 10,481 |

The `pivot_trace_column_row_base85` payloads, their little-endian `uint32`
pair encoding, counts, and SHA-256 hashes in
`results/integrability/extensive_charges.json` are the two-prime finite-field
certificates.  Per-case timings and the conservative process-lifetime RSS
observations are recorded in that artifact; the explicit 7 GiB wall was never
reached by either completed transverse case.

**[THEOREM — finite-box generic-ratio corollary].** For each completed
nontrivial box in the preceding table, the quotient dimension is one for all
but finitely many nonzero complex ratios \(g=b/a\).

*Proof.* In the ordered-Pauli basis, \(M_R(g)=M_{R,X}+gM_{R,ZZ}\) has entries
in \(\mathbb Q[g]\).  The certified rank at \(g=1\) exhibits a maximal-size
minor that is nonzero at \(1\), so that minor is a nonzero polynomial.  Away
from its finite zero set, the rank is at least the rank displayed above.
For every nonzero \(g\), the identity and
\(X_0+g\sum_iZ_0Z_{e_i}\) provide the two-dimensional kernel used in (9), so
the rank cannot be larger.  Therefore the generic quotient is exactly the
Hamiltonian class. \(\square\)

The theorem is finite-radius only.  It does not name the exceptional ratios,
and it does not turn a finite computation into an all-range assertion.

## 5. Interpretation and boundaries

**[COMPUTATION]** At \(d=3,R=1\), the exact finite quotient has dimension
one.  That one class is the already present Hamiltonian density (8), modulo
identity and finite divergences.  No additional non-Hamiltonian three-
dimensional density was found.  Therefore the requested independent periodic
finite-box verification for a new 3D density is not triggered.

**[UNRESOLVED]** No all-range conclusion follows.  In particular, this work
does not rule out a density supported in a larger box, a quasilocal density,
a non-translation-covariant charge, or another mechanism of integrability.

**[UNRESOLVED]** The calculation directly evaluates only \(a=b=1\).  The
generic-ratio corollary leaves a finite, unidentified exceptional set of
nonzero ratios for each finite box; it is not a uniform \(ab\ne0\) theorem.

**[UNRESOLVED — input-size preflight, not an observed wall]** In
\(\mathbb Z^2\), \(B_3\) has 16 sites and
\(4^{16}=4,294,967,296\) density coefficients before currents.  This case was
not launched.  The count is an input-size projection under the stated run
budget, not a demonstrated 7 GiB / 240-second lower bound; every
\(\mathbb Z^2\) radius \(R\ge3\) remains unresolved.

**[UNRESOLVED — input-size preflight, not an observed wall]** In the stated
convention \(B_2\subset\mathbb Z^3\) has 27 sites and
\(4^{27}=18,014,398,509,481,984\) density coefficients before currents.  This
case was not launched.  The count is an input-size projection, not an observed
resource wall; every \(\mathbb Z^3\) radius \(R\ge2\) remains unresolved here.
