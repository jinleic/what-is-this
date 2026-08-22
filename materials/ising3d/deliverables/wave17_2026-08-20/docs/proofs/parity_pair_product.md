# Ordering-free parity-pair product obstruction for the open `2 x 3` layer

**[THEOREM] Scope and status.** This note proves an ordering-free point
obstruction at `t = tanh(K*/2) = 1/3` and an all-but-finitely-many-`t`
upgrade, both for the one open `2 x 3` layer.  It uses neither eigenvalue
ordering nor inertia windows.  The producer is
`experiments/e119_parity_pair_product.py`, its exact artifact is
`results/spectral/parity_pair_product.json`, and the independent standalone
verifier is `tests/test_parity_pair_product.py`.

**[LEMMA] Physical-parity hypothesis.** Let
\(
P=\prod_{v=1}^{6}X_v
\), and suppose that its two physical sectors are the parity halves of one
six-mode subset-product spectrum:
\[
 \lambda_\epsilon=a\prod_{i=1}^{6}u_i^{\epsilon_i},
 \qquad \epsilon\in\{0,1\}^{6}.
\tag{1}
\]
Under one allowed assignment, `P=+1` is the even-cardinality half and `P=-1`
is the odd-cardinality half; under the other assignment the labels are
swapped.  The argument below permits arbitrary parameters in an algebraic
closure, including coincidences among the displayed values, because it works
with labelled eigenvalue slots.

## 1. Exact exponent-class counts

Put
\[
 H_p=\{\epsilon\in\{0,1\}^m:|\epsilon|\equiv p\pmod2\},
 \qquad p\in\{0,1\}.
\]
For \(f\in\{0,1,2\}^m\), write
\[
 r(f)=\#\{i:f_i=1\},\qquad s(f)=\#\{i:f_i=2\}.
\]

**[LEMMA] (1 — attainable parity-pair exponent classes).** For
\(\epsilon\in H_p\), \(\delta\in H_q\), and \(f=\epsilon+\delta\):

1. \(r(f)\equiv p+q\pmod2\).
2. If \(r(f)>0\), every \(f\in\{0,1,2\}^m\) satisfying that congruence is
   attained by some such \((\epsilon,\delta)\).
3. If `p=q` and the two labelled slots are required to be distinct, the
   attainable classes are exactly those with `r(f)` positive and even, hence
   `r(f) >= 2`.
4. If `p!=q`, the slots are automatically distinct and the attainable classes
   are exactly those with `r(f)` odd, hence `r(f) >= 1`.

*Proof.* A coordinate where \(f_i=0\) forces
\((\epsilon_i,\delta_i)=(0,0)\), a coordinate where \(f_i=2\) forces
\((1,1)\), and a coordinate where \(f_i=1\) permits either `(1,0)` or
`(0,1)`.  Thus
\[
 |\epsilon|+|\delta|=r(f)+2s(f),
\]
which proves (1).  For `r(f)>0`, choose a subset \(A\) of the `r(f)`
one-coordinates to receive `(1,0)` and give the complement `(0,1)`.  Then
\[
 |\epsilon|=s(f)+|A|,
 \qquad |\delta|=s(f)+r(f)-|A|.
\]
Among the subsets of a nonempty set there are subsets of either cardinality
parity, so choose \(A\) with
\(|A|\equiv p-s(f)\pmod2\).  The displayed congruence for `r(f)` then gives
\(|\delta|\equiv q\pmod2\).  Conversely, the forced-coordinate description
shows that no other `f` can occur.  When `p=q`, `r(f)` is even; `r(f)=0` means
\(\epsilon=\delta\), so distinct labelled slots force `r(f)>=2`.  When
`p!=q`, `r(f)` is odd and the two slots cannot coincide.  []

**[LEMMA] (2 — within-even and within-odd counts).** The products of two
distinct slots within `H_0` have exactly
\[
 A_{00}(m)=\sum_{\substack{2\le r\le m\\r\ \mathrm{even}}}
 {m\choose r}2^{m-r}
\tag{2}
\]
attainable exponent classes.  The same exact count holds within `H_1`:
\(A_{11}(m)=A_{00}(m)\).

*Proof.* By Lemma 1, choose the `r` coordinates equal to one and then assign
independently either zero or two to each of the other `m-r` coordinates.  This
gives \({m\choose r}2^{m-r}\) vectors for every permitted even `r>=2`.
The cases indexed by `r` are disjoint because `r(f)` is intrinsic to `f`.  The
proof does not depend on whether the common sector parity is zero or one,
because the orientation of the one-coordinates supplies the required parity
of \(\epsilon\).  []

**[LEMMA] (3 — cross-even/odd count).** The products of one `H_0` slot and one
`H_1` slot have exactly
\[
 A_{01}(m)=\sum_{\substack{1\le r\le m\\r\ \mathrm{odd}}}
 {m\choose r}2^{m-r}
\tag{3}
\]
attainable exponent classes.

*Proof.* Lemma 1 makes the permitted values exactly the odd positive values
of `r(f)`.  The same choice of the `r` one-coordinates and the independent
zero/two assignments on the complement yields the summand in (3), with no
overlap between different `r`.  []

**[COMPUTATION] Six-mode arithmetic.** At `m=6`, the three exact maxima are
\[
\begin{aligned}
 A_{00}(6)=A_{11}(6)
 &= {6\choose2}2^4+{6\choose4}2^2+{6\choose6}=240+60+1=301,\\
 A_{01}(6)
 &= {6\choose1}2^5+{6\choose3}2^3+{6\choose5}2
 =192+160+12=364.
\end{aligned}
\tag{4}
\]
The standalone verifier independently enumerates all labelled slots and all
exponent sums for every `m=2,...,6`, recovering (2)--(4) without importing the
producer.

## 2. Forced gcd thresholds

For a sector with eigenvalue slots \(\{\lambda_\epsilon\}\), define the
monic pair polynomials
\[
\begin{aligned}
 C_{00}(z)&=\prod_{\substack{\epsilon,\delta\in H_0\\\epsilon<\delta}}
             (z-\lambda_\epsilon\lambda_\delta),\\
 C_{11}(z)&=\prod_{\substack{\epsilon,\delta\in H_1\\\epsilon<\delta}}
             (z-\lambda_\epsilon\lambda_\delta),\\
 C_{01}(z)&=\prod_{\epsilon\in H_0,\ \delta\in H_1}
             (z-\lambda_\epsilon\lambda_\delta).
\end{aligned}
\tag{5}
\]
The first two use distinct occurrences even when their numerical eigenvalues
coincide.  Their degrees are \({32\choose2}=496\); the cross polynomial has
degree \(32\cdot32=1024\).

**[LEMMA] (4 — collision degree in characteristic zero).** If a monic
characteristic-zero polynomial of degree `N` has `d` distinct roots in an
algebraic closure, then
\[
 \deg\gcd(C,C')=N-d.
\tag{6}
\]

*Proof.* Factor \(C=\prod_\nu(z-\nu)^{m_\nu}\) over an algebraic closure.
Characteristic zero makes the local derivative multiplicity at \(\nu\)
equal to \(m_\nu-1\), so the gcd is
\(\prod_\nu(z-\nu)^{m_\nu-1}\).  Its degree is
\(\sum_\nu(m_\nu-1)=N-d\).  []

**[THEOREM] (1 — parity-pair necessary thresholds).** If (1) realizes the true
even and odd halves, then
\[
\begin{aligned}
 \deg\gcd(C_{00},C_{00}')&\ge496-301=195,\\
 \deg\gcd(C_{11},C_{11}')&\ge496-301=195,\\
 \deg\gcd(C_{01},C_{01}')&\ge1024-364=660.
\end{aligned}
\tag{7}
\]

*Proof.* Every root in (5) is `a^2` times a monomial indexed by one of the
classes counted in Lemmas 2--3.  Accidental relations among the `u_i`, zero
parameters, or repeated eigenvalues can merge roots but cannot create an
additional exponent class.  Thus the three distinct-root counts are at most
301, 301, and 364.  Lemma 4 gives (7).  []

**[LEMMA] (5 — both physical assignments).** Let `C_{++}`, `C_{--}`, and
`C_{+-}` denote the polynomials made from the physical `P=+1` sector, the
physical `P=-1` sector, and one slot from each, respectively.  The complete
necessary conditions are:

| physical-to-fermionic assignment | `C_{++}` | `C_{--}` | `C_{+-}` |
|---|---:|---:|---:|
| `P=+1 -> even`, `P=-1 -> odd` | EE `>=195` | OO `>=195` | EO `>=660` |
| `P=+1 -> odd`, `P=-1 -> even` | OO `>=195` | EE `>=195` | EO `>=660` |

*Proof.* Substitute the physical labels into (5) and then apply Theorem 1.
The two within-half numerical floors happen to agree, but their EE/OO roles
are retained in the table; the two assignments are logically distinct and
both must be rejected.  []

## 3. Integral finite-field certificate

**[LEMMA] (6 — monic integral sector models).** Let \(R_+\) and \(R_-\) be the
rational restrictions of the symmetrized transfer representative to the two
physical sectors.  If positive integers \(D_\pm\) clear their entry
denominators and \(M_\pm=D_\pm R_\pm\), then the within-sector pair
polynomials are the monic integral characteristic polynomials of
\(\bigwedge^2M_\pm\), and the cross-sector pair polynomial is the monic
integral characteristic polynomial of \(M_+\otimes M_-\).

*Proof.* Over an algebraic closure, triangularize `M_+` and `M_-`.  The
diagonals of their exterior squares are exactly the products of distinct
same-sector eigenvalue slots; the diagonal of their tensor product is exactly
the cross product of every `+` and `-` slot.  Exterior powers and tensor
products of integer matrices are integer matrices, so all three
characteristic polynomials are monic integral.  Scaling by the nonzero
\(D_\pm\) preserves equality relations among their respective roots.  []

**[LEMMA] (7 — modular gcd is a one-sided certificate).** For any monic
integral polynomial `C`, reduction modulo every prime `p` satisfies
\[
 \deg\gcd_{\mathbb Q}(C,C')\le
 \deg\gcd_{\mathbb F_p}(\bar C,\bar C').
\tag{8}
\]

*Proof.* The monic rational gcd is monic integral by Gauss's lemma and divides
both integral polynomials.  Its leading coefficient remains one after every
reduction, so it remains a common divisor with unchanged degree.  []

**[COMPUTATION] Exact rational input and scales.** For the open `2 x 3` bond
set at `t=1/3`, the common symmetrized model has
\(e^{2K}=5/3\), bond parity one, and commutes entrywise with the fixed-point
free complement permutation representing `P`.  In the unnormalised bases
\(b_a^\pm=e_a\pm e_{\bar a}\),
\[
 R_\pm=\tfrac12B_\pm^TRB_\pm,
 \qquad B_\pm^TB_\pm=2I.
\tag{9}
\]
The exact denominator LCMs are equal:
\[
 D_+=D_-=8{,}968{,}066{,}875=3^{15}5^4.
\tag{10}
\]
The script asserts integrality entry by entry after multiplication by these
LCMs.

**[COMPUTATION] Modular construction.** At each trial-division-certified prime
\(p\in\{1{,}000{,}003,2{,}000{,}003\}\), the script reconstructs a sector
characteristic polynomial from its power traces using Newton identities,
extends traces by Cayley--Hamilton, and checks selected recurrence traces with
independent binary matrix powering.  For a within sector it uses
\[
 q_m=\frac{(\operatorname{tr}M^m)^2-\operatorname{tr}M^{2m}}2,
 \qquad 1\le m\le496,
\tag{11}
\]
and for the cross polynomial it uses
\[
 q_m=\operatorname{tr}(M_+^m)\operatorname{tr}(M_-^m),
 \qquad 1\le m\le1024.
\tag{12}
\]
Newton reconstruction of the monic degree-496 or degree-1024 polynomial is
followed by an explicit Euclidean gcd with its formal derivative.  The
artifact records monicity, zero remainders on division by the returned gcd,
Cayley--Hamilton checks, trace spots, and coefficient digests.  The exact
integer guard at the larger prime is
\[
 1024(2{,}000{,}002)^2=4{,}096{,}008{,}192{,}004{,}096<2^{63},
\]
so every used `int64` modular dot/matrix product is bounded before reduction.

**[COMPUTATION] Degree table.** The two entries in every row are the gcd
degrees at `p=1000003` and `p=2000003`, respectively.  They are exact
finite-field values, not asserted characteristic-zero equalities.

| case | within `P=+1` | within `P=-1` | cross `P=+1 x P=-1` |
|---|---:|---:|---:|
| true synthetic six-mode Gaussian, actual E/O split | `195, 195` | `195, 195` | `660, 660` |
| open six-site chain physical sectors | `195, 195` | `195, 195` | `660, 660` |
| open `2 x 3` layer physical sectors | **`177, 177`** | **`15, 15`** | **`192, 192`** |

**[COMPUTATION] Controls.** Both controls meet all three thresholds in (7)
exactly at both certificate primes.  This validates that the three separate
slot conventions—EE, OO, and EO—are implemented rather than silently
collapsed into the full-spectrum pair polynomial.  The chain row is a
consistency control, not an independent all-coupling Gaussianity theorem.

## 4. Conclusion

**[THEOREM] (PP — ordering-free physical-parity no-go).** At
\(t=\tfrac13\), the open `2 x 3` Ising-layer transfer operator cannot have
its `P=+1` and `P=-1` spectra equal, in either assignment, to the even- and
odd-subset-product halves of one six-mode Gaussian spectrum (1).

*Proof.* In the assignment `P=+1 -> even`, `P=-1 -> odd`, the exact modular
values `177<195`, `15<195`, and `192<660` each contradict the corresponding
necessary condition in Lemma 5 after applying Lemma 7.  In the swapped
assignment, the same physical polynomials occupy the OO, EE, and EO slots,
respectively, so the values `15<195`, `177<195`, and `192<660` again
contradict Lemma 5.  Thus both possible identifications of physical spin flip
with fermionic parity are impossible.  []

## 5. Generic-coupling upgrade

**[THEOREM] (GP — headline all-but-finitely-many-couplings no-go).** There is
a nonzero polynomial \(\rho(t)\in\mathbb Q[t]\) such that, for every
\[
 t\in\mathbb C\setminus\bigl(\{0,i,-i\}\cup\{t:\rho(t)=0\}\bigr),
\tag{13}
\]
the open `2 x 3` layer cannot realize either physical-to-fermionic parity
assignment in (1).  The displayed exceptional set has at most
\[
 25{,}534{,}083
\tag{14}
\]
complex parameter values under the explicit crude bound proved below.  Thus
the physical-parity six-mode loophole is excluded for all but finitely many
couplings of this one layer.  The exceptional polynomial is certified to
exist by the algebraic argument below; it is not expanded or root-isolated in
this work.

**[LEMMA] (8 — polynomial physical `P=+1` family).** Write
\[
 q(t)=\frac{1+t^2}{2t},\qquad C(t)=8t^3(1+t^2)^4.
\tag{15}
\]
For the open `2 x 3` bonds, the exponents of \(q(t)\) in the rational
construction of `R` are exactly in \(\{-4,-3,\ldots,3\}\).  Consequently
\[
 M_+(t):=C(t)R_+(t)\in\mathbb Q[t]^{32\times32},
\tag{16}
\]
and every entry has `t`-degree at most 26.

*Proof.* For \(0\le m\le3\),
\[
 C(t)q(t)^m=2^{3-m}t^{3-m}(1+t^2)^{4+m},
\tag{17}
\]
whereas for \(m=-r\), \(1\le r\le4\),
\[
 C(t)q(t)^{-r}=2^{3+r}t^{3+r}(1+t^2)^{4-r}.
\tag{18}
\]
These are polynomials of degree at most 14 and 10, respectively.  Each
kernel product contributes \(t^{\operatorname{Ham}(a,k)+\operatorname{Ham}(b,k)}\)
of degree at most 12, so every entry of \(C(t)R(t)\) has degree at most
\(14+12=26\).  The fixed physical-parity projection is
\(R_+=\tfrac12B_+^TRB_+\), with \(B_+\) independent of `t`, hence preserves
polynomiality and the degree bound over \(\mathbb Q[t]\).  []

**[COMPUTATION] Symbolic construction check.** The producer constructs the
polynomial matrix (16) directly with integer coefficient lists, without
SymPy or floating-point arithmetic.  Its actual maximum entry degree is 26,
all sector coefficients are integral, and exact substitution
\(t=\tfrac13\) equals \(C(\tfrac13)R_+(\tfrac13)\) entry by entry.

Define
\[
 W(t)=\bigwedge^2 M_+(t),\qquad
 A(z,t)=\det(zI-W(t)),\qquad B(z,t)=\partial_z A(z,t).
\tag{19}
\]

**[LEMMA] (9 — bidegree bounds).** The polynomial \(A\) is monic of
`z`-degree 496 and belongs to \(\mathbb Q[t][z]\).  If
\[
 A(z,t)=z^{496}+a_1(t)z^{495}+\cdots+a_{496}(t),
\tag{20}
\]
then \(\deg_t a_j\le52j\).  The same slope bound holds for the coefficients
of \(B\), relative to its leading `z` coefficient.

*Proof.* An entry of \(\bigwedge^2M_+\) is a two-by-two minor of \(M_+\), so
has `t`-degree at most \(2\cdot26=52\).  The coefficient \(a_j\) is a sum of
`j`-by-`j` principal minors of this 496-dimensional matrix, giving
\(\deg_ta_j\le52j\).  Differentiating in `z` leaves all `t` degrees
unchanged.  []

**[LEMMA] (10 — specialization cannot lower the generic gcd).** Let `G` be
the monic gcd of \(A\) and \(B\) in \(\mathbb Q(t)[z]\), and put
\(g_{\rm gen}=\deg_zG\).  Then
\[
 \deg_z\gcd\bigl(A(z,t_0),B(z,t_0)\bigr)\ge g_{\rm gen}
\tag{21}
\]
for every \(t_0\in\mathbb C\).

*Proof.* Since \(A\) is monic and \(\mathbb Q[t]\) is a UFD, Gauss's lemma
puts both the monic factor \(G\) and \(A/G\) in \(\mathbb Q[t][z]\).  Division
of \(B\) by the monic \(G\) inside \(\mathbb Q[t][z]\) has zero remainder
because it has zero remainder over \(\mathbb Q(t)[z]\); hence
\(J=B/G\in\mathbb Q[t][z]\).  Since \(B=\partial_zA\) and \(A\) is monic of
`z`-degree 496, \(B\) has constant leading coefficient 496; division by the
monic \(G\) leaves \(J\) with that same nonzero leading coefficient (so `J`
is not asserted monic).  After every specialization, \(G(z,t_0)\) remains
monic of degree \(g_{\rm gen}\) and divides both specialized polynomials.
[]

**[COMPUTATION] Generic-degree input from the point certificate.** At
\(t=\tfrac13\), the point computation obtains modular gcd degree 177 for the
within-`P=+1` degree-496 polynomial.  Lemma 7 makes 177 an upper bound on the
characteristic-zero specialized gcd degree.  The integer model used there is
a nonzero scalar multiple of (16) at \(t=\tfrac13\), so it has the same pair
collision degree.  Lemma 10 therefore gives
\[
 g_{\rm gen}\le177.
\tag{22}
\]

**[LEMMA] (11 — generic equality away from one resultant).** Let
\[
 H=A/G,\qquad J=B/G,\qquad
 \rho(t)=\operatorname{Res}_z(H,J).
\tag{23}
\]
The factor \(H=A/G\) is monic, while \(J=B/G\) lies in
\(\mathbb Q[t][z]\) with constant leading `z` coefficient 496; moreover
\(H\) and \(J\) are coprime over \(\mathbb Q(t)\).  Then
\(\rho\in\mathbb Q[t]\) is nonzero, and at every \(t_0\) with
\(\rho(t_0)\ne0\),
\[
 \deg_z\gcd\bigl(A(z,t_0),B(z,t_0)\bigr)=g_{\rm gen}.
\tag{24}
\]

*Proof.* Lemma 10 places `H` and `J` in \(\mathbb Q[t][z]\), with `H` monic
and the stated nonzero leading coefficient for `J`.  They are coprime over
\(\mathbb Q(t)\) by the definition of `G`, so their resultant is nonzero.  A
nonzero specialized resultant makes \(H(z,t_0)\) and \(J(z,t_0)\) coprime.
Since \(G(z,t_0)\) divides both specialized original polynomials and is monic
of fixed degree, it is exactly their gcd.  []

**[LEMMA] (12 — effective resultant-size bound).** The resultant in (23)
satisfies
\[
 \deg_t\rho\le104\cdot496\cdot495=25{,}534{,}080.
\tag{25}
\]

*Proof.* A monic factor of a polynomial whose `j`th coefficient has
`t`-degree at most \(52j\) obeys the same coefficient slope: substitute
\(z=t^{52}y\), divide by the top power of `t`, and apply Gauss's lemma to the
resulting monic factor over \(\mathbb Q[t^{-1}][y]\).  Thus `H` has that
slope.  Monic long division of `B` by `G` then gives the same slope for `J`.
If \(h=\deg_zH\), then \(\deg_zJ=h-1\).  In the Sylvester determinant for
\(\operatorname{Res}(H,J)\), each term uses \(h-1\) coefficients of `H` and
\(h\) coefficients of `J`, so its `t`-degree is at most
\[
 (h-1)(52h)+h(52(h-1))=104h(h-1)
\le104\cdot496\cdot495.
\]
[]

**[THEOREM] (GP, proof).** For \(C(t)\ne0\) and \(\rho(t)\ne0\), Lemmas
10--11 and (22) give a within-`P=+1` gcd degree at most 177.  Lemma 4 then
gives at least \(496-177=319\) distinct `P=+1` pair products.  Under either
physical assignment, that sector would be either the even or the odd
six-mode half, for which Lemma 2 allows at most 301 distinct within-half pair
products.  Since \(319>301\), both assignments are impossible.  Equation
(25) bounds the root set of \(\rho\) by 25,534,080 values; adding the three
zeros of \(C(t)\) proves (14).  []

**[COMPUTATION] Deterministic rational-grid consistency check.** Each listed
point is rebuilt from the raw exact rational transfer matrix, split by the
physical `P`, integralized independently, and tested at
\(p=1{,}000{,}003\).  Every row has monic degree 496 and a returned gcd that
divides both the pair polynomial and its derivative.

| `t` | modular `P=+1` within-sector gcd degree |
|---:|---:|
| `1/4` | 177 |
| `1/3` | 177 |
| `2/5` | 177 |
| `1/2` | 177 |
| `2/3` | 177 |

**[COMPUTATION] Grid interpretation.** Every table row is itself an
independent pointwise no-go because `177<195`; the grid is a deterministic
consistency check on the generic argument, not a replacement for the
nonzero-resultant proof.

**[UNRESOLVED] Exact characteristic-zero collision counts and exceptional
parameters.** Agreement of modular degrees does not prove that a
characteristic-zero gcd degree is equal to `177`, `15`, or `192`; modular
reduction can add collisions.  Nor is \(\rho\) expanded or its finite root
set identified here.  The proved generic conclusion needs only the one-sided
upper bound in (22).

**[UNRESOLVED] Deliberately narrow exclusions.** The theorem does not exclude
a restriction of a Gaussian on a **larger** space, an arbitrary submultiset or
non-parity-invariant restriction of such a Gaussian, a fermionic grading not
identified with physical `P`, or a different layer.  It also does not decide
the finitely many exceptional `t` values in (13), including the denominator
roots, and is not an all-size or thermodynamic-limit statement.  It is an
all-but-finitely-many-parameter result only for the stated finite six-mode,
physical-parity-pair model of this one layer.
