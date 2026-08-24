# [THEOREM] Exact Galois spectra at a rational Ising specialization

## [LEMMA] Normalization and scope

[COMPUTATION] Fix

a) \(t=1/3\), and

b) \(q=(1+t^2)/(2t)=5/3\).

[LEMMA] For a layer graph \(G\) on \(n\) sites, the canonical rational representative reused from `e38_gaussianity_certificate.build_R` is

\[
R_G=P_t\,\operatorname{diag}\!\left(q^{(b_\sigma-\epsilon)/2}\right)P_t,
\qquad (P_t)_{\sigma\tau}=t^{d_H(\sigma,\tau)}.
\]

[COMPUTATION] The producer clears the actual entry-denominator lcm of \(R_G\), divides by the gcd of all resulting integer entries, and calls the primitive integer matrix \(B_G\).

[LEMMA] Thus \(B_G=c_G R_G\) for a nonzero rational \(c_G\). If a degree-six factor of \(\chi_{B_G}(x)\) is \(f(x)\), the displayed polynomial is

\[
g(y)=s^{-6}f(sy)
\]

for a positive integer \(s\). Its roots are the roots of \(f\) divided by \(s\), so \(f\) and \(g\) have the same splitting field and Galois group over \(\mathbb Q\).

[UNRESOLVED] The theorem below is about this rational normalization at the single exact point \(t=1/3\). It is not an all-parameter, all-size, thermodynamic-limit, or critical-coupling theorem.

## [COMPUTATION] Exact characteristic factorizations

[COMPUTATION] In the pattern column, \(d^e\) means one irreducible degree-\(d\) factor occurring with characteristic-polynomial exponent \(e\); an omitted exponent is one.

| tag | layer graph | physical role | dimension | exact irreducible degree/exponent pattern over \(\mathbb Q\) |
|---|---|---|---:|---|
| [COMPUTATION] | one site, no transverse edge | 1D control | 2 | \(1,1\) |
| [COMPUTATION][EXTERNAL] | periodic \(C_4\) chain | finite-width 2D integrable control | 16 | \(1^2,1^2,1^4,2,2,4\) |
| [COMPUTATION] | open \(2\times2\) square | 3D layer | 16 | \(1^2,1^2,1^4,2,2,4\) |
| [COMPUTATION][EXTERNAL] | open \(P_6\) chain | finite-width 2D/free-fermion control | 64 | \(1,1,6,6,15,15,20\) |
| [COMPUTATION] | open \(2\times3\) rectangle | 3D layer | 64 | \(1,1,6,6,6,6^2,9,9,14\) |

[COMPUTATION] Every coefficient of every factor, its exponent, its exact hash, and the exact reconstruction check is stored in `results/spectral/galois_spectrum.json`.

[THEOREM] The open \(2\times2\) square is not a new spectral control: the site map

\[
(0,1,2,3)\longmapsto(0,1,3,2)
\]

maps the cycle edges onto the square edges. The corresponding transfer matrices are permutation-similar, and their exact factor lists agree.

[THEOREM] All irreducible factors in this \(C_4\)/open-square case have degree at most four. Each factor therefore has a solvable Galois group, and the Galois group of the compositum of their splitting fields is solvable. This small control does not persist at width six.

## [THEOREM] The open \(2\times3\) layer has an \(S_6\) factor

[COMPUTATION] The lexicographically first irreducible sextic factor of the primitive open \(2\times3\) matrix is divided by the rational root scale \(s=57600\). The resulting explicit factor is

\[
\begin{aligned}
g_{3D}(x)={}&x^6-6613x^5+6886615x^4-2793624000x^3\\
&+498389850000x^2-34838100000000x\\
&+437400000000000.
\end{aligned}
\]

[COMPUTATION] Its exact discriminant is

\[
16083546436163715783952763789775981387854676312890625000000000000000000000000000000
\]

and its certified prime factorization is

\[
2^{30}3^{18}5^{38}11^2\cdot103\cdot120750823\cdot1343479977043\cdot52564834937707.
\]

[LEMMA] The odd exponents prove directly that the discriminant is not a square.

[COMPUTATION] The following are exact factorizations in the indicated finite fields. The displayed discriminant residues are nonzero, so all three primes are good.

| tag | prime | discriminant mod \(p\) | exact factorization mod \(p\) | cycle type |
|---|---:|---:|---|---|
| [COMPUTATION] | 13 | 7 | \(x^6+4x^5+8x^4+2x^3+3x^2+2x+11\), irreducible | \((6)\) |
| [COMPUTATION] | 53 | 6 | \((x+22)(x^5+43x^4+15x^3+52x^2+47x+3)\) | \((1,5)\) |
| [COMPUTATION] | 197 | 106 | \((x+42)(x+77)(x+139)(x+174)(x^2+47x+22)\) | \((1,1,1,1,2)\) |

[LEMMA] Irreducibility modulo 13 implies irreducibility of the monic integer polynomial over \(\mathbb Q\).

[LEMMA] At a prime not dividing the discriminant, the degrees of the irreducible modular factors are the cycle lengths of a Frobenius element in \(\operatorname{Gal}(g_{3D}/\mathbb Q)\).

[THEOREM] Let \(G\leq S_6\) be this Galois group. The irreducible reduction makes \(G\) transitive. The \((1,5)\) witness gives a 5-cycle. A transitive degree-six group containing a 5-cycle is primitive: a nontrivial block would have size two or three, while the 5-cycle and its fixed point force the block containing that point to have size one or six. The \((1,1,1,1,2)\) witness gives a transposition. Conjugates of that transposition are the edges of a \(G\)-invariant graph; primitivity makes the graph connected, and edge transpositions of a connected graph generate \(S_6\). Therefore

\[
\operatorname{Gal}(g_{3D}/\mathbb Q)=S_6.
\]

[THEOREM] Since \(S_6\) is not solvable, this sextic is not solvable by radicals over \(\mathbb Q\).

## [THEOREM] Decisive open-chain control counterexample

[COMPUTATION] For the open \(P_6\) layer at the same \(t=1/3\), the corresponding deterministic sextic and rational root scale \(s=4800\) are

\[
\begin{aligned}
g_{\mathrm{open}}(x)={}&x^6-782x^5+204493x^4-20321664x^3\\
&+615904560x^2-7138368000x+27993600000.
\end{aligned}
\]

[COMPUTATION] Its exact discriminant is

\[
171006597737218576740614323521392515712286720000000000000000
\]

with certified prime factorization

\[
2^{38}3^{29}5^{16}\cdot59406856868547588648571.
\]

[LEMMA] This discriminant is nonsquare.

[COMPUTATION] The independent good-prime cycle witnesses are:

| tag | prime | discriminant mod \(p\) | exact factorization mod \(p\) | cycle type |
|---|---:|---:|---|---|
| [COMPUTATION] | 7 | 5 | \(x^6+2x^5+2x^4+x^3+5x^2+6x+2\), irreducible | \((6)\) |
| [COMPUTATION] | 43 | 23 | \((x+9)(x^5+26x^4+9x^3+12x^2+15x+31)\) | \((1,5)\) |
| [COMPUTATION] | 587 | 473 | \((x+202)(x+444)(x+464)(x+533)(x^2+510x+338)\) | \((1,1,1,1,2)\) |

[THEOREM] The same transitive-plus-5-cycle-plus-transposition argument gives

\[
\operatorname{Gal}(g_{\mathrm{open}}/\mathbb Q)=S_6.
\]

[EXTERNAL] The \(P_6\) layer is the standard open transverse-chain transfer operator for a finite-width 2D Ising strip and is the conventional free-fermion/integrable control used elsewhere in this repository.

[THEOREM][EXTERNAL] Consequently, the proposed contrast “a nonsolvable characteristic factor occurs for the 3D layer but not for the 2D-integrable/open-chain control” is false at this exact specialization: both have an \(S_6\) factor.

[UNRESOLVED] This result must not be inverted into a nonintegrability or non-Gaussianity claim. A free-fermion reduction can still leave a finite one-particle secular polynomial with a nonsolvable Galois group. Non-solvability by radicals of a many-body characteristic factor is therefore not, by itself, a Gaussianity test.

## [LEMMA] Specialization direction

[LEMMA] The finite-field direction used above is one-way and local to each fixed polynomial: at a good prime, modular factor degrees exhibit elements of the characteristic-zero specialized Galois group. A modular factorization does not by itself identify the entire Galois group; the group-theoretic argument using all three witnessed cycle types does that here.

[LEMMA] For a separately constructed separable family \(F(t,x)\in\mathbb Q(t)[x]\), a good specialization group embeds into the generic Galois group. Thus, if the displayed sextic were first tracked as one degree-six factor in \(\mathbb Q(t)[x]\) and \(t=1/3\) were proved good for that family, the specialized \(S_6\) would force the generic degree-six group to be \(S_6\).

[UNRESOLVED] No such symbolic factor family over \(\mathbb Q(t)\) is constructed here. The exact \(t=1/3\) result therefore proves neither an all-\(t\) statement nor even a generic-\(t\) statement for a named family factor.

[UNRESOLVED] The benchmark value \(K_c=0.221654626\) was not used to choose \(t\), select a factor, search for a prime, fit a polynomial, or validate any conclusion.

## [COMPUTATION] Independent replay

[COMPUTATION] The standalone verifier `tests/test_galois_spectrum.py` imports none of `e203`, `e204`, or `e205`. It rebuilds the integer transfer matrices directly from tensor-entry integers, recomputes all characteristic factorizations, and checks the stored matrix and polynomial digests.

[COMPUTATION] Its finite-field audit does not trust factor degrees alone: it multiplies the stored modular factors with an independent coefficient-list implementation, checks squarefreeness by Euclid, and applies the finite-field irreducibility criterion using modular powers of \(x\).

[UNRESOLVED] Per the execution contract for this front, the verifier is written but not run by the producer agent; it is reserved for the lead's independent replay.
