# Streamed symmetry-reduced circle evaluations and the remaining $5^3$ zero-count wall

## [COMPUTATION] Result and scope

This artifact implements an exact, streamed alternative to the field-graded
$2^{25}$ layer array.  At one rational circle coordinate

\[
t=\frac{z+z^{-1}}2\in(-1,1),
\]

it evaluates the reduced Lee--Yang polynomial $Q(t)$ directly, rather than
materializing its 126 field coefficients.  The implementation uses exact
residues in

\[
\mathbb F_p[z]/(z^2-2tz+1),
\]

CRT reconstruction with a proved signed bound, spatial cross-section
symmetries, spin-flip conjugation, and bounded pair strips in every
interlayer butterfly.

It is a valid compressed exact evaluator, but it does **not** produce an
exact $5\times5\times5$ zero or a Sturm root count.  No $5^3$ value is
claimed in this note or in `results/lee_yang/zeros_53.json`.  The fixed-value
method removes the old tens-of-GiB *storage* obstacle, but reconstructing the
degree-62 $Q$ by enough exact evaluations is projected to be slower than the
already-recorded graded transfer and still requires subsequent exact Sturm
work.  Thus this is the requested validated fallback, not a resolution of the
isotropic Lee--Yang edge obstruction.

The requested historical path `proofs/lee_yang_zeros_l5.md` is absent in the
checkout.  The applicable prior artifact is
`proofs/lee_yang_zeros.md` with `results/lee_yang/zeros_l5.json`.

## [LEMMA] Spatial-plus-spin quotient at fixed rational $t$

Let $S=\{0,1\}^{n_s}$ be one open cross-section, let $G$ be its spatial
automorphism group, and let $v_j(s;z)$ be the partial partition function of
$j$ layers whose final state is $s$.  The in-layer weight is
$w_s(z)=x^{b(s)}z^{d(s)}$, with the rational bond fugacity cleared modulo a
prime, and the interlayer kernel is

\[
K(s,r)=d^{n_s-h(s,r)}n^{h(s,r)}.
\]

Here $x=n/d$, $b(s)$ is the number of broken in-layer bonds, $d(s)$ is the
number of down spins, and $h$ is Hamming distance.  For every $g\in G$,

\[
v_j(gs;z)=v_j(s;z),\qquad
v_j(\bar s;z)=z^{jn_s}v_j(s;z^{-1}). \tag{1}
\]

**Proof.**  At $j=1$, both statements follow from invariance of $b$ and $d$
under $G$, and from $b(\bar s)=b(s)$ and $d(\bar s)=n_s-d(s)$.  If they hold
at layer $j$, then $K(gs,gr)=K(s,r)$ and
$K(\bar s,\bar r)=K(s,r)$.  Applying the transfer sum and multiplying by
$w_s(z)$ preserves both identities, giving the $j+1$ case.  Induction proves
(1). $\square$

At rational $t$, write each stored value as $a+bz$ in the displayed
quadratic algebra.  For a state in the complementary spatial orbit, (1) is
recovered by the exact semilinear map

\[
z^{jn_s}\overline{(a+bz)},\qquad \bar z=z^{-1}=2t-z.
\]

Thus one pair of base-field residues represents an entire $G\times\{1,
\mathrm{flip}\}$ orbit.  This is not an invalid same-grade spin quotient:
the factor $z^{jn_s}$ and algebra conjugation are explicitly retained.
Because the butterfly coefficients are in the base field, its two algebra
components may be streamed separately through the same exact factorized
butterflies.

For even $N=2m$ the final scalar is

\[
Q(t)=z^{-m}P(z),
\]

and for odd $N=2m+1$ it is

\[
Q(t)=z^{-m}\frac{P(z)}{z+1}.
\]

The final $z$ coefficient vanishes exactly by spin reversal.  This invariant
is checked after every modular pass.

## [LEMMA] Signed CRT bound for a circle evaluation

For $t=a/b\in[-1,1]$ in lowest terms and $B$ open-box bonds, define

\[
U=2^N\max(n,d)^B.
\]

Then every cleared field coefficient is nonnegative and their sum is at most
$U$.  Dividing an odd palindromic polynomial by $z+1$ produces alternating
partial sums, each with absolute value at most $U$.  Writing the resulting
palindrome in the Chebyshev basis and using $|T_j(t)|\leq1$ gives

\[
\left|b^mQ(a/b)\right|\leq (N+1)U b^m. \tag{2}
\]

Consequently, CRT reconstruction of the signed integer on the left is unique
once the product of distinct admissible moduli exceeds twice the right-hand
side.  The code uses (2), centers the CRT residue, and rejects a reconstruction
outside this exact bound.

## [COMPUTATION] Exact independent reproduction certificates

At $x=2/3$, the new code exactly reproduces $Q(t)$ from independent complete
field-polynomial data at two rational points for each shape:

| shape | $t$ | exact $Q(t)$ | CRT primes | butterfly strips | wall (s) |
|---|---:|---:|---:|---:|---:|
| $4\times4\times4$ | $0$ | $-40319076205794220175421543617960442199318722643939536904540936745854$ | 10 | 960 | 0.2602220829867292 |
| $4\times4\times4$ | $1/2$ | $-7789132650221794414263456295247747100430473906411092662630309836953089$ | 12 | 1152 | 0.2911415419948753 |
| $5\times4\times4$ | $0$ | $185927737984183974446086767901155045972663846449223331664945221900263268910379755208386$ | 13 | 1560 | 6.889724416003446 |
| $5\times4\times4$ | $1/2$ | $-126762842231443377543258365089632596233781543365159619135870010223242755187846025646572129$ | 14 | 1680 | 7.481306665999 |

The $4^3$ reference is rebuilt through the stable
`rational_field_polynomial_transfer` API.  The $5\times4\times4$ reference
is the independently produced, axis-permuted $4\times4\times5$ coefficient
vector stored in `zeros_l5.json`; the equality tested is an exact Fraction
identity, not a floating comparison.  The standalone test separately
brute-forces the $4\times4$ spatial-plus-spin orbit count, independently
computes $Q(t)$ from the full coefficient vector, checks a forward versus
reverse butterfly-strip order modulo an independent CRT prime, exercises the
odd-$N$ $z+1$ removal identity against an independent $3^3$ polynomial, and
repeats that odd-size check through a temporary disk-backed vector path.

## [COMPUTATION] Exact Burnside counts and refined $5^3$ layout

For the $5\times5$ square cross-section, the spatial group is $D_8$.  Burnside
counting gives

\[
\frac{2^{25}+4\,2^{15}+2^{13}+2\,2^7}{8}
=4{,}211{,}744
\]

spatial orbits.  Every $D_8$ operation fixes the central site, so a
spatial-complement operation has no fixed binary state.  Including spin flip
therefore gives exactly

\[
\frac{33{,}693{,}952}{16}=2{,}105{,}872
\]

stored algebra orbits, versus $2^{25}=33{,}554{,}432$ raw states.

With two raw uint64 component buffers, a uint32 state-to-orbit map plus a
boolean flip map, ten long-lived uint64 orbit-component vectors, and the
persistently retained representative/count/statistics/flip-index metadata,
the explicit implementation subtotal is

\[
536{,}870{,}912+167{,}772{,}160+168{,}469{,}760+113{,}438{,}048
=986{,}550{,}880\ \text{bytes}.
\]

The ten component vectors are the two layer factors, two derived factor
vectors, current pair, transformed pair, and next pair (or the final flipped
pair).  The metadata term includes the uint32 representative and flip-index
arrays, two uint64 multiplicity arrays, and uint8 broken/down arrays.  This
still excludes NumPy expression temporaries, allocator overhead, transient
symmetry-map construction arrays, and optional disk-page effects; it is not a
measured $5^3$ peak.

The validated $4^3$, $t=0$ run used 10 CRT primes, 96 exact butterfly strips
per prime, and 0.2602220829867292 seconds total.  Counting scalar butterfly
slots exactly gives

\[
\frac{2(5-1)25\,2^{25}}{2(4-1)16\,2^{16}}
=\frac{3200}{3}
\]

times the $4^3$ fixed-evaluation work.  The resulting **linear extrapolation
only** is 582.8974658902735 seconds for $Q(0)$ (21 CRT primes) and
638.4115102607757 seconds for $Q(1/2)$ (23 primes).  It excludes $5\times5$
symmetry-map construction, memory-bandwidth/cache effects, disk traffic, and
contention; it is not promoted to a guaranteed runtime.

A full coefficient reconstruction through this evaluator would require 63
distinct values for the degree-62 reduced polynomial.  On the concrete grid
$j/32$, $-31\leq j\leq31$, the worst signed bound has 918 bits and needs 31
of the existing 30-bit CRT primes.  The same linear arithmetic projects
54209.46432779543 seconds (about 15 hours) for the 63 values, before map
construction, interpolation, or Sturm.  This quantified route is therefore
not a competitive exact zero-count procedure relative to the recorded
3.2-hour graded $x=2/3$ coefficient-wall estimate.

## [UNRESOLVED] Correction-amplitude bridge alternative

No monotone majorant relating the slab first zeros of $5\times4\times4$ or
$5\times5\times4$ to the isotropic $5^3$ first zero is established here.
The exact fourth slab points in `proofs/lee_yang_zeros.md` still allow both
$p=2$ and $p=4$ analytic correction models when the edge is allowed to move.
Independently, `proofs/lee_yang_prg.md` gives exact, globally decreasing
isotropic continuations through the $L=2,3,4$ cube data with the common
positive edge $1/50$ and exponents separated by $3/4$.

The stored finite-volume data alone therefore do not supply a correction
amplitude bound.  A bridge would require a new proved inequality for complex
field zeros under transverse layer addition, bond addition, or box inclusion;
real-field correlation monotonicity is not such an inequality.  No such lemma
is asserted.

## [COMPUTATION] Reproducibility

```sh
timeout 7200 .venv/bin/python experiments/e98_lee_yang_53.py
timeout 1800 .venv/bin/python tests/test_lee_yang_53.py
```

The producer writes `results/lee_yang/zeros_53.json`.  The standalone test
prints `PASS` after all exact checks.  Neither command uses a three-dimensional
critical benchmark.
