# Exact Lee--Yang zeros for the 5x4x4 and 5x5x4 boxes

## [COMPUTATION] Result and scope

At the two frozen high-temperature protocol couplings $x=e^{-2K}=2/3$ and
$7/10$, `experiments/e88_lee_yang_zeros.py` computes the exact integer field
polynomials for the open all-direction boxes

\[
5\times4\times4 \quad\text{and}\quad 5\times5\times4,
\]

including every coefficient and exact isolating certificates for the first
Lee--Yang zero.  The producer records the full coefficients, CRT provenance,
and certificate data in `results/lee_yang/zeros_l5.json`; the clean-room
verifier is `tests/test_lee_yang_zeros.py`.

These are finite-box results.  They neither construct an exact $5^3$ cube
polynomial nor imply an all-size theorem or a thermodynamic edge exponent.

| coupling $x$ | box | transfer orientation | $\theta_1$ certified center | degree of $P(z)$ | roots of $Q(t)$ in $(-1,1)$ |
|---|---:|---:|---:|---:|---:|
| $2/3$ | $5\times4\times4$ | $4\times4\times5$ | $0.287184747933263041947248784497375583444167444708625700042242\ldots$ | 80 | 40 |
| $2/3$ | $5\times5\times4$ | $4\times5\times5$ | $0.252966574835448745784246982233424911589498615726051779125443\ldots$ | 100 | 50 |
| $7/10$ | $5\times4\times4$ | $4\times4\times5$ | $0.384639014895686190816470540431531831860008122307113075774992\ldots$ | 80 | 40 |
| $7/10$ | $5\times5\times4$ | $4\times5\times5$ | $0.348821331875047834529437300180295081015490745394908718352315\ldots$ | 100 | 50 |

Every displayed center lies in a rational angle interval of radius $10^{-60}$.
The corresponding first-$t$ interval has exact width

\[
2^{-260}=
\frac{1}{1852673427797059126777135760139006525652319754650249024631321344126610074238976}.
\]

All reduced roots are certified simple in every four headline calculation.
The CRT pass counts are respectively $13,16$ at $x=2/3$ and $24,30$ at
$x=7/10$, in the row order above.

## [LEMMA] Exact field-polynomial and circle-restriction certificate

For $x=n/d$, let $b(\sigma)$ be the number of broken bonds and let
$k(\sigma)$ be the number of down spins.  The producer stores the exact
cleared field polynomial

\[
P(z)=\sum_{k=0}^{N} A_kz^k,
\qquad
A_k=\sum_{k(\sigma)=k}n^{b(\sigma)}d^{|B|-b(\sigma)}\in\mathbb Z_{\ge0}.
\]

Spin reversal gives $A_k=A_{N-k}$.  For even $N=2m$, define the exact integer
polynomial $Q\in\mathbb Z[t]$ by

\[
z^{-m}P(z)=Q\!\left(\frac{z+z^{-1}}2\right).
\]

(For odd $N$, the same construction first removes the exact factor $z+1$.)
If a root $t\in(-1,1)$ is isolated by the exact polynomial $Q$, then
$z=e^{\pm i\arccos t}$ are the associated field roots.  This is the
circle-restriction identity used in the artifact and standalone verifier.

## [COMPUTATION] Exact zero isolation and finite unit-circle checks

For every produced box the code constructs $Q$ over `ZZ`, uses exact Sturm
root counting on $(-1,1)$, and uses exact rational real-root isolation for
the largest $t$ root.  It then certifies the angle enclosure with alternating
Taylor bounds for cosine: if $[t_-,t_+]$ is the isolated interval and
$[\theta_-,\theta_+]$ is stored, the certificate proves

\[
\cos(\theta_+)<t_-<t_*<t_+<\cos(\theta_-).
\]

For the headline boxes, the number of real reduced roots equals the reduced
degree ($40$ or $50$), all lie in $(-1,1)$, and all are simple.  Consequently
the $80$ or $100$ roots reconstructed from the reduced roots account for the
entire field polynomial and lie on $|z|=1$.  This is an exact **finite-box**
unit-circle check, not an all-size theorem.

## [COMPUTATION] Independent controls

The artifact records all of the following controls as passing.

1. The new $3\times3\times3$ transfer polynomial agrees coefficient-by-
   coefficient with the independent brute-force field density of states at
   both couplings.
2. The new $4\times4\times4$ polynomial agrees coefficient-by-coefficient
   with the frozen `scaling.json` cube data at both couplings.
3. For every new box, $\sum_kA_k$ agrees exactly with a separate scalar
   zero-field modular transfer at $z=1$.  Small admissible boxes also agree
   with `box_broken_bond_poly` evaluated at the same rational fugacity.
4. The complex128 fixed-field recurrence is only a diagnostic control; on
   every exact shape it agrees with the certified first angle to less than
   $5.2\times10^{-12}$.  It is not used for a $5^3$ claim.
5. `tests/test_lee_yang_zeros.py` does not import the producer.  It rebuilds
   the $3^3$ field polynomial from brute force, reconstructs the circle
   polynomial, reruns exact Sturm counts for both headline boxes at $x=2/3$,
   checks the cosine enclosures, independently checks a zero-field
   transfer-matrix identity, re-derives the four-size divided-difference
   edge interval, and checks the resource arithmetic.  It prints `PASS`.

## [LEMMA] Four-size analytic-correction test on exact slab ladders

For a fixed transverse cross-section and transfer lengths $L=2,3,4,5$, test

\[
\Theta_p(L)=e+L^{-p}\left(c_0+\frac{c_1}{L}+\frac{c_2}{L^2}\right).
\]

Let $u_L=1/L$ and

\[
w_L=\left(\prod_{M\ne L}(u_L-u_M)\right)^{-1}.
\]

The four values lie on a quadratic correction polynomial exactly only if the
third divided difference vanishes.  Since that condition is linear in the
edge,

\[
e_*(p)=
\frac{\sum_L w_LL^p\theta_L}{\sum_Lw_LL^p}.
\]

All operations in this calculation use exact rational intervals for the four
angles.  With $e=e_*(p)$, interpolation of the first three sizes gives
intervals for $c_0,c_1,c_2$.  The certificate accepts a model only when

\[
c_0>0,\quad c_2>0,\quad
\operatorname{disc}C<0,\quad
\operatorname{disc}\!\left(
 p c_0L^2+(p+1)c_1L+(p+2)c_2
\right)<0.
\]

These conditions prove $C(u)>0$ and $\Theta_p'(L)<0$ for every real $L>0$.

## [COMPUTATION] Continued degeneracy on two four-point slab ladders

Both adversarial powers remain admissible on every exact slab sequence once
its edge is allowed to be determined by the fourth point.  The table gives
the certified edge centers; their radii are inherited from the exact
$10^{-60}$ angle intervals and are serialized exactly in the artifact.

| $x$ | ladder | $e_*(p=2)$ | $e_*(p=4)$ | verdict |
|---|---:|---:|---:|---|
| $2/3$ | $4\times4\times L$ | $0.200597192754104906436817718830\ldots$ | $0.242767129179769430177209875249\ldots$ | both admissible |
| $2/3$ | $4\times5\times L$ | $0.172012041243902093830500768009\ldots$ | $0.211321346108943105354649413802\ldots$ | both admissible |
| $7/10$ | $4\times4\times L$ | $0.295523908092866598226145962152\ldots$ | $0.338745398315739526823964551806\ldots$ | both admissible |
| $7/10$ | $4\times5\times L$ | $0.264381894900782026819777762580\ldots$ | $0.305254473651059781600984614704\ldots$ | both admissible |

For all eight admissible cases, the artifact contains exact negative upper
bounds for both discriminants.  Thus a fourth **anisotropic** exact datum
does not itself remove the two-power ambiguity in this analytic-correction
class.

## [FALSIFIED] Fixed $e=1/50$ slab extensions of the wave-8 models

At each of the four slab/coupling combinations and for both $p=2,4$, the
exact third divided-difference interval at the wave-8 fixed edge $e=1/50$
excludes zero.  Therefore the **fixed-edge correction class** used in wave 8
(with its three correction coefficients reinterpolated on the slab data)
cannot extend to the fourth point of any of these slab ladders.

This is **not** an exclusion of either $p=2$ or $p=4$ for isotropic cubes:
allowing the edge to move restores certified globally positive, decreasing
models on both slabs, and slab geometry does not obey the isotropic relation
$\sigma=3/p-1$.

## [COMPUTATION] Exact resource wall for the $5\times5\times5$ cube

The open isotropic cube has $N=125$, $|B|=300$, and the implemented
axis-parallel transfer has a $5\times5$ cross-section.  Its exact
field-graded state array has $2^{25}$ state rows.  Even after using the
spin-reversal palindrome
to retain only grades $0,\ldots,62$, one final grade buffer requires

\[
2^{25}\cdot63\cdot8
=16{,}911{,}433{,}728\ \text{bytes}
=15.75\ \text{GiB}.
\]

The full 126-grade buffer is $33{,}822{,}867{,}456$ bytes ($31.5$ GiB).
The current exact butterfly path needs up to three simultaneous buffers;
the truncated conservative peak is $50{,}734{,}301{,}184$ bytes
($47.25$ GiB), while the full-width conservative peak is
$101{,}468{,}602{,}368$ bytes ($94.5$ GiB).  The host has 96 GiB shared
memory and was shared by concurrent agents.

A measured one-prime exact pass on the same $5\times4$ cross-section path
took $10.94$ seconds, corresponding to the recorded nominal rate
$1.257530960\times10^9$ element-operations/second.  Linear extrapolation of
the same code path gives the following lower resource requirements for the
palindrome-truncated $5^3$ computation:

| coupling | exact coefficient-bound bits | 30-bit CRT requirement | estimated compute time |
|---|---:|---:|---:|
| $2/3$ | 601 | 21 stored primes suffice | $3.1597$ hours |
| $7/10$ | 1122 | at least 38 30-bit moduli | $5.7175$ hours |

The checked-in CRT table has only 32 primes, so it cannot itself reconstruct
the $x=7/10$ cube coefficients even if memory were available.  These times
exclude contention and any allocator/paging penalty.  No $2^{25}$ run was
attempted, no $5^3$ polynomial was written, and no non-exact fixed-field
surrogate is promoted to an input certificate.

## [UNRESOLVED] Frozen PRG verdict for isotropic cubes

The frozen protocol in `proofs/lee_yang_prg.md` requires raw isotropic sizes
$L=4,\ldots,15$.  Exact isotropic cube data remain available only at
$L=2,3,4$; the exact $L=5$ cube datum is blocked as above.  The protocol is
therefore not instantiated, shortened, refit, shifted, or retuned.

Accordingly, neither the $p=2$ ($\sigma=1/2$) nor the $p=4$
($\sigma=-1/4$) wave-8 adversarial model is selected or excluded for
**isotropic cubes**.  The exact box and slab results narrow the computational
frontier and falsify only the fixed-edge slab extensions; they do not resolve
the three-dimensional Lee--Yang edge exponent.

## [COMPUTATION] Reproducibility

Run:

```sh
timeout 7200 .venv/bin/python experiments/e88_lee_yang_zeros.py
timeout 1800 .venv/bin/python tests/test_lee_yang_zeros.py
```

The producer prints one `PASS`/`FAIL` line for every check and writes
`results/lee_yang/zeros_l5.json`.  No external numerical benchmark or
external source is used in selection, fitting, or any certificate.
