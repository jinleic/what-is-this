# Falsification of a proposed polynomial spectral obstruction

Artifacts: `experiments/e69_spectral_polynomial.py`, `results/spectral/polynomial_identity.json`, and `tests/test_spectral_polynomial.py`.

## Result

**[FALSIFIED]** The candidate polynomial relation developed in this front is
not a necessary identity for six positive Gaussian modes. It also arose from an
incorrect squared-trace normalization. The open 2x3 layer requires six modes,
so this candidate gives no obstruction, exceptional-root set, or continuum
Gaussianity theorem. The requested symbolic 2x3 application remains
**[UNRESOLVED]**.

All counterexamples below use exact integer algebra. No floating-point
calculation and no critical-coupling benchmark is used.

## Basis-free starting point

For a positive $n$-mode Gaussian subset-product spectrum

$$
\left\{c\prod_{i=1}^n r_i^{\epsilon_i}:\epsilon_i\in\{0,1\}\right\},
\qquad c,r_i>0,
$$

put $P_k=\operatorname{Tr}(V^k)$ and
$h=\det(V)^{1/2^n}>0$. Direct multiplication gives the correct basis-free
identity

$$
S_k=\frac{P_k}{h^k}
   =\prod_{i=1}^n\left(r_i^{k/2}+r_i^{-k/2}\right).
$$

This formula avoids all eigenvalue ordering, so crossings and degeneracies are
irrelevant. It is a sound starting point, but the subsequent attempted
elimination was not.

## Exact normalization error

Set $z=r+r^{-1}+2$. The square of one centered mode factor is

$$
(r^{k/2}+r^{-k/2})^2=r^k+r^{-k}+2=:g_k(z).
$$

The $g_k$ are Chebyshev polynomials in $z-2$. In particular,

$$
g_1(z)=z,\qquad g_2(z)=(z-2)^2,
$$

not $g_2(z)=z-2$ as used in the invalidated derivation. The discrepancy

$$
(z-2)^2-(z-2)=z^2-5z+6
$$

is a nonzero polynomial. Therefore the proposed norm equations were already
incorrect at order two.

The standalone test independently derives a correct one-mode identity:
if $Q_k=S_k^2$, then

$$
Q_2-(Q_1-2)^2=0,
$$

and verifies it exactly at $r=4,9,25$.

## Exact six-mode counterexample

There is a second, independent failure. Even a correctly normalized identity
obtained by eliminating exactly three elementary generators is not necessary
for six modes. Pairing two centered factors is invalid:

$$
(a^{k/2}+a^{-k/2})(b^{k/2}+b^{-k/2})
$$

contains four monomials and is not one factor
$R^{k/2}+R^{-k/2}$ for a fixed $R$ independent of $k$.
A six-mode norm has six independent elementary generators, not three.

Take the valid positive six-mode Gaussian control

$$
r_1=r_2=\cdots=r_6=1.
$$

Then $S_k=2^6=64$ for every $k$. The ratio variables used in the invalidated
attempt become

$$
(Q_2,Q_3,Q_4,Q_5)=(64,1,64,1).
$$

Direct substitution into the proposed polynomial $F$ gives

$$
F(64,1,64,1)=64(64-8)(64-1)(64+1)=14{,}676{,}480\ne0.
$$

A necessary six-mode identity must accept every positive six-mode Gaussian
control, so this one exact counterexample falsifies the proposed necessity.
Moreover, even if the variables are corrected to the genuine norm convention
$Q_k=S_k^2=4096$, exact substitution gives the nonzero residual
$-280{,}788{,}234{,}870{,}780$. Thus both the attempted normalization and the
three-to-six-mode generalization fail independently.

## Consequences for the open 2x3 layer

The layer has 64 eigenvalues and hence the relevant unprojected Gaussian ansatz
has six modes. Because the candidate relation is not necessary for six modes:

1. substituting exact layer traces into it cannot produce a valid Gaussian
   obstruction;
2. any degree-370 numerator, factorization, or Sturm count obtained from that
   substitution has no implication for Gaussianity;
3. no all-but-finitely-many theorem and no no-go on $[1/5,1/2]$ follows.

Those invalid headline fields have been removed from the JSON artifact. It is
marked `invalidated: true`, tagged `[FALSIFIED]`, and has status `UNRESOLVED`.
The artifact instead stores both exact counterexamples and the normalization
discrepancy.

## Correct frontier

A legitimate version of the proposed strategy must use

$$
g_k(z)=r^k+r^{-k}+2
$$

and retain all six elementary generators of a degree-six polynomial
$E(z)=\prod_{i=1}^6(z-z_i)$. It must eliminate those six generators from at
least seven exact norm equations, or derive another independently valid
six-mode invariant. That elimination was not completed here.

The honest quantified result is therefore negative: one explicit degree-four
candidate fails a positive six-mode Gaussian control in two exact conventions,
and its underlying norm construction is wrong already at $k=2$. This front
proves no new statement about the 2x3 layer or the 3D Ising model.
