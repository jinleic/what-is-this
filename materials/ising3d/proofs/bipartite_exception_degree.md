# A graph-uniform degree bound for full-spectrum exceptional couplings

Artifacts:

- `experiments/e245_bipartite_exception_degree.py`
- `results/spectral/bipartite_exception_degree.json`
- `tests/test_bipartite_exception_degree.py` (clean-room; no producer import)

Reproduce from the repository root:

```sh
.venv/bin/python experiments/e245_bipartite_exception_degree.py
.venv/bin/python tests/test_bipartite_exception_degree.py
```

## 1. Outcome and exact scope

Let \(G\) be a finite connected simple nonpath graph with \(n\) vertices and \(m\)
edges, and put \(N=2^n\). Wave 22 proved that the physical couplings at which the
complete positive-definite layer spectrum is one full nonzero \(n\)-mode subset-product
multiset form a finite set. The result here makes that qualitative statement effective:

\[
 \boxed{\#\mathcal E_G\le
 \left[2^n(2n+2m)\right]^{n+2}.}
 \tag{1}
\]

Because a simple graph has \(m\le n(n-1)/2\), this gives a bound uniform over
connected simple nonpath graphs of order \(n\):

\[
 \boxed{\#\mathcal E_G\le
 \left[2^n(n^2+n)\right]^{n+2}.}
 \tag{2}
\]

The count in the proof is deliberately larger: it counts the complete complex
coefficient-incidence image, including singular vertical boundary fibers. The physical
exceptional set in \(0<t<1\) is a subset. Thus no reality assumption is hidden in the
degree argument.

This is the first graph-uniform degree control on the finite exceptional schemes of
connected simple nonpaths. Paths are excluded: their full-spectrum incidence is cofinite.
It does **not** prove that the nonpath exceptional schemes are empty. The bound is
intentionally conservative and far too large for root enumeration.

## 2. Polynomial transfer representative

Write

\[
 q(t)=\frac{1+t^2}{2t},\qquad
 P_t(x,y)=t^{|x\mathbin\triangle y|}.
\]

For a spin configuration \(\sigma\), let \(a_\sigma=m-\operatorname{cut}_G(\sigma)\)
be the number of aligned edges. The uncentered symmetric representative

\[
 R_{\rm unc}(t)=P_t\,\operatorname{diag}(q(t)^{a_\sigma})\,P_t
 \tag{3}
\]

differs from the determinant-centered physical representative only by the scalar
\(q(t)^{\lceil m/2\rceil}\). Scalars do not change the subset-product property: they
are absorbed into the common slot \(c\).

Multiply (3) by \((2t)^m\). Every diagonal term becomes

\[
 (2t)^m q(t)^a=(2t)^{m-a}(1+t^2)^a,
 \qquad 0\le a\le m.
 \tag{4}
\]

It is a polynomial of degree \(m+a\le2m\). Each of the two Hamming kernels has
entry degree at most \(n\), so

\[
 \widetilde R_G(t):=(2t)^mR_{\rm unc}(t)\in\mathbb Z[t]^{N\times N},
 \qquad \deg_t(\widetilde R_{xy})\le\delta:=2n+2m.
 \tag{5}
\]

If

\[
 \det(xI-\widetilde R)=x^N+a_1(t)x^{N-1}+\cdots+a_N(t),
\]

then every \(k\)-th characteristic coefficient is a sum of \(k\)-minors, hence

\[
 \deg_t a_k\le k\delta.
 \tag{6}
\]

### Determinant boundary

The Hamming matrix is the \(n\)-fold Kronecker power of
\(\begin{psmallmatrix}1&t\\t&1\end{psmallmatrix}\), so

\[
 \det(P_t)^2=(1-t^2)^{nN}.
\]

Every edge is aligned in exactly \(N/2\) spin configurations, giving
\(\sum_\sigma a_\sigma=mN/2\). Therefore

\[
 \det\widetilde R_G(t)=
 (2t)^{mN/2}(1+t^2)^{mN/2}(1-t^2)^{nN}.
 \tag{7}
\]

Its roots are exactly \(0,\pm i,\pm1\), all outside the physical open interval.
Equation (7) also proves that the determinant is not the zero polynomial.

## 3. The affine coefficient-incidence family

A full \(n\)-mode subset-product polynomial has slots

\[
 \lambda_S=c\prod_{i\in S}u_i,
 \qquad S\subseteq[n].
\]

Let \(E_k(u)\) denote the \(k\)-th elementary symmetric polynomial in those \(N\)
Boolean monomials. Characteristic-polynomial equality is the affine system

\[
 F_k(t,c,u):=a_k(t)-(-1)^k c^kE_k(u)=0,
 \qquad k=1,\ldots,N.
 \tag{8}
\]

The mode term has total degree at most \(k(n+1)\). Since \(\delta\ge n+1\) for
the graphs under consideration, every equation in (8) has total degree at most

\[
 d=N\delta=2^n(2n+2m).
 \tag{9}
\]

The ambient affine space has exactly \(M=n+2\) coordinates
\((t,c,u_1,\ldots,u_n)\).

## 4. Why the full affine projection is finite

Let \(V\subset\mathbb A^{n+2}_{\overline{\mathbb Q}}\) be the common zero locus of
all equations (8), without saturation. On the open set \(a_N(t)\ne0\), the constant
coefficient equation forces \(c\ne0\) and every \(u_i\ne0\). Thus every point there
belongs to the nonzero-mode incidence family used by the Wave-22 theorem in
`e236_token_splitting_obstruction.py`.

That theorem combines the high-temperature token-graph obstruction with
[Chevalley's theorem](https://stacks.math.columbia.edu/tag/054K): the nonzero-mode
\(t\)-projection is constructible in the affine line, excludes a real open interval,
and is therefore finite. Points of the unsaturated variety with \(c=0\) or some
\(u_i=0\) can occur only where \(a_N(t)=0\), namely at the five roots in (7).
Consequently the full coordinate image \(\pi_t(V)\) is finite. No regular-sequence or
zero-dimensional-fiber assumption is needed; vertical fibers may have positive
dimension.

At a physical \(0<t<1\), \(\widetilde R_G(t)\) is a positive scalar times a
positive-definite matrix. In any complex solution of (8), both \(c\) and each
\(cu_i\) occur among its eigenvalue slots. They are therefore positive real, and so
are the \(u_i\). This recovers the physical-positive interpretation without imposing
positivity algebraically.

## 5. Affine degree gives the count

We use the affine degree bound introduced by Heintz:

> If an affine algebraic set in \(\mathbb A^M\) is a common zero locus of
> polynomials of degree at most \(d\), then its degree is at most \(d^M\).

A recent direct theorem with no dependence on the number of equations is Theorem 2.1
of Bierstone, Grigoriev, Milman, and Włodarczyk,
[*Effective resolution of singularities*, arXiv:2511.07639v1
(2025)](https://arxiv.org/abs/2511.07639v1): homogeneous equations of degree at
most \(d\) in projective \(M\)-space have common zero locus of degree at most
\(d^M\). Homogenizing (8), the affine components are exactly the projective
components meeting the affine chart, so their total degree satisfies the same
bound. Fact 19(g) of
[Helfgott, *Growth in finite simple groups of Lie type of bounded rank*,
arXiv:1005.1858](https://arxiv.org/abs/1005.1858) gives the affine formulation.
The original source is Joos Heintz, *Definability and fast quantifier
elimination in algebraically closed fields*, Theoretical Computer Science 24
(1983), 239--277,
[doi:10.1016/0304-3975(83)90002-6](https://doi.org/10.1016/0304-3975(83)90002-6).

Applying it to (8) gives

\[
 \deg V\le d^{n+2}.
 \tag{10}
\]

Because \(\pi_t(V)\) is finite, every irreducible component of \(V\) has one-point
image: the image closure of an irreducible set is irreducible, and an irreducible
finite subset of \(\mathbb A^1\) is one point. Hence

\[
 \#\pi_t(V)
 \le \#\{\text{irreducible components of }V\}
 \le \deg V
 \le d^{n+2}.
\]

This proves (1). Substituting \(m\le n(n-1)/2\) into \(\delta=2n+2m\) proves (2).
For the physical grid families, with \(L\ge2\) so neither graph is a path,

\[
 \delta_{2\times L}=10L-4,
 \qquad
 \delta_{3\times L}=16L-6.
\]

## 6. Exact finite controls

The producer independently enumerates aligned-edge exponents for a path, a cycle,
a claw, the open \(2\times3\) grid, and the open \(3\times3\) grid. It checks (4)--(10),
the determinant exponent sum, every Boolean-mode coefficient degree, and the
simple-graph uniformization.

For the disposed open \(2\times3\) layer, \((n,m,N)=(6,7,64)\),

\[
 \delta=26,
 \qquad d=1664,
 \qquad d^8=58779593138084256579321856.
\]

This number is only a scale control. No degree-971 elimination or branch computation
is rerun.

## 7. Limitations

Proved:

- a finite explicit complex-candidate bound for every fixed connected simple nonpath
  graph;
- the uniform connected-simple-nonpath \(n\)-vertex bound (2);
- the same bound for physical full-spectrum exceptional couplings on every bipartite
  rectangular layer with both side lengths at least two.

Not proved:

- emptiness beyond the already disposed open \(2\times3\) layer;
- a practical root count, a root list, or a graph-uniform root-free interval;
- parity-sector unions, larger auxiliary spectra, or any other notion weaker than one
  full nonzero \(n\)-mode spectrum;
- integrability, a thermodynamic solution, a critical point, or a moved critical
  endpoint.
