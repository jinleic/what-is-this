# Simon–Lieb finite-box lower bound for the simple-cubic Ising model

## Result

For the nearest-neighbour ferromagnetic Ising model on \(\mathbb Z^3\), with
\(K=\beta J\) and \(v=\tanh K\), the computation in
`experiments/e12_simon_lieb.py` proves

\[
\boxed{K_c \;>\; 0.20517823158099.}
\]

More precisely, for the free box \(B=\{0,1,2,3\}^2\times\{0,\ldots,31\}\),
with the origin placed at `(1,1,15)`, the exact rational endpoint

\[
v_- = \frac{20234669789378}{10^{14}}
\]

satisfies the Simon–Lieb criterion \(\kappa_B(v_-)<1\). Directed interval
arithmetic gives

\[
\operatorname{atanh}(v_-)
\in [
0.2051782315809953090983856061486869898893488540691077937478209740588918613438798145708,
\]
\[
0.2051782315809953090983856061486869898893488540691077937478209740588918613438798169428
].
\]

The displayed 14-place bound is the lower endpoint rounded **down**. The
published numerical estimate \(0.221654626\) was used only as a falsification
check, never as an input or fit target.

## 1. The Simon–Lieb inequality and its hypotheses

Let \(G=(V,E)\) be a locally finite graph and consider the zero-field Ising
Hamiltonian with pair couplings \(J_{xy}=J_{yx}\ge 0\). Assume the sums that
occur below converge (automatic for nearest-neighbour \(\mathbb Z^3\)). Let
\(\langle\cdot\rangle_\beta^+\) denote the infinite-volume plus Gibbs state.
For a finite \(S\subset V\), let
\(\langle\cdot\rangle_{S,\beta,0}\) denote the finite Ising model containing
only the spins and couplings internal to \(S\), at zero field—equivalently,
**free boundary conditions**.

For \(o\in S\) and \(z\notin S\), the modified Simon–Lieb inequality is

\[
\langle\sigma_o\sigma_z\rangle_\beta^+
\le
\sum_{\substack{x\in S\\y\notin S}}
\tanh(\beta J_{xy})
\langle\sigma_o\sigma_x\rangle_{S,\beta,0}
\langle\sigma_y\sigma_z\rangle_\beta^+ .
\tag{SL}
\]

Only pairs with \(J_{xy}>0\) contribute. This statement requires
ferromagnetic couplings, zero field for the restricted correlation, and a
finite set containing the first endpoint; no periodic boundary condition is
inserted on \(S\).

The historical sources are:

1. B. Simon, “Correlation inequalities and the decay of correlations in
   ferromagnets,” *Communications in Mathematical Physics* **77** (1980),
   111–126, [doi:10.1007/BF01982711](https://doi.org/10.1007/BF01982711).
2. E. H. Lieb, “A refinement of Simon’s correlation inequality,”
   *Communications in Mathematical Physics* **77** (1980), 127–135,
   [doi:10.1007/BF01982712](https://doi.org/10.1007/BF01982712). Lieb’s
   refinement is what permits the first correlation to be evaluated in the
   finite inside system with free boundary conditions.

For the exact crossing-edge form (SL), including its coefficient, see Lemma
2.7 (“Modified Simon’s inequality”) of H. Duminil-Copin and V. Tassion, “A
new proof of the sharpness of the phase transition for Bernoulli percolation
and the Ising model,” *Communications in Mathematical Physics* **343**
(2016), 725–745,
[doi:10.1007/s00220-015-2480-z](https://doi.org/10.1007/s00220-015-2480-z),
[arXiv:1502.03050](https://arxiv.org/abs/1502.03050).

### Why the coefficient is \(\tanh K\), not \(K\)

The backbone proof of (SL) isolates the first edge \(\{x,y\}\) that leaves
\(S\). Its factor is the zero-field correlation of the isolated two-spin
system:

\[
\langle\sigma_x\sigma_y\rangle_{\{x,y\},\beta,0}
=
\frac{e^{\beta J_{xy}}-e^{-\beta J_{xy}}}
     {e^{\beta J_{xy}}+e^{-\beta J_{xy}}}
=	anh(\beta J_{xy}).
\]

Thus \(\tanh K\) is the verified coefficient used here. Some presentations
replace it by \(\beta J_{xy}\). That remains valid but is weaker, because
\(\tanh t\le t\) for \(t\ge0\); it would produce a smaller lower bound. No
such weakening is needed in this computation.

## 2. From (SL) to a finite-box criterion

For unit nearest-neighbour coupling on \(\mathbb Z^3\), define

\[
q_B(x)=\#\{y\notin B:y\sim x\}=6-\deg_B(x)
\]

and

\[
\kappa_B(K)
=
\tanh K\sum_{x\in B}q_B(x)
\langle\sigma_o\sigma_x\rangle_{B,K,0}.
\tag{1}
\]

Translation invariance turns (SL) into a convolution inequality whose total
kernel mass is exactly \(\kappa_B(K)\). If \(\kappa_B(K)<1\), iterating the
inequality each time a correlation leaves a translate of \(B\) gives
exponential decay with distance. Equivalently, summing the iterated
inequality gives finite susceptibility; the standard estimate in the
modified Simon argument is of the form

\[
\chi(K)\le \frac{|B|}{1-\kappa_B(K)}<\infty.
\]

Therefore the plus magnetization vanishes. Moreover, (1) is continuous in
\(K\) for finite \(B\), so a strict inequality persists on an interval
\([K,K+\epsilon]\). Consequently \(K+\epsilon\le K_c\), and hence

\[
\kappa_B(K)<1\quad\Longrightarrow\quad K<K_c.
\tag{2}
\]

For the singleton \(B=\{o\}\), the free correlation is
\(\langle\sigma_o^2\rangle=1\) and \(q_B(o)=6\), so (1) reduces exactly to

\[
6\tanh K<1,
\qquad
K<\operatorname{atanh}(1/6)=0.168236118310606\ldots .
\]

This is an exact normalization check on both the coordination number and the
edge factor.

## 3. Exact finite-volume correlations

Expand every internal bond as

\[
e^{K\sigma_i\sigma_j}=\cosh K\,(1+v\sigma_i\sigma_j).
\]

Let

\[
P_B(v)=\sum_{F\subset E(B):\,\partial F=\varnothing}v^{|F|},
\qquad
P_{o,x}(v)=\sum_{F\subset E(B):\,\partial F=\{o,x\}}v^{|F|}.
\]

The spin sum forces even degree at every vertex in the partition function and
odd degree precisely at \(o,x\) in the correlation numerator. Thus, exactly,

\[
\langle\sigma_o\sigma_x\rangle_{B,K,0}
=\frac{P_{o,x}(v)}{P_B(v)}.
\]

For \(x=o\), the source set is empty and \(P_{o,o}=P_B\), as required by
\(\sigma_o^2=1\). Define the non-negative integer polynomial

\[
Q_B(v)=\sum_{x\in B}q_B(x)P_{o,x}(v).
\]

Then the implemented criterion is the exact rational function

\[
\boxed{\kappa_B(v)=\frac{vQ_B(v)}{P_B(v)}}.
\tag{3}
\]

Ferromagnetic correlation monotonicity and monotonicity of \(v=\tanh K\)
show that \(\kappa_B\) is increasing. It starts at zero and tends to the
number of crossing boundary bonds as \(v\uparrow1\), so the reported root is
unique.

## 4. Two independent exact computations

### 4.1 Exhaustive enumeration for the \(3^3\) cube

For the 27-site cube, the code fixes \(\sigma_o=+1\) and enumerates all
\(2^{26}=67,108,864\) remaining assignments in chunks. Global spin flip
accounts exactly for the omitted half. For each assignment it records the
number \(b\) of satisfied internal bonds and
\(\sigma_o\sum_xq_B(x)\sigma_x\). If \(H_b\) and \(G_b\) are the resulting
integer histograms and \(E=54\), then

\[
P_B(v)=2^{-26}\sum_b H_b(1+v)^b(1-v)^{E-b},
\]

\[
Q_B(v)=2^{-26}\sum_b G_b(1+v)^b(1-v)^{E-b}.
\]

Every coefficient division by \(2^{26}\) is checked to be exact. The output
passes

\[
P_B(0)=1,
\quad P_B(1)=2^{E-N+1}=2^{28},
\quad Q_B(1)=54\,2^{28},
\]

and all coefficients of both polynomials are non-negative integers. The full
coefficient lists and SHA-256 fingerprints are in
`results/bounds/simon_lieb_bounds.json`.

As an independent check, these polynomials evaluated at \(v=1/5\) agree
exactly with the parity-transfer computation described next.

### 4.2 Parity transfer for side four and long prisms

A layer has \(m=n_xn_y\) sites and a state is an \(m\)-bit mask giving the
parities entering that layer. For each intralayer edge with endpoint mask
\(e\), selecting or omitting that edge applies

\[
f(s)\longleftarrow f(s)+v f(s\mathbin\oplus e).
\tag{4}
\]

If the prescribed source mask in the layer is \(t\), the outgoing vertical
edge mask \(r\) is uniquely fixed by parity, and the layer transition is

\[
f_{\mathrm{next}}(r)=v^{|r|}f(r\mathbin\oplus t).
\tag{5}
\]

Three vector streams compute \(P\) and the weighted two-source sum \(Q\) in
one pass: a source-free stream, a stream with the fixed source \(o\), and a
marked stream summing over the possible second source \(x\) with weight
\(q_B(x)\). This costs \(O(n_z(m+|E_{\rm layer}|)2^m)\) operations and only
\(O(2^m)\) states. The largest calculation here has \(m=16\), hence 65,536
states, even for the \(4\times4\times32\) prism.

At a rational \(v=p/q\), (4) is replaced by the exact integer update

\[
f(s)\longleftarrow qf(s)+pf(s\mathbin\oplus e),
\]

and the vertical weight is \(p^{|r|}q^{m-|r|}\). Therefore the returned
integers are exactly \(q^E P(p/q)\) and \(q^E Q(p/q)\), with no floating-point
operation. A two-stream variant directly evaluates an integer with the sign
of \(pQ(p/q)-qP(p/q)\). This exact sign certifies the boxes through \(4^3\).

## 5. Numerical-error audit for the long prisms

For \(4\times4\times8\), \(4\times4\times16\), and
\(4\times4\times32\), the same recurrence was evaluated independently with
`mpmath` at 50 and 80 decimal digits. All recurrence terms are non-negative,
so there is no cancellation in \(P\) or \(Q\).

For a box with \(E\) internal edges and \(N\) sites, the implementation
bounds the arithmetic depth by

\[
n_{\rm op}=32(E+N+1).
\]

It deliberately uses the inflated unit roundoff
\(u=10^{-(\mathrm{dps}-5)}\), five decimal orders larger than the working
unit, forms \(\gamma=n_{\rm op}u/(1-n_{\rm op}u)\), and assigns the relative
allowance \(8\gamma\) to the final \(vQ/P\). This covers input rounding,
power construction, both positive transfer streams, and the final
multiplication and division. The 80-digit value is enlarged by this allowance
before testing it against one.

The audit values at the guarded rational endpoints are:

| box | exact decimal \(v_-\) | \(\kappa\), 50 dps | \(|\kappa_{50}-\kappa_{80}|\) | roundoff-enlarged margin \(1-\kappa\) |
|---|---:|---:|---:|---:|
| \(4\times4\times8\) | 0.20184064196866 | 0.999999999999692092108692011250831020878818768733109… | \(1.5343\times10^{-50}\) | \(3.079078913079887\times10^{-13}\) |
| \(4\times4\times16\) | 0.20232906509514 | 0.999999999999583651305851731955535132945536558499338… | \(5.5704\times10^{-51}\) | \(4.163486941482680\times10^{-13}\) |
| \(4\times4\times32\) | 0.20234669789378 | 0.999999999999536548929913873294324679678507417138722… | \(3.2245\times10^{-50}\) | \(4.634510700861267\times10^{-13}\) |

For the final box, the conservative relative allowances are
\(4.54912\times10^{-40}\) at 50 dps and \(4.54912\times10^{-70}\) at 80
dps. They are more than 27 orders of magnitude smaller than the certified
margin. Finally, `mpmath.iv` at 80 dps encloses
\(\operatorname{atanh}(v_-)\), and the reported decimal is rounded down from
that interval’s lower endpoint.

## 6. Certified sequence

The numerical root column locates \(\kappa_B=1\); only the guarded rational
endpoint and the last column are used as rigorous bounds.

| free box \(B\) | sites | internal bonds | crossing bonds | numerical root \(K_B\) | certified \(K<K_c\) endpoint |
|---|---:|---:|---:|---:|---:|
| \(1\times1\times1\) | 1 | 0 | 6 | 0.16823611831060645 | 0.16823611831058 |
| \(2\times2\times2\) | 8 | 12 | 24 | 0.18346372155215190 | 0.18346372155214 |
| \(3\times3\times3\) | 27 | 54 | 54 | 0.19740388284055022 | 0.19740388284053 |
| \(4\times4\times4\) | 64 | 144 | 96 | 0.20193444723608650 | 0.20193444723606 |
| \(4\times4\times8\) | 128 | 304 | 160 | 0.20465062730735270 | 0.20465062730734 |
| \(4\times4\times16\) | 256 | 624 | 288 | 0.20515984606810345 | 0.20515984606808 |
| \(4\times4\times32\) | 512 | 1264 | 544 | 0.20517823158101176 | **0.20517823158099** |

Every certified endpoint in this chosen sequence is non-decreasing. This is a
checked fact for the displayed boxes, not a claim that arbitrary nested
finite sets must produce monotone Simon–Lieb roots. The long non-cubic boxes
were included because they preserve a 16-site transfer cross-section while
moving the two end faces farther from the origin.

## 7. Comparison with the other lower bound and with the benchmark

The independently computed SAW/connective-constant result in
`results/bounds/kc_bounds.json` is

\[
K_c\ge 0.2074277114992039908436804465100887760027.
\]

It is stronger than the best Simon–Lieb number here by

\[
0.002249479918214003504317588522098958492279\ldots .
\]

Both rigorous lower bounds remain below the published simple-cubic numerical
estimate \(0.221654626\). The Simon–Lieb bound is lower by
\(0.01647639441901\ldots\), so it passes the required wrong-side
falsification test. The benchmark played no role in the construction of the
boxes, the transfer recurrence, or the rounded endpoints.

## Reproduction

From the repository root:

```text
.venv/bin/python tests/test_simon_lieb.py
.venv/bin/python experiments/e12_simon_lieb.py
```

The test prints `ALL Simon--Lieb checks: PASS`. The experiment writes the
machine-readable certificate to
`results/bounds/simon_lieb_bounds.json`, including provenance, exact
polynomials for \(3^3\), integer-sign fingerprints, both multiprecision
values, forward-error allowances, interval endpoints, and every PASS check.
