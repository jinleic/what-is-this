# Bipartite log-spectrum moment obstruction

## 1. Scope and notation

Let `n >= 1`, `N = 2^n`, and let a positive spectrum mean a multiset of
`N` strictly positive real numbers, counted with algebraic multiplicity.  Its
empirical log random variable is

\[
L=\log \lambda_I,\qquad \Pr(I=i)=N^{-1},
\]

and its determinant-centered version is

\[
X=L-\mathbb E L
  =\log\lambda_I-\frac1N\log\!\left(\prod_{i=1}^N\lambda_i\right).
\]

The centered empirical cumulants are

\[
\kappa_r=\left.\frac{d^r}{dz^r}
  \log\mathbb E(e^{zX})\right|_{z=0}.
\]

The symbol `p_j` below **does not mean a cumulant or a spectral power trace**.
It means the power sum of squared mode log-energies,

\[
p_j=\sum_{i=1}^n y_i^j,\qquad y_i=E_i^2\ge 0.
\]

[UNRESOLVED] The finite theorem proved here applies deliberately to the open
`2x2` and `2x3` bipartite layer graphs on one explicit rational interval of the
physical isotropic parameter; the `P4` and `P6` chains are controls, not
applications of the theorem.  It makes no all-coupling, all-size, `2x4`,
trace-reciprocity, isotropic-localization, varying-conjugator, or pair-polynomial
sequence claim.

## 2. Determinant centering gives a Rademacher sum

[LEMMA] Suppose the positive multiset is a full `n`-mode subset-product
multiset,

\[
\lambda_S=a\prod_{i\in S}u_i,
\qquad S\subseteq\{1,\ldots,n\}.
\]

Then `a>0` and every `u_i>0`: the empty-subset slot is `a`, and the corresponding
singleton slot is `a u_i`.  Multiplying all `N` slots gives

\[
\prod_S\lambda_S
 =a^N\prod_{i=1}^n u_i^{N/2}.
\]

Thus the positive geometric mean is

\[
g=\left(\prod_S\lambda_S\right)^{1/N}
  =a\prod_{i=1}^n u_i^{1/2}.
\]

For a uniformly random subset `S`, put

\[
\epsilon_i=2\mathbf 1_{i\in S}-1,
\qquad E_i=\frac12\log u_i.
\]

The `epsilon_i` are independent uniform signs and

\[
X_S=\log(\lambda_S/g)=\sum_{i=1}^n\epsilon_i E_i. \tag{2.1}
\]

This proves every normalization and uses no ordering, distinctness, nonzero-mode,
or simplicity assumption.  A factor `u_i=1` merely gives `E_i=0`.  Equation
(2.1) also pairs `S` with its complement, so the centered spectrum is reciprocal:
`lambda_S lambda_{S^c}=g^2`, and all odd centered cumulants vanish.

## 3. Exact cumulant normalization

[LEMMA] From (2.1),

\[
\log\mathbb E(e^{zX})
 =\sum_{i=1}^n\log\cosh(zE_i).
\]

The Taylor expansion

\[
\log\cosh w
 =\frac{w^2}{2}-\frac{w^4}{12}+\frac{w^6}{45}
  -\frac{17w^8}{2520}+O(w^{10})
\]

therefore gives

\[
\boxed{
 p_1=\kappa_2,
 \quad p_2=-\frac{\kappa_4}{2},
 \quad p_3=\frac{\kappa_6}{16},
 \quad p_4=-\frac{\kappa_8}{272}.}
 \tag{3.1}
\]

In particular

\[
\kappa_2\ge0,\qquad \kappa_4\le0,
\qquad \kappa_6\ge0,\qquad \kappa_8\le0. \tag{3.2}
\]

These signs alone are weak; both obstructed grids below satisfy them.

For reproducibility, the cumulants are evaluated from centered empirical moments
`m_r=E[X^r]` using

\[
\begin{aligned}
\kappa_2={}&m_2,\\
\kappa_4={}&m_4-3m_2^2,\\
\kappa_6={}&m_6-15m_4m_2+30m_2^3,\\
\kappa_8={}&m_8-28m_6m_2-35m_4^2+420m_4m_2^2-630m_2^4.
\end{aligned} \tag{3.3}
\]

## 4. Newton and Hankel necessities

### 4.1 Newton inequalities

[LEMMA] Let `e_k=e_k(y_1,...,y_n)` with `e_0=1`.  Newton's identities read

\[
k e_k=\sum_{j=1}^k(-1)^{j-1}e_{k-j}p_j. \tag{4.1}
\]

Because every `y_i>=0`, every `e_k>=0`, and because there are exactly `n`
variables, `e_k=0` for `k>n`.  The first nontrivial exact necessities are

\[
\begin{aligned}
p_1^2-p_2&=2e_2\ge0,\\
p_1^3-3p_1p_2+2p_3&=6e_3\ge0,\\
p_1^4-6p_1^2p_2+3p_2^2+8p_1p_3-6p_4&=24e_4\ge0.
\end{aligned} \tag{4.2}
\]

The producer evaluates all three numerators.  They remain strictly positive for
both controls and both grids, so Newton positivity does not furnish the decisive
sign in these cases.

### 4.2 Stieltjes/Hankel positivity

[LEMMA] Set `p_0=n`.  For every finite real vector `c=(c_0,...,c_r)`,

\[
\begin{aligned}
\sum_{a,b=0}^r c_a c_b p_{a+b}
 &=\sum_{i=1}^n\left(\sum_{a=0}^r c_a y_i^a\right)^2\ge0,\\
\sum_{a,b=0}^r c_a c_b p_{a+b+1}
 &=\sum_{i=1}^n y_i\left(\sum_{a=0}^r c_a y_i^a\right)^2\ge0.
\end{aligned} \tag{4.3}
\]

Hence both Hankel matrices `[p_(a+b)]` and their one-step localizations
`[p_(a+b+1)]` are positive semidefinite.  Their rank is at most `n`, because
each is a sum of `n` rank-one matrices.  Three low-order `2x2` minors have the
explicit sum-of-squares forms

\[
\begin{aligned}
n p_2-p_1^2
 &=\sum_{i<j}(y_i-y_j)^2\ge0,\\
p_1p_3-p_2^2
 &=\sum_{i<j}y_i y_j(y_i-y_j)^2\ge0,\\
p_2p_4-p_3^2
 &=\sum_{i<j}y_i^2y_j^2(y_i-y_j)^2\ge0.
\end{aligned} \tag{4.4}
\]

Combining (3.1) and (4.4) yields

\[
\boxed{8\kappa_2\kappa_6-17\kappa_4^2
      =8p_1p_3-17p_2^2}
 =\sum_{i<j}y_i y_j(8y_i-17y_j)^2, \tag{4.5}
\]

which is not a sign invariant and is not used.  The converted Hankel tests use
the Gram form, not products of cumulants with unrelated signs:

\[
\boxed{J_1:=\kappa_2\kappa_6-4\kappa_4^2
      =16(p_1p_3-p_2^2)\ge0,} \tag{4.6}
\]

\[
\boxed{J_2:=\kappa_4\kappa_8-68\kappa_6^2
      =544(p_2p_4-p_3^2)\ge0.} \tag{4.7}
\]

Only `J_1`, which needs cumulants through order six, is required for the theorem
below.  `J_2` is evaluated as a redundant stronger-order obstruction.  Both raw
cumulant forms are combinations of cumulants of the **same** order-2r sign
parity, so they have the stated necessary sign.

## 5. Why the physical bipartite spectra are positive reciprocal

Let

\[
A=\sum_{v\in V}X_v,\qquad
B=\sum_{uv\in E}Z_uZ_v,
\]

and use the symmetric physical layer operator

\[
S=e^{\alpha A/2}e^{\beta B}e^{\alpha A/2}
\]

at real physical couplings.  It is positive definite.  Also `tr A=tr B=0`, so
`det S=1`.

[LEMMA] If the layer graph is bipartite with one color class `C`, define

\[
U_C=\prod_{v\in C}X_v,\qquad Z_V=\prod_{v\in V}Z_v,
\qquad J=Z_VU_C.
\]

Every edge has exactly one endpoint in `C`, so conjugation by `U_C` fixes `A`
and negates `B`.  Conjugation by `Z_V` negates `A` and fixes `B`.  Therefore

\[
JSJ^{-1}=e^{-\alpha A/2}e^{-\beta B}e^{-\alpha A/2}=S^{-1}. \tag{5.1}
\]

Thus the physical spectrum is positive and reciprocal, including
multiplicities.  This exact graph argument, not a floating pairing test, handles
all four finite cases used here.

For exact arithmetic set

\[
t=\tanh(\alpha/2),\qquad q=e^{2\beta}=\frac{1+t^2}{2t}.
\]

Up to a positive scalar, which determinant centering removes,

\[
R(t)=P_t\,\operatorname{diag}\!\left(q^{(b_s-\varepsilon)/2}\right)P_t,
\quad (P_t)_{rs}=t^{\operatorname{Ham}(r,s)}, \tag{5.2}
\]

where `b_s` is the in-layer bond energy and all `b_s` have the bond-count parity
`epsilon`.  At the rational physical point

\[
t_0=\frac13,\qquad q(t_0)=\frac53, \tag{5.3}
\]

all entries of `R` are rational.  The comparison value `K_c` is neither used nor
consulted to select or tune (5.3).

## 6. Full-spectrum certificate; floating point carries no claim

For each finite graph, let `D` be the exact least common multiple of all entry
denominators of `R(t_0)` and set `M=DR`, an integer symmetric matrix.  The
producer uses float64 diagonalization only to propose a matrix `V` and sorted
diagonal `d`; it stores the integer-rounded proposal at scale `10^14`.
Everything after that proposal is exact.

Define, over the rationals,

\[
G=V^TV,\qquad C=V^TMV.
\]

Exact integer products and exact row sums establish

\[
\|G-I\|_2\le\|G-I\|_\infty\le f<1,
\qquad
\|C-\operatorname{diag}(d)\|_2
 \le\|C-\operatorname{diag}(d)\|_\infty\le e. \tag{6.1}
\]

The first inequality proves `G` positive definite and hence `V` invertible.
The eigenvalues of `M` are exactly the generalized eigenvalues of `(C,G)`.
Weyl's inequality for `C` and the generalized Courant--Fischer principle give,
for each ordered eigenvalue, including every multiplicity,

\[
\frac{d_i-e}{1+f}\le\lambda_i(M)\le\frac{d_i+e}{1-f}. \tag{6.2}
\]

Division by `D` encloses every ordered eigenvalue of `R`.  No spectral simplicity
or separation is assumed.  In particular, the exact certificate remains valid
for the symmetry-forced degeneracies in the rectangular grids.

The exact integral-matrix anchors are:

| case | `N` | `D` | SHA-256 of row-major integer `M` |
|---|---:|---:|---|
| open chain `P4` | 16 | 492075 | `71101318dcb0eade93f96cb3ac02d351cc7a56ef95cb12c3d9fc51d43a9109d2` |
| open chain `P6` | 64 | 597871125 | `ec9c4c5001b3a4934628a10fbbde722261a1d4997488addd8661733fd3043944` |
| open grid `2x2` | 16 | 1476225 | `90af9b82f62591144b99a41463a2ce69d1c63fe6e3339f2d971fc10fb67cb8b3` |
| open grid `2x3` | 64 | 8968066875 | `ca37564f5c8e4c2c1efb6f08497bc5d625d6b667ad9fced3ca016cfe28f7c9ee` |

The independent verifier reconstructs (5.2) directly from its entry sum,
recomputes `D` and every integer in (6.1), and then derives (6.2) without
importing the producer.

## 7. Certified logarithms and cumulants

Each rational eigenvalue enclosure from (6.2) is converted to a logarithm
interval without a library transcendental.  For a positive rational `x`, choose
an integer `k` so that `x=2^k y` with `1<=y<2`, and set

\[
z=\frac{y-1}{y+1}\in[0,1/3].
\]

Then

\[
\log y=2\sum_{j=0}^{m-1}\frac{z^{2j+1}}{2j+1}+T_m,
\qquad
0\le T_m\le
\frac{2z^{2m+1}}{(2m+1)(1-z^2)}
\le\frac{9z^{2m+1}}{4(2m+1)}. \tag{7.1}
\]

The producer uses `m=48`.  Every rational input is rounded outward to a common
`2^-160` grid; every addition, multiplication, integer power, rational scaling,
and tail addition is then rounded outward with integer floor/ceiling division.
`log 2` is certified by the same formula at `z=1/3`.  Monotonicity gives
`log[lambda_lo,lambda_hi]=[log(lambda_lo),log(lambda_hi)]`.  The mean log,
centered log slots, moments, (3.3), (3.1), and (4.2)--(4.7) are all propagated by
the same outward interval arithmetic.

The verifier uses a separately written implementation with `2^-192` endpoints
and `m=64`, starting from the independently reconstructed eigenvalue enclosures.
It ignores every approximate decimal in the artifact.

## 8. Exact point evaluation

[COMPUTATION] At `t=1/3`, the decisive normalized Hankel quantity is:

| case | certified interval for `p1*p3-p2^2` | result |
|---|---:|---|
| open chain `P4` | `[0.1905423981203252, 0.19054240164388778]` | strict non-violation |
| open chain `P6` | `[0.46926576632140665, 0.4692664689158749]` | strict non-violation |
| open grid `2x2` | `[-0.06403871112577236, -0.06403870542584279]` | strict violation |
| open grid `2x3` | `[-0.16346806229912153, -0.16346717290955828]` | strict violation |

The displayed decimals are only readable renderings.  For example, the exact
`2x3` interval is

\[
\frac1{2^{160}}
[-238908840701476159867803973525457574463607181484,
 -238907540857173253730947686396928296311422392346], \tag{8.1}
\]

whose upper integer is strictly negative.  This is the reproducible sign
certificate.

For comparison with raw cumulants, the `2x3` enclosures are

\[
\begin{aligned}
\kappa_2&\in[3.4672496186767248,3.4672496198096954],\\
\kappa_4&\in[-4.505539043018897,-4.505538994894833],\\
\kappa_6&\in[22.66466159591826,22.664665192416123].
\end{aligned}
\]

Substitution in (4.6), performed by interval arithmetic rather than rounded
decimals, makes `J_1=16(p1*p3-p2^2)` strictly negative.  The higher-order minor
`p2*p4-p3^2` is also strictly negative for both grids, while all four normalized
power sums and all three displayed Newton numerators are strictly positive.
Therefore the failure is specifically a Hankel/Stieltjes moment failure, not a
mere cumulant-sign failure.

Satisfying finitely many necessities does not prove Gaussianity.  Accordingly,
the two chain rows are controls and non-violations only.

## 9. Uniform rational physical interval

Point strictness implies an unspecified neighborhood by continuity, but the
artifact proves an explicit one.  Take

\[
I=\left[\frac{999999997}{3000000000},
         \frac{1000000003}{3000000000}\right]
 =\left[\frac13-10^{-9},\frac13+10^{-9}\right]. \tag{9.1}
\]

On `0<t<1`, `q(t)=(t+t^{-1})/2` is decreasing.  For every entry in (5.2), the
producer combines exact positive rational intervals for `t^h` and `q(t)^m`.
The resulting natural intervals enclose every entry of `R(t)` for every `t in I`,
including negative integer `m`.  If `eta` is the maximum exact row sum of the
entrywise deviations from `R(t_0)`, symmetry gives

\[
\|R(t)-R(t_0)\|_2\le\|R(t)-R(t_0)\|_\infty\le\eta. \tag{9.2}
\]

Weyl enlarges every point enclosure by `[-eta,+eta]`; all enlarged lower
endpoints remain positive.  Repeating the exact logarithm/cumulant interval
calculation gives:

| case | exact-method interval for `p1*p3-p2^2` over every `t in I` |
|---|---:|
| open grid `2x2` | `[-0.06469247890703354, -0.06338493944082486]` |
| open grid `2x3` | `[-0.21674380000430868, -0.11019377850467191]` |

Both upper endpoints are strictly negative.  The artifact stores their exact
`2^-160` numerators; in particular the `2x3` upper numerator is

```
-161048387708256849652060785569905921519021276529
```

and is not a floating-point comparison.

[THEOREM] For every physical `t` in the rational interval (9.1), the positive
reciprocal spectrum of each open `2x2` and open `2x3` bipartite Ising layer
operator is not a full `n`-mode subset-product spectrum.  Indeed, every such
subset-product spectrum must satisfy (4.6), whereas each certified upper bound
is strictly negative.

## 10. Artifacts and honest closure

The durable files are:

- producer: `experiments/e228_bipartite_log_moments.py`;
- machine certificate: `results/spectral/bipartite_log_moments.json`;
- independent verifier: `tests/test_bipartite_log_moments.py`.

[COMPUTATION] The producer records `17/17` established checks.  The recorded
run's peak RSS is about `37.9 MB`; process CPU is recorded in the artifact
metadata and is O(2 s).  The exact
matrix, basis-error, eigenvalue, dyadic logarithm, cumulant, and uniform
perturbation data are retained in the JSON rather than inferred from the tables
above.

[UNRESOLVED] No conclusion is drawn for `2x4`, for every coupling, or for an
all-size bipartite family.  This route supplies a decisive finite obstruction
and an explicit local coupling interval; it does not turn finite Hankel signs
into a monotonicity or extension theorem.
