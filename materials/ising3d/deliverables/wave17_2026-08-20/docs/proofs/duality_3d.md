# Exact finite-volume duality of the 3D Ising and \(\mathbb Z_2\) gauge models

## Statement and conventions

This note uses \(K_g=\beta J_g\) for the gauge coupling and \(K_*\) for the dual Ising coupling.  Put

\[
 v_g=\tanh K_g,
 \qquad x_*=e^{-2K_*},
 \qquad x_*=v_g,
 \qquad
 A(K_g)=\cosh K_g\,e^{-K_*}
        =\sqrt{\frac{\sinh(2K_g)}{2}}.
\]

Equivalently,

\[
 e^{-2K_*}=\tanh K_g
 \quad\Longleftrightarrow\quad
 \sinh(2K_*)\sinh(2K_g)=1.
\]

Let \(X\) be a finite rectangular 3D cubical cell complex.  All chains below have coefficients in
\(\mathbb F_2\):

\[
 0\longrightarrow C_3(X)\mathop{\longrightarrow}^{\partial_3}C_2(X)
 \mathop{\longrightarrow}^{\partial_2}C_1(X)
 \mathop{\longrightarrow}^{\partial_1}C_0(X)\longrightarrow0.
\]

Write \(n_i=\dim C_i\), \(b_i=\dim H_i(X;\mathbb F_2)\), and \(|s|\) for the Hamming weight of a chain.  Links are distinct 1-cells and plaquettes are distinct 2-cells.  Thus a periodic direction of length two has two parallel links, as required by the repository convention.  Gauge variables \(U_l\in\{\pm1\}\) live on links, and

\[
 Z_G(X;K_g)=\sum_{\{U_l\}}
 \exp\!\left(K_g\sum_{p\in C_2}\prod_{l\in\partial p}U_l\right)
\]

is the **raw** gauge partition function: gauge-equivalent link assignments are all included.

Choose one closed surface \(h\in\ker\partial_2\) from each class in

\[
 H_2(X;\mathbb F_2)=\ker\partial_2/\operatorname{im}\partial_3.
\]

For that class, put a dual Ising spin \(\sigma_c=(-1)^{a_c}\) on every cube and define

\[
 Z_I^{[h]}(X^*;K_*)
 =\sum_{a\in C_3}
   \exp\!\left[K_*\sum_{p\in C_2}(-1)^{h_p+(\partial_3a)_p}\right].
\]

On a closed torus every plaquette is incident on two cubes, so this is the ordinary nearest-neighbor Ising model with a representative antiperiodic seam.  On an open boundary a boundary plaquette is incident on one cube; its term couples that cube spin to a fixed exterior \(+\) spin.  This fixed exterior condition is not optional: it is what makes a domain wall a closed surface in the original finite box.

### Theorem (all finite-volume constants)

For every such finite complex and every positive coupling related as above,

\[
 \boxed{
 Z_G(X;K_g)
 =2^{\,n_1-b_3}\,A(K_g)^{n_2}
   \sum_{[h]\in H_2(X;\mathbb F_2)} Z_I^{[h]}(X^*;K_*) .
 }
 \tag{1}
\]

Equivalently,

\[
 Z_G(X;K_g)
 =2^{\,n_1-b_3-n_2/2}
  [\sinh(2K_g)]^{n_2/2}
  \sum_{[h]\in H_2} Z_I^{[h]}(X^*;K_*).
 \tag{2}
\]

The form (1) is preferable if \(n_2\) is odd because every positive square root is explicit.  If the gauge partition function is normalized by the volume of a local gauge orbit, then

\[
 Z_G^{\rm orb}=2^{-\operatorname{rank}\partial_1}Z_G,
\]

so the exponent of two in (1) becomes
\(n_1-b_3-\operatorname{rank}\partial_1\).  For a connected complex,
\(\operatorname{rank}\partial_1=n_0-1\).

## Direct derivation

### 1. Gauge high-temperature expansion gives closed surfaces

Write \(U_l=(-1)^{u_l}\), \(u\in C_1\).  For a plaquette \(p\),

\[
 \prod_{l\in\partial p}U_l=(-1)^{\langle u,\partial_2p\rangle}.
\]

Using \(e^{K_gF}=\cosh K_g(1+v_gF)\) for \(F=\pm1\), expand once for every plaquette:

\[
\begin{aligned}
 Z_G
 &= (\cosh K_g)^{n_2}
    \sum_{u\in C_1}\prod_{p\in C_2}
    \left[1+v_g(-1)^{\langle u,\partial_2p\rangle}\right]\\
 &= (\cosh K_g)^{n_2}
    \sum_{s\in C_2}v_g^{|s|}
    \sum_{u\in C_1}(-1)^{\langle u,\partial_2s\rangle}.
\end{aligned}
\]

Character orthogonality on \(C_1\) gives

\[
 \sum_{u\in C_1}(-1)^{\langle u,\partial_2s\rangle}
 =\begin{cases}
   2^{n_1},&\partial_2s=0,\\
   0,&\partial_2s\ne0.
  \end{cases}
\]

Therefore, with the exact closed-surface polynomial

\[
 S_X(v)=\sum_{s\in\ker\partial_2}v^{|s|},
\]

one has

\[
 \boxed{Z_G=2^{n_1}(\cosh K_g)^{n_2}S_X(v_g).}
 \tag{3}
\]

The factor \(2^{n_1}\) is the sum over **all raw link variables**, not the gauge-orbit volume.

### 2. Dual Ising low-temperature expansion and its multiplicity

For fixed \([h]\), the set of unsatisfied signed dual Ising interactions is exactly

\[
 s=h+\partial_3a.
\]

It is closed because \(\partial_2h=0\) and \(\partial_2\partial_3=0\).  Consequently

\[
 Z_I^{[h]}(K_*)
 =e^{K_*n_2}Q_{[h]}(x_*),
 \qquad
 Q_{[h]}(x)=\sum_{a\in C_3}x^{|h+\partial_3a|}.
 \tag{4}
\]

The map

\[
 ([h],a)\longmapsto h+\partial_3a\in\ker\partial_2
\]

is onto.  A given closed surface determines its homology class uniquely, and the solutions for
\(a\) form a coset of \(\ker\partial_3\).  Every surface is therefore counted exactly

\[
 |\ker\partial_3|=2^{b_3}
\]

times.  This proves the exact **integer polynomial identity**

\[
 \boxed{
 \sum_{[h]\in H_2}Q_{[h]}(x)=2^{b_3}S_X(x)
 }
 \tag{5}
\]

coefficient by coefficient.  Combining (3)--(5) at \(x_*=v_g\) gives

\[
 Z_G
 =2^{n_1-b_3}(\cosh K_g)^{n_2}e^{-K_*n_2}
  \sum_{[h]}Z_I^{[h]},
\]

which is (1).  No citation or thermodynamic-limit argument enters this derivation.

### 3. The exact DOS conversion used in the computation

If \(g_G(b)\) is the number of raw link assignments with exactly \(b\) satisfied plaquettes, then

\[
 Z_G(K_g)=e^{-K_gn_2}\sum_{b=0}^{n_2}g_G(b)e^{2K_gb}.
\]

Expanding the same configurations in \(v=\tanh K_g\) gives

\[
 S_X(v)=2^{-n_1}\sum_{b=0}^{n_2}g_G(b)(1+v)^b(1-v)^{n_2-b}.
 \tag{6}
\]

Every coefficient on the right was computed as a Python integer and checked to be divisible by
\(2^{n_1}\).  The resulting polynomial was compared with a separate nullspace enumeration of
\(\ker\partial_2\).  It was then compared with the independently enumerated cube-spin polynomials
through (5).  Thus the check does not merely evaluate the two partition functions at selected floating-point couplings.

## Boundary conditions, topological sectors, and the factors that are easy to miss

### Open box

For an \(L_x\times L_y\times L_z\) box of gauge vertices,

\[
\begin{aligned}
 n_0&=L_xL_yL_z,\\
 n_1&=(L_x-1)L_yL_z+L_x(L_y-1)L_z+L_xL_y(L_z-1),\\
 n_2&=(L_x-1)(L_y-1)L_z+(L_x-1)L_y(L_z-1)
      +L_x(L_y-1)(L_z-1),\\
 n_3&=(L_x-1)(L_y-1)(L_z-1).
\end{aligned}
\]

The box is contractible: \((b_0,b_1,b_2,b_3)=(1,0,0,0)\).  Equation (1) contains one dual sector and no global-spin factor.  The dual model has \(n_3\) spins and fixed \(+\) exterior couplings across every boundary plaquette.  Replacing those by free Ising boundary conditions gives the wrong finite-volume identity.

### Three-torus

For a fully periodic \(L_x\times L_y\times L_z\) lattice with \(N=L_xL_yL_z\),

\[
 (n_0,n_1,n_2,n_3)=(N,3N,3N,N),
 \qquad (b_0,b_1,b_2,b_3)=(1,3,3,1).
\]

There are \(2^{b_2}=8\) dual Ising sectors.  They may be represented by adding any combination of the three noncontractible coordinate sheets; on the dual lattice these are the periodic/antiperiodic twists in the three directions.  The kernel \(\ker\partial_3\) has size two because flipping all dual spins changes no domain wall.  Hence

\[
 \boxed{
 Z_G^{T^3}(K_g)=2^{3N-1}A(K_g)^{3N}
 \sum_{\tau\in\mathbb F_2^3}Z_I^\tau(K_*).
 }
 \tag{7}
\]

The local gauge orbit has size \(2^{N-1}\).  Flat gauge fields have multiplicity

\[
 2^{n_1-\operatorname{rank}\partial_2}
 =2^{N+2}=2^{N-1}\times2^3.
 \tag{8}
\]

The last \(2^3\) is the set of gauge-inequivalent Wilson/Polyakov (torelon) holonomies in
\(H_1(T^3;\mathbb F_2)\).  It is already present in the raw link sum; multiplying (7) by another eight would double-count it.  Conversely, replacing the raw link sum by one representative per gauge orbit requires division by exactly \(2^{N-1}\), not by \(2^N\), because the global vertex transformation acts trivially.

### Mixed boundaries

If exactly \(p\) directions are periodic and the others are open, the topology is
\(T^p\times I^{3-p}\), so

\[
 b_k=\binom pk.
\]

The counts of gauge holonomies and dual surface sectors are therefore generally different:
\(2^{b_1}=2^p\) flat gauge orbits but \(2^{b_2}=2^{\binom p2}\) dual sectors.  The exact
\(2\times2\times2\) two-periodic/one-open check has four flat gauge orbits but only two dual Ising sectors.  Treating these as a single unnamed “topological factor” fails.

## Exact finite-volume results

The command

```text
.venv/bin/python experiments/e09_duality.py
```

produced `results/duality/duality.json`.  The following are exact integer results; omitted coefficients are zero.

| Gauge lattice | \((n_0,n_1,n_2,n_3)\) | \((b_1,b_2,b_3)\) | \(S_X(v)\) | Verified dual relation |
|---|---:|---:|---|---|
| \(2\times2\times2\) open | \((8,12,6,1)\) | \((0,0,0)\) | \(1+v^6\) | \(Q_0=S_X\) |
| \(3\times2\times2\) open | \((12,20,11,2)\) | \((0,0,0)\) | \(1+2v^6+v^{10}\) | \(Q_0=S_X\) |
| \(2\times2\times2\), POO | \((8,16,10,2)\) | \((1,0,0)\) | \(1+2v^6+v^8\) | \(Q_0=S_X\) |
| \(2\times2\times2\), PPO | \((8,20,16,4)\) | \((2,1,0)\) | \(1+2v^4+4v^6+17v^8+4v^{10}+4v^{12}\) | \(Q_0+Q_1=S_X\) |
| \(2\times2\times2\) torus | \((8,24,24,8)\) | \((3,3,1)\) | \(1+6v^4+8v^6+111v^8+120v^{10}+532v^{12}+120v^{14}+111v^{16}+8v^{18}+6v^{20}+v^{24}\) | \(\sum_{\tau=0}^7Q_\tau=2S_X\) |

For the \(2^3\) torus the raw gauge DOS, indexed by the number of satisfied plaquettes, is

\[
[1024,0,0,0,43008,0,229376,0,1588224,0,3440640,0,6172672,
0,3440640,0,1588224,0,229376,0,43008,0,0,0,1024].
\]

It sums to \(2^{24}\), as a 24-link brute-force enumeration must.  The all-satisfied coefficient is
\(1024=2^7\times8\), independently exposing the gauge-orbit and Polyakov factors in (8).

A useful failed-naive check is the torus trivial Ising sector alone:

\[
 Q_0(x)=2+16x^6+30x^8+48x^{10}+64x^{12}
        +48x^{14}+30x^{16}+16x^{18}+2x^{24}.
\]

It cannot equal \(S_X\) or \(2S_X\): for example, its \(x^4\) coefficient vanishes whereas
\([x^4]2S_X=12\).  Only the sum of all eight twists gives (5).  This is the finite-volume correction that a schematic “closed surfaces equal domain walls” statement misses.

## Two-dimensional self-duality control

For a square \(N=L_xL_y\) torus, let \(\alpha=(\alpha_x,\alpha_y)\in\mathbb F_2^2\) label the signs of the two wrap seams.  Let \(Q_\beta(x)\) be the low-temperature broken-bond polynomial in sector \(\beta\), and let \(R_\alpha(v)\) be the signed high-temperature even-subgraph polynomial in sector \(\alpha\).  With sector order
\((00,10,01,11)\), exact Kramers--Wannier duality is the four-sector transform

\[
 2R_\alpha(v)=\sum_\beta H_{\alpha\beta}Q_\beta(v),
 \qquad
 H=\begin{pmatrix}
 1&1&1&1\\
 1&1&-1&-1\\
 1&-1&1&-1\\
 1&-1&-1&1
 \end{pmatrix}.
 \tag{9}
\]

The crosswise signs are
\(H_{\alpha\beta}=(-1)^{\alpha_x\beta_y+\alpha_y\beta_x}\), because a seam and a dual winding cycle intersect in the perpendicular direction.  In partition-function form,

\[
 \boxed{
 Z_\alpha(K)=\frac{[\sinh(2K)]^N}{2}
 \sum_\beta H_{\alpha\beta}Z_\beta(K_*),
 \qquad e^{-2K_*}=\tanh K.
 }
 \tag{10}
\]

The experiment verified every coefficient in all four rows of (9), separately on the
\(4\times4\) and \(5\times5\) tori.  Each sector polynomial sums to \(2^N\), and the untwisted
high-temperature polynomial obeys \(R_{00}(1)=2^{2N-N+1}=2^{N+1}\).  The untwisted
broken-bond polynomial also agrees with the repository's independent transfer matrix on both
sizes.  On the \(4\times4\) torus the calculation reproduces the known 33-coefficient
even-subgraph polynomial exactly.  The \(5\times5\) torus contains odd-length noncontractible
cycles, so it additionally detects the erroneous assumption that every torus polynomial has only
even powers.

At the self-dual coupling,

\[
 K=K_*=\frac12\log(1+\sqrt2),\qquad \sinh(2K)=1,
\]

(10) still mixes four sectors; even in 2D it is incorrect at finite volume to write
\(Z_{PP}(K)=Z_{PP}(K_*)\) without the sector transform.  What makes the 2D argument close is that the dual objects are again ordinary 0-form Ising spins on a square lattice, so the same model occurs on both sides.

## Does 3D duality close an equation for the Ising free energy?

**No.**  On the periodic cubic lattice, divide the logarithm of (7) by \(N\) and take the thermodynamic limit.  A finite number of twist sectors has the same bulk free-energy density, so for the raw gauge normalization

\[
 \phi_G(K_g)
 =\phi_I(K_*)+\frac32\log\!\left(2\sinh(2K_g)\right),
 \qquad \sinh(2K_*)\sinh(2K_g)=1.
 \tag{11}
\]

This is one equation relating two different unknown functions, \(\phi_G\) and \(\phi_I\).  The gauge constraint \(\partial_2s=0\) is precisely the condition already used to obtain (3); imposing it again supplies no second equation.

Several apparent ways to close the system do not work:

1. **Set \(K_g=K_*\).**  This would give \(\sinh(2K)=1\), but it equates the coupling coordinates of two different theories.  A 3D 0-form Ising model is not the same statistical model as a 3D 1-form gauge theory, so this is not a symmetry and does not determine the 3D Ising critical point.
2. **Dualize the gauge theory again.**  The coupling map \(D\), defined by \(e^{-2D(K)}=\tanh K\), obeys \(D(D(K))=K\).  Solving (11) for \(\phi_I\) gives
   \[
   \phi_I(K)=\phi_G(DK)+\frac32\log\!\left(\frac{\sinh(2K)}{2}\right).
   \]
   Substitution of (11) at \(DK\) cancels the two logarithms and returns the identity \(\phi_I(K)=\phi_I(K)\).  The second duality is the inverse of the first, not an independent relation.
3. **Use local gauge symmetry or the flat-holonomy sectors.**  Local gauge transformations explain the orbit factor, and \(H_1\) labels Polyakov sectors.  Neither changes \(K_g\) or relates \(\phi_G(K_g)\) to \(\phi_G\) at another coupling.
4. **Invoke electric--magnetic duality as an extra map.**  For \(\mathbb Z_2\) lattice models in dimension \(d\), the local Fourier/cell duality changes a \(p\)-form model into a \((d-p-2)\)-form model.  Thus 2D spins (\(p=0\)) are self-dual, and 4D link gauge theory (\(p=1\)) is self-dual.  In 3D, however, \(p=0\leftrightarrow1\): Ising and gauge exchange, and composition merely returns to the starting theory.

The exact consequence for critical couplings is only

\[
 \sinh(2K_{I,c})\sinh(2K_{G,c})=1,
\]

provided the corresponding transitions are mapped into each other.  Without an independent exact value or an independent coupling-changing symmetry of either model, neither critical coupling is fixed.  The obstruction is therefore structural, not a missing gauge-orbit or topological factor: after all finite-volume factors are restored, duality still exchanges two different model classes.

## Reproduction

```text
.venv/bin/python tests/test_duality.py
.venv/bin/python experiments/e09_duality.py
```

The standalone test checks five 3D lattices (including the \(3\times2\times2\) open box and the
24-link \(2^3\) torus) and both 2D controls with exact integer equality, then prints `PASS`.
