# 3D-Ising literature matrix

Status counts in `sources/manifest.yaml`: **26 full-text**, **1 abstract/partial-text**, and **13 metadata-only** records (40 total). The 26 full-text records correspond to the 26 locally hashed PDFs. `reproduced_locally` is false throughout: this file records literature extraction, not independent reproduction.

## Matrix

| key | year | obtained | key result used | exact location | reproduced_locally |
|---|---:|---|---|---|---|
| onsager1944 | 1944 | metadata_only | Exact 2D free energy; $\sinh 2K_x\sinh 2K_y=1$ | pp. 117-149; equations unobtained | false |
| kaufman1949 | 1949 | abstract_only | Spinor reduction and four torus sectors | UCLA Eqs. (1)-(9); later formula unobtained | false |
| leeyang1952 | 1952 | metadata_only | Ferromagnetic field zeros lie on $|z|=1$ | pp. 410-419; theorem number unobtained | false |
| yang1952 | 1952 | metadata_only | $M=[1-\sinh^{-4}(2K)]^{1/8}$ | pp. 808-816; equation number unobtained | false |
| kacward1952 | 1952 | metadata_only | Planar Ising determinant/Pfaffian | pp. 1332-1337; equations unobtained | false |
| schultzmattisl1964 | 1964 | metadata_only | 2D transfer matrix as quadratic fermions | pp. 856-871; equations unobtained | false |
| istrail2000 | 2000 | full | NP-completeness on non-planar lattice families | main theorems / extended abstract | false |
| zamolodchikov1980 | 1980 | metadata_only | Tetrahedron equation | JETP 52, 325-336; equation number unobtained | false |
| zamolodchikov1981 | 1981 | metadata_only | Tetrahedron consistency for string scattering | CMP 79, 489-505 | false |
| baxter1983 | 1983 | metadata_only | Analysis of Zamolodchikov solution | CMP 88, 185-205 | false |
| baxter1986 | 1986 | metadata_only | Yang-Baxter/tetrahedron model relation | Physica D 18, 321-347 | false |
| ferrenberg2018 | 2018 | full | $K_c=0.221654626(5)$ and exponents | abstract; Eqs. (39), (46), (47); Table X | false |
| campostrini2002 | 2002 | full | 25th-order HT data and improved exponents | Eq. (2.4), Sec. II, Table 2 | false |
| butera2011 | 2011 | full | Bivariate free energy through $K^{24}$ | Sec. II; Eqs. (5), (8), (9), (13), (14) | false |
| kos2016 | 2016 | full | Precision 3D-Ising CFT island | Eqs. (3.1)-(3.3) | false |
| arisuefujiwara2002 | 2002 | full | Exact free-energy coefficients through $t^{46}$ | Table 1 | false |
| arisuefujiwara2003 | 2002 | full | Free energy to $\beta^{46}$, susceptibility to $\beta^{32}$ | abstract and algorithm sections | false |
| arisuetabata1994 | 1994 | full | LT second-moment series through $u^{26}$ | abstract and tables | false |
| bhanot1990 | 1990 | metadata_only | DOS/histogram reference | metadata only | false |
| pearson1982 | 1982 | metadata_only | Exact periodic $4^3$ partition function exists | PRB 26, 6285; table unobtained | false |
| valani2011 | 2011 | full | Exact finite transfer-matrix polynomials | Eqs. (C.0.1)-(C.0.2), Table E.1 | false |
| chapman2020 | 2020 | full | Injective free-fermion map iff frustration graph is a line graph | Def. 1; Thm. 1; Eqs. (25), (78) | false |
| elman2021 | 2021 | full | ECF implies explicit free-fermion solution | Def. 3; Thms. 1-2; Eq. (5) | false |
| chapman2023 | 2023 | full | Claw-free plus simplicial clique is sufficient | abstract and main theorem | false |
| wiersema2023 | 2023 | full | 1D 2-local DLA classification; TFIM $\mathfrak{so}$ algebras | Thms. IV.1-IV.2; Lemmas C.6, C.8 | false |
| wosiek1982 | 1982 | metadata_only | Higher-D local fermion representation needs constraints | pp. 543-551; equations unobtained | false |
| ball2005 | 2005 | full | Local spin description with conserved auxiliaries | construction in arXiv paper | false |
| verstraete2005 | 2005 | full | Auxiliary Majoranas preserve locality | mapping construction | false |
| chen2018 | 2018 | full | Exact 2D bosonization to a gauge theory | duality/Gauss-law construction | false |
| chen2019 | 2019 | full | 3D fermions map to a 2-form $\mathbb Z_2$ gauge theory | main mapping/Gauss law | false |
| wegner1971 | 1971 | full | Generalized Ising/gauge duality | reprint arXiv:1411.5815 | false |
| balian1975 | 1975 | metadata_only | Gauge-invariant Ising dual setting | PRD 11, 2098-2103 | false |
| zhang2007 | 2007 | full | Conjectured golden-ratio $K_c$ and fourfold integral | Eqs. (3.37), (3.54), (3.60), (3.62) | false |
| wu2008comment | 2008 | full | HT/LT series disproof | pp. 2-5 | false |
| zhang2008response | 2008 | full | Defense of two conjectures | abstract and response body | false |
| wu2008rejoinder | 2008 | full | Eleven fitted terms are not predictions | entire two-page rejoinder | false |
| perk2012 | 2012 | full | Unknown weights make integral non-solution | Sec. 3.1 | false |
| zhang2021 | 2021 | full | Claimed operator-algebra solution gives $H_c=0.30468893$ | Eqs. (15)-(17), (32), (44)-(45) | false |
| perk2022comment | 2023 | full | Boundary/operator/rotation tests refute claim | Eqs. (1)-(3) and numerical test | false |
| zhang2023reply | 2023 | full | Reply and appended counter-reply | pp. 1-2 | false |

## Source-by-source implementation notes

## onsager1944 — Crystal Statistics I

For anisotropic square-lattice couplings $K_x=\beta J_x$ and $K_y=\beta J_y$, the critical manifold is

$$\sinh(2K_x)\sinh(2K_y)=1.$$

This is a two-dimensional thermodynamic-limit control result. The APS PDF returned HTTP 403, so the original displayed free-energy equation and number are deliberately not transcribed here.

## kaufman1949 — Crystal Statistics II

The accessible UCLA transcription establishes $Z=\operatorname{tr}(V^m)$ (Eq. (3)), the factors $V_1,V_2$ (Eq. (4)), the Pauli/quaternion generators (Eqs. (5)-(6)), the dual coupling (Eq. (7)), and the scalar-normalized transfer matrix (Eqs. (8)-(9)). The later finite-torus formula is image-only and was not recovered.

A commonly used implementation organization (recorded here as a **cross-check template, not a verbatim Kaufman quotation**) defines a signed set $\gamma_r$ from

$$\cosh\gamma_r=\cosh(2K_x^*)\cosh(2K_y)-\sinh(2K_x^*)\sinh(2K_y)\cos(\pi r/n),$$
$$\sinh(2K_x)\sinh(2K_x^*)=1,$$

and forms four products from $2\cosh(m\gamma_{2r+1}/2)$, $2\sinh(m\gamma_{2r+1}/2)$, $2\cosh(m\gamma_{2r}/2)$, and $2\sinh(m\gamma_{2r}/2)$. The sector signs and the sign convention for $\gamma_0$ must be checked against a scan of the original formula before this can be used as an oracle. No equation number or claimed verbatim formula is fabricated.

## leeyang1952 — Lee-Yang circle theorem

For a finite Ising ferromagnet with pair interactions $J_{ij}\ge0$, write the field dependence as a polynomial in

$$z=e^{-2\beta h}.$$

The circle theorem places every zero on $|z|=1$. Therefore, for real $h\ne0$ (positive real $z\ne1$), the finite-volume partition function is nonzero; thermodynamic singularities driven by these zeros can reach the real-field axis only at $h=0$. The exact theorem numbering was not obtained from APS.

## yang1952 — spontaneous magnetization

For the isotropic square lattice in zero field,

$$M(K)=\begin{cases}\left[1-\sinh^{-4}(2K)\right]^{1/8},&K>K_c,\\0,&K\le K_c.\end{cases}$$

This is the decisive dimensional-reduction check used by Wu et al.: a proposed 3D expression that reduces to 2D must recover exponent $1/8$, not $3/8$.

## kacward1952 — planar determinant

For a planar graph, introduce directed edges and the non-backtracking matrix $T$ whose transition $e\to e'$ carries $\tanh(\beta J_e)$ and the half-turning-angle phase, with immediate reversal forbidden. The implementable structural identity is

$$Z=2^{|V|}\prod_{e\in E}\cosh(\beta J_e)\sqrt{\det(I-T)}.$$

Planarity controls the cancellation of intersecting loops. Applying this determinant unchanged to a non-planar 3D lattice is not justified.

## schultzmattisl1964 — many-fermion formulation

The row transfer matrix is Jordan-Wigner transformed to Majorana/fermion operators. Solvability follows because the exponent is quadratic, so diagonalization reduces to a one-particle orthogonal/Bogoliubov problem. The implementation invariant to test in any proposed 3D reduction is closure in the quadratic Majorana algebra; products of four or more independent Majoranas are interactions, not free hopping.

## istrail2000 — complexity boundary

The result concerns exact finite-instance partition-function computation over non-planar lattice families and establishes NP-completeness in that setting. It does **not** state that the translationally invariant, homogeneous, zero-field simple-cubic thermodynamic-limit free energy is NP-complete, nor does it logically forbid a special closed form. Use it as a scope constraint against claims of a generic non-planar Pfaffian algorithm.

## zamolodchikov1980 — tetrahedron equation

In a common tensor-index convention the consistency equation has the form

$$R_{123}R_{145}R_{246}R_{356}=R_{356}R_{246}R_{145}R_{123}.$$

The placement convention varies, and the original PDF was not obtained; treat this as the algebraic shape, not a certified original equation transcription. A credible 3D integrability claim must specify the local $R$ operator and verify this identity, not merely invoke an extra dimension.

## zamolodchikov1981 — straight-string scattering

The paper applies tetrahedron consistency to factorized scattering of straight strings in $2+1$ dimensions. Its relevance is methodological: 3D consistency is an explicit operator identity. It is not a solution of the nearest-neighbor simple-cubic Ising partition function.

## baxter1983 — tetrahedron solution analysis

Baxter analyzes Zamolodchikov's explicit solution and the associated commuting-transfer structure. No formula from this metadata-only source is used as a numerical Ising oracle.

## baxter1986 — Yang-Baxter and Zamolodchikov model

This paper relates ordinary Yang-Baxter reductions to the Zamolodchikov model. It supplies a comparison class of genuine 3D integrability, but no identification with the standard 3D Ising model was found.

## ferrenberg2018 — Monte Carlo benchmark

The primary numerical targets are

$$K_c=0.221654626(5),\qquad \nu=0.629912(86),$$
$$\gamma=1.23708(33),\qquad \beta=0.32630(22),\qquad \gamma/\nu=1.96390(45).$$

Parentheses are one-standard-uncertainty digits. Locations: $\nu$ Eq. (39), $\gamma$ Eq. (46), $\beta$ Eq. (47), and $\gamma/\nu$ in Sec. III C. Table X compares HT-series and bootstrap values. A candidate exact solution should predict $K_c$ within these uncertainties before any exponent fitting.

## campostrini2002 — 25th-order high-temperature series

Improved-model estimates in the abstract are

$$\gamma=1.2373(2),\ \nu=0.63012(16),\ \alpha=0.1096(5),\ \eta=0.03639(15),$$
$$\beta=0.32653(10),\ \delta=4.7893(8),\ \Delta=0.52(3)\ \text{(standard model)}.$$

For $v=\tanh\beta$, Eq. (2.4) prints, among other tails,

$$\chi(v)=\cdots+18554916271112254v^{24}+85923704942057238v^{25}+O(v^{26}),$$
$$m_2(v)=\cdots+977496788431483776v^{24}+47677378698515169334v^{25}+O(v^{26}).$$

These exact integer coefficients are direct symbolic regression tests; do not compare against a series in $\beta$ without converting variables.

## butera2011 — bivariate equation of state

The paper extends the high-temperature/low-field bivariate free energy through $K^{24}$ for simple-cubic and body-centered-cubic lattices and obtains higher field derivatives through the same order. The source explicitly defers full numerical coefficient tables to a separate publication. Thus this record contributes the scaling and equation-of-state structure (Eqs. (5), (8), (9), (13), and Eq. (14) onward), not fabricated tables.

## kos2016 — conformal-bootstrap island

The 3D-Ising CFT data are

$$(\Delta_\sigma,\Delta_\epsilon,\lambda_{\sigma\sigma\epsilon},\lambda_{\epsilon\epsilon\epsilon})
=(0.5181489(10),1.412625(10),1.0518537(41),1.532435(19)).$$

Use $\eta=2\Delta_\sigma-1$ and $\nu=1/(3-\Delta_\epsilon)$, giving $\eta=0.0362978(20)$ and $\nu=0.629971(4)$. These are universality-class checks, not direct finite-lattice partition coefficients.

## arisuefujiwara2002 — exact HT free-energy coefficients

The source defines $t=\tanh\beta$ and prints the following Table 1 coefficients $a_n$ (the displayed normalization line in the source should be checked because it prints `3 cosh(beta)` where the standard reduced free energy uses logarithms):

| $n$ | $a_n$ | $n$ | $a_n$ |
|---:|---:|---:|---:|
| 2 | 0 | 4 | 3 |
| 6 | 22 | 8 | $375/2$ |
| 10 | 1980 | 12 | 24044 |
| 14 | 319170 | 16 | $18059031/4$ |
| 18 | $201010408/3$ | 20 | $5162283633/5$ |
| 22 | 16397040750 | 24 | 266958797382 |
| 26 | 4437596650548 | 28 | $525549581866326/7$ |
| 30 | $6448284363491202/5$ | 32 | $179577198475709847/8$ |
| 34 | 395251648062268272 | 36 | $21093662188820520521/3$ |
| 38 | 126225408651399082182 | 40 | $4569217533196761997785/2$ |
| 42 | $291591287110968623857940/7$ | 44 | $8410722262379235048686604/11$ |
| 46 | 14120314204713719766888210 |  |  |

The same section reports $\alpha=0.1045(1)$ at $\beta_c=0.22165459(10)$ or $\alpha=0.1077(2)$ at $\beta_c=0.2216595(15)$, illustrating sensitivity to the imposed critical point.

## arisuefujiwara2003 — short algorithm presentation

This companion paper states that the finite-lattice method reaches $\beta^{46}$ for free energy and $\beta^{32}$ for susceptibility. Use the coefficient table in `arisuefujiwara2002` as the exact data source; this record documents the algorithm and susceptibility reach.

## arisuetabata1994 — low-temperature finite-lattice series

The paper computes the second moment of the 3D Ising correlation function through $u^{26}$ and the ASOS free energy through $u^{23}$. Its variable and observable must be preserved when comparing coefficients; it is not the spontaneous-magnetization series quoted by Wu et al.

## bhanot1990 — density of states

This is a histogram/DOS benchmark cited by later series work. No coefficient table or exact $4^3$ polynomial was obtained, so the manifest deliberately marks it metadata-only.

## pearson1982 — periodic $4\times4\times4$

Pearson reports an exact periodic $4^3$ partition function using symmetry to reduce enumeration from $2^{64}$ to $2^{32}$. The actual DOS/polynomial table was not legally retrieved in this run. It must not be substituted by Valani's different $5\times5\times10'$ boundary-condition table.

## valani2011 — finite transfer-matrix data

The thesis convention is

$$Z=a_Nx^N+a_{N-1}x^{N-1}+\cdots+a_0,\qquad x=e^\beta \tag{C.0.1}$$

and the layered computation is

$$Z=\operatorname{Tr}(T^{N_z}).\tag{C.0.2}$$

Table E.1 prints an exact $5\times5\times10'$ coefficient sequence. Before using it, recover the prime/boundary convention from the thesis and verify $\sum_i a_i=2^{N_xN_yN_z}$ as prescribed in Appendix C.

## chapman2020 — line-graph criterion

For $H=\sum_{j\in V}b_j\sigma^j$, the frustration graph $G(H)$ has a vertex for each nonzero Pauli term and an edge exactly when two terms anticommute. Definition 1 gives

$$L(R)=(E,F),\qquad F=\{(e_1,e_2):e_1,e_2\in E,\ |e_1\cap e_2|=1\}.\tag{24}$$

Theorem 1 (restated as Eq. (78)) says an **injective** map satisfying Eqs. (20)-(21), from each Pauli generator to a Majorana bilinear, exists iff

$$G(H)\simeq L(R).$$

This is necessary and sufficient only for that injective generator-to-generator mapping. The paper explicitly notes possible non-injective mappings after fixing stabilizer sectors.

## elman2021 — ECF sufficiency

A claw is $K_{1,3}$. A hole is an induced cycle, so an even hole is an induced $C_{2k}$ for $k\ge2$. Definition 3 calls $G$ ECF when neither a claw nor an even hole occurs as an induced subgraph.

Theorem 1: every ECF Hamiltonian has a free spectrum, with single-particle energies satisfying

$$P_G(-1/\epsilon_j^2)=0,$$

where $P_G$ is the vertex-weighted independence polynomial. Theorem 2: every ECF Hamiltonian is free-fermion-solvable via Eq. (2), with modes constructed by Definition 5's incognito-mode transfer operator. This is a **sufficient**, not necessary, criterion; the paper itself exhibits special non-ECF equal-coupling free spectra.

## chapman2023 — unified graph criterion

The sufficient condition is: the frustration graph is claw-free and contains a simplicial clique. Cycle symmetries commute, allowing the Fendley construction to be applied separately in their common eigenspaces. This contains both line-graph and ECF families but is not claimed here as a necessary criterion for every conceivable exact solution.

## wiersema2023 — dynamical Lie algebras

For the open-chain TFIM-equivalent family, generators can be taken as the orbit of `{XX, XZ}` (equivalently `{ZZ, IX}` after local relabeling). Supplemental Eq. (C19) yields a path frustration graph of $2n-2$ generators, and

$$\mathfrak a_8(n)\cong\mathfrak{so}(2n-1),\qquad \dim=(n-1)(2n-1).$$

For periodic boundaries, Lemma C.8 gives

$$\mathfrak a_8^\circ(n)\cong\mathfrak{so}(2n)^{\oplus2},\qquad \dim=2n(2n-1).$$

This classification is one-dimensional. The paper explicitly lists extension to 2D/3D topologies as future work; it does not establish a 3D-Ising free-fermion algebra.

## wosiek1982 — local higher-dimensional fermions

The construction maps higher-dimensional fermions to locally interacting spins with constraints. The key implementation lesson is that locality alone does not make the spin Hamiltonian unconstrained or quadratic; track the physical/gauge subspace explicitly. Full equations were not obtained.

## ball2005 — fermions without fermion fields

Ball represents arbitrary hopping Hamiltonians using locally identified operators that commute on different sites. Extra conserved degrees of freedom appear. For a proposed Ising fermionization, count these sectors and show which one reproduces the original spin Hilbert space.

## verstraete2005 — auxiliary Majoranas

Auxiliary Majorana modes convert Jordan-Wigner strings in $d>1$ into local interactions. The enlarged system comes with a gauge/stabilizer sector that must be fixed. Dropping the auxiliary constraint changes the spectrum and partition function.

## chen2018 — exact 2D bosonization

A 2D fermionic lattice model maps to a lattice gauge theory. On simply connected space the map is an equivalence after imposing the Gauss law. On nontrivial topology, global sectors also matter. This is a precise counterexample to the idea that higher-dimensional local fermionization should produce an unconstrained ordinary spin model.

## chen2019 — 3D bosonization

In three spatial dimensions, fermions map locally to a two-form $\mathbb Z_2$ gauge field with an unusual Gauss law. Any use as an Ising tool would have to derive that gauge theory and its constraints; it does not directly diagonalize the 3D nearest-neighbor Ising Hamiltonian.

## wegner1971 — generalized Ising duality

Wegner's generalized interactions produce dual models and include gauge-invariant examples with nonlocal order diagnostics. The ordinary 3D nearest-neighbor Ising model is dual to a different gauge-type interaction structure, not to itself, so duality alone does not fix $K_c$ by a self-dual equation.

## balian1975 — gauge-invariant Ising model

This source studies a $\mathbb Z_2$ lattice-gauge Ising model, not the ordinary nearest-neighbor spin model. Keep the degrees of freedom (links versus sites), gauge redundancy, and plaquette interactions distinct in any duality computation.

## zhang2007 — first claimed exact solution

The claim depends explicitly on two conjectures: an extra rotation in a curled-up fourth dimension and temperature-dependent eigenvector weights. For the isotropic model it claims

$$\sinh(2K_c)\sinh(6K_c)=1,$$
$$x_c=e^{-2K_c}=\frac{\sqrt5-1}{2},\qquad K_c=-\frac12\log x_c=0.24060591\ldots.$$

It also claims $(\alpha,\beta,\gamma,\delta,\eta,\nu)=(0,3/8,5/4,13/3,1/8,2/3)$. Eqs. (3.37)/(3.62) are fourfold logarithmic integrals containing $w_x,w_y,w_z$. Because those functions are not derived, the expression is not closed or directly predictive.

## wu2008comment — direct series falsification

Wu, McCoy, Fisher, and Chayes identify two decisive tests:

1. Zhang fits 11 unknown coefficients in $w_y=w_z$ to 11 already-known HT coefficients, then sets $w_y=w_z=0$ for finite temperature. The second choice no longer reproduces the exact HT series.
2. With $x=e^{-2K}$, the claimed magnetization gives
   $$1-6x^8-12x^{10}-18x^{12}-\cdots,$$
   while the exact simple-cubic LT expansion begins
   $$1-2x^6-12x^{10}+14x^{12}-\cdots.$$

The proposed formula also reduces to exponent $3/8$ in 2D, contradicting Yang's exact $1/8$.

## zhang2008response — response to the series tests

The response describes Conjecture 1 as smoothing knot crossings and the weights of Conjecture 2 as a topological phase. It explicitly says the conjectures are not rigorous. It argues that standard HT/LT expansions need not govern finite temperatures, but supplies no independent next coefficient that would distinguish prediction from fit.

## wu2008rejoinder — fit versus prediction

The rejoinder emphasizes that choosing 11 free series coefficients to match 11 known terms is interpolation, not a prediction, and that the actual conjectured weights used afterward fail both rigorous HT and LT expansions. Bibliographic correction: `0812.0264` is unrelated cosmology; the rejoinder is arXiv:0812.2837.

## perk2012 — unknown-weight objection

Perk states that Zhang's Eq. (49) is an integral transform of unknown functions $w_x,w_y,w_z$. Fitting Eq. (A2) reproduces only supplied data; the alternative $(w_x,w_y,w_z)=(1,0,0)$ has a wrong first nontrivial HT term. A coding interface containing arbitrary weight functions is therefore not an exact solution until an independent rule determines them.

## zhang2021 — second claimed exact solution

The claimed critical condition is

$$\sinh(2H)\sinh(2H_1+2H_2)=1.\tag{32}$$

For $H_1=H_2=H$ it gives $H_c=0.30468893$, already far from the numerical benchmark. The paper's screw-boundary HT series is

$$\lambda_\infty=2\cosh^3H\left(1-3k^4-62k^6-1036k^8-20838k^{10}-\cdots\right),\tag{44}$$

whereas it prints the periodic result

$$\lambda^p_\infty=2\cosh^3H\left(1+3k^4+22k^6+192k^8+2046k^{10}+\cdots\right).\tag{45}$$

A bulk free energy cannot retain this boundary-condition discrepancy in the thermodynamic limit.

## perk2022comment — boundary, operator, and rotation tests

Peierls-Bogolyubov gives

$$|F[H_2]-F[H_1]|\le\|H_2-H_1\|.$$

Moving $O(ln)$ boundary bonds and dividing by $lmn$ yields

$$|f_{\rm screw}-f_{\rm periodic}|\le\frac{2|J_1|}{m}\to0.$$

The claimed equality $A_{p,1}=\bar A_{p,1}$ is invalid: commuting operators can share an eigenvalue set and eigenvectors while assigning those eigenvalues differently. Explicitly, $\sigma_p^z\sigma_{p+m}^z$ is not equal to $\sigma_p^z\sigma_{p+1}^x\cdots\sigma_{p+m-1}^x\sigma_{p+m}^z$.

Finally, permuting $(J_1,J_2,J)=(0.5,1.0,1.5)$ in Eq. (32) gives three alleged values $\beta_c=0.3503982204$, $0.3046889317$, and $0.2937911957$. Rotational invariance requires one value. This is a compact executable falsification test.

## zhang2023reply — reply and appended response

The reply introduces primed spin variables to defend the disputed sum and accepts that bulk physics is boundary-condition independent. It claims $\alpha=0$, $\nu=2/3$, and $\mu=4/3$. Perk's appended response points out that the primed operator is not the original physical Ising spin, so substituting it does not establish the required operator equality; accepting boundary independence also removes the proposed explanation for Eqs. (44)-(45).
