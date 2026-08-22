# Theorem J-E′′: the demote fixed-point structure, and the catalogue-wide all-or-nothing corollary

**Status (2026-08-21):** the general theorem (fixed set = stable image
$I^\infty M$ = idempotent image; trichotomy of cases) is **proved**, and
**independently audited** (codex CLI read-only review of the premise link and
the stable-image lemma, 2026-08-21: all algebraic steps confirmed; three
presentational defects found and applied here: local-unit phrasing of the
local step, the nonzero-class fraction convention, and the descent/stop-rule
wording). The
*all-or-nothing* sharpening — demote fraction $\in\{0,1\}$ — is a
**machine-exact corollary for the 202 published catalogue parents**, not a
general law: a mixed parent is a coherent algebraic object, and whether any
BB parent realises the mixed case is left **open** (problem J-A below).
Machine confirmation on all 202 parents uses the ideal-power route
(EXP-052); an independent exact enumeration over the 196 parents with
$k_P \le 20$ (EXP-050: 2,693,100 nonzero quotient classes) agrees with the
algebraic predictions exactly.

## Statement

Let $R = \mathrm{GF}(2)[x,y]/(x^\ell-1,\,y^m-1)$ (finite commutative, hence
Artinian: $R=\prod_{\beta} R_\beta$, a product of local Artinian rings — note
that the coarse primary blocks $\mathrm{GF}(2)[x]/(p^a)\otimes
\mathrm{GF}(2)[y]/(q^b)$ are NOT in general local, e.g.
$\mathrm{GF}(4)\otimes\mathrm{GF}(4)\cong\mathrm{GF}(4)\times\mathrm{GF}(4)$
on the $6\times6$ lattice; the proof works at local granularity, but the
computations never need the decomposition). For a CSS BB parent write

$$M = \ker H_X / S_Z \quad (\dim_{\mathrm{GF}(2)} M = k_P),\qquad
I = L_{\rm pre} = \{\lambda \in R : \lambda A = \lambda B = 0\} \subset R.$$

By Theorem J-E the single-row syzygy class of $y\in M$ demotes iff
$y\notin Iy$.

**Lemma (fixed set).** $S := \{y\in M : y\in Iy\} = e_I\,M = I^\infty M$,
where $e_I = \sum_{\beta:\,I_\beta = R_\beta} e_\beta$ is the idempotent of
$I$'s full-factor support.

**Proof.** Finite commutative $R$ decomposes as $\prod R_\beta$; every ideal
and every finitely generated module splits componentwise: $I=\prod_\beta
I_\beta$ (for $a_\beta\in I_\beta$ pick $j\in I$ with $j_\beta=a_\beta$; then
$e_\beta j\in I$ has only that component) and $M=\bigoplus_\beta e_\beta M$
(a quotient module decomposes exactly like a submodule). In a local Artinian
ring $R_\beta$ with maximal ideal $\mathfrak m_\beta$, if
$I_\beta \subsetneq R_\beta$ then $I_\beta\subseteq\mathfrak m_\beta$, and
$y_\beta = i\,y_\beta$ for $i\in I_\beta$ gives $(1-i)y_\beta = 0$ **in
$R_\beta$**, where $1-i$ is a *local* unit, hence $y_\beta=0$ (this is the
only "Nakayama" needed). If $I_\beta = R_\beta$ the constraint is trivial.
Since $I$ contains $e_\beta$ for every full factor
($I_\beta = e_\beta I = R_\beta = e_\beta R \ni e_\beta$),
$y\in Iy \iff y = e_I y$. Hence $S = e_I M$. On the other hand
$e_I M = \prod_{\beta\ \mathrm{full}} M_\beta = I^\infty M$: the descending
chain $F_{r+1}=IF_r$, $F_0=M$, stabilizes at
$\prod_{\mathrm{full}}M_\beta \times \prod_{\mathrm{proper}}(I_\beta^\infty
M_\beta=0)$, since each proper $I_\beta\subseteq\mathfrak m_\beta$ is
nilpotent and there are finitely many factors. $\square$

**Theorem J-E′′ (fixed-point trichotomy).** Assume $k_P\neq 0$ (i.e.
$M\neq 0$; at $M=0$ cases (i) and (ii) coincide vacuously and the demote-set
is empty). Then the demote-set is $D = M\setminus e_I M$, and exactly one of:

1. **immune**: $e_I M = M \iff IM = M \iff 1\in I+\mathrm{Ann}_R(M)$
   (equivalent to Theorem J-E);
2. **demote-full**: $e_I M = 0 \iff I$ is nilpotent on $M$;
   every nonzero quotient class demotes;
3. **mixed**: $0\neq e_I M \neq M$: exactly the classes supported on $e_I M$
   fail to demote; the demote fraction **among nonzero quotient classes** is
   $1-\frac{2^{\dim e_I M}-1}{2^{k_P}-1}$ (among all elements of $M$ it would
   be $1-2^{\dim S}/2^{k_P}$; enumeration counts nonzero classes).

The chain $F_{r+1}=IF_r$ in EXP-047's quotient coordinates computes $\dim S$
by pure GF(2) linear algebra (each transition is one row-space computation
under the shift action). Containment $F_{r+1}\subseteq F_r$ is automatic
(ideal powers nest), so rank equality already certifies stabilization
mathematically; EXP-052 additionally *asserts* the membership containment
$F_{r+1}\subseteq F_r$ at every step as an implementation tripwire (it would
catch a wrongly coded action that silently violates the monotone law).

## Machine evidence

| scope | parents | method | immune | demote-full | mixed |
|---|---|---|---|---|---|
| catalogue | 202/202 | ideal-power chain (EXP-052) | 10 | 192 | **0** |
| $k_P\le 20$ | 196/202 | exact enumeration (EXP-050) | 10 | 186 | **0** |

Cross-validation: the algebraic immune set coincides with EXP-047's
`immunity_exact` list (10/10); the predicted integer count
$2^{k_P}-2^{\dim S}$ equals the enumerated `classes_demote` on all 196
enumerated parents, exactly (integer equality, not float comparison).

**The six previously unreachable parents are decided** (enumeration needs
$2^{24}$–$2^{60}$ classes; algebra needs up to three chain steps):
`phase2_90` ($k_P{=}24$), `30_6_0013` (24), `30_6_0039` (40),
`30_6_0007/0008/0012` (60): **all demote-full** (`30_6_0007/0008/0012` have
chain $60\to 20\to 0$; the others descend in one step).

Structural sharpening, exact on the catalogue (corollary of EXP-052, not a
general theorem):

- every demoting catalogue parent has $I^2M = 0$ (nilpotency index $\le 2$),
  and 165 of the 192 have $IM=0$ already;
- chain lengths are $\le 3$ everywhere: shapes observed are
  $(k_P,k_P)$-immune $\times 10$, $(k_P,0)$ $\times 165$, and
  $(k_P, u, 0)$ with $u\in\{4,8,20\}$ $\times 27$;
- the flagship `12_6_0193` has chain $(12,0)$: its 4095/4095 enumeration is
  now one rank test.

## Open problem J-A (the only remaining layer)

Whether a mixed BB parent exists at all. Mixed requires $M$-support on both
an $I$-full factor (a local factor where the pair $(A,B)$ vanishes
identically) and a degenerate-active factor (nonzero pair, deficient ranks).
One *heuristic* (unproved): for weight-$\le3$ pairs the vanishing condition
may saturate the polynomial — three monomials forming one small quotient-ring
relation — which would make the mixed shape expensive; this is a guess, not a
lemma. Constructing one mixed parent, or proving its impossibility for the
trinomial family, settles the last layer.

## Consequences for the no-go paper

- Demotion immunity on the catalogue is decided by $I^\infty M = M$ (one
  chain), and the demote-fraction is $\in\{0,1\}$ **for every catalogue
  parent** (machine-exact corollary) — the universal statement stands or
  falls with J-A.
- The X-side classification of the paper (192 demotion-realized / 8 certified
  X-monotone / 2 X-U) stands, with the collapse mechanism now certified
  without enumeration anywhere; the collapse *fractions* are catalogue facts.
