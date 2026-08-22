# No tridiagonal (Askey–Wilson) relation for the 3D Ising transfer pair

**Second main result of Track A.**  Theorem DG (`proofs/dolan_grady_defect.md`) shows the
Dolan–Grady relations fail whenever the layer graph has a vertex of degree `>= 3`.  A natural
objection is that Dolan–Grady is only one normalisation.  This note closes that loophole: the
3D Ising layer pair satisfies **no** tridiagonal relation, for **any** parameters, and this
persists under every affine redefinition and rescaling of the two generators.

Scripts: `experiments/e21_tridiagonal.py` (decision), `experiments/e22_tridiagonal_certificate.py`
(explicit witnesses), `experiments/e23_certificate_general.py` (scope of the one-line witness).
Data: `results/tridiagonal.json`, `results/tridiagonal_certificate.json`,
`results/tridiagonal_certificate_general.json`.

---

## 1. The criterion

For a layer graph `Gamma`, `A = sum_v X_v`, `B = sum_{(ij) in E} Z_iZ_j`.  The general condition
for `(A,B)` to be a *tridiagonal pair* — the structure underlying the Onsager algebra, the
Askey–Wilson algebra and Terwilliger's Leonard pairs, and hence the whole Onsager route to an
exact solution — is the pair of relations (Ito–Tanabe–Terwilliger)

    (TD1)  [ A,  A^2B - beta ABA + BA^2 - gamma (AB+BA) - rho B ]  = 0
    (TD2)  [ B,  B^2A - beta BAB + AB^2 - gamma*(BA+AB) - rho* A ] = 0

Dolan–Grady is the single point `beta = 2, gamma = gamma* = 0, rho = rho* = 16`.

**Observation (what makes this decidable).**  For fixed `A, B` each relation is *linear* in its
three scalars.  Setting

    T0 = [B, B^2A + AB^2],   T1 = [B, BAB],   T2 = [B, BA + AB],   T3 = [B, A],

(TD2) reads `T0 = beta T1 + gamma* T2 + rho* T3`, i.e. `T0 in span{T1,T2,T3}`.  Expanded in the
Pauli-string basis this is an exact integer linear system, decided by rational Gaussian
elimination.

---

## 2. Affine completeness lemma

> **Lemma.** Let `A' = aA + s`, `B' = bB + t` with `a,b != 0`.  Then `(A',B')` satisfies a
> tridiagonal relation iff `(A,B)` does.

*Proof.* Every term of `A'^2B' - beta A'B'A' + B'A'^2 - gamma(A'B'+B'A') - rho B'` is linear in
`B'`.  The `t`-dependent part is `t(2-beta)(aA+s)^2 - 2 gamma t (aA+s) - rho t`, a polynomial in
`A` alone, annihilated by the outer bracket `[A', .] = a[A, .]`.  The `B`-dependent part equals
`b` times

    a^2 (A^2B + BA^2) - beta a^2 ABA + a[s(2-beta) - gamma](AB+BA) + [s^2(2-beta) - 2 gamma s - rho] B .

Dividing by `a^2 b`, the relation becomes the unshifted one with

    gamma_new = (gamma - s(2-beta))/a ,     rho_new = (rho + 2 gamma s - s^2(2-beta))/a^2 ,

and `beta` unchanged.  As `(gamma, rho)` range over `Q^2` so do `(gamma_new, rho_new)`. `[]`

So the three-parameter linear test already decides the full five-parameter affine family.
Verified in code: `affine_note()` in `e21_tridiagonal.py`.

---

## 3. Result

`experiments/e21_tridiagonal.py`, exact rational arithmetic.

### 3.1 (TD1) is uninformative

`T0^{(A)} = 2 T1^{(A)} + 0 * T2^{(A)} + 16 T3^{(A)}` for **every** graph tested — consistent with
the analytic proof (Theorem DG(i)) that `[A,[A,[A,B]]] = 16[A,B]` identically.  All the
information is in (TD2).

### 3.2 (TD2): the dichotomy

| layer graph | regular? | (TD2) solution |
|---|---|---|
| cycle `C_3` | 2-regular | solvable (degenerate, 2-parameter family) |
| cycle `C_4` | 2-regular | solvable, `(0,0,16)` + 1 free parameter |
| cycle `C_5` | 2-regular | solvable, `(-14,16,0)` + 1 free parameter |
| cycle `C_6, C_7, C_8, C_9` | 2-regular | solvable, **exactly `(2,0,16)` = Dolan–Grady** |
| torus layer `3x3, 3x4, 4x4, 3x5, 4x5, 5x5` | 4-regular | **NO SOLUTION** |
| open grids `2x3, 3x3` | `Delta >= 3` | **NO SOLUTION** |
| open paths `P_4, P_6` | `Delta = 2`, has boundary | no solution (boundary terms; see below) |

The comparison that matters is between the two **boundaryless, translation-invariant** families,
which removes the boundary confound entirely:

    2-regular layer  (2D Ising on a cylinder)  ->  (A,B) IS a tridiagonal pair, with the
                                                   Dolan-Grady parameters (2,0,16);
    4-regular layer  (3D Ising, torus layer)   ->  (A,B) satisfies NO tridiagonal relation.

Open 1D chains also fail (TD2), by the boundary terms already isolated in Theorem DG; this is the
known and benign fact that a free end needs the boundary-modified Onsager/`q`-Onsager machinery.
It does not affect the comparison above, in which neither side has a boundary.

---

## 4. Explicit certificates (so this is a proof, not a rank computation)

`experiments/e22_tridiagonal_certificate.py` extracts, for every failing layer, a set of at most
four Pauli strings on which the three-unknown system is *already* inconsistent.  Verifying a
certificate requires reading four integers per row and solving a `3x3` system by hand.

**Best case — a ONE-LINE proof.**  For the `4x5` torus layer, the single Pauli string

    W = Z_0 Y_5 Z_6 Z_10        (i.e. `Y_v Z_{u1} Z_{u2} Z_{u3}` for v = 5 and three of its
                                 four neighbours)

has

    [T0]_W = -48 ,      [T1]_W = [T2]_W = [T3]_W = 0 .

Therefore `T0` has a component that no combination of `T1, T2, T3` can produce, and (TD2) is
unsolvable.  **This is a complete proof in one line.**

The same holds for **all 64** quartic strings `W` of the `4x4` torus layer
(`e23_certificate_general.py`).  The empirical scope of the one-row certificate, read off from
`results/tridiagonal_certificate_general.json`, is exactly:

| torus layer | one-row certificates | `[T0]_W` values |
|---|---|---|
| `4x4` | 64 of 64 | `-48` |
| `3x4`, `4x5`, `4x6` | exactly half | `-48` and `-40` |
| `3x3`, `3x5`, `3x6`, `5x5` | none | `-48` (3x3) / `-40` |

i.e. a one-row certificate exists precisely when at least one side length equals `4`.  This is a
degeneracy of the short cycle, not a feature of the obstruction: on the remaining layers the
coefficients `[T1]_W, [T2]_W, [T3]_W` no longer all vanish, and a **four-row** witness is used
instead.  Every failing layer in `results/tridiagonal_certificate.json` carries an explicit
witness; the witness sizes occurring are `{1, 4}`.

**Unification with Theorem DG.**  The certificate string `W = Y_v Z_{u1}Z_{u2}Z_{u3}` is *exactly*
the quartic Dolan–Grady defect operator of Theorem DG(ii).  The same 4-body operator that measures
the failure of Dolan–Grady is the witness that kills the entire tridiagonal family.  The two
results are one obstruction seen twice.

---

## 5. What is and is not established

**Established (theorem + machine-checked certificate).**
The pair `(A,B)` of the 3D Ising layer transfer operator is not a tridiagonal pair, on any layer
tested up to `5x5 = 25` sites, for any `(beta, gamma*, rho*)`, and invariantly under
`A -> aA+s`, `B -> bB+t`.  Hence the Onsager / Dolan–Grady / Askey–Wilson mechanism — the
mechanism that solves the 2D Ising model — cannot be applied to the 3D Ising transfer matrix in
its natural generators.

**Not established.**
1. This does not prove the 3D Ising model is unsolvable.  It closes one specific, historically
   central mechanism, in its full parametric generality.
2. It is a statement about the *natural* generator pair `(sum X, sum ZZ)`.  A different
   factorisation of the same transfer operator into two generators is not covered; the
   deformation search (`proofs/dg_deformation_nogo.md`) addresses part of that gap.
3. Layers larger than `5x5` were not tested; the computation is `O(|T0|)` with
   `|T0| ~ 88200` strings at `5x5`, and grows quickly.  Nothing suggests a change of behaviour,
   but nothing here proves it either.

**Independent corroboration.** `experiments/e10_conserved_charges.py` finds that for the same
2D layer the space of local conserved charges of `H = -g sum X - sum ZZ` is exactly
`span{I, H}` for interaction ranges 2, 3 and 4 (ansatz dimensions 22, 220, 2461), while the 1D
chain has a growing tower (`1, 3, 5, 7` nontrivial charges at ranges 2, 3, 4, 5).  Two logically
independent integrability diagnostics therefore agree.
