# P-adic lifting of the full scalar Kac--Ward construction system

Artifacts: `experiments/e42_kw_hensel.py`, `results/kac_ward/hensel_lift.json`, and `tests/test_kw_hensel.py`.

## 1. Exact input and scope

**[COMPUTATION]** The experiment loads, rather than regenerates, the six stored integer polynomials in the 25 variables of the generic five-edge gauge chart from `results/kac_ward/full_family.json`. It reparses the expressions over `ZZ` and recomputes the canonical SymPy hash

```text
45ac8bdd54b7a309b197f1905d81c0ed255573e4dbc8d2d0b0465837b9c9cded
```

before doing any lift. The six equations are the `k=4,6,8` trace equations for the free `3x3x2` and `2x2x3` construction boxes.

There are six equations and 25 variables. Thus “nonsingular” below means full **row** rank six, not a 25-by-25 isolated point. At each full-row-rank witness the computation chooses the first six pivot columns

```text
u_px_px, u_mx_mx, u_mx_py, u_mx_my, u_mx_pz, u_mx_mz
```

as Newton variables. The remaining 19 coordinates are held exactly at their least nonnegative finite-field representatives. This defines a reproducible local p-adic section; it is not a parametrization or classification of the positive-dimensional construction variety.

## 2. Hensel and reconstruction certificates

**[LEMMA]** Suppose `F(x)=0 mod p^m` and a six-by-six Jacobian minor `A(x)` has determinant nonzero modulo `p`. For `m' <= 2m`, solve

```text
A(x) y = -F(x)/p^m  (mod p^(m'-m))
```

and replace the six pivot coordinates by `x+p^m y`. Taylor expansion shows that the new point solves `F=0 mod p^m'`; the terms quadratic in the correction are divisible by `p^(2m)`. Since the minor determinant is a p-adic unit, every modular linear solve is valid. The experiment uses `m'=min(2m,N)`, so precision doubles until the requested exponent is reached.

**[LEMMA]** For modulus `M`, set

```text
H = floor(sqrt((M-1)/2)).
```

Then `2 H^2 < M`. A reduced fraction `a/b` with `|a| <= H`, `1 <= b <= H`, and `gcd(b,M)=1` is unique if `a = r b mod M`: two such fractions would give a nonzero integer determinant of absolute value at most `2H^2` divisible by `M`. The extended-Euclidean rational-reconstruction algorithm is complete under this bound.

Consequently, a reconstruction failure for one coordinate certifies that no fraction within this numerator/denominator box represents that coordinate residue. If at least one of 25 coordinates fails, there is no rational 25-vector of coordinate height at most `H` in that **selected residue class**. It does not exclude rational points in other residue classes or other local sections.

## 3. Per-prime results

**[COMPUTATION]** All residues and ranks in this table are exact.

| prime | Jacobian row rank | lift result | final modulus | decimal digits below modulus | `H` | reconstructed coordinates |
|---:|---:|---|---|---:|---:|---:|
| 5 | 6 | Newton/Hensel lift through `5^86` | `1292469707114105741986576081359316958696581423282623291015625` | 60 | `803887338846092830409411774445` | 23/25 |
| 7 | 5 | no lift of this witness even to `7^2` | `49` attempted | — | — | — |
| 11 | 6 | Newton/Hensel lift through `11^58` | `2516377186292711566730985912068419625116019959228909823321881` | 60 | `1121690061089227884474199393243` | 22/25 |

For `p=5`, the reconstruction failures are `u_mx_mx` and `u_mx_py`. The other 23 coordinates reconstruct, but 19 of those are deliberately held integer coordinates; among the moving coordinates, some reconstructions have numerator and denominator near the height ceiling. The moving coordinate `u_mx_pz` happens to reconstruct as `-11`, but it does not reconstruct to the same value on the `p=11` section.

For `p=11`, the reconstruction failures are `u_px_px`, `u_mx_mx`, and `u_mx_my`. Again, all 19 held coordinates reconstruct trivially, while the successful moving-coordinate fractions are not cross-prime consistent with the `p=5` lift.

**[COMPUTATION]** The stored `F_7` point has rank five. Evaluating the six integer equations at its least nonnegative representatives modulo 49 gives

```text
F(x) mod 49 = (42, 7, 28, 42, 42, 21),
-b := -F(x)/7 mod 7 = (1, 6, 3, 1, 1, 4).
```

Here `b` is written in the convention used by the linearized equation `J t = -F(x)/7`. The exact left-kernel vector

```text
lambda = (4, 0, 0, 1, 0, 0)
```

satisfies `lambda J = 0 mod 7`, while its pairing with the displayed right-hand side is `5 mod 7`. Therefore the linearized correction equation is inconsistent.

**[THEOREM]** This particular stored `F_7` point has no lift to a solution modulo 49. This theorem is about one finite-field point only. Other `F_7` points, and hence the branch as a whole, remain undecided.

## 4. Cross-prime partial reconstruction

**[COMPUTATION]** Comparing the selected `p=5` and `p=11` sections coordinate by coordinate gives two consistent reconstructed values:

```text
u_pz_mx = 2,  u_pz_my = 2.
```

Both coordinates are in the deliberately held block at both primes, so this agreement is an artifact-compatible partial coincidence rather than evidence for a shared rational solution. Nineteen coordinates reconstruct at both primes but to different rationals; four coordinates are unavailable because reconstruction fails at one or both primes. No moving coordinate reconstructs consistently across the two primes.

**[NOT PROVED]** These two p-adic points need not be reductions of one characteristic-zero point. Cross-prime inconsistency of selected sections does not rule out another rational or algebraic point.

## 5. Zero-pattern information and the bounded Groebner retry

The predecessor artifact has no solution counts across all 32 patterns. Its sampled Jacobian ranks are ranks at random points, not evidence that those points solve the equations. The only stored solution-pattern statement available for `p=5,7,11` is that pattern `11111` is nonempty: each stored generic-chart witness has all five gauge-tree entries normalized to one. The `p=7` nonlift result above does not make that finite-field pattern empty.

Patterns `00000` and `00001` tie for the smallest active-variable count, 12. The experiment deterministically chooses lexicographically first `00000`, orders its four globally linear active variables in a lexicographic elimination block followed by the other eight in a grevlex block, and runs exact SymPy Groebner computations with a hard 2400-second limit each.

**[COMPUTATION]** Both bounded computations completed in under one second of worker time and returned the reduced basis `[1]`, over `F_101` and over `F_32003`. The stored canonical basis hash is

```text
6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b
```

for each prime.

The modular output led to a much shorter exact certificate. After setting all five selected tree weights to zero, let `F_(B,k)` denote the stored trace equation for box `B` and order `k`. Direct symbolic expansion gives

```text
2 F_(3x3x2,4) - 3 F_(2x2x3,4) = 56.
```

The standalone test re-derives this identity from `full_weight_finite_system(..., gauge_fix=False)` rather than reading it from the result JSON.

**[THEOREM]** Tree-zero pattern `00000` has no common solution over `Q`, and no solution over any field whose characteristic does not divide 56. Indeed, a common zero would make the left side zero while the right side is nonzero. This excludes exactly one of the 32 recorded tree patterns; it does not exclude pattern `00001`, any of the other 30 patterns, or support charts not represented by this selected tree stratification.

For clarity about specialization: a basis `[1]` modulo `p` proves emptiness only for that reduced affine system over the algebraic closure of `F_p` and excludes a p-integral lift on that chart. Modular emptiness alone does **not** imply emptiness over `Q`, because a rational point may have a denominator divisible by `p`. Conversely, characteristic-zero emptiness only forces `[1]` after reduction at primes not dividing a cleared-denominator certificate. Here the independent integer identity with constant 56, not the two modular bases by themselves, supplies the characteristic-zero theorem.

No fresh `p=13` or `p=17` search was run. Uniform random search for six equations has heuristic hit costs `13^6=4,826,809` and `17^6=24,137,569`, respectively, and the stored artifact supplies no triangular finite-field solver. A tiny unsuccessful random sample would have had no honest negative interpretation, so it was not performed.

## 6. Exact-candidate and epistemic status

**[COMPUTATION]** Neither full-rank lift reconstructs all 25 coordinates. Therefore no rational construction candidate was produced and the exact `3x3x3` holdout pipeline was not triggered. The experiment contains a guarded path that, if all coordinates reconstruct, first substitutes the candidate exactly into all six construction equations and then recomputes exact `k=4,6,8` targets on independent `3x3x3` and `4x2x3` boxes before allowing a `THEOREM-CANDIDATE` label.

**[THEOREM]** The exact constant certificate excludes branch `00000` over characteristic zero.

**[COMPUTATION]** The selected full-row-rank `p=5` and `p=11` points lift beyond 60 decimal digits, but each selected p-adic section fails all-coordinate rational reconstruction at the certified height shown above. The stored `p=7` point has the exact first-order nonlift certificate above.

**[UNRESOLVED]** None of these facts decides the remaining construction variety, the other zero patterns, or the full translation-invariant scalar Kac--Ward family. Reconstruction failure is evidence against low-height rational points only in two explicitly selected p-adic residue classes. The family remains **UNRESOLVED**.

## 7. Reproduction

```text
.venv/bin/python experiments/e42_kw_hensel.py
.venv/bin/python tests/test_kw_hensel.py
```

The experiment and standalone test both print a final `PASS`.
