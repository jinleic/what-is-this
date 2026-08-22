# Exact gauge-tree branch map for the scalar Kac--Ward family

Artifacts: `experiments/e50_kw_branches.py`, `results/kac_ward/branches.json`, and `tests/test_kw_branches.py`.

## 1. Exact branch reconstruction

**[COMPUTATION]** The experiment reconstructs the wave-4/wave-6 tree in the fixed order

```text
+x->+y, +x->-y, +x->+z, +x->-z, +y->-x
```

with `0=vanishes` and `1=normalised to one`. Lexicographic binary enumeration gives exactly `2^5=32` labels from `00000` through `11111`. For every label, fresh substitution into the ungauged integer equations reproduces the stored zero pattern, nonzero pattern, active-variable count, and globally-linear-variable count.

**[UNRESOLVED]** These 32 labels are the strata induced by this selected five-edge direction-state tree. They are not a classification of all support graphs among the original 30 weights, so no claim about an omitted support chart follows from this inventory.

## 2. Seventeen exact characteristic-zero emptiness certificates

**[COMPUTATION]** The constant search uses the `k=4,6,8` identities on the ten oriented boxes

```text
3x3x2, 2x2x3, 3x2x3, 2x2x2, 3x2x2,
2x3x2, 4x2x2, 2x4x2, 2x2x4, 2x3x3.
```

For each branch it expands all 30 specialised polynomials over `ZZ`, forms the exact coefficient matrix of every nonconstant monomial, and computes its nullspace over `Q`. Thus failure of this search is complete only for constant-coefficient linear combinations of these 30 constraints; it does not exclude polynomial-multiplier Nullstellensatz certificates.

**[THEOREM]** On the twelve branches

```text
00000, 00001, 00010, 00011, 00100, 00101,
00110, 00111, 10000, 10010, 10100, 10110
```

direct expansion gives the same integer identity

```text
2 F_(3x3x2,4) - 3 F_(2x2x3,4) = 56.
```

At a common zero its left side would vanish, whereas `56 != 0` in `Q`; every listed branch is therefore `EMPTY_OVER_Q`. This extends the wave-6 certificate from one branch to twelve branches.

**[THEOREM]** On the five branches

```text
01000, 01001, 10001, 11000, 11001
```

direct expansion gives the independent oriented-box identity

```text
F_(3x3x2,4) - 2 F_(3x2x2,4) = -16.
```

The same common-zero contradiction proves each listed branch `EMPTY_OVER_Q`.

**[THEOREM]** Exactly 17 of the 32 named branches are therefore empty over `Q` by explicit primitive integer constant certificates. The result JSON stores the full 30-entry coefficient vector, nonzero terms, and constant separately for every branch rather than inferring one branch from another by symmetry.

## 3. Exact per-branch ledger

| branch | status | strongest exact evidence |
|---|---|---|
| `00000` | **[THEOREM] EMPTY_OVER_Q** | `2 F_(3x3x2,4)-3 F_(2x2x3,4)=56` |
| `00001` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `00010` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `00011` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `00100` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `00101` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `00110` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `00111` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `01000` | **[THEOREM] EMPTY_OVER_Q** | `F_(3x3x2,4)-2 F_(3x2x2,4)=-16` |
| `01001` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `01010` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `01011` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `01100` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `01101` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `01110` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `01111` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `10000` | **[THEOREM] EMPTY_OVER_Q** | `2 F_(3x3x2,4)-3 F_(2x2x3,4)=56` |
| `10001` | **[THEOREM] EMPTY_OVER_Q** | `F_(3x3x2,4)-2 F_(3x2x2,4)=-16` |
| `10010` | **[THEOREM] EMPTY_OVER_Q** | `2 F_(3x3x2,4)-3 F_(2x2x3,4)=56` |
| `10011` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `10100` | **[THEOREM] EMPTY_OVER_Q** | `2 F_(3x3x2,4)-3 F_(2x2x3,4)=56` |
| `10101` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `10110` | **[THEOREM] EMPTY_OVER_Q** | `2 F_(3x3x2,4)-3 F_(2x2x3,4)=56` |
| `10111` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `11000` | **[THEOREM] EMPTY_OVER_Q** | `F_(3x3x2,4)-2 F_(3x2x2,4)=-16` |
| `11001` | **[THEOREM] EMPTY_OVER_Q** | same freshly specialised identity |
| `11010` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `11011` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `11100` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `11101` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `11110` | **[UNRESOLVED]** | no constant certificate; exact-`Q` Gröbner timeout at 60 s |
| `11111` | **[UNRESOLVED]** | exact `F_3,F_5,F_7,F_11` points plus the p-adic evidence below; exact-`Q` Gröbner timeout at 60 s |

**[COMPUTATION]** Each of the 15 unresolved systems received a separate exact characteristic-zero Gröbner attempt on the branch-specialised six-equation construction ideal. Exact row reduction first removed rationally dependent input polynomials; variables linear in every retained equation formed a lexicographic elimination block, followed by a grevlex block. Every worker reached its hard 60-second bound. A timeout is recorded as `UNRESOLVED`, never as evidence for a rational point.

## 4. P-adic evidence on branch `11111`

**[COMPUTATION]** The only branch with stored finite-field witnesses is `11111`; the predecessor artifact supplies exact construction-system points over `F_3`, `F_5`, `F_7`, and `F_11`. All four were freshly reverified. The successful selected construction lifts and exact continued-fraction reconstruction results are:

| prime | exact lift | height bound `H` with `2H^2<M` | reconstructed coordinates | failed coordinates |
|---:|---|---:|---:|---|
| 3 | **[COMPUTATION]** singular selected path through `3^126` | `809327037926455399226508505424` | 8/25 | `u_mx_mx,u_mx_pz,u_mx_mz,u_py_px,u_py_pz,u_py_mz,u_my_px,u_my_my,u_my_pz,u_my_mz,u_pz_px,u_pz_mx,u_pz_my,u_mz_px,u_mz_py,u_mz_my,u_mz_mz` |
| 5 | **[COMPUTATION]** full-row-rank Newton lift through `5^86` | `803887338846092830409411774445` | 23/25 | `u_mx_mx,u_mx_py` |
| 7 | **[THEOREM]** this witness has no lift modulo `49` | — | — | exact left-kernel obstruction |
| 11 | **[COMPUTATION]** full-row-rank Newton lift through `11^60` | `12338590671981506729216193325676` | 23/25 | `u_mx_my,u_mx_pz` |

**[LEMMA]** For a residue `r mod M`, the extended-Euclidean reconstruction used here is complete for reduced `a/b` satisfying `|a|<=H`, `1<=b<=H`, `gcd(b,M)=1`, and `2H^2<M`. Consequently every recorded coordinate failure excludes such a fraction for that selected p-adic residue. It does not exclude a different local section or finite-field point.

**[COMPUTATION]** Every coordinate of every successful lift was attempted. No lift reconstructed all 25 coordinates, so no rational vector was available for exact all-coordinate substitution and the stored rational-candidate list is empty.

**[COMPUTATION]** Fresh exact evaluation on the `3x3x3` `k=4,6,8` identities falsifies all three selected construction lifts at finite p-adic precision. Their capped valuation triples are

```text
p=3:  (126,1,1)
p=5:  (86,0,3)
p=11: (60,0,0).
```

A valuation equal to the lift exponent means that component is zero modulo the full modulus; a smaller valuation is an exact nonzero residual. These falsifications concern the selected construction lifts only.

**[THEOREM]** For the stored `F_7` construction point, `lambda=(4,0,0,1,0,0)` annihilates the six-by-25 Jacobian modulo 7 while pairing to `5 mod 7` with the linearised right-hand side. Therefore this point has no construction-system lift modulo 49; other `F_7` points remain unresolved.

**[THEOREM]** The stored `F_3` point passes the `3x3x3` identities modulo 3, but it has no simultaneous lift modulo 9 after those three identities are adjoined to the six construction equations. For the nine-equation system,

```text
F(x) mod 9 = (6,0,3,6,0,6,0,3,3),
-F(x)/3 mod 3 = (1,0,2,1,0,1,0,2,2),
lambda = (0,0,0,0,0,0,0,1,0).
```

Exact recomputation gives `lambda J=0 mod 3` and `lambda*(-F/3)=2 mod 3`. This excludes this one combined residue point, not every `F_3` point.

## 5. Scope and conclusion

**[THEOREM]** The exact result is a strict partial no-go: 17 named gauge-tree branches are `EMPTY_OVER_Q`.

**[UNRESOLVED]** The other 15 named branches have neither an emptiness certificate nor a verified rational point. In particular, finite-field solutions, selected p-adic lifts, height bounds, exact holdout failures, and bounded Gröbner timeouts do not justify a full-family theorem.

**[UNRESOLVED]** Because not all 32 branches are empty, no global claim that the full 25-coordinate `k<=8` family admits no rational realization is made.

## 6. Reproduction

**[COMPUTATION]** Run only the dedicated experiment and independent standalone test:

```text
.venv/bin/python experiments/e50_kw_branches.py
.venv/bin/python tests/test_kw_branches.py
```

**[COMPUTATION]** Both commands print a final `PASS`; the test independently regenerates every stored empty-branch identity from fresh lattice constants, checks every high-precision construction residue and rational-reconstruction outcome, and re-derives the `F_7 mod 49` and combined `F_3 mod 9` left-kernel certificates.
