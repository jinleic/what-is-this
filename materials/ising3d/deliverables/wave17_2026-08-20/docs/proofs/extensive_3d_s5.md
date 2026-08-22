# The `s <= 5` extensive-charge class in the `Z^3` radius-2 box is decided

Artifacts: `experiments/e135_extensive_3d_s5.py`,
`tests/test_extensive_3d_s5.py`, `results/integrability/extensive_3d_s5.json`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e135_extensive_3d_s5.py
PYTHONPATH=src .venv/bin/python tests/test_extensive_3d_s5.py
```

**Status convention.** `[THEOREM]` and `[LEMMA]` are proved here; `[COMPUTATION]` is the exact
finite certificate; `[UNRESOLVED]` is what is deliberately not decided.

---

## 1. Setting

`H = sum_x X_x + sum_(<xy>) Z_x Z_y` on `Z^3`, `D = (1/2)[H, .]`,
`pi` the translation-orbit sum. The **support-size stratification** was proved in
`proofs/extensive_3d_r2.md`: a bond commutator
`(1/2)[Z_u Z_v, X^a Z^b] = -[a_u xor a_v = 1] X^a Z^(b+e_u+e_v)` is nonzero only when exactly one
endpoint already carries an `X`, so the support of a density word grows by at most one site and
the `s`-class system is closed. The class here is

```
C(3,2,5): rational span of ordered-Pauli words supported in B_2 = {0,1,2}^3
          with word support size <= 5.
```

The trivial classes are the identity density and the finite lattice divergences
`sum_i (g_i - tau_(e_i) g_i)` (Lemmas A/B of `proofs/extensive_3d_r2.md`), so the
nontrivial quotient is

```
Q(C) = dim ker(pi D |_C) / (ker(pi |_C) + Q I) = rank S_C - rank M_C - 1
```

with `S_C = pi |_C` and `M_C = pi D |_C`.

## 2. The certificate

**[COMPUTATION — sizes]** `|C(3,2,5)| = 21121156` columns by the exact closed form
(the wave-13 stock: `4^N - 4^(N-1) - ...`), `rank S_C = 14757412`, and the
Hamiltonian density `h = X_0 + sum_i Z_0 Z_(e_i)` is inside the class, so the
analytic upper bound `rank_Q M_C <= rank S_C - 2` holds (the kernel contains `I` and `h`,
linearly independent in `C / ker pi`).

**[COMPUTATION — two-prime exact rank]** Column-oriented elimination over `F_p` with compact
CSR pivot storage gives, at `p = 2147483647` and `p = 2147483629`,

```
rank_Fp M_C = 14757410    (= rank S_C - 2)
```

and since modular rank is an upper bound for rational rank (`rank_Fp <= rank_Q`), the
sandwich forces

```
[THEOREM]   rank_Q M_C = 14757410,     Q(C(3,2,5)) = 1.
```

The two primes agree exactly on the lead stream (identical sha256 of the field-independent
lead stream is a strong cross-check), `p1` took `1446.7 s` process time, `p2`
`1460.16 s`, peak RSS `2292.7 MiB` (under the 6 GiB cap).

**[COMPUTATION — greedy gap, recorded honestly]** The certified greedy *lower-bound* certificate
stalls at `744881 / 822332 = 0.905815` of the exact rank on the calibration class
`C(3,2,4)`. The exact column echelon, not the greedy one, is the certificate; the
greedy value is recorded as a `[COMPUTATION]` gap, not a bound that was used.

**[COMPUTATION — orbit-collision scan]** The elimination's direct triangular-minor claim was
independently checked on the calibration class: `1,503,765` images, `26,471,880` terms,
`ZERO` collisions, `ZERO` zero-sum orbits — so the recorded direct certificates are not
collision-collapsed.

## 3. The regression that guards it

The certified `s <= 4` values are reproduced by the new pipeline BEFORE the `s <= 5` run:
`columns = 1503766`, `rank S = 822334`, `rank M = 822332` at both primes, quotient `1`.
This is the same class (with the sandwich) whose committed 5 GiB wall was recorded in the
wave-13 note as `s <= 5` — the current `s <= 5` result is therefore a strict advancement
one class up the stratification, not a resizing.

## 4. What this does and does not decide

**[THEOREM — understand these traps.]** The artifact's checks include:
`D(X_0)` has exactly six bond terms, all `-1`; `D(Z_0) = + X_0 Z_0` (sign
convention fixed); `pi D h = 0` on orbit sums; the closed-form `rank_S` matches the
injectivity claim; and the `[UNRESOLVED]` that the FULL `Z^3` radius-2 box is **still not
decided** — concretely, any density whose words need `>= 6` non-identity sites stays open.

**[UNRESOLVED]** `C(3,2,6)` is, by the same closed form, `>= 4^N ...` (the full box
inventory) and was not attempted; the `s >= 6` classes and, a fortiori, the full
radius-2 box, are untouched. No all-size charge theorem follows from this finite class
statement.

## 5. Reading the literal table

| class | columns | rank S | rank_Fp M at both primes | analytic upper bound | certified quotient |
|---|---:|---:|---:|---:|---:|
| `s<=4` (regression, both primes) | 1 503 766 | 822 334 | 822 332 | 822 332 | 1 |
| **`s<=5` (headline)** | **21 121 156** | **14 757 412** | **14 757 410** | **14 757 410** | **1** |

The `s<=5` row is the certification. No row above changes the wave-13 statement: the
FULL `Z^3` `R=2` box is not decided.

## 6. Scope, resources, and provenance

The only gate in the pipeline is `time.process_time()`; RSS was measured with `resource`.
The producer's actual payload contains the full lead streams (not the launch-time validation
placeholder that was briefly mislabelled in a scratch file); the final artifact is the one
whose schema the test gates pass and whose checks all pass (6/6).
