# The duality pairing of the ladder sector-saturation route: the pairing matrix and its kernel

Artifacts: `experiments/e132_sector_saturation_pairing.py`,
`tests/test_sector_saturation_pairing.py`,
`results/algebra_growth/sector_saturation_pairing.json`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e132_sector_saturation_pairing.py
PYTHONPATH=src .venv/bin/python tests/test_sector_saturation_pairing.py
```

**Status convention.** `[THEOREM]` is proved here and machine-certified over `Q`;
`[COMPUTATION]` is the exact finite calculation recorded in the artifact;
`[COND]` is conditional; `[UNRESOLVED]` is what is deliberately not decided.

---

## 1. Setting

Same base as `proofs/ladder_alll_proof.md` and the wave-14 sector-saturation
fronts.  Let `L >= 3`, let the `2xL` ladder have sites `2i` (top of
rung `i`) and `2i+1` (bottom), and let the generators

```
A_L = sum_v X_v,   B_L = sum_{(u,v) in E} Z_u Z_v,   psi = |+>^(2L)
```

act on the configuration space `(C^2)^(2L)` in the `X`-eigenbasis
(`|x>` with `X_v|x> = (-1)^{x_v}|x>`, occupation bit `x_v`).  The
leg-swap `tau` and rung-reversal `rho` both fix `A_L`, `B_L`, `psi`;
every vector of the stabiliser module `U(A,B) psi = span{W(A,B) psi}`
is even-particle and `<tau,rho>`-invariant, so

```
U(A,B) psi  ⊆  K_L = {even-particle, <tau,rho>-invariant states},
dim K_L = 2^(2L-3) + 3*2^(L-2) = (2^(2L-1) + 3*2^L)/4.
```

This note works with the **orbit-sum basis** `{e_o = sum_{x in o} |x>}`
indexed by the `<tau,rho>`-orbits `o` of even configurations, `o` its
least configuration; `|o| in {1,2,4}`.

**The duality pairing.**  For every even configuration `x`, the delta-functional
`delta_x(v) = v_x` on `K_L` depends on `x` only through its orbit, so the
*distinct* delta-functionals of the sector are indexed by the orbits.  The
**pairing matrix** of the front is

```
P_L[o, j] = <delta_o, v_j>_conf = |o| * v_j[rep(o)],
```

against the explicit stabiliser-module columns `{v_j}`, i.e. the closure of
`psi` under the two symmetric operators `A_L = 2L - 2N` (diagonal in the
basis) and `B_L` (integer in the orbit-sum basis).  `P_L` is a `K_L x m`
integer matrix.  The kernel of `P_L` is the **dual annihilator**

```
W_L = {c in Q^{K_L} : sum_o |o| c_o v_o = 0 for all v in U(A,B) psi}
    = (U(A,B) psi)^perp   (orthogonal complement inside K_L in <.,.>_conf).
```

**[COMPUTATION]** Ranking `P_L` with rows only at distinct configurations
(`2^(2L-1)` rows) versus rows at the `K_L` orbit representatives gives the
same rank for `L = 3..6` — an invariant vector is constant on orbits, so the
surplus delta-rows are duplicates (`C_distinct_delta_rank_*`).

---

## 2. [THEOREM] Identity of the kernel with the certified W blocks

**Theorem (annihilator recursion).**  `W = (U(A,B) psi)^perp` is A- and
B-invariant and vacuum-orthogonal; conversely every A-,B-invariant subspace
of `K_L` orthogonal to `psi` lies in `(U(A,B) psi)^perp`.  Indeed `A_L,
B_L` are symmetric for `<.,.>_conf` and preserve `K_L`, so for `w in
(U psi)^perp` and any word `W(A,B)`, `<A w, W psi> = <w, A W psi> = 0`
(and likewise `B`); and a superspace symmetric argument gives the converse.
Hence `W` is exactly the largest A,B-invariant subspace of `K_L` orthogonal
to `psi`, i.e. **the pairing kernel is exactly the annihilator of the Krylov
span — the object whose dimensions are the wave-14 `W` blocks.**

**[COMPUTATION + THEOREM-CERTIFIED]  Exact cyclic dimensions and kernels.**

A two-prime modular closure over `F_p`, `p1 = 2147483647`,
`p2 = 2147483629`, gives `rank_{F_p} P_L = 14, 42, 142, 494`
(`L = 3..6`), identical at both primes.  Weighting the pivot rows by `|o|
mod p` and completing the free coordinates produces integer lifts, each required to
have all |coefficients| `<= 128`; these are then validated **over Q**: each
lift is vacuum-orthogonal, homogeneous in one particle sector, linearly
independent over Q, and the whole family is A- and B-invariant.  By the
theorem, `W_L` (the family) is contained in `(U psi)^perp`, so
`dim_Q U psi <= K_L - dim W_L`; the two-prime modular rank is the reverse
bound `dim_Q U psi >= rank_{F_p} = K_L - dim W_L`, hence both are **equalities over Q**:

```
L        dim K_L     dim(U(A,B) psi) = K_L - dim W_L     W sectors              W
3        14         14/14                                    {}                   0
4        44         42/44                                    {4:2}                2
5        152        142/152                                  {4:5, 6:5}         10
6        560        494/560                                  {4:15, 6:36, 8:15}  66
```

These are **exactly** the certified one-rung-TENSOR-front numbers
(`results/algebra_growth/sector_saturation_tensor.json`): cyclic dims
`14/14, 42/44, 142/152, 494/560` and the 4B-invariant annihilator
`W` dims `0, 2, 10, 66`.  The pairing route therefore lands on the same
`W` — the object the tensor front certified to obstruct cyclicity — by the
dual characterisation rather than by tensor-cut annihilator lifting.

---

## 3. What the pairing picture shows

**(a) Why the duality pairing is NOT cyclic.**  `SectorSaturationCyclicity(L)`
(the missing lemma named in `proofs/ladder_alll_proof.md`, viz.
`U(A,B) psi = K_L`) is equivalent to `W_L = 0`, i.e. to **non-degeneracy**
of the pairing.  The certified non-degenerate kernel `dim W_L = 0, 2, 10,
66` (L = 3..6) is the precise statement that the pairing degenerates from L = 4
on, and that the degeneracy family is the certified `W`.  The pairing-kernel
members are the orthogonal-complement witnesses: integral orbit-sum weight
families `{4:2}`, `{4:5, 6:5}`, `{4:15, 6:36, 8:15}`.

**(b) The evaluation bound survives unchanged.**  `dim g_L >= dim(U(A,B) psi)
>= K_L - dim W_L`:

```
L = 3..6:  dim g_L >= 14, 42, 142, 494,
```

the same certified evaluation dimensions as the tensor front (and, a fortiori,
`dim g_L >= 2^L` for these L).  The strongest `2^(2L-3)`-level target for
`dim g_L` is **not** certified by W (`494 < 512` at L = 6).

**(c) Regressions.**  Both mandatory regressions are embedded in the test and the
artifact: the cyclic dims `14, 42, 142, 494` (two primes) and the
distinct-delta-functional ranks at `2^(2L-1)` rows (`C_cyclic_dims_*`,
`C_distinct_delta_rank_*`, `C_exact_annihilator_*`), identical to the
sibling's certified values.

---

## 4. Scope and honest negatives

* All exact claims are for `L = 3, 4, 5, 6`; `L >= 7` is not
  computed here (the sibling front certifies `W_7 = 364`; its `W` and this
  note's kernel are the same object, but this front records `L = 3..6` only).
* The per-sector splits `{}`, `{4:2}`, `{4:5,6:5}`, `{4:15,6:36,8:15}`
  are certified per L here; the sibling did not tabulate per-sector splits, so
  the two fronts agree on the aggregate `W` dims and these per-sector families are
  additional (correct, consistent) detail.
* No statement about `dim g_L` itself is made besides `dim g_L >= K_L - dim W_L`;
  the asymptotics of `dim W_L / K_L` (0, 4.5%, 6.6%, 11.8%) stay
  `[UNRESOLVED]` exactly as in the tensor front.
* The lift bound `128` is an implementation contract of the certificate, not a claim
  about `W` (validation, not bound-search, is what certifies `W` exactly).
