# Characteristic-zero structure of the 2x4-layer 18-local-term DLA

## Object and statement

Let `g` be the dynamical Lie algebra generated over `Q` by the 18 local terms of
the open `2 x 4` grid: the eight one-site `X_i` and the ten bond terms
`Z_u Z_v` of `grid_graph(2,4)` (6 horizontal: 2 rows x 3; 4 vertical).
This is the direct analogue, one layer up, of Theorem C16's object (the 2x3
layer's 13 local terms, `experiments/e143_char0_quotient.py`).

*Hygiene.* This is **not** the two-generator algebra `<A,B>` on the same layer:
that object has dim 2952 over Q (Theorem OA-2x4, report sec. 3.8,
`proofs/char0_complete_2x4.md`), centre of dimension 1, and type
`D21 + 2B13 + A23 + A27 + 2A3`.  It is a *different* algebra and is fully
settled elsewhere; this note's disambiguation (H396 class) is recorded below.

**Theorem 1.**  `dim_Q g = 16256`; `g` is perfect with a nondegenerate Killing
form that is **diagonal** in the Pauli-string basis, of exact signature
`(8192, 8064)`; the natural module `V = Q^256` splits under `P = prod X` as
`V+ (+) V-` (128/128) with both summands `g`-modules; on each half there is an
exact nondegenerate **symmetric** rational form `F_+ = F_-` (det 1) such that

```
pi_+(g) = so(F_+, Q),        pi_-(g) = so(F_-, Q),       dim = 8128 each,
g = ker(pi_+)  (+)  ker(pi_-),      each factor an ideal of dim 8128
```

and consequently

```
g  =  so_128(Q)  (+)  so_128(Q)      (split D_64, dim 8128 each, total 16256).
```

The predicted `sp_128 (+) sp_128` (dim 8256 + 8256 = 16512) is **refuted**:
the invariant form is symmetric, not alternating, exactly as in the C16 control
`grid 2x2` (`so_8 (+) so_8`); the 2x3 layer is the alternating case
(`sp_32 (+) sp_32`).

**Corollary 2** (quotient dimensions).  The quotient dimensions of `g` by its
ideals are exactly `{0, 8128, 16256}`; in particular `sp_128` (8256), the
2952-dim OA-2x4 algebra, and any dimension outside this set are not quotients of
`g`.

## 1. Closure, Killing form, derived series

**[COMPUTATION]** The 18 generator strings close under the anticommuting-xor rule
(`v,w -> v xor w` when `<v,w> = 1`) to a set `S` of size `d = 16256`
(16,256 of the 65,536 Pauli strings; the two engines -- the independent BFS of
`e148` and the repository engine `ising.clifford.dla_pauli_closure` --
agree set-for-set; `results/dla_table.json` already records 16256 for
`grid_open 2x4`).

**Lemma 1 (diagonal Killing form).**  For the string basis `{Q_v : v in S}`,
`Tr(ad Q_v ad Q_w) = 0` for `v != w`: `ad Q_v` shifts string indices by
`v xor w`, and the trace vanishes unless the shift is trivial.  Hence the Killing
form is diagonal with
`kappa_vv = sum_{j in S} c(v, j) c(v, j xor v)`, `c(v, j) in {0, +-2}` the
Pauli structure constants.

**[COMPUTATION]** The vectorised row-wise evaluation gives `kappa_vv != 0` for
all `v` (min abs 32256), exact signature `(8192, 8064)`, and every element of
`S` is an anticommuting pair sum (`[g,g]` has full support = `g`, i.e. `g` is
perfect).  By Cartan's criterion `[EXTERNAL]` with nondegenerate diagonal Killing
form, `g` is semisimple; the radical is `{x : kappa(x,[g,g]) = 0} = {0}`.

## 2. Parity split and the half representations

**[LEMMA]** Every string of `S` has even Z-weight (generators do, and xor keeps
parity), so `P = prod_i X_i` commutes with `g`; `V = V+ (+) V-` with
`dim V+ = dim V- = 128` is a `g`-module decomposition.  The half-blocks are
built from the exact integer `pauli_matrix` and the exact integer pair-basis change
of basis; `block_float` agrees with `e143.parity_blocks` matrix-for-matrix on
every string of every control case (exhaustive).

## 3. The invariant form is symmetric (so, not sp)

**[COMPUTATION]** The linear system `X^T B + B X = 0` for the 18 generator
blocks of each half is set up as exact equations and reduced over `F_p` by sparse
elimination (per-generator equations kept separate; pivots normalized; kernel by
back-substitution).  Reproduces `e143.solve_invariant_form` on all controls.
For `2x4` on each half: `kernel dim = 1`, rational reconstruction gives an integer
matrix `F` with `det F = 1`, `max abs entry 1`, `F^T = F` (SYMMETRIC), and
exact verification `X^T F + F X = 0` for all 18 `X`.

**Lemma 2.**  Every element of `pi_+(g)` (resp. `pi_-(g)`) is `F`-skew:
the `F`-skew matrices form a Lie subalgebra and the 18 generators (hence the
whole algebra they generate) lie in it.  Therefore
`pi_+-(g) <= so(F, Q)`, `dim pi_+-(g) <= 128*127/2 = 8128`.

## 4. Faithfulness bound and the equality

**[THEOREM]** Distinct Pauli strings are linearly independent matrices over `Q`, so
the natural representation `rho = pi_+ (+) pi_- : g -> gl(V+ (+) V-)` is
injective and `dim rho(g) = dim g = 16256`.  Hence
`r_+ + r_- >= 16256` where `r_+- = dim pi_+-(g)`, and by Lemma 2
`r_+- <= 8128`, forcing `r_+ = r_- = 8128` and `rho(g) = so(F_+) (+) so(F_-)`
(same dimension as the container).  Since `rho` is injective,

```
g  =  ker(pi_+) (+) ker(pi_-)   ~=   so_128(Q) (+) so_128(Q),
ker(pi_+) ~= so(F_-),   ker(pi_-) ~= so(F_+),
```

with each factor and each `pi_+-` an isomorphism onto `so(F_+-)`.
The reconstructed integer form is a sum of 64 hyperbolic planes (128 nonzero
off-diagonal entries, max abs 1, signature (64,64), and diag(1,-1) ~= [[0,1],[1,0]]
over Q), so `so(F) = so_128(Q)` is the split `D_64`.  The Killing signature
`(8192, 8064) = (4096, 4032) + (4096, 4032)` is the split-`so(128)` pair
(D_64: (D+r)/2 x (D-r)/2 with D = 8128, r = 64), consistent with Theorem 1.

## 5. Half-image ranks: Walsh-domain certificate (corroboration)

The rank identity `r_+- = 8128` is, in addition, certified directly by the
Walsh-domain structure lemma:

**Lemma 3 (regular-abelian monomial family).** Each even-Z block is a signed
permutation in the pair basis.  The permutation parts act regularly as an abelian
group of order 128 on the `(d, m)` index set; after conjugation by the signature
Walsh matrix `W` of size 128 every block is a character-scaled matching:
`W^{-1} M W` has, in every row, exactly one nonzero entry, at column
`y xor d` with value `+-128 * chi_m(y)`.  The matchings for distinct `d` are
disjoint, and distinct characters are independent over any field of char != 2;
hence `r = #{distinct (d, m)}` over the basis strings.

**[COMPUTATION]** `read_dz` verifies the pure matching + character row-profile
for **all 16,256 strings** of both halves (0 failures), the count is
`r_+ = r_- = 8128`, and it agrees with the dense mod-p row rank (e143
`rank_mod_p_cert`, both primes) on every control case; per-d integer minors of
the conjugated stack are certified nonzero mod `p1 = 2147483647` and
`p2 = 2147483629` (127 matching offsets, sum 8128 rows per half per prime).

## 6. Verification table

| case             | d     | half | form      | r+ = r- | type            | sig           | claim                      |
|------------------|-------|------|-----------|---------|-----------------|---------------|----------------------------|
| chain n=2        | 6     | 2    | alternating| 3      | sp(2) = sl_2    | (4,2)         | sl_2 (+) sl_2              |
| chain n=3        | 15    | 4    | none      | 15      | sl(4) = so(6)   | (9,6)         | sl_4 (simple)              |
| grid 2x2         | 56    | 8    | symmetric| 28      | so(8) [D_4]     | (32,24)       | so_8 (+) so_8             |
| grid 2x3         | 1056  | 32   | alternating| 528   | sp(32) [C_16]   | (544,512)     | sp_32 (+) sp_32 (C16)     |
| grid 2x4 (main)  | 16256 | 128  | symmetric| 8128   | so(128) [D_64]  | (8192,8064)   | so_128 (+) so_128         |

Every row of the table is re-derived by the stand-alone test
`tests/test_char0_2x4.py` with an independent implementation.

## 7. Disambiguation (H396 class)

- **vs `<A,B>` on 2x4 (Theorem OA-2x4).**  `Lie_Q<sum X_i, sum ZZ_e>` has
  dim 2952, centre 1, `D21 + 2B13 + A23 + A27 + 2A3`.  The 18-local-term
  DLA has dim 16256, centre 0.  Different objects, both char-0 certified.
- **vs the predicted sp pattern.**  `sp_128 (+) sp_128` would have dim
  16,512 = 8256 + 8256 > 16256 = dim g.  The half-form is symmetric, giving
  `so_128 (+) so_128`.  The 2x2 row of the same table (symmetric)
  and the 2x3 row (alternating) show the form type is not fixed by the pattern.

## 8. Resource record

`experiments/e148_char0_2x4.py` reuses e143's exact machinery and adds
vectorised/sparse generalisations.  Observed wall clocks on this machine
(p = 2147483647/2147483629): closure 71 s, pair table 1.4 s, Killing 2.5 s,
perfect 1.7 s, form solve 1.4 s, Walsh ranks 2 x 20 s, conjugated minors
~190 s total; total end-to-end in `results/algebra/char0_2x4.json`.
[Dense 16256 x 16384 block stacks would need ~2.1 GB/prime and are avoided by
the sparse route.]
