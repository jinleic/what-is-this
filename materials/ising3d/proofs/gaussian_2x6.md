# Gaussian purity of the 2x6 open Ising layer at t = 1/3 — (P,rho)-block certificate

**Status: e154 run in flight (`results/spectral/gaussian_2x6_work/`); this note is finalized
when the decisive-block numbers land. Sections §2-§4 (machinery + calibration) are proven;
§5 records the decisive outcome. Classification tags: [THEOREM] (ceiling + pipeline-correct case),
[COMPUTATION] (exact modular certificates).**

## 1. Question and invariant

For the open `2x6` layer (`n = 12`, 16 bonds), at `t = tanh(K*/2) = 1/3` (`exp(2K) = 5/3`,
the canonical curve point of Theorem S era), is the physical `P = prod X = +1` sector of the
bond-operator's transfer representative `R = P_t diag(q^e) P_t` (an exact rational symmetric
matrix up to scalar) the even half of a 12-mode Gaussian subset-product multiset (i.e. is the
layer's spectrum free-fermionic), in **either assignment of physical `P` to fermion parity**?

**Theorem (ceiling, `proofs/parity_pair_product.md`, m = 12).** Any *subset* `U` of one
parity half of an `m`-mode Gaussian subset-product multiset has at most
`A_same(m) = sum_{2<=r<=m, r even} C(m,r) 2^{m-r}` distinct within-`U` pair products
(y_i y_j). Both assignments of physical P to fermion parity test *within*-half products, so
both have the same ceiling (A_cross applies only to P+xP- cross products). For the
`(P+,rho+) = 1056`-dimensional block: `C(1056,2) = 557040`, kill-cap
`557040 - A_same(12) = 295415`.

**Certificate.** Compute the characteristic polynomial `C(z)` of the pair-product
(multiplicative) system on the block over two independent primes (`1000003`, `2000003`,
both `> 557040`); an exact gcd / early-exit with degree `< 295415` proves the block is NOT
a Gaussian half-subset, in both primes, in both parity assignments.

## 2. Route

1. Exact integral model `W` of the layer (e136 split), `M = W mod p`.
2. Joint `(prod-X, rho)` character basis `B` over the abelian involution group; basis rows
   integer orbit-sums; `gram` denominators in `{1,2,4,8}`; the restriction
   `R|_block = diag(gram)^{-1} B^T M B`; certificate on `Y = 8 * that` (integral).
3. Power traces (float64-BLAS with int64 round-trip first-product audit; exact 11-bit
   split-fallback) -> Newton charpoly -> Cayley-Hamilton trace recurrence to `2*Npairs`
   -> pair power sums `q_m = (s_m^2 - s_{2m})/2` -> Newton pair polynomial.
4. Euclid `gcd(C, C')` with early exit at the cap.

[THEOREM-grade honesty]: all linear-algebraic identities used (Newton identities, CH
recurrence, pair-power identity, Euclid correctness) are **proved** in
`proofs/parity_pair_product.md` and rechecked by machine controls below.

## 3. Controls (all PASS, two primes)

- `chain n=4`: gcd == 55 (the Gaussian floor); `chain n=6`: == 1351; both exact modes.
- `2x3 layer t=1/3` stored-digest replay: charpoly + pair-poly sha256 match
  `results/spectral/pair_product_scale.json` on both primes (gcd == 385).
- `chain n=10` at the Gaussian scale: gcd == 465751 (pairs 523,776).
- **Structural guard**: for the `2x3` layer at `t=1/3`, for ALL FOUR `(flip,rho)`
  characters, `cp(block) | cp(full)` divides exactly (characteristic polynomial of the
  restriction divides that of the full matrix).
- **Pilot (geometry-corrected)**: the 1D chain with the true 1D reflection, same
  character classes and same scale: `pilot_ch12_pp_gcd_forced` = **exact gcd 442067**
  on BOTH primes, satisfying the `>= 295415` forced floor with margin — validating the
  entire chain at the decisive scale.

## 4. Bugs retired on the way (each caught by a contradiction, never landed)

1. `orbit_character_data` stabilizer test used per-generator checks instead of the
   character product (dims 4064 != 4096 artifact).
2. `next_prime` could return EVEN numbers (`next_prime(2098176) == 2098178`).
3. `poly_divmod` fixed-index padding misaligned after zero quotient coefficients
   (unit `a = q.b + r` failures; rewritten as re-trimming schoolbook, verified on
   500 random divisions and manual schoolbook cases).
4. Pilot's kappa geometry: the 2xL-layer (row-pair) reversal applied to a 1D-chain
   model — non-invariant blocks; fixed to true 1D reflection.
Every landed stage file is content-keyed `(code-hash, schema, model, char, prime, label)`,
written atomically; Euclid resume validates the input digest of its starting remainder pair.

## 5. Decisive block — verdict: NO-GO certified on both primes **[THEOREM]**

Recorded (`results/gaussian_2x6.json`):

| prime | mode | cert_upper (deg gcd bound) | cap | steps |
|---|---|---|---|---|
| 1000003 | exit | 295414 | 295415 | 261625 |
| 2000003 | exit | 295414 | 295415 | 261625 |

The `(P+ = prod X, rho+)` block of the open `2x6` layer's transfer operator at `t=1/3`
satisfies `deg gcd(C, C') <= 295414 < 295415 = C(1056,2) - A_same(12)`, i.e. the block's
within-half pair-product multiset has **at least A_same(12) + 1 = 261626 distinct values** —
at least one more than ANY parity-half of a 12-mode Gaussian can carry. Therefore (both
primes, both assignments of physical P to fermion parity, within-half products in both
cases):

**THEOREM (wave 17).** The physical `P = +1` sector of the open `2x6` Ising layer at
`t = 1/3` is not the even/odd half of any 12-mode Gaussian subset-product multiset. The
largest layer previously certified was `2x5`.

*Precision discipline.* The early-exit proves `deg gcd <= 295414` (a one-sided bound);
the excess distinct-product count is therefore at least 261626, NOT computed exactly —
the Euclid's own step count (261625) is unrelated to the distinct-root count and must not
be interpreted as a certificate.
Controls/zero-failure: pilot at the decisive scale exact `442067` on both primes; every
structural judge in the producer passed; `failures: []`; the exit states and all inputs are
content-keyed, atomically stored under `results/spectral/gaussian_2x6_work/`. An independent
list-schoolbook replay reproduces the exit at deg 295414 on both primes with bit-identical
remainders, and the modular reduction step of the certificate uses the established
two-prime protocol (integral exact content + per-field ceiling counting).

## Artifacts

`results/spectral/gaussian_2x6_work/` (content-keyed stage states),
`results/gaussian_2x6.json` (final), `experiments/e154_gaussian_2x6.py`,
progress journal `results/gaussian_2x6_progress.jsonl`.
