# Ladder annihilator at L = 9: exact W_9, the mirror theorem, and the annihilator-route horizon

Artifacts: `experiments/e151_ladder_w9.py`, `tests/test_ladder_w9.py`,
`results/ladder/w9_saturation.json`, `results/ladder/w9_basis.json`,
`results/ladder/w9_progress.jsonl`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py regress
PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py closure 999983 --save-rows
PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py closure 2147483647
PYTHONPATH=src .venv/bin/python experiments/e151_ladder_w9.py assemble
PYTHONPATH=src .venv/bin/python tests/test_ladder_w9.py
```

**Status convention.** `[THEOREM]` is proved here; `[COMPUTATION]` is the
exact finite calculation recorded in the artifact (modular ranks at the two
listed primes, exact-Q kernel validation); `[CONJECTURE]` is the unproved
growth model.  Claim tags: `EXACT-Q` / `MODULAR-TWO-PRIME-LOWER-BOUND` /
`[COMPUTATION]+[CONJECTURE]` as tabulated in the artifact scope.

---

## 1. Setting and the wave-15 question

The `2xL` ladder, `K_L = {even-particle, <tau,rho>-invariant states}` in the
`X`-eigenbasis, the vacuum `psi`, the operators
`A_L = sum_v X_v = 2L - 2N`, `B_L = sum_{edges} Z_u Z_v` — notation identical
to `proofs/sector_saturation_pairing.md` and `proofs/ladder_alll_proof.md`.
The one-rung state-space certificate for the ladder evaluation bound is

```
dim g_L >= dim U(A,B) psi = dim K_L - dim W_L,
```

where `W_L = (U(A,B) psi)^perp` is the pairing kernel — the annihilator of
the Krylov span.  The certified sequence stood at

```
L        3   4   5    6    7     8
W_L      0   2   10   66   364   1822
W_L/K_L  0   4.5% 6.6% 11.8% 17.0% 21.7%
```

with `W_L` growing about x5 - x6.6 per step against a x4 `K_L` base.  The
wave-16 question: compute `W_9` exactly and determine whether the trend — if
it persisted — would make the annihilator certificate `K_L - W_L >= 2^L` go
vacuous at finite `L`.

---

## 2. [THEOREM] Sector decomposition (sector-pure closure)

On the even-parity invariant space `K_L`, `A_L` is diagonal with eigenvalue
`2L - 2k` on the `k`-particle sector `K_L^{(k)}` (`k` even, `0 <= k <= 2L`),
and the eigenvalues are distinct for distinct `k` (over Q and over every
field of characteristic `char > 4L`).  Hence each sector projector `P_k` is a
complex polynomial in `A_L` (Lagrange interpolation), and

```
U(A,B) psi = <psi> closed under A, B
           = <psi> closed under { P_k B : k even }        (sector-pure closure)
span{U}  = direct sum over k of U^{(k)},  U^{(k)} = U cap K^{(k)}.
```

So closure rows may be kept sector-pure and `A`-images of sector-pure rows
(scalar multiples) never need to be generated.

---

## 3. [THEOREM] Particle-hole mirror (half the sectors for free)

Let `P = prod_v Z_v` be the global `X`-eigenbasis bit complement.

(a) `[P, B_L] = 0`, `P A_L P = -A_L`, `P` maps `K_L^{(k)}` bijectively onto
    `K_L^{(2L-k)}`, mapping each `<tau,rho>`-orbit onto an orbit of the same
    size, and `P psi = |all-ones>`.

(b) `psi`'s cyclic module contains `|all-ones>`: with the hard-core
    decomposition `B = D + F + D-dagger` of `proofs/ladder_alll_proof.md`
    (Theorem 1 there: `D = (1/32)[A,[A,B]] - (1/8)[A,B]` lies in the
    enveloping algebra), the element `D^L psi` has particle number `2L`, and
    its coefficient on the fully occupied configuration is the number of
    ordered perfect matchings of the ladder by its edges — at least `L!`
    (fill the rungs in any order).  Hence
    `K_L^{(2L)} = < |all-ones> > subseteq U`.

(c) `U` is `P`-invariant: `P U = alg(PA,PB) P psi = alg(A,B) P psi subseteq U`
    by (b) and module invariance; `P` is a bijection, so `P U = U`, and
    therefore

```
U^{(2L-k)} = P (U^{(k)}),    W^{(2L-k)} = P (W^{(k)}),
```

    where `W = U^perp` and `P` is orthogonal for the configuration pairing.
    In particular `dim W^{(k)} = dim W^{(2L-k)}`: the experimentally observed
    palindromes of the wave-14/15 kernel sector splits
    (`{70,448,786,448,70}` at L = 8 etc.) are theorems.

(d) *Mirror-injected closure.*  Let `low = {k <= L}` and inject: from a stored
    low-sector row, every `k > L` image component is instead added as its
    `P`-image in sector `2L - k`.  The fixpoint span `S` equals `Pi_low U`:
    `S subseteq Pi_low U` because every direct or injected image of a
    `U`-element is in `U` intersector-low by (c); `Pi_low U subseteq S` by
    induction on generating words, using `P_s B (P w) = P (P_{2L-s} B w)` and
    the strengthened hypothesis `P Pi_high V subseteq S`.

Consequences used here: the closure only stores and processes sectors
`k <= L`; the modular rank is
`rank U = sum_{k<=L} (2 - delta_{k,L}) rank U^{(k)}`; kernel lifts for
`k > L` are the `P`-images of the low-sector lifts; the exact-Q validation of
a mirrored basis needs only the mirror invariance of the pairing plus the
standard validation of one half (both halves are validated in the artifact
regardless).

---

## 4. [COMPUTATION] L = 3..8 regression of the new engines

The sector-pure mirror closure + low-sector triangular lift + `P`-mirroring,
triangular back-substitution at `q = 999983` with
`e142_sector_saturation_l8.kernel_lifts_dense` and exact-Q validation with
`e142_sector_saturation_l8.validate_W_over_Q`:

```
L        closure rank       cyclic/K dim      W       W sectors
3        14                 14/14             0       {}
4        42                 42/44             2       {4:2}
5        142                142/152           10      {4:5, 6:5}
6        494                494/560           66      {4:15, 6:36, 8:15}
7        1780               1780/2144         364     {4:35, 6:147, 8:147, 10:35}
8        6562               6562/8384         1822    {4:70, 6:448, 8:786, 10:448, 12:70}
```

agrees with the wave-14/15 certificates on every entry (cyclic dims, kernel
dims, per-sector splits, sandwich identity, full Q-invariance).
`results/ladder/_w9_work/regression.json`; replayed independently by the
standalone test.

---

## 5. [COMPUTATION] The L = 9 certificate

`dim K_9 = (2^17 + 3*2^9)/4 = 33152`, confirmed three independent ways
(Burnside formula, direct `<tau,rho>`-orbit enumeration on `2^17` even
configurations, per-sector Burnside sum; artifact `K9_dimension_record`).

The mirror closure at L = 9 (sectors `0, 2, 4, 6, 8` stored, all
`{A,4B}`-properties as in section 3) gives, at two primes:

* `rank_{F_999983} U = 24566` (12,283 processed low rows, 62,057,796
  reduction steps, 2504 s cpu; peak RSS 1.36 GB);
* `rank_{F_2147483647} U = 24566` (identical row/step counts within 67 steps,
  2612 s cpu; peak RSS 1.37 GB);

with low-sector ranks `{0:1, 2:45, 4:666, 6:3570, 8:8001}` agreeing exactly
at both primes.

The per-sector kernel lifts at `q = 999983` (`bound = 1<<20`, observed
`|c| <= 4`), plus their `P`-images, form `4,293 + 4,293 = 8,586` integer
orbit-sum vectors that validate over Q: vacuum-orthogonal, homogeneous in one
particle sector each, linearly independent, and closed under `A_L` and
`4 B_L`.  By the annihilator-recursion theorem of
`proofs/sector_saturation_pairing.md` (section 2) this family lies in
`(U psi)^perp`, so `dim U <= K_9 - 8586`; with the modular ranks at both
primes:

```
dim_Q U(A_9, B_9) psi = K_9 - W_9 = 33152 - 8586 = 24566   exactly,
W_9 = 8586,  sector split {4:126, 6:1134, 8:3033, 10:3033, 12:1134, 14:126}
             (palindrome under k -> 18-k, as required by section 3),
dim g_9 >= K_9 - W_9 = 24566 >= 512 = 2^9   (certificate form preserved).
```

Per-sector deficits (`dim K_9^{(k)} - rank U_9^{(k)}`):

```
sector k       0   2   4    6     8     10    12    14   16  18
dim K_9^(k)    1   45  792  4704  11034 11034 4704  792  45  1
deficit        0   0   126  1134  3033  3033  1134  126  0   0
```

Claim tags: `W_9 = 8586`, the cyclic dim `24,566` and the per-sector
deficits are `EXACT-Q`; the two ranks are
`MODULAR-TWO-PRIME-LOWER-BOUND`.

---

## 6. [COMPUTATION]+[CONJECTURE] The growth model and the horizon

Exact points `(L, W_L)` for `L = 3..9`:

```
W = (0, 2, 10, 66, 364, 1822, 8586),
consecutive ratios 5.0, 6.6, 5.515151..., 5.005494..., 4.712404...,
W_L/K_L 0, 4.55%, 6.58%, 11.79%, 16.98%, 21.73%, 25.90%.
```

*Model.* `W_L = C r^L`, least squares on `log W` over `L = 4..9` (mpmath
`dps = 50`): `r_hat = 5.41894...`, maximal relative log-space residual
`11.7%`.  Two-point anchored variants `W_L = W_9 r^(L-9)` for
`r in {5.0, 6.6}` (the prior ratio band) and for the new ratio
`r = 8586/1822 = 4.712404...` bound the plausible family.  `W_9 = 8586`
deviates from the last-prior-ratio prediction `1822 * 5.005494... = 9120.01...`
by `5.86%` — below the 20% verdict-change threshold recorded in the brief,
so the refutation trend persists unchanged in classification.

*Horizon.*  Solving `C r^L = K_L = (2^{2L-1} + 3*2^L)/4`:

* LSQ model (`r_hat = 5.41894`): continuous crossing `L* = 13.057...`;
  first integer with model `W_L >= K_L` is `L0 = 14`;
* anchored band: `r = 5.0` gives `L* = 15.003...` (first integer 16);
  `r = 6.6` gives `L* = 11.678...` (12);
  `r = 4.712404...` (new ratio) gives `L* = 17.172...` (18).

Under the same models, the first `L` where the certificate `K_L - W_L < 2^L`
would go vacuous is `L = 14` (LSQ) and `L = 18` (anchored at the new ratio).

**Verdict.** `W_9 = 8586` is exact; the refutation continues to grow
(`W_9 = 4.71 * W_8`, within 6% of the last-ratio x5.005 trend).  The
`[CONJECTURE]` horizon: if the fitted geometric law persisted, the
annihilator certificate `K_L - W_L >= 2^L` would go vacuous at
`L0 = 14` (LSQ; bank range `12 <= L0 <= 18` from the anchored models), and
`W_L >= K_L` (empty-annihilator) at `L_hat = 14`, i.e. `L0` and the empty-
kernel crossing coincide under the LSQ fit.  Seven points are not a proof of
a law; the tabulated per-`L` model values in the artifact
(`extrapolation_table_L10_L20`) make the gap between data and model visible
at each step.

---

## 7. Scope and honest negatives

* The proof of sector decomposition (section 2) and of the mirror lemma
  (section 3) is complete; the L = 9 numbers rest on the certified sandwich
  of section 5, exactly as in e142: one prime suffices for the lower rank,
  the second prime is redundancy against implementation error.
* The growth model is fitted to seven exact points; it is
  `[COMPUTATION]+[CONJECTURE]`, no law is proved, and the crossing `L0` is a
  conditional estimate, not a theorem.
* `dim g_9 >= K_9 - W_9` remains an evaluation bound via the stabiliser
  module; no statement about `dim g_9` itself.
* The sibling repository result `ladder_nonstat.json` (literal-Pauli
  512-certificate at L = 9) is unaffected; this note strengthens the same
  `2^9` bound to `24566`.
* The engine's measured cost at L = 9: ~2.5 - 2.6 K seconds cpu per prime,
  62.1M reduction steps, 12,283 processed low rows, max echelon row nnz
  5,543, peak RSS 1.4 GB; the exact-Q lift + validation of all 8,586 kernel
  vectors took 283 s.  The watchdog `w9_progress.jsonl` records block-level
  checkpoints throughout (a kill would lose at most one block).
