# The one-rung induction route, and the falsification of SectorSaturationCyclicity

Artifacts: `experiments/e131_sector_saturation_tensor.py`,
`tests/test_sector_saturation_tensor.py`, `results/algebra_growth/sector_saturation_tensor.json`.

Reproduce:

```
PYTHONPATH=src .venv/bin/python experiments/e131_sector_saturation_tensor.py
PYTHONPATH=src .venv/bin/python tests/test_sector_saturation_tensor.py
```

**Status convention.** `[THEOREM]` is proved here and machine-certified over `Q`;
`[COMPUTATION]` is the exact finite calculations recorded in the artifact;
`[COND]` is conditional; `[UNRESOLVED]` is what is deliberately not decided.

---

## 1. Setting

In the `X`-eigenbasis the ladder generators decompose (`proofs/ladder_alll_proof.md`,
verified term by term here) as

```
A_L = 2L - 2N,   B_L = D + F + Ddag,
```

with `D` the hard-core-boson pair creation on a bond, `F = B - (1/16)[A,[A,B]]` the
hopping operator, and `Ddag` the pair annihilation, at bracket depth 3, with
`ad_A` eigenvalues `-4, 0, +4`. Let `psi = |+>^(2L)` be the vacuum and let

```
K_L = the (row-swap, rung-reversal)-invariant 4B-invariant space
```

of the `2xL` ladder configuration Hilbert space, restricted to **even particle
number**. (The parity-inclusive counterpart is twice as large and contradicts
the `T1` formula; the artifact's own scope string is the correct one.) This note
was written to prove that the one-rung tensor induction closes `g_L psi = K_L`;
it reaches the opposite certified conclusion, and that is its result.

## 2. The tensor factorisation of the invariant space

**[THEOREM T1 — Burnside.**] The invariant space dimensions are exactly

```
dim K_L = (2^(2L-1) + 3 * 2^L) / 4,
```

independently computed from the bipartite orbits of the global row-swap and the
rung-reversal reflections (check `T1_Burnside_orbit_dimensions`).

**[THEOREM T2 — tensor extension, exact.**] Writing the `2x(L+1)` configuration
space as `(2xL) tensor` the four states of the new rung, and
`Pi_{L+1}` for the invariant projection of the `2x(L+1)` reflection group, the exact
ranks `Pi_{L+1}(K_L tensor span{00, 11})` are computed by sparse `Fraction` elimination
with pivot-complement direct sums, per `L`:

```
L = 3 -> 4 : rank E = 25 + C = 19,   total 44 = K_4
L = 4 -> 5 : rank E = 88 + C = 64,   total 152 = K_5
L = 8 -> 9 : rank E = 16768 + C = 16384,  total 33152 = K_9
```

(the `L = 8 -> 9` cut factorisation is in the artifact although the cyclic/W
certificates stop at `L = 7`; see §7). The intermediate rows `L = 5 -> 6` and
`L = 6 -> 7` follow the same structure and are recorded row-by-row in the
artifact rather than displayed here.

with the correction representatives and the block histograms stored in the artifact
(`T2_exact_tensor_extension_rank`).

**[LEMMA — the cut is NOT `K`-preserving.**] The rung-hopping operator
`F_new` on the two new bonds maps `K_L tensor span{00}` outside the invariant space:
an explicit `rho`-breaking orbit witness at `4 -> 5` is recorded in the artifact
(`T6_local_cut_not_K_endomorphism`). The invariant space of `K_(L+1)` therefore
cannot be built from the new rung's local occupancy alone: the induction requires
the genuinely new stabilisers, and that is where it breaks.

## 3. Certified annihilator blocks inside `K_L`

**[THEOREM T3 — exact W certificates.**] For every `L` below, there is an
explicitly exhibited subspace `W_L subset K_L` that is *orthogonal to*
`g_L.psi` and 4B-invariant, with

```
dim W_3 = 0,   dim W_4 = 2,   dim W_5 = 10,   dim W_6 = 66,   dim W_7 = 364
```

and the orthogonal complement of `W_L` in `K_L` exactly spans the `g_L.psi` cyclic
subspace. The certificate is formed by sparse `Fraction` elimination (T3), and it
is confirmed independently by a two-sided sandwich: the modular cyclic-closure
rank gives `dim g_L.psi >= K_L - dim span_Q(W)`, and the orthogonal-complement
spanning set gives `dim g_L.psi <= K_L - W_L` over `Q`; one prime suffices because
the complementary side is exact over `Q`.

## 4. The decisive numbers

**[THEOREM — SectorSaturationCyclicity as a universal hypothesis is FALSE; certified
falsity at `L = 4, 5, 6, 7`.]**

*Quantifier, precisely.* The cyclic/W certificates below are per-`L` at
`L = 4, 5, 6, 7`. No falsity (and no cyclicity) at any `L >= 8` is claimed
here. This is a refutation of the *universally quantified hypothesis* in the
wave-13 conditional theorem (which requires `g_Lψ = K_L` for every `L`), and a
sharpening of the certified bound `dim g_L ≥ K_L − W_L` at `L ≤ 7`; it is
**not** a claim of falsity at every `L >= 4`.

With

```
dim U(A,B) psi = K_L - dim W_L   (exact, certified over Q and one prime):
   L = 3 : 14/14      (saturation, consistent with `dim g_3 = 263`)
   L = 4 : 42/44      (W = 2)      <- saturation FAILS here, first place it can
   L = 5 : 142/152    (W = 10)
   L = 6 : 494/560    (W = 66)
   L = 7 : 1780/2144  (W = 364)
```

## 5. What this means

The conditional theorem in `proofs/ladder_alll_proof.md` (Section `dim K_L = 2^(2L-3)+3*2^(L-2)`
and `dim(g_L) >= 2^(2L-3) >= 2^L` IF `g_L psi = K_L`) is now **proved in the form
`dim(g_L) >= dim K_L - dim W_L`**, with `dim K_L` canonical and `dim W_L` certified:

```
liminf (dim g_L)  >=  K_L - W_L,  certified at L = 3..7 as
 14, 42, 142, 494, 1780.
```

Observe: `42 >= 2^4 = 16`, `142 >= 2^5 = 32`, `494 >= 2^6 = 64`, `1780 >= 2^7 = 128` —
the standing conjecture's `2^L` certificates at L=3..7 are recovered as trivial
consequences of much stronger certified values. Observe: `42 >= 2^4 = 16`, `142 >= 2^5 = 32`, `494 >= 2^6 = 64`, `1780 >= 2^7 = 128` —
the standing conjecture's `2^L` certificates at L=3..7 are recovered as trivial
consequences of much stronger certified values.

**[UNRESOLVED]** The asymptotic ratio `W_L / K_L`: 0, 4.5%, 6.6%, 11.8%, 17.0% at
`L = 3..7`. Whether it stays below `1 - 2^L/K_L` for every `L` — which would rescue
`dim g_L >= 2^L` unconditionally — is unproved. The certified data do not exclude
either outcome, and this note does not attempt to forecast it.

## 6. What the wider program inherits

**(a)** The ladder-growth question should now be asked with `W_L`, not
`g_L psi = K_L`: the invariant space is **not** cyclic, so the promising object is
the pair `(extension tensor step, W_L family)` with the orbit witness structure from §3.

(m) The free-fermion no-go via mode counting (`m` modes force the Lie algebra into
`so(2m)`, so a representation on `2^L` qubits would need `m(2m-1) >= dim g_L`):
the certified dimension bound `1780` at `L = 7` forces `m(2m-1) >= 1780`, hence
`m >= 31` modes, against only `2L = 14` qubits. Any stronger certified `dim g_L`
lower bound strengthens this arithmetically, not structurally; it is recorded
here as `m >= 31`, the exact `2^L`-style lower bound, not as a mode-space
identity.

**(c)** No statement is made about `dim g_L`: only its evaluation lower bound
`dim g_L >= K_L - W_L` is certified. `dim g_L >= dim(g_L psi)` could in principle gain
more from other vectors, and the unrestricted `dim g_L` for `L >= 4` is **not determined
here**; it stays in the wave-10/13 tabulation (263 at `L = 3`, 2952 at `L = 4`).

## 7. Scope and honest negatives

* No exact `dim g_L` is claimed for `L >= 4`.
* No asymptotic statement about `W_L` is claimed: the five values 0, 2, 10, 66, 364
  are the full certificate.
* The one-rung induction is certified to FAIL in the way specified by §3 (the new
  hopping operator does not preserve the invariant space); nothing else about tensor
  extensions of this basis is conjectured here.
* The cyclic/W certificates were produced within `process_time` budgets per `L`
  (30 s for `L <= 5`, 60 s for `L = 6`, 300 s for `L = 7`) and stop at `L = 7`.
  The tensor-cut factorisations additionally cover `L = 8 -> 9`; no further `L`
  is claimed.
