# `shadows/` — contractive-shadow exact census

**Status: ACTIVE (gates A/B closed, PROVED; gate C closed with the k=4
361/353 weight improvement; k-extension closed EDGE-VANISHES; every
pinned gate-C family is now classified for all k≥1 by exact algebra:
only pure `cx_chain`/`cx_mesh_lex` at k=4 strictly improve on U_ct.
No further finite-k campaign.
Noisy-variance remains FROZEN-CERTIFIED — exact 132-cell no-flip
regression).** Exact GF(2)/rational computation and symbolic derivation;
no Monte Carlo, no RNG.

## Reference

Wu, Wang, Yao, Zhai, You, Zhang, *Contractive unitary and classical shadow
tomography*, npj Quantum Information 12, 86 (2026), DOI
10.1038/s41534-026-01227-w (= arXiv:2412.01850, v1 2024-11-28) —
deterministic commuting Clifford U_ct = ∏_{i<j} exp(iπ/4 Z_i Z_j), claimed
shadow norm ≈ 2×1.8^k for size-k Paulis vs the ~2^k random-Clifford
barrier; closed-form Pauli weight
w(k) = ½[3^-k + (-1)^k 9^-k] + ½[(5/9)^k − 9^-k].

Paper read first-hand 2026-08-29 (full text incl. Methods, Table 2
verbatim from the tables/2 landing page). Follow-up paper arXiv:2608.18935
(Hingane & Koh, "Real Classical Shadows with Noise", 2026-08-19) read
first-hand 2026-09-03 from the arXiv HTML v1 — Eqs. (28), (29), (33),
(36), (37), (39), (44) transcribed verbatim in pre_statement.md
Revision 3 — see Follow-up below.

## Gates — CLOSED, all three, campaign
`campaigns/2026-08-29T232134Z_4f1a9cab_e274dc65db35/` (194 s wall,
`nice -n 10`, one core; artifact sha256-12 `e274dc65db35` over
gates_a+b+c.json)

1. **Gate A (census) — CONFIRMED, [PROVED].** Exhaustive enumeration of all
   65,536 GF(2) 4×4 maps → exactly 720 symplectic (hard-checked = |Sp(4,2)|)
   → 720 × 16 = 11,520 mod-phase Clifford labels [DERIVED count; the paper
   never states it]. Paper's Methods lemma (max 4 of the 9 size-2 Paulis
   contractible to size-1) CONFIRMED: histogram over the 720 class actions
   is **bimodal — 72 actions contract 0 of 9, 648 contract exactly 4,
   nothing in between, max = 4**. Anchor exp(iπ/4 Z Z) reproduces the paper's
   exact 4-contracted / 5-unchanged pattern (Methods list reproduced
   elementwise). All 9 achieving sets = (row ∪ column)\{cell} of the paper's
   Table 2 [DERIVED, hand-checked]. Sharpening: no Clifford is
   "partially contractive" — every two-qubit Clifford contracts either 0 or
   4 [PROVED, same census].
2. **Gate B (re-derive) — PASS, [PROVED for k≤10 / integer-exact k≤16].**
   Eq. (3) verified per-string exhaustively for all 3^k strings, k≤10 (both
   conjugation directions — the class-level action of the π/4 roots is
   direction-invariant); size histogram == Eq. (5); exact mesh weight ==
   Eq. (4) for k≤10 by Fractions and for k≤16 via the integer Eq. (5)-sum;
   Eq. (7) identity-insertion verified 27/27 (k̃+q ≤ 7) with per-string
   Eq. (6) checks; random-Clifford derivative 1/(2^k+1) re-derived exactly
   16/16. Shadow norm ∥O∥²/(2·1.8^k) rises monotonically 0.735 (k=2) →
   0.99972 (k=16) [NUMERICAL tag only].
   Frozen odd-N_XY mode at k≤10: m = 2⌊k/3⌋+1 [DERIVED, exact-count
   maximizer] — consistent with the paper's "peak near m/k ≈ 2/3".
3. **Gate C (search) — IMPROVEMENT-FOUND, [PROVED within pinned scope].**
   11,919 + 1,701 exact family composites at k ≤ 8 (all single-letter
   meshes X/Y/Z, all near-perfect matchings with lettered pair rotations,
   CZ mesh, CX chain, CX mesh, U_ct∘CZ). Result: **the pure `cx_chain`
   and `cx_mesh_lex` maps have weight 361/6561 versus 353/6561 for
   `U_ct` at k=4, giving shadow-norm ratio 353/361 (~2.2% smaller; exact
   Fractions)**. The literal `cx_chain`-after-`U_ct` composite is a
   different, excluded map with weight 1/81 at k=4. Mesh_X/mesh_Y are
   weight-identical to U_ct (letter symmetry, exact 1/1 ratios all k).
   Scope-limited per pre_statement.md Revision 1 — NOT an obstruction
   proof, NOT claimed to beat 2×1.8^k asymptotically (single k point).

Frozen gate text + exact acceptance criteria: `pre_statement.md`
(2026-08-29; Revision 1 appended there documents the gate-C comparison
handedness fix — code and statement agree on w_M > w_ct).

## How to run

```
nice -n 10 physics/.venv/bin/python physics/shadows/src/shadows_exact.py   # selftest (~0.3 s)
nice -n 10 physics/.venv/bin/python physics/shadows/src/run_campaign.py    # full gates (~3.3 min)
nice -n 10 physics/.venv/bin/python physics/shadows/src/shadows_noisy_variance.py --check   # 132-cell noisy-variance recompute (~0.2 s, writes nothing)
```

Exactness argument: Paulis are GF(2) bit-label ints; Clifford class maps are
symplectic 4-tuples verified via the symplectic Gram matrix; weights are
`fractions.Fraction` counts over the full 3^k token space. No RNG anywhere;
floats only under artifact keys tagged `numerical*`. Exhaustive ⇒ universal.

## Follow-up

### k-extension — RAN 2026-08-30, closed EDGE-VANISHES; no further k campaign (decided 2026-09-03)

`campaigns/2026-08-30T114124Z_ba4f8430_65d2435684e3/` (pre_statement.md
Revision 2, appended before any k > 8 run; carrier `cx_chain`, regression
gate PASS, 23.8 s single core nice-10, exact Fractions, no RNG). Exact
per-k weight differences w(cx_chain) − w(U_ct):

| k | diff | verdict |
|---|---|---|
| 9  | −54080/129140163      | EDGE-VANISHES |
| 10 | −345760/1162261467    | EDGE-VANISHES |
| 11 | −698240/3486784401    | EDGE-VANISHES |
| 12 | −36680152/282429536481 | EDGE-VANISHES |

The k=4 edge (norm ratio 353/361) does not recur in the enumerated range:
`cx_chain` ties `U_ct` at k=5 and has lower weight at every k=6..12. Across
the extension rows k=9..12, the absolute weight gap
|w(`cx_chain`)−w(`U_ct`)| narrows monotonically, while
w(`cx_chain`)/w(`U_ct`) decreases monotonically, so the relative and
shadow-norm disadvantage widens. The largest k explicitly enumerated is 12
(`largest_k_completed: 12`). **Decision: DO NOT run a further k-extension
campaign.** Further finite-k rows cannot change the already-closed k≤8
Gate-C verdict. The analytic proof below now settles the k→∞ law without
further enumeration. In the measured
range k=4..12, `cz_mesh` matches `U_ct`, `cx_mesh_lex` matches `cx_chain`,
and `uct_then_cz_mesh` remains below `U_ct`; family (b) was searched only
through k=8. Incremental extension-row times were 0.468, 1.535, 5.042, and
16.356 s, growing by about 3.25× per k (k=12 accounted for ≈16.4 s of the
23.8 s run). Any k>12 run still requires a new pre_statement revision first.

### Analytic closure — pure cx_chain for every k

**[DERIVED; PROVED within the pinned model, 2026-09-05].** For the open,
sequential CNOT staircase in `src/ext_k.py::build_layers`, with independent
uniform full-support X/Y/Z input letters and weight
$w=\mathbb{E}[3^{-\mathrm{output\ size}}]$, the only strict improvement
over $U_{\rm ct}$ is **k=4**. There are ties at **k=1,2,3,5** and strict
losses for **every integer k≥6**. This is an algebraic result, not a
new finite-k campaign or an extrapolation from the stored rows.

**Transfer construction.** The pinned ascending CNOT order gives
$x'_i=P_i=x_0\oplus\cdots\oplus x_i$ and
$z'_i=z_i\oplus z_{i+1}$, except $z'_{k-1}=z_{k-1}$.
An interior output is identity exactly when $P_i=0$ and $z_i=z_{i+1}$.
Carry the four states $(P_i,z_i)=(0,0),(0,1),(1,0),(1,1)$.
For the next input letter $(x,z)$, move to $(P_i\oplus x,z)$ with its
probability $1/3$, multiplying by 1 for an identity output or $1/3$
otherwise. The initial and final boundaries give

```text
M = 9T = [[0,1,3,1], [0,3,1,3], [1,1,0,1], [1,1,0,1]]
u = (0,1,1,1)/3
v = (1,1/3,1/3,1/3)^T
w_cx(k) = u (M/9)^(k-1) v
```

The last two components of $v$ are equal and remain equal under $M$,
so reduce to

```text
R = [[0,1,4], [0,3,4], [1,1,1]]
f = (0,1,2),  g = (3,1,1)^T
A_k = 9^k w_cx(k) = f R^(k-1) g
fR = (2,5,6),  fR^2 = (6,23,34)
f(R^2 - 3R - 8I) = 0
```

Consequently **$A_1=3$, $A_2=17$ and
$A_{k+2}=3A_{k+1}+8A_k$ for every k≥1**. The recurrence follows from
the transfer identity, not a fit. With
$\lambda_\pm=(3\pm\sqrt{41})/2$, its closed form is

$$
w_{\rm cx}(k)=
\frac{\lambda_+^{k+1}-\lambda_-^{k+1}}{\sqrt{41}\,9^k}.
$$

**All-k comparison.** The pinned Eq. (4) gives
$N_k=9^k w_{\rm ct}(k)=(5^k+3^k+(-1)^k-1)/2$. Set $D_k=N_k-A_k$.
Then

$$
D_{k+2}=3D_{k+1}+8D_k+F_k,\qquad
F_k=5^k-4\,3^k-2(-1)^k+5.
$$

For k≥3, $5^k>4\,3^k$ and $5-2(-1)^k\ge3$, hence $F_k>0$.
The initial values are
$(D_1,\ldots,D_7)=(0,0,0,-8,0,240,2880)$.
Induction from $(D_6,D_7)>0$ proves $D_k>0$ for every k≥6.
Because larger weight means smaller squared shadow norm, this establishes
the complete classification above; a weight tie need not be a size-
distribution identity.

**Asymptotics.** The dominant coefficient is nonzero,
$|\lambda_-|<\lambda_+$, and $\lambda_+<5$ because $\sqrt{41}<7$.
Therefore

$$
\frac{w_{\rm cx}(k)}{w_{\rm ct}(k)}
\sim\left(1+\frac3{\sqrt{41}}\right)
\left(\frac{3+\sqrt{41}}{10}\right)^k
\longrightarrow0.
$$

The reciprocal shadow-norm ratio diverges. The chain's squared-shadow-
norm exponential base is $18/(3+\sqrt{41})$, strictly between $9/5$ and
2. The bounded k=4 edge is not an asymptotic advantage.

**Verification and provenance.** Independent reviews checked the transfer
construction and the all-k induction. A four-state integer-polynomial
check, without token-string enumeration, matched all nine stored k=4..12
weight rows and every coefficient in all four stored k=9..12 size
histograms. The [verification record](../../data/shadows-cx-chain-analytic-20260905.json)
contains the matrices, recurrence, exact differences and input byte hashes.
The legacy extension manifest's `65d2435684e3` is a preregistration/source
identifier computed by `ext_k.py::artifact_id`, not a raw digest of
`gates_ext.json`; no legacy raw-output SHA ledger is claimed.

No campaign or source file was changed, and no k>12 row was generated.
The proof is only for the pinned pure open chain: it does not establish
claims for inserted identities, noise, rings, other gate orders, or
family (b). The historical finite-k verdicts and the separate noisy-
variance FROZEN-CERTIFIED result remain unchanged.

### Analytic closure — complete pinned-family comparison

**[DERIVED / PROVED within the pinned model, 2026-09-05.]** Together with
the chain recurrence above, the following identities classify every
gate-C family for every integer k≥1. The token law is uniform iid X/Y/Z
on all k sites, and larger `w=E[3^-size]` means smaller shadow norm.
This is a statement about those families, not arbitrary Clifford circuits.

| Pinned family or layer | Exact weight or distribution | Comparison with U_ct |
|---|---|---|
| Uniform-letter complete meshes X/Y/Z | Same size histogram as U_ct | Tie at every k |
| Near-perfect matchings of same-letter pair rotations | `(17/81)^m * 3^-r`, k=2m+r, r∈{0,1} | Ties k=1,2; strict losses for every k≥3 |
| Pure ascending `cx_chain` and control-major lex `cx_mesh_lex` | Same histogram; chain recurrence above | Ties k=1,2,3,5; sole strict win k=4; strict losses for every k≥6 |
| `cz_mesh` | Same per-token output size as U_ct | Tie at every k |
| `uct_then_cz_mesh` | `3^-k` | Tie k=1; strict losses for every k≥2 |

**Disjoint-pair factorization.** For a pair carrying letter L, its class
map is the transvection `v → v xor B(v,L_i L_j)*(L_i L_j)`.
Exactly one input letter equals L in four of the nine full-support pair
tokens; those four contract to size1. The other five commute with the
generator and stay size2. Thus the pair's mean weight is `17/81`.
Disjoint pairs have independent inputs and additive output sizes; an
uncovered site stays size1. Every matching and every assignment of its
pair letters therefore has histogram polynomial
`(4t+5t^2)^m*(3t)^r` and weight `(17/81)^m*3^-r`.

For an exact comparison, put `N_k=9^k*w_ct(k)` as above and

```
E_m = ((17+8)^m + (17-8)^m)/2 - 17^m
    = sum_{2<=j<=m, j even} binom(m,j)*17^(m-j)*8^j.
```

Then `E_0=E_1=0` and `E_m>0` for m≥2. At even k=2m the
numerator gap `N_k-17^m` is `E_m`; at odd k=2m+1 the gap
`N_k-3*17^m` is `3*E_m+25^m-1`, positive for m≥1.
This proves the matching row without another family search.

**Layer identities.** All equations below are class-level GF(2) equations
with sites numbered 0..k−1; xor is addition.

- The lex all-pairs CNOT mesh has
  `x'_i=x_i xor x_(i-1)` (`x_(-1)=0`) and
  `z'_i=xor_{j>=i} z_j`. These telescope against the chain's prefix-x,
  adjacent-z equations, so the lex mesh is the chain's inverse.
  More importantly, it is `C * cx_chain * C`, where C is site reflection
  followed by H on every site. C permutes the full-support input tokens
  and preserves output size; hence the two full histograms agree.
  Inversion alone would not establish that histogram claim.
- Put `q=xor_i x_i`. U_ct keeps x and sends
  `z_i → z_i xor (k mod 2)*x_i xor q`; the CZ mesh keeps x and sends
  `z_i → z_i xor x_i xor q`. They differ only by a local S layer
  (`z_i → z_i xor x_i`) for even k, so their per-token sizes agree.
  Composing them cancels q and leaves
  `z_i → z_i xor ((k+1) mod 2)*x_i`: identity for odd k, local S for
  even k. Thus `uct_then_cz_mesh` preserves full support exactly.
  Its gap is `N_k-3^k=(5^k-3^k+(-1)^k-1)/2`, zero only at k=1
  and positive for every k≥2.
- A uniform-letter mesh is `R * U_ct * R^-1`, with the same one-site
  Clifford sending Z to the chosen letter at every site. The inner map
  permutes the input token ensemble and the outer map preserves size,
  proving equality of the full histograms for X/Y/Z meshes.

Two independent symbolic reviews validated the matching and layer proofs.
Main then checked all eight layer identities on complete GF(2) bases at
k=1..12 (96 exact map checks), the three pair transvections, and matching
direct sums at k=1..8. These checks took less than one second combined;
they did not enumerate token strings or matching families. The
[verification record](../../data/shadows-pinned-family-analytic-20260905.json)
pins the reviewed source bytes and records the exercised checks.

No campaign, scientific source, frozen payload, or terminal verdict was
changed. Largest empirical k remains12; exhaustive matching-family
enumeration still ends at k8. The k1 statements use the empty-gate
identity, not a new experimental row. No claim covers other CNOT orders,
overlapping or repeated pairs, inserted identities, rings, or noise.

### Noisy-variance gate — RAN 2026-09-04, closed FROZEN-CERTIFIED

`campaigns/20260904T034040Z_bed4a462_be2fd1106eed/` (gate
`noisy-variance-r7`, pre_statement.md Revision 7; closed
2026-09-04T03:41:22Z; prereg sha256 `8f9d51ba…9ad0fbf2`, module snapshot
`shadows_noisy_variance.snapshot.py` sha256 `156b09c0…59ada5a0`,
`results.json` sha256 `fcb733b7…ddc82d2`, `summary.json` sha256
`36885050…4675227`; all six frozen artifacts re-hash clean against the
dir's `sha256s.txt`). Closed-form `fractions.Fraction` arithmetic over the
frozen grid — no sampling, no RNG, single core.

Grid, frozen before the run: d ∈ {4, 8, 16} × two channels (depolarizing
60 cells + amplitude damping 72 cells = 132 unique cells) × three anchors
(GHZ projector, maximally-mixed Pauli-Z, GHZ Pauli-X). All 31
preregistered regression identities hold: the Eq. (37) second-moment
ratio equals the direct per-cell ratio in all 132 rows, and the Eq. (44)
variance ratio is ≥ the second-moment ratio in all 132 rows. Result:

- **Second-moment ratio** (Eq. 37): min **1955/1539** at depolarizing
  d=16, p=1, GHZ projector; max **17/9** at depolarizing d=16, p=1,
  maximally-mixed Pauli-Z.
- **Variance ratio** (Eq. 44): min **6815/4671** at depolarizing d=16,
  p=3/4, GHZ projector; max **2** at depolarizing d=4, p=1, GHZ Pauli-X.
- **No flip anywhere**: `flip_cells: []`, no ratio ≤ 1 in any of the 132
  cells.

Scope, frozen before the run (Revision 3): the paper proves the no-flip
result analytically (Corollary 3.11 / Eq. (44), with R ≥ ρ_S(d) > 1 for
every d ≥ 2 on the invertible domain 1 < β ≤ d). This campaign is
**regression-grade exact confirmation of the formula chain on the frozen
grid plus margin quantification — not a discovery search and not a proof
beyond the grid**, and it is not evidence on whether the earlier pure
`cx_chain` k=4 bounded edge persists. A ratio ≤ 1 in any cell would have
indicated instrument error first.

Provenance: an earlier r4-instrument run
(`campaigns/20260904T031748Z_cbbe7c0b_fbc8a23e2ed9/`) was closed
SUPERSEDED 2026-09-04 before any result was used — its cardinality checks
could not reject a same-sized wrong-p ladder and its source snapshot was
read only after evaluation; Revisions 5–7 hardened the instrument without
changing any grid value, formula, threshold, or scientific branch.

### Control-plane note

`state.json` now carries `latest_campaign:
20260904T034040Z_bed4a462_be2fd1106eed` / `FROZEN-CERTIFIED`, imported
from that run's `status.json` — the first conforming control-plane run
dir in this target. The five pre-control-plane run dirs (2026-08-29/30)
remain name-non-conforming with no `status.json`, so `campaign.py state
refresh` still cannot import a verdict from them; the evidence is intact
and hashed inside each dir (`SUMMARY.json`, `MANIFEST.txt`); nothing is
edited inside a frozen dir to fix a cosmetic field.

## Layout

- `campaigns/2026-08-29T232134Z_4f1a9cab_e274dc65db35/` — frozen snapshot
  (gates_a/b/c.json + SUMMARY.json + MANIFEST.txt).
- `campaigns/2026-08-30T114124Z_ba4f8430_65d2435684e3/` — frozen
  k-extension snapshot (k=9..12 EDGE-VANISHES).
- `campaigns/20260904T034040Z_bed4a462_be2fd1106eed/` — frozen
  noisy-variance snapshot (132-cell no-flip regression,
  FROZEN-CERTIFIED).
- `pre_statement.md` — frozen gates, conventions, Revisions 1–7.
- `src/shadows_exact.py` — GF(2) symplectic library + all gate
  implementations + selftest.
- `src/run_campaign.py` — campaign runner (freezes artifacts, hashes them).
- `src/shadows_noisy_variance.py` — noisy-variance runner (closed-form
  Fractions; snapshots module + prereg into the run dir).
- `scratch/` — non-authoritative exploration (currently empty).
