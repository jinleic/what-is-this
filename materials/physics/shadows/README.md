# `shadows/` — contractive-shadow exact census

**Status: ACTIVE (gates A/B closed 2026-08-29, PROVED; gate C closed —
exact 361/353 weight improvement over U_ct found at k=4; follow-up
campaign pending).** Exact stabilizer/Clifford arithmetic censuses of
classical-shadow measurement ensembles. Everything here is exact
GF(2)/rational computation — no Monte Carlo, no RNG, single core.

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
(Hingane & Koh, "Real Classical Shadows with Noise", 2026-08-19) read at
abstract level — see Follow-up below.

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
   CZ mesh, CX chain, CX mesh, U_ct∘CZ). Result: **the CX chain /
   CX-mesh layers compose with U_ct to weight 361/6561 vs U_ct's
   353/6561 at k=4 — shadow norm ratio 353/361 (~2.2% smaller than the
   paper's optimum family at that k; exact Fractions)**. Mesh_X/mesh_Y are
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
```

Exactness argument: Paulis are GF(2) bit-label ints; Clifford class maps are
symplectic 4-tuples verified via the symplectic Gram matrix; weights are
`fractions.Fraction` counts over the full 3^k token space. No RNG anywhere;
floats only under artifact keys tagged `numerical*`. Exhaustive ⇒ universal.

## Follow-up (next campaign, not started)

- Channel-resolved exact variance table for the real-Clifford ensemble
  (arXiv:2608.18935) at d ∈ {4,8,16} under depolarizing + amplitude
  damping; search for an advantage-flipping channel (variance ratio < 1).
- Extend gate C: family (b) factorization shortcut (block-constant matchings)
  to push k past 8; per-k norm table for the U_ct∘CX-chain family to see
  whether the 353/361 edge decays with k (decides if it is a k=4
  finite-size effect). NOTE: gate C is descriptor-frozen at k ≤ 8;
  changing k requires a new revision entry in pre_statement.md.

## Layout

- `campaigns/2026-08-29T232134Z_4f1a9cab_e274dc65db35/` — frozen snapshot
  (gates_a/b/c.json + SUMMARY.json + MANIFEST.txt).
- `pre_statement.md` — frozen gates, conventions, Revision 1.
- `src/shadows_exact.py` — GF(2) symplectic library + all gate
  implementations + selftest.
- `src/run_campaign.py` — campaign runner (freezes artifacts, hashes them).
- `scratch/` — non-authoritative exploration (currently empty).
