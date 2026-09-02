# PRE-STATEMENT — exact vectorized-RREF replay of the fixed direct-route toy census

Created first, before this campaign's source, manifest, state, controls, or census. Parent campaign: `../2026-08-31T08-44Z_DA168C79/`.

## Why a successor is needed

The parent source is mathematically frozen and remains the sole released low-priority process. Its exact scalar Python Gauss–Jordan kernel repeatedly eliminates hundreds of rows across hundreds of columns. The first immutable production cell, `state/cell_8.json`, took minutes after controls and already contains the complete registered P1–P4 evidence. At campaign creation the parent is still running and must not be interrupted, edited, or overlapped by this successor.

This campaign changes only the execution of exact finite-field row operations and removal of repeated per-cell matrix/kernel construction. It changes no field, matrix entry, branch, quantifier, predicate, output semantics, or scope. Release is forbidden until the parent process has exited or reached its registered two-hour cap.

## Pre-campaign development disclosure

Before this pre-statement, Main prototyped a NumPy table-indexed Gauss–Jordan routine in a persistent in-memory kernel and compared its complete return tuple `(rows,pivots,feasible,particular)` against the parent's scalar routine on 300 deterministic random small affine systems. All 300 matched. That probe wrote no campaign file and is **DEVELOPMENT EVIDENCE ONLY**; it is not imported as a control or result here.

## Frozen parent anchors

The successor must hard-hash the parent source ledger and every copied source before import. Two already-complete parent artifacts are immutable semantic anchors:

- `state/controls_result.json` SHA-256 `41652388f5ce7b7c8db82648e06be56bf21c7a0fe0e56e9b4020b1f61ccc8570`;
- `state/cell_8.json` SHA-256 `b7c4d97607f43385f14b4d9b1e1cf43734b809727214176b4f5e1881c075edaa`.

The beta-8 anchor records P1 PASS on all five U3 and five U4 branches, P2 FAIL (33 accepted labels rather than the expected five), P3 FAIL (five mapped branches and 28 unmapped labels), and P4 FAIL (rank 36, nullity 244, both full and block admission false). These are parent-source results, not yet successor results.

## Exact implementation delta

1. Retain the parent's scalar `rref_affine` as `rref_affine_reference` for controls only.
2. Implement `rref_affine` with the same first-nonzero pivot order and the same full-row Gauss–Jordan operations, but store the augmented matrix as `uint16` and evaluate field multiplication solely through the frozen `32 x 32` `ef.MUL` table. Pivot normalization is `MUL[row,INV[pivot]]`; elimination is exact XOR with `MUL[factor,pivot_row]`. No floating-point operation participates.
3. Build `E` and `ker(E)` once per public record. Pass that exact pair into P4, the public label scans, and all ten P1 branch checks. Rebuilding remains available only in a control. Cached values are immutable within a cell.
4. Preserve the parent's witness-bearing output schema, deterministic ordering, write-once/fsync/rename perimeter, source/control hashes, negative-cell semantics, and cell-level P1 aggregation.

## Mandatory controls before production

All controls run under the eventual release gate and before the first successor cell.

1. Re-run every parent semantic/counterfactual control under the fast kernel.
2. Compare scalar and vectorized complete RREF return tuples on at least 1,000 deterministic systems covering zero rows, zero matrices, rectangular matrices, rank deficiency, consistent affine systems, and inconsistent augmented pivots.
3. Compare `rank_of` and `nullspace_of` outputs exactly, including basis order, on the same systems.
4. Rebuild beta 8 from the primary binary-Goppa definition and compare the complete semantic payload to the hard-hashed parent cell after removing only provenance fields that must differ (`header`, `t_utc`). Every P1 witness coefficient/vector, accepted-label order, P2/P3/P4 field, and verdict must match exactly.
5. Rebuild beta 8 twice, once with per-cell caching and once with forced uncached construction; canonical semantic payloads must be byte-identical.
6. Plant a wrong multiplication table and require at least one deterministic RREF control mismatch.
7. Verify all five source hashes before import, normal Python assertions, niceness at least 10, and all five thread variables exactly one.

Any mismatch is **INSTRUMENT FAILURE**. No successor census number may then be reported.

## Domain, order, stopping, and evidence

The production population is exactly the parent's 24 monic degree-one Goppa polynomials `G_beta=Z+beta`, `beta=8,9,...,31`, over `F_32`, fixed ordered support `(0,...,7)`, fixed hold-outs `(0,1,2,3)`, and fixed tuple `(m,n,t,k,ell,n_ell,k_ell,D_ell,d,s,h)=(5,8,1,3,0,8,3,5,7,5,4)`. It evaluates every registered P1_U3/P1_U4/P2/P3/P4 predicate in that order and never stops on a mathematical failure. Cells are atomic and resumable only from exact source/control/header hashes. An internal two-hour wall cap preserves every completed byte.

Exact finite-field identities, exhaustive counts, and equality with the parent anchor are eligible for **MACHINE-VERIFIED**. Wall time and speedup are **COMPUTATIONAL-EVIDENCE** only and prove no mathematical claim. This toy domain is outside current 1786 Table 1 and Assumption 1's five Classic McEliece cells; even a 24/24 result is only a complete theorem for this one toy population. A failed predicate does not refute either paper.

## Release hold

Status: **PRE-REGISTERED / CODE NOT YET FROZEN / COMPUTE HELD BEHIND PARENT**. No command in this campaign may execute matrix construction, controls, or census until Main verifies the parent process has exited. Static source construction, source-text compilation, and hash/refusal checks are permitted while held.

## Static freeze while parent remains active

The five successor sources now compile from text. The unreleased authoritative
runner refused before source import, state creation, controls, or census.
No successor `state/` exists. Parent process ownership is unchanged.

Frozen source SHA-256 values:

- `toy_census.py`:
  `3fd89fd7d4af94a77d82835796c5d560a574cf366f85e5f80622e8a31c467436`;
- `toy_static_guards.py`:
  `9a433a2cfb75eeb2a09ba212a3afc8444ebb2d62d04e274aa83b3e6e432db673`;
- `run_toy_census.py`:
  `90ea138efa73df089b4cd9900b22b8cdc09a67e4ab4baf6a32db816b3a92737c`;
- `fastfield_frozen.py`:
  `209fe885f3d7837650e03fb63350d73a20492b96a17467969ac4ccf3ce62e560`;
- `gfield_frozen.py`:
  `69c6920e7775f80fd14d8d05166e5da6155c8fa778948f324d96bbaaef1fa6ae`.

The five-entry static-ledger SHA-256 is
`0da4a55800d75580bfb721cfcd30fc10ef63ba09a8bac4881af8de789196be4d`.
Status remains **COMPUTE HELD UNTIL PARENT TERMINATION**.

## Pre-release resume-perimeter refreeze

Static review found that the copied parent runner would rerun controls on every
resume and then refuse to replace the already-complete
`controls_result.json`. The successor now reuses a complete control artifact
only after matching its frozen source-hash map and every persistent plant
hash. A partial or mismatched control state remains preserved and aborts.
Cell, manifest, and checksum write-once semantics are unchanged.

Refrozen runner SHA-256:
`1faac0381530322659e2dba82cf9f35a1bc990046b0e9a0ced77c07539171d73`.
Superseding five-entry source-ledger SHA-256:
`368d1db28855f45581feec3423c06e8578e43833e8b3f6ef6175ba90ed5830ca`.
The runner compiles from text; no control or census arithmetic ran, and the
parent remains the sole active process.
