# Exact XOR-circuit headroom

## Claim and scope

Can an exact binary-XOR circuit compute the AES MixColumns 32-by-32 binary
linear map with at most 87 XOR gates? The strongest primary result located on
2026-09-06 is Jérémy Jean's **88-XOR construction**,
[ePrint 2026/1481](https://eprint.iacr.org/2026/1481), following the 89-XOR
construction of [Sun, Yang and Li](https://eprint.iacr.org/2025/1493).
A bounded search did not locate a subsequent smaller construction; that is not
exhaustive priority clearance. Jean explicitly credits Codex assistance.

The metric is unrestricted binary XOR count with free wires/fanout and free
renaming, not XOR depth, in-place operations, silicon area, or cryptographic
security. A smaller independently checked circuit would improve an externally
published construction. An algorithmic AI-superiority claim is a different task.

### Verified bit/byte interface

The two-page ePrint names `x0..x31` and `y0..y31` without an explicit
byte/bit-order definition. Our verified interface is
`x_(8*j+b) = coefficient of t^b in s[j,c]` for one AES column; likewise for
`y`. Thus bits are LSB-first within each byte, and the local integer places
byte `j` at bits `8*j..8*j+7`. This is an explicit interface, not attributed
author intent or FIPS bit-sequence numbering: FIPS sequence bit
`r_(8*j+7-b)` corresponds to this circuit's `x_(8*j+b)`.

[FIPS 197-upd1](https://doi.org/10.6028/NIST.FIPS.197-upd1) sections 3.2–3.4
and equations 4.1, 4.3, 5.7–5.8 define the byte polynomial, modulus and
MixColumns matrix. An independent carryless-polynomial/long-division replay
checked all 32 basis inputs under 16 normal/reversed-byte and LSB/MSB-order
input/output combinations; only the stated interface matched within that
census. The NIST Appendix B column `d4 bf 5d 30 -> 04 66 81 e5` also passed.
See the [exact convention receipt](../../data/frontier-scout-20260906/aes-convention-adjudication.json).
These are post-assay functional checks, not new candidate evaluations.

## Candidate gate

**Initial rank: A for a cheap assay, not permission for a large campaign.**
Importance is concrete: MixColumns is a standardized, widely used linear map.
Inputs are obtainable: the two-page primary contains all 88 assignments.
Verification is exact: compare the full linear map against the AES field
multiplication definition, not just the author's circuit or random vectors.
Headroom is **unknown**, not inferred from a loose global lower bound.

The first representation is deliberately modest: re-synthesize short windows
of the current circuit over the exact GF(2) forms available at their boundaries.
Classical exact search, not an LLM, is the first control. This directly tests
whether that cheap neighborhood has already been exhausted. Failure prevents
wasting program-search calls there; it is not a new global optimality theorem.

## Frozen first assay

Before evaluation, `scratch/headroom-preregistration.json` pins the exact primary
PDF, transcribed 88-gate input, source hashes, and all 258 contiguous windows of
widths 2, 3 and 4 in published order. Boundary inputs include every earlier wire;
required outputs include all wires consumed later and every final output.
Every replacement uses at most one fewer gate. Exact search covers at most
three binary XORs, with a complete finite lower-bound/search argument.

The independent checker uses the AES polynomial x^8+x^4+x^3+x+1 and the
MixColumns coefficient matrix. It checks all 32 basis vectors and the FIPS
`db 13 53 45 -> 8e 4d a1 bc` example, and rejects a deliberately corrupted
circuit. No floating point, private training instances, or source-text-based
functional claims are used.

Limits: one CPU core, nice 10, 60 CPU seconds, 90 wall seconds, no automatic
retries, no extension of windows after seeing outcomes. A timeout is
inconclusive; no partially explored neighborhood is called optimal. Primary
metric is integer gate savings. Global optimum is unknown, so global gap
closure is null; a zero denominator is never called 100% closure.

Success means a strictly smaller circuit whose full map independently checks.
Negative means no saving in the completely explored stated neighborhood only.
Both outcomes stop this assay. A positive result needs refreshed priority
comparison before a record claim; an algorithmic follow-on must freeze public
matrices, development/holdout split, candidate-evaluation and compute budgets,
and compare against contemporary non-AI local optimization and LCB-BP controls.

## Failure and publication boundary

The public circuit may already be irreducible under these local substitutions.
A breakthrough may require nonlocal sharing or a different circuit ordering.
The exact verifier does not establish an efficient search strategy. A local
irreducibility result alone is unlikely to merit a research paper. A genuinely
smaller AES circuit, a substantive new lower bound, or a transferable method
beating current public comparators is the publication-level target.

Source acquisition and selection evidence:
[`data/frontier-scout-20260906`](../../data/frontier-scout-20260906/).
Campaigns must be minted and frozen with `../../scripts/campaign.py`.

## Current state (agent Main, 2026-09-06)

**C: do not extend this local-search lane.** The
[`first assay`](campaigns/20260906T021828Z_57f63dca_6cc09db616c6/)
closed **FROZEN-INCONCLUSIVE**, not certified irreducibility.
The 88-gate public circuit reproduced exactly; the normal-form search found
no saving in all 258 registered windows (0.015378 CPU seconds).
Independent byte-field replay passed every basis vector and the worked byte example;
the deliberately corrupted operand was rejected.

The independent search checked or lower-bounded 239 windows, but **19 Z3
queries returned UNKNOWN** at the frozen 500 ms limit (14.424723 CPU seconds
for verification). No timeout was relabeled UNSAT, no limit was expanded,
and no 87-gate circuit, global optimum, AI advantage or completed local
irreducibility certificate is claimed. Move to the next candidate.

For fresh replay from the workspace root, use the parent-supervised command:

```sh
physics/.venv/bin/python -I -B cs/xor/replay.py \
  cs/xor/campaigns/20260906T021828Z_57f63dca_6cc09db616c6 \
  cs/xor/scratch/NEW_REPLAY_DIRECTORY
```

The output path must not exist. The wrapper runs the unchanged archived
scripts with one thread and their original 60-CPU/90-wall-second limits;
an outer 95-second per-stage deadline owns interruption receipts.
SIGALRM, SIGXCPU, nonzero exits and the outer deadline produce explicit
INCONCLUSIVE receipts, never UNSAT or success. Partial receipts are retained.
The verifier requires Z3 5.1.0 in `physics/.venv`; exit 2 means inconclusive.
Do not run the raw scripts against the frozen directory. The archived
`known_fips_vector` field names the `db 13 53 45` worked example; the explicitly
sourced NIST Appendix B check is in the post-assay convention receipt above.
