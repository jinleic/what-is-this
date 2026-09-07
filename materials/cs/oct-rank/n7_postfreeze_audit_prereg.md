# N7 postfreeze premise audit

Gate: `n7-postfreeze-premise-audit`. This is an independent verification of a
known claim, not a blinded prediction or a new search population.

## Claim and scope

The legacy N7 record claims exact rank 156 for the 192-by-247 complex
Jacobian at the twelve-term N3 witness with a zero thirteenth slot. Its
field-independent tangent-containment proof then bounds every single-term
split by 156 + 10 = 166 < 192. This excludes full-rank square charts in that
split family; it does not determine real tensor rank or exclude non-gauge
multi-term deformations.

The legacy N7 directory remains untouched. Independent review found missing
producer freeze/close artifacts, stale pilot-hash attribution, unpublished
independent-check evidence, a merely checkpointed RSS limit, and a prereg
hash that was recorded but not enforced. This new producer-generated run
will preserve these disclosures rather than retroactively repair history.

## Fixed input and independent method

Read only the twelve terms from the frozen N3 witness
`campaigns/20260904T034332Z_f5a61843_0333458ed434/witness_tf_c12_qi.json`,
SHA-256 `251e00577b43b400aa97d306bb923a98eb79ca7850ac67652d45288c2e570037`.
Decode `a`, `b_orig`, and `c_orig` as exact pairs of rationals; add a zero
thirteenth slot. No N6/N7 code may be imported.

1. Enforce input, root-prereg and run-copy hashes before matrix work.
2. Rebuild the two quaternion left-multiplication blocks independently;
   check all 192 tensor entries exactly and reject a one-entry corruption.
3. Check the realification sign convention on the complex rank-one matrix
   `[[1,i],[i,-1]]` (real rank 2); changing its final entry to zero must
   give real rank 4.
4. Construct each CP derivative directly in row order `(p,j,k)` and column
   order A-major, B-major, C-major, with term index innermost.
5. For `J=X+iY`, compute exact rational rank of `[[X,-Y],[Y,X]]` using
   python-flint 0.9.0. Hash the canonical rational matrix row serialization.
   An odd rank is an instrument failure. Divide an even rank by two.
6. Persist exact ranks, all control outcomes, input/driver/prereg hashes,
   runtime and the implication bound. PASS requires rank 312 (complex156).
   Any mismatch is an explicit audit failure, never silently reconciled.

## Lifecycle and local resources

Commit this statement, initialize with `scripts/campaign.py`, copy it
byte-identically, bind the standalone driver before launch, and obtain an
independent source review. Use the owner resource guard: at most two
single-thread jobs, sampled total-CPU pause at40%/resume below30%, 100GiB disk
reserve, sampled2GiB group-RSS stop, 64MiB per-file ceiling. This audit adds
at most8MiB and has a180-second outer wall ceiling. No scheduler bypass,
package installation, bulk copy, or deletion is authorized. Resource stops
produce an inconclusive run with partial evidence, not an audit PASS.

After a complete independent audit, freeze through the producer and close
exactly once: `FROZEN-CERTIFIED` for PASS, `FROZEN-NEGATIVE` for an exact
contradiction, `FROZEN-INCONCLUSIVE` for incomplete/resource-stopped work.
Only then append authoritative results and refresh generated state. A PASS
supports only the stated derivative/split-family theorem, not a new tensor
rank bound. No frozen predecessor bytes are modified.
