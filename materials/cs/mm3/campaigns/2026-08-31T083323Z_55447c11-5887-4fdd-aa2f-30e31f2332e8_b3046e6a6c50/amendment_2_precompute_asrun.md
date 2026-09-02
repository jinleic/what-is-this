# Amendment 2 — precompute as-run provenance seal

**Additive amendment committed before target compute. It changes no domain, action, order, cap, arithmetic, or adjudication rule. Source, pre-statement, and amendment 1 are untouched. Status: COMPUTE PENDING.**

Amendment 1 deferred byte-freezing until after the run. This amendment corrects that plan: the exact source and every local code/data dependency were byte-copied into `asrun/` before compute. The source hash remains `b3046e6a6c506ffd5c78871a9f057770a552968314ac6fd225ad60e2b1d4e77c`, matching the campaign suffix.

Because the exact source contains absolute live-tree import paths, `asrun/frozen_launcher.py.asrun` preloads every canonical local module name from the campaign copies, executes the campaign source copy with `runpy`, and patches only its `SRC`/`SCRATCH` data roots to the campaign copies before calling the unchanged `run_campaign`. `CAMPAIGNS` remains the live target campaign parent so the source's path guard still authenticates this campaign. No local dependency is imported from the live tree by the frozen run.

Frozen command: `asrun/frozen_command.sh.asrun`. It invokes `nice -n 10` and pins OMP, OpenBLAS, MKL, vecLib, and NumExpr to one thread.

| as-run dependency | SHA-256 |
|---|---|
| `asrun/frozen_command.sh.asrun` | `9d60d64c8443c4b36e2c2c1eea27396ce7db7b90cf1edddbdb9b70ba69388648` |
| `asrun/frozen_launcher.py.asrun` | `852798a4a1656a1bb4da39e40374242e5e30fb8fc8dfeb08995884a201f8eab3` |
| `asrun/scratch/cr58_cn122_ZT_reduced.json` | `2821574f8182f0abc97f1b31168e5077a3d6c63d756dbe1f06f24075b0d9926e` |
| `asrun/scratch/mws59_layout.txt` | `7bc075d41f6cd2ef71a7c5d5d76d6b3d669e0845728535b857c15ef1d735208f` |
| `asrun/scratch/sun56_verify.py.asrun` | `3b7fb40aae051fcc2236dd151ad60f802fb68244df296686667d9a313fb18033` |
| `asrun/src/gate_a.py.asrun` | `1ddb004b342c64ddcfefea9ee7b66303fbe07f22951303c485aaf559f306a260` |
| `asrun/src/gate_b_floor.py.asrun` | `c8caffa8fb48de6d3542a5439f84c35e4d74f1b1cab6710a154842ce20002c1a` |
| `asrun/src/gatec_decomps.py.asrun` | `e96e4f52053be7f02e23b3fd694c68845711da2c429ba887a8b2daa8e8229eed` |
| `asrun/src/gatec_sweep.py.asrun` | `0dcb5844f3065fa445d0f079d8e9055cf6f9d8158d9f78fb72bca50f7e4b4ed4` |
| `asrun/src/new_decomp_offdiag.py.asrun` | `b3046e6a6c506ffd5c78871a9f057770a552968314ac6fd225ad60e2b1d4e77c` |
| `asrun/src/stapleton60_data.py.asrun` | `70da1bfb64f8c51637e540eb19873449c712cb80b9d5f272de10020cbe5ee563` |
| `asrun/src/tensor_data.py.asrun` | `b55b18a8eb95fb950973b944a2e8c703bbba78d020a83a6ef3a32c36f9f8a793` |
| `asrun/src/verify_anchors.py.asrun` | `0ffcd100a4c08200f8c2d3aec04cf0a5e53b16c35daf99f596ac64daa91b3239` |

`checksums_precompute.sha256` seals the complete precompute campaign contents. It remains immutable after compute; final artifacts will receive a separate final `checksums.sha256`. Staging performed byte reads/copies/hashes only: no Python target import, anchor, control, census, test, solver, or checksum-verification run.

## Owner release-perimeter correction

The first frozen launcher relied on the external checksum ledger and only
refused a completed `campaign_result.json`. It did not itself verify the
campaign-local dependencies, gate release, assert resources, or protect the
other result paths from a partial-run overwrite. Its launcher and command
hashes above are revoked before target compute.

The accepted launcher now hard-verifies all eleven campaign-local code/data
dependencies before import, requires `MM3_CAMPAIGN_RELEASED=1`, `__debug__`,
nice at least 10, and all five thread variables equal to one. It refuses if
any reserved result path, final checksum, run log, or campaign-level `.tmp`
already exists. Therefore a partial run is preserved and cannot be silently
restarted over. The frozen command sets the release gate and disables bytecode
writes. `--check-only` exits after this perimeter without importing a target
module.

Static perimeter evidence: an unreleased direct launcher invocation refused
before module import; the released, resource-pinned `--check-only` invocation
verified every frozen dependency and every reserved output, then exited before
target import or compute.

- `asrun/frozen_command.sh.asrun`:
  `55c23029adc490179f934fa87b6225127457901f5c1e40e2699d8f7b7ca952cf`
- `asrun/frozen_launcher.py.asrun`:
  `f94f57c5e723923abd6dcac28b07e758d19b8de817805d6134cbc6381a6ad14d`

The mathematical source remains byte-identical at `b3046e6a6c50…`; domain,
order, arithmetic, and adjudication are unchanged. Status remains
**COMPUTE PENDING**.

This owner correction supersedes line 27's earlier “immutable” seal: the
precompute checksum ledger was refrozen once, before target compute, to bind
the hardened launcher, command, manifest, and this additive correction. It is
immutable from this corrected freeze forward.
