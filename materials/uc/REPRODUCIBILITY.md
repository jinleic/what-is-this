# Campaign I reproducibility and clean-room replay

## Status vocabulary and scope

- **MACHINE-VERIFIED** means that the named executable check completed successfully on the exact bytes identified below.
- **HUMAN-AUDITED** means that the mathematical argument is meant to be checked from the displayed derivation; it is not established by a computation in this directory.
- **COMPUTATIONAL-EVIDENCE** means that an archived or newly measured computation is useful evidence but is not, by itself, an accepted certificate replay.
- **CITED-DEPENDENCY** means that the step comes from the cited source and is not reproved here.
- **OPEN** means that the check has not been completed by the clean-room harness described here.
- **FAILED** means that a check rejected its inputs, environment, execution, or output. A failed or interrupted slice contributes no certificate claim.

**MACHINE-VERIFIED (identity and structure):** `uc/verification/structural_trace_audit.py` independently streamed and hashed all 488,465,854 trace bytes, accepted only documented opcodes, verified DFS topology/termination, and reproduced all committed structural tallies. The hardened repository-local report is `verification/results/structural-traces.json`, canonical SHA-256 `7270e5418be72b09e607457438bc9f8109154c31b491abca9546f4ebee3cf5a3`.

**MACHINE-VERIFIED (complete hardened clean-room replay):** all eight fresh read-only staged slice reports passed external pins, environment/module hashes, exact stdout/tallies, and before/after source/stage checks. `aggregate_replay_reports.py` accepted 488,465,854 total nodes and zero residuals. Composite: `verification/results/clean-room-arithmetic/composite.json`, canonical SHA-256 `57ca3f52c054bd9bccfefc349f440674c74fe7ef671d06c103ef3f52dc33d53d`.

**MACHINE-VERIFIED (complete fresh direct replay):** eight isolated `python -I -B` processes replayed every trace from the externally pinned raw launch/result/trace inputs. All exited 0; all stdout hashes and tallies matched; residual was zero; total processed was 488,465,854. Independent structural audits before and after proved the input maps unchanged. Report: `verification/results/direct-isolated-replay/summary.json`, canonical SHA-256 `6ef126ced337b035e5c5f22a130b1be3697666e60f16586b75540ef6fc734b16`.

**FAILED (superseded launch):** an earlier eight-process clean-room batch was stopped after roughly four minutes with exit 143 because the verifier was hardened during execution; its before/after verifier hashes could not match. It produced only `START` events and no accepted `latest.json`. The final batch was restarted from the hardened bytes.

**COMPUTATIONAL-EVIDENCE (archived, not accepted as fresh):** the campaign-local collector/replay marker files are excluded from staging and were never used by either fresh replay.

**MACHINE-VERIFIED (separate bridge check):** `uc/verification/cambie_bridge_strictness.py` independently checks the self-contained point-mass value at 256 Arb bits. This is separate from the Campaign I certificate replay; see [Cambie bridge strictness](#cambie-bridge-strictness).

## Frozen target and external pins

**MACHINE-VERIFIED:** the replay target is only:

```text
uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8
```

**MACHINE-VERIFIED:** `uc/verification/campaign-lock.json` is the checked external pin manifest. Both verifiers hard-code its raw SHA-256 and reject a replacement, so Campaign I identity is not inferred merely from self-consistent values supplied by `launch.json` or the result records.

| Pinned object | Exact value |
|---|---|
| External pin manifest raw SHA-256 | `d254a23a4df7cb3a4c6aae2453883324fba2381c1890c4f011d61c7971b84b48` |
| External pin manifest size | `12590` bytes |
| Campaign ID | `20260818T212601Z_425f109c15b64a6198785c6cebbbdaab` |
| Raw `launch.json` SHA-256 | `63bb6fbe1cfbcb6b7d4c50d7b66882255e5e778ce604c92ecd1abdfa1ea8337e` |
| Raw `launch.json` size | `4339` bytes |
| Canonical launch object SHA-256 | `6414affcaa54bd5e37e9f7bd351c6ed4cf675167b10b55d86cce1436ebdf65cc` |
| Combined executable `code_sha256` | `2f23a58ebdb8b14275284e2ffefca7862373f7baf5d354a06e607cff8d78ce05` |
| Frozen collector SHA-256 | `95d09322a7821f41850ef1ce77345e5e93eb4fa570cc0e7f58f817716cf681ec` |
| Frozen replayer SHA-256 | `4bfc2568e6d5b5364516f1acf734eb049d0f670bbec6cc6c95728466146ebeef` |
| Trace schema | `cert3-trace-v1` |
| Result schema | `cert3-result-v2` |
| Campaign schema | `cert3-campaign-v2` |

**MACHINE-VERIFIED:** the exact 15-file snapshot inventory and raw source hashes are:

| Snapshot file | SHA-256 | Bytes |
|---|---|---:|
| `arbcore.py` | `42f7f7afdd69984a25ca3e6b1db43b6b3792db43c925f4e5eafe929fdb05e5bd` | 4,013 |
| `bound_kkt.py` | `ece39bd496a7fe88bd63d1c337cf242fcce6264c28bc94132ff266e94f35b4f4` | 33,053 |
| `bridge_uc.py` | `49af3194090798a55609fe17b926ca95ed552c87adb06ea55f8e7a553d61cd54` | 9,999 |
| `cert2.py` | `b6cc756d07ed360af4eb3913cc899cb73c3288af72ae07345cf81bfe20e6c1fa` | 26,364 |
| `cert3.py` | `984a568fbd70ce4860f981e4f7e105b7028ec3c3869ecc7560b8b4ac08a06b3f` | 70,700 |
| `cert3_collect.py` | `95d09322a7821f41850ef1ce77345e5e93eb4fa570cc0e7f58f817716cf681ec` | 24,056 |
| `cert3_par.py` | `a741aa167a79fb1058d89373db15e19d4124576a7eb7fb76c57ee3f70d3709bc` | 29,848 |
| `cert3_replay.py` | `4bfc2568e6d5b5364516f1acf734eb049d0f670bbec6cc6c95728466146ebeef` | 14,770 |
| `decomposition.py` | `622d297cc61e0cf9adb67b78c09afaac1b32d15deb0661ee9cee7bb441735515` | 16,219 |
| `diag_exhaust.py` | `23d8ec78fc2bc1ac29262577705551eb25a4a98c24ffb6062d656f39a7e00016` | 24,222 |
| `entropy.py` | `2fbddcd8b2458d652ffd0d4268f9c868353967c6aacd7a66183556c65ece88a5` | 2,459 |
| `lemma_rh_proof.py` | `4236dcfd159ea091638d4d7e3fe4f1a804c8d551c3490a7e311f56f6be7a9fbd` | 14,602 |
| `margin_lemma.py` | `d11d2eb57269abe849acf04ce8acb4a7d908b2d210370cb48a2140cdcd895ace` | 18,891 |
| `reduction.py` | `b2af27e52b28339258f384f98ffd451b54b4a2cef2aae995fa0993880cfecc65` | 23,274 |
| `thmB3_proof.py` | `c459d70d8e06a855031c52fb76f60db3510c099bfa55f6901fc97b594e5a214f` | 24,473 |

**MACHINE-VERIFIED:** an extra or missing snapshot entry, including `__pycache__`, a `.pyc`, a native module, or a same-named planted dependency, is a hard failure. Extra or missing `result_*.json` or `trace_*.bin` artifacts are also a hard failure. Other campaign files are classified as non-inputs, listed in each report, and never copied into the stage.

## Environment contract and capture

**MACHINE-VERIFIED:** the launch records, and the clean-room harness requires, this exact version/platform tuple:

| Component | Required value |
|---|---|
| Python | `3.14.3` |
| Implementation | `CPython` |
| Platform string | `macOS-26.5.2-arm64-arm-64bit-Mach-O` |
| python-flint | `0.9.0` |
| mpmath | `1.3.0` |
| NumPy | `2.5.2` |
| SciPy | `1.18.0` |
| SymPy | `1.14.0` |

**MACHINE-VERIFIED:** the harness now requires these exact top-level dependency module-file hashes. They narrow accidental environment drift but do not pin transitive Python modules or native-library bytes.

| Module file | Observed SHA-256 |
|---|---|
| `flint/__init__.py` | `2e5f8f1768d14eccd7961353c635195bd557f297edff8b8de09e2d211f03ec2d` |
| `mpmath/__init__.py` | `b241584d2c1fc0304b0a1015ea923749d7b0800411dd406dcab7c82bf25d9fe8` |
| `numpy/__init__.py` | `09295a80660f17925ae23765ce8cbd7ff7ceae968d5f2f89349f1cb74c0b9e11` |
| `scipy/__init__.py` | `39ccae300a4739cc53719bd3b39ce3f8f66e736c6eb672006e338a4c5cc3dc76` |
| `sympy/__init__.py` | `4e9476348ba105feab28d82f5bcf6cdba2e3e84de6e059bbfe7a13728c0a4ab0` |

**MACHINE-VERIFIED with explicit platform trust:** the harness requires `-I -B`; the child uses `-I -B -s -E`, an empty working directory and home, no `PYTHONPATH` or user site, one-thread numerical settings, exact versions, and the five module hashes above. Python still processes the system site because `-S` is not used; system-site startup behavior and transitive/native dependency bytes remain explicit trust premises rather than a sealed-environment guarantee.

## What the clean-room harness accepts

**MACHINE-VERIFIED (preflight behavior):** `uc/verification/clean_room_replay.py` performs these steps for one slice:

1. It verifies the hard-coded raw hash and size of `campaign-lock.json`, then verifies the independent Campaign I code, launch, collector, snapshot, result, and trace pins in that manifest.
2. It rejects symlink/non-regular required files; verifies the exact snapshot, result, and trace inventories; validates canonical launch/result/worker object hashes; and checks the exact rational partition, completion fields, zero residual/budget/stack state, and committed tallies.
3. It hashes every required source input before staging. A slice input set consists of raw `launch.json`, all 15 frozen source files, that slice's raw result record, and that slice's raw trace.
4. It creates a fresh temporary copy containing only those inputs, makes its files `0444` and directories `0555`, and hashes the staged copy before execution.
5. It invokes the frozen `snapshot/cert3_replay.py` directly from raw staged `launch.json`, result JSON, and trace bytes. No existing replay marker or collector output is an argument or staged file.
6. It accepts exactly exit code 0, empty stderr, and exactly one JSON stdout line whose 11 replay tallies equal the externally pinned expected values.
7. It re-lists and re-hashes the stage, re-hashes the original source inputs, and re-captures the environment. Any before/after difference makes the report `FAILED`.
8. It atomically writes a per-slice report and appends `START`/`FINISH` events to `slice-N/events.jsonl`. A report includes source/stage before/after hashes, dependency capture, actual elapsed time, process termination, exit code, stdout/stderr text and SHA-256, resource settings, result tallies, and a canonical `report_sha256` excluding that field.

**MACHINE-VERIFIED (cache exclusion):** `collect_output.txt`, `replay_slice*.txt`, `.replay_slice*_OK`, `.worker_slice*.json`, bytecode, and arbitrary cache files are never evidence for a fresh replay. The harness does not parse them. This is intentional even when their content agrees with expected tallies.

**FAILED-CLOSED resume policy:** repository reports have only public unkeyed hashes, so they cannot authenticate that a replay actually ran. `--resume` is disabled; every proof-acceptance invocation performs a fresh replay. A protected external attestation system may add resumability outside this repository.

## Exact low-CPU commands

Run from the `math` directory. `OUT` must be outside the frozen campaign; keep it on durable storage for logs and reports.

```bash
PY="$PWD/.venv/bin/python"
CAMPAIGN="$PWD/uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8"
OUT="/path/to/durable/uc-campaign-i-replay"
STAGE="${TMPDIR:-/tmp}/uc-campaign-i-stage"
mkdir -p "$OUT" "$STAGE"
```

### Fast preflight without arithmetic

**OPEN:** `PREFLIGHT_PASS` establishes only identity, staging, environment, and before/after hashing. It does not establish a discharge rule or a certificate.

```bash
"$PY" -I -B uc/verification/clean_room_replay.py \
  --campaign "$CAMPAIGN" \
  --slice 0 \
  --preflight-only \
  --output-dir "$OUT" \
  --staging-root "$STAGE" \
  --nice 10
```

### Independent structural audit of all traces

**MACHINE-VERIFIED:** this command is single-process, single-core, nice-adjusted, and does not import the proof arithmetic.

```bash
"$PY" -I -B uc/verification/structural_trace_audit.py \
  --campaign "$CAMPAIGN" \
  --output "$OUT/structural-traces.json" \
  --nice 10 \
  --cpu-limit-seconds 600
```

**MACHINE-VERIFIED:** the hardened all-slice run completed in 54.46924633299932 measured seconds, with report object SHA-256 `7270e5418be72b09e607457438bc9f8109154c31b491abca9546f4ebee3cf5a3`. This is evidence for that run only.

### Sequential fresh arithmetic replay

This sequence launches no parallel slices. Every invocation performs fresh arithmetic; the harness forces common numerical thread pools to one thread and applies `nice +10`.

```bash
for SLICE in 0 1 2 3 4 5 6 7; do
  "$PY" -I -B uc/verification/clean_room_replay.py \
    --campaign "$CAMPAIGN" \
    --slice "$SLICE" \
    --output-dir "$OUT" \
    --staging-root "$STAGE" \
    --nice 10 \
    --cpu-limit-seconds 43200 \
    --wall-timeout-seconds 43200
done
```

After all eight fresh reports exist, validate and aggregate them:

```bash
"$PY" -I -B uc/verification/aggregate_replay_reports.py \
  --reports "$OUT" \
  --output "$OUT/composite.json"
```

The aggregator recomputes every report object hash, external input map,
environment/tool match, stream hash, parsed stdout tally, implementation hash,
slice identity, and the 488,465,854-node total.  Its output remains a public
self-consistency record, not an authenticated execution attestation.

**COMPUTATIONAL-EVIDENCE:** `43,200` seconds is an explicit operational cap, not a claimed runtime. It is larger than the slowest archived replay slice (`25,329.8` seconds) but can still reject a slower machine. Set either limit to `0` to disable that limit; doing so does not enable parallelism or change the one-thread settings.

**MACHINE-VERIFIED (resource behavior):** CPU-limit expiry, wall timeout, signal termination, nonzero exit, stderr, missing stdout, malformed JSON, or a tally mismatch produces a nonzero harness exit and a `FAILED` report. A two-second bounded smoke exercise reached the isolated frozen replayer, terminated it deliberately, re-hashed source and stage inputs, and correctly reported `wall-timeout`/exit `-15` as `FAILED`; this is a failure-path test, not a replay result.

## Expected slice inputs and output hashes

**MACHINE-VERIFIED (byte identity and structure):** the independent structural run re-hashed every trace and raw result below and reproduced the expected processed/split topology. `Raw result SHA-256` is the byte hash of the JSON file; `record_sha256` is the canonical JSON object hash excluding its own field. `Expected stdout SHA-256` hashes the exact successful frozen-replayer JSON line including its final newline. Successful stderr must be empty, whose SHA-256 is always `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.

| Slice | Exact root in `w` | Raw result SHA-256 | `record_sha256` | Trace bytes | Trace SHA-256 | Expected stdout SHA-256 |
|---:|---|---|---|---:|---|---|
| 0 | `[1/2,9/16]` | `1563f7310eefc2677b53d74448e31fc01c76a5d558ce545a2fb101ac1c043e30` | `8a1bbb4d960b43128fd5884eec2054bbbbc267cba9b940a21b8d97ae32010133` | 21,135,611 | `4f55bcc261c541ab1adadff85f2b1148ce87fe5d34a9782f5cc14d1ba6588240` | `b0d61e1336654f973610485f44eba24c8f53a68c13ad0cc01345f8d4239e65b0` |
| 1 | `[9/16,5/8]` | `2692a0138d6a27806a34978d42a3073e073635239abe99639610a92427832b0a` | `833ba059a088bfe67e2b1125485f727de21824623a0d307fadd166fbd7d95a89` | 20,827,065 | `4f267f19cd0a99e7ed9dbb423b35fd1b4f2c9bf86f170e5527b2c31dd3870adb` | `92c54a2e7ab1174c11f9f4243ab55e3392798f064fd8d4b4a7cdd6aa0c3312a6` |
| 2 | `[5/8,11/16]` | `98575a068643f3d31e64876f7f3f4dcae7c45e5fcacd323551d3661c0669045d` | `cde3056e999d2eed17b1da7acfb4930bde7531ddd7bec7051acd237b40d27658` | 23,486,265 | `e31cd8fdcec4db33a2eac19a153ef6d1b70772bb53416c7255f4b9ae088ec0ae` | `9690ebe1efce97db953d29437c29a07a33c9052d0caf12e344cfd9372d675141` |
| 3 | `[11/16,3/4]` | `6ed4deef8cf9305e86db8529b2446ac9ed884ecba746a10afa24de236f90f493` | `41d1b64aacfe84fc8cd087dea77494fc3db3a052421f10c89d179d8dc87622f2` | 34,186,855 | `d6cfc79989e17d69bfd449ec2a0dbcbf75f2f859a531a19114a11748a3f357e5` | `6c0c74133fb96ed64c6645457ef30caed7fde0e1a5060127d1021bd7cbdad43b` |
| 4 | `[3/4,13/16]` | `572d73f4ffef9d8df556bee2784d71f6075e2179f49499e04ca581e69d3bba12` | `b29baa18da4cc596bbc2b176211bae9e251e8589c4450ead694215ecb1b2b11e` | 71,405,987 | `c1434d19d3a7ce44224c295fc0eef0e56a14fd99351d9604b4edc84b0ff44f70` | `2d1f31a8950d15a12c5507021ed1f496503c15651ad20dfeea9f8e6d8ff62e52` |
| 5 | `[13/16,7/8]` | `9e31c7f73486874893722c946837f2e9bce1c8db53483249fc236860d32ca759` | `c725d9b10bff5a043fa0bd99c5e3647b0cc3f9bae85ea7c16052645f79a13a23` | 117,709,435 | `ae584e405a97229755ddb81dc77cfb5113938854e7566a6226c5ea5dee392f3c` | `0c7d9aaa36ead48d7135d0d07e6db4e2b56cf08726a1cfa52bec0f09fd8ca7e6` |
| 6 | `[7/8,15/16]` | `894c2c1a412395ca4feffcb52d6473f47ed12fbb5dd5cd635968233e16b7f1dd` | `4bdc2d23b1ccfea0dd0f5f9e790fe3b3f0a95c34f54c1cc67fd43febc00fa4cd` | 121,456,775 | `4a483eda93ebfbacc033f95f363b55af795c6e8664f67504849b943f056d03ef` | `cf459a39ab8c13cd6f2b80df203179c8e6067b71589ab3562ed155852879bf1e` |
| 7 | `[15/16,1]` | `ed01f108e7642a4c85c3b09d93840a46cc51cf83f04596e0fdda5aefe325c603` | `1f72efcd57f1a28ea15891ab41674183d16ce92d4ee112c57c5c3c82d4321f32` | 78,257,861 | `302fc8b877bf2331c10b7d67a717ca18a513f748cf9cf6929da8db6108bf5156` | `f4407f29f33837d3b04767d80895514ea7b9debfb27c274f44b2abfc53b71ea3` |

**MACHINE-VERIFIED (structural tallies):** these are exactly the 11 tallies encoded by trace opcodes and expected from a successful frozen replay.

| Slice | Processed | Split | Infeasible | Corner | Ratio | Center | Center mixed | Center mixed swap | Center `w` | Face | Residual |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 21,135,611 | 10,567,805 | 766,564 | 3,218,757 | 356,182 | 4,166,369 | 1,059,513 | 831,310 | 27,088 | 142,023 | 0 |
| 1 | 20,827,065 | 10,413,532 | 835,153 | 2,837,057 | 446,320 | 4,074,579 | 1,107,712 | 858,664 | 23,434 | 230,614 | 0 |
| 2 | 23,486,265 | 11,743,132 | 1,159,349 | 3,338,422 | 694,759 | 4,198,151 | 1,025,766 | 827,217 | 29,236 | 470,233 | 0 |
| 3 | 34,186,855 | 17,093,427 | 2,637,747 | 4,699,405 | 777,097 | 5,261,749 | 999,230 | 592,501 | 14,975 | 2,110,724 | 0 |
| 4 | 71,405,987 | 35,702,993 | 6,969,122 | 9,899,146 | 84,904 | 11,895,793 | 1,269,484 | 95,757 | 1,427 | 5,487,361 | 0 |
| 5 | 117,709,435 | 58,854,717 | 12,719,058 | 12,155,506 | 27,261 | 22,256,705 | 1,786,704 | 26,246 | 531 | 9,882,707 | 0 |
| 6 | 121,456,775 | 60,728,387 | 13,349,700 | 13,820,378 | 23,348 | 20,315,341 | 2,937,116 | 16,355 | 313 | 10,265,837 | 0 |
| 7 | 78,257,861 | 39,128,930 | 5,057,436 | 16,794,114 | 13,108 | 11,466,023 | 4,150,412 | 2,473 | 80 | 1,645,285 | 0 |

**MACHINE-VERIFIED (topology):** for every slice, terminal events equal split events plus one, no event occurs after the DFS pending count reaches zero, and EOF leaves zero pending nodes. Maximum pending-node counts were `31, 35, 36, 38, 43, 42, 42, 38` for slices 0 through 7. Opcode 9 (`residual`) occurs zero times.

## What is and is not independent

| Check | Status | Arithmetic independence |
|---|---|---|
| External pin manifest, regular-file/inventory checks, raw hashes, canonical JSON hashes, exact rational roots, before/after comparison | **MACHINE-VERIFIED** | Standard-library implementation independent of campaign code |
| One-byte opcode inventory, split-coordinate labels, DFS pending-node topology, termination, and trace-derived tallies | **MACHINE-VERIFIED** | `structural_trace_audit.py` imports no campaign or numerical module |
| A terminal opcode actually clears its reconstructed box | **MACHINE-VERIFIED** by the two same-source replays and the complete secure byte-zero replay | Secure replay uses an independent Python formula/control implementation and interval AD; it shares the trace and Arb primitives |
| Frozen `cert3_collect.py` provenance checks | **COMPUTATIONAL-EVIDENCE** in archived output only | Its arithmetic phase calls the same frozen `cert3_replay.replay_trace`; it is not a second arithmetic implementation |
| Structural trace audit | **MACHINE-VERIFIED** | Structural only; it cannot validate interval contraction, corner, ratio, centered, pin, or face predicates |
| Cambie point-mass strictness value | **MACHINE-VERIFIED** | Self-contained formula, independent of Campaign I source modules; still depends on python-flint/Arb |
| Entropy-to-union-closed implication | **HUMAN-AUDITED** | Reconstructed in `uc/AUDIT.md` and the revised paper; Cambie is the source |

**HUMAN-AUDITED:** exact dyadic roots matter to the structural interpretation: the first four coordinates are `[0,1]`, and the fifth is the pinned sixteenth-width interval for the slice. Each split opcode 16 through 20 selects one coordinate; exact midpoint bisection preserves dyadic rational endpoints. The structural auditor verifies the root partition and split opcode range but deliberately does not reproduce floating contractors or interval discharge arithmetic.

## Standalone byte-zero arithmetic implementation

The new
[`independent_arithmetic_replay_secure.py`](verification/independent_arithmetic_replay_secure.py)
imports no frozen campaign arithmetic.  It independently transcribes the
entropy, \(s^\star\), \(\operatorname{rh}\), and \(\rho\) interval primitives;
recomputes the two 20,000-cell global covers; derives centered derivatives
through interval automatic differentiation rather than copied expanded
formulas; reconstructs the exact decimal/Fraction contractor and DFS; and
rediscovers every nested face proof because opcode 8 stores no inner face
trace.  It still trusts the authenticated outer trace partition and the same
python-flint/Arb primitive library.

The exact 80-bit contractor threshold omitted from the historical trace is
pinned and asserted as
\[
\frac{115472366449356389981023}{2^{78}},
\]
the upper endpoint of the Arb ball for `0.3820660112501052`.  Center rules use
the two authenticated nonnegative shifts `0.0` and
`1.3682158100251758`; face rules enumerate all seven launch values.

An adversarial review found that the first standalone implementation's public
checkpoint HMAC was only a corruption checksum: a forged EOF checkpoint could
have supplied public expected tallies.  Its direct byte-zero launch remains
diagnostic evidence only.  The secure implementation closes this acceptance
path:

- checkpoint resume is allowed only for finite `PARTIAL_AUDIT` work and can
  never produce `PASS`;
- every accepted report states `initial_trace_offset: 0` and
  `resume_used: false`;
- lock, launch, result, trace, verifier, and 343-file runtime-seal bytes are
  checked before and after arithmetic;
- the executable aggregator pins the exact lock, campaign, verifier, runtime,
  target, precision, eight adjacent roots, offsets, tallies, and topology;
- aggregation requires an out-of-band raw SHA-256 for a report lock that fixes
  the eight raw reports, launch manifest, and live-process runtime attestation.

The secure source SHA-256 is
`e67058df1347545c39b6c2d27c53e758d68f4ff338e42225acff51836aa91a9a`.
The runtime seal contains 343 regular files, has raw SHA-256
`58c49334a92b2d9e8b02b672d9e67da01b772a626b0ec29e88347d1cc7f92a0c`,
and canonical report SHA-256
`2bf686b6d0ba496383c2dd2afe51e4bffb34ee7b9e18f3bd49e7f68f4bb5184c`.
A mid-run `ps`/`lsof` attestation verified all eight exact no-resume commands
and hashed 64 loaded non-system images per worker: 62 matched the prelaunch
seal and the two additional Python framework images were explicitly pinned.
Its canonical SHA-256 is
`f54adf0e02a0ebb306107c63a974a813815d29485d6507cf33bd2d785a6bc707`.

Run the secure tests from `math/uc/`:

```sh
../.venv/bin/python -I -B verification/test_independent_arithmetic_replay_secure.py
../.venv/bin/python -I -B verification/test_independent_report_acceptance.py
```

Run one slice from trace byte zero:

```sh
../.venv/bin/python -I -B verification/independent_arithmetic_replay_secure.py \
  --runtime-seal verification/results/runtime-environment-seal-secure.json \
  --lock verification/campaign-lock.json \
  --campaign campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8 \
  --slice 0 \
  --output verification/results/independent-arithmetic/secure-full/slice-0.json \
  --checkpoint verification/results/independent-arithmetic/secure-full/checkpoint-0.json
```

**MACHINE-VERIFIED — independently reproduced subject to the trust boundary
below.**  All eight secure workers started at trace byte zero, used no resume
or event limit, reached authenticated EOF with an empty DFS stack, reproduced
every locked terminal tally, revalidated their 343-file runtime seal, and
exited `0`.  The aggregate proves:

- processed events: `488465854`;
- splits: `244232923`;
- terminal leaves: `244232931`;
- residuals: `0`;
- every opcode 0--8 occurs in the accepted reports; opcode 9 occurs zero times.

| Slice | Processed | Elapsed seconds | Canonical report SHA-256 |
|---:|---:|---:|---|
| 0 | 21,135,611 | 6,272.135 | `4d9e46f66124bfb580db36bdc4d391c0a7e82b4a5e11db311bc838dbc0b0b4d6` |
| 1 | 20,827,065 | 6,309.982 | `23b38b2a581bf460cee02e0e2a17ee3c978434e700fd3047e1d42febc6928a97` |
| 2 | 23,486,265 | 6,708.251 | `c2b4385a4480d4f8ed76276921448e3faab94d00942ce2f0c604e987f25709a2` |
| 3 | 34,186,855 | 9,677.585 | `87bccbaf357352e841ae8d2cbabe95d6a5ec116bf142b84e97fd023ba5aa6001` |
| 4 | 71,405,987 | 21,249.065 | `9c49621cdb522de19821bab4dfc0aef9f97e3baa4289845220529ef5be8a559c` |
| 5 | 117,709,435 | 38,905.739 | `97c7f81459d70efc5ec0d11eab3a2c9f86881c923f5b1649894f58e05191b708` |
| 6 | 121,456,775 | 35,566.039 | `09e5bedcf8b23a4eb88c01a2410737c0c32233fafce59c3a62ae94bfe30bd294` |
| 7 | 78,257,861 | 17,728.386 | `7f3a117298c88d477533b0c458b13b183891ad01bf7d01110a5d49a5bc81b91c` |
| **Sum** | **488,465,854** | **142,417.181** | |

The immutable launch record is
[`launch-manifest.json`](verification/results/independent-arithmetic/secure-full/launch-manifest.json).
The live eight-process/image attestation has raw SHA-256
`19c8cb83abd9fba16c2bcae63302e02bcd4602db842a8b3675085902801e126c`
and canonical SHA-256
`f54adf0e02a0ebb306107c63a974a813815d29485d6507cf33bd2d785a6bc707`.
Complete checkpoint/exit logs are in
[`secure-full/logs/`](verification/results/independent-arithmetic/secure-full/logs/);
their `sha256sum -c` manifest has raw SHA-256
`94d9fd6c029e57b87c4425200d8974ef0b962ee0956bbf79ad491fcc820ea880`.

The externally trusted report lock has raw SHA-256
`af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730`
and canonical SHA-256
`7fe63ab9d3c420fd60b648e40aa342bf201b53a03c3b1b46f9f89067534a5c84`.
It fixes all eight raw reports, the launch manifest, and the live runtime
attestation.  The accepted
[`composite.json`](verification/results/independent-arithmetic/secure-full/composite.json)
has raw SHA-256
`92edd4cdf607f65d7a7479bd41b025385b54847a80047d716d8b97f61e5c557a`
and canonical SHA-256
`4acbd3b935bda7e51ed387e42e0598debf22f4a85b7a24ad976c3eaca4f23243`.

To repeat the acceptance step after independently fixing the report-lock raw
digest:

```sh
../.venv/bin/python -I -B verification/aggregate_independent_reports.py \
  --reports verification/results/independent-arithmetic/secure-full \
  --lock verification/campaign-lock.json \
  --report-lock verification/results/independent-arithmetic/secure-full/report-lock.json \
  --report-lock-sha256 af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730 \
  --output /path/to/new-composite.json
```

The accepted trust boundary is explicit: Arb/python-flint mathematical
correctness; custody of the out-of-band report-lock digest; the authenticated
outer trace partition; and a nonmalicious OS, kernel, interpreter, `ps`,
`lsof`, file system, and SHA-256 implementation.  A hostile same-inode
rewrite race or malicious system stack is not excluded.  Thus this closes the
previous same-Python common-mode arithmetic gap, but it is not hardware-backed
execution attestation and does not machine-check the human analytic chain.

## Runtime provenance, not a promised runtime

**COMPUTATIONAL-EVIDENCE:** the committed result records contain producer `elapsed_ms` fields. Their sum is 468,785.706 seconds (`130 h 13 m 5.706 s`) across slices; this is the sum of worker-reported generation durations, not a campaign wall-clock measurement and not a replay estimate.

**COMPUTATIONAL-EVIDENCE:** archived `collect_output.txt` has raw SHA-256 `c521b67208d5169098e0f091a8cf56ba49f6062dcc44ef239a8c705be62512c0`. It reports these prior replay durations. The file is used only as the provenance of the estimate in this table and is not accepted by either new verifier.

| Slice | Producer result `elapsed_ms` | Archived replay seconds |
|---:|---:|---:|
| 0 | 10,365,738 | 2,436.6 |
| 1 | 12,197,641 | 2,442.1 |
| 2 | 17,214,494 | 2,721.5 |
| 3 | 35,005,187 | 4,326.3 |
| 4 | 87,220,753 | 9,746.8 |
| 5 | 138,200,788 | 25,329.8 |
| 6 | 128,286,977 | 21,622.3 |
| 7 | 40,294,128 | 9,115.2 |
| **Sum** | **468,785,706** | **77,740.6** |

**COMPUTATIONAL-EVIDENCE:** the archived replay sum is `77,740.6` seconds (`21 h 35 m 40.6 s`) if the slice durations are added sequentially. It is a historical same-campaign observation, not a guarantee for another host, interpreter build, thermal state, or resource limit. No runtime in this document is extrapolated from an invented throughput.

**MACHINE-VERIFIED fresh-run timing:** the hardened staged slice arithmetic times sum to `58,266.973938498995` seconds; the slowest slice was `16,296.533956208004` seconds. These are observed per-process timings, not guarantees.

## Machine-readable reports and hash interpretation

**MACHINE-VERIFIED:** `uc/verification/report.schema.json` constrains report shape and PASS invariants. Schema validation alone is never acceptance: an executable consumer must also recompute the canonical report hash and compare external pins, per-slice stdout/tallies, input maps, and slice uniqueness.

**MACHINE-VERIFIED:** report hashes have three distinct meanings:

- `inputs.*.*.sha256` is the raw byte hash of an input, recorded before and after.
- `command.stdout_sha256` and `command.stderr_sha256` hash the exact subprocess byte streams.
- `report_sha256` hashes canonical ASCII JSON with sorted keys and compact separators after omitting only `report_sha256`. It intentionally differs from the raw hash of the pretty-printed report file and changes when timestamps, elapsed time, paths, or environment evidence change.

**MACHINE-VERIFIED fresh-run record:** `events.jsonl` is operational logging and `slice-N/latest.json` records the most recent fresh attempt. Its unkeyed hash proves self-consistency, not authenticity; an independent recipient must run the replay rather than trust a copied report.

## Success and failure interpretation

- **MACHINE-VERIFIED:** a fresh slice arithmetic replay succeeds only with `claim_status: "MACHINE-VERIFIED"`, `outcome: "PASS"`, no errors, child termination `exit`, integer exit code `0`, empty stderr, the exact expected stdout hash/tallies, matching environment captures, and identical source/stage before/after hashes.
- **OPEN:** `claim_status: "OPEN"`, `outcome: "PREFLIGHT_PASS"` proves no interval discharge and cannot be counted as a replayed slice.
- **FAILED:** `--resume` is deliberately rejected because a repository-contained unkeyed report cannot authenticate execution.
- **FAILED:** `FAIL`, `FAIL_CLOSED`, timeout, interruption, nonzero exit, stderr, unknown/missing output, missing input, symlink, inventory difference, hash difference, environment drift, malformed JSON, unknown opcode, invalid topology, nonempty DFS, residual, or tally difference rejects the slice.
- **MACHINE-VERIFIED composite:** all eight unique fresh reports pass; the aggregator accepted them in `clean-room-arithmetic/composite.json` with report SHA-256 `57ca3f52c054bd9bccfefc349f440674c74fe7ef671d06c103ef3f52dc33d53d`.

## Cambie bridge strictness

**MACHINE-VERIFIED:** run the check independently of certificate replay:

```bash
"$PY" -I -B uc/verification/cambie_bridge_strictness.py \
  --output "$OUT/cambie-bridge-strictness.json"
```

**MACHINE-VERIFIED:** at 256 Arb bits, with exact rationals

```text
t     = 955165028125263 / 2500000000000000
alpha = 356069 / 10000000
```

and binary entropy `h`, the standalone script certifies

```text
f(t) = (1-alpha) h(2t-t^2) + alpha - h(t)

0.001292783842665984988348621200827023977631257987517743548011480508341589953
  < f(t) <
0.001292783842665984988348621200827023977631257987517743548011480508341589955.
```

**MACHINE-VERIFIED:** the observed Arb ball was

```text
[0.00129278384266598498834862120082702397763125798751774354801148050834158995427309831484301405 +/- 2.31e-76]
```

The hardened repository-local rerun has canonical object SHA-256 `68f3157157d30136f5bbe69332252594a9979f642ed9038cee52fbe8b825b67b`.

**HUMAN-AUDITED:** the bridge use of this positive endpoint also requires the derivative argument `f'(u) < 0` on `[psi,t]` and the elementary `u <= psi` case. The script records those as human-audited dependencies and does not pretend that evaluating one point proves monotonicity.

**HUMAN-AUDITED:** the sequential entropy-to-union-closed proof is in
`uc/AUDIT.md` and `uc/paper/main.tex`; Cambie's Question 2 / Section 4 is its
source. The standalone Arb computation checks only the strict endpoint.

## Exhaustive entropy-bridge control

**MACHINE-VERIFIED (exhaustive finite):** the sequential coupling behind
`cor:uc` is rebuilt from the corollary statement alone and exercised on every
nonempty family of subsets of `[n]` for `n <= 4`:

```bash
"$PY" -I -B uc/verification/entropy_bridge_exhaustive.py \
  --max-coordinates 4 \
  --sample 5:5000 --sample 6:2000 --sample 7:500 --sample 8:150 \
  --sample-closed 5:3000 --sample-closed 6:1500 \
  --sample-closed 7:600 --sample-closed 8:200 \
  --output "$OUT/entropy-bridge-exhaustive.json"
```

Single core, 110.0 s observed. The verifier imports no campaign snapshot, no
certificate module, and not `uc/bridge_uc.py`. Exhaustive coverage:

| n | families | coupled prefix states | union-closed families | smallest certified positive iid slack | smallest certified positive coupled slack |
|---:|---:|---:|---:|---:|---:|
| 1 | 3 | 3 | 3 | none (all ties) | none (all ties) |
| 2 | 15 | 43 | 13 | 4.337e-2 | 1.258e-1 |
| 3 | 255 | 2,083 | 121 | 2.878e-3 | 9.910e-3 |
| 4 | 65,535 | 1,629,751 | 4,959 | 1.278e-3 | 5.857e-4 |

Seeded deterministic samples beyond the exhaustive range (seed 20260827):

| n | uniform families | prefix states | random union-closed families | largest closed family | smallest max frequency |
|---:|---:|---:|---:|---:|---:|
| 5 | 5,000 | 391,840 | 3,000 | 32 | 1/2 |
| 6 | 2,000 | 510,708 | 1,500 | 58 | 1/2 |
| 7 | 500 | 391,018 | 600 | 103 | 4/7 |
| 8 | 150 | 388,571 | 200 | 146 | 20/31 |

Exact rational arithmetic checks the four coupling masses, both Bernoulli
marginals, `P(A_i or C_i | prefixes) = s*`, prefix-by-prefix uniformity of `A`
and `C`, `law(P_i) = law(R_i)` with mean equal to the frequency of element
`i`, and the chain rule `H(A) = sum_i E h(P_i)`. 256-bit Arb certifies both
data-processing inequalities; every enclosure is either provably nonnegative
or a structural tie inside `+/-2^-200`. Every union-closed family reached,
enumerated or sampled, has an element of frequency at least `t_cert`, the
extreme case being exactly `1/2`.

Report: [`results/entropy-bridge-exhaustive.json`](verification/results/entropy-bridge-exhaustive.json),
canonical object SHA-256
`6681f8faf13d9344f8e7bf3a5d7785ed289f1e82c5e9badd1e7967c03e3ebfef`.
That digest identifies the stored artifact, which carries `finished_utc` and
`elapsed_seconds`; a rerun reproduces every count and verdict but not the
digest. The coverage counts and the seeded sample identity are the
reproducible quantities.

**MACHINE-VERIFIED (fail-closed):**
`uc/verification/test_entropy_bridge_exhaustive.py` passes 4/4 and shows the
control has teeth: dropping the `max` branch of `s*`, replacing the prefix
conditional probability by a constant, or inflating the union bit each makes
it reject. Replacing the `1/2` clip by another admissible cap does not, since
that coupling is still valid; the identification with Cambie's `s*` is the
separate three-case proof in `uc/bridge_uc.py`.

**OPEN:** this is a finite control. The universal statement is the induction
in `uc/paper/main.tex`, and the control does not touch Theorem `thm:main`, the
certificate, or the strictness estimate at `t_cert`.
