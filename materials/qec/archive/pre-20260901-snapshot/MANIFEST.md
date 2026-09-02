# Historical snapshot: divergent members of `math/qec.zip`

- Created: 2026-09-01
- Source container: `2026/2026-09-01/superseded-bundles/math/qec.zip` (SHA-256 `4e54a41d8ccdc603ad747f660315cc2e4dfcb6856a2c81c03b8e4d1337b57e08`), 281,927,858 bytes.
- Reason: `math/qec/.git` contains ZERO commits, so git conserves no history for this tree.
  The zip was the only prior snapshot. Rather than retain 281 MB of overwhelmingly
  redundant bytes, the divergent members are conserved here and the container is retired.

## Proof the live tree loses nothing

Compared all 5,119 real members (excluding `__MACOSX/` resource forks, `.git/`, `.venv/`,
`__pycache__/`, `.pytest_cache/`) against the live `math/qec/` tree by SHA-256:

| outcome | count |
|---|---|
| identical to live | 5,066 |
| path missing from live | **0** |
| differing content | 53 |

The 53 differing members were each shown to be an OLDER state of a live file:

- `experiments/exp023_light_gensets.py`, `exp027_delta_audit.py`, `exp030_surface_baseline.py`
  differ by exactly -4 bytes each: the 2026-09-01 alias migration rewrote the string
  `qec-codesign` (12 chars) to `math/qec` (8 chars). Live is correct.
- `checkpoints/verified_results.json` is 72,542 bytes SMALLER live, which looked like loss.
  It is not: live holds 57 result records vs the zip's 37, with **0** record ids lost,
  **0** dropped fields on surviving records, and **0** bytes of field-value shrinkage.
  The size drop is JSON compaction; live is a strict superset.
- The remaining 49 members are all LARGER live (accumulated notes, checkpoints, reports,
  partial-run logs).

## Contents

| path | bytes | sha256 |
|---|---|---|
| `README.md` | 18,190 | `5d3b6dea3df83518bd0a23a0acb281e5f062229a24d1016a2c3b9386feeae402` |
| `checkpoints/current_state.md` | 23,010 | `6e8f4c3eef4fde57717707daa17278480ac5988614ae57047dfbacc3044a58bc` |
| `checkpoints/next_actions.md` | 9,803 | `dc480106bcd87f899b964f93f5652c239756572bf53b0e7a220fb53d47392be9` |
| `checkpoints/verified_results.json` | 379,783 | `d3aa4f9969ac369b7a832439616f217831d9e092a3d354120541f1c4e7d5bb82` |
| `experiments/exp023_light_gensets.py` | 84,841 | `39d3f757f9857831481a279dbae8ef1de70281073432c52cf104207617a9552f` |
| `experiments/exp026_decoder_codesign.py` | 53,300 | `e066e0998840f3d4f9b08f80c770e9fd2b28b147b92954f34624c3c7bc9642ee` |
| `experiments/exp027_delta_audit.py` | 57,273 | `a60b11ed00528195af423ac40d039ee15eb88cc23abac0bd5868f9cc2ae03ac5` |
| `experiments/exp030_surface_baseline.py` | 40,258 | `b271658058f6d3f9f3eae318314c252841313a69ae5dca9b52b4e82757bcc489` |
| `notes/failed_routes.md` | 43,653 | `d37f8b70a69300a49c6c6a1e10d0b0f20f09ebdf7f38a1cd6f9854db61d3eff8` |
| `notes/novelty_matrix.md` | 9,288 | `32eeb5c56de43e8faaae5dcf0fcb5af9df28362856d1a92908147112eac8dac8` |
| `notes/open_status.md` | 9,576 | `cdda08b511e03968518d212904d6815da33b22513b77db6e7a19d89c769d137f` |
| `pyproject.toml` | 823 | `2761202893e4430e2e0bc140b5a08872fc8d8f8550661cab8180c2491a388584` |
| `reports/arxiv_metadata.json` | 11,366 | `7957e925d82f8a72840b6bb2ee2b88492eb2055836ac8079ad9617ccca83f9ce` |
| `reports/arxiv_package.md` | 16,444 | `88bccc7244564b4a95ffb508f29913b9b9c34ca6496992715137071e3fe28ce7` |
| `reports/claims_matrix.md` | 35,288 | `f9649c02384c57d7243952bc0d73de70e412327a28251f6f61944ede9e17d1fa` |
| `reports/paper_pbb_nogo.aux` | 3,981 | `57eb6e7189b7b799eb9e46797296e9e510179a07edaab5990349a2fc50b8b37e` |
| `reports/paper_pbb_nogo.bbl` | 1,978 | `202ca3aa9b2fca1ac480f95b6aedaea490bdf0ea2a6ea5e43f4cc7b13a292a93` |
| `reports/paper_pbb_nogo.bib` | 2,772 | `9dd260e89646141bb9183a8319a252392e888d870b8d3a24ab8c080de20e34e1` |
| `reports/paper_pbb_nogo.blg` | 1,213 | `cd37b93a21bbc8bba61209c6f60a6bb30d9297769cdd0ee206839fca00c705e6` |
| `reports/paper_pbb_nogo.log` | 21,186 | `3ca2d73799fc837803c3ff01463fd6f5e0863f5789f41ddf820d7ba683b5082c` |
| `reports/paper_pbb_nogo.md` | 35,628 | `cdbdc86b8a8ab7c7c01e89da2d7304512dbe6a778047bc9e26c402eee4398669` |
| `reports/paper_pbb_nogo.out` | 3,240 | `24b738d5afe41547581d896e330557e6fc80684713ccc5b67fd2ba7919a9c67a` |
| `reports/paper_pbb_nogo.pdf` | 626,136 | `2fa06a95f1f2c517d44fe482dcdc4a3c74d9c0d3412f7103cede1ea1b16e46e2` |
| `reports/paper_pbb_nogo.tex` | 39,866 | `c4995c815a2afcbf7cf98e8984ee9156d5ea8931ddb2b34bc92abfed84a3e15b` |
| `results/partial_runs/exp037/row_0051.json` | 3,293 | `c0f358d3d61072a68b55b61fbb33c6006ebb4c2e3888c2b1a070a01ea42a0bc5` |
| `results/partial_runs/exp037/row_0052.json` | 3,292 | `14594d1e4cfb0e703b6b4dfcfb42a99af7e49fc068967fc7751908f6849c153f` |
| `results/partial_runs/exp037_envelope_classification.json` | 125,012 | `8895dc55caa95937423037b98b343db92d3d157f2c580779fabc76dfceaecf4f` |
| `results/partial_runs/exp039/shard_0.log` | 342 | `9c5a897f12dce1ffa1e59dd518d17c1b58e45422e2116e5d61be619e937918d4` |
| `results/partial_runs/exp039/shard_1.log` | 240 | `8c209265adc431fcdb0a6b786fc170bfea465fd776da3f0a4ca553137335feb8` |
| `results/partial_runs/exp039/shard_2.log` | 141 | `3e108ced2934da0463ca6930cbbf7d64c178aec51163fb4eed746e5fe41cfa53` |
| `results/partial_runs/exp039/shard_3.log` | 444 | `ad09b81134a985c8455f38200cdffb1ef329ac2dce8a4907a37e1318be5779ad` |
| `results/partial_runs/exp039_nogo_module.json` | 66,430 | `e8ff96f5403fb946edca83599bf2d1027376c16ee2e2cfc024bf31ed10326842` |
| `results/partial_runs/exp040/mc.lock` | 5 | `d74acceddbb2edbf0f315323cdc61b44616128742edb2927bf17bbf9a923f96e` |
| `results/partial_runs/exp040/mc_driver.log` | 992 | `eb0d578628b06702a224d7e05aeda3651617371423669112212754198f2b8934` |
| `results/partial_runs/exp040/mc_partial.jsonl` | 23,084 | `d62d40a9b0bbc54acd4fd1167728e6f507272a600d4ff3c6f2d62f2d2a3c5d13` |
| `sources/bibliography.bib` | 5,230 | `1ebe0032c4352a1d1edfa6f48f63746f649f49e3764d8d7a619a4dbcebc47623` |
| `sources/manifest.yaml` | 6,917 | `f62b5017869894a89a8da51f2d7418fddf45a0f118f372f906e716ba2d18d1a2` |
| `src/qec_research/distance/sat_decide.py` | 20,934 | `eec17e07f3c45972cf0af1cfc33f1e29570d7e1ffb32c0937a38f9f4cabfebcc` |
| `tests/test_exp036_delta_closure.py` | 10,884 | `571d29d4979d3bef98ba4c45a25d9240247153d562e68a41aef36be3501831d1` |
| `tests/test_paper_claims.py` | 12,099 | `1118f5aceb18d17dd5b8f158584995be526555a6a67cdfd52450a7225b6239e1` |
| `third_party/qLDPC/docs/source/examples/basics.ipynb` | 30 | `c17688d8b7a1655e8d8a23daeed7dded7b3acd1a8764de4d5e16d2e9b3702181` |
| `third_party/qLDPC/docs/source/examples/bivariate_bicycle_codes.ipynb` | 47 | `c3e0ef4b4be4b2cfe1c952797052c3bf13f2eb03479f26b1685d095ad3442187` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/1_code_capacity.ipynb` | 62 | `d53f9dfe05116398c7a6dbd28abb71b316440e56acc24b18fbebb4b5e3963c1f` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/2_quantum_memory_x_or_z.ipynb` | 70 | `9d43b4733077c87bb8651944aefb1985c31eb7f4115e9a95a8825aa04d0196a8` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/3_quantum_memory_combined.ipynb` | 72 | `d887882d427b3975b932c41cc3a876c996b32752da204554d8e78de08b791ca7` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/4_sliding_window_decoding.ipynb` | 72 | `c754467446cd0c146de8d7535e5b7c8b727d586298c5e91067e3ea6c9017e766` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/5_state_preparation.ipynb` | 66 | `175c9f81b4cca1bb1471afa99156351cc6948fb8adfed3f4ee7fac22f31fd4a0` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/6_knill_qec.ipynb` | 58 | `983aae6931d66e0234c48435fdf1a2a3585d53bb69aac2e8a79811ee3bff1d8b` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/misc/alpha_syndrome.ipynb` | 69 | `2dd95b4aab78e20460e589910391a07c6fbec7217900eeac1b664bd61025d4e5` |
| `third_party/qLDPC/docs/source/examples/logical_error_rates/misc/decoding_tqec_circuits.ipynb` | 77 | `479716c3caffd52a3bbe293c5f1243c9a59507a0ca67db3fed108e7feb841a83` |
| `third_party/qLDPC/docs/source/examples/noise_models.ipynb` | 36 | `73c1042cd87c62f4604e4421f6676a57e895b4a519631d80721cbbbc2b744379` |
| `third_party/qLDPC/docs/source/examples/transversal_gates.ipynb` | 41 | `6c32f29d1a00f95629897eb7f8cd8db5c76cf0f3e9ab1a8da91c0618cd005f77` |
| `uv.lock` | 206,223 | `0a9ad8571fc7e8aa291965225215c05a6ba85a705993559685a0690b411a877d` |

## Known broken links in the archived README (intentional)

The stored `README.md` (hash above, never edited) references two sibling
targets that the conservation policy deliberately did NOT copy into this
snapshot, because the live tree already holds them unchanged:

- `reports/technical_report.md` — live at `math/qec/reports/technical_report.md`
- `proofs/pbb_structure.md` — live at `math/qec/proofs/pbb_structure.md`

The README's links are therefore intentionally dangling inside the snapshot
and are reported by `scripts/knowledge-system-audit.sh` as archival warnings,
not actionable rot. Restore readers should resolve them against the live
`math/qec/` tree per the Restore section above.

## Restore

```bash
# these files are the PRIOR state; the live tree is newer. To inspect one:
diff math/qec/archive/pre-20260901-snapshot/<path> math/qec/<path>
```
