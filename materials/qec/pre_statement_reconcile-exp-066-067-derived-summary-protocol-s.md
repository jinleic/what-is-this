# FREEZE — EXP-066/067 derived-summary and package dependency reconciliation

Cycle: `20260909T103627Z_e1b43a`. Actor: GPT-6 Astra, directly in the coordinator. Independent reviewer: read-only `institute-opus-reviewer` (Claude Opus 5); no self-review.

## Funded question

Can the EXP-066/067 derived summary protocol, shard and global-screen bindings, and arXiv/package/checkpoint descriptors be reconciled against certified EXP-055 run `20260909T095538Z_4809c2f8_816c74c3bcbd`, with all scientific evidence and scope unchanged, before any package or repository-consistency publication?

The source run's RESULTS, AUDIT, publication-verification.json and downstream-publication-holds.json are the prior evidence authority. The funded objective is the brief's complete next_action, not new Liang or residual screening. The initial discrepancy is already established by the source run; do not repeat its scientific searches or its witness-wide reconciliation campaign.

## Falsifiable acceptance criteria

1. **C1 — Certified inputs remain exact.** Verify all 25 current EXP-055 artifacts against the source run's installed hashes. Preserve its 24 shards, n<=234 aggregate, exact-reference certificates, preimage archives and frozen run bytes. Preserve all n=270 scientific records and threshold 20 / row-a attribution; do not execute the 21 residual k=8 searches or any Liang search.
2. **C2 — Complete summary dependency binding.** Regenerate only the derived EXP-066 summary via EXP-067's existing `_refresh_exp066_summary` path using the current certified inputs. Run the existing strict EXP-066 summary validator in a fresh process. Its current protocol, both shard hashes, record identities and global-screen hash must agree with the actual bytes. All scientific summary fields, scope, automorphism partition, counts and verdicts must equal the predecessor; only generation/provenance bindings may change. Retain its exact preimage in this new run.
3. **C3 — Current versus historical attestation.** Keep the old empty-survivor certificate byte-unchanged and explicitly historical. If a current attestation is produced, invoke the existing fail-closed screen-certification producer with its output redirected into this new run; assert the already-closed n<=234 scope and empty survivor list before invocation, prohibit solver calls, and validate the resulting exact aggregate hash. Do not imply n=270 closure.
4. **C4 — Package/checkpoint dependency closure.** Reconcile `reports/arxiv_metadata.json`, `reports/arxiv_package.md` and affected `checkpoints/verified_results.json` evidence descriptors to the current summary, aggregate, current attestation and source run. Check actual file sizes/hashes and local links for the affected package inventory; retain original descriptor bytes in the new run. Historical scientific statements remain dated historical statements, not silently promoted to current claims. Preserve paper source/build bytes, prior theorem claims, authorship placeholders, `do_not_post` and no-external-publication status. Distinguish this dependency-consistency gate from full manuscript validation or a repository-wide green test suite.
5. **C5 — Fresh audit and independent review.** Run the exact verification command below afresh after candidate construction, exercise focused existing EXP-066/067 consumer checks, and verify adversarial stale protocol/shard/global bindings and wrong aggregate scope are refused. Verify installed bytes equal the reviewed candidates, unchanged inputs remain pinned, and the required independent Opus review is complete before choosing the terminal verdict. Record limitations, command receipts and measured wall usage. No publication or external communication is authorized.

## Scope and execution bounds

Allowed mutations: this preregistration before init; the minted run's evidence and closing artifacts before freeze; derived `results/processed/exp066_n234_frontier.json`; the three descriptor files in C4; narrowly necessary EXP-066/067 producer/validator fixes only if a demonstrated defect prevents C2 (with fail-before/pass-after evidence and all callers migrated); the brief's owner state fields and math ledgers after close; the cycle usage record. No existing certificate, screen shard, aggregate, preimage archive, frozen run, source preregistration, institute budget/lifecycle/priority block or monitor file is hand-edited. No service install/start, no new cycle, no science-worker delegation.

Heavy deterministic verification runs on mini-pro via local `scripts/institute_monitor.py run`, in the existing path-identical QEC tree, one lane, `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`, nice 15. Stable scientific objective id: `exp066-067-dependency-reconciliation`, unchanged across attempts/hosts/runs. Fixed UTC deadline: **2026-09-09T12:40:00Z**, including replay and at least 60 seconds termination margin. Stall bound: 1800 seconds. Projection: at most 30 minutes of deterministic validation plus descriptor work/review within a two-hour campaign window; no solver search, full project suite, paper rebuild or multi-day work. If that cannot fit, stop under the original objective; do not relabel retries to extend the window. Existing compute, if any, must be identified and adopted, never relaunched.

## Exact verification command

From `math/qec`, with `<RUN>` substituted only by `campaign.py init`'s minted id:

```text
.venv/bin/python -B campaigns/<RUN>/reconcile_dependencies.py audit
```

Run that command on mini-pro through the prescribed monitor, with original deadline and run binding. The run-local verifier owns C1-C4 checks and negative controls and emits `audit.json`, `audit.log` and `progress.jsonl`. Focused existing consumer checks: `.venv/bin/python -B -m pytest tests/test_exp066_n234_frontier.py tests/test_exp067_n234_connected_cluster.py -q -p no:cacheprovider`, also through the same monitor objective. It is not a project-wide suite. Candidate construction uses the same verifier's `produce` command; no rerun of `exp066 run --force`, `promote-reference`, or native Liang solver is allowed.

## Verdict mapping

- `FROZEN-CERTIFIED`: all C1-C5 pass, reconciliation installed exactly as independently reviewed; certification is dependency reconciliation only, with unchanged scientific claims.
- `FROZEN-NEGATIVE`: reproducible scientific evidence contradicts the funded consistency proposition (for example a genuinely incompatible scientific summary under the unchanged certified inputs), and the independent review confirms the obstruction. Ordinary stale metadata repaired here is not a negative scientific result.
- `FROZEN-INCONCLUSIVE`: valid bounded verification yields insufficient scientific evidence to settle the original question; specify the unclosed obligation and do not claim consistency.
- `CRASHED`: operational/tool/execution failure prevents completion, preserving exact run/evidence for recovery. `REJECTED`: administrative invalidity before meaningful science. Neither maps to a scientific negative.

Freeze closing evidence only after the audit and independent review. Close with the honest mapped verdict; then write measured usage outside the sealed run, append the math/PROGRESS.md and target README entries with exact-text deduplication, update owner semantic fields only, and run campaign.py state refresh locally.
