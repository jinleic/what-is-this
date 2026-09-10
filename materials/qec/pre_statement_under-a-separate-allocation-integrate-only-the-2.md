# FREEZE — integrate exactly 21 reviewed n270 k8 witnesses

Cycle: `20260910T004403Z_37f99d`. Actor: direct GPT-6 Astra Main. Required independent reviewer: `institute-opus-reviewer` (Claude Opus 5, read-only), never Astra self-review. One campaign; no new cycle.

## Gate question

Can exactly the 17+4 reviewed weight20 logical witnesses in the closed source campaign `campaigns/20260910T001649Z_faeb16bd_a50d496cd443/search.json` be installed in the canonical EXP-055 15x9/27x5 shards using the existing EXP-069 finalization and EXP-055 strict validation conventions, with correct provenance/archives and no changes to other records, threshold/source, or closed n<=234 aggregate?

Source search SHA-256: `4db51fd1451372d60ef4e1e545a50d3d8205fb0625e2d086677f329e6ef4797d`. Exact source inputs are the full `input.record` objects, not run-id or record-index matches alone. Selection: 15x9 indices 1891,1909,1923,2281,2313,2316,2430,2592,2849,3149,3297,3477,3607,4261,4704,4969,4995; 27x5 indices 282,471,532,577. Current canonical records equal all 21 full source inputs at admission. No integrity_hold is present in the read target state.

## Falsifiable acceptance criteria

1. Pin source-run frozen hashes, search bytes, full preimage shards, science modules and required certificate/dependency bytes before execution. Verify source status FROZEN-CERTIFIED and immutable inputs. Exact counts: 5024+688 records; exactly the selected 21 currently undecided k8 records, each threshold20 / `EXP-070 [[270,8,20]] Liang row a`. Admission mismatch stops installation.
2. Reverify all21 imported supports from matrices in a fresh monitored process: integer unique in-range coordinates, physical weight20, kernel membership, non-stabilizer membership, X/Z transport and both existing NumPy/bitset paths. Compare the freshly rebuilt verification with the exact source verification. These are existing independent algorithm paths, not an independent external software implementation.
3. Apply `exp069_n270_frontier.finalize_deep_record` only to the selected original record and its exact source `deep` result. No new deepening/search. Recompute shard derived aggregate fields through the existing `_screen_lattice_payload` conventions, retaining other preexisting metadata. New provenance must distinguish historical 2000-try pipeline metadata from the imported 20000-request results and distinguish prior search time from this integration's time; actual prior completed order counts are unknown. Preserve all 5691 nonselected record objects verbatim semantically, survivor/no-reference lists, code identities, orbit multiplicities, thresholds/sources, and solver-call accounting. The only scientific transitions are 21 undecided -> dominated_by_deep_reduction.
4. Preserve every old archive/certificate. Archive both complete immediate shard preimages using the existing validator-version/content-addressed archive layout; bind the integration to those hashes and the source search. Keep the existing reference-rebind binding historical and unchanged; verify both its archive and the immediate integration archive. Strict validator-v11 aggregate/record checks must pass on both final shards. Historical archive checks may use only the existing explicit stale-threshold mode; never weaken current shard validation.
5. Fresh audit command below must pass again in a separate process on the exact candidate bytes. Exercise adversarial controls rejecting wrong selected identity, wrong threshold/source, malformed/corrupted supports, unselected-record mutation, stale protocol/binding, stale aggregate and premature n270 closure. Controls are run-local; no new permanent test suite or source special case.
6. Check actual derived bindings: the current EXP-055 aggregate retains the same explicit 22-lattice n<=234 list, 4862 classes / 150581 pairs, 4658 dominated and 204 no-reference; no survivors/undecided and no n270 inclusion. Its bytes, EXP-066 summary, arXiv metadata/package, checkpoint descriptors, all other shards and paper bytes remain unchanged. Current descriptors do not name the two n270 shards; verify current referenced file hashes/identities and document any historical evidence separately. No publication-readiness inference.
7. Install only the reviewed/fresh-audited candidates (two canonical shards and two additive preimage archives), verify installed hashes equal reviewed output hashes and recheck protected local inputs. Require completed independent Opus review with no unresolved blocking issue before choosing certification. Report final zero k8 undecided, still147 k12 undecided, and210 unpromoted k20/k24 survivors; global n270 closure remains open.

## Scope and resource boundary

No resolve_screen, screen_initial, new search, reduced_witness_bound, Liang rerun, solver invocation, lower-bound proof, k12 processing, survivor promotion, exact-distance claim, global n270 closure, publication or source algorithm modification. Data integration and run-local verification only. Validation of unchanged certificate bytes is allowed; executing their solvers is not. The historical n<=234 aggregate and publication holds must remain unchanged. Canonical installation is withheld until validation and review; full preimages supply exact recovery evidence, with no deletion of prior artifacts.

Stable scientific-objective task id: `n270-k8-residual21-canonical-integration`. One lane on freshly identity-verified `mini-pro.local`, nice15, one BLAS/OpenMP thread. Original task start `2026-09-10T00:49:00Z`; fixed objective deadline `2026-09-10T01:49:00Z`; monitor termination deadline `2026-09-10T01:48:00Z`, leaving60 seconds. Every attempt/replay retains that task/start/deadline. Launch only through local `scripts/institute_monitor.py run` with explicit stall bound1200 seconds, run binding and real stage/progress outputs. No service or daemon installation.

Projection: preparation/deployment300s + produce/strict validation900s + independent fresh audit900s + review/installation/closure900s + termination reserve60s =3060s, below3600s and72h. Inspect actual progress and do not launch another verification that cannot fit the remaining original window. Heavy validation runs on the mini, all campaign/state/ledger commands run locally. Full actor wall usage is measured from cycle dispatch `2026-09-10T00:44:03Z`; tokens null if unavailable. The available target allocation is not enlarged.

## Exact verification command

After campaign init, `<run>` is only the minted id, not a free parameter. The preserved run-local driver takes its pinned `inputs.json`; stage output and progress paths are inside that run on the mini:

```text
/Users/jinleic/jinleic-workspace/math/qec/.venv/bin/python -B /Users/jinleic/jinleic-workspace/math/qec/campaigns/<run>/integrate_witnesses.py audit --output audit-fresh.json
```

Execute via local `python3 scripts/institute_monitor.py run --domain math --target qec --run <run> --host mini-pro.local --task n270-k8-residual21-canonical-integration --task-started-utc 2026-09-10T00:49:00Z --deadline-utc 2026-09-10T01:48:00Z --remote-cwd /Users/jinleic/jinleic-workspace/math/qec/campaigns/<run> --output progress.jsonl --output audit-fresh.json --stall-seconds 1200 -- /usr/bin/nice -n 15 /usr/bin/env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 <exact verification command>`. Fetch outputs and completion receipts and verify hashes before audit/close. Exit0 alone is insufficient.

## Verdict mapping and closing contract

- FROZEN-CERTIFIED: every acceptance criterion passes and independent Opus review has no unresolved blocker; exactly21 imported canonical domination witnesses only.
- FROZEN-NEGATIVE: a reproducible scientific counterexample refutes the frozen witness/identity gate with intact execution and evidence; never an infrastructure or budget failure.
- FROZEN-INCONCLUSIVE: scientifically valid partial evidence within the original gate, but remaining scientific obligations prevent certification. No narrowed substitute gate or partial installation claimed as success.
- CRASHED: runtime/tool/model/compute failure or deadline prevents completing the contract. Administrative REJECTED is allowed for an admission/authorization obstruction. Preserve exact run/evidence for recovery; never start another run.

Write PLAN, command/outcome ITERATE, fresh skeptical AUDIT, RESULTS and provenance before campaign freeze/close. No writes into the run after sealing. Main is the sole domain owner/coordinator: exact-full-entry dedupe into math/PROGRESS.md and target README, never math/RESULTS.md. Update only headline/next_action and current_gate if actually moved; campaign.py owns mechanics, institute owns lifecycle/budget/priority/status. Write measured usage outside the run and refresh target state locally; verify closed status and usage on disk.
