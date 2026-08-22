# Certification closeout checklist — campaign I

Run through this **the moment** `collect_output.txt` in
`uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/`
prints `COMPOSITE CERTIFICATE` — not before.

## 0. Verify the certificate itself (no shortcuts)

```sh
cd /Users/jinleic/jinleic-workspace/math
./.venv/bin/python -B uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/snapshot/cert3_collect.py \
    uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/launch.json
# Must print: ALL MACHINE CHECKS PASS … COMPOSITE CERTIFICATE
```

Record verifier parity (determinism): slice 0 trace SHA-256 must equal
campaign H's slice-0 trace SHA-256 and campaign E/F's (`9eed0cc1382f5485…`).
Add all 8 SHAs to the ledger below.

## 1. Ledger (PROGRESS.md)

New entry `2026-08-1X — CAMPAIGN I CERTIFIED\(t=\psi+10^{-4}\)`:
- 8-slice table (boxes, stack=0, residual=0, budget_time=0 each)
- all 8 trace SHA-256 prefixes; total boxes (sum over slices)
- replay evidence (iwatch markers, 11 tallies, tree invariant); collector
  output block pasted verbatim
- statement of the implied theorem with the Cambie Q2 conditional spelled
  out (t rational-exact proof)
- note that ψ=0.3819660112501051 was the largest previous constant with a
  complete proof, so this certificate beats it (librarian: re-check arXiv
  for any later full-proof result before claiming)

## 2. `uc/campaigns/README.md`

Replace the I row's `pending` with:
`CERTIFIED — collector replayed all 8 traces; COMPOSITE CERTIFICATE printed`.

## 3. `uc/README.md`

Section "The ladder, and what is actually proved": move our result to a new
row above ψ with status **proved (machine-checked + replayed)**.

## 4. `uc/PROOF.md`

Fill `[I-FINAL]` placeholders with I's collector output. State H/I slices-7
byte-identity conclusion (determinism).

## 5. `uc/ANNOUNCEMENT.md`

Replace `[H-FINAL]`/`[I-FINAL]` with the final numbers; timestamp it.

## 6. `uc/paper/main.tex`

Fill `[I-FINAL]`; add the collector output block (verbatim); re-run
`pdflatex` twice; zero errors / zero undefined references required.
Then write `uc/paper/main.md` (Markdown mirror) and archive
`uc/paper/main.pdf` alongside.

## 7. Final readiness

- `hi-precision` block: rerun `uc/...` sanity (already PROVED pieces stand).
- One-sentence English statement of victory in the final report to the user
  (claim shape: "largest explicit constant in the union-closed sets
  conjecture with a complete, independently replayable machine-checked
  proof", with the Cambie-implication caveat).
