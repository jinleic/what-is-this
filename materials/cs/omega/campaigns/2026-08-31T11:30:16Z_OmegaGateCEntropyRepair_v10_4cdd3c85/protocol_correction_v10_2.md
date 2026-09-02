# v10.1 audit round-2 fix (2026-08-31, static only)

Status: ZERO real point compute. Only py_compile + prepare-only smoke +
mock-self-test.

## Frozen artifact provenance corrections (Main round-2)

1. `HERE` now resolves to `dirname(abspath(__file__))`, so the v10
   artifact loads THIS campaign's checksummed dependencies (parent v9
   campaign is no longer referenced).
2. The three expected dependency SHA-256 digests are hard-coded and
   verified against the actual bytes BEFORE any import; mismatch
   hard-aborts.  Removed:
   - the tautological `expected = _sha256(path)` check,
   - the needless live-source equality dependency.
3. Guards metadata made exact:
   - `branch_equation` now records the FULL mass-aware generalized dual
     `U_m(y) = m*logsumexp(B^T y) - y^T b - m*log(m)` with KKT candidate
     `q(y) = m*softmax(B^T y)` (replaced the obsolete simplex
     `logsumexp(B^T y) - y^T b` record).
   - `v10_feature_matrix_fix` now records the exact compact-lambda +
     transpose construction (Part uses `m[:, lo:hi+1].T` per dim with
     lo/hi from `lam_low`/`lam_high`; GlobalStage uses full `m.T`;
     omitted Part columns asserted structurally zero; rows-per-region ==
     lam width; columns == distribution length), replacing the plain
     `vstack(m.T)` description.
4. `__pycache__/` preserved (not deleted per directive) and explicitly
   classified as static py_compile output, EXCLUDED from the checksum
   ledger (marked unchecksummed/non-evidence in manifest.json).

## Artifact inspection (performed on the frozen artifact bytes)

- HERE = dirname(abspath(__file__)) confirmed in asrun_v10.
- EXPECTED_DEPENDENCY_SHA256 dict matches Main's three digests verbatim
  (67ca5f..., 31657c..., cc423a...) and is consumed before loading.
- prepare-only smoke (on the frozen artifact) confirmed the new
  `branch_equation` string, the new `v10_feature_matrix_fix` string,
  hard-coded dependency digests, structural guard 63 Parts + glob with
  B (7,4)-(27,45), and no live-source reference.
- Mock self-test still PASSES (no real point evaluation).
