# How to verify this bundle

Two modes. Neither ever skips a check silently.

## A. Standalone, from the extracted archive (seconds)

The archive carries sources, artifacts and the paper, but **not** the
8,500,211-graph `.g6` catalogs (≈200 MB). Provenance re-hashing therefore
cannot run here, and the verifier says so out loud:

```bash
python r55/src/verify_mixed_independent.py --no-input-hashes \
    r55/data/engstrom_identity.json
```

Prints `inputs NOT re-hashed (--no-input-hashes): provenance of 64 recorded
files is UNVERIFIED in this run`, then still proves, independently of every
producer module:

- schema, campaign id, disposition, key sets, trust-root sources;
- the 656-graph and n=49 replay records;
- all 53 catalog classes: histograms sum to their graph counts, the recorded
  catalog windows equal the histogram extremes, the outer windows contain them,
  and the total is exactly 8,500,211;
- all **3,215** F intervals, recomputed from the artifact's own motif windows
  with independently written affine-interval arithmetic and the verbatim n=45
  coefficients;
- the (d,a,b) → (deficiency, balance, g, h) projection against the frozen
  `higher_identity_m4.json`;
- every certificate: sign gates, per-state affine inequality, least-α
  canonicality, exact 45α;
- every primal witness: canonical nonneg weights, Σw = 45, Σw·balance = 0, and
  Σw·g_lo ≤ 0 ≤ Σw·g_hi, Σw·h_lo ≤ 0 ≤ Σw·h_hi, Σw·f_lo ≤ 0 ≤ Σw·f_hi;
- the derived route statuses and the terminal disposition.

Expected last line: `MIXED EVIDENCE VERIFIED: MIXED_NO_CUT_IN_FROZEN_BASIS`.

## B. Full, inside the repository working tree (seconds + optional 75 min)

```bash
cd <repo>/math
./.venv/bin/python r55/src/verify_mixed_independent.py \
    r55/data/engstrom_identity.json          # adds: re-hash all 64 inputs
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*mixed*.py' -v
./.venv/bin/python r55/src/mixed_deficiency_cone.py \
    --output r55/data/engstrom_identity.json # optional full rebuild, ~75 min
```

## Manifest

Every archive entry is listed in `MANIFEST.md` with its SHA-256 and byte size.
Check one with `shasum -a 256 <path>`.

## What is NOT verified anywhere yet

The plan's Task 4 disjoint checker is incomplete. Specifically: no second,
independently written motif kernel has re-swept the 8.5 M graphs; the g and h
endpoints are anchored by exact equality against the frozen m=4 artifact rather
than recomputed from the R1–R5 relation system; and the envelope-LP windows for
the 77 catalog-missing stratum classes are not independently re-derived (the
m=4 campaign scoped this identically). Treat the disposition as
producer-verified, independently re-verified at the search and F-interval
layers, and cross-anchored — not as a fully disjoint replication.
