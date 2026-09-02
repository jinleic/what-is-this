# Pre-statement — targeted cached refinement of the unresolved direct-D.4 midband

Campaign: `20260831T145732Z_kg_direct_d4_midband_panelcache`  
Status: **COMPUTE PENDING** — written before any import, static probe, envelope
evaluation, cache build, or band run in this campaign.

## Parent evidence and fixed target

The frozen restart result
`20260831T082425Z_kg_direct_d4_restart/logs/band_1p75_3p5_result.json`
(`324a44396ca0ea388e327059028b935dc37e3ca95471a88041a2e30454241347`)
certified 25 of 36 tiles. The first 24 form a continuous prefix through the
exact ratio-tiling boundary

`671092316289340128056086075422694163338407 /
238418579101562500000000000000000000000000`

(approximately `2.814765186581644456...`). One terminal tile from the exact
boundary

`4074311253101677331106625434260097734919149911558053987299757 /
1164153218269348144531250000000000000000000000000000000000000`

(approximately `3.499806717159296539...`) through `3.5` also passed. The 11
interior tiles between those boundaries stopped at depth 14 and 32,768 panels.
Their observed envelope lower margins ranged from
`-1.63743914468701138232158863988e-6` to
`-1.07096659996921571889512310197e-4`. Those signs are failures to certify an
upper bound, not lower bounds on `J-d` and not paper refutations.

This campaign covers the deliberately widened decimal interval
`[2.81476518658164, 3.49980671715930]`. Exact-rational import-time assertions
prove its left endpoint lies below the parent's exact open boundary and its
right endpoint lies above the parent's exact terminal boundary. Thus a full
PASS covers the complete unresolved interior with overlap into both already
certified neighboring tiles; endpoint rounding cannot leave a gap.

## Frozen instrument

`code/gate_direct_cached.py` and
`asrun/gate_direct_cached.py.asrun` are byte-identical, SHA-256
`ac01313afd3ec96b5ceaf74849c939b5cf1386a8716271b49223c200c53faa51`.
They are a domain-only refreeze of the owner-hardened panel-cache instrument
from campaign `20260831T111226Z_kg_direct_d4_highband_panelcache`; the inherited
`core.py` copies are byte-identical, SHA-256
`88c3cac0433105df58e9fbe2f597d0d21d94d2a67c950d09627233e70e5f4037`.

The machinery remains:

- direct Lemma D.4 monotone combined-hinge envelope
  `f(E,W;t)=(E+W-t)_+ + (|E-W|-t)_+`;
- full even unit disk via the exact 12x12 root cover (132 intersecting boxes),
  all compatible odd coefficients, and binary longer-axis subdivision;
- Arb precision 256, `S_MAX=8`, exact Gaussian panel masses, closed tail;
- fixed ratio `1.02`, panel ladder `1024 -> 4096 -> 16384 -> 32768` because
  every left endpoint is below 4;
- subdivision floor `1/16384`, depth cap 40, panel cap 32768, budget 260000;
- acceptance only when `d(c_U).lower()-U(B,c_L).upper()>0` on every leaf;
- exact panel-only Arb cache keyed by `(precision bits, panel count)`;
- exclusive/durable JSON checkpoints and refusal of any pre-existing result,
  progress, or temporary path.

Only the fixed domain, release environment name, output tag, provenance, and
reported command differ from the high-band cached runner. This campaign does
not alter the inequality or search another construction.

## Release boundary and decision

The source requires `KG_MIDBAND_RELEASED=1`, normal Python assertions,
niceness at least 10, and all five named numerical thread variables exactly
one before `run_band`. It also requires the exact fixed domain, budget, tag,
and campaign-local `logs/` directory.

Exactly one full-band run is registered. A PASS requires every generated tile
and every disk leaf to pass. Any open tile, budget stop, depth/floor stop, or
panel-cap stop leaves the interval **FAILURE TO CERTIFY / OPEN**. No second run,
sub-band probe, cell probe, or parameter search is authorized by this
pre-statement.

Release order is explicit: this run remains held until Main has completed the
untouched baseline `[3.5,4.083]` and the separately frozen high-band cached
rerun `[4.083,6.0]`.

## Registered command

```sh
cd cs/kg/campaigns/20260831T145732Z_kg_direct_d4_midband_panelcache/code
set -o noclobber
nice -n10 env KG_MIDBAND_RELEASED=1 PYTHONDONTWRITEBYTECODE=1 \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 \
  ../../../../.venv/bin/python gate_direct_cached.py \
  --c-lo 2.81476518658164 --c-hi 3.49980671715930 \
  --budget 260000 --out-dir ../logs \
  --tag band_2p814765_3p499807_panelcache \
  > ../logs/band_2p814765_3p499807_panelcache_stdout.log 2>&1
```

Reserved machine outputs are
`logs/band_2p814765_3p499807_panelcache_progress.json`,
`logs/band_2p814765_3p499807_panelcache_result.json`, and
`logs/band_2p814765_3p499807_panelcache_stdout.log`.

## Owner static refreeze

The first import-only release-refusal probe exposed a staging defect before
the release gate: python-flint `fmpq` does not parse decimal strings. It
aborted at the exact-domain assertion, before `run_band`, cache construction,
any envelope evaluation, or any output. The domain decimals are now encoded
as explicit numerator/denominator `fmpq` constants. This changes no endpoint.

The refrozen code/as-run SHA-256, superseding the earlier hash in this file, is
`319b0f28554e797978827d12c5d2ae71db2acfd1e32282863f4918eb279faa9a`.
Runner twins and core twins compare byte-identical. Source-text compilation
passes. The exact registered command without `KG_MIDBAND_RELEASED=1` now
reaches and passes the intended release refusal without calling `run_band`.
An import-only tiling check proves that 12 ratio-1.02 tiles cover the widened
domain, the panel ladder is `[1024,4096,16384,32768]`, and the cache remains
empty. The initial assumption of 11 generated tiles in the static probe was
wrong: the parent has 11 unresolved tiles, while restarting the ratio tiling
at the deliberately widened left endpoint generates 12. No logs or machine
result exist. Status remains **COMPUTE PENDING** in the registered release
order.

## Released production result

Main released the exact registered command only after the baseline
`[3.5,4.083]` and cached `[4.083,6]` runs completed and were inspected. This
run was the sole low-priority process with all five thread pools fixed to one.

The widened domain produced 12 ratio-1.02 tiles. Six certified and six
remained OPEN after 56,451 envelope evaluations and
`1314.9180881977081` recorded seconds. All six open leaves reached depth 22
and 32,768 panels. The smallest certified lower margin was
`+3.51851846436958939716186293650e-9`; the closest open upper-envelope lower
margin was `-4.52960421419799285256261953577e-8`, and the most negative was
`-6.84818667525655502520446653037e-6`.

The first five adjacent tiles certify the continuous interval from
`2.81476518658164` through
`3.107728208020454953573248`. Together with the parent campaign's certified
prefix, this extends the continuous `[1.75,3.5]` certification endpoint from
`2.81476518658164445646615405049` to
`3.107728208020454953573248`. A sixth tiny terminal tile
`[3.49980671715929099891398176911,3.49980671715930]` also certifies and
overlaps the parent's already certified terminal tile through `3.5`. The
remaining contiguous interior gap is therefore
`[3.107728208020454953573248,3.49980671715929653996950394581]`.

This is a strict extension of the certified subdomain but still **PARTIAL /
FAILURE TO CERTIFY** for the full band. Negative open margins remain failures
of the upper envelope at registered caps, not lower bounds on `J-d` or paper
refutations. Machine artifact SHA-256 values:

- result:
  `f97604585aba167e583e8d3c8859e25442db5e09c7a482516715af02d540592d`;
- progress:
  `0c5087efdaf7ace87c110334d0a2b59fd7eb4d28a2faec39cc3da797e8fa4d9d`;
- combined stdout/stderr:
  `08ec1472318a07ed94f23f4e8326770efce62c0bf0b5e188b80be9ce975910c2`.
