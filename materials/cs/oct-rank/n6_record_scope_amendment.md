# N6 precompute amendment — collateral record scope and exact ladder

Run `20260904T052006Z_e5f2adce_0e7798f95678`; gate
`n6-qi-real13-descent`. Bound after the corrected historical formula audit and
before any N6 continuation, plant, target, or certificate computation.

## Exact N6 radius ladder

N6 preregistration wrote the copied Route-F ladder compactly as
`1,10^-1,...,10^-12`. To remove every ellipsis/heritage ambiguity, N6 binds
exactly these 15 rational radii, in this order:

```
1,
1/10,
1/100,
1/1000,
1/10000,
1/100000,
1/1000000,
1/10000000,
1/100000000,
3/1000000000,
1/1000000000,
3/10000000000,
1/10000000000,
1/100000000000,
1/1000000000000.
```

The driver, positive control, corrupted same-center control, target
certificate, and independent replay must all use that exact full list. Unlike
the historical Route-F preregistration, whose registered ladder started at
`1/100`, N6's committed preregistration explicitly started at `1`; this
amendment makes the N6 list fully explicit rather than changing it.

## Historical collateral corrections not inherited by N6

1. Historical Route-F/N1 source evaluated 15 radii, but the Route-F
   preregistration registered only the 13-rung suffix beginning at `1/100`.
   The extra `1` and `1/10` rows were unregistered. They were non-decisive;
   every first containment/exclusion rung used in the corrected classification
   lies in the registered 13-rung suffix. The corrected formula audit evaluated
   all 15 to inventory actual execution, not to retroactively register the
   first two.
2. N1 described PWRNG as using the “best N1 candidate,” although controls ran
   before the N1 sweep. Code actually loaded Route-F's frozen selected-CP
   candidate CSV, whose hash was omitted from N1 `seed_input_hashes`. PWRNG is
   therefore an unregistered-center, box-local diagnostic, not a valid N1
   preregistered control. Corrected replay establishes only the diagnostic's
   actual fixed-box classification. N6 instead uses the explicitly bound
   planted center and its same-center one-entry corruption from
   `n6_precompute_amendment.md`; it consumes no N1 PWRNG semantics.
3. N1 `VERDICT.md` calls R18 the best seed, while `n1_results.json` identifies
   R22. The numeric minimum in the root results record is unaffected. N6
   consumes neither seed.
4. Route-F polished CP and EXT but called `certify` only on selected CP. The
   replay's decimal center is a second rational parsing of that same selected
   CP CSV, not EXT. Claims of certified boxes around “both” polished candidates
   are retracted: EXT received no Krawczyk adjudication and its center was not
   serialized. N6 consumes only its own candidate and separately serialized
   continuation/certification charts.

The automated corrected-audit string `PASS-CLAIMS-UNCHANGED` is interpreted
narrowly: the five actually executed Route-F classifications and their shared
N1 formula classifications remain the same under the corrected Taylor bound.
It does not excuse the registration, input-provenance, best-seed-label, or
both-candidate scope defects above. Frozen files remain unchanged; README and
provenance carry the append-only correction.
