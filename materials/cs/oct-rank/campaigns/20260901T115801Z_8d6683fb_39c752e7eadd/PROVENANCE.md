# PROVENANCE — N2 campaign (run 20260901T115801Z_8d6683fb_39c752e7eadd)

- Minted by `python3 scripts/campaign.py init --gate n2_groebner_infeasibility
  --prereg n2_prereg.md` from cs/oct-rank, agent OctRankNext, AFTER the
  path-scoped prereg commit `028f5880ddde433653ef035dc44828039604193e`
  (message: "oct-rank: preregister N2 Groebner/elimination campaign …");
  prereg sha256 46e596bc80b91c290694bb639c1d50f4bda2a99050c31243b0a7353a3cf423a4,
  identical across workspace file / campaign.py hash / this dir's
  pre_statement.md; byte-identical copy made AFTER init, BEFORE any
  build/compute, per Main's standing order).
- N1 relationship: N1 (20260901T103338Z_06de66a3_d0b7bc2dc6a5,
  FROZEN-INCONCLUSIVE, FAILURE TO CERTIFY) did NOT make N2 logically
  unnecessary — documented in prereg section 1 (the N1 implication clause
  fires only on a witness; none was produced).
- Environment probes before commit (disclosed in the prereg Header):
  sympy groebner sizing only; NO msolve/flint build artifacts were created
  for this campaign before init. msolve 0.8.0 is built INSIDE this run dir
  (BUILD_LOG), from the source with sha256
  319ba0de67dca967dea40cb6e4dacf44eab387c2f0c2416b71842c11affbadfd.
