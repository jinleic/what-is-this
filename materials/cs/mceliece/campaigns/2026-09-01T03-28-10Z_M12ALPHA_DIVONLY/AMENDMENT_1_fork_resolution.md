# AMENDMENT 1 — fork resolution by measurement (recorded BEFORE the m=12 verdict compute)

Per pre_statement.md section 1 (fork) and section 4 (calibration + budget):

**What was measured (calib.json, in-process CPU only):**

1. m=11 anchor (11,2048,48,6211) full unmodified guard chain, including the
   FULL alpha grid: cpu 283.37 s / wall 283.62 s, vs frozen 285.91 s —
   agreement 0.9 percent. All guards true on this fresh build.
2. m=12 instance build (12,3488,64,16384): cpu 1783.47 s / wall 1785.41 s.
   Guards measured during build: beta_delta_nonzero true (pair (0,1)),
   gamma_k true, eps_gcd_const true, eps_maxdeg_is_D true.
3. m=12 alpha inner cost, unmodified instrument (`instance.py::_guards`
   alpha block op, vectorized over all k=2720 coordinates): measured
   0.1504 cpu s per support point (probes at support indices 0..4; value
   relation TRUE on all probes, which is expected for the true F).
   Projection: 0.1504 x 3488 = 524.6 cpu s for the FULL alpha grid.

**Resolution:** projected full m=12 alpha measurement = 524.6 cpu s, versus
the registered 10,800 cpu-s (3 h) verdict budget. The measurement path is
AFFORDABLE with >20x headroom even after adding the measured build cost
(1783.5 s) and the delta chain (order 100 s, O(n D)). Per pre_statement
section 4: "If the extrapolation instead fits the budget, the path stays
MEASUREMENT and no amendment is filed" — this file is the fork-resolution
RECORD (the pre-statement requires the resolution be recorded); the alpha
verdict at m=12 will be produced by the UNMODIFIED m<=11 instrument over
all n = 3488 support points, hence alpha at m=12 is MEASURED, not derived,
and carries label MEASURED (MACHINE-VERIFIED arithmetic), NOT
CITED-DEPENDENCY. No m<=11 label changes.

**Contingency (pre-declared):** if the verdict run hits the registered
budget ceiling mid-alpha (it should not, 20x headroom), the run STOPS at
the ceiling and freezes whatever prefix of the alpha grid completed; the
bounded statement then names the exact completed point count and carries
MEASURED-ON-PREFIX labeling, never a silent extrapolation.

— MceliecelM12, 2026-09-01, before the m=12 verdict compute.
