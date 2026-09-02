# N2 pre-statement — upward Gröbner/elimination attack on rank((L_1, L_i, L_j))

Campaign: to be minted by `python3 scripts/campaign.py init --gate
n2_groebner_infeasibility --prereg n2_prereg.md` from
/Users/jinleic/jinleic-workspace/cs/oct-rank AFTER this file is committed
(path-scoped commit; Main's strict order). Owner agent: `OctRankNext`.
No N2-parameter-dependent computation has been run. Tool INSPECTION only was
performed before this commit (allowed): (a) venv Groebner engine sizing
(sympy groebner: 3-var smoke 0.02 s; a generic 6-var quadratic system took
114.6 s in lex mod 32003 and 0.09 s in grevlex — orders of magnitude short
of 247-variable elimination); (b) cmake 4.2.3 / Apple clang 17 confirmed
present; (c) the exact msolve version/build is pinned below and built only
AFTER this commit. Disclosed precommit probes (no evidence value): the
disclosed msolve-0.8.0 source download (hash below) and the abandoned
flint-3.6.0 /tmp build per Main's steering — nothing from either entered
any campaign path.

## 0. Question and standing state (frozen facts, re-read first-hand)

Decide UPWARD: rank(T) >= 14 for T = (L_1, L_i, L_j) on the Cayley-Dickson
octonions, equivalently for the rank-equivalent certification target
TF = blockdiag(tau, tau) (two frozen entrywise diagonal-similitude
identities, rf_krawczyk.py / s3_routeF3_global anchors).

Equivalently: certify that the rank-13 CP system of TF has NO real solution:
  g(x) = CP13(x) - TF = 0,  CP13[p,b,c] = sum_r a[p,r] b[b,r] c[c,r],
  192 equations, 247 unknowns (a: 3x13, b: 8x13, c: 8x13),
all coefficients integers (T entries in {-1,0,1}).

Combined with the frozen floor rank(T) >= 13 [MACHINE-VERIFIED chain:
1 peel + pencil 12; Route A prov closed the only live >13 route at 12 < 13],
a certified no-solution statement would give rank(T) = 14 EXACTLY.

Frozen context consumed:
- N1 campaign (20260901T103338Z_06de66a3_d0b7bc2dc6a5,
  FROZEN-INCONCLUSIVE): 40 seeds, 0 admitted; verdict FAILURE TO CERTIFY —
  no feasibility statement either direction.
- Certified-local Route-F facts: the two frozen candidates have certified
  no-root boxes (rho <= 1e-6) — local only, never global; nothing here
  re-litigates them.
- Route E (S3 campaign, 2026-08-30) BUDGET-FAILed on 312 unknowns by
  pre-registration; only a necessary-only linear diagnostic ran. N2 is the
  first FULL-system exact elimination attempt in this repository.

## 1. Logical relationship to N1 (assignment item 2, documented)

N1 attacked DOWNWARD (produce a 13-witness). Its failure produces no
evidence about feasibility. N2 attacks UPWARD (certify 13-infeasibility =>
rank 14). N1 did NOT make N2 logically unnecessary; the N1 prereg's
implication clause (fires only on a witness) is recorded there and did not
trigger. Accordingly N2 opens as a distinct lifecycle campaign. If N2
certifies infeasibility, the 13-witness question becomes moot
(rank = 14 EXACTLY); if N2 budget-fails, {13, 14} stays OPEN exactly as
before. Published window 18 <= R_R(T_O) <= 25 untouched by every outcome;
no published work is refuted.

## 2. Gauge and soundness policy (fixed NOW; the subtle point)

The CP parametrization carries a 13-dim scaling torus
(a,b,c) ~ (aD^-1, Db, Dc), D diagonal. A NEGATIVE certificate must cover
ALL decompositions, so the gauged system must be equivalent to the original
on the relevant solution set. THE SOUND REDUCTION, fixed here:

- Every term of any rank-13 decomposition is NONZERO (a term with
  a[:,s]==0 or b[:,s]==0 or c[:,s]==0 contributes nothing and could be
  dropped, contradicting the frozen floor 13; note this uses the frozen
  floor as an INPUT — recorded).
- Sound gauge elimination (THE fixed formulation; the subtle point was
  worked out BEFORE commit): the naive normalization c[0,s] = a[0,s] is
  NOT WLOG (scaling sends a[0,s] -> a[0,s]/d_s and c[0,s] -> c[0,s]*d_s,
  so one-sided zeros can never meet). The sound reduction appends, per
  term s, the TWO quadratic equations in 2 fresh unknowns (d_s, g_s):
      a[0,s] * d_s - c[0,s] * g_s = 1,
      c[0,s] * d_s + a[0,s] * g_s = 0.
  As a linear system in (d_s, g_s) it has the 2x2 coefficient matrix with
  determinant a[0,s]^2 + c[0,s]^2 and RHS (1, 0); over R it is solvable
  IFF a[0,s]^2 + c[0,s]^2 != 0, i.e. (a[0,s], c[0,s]) != (0, 0). Over
  F_p the same linear-algebra criterion holds: solvable iff
  a[0,s]^2 + c[0,s]^2 != 0 in F_p (for p = 32003, 2147483647 — both
  3 mod 4 — this is equivalent to (a[0,s], c[0,s]) != (0,0); for
  p = 65537 = 1 mod 4 there are nonzero isotropic pairs, so a {1}
  basis at 65537 certifies infeasibility of a SUPERSET of the
  all-nonsingular family and therefore still implies: no solution in
  which every term avoids BOTH the singular pattern AND the isotropic
  pattern — recorded per-prime with this exact wording).
  SEMANTICS (sound by construction): the extended system has a solution
  IFF the original 247-var system has a solution in which EVERY term s
  satisfies (a[0,s], c[0,s]) != (0, 0). The residual uncovered family —
  solutions with a[0,s] = c[0,s] = 0 for some s — is handled by the
  PRE-REGISTERED SECOND SYSTEM E1: fix a[0,1] = 0 AND c[0,1] = 0 (two
  linear equations; immaterial which term index — column permutation of
  the 13 terms is a solution-set bijection, recorded as the second hand
  lemma), plus the same 2-equation pair for terms s = 2..13 (24 fresh
  unknowns, 24 equations). E1 covers exactly the solutions with term 1 in
  the uncovered family and all other terms non-singular at (a[0,s],
  c[0,s]); iterating for k >= 2 zero-terms is NOT needed: E1 already
  covers any solution having AT LEAST ONE (a[0,s], c[0,s]) = (0,0) term —
  relabel that term to index 1 (permutation lemma). Joint semantics
  (subject to the section-3 characteristic-zero protocol, which governs):
  E and E1 both certified infeasible by the section-3 protocol =>
  rank >= 14; with the frozen floor 13 replayed in-run: rank(T) = 14
  EXACTLY [MACHINE-VERIFIED (i) / per-protocol strength (ii)].
  Sizes: E: 273 variables (247 + 26), 218 equations (192 + 26).
  E1: 275 variables, 220 equations.
- Ordering/equation bookkeeping: elimination proceeds on the extended
  system; infeasibility certificates are Groebner bases, not hand logs.

## 3. Modular primes and reconstruction policy

Primes (fixed list, tried in order): 32003, 65537, 2147483647.
Characteristic-zero policy (the standard reduction hazard discharged
BEFORE commit): a Q-point of any of these systems does NOT necessarily
reduce to an F_p-point at a FIXED prime (p may divide the cleared
denominator). The pre-registered decision protocol therefore makes the
rank-14 claim ONLY in one of the two certificate-strong ways:
  (i) msolve 0.8.0's Q-mode (rational Groebner with certified
      denominators) certifies {1} over Q within budget; or
  (ii) the F_p runs return {1} at ALL THREE fixed primes AND the
      "primitive-vector" argument below is checked in-run exactly:
      a hypothetical Q-point, cleared to a PRIMITIVE integer vector
      (gcd of all coordinates = 1; the appended unknowns scale along
      since the appended pairs are equations in differences of products
      with RHS — after homogenizing the RHS by one extra variable h and
      the equation h*e_s - 1 replaced by h*e_s - z^2 with z the
      homogenizing variable... the exact homogenized form is FIXED as:
      dehomogenize by h; a primitive Z-point of the dehomogenized
      system with h=0 corresponds to solutions "at infinity" which are
      separately excluded by certifying {1} on the homogenized system's
      h=0 slice), reduces to an F_p-point for EVERY prime NOT dividing
      the gcd — primitivity guarantees AN INFINITE set of such primes,
      but no specific prime is guaranteed good; HENCE clause (ii)
      alone cannot be airtight at fixed primes.
  ADOPTED PROTOCOL (final): the campaign attempts (i) first at tau
  scale (calibration), then for the main system; if the Q-mode exceeds
  budget, the F_p runs are certified and frozen as mod-p infeasibility
  certificates [COMPUTATIONAL-EVIDENCE, each {1} basis a checkable exact
  artifact over F_p], and the honest campaign outcome is: the {13,14}
  question stays OPEN with the strongest exact statements being the
  frozen mod-p certificates. This meets the assignment requirement:
  every negative claim either has a checkable exact certificate ({1}
  basis over F_p at a fixed prime, or over Q) or is labeled OPEN.
  NO claim of rank >= 14 is made from mod-p runs alone.

## 4. Tool, version, and build (post-commit only)

- Solver: msolve 0.8.0 (F4-style GB over F_p; the binary and library are
  built from the source tarball sha256
  319ba0de67dca967dea40cb6e4dacf44eab387c2f0c2416b71842c11affbadfd — the
  same artifact Sage 10.7 pins), linked against GMP 6.3.0 (homebrew) and
  the flint 3.6.0 ALREADY BUILT AND DISCLOSED pre-commit; the binary build
  is repeated from scratch INSIDE the campaign to keep the artifact chain
  clean. Build recorded in BUILD_LOG inside the run dir.
- Justification vs existing backends (why the venv cannot serve): sympy
  groebner measured (inspection above) ~2 min for SIX generic quadratic
  variables; the extended system has 273 variables with degree-3
  equations; flint provides no multivariate GB; python-flint fmpq
  Groebner equivalents do not exist at this size. msolve is the only
  in-policy exact solver with published results at >=100 variables.

## 5. Stages and budgets (fixed NOW; every stage outcome reported)

S1 (calibration, cap 30 min CPU each): run the pipeline on
  (a) tau r=7 CP system (48 eqs, 77 unknowns; KNOWN feasible):
      expected NON-{1} (feasible) outcome;
  (b) tau r=6 CP system (KNOWN infeasible — frozen rank(tau)=7): the
      pipeline MUST certify {1} (or equivalent) within cap;
  (c) tau r=6 with ONE target entry perturbed by +1 (smoke only; outcome
      recorded, not adjudicated — same semantics as K-3 below).
  If (b) fails to certify infeasibility within cap: BUDGET-FAIL for the
  whole campaign (the instrument cannot scale even to the known case; no
  claims). If (a) yields {1}: FALSE INFEASIBILITY — quarantine;
  verdict INVALID INSTRUMENT; nothing else reported.
S2 (the real attack, cap 6 h CPU per prime, 24 h campaign CPU cap total,
  including S1 and the build): extended E-system (273 vars) and E1-system
  (275 vars) over F_p for the three fixed primes in order; clause-(i)
  Q-mode attempt first within the same cap. The rank-14 verdict is made
  ONLY per the section-3 protocol ((i) Q-mode {1}, or (ii) with its
  stated strength limits — mod-p-only results are frozen as
  COMPUTATIONAL-EVIDENCE certificates and {13,14} stays OPEN). Per-prime
  records: GB degree profile, timeout status, msolve log.
  Exhaustion of all caps without a protocol-grade certificate:
  BUDGET-FAIL; {13,14} stays OPEN.

## 6. Controls (through the IDENTICAL pipeline; run BEFORE S2, after build)

- K-1 (accept known-feasible): tau r=7 system — must NOT certify {1}.
- K-2 (reject known-infeasible): tau r=6 system — must certify {1}.
- K-3 (one-coordinate perturbation smoke control): the tau r=7 CP system
  with ONE target entry (the I-slice (1,1) entry) perturbed by +1, at
  r = 7. Fixed semantics: SMOKE ONLY — must run to completion without a
  crash; its outcome (either direction) is RECORDED, not adjudicated: a
  rank statement for the perturbed tensor is NOT a machine-verified frozen
  fact, so neither outcome is a calibrated signal. The load-bearing
  calibrated controls are K-1/K-2/K-4.
- K-4 (gauge soundness): the extended tau r=7 system (26 added vars/eqs)
  must NOT certify {1} (feasibility survives the extension), and the
  extended tau r=6 system must STILL certify {1} (infeasibility survives
  the extension) — the soundness lemma live-tested at tau scale.

Quarantine semantics: any K-1 false-infeasibility or K-4 extension
breakage => verdict INVALID INSTRUMENT; main S2 outcomes not reported as
evidence.

## 7. Adjudication vocabulary (fixed now)

- "rank(T) = 14 EXACTLY": ONLY via the section-3 protocol: clause (i)
  Q-mode certified {1} over Q, or clause (ii) with its stated strength
  limits; in BOTH cases with K-1..K-4 clean + frozen floor 13 replayed
  in-run (the S3/lower13 chain replay script, exit 0). Clause-(ii)
  success labels the rank statement with exactly the strength its
  certificates carry (MACHINE-VERIFIED for Q-mode; per-prime
  certificate strength otherwise, and OPEN where the protocol says so).
- "BUDGET-FAIL": instrument cannot scale (S1/S2 caps) or no
  protocol-grade certificate; no rank claims; both frozen sentences
  verbatim in the verdict: "Absence of a 13-witness is NOT evidence for
  14; absence of an impossibility argument is NOT evidence for 13."
- status.json mapping: protocol-grade rank-14 success ->
  FROZEN-CERTIFIED; S1/S2 budget exhaustion with clean controls ->
  FROZEN-INCONCLUSIVE (BUDGET-FAIL);
  control break -> FROZEN-INCONCLUSIVE (INVALID INSTRUMENT); crash ->
  FROZEN-INCONCLUSIVE (CRASHED).

## 8. Reproducibility of a negative resource result (assignment item 3)

If BUDGET-FAIL: the run dir records (a) the EXACT msolve command lines with
all flags, (b) per-prime stage: GB degree profile / max degree reached,
time, memory, exit status, (c) the built binary's sha256, (d) the full
input system files (the extended system serialized in msolve input format
with checksums), (e) the CPU/wall accounting. Any future re-run can replay
the identical commands against the identical input files.

## 9. Lifecycle and freeze plan

Commit this file (path-scoped) -> init --gate n2_groebner_infeasibility ->
byte-identical copy into the minted run dir + PROVENANCE.md (source commit
+ sha256) -> build msolve inside the run dir (BUILD_LOG) -> controls
(S1: K-1..K-4 through tau) -> S2 primes -> freeze -> close with the
section-7 mapping. Pristine discipline: PYTHONDONTWRITEBYTECODE=1 in every
script; no __pycache__ in the run dir at freeze.

## 10. Rule-7 prospective scope

Covered: T = (L_1, L_i, L_j) (and its rank-equivalent TF) at exactly
rank 13, via the extended-gauge 273-variable integer system, through
msolve 0.8.0 at the three fixed primes, with the tau-scale controls; exact
integer/F_p arithmetic end to end; CPU caps as fixed.
NOT covered, hence not decided by anything here: any rank other than 13;
T_O or any other slice family; border rank; complex-rank statements beyond
the recorded char-0 lifting argument; Koiran commuting-extension systems;
numerical-certification instruments (N1's domain); any positive claim
about rank 13 feasibility; the published window
18 <= R_R(T_O) <= 25 (untouched).
