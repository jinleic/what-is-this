# Preregistration: Input-I gap closure

Cycle: `20260908T155849Z_f7f167`. Target: `math/e389`. One campaign run only.
Orchestrator and all research/review agents: user-selected GPT-6 Astra.

## Gate

Determine whether the variable interval and nonprimitive-frequency gaps in
`LOCALIZATION.md` section 2 / section 10 item 5 can actually be closed from the
cited inverse-prime sum bounds, with the quantifiers needed by `THEOREMS.md`
section 15.6. Do not assume that the existing description "finite work" is true.
The proposed frontier is theta < 17/32, not a solution of Erdős #389.

## Falsifiable acceptance criteria

1. Reconstruct the missing two gap statements from equations (40)--(43).
2. For the moving interval, give a valid uniform reduction with explicit prime
   subinterval lengths, coefficients, endpoint treatment, and error loss, or
   exhibit the precise inequality/hypothesis that fails. A fixed-arc theorem
   cannot silently be applied to prime-dependent arcs.
3. For nonprimitive frequencies, track gcd(h,q), the effective modulus, prime
   exclusions, the zero/major-frequency contribution, and theorem range. Treat
   noncoprime classes and the reduced-modulus parameter d explicitly. Do not
   assume uniform power cancellation for every nonzero h.
4. Check the proposed exponent only after both reductions, and state every
   restriction on m, M, a, d, and N. If only a restricted theorem survives, do
   not promote it to every class/every modulus.
5. Give a separate skeptical review, exact cited source statements, and a
   bounded arithmetic sanity check on mini-0. Finite checks validate only the
   checked identities/counterexamples, never an asymptotic theorem.

## Verdict mapping

- FROZEN-CERTIFIED: both gaps have complete arguments, all cited hypotheses
  checked at primary sources, and the advertised uniform frontier follows.
- FROZEN-NEGATIVE: a rigorous obstruction or counterexample invalidates a
  required claim or the advertised black-box reduction. This does not refute
  Erdős #389 or rule out every alternative method.
- FROZEN-INCONCLUSIVE: partial progress without a complete positive proof or
  decisive obstruction. Preserve the exact remaining lemma.
- CRASHED/REJECTED: execution/integrity failure, not a mathematical verdict.

## Scope, ownership, and bounds

Main opens the campaign, stages inputs, runs finite verification, integrates
accepted documentation, and alone freezes/closes/settles it. One OMP topic
owner investigates the moving-interval reduction and drafts the combined
proof; an independent Astra reviewer investigates nonprimitive frequencies
and adversarial classes in a disjoint artifact. No concurrent validation.
A bounded native Humanize2 Astra/Astra review follows the draft, with local
orchestration only and no local scientific execution.

All scientific computation, including tiny checks, must run only on
`ssh jinleic@mini-0.local`. The host's /usr/bin/python3 is an unusable developer
tools stub; /usr/bin/perl 5.34.1 is verified and requires no installation.
Use a new per-run directory `/Users/jinleic/institute/e389/<run_id>`; copy only
this gate's public inputs and scripts. No credentials, whole-workspace copies,
package installation, large sweeps, or heavy compute. Maximum finite-check
command duration: 60 seconds; total planned scientific CPU work: below three
minutes. Stop for approval rather than extend these limits.

The reproducible finite verification command, after substituting the minted
run id once, is:

```sh
ssh -o BatchMode=yes jinleic@mini-0.local /usr/bin/perl /Users/jinleic/institute/e389/<run_id>/verify_input_i.pl
```

The verifier must exercise elementary identities, reduced-frequency arithmetic,
and interval containment/endpoints on explicitly recorded bounded ranges.
Keep its source and actual remote output in the campaign before freeze. Keep
usage outside the sealed run and measure elapsed wall time rather than guess.

No unrelated targets, old sealed campaigns, budgets, lifecycle policy, or
integrity holds may be changed. Main owns the math PROGRESS entry and target
README changes. math/RESULTS.md is outside this target's ledger contract.
