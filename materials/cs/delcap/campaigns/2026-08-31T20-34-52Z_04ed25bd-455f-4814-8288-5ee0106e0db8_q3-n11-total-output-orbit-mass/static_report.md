# q=3,n=11 static release report

## Verdict

**PASS / ZERO BA ITERATIONS / ZERO CERTIFICATE ROWS.** The final frozen runner SHA-256 is `08ecbc493e1c5c674f04ea4a0d0c9c6efd9ac2254fe9d6513db90933b3cdf114`; the manifest SHA-256 is `198cb863a104857dfd7b067e2035ebca9b08c6c0400e043c01508d95666a183a`. All seven static checksum entries pass.

The current-source static guard completed at 2026-08-31T21:10:33Z in `7.517722 s` wall / `7.512062 s` CPU with maximum recorded RSS `944,291,840` bytes (`COMPUTATIONAL-EVIDENCE`). It wrote only `static_guard.json` and its captured stdout; resume validation remained exactly zero rows.

## Machine controls

`MACHINE-VERIFIED` controls:

- exact row box: four cells, only `q=3,n=11`, in deletion-probability order `1/2,1/5,1/10,1/20`;
- independent Burnside counts: 14,884 input orbits and 22,450 output orbits;
- exact full structure: census tuple `(14884,22450,6148309,3573542,9721851)`, input/output orbit-size sums `177147` and `265720`;
- 6,480 exhaustive small-instance joint value/reversal equivariance cases;
- position plant breaks with subsequence counts 1 versus 2;
- planted output-orbit split has five failed generator equalities, is detected, and is rejected before KL;
- exact orbit-total uniform control passes 1,092 expanded generator checks;
- zero-support dual returns infinity;
- negative-width plant is rejected;
- immutable resume validates zero rows.

## Independent review and owner hardening

A read-only `claude-fable-5` release audit returned `PASS` on runner hash `0fb61277...`, finding no correctness, provenance, or release-safety blocker. It identified four optional hardenings. Before static execution, owner applied three:

1. enforce `PYTHONDONTWRITEBYTECODE=1` before any campaign dependency import;
2. assert Darwin before interpreting `ru_maxrss` as bytes;
3. forbid standalone table rendering until a complete terminal summary exists.

The fourth suggestion—derive display width from the already-outward binary rationals—was not applied: the certificate width is already formed from the Arb upper/lower endpoints, moved outward by one binary64 ULP, converted to an exact rational, and machine-asserted outside the exact endpoint difference. The final post-hardening source was compiled from text, rebound in the manifest and static ledger, revalidated at zero rows, and then exercised by the passing static guard above.

The release key, niceness, thread-pool, assertions, no-bytecode, static-prerequisite, terminal-artifact, resource, atomic-write, and immutable-resume guards are active. On the reviewed pre-static candidate, unreleased production, released production without niceness 10, and released production before static acceptance each refused before BA or row output. The final hash changed only to add the three hardenings above; its new pre-import no-bytecode refusal was re-smoked directly.

## Release boundary

Production remains limited to the exact command frozen in `pre_statement.md`. The first `d=1/2` production row is the startup anchor and canonical row zero; it will not be duplicated by a separate anchor run. No result claim exists at this stage.
