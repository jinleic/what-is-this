# AMENDMENT 2 — counterfactual plant suite, corrected BEFORE the verdict compute

Pre-statement section 3 registered two plants. Working the CF-2 arithmetic
through BEFORE running it exposed a defect in the pre-registered
expectation (not in the instrument), and this amendment fixes the suite
per rule 16 (declared before the verdict compute; the plants have not run
yet; nothing result-bearing has been computed past the calibration
already recorded in AMENDMENT 1).

**Defect found in the CF-2 registration.** The registered corruption
"F -> F + Pi*Z^0" (adding the full-support vanishing polynomial) is — by
direct computation at any support point a:

    Pi(a) = 0,  Pi'(a) = prod_{l != i}(a_i + a_l)  (nonzero),
    (F+Pi)(a) = F(a),  (F+Pi)'(a) = F'(a) + Pi'(a).

    LHS = Pi(a)(F+Pi)'(a) + Pi'(a)(F+Pi)(a)
        = 0*(F'(a) + Pi'(a)) + Pi'(a)F(a) = Pi'(a)F(a);
    RHS = G(a)^2 (F+Pi)(a)^2 = G(a)^2 F(a)^2.

So the value-level alpha check at support points sees LHS = Pi'(a)F(a)
and RHS = G(a)^2 F(a)^2, which are EQUAL for the TRUE F by the identity
itself (the identity says Pi F' + Pi' F = G^2 F^2 on F, and the Pi F'
term vanishes at support points). The plant is therefore INVISIBLE to the
value-level check — exactly the failure mode ADDENDUM 2 clause (i)
warns about verbatim: "f and f + Pi agree at every support point yet
differ, so without deg <= D the agreement-at-n-points argument
collapses." My original registration expected a VALUE-level rejection;
the correct rejection for this plant is the DEGREE GATE (deg(F+Pi) = D+1
> D, so deg_ok = false and ABORT fires producing NO delta verdict). The
pre-statement's own text actually says this ("Required: the degree gate
flags deg_ok = false and NO delta verdict is produced"); the incorrect
part was expecting additionally that the identity FAILS on some support
point value-wise. It does not and cannot. Corrected expectation: the
plant is accepted as REJECTED iff (deg_ok == false) AND (the ABORT
semantics visibly produce NO delta verdict). This is itself a
demonstration that the guard is load-bearing, which is the point of the
plant.

**Additional plant registered now (CF-3), giving the suite an
alpha-level rejection:** F -> F with row 0 replaced by 7*f_0 (scalar
multiply row 0 by the field element 7). This preserves delta
(lam*(7f_0)(a) = 7*Y[0,a], no longer binary — so it is ALSO expected to
trip the delta binary check; to make it a PURE alpha plant we instead
use: F -> F with row 0 replaced by f_0 + G*(Z^{D-2t-1}) — no; simpler
and cleanest: row 0 -> f_0 XOR f_1 (coordinate-wise XOR of two rows,
both of degree <= D). Then lam*F'(a) entries stay binary (sum of two
binary columns is binary), delta still passes pointwise — but the
IDENTITY Pi F' + Pi' F = G^2 F^2 is F_2-LINEAR in the coordinate
selector ONLY through the SQUARED RHS being quadratic: a single row XOR
breaks the identity exactly when the cross term f_0*f_1 fails to cancel.
Concretely the corruption must be checked to fail; if it happens to
pass (possible: some pairs are compatible), the plant is re-drawn with
the NEXT row pair (row 1 XOR row 2, then row 0 XOR row 2, ...) in the
fixed order until one FAILS or the list of three pairs is exhausted;
the attempt order is fixed here in advance and every attempt is logged.
Expected: at least one of the three pairs fails the alpha check value-
wise at some support point; that failure is the alpha-level rejection
certificate. If (adversarially) all three pass, THAT is itself a
finding: the alpha identity is not F_2-coordinate-sensitive at row level
on this instance; report it and rely on CF-2's degree-gate rejection as
the plant evidence of record.

Summary of the corrected suite:

- CF-1 (delta plant, unchanged): duplicate-support m=11 build; expect
  the build's distinct-support assert to fire OR lagrange_unit_ok=false.
- CF-2 (degree-gate plant, expectation corrected above): F -> F + Pi;
  expect deg_ok=false, ABORT recorded, NO delta verdict.
- CF-3 (alpha plant, NEW, fixed attempt order (0,1), (1,2), (0,2)):
  row-XOR corruption; expect >= 1 pair to fail the alpha value check;
  log every attempt.

— MceliecelM12, 2026-09-01, before the m=12 verdict compute.
