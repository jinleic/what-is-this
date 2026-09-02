# AMENDMENT 3 — CF-3 plant mechanism corrected; row-XOR proven inert a priori

Recorded BEFORE running any plant (no plant has executed yet).

**CF-3 as registered in AMENDMENT 2 (row-XOR corruption) is mathematically
incapable of failing the alpha check.** One-line proof: in char 2 the
identity's both sides are F_2-ADDITIVE in each coordinate separately —
LHS(f_p + f_q) = Pi(f_p'+f_q') + Pi'(f_p+f_q) = LHS(f_p) + LHS(f_q), and
RHS(f_p + f_q) = G^2 (f_p+f_q)^2 = G^2 (f_p^2 + f_q^2) = RHS(f_p)+RHS(f_q)
by the Frobenius additivity (a+b)^2 = a^2 + b^2. So for rows built from a
passing F, EVERY row-XOR combination passes alpha identically. The
registered expectation ("at least one of the three pairs fails the alpha
check") is WRONG for ANY pair; running the registered suite would have
produced three vacuous "attempts" and either a false finding or a wasted
fallback. This is exactly the class of pre-run error rule 16 exists for.

**Corrected CF-3 (registered before any run):** replace row 0 of the m=11
anchor instance by 3 * f_0 (scalar multiplication by the F_4096 element 3,
c ∉ F_2). LHS(c f) = c LHS(f) = c G^2 f^2, RHS = G^2 c^2 f^2; equal iff
c^2 = c iff c in F_2. Since 3^2 = 5 != 3 in F_4096 (3*3 in GF(2^12): 3 =
x+1, squared = x^2+1 = 5; check in-run), the identity FAILS at every
support point where G(a)^2 f_0(a)^2 != 0 and the value relation is
exercised, i.e. wherever Y[0,a] = 1. Expected rejection certificate:
alpha value-check failure count = #{a in support : Y[0,a] = 1} recorded
exactly, first failing point recorded with lhs/rhs values. NOTE: this
corruption ALSO trips the delta binary check (lam*3f_0(a) need not be
binary) — a simultaneous delta trip is RECORDED as part of the
certificate; the plant establishes the alpha-level rejection direction,
which is its purpose. deg stays <= D (scalar multiply), so the degree
gate is NOT the rejector here — the alpha value check itself rejects.

**CF-1 and CF-2 run exactly as amended by AMENDMENT 2** (CF-1:
duplicate-support build, expect the distinct-support guard to fire;
CF-2: f_0 -> f_0 + Pi, expect deg_ok = false and ABORT with NO delta
verdict; NO value-level failure expected for CF-2 — see AMENDMENT 2's
correction).

— MceliecelM12, 2026-09-01, before any plant run.
