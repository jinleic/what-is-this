#!/usr/bin/env python3
"""Shiu (J. reine angew. Math. 313 (1980) 161-170) -- effective per-class
margin at the e389 Input-II corner u = 2.

Source audit (this session): full text from GDZ page images, printed
pp. 161-169 = GDZ iiif PPN243919689_0313:00000165..173.  Wright's tex quotes
Shiu's Lemma 2 as (tex l. 361): Phi <= y/(phi(k) log z) + z^2 -- the printed
Lemma 2 (p. 164) is an honest INEQUALITY with constant 1:

    Phi(x,y,z;k,a) <= y/(phi(k) log z) + z^2      (z >= 2, k < y <= x)

Notation map (verified on the scans): Shiu's (alpha, beta) = (modulus range
k < y^(1-alpha), length range x^beta < y <= x) correspond to (alpha, kappa)
in Wright's tex.  For f = 1_S: A1 = A2 = 1 exactly (f(p^l) = 1, f(n) <= 1).

Corner: x = 2N, y = Q = sqrt(4N) = 2 sqrt N, k = M, u = log x/log y = 2.
Shiu's hypothesis (2.2) (k < y^(1-alpha), any alpha > 0) implies M < y:
the theorem covers exactly M in [sqrt(2N), 2 sqrt N), silent above 2 sqrt N.
Within the band, (5.1) z = y^(alpha/10), and since alpha <= alpha_max =
ln(y/M)/ln y, the sieve half-level satisfies ln z' = ln(y/M)/20 -- the depth
collapses across the band and hits (ln 2)/40 at the frontier M = sqrt(2N).

The per-class (45)-meter.  Class-I bound per window with the b-sum
B(z,M) = sum_{b <= z, (b,M)=1} f(b)/b:

    survivors(window) <= y B(z,M)/(phi(M) ln z)  +  z'^(3/2) M/y (+ edges)

amortized over N/y windows (Hildebrand subadditivity, as in the Wright
audit); the junk is additive/nonnegative (only helps a floor for R), no
o(1) claim is made:

  R(alpha) := (amortized bound)/(N/M) = B(z,M)(M/phi(M)) / ln z'
            = 2 B(z,M)(M/phi(M)) / ln z.

In the band z <= 2^(1/20) < 2 so ONLY b = 1 contributes: B = 1 EXACTLY.
Hence in-band

    R_min = (M/phi(M)) / ln z'_max = (40/ln 2)(M/phi(M)) / log_2(4N/M^2),

minimized at alpha = alpha_max; frontier M = sqrt(2N): log_2(4N/M^2) = 1,

    R_min >= (40/ln 2)(M/phi(M)) = 57.7078016...(M/phi(M)) >= 57.708.

Contrast (balanced sieve, junk z'^2 balanced against the main term):
ln z' <= (1/2) ln(y/M) + O(ln ln), R >= 2(M/phi)/ln(y/M) >= (4/ln 2)(M/phi)
= 5.7708(M/phi) at the frontier.

Universal deep-depth statement (M depth -> infinity, exact Mobius identity
B(z,M) = sum_{d|M} mu(d)/d H_{floor(z/d)}, upper-friendly Mertens: with the
(phi(M)/M) average density and + O(ln M) divisor correction,

    B ~ (phi/M)(ln z + gamma + sum_{p|M} ln p/(p-1)),

so the amortized cascade has R -> 2 (one factor phi/M from the density, one
ln z'/ln z' = 1 from B ~ density * ln z; the factor 2 = one from the
amortized y/(phi ln z) main term, one from B's density factor phi/M).  The
cascade NEVER proves R < 1 at any depth: its realized constant is 2 at
unbounded depth, not 1.  (Earlier "R -> 1^+" was wrong: it ignored the b-sum
density.  Corrected here.)

e_req-analog: (45) needs R < 1 <=> ln z' > M/phi(M): Lemma-2's per-class
shape supplies at best ln z' ~ (1/2) ln(y/M) (balanced), so the deficit is
ln(y/M) < 40 (M/phi) at Shiu's depth (frontier: 57.708(M/phi) shortfall)
and ln(y/M) < 4(M/phi) at balanced depth (frontier: 5.7708(M/phi)).

Reconciliation with the ledger: "Shiu margin e^{B1} = 1.299 per rho(2) unit
(4.23x truth)" normalizes the printed exp-form; two corrections: (a) the
proof-realized idealized-depth constant is 2(M/phi(M)) (theorem-form is
M/phi-fold looser than its own proof), (b) at every finite N in the band
the realized margin is >= 57.708(M/phi(M)) because (2.2) collapses the
sieve DEPTH, not the constant.

Korolev 2019 cross-check: NO / orthogonal (stdout).  Three-way split:
(i) printed-literal, (ii) Arb-certified derivations, (iii) EMPTY for
f = 1_S (the sole external link (4.1)-PNT is effective, HUMAN-AUDITED).

Reproduction: ./.venv/bin/python -I -B e389/shiu_effective_margin.py from
math/.  Byte-stable stdout; every [CERT]/[EXACT] decimal is the common part
of matching Arb lower/upper balls at prec >= 300; [PUB] flags published-
bound dependencies; [NOTE] lines are static prose.
"""

import sys

import flint
from flint import arb, fmpq, fmpz

flint.ctx.prec = 300


def dec(x: arb, digits: int) -> str:
    """Certified decimal expansion (common part of lower/upper balls)."""
    m = fmpz(10) ** digits
    while True:
        lo = (x.lower() * m).floor().unique_fmpz()
        hi = (x.upper() * m).floor().unique_fmpz()
        if lo == hi:
            neg = lo < 0
            a = -lo if neg else lo
            s = str(a).rjust(digits + 1, "0")
            out = s[:-digits] + "." + s[-digits:]
            return ("-" + out) if neg else out
        digits -= 1
        if digits < 8:
            raise AssertionError("dec: balls do not agree to 8 digits")


def phi_exact(n: int) -> int:
    """Euler phi by trial division (exemplar moduli only, manageable)."""
    m, result = n, n
    p = 2
    while p * p <= m:
        if m % p == 0:
            result -= result // p
            while m % p == 0:
                m //= p
        p += 1
    if m > 1:
        result -= result // m
    return result


def moebius(n: int) -> int:
    """Mobius mu by trial division (n <= 40 here)."""
    if n == 1:
        return 1
    m, c = n, 0
    p = 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            c += 1
        p += 1
    if m > 1:
        c += 1
    return 1 if c % 2 == 0 else -1


assert [moebius(d) for d in [1, 2, 3, 4, 6, 8, 30]] == [1, -1, -1, 0, 1, 0, -1]

# ------------------------------------------------------------------ constants
B1 = arb("0.2614972128476427837554268386086958590516")   # Meissel-Mertens
rho2 = 1 - arb(2).log()
ln2 = arb(2).log()
gamma_e = arb.const_euler()
e_B1 = B1.exp()

W = 78
out = []
out.append("=" * W)
out.append("Shiu 1980 (Crelle 313) -- effective per-class margin at the e389")
out.append("Input-II corner u = 2  (x = 2N, y = Q = sqrt(4N), k = M)")
out.append("=" * W)
out.append("")
out.append("[EXACT] source constants (printed pp. 161-169, GDZ images 165-173):")
out.append("  Lemma 2 (p.164): Phi(x,y,z;k,a) <= y/(phi(k) log z) + z^2   [const 1]")
out.append("  Lemma 1 (p.164): Psi(x, log x log log x) <= exp(3 log x/(log log x)^(1/2))")
out.append("  (4.1): sum_{p<=y} 1/log p <= 2y/log^2 y   [PUB: Rosser-Schoenfeld/Dusart]")
out.append("  Lemma 4 (pp.165-166): ... exp(-r log r/10 + 2 A1 r^(1/4)), r <= log z/log log z")
out.append("  (5.1) z = y^(alpha/10);  classes I-IV;  assembly p.169: z^2 < y/(phi(k) log z)")
out.append("  Theorem 1 (p.163) hypotheses (2.2): k < y^(1-alpha), x^beta < y <= x")
out.append("  (Shiu's beta = Wright's kappa; notational map verified on the scans)")
out.append("")
out.append("[EXACT] f = 1_S: A1 = A2 = 1 (f(p^l) = 1, f(n) <= 1)")
out.append("[EXACT] Lambda3 = sum_p sum_{l>=2} f(p^l)/p^l <= sum_{n>=2} 1/(n(n-1))")
out.append("        = lim_K (1 - 1/K) = 1   (exact telescoping)      -> Lambda3 <= 1")
K = 1000
tel = sum(fmpq(1, n * (n - 1)) for n in range(2, K + 1))
assert tel == 1 - fmpq(1, K)
out.append("[EXACT]   partial-sum identity to K=" + str(K) + ": 1 - 1/K  (verified exactly)")
Nz = 20000
zeta32 = arb(fmpq(0))
for n in range(1, Nz + 1):
    zeta32 += arb(1) / (arb(n) ** fmpq(3, 2))
tail = 2 / (arb(Nz) - 1).sqrt()
A3 = 2 * (zeta32 + tail - 1)
out.append("[CERT] A3 ((2.3)-constant for 1_S) <= 2(zeta(3/2)-1) = " + dec(A3, 12))
out.append("")
out.append("[CERT] direct b-sum B(z,M): band exact value z < 2 => B = 1 [EXACT];")
out.append("[CERT]   outside band, TRUE value via Mobius identity (see exemplars);")
out.append("[CERT]   gamma (Euler-Mascheroni)               = " + dec(gamma_e, 20))
KIV = arb(fmpq(0))
r = 2
while True:
    expo = (-arb(r) * arb(r).log()) * fmpq(1, 10)
    term = (arb(r + 1) * expo.exp()).upper()
    if term < arb(10) ** -30 and r > 100:
        break
    KIV += arb(r + 1) * expo.exp()
    r += 1
    if r > 5000:
        raise AssertionError("KIV series did not converge by r=5000")
out.append("[CERT] class-IV r-sum, Lemma-4 form (A5 = 1):  sum_{r>=2} (r+1) e^(-r ln r/10)")
out.append("[CERT]   < " + dec(KIV.upper(), 12) + "  (direct; r <= " + str(r - 1) + " + certified tail)")
out.append("[CERT] class-IV direct b-sum form (no Lemma 4): per-r <= (1/2) ln 2 = "
           + dec(ln2 / 2, 12))
out.append("")
out.append("=" * W)
out.append("THE PER-CLASS (45)-METER")
out.append("=" * W)
out.append("")
out.append("[EXACT] (2.2) k < y^(1-alpha)  =>  alpha_max(M) = ln(y/M)/ln y  (> 0 iff M < y)")
out.append("[EXACT]   Shiu's theorem covers exactly M in [sqrt(2N), y = 2 sqrt N);")
out.append("[EXACT]   above 2 sqrt N the theorem is silent (hypotheses fail).")
out.append("[EXACT] (5.1) z = y^(alpha/10);  ln z' = ln(z^(1/2)) = alpha ln(y)/20;")
out.append("[EXACT]   at alpha = alpha_max:  ln z' = ln(y/M)/20  [EXACT identity]")
out.append("[EXACT] block bound: survivors(window) <= y B(z,M)/(phi(M) ln z)")
out.append("        + z'^(3/2) M/y + window-edge floors")
out.append("[EXACT] amortized over N/y windows: junk additive, >= 0, only helps R;")
out.append("[EXACT]   no o(1) decay claim (would be false at finite N); at band")
out.append("[EXACT]   depths z'^(3/2) M/y <= (z')^2 = o(y) absorbed into additive junk")
out.append("[EXACT]   R(alpha) := (amortized bound)/(N/M) = B(z,M) (M/phi(M)) / ln z'")
out.append("            = 2 B(z,M) (M/phi(M)) / ln z")
out.append("[EXACT] R minimized at alpha = alpha_max; in-band B = 1 exactly, so")
out.append("[EXACT]   R_min = (M/phi(M)) / ln z'_max = (40/ln 2)(M/phi(M)) / log_2(4N/M^2)")
out.append("")
out.append("[CERT] frontier (M = sqrt(2N)): log_2(4N/M^2) = 1,  ln z' = (ln 2)/40 = "
           + dec(ln2 / 40, 12))
out.append("[CERT]   R_min >= (40/ln 2)(M/phi) = " + dec(40 / ln2, 12) + " (M/phi(M)) >= 57.708")
out.append("[CERT] universal (any depth): R >= (M/phi(M)) / ln z' with ln z' > 0;")
out.append("[CERT]   R -> infinity as M -> 2 sqrt N (alpha_max -> 0): vacuous at the top.")
out.append("[CERT] balanced-sieve contrast (junk z'^2 balanced against main term):")
out.append("[CERT]   ln z' <= (1/2) ln(y/M) + O(ln ln), so")
out.append("[CERT]   R >= 2 (M/phi) / ln(y/M) >= (4/ln 2)(M/phi) = " + dec(4 / ln2, 12)
           + " (M/phi) >= 5.7708")
out.append("[CERT] idealized unbounded depth (Mobius-asymptotic B above): R -> 2:")
out.append("[CERT]   the cascade NEVER proves R < 1 at any depth (its own b-sum")
out.append("[CERT]   density phi/M contributes the second factor 2).")
out.append("[CERT] e_req-analog: (45) needs R < 1  <=>  ln z' > M/phi(M); Lemma-2 shape")
out.append("[CERT]   supplies ln z' = ln(y/M)/20 (Shiu depth) or ~ (1/2) ln(y/M)")
out.append("[CERT]   (balanced); deficit: ln(y/M) < 40 (M/phi) resp. 4 (M/phi);")
out.append("[CERT]   frontier shortfall = (40/ln 2)(M/phi) >= 57.708 (M/phi(M)).")
out.append("[EXACT] above M = y = 2 sqrt N: (2.2) fails for all alpha > 0: no theorem.")
out.append("")
out.append("=" * W)
out.append("EXEMPLARS (exact phi and Mobius; R_min = (M/phi)/ln z'_max, in-band)")
out.append("=" * W)
out.append("")
exemplars = [
    (10 ** 12, 10 ** 6),
    (10 ** 12, 1414214),          # ~sqrt(2N): the band frontier
    (10 ** 12, 1999999),          # just under y = 2*10^6 (alpha_max -> 0)
    (10 ** 300, 10 ** 40),
    (10 ** 300, 10 ** 70),
]
for (N, M) in exemplars:
    ph = phi_exact(M)
    ratio = fmpq(M, ph)
    ln4N = arb(4 * N).log()
    ln_y = ln4N / 2
    alpha_max = 1 - arb(M).log() / ln_y
    ln_zmax = (ln4N - 2 * arb(M).log()) / 40   # ln z' = ln(y/M)/20  [EXACT]
    alt = (ln_y - arb(M).log()) / 20
    assert arb((ln_zmax - alt).upper()) < arb(10) ** -30
    if ln_zmax.lower() < ln2:                  # band: z' < 2, B = 1 exactly
        Bs = arb(1)
        B_note = "B = 1 exactly (band, z' < 2; only b = 1 contributes)"
        in_band = True
    else:                       # contrast: z' >= 2, TRUE b-sum via Mobius
        in_band = False
        zfl = ln_zmax.lower().exp().floor().unique_fmpz()
        divs = [d for d in range(1, 41) if M % d == 0]
        B_lo = arb(0)
        B_hi = arb(0)
        for d in divs:
            nd = zfl // d
            Hlo = (arb(nd + 1)).log()          # H_n >= ln(n+1)
            Hhi = (arb(nd)).log() + arb(1)     # H_n <= ln n + 1
            mu = moebius(d)
            if mu == 1:
                B_lo += Hlo / fmpq(d)
                B_hi += Hhi / fmpq(d)
            elif mu == -1:
                B_lo += -Hhi / fmpq(d)
                B_hi += -Hlo / fmpq(d)
        Bs = B_lo
        B_note = "B via Mobius identity (exact divisors d|M <= 40, H-brackets)"
    R_min = Bs * arb(ratio) / ln_zmax.lower() * (1 - arb(10) ** -30)
    tag = "IN BAND" if in_band else "BELOW BAND (survey only)"
    out.append("(N = 10^" + str(len(str(N)) - 1) + ", M = " + str(M) + "):  " + tag)
    out.append("  phi(M) = " + str(ph) + "   M/phi = " + str(ratio) + " = "
               + dec(arb(ratio), 12))
    out.append("  alpha_max = " + dec(alpha_max, 12) + "   ln z'_max = ln(y/M)/20 = "
               + dec(ln_zmax, 12))
    out.append("  " + B_note)
    out.append("  R_min >= " + dec(R_min, 9))
    out.append("")
out.append("[CERT] in-band rows: B = 1 (exact), lower wall of ln z', exact M/phi:")
out.append("[CERT]   honest lower bounds; junk would only increase R.")
out.append("[CERT] BELOW-BAND rows are survey-only (M < sqrt(2N) is outside the")
out.append("[CERT]   (45)-band: nothing there bears on (45)): B computed from the")
out.append("[CERT]   exact Mobius identity B(z,M) = sum_{d|M} mu(d)/d H_{floor(z/d)},")
out.append("[CERT]   asymptotics: B = (phi/M)(ln z + gamma + sum_{p|M} ln p/(p-1)) + o(1).")
out.append("")
out.append("=" * W)
out.append("RECONCILIATION WITH THE LEDGER FIGURE (section 3 row)")
out.append("=" * W)
out.append("")
out.append("[CERT] recorded normalization: e^{B1}/rho(2) = " + dec(e_B1 / rho2, 12)
           + "  (= 4.23x)")
out.append("[CERT]   B1 (Meissel-Mertens)               = " + dec(B1, 20))
out.append("[CERT]   rho(2) = 1 - ln 2                   = " + dec(rho2, 20))
out.append("[CERT] (a) the printed exp-form charges exp(B1)-type truncated Mertens;")
out.append("[CERT]     the proof recovers the (phi/M) b-sum density, so the")
out.append("[CERT]     proof-realized idealized-depth constant is 2 (M/phi(M)) --")
out.append("[CERT]     the theorem-form is (M/phi(M))-fold looser than its own proof.")
out.append("[CERT] (b) at every finite N in the band the realized margin is")
out.append("[CERT]     >= (40/ln 2)(M/phi(M)) = 57.708 (M/phi(M)).  Moreover (2.2)")
out.append("[CERT]     k < y^(1-alpha) forces M < y = 2 sqrt N at any alpha > 0: the")
out.append("[CERT]     (45) band M in [sqrt(2N), (2N)^(2/3)] loses its upper half")
out.append("[CERT]     to coverage failure (no theorem above 2 sqrt N) and its lower")
out.append("[CERT]     half to the level failure (57.708 (M/phi) depth collapse).")
out.append("")
out.append("=" * W)
out.append("KOROLEV 2019 CROSS-CHECK (assignment item 4)")
out.append("=" * W)
out.append("")
out.append("[NOTE] arXiv:1911.09981 (abstract checked 2026-08-29): estimates for")
out.append("[NOTE] Kloosterman sums over primes to COMPOSITE modulus q, non-trivial in")
out.append("[NOTE] the range q^(3/4+eps) <= X << q^(3/2) (the ledger section-2 row")
out.append("[NOTE] records eq. (7) with X >= q^(7/10+eps) -- a premise-level nuance")
out.append("[NOTE] flagged to Main, not re-audited here).")
out.append("[NOTE] VERDICT: NO / orthogonal.  Korolev's (X, q) parameters live on the")
out.append("[NOTE] Input-I side (inverse-prime sums; prime-side eviction).  His")
out.append("[NOTE] modulus premise carries no per-class sieve depth for")
out.append("[NOTE] smooth-indicator progressions; no bridge lets his saving structure")
out.append("[NOTE] dodge the Q = y depth collapse.  Serves the prime-side (Weyl-sum,")
out.append("[NOTE] input III) route only.")
out.append("")
out.append("=" * W)
out.append("THREE-WAY SPLIT SUMMARY (assignment item 1)")
out.append("=" * W)
out.append("")
out.append("(i) printed-literal in source: Lemma 2's constant 1 and z^2 term;")
out.append("    Lemma 1's 3 and 2 (incl. (4.1)'s 2y/log^2 y); Lemma 4's -(1/10) r ln r")
out.append("    with the 2 A1 r^(1/4) terms and r-border log z/(4 log log z); (5.1)")
out.append("    z = y^(alpha/10); (5.3)-(5.8) factors 2 (log z^(1/2) split), z^(-1/4),")
out.append("    z^(-1/8), f(n) <= n^(alpha beta/80), r0 = [log z/log(log x log log x)],")
out.append("    A5 = A1^(20/(alpha beta)); p.169 assembly chain z^2 < y/(phi(k) log z).")
out.append("(ii) derivable from stated inequalities, Arb-certified here:")
out.append("    A1 = A2 = 1 for f = 1_S; Lambda3 <= 1 (exact telescoping);")
out.append("    A3 <= 2(zeta(3/2)-1) = " + dec(A3, 9) + "; band b-sum B = 1 (exact);")
out.append("    class-IV r-sum both forms (" + dec(KIV.upper(), 9) + " resp. <= (1/2) ln 2);")
out.append("    alpha_max, ln z'_max = ln(y/M)/20; the R-meter with floors")
out.append("    57.708 (M/phi) [Shiu depth], 5.7708 (M/phi) [balanced], 2 [idealized].")
out.append("(iii) genuinely ineffective: NONE for f = 1_S.  Every << and O(1) in the")
out.append("    f = 1_S-relevant chain realizes explicitly; the only external link is")
out.append("    (4.1)'s PNT form, effective via published explicit bounds")
out.append("    (Rosser-Schoenfeld/Dusart; labeled HUMAN-AUDITED, not re-derived).")
out.append("")
out.append("OVERALL VERDICT (assignment item 3): NO -- the effective Shiu margin never")
out.append("enters (45) with constant < 1: R >= 57.708 (M/phi(M)) >= 57.708 at Shiu's")
out.append("own depth across the applicable band M in [sqrt(2N), 2 sqrt N); >= 5.7708")
out.append("(M/phi(M)) even at optimal Selberg balance; -> 2 (never < 1) at idealized")
out.append("depth; no theorem at all above M = 2 sqrt N.  The binding obstruction is the")
out.append("per-class Selberg sieve LEVEL (Lemma 2's 1/ln z' per-class shape with depth")
out.append("ln z' = ln(y/M)/20 collapsing across the band), not any ineffective constant")
out.append("and not the Meissel-Mertens normalization.")
out.append("=" * W)
sys.stdout.write("\n".join(out) + "\n")

# ---------------------------------------------------------------- assertions
assert (40 / ln2).lower() > arb(57), "frontier floor below 57"
assert (4 / ln2).lower() > arb(5), "balanced floor below 5"
if not exemplars:
    raise AssertionError("no exemplars")
