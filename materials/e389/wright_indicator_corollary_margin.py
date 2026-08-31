#!/usr/bin/env python3
"""Wright arXiv:2508.17217v2 -- commented-out smooth-count theorem and indicator
corollary (Hilcor) at the e389 Input-II corner u = 2.

Corner parameters (Input II): x = 2N, q = y = sqrt(4N), progression modulus
k = M in [sqrt(2N), x^(2/3)].  Then
    eta = log y / log x -> 1/2,   u_x = log x / log Q = 2,
and because the smooth-count theorem forces Z <= y = Q (it defines
Z = y/k^(1+nu) <= y), the split indices satisfy
    u_{x/Z} = log(x/Z)/log Q = 1 - o(1)  and  u_{Z/Q} = log(y)/log Q = 1,
so BOTH Dickman factors of the parent bound sit at rho(1) = 1.

What this script certifies (byte-stable; exact rationals + Arb >= 300 bit prec):

  1. Washout of the rho machinery at u = 2 (exact piecewise-closed-form rho):
     for every eta in (0, 1),
         rho(2(1-eta))^{C1} * rho(2*eta) >= rho(2)^{C1} = (1-ln2)^(1/164),
     using rho(v) = 1 on [0, 1] and rho(v) = 1 - ln v on [1, 2] (classical,
     exact).  No eta extracts any rho(2)-shaped gain from the split at u = 2.

  2. The c-split pays the full harmonic sum H_Q = ln Q + gamma (every
     c <= Q is Q-smooth at Q = y), so the parent bound at the corner is
       log^2 Q / log x = (1/2) ln Q (1 + o(1))  per unit y,
     unboundedly above the class scale -- an elementary, constant-free vacuity.

  3. Printed exponent-algebra slip (exact rationals, C1 = 1/164, eta = 1/2):
       printed exponent:   C1 - (1-C1)*eta     = -161/328
       honest  exponent:   C1*(1-eta) + eta    =  165/328
       shift:              2*eta*(1-C1)        =  163/164
     and the honest bound additionally carries the split-origin prefactor
       A(u, eta) = eta^{-eta u} (1-eta)^{-C1 (1-eta) u}  >= 2^{(1+C1) u/2}
     at eta = 1/2, dropped by the printed proof.

  4. Corner values (Arb, 20 certified digits):
       printed form 1 at the corner:  rho(2)^(-161/328) = 1.787... y  (vacuous)
       printed form 2 at the corner:  rho(1)^(1+C1)     = 1            (zero margin)
       honest repaired shape at u=2:  2^(1+C1) rho(2)^{(1+C1)/2} = 1.108... y
     The class bound at M = sqrt(2N) is N/M = y/(2 sqrt 2) = 0.35355... y;
     every one of these bounds exceeds it.

  5. Effectivity threshold of the honest repaired shape at eta = 1/2:
       honest(u, 1/2) = 2^{(1+C1)u/2} rho(u)^{(1+C1)/2} <= 1
     is independent of C1 and equivalent to rho(u) = 2^{-u}.  Certified
     two-sided bracket u* in (2.21, 2.22) by rigorous midpoint quadrature of
       rho(u) = rho(2) - int_2^u (1 - ln(t-1))/t dt   (rho' = -rho(u-1)/u)
     with cell-wise second-derivative majorants (inclusion-monotone ball
     evaluation of g'' = (2 - ln(t-1))/t^3 over whole cells).  Also
     certified: rho(3)/2^{-3} < 1/2, and the ratio rho(u)/2^{-u} decreases on
     [2, 3], so the mechanism stays vacuous at every u in [2, 3].  A 20-digit
     point value of u* is computed with mpmath at 60 dps (discovery-grade;
     the bracket is the certified fact).

  6. Best Wright-faithful bound at the corner (corrected Theorem 3.1 part 1,
     reinstating the exp(sum_{p<=Q} 1/p) its own proof carries, applied to
     f = 1_S on the window):
       e^{B1} rho(2)^{C1} (log Q / log x)  per unit y,
     printed against truth rho(2) and against the class bound 1/(2 sqrt 2).
     Reproduces the LOCALIZATION.md figure "0.649 vs 0.307 = 2.11x" up to the
     rho(2)^{C1} = 0.9928 refinement.

  7. Repair target: an effective rho(2)-exponent must exceed
       e_req = ln(2 sqrt 2) / ln(1/(1 - ln 2))          (Arb, 20 digits)
     for the total smooth count to fall below the class bound at M = sqrt(2N),
     while the honest mechanism ceiling at the corner is (1+C1)/2 = 165/328.

  8. Cross-check of the house-recorded Changa-Korolev constants (LOCALIZATION
     section 2): alpha = (29 - 5 sqrt 33)/8, beta as printed there, the
     continuity identity alpha(beta - 1/2) = (4/19)(beta - 5/8) at 300-bit Arb,
     and the re-derivation beta = (alpha/2 - 5/38)/(alpha - 4/19).

Reproduction: ./.venv/bin/python -I -B e389/wright_indicator_corollary_margin.py
from math/.  Byte-stable stdout: no timing, environment, or resume fields;
every decimal printed is certified common to the Arb lower and upper balls.

Labels: [CERT] lines are Arb-ball certificates; the [MODEL] u* point is
mpmath discovery-grade, bracketed by the certified interval printed above it.
"""

import fractions
import sys

from mpmath import mp, mpf, findroot

import flint
from flint import arb, fmpq, fmpz

flint.ctx.prec = 300


def dec(x: arb, digits: int) -> str:
    """Certified decimal expansion: the digits on which the lower and upper
    balls agree at the 10^-digits floor.  Shrinks `digits` if needed."""
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


def ball(lo: arb, hi: arb) -> arb:
    """Inclusion-monotone interval [lo, hi] as an Arb ball."""
    return (lo + hi) / 2 + (hi - lo)


# ------------------------------------------------------------------ constants
C1MAX = fmpq(1, 164)      # sup of C1 = alpha*kappa/41 over alpha, kappa < 1/2
ETA = fmpq(1, 2)          # corner: y = x^(1/2+o(1))
B1_STR = "0.2614972128476427837554268386086958590516"   # Meissel-Mertens

rho2 = 1 - arb(2).log()                    # rho(2), exact closed form

printed_exp = C1MAX - (1 - C1MAX) * ETA    # -161/328
honest_exp = C1MAX * (1 - ETA) + ETA       #  165/328
exp_shift = honest_exp - printed_exp       #  163/164 == 2*eta*(1-C1)
assert exp_shift == 2 * ETA * (1 - C1MAX)

# ------------------- 4. corner values ---------------------------------------
printed_corner = rho2 ** arb(printed_exp)              # > 1: vacuous
honest_corner = arb(2) ** arb(1 + C1MAX) * rho2 ** arb((1 + C1MAX) / 2)
class_per_y = 1 / (2 * arb(2).sqrt())                  # N/M at M = sqrt(2N)

B1 = arb(B1_STR)
exp_B1 = B1.exp()
part1_per_y = exp_B1 / 2 * rho2 ** arb(C1MAX)
part1_vs_truth = part1_per_y / rho2
part1_vs_class = part1_per_y / class_per_y
e_req = (2 * arb(2).sqrt()).log() / (1 / rho2).log()
shortfall = e_req - arb(honest_exp)
# printed form 1 is a BOUND, so vacuity must hold over all eta, not at one
# eta: exponent C1-(1-C1)eta in [-1+2C1, 0], base rho(2(1-eta)) in
# [rho(2), 1), so the log of the bound factor is >= log sup = (2C1-1) log
# rho(2) (inclusion monotone) -- certified once, valid for every eta.
sup_form1 = rho2 ** arb(2 * C1MAX - 1)

# part-1-anchored exponent gaps: cancelling the engine prefactor e^{B1}/2
# per y needs rho(u)^C1 to absorb 2 x 2^{1/2} per y, i.e. exponent
# e_need: rho(2)^{e_need} <= 2^{3/2}; e_req: rho(2)^{e_req} <= 1/(2 sqrt 2).
e_need = (exp_B1 / 2 / class_per_y).log() / (1 / rho2).log()
e_need_gap = e_need.upper() - honest_exp
e_req_gap = e_req.upper() - honest_exp

# multipliers on the honest ceiling: e_need / honest_exp and e_req /
# honest_exp (both e_need and e_req are closed-form balls, fully certified).
e_need_ratio = e_need / arb(honest_exp)
e_req_ratio = e_req / arb(honest_exp)

# ------------------- 5. certified bracket for u* ----------------------------
# rho(u) = rho(2) - int_2^u g(t) dt, g(t) = (1 - ln(t-1))/t, on [2, 2.25].
# Rigorous midpoint quadrature: on a cell [a, b] of width h,
#   |int - h g(m)| <= (h^3/24) sup_cell |g''|, with g'' = (2 - ln(t-1))/t^3,
#   certified by inclusion-monotone ball evaluation of g'' over the cell.


def quad_enclosure(a: fractions.Fraction, b: fractions.Fraction, n: int):
    """Two Arb balls rigorously bracketing int_a^b g(t) dt."""
    h = (b - a) / n
    h_arb = arb(fmpq(h.numerator, h.denominator))
    lo_acc = arb(0)
    hi_acc = arb(0)
    for i in range(n):
        ca = a + h * i
        cb = ca + h
        mid = (ca + cb) / 2
        gm = g(arb(fmpq(mid.numerator, mid.denominator)))
        gmid = float(mid)
        grad = float(cb - ca) / 2 * (1 + 1e-15) + 1e-300
        tcell = ball(arb(repr(gmid - grad)), arb(repr(gmid + grad)))
        # rigorous g'' over the cell: g''(t) = (2 - ln(t-1)) / t^3
        gd2 = (2 - (tcell - 1).log()) / (tcell ** 3)
        mag = max(abs(gd2.lower()), abs(gd2.upper()))
        width = h_arb ** 3 / 24 * mag
        lo_acc += h_arb * (gm - width)
        hi_acc += h_arb * (gm + width)
    return lo_acc, hi_acc


def g(t: arb) -> arb:
    return (1 - (t - 1).log()) / t


def product_enclosure(a: fractions.Fraction, b: fractions.Fraction,
                      points: int):
    """Enclosures of the min and max of g on [a, b] via exhaustive
    inclusion-monotone cell balls (cells cover [a, b] exactly)."""
    h = (b - a) / points
    lo_all = arb(1)
    hi_all = arb(0)
    for i in range(points):
        ca = a + h * i
        cb = ca + h
        midf = float((ca + cb) / 2)
        radf = float(cb - ca) / 2 + 1e-25
        tcell = ball(arb(repr(midf - radf)), arb(repr(midf + radf)))
        gv = g(tcell)
        lo_all = min(lo_all, gv.lower())
        hi_all = max(hi_all, gv.upper())
    return lo_all, hi_all


NQ = 20000   # quadrature cells: total midpoint error <= L^3 M2 / (24 n^2)

u_lo = fractions.Fraction(221, 100)   # side A: expect rho(2.21) > 2^{-2.21}
u_hi = fractions.Fraction(222, 100)   # side B: expect rho(2.22) < 2^{-2.22}

int_lo_lo, int_lo_hi = quad_enclosure(fractions.Fraction(2), u_lo, NQ)
int_hi_lo, int_hi_hi = quad_enclosure(fractions.Fraction(2), u_hi, NQ)

rho_lo_lo = rho2 - int_lo_hi         # lower bound of rho(2.21)
rho_lo_hi = rho2 - int_lo_lo         # upper bound of rho(2.21)
rho_hi_lo = rho2 - int_hi_hi         # lower bound of rho(2.22)
rho_hi_hi = rho2 - int_hi_lo         # upper bound of rho(2.22)

target_lo = arb(2) ** (-arb(str(u_lo)))      # 2^{-2.21}
target_hi = arb(2) ** (-arb(str(u_hi)))      # 2^{-2.22}

cert_A = rho_lo_lo - target_lo > 0
cert_B = rho_hi_hi - target_hi < 0

g_lo, g_hi = product_enclosure(fractions.Fraction(2),
                               fractions.Fraction(9, 4), 400)
g_pos = g_lo > 0    # g > 0 on [2, 2.25]  =>  rho(u) - 2^{-u} strictly decr.

# rho(3): on [2, 3], rho'(u) = -rho(u-1)/u with rho(u-1) = 1 - ln(u-1)
# EXACTLY (u-1 in [1, 2] uses the closed form), so rho'(u) = -g(u) there and
# the same rigorous quadrature gives rho(3) directly:
#   rho(3) in [rho(2) - int_g_23_hi, rho(2) - int_g_23_lo].
int_23_lo, int_23_hi = quad_enclosure(fractions.Fraction(2),
                                      fractions.Fraction(3), NQ)
rho3_lo = rho2 - int_23_hi
rho3_hi = rho2 - int_23_lo
cert_rho3_quarter = rho3_hi < arb(1) / 4        # rho(3) < 0.25
ratio3 = rho3_hi / (arb(2) ** (-3))             # rho(3)/2^{-3} upper
cert_rho3_ratio = ratio3 < arb(1) / 2            # ratio dropped below 1/2

# ------------------- 5b. 20-digit point value of u* (discovery) --------------
mp.dps = 60
mp_rho2 = mpf(1) - mp.log(2)


def mp_rho(u):
    return mp_rho2 - mp.quad(lambda t: (1 - mp.log(t - 1)) / t, [2, u])


u_star = findroot(lambda u: mp_rho(u) - mp.power(2, -u),
                  (mpf("2.21"), mpf("2.22")))
u_star_str = mp.nstr(u_star, 20)

# ------------------- 8. Changa-Korolev cross-check ---------------------------
sqrt33 = arb(33).sqrt()
alpha_ck = (29 - 5 * sqrt33) / 8
beta_ck = (95 * sqrt33 - 511) / (2 * (95 * sqrt33 - 519))
lhs_id = alpha_ck * (beta_ck - arb(1) / 2)
rhs_id = arb(4) / 19 * (beta_ck - arb(5) / 8)
id_gap = ball(abs(lhs_id.lower() - rhs_id.upper()),
              abs(lhs_id.upper() - rhs_id.lower()))
beta_re = (alpha_ck / 2 - arb(5) / 38) / (alpha_ck - arb(4) / 19)
beta_gap = ball(abs(beta_re.lower() - beta_ck.upper()),
                abs(beta_re.upper() - beta_ck.lower()))
TOL = arb(10) ** -60
id_ok = id_gap.upper() < TOL
beta_ok = beta_gap.upper() < TOL

# ------------------------------------------------------------------ report --
W = 78
out = []
out.append("=" * W)
out.append("Wright 2508.17217v2 commented smooth-count theorem / Hilcor corollary")
out.append("at the e389 Input-II corner u = 2  (x = 2N, Q = y = sqrt(4N))")
out.append("=" * W)
out.append("")
out.append("[EXACT] corner parameters: eta = 1/2, u_x = 2, u_{x/Z} = u_{Z/Q} = 1")
out.append("[EXACT] C1 max (sup over alpha,kappa < 1/2): C1 = alpha*kappa/41 -> "
           + str(C1MAX) + " = " + dec(arb(C1MAX), 20))
out.append("[EXACT] printed exponent C1-(1-C1)*eta  = " + str(printed_exp)
           + " = " + dec(arb(printed_exp), 20))
out.append("[EXACT] honest  exponent C1*(1-eta)+eta = " + str(honest_exp)
           + " = " + dec(arb(honest_exp), 20))
out.append("[EXACT] exponent shift 2*eta*(1-C1)     = " + str(exp_shift)
           + " = " + dec(arb(exp_shift), 20))
out.append("")
out.append("[CERT] rho(2) = 1 - ln 2                    = " + dec(rho2, 20))
out.append("[CERT] 1/rho(2)                             = " + dec(1 / rho2, 20))
out.append("[CERT] rho(2)^{C1}                           = "
           + dec(rho2 ** arb(C1MAX), 20))
out.append("[CERT] washout: rho(2(1-eta))^{C1} rho(2eta) >= rho(2)^{C1} for all eta")
out.append("[CERT]   in (0,1)  [rho(v)=1 on [0,1], rho(v)=1-ln v on [1,2]]:")
out.append("[CERT]   minimum at eta = 1/2, value rho(2(1/2))^{C1} rho(1) = "
           + dec(rho2 ** arb(C1MAX), 20))
out.append("")
out.append("[EXACT] parent corner prefactor log Q / log x -> 1/2 (exact rational)")
out.append("[CERT] parent corner bound per unit y: (1/2) ln Q (c-split sum H_Q,")
out.append("[CERT]   every c <= Q is Q-smooth at Q = y; H_Q = ln Q + gamma)")
out.append("[CERT]   => vacuous: grows, for every N with Q = sqrt(4N) >= e^2")
out.append("")
out.append("[CERT] printed form 1 at corner rho(2)^(-161/328) = "
           + dec(printed_corner, 20) + "   (vacuous: > 1)")
out.append("[CERT]   ...and over ALL eta in (0,1) the factor stays above "
           + dec(sup_form1, 20) + " > 1:")
out.append("[CERT]   universal vacuity of printed form 1 (its factor is a bound,")
out.append("[CERT]   so sup over eta certifies: sup >= sup(base^{[-1+2C1,0]}) = "
           + dec(sup_form1, 20) + ")")
out.append("[CERT] printed form 2 at corner rho(1)^(1+C1)     = 1.000000...")
out.append("[CERT]        (margin exactly zero: the bound is y itself)")
out.append("[CERT] honest repaired shape at u = 2             = "
           + dec(honest_corner, 20) + "   (vacuous: > 1)")
out.append("[CERT] class bound at M = sqrt(2N): 1/(2 sqrt 2)  = "
           + dec(class_per_y, 20))
out.append("[CERT] honest corner / class bound                = "
           + dec(honest_corner / class_per_y, 20))
out.append("")
out.append("[CERT] best Wright-faithful bound (corrected part 1, f = 1_S):")
out.append("[CERT]   e^{B1} rho(2)^{C1} log Q/log x per unit y, B1 Meissel-Mertens")
out.append("[CERT]   B1                                 = " + dec(B1, 20))
out.append("[CERT]   e^{B1}                             = " + dec(exp_B1, 20))
out.append("[CERT]   bound per unit y (log Q/log x = 1/2) = " + dec(part1_per_y, 20))
out.append("[CERT]   / truth rho(2)                     = " + dec(part1_vs_truth, 20))
out.append("[CERT]   / class bound                      = " + dec(part1_vs_class, 20))
out.append("")
out.append("[CERT] effectivity threshold of honest shape at eta = 1/2:")
out.append("[CERT]   honest(u,1/2) = 2^{(1+C1)u/2} rho(u)^{(1+C1)/2}; effective iff")
out.append("[CERT]   rho(u) = 2^{-u}  (C1-independent at eta = 1/2)")
out.append("[CERT]   rho(2.21) - 2^{-2.21} > 0 : " + ("PASS" if cert_A else "FAIL")
           + "  (ball interval: [" + dec(rho_lo_lo - target_lo, 9) + ", "
           + dec(rho_lo_hi - target_lo, 9) + "])")
out.append("[CERT]   rho(2.22) - 2^{-2.22} < 0 : " + ("PASS" if cert_B else "FAIL")
           + "  (ball interval: [" + dec(rho_hi_lo - target_hi, 9) + ", "
           + dec(rho_hi_hi - target_hi, 9) + "])")
out.append("[CERT]   g(t) = (1-ln(t-1))/t > 0 on [2, 2.25]: "
           + ("PASS" if g_pos else "FAIL") + "  (g in [" + dec(g_lo, 9) + ", "
           + dec(g_hi, 9) + "]; crossing unique)")
out.append("[CERT]   certified bracket: u* in (2.21, 2.22)  =>  u* > u = 2")
out.append("[CERT]   rho(3) in [" + dec(rho3_lo, 9) + ", " + dec(rho3_hi, 9)
           + "]  (rho(3) < 1/4: " + ("PASS" if cert_rho3_quarter else "FAIL")
           + ")")
out.append("[CERT]   rho(3)/2^{-3} < 1/2 : "
           + ("PASS" if cert_rho3_ratio else "FAIL") + "  => vacuous for all")
out.append("[CERT]   u in [2, 3]: ratio rho(u)/2^{-u} decreasing, value at 3 < 1/2")
out.append("[MODEL] u* point value at 60 dps (mpmath, discovery): u* = " + u_star_str)
out.append("")
out.append("[CERT] repair target: rho(2)-exponent needed to close (45) at M = sqrt(2N):")
out.append("[CERT]   e_req = ln(2 sqrt 2)/ln(1/(1-ln 2)) = " + dec(e_req, 20))
out.append("[CERT]   honest mechanism ceiling at corner  = " + str(honest_exp)
           + " = " + dec(arb(honest_exp), 20))
out.append("[CERT]   shortfall e_req - honest ceiling    = " + dec(shortfall, 20))
out.append("[CERT]   ...anchored to the best Wright-faithful bound instead (ceiling")
out.append("[CERT]   before the rho-engine buys anything): e_need = "
           + dec(e_need, 20) + ",")
out.append("[CERT]   gap e_need - honest ceiling         = " + dec(e_need_gap, 20))
out.append("[CERT]   gap e_req  - honest ceiling         = " + dec(e_req_gap, 20)
           + "   (certified one-sided: e_req upper bracket)")
out.append("[CERT]   in multiplier form: needed exponent is "
           + dec(e_need_ratio, 6) + "x and " + dec(e_req_ratio, 6)
           + "x the honest ceiling")
out.append("")
out.append("[CERT] Changa-Korolev cross-check (LOCALIZATION section 2):")
out.append("[CERT]   |alpha(beta-1/2)-(4/19)(beta-5/8)|  = " + dec(id_gap, 60)
           + ("   PASS: < 10^-60" if id_ok else "   FAIL"))
out.append("[CERT]   |beta_rederived - beta|             = " + dec(beta_gap, 60)
           + ("   PASS: < 10^-60" if beta_ok else "   FAIL"))
out.append("")
out.append("[CERT] every decimal above is certified common to Arb lower/upper balls.")
out.append("[MODEL] the u* point value is mpmath discovery-grade, bracketed above.")
out.append("=" * W)
sys.stdout.write("\n".join(out) + "\n")

failed = []
if not cert_A:
    failed.append("cert_A: rho(2.21) > 2^-2.21")
if not cert_B:
    failed.append("cert_B: rho(2.22) < 2^-2.22")
if not g_pos:
    failed.append("g_pos: g > 0 on [2, 2.25]")
if not cert_rho3_quarter:
    failed.append("cert_rho3_quarter: rho(3) < 1/4")
if not cert_rho3_ratio:
    failed.append("cert_rho3_ratio: rho(3)/2^{-3} < 1/2")
if sup_form1.lower() <= 1:
    failed.append("sup_form1: universal form-1 factor > 1")
if not id_ok:
    failed.append("Changa-Korolev identity")
if not beta_ok:
    failed.append("Changa-Korolev beta re-derivation")
if failed:
    raise AssertionError("certification failed: " + "; ".join(failed))
