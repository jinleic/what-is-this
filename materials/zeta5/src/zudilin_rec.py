"""Baseline for the zeta(5) campaign: exact verification of Zudilin's recursion.

Primary source (read first-hand 2026-08-15):
  W. Zudilin, "A third-order Apery-like recursion for zeta(5)",
  arXiv:math/0206178v2 (2002).

The PDF-to-text conversion detached several exponents; the reconstruction used
here is pinned four independent ways, all checked below:
  * all four coefficient polynomials of (1) must have degree 9, which forces
    (n+1)^6 on the x_{n+1} term and (n-1)^4 on the x_{n-2} term;
  * the leading coefficients must reproduce the printed characteristic
    polynomial  mu^3 + 2368 mu^2 - 752 mu - 16  (e.g. 2*48802112/41218 = 2368);
  * the n=1 instance of (1) reproduces the printed q_2 = -17934 exactly;
  * the generated convergents p_n/q_n equal the paper's printed fractions
    for every n <= 7.

Checks (exact unless labelled NUMERICAL):
  table        p_n/q_n equals the printed fraction, n = 0..7; the printed
               |zeta(5) - p_n/q_n| bounds hold at n = 3..7, 10, 20, 50.
  integrality  paper's (6):  q_n in Z,  2 D_n^5 p_n in Z,  2 D_n^3 pt_n in Z,
               1 <= n <= N, where D_n = lcm(1..n)  (empirical for each n; the
               paper itself marks (6) as observed, with only the weaker (14)
               proved).
  signs        paper's (3): l_n > 0, lt_n < 0; and (-1)^(n-1) q_n, p_n, pt_n > 0.
  roots        polyroots of mu^3+2368mu^2-752mu-16 match the printed decimals
               (5); mu_1 = lam_1*lam_2, mu_2 = lam_1*lam_3, mu_3 = lam_2*lam_3
               for the auxiliary polynomial lam^3 - 188 lam^2 - 2368 lam + 4.
  scoping      the n=1 instance of (1) holds for {q_n} but NOT for {p_n}: the
               recursion is imposed for n >= 2 only.
  rates        NUMERICAL: (1/n) log|l_n| -> log|mu_2| = -1.08607936...,
               (1/n) log|q_n| -> log|mu_3|, sampled at n = N.
  deficit      irrationality of zeta(5) via these forms would need
               lim (1/n) log|l_n| < -5 (with (6); -7 with only the proved
               (14)-type inclusions).  This construction achieves
               log|mu_2| ~ -1.086: deficit ~ 3.914.  Recorded, not asserted
               away -- the deficit IS the baseline result.

Writes ../data/BASELINE.json.  Usage:
  cd math && ./.venv/bin/python zeta5/src/zudilin_rec.py [N]   # default N=300
"""

from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

ARXIV = "arXiv:math/0206178v2"


# --- recursion (1), transcribed -------------------------------------------

def a0(n: int) -> int:
    return 41218 * n**3 - 48459 * n**2 + 20010 * n - 2871


def a1(n: int) -> int:
    return 2 * (48802112 * n**9 + 89030880 * n**8 + 36002654 * n**7
                - 24317344 * n**6 - 19538418 * n**5 + 1311365 * n**4
                + 3790503 * n**3 + 460056 * n**2 - 271701 * n - 60291)


def a2(n: int) -> int:
    return (3874492 * n**8 - 2617900 * n**7 - 3144314 * n**6
            + 2947148 * n**5 + 647130 * n**4 - 1182926 * n**3
            + 115771 * n**2 + 170716 * n - 44541)


def coeffs(n: int) -> tuple[int, int, int, int]:
    """(c3, c2, c1, c0) with  c3*x_{n+1} + c2*x_n + c1*x_{n-1} + c0*x_{n-2} = 0."""
    return ((n + 1)**6 * a0(n),
            a1(n),
            -4 * (2 * n - 1) * a2(n),
            -4 * (n - 1)**4 * (2 * n - 1) * (2 * n - 3) * a0(n + 1))


# initial data as printed (q_2 is redundant: the n=1 instance regenerates it)
Q_INIT = (Fraction(-1), Fraction(42), Fraction(-17934))
P_INIT = (Fraction(0), Fraction(87, 2), Fraction(-1190161, 64))
PT_INIT = (Fraction(0), Fraction(101, 2), Fraction(-344923, 16))


def generate(N: int) -> tuple[list[Fraction], list[Fraction], list[Fraction]]:
    """q_n, p_n, pt_n for 0 <= n <= N via recursion (1) imposed for n >= 2."""
    q, p, pt = list(Q_INIT), list(P_INIT), list(PT_INIT)
    for n in range(2, N):
        c3, c2, c1, c0 = coeffs(n)
        for s in (q, p, pt):
            s.append(Fraction(-(c2 * s[n] + c1 * s[n - 1] + c0 * s[n - 2]), c3))
    return q, p, pt


# --- paper's table (page 2), exact fractions and printed error bounds -----

TABLE_FRACTIONS = {
    0: Fraction(0),
    1: Fraction(29, 28),
    2: Fraction(24289, 23424),
    3: Fraction(7682021239, 7408444032),
    4: Fraction(24943788950905, 24055474286592),
    5: Fraction(81875586674776013003, 78959779279372800000),
    6: Fraction(282653756112686336975107, 272587704119854963200000),
    7: Fraction(215903781003833520407770175189,
                208214873150908926517286400000),
}
TABLE_BOUNDS = {3: "2.80e-11", 4: "4.13e-15", 5: "6.02e-19", 6: "8.71e-23",
                7: "1.26e-26", 10: "3.71e-38", 20: "1.32e-76", 50: "5.52e-192"}

# printed root decimals (5) and both characteristic polynomials; the paper
# prints truncations ("..."), so agreement is checked numerically to the
# number of printed fractional digits
MU_POLY = [1, 2368, -752, -16]          # mu^3 + 2368 mu^2 - 752 mu - 16
LAM_POLY = [1, -188, -2368, 4]          # lam^3 - 188 lam^2 - 2368 lam + 4
MU_PRINTED = ("-0.02001512", "0.33753726", "-2368.31752213")
LOG_MU2_PRINTED = "-1.08607936"


def matches_printed(value: mp.mpf, printed: str) -> bool:
    frac_digits = len(printed.split(".")[1])
    return abs(value - mp.mpf(printed)) < 2 * mp.mpf(10) ** (-frac_digits)


def main() -> int:
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    if N < 60:
        sys.exit("N >= 60 required (table checks reach n=50)")
    t0 = time.time()
    failures: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        tag = "PASS" if ok else "FAIL"
        print(f"[{tag}] {name}" + (f" -- {detail}" if detail else ""))
        if not ok:
            failures.append(name)

    # precision: |q_N| ~ 2368^N, |l_N| ~ 0.3375^N; margin on top of both
    mp.mp.dps = int(3.86 * N) + 60
    zeta5 = mp.zeta(5)
    zeta3 = mp.zeta(3)

    q, p, pt = generate(N)

    def mpf_of(x: Fraction) -> mp.mpf:
        return mp.mpf(x.numerator) / mp.mpf(x.denominator)

    ell = [mpf_of(q[n]) * zeta5 - mpf_of(p[n]) for n in range(N + 1)]
    ellt = [mpf_of(q[n]) * zeta3 - mpf_of(pt[n]) for n in range(N + 1)]

    # 1. table: exact convergents and printed error bounds
    check("table: p_n/q_n exact, n=0..7",
          all(p[n] / q[n] == TABLE_FRACTIONS[n] for n in TABLE_FRACTIONS))
    err1 = abs(zeta5 - mpf_of(p[1] / q[1]))
    err2 = abs(zeta5 - mpf_of(p[2] / q[2]))
    check("table: printed decimals at n=1,2",
          int(err1 * 10**9) == 1213469 and int(err2 * 10**9) == 182,
          f"|err_1|={mp.nstr(err1, 10)} |err_2|={mp.nstr(err2, 10)}")
    check("table: printed bounds at n=3..7,10,20,50",
          all(abs(zeta5 - mpf_of(p[n] / q[n])) < mp.mpf(b)
              for n, b in TABLE_BOUNDS.items()))

    # 2. integrality (6), each n up to N (empirical, as in the paper)
    D = 1
    ok6 = True
    for n in range(1, N + 1):
        D = math.lcm(D, n)
        if not (q[n].denominator == 1
                and (2 * D**5 * p[n]).denominator == 1
                and (2 * D**3 * pt[n]).denominator == 1):
            ok6 = False
            break
    check(f"integrality (6): q_n in Z, 2D^5 p_n in Z, 2D^3 pt_n in Z, n<={N}", ok6)

    # 3. signs (3) and alternation
    check("signs (3): l_n>0 and lt_n<0, 1<=n<=N",
          all(ell[n] > 0 and ellt[n] < 0 for n in range(1, N + 1)))
    check("alternation: (-1)^(n-1) * (q_n, p_n, pt_n) > 0, 1<=n<=N",
          all((-1)**(n - 1) * s[n] > 0
              for n in range(1, N + 1) for s in (q, p, pt)))

    # 4. characteristic roots: printed decimals and the mu = lam*lam pairing
    mu = sorted(mp.polyroots(MU_POLY), key=abs)
    lam = sorted(mp.polyroots(LAM_POLY), key=abs)
    mu = [mp.mpf(x.real) for x in mu]     # all real
    lam = [mp.mpf(x.real) for x in lam]
    check("roots: printed decimals (5)",
          all(matches_printed(m, s) for m, s in zip(mu, MU_PRINTED)),
          " ".join(mp.nstr(m, 12) for m in mu))
    pairing = (lam[0] * lam[1], lam[0] * lam[2], lam[1] * lam[2])
    check("roots: mu_i = (lam lam) pairing across (1) and (8)/(9)",
          all(abs(m - pr) < mp.mpf(10) ** (-(mp.mp.dps // 2))
              for m, pr in zip(mu, pairing)))

    # 5. scoping: n=1 instance holds for q, not for p
    c3, c2, c1, c0 = coeffs(1)
    assert c0 == 0  # (n-1)^4 kills the x_{-1} term
    check("scoping: n=1 instance of (1) holds for {q_n}",
          c3 * q[2] + c2 * q[1] + c1 * q[0] == 0)
    check("scoping: n=1 instance of (1) FAILS for {p_n} (so (1) is n>=2 only)",
          c3 * p[2] + c2 * p[1] + c1 * p[0] != 0)

    # 6. rates (4)-(5)  [NUMERICAL at n=N] and the deficit.
    # Estimator: log|x_N / x_{N-1}|.  With x_n ~ C mu^n n^theta (Poincare /
    # Birkhoff-Trjitzinsky) this is log|mu| + theta/N + o(1/N), so the
    # tolerance 6/N covers any |theta| <= 6; the naive (1/N)log|x_N| carries
    # a Theta(log N / N) correction and is reported as INFO only.
    log_mu2 = mp.log(abs(mu[1]))
    log_mu3 = mp.log(abs(mu[2]))
    tol = mp.mpf(6) / N
    r_ell = mp.log(abs(ell[N] / ell[N - 1]))
    r_ellt = mp.log(abs(ellt[N] / ellt[N - 1]))
    r_q = mp.log(abs(mpf_of(q[N] / q[N - 1])))
    check("rates (4): log|l_N/l_{N-1}| ~ log|mu_2|  [NUMERICAL]",
          abs(r_ell - log_mu2) < tol and abs(r_ellt - log_mu2) < tol,
          f"l: {mp.nstr(r_ell, 10)}  lt: {mp.nstr(r_ellt, 10)} "
          f" target {mp.nstr(log_mu2, 10)}  tol {mp.nstr(tol, 3)}")
    check("rates (5): log|q_N/q_{N-1}| ~ log|mu_3|  [NUMERICAL]",
          abs(r_q - log_mu3) < tol,
          f"q: {mp.nstr(r_q, 10)}  target {mp.nstr(log_mu3, 10)}")
    check("rates: printed log|mu_2| decimal",
          matches_printed(log_mu2, LOG_MU2_PRINTED))
    print(f"[INFO] cumulative (1/N)log|l_N| = "
          f"{mp.nstr(mp.log(abs(ell[N])) / N, 10)}, (1/N)log|q_N| = "
          f"{mp.nstr(mp.log(abs(mpf_of(q[N]))) / N, 10)} "
          f"(carry Theta(log N/N) corrections)")

    # deficit: with (6), need lim (1/n)log|l_n| < -5; have log|mu_2|
    deficit6 = 5 + log_mu2
    deficit14 = 7 + log_mu2
    needed_mu2 = mp.exp(-5)
    print(f"[INFO] deficit vs (6):  5 + log|mu_2| = {mp.nstr(deficit6, 10)}"
          f"  (need < 0; equivalently |mu_2| < e^-5 = {mp.nstr(needed_mu2, 6)},"
          f" have {mp.nstr(abs(mu[1]), 9)})")
    print(f"[INFO] deficit vs proved (14)-type: 7 + log|mu_2| ="
          f" {mp.nstr(deficit14, 10)}")
    print(f"[INFO] Apery zeta(3) margin, for scale: needed < -3,"
          f" achieved -4log(1+sqrt2) = {mp.nstr(-4 * mp.log(1 + mp.sqrt(2)), 9)}")

    wall = time.time() - t0
    record = {
        "date": "2026-08-15",
        "source": ARXIV,
        "N": N,
        "dps": mp.mp.dps,
        "mu": [mp.nstr(m, 20) for m in mu],
        "log_mu2": mp.nstr(log_mu2, 20),
        "log_mu3": mp.nstr(log_mu3, 20),
        "rate_ratio_ell_at_N": mp.nstr(r_ell, 20),
        "rate_ratio_ellt_at_N": mp.nstr(r_ellt, 20),
        "rate_ratio_q_at_N": mp.nstr(r_q, 20),
        "rate_cumulative_ell_at_N": mp.nstr(mp.log(abs(ell[N])) / N, 20),
        "rate_cumulative_q_at_N": mp.nstr(mp.log(abs(mpf_of(q[N]))) / N, 20),
        "deficit_vs_inclusions_6": mp.nstr(deficit6, 20),
        "deficit_vs_inclusions_14": mp.nstr(deficit14, 20),
        "needed_abs_mu2": mp.nstr(needed_mu2, 20),
        "integrality_6_verified_to": N,
        "failures": failures,
        "wall_seconds": round(wall, 2),
    }
    out = Path(__file__).resolve().parent.parent / "data"
    out.mkdir(exist_ok=True)
    (out / "BASELINE.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"[INFO] wrote {out / 'BASELINE.json'}  ({wall:.1f}s)")

    if failures:
        print(f"FAILURES: {failures}")
        return 1
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
