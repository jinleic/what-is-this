"""Cross-check the two disputed literature rows against the PUBLISHED formula.

Panteleev-Kalachev Prop. 1 / Wang-Mueller Eq. (11):  for a coprime lattice the
BB code is a generalised bicycle code over the cyclic ring GF(2)[pi]/(pi^N - 1),
N = ell*m, and

    k = 2 * deg gcd( a(pi), b(pi), pi^N - 1 ).

Computed here with plain GF(2) polynomial arithmetic on bitmask integers, which
shares no code with the ring/annihilator route used by EXP-055.  If the two
routes disagree, the transcription of the paper row is wrong; if they agree with
each other but not with the printed k, the printed row is the outlier.
"""


def deg(p: int) -> int:
    return p.bit_length() - 1


def polymod(a: int, b: int) -> int:
    db = deg(b)
    while a and deg(a) >= db:
        a ^= b << (deg(a) - db)
    return a


def polygcd(a: int, b: int) -> int:
    while b:
        a, b = b, polymod(a, b)
    return a


def from_exponents(exps) -> int:
    v = 0
    for e in exps:
        v ^= 1 << e
    return v


CASES = [
    # (ell, m, pi_A, pi_B, k_published, d_published, source)
    (5, 9, [0, 1, 4], [0, 8, 34], 4, 12, "2408.10001v4 App.C Table 4"),
    (7, 11, [0, 4, 31], [0, 19, 53], 6, 16, "2408.10001v4 App.C Table 4"),
    # controls: rows where our ring route already agreed with the paper
    (3, 5, [0, 1, 2], [1, 3, 8], 4, 6, "2408.10001v4 Table 2 (control)"),
    (3, 7, [0, 2, 3], [1, 3, 11], 6, 6, "2408.10001v4 Table 2 (control)"),
    (5, 7, [0, 1, 5], [0, 1, 12], 6, 8, "2408.10001v4 Table 2 (control)"),
    (7, 9, [0, 1, 58], [3, 16, 44], 12, 10, "2408.10001v4 Table 2 (control)"),
]

for ell, m, ea, eb, k_pub, d_pub, src in CASES:
    N = ell * m
    mod = (1 << N) | 1                      # pi^N - 1 = pi^N + 1 over GF(2)
    a = from_exponents([e % N for e in ea])
    b = from_exponents([e % N for e in eb])
    g = polygcd(polygcd(a, b), mod)
    k_formula = 2 * deg(g)
    print(f"({ell},{m}) N={N:3d} k_published={k_pub:2d} k_gcd_formula={k_formula:2d} "
          f"{'AGREE' if k_formula == k_pub else 'DISAGREE'}   {src}")
