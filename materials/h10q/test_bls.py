"""BLS-relaxed Pocklington: pinned-regime unit check.

Builds n = F*s+1 with certified F-part (2, q1, q2) and a brent-hard
25-26 digit semiprime tail s = s1*s2, such that F^2 <= n < F^3; the old
code raises PrimalityBound, the new code must prove primality for prime n
and return False for composite n.  NOT part of the suite; manual use.
"""
import time
from h10q import _is_prime, _mr_screen, PrimalityBound


def next_prime(x):
    x |= 1
    while not _is_prime(x):
        x += 2
    return x


q1, q2 = next_prime(1234577), next_prime(7654321)
s1 = next_prime(10**12 + 39)
Fpart = 2*q1*q2

found = 0
s2 = next_prime(2*10**12 + 61)
t0 = time.time()
tested = 0
while found < 3 and time.time() - t0 < 600:
    if _mr_screen(s2):
        s = s1*s2
        n = Fpart*s + 1
    else:
        s2 += 2
        continue
    s2 += 2
    if not _mr_screen(n):
        continue
    s = s1*(s2 - 2)
    tested += 1
    assert (n - 1) == Fpart*s and Fpart**2 <= n <= Fpart**3
    t1 = time.time()
    try:
        v = _is_prime(n)
    except PrimalityBound:
        print(f"  PrimalityBound: {n}")
        continue
    print(f"  _is_prime({n}) = {v} ({time.time()-t1:.2f}s)")
    assert v is True
    # composite cousin: n*s2 must be refused-or-False but NEVER True
    try:
        vc = _is_prime(n*s2)
        assert vc is False, "FALSE PRIME PROOF ON COMPOSITE"
    except PrimalityBound:
        pass
    found += 1
print(f"BLS pinned-regime: {found} prime proofs with composite-cousin guard, "
      f"{time.time()-t0:.1f}s total")
assert found >= 1, "no pinned-regime vector found in time budget"
print("BLS unit check OK")
