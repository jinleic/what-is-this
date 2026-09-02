#!/usr/bin/env python3
"""Anchor verification: all five decompositions must pass 729/729 Brent over Z
BEFORE the sweep starts. Any miss aborts with exit 1."""
import sys
from flint import fmpz

sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/src")
sys.path.insert(0, "/Users/jinleic/jinleic-workspace/cs/mm3/scratch")

from gatec_decomps import LOADERS, META, R


def brent_failures(U, V, W):
    fails = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for ip in range(3):
                    for jp in range(3):
                        for kp in range(3):
                            s = fmpz(0)
                            for r in range(R):
                                ur, vr, wr = U[r][3 * i + k], V[r][3 * kp + j], W[r][3 * ip + jp]
                                if ur and vr and wr:
                                    s += ur * vr * wr
                            expect = 1 if (i == ip and j == jp and k == kp) else 0
                            if int(s) != expect:
                                fails += 1
    return fails


def main():
    all_ok = True
    for name, loader in LOADERS.items():
        U, V, W = loader()
        tern = all(x in (-1, 0, 1) for blk in (U, V, W) for row in blk for x in row)
        f = brent_failures(U, V, W)
        ok = (f == 0) and tern
        all_ok &= ok
        print(f"{name:13s} brent={729 - f}/729 ternary={tern} anchor={META[name]['anchor']:3d} "
              f"{'OK' if ok else 'FAIL'}")
    if not all_ok:
        sys.exit(1)
    print("ALL ANCHORS PASS — sweep counter may proceed")


if __name__ == "__main__":
    main()
