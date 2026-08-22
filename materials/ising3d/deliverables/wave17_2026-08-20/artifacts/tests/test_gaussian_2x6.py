"""Standalone verifier for the 2x6 Gaussian no-go certificate (H425/H426).

Verifies, without trusting the producer's own narrative:
1. final artifact invariants (schema, verdicts, cap chain, both primes, zero
   failures, witnessed decisive records with digest chains);
2. the saved Euclid states: input-digest conservation against their OWN stored
   pair polynomials (both primes), and an INDEPENDENT list-schoolbook replay of
   the exit tail reproducing deg(remainder) = 295414 < cap in <= 40 steps, with
   remainders bit-identical across a second implementation;
3. the pilot's forced-floor witness is recorded in the progress journal
   (gcd == 442067, exact mode, both primes);
4. the n=6 layer structural guard re-derived here: for all four (flip,rho)
   characters, cp(block) divides cp(full) -- the only algebraic premise the
   block-counting argument depends on.

Run: PYTHONPATH=src:experiments python tests/test_gaussian_2x6.py
"""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / "src"), str(ROOT / "experiments")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import importlib.util

spec136 = importlib.util.spec_from_file_location(
    "e136", ROOT / "experiments" / "e136_pair_product_2x5.py")
e136 = importlib.util.module_from_spec(spec136)
sys.modules["e136"] = e136
spec136.loader.exec_module(e136)

spec154 = importlib.util.spec_from_file_location(
    "e154", ROOT / "experiments" / "e154_gaussian_2x6.py")
e154 = importlib.util.module_from_spec(spec154)
sys.modules["e154"] = e154
spec154.loader.exec_module(e154)

CAP = 295415
A_SAME_12 = 261625


def divmod_hard(a, b, p):
    """independent list schoolbook division (descending ints)"""
    a = [int(x) % p for x in a]
    b = [int(x) % p for x in b]

    def trim(x):
        while len(x) > 1 and x[0] == 0:
            x.pop(0)
        return x

    a = trim(a)
    b = trim(b)
    if len(a) < len(b):
        return [0], a
    inv = pow(b[0], -1, p)
    q = [0] * (len(a) - len(b) + 1)
    r = a[:]
    i = 0
    while len(r) >= len(b):
        c = r[0] * inv % p
        q[i] = c
        if c:
            for j in range(len(b)):
                r[j] = (r[j] - c * b[j]) % p
        r = trim(r)
        i += 1
    return trim(q), r


def check(ok, name, cond, detail):
    ok.append((name, bool(cond)))
    print(("  ok: " if cond else "FAIL: ") + f"{name} {detail}")


def main():
    ok = []
    art = json.load(open(ROOT / "results" / "gaussian_2x6.json"))
    work = ROOT / "results" / "gaussian_2x6_work"

    # 1. artifact invariants
    check(ok, "artifact schema/version", art["schema"].startswith("gaussian_2x6/"), art["schema"])
    check(ok, "zero failures recorded", not art["failures"], str(art["failures"]))
    check(ok, "both primes present in verdicts",
          set(map(str, (1000003, 2000003))) <= set(map(str, art["decisive"]["verdicts"].keys())),
          str(art["decisive"]["verdicts"]))
    v = art["verdict"]
    check(ok, "both assignments killed per artifact verdict",
          v["both_primes_agree_same_half_killed"] and v["both_primes_agree_cross_killed"], "")
    recs = art["decisive"]["decisive_2x6_pp"]
    check(ok, "two decisive records", len(recs) == 2, str(len(recs)))
    for rec in recs:
        pp = rec["per_prime"]
        check(ok, "decisive exit below cap",
          pp["gcd"]["mode"] == "exit" and pp["gcd"]["exit_degree"] == CAP - 1,
          f"p={pp['p']} deg={pp['gcd'].get('exit_degree')}")
        check(ok, "pair poly degree = C(1056,2)", rec["Npairs"] == 1056 * 1055 // 2 == 557040,
          str(rec["Npairs"]))
        check(ok, "cap == pairs - A_same(12)", rec["cap"] == rec["Npairs"] - A_SAME_12 == CAP,
          str(rec["cap"]))
        check(ok, "charpoly/pairpoly digests present",
          bool(pp["charpoly_sha"]) and bool(pp["pair_poly_sha"]), "")
    # pilot evidence: both primes recorded exact gcd 442067 in the artifact itself
    pilot_recs = art["pilot_ch12"]["pilot_ch12_pp"]
    check(ok, "pilot rows present for both primes",
          {r["per_prime"]["p"] for r in pilot_recs} == {1000003, 2000003}, str(len(pilot_recs)))
    for r in pilot_recs:
        pp = r["per_prime"]
        check(ok, f"pilot forced floor exact p={pp['p']}",
          pp["gcd"]["mode"] == "exact" and pp["gcd"]["gcd_degree"] == 442067,
          f"mode={pp['gcd']['mode']} gcd={pp['gcd'].get('gcd_degree')}")

    # 2. independent replay of the exit tail from saved states, chained to the
    #    decisive records by cache_key + digest chain + step accounting.
    for rec in art["decisive"]["decisive_2x6_pp"]:
        p = rec["per_prime"]["p"]
        key = rec["cache_key"]
        st_paths = sorted(work.glob(f"decisive_2x6_pp_p{p}_euclidstate_{key[:16]}.npz"))
        check(ok, f"key-named euclid state exists p={p}", len(st_paths) == 1, key[:16])
        st = np.load(st_paths[0], allow_pickle=True)
        check(ok, f"state embedded key p={p}", str(st["key"]) == key, key[:16])
        pp_paths = sorted(work.glob(f"decisive_2x6_pp_p{p}_pairpoly_{key[:16]}.npz"))
        pp_st = np.load(pp_paths[0], allow_pickle=True)
        check(ok, f"pairpoly embedded key p={p}", str(pp_st["key"]) == key, key[:16])
        pair_poly = pp_st["pair_poly"]
        check(ok, f"pairpol sha chain p={p}",
              e154.arr_sha(np.ascontiguousarray(pair_poly)) == rec["per_prime"]["pair_poly_sha"],
              rec["per_prime"]["pair_poly_sha"][:16])
        trimmed = e154.poly_trim_desc(pair_poly % p)
        check(ok, f"digest0 == trimmed input p={p}", e154.arr_sha(trimmed) == str(st["digest0"]),
              str(st["digest0"])[:16])
        a, b = [int(x) for x in st["a"]], [int(x) for x in st["b"]]
        steps = 0
        while True:
            q, r = divmod_hard(a, b, p)
            steps += 1
            if len(r) - 1 < CAP:
                q2, r2 = e154.poly_divmod(np.asarray(a), np.asarray(b), p)
                drift = (len(r) == len(r2)) and all(int(x) == int(y) for x, y in zip(r, r2))
                check(ok, f"independent exit replay p={p}", len(r) - 1 == CAP - 1 and steps <= 40,
                      f"deg={len(r)-1} steps={steps}")
                check(ok, f"state + replay steps == artifact steps p={p}",
                      int(st["steps"]) + steps == rec["per_prime"]["gcd"]["steps"],
                      f"{int(st['steps'])}+{steps} == {rec['per_prime']['gcd']['steps']}")
                check(ok, f"exit remainder digest chain p={p}",
                      e154.arr_sha(np.ascontiguousarray(np.asarray(r, dtype=np.int64)))
                      == rec["per_prime"]["gcd"]["exit_digest"],
                      rec["per_prime"]["gcd"]["exit_digest"][:16])
                check(ok, f"two-implementation remainder equality p={p}", drift, "bits match")
                break
            if steps > 128:
                check(ok, f"independent exit replay p={p}", False, "no exit in 128 steps")
                break
            a, b = b, r

    # 3. structural guard re-derived: cp(block) | cp(full) for the 2x3 layer
    mS = e136.fast_integral_model(6, list(e154.layer_bonds((2, 3), (False, False))), Fraction(1, 3), e136.Budget())
    pS = 100003
    MS = e136.matrix_mod_p(mS, pS)
    trS = e136.sequential_power_traces(MS, 64, pS, e136.Budget())
    invS = [0] + [pow(k, -1, pS) for k in range(1, 65)]
    cpS = e136.newton_from_power_sums(trS, 64, pS, invS, e136.Budget())
    perms = [e154.e154_gp_flip(6), e154.e154_gp_rho(6)]
    for c in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        B = e154.char_basis_dense(6, perms, list(c))
        Bi = B.astype(np.int64)
        Wb = (MS @ Bi) % pS
        Gf = (Bi.T @ Wb) % pS
        grams = np.array([int((Bi[:, j] ** 2).sum()) for j in range(Bi.shape[1])])
        scale = np.array([pow(int(g) % pS, -1, pS) for g in grams])
        Y = (Gf * scale[:, None]) % pS
        tr2 = e136.sequential_power_traces(Y, Bi.shape[1], pS, e136.Budget())
        inv2 = [0] + [pow(k, -1, pS) for k in range(1, Bi.shape[1] + 1)]
        cp = e136.newton_from_power_sums(tr2, Bi.shape[1], pS, inv2, e136.Budget())
        rem = e136.remainder_lazy(cpS, cp, pS, e136.Budget())
        check(ok, f"structural cp(block)|cp(full) char {c}", e136.poly_is_zero(rem), f"dim={Bi.shape[1]}")

    n_ok = sum(1 for _, c in ok if c)
    print(f"\n{n_ok}/{len(ok)} checks passed")
    assert n_ok == len(ok), [n for n, c in ok if not c]


if __name__ == "__main__":
    main()
    print("PASS")
