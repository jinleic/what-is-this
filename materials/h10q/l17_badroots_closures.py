#!/usr/bin/env python3
"""Bad signed root tables for the closure-authority progressions.

The p-adic engine is the checked-in l17_ratemodel_padic.py in this
directory.  Here it is run with every base prime p<=10,000 plus the >10,000
primes that actually occur in l17_sieve.jsonl for that exact cell.  The
latter targeted scan is needed because the sieve cutoff is 100,000.  Each
root coordinate is k mod p for the closure record's own (eps,f,q1,N), never
a legacy cohort row.
"""
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = os.environ.get("H10Q_ROOT", str(HERE))
ENGINE_PATH = os.environ.get("L17_PADIC_ENGINE", str(HERE / "l17_ratemodel_padic.py"))
AUTHORITY = Path(ROOT) / "data/l13h_all_closures.json"
SIEVE = Path(ROOT) / "data/l17_sieve.jsonl"
OUT = Path(os.environ.get("L17_OUT", str(Path(ROOT) / "data/l17_badroots_closures.jsonl")))
BASE_LIMIT = int(os.environ.get("L17_P_LIMIT", "10000"))
SIEVE_LIMIT = 100000

spec = importlib.util.spec_from_file_location("l17_padic_engine", ENGINE_PATH)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)


def cell_key(cell):
    return (cell[0], tuple(cell[1]))


def load_sieve_observations():
    observed = []
    by_cell = {}
    with SIEVE.open() as fh:
        for line in fh:
            row = json.loads(line)
            if "cell" not in row or "members" not in row:
                continue
            cell = cell_key(row["cell"])
            by_cell.setdefault(cell, set())
            for member in row["members"]:
                k = int(member["k"])
                for p, symbol in member.get("em_small", []):
                    if int(symbol) == -1:
                        p = int(p)
                        observed.append((cell, k, p))
                        by_cell[cell].add(p)
    return observed, by_cell


def compact_row(src, sf, table, smeta, observed_primes):
    lists = smeta["bad_root_lists"]
    ftable = [[q["p"], q["local_mass_num"], q["local_mass_den"]]
              for q in lists if not q["excluded"]]
    excluded_bound = sum(
        float(root.get("residual_mass_upper_bound_float", 0.0))
        for q in lists for root in q["roots"]
        if q["excluded"] and root.get("status") == "excluded-unresolved-singular-branch"
    )
    return {
        "type": "class",
        "cell": src["cell"],
        "family": src["family"],
        "class_index": src.get("class_index"),
        "params": {
            "a": src["a"], "eps": src["eps"], "f": src["f"],
            "q1": src["q1"], "N": src["N"],
        },
        "root_coordinate": "k_mod_p",
        "observed_k_zero": src.get("k_zero"),
        "small_factor_padic": sf,
        "f_p_table": ftable,
        "bad_root_lists": lists,
        "n_bad_prime_lists": len(lists),
        "n_bad_root_residues": sum(len(q["roots"]) for q in lists),
        "n_simple_roots": smeta["n_simple_roots"],
        "n_step_roots": smeta["n_step_roots"],
        "n_bad_simple_roots": smeta["n_bad_simple_roots"],
        "n_bad_step_roots": smeta["n_bad_step_roots"],
        "n_excluded_pairs": smeta["n_excluded_pairs"],
        "partial": smeta["partial"],
        "excluded_mass_upper_bound": excluded_bound,
        "sieve_observed_prime_count": len(observed_primes),
        "targeted_sieve_primes_above_base": sorted(p for p in observed_primes if p > BASE_LIMIT),
        "scan_prime_count": len(table),
    }


def main():
    authority = json.loads(AUTHORITY.read_text())
    records = authority["records"]
    class_limit = int(os.environ.get("L17_CLASS_LIMIT", str(len(records))))
    records = records[:class_limit]
    observed, observed_by_cell = load_sieve_observations()
    allowed_cells = {cell_key(r["cell"]) for r in records}
    observed = [x for x in observed if x[0] in allowed_cells]
    observed_by_cell = {c: ps for c, ps in observed_by_cell.items() if c in allowed_cells}
    base_primes = tuple(engine.primerange(2, BASE_LIMIT + 1))
    all_sieve_primes = set(engine.primerange(2, SIEVE_LIMIT + 1))
    bad_target_claims = sorted({p for _, _, p in observed if p > BASE_LIMIT})
    assert set(bad_target_claims) <= all_sieve_primes

    rows = []
    for i, src in enumerate(records, 1):
        cell = cell_key(src["cell"])
        targets = sorted(p for p in observed_by_cell.get(cell, set())
                         if p > BASE_LIMIT)
        scan_primes = tuple(base_primes) + tuple(targets)
        # The imported engine's only global scan input is PRIMES.
        engine.PRIMES = scan_primes
        engine.LIMIT = BASE_LIMIT
        sf, table, smeta = engine._class_small_factor(src)
        row = compact_row(src, sf, table, smeta, observed_by_cell.get(cell, set()))
        rows.append(row)
        print(f"roots {i}/{len(records)} {src['family']} {src['cell']} "
              f"scan={len(scan_primes)} badp={len(row['bad_root_lists'])}", flush=True)
    # Reconcile every observed (cell,k,p) odd-minus-one occurrence against
    # a listed bad-signed residue r == k mod p in that exact progression.
    listed = set()
    for row in rows:
        cell = cell_key(row["cell"])
        for q in row["bad_root_lists"]:
            p = int(q["p"])
            for root in q["roots"]:
                listed.add((cell, p, int(root["r"])))
    covered = [x for x in observed if (x[0], x[2], x[1] % x[2]) in listed]
    missing = [x for x in observed if (x[0], x[2], x[1] % x[2]) not in listed]
    print(f"RECONCILIATION covered={len(covered)}/{len(observed)} "
          f"missing={len(missing)}", flush=True)
    if missing:
        print("MISSING_SAMPLE", missing[:20], flush=True)
        raise AssertionError(("reconciliation", len(missing), missing[:5]))

    total_lists = sum(len(r["bad_root_lists"]) for r in rows)
    total_roots = sum(r["n_bad_root_residues"] for r in rows)
    total_excluded = sum(r["n_excluded_pairs"] for r in rows)
    meta = {
        "type": "meta",
        "artifact": "l17_badroots_closures",
        "authority": "data/l13h_all_closures.json",
        "sieve": "data/l17_sieve.jsonl",
        "root_coordinate": "k_mod_p",
        "n_classes": len(rows),
        "base_prime_cap": BASE_LIMIT,
        "sieve_prime_cap": SIEVE_LIMIT,
        "base_prime_count": len(base_primes),
        "targeted_sieve_primes_above_base": bad_target_claims,
        "target_scan_mode": "all p<=10000 per class plus observed sieve p>10000 for that cell",
        "p_adic_mass": "simple bad root 1/(p+1); singular roots recursive-lifted or exclusion-labeled",
        "reconciliation": {
            "observed_odd_minus_one_occurrences": len(observed),
            "covered": len(covered),
            "missing": len(missing),
            "coverage": 1.0 if not observed else len(covered) / len(observed),
        },
    }
    summary = {
        "type": "summary",
        "n_classes": len(rows),
        "n_bad_prime_lists": total_lists,
        "n_bad_root_residues": total_roots,
        "n_excluded_pairs": total_excluded,
        "reconciliation": meta["reconciliation"],
        "root_coordinate": "k_mod_p",
    }
    with OUT.open("w") as fh:
        fh.write(json.dumps(meta, separators=(",", ":")) + "\n")
        for row in rows:
            fh.write(json.dumps(row, separators=(",", ":")) + "\n")
        fh.write(json.dumps(summary, separators=(",", ":")) + "\n")
    print(json.dumps({"out": str(OUT), "classes": len(rows), "lists": total_lists,
                      "roots": total_roots, "excluded": total_excluded,
                      "reconciliation": meta["reconciliation"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
