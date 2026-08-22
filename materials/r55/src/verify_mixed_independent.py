#!/usr/bin/env python3
"""Independent re-verification of the mixed n=45 artifact's search layer.

This script deliberately imports NOTHING from the producers
(`mixed_subgraph_identities`, `mixed_deficiency_cone`, `search_mixed_cuts`,
`m3_deficiency_cone`, `m4_deficiency_cone`).  It uses only the Python
standard library, so every number it reports is recomputed from the
artifact bytes by a second, independent implementation.

Scope (deliberately stated, mirroring the m4 checker's scoping discipline):

  VERIFIED HERE
    1. schema/campaign/disposition/key sets and record shapes;
    2. SHA-256, byte size, and graph counts of every recorded input;
    3. exact reconstruction of every per-state F interval from the
       artifact's OWN motif windows, using this file's own affine-interval
       arithmetic and the plan's verbatim n=45 coefficients;
    4. the (d,a,b,g,h) projection against the frozen m4 artifact;
    5. catalog windows against their own histograms, and outer windows
       containing every catalog window;
    6. every certificate's sign gates, per-state affine inequality,
       least-alpha canonicality, and the exact 45*alpha bound;
    7. every primal witness's weights, the seven aggregate relations
       (including sum f_lo <= 0 <= sum f_hi), and its exact value;
    8. the derived route statuses and the terminal disposition.

  NOT VERIFIED HERE (remaining tier, Task 4 of the plan)
    the 8,500,211-graph catalog re-sweep with an independently written
    motif kernel, and the independent recomputation of the m3 g row and
    m4 h row endpoints from the R1-R5 relation system.

Exit status is 0 only when every check above passes.
"""

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

N = 45
EXPECTED_STATES = 3215
EXPECTED_CATALOG_GRAPHS = 8_500_211
SCHEMA = 3
CAMPAIGN = "higher_order_identity_positive_deficiency_mixed"
DISPOSITIONS = (
    "MIXED_ANALYZED", "MIXED_ACCEPTED_CUT", "MIXED_NO_CUT_IN_FROZEN_BASIS",
    "MIXED_CERTIFICATION_UNRESOLVED",
)
STREAM_KEYS = (
    "q", "t", "p3", "pc3", "i4", "k4", "diamond", "k3k1", "c4",
    "claw", "paw", "p4", "two_k2", "k2_2k1", "p3_k1",
)
# Acceptance edges, frozen independently of the producer.
EDGES = {
    "total_deficiency": (Fraction(315), True),
    "degree20_count": (Fraction(1), False),
    "deficiency_ge8_count": (Fraction(1), False),
}
COEFFS = ("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta",
          "theta")
# Published R(4,5) maximum edge counts at the five feasible degrees.
E45 = {20: 100, 21: 107, 22: 114, 23: 122, 24: 132}


def fail(message):
    print(f"MIXED EVIDENCE REJECTED: {message}", file=sys.stderr)
    raise SystemExit(1)


def c2(value):
    return value * (value - 1) // 2


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def objective_value(name, state):
    if name == "total_deficiency":
        return Fraction(state["deficiency"])
    if name == "degree20_count":
        return Fraction(1 if state["d"] == 20 else 0)
    if name == "deficiency_ge8_count":
        return Fraction(1 if state["deficiency"] >= 8 else 0)
    fail(f"unregistered objective {name!r}")


def canonical(text, label):
    if not isinstance(text, str):
        fail(f"{label} is not a string: {text!r}")
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        fail(f"{label} is not rational: {text!r}")
    if str(value) != text:
        fail(f"{label} {text!r} is not canonical")
    return value


def affine(coefficient, window):
    """This file's own sign-aware affine image of an integer interval."""
    low, high = window
    if low > high:
        fail(f"inverted window {window}")
    products = (coefficient * low, coefficient * high)
    return min(products), max(products)


def f_interval(state, windows):
    """Recompute one state's F interval from the artifact's own windows.

    Verbatim n=45 coefficients of the plan's Task 2 expansion.  The x
    stratum is the neighbourhood class, the z stratum is the complement of
    the anti-neighbourhood; e_y = C(44-d,2) - e_z.
    """
    d = state["d"]
    m = 44 - d
    e_x = E45[d] - state["a"]
    e_z = E45[m] - state["b"]
    e_y = c2(m) - e_z
    x_cls = f"{d},{e_x}"
    z_cls = f"{m},{e_z}"
    exact = (
        1890 * d - 2109 * d ** 2 + 135 * d ** 3 - 2 * d ** 4
        + 4124 * e_x - 12 * e_x ** 2 - 528 * d * e_x + 12 * d ** 2 * e_x
        + 4 * e_y ** 2 - 86 * d * e_y + 4 * d ** 2 * e_y
    )
    terms = (
        (516, "t", x_cls),
        (72, "c4", x_cls),
        (24, "claw", x_cls),
        (24, "p4", x_cls),
        (24, "paw", x_cls),
        (564 - 24 * d, "p3", x_cls),
        (32, "diamond", x_cls),
        (74 + 4 * d, "pc3", z_cls),
        (-12, "k3k1", z_cls),
        (-8, "two_k2", z_cls),
        (-8, "p3_k1", z_cls),
        (-24, "k2_2k1", z_cls),
    )
    low = high = 0
    for coefficient, key, cls in terms:
        if cls not in windows:
            fail(f"state {state['d'], state['a'], state['b']}: no windows "
                 f"for stratum {cls}")
        piece = affine(coefficient, windows[cls][key])
        low += piece[0]
        high += piece[1]
    return exact + low, exact + high


def main(argv):
    arguments = [value for value in argv[1:] if not value.startswith("--")]
    flags = {value for value in argv[1:] if value.startswith("--")}
    unknown = flags - {"--no-input-hashes"}
    if len(arguments) != 1 or unknown:
        print("usage: verify_mixed_independent.py [--no-input-hashes] "
              "<engstrom_identity.json>", file=sys.stderr)
        return 2
    skip_input_hashes = "--no-input-hashes" in flags
    target = Path(arguments[0]).resolve()
    math_root = target.parents[2]
    document = json.loads(target.read_text(encoding="ascii"))

    # ---- stage 1: schema, campaign, disposition, key sets ----------------
    if document.get("schema_version") != SCHEMA:
        fail("schema_version drifted")
    if document.get("campaign_id") != CAMPAIGN:
        fail("campaign_id drifted")
    disposition = document.get("disposition")
    if disposition not in DISPOSITIONS:
        fail(f"unknown disposition {disposition!r}")
    for key in ("inputs", "states", "searches", "trust_roots",
                "catalog_motif_histograms", "catalog_motif_windows",
                "outer_motif_windows", "replay", "n49_replay"):
        if key not in document:
            fail(f"missing top-level key {key!r}")
    roots = document["trust_roots"]
    if not isinstance(roots, list) or len(roots) != 8:
        fail("expected exactly eight trust roots")
    allowed = (
        "https://arxiv.org/abs/1002.4304",
        "https://arxiv.org/abs/2409.15709",
        "https://arxiv.org/abs/1703.08768",
        "https://users.cecs.anu.edu.au/~bdm/papers/r55.pdf",
        "https://users.cecs.anu.edu.au/~bdm/data/ramsey.html",
    )
    for root in roots:
        if set(root) != {"statement", "url"} or not root["statement"]:
            fail("malformed trust root")
        if root["url"] not in allowed:
            fail(f"trust root cites an unrecognized source: {root['url']!r}")

    # ---- stage 2: replay records ----------------------------------------
    if document["replay"] != {"published_graphs": 328, "complements": 328,
                              "residual_zero": 656}:
        fail("replay record drifted")
    if document["n49_replay"] != {"extremal_graphs": 2, "combos": 4,
                                  "row_values": [0, 144, 288, 432],
                                  "zero_corners": 1}:
        fail("n49 replay record drifted")

    # ---- stage 3: input hashes and sizes --------------------------------
    # Refusing to hash is a loud, explicit scope reduction, never a silent
    # skip: the caller must ask for it and the notice names what is unproved.
    if skip_input_hashes:
        for record in document["inputs"]:
            if set(record) != {"relative_path", "sha256", "bytes",
                               "graph_count"}:
                fail("input record keys drifted")
        print(f"inputs NOT re-hashed (--no-input-hashes): provenance of "
              f"{len(document['inputs'])} recorded files is UNVERIFIED in "
              f"this run; rerun inside the repository working tree to prove "
              f"it")
    else:
        for record in document["inputs"]:
            if set(record) != {"relative_path", "sha256", "bytes",
                               "graph_count"}:
                fail("input record keys drifted")
            path = math_root / record["relative_path"]
            if not path.is_file():
                fail(f"missing input {record['relative_path']}")
            if path.stat().st_size != record["bytes"]:
                fail(f"{record['relative_path']}: byte size drifted")
            if sha256(path) != record["sha256"]:
                fail(f"{record['relative_path']}: sha256 drifted")
        print(f"inputs re-hashed: {len(document['inputs'])} files")

    # ---- stage 4: catalog histograms, windows, outer containment --------
    catalog_windows = {}
    total_graphs = 0
    for record in document["catalog_motif_histograms"]:
        cls = f"{record['order']},{record['edges']}"
        total_graphs += record["graph_count"]
        derived = {}
        for key, entries in record["histograms"].items():
            values = [entry["value"] for entry in entries]
            counts = [entry["count"] for entry in entries]
            if values != sorted(values) or any(c <= 0 for c in counts):
                fail(f"{cls} {key}: malformed histogram")
            if sum(counts) != record["graph_count"]:
                fail(f"{cls} {key}: histogram count drift")
            derived[key] = (values[0], values[-1])
        catalog_windows[cls] = derived
    if total_graphs != EXPECTED_CATALOG_GRAPHS:
        fail(f"catalog total {total_graphs}, expected "
             f"{EXPECTED_CATALOG_GRAPHS}")
    for record in document["catalog_motif_windows"]:
        cls = f"{record['order']},{record['edges']}"
        claimed = {k: tuple(v) for k, v in record["windows"].items()}
        if claimed != catalog_windows.get(cls):
            fail(f"{cls}: claimed catalog windows differ from histograms")
    outer = {}
    for record in document["outer_motif_windows"]:
        cls = f"{record['order']},{record['edges']}"
        outer[cls] = {k: tuple(v) for k, v in record["windows"].items()}
        if set(outer[cls]) != set(STREAM_KEYS):
            fail(f"{cls}: outer window keys drifted")
    for cls, windows in catalog_windows.items():
        if cls not in outer:
            continue
        for key in STREAM_KEYS:
            lo, hi = outer[cls][key]
            clo, chi = windows[key]
            if lo > clo or chi > hi:
                fail(f"{cls} {key}: outer window does not contain catalog")
    print(f"catalog re-derived: {total_graphs} graphs across "
          f"{len(catalog_windows)} classes; outer windows contain all")

    # ---- stage 5: independent F endpoints, and the m4 g/h projection ----
    states = document["states"]
    if len(states) != EXPECTED_STATES:
        fail(f"{len(states)} states, expected {EXPECTED_STATES}")
    frozen_path = target.parent / "higher_identity_m4.json"
    frozen = {
        (s["d"], s["a"], s["b"]): s
        for s in json.loads(frozen_path.read_text(encoding="ascii"))["states"]
    }
    if len(frozen) != EXPECTED_STATES:
        fail("frozen m4 projection state count drifted")
    for state in states:
        key = (state["d"], state["a"], state["b"])
        record = frozen.get(key)
        if record is None:
            fail(f"state {key} is absent from the frozen m4 projection")
        for field in ("deficiency", "excess_balance", "g_lo", "g_hi",
                      "h_lo", "h_hi"):
            if state[field] != record[field]:
                fail(f"state {key}: {field} disagrees with frozen m4")
        if state["deficiency"] != state["a"] + state["b"]:
            fail(f"state {key}: deficiency is not a + b")
        low, high = f_interval(state, outer)
        if (low, high) != (state["f_lo"], state["f_hi"]):
            fail(f"state {key}: recomputed F interval ({low}, {high}) is not "
                 f"the recorded ({state['f_lo']}, {state['f_hi']})")
    print(f"states re-derived: {len(states)} F intervals recomputed from the "
          f"artifact's own windows; g/h/balance match the frozen m4 cone")

    # ---- stage 6: certificates, witnesses, statuses, disposition -------
    searches = document["searches"]
    if not searches:
        print("searches empty: artifact is still at MIXED_ANALYZED")
        print(f"MIXED EVIDENCE VERIFIED (analysis tier): {disposition}")
        return 0
    if len(searches) != 4:
        fail("expected exactly four search records")
    statuses = []
    for record in searches:
        name = record["objective_id"]
        if name == "required_local_family":
            if record["route_status"] != "UNAVAILABLE_NO_COVER_CERTIFICATE":
                fail("the unavailable route carries a status it cannot have")
            if any(record[field] is not None for field in
                   ("certificate", "exact_upper_bound", "primal_witness",
                    "exact_witness_value")):
                fail("the unavailable route carries evidence")
            continue
        limit, inclusive = EDGES[name]
        bound = None
        if record["certificate"] is not None:
            cert = record["certificate"]
            if set(cert) != set(COEFFS) | {"objective_id"}:
                fail(f"{name}: certificate arity drifted")
            values = {k: canonical(cert[k], f"{name}.{k}") for k in COEFFS}
            for gate in ("gamma", "delta", "epsilon", "zeta", "eta", "theta"):
                if values[gate] < 0:
                    fail(f"{name}: {gate} is negative")
            required = None
            for state in states:
                need = (
                    objective_value(name, state)
                    - values["beta"] * state["excess_balance"]
                    - values["gamma"] * state["g_lo"]
                    + values["delta"] * state["g_hi"]
                    - values["epsilon"] * state["h_lo"]
                    + values["zeta"] * state["h_hi"]
                    - values["eta"] * state["f_lo"]
                    + values["theta"] * state["f_hi"]
                )
                required = need if required is None else max(required, need)
            if values["alpha"] < required:
                fail(f"{name}: alpha {values['alpha']} violates a state")
            if values["alpha"] != required:
                fail(f"{name}: alpha {values['alpha']} is not the least "
                     f"admissible {required}")
            bound = N * values["alpha"]
            if canonical(record["exact_upper_bound"],
                         f"{name}.bound") != bound:
                fail(f"{name}: recorded bound is not the recomputed {bound}")
        value = None
        if record["primal_witness"] is not None:
            total = balance = _ = Fraction(0)
            aggregates = {k: Fraction(0) for k in
                          ("g_lo", "g_hi", "h_lo", "h_hi", "f_lo", "f_hi")}
            value = Fraction(0)
            previous = None
            for pair in record["primal_witness"]["weights"]:
                index = pair["state_index"]
                if previous is not None and index <= previous:
                    fail(f"{name}: witness indices do not increase")
                previous = index
                if not 0 <= index < len(states):
                    fail(f"{name}: witness index out of range")
                weight = canonical(pair["weight"], f"{name}.weight")
                if weight < 0:
                    fail(f"{name}: negative witness weight")
                state = states[index]
                total += weight
                balance += weight * state["excess_balance"]
                for field in aggregates:
                    aggregates[field] += weight * state[field]
                value += weight * objective_value(name, state)
            if total != N:
                fail(f"{name}: witness weights sum to {total}, not {N}")
            if balance != 0:
                fail(f"{name}: witness aggregate balance is {balance}")
            for lo_field, hi_field in (("g_lo", "g_hi"), ("h_lo", "h_hi"),
                                       ("f_lo", "f_hi")):
                if aggregates[lo_field] > 0:
                    fail(f"{name}: aggregate {lo_field} exceeds zero")
                if aggregates[hi_field] < 0:
                    fail(f"{name}: aggregate {hi_field} is below zero")
            if canonical(record["exact_witness_value"],
                         f"{name}.value") != value:
                fail(f"{name}: recorded witness value is not the recomputed "
                     f"{value}")
        accepts = (
            bound is not None
            and (bound <= limit if inclusive else bound < limit)
        )
        rejects = (
            value is not None
            and not (value <= limit if inclusive else value < limit)
        )
        if accepts and rejects:
            fail(f"{name}: contradictory evidence")
        derived = ("ACCEPTED_EXACT_CUT" if accepts else
                   "REJECTED_BY_EXACT_WITNESS" if rejects else
                   "CERTIFICATION_UNRESOLVED")
        if record["route_status"] != derived:
            fail(f"{name}: recorded status {record['route_status']!r} is not "
                 f"the derived {derived}")
        statuses.append(derived)
        print(f"route {name}: bound={bound} witness={value} status={derived}")
    if len(statuses) != 3:
        fail("expected three instantiable routes")
    if any(s == "ACCEPTED_EXACT_CUT" for s in statuses):
        terminal = "MIXED_ACCEPTED_CUT"
    elif all(s == "REJECTED_BY_EXACT_WITNESS" for s in statuses):
        terminal = "MIXED_NO_CUT_IN_FROZEN_BASIS"
    else:
        terminal = "MIXED_CERTIFICATION_UNRESOLVED"
    if terminal != disposition:
        fail(f"disposition {disposition!r} is not the derived {terminal}")
    print(f"MIXED EVIDENCE VERIFIED: {disposition}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
