"""Machine checks that reports/paper_pbb_nogo.md matches its artifacts.

The artifacts are the single source of truth; the paper is a derived
presentation of them.  A paper number that no longer re-derives is a defect in
the paper, so every headline figure it asserts is parsed back out of the prose
here and compared against the JSON it claims to come from.

This is the guard that catches the failure mode we actually hit: a long-running
sweep advances, the artifact grows, and the prose silently keeps a stale count.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "reports" / "paper_pbb_nogo.md"
NOGO = ROOT / "results" / "partial_runs" / "exp039_nogo_module.json"
ENVELOPE = ROOT / "results" / "processed" / "exp036_envelope_check.json"
TRICHOTOMY = ROOT / "results" / "processed" / "exp052_ideal_power_trichotomy.json"
IDEAL = ROOT / "results" / "processed" / "exp053_ideal_classification.json"
CENSUS = ROOT / "results" / "processed" / "exp054_mixed_census.json"
ODDSWEEP = ROOT / "results" / "processed" / "exp055_odd_lattice_sweep.json"
ODDLIT = ROOT / "results" / "processed" / "exp055_literature_validation.json"
ODDSCREEN = ROOT / "results" / "processed" / "exp055_odd_lattice_screen.json"
ODDDISC = ROOT / "results" / "certificates" / "exp055_discovered_references.json"


@pytest.fixture(scope="module")
def paper() -> str:
    if not PAPER.exists():
        pytest.skip("paper not present")
    return PAPER.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def nogo() -> dict:
    if not NOGO.exists():
        pytest.skip("EXP-039 artifact not present")
    return json.loads(NOGO.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def trichotomy() -> dict:
    assert TRICHOTOMY.exists(), f"required artifact missing: {TRICHOTOMY}"
    d = json.loads(TRICHOTOMY.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp052-ideal-power-trichotomy-v1-assembled", d.get("schema")
    return d


@pytest.fixture(scope="module")
def ideal() -> dict:
    assert IDEAL.exists(), f"required artifact missing: {IDEAL}"
    d = json.loads(IDEAL.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp053-ideal-classification-v1", d.get("schema")
    return d


@pytest.fixture(scope="module")
def oddsweep() -> dict:
    assert ODDSWEEP.exists(), f"required artifact missing: {ODDSWEEP}"
    d = json.loads(ODDSWEEP.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-odd-lattice-sweep-v1", d.get("schema")
    return d


@pytest.fixture(scope="module")
def oddlit() -> dict:
    assert ODDLIT.exists(), f"required artifact missing: {ODDLIT}"
    d = json.loads(ODDLIT.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-literature-v3", d.get("schema")
    return d


@pytest.fixture(scope="module")
def oddscreen() -> dict:
    assert ODDSCREEN.exists(), f"required artifact missing: {ODDSCREEN}"
    d = json.loads(ODDSCREEN.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-odd-lattice-screen-v2", d.get("schema")
    return d


@pytest.fixture(scope="module")
def odddisc() -> dict:
    assert ODDDISC.exists(), f"required artifact missing: {ODDDISC}"
    d = json.loads(ODDDISC.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp055-discovered-references-v1", d.get("schema")
    return d


def test_section_62_numbers_match_artifacts(paper: str, oddsweep: dict,
                                            oddlit: dict) -> None:
    """Section 6.2's counts must re-derive from EXP-055, not drift with prose."""
    v = oddsweep["verdict"]
    assert re.search(rf"all \${v['lattices_swept']}\$ odd lattices", paper)
    # 4.23e9 pairs, quoted to three significant figures
    assert f"${v['pairs_total'] / 1e9:.2f}" in f"${4.23:.2f}"
    assert re.search(r"4\.23\\times10\^\{9\}", paper)
    assert abs(v["pairs_total"] / 1e9 - 4.23) < 0.005, v["pairs_total"]
    assert re.search(rf"\(\${v['idempotence_tested']}\$ idempotence tests", paper)
    assert v["k_mismatches"] == 0 and v["idempotence_violations"] == 0
    # literature battery
    assert re.search(rf"we sourced \${oddlit['instances']}\$ instances", paper)
    assert oddlit["pole_isomorphism_all"] is True
    assert re.search(rf"passes\s+\*\*{oddlit['instances']}/{oddlit['instances']}\*\*",
                     paper)
    assert re.search(rf"printed \$k\$ on \${oddlit['reproduced_instances']}\$",
                     paper)
    assert oddlit["reported_ceiling_sanity_holds"] is True
    assert oddlit["exact_ceiling_never_violated"] is True
    assert re.search(r"Eight rows\s+are independently exact-certified", paper)
    assert re.search(r"BP-OSD `distance_upperbound`", paper)
    triplet = (oddlit["exact_slack_min"], oddlit["exact_slack_median"],
               oddlit["exact_slack_max"])
    assert re.search(rf"slack min/median/max "
                     rf"\${triplet[0]}/{triplet[1]}/{triplet[2]}\$", paper), triplet


def test_section_62_screen_and_discoveries_match_artifacts(
        paper: str, oddscreen: dict, odddisc: dict) -> None:
    v, scope = oddscreen["verdict"], oddscreen["scope"]
    assert v["complete"] and v["all_referenced_decided"]
    assert v["all_referenced_dominated"]
    assert re.search(rf"Exact \$n={scope['n_max']}\$ closure", paper)
    assert f"**{v['candidates_after_symmetry']:,}**" in paper
    assert f"**{v['orbits_represented']:,}**" in paper
    assert f"**{v['dominated']:,}/{v['with_reference']:,}**" in paper
    assert (
        f"{v['verdicts']['dominated_by_witness']:,} "
        "direct physical witnesses"
    ) in paper
    assert re.search(
        rf"{v['verdicts']['dominated_by_automorphism_transport']:,}\s+"
        rf"full-automorphism\s+witness transports",
        paper,
    )
    assert (
        f"{v['verdicts']['dominated_by_cdcl_witness']:,} "
        "bounded CDCL witnesses"
    ) in paper
    assert re.search(
        rf"{v['verdicts']['dominated']:,}\s+exact CP-SAT\s+fallbacks",
        paper,
    )
    assert re.search(rf"{v['no_reference']:,}\s+are `no_reference`", paper)
    assert v["survivors"] == 0 and v["undecided"] == 0
    assert "zero survivors or undecided" in paper

    params = {(r["n"], r["k"], r["d"]) for r in odddisc["records"]}
    assert params == {(30, 8, 4), (54, 8, 6), (126, 12, 10)}
    for n, k, d in params:
        assert f"[[{n},{k},{d}]]" in paper


@pytest.fixture(scope="module")
def census() -> dict:
    assert CENSUS.exists(), f"required artifact missing: {CENSUS}"
    d = json.loads(CENSUS.read_text(encoding="utf-8"))
    assert d.get("schema") == "exp054-mixed-census-v1", d.get("schema")
    return d


def test_census_pair_count_in_prose_matches_artifact(paper: str, census: dict) -> None:
    """The 653,022,021-pair census figure must re-derive from EXP-054."""
    quoted = {int(s.replace("{,}", "").replace(",", ""))
              for s in re.findall(r"\$?([0-9]{3}\{,\}[0-9]{3}\{,\}[0-9]{3})\$?", paper)}
    assert census["verdict"]["pairs_total"] in quoted, sorted(quoted)
    assert census["verdict"]["lattices_with_mixed"] == []
    assert census["verdict"]["crosscheck_mismatches"] == 0
    assert f"${census['verdict']['lattices_scanned']}$ lattices" in paper


def test_ideal_invariance_prose_matches_artifact(paper: str, ideal: dict, trichotomy: dict) -> None:
    """Section 6.1's counts are the EXP-053 verdict, and dim S = 2 dim I^infty
    is re-derived here against EXP-052's independent module chains."""
    v = ideal["verdict"]
    assert v["cases"]["demote_full"] == 192 and v["cases"]["immune"] == 10
    assert v["cases"]["mixed"] == 0
    assert v["ideal_route_matches_exp052"] and v["coset_criterion_exact"]
    assert v["nilpotency_index_2_on_all_demoting"] and v["odd_corollary_holds"]
    ref = trichotomy["records"]
    for r in ideal["catalogue"]["records"]:
        assert r["dim_S_pred"] == 2 * r["dim_I_infty"] == ref[r["fingerprint"]]["dim_S"]
    # prose invariants that must not drift
    assert re.search(r"\$I\^2=0\$\s+on all \$192\$\s+demoting parents", paper)
    assert re.search(r"\$I\^2=I\$\s+on all \$10\$\s+immune ones", paper)
    assert re.search(r"\$\[\[90,8,10\]\]\$ code on the \$\(15,3\)\$\s+lattice", paper)


def test_dual_ideal_prose_matches_artifact(paper: str, ideal: dict) -> None:
    """Section 6.1's dual-ideal audit numbers must re-derive from EXP-053."""
    da = ideal["duality_audit"]
    assert da["audit_holds"] and da["annihilators_interchangeable"] == []
    bar_invariant = sum(bool(r["ideal_AB_is_bar_invariant"]) for r in da["records"])
    assert re.search(rf"bar-invariant on only \${bar_invariant}\$", paper), bar_invariant
    imm = [r for r in da["records"] if r["case"] == "immune"]
    # the wrong route returns 2(dim R - dim I) on every immune parent
    for r in imm:
        assert r["wrong_route_dim_S"] == 2 * (r["dim_R"] - r["dim_I"]), r["label"]
    quoted = sorted({r["wrong_route_dim_S"] for r in imm})
    assert quoted == [100, 172, 344], quoted
    for v in quoted:
        assert f"${v}$" in paper, v
    assert re.search(r"2\(\\ell m-k_P/2\)", paper)
    assert re.search(r"\\operatorname\{Ann\}_R\(M\)=\(\\bar a,\\bar b\)", paper)


def test_mixed_witness_numbers_in_prose(paper: str, ideal: dict) -> None:
    w = ideal["mixed_witness"]
    assert (w["case"], w["k_parent"], w["dim_S_pred"]) == ("mixed", 8, 4)
    assert w["census"]["classes_demote"] == 240 and w["census"]["classes_total"] == 255
    assert w["agree_ideal_vs_module"] and w["agree_ideal_vs_census"]
    assert "$240$ of $255$" in paper
    assert "$2^8-2^4=240$" in paper


def _one(pattern: str, text: str) -> str:
    """Exactly one match, or the claim is ambiguous and the test must fail."""

    found = re.findall(pattern, text)
    assert len(found) == 1, f"pattern {pattern!r} matched {len(found)} times"
    return found[0]


def test_summary_table_matches_nogo_artifact(paper: str, nogo: dict) -> None:
    summary = nogo["counts"]
    table = {
        "distinct parents": r"\| distinct parents \| \$(\d+)\$ \|",
        "certified": r"parents with \*\*exact\*\* \$T\$[^|]*\| \$(\d+)\$ \|",
        "closed": r"parents \*\*family-closed\*\*[^|]*\| \$\\mathbf\{(\d+)\}\$ \|",
        "capped": r"catalogue rows capped \*a priori\*[^|]*\| \$\\mathbf\{(\d+)\}\$ of \$368\$ \|",
        "permitted": r"catalogue rows Theorem H permits[^|]*\| \$(\d+)\$ \|",
    }
    assert int(_one(table["distinct parents"], paper)) == summary["distinct_parents_total"]
    assert int(_one(table["certified"], paper)) == summary["parents_T_exact"]
    assert int(_one(table["closed"], paper)) == summary["parents_family_closed"]
    assert int(_one(table["capped"], paper)) == summary["catalogue_rows_capped_a_priori"]
    assert int(_one(table["permitted"], paper)) == summary["catalogue_rows_theorem_permits"]


def test_abstract_counts_agree_with_summary_table(paper: str, nogo: dict) -> None:
    # The abstract restates the two headline counts in prose; both must match.
    summary = nogo["counts"]
    assert int(_one(r"exactly for \$(\d+)\$ of the \$202\$ distinct parents", paper)) == (
        summary["parents_T_exact"]
    )
    assert int(_one(r"\*\*\$(\d+)\$ parents are family-closed\*\*", paper)) == (
        summary["parents_family_closed"]
    )
    assert int(_one(r"settles \$(\d+)\$ of the \$368\$ catalogue rows", paper)) == (
        summary["catalogue_rows_capped_a_priori"]
    )


def test_per_length_table_matches_artifact(paper: str, nogo: dict) -> None:
    # The per-length table lists parents / certified / closed / rows / capped.
    # Parents and rows are parent-side facts recomputed from the experiment;
    # certified, closed and capped must match the artifact's by_n ledger.
    import collections
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    parents = module.distinct_parents(module.E27.load_catalogue())

    tot: collections.Counter = collections.Counter()
    rws: collections.Counter = collections.Counter()
    for entry in parents.values():
        n = entry["n"]
        tot[n] += 1
        rws[n] += len(entry["members"])

    for n_raw, stats in nogo["by_n"].items():
        n = int(n_raw)
        expected = (
            rf"\n\| \${n}\$ \| \${tot[n]}\$ \| "
            rf"\${stats['T_exact']}/{tot[n]}\$ \| "
            rf"\${stats['family_closed']}\$ \| "
            rf"\${rws[n]}\$ \| "
        )
        assert re.search(expected, "\n" + paper), (
            f"n={n}: paper table does not state {expected!r}"
        )

    # The bold n<=144 subtotal row is the sum and is what the prose quotes.
    le144 = [n for n in tot if n <= 144]
    parents_total = sum(tot[n] for n in le144)
    certified_total = sum(int(nogo["by_n"][str(n)]["T_exact"]) for n in le144)
    closed_total = sum(int(nogo["by_n"][str(n)]["family_closed"]) for n in le144)
    rows_total = sum(rws[n] for n in le144)
    subtotal = _one(
        r"\n\| \*\*\$\\le 144\$\*\* \| \*\*\$(\d+)\$\*\* \| \*\*\$(\d+)/(\d+)\$\*\* \| "
        r"\*\*\$(\d+)\$\*\* \| \*\*\$(\d+)\$\*\*", "\n" + paper
    )
    assert int(subtotal[0]) == parents_total
    assert int(subtotal[1]) == certified_total == int(subtotal[2]) == parents_total
    assert int(subtotal[3]) == closed_total
    assert int(subtotal[4]) == rows_total


def test_confinement_table_matches_recomputation(paper: str, nogo: dict) -> None:
    # Section 5.1 asserts a row count per (k_Q/k_P, T/k_P) cell and that every
    # cell above k_Q/k_P = 1/2 is capped.  Recompute both from the artifact.
    from fractions import Fraction

    rows = nogo["rows_capped_a_priori"] + nogo["rows_theorem_permits"]
    cross: dict[tuple[str, str], int] = {}
    for item in rows:
        kq, kp, T = int(item["k_pbb"]), int(item["k_parent"]), int(item["T"])
        key = (str(Fraction(kq, kp)), str(Fraction(T, kp)))
        cross[key] = cross.get(key, 0) + 1

    total = int(_one(r"against their \$(\d+)\$ catalogue rows", paper))
    assert total == len(rows)

    for (kq_kp, t_kp), count in cross.items():
        # Cell values appear as bare integers in the row for this k_Q/k_P.
        row = _one(rf"\n\| \${re.escape(kq_kp)}\$ \|([^\n]*)\n", paper)
        cells = re.findall(r"\$?\\?m?a?t?h?b?f?\{?(\d+)\}?\$?", row)
        assert str(count) in cells, (
            f"cell (k_Q/k_P={kq_kp}, T/k_P={t_kp}) = {count} absent from row {row!r}"
        )

    # The structural claim itself, not just its tabulation.
    for item in rows:
        kq, kp, T = int(item["k_pbb"]), int(item["k_parent"]), int(item["T"])
        if Fraction(kq, kp) > Fraction(1, 2):
            assert kq > kp - T, f"{item['label']}: keeps >1/2 of k but is not capped"


def test_dimension_upper_law_never_exceeds_T(paper: str, nogo: dict) -> None:
    # Conjecture A's second evidence layer: on every classified row the
    # dressing dimension k_P - k_Q (Theorem G(ii), verified 368/368) never
    # exceeds T(P), and the rows at equality are exactly the rows Theorem H
    # permits.  A violation of either half is a structural discovery, not a
    # bookkeeping slip, and must fail loudly here.
    capped = nogo["rows_capped_a_priori"]
    permitted = nogo["rows_theorem_permits"]
    for item in capped + permitted:
        dim_delta = int(item["k_parent"]) - int(item["k_pbb"])
        T = int(item["T"])
        assert dim_delta <= T, (
            f"{item['label']}: dim(Delta_bar)={dim_delta} exceeds T={T} -- "
            "the catalogue upper law is broken"
        )
    # The paper states the upper law twice (evidence (ii) and the spread note);
    # require it verbatim at least once rather than uniquely.
    assert len(re.findall(r"zero exceptions", paper)) >= 1
    assert f"${len(capped) + len(permitted)}$ classified rows" in paper


def test_gate_claims_match_artifact(paper: str) -> None:
    # Section 5.2 claims 7/7 checked, 7/7 non-vacuous, zero violations, slack 0.
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    gate = module.falsification_gate(
        module.distinct_parents(module.E27.load_catalogue()),
        module.load_certificates(),
    )
    checked = int(_one(r"\$(\d+)/7\$ checked", paper))
    assert checked == gate["reversals_checked"] == 7
    assert gate["reversals_non_vacuous"] == 7
    assert not gate["theorem_refuted"] and not gate["violations"]
    assert gate["all_saturate_bound"]
    assert all(c["slack"] == 0 for c in gate["checks"])
    assert "slack exactly $0$ in every" in paper


def test_envelope_table_matches_envelope_artifact(paper: str) -> None:
    if not ENVELOPE.exists():
        pytest.skip("envelope artifact not present")
    env = json.loads(ENVELOPE.read_text(encoding="utf-8"))
    assert env["schema"] == "exp036-envelope-check-v2"
    assert f'schema\n`{env["schema"]}`' in paper
    assert env["all_reversals_css_dominated_at_equal_n"]
    assert int(_one(r"against \$(\d+)\$ certified-exact CSS candidates", paper)) == (
        env["num_exact_css_candidates"]
    )
    # Every reversal's kd^2/n pair and dominator must appear as the paper states.
    for check in env["checks"]:
        dom = check["dominator"]
        pbb = rf"\$\[\[{check['n']},{check['k_pbb']},{check['d_pbb_upper']}\]\]\$"
        css = rf"\$\[\[{check['n']},{dom['k_css']},{dom['d_css_exact']}\]\]\$"
        assert re.search(pbb, paper), f"{check['label']}: PBB triple {pbb} absent"
        assert re.search(css, paper), f"{check['label']}: dominator {css} absent"
        assert f"${check['kd2_over_n_pbb']:.2f}$" in paper
        assert f"${check['kd2_over_n_css']:.2f}$" in paper
    factors = [c["kd2_over_n_css"] / c["kd2_over_n_pbb"] for c in env["checks"]]
    span = _one(r"a \$(\d+\.\d+)\$–\$(\d+\.\d+)\\times\$ deficit", paper)
    lo, hi = float(span[0]), float(span[1])
    # Quoted to one decimal; require the true values to round-trip within the
    # stated endpoints (rounding pushes the endpoints outward, never inward).
    assert lo <= min(factors) and round(max(factors), 1) <= hi, (
        f"stated {lo}-{hi}x does not bracket observed "
        f"{min(factors):.2f}-{max(factors):.2f}x"
    )


def test_gross_parent_row_is_the_literal_gross_code(paper: str) -> None:
    # The headline corollary names the Gross code; the row must be the real one.
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    E27 = module.E27
    rows = E27.load_catalogue()
    target = None
    for index, row in enumerate(rows):
        if row["A_terms"] == [[3, 0], [0, 1], [0, 2]] and row["B_terms"] == [
            [0, 3],
            [1, 0],
            [2, 0],
        ]:
            target = (index, row)
            break
    assert target is not None, "Gross polynomials absent from the catalogue"
    index, row = target
    assert (int(row["ell"]), int(row["m"])) == (12, 6)
    label = E27.catalogue_label(row, index)
    assert f"| `{label}` |" in paper

    _, HX, HZ = E27.parent_matrices(row)
    fingerprint = E27.matrix_fingerprint(HX, HZ)
    certs = module.load_certificates()
    cert = certs.get(fingerprint)
    if cert is None:
        pytest.skip("Gross parent not yet certified in this sweep")
    assert cert["k_parent"] == 12 and cert["d_z_parent"] == 12
    assert cert["T"] == 12 and cert["T_is_exact"] and cert["family_closed"]
    assert f"| `{fingerprint[:16]}` |" in paper


def test_trichotomy_claims_match_exp052_artifact(paper: str, trichotomy: dict) -> None:
    # EXp-052 headline: 202 parents classified by the ideal-power route,
    # 192 demote-full / 10 immune / 0 mixed, chains <= 3, 196 exact crosschecks.
    assert int(_one(r"on all \$(\d+)\$\s+catalogue parents", paper)) == trichotomy["parents_total"]
    pair = _one(r"\(\$(\d+)/(\d+)\$, zero mixed\)", paper)
    dm, im = int(pair[0]), int(pair[1])
    cases = trichotomy["cases"]
    assert dm == cases["demote_full"]
    assert im == cases["immune"]
    assert trichotomy["mixed_labels"] == []
    stated_chain = int(_one(r"chain length \$\\le (\d+)\$", paper))
    assert stated_chain == trichotomy["max_chain_steps"]
    assert int(_one(r"over the \$(\d+)\$ parents with \$k_P\\le 20\$", paper)) == (
        trichotomy["crosschecked_fractions"]
    )
    assert trichotomy["crosschecked_fractions"] == 196
    assert trichotomy["fraction_mismatches_vs_exp050"] == []
    assert trichotomy["immune_mismatch_vs_exp047"] == []


def test_no_stale_test_count(paper: str) -> None:
    # The suite count is quoted in the provenance section; keep it honest by
    # requiring it to be at least as large as the tests that exist.
    stated = int(_one(r"Full suite: \$(\d+)\$ passing tests", paper))
    collected = sum(
        len(re.findall(r"^def test_", path.read_text(encoding="utf-8"), re.M))
        for path in sorted((ROOT / "tests").glob("test_*.py"))
    )
    assert stated >= collected, (
        f"paper claims {stated} passing tests but {collected} test functions exist"
    )
