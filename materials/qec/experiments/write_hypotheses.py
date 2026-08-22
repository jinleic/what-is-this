"""Regenerate notes/hypotheses.csv with guaranteed-correct CSV quoting."""
from __future__ import annotations

import csv
from pathlib import Path

COLS = ["id", "date", "research_track", "precise_hypothesis", "motivation",
        "mathematical_predictions", "empirical_predictions",
        "cheapest_decisive_test", "development_data", "held_out_data",
        "evidence_for", "evidence_against", "status", "artifact_paths",
        "next_action"]

R = [
 dict(id="H001", research_track="C",
  precise_hypothesis="The PBB commutation condition is symmetry of M = A C^T + B D^T over GF(2), not M = 0.",
  motivation="An internal summary of the seed paper stated M = 0 while the reference implementation tests symmetry. Only one can be right, and the primary source must decide it.",
  mathematical_predictions="H_x H_z^T = [[A C^T + B D^T ; AB+BA],[0;0]]. AB=BA because R is commutative, so validity <=> M + M^T = 0 <=> M symmetric. M=0 is strictly stronger and would exclude valid codes.",
  empirical_predictions="All catalogue codes satisfy symmetry; some violate M = 0, so a source stating M = 0 would be wrong.",
  cheapest_decisive_test="Derive the symplectic Gram matrix analytically; then evaluate both predicates on all 368 catalogue codes using an independent construction.",
  development_data="qcode-discovery campaign7 catalogue (368 codes)",
  held_out_data="7 Bravyi BB instances (C=D=0)",
  evidence_for="Analytic derivation (Lemma 0). Independent rebuild of all 368 codes: zero commutation failures under the symmetry test. arXiv:2606.02418v1 Sec III.2 states the symmetry condition verbatim, so the paper agrees with our derivation.",
  evidence_against="none", status="CONFIRMED",
  artifact_paths="proofs/pbb_structure.md; src/qec_research/codes/bicycle.py",
  next_action="Done: primary source quoted. No discrepancy exists; our own paraphrase was at fault (see failed_routes FR-006)."),
 dict(id="H002", research_track="C",
  precise_hypothesis="k(PBB) = k(BB(A,B)) - delta, with delta = dim of {u[C D] : u H_X = 0} modulo rowspace(H_Z).",
  motivation="If true the perturbation can only remove logical qubits, never add them, which bounds the whole PBB design space.",
  mathematical_predictions="rank(H_PBB) = rank(H_X) + rank(H_Z) + delta, via the projection (u,w)H -> u H_X whose kernel is exactly the pure-Z subgroup S_Z.",
  empirical_predictions="Catalogue k equals k_BB - delta for every entry, and agrees with a direct symplectic rank.",
  cheapest_decisive_test="Compute all three quantities independently for all 368 catalogue codes.",
  development_data="368 catalogue codes", held_out_data="7 Bravyi BB instances",
  evidence_for="368/368 exact agreement. delta distribution {0:213; 2:115; 4:34; 6:2; 8:2; 10:2}.",
  evidence_against="none", status="CONFIRMED",
  artifact_paths="src/qec_research/codes/pbb_theory.py; results/processed/exp006_parent_domination.json",
  next_action="Use as the basis of Theorem 1."),
 dict(id="H003", research_track="C",
  precise_hypothesis="Six catalogue PBB [[144,12,12]] codes marked EXACT possess a nontrivial weight-6 logical operator, so their distance claim is false.",
  motivation="An early audit found weight-6 vectors in ker(H_X) apparently outside the Z-stabilizer span.",
  mathematical_predictions="If true then d <= 6 for those codes, contradicting the reported d = 12.",
  empirical_predictions="CP-SAT returns a weight-6 witness that passes an independent is_logical check.",
  cheapest_decisive_test="Re-test each witness with StabilizerCode.is_logical: symplectic product against every row plus row-space membership.",
  development_data="6 catalogue codes", held_out_data="none", evidence_for="none",
  evidence_against="Every witness FAILED the independent check - each lies in the stabilizer group. Root cause: a pure-Z candidate must be paired against the x-part of a logical basis element, not against a Z-type representative.",
  status="FALSIFIED",
  artifact_paths="notes/failed_routes.md FR-001; src/qec_research/distance/sectors.py",
  next_action="Claim retracted. Corrected solver written and controlled on [[72,12,6]] -> 6 and [[108,8,10]] -> 10."),
 dict(id="H004", research_track="C",
  precise_hypothesis="For delta = 0 the nontrivial pure-Z logicals of the PBB code are identical to those of its parent BB code, hence d(PBB) <= d_Z(BB(A,B)).",
  motivation="Follows from H002 plus the fact that the pure-Z centraliser is ker(H_X) regardless of C and D.",
  mathematical_predictions="Pure-Z centralisers coincide (Prop 1); pure-Z stabilizers coincide when delta = 0 (Prop 2); so the minimum pure-Z logical weights are EQUAL, not merely ordered.",
  empirical_predictions="min pure-Z logical weight of PBB 12_6_0187 equals that of its parent exactly.",
  cheapest_decisive_test="Two certified CP-SAT minimisations, both required to reach status OPTIMAL.",
  development_data="PBB 12_6_0187 and its parent",
  held_out_data="controls [[72,12,6]] and [[108,8,10]]",
  evidence_for="Both returned 12 with CP-SAT status OPTIMAL (1086 s and 272 s). Controls reproduced the known exact distances 6 and 10.",
  evidence_against="none", status="CONFIRMED",
  artifact_paths="proofs/pbb_structure.md Theorem 1",
  next_action="Extend to Theorem 2 (CSS shadow) covering delta > 0."),
 dict(id="H005", research_track="C",
  precise_hypothesis="Every bivariate-bicycle code satisfies d_X = d_Z.",
  motivation="Needed to turn Theorem 1's bound on d_Z into a bound on the parent's full distance.",
  mathematical_predictions="The qubit permutation sigma mapping (L,g) to (R,-g) carries rowspace(H_X) onto rowspace(H_Z): swapping halves gives [B A] and lattice inversion transposes each circulant.",
  empirical_predictions="Row-space equality after applying sigma, for every BB code tested.",
  cheapest_decisive_test="Rank identities rank([sigma(H_X) ; H_Z]) = rank(H_Z) and its mirror.",
  development_data="7 Bravyi instances", held_out_data="368 catalogue parents",
  evidence_for="375/375 codes pass both rank tests. Independently exact CP-SAT gave d_X = d_Z for [[72,12,6]] [[90,8,10]] [[108,8,10]] and [[144,12,12]].",
  evidence_against="none", status="CONFIRMED",
  artifact_paths="proofs/pbb_structure.md Prop 4",
  next_action="Combine with H004 into Corollary 1 (parent domination)."),
 dict(id="H006", research_track="A",
  precise_hypothesis="A conflict-free (edge-coloured) gate schedule suffices for a one-ancilla syndrome-extraction circuit.",
  motivation="Standard scheduling intuition: never use a qubit twice in a layer.",
  mathematical_predictions="Konig colouring attains Delta = max(max check weight, max qubit degree) layers, which would be optimal.",
  empirical_predictions="A noiseless circuit built from a Konig colouring fires zero detectors.",
  cheapest_decisive_test="Build the Steane circuit from a Konig colouring and sample 200 noiseless shots.",
  development_data="Steane [[7,1,3]]", held_out_data="[[5,1,3]] and Shor [[9,1,3]]",
  evidence_for="none",
  evidence_against="1402 detector firings in 200 noiseless shots; 5 pairs violate the anticommuting-overlap parity rule.",
  status="FALSIFIED", artifact_paths="notes/failed_routes.md FR-002",
  next_action="Replaced by Lemma C1 and a parity-enforcing CP-SAT scheduler."),
 dict(id="H007", research_track="A",
  precise_hypothesis="Correct one-ancilla scheduling requires the count of j in J(a,b) with a acting before b to be even for every pair, where J(a,b) is the anticommuting shared support.",
  motivation="C-P and C-Q on a shared target differ by a CZ between the two ancillas when P and Q anticommute, and CZ squared is the identity.",
  mathematical_predictions="The criterion is exactly necessary and sufficient for the circuit to measure the intended generators.",
  empirical_predictions="A scheduler enforcing it gives zero noiseless detector firings on CSS and non-CSS codes alike.",
  cheapest_decisive_test="CP-SAT with reified order variables and an XOR parity constraint; test on Steane, Shor and the non-CSS [[5,1,3]].",
  development_data="Steane; [[5,1,3]]; Shor",
  held_out_data="BB [[72,12,6]] through [[288,12,18]]; PBB [[144,12,12]]",
  evidence_for="Zero detector firings on all toy codes and on every BB and PBB circuit built, 4000 noiseless shots each.",
  evidence_against="none", status="CONFIRMED",
  artifact_paths="proofs/pbb_structure.md Lemma C1; src/qec_research/circuits/scheduling.py",
  next_action="Adopt as the definition of a valid schedule throughout."),
 dict(id="H008", research_track="A",
  precise_hypothesis="The parity criterion alone forces syndrome depth 7 for weight-6 BB codes within the translation-invariant (orbit) schedule class.",
  motivation="Bravyi et al. report a depth-7 cycle while the combinatorial lower bound is only 6; if the parity rule closes that gap in the searched class, depth 7 is derived rather than assumed there.",
  mathematical_predictions="TI-model T = 6 infeasible and T = 7 achievable for every weight-6 BB code.",
  empirical_predictions="CP-SAT on the orbit model reports INFEASIBLE at T=6 and returns a verified schedule at T=7.",
  cheapest_decisive_test="Run the orbit-model scheduler at T=6 then T=7 for four BB instances.",
  development_data="[[72,12,6]]; [[108,8,10]]",
  held_out_data="[[144,12,12]]; [[288,12,18]]",
  evidence_for="All four: TI T=6 proven INFEASIBLE (not a timeout); T=7 valid with zero noiseless detector firings. Scope caught by the 2026-08-13 audit: the exp004 run is orbit-restricted, so no unrestricted depth-6 exclusion is ours; ASC (arXiv:2603.21499) certifies that exclusion externally.",
  evidence_against="none", status="CONFIRMED",
  artifact_paths="proofs/pbb_structure.md Theorem C3; results/raw/exp004_circuit_cost_survey.json",
  next_action="Report as a TI-class derivation of the published depth-7 cycle; cite ASC for the unrestricted exclusion."),
 dict(id="H009", research_track="A",
  precise_hypothesis="83 of 318 catalogue PBB codes admit no one-ancilla syndrome-extraction schedule at any depth.",
  motivation="CP-SAT returned INFEASIBLE at every probed depth for those codes.",
  mathematical_predictions="If the obstruction is intrinsic, the unrestricted model is also infeasible.",
  empirical_predictions="The general model reproduces INFEASIBLE on the smallest such codes.",
  cheapest_decisive_test="Re-run the three smallest orbit-unschedulable codes on the unrestricted model.",
  development_data="3 smallest orbit-unschedulable codes", held_out_data="none",
  evidence_for="none",
  evidence_against="All three ARE schedulable in the general model with status OPTIMAL: phase2_26 at depth 7; phase2_21 and phase2_44 at depth 9. Each verified with zero detector firings in 4000 noiseless shots.",
  status="FALSIFIED",
  artifact_paths="notes/failed_routes.md FR-005; results/processed/exp009_general_schedule_nogo.json",
  next_action="Restated as Proposition C5 scoped to translation-invariant schedules; depth separation re-derived from Proposition C2."),
 dict(id="H010", research_track="A",
  precise_hypothesis=("Every catalogue PBB [[144,12,12]] needs strictly more syndrome "
                      "depth than the CSS Gross code, for ANY one-ancilla schedule "
                      "measuring ANY generating set of its stabiliser group: the "
                      "separation is basis-independent, a property of the code."),
  motivation=("Proposition C2 bounds depth below by the weight of the measured "
              "generators; EXP-023 removes the generating-set dependence by bounding "
              "the light span V_7 of the whole group."),
  mathematical_predictions=("For every member, rank(V_6)=rank(V_7)=66 < 132 = n-k, so every "
                            "generating set contains a generator of symplectic weight >= 8; "
                            "hence depth >= 8 > 7 = Gross for every one-ancilla schedule."),
  empirical_predictions=("CP-SAT theorem gate 'exists weight<=7 element with nonzero X-part' "
                         "returns INFEASIBLE for all 14; pure-Z subgroup rank is exactly 66; "
                         "Gross control has rank(V_5)=0 and weight-6 rows spanning 132."),
  cheapest_decisive_test=("Certified-complete enumeration of weight<=7 X-codewords plus exact "
                          "GF(2) coset-budget tests; independent monolithic encoding as cross-check."),
  development_data="14 PBB [[144,12,12]] codes", held_out_data="CSS Gross code (control)",
  evidence_for=("EXP-023: all 14 gates INFEASIBLE with completeness certificates; r_Zsub=66 by "
                "three independent rank routes; positive SAT controls pass; every witness and "
                "rank re-verified in physical 288-bit space; w* exactly 8 for 4 members, "
                "in [8,9] for three and [8,10] for seven (does not affect the depth bound). "
                "Gross w*=6 exact. The independent monolithic CP-SAT encoding was run on "
                "benchmark member 12_6_0193 only, where it agreed."),
  evidence_against="none",
  status="CONFIRMED",
  artifact_paths=("proofs/pbb_structure.md Theorems C4 + C6; results/processed/exp023_light_gensets.json; "
                  "results/raw/exp023_independent_verification.json; tests/test_exp023_light_gensets.py"),
  next_action=("Feed the basis-independent separation into the paper draft; remaining open piece "
               "is the multi-ancilla model (EXP-024), not the generating set.")),
dict(id="H011", research_track="A,D",
  precise_hypothesis="Under matched circuit-level noise, rounds, decoder and observables, the non-CSS PBB [[144,12,12]] has a strictly higher logical error rate than the CSS Gross [[144,12,12]].",
  motivation="They share n, k and d, so any difference isolates the circuit and decoding cost of mixed stabilizers.",
  mathematical_predictions="Extra gate count (1008 vs 864), extra depth (8 vs 7) and mixed-check hook errors all push the PBB logical error rate up.",
  empirical_predictions="Higher schedule-averaged LER for the PBB code, with the circuit held fixed across noise points.",
  cheapest_decisive_test="Exhaustively enumerate the minimum-depth schedules, sample uniformly from the complete set, pin and persist each slot map, and compare the LER distributions.",
  development_data="EXP-007 pilot at p=0.002 (schedule NOT controlled; downgraded)",
  held_out_data="none. EXP-016 was designed AFTER seeing the exploratory EXP-007 result and after the schedule confound was discovered, so it is a corrected confirmatory measurement, not held-out validation. Genuine held-out checks would be a second p value and a second PBB code, neither of which was run under schedule control.",
  evidence_for="Minimum-depth TRANSLATION-INVARIANT schedules exhaustively enumerated (8496 CSS at depth 7; 9968 PBB at depth 8; CP-SAT status OPTIMAL). Schedule-averaged LER 8.90e-3 (CSS) vs 3.16e-2 (PBB) = 3.55x. Welch t = 3.83, p = 0.014. Mann-Whitney one-sided U = 25 of a maximum 25, p = 0.0040: every sampled PBB schedule is worse than every sampled CSS schedule.",
  evidence_against="The stricter pointwise test FAILS: with Bonferroni-corrected simultaneous intervals at familywise 95 percent the smallest PBB lower bound 1.30e-2 overlaps the largest CSS upper bound 1.91e-2. Also, schedule choice alone moves LER by 2.8x within each code, comparable to the between-code effect, so neither code was schedule-optimised.",
  status="CONFIRMED",
  artifact_paths="results/raw/exp016_schedule_controlled.json; notes/failed_routes.md FR-007",
  next_action="More schedules per code would be needed for pointwise familywise separation; between-schedule variance dominates, so more schedules beats more shots."),
 dict(id="H012", research_track="D",
  precise_hypothesis="Decoding the non-CSS PBB circuit is strictly more expensive than decoding the CSS Gross circuit at matched decoder settings.",
  motivation="Mixed checks densify the detector error model and destroy the CSS block structure that BP exploits.",
  mathematical_predictions="More DEM error mechanisms and higher BP+OSD latency percentiles for the PBB code.",
  empirical_predictions="Higher p50/p95/p99 per-shot decode latency at identical settings.",
  cheapest_decisive_test="Compare DEM size and latency percentiles under identical decoder settings on an idle machine.",
  development_data="quarantined shared-load pilot (FR-004)",
  held_out_data="none. EXP-013 is a corrected re-measurement of the same quantity, not held-out validation. A genuine held-out check would be a second p value and a second code pair; neither was run.",
  evidence_for="EXP-013 on an idle machine (load 1.83 at start, 2.62 at end), OMP_NUM_THREADS=1, single worker, warm-up discarded, the two codes interleaved in alternating blocks. p50 118.2 ms (CSS) vs 1220.2 ms (PBB) = 10.3x; p95 579.2 vs 9614.2 = 16.6x; p99 2159.0 vs 13629.1 = 6.3x; throughput 4.89 vs 0.43 shots/s/core = 11.4x. The DEM size ratio is only 1.24x and the detection-event ratio 1.37x, so the slowdown is superlinear in both.",
  evidence_against="Load was low but not exactly zero; measured at a single p and a single round count. Our proposed CAUSE was refuted: EXP-017 instrumented BpOsdDecoder.converge and found BP converges on 0.0 percent of shots for BOTH codes, so differential BP convergence explains nothing and the mechanism of the slowdown remains unresolved. Also, because BP never converges at these settings the comparison is effectively OSD-0 on the two DEMs, which is matched but is not the artifact decoder configuration.",
  status="CONFIRMED", artifact_paths="results/raw/exp013_isolated_latency.json",
  next_action="Optional: repeat at a second p to confirm the ratio is not p-specific."),
 dict(id="H013", research_track="B",
  precise_hypothesis=("On the held-out lattices (ell,m) = (12,12) and (15,12), every non-CSS PBB code "
                      "built from a canonical BB parent satisfies Prop. 3 (k_PBB = k_BB - delta with "
                      "delta >= 0), Prop. 4 on its parent (rank H_X = rank H_Z), and the Theorem-2 "
                      "shadow identity k(Q') = k(Q); consequently every delta = 0 code on these "
                      "lattices is weakly dominated in [[n,k,d]] by its parent CSS BB code."),
  motivation=("Every proposition in this repository was derived and verified on the 368-code published "
              "catalogue. These two lattices appear in the brief as underexplored for PBB and were "
              "never used to build or tune the theory, so they are genuine held-out data (GATE 10). "
              "Track B asks whether PBB codes there can exceed the CSS rate-distance envelope."),
  mathematical_predictions=("Zero violations of Prop. 3, Prop. 4 and the shadow-k identity over every "
                            "enumerated perturbation. delta >= 0 always. For delta = 0 codes Corollary 1 "
                            "then settles Track B with no distance computation at all. NOTE the theory "
                            "predicts nothing for delta > 0: Theorem 2 is a pure-Z ceiling and does not "
                            "bound d_X of the shadow, and two catalogue codes already escape it."),
  empirical_predictions=("Exhaustive k enumeration returns a dimension envelope only. Among all (C,D) of "
                         "weight <= 4 in the exact commutation null space, structural checks pass "
                         "100 percent. Max PBB check weight strictly exceeds the parent's |A|+|B| = 6 "
                         "whenever (C,D) != 0, so the circuit cost result carries over."),
  cheapest_decisive_test=("Enumerate the commutation condition as an exact GF(2) null space, enumerate "
                          "ALL its members of weight <= 4 by CP-SAT all-solution search (validated "
                          "against brute force on the known (6,6) lattice: 12 = 12), and run the "
                          "structural predicates on every one. A single violation falsifies."),
  development_data="the 368-code published catalogue on (6,6), (9,6), (12,6), (30,6), (18,6)",
  held_out_data="(12,12) and (15,12); no code on either lattice informed any proposition",
  evidence_for=("EXP-018 + EXP-019 on both held-out lattices. Structural: "
                 "896 non-CSS PBB codes on (12,12) and 1200 on (15,12) (weight<=2 of "
                 "6 representative parents each), 0 violations of Prop. 3, Prop. 4, "
                 "the Theorem-2 shadow identity, or P2 (delta=0 => parent's Z-logicals). "
                 "Distance: 400+400 children of 10+10 QEC-relevant parents (8<=k<=24), "
                 "weight<=2 perturbations enumerated COMPLETELY on every parent "
                 "(CP-SAT all-solution search; the two parents reporting 0 members were "
                 "confirmed by independent brute force), every one weakly dominated by "
                 "its parent CSS BB code with a certified witness. Parent distances exact: "
                 "d=6 at k=24, d=8 at k=16 on (15,12); d in {4,6} on (12,12)."),
  evidence_against="none observed",
  status="CONFIRMED",
  artifact_paths="results/processed/exp018_heldout_12x12.json; results/processed/exp018_heldout_15x12.json",
  next_action=("If structural checks pass, Track B is settled for delta = 0 by Corollary 1 and the "
               "residual question is exactly the delta > 0 case, which needs certified distances "
               "(stage 4) and is the same gap left open on the catalogue.")),
 dict(id="H014", date="2026-08-13", research_track="A",
  precise_hypothesis=("For unrestricted one-ancilla measurement of a fixed commuting "
                      "stabilizer generating set, the exact Lemma-C1 incidence-colouring "
                      "CSP predicts the minimum depth, and that minimum is always either "
                      "w_max or w_max+1."),
  motivation=("EXP-031's exact CSP matched every certified minimum but its "
              "translation-invariant PBB tier exposed that schedule-class emptiness can "
              "masquerade as a larger depth penalty. EXP-033 removes that restriction."),
  mathematical_predictions=("Whenever the depth-w_max CSP is satisfiable, depth w_max is "
                            "exact; otherwise depth w_max+1 is feasible. The claim is only "
                            "for unrestricted edge variables, not orbit-constrained schedules."),
  empirical_predictions=("The three TI-refuting weight-9 PBB instances admit verified "
                         "unrestricted schedules at depth 9, and the full unrestricted test "
                         "domain contains no two-value-law counterexample."),
  cheapest_decisive_test=("Run unrestricted CP-SAT at w_max and w_max+1 on the three "
                          "certified TI refutations; a certified INFEASIBLE result at both "
                          "depths falsifies the law, while a valid depth-w_max witness settles "
                          "each instance exactly by the check-weight lower bound."),
  development_data="150 random commuting stabilizer sets in EXP-031",
  held_out_data="PBB 12_6_0188, 12_6_0190, and catalogue row 41 selected only because they refuted the TI-scoped prediction",
  evidence_for=("EXP-031: exact predictor on all 208 decided instances, including 150 "
                "unrestricted random sets. EXP-033: all three adversarial PBB instances "
                "return valid unrestricted depth-9 schedules (OPTIMAL); w_max=9 proves "
                "each minimum exactly."),
  evidence_against=("The same three instances have no translation-invariant schedule "
                    "through depth 13, so the universal statement is false without the "
                    "unrestricted-class qualifier. The weak closed-form anticommuting-overlap "
                    "proxy has 157 false positives and is not the exact CSP."),
  status="PROVISIONALLY SUPPORTED",
  artifact_paths=("results/processed/exp031_depth_criterion.json; "
                  "results/processed/exp033_unrestricted_three.json; "
                  "proofs/pbb_structure.md Proposition C7"),
  next_action=("Prove the unrestricted w_max+1 sufficiency construction or search larger "
               "adversarial commuting sets; three PBB witnesses and 150 random instances "
               "are evidence, not a theorem.")),
]

if __name__ == "__main__":
    out = Path(__file__).resolve().parents[1] / "notes" / "hypotheses.csv"
    for r in R:
        r.setdefault("date", "2026-08-11")
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, quoting=csv.QUOTE_ALL)
        w.writeheader()
        for r in R:
            w.writerow({c: r.get(c, "") for c in COLS})
    import csv as _c
    rows = list(_c.DictReader(out.open()))
    assert len(rows) == len(R), (len(rows), len(R))
    for row in rows:
        assert len(row) == len(COLS) and None not in row, row["id"]
        assert row["status"] in {"CONFIRMED", "FALSIFIED", "PROVISIONALLY SUPPORTED",
                             "PREREGISTERED", "DOWNGRADED",
                             "CONFIRMED (scoped to published generators)",
                                 "UNRESOLVED"}, (row["id"], row["status"])
    print(f"wrote {out} : {len(rows)} hypotheses, {len(COLS)} columns, all well-formed")
    from collections import Counter
    print(Counter(r["status"] for r in rows))
