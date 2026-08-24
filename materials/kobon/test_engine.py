import unittest
from itertools import combinations

import engine


class FaceMultiplicityCountsTest(unittest.TestCase):
    def test_line_through_vertex_crossing_interior_counts(self):
        # Triangle from y=0, y=x, y=-x+2.  y=x/2 passes through (0,0)
        # and then enters the open triangle before crossing its opposite side.
        ms = [0, 1, -1, engine.Fr(1, 2)]
        bs = [0, 0, 2, 0]
        crossings, bounded_faces = engine.face_multiplicity_counts(
            4, ms, bs, [(0, 1, 2)])
        self.assertEqual((crossings, bounded_faces), (1, 2))

    def test_line_through_vertex_tangent_does_not_count(self):
        # y=2x passes through the same vertex but stays outside the triangle.
        ms = [0, 1, -1, 2]
        bs = [0, 0, 2, 0]
        crossings, bounded_faces = engine.face_multiplicity_counts(
            4, ms, bs, [(0, 1, 2)])
        self.assertEqual((crossings, bounded_faces), (0, 2))

    def test_embedded_certificates_satisfy_face_multiplicity(self):
        for n, lines, selected in (
            (10, engine.N10_LOWER_BOUND_LINES,
             engine.N10_LOWER_BOUND_TRIANGLES),
            (12, engine.N12_LOWER_BOUND_LINES,
             engine.N12_LOWER_BOUND_TRIANGLES),
        ):
            ms, bs = zip(*lines)
            ok, reason = engine.verify_selection(
                n, list(ms), list(bs), selected, minimum=len(selected))
            self.assertTrue(ok, reason)
            crossing_count, bounded_faces = engine.face_multiplicity_counts(
                n, list(ms), list(bs), selected)
            self.assertLessEqual(len(selected) + crossing_count, bounded_faces)


class CrossingFaceBoundCnfTest(unittest.TestCase):
    def test_face_bound_counts_supplied_crossing_literals(self):
        n, target = 4, 2
        cnf = engine.CNF()
        pool = engine.IDPool()
        crossing_lits = [pool.id(("TEST_CROSS", i)) for i in range(2)]
        engine.add_face_bound(
            cnf, pool, n, target, crossing_lits=crossing_lits)

        no_degeneracy = [
            pool.id(("P", i, j)) for i, j in combinations(range(n), 2)
        ] + [
            -pool.id(("C",) + t) for t in combinations(range(n), 3)
        ]
        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(
                solver.solve(assumptions=no_degeneracy
                             + [crossing_lits[0], -crossing_lits[1]]))
            self.assertFalse(
                solver.solve(assumptions=no_degeneracy + crossing_lits))

    def test_crossing_indicator_is_exact_for_selected_triangle(self):
        n = 4
        cnf = engine.CNF()
        pool = engine.IDPool()
        for t in combinations(range(n), 3):
            pool.id(("S",) + t)
            m = next(r for r in range(n) if r not in t)
            for p, q in combinations(t, 2):
                key = tuple(sorted((p, q))) + (m,)
                pool.id(("SA",) + key)
                pool.id(("SB",) + key)

        crossing_lits = engine.add_triangle_crossing_indicators(cnf, pool, n)
        t, m = (0, 1, 2), 3
        selected = pool.id(("S",) + t)
        crossing = pool.id(("TC",) + t + (m,))
        above = [pool.id(("SA",) + pq + (m,))
                 for pq in combinations(t, 2)]
        below = [pool.id(("SB",) + pq + (m,))
                 for pq in combinations(t, 2)]
        self.assertIn(crossing, crossing_lits)

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(
                assumptions=[selected, above[0], below[1], -crossing]))
            self.assertFalse(solver.solve(assumptions=[crossing, -selected]))
            self.assertFalse(solver.solve(
                assumptions=[crossing] + [-v for v in above]))
            self.assertFalse(solver.solve(
                assumptions=[crossing] + [-v for v in below]))
            self.assertTrue(solver.solve(
                assumptions=[selected, above[0], below[1], crossing]))
            self.assertTrue(solver.solve(assumptions=[-selected, -crossing]))

    def test_crossed_triangle_indicator_is_forced_by_straddling_line(self):
        n = 4
        cnf = engine.CNF()
        pool = engine.IDPool()
        for t in combinations(range(n), 3):
            pool.id(("S",) + t)
            m = next(r for r in range(n) if r not in t)
            for p, q in combinations(t, 2):
                key = tuple(sorted((p, q))) + (m,)
                pool.id(("SA",) + key)
                pool.id(("SB",) + key)

        crossed_lits = engine.add_crossed_triangle_indicators(cnf, pool, n)
        t, m = (0, 1, 2), 3
        selected = pool.id(("S",) + t)
        crossed = pool.id(("TX",) + t)
        above = pool.id(("SA", 0, 1, m))
        below = pool.id(("SB", 0, 2, m))
        self.assertIn(crossed, crossed_lits)
        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(
                assumptions=[selected, above, below, -crossed]))
            self.assertFalse(solver.solve(assumptions=[crossed, -selected]))
            self.assertTrue(solver.solve(
                assumptions=[selected, above, below, crossed]))
            self.assertTrue(solver.solve(assumptions=[-selected, -crossed]))

    def test_compact_no_concurrency_flag_is_forced_by_mixed_sides(self):
        n = 4
        cnf = engine.CNF()
        pool = engine.IDPool()
        for t in combinations(range(n), 3):
            pool.id(("S",) + t)
            m = next(r for r in range(n) if r not in t)
            for p, q in combinations(t, 2):
                pool.id(("SA",) + tuple(sorted((p, q))) + (m,))

        crossed_lits = (
            engine.add_no_concurrency_crossed_triangle_indicators(
                cnf, pool, n))
        t, m = (0, 1, 2), 3
        selected = pool.id(("S",) + t)
        crossed = pool.id(("TX",) + t)
        sides = [pool.id(("SA",) + pq + (m,))
                 for pq in combinations(t, 2)]
        self.assertIn(crossed, crossed_lits)
        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(assumptions=[
                selected, sides[0], -sides[1], -crossed]))
            self.assertTrue(solver.solve(assumptions=[
                selected, sides[0], -sides[1], crossed]))
            self.assertFalse(solver.solve(assumptions=[crossed, -selected]))


    def test_default_face_bound_omits_experimental_crossing_penalty(self):
        cnf, pool = engine.build_model(4, 2)
        engine.add_face_bound(cnf, pool, 4, 2)
        self.assertNotIn(("TX", 0, 1, 2), pool.obj2id)
        self.assertNotIn(("TC", 0, 1, 2, 3), pool.obj2id)




class SubarrangementFaceBoundCnfTest(unittest.TestCase):
    def test_caps_faces_supported_on_every_eleven_line_subset(self):
        n = 12
        cnf = engine.CNF()
        pool = engine.IDPool()
        for triple in combinations(range(n), 3):
            pool.id(("S",) + triple)

        engine.add_subarrangement_face_bounds(cnf, pool, n, {11: 32})
        inside = [
            pool.id(("S",) + triple)
            for triple in combinations(range(11), 3)
        ]

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(solver.solve(assumptions=inside[:32]))
            self.assertFalse(solver.solve(assumptions=inside[:33]))

    def test_exact_target_form_enforces_the_deletion_lower_bound(self):
        n = 5
        cnf = engine.CNF()
        pool = engine.IDPool()
        triples = list(combinations(range(n), 3))
        selected = [pool.id(("S",) + triple) for triple in triples]

        engine.add_subarrangement_face_bounds(
            cnf, pool, n, {4: 2}, exact_target=3)
        inside = [
            pool.id(("S",) + triple)
            for triple in combinations(range(4), 3)
        ]
        outside = next(
            pool.id(("S",) + triple) for triple in triples if 4 in triple)

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            all_false = [-literal for literal in selected]
            three_inside = all_false.copy()
            for literal in inside[:3]:
                three_inside[selected.index(literal)] = literal
            self.assertFalse(solver.solve(assumptions=three_inside))

            two_inside_one_outside = all_false.copy()
            for literal in inside[:2] + [outside]:
                two_inside_one_outside[selected.index(literal)] = literal
            self.assertTrue(solver.solve(
                assumptions=two_inside_one_outside))


class FaceVertexSectorBoundCnfTest(unittest.TestCase):
    def test_simple_vertex_allows_four_faces_but_not_five(self):
        n = 7
        cnf = engine.CNF()
        pool = engine.IDPool()
        pair_faces = [
            pool.id(("S", 0, 1, third)) for third in range(2, n)
        ]

        engine.add_selected_face_sector_bounds(cnf, pool, n)
        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(solver.solve(assumptions=pair_faces[:4]))
            self.assertFalse(solver.solve(assumptions=pair_faces))

    def test_multipoint_allows_at_most_two_faces_on_a_line_pair(self):
        n = 7
        cnf = engine.CNF()
        pool = engine.IDPool()
        pair_faces = [
            pool.id(("S", 0, 1, third)) for third in range(3, 6)
        ]
        concurrent = pool.id(("C", 0, 1, 2))

        engine.add_selected_face_sector_bounds(cnf, pool, n)
        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(solver.solve(
                assumptions=[concurrent] + pair_faces[:2]))
            self.assertFalse(solver.solve(
                assumptions=[concurrent] + pair_faces))

    def test_nonadjacent_multipoint_lines_support_no_face(self):
        n = 7
        cnf = engine.CNF()
        pool = engine.IDPool()
        inside = pool.id(("C", 1, 2, 4))
        outside = pool.id(("C", 0, 1, 4))
        selected = pool.id(("S", 1, 4, 6))

        engine.add_selected_face_sector_bounds(cnf, pool, n)
        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(
                assumptions=[inside, outside, selected]))


class FaceSelectionCnfTest(unittest.TestCase):
    def test_selected_face_rejects_an_outside_straddling_line(self):
        cnf, pool = engine.build_model(4, 1, rules=("FACE",))
        selected = pool.id(("S", 0, 1, 2))
        above = pool.id(("SA", 0, 1, 3))
        below = pool.id(("SB", 0, 2, 3))

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(
                assumptions=[selected, above, below]))

    def test_selected_face_rejects_an_intersection_inside_a_side(self):
        cnf, pool = engine.build_model(4, 1, rules=("FACE",))
        selected = pool.id(("S", 0, 1, 2))
        between = pool.id(("B", 0, 3, 1, 2))

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(
                assumptions=[selected, between]))

class NoConcurrencyCapacityCnfTest(unittest.TestCase):
    def test_global_target_implies_per_line_lower_bounds(self):
        n, target = 5, 5
        cnf = engine.CNF()
        pool = engine.IDPool()
        triples = list(combinations(range(n), 3))
        for t in triples:
            pool.id(("S",) + t)

        engine.add_no_concurrency_capacities(
            cnf, pool, n, [0] * n, target=target)
        on_line_zero = [pool.id(("S",) + t) for t in triples if 0 in t]
        cycle = {(0, 1, 2), (1, 2, 3), (2, 3, 4),
                 (0, 3, 4), (0, 1, 4)}
        cycle_model = [pool.id(("S",) + t)
                       if t in cycle else -pool.id(("S",) + t)
                       for t in triples]

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(
                assumptions=on_line_zero[:2]
                + [-v for v in on_line_zero[2:]]))
            self.assertTrue(solver.solve(assumptions=cycle_model))

    def test_target_is_required_for_coupled_lower_bounds(self):
        n = 5
        cnf = engine.CNF()
        pool = engine.IDPool()
        triples = list(combinations(range(n), 3))
        for t in triples:
            pool.id(("S",) + t)

        engine.add_no_concurrency_capacities(cnf, pool, n, [0] * n)
        on_line_zero = [pool.id(("S",) + t) for t in triples if 0 in t]
        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(solver.solve(
                assumptions=on_line_zero[:2]
                + [-v for v in on_line_zero[2:]]))


class NoConcurrencyCrossingBudgetCnfTest(unittest.TestCase):
    def test_budget_caps_supplied_crossing_incidences(self):
        cnf = engine.CNF()
        pool = engine.IDPool()
        crossings = [pool.id(("TEST_CROSS", i)) for i in range(4)]
        engine.add_no_concurrency_crossing_budget(
            cnf, pool, 4, 2, [0] * 4, crossing_lits=crossings)

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(solver.solve(
                assumptions=crossings[:2] + [-v for v in crossings[2:]]))
            self.assertFalse(solver.solve(assumptions=crossings[:3]))

    def test_per_line_budget_couples_sides_and_crossings(self):
        cnf, pool = engine.build_model(4, 1)
        engine.add_no_concurrency_crossing_budget(
            cnf, pool, 4, 1, [0] * 4, per_line=True)
        sides_on_zero = [
            pool.id(("S",) + t)
            for t in combinations(range(4), 3) if 0 in t
        ]
        crossing_zero = pool.id(("TC", 1, 2, 3, 0))

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertFalse(solver.solve(
                assumptions=sides_on_zero[:2] + [crossing_zero]))


class MultipointCapacityCnfTest(unittest.TestCase):
    def test_triple_point_adds_only_one_side_slot_on_each_incident_line(self):
        n = 5
        cnf, pool = engine.build_model(n, 1, rules=())
        try:
            engine.add_face_bound(cnf, pool, n, 1, per_line=True)
        except TypeError:
            self.fail("add_face_bound does not expose multipoint capacities")

        selected = [
            pool.id(("S",) + t)
            for t in combinations(range(n), 3)
            if 0 in t and t != (0, 1, 2)
        ]
        all_cross = [
            pool.id(("P", i, j))
            for i, j in combinations(range(n), 2)
        ]
        concurrency = [
            pool.id(("C",) + t) if t == (0, 1, 2)
            else -pool.id(("C",) + t)
            for t in combinations(range(n), 3)
        ]

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(solver.solve(
                assumptions=selected[:4] + [-selected[4]]
                + all_cross + concurrency))
            self.assertFalse(solver.solve(
                assumptions=selected + all_cross + concurrency))

    @staticmethod
    def _tc_family(pool, n, triples):
        return [
            pool.id(("TC",) + t + (r,))
            for t in triples for r in range(n) if r not in t
        ]

    def test_per_line_rejects_non_exact_crossing_literals(self):
        cnf = engine.CNF()
        pool = engine.IDPool()
        n = 4
        triples = list(combinations(range(n), 3))
        for t in triples:
            pool.id(("S",) + t)
        fake = self._tc_family(pool, n, triples)
        fake[0] = pool.id(("TX",) + triples[0])

        with self.assertRaises(ValueError):
            engine.add_face_bound(
                cnf, pool, n, 1, crossing_lits=fake, per_line=True)

    def test_triple_point_couples_sides_and_crossings_per_line(self):
        n = 5
        cnf, pool = engine.build_model(n, 1, rules=())
        triples = list(combinations(range(n), 3))
        crossing_lits = self._tc_family(pool, n, triples)
        try:
            engine.add_face_bound(
                cnf, pool, n, 1,
                crossing_lits=crossing_lits, per_line=True)
        except ValueError as exc:
            self.fail(f"per_line rejected the exact TC family: {exc}")

        sides = [
            pool.id(("S",) + t)
            for t in ((0, 1, 3), (0, 1, 4), (0, 2, 4))
        ]
        incidences = [
            pool.id(("TC",) + t + (0,))
            for t in ((1, 2, 3), (1, 2, 4))
        ]
        all_cross = [
            pool.id(("P", i, j))
            for i, j in combinations(range(n), 2)
        ]
        concurrency = [
            pool.id(("C",) + t) if t == (0, 1, 2)
            else -pool.id(("C",) + t)
            for t in triples
        ]

        with engine.Solver(name="cadical195",
                           bootstrap_with=cnf.clauses) as solver:
            self.assertTrue(solver.solve(
                assumptions=sides + [incidences[0], -incidences[1]]
                + all_cross + concurrency))
            self.assertFalse(solver.solve(
                assumptions=sides + incidences + all_cross + concurrency))


if __name__ == "__main__":
    unittest.main()
