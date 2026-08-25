"""EXACT random-insertion identity behind projection contraction.

Let P be a distribution on n binary coordinates, fix coordinate j, and let
P^- be its projection after deleting j.  For a permutation sigma of the other
n-1 coordinates, write

  u_k = F(law(X_{sigma_k} | earlier sigma coordinates)) under P^-,
  r_k = F(law(X_{sigma_k} | X_j, earlier sigma coordinates)) under P,
  v_r = F(law(X_j | first r coordinates of sigma)).

Insert j into one of the n slots of sigma uniformly.  Coordinate sigma_k sees
j after it in n-k slots and before it in k slots.  Therefore

  A_n(P) = E_sigma [ (sum_k ((n-k)u_k+k r_k) + sum_r v_r) / n ].

Since A_{n-1}(P^-)=E_sigma sum_k u_k, the now-refuted 2/3 contraction candidate
was exactly

  E_sigma [ (n/3)sum_k u_k + sum_k k(r_k-u_k) + sum_r v_r ] >= 0.

The identity remains the exact accounting equation for any repaired potential:
negative refinement increments, insertion-coordinate laws, and the projected
reserve cannot be analyzed independently.  No convexity or union-closure
assumption enters the identity itself.

The executable control evaluates the complete-[5] contraction-ratio extremizer
P={(000):2/3,(011):1/6,(110):1/6}, deleting its middle (OR) coordinate.

Run: ./.venv/bin/python uc/shapley_insertion.py
"""

from fractions import Fraction
from itertools import permutations

from shapley_entropy import functional
from shapley_n5_projections import conditional_law, direct_value


def delete_coordinate(states, probabilities, coordinate, dimension):
    kept = tuple(j for j in range(dimension) if j != coordinate)
    mass = {}
    for state, probability in zip(states, probabilities):
        image = sum(((state >> old) & 1) << new
                    for new, old in enumerate(kept))
        mass[image] = mass.get(image, Fraction()) + probability
    projected_states = tuple(sorted(mass))
    return projected_states, tuple(mass[state] for state in projected_states)


def insertion_decomposition(states, probabilities, dimension, inserted):
    others = tuple(j for j in range(dimension) if j != inserted)
    child_states, child_probabilities = delete_coordinate(
        states, probabilities, inserted, dimension)
    child_index = {old: new for new, old in enumerate(others)}

    base = 0.0
    refinement = 0.0
    insertion = 0.0
    order_count = 0
    for order in permutations(others):
        order_count += 1
        unrefined = []
        for position, coordinate in enumerate(order, 1):
            predecessors = order[:position - 1]
            child_predecessors = tuple(child_index[j] for j in predecessors)
            child_coordinate = child_index[coordinate]
            u_value = functional(conditional_law(
                child_states, child_probabilities,
                child_coordinate, child_predecessors))[0]
            r_value = functional(conditional_law(
                states, probabilities, coordinate,
                predecessors + (inserted,)))[0]
            unrefined.append(u_value)
            refinement += position * (r_value - u_value)

        base += dimension * sum(unrefined) / 3.0
        for position in range(dimension):
            insertion += functional(conditional_law(
                states, probabilities, inserted, order[:position]))[0]

    return {
        "base": base / order_count,
        "refinement": refinement / order_count,
        "insertion": insertion / order_count,
        "parent": direct_value(states, probabilities, dimension),
        "child": direct_value(
            child_states, child_probabilities, dimension - 1),
        "child_distribution": (child_states, child_probabilities),
    }


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    states = (0, 3, 6)
    probabilities = (Fraction(2, 3), Fraction(1, 6), Fraction(1, 6))
    terms = insertion_decomposition(states, probabilities, 3, 1)
    left = terms["base"] + terms["refinement"] + terms["insertion"]
    right = 3 * (terms["parent"] - (2.0 / 3.0) * terms["child"])
    assert abs(left - right) < 2e-12
    assert right > 0

    print("PROVED [permutation slot counting identity]")
    print("NUMERICAL CONTROL [float64 HiGHS local functionals]")
    print("  parent distribution =", (states, probabilities))
    print("  projected child =", terms["child_distribution"])
    print("  A(parent) = %.15f" % terms["parent"])
    print("  A(child)  = %.15f" % terms["child"])
    print("  (n/3) A(child) term = %+.15f" % terms["base"])
    print("  weighted refinement term = %+.15f" % terms["refinement"])
    print("  inserted-coordinate term = %+.15f" % terms["insertion"])
    print("  sum = n[A(parent)-(2/3)A(child)] = %+.15f" % left)
    print()
    print("FRONTIER")
    print("  Refinement plus insertion is negative; the child reserve is essential.")
    print("  The aggregate is not universally nonnegative (certified elsewhere).")
    print("  Reuse the identity for nonlinear-potential restriction inequalities.")
