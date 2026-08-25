// Exact small-size cap-2/5 n=6 family enumerator.
//
// Emits one canonical uint64 family mask per S_6 coordinate-permutation orbit
// for every nontrivial uniform simple family of size 3 through 7 satisfying
//   max_i count(i) <= floor(2m/5),
//   total incidence >= ceil(m log_2(m)/2).
//
// Enumeration and cap/incidence filters use integer arithmetic.  Once one
// labeled member of an orbit is reached, all of its coordinate permutations
// are inserted into a seen set; this canonicalizes once per orbit rather than
// once per labeled family.  The final seen-set size must equal the independently
// accumulated labeled count.
//
// Build:
//   c++ -O3 -std=c++20 uc/shapley_n6_small_enumerate.cpp \
//       -o /tmp/shapley_n6_small_enum

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_set>
#include <vector>

namespace {

struct Spec {
    int family_size;
    int coordinate_cap;
    int incidence_need;
};

constexpr std::array<Spec, 5> SPECS{{
    {3, 1, 3},
    {4, 1, 4},
    {5, 2, 6},
    {6, 2, 8},
    {7, 2, 10},
}};

using Permutation = std::array<int, 6>;
using PermutedSubsets = std::array<std::uint8_t, 64>;

std::vector<Permutation> make_permutations() {
    std::vector<Permutation> result;
    Permutation permutation{{0, 1, 2, 3, 4, 5}};
    do {
        result.push_back(permutation);
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

const std::vector<Permutation> PERMUTATIONS = make_permutations();

std::vector<PermutedSubsets> make_permuted_subsets() {
    std::vector<PermutedSubsets> result;
    result.reserve(PERMUTATIONS.size());
    for (const Permutation& permutation : PERMUTATIONS) {
        PermutedSubsets images{};
        for (int subset = 0; subset < 64; ++subset) {
            int image = 0;
            for (int old_coordinate = 0; old_coordinate < 6;
                 ++old_coordinate) {
                if ((subset >> old_coordinate) & 1) {
                    image |= 1 << permutation[old_coordinate];
                }
            }
            images[subset] = static_cast<std::uint8_t>(image);
        }
        result.push_back(images);
    }
    return result;
}

const std::vector<PermutedSubsets> PERMUTED_SUBSETS =
    make_permuted_subsets();

class Enumerator {
  public:
    explicit Enumerator(Spec spec) : spec_(spec) {
        seen_.reserve(1 << 20);
    }

    bool run() {
        choose(0, 0, 0, 0);
        return seen_.size() == labeled_count_;
    }

    std::uint64_t labeled_count() const { return labeled_count_; }
    std::uint64_t orbit_count() const { return orbit_count_; }

  private:
    void emit_orbit(std::uint64_t family_mask) {
        if (seen_.contains(family_mask)) {
            return;
        }
        std::uint64_t canonical = family_mask;
        for (const PermutedSubsets& images : PERMUTED_SUBSETS) {
            std::uint64_t image = 0;
            for (int index = 0; index < spec_.family_size; ++index) {
                image |= std::uint64_t{1} << images[selected_[index]];
            }
            canonical = std::min(canonical, image);
            seen_.insert(image);
        }
        std::cout << std::hex << std::setw(16) << std::setfill('0')
                  << canonical << '\n';
        ++orbit_count_;
    }

    void choose(int next_subset, int chosen, int incidence,
                std::uint64_t family_mask) {
        const int remaining = spec_.family_size - chosen;
        if (remaining == 0) {
            if (incidence < spec_.incidence_need) {
                return;
            }
            ++labeled_count_;
            emit_orbit(family_mask);
            return;
        }
        if (64 - next_subset < remaining) {
            return;
        }

        int remaining_coordinate_capacity = 0;
        for (int count : coordinate_counts_) {
            remaining_coordinate_capacity += spec_.coordinate_cap - count;
        }
        const int incidence_upper = incidence
            + std::min(6 * remaining, remaining_coordinate_capacity);
        if (incidence_upper < spec_.incidence_need) {
            return;
        }

        const int last_start = 64 - remaining;
        for (int subset = next_subset; subset <= last_start; ++subset) {
            bool admissible = true;
            int subset_size = 0;
            for (int coordinate = 0; coordinate < 6; ++coordinate) {
                if ((subset >> coordinate) & 1) {
                    ++subset_size;
                    if (coordinate_counts_[coordinate]
                            == spec_.coordinate_cap) {
                        admissible = false;
                        break;
                    }
                }
            }
            if (!admissible) {
                continue;
            }
            for (int coordinate = 0; coordinate < 6; ++coordinate) {
                coordinate_counts_[coordinate] +=
                    (subset >> coordinate) & 1;
            }
            selected_[chosen] = subset;
            choose(subset + 1, chosen + 1, incidence + subset_size,
                   family_mask | (std::uint64_t{1} << subset));
            for (int coordinate = 0; coordinate < 6; ++coordinate) {
                coordinate_counts_[coordinate] -=
                    (subset >> coordinate) & 1;
            }
        }
    }

    Spec spec_;
    std::array<int, 6> coordinate_counts_{{0, 0, 0, 0, 0, 0}};
    std::array<int, 7> selected_{{0, 0, 0, 0, 0, 0, 0}};
    std::unordered_set<std::uint64_t> seen_;
    std::uint64_t labeled_count_ = 0;
    std::uint64_t orbit_count_ = 0;
};

}  // namespace

int main() {
    if (PERMUTATIONS.size() != 720
            || PERMUTED_SUBSETS.size() != 720) {
        std::cerr << "permutation table mismatch\n";
        return 2;
    }

    std::uint64_t labeled_total = 0;
    std::uint64_t orbit_total = 0;
    for (const Spec spec : SPECS) {
        const int expected_cap = 2 * spec.family_size / 5;
        const int expected_need = static_cast<int>(std::ceil(
            spec.family_size * std::log2(spec.family_size) / 2.0 - 1e-12));
        if (expected_cap != spec.coordinate_cap
                || expected_need != spec.incidence_need) {
            std::cerr << "threshold table mismatch\n";
            return 3;
        }
        Enumerator enumerator(spec);
        if (!enumerator.run()) {
            std::cerr << "orbit coverage mismatch at m="
                      << spec.family_size << '\n';
            return 4;
        }
        labeled_total += enumerator.labeled_count();
        orbit_total += enumerator.orbit_count();
        std::cerr << std::dec
                  << "m=" << spec.family_size
                  << " labeled=" << enumerator.labeled_count()
                  << " orbits=" << enumerator.orbit_count() << '\n';
    }
    std::cerr << "total_labeled=" << labeled_total
              << " total_orbits=" << orbit_total << '\n';
    return std::cout.good() ? 0 : 5;
}
