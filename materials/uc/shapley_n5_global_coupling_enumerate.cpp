// Exact cap-2/5 n=5 family enumerator for the global-coupling experiment.
//
// Emits one canonical uint32 family mask per S_5 coordinate-permutation orbit.
// Accepted families are exactly the nontrivial uniform simple families with
//
//   max_i frequency(i) <= 2/5,
//   2 * average set size >= log_2 |F|.
//
// Family enumeration, cap/incidence filtering, and orbit canonicalization are
// exact integer operations.  The stderr counts include all labeled families
// and canonical representatives separately.
//
// Build:
//   c++ -O3 -std=c++17 uc/shapley_n5_global_coupling_enumerate.cpp \
//       -o /tmp/shapley_n5_global_coupling_enum

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

struct Spec {
    int family_size;
    int coordinate_cap;
    int incidence_need;
};

constexpr std::array<Spec, 11> SPECS{{
    {3, 1, 3},
    {4, 1, 4},
    {5, 2, 6},
    {6, 2, 8},
    {7, 2, 10},
    {8, 3, 12},
    {9, 3, 15},
    {10, 4, 17},
    {11, 4, 20},
    {13, 5, 25},
    {15, 6, 30},
}};

using Permutation = std::array<int, 5>;

std::vector<Permutation> permutations() {
    std::vector<Permutation> result;
    Permutation permutation{{0, 1, 2, 3, 4}};
    do {
        result.push_back(permutation);
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

const std::vector<Permutation> PERMUTATIONS = permutations();

int permute_subset(int subset, const Permutation& permutation) {
    int image = 0;
    for (int old_coordinate = 0; old_coordinate < 5; ++old_coordinate) {
        if ((subset >> old_coordinate) & 1) {
            image |= 1 << permutation[old_coordinate];
        }
    }
    return image;
}

std::uint32_t permute_family(std::uint32_t family_mask,
                             const Permutation& permutation) {
    std::uint32_t image = 0;
    for (int subset = 0; subset < 32; ++subset) {
        if ((family_mask >> subset) & 1U) {
            image |= std::uint32_t{1} << permute_subset(subset, permutation);
        }
    }
    return image;
}

bool is_canonical(std::uint32_t family_mask) {
    for (const Permutation& permutation : PERMUTATIONS) {
        if (permute_family(family_mask, permutation) < family_mask) {
            return false;
        }
    }
    return true;
}

class Enumerator {
  public:
    explicit Enumerator(Spec spec) : spec_(spec) {}

    void run() {
        choose(0, 0, 0, 0U);
    }

    std::uint64_t count() const { return count_; }
    std::uint64_t canonical_count() const { return canonical_count_; }

  private:
    void choose(int next_subset, int chosen, int incidence,
                std::uint32_t family_mask) {
        const int remaining = spec_.family_size - chosen;
        if (remaining == 0) {
            if (incidence < spec_.incidence_need) {
                return;
            }
            ++count_;
            if (is_canonical(family_mask)) {
                std::cout << std::hex << std::setw(8) << std::setfill('0')
                          << family_mask << '\n';
                ++canonical_count_;
            }
            return;
        }
        if (32 - next_subset < remaining) {
            return;
        }

        int remaining_coordinate_capacity = 0;
        for (int count : coordinate_counts_) {
            remaining_coordinate_capacity += spec_.coordinate_cap - count;
        }
        const int incidence_upper = incidence
            + std::min(5 * remaining, remaining_coordinate_capacity);
        if (incidence_upper < spec_.incidence_need) {
            return;
        }

        const int last_start = 32 - remaining;
        for (int subset = next_subset; subset <= last_start; ++subset) {
            bool admissible = true;
            int subset_size = 0;
            for (int coordinate = 0; coordinate < 5; ++coordinate) {
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
            for (int coordinate = 0; coordinate < 5; ++coordinate) {
                coordinate_counts_[coordinate] += (subset >> coordinate) & 1;
            }
            choose(subset + 1, chosen + 1, incidence + subset_size,
                   family_mask | (std::uint32_t{1} << subset));
            for (int coordinate = 0; coordinate < 5; ++coordinate) {
                coordinate_counts_[coordinate] -= (subset >> coordinate) & 1;
            }
        }
    }

    Spec spec_;
    std::array<int, 5> coordinate_counts_{{0, 0, 0, 0, 0}};
    std::uint64_t count_ = 0;
    std::uint64_t canonical_count_ = 0;
};

}  // namespace

int main() {
    if (PERMUTATIONS.size() != 120) {
        std::cerr << "permutation table mismatch\n";
        return 2;
    }
    std::uint64_t total = 0;
    std::uint64_t canonical_total = 0;
    for (const Spec spec : SPECS) {
        const int expected_cap = 2 * spec.family_size / 5;
        const int expected_need = static_cast<int>(
            std::ceil(spec.family_size * std::log2(spec.family_size) / 2.0
                      - 1e-12));
        if (expected_cap != spec.coordinate_cap
                || expected_need != spec.incidence_need) {
            std::cerr << "threshold table mismatch\n";
            return 3;
        }
        Enumerator enumerator(spec);
        enumerator.run();
        total += enumerator.count();
        canonical_total += enumerator.canonical_count();
        std::cerr << "m=" << spec.family_size
                  << " count=" << enumerator.count()
                  << " canonical=" << enumerator.canonical_count() << '\n';
    }
    std::cerr << "total=" << total
              << " canonical_total=" << canonical_total << '\n';
    return std::cout.good() ? 0 : 4;
}
