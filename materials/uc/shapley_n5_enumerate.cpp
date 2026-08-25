// Exact constraint enumerator for the n=5 direct-Shapley experiment.
//
// Emits one 8-hex-digit uint32 family mask per line. Bit s says that subset
// s of [5] belongs to the family. The accepted families are exactly the
// nontrivial families satisfying
//
//   max_i frequency(i) <= 0.38261,
//   2 * average set size >= log_2 |F|.
//
// The only feasible family sizes and integer incidence thresholds are listed
// below. They follow directly by flooring 0.38261*m and ceiling
// m*log_2(m)/2; the runtime assertions independently check the table.
//
// Build/run:
//   c++ -O3 -std=c++17 uc/shapley_n5_enumerate.cpp -o /tmp/shapley_n5_enum
//   /tmp/shapley_n5_enum > /tmp/shapley_n5_families.txt

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

constexpr std::array<Spec, 7> SPECS{{
    {3, 1, 3},
    {4, 1, 4},
    {6, 2, 8},
    {7, 2, 10},
    {8, 3, 12},
    {9, 3, 15},
    {11, 4, 20},
}};

class Enumerator {
  public:
    explicit Enumerator(Spec spec) : spec_(spec) {}

    std::uint64_t run() {
        choose(0, 0, 0, 0U);
        return count_;
    }

  private:
    void choose(int next_subset, int chosen, int incidence,
                std::uint32_t family_mask) {
        const int remaining = spec_.family_size - chosen;
        if (remaining == 0) {
            if (incidence >= spec_.incidence_need) {
                std::cout << std::hex << std::setw(8) << std::setfill('0')
                          << family_mask << '\n';
                ++count_;
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
};

}  // namespace

int main() {
    std::uint64_t total = 0;
    for (const Spec spec : SPECS) {
        const int expected_cap = static_cast<int>(
            std::floor(38261.0 * spec.family_size / 100000.0));
        const int expected_need = static_cast<int>(
            std::ceil(spec.family_size * std::log2(spec.family_size) / 2.0
                      - 1e-12));
        if (expected_cap != spec.coordinate_cap
                || expected_need != spec.incidence_need) {
            std::cerr << "threshold table mismatch\n";
            return 2;
        }
        Enumerator enumerator(spec);
        const std::uint64_t count = enumerator.run();
        total += count;
        std::cerr << "m=" << spec.family_size << " count=" << count << '\n';
    }
    std::cerr << "total=" << total << '\n';
    return std::cout.good() ? 0 : 3;
}
