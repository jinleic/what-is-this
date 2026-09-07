#include "union_closed_interval_common.hpp"

int main() {
    std::cout << std::setprecision(17);
    const bool ok = verify_constants() && scalar_checks();
    std::cout << (ok ? "PROPOSITION A.4 CERTIFICATE PASSED\n" : "PROPOSITION A.4 CERTIFICATE FAILED\n");
    return ok ? 0 : 1;
}
