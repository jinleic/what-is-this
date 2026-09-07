#include "union_closed_interval_common.hpp"

int main() {
    std::cout << std::setprecision(17);
    const bool ok = verify_constants() && certify_combined();
    std::cout << (ok ? "PROPOSITION A.6 CERTIFICATE PASSED\n" : "PROPOSITION A.6 CERTIFICATE FAILED\n");
    return ok ? 0 : 1;
}
