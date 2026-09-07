#include "union_closed_interval_common.hpp"

int main() {
    std::cout << std::setprecision(17);
    const bool ok = verify_constants() && certify_local_hessian_and_stationary_point();
    std::cout << (ok ? "PROPOSITION A.5 CERTIFICATE PASSED\n" : "PROPOSITION A.5 CERTIFICATE FAILED\n");
    return ok ? 0 : 1;
}
