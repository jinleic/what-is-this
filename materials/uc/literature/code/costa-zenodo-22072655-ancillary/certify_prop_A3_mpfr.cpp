#include "union_closed_interval_common.hpp"

static bool verify_xstar_bracket() {
    auto poly=[](const I&x){return x*x*x*x-I(2.0)*x*x*x+I(3.0)*x*x-I(1.0);};
    I pl=poly(decimal_point(0.69078759392498));
    I pu=poly(decimal_point(0.69078759392500));
    I dp=I(4.0)*XSTAR*XSTAR*XSTAR-I(6.0)*XSTAR*XSTAR+I(6.0)*XSTAR;
    const bool left=pl.upper()<0;
    const bool right=pu.lower()>0;
    const bool deriv=dp.lower()>0;
    std::cout << "A.3 root-bracket checks: left=" << left
              << ", right=" << right << ", derivative=" << deriv << "\n";
    return left && right && deriv;
}

int main() {
    std::cout << std::setprecision(17);
    const bool ok = verify_xstar_bracket() && certify_minorant();
    std::cout << (ok ? "PROPOSITION A.3 CERTIFICATE PASSED\n" : "PROPOSITION A.3 CERTIFICATE FAILED\n");
    return ok ? 0 : 1;
}
