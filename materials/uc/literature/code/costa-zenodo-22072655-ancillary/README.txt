BEYOND LIU'S 0.382709 THRESHOLD FOR THE UNION-CLOSED SETS CONJECTURE
====================================================================

Main results
------------

Let

  c_L = 1 - x_* H(x_*) / H(x_*^2),

where x_* is the unique root in (0,1) of

  x^4 - 2 x^3 + 3 x^2 - 1 = 0.

The paper proves

  c_UC >= c_L > 0.382709087918735

and then the strict improvement

  c_UC >= c_L + 2e-5 > 0.382729087918735.

The first result combines the independent construction with Liu's Example 5.
The second combines that two-construction estimate with Liu's Example 4.
Proposition A.3 supplies the Example 4 product lower bound. Propositions A.4--A.6
verify the stronger pointwise inequality at mean

  m_cert = m_* - 2.05e-5.

The theorem is stated at mean m_* - 2e-5, leaving the factor

  (m_* - 2e-5) / (m_* - 2.05e-5) > 1,

as required by Liu's criterion.

Document structure
------------------

The main paper contains the mathematical argument and six finite computational
propositions. Propositions A.1 and A.3--A.6 are interval certificates;
Proposition A.2 is an exact rational Bernstein-coefficient calculation.
Operational details, ancillary filenames, build commands, subdivision
statistics, and reproduction data are kept in the computational supplement
and this README. The main paper and the computational supplement are
distributed separately; this archive contains only the ancillary verification
files listed below.

Sections 1--6 of the computational supplement correspond, in order, to
Propositions A.1--A.6. Each proposition has its own verification source and output file. The C++ drivers for A.3--A.6 share the implementation in
union_closed_interval_common.hpp; this header is common code, not a separate
certificate. Each of these drivers independently verifies the common constant
intervals it uses before running its proposition-specific checks. The
ancillary/arithmetic preamble and reproduction record are unnumbered.

Files
-----

  certify_prop_A1_mpfr.cpp
      C++ certificate driver for Proposition A.1.

  certify_prop_A1_mpfr_output.txt
      Complete output of the Proposition A.1 certificate.

  verify_prop_A2_bernstein.py
      Standard-library exact rational checker for Proposition A.2.

  verify_prop_A2_bernstein_output.txt
      Complete output of the Proposition A.2 checker.

  certify_prop_A3_mpfr.cpp
      C++ certificate driver for Proposition A.3.

  certify_prop_A3_mpfr_output.txt
      Complete output of the Proposition A.3 certificate.

  certify_prop_A4_mpfr.cpp
      C++ certificate driver for Proposition A.4.

  certify_prop_A4_mpfr_output.txt
      Complete output of the Proposition A.4 certificate.

  certify_prop_A5_mpfr.cpp
      C++ certificate driver for Proposition A.5.

  certify_prop_A5_mpfr_output.txt
      Complete output of the Proposition A.5 certificate.

  certify_prop_A6_mpfr.cpp
      C++ certificate driver for Proposition A.6.

  certify_prop_A6_mpfr_output.txt
      Complete output of the Proposition A.6 certificate.

  union_closed_interval_common.hpp
      Shared interval-arithmetic code used by the A.3--A.6 drivers.

  README.txt
      Package map and reproduction guide.

Compilation and verification
----------------------------

Proposition-specific interval certificates (with MPFR/GMP development
libraries installed):

  for A in A1 A3 A4 A5 A6; do
    g++ -O2 -std=c++17 -frounding-math -fno-fast-math -ffp-contract=off \
        certify_prop_${A}_mpfr.cpp -lmpfr -lgmp -o certify_${A}
    ./certify_${A}
  done

Exact Bernstein checker for Proposition A.2:

  python3 verify_prop_A2_bernstein.py

Successful output ends with the proposition-specific messages

  PROPOSITION A.1 CERTIFICATE PASSED
  PROPOSITION A.2 CERTIFICATE PASSED
  PROPOSITION A.3 CERTIFICATE PASSED
  PROPOSITION A.4 CERTIFICATE PASSED
  PROPOSITION A.5 CERTIFICATE PASSED
  PROPOSITION A.6 CERTIFICATE PASSED

The five C++ drivers use outward binary64 interval arithmetic and MPFR at
192-bit precision for transcendental endpoint evaluations. Decimal domain endpoints
from the paper are placed between adjacent binary64 values before subdivision,
and printed lower/upper bounds are moved one additional binary64 step outward
before decimal conversion. The A.1 driver also contains a separate 256-bit
directed MPFR check certifying

  c_L > 0.382709087918735.

The A.2 checker uses only fractions.Fraction and performs the entire polynomial
and power-to-Bernstein conversion over the rational numbers.

Reproducibility data from the supplied interval runs
-----------------------------------------------------

MPFR version: 4.2.2

Proposition A.1 (two-construction pointwise inequality):
  local Hessian determinant lower bound: 0.050114365218749776
  compact terminal cells: 238 + 363 + 245 + 722 + 671
  unresolved cells: 0

Proposition A.3 (Example 4 product lower bound):
  terminal cells: 1399
  maximum depth: 17
  unresolved cells: 0
  minimum certified lower bound: 2.3311935604120031e-08

Proposition A.4 (scalar checks):
  c_epsilon lower bound: 0.53084624205541608
  small-t lower bound: 0.035919173514281852
  boundary lower bound: 6.4360829338080744e-06
  large-t lower bound: 0.024386055203189365

Propositions A.5--A.6 (three-construction pointwise inequality):
  local Hessian grid: 40 x 20
  local stationary-point value lower bound: 6.7220206068085965e-07
  former equality point lower bound: 2.9270160742389521e-05
  A.6 repeated boundary scalar lower bound: 6.4360829338080744e-06
  complementary terminal cells: 4557
  maximum depth: 12
  unresolved cells: 0
  minimum certified lower bound: 5.3640462249763984e-09

Proposition A.2 (exact Bernstein check):
  bidegree: (3,3)
  number of coefficients: 16
  all coefficients strictly positive
