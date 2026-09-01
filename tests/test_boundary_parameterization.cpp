#include "test_support.hpp"

#include <stdexcept>
#include <string>

#include <meow/cumes/boundary_parameterization.hpp>

int main() {
    using cumes_meow_example::BoundaryFamily;
    using cumes_meow_example::StellaratorSymmetricBoundaryParameterization;
    using meow::test::check;

    cumes::ProblemSpec problem;
    problem.rbc = {{0, 0, 1.0}, {0, 1, 0.1}, {1, -1, 0.2}, {5, 5, 0.9}};
    problem.zbs = {{0, 0, 0.0}, {1, 0, -0.3}};

    StellaratorSymmetricBoundaryParameterization modes4(4);
    check(modes4.size() == 80,
          "max_mode=4 has the same 80 active DOFs as SIMSOPT");
    check(modes4.name(0) == "rbc(0,1)", "RBC m=0 positive-n modes come first");
    check(modes4.name(4) == "rbc(1,-4)",
          "RBC m>=1 modes include signed toroidal indices");
    check(modes4.name(40) == "zbs(0,1)", "ZBS modes follow all RBC modes");

    const meow::Vector initial = modes4.values(problem);
    check(initial[0] == 0.1, "existing coefficients are extracted");
    check(initial[4] == 0.0, "missing coefficients are initialized to zero");
    check(initial[7] == 0.2, "signed-n coefficients use cuMES (m,n) order");

    meow::Vector changed = initial;
    changed[0] = 0.125;
    changed[4] = -0.05;
    const cumes::ProblemSpec applied = modes4.apply(problem, changed);
    const meow::Vector round_trip = modes4.values(applied);
    check((round_trip - changed).norm() == 0.0,
          "absolute boundary vector round trips through ProblemSpec");
    check(applied.rbc.front().value == 1.0, "fixed rbc(0,0) is not changed");
    check(applied.rbc[3].value == 0.9,
          "modes outside the active range remain fixed");

    cumes::ProblemSpec predictor_problem;
    predictor_problem.ntor = 3;
    predictor_problem.raxis_c = {9.0};
    predictor_problem.zaxis_s = {8.0};
    predictor_problem.rbc = {{0, 0, 1.0},  {0, 1, 0.12}, {0, 3, -0.02},
                             {0, -1, 7.0}, {1, 2, 6.0},  {0, 4, 5.0}};
    predictor_problem.zbs = {{0, 0, 0.0}, {0, 2, -0.03}, {2, 1, 4.0}};
    cumes_meow_example::refresh_axis_predictor_from_boundary_centerline(
        predictor_problem);
    check(predictor_problem.raxis_c ==
              std::vector<double>({1.0, 0.12, 0.0, -0.02}),
          "R-axis predictor follows nonnegative m=0 boundary harmonics");
    check(predictor_problem.zaxis_s ==
              std::vector<double>({0.0, 0.0, -0.03, 0.0}),
          "Z-axis predictor zero-fills missing centerline harmonics");
    check(predictor_problem.has_raxis_c && predictor_problem.has_zaxis_s,
          "refreshed axis predictor is serialized explicitly");

    StellaratorSymmetricBoundaryParameterization modes5(5);
    check(modes5.size() == 120,
          "max_mode=5 has the same 120 active DOFs as SIMSOPT");

    bool rejected = false;
    try {
        [[maybe_unused]] StellaratorSymmetricBoundaryParameterization invalid(
            0);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "zero max_mode is rejected");

    return meow::test::summary();
}
