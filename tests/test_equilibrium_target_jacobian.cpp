#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <cumes/config/json_reader.hpp>
#include <cumes/config/validated_problem.hpp>
#include <cumes/solver/equilibrium_linearization.hpp>
#include <cumes/solver/equilibrium_solver.hpp>
#include <meow/cumes/boundary_parameterization.hpp>
#include <meow/cumes/landreman_workflow.hpp>
#include <meow/cumes/quasisymmetry_target.hpp>
#include <meow/cumes/quasisymmetry_target_jvp.hpp>

namespace {

using cumes_meow_example::BoundaryFamily;
using cumes_meow_example::LandremanCase;
using cumes_meow_example::LandremanSelection;
using cumes_meow_example::LandremanWorkflow;
using cumes_meow_example::StellaratorSymmetricBoundaryParameterization;

void set_stage_resolution(cumes::ProblemSpec& problem,
                          const cumes_meow_example::LandremanStage& stage) {
    problem.mpol = stage.mpol;
    problem.ntor = stage.ntor;
    for (auto& request : problem.stages) {
        request.max_iterations =
            std::max(request.max_iterations,
                     static_cast<std::size_t>(stage.minimum_niter));
        request.tolerance = std::max(request.tolerance, stage.tolerance_floor);
    }
    problem.raxis_c.resize(static_cast<std::size_t>(problem.ntor + 1), 0.0);
    problem.zaxis_s.resize(static_cast<std::size_t>(problem.ntor + 1), 0.0);
}

cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec target_spec(int nfp) {
    cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec spec;
    spec.helicity_m = 1;
    spec.helicity_n = -nfp;
    for (int index = 0; index <= 10; ++index) {
        spec.normalized_toroidal_flux_surfaces.push_back(index / 10.0);
        spec.surface_weights.push_back(1.0);
    }
    return spec;
}

double vector_norm(const std::vector<double>& values) {
    double squared = 0.0;
    for (double value : values) squared += value * value;
    return std::sqrt(squared);
}

}  // namespace

int main() {
    using meow::test::check;

    const LandremanSelection selection{LandremanCase::QH,
                                       LandremanWorkflow::CONSTRUCTION};
    const std::string input_path =
        std::string(MEOW_SOURCE_DIR) + "/examples/landreman/qh_analytic.json";
    cumes::ParsedProblem parsed = cumes::read_problem_spec(input_path, {});
    check(parsed.report.ok(), "analytic QH input parses");
    cumes::ProblemSpec baseline = std::move(parsed.spec);
    cumes_meow_example::refresh_axis_predictor_from_boundary_centerline(
        baseline);
    set_stage_resolution(
        baseline, cumes_meow_example::landreman_stages(selection).front());
    const auto validated = cumes::validate(baseline, {});
    check(validated.has_value(), "mode-1 analytic QH input validates");
    if (!validated.has_value()) { return meow::test::summary(); }

    cumes::SolveRequest request;
    request.radial_transfer = cumes::RadialTransferPolicy::CATMULL_ROM;
    cumes::EquilibriumSolver solver;
    const cumes::SolveOutcome primal = solver.solve(validated.value(), request);
    check(primal.converged && primal.has_complete_equilibrium(),
          "mode-1 analytic QH equilibrium converges");
    if (!primal.converged || !primal.has_complete_equilibrium()) {
        return meow::test::summary();
    }

    StellaratorSymmetricBoundaryParameterization boundary(1);
    const meow::Vector center = boundary.values(baseline);
    const auto& degrees = boundary.degrees_of_freedom();
    const auto selected =
        std::find_if(degrees.begin(), degrees.end(), [](const auto& degree) {
            return degree.family == BoundaryFamily::RBC && degree.m == 0 &&
                   degree.n == 1;
        });
    check(selected != degrees.end(), "QH boundary direction exists");
    const std::size_t column =
        static_cast<std::size_t>(selected - degrees.begin());

    cumes::EquilibriumLinearization linearization(validated.value(),
                                                  primal.equilibrium);
    cumes::TangentLinearOptions tangent_options;
    tangent_options.max_iterations = 1000;
    tangent_options.restart = 300;
    tangent_options.relative_tolerance = 1.0e-6;
    tangent_options.absolute_tolerance = 1.0e-11;
    const auto spectral = linearization.solve_boundary_tangent(
        boundary.tangent(validated.value(), column), tangent_options);
    check(spectral.converged, "QH equilibrium tangent solve converges");
    if (!spectral.converged) { return meow::test::summary(); }
    const cumes::EquilibriumTangent tangent = linearization.materialize_tangent(
        spectral.state_tangent, primal.equilibrium, primal.profiles);
    const auto spec = target_spec(primal.report.input_params.nfp);
    const std::vector<double> analytic =
        cumes_meow_example::calculate_qh_target_jvp(
            primal.equilibrium, primal.profiles, tangent,
            primal.report.input_params, spec);

    constexpr double epsilon = 1.0e-3;
    const auto solve_perturbed = [&](double sign) {
        meow::Vector values = center;
        values[static_cast<Eigen::Index>(column)] += sign * epsilon;
        cumes::ProblemSpec problem = boundary.apply(baseline, values);
        problem.stages = {problem.stages.back()};
        cumes_meow_example::track_axis_predictor_from_accepted_boundary(
            problem, baseline);
        const auto perturbed = cumes::validate(std::move(problem), {});
        if (!perturbed.has_value()) {
            throw std::runtime_error("perturbed QH boundary did not validate");
        }
        cumes::SolveRequest hot_request = request;
        hot_request.restart = std::cref(primal.equilibrium);
        return solver.solve(perturbed.value(), hot_request);
    };
    const cumes::SolveOutcome plus = solve_perturbed(1.0);
    const cumes::SolveOutcome minus = solve_perturbed(-1.0);
    check(plus.converged && minus.converged,
          "centered QH target oracle equilibria converge");
    if (!plus.converged || !minus.converged) { return meow::test::summary(); }
    const auto plus_target = cumes_meow_example::calculate_qh_target(
        plus.equilibrium, plus.profiles, plus.report.input_params, spec,
        cumes_meow_example::landreman_target_aspect(LandremanCase::QH));
    const auto minus_target = cumes_meow_example::calculate_qh_target(
        minus.equilibrium, minus.profiles, minus.report.input_params, spec,
        cumes_meow_example::landreman_target_aspect(LandremanCase::QH));
    check(analytic.size() == plus_target.residuals.size() &&
              analytic.size() == minus_target.residuals.size(),
          "analytic and nonlinear QH target tangents have equal extent");

    std::vector<double> finite_difference(analytic.size());
    std::vector<double> difference(analytic.size());
    for (std::size_t index = 0; index < analytic.size(); ++index) {
        finite_difference[index] =
            (plus_target.residuals[index] - minus_target.residuals[index]) /
            (2.0 * epsilon);
        difference[index] = analytic[index] - finite_difference[index];
    }
    const double relative_error =
        vector_norm(difference) / vector_norm(finite_difference);
    const double objective_fd =
        (plus_target.value - minus_target.value) / (2.0 * epsilon);
    const auto primal_target = cumes_meow_example::calculate_qh_target(
        primal.equilibrium, primal.profiles, primal.report.input_params, spec,
        cumes_meow_example::landreman_target_aspect(LandremanCase::QH));
    double objective_analytic = 0.0;
    for (std::size_t index = 0; index < analytic.size(); ++index) {
        objective_analytic +=
            2.0 * primal_target.residuals[index] * analytic[index];
    }
    const double objective_relative_error =
        std::abs(objective_analytic - objective_fd) /
        std::max(std::abs(objective_fd), 1.0e-12);
    std::cout << "QH target tangent column=" << boundary.name(column)
              << " GMRES_iterations=" << spectral.iterations
              << " residual_relative_error=" << relative_error
              << " objective_relative_error=" << objective_relative_error
              << " objective_analytic=" << objective_analytic
              << " objective_fd=" << objective_fd << '\n';
    check(std::isfinite(relative_error) &&
              std::isfinite(objective_relative_error),
          "QH equilibrium target tangent comparison is finite");
    check(relative_error < 6.0e-2,
          "QH target residual tangent agrees with the nonlinear oracle");
    check(objective_relative_error < 1.0e-2,
          "QH objective tangent agrees with the nonlinear oracle");

    return meow::test::summary();
}
