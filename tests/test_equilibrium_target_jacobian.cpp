#include "test_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/SVD>
#include <cumes/config/json_reader.hpp>
#include <cumes/config/validated_problem.hpp>
#include <cumes/solver/equilibrium_linearization.hpp>
#include <cumes/solver/equilibrium_solver.hpp>
#include <meow/cumes/boundary_parameterization.hpp>
#include <meow/cumes/landreman_workflow.hpp>
#include <meow/cumes/magnetic_gradient_target.hpp>
#include <meow/cumes/plasma_size_target_jvp.hpp>
#include <meow/cumes/quasisymmetry_target.hpp>
#include <meow/cumes/quasisymmetry_target_jvp.hpp>

namespace {

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

double max_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) result = std::max(result, std::abs(value));
    return result;
}

struct VectorComparison {
    double reference_norm = 0.0;
    double error_norm = 0.0;
    double relative_error = 0.0;
};

VectorComparison compare_vectors(const std::vector<double>& actual,
                                 const std::vector<double>& reference) {
    if (actual.size() != reference.size()) {
        throw std::runtime_error("comparison vector shape mismatch");
    }
    std::vector<double> error(actual.size());
    for (std::size_t index = 0; index < error.size(); ++index) {
        error[index] = actual[index] - reference[index];
    }
    VectorComparison result;
    result.reference_norm = vector_norm(reference);
    result.error_norm = vector_norm(error);
    result.relative_error =
        result.error_norm / std::max(result.reference_norm, 1.0e-14);
    return result;
}

void print_comparison(const std::string& column,
                      const std::string& quantity,
                      const std::vector<double>& actual,
                      const std::vector<double>& reference) {
    const VectorComparison comparison = compare_vectors(actual, reference);
    std::cout << "QH field tangent column=" << column
              << " quantity=" << quantity
              << " reference_norm=" << comparison.reference_norm
              << " error_norm=" << comparison.error_norm
              << " relative_error=" << comparison.relative_error << '\n';
}

std::vector<double> difference(const std::vector<double>& perturbed,
                               const std::vector<double>& primal,
                               double step) {
    if (perturbed.size() != primal.size()) {
        throw std::runtime_error("finite-difference vector shape mismatch");
    }
    std::vector<double> result(primal.size());
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] = (perturbed[index] - primal[index]) / step;
    }
    return result;
}

double relative_difference(const std::vector<double>& actual,
                           const std::vector<double>& reference) {
    std::vector<double> delta(actual.size());
    if (actual.size() != reference.size()) {
        throw std::runtime_error("relative-error vector shape mismatch");
    }
    for (std::size_t index = 0; index < delta.size(); ++index) {
        delta[index] = actual[index] - reference[index];
    }
    return vector_norm(delta) / std::max(vector_norm(reference), 1.0e-14);
}

cumes::EquilibriumTangent finite_difference_tangent(
    const cumes::SolveOutcome& primal,
    const cumes::SolveOutcome& perturbed,
    double step) {
    cumes::EquilibriumTangent result = cumes::EquilibriumTangent::zero_like(
        primal.equilibrium, primal.profiles);
    for (std::size_t family = 0; family < cumes::EquilibriumSnapshot::COUNT;
         ++family) {
        result.equilibrium.families[family] =
            difference(perturbed.equilibrium.families[family],
                       primal.equilibrium.families[family], step);
    }
    for (std::size_t field = 0;
         field < cumes::EquilibriumSnapshot::HALF_FIELD_COUNT; ++field) {
        result.equilibrium.half_fields[field] =
            difference(perturbed.equilibrium.half_fields[field],
                       primal.equilibrium.half_fields[field], step);
    }
    for (std::size_t field = 0;
         field < cumes::EquilibriumSnapshot::FULL_FIELD_COUNT; ++field) {
        result.equilibrium.full_fields[field] =
            difference(perturbed.equilibrium.full_fields[field],
                       primal.equilibrium.full_fields[field], step);
    }
    result.profiles.toroidal_flux_derivative =
        difference(perturbed.profiles.toroidal_flux_derivative,
                   primal.profiles.toroidal_flux_derivative, step);
    result.profiles.poloidal_flux_derivative =
        difference(perturbed.profiles.poloidal_flux_derivative,
                   primal.profiles.poloidal_flux_derivative, step);
    result.profiles.rotational_transform =
        difference(perturbed.profiles.rotational_transform,
                   primal.profiles.rotational_transform, step);
    result.profiles.poloidal_covariant_field =
        difference(perturbed.profiles.poloidal_covariant_field,
                   primal.profiles.poloidal_covariant_field, step);
    result.profiles.toroidal_covariant_field =
        difference(perturbed.profiles.toroidal_covariant_field,
                   primal.profiles.toroidal_covariant_field, step);
    return result;
}

std::vector<double> spectral_state(const cumes::EquilibriumTangent& tangent) {
    const std::size_t family_size = tangent.equilibrium.family_size();
    std::vector<double> result(cumes::EquilibriumSnapshot::COUNT * family_size);
    for (std::size_t family = 0; family < cumes::EquilibriumSnapshot::COUNT;
         ++family) {
        std::copy(tangent.equilibrium.families[family].begin(),
                  tangent.equilibrium.families[family].end(),
                  result.begin() + family * family_size);
    }
    return result;
}

cumes::EquilibriumSnapshot add_scaled_tangent(
    const cumes::EquilibriumSnapshot& primal,
    const cumes::EquilibriumTangent& tangent,
    double scale) {
    cumes::EquilibriumSnapshot result = primal;
    for (std::size_t field = 0;
         field < cumes::EquilibriumSnapshot::HALF_FIELD_COUNT; ++field) {
        for (std::size_t index = 0; index < result.half_fields[field].size();
             ++index) {
            result.half_fields[field][index] +=
                scale * tangent.equilibrium.half_fields[field][index];
        }
    }
    return result;
}

cumes::EquilibriumProfiles add_scaled_tangent(
    const cumes::EquilibriumProfiles& primal,
    const cumes::EquilibriumTangent& tangent,
    double scale) {
    cumes::EquilibriumProfiles result = primal;
    const auto add = [scale](std::vector<double>& values,
                             const std::vector<double>& direction) {
        for (std::size_t index = 0; index < values.size(); ++index) {
            values[index] += scale * direction[index];
        }
    };
    add(result.toroidal_flux_derivative,
        tangent.profiles.toroidal_flux_derivative);
    add(result.poloidal_flux_derivative,
        tangent.profiles.poloidal_flux_derivative);
    add(result.rotational_transform, tangent.profiles.rotational_transform);
    add(result.poloidal_covariant_field,
        tangent.profiles.poloidal_covariant_field);
    add(result.toroidal_covariant_field,
        tangent.profiles.toroidal_covariant_field);
    return result;
}

cumes_meow_example::MagneticGradientFields magnetic_observable_tangent(
    const cumes::SolveOutcome& primal,
    const cumes::EquilibriumTangent& tangent,
    int nfp) {
    constexpr double EPSILON = 1.0e-6;
    const auto plus_equilibrium =
        add_scaled_tangent(primal.equilibrium, tangent, EPSILON);
    const auto minus_equilibrium =
        add_scaled_tangent(primal.equilibrium, tangent, -EPSILON);
    const auto plus_profiles =
        add_scaled_tangent(primal.profiles, tangent, EPSILON);
    const auto minus_profiles =
        add_scaled_tangent(primal.profiles, tangent, -EPSILON);
    const auto plus = cumes_meow_example::calculate_magnetic_gradient_fields(
        plus_equilibrium, plus_profiles, nfp);
    const auto minus = cumes_meow_example::calculate_magnetic_gradient_fields(
        minus_equilibrium, minus_profiles, nfp);
    cumes_meow_example::MagneticGradientFields result;
    result.field_strength =
        difference(plus.field_strength, minus.field_strength, 2.0 * EPSILON);
    result.b_dot_grad_b =
        difference(plus.b_dot_grad_b, minus.b_dot_grad_b, 2.0 * EPSILON);
    result.b_cross_grad_s_dot_grad_b =
        difference(plus.b_cross_grad_s_dot_grad_b,
                   minus.b_cross_grad_s_dot_grad_b, 2.0 * EPSILON);
    result.b_cross_grad_toroidal_flux_dot_grad_b =
        difference(plus.b_cross_grad_toroidal_flux_dot_grad_b,
                   minus.b_cross_grad_toroidal_flux_dot_grad_b, 2.0 * EPSILON);
    result.b_cross_grad_psi_p_dot_grad_b =
        difference(plus.b_cross_grad_psi_p_dot_grad_b,
                   minus.b_cross_grad_psi_p_dot_grad_b, 2.0 * EPSILON);
    return result;
}

double objective_derivative(const std::vector<double>& residual,
                            const std::vector<double>& tangent) {
    if (residual.size() != tangent.size()) {
        throw std::runtime_error("objective derivative shape mismatch");
    }
    double result = 0.0;
    for (std::size_t index = 0; index < residual.size(); ++index) {
        result += 2.0 * residual[index] * tangent[index];
    }
    return result;
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
    cumes::EquilibriumLinearization linearization(validated.value(),
                                                  primal.equilibrium);
    cumes::TangentLinearOptions tangent_options;
    tangent_options.max_iterations = 1000;
    tangent_options.restart = 300;
    tangent_options.relative_tolerance = 5.0e-6;
    tangent_options.absolute_tolerance = 1.0e-11;
    const auto spec = target_spec(primal.report.input_params.nfp);
    const auto primal_target = cumes_meow_example::calculate_qh_target(
        primal.equilibrium, primal.profiles, primal.report.input_params, spec,
        cumes_meow_example::landreman_target_aspect(LandremanCase::QH));

    const auto check_zero_basis_lambda = [&](std::size_t component,
                                             const std::string& name) {
        std::vector<double> state(linearization.state_size(), 0.0);
        const std::size_t family_size = primal.equilibrium.family_size();
        for (int surface = 0; surface < primal.equilibrium.ns; ++surface) {
            state[component * family_size + static_cast<std::size_t>(surface)] =
                1.0 + static_cast<double>(surface) /
                          static_cast<double>(primal.equilibrium.ns);
        }
        const cumes::ResidualJvp residual = linearization.residual_jvp(state);
        const cumes::EquilibriumTangent lambda_tangent =
            linearization.materialize_tangent(state, primal.equilibrium,
                                              primal.profiles);
        double public_field_max = 0.0;
        for (const auto& field : lambda_tangent.equilibrium.half_fields) {
            public_field_max = std::max(public_field_max, max_abs(field));
        }
        double profile_max = 0.0;
        for (const auto* profile :
             {&lambda_tangent.profiles.toroidal_flux_derivative,
              &lambda_tangent.profiles.poloidal_flux_derivative,
              &lambda_tangent.profiles.rotational_transform,
              &lambda_tangent.profiles.poloidal_covariant_field,
              &lambda_tangent.profiles.toroidal_covariant_field}) {
            profile_max = std::max(profile_max, max_abs(*profile));
        }
        const std::vector<double> target_tangent =
            cumes_meow_example::calculate_qh_target_jvp(
                primal.equilibrium, primal.profiles, lambda_tangent,
                primal.report.input_params, spec);
        const double aspect_tangent =
            cumes_meow_example::calculate_plasma_size_jvp(
                primal.equilibrium, lambda_tangent, primal.report.input_params)
                .aspect_ratio;
        std::cout << "QH zero-basis lambda component=" << name
                  << " residual_jvp_max=" << max_abs(residual.tangent)
                  << " public_field_max=" << public_field_max
                  << " profile_max=" << profile_max
                  << " aspect_tangent=" << aspect_tangent
                  << " target_jvp_max=" << max_abs(target_tangent) << '\n';
        check(max_abs(residual.tangent) < 1.0e-13 &&
                  public_field_max < 1.0e-13 && profile_max < 1.0e-13 &&
                  std::abs(aspect_tangent) < 1.0e-13 &&
                  max_abs(target_tangent) < 1.0e-13,
              "zero-basis " + name +
                  " is absent from the residual and target bridge");
    };
    check_zero_basis_lambda(cumes::EquilibriumSnapshot::LMNSC, "LMNSC(0,0)");
    check_zero_basis_lambda(cumes::EquilibriumSnapshot::LMNCS, "LMNCS(0,0)");

    const auto finite_difference_policy =
        cumes_meow_example::landreman_finite_difference_policy(selection);
    double worst_target_chain_error = 0.0;
    double worst_analytic_error = 0.0;
    double worst_objective_error = 0.0;
    double worst_hybrid_policy_error = 0.0;
    std::vector<std::vector<double>> branch_differences;
    std::vector<std::string> branch_names;

    for (std::size_t column = 0; column < boundary.size(); ++column) {
        const cumes::BoundaryTangent boundary_tangent =
            boundary.tangent(validated.value(), column);
        const auto spectral = linearization.solve_boundary_tangent(
            boundary_tangent, tangent_options);
        std::cout << "QH equilibrium tangent column=" << boundary.name(column)
                  << " converged=" << spectral.converged
                  << " iterations=" << spectral.iterations
                  << " initial_residual=" << spectral.initial_residual
                  << " final_residual=" << spectral.final_residual
                  << " relative_residual="
                  << spectral.final_residual / spectral.initial_residual
                  << '\n';
        check(spectral.converged && std::isfinite(spectral.final_residual) &&
                  spectral.state_tangent.size() == linearization.state_size(),
              "QH equilibrium tangent solve reaches the production tolerance "
              "for " +
                  boundary.name(column));
        const cumes::EquilibriumTangent tangent =
            linearization.materialize_tangent(
                spectral.state_tangent, primal.equilibrium, primal.profiles);
        const std::vector<double> analytic =
            cumes_meow_example::calculate_qh_target_jvp(
                primal.equilibrium, primal.profiles, tangent,
                primal.report.input_params, spec);

        const double step =
            std::max(finite_difference_policy.relative_step *
                         std::abs(center[static_cast<Eigen::Index>(column)]),
                     1.0e-4);
        meow::Vector values = center;
        values[static_cast<Eigen::Index>(column)] += step;
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
        const cumes::SolveOutcome plus =
            solver.solve(perturbed.value(), hot_request);
        check(plus.converged,
              "forward QH target oracle equilibrium converges for " +
                  boundary.name(column));
        if (!plus.converged) { continue; }
        std::cout << "QH hot-restart column=" << boundary.name(column)
                  << " nonlinear_iterations=" << plus.total_iterations << '\n';
        const auto plus_target = cumes_meow_example::calculate_qh_target(
            plus.equilibrium, plus.profiles, plus.report.input_params, spec,
            cumes_meow_example::landreman_target_aspect(LandremanCase::QH));
        const std::vector<double> finite_difference =
            difference(plus_target.residuals, primal_target.residuals, step);
        const cumes::EquilibriumTangent nonlinear_tangent =
            finite_difference_tangent(primal, plus, step);
        const std::vector<double> nonlinear_state =
            spectral_state(nonlinear_tangent);
        std::vector<double> branch_difference(spectral.state_tangent.size());
        for (std::size_t index = 0; index < branch_difference.size(); ++index) {
            branch_difference[index] =
                spectral.state_tangent[index] - nonlinear_state[index];
        }
        const std::size_t family_size = primal.equilibrium.family_size();
        double lambda_energy = 0.0;
        double m1_mixed_rz_energy = 0.0;
        double other_rz_energy = 0.0;
        for (std::size_t family = 0; family < cumes::EquilibriumSnapshot::COUNT;
             ++family) {
            for (int mode = 0; mode < primal.equilibrium.mnmax; ++mode) {
                const int m = mode / (primal.report.input_params.ntor + 1);
                for (int surface = 0; surface < primal.equilibrium.ns;
                     ++surface) {
                    const std::size_t index =
                        family * family_size +
                        static_cast<std::size_t>(mode) * primal.equilibrium.ns +
                        surface;
                    const double energy =
                        branch_difference[index] * branch_difference[index];
                    if (family == cumes::EquilibriumSnapshot::LMNSC ||
                        family == cumes::EquilibriumSnapshot::LMNCS) {
                        lambda_energy += energy;
                    } else if (m == 1 &&
                               (family == cumes::EquilibriumSnapshot::RMNSS ||
                                family == cumes::EquilibriumSnapshot::ZMNCS)) {
                        m1_mixed_rz_energy += energy;
                    } else {
                        other_rz_energy += energy;
                    }
                }
            }
        }
        const double total_branch_energy =
            lambda_energy + m1_mixed_rz_energy + other_rz_energy;
        std::cout << "QH branch difference column=" << boundary.name(column)
                  << " norm=" << std::sqrt(total_branch_energy)
                  << " implicit_norm=" << vector_norm(spectral.state_tangent)
                  << " nonlinear_norm=" << vector_norm(nonlinear_state)
                  << " lambda_fraction="
                  << lambda_energy / std::max(total_branch_energy, 1.0e-30)
                  << " m1_mixed_rz_fraction="
                  << m1_mixed_rz_energy / std::max(total_branch_energy, 1.0e-30)
                  << " other_rz_fraction="
                  << other_rz_energy / std::max(total_branch_energy, 1.0e-30)
                  << '\n';
        branch_differences.push_back(std::move(branch_difference));
        branch_names.push_back(boundary.name(column));
        const cumes::ResidualJvp boundary_residual =
            linearization.boundary_residual_jvp(boundary_tangent);
        const cumes::ResidualJvp implicit_combined =
            linearization.residual_jvp(spectral.state_tangent);
        const cumes::ResidualJvp nonlinear_combined =
            linearization.residual_jvp(spectral_state(nonlinear_tangent));
        std::cout << "QH equilibrium response column=" << boundary.name(column)
                  << " boundary_residual_norm="
                  << vector_norm(boundary_residual.tangent)
                  << " implicit_combined_residual_norm="
                  << vector_norm(implicit_combined.tangent)
                  << " nonlinear_combined_residual_norm="
                  << vector_norm(nonlinear_combined.tangent) << '\n';
        const std::vector<double> target_chain =
            cumes_meow_example::calculate_qh_target_jvp(
                primal.equilibrium, primal.profiles, nonlinear_tangent,
                primal.report.input_params, spec);
        std::vector<double> implicit_rz_fd_lambda = spectral.state_tangent;
        std::vector<double> fd_rz_implicit_lambda = nonlinear_state;
        for (const std::size_t lambda_family :
             {static_cast<std::size_t>(cumes::EquilibriumSnapshot::LMNSC),
              static_cast<std::size_t>(cumes::EquilibriumSnapshot::LMNCS)}) {
            const std::size_t begin = lambda_family * family_size;
            const std::size_t end = begin + family_size;
            std::copy(fd_rz_implicit_lambda.begin() + begin,
                      fd_rz_implicit_lambda.begin() + end,
                      implicit_rz_fd_lambda.begin() + begin);
            std::copy(spectral.state_tangent.begin() + begin,
                      spectral.state_tangent.begin() + end,
                      fd_rz_implicit_lambda.begin() + begin);
        }
        const cumes::EquilibriumTangent implicit_rz_fd_lambda_tangent =
            linearization.materialize_tangent(
                implicit_rz_fd_lambda, primal.equilibrium, primal.profiles);
        const cumes::EquilibriumTangent fd_rz_implicit_lambda_tangent =
            linearization.materialize_tangent(
                fd_rz_implicit_lambda, primal.equilibrium, primal.profiles);
        const std::vector<double> implicit_rz_fd_lambda_target =
            cumes_meow_example::calculate_qh_target_jvp(
                primal.equilibrium, primal.profiles,
                implicit_rz_fd_lambda_tangent, primal.report.input_params,
                spec);
        const std::vector<double> fd_rz_implicit_lambda_target =
            cumes_meow_example::calculate_qh_target_jvp(
                primal.equilibrium, primal.profiles,
                fd_rz_implicit_lambda_tangent, primal.report.input_params,
                spec);
        std::vector<double> rz_response_difference(target_chain.size());
        std::vector<double> lambda_response_difference(target_chain.size());
        std::vector<double> total_response_difference(target_chain.size());
        double cross_inner_product = 0.0;
        for (std::size_t row = 0; row < target_chain.size(); ++row) {
            rz_response_difference[row] =
                implicit_rz_fd_lambda_target[row] - target_chain[row];
            lambda_response_difference[row] =
                fd_rz_implicit_lambda_target[row] - target_chain[row];
            total_response_difference[row] = analytic[row] - target_chain[row];
            cross_inner_product +=
                rz_response_difference[row] * lambda_response_difference[row];
        }
        const double rz_difference_norm = vector_norm(rz_response_difference);
        const double lambda_difference_norm =
            vector_norm(lambda_response_difference);
        const double total_difference_norm =
            vector_norm(total_response_difference);
        const double cancellation_cosine =
            cross_inner_product /
            std::max(rz_difference_norm * lambda_difference_norm, 1.0e-30);
        const double cancellation_ratio =
            total_difference_norm /
            std::max(rz_difference_norm + lambda_difference_norm, 1.0e-30);
        std::cout << "QH lambda attribution column=" << boundary.name(column)
                  << " implicit_all_error="
                  << relative_difference(analytic, finite_difference)
                  << " implicit_rz_fd_lambda_error="
                  << relative_difference(implicit_rz_fd_lambda_target,
                                         finite_difference)
                  << " fd_rz_implicit_lambda_error="
                  << relative_difference(fd_rz_implicit_lambda_target,
                                         finite_difference)
                  << " fd_all_chain_error="
                  << relative_difference(target_chain, finite_difference)
                  << " rz_difference_norm=" << rz_difference_norm
                  << " lambda_difference_norm=" << lambda_difference_norm
                  << " total_difference_norm=" << total_difference_norm
                  << " cancellation_cosine=" << cancellation_cosine
                  << " cancellation_ratio=" << cancellation_ratio << '\n';
        const std::string column_name = boundary.name(column);
        constexpr std::array<const char*, cumes::EquilibriumSnapshot::COUNT>
            FAMILY_NAMES = {"RMNCC", "ZMNSC", "LMNSC",
                            "RMNSS", "ZMNCS", "LMNCS"};
        for (std::size_t family = 0; family < cumes::EquilibriumSnapshot::COUNT;
             ++family) {
            print_comparison(column_name, FAMILY_NAMES[family],
                             tangent.equilibrium.families[family],
                             nonlinear_tangent.equilibrium.families[family]);
        }
        constexpr std::array<const char*,
                             cumes::EquilibriumSnapshot::HALF_FIELD_COUNT>
            HALF_FIELD_NAMES = {"SQRTG", "BSUPS", "BSUPU", "BSUPV",
                                "BSUBS", "BSUBU", "BSUBV"};
        for (std::size_t field = 0;
             field < cumes::EquilibriumSnapshot::HALF_FIELD_COUNT; ++field) {
            print_comparison(column_name, HALF_FIELD_NAMES[field],
                             tangent.equilibrium.half_fields[field],
                             nonlinear_tangent.equilibrium.half_fields[field]);
        }
        print_comparison(column_name, "TOROIDAL_FLUX_PRIME",
                         tangent.profiles.toroidal_flux_derivative,
                         nonlinear_tangent.profiles.toroidal_flux_derivative);
        print_comparison(column_name, "POLOIDAL_FLUX_PRIME",
                         tangent.profiles.poloidal_flux_derivative,
                         nonlinear_tangent.profiles.poloidal_flux_derivative);
        print_comparison(column_name, "IOTA",
                         tangent.profiles.rotational_transform,
                         nonlinear_tangent.profiles.rotational_transform);
        print_comparison(column_name, "I_COVARIANT",
                         tangent.profiles.poloidal_covariant_field,
                         nonlinear_tangent.profiles.poloidal_covariant_field);
        print_comparison(column_name, "G_COVARIANT",
                         tangent.profiles.toroidal_covariant_field,
                         nonlinear_tangent.profiles.toroidal_covariant_field);

        const auto analytic_magnetic = magnetic_observable_tangent(
            primal, tangent, primal.report.input_params.nfp);
        const auto primal_magnetic =
            cumes_meow_example::calculate_magnetic_gradient_fields(
                primal.equilibrium, primal.profiles,
                primal.report.input_params.nfp);
        const auto plus_magnetic =
            cumes_meow_example::calculate_magnetic_gradient_fields(
                plus.equilibrium, plus.profiles, plus.report.input_params.nfp);
        print_comparison(column_name, "B", analytic_magnetic.field_strength,
                         difference(plus_magnetic.field_strength,
                                    primal_magnetic.field_strength, step));
        print_comparison(column_name, "B_DOT_GRAD_B",
                         analytic_magnetic.b_dot_grad_b,
                         difference(plus_magnetic.b_dot_grad_b,
                                    primal_magnetic.b_dot_grad_b, step));
        print_comparison(
            column_name, "B_CROSS_GRAD_S_DOT_GRAD_B",
            analytic_magnetic.b_cross_grad_s_dot_grad_b,
            difference(plus_magnetic.b_cross_grad_s_dot_grad_b,
                       primal_magnetic.b_cross_grad_s_dot_grad_b, step));
        print_comparison(
            column_name, "B_CROSS_GRAD_PSI_P_DOT_GRAD_B",
            analytic_magnetic.b_cross_grad_psi_p_dot_grad_b,
            difference(plus_magnetic.b_cross_grad_psi_p_dot_grad_b,
                       primal_magnetic.b_cross_grad_psi_p_dot_grad_b, step));
        const double analytic_aspect =
            cumes_meow_example::calculate_plasma_size_jvp(
                primal.equilibrium, tangent, primal.report.input_params)
                .aspect_ratio;
        const double finite_difference_aspect =
            (plus_target.plasma_size.aspect_ratio -
             primal_target.plasma_size.aspect_ratio) /
            step;
        print_comparison(column_name, "ASPECT", {analytic_aspect},
                         {finite_difference_aspect});
        print_comparison(column_name, "TARGET_RESIDUAL", analytic,
                         finite_difference);
        const double target_chain_error =
            relative_difference(target_chain, finite_difference);
        const double analytic_error =
            relative_difference(analytic, finite_difference);
        const double objective_fd =
            (plus_target.value - primal_target.value) / step;
        const double objective_analytic =
            objective_derivative(primal_target.residuals, analytic);
        const double objective_error =
            std::abs(objective_analytic - objective_fd) /
            std::max(std::abs(objective_fd), 1.0e-12);
        const bool uses_blackbox =
            cumes_meow_example::requires_blackbox_finite_difference(
                boundary.degrees_of_freedom()[column]);
        const double hybrid_policy_error =
            uses_blackbox ? 0.0 : objective_error;
        worst_target_chain_error =
            std::max(worst_target_chain_error, target_chain_error);
        worst_analytic_error = std::max(worst_analytic_error, analytic_error);
        worst_objective_error =
            std::max(worst_objective_error, objective_error);
        worst_hybrid_policy_error =
            std::max(worst_hybrid_policy_error, hybrid_policy_error);
        std::cout << "QH target tangent column=" << boundary.name(column)
                  << " step=" << step
                  << " GMRES_iterations=" << spectral.iterations
                  << " GMRES_relative_residual="
                  << spectral.final_residual / spectral.initial_residual
                  << " target_chain_error=" << target_chain_error
                  << " analytic_residual_error=" << analytic_error
                  << " objective_error=" << objective_error
                  << " hybrid_policy_error=" << hybrid_policy_error
                  << " blackbox=" << uses_blackbox
                  << " objective_analytic=" << objective_analytic
                  << " objective_fd=" << objective_fd << '\n';
        check(std::isfinite(target_chain_error) &&
                  std::isfinite(analytic_error) &&
                  std::isfinite(objective_error),
              "QH target tangent diagnostics are finite for " +
                  boundary.name(column));
    }
    std::cout << "QH target tangent worst_target_chain_error="
              << worst_target_chain_error
              << " worst_analytic_residual_error=" << worst_analytic_error
              << " worst_objective_error=" << worst_objective_error
              << " worst_hybrid_policy_error=" << worst_hybrid_policy_error
              << '\n';
    if (!branch_differences.empty()) {
        Eigen::MatrixXd normalized_branch_matrix(
            static_cast<Eigen::Index>(branch_differences.front().size()),
            static_cast<Eigen::Index>(branch_differences.size()));
        for (std::size_t column = 0; column < branch_differences.size();
             ++column) {
            const double scale =
                std::max(vector_norm(branch_differences[column]), 1.0e-30);
            for (std::size_t row = 0; row < branch_differences[column].size();
                 ++row) {
                normalized_branch_matrix(static_cast<Eigen::Index>(row),
                                         static_cast<Eigen::Index>(column)) =
                    branch_differences[column][row] / scale;
            }
        }
        const Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(
            normalized_branch_matrix,
            Eigen::ComputeThinU | Eigen::ComputeThinV);
        std::cout << "QH normalized branch-difference singular values=";
        for (Eigen::Index index = 0;
             index < decomposition.singularValues().size(); ++index) {
            if (index != 0) std::cout << ',';
            std::cout << decomposition.singularValues()[index];
        }
        std::cout << '\n';
        const Eigen::MatrixXd gram =
            normalized_branch_matrix.transpose() * normalized_branch_matrix;
        for (std::size_t row = 0; row < branch_names.size(); ++row) {
            std::cout << "QH branch correlation row=" << branch_names[row]
                      << " values=";
            for (std::size_t column = 0; column < branch_names.size();
                 ++column) {
                if (column != 0) std::cout << ',';
                std::cout << gram(static_cast<Eigen::Index>(row),
                                  static_cast<Eigen::Index>(column));
            }
            std::cout << '\n';
        }
    }
    check(worst_target_chain_error < 2.0e-2,
          "QH target chain rule agrees with the nonlinear oracle");
    check(worst_hybrid_policy_error < 1.0e-2,
          "hybrid QH objective Jacobian agrees with the nonlinear oracle");

    return meow::test::summary();
}
