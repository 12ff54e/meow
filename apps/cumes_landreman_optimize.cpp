// Execute a JSON-described Landreman-Paul relaxation rundown.
#include "clap.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <filesystem>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cumes/config/json_reader.hpp>
#include <cumes/config/json_writer.hpp>
#include <cumes/config/validated_problem.hpp>
#include <cumes/io/output_spec.hpp>
#include <cumes/io/writer.hpp>
#include <cumes/solver/equilibrium_linearization.hpp>
#include <cumes/solver/equilibrium_solver.hpp>
#include <meow/config/relaxation_rundown.hpp>
#include <meow/cumes/boundary_parameterization.hpp>
#include <meow/cumes/landreman_workflow.hpp>
#include <meow/cumes/quasisymmetry_target.hpp>
#include <meow/cumes/quasisymmetry_target_jvp.hpp>
#include <meow/trf.hpp>

namespace {

using meow::config::RelaxationRundown;

enum class JacobianMethod {
    ANALYTIC,
    BROYDEN,
    AGGRESSIVE_BROYDEN,
    PARALLEL_AGGRESSIVE_BROYDEN,
    FOUR_WORKER_AGGRESSIVE_BROYDEN,
    EXTENDED_FOUR_WORKER_BROYDEN,
    JACOBIAN_SCALED,
    TWO_ACCURACY,
    TWO_ACCURACY_BROYDEN,
    GEOMETRY_RESTART_CHECK,
    GEOMETRY_RESTART_FINITE_DIFFERENCE,
    PARALLEL_FINITE_DIFFERENCE_CHECK,
    RELAXED_PARALLEL_FINITE_DIFFERENCE_CHECK,
    PARALLEL_WORKER_COUNT_CHECK,
    PARALLEL_FINITE_DIFFERENCE,
    HOT_FINITE_DIFFERENCE,
    WARM_FINITE_DIFFERENCE,
    FINITE_DIFFERENCE
};

JacobianMethod parse_jacobian_method(const std::string& name) {
    if (name == "analytic") { return JacobianMethod::ANALYTIC; }
    if (name == "broyden") { return JacobianMethod::BROYDEN; }
    if (name == "aggressive-broyden") {
        return JacobianMethod::AGGRESSIVE_BROYDEN;
    }
    if (name == "parallel-aggressive-broyden") {
        return JacobianMethod::PARALLEL_AGGRESSIVE_BROYDEN;
    }
    if (name == "four-worker-aggressive-broyden") {
        return JacobianMethod::FOUR_WORKER_AGGRESSIVE_BROYDEN;
    }
    if (name == "extended-four-worker-broyden") {
        return JacobianMethod::EXTENDED_FOUR_WORKER_BROYDEN;
    }
    if (name == "jacobian-scaled") { return JacobianMethod::JACOBIAN_SCALED; }
    if (name == "two-accuracy") { return JacobianMethod::TWO_ACCURACY; }
    if (name == "two-accuracy-broyden") {
        return JacobianMethod::TWO_ACCURACY_BROYDEN;
    }
    if (name == "geometry-restart-check") {
        return JacobianMethod::GEOMETRY_RESTART_CHECK;
    }
    if (name == "geometry-restart-finite-difference") {
        return JacobianMethod::GEOMETRY_RESTART_FINITE_DIFFERENCE;
    }
    if (name == "parallel-finite-difference-check") {
        return JacobianMethod::PARALLEL_FINITE_DIFFERENCE_CHECK;
    }
    if (name == "relaxed-parallel-finite-difference-check") {
        return JacobianMethod::RELAXED_PARALLEL_FINITE_DIFFERENCE_CHECK;
    }
    if (name == "parallel-worker-count-check") {
        return JacobianMethod::PARALLEL_WORKER_COUNT_CHECK;
    }
    if (name == "parallel-finite-difference") {
        return JacobianMethod::PARALLEL_FINITE_DIFFERENCE;
    }
    if (name == "hot-finite-difference") {
        return JacobianMethod::HOT_FINITE_DIFFERENCE;
    }
    if (name == "warm-finite-difference") {
        return JacobianMethod::WARM_FINITE_DIFFERENCE;
    }
    if (name == "finite-difference") {
        return JacobianMethod::FINITE_DIFFERENCE;
    }
    throw std::invalid_argument(
        "JACOBIAN_METHOD must be analytic, broyden, aggressive-broyden, "
        "parallel-aggressive-broyden, "
        "four-worker-aggressive-broyden, "
        "extended-four-worker-broyden, "
        "jacobian-scaled, "
        "two-accuracy, two-accuracy-broyden, geometry-restart-check, "
        "geometry-restart-finite-difference, "
        "parallel-finite-difference-check, "
        "relaxed-parallel-finite-difference-check, "
        "parallel-worker-count-check, "
        "parallel-finite-difference, "
        "hot-finite-difference, warm-finite-difference, or "
        "finite-difference");
}

const char* jacobian_method_name(JacobianMethod method) {
    switch (method) {
        case JacobianMethod::ANALYTIC:
            return "analytic";
        case JacobianMethod::BROYDEN:
            return "broyden";
        case JacobianMethod::AGGRESSIVE_BROYDEN:
            return "aggressive-broyden";
        case JacobianMethod::PARALLEL_AGGRESSIVE_BROYDEN:
            return "parallel-aggressive-broyden";
        case JacobianMethod::FOUR_WORKER_AGGRESSIVE_BROYDEN:
            return "four-worker-aggressive-broyden";
        case JacobianMethod::EXTENDED_FOUR_WORKER_BROYDEN:
            return "extended-four-worker-broyden";
        case JacobianMethod::JACOBIAN_SCALED:
            return "jacobian-scaled";
        case JacobianMethod::TWO_ACCURACY:
            return "two-accuracy";
        case JacobianMethod::TWO_ACCURACY_BROYDEN:
            return "two-accuracy-broyden";
        case JacobianMethod::GEOMETRY_RESTART_CHECK:
            return "geometry-restart-check";
        case JacobianMethod::GEOMETRY_RESTART_FINITE_DIFFERENCE:
            return "geometry-restart-finite-difference";
        case JacobianMethod::PARALLEL_FINITE_DIFFERENCE_CHECK:
            return "parallel-finite-difference-check";
        case JacobianMethod::RELAXED_PARALLEL_FINITE_DIFFERENCE_CHECK:
            return "relaxed-parallel-finite-difference-check";
        case JacobianMethod::PARALLEL_WORKER_COUNT_CHECK:
            return "parallel-worker-count-check";
        case JacobianMethod::PARALLEL_FINITE_DIFFERENCE:
            return "parallel-finite-difference";
        case JacobianMethod::HOT_FINITE_DIFFERENCE:
            return "hot-finite-difference";
        case JacobianMethod::WARM_FINITE_DIFFERENCE:
            return "warm-finite-difference";
        case JacobianMethod::FINITE_DIFFERENCE:
            return "finite-difference";
    }
    return "unknown";
}

std::size_t parse_positive_size(std::string_view text,
                                std::string_view option_name,
                                bool eigen_index_bound = false) {
    std::size_t result = 0;
    const char* const begin = text.data();
    const char* const end = begin + text.size();
    const auto [parsed_end, error] = std::from_chars(begin, end, result);
    if (text.empty() || parsed_end != end || error != std::errc{} ||
        result == 0 ||
        (eigen_index_bound &&
         result > static_cast<std::size_t>(
                      std::numeric_limits<Eigen::Index>::max()))) {
        throw std::invalid_argument(std::string(option_name) +
                                    " must be a positive integer");
    }
    return result;
}

std::string first_error(const cumes::ValidationReport& report,
                        const std::string& fallback) {
    const auto errors = report.errors();
    return errors.empty() ? fallback : errors.front();
}

void write_problem(const std::string& path,
                   cumes::ProblemSpec problem,
                   const cumes::SolverOptions& validation_options) {
    cumes::ValidationResult validated =
        cumes::validate(std::move(problem), validation_options);
    if (!validated.has_value()) {
        throw std::runtime_error(
            first_error(validated.error(), "output validation failed"));
    }
    cumes::write_problem_spec(path, validated.value().spec());
}

class LandremanResidual {
   public:
    LandremanResidual(
        cumes::ProblemSpec baseline,
        cumes_meow_example::StellaratorSymmetricBoundaryParameterization
            boundary,
        RelaxationRundown rundown,
        cumes::SolverOptions validation_options,
        const cumes::EquilibriumSnapshot* initial_restart = nullptr,
        bool warm_start_trials = false)
        : baseline_(std::move(baseline)),
          boundary_(std::move(boundary)),
          rundown_(std::move(rundown)),
          validation_options_(validation_options),
          initial_restart_(initial_restart),
          warm_start_trials_(warm_start_trials) {}

    meow::Vector operator()(const meow::Vector& x) {
        if (cached_x_.has_value() && cached_x_->size() == x.size() &&
            (cached_x_->array() == x.array()).all()) {
            return cached_residual_;
        }

        cumes::ProblemSpec problem = trial_problem(x);
        cumes::ValidationResult validated =
            cumes::validate(std::move(problem), validation_options_);
        if (!validated.has_value()) {
            throw std::runtime_error(first_error(
                validated.error(), "cuMES boundary validation failed"));
        }

        cumes::SolveRequest request;
        const bool used_initial_restart = initial_restart_ != nullptr;
        const cumes::EquilibriumSnapshot* restart = initial_restart_;
        if (restart == nullptr && warm_start_trials_ &&
            accepted_equilibrium_.has_value()) {
            restart = &*accepted_equilibrium_;
        }
        const bool used_restart = restart != nullptr;
        if (used_restart) { request.restart = std::cref(*restart); }
        if (rundown_.initialization.radial_transfer ==
            meow::config::RadialTransfer::CATMULL_ROM) {
            request.radial_transfer = cumes::RadialTransferPolicy::CATMULL_ROM;
        }
        std::optional<cumes::ValidationResult> restart_validated;
        if (used_restart) {
            cumes::ProblemSpec restart_problem = validated.value().spec();
            restart_problem.stages = {restart_problem.stages.back()};
            restart_validated.emplace(cumes::validate(
                std::move(restart_problem), validation_options_));
            if (!restart_validated->has_value()) {
                throw std::runtime_error(
                    first_error(restart_validated->error(),
                                "stage-restart boundary validation failed"));
            }
        }
        cumes::SolveOutcome solved = solver_.solve(
            used_restart ? restart_validated->value() : validated.value(),
            request);
        initial_restart_ = nullptr;
        if (used_restart) { ++restart_solve_attempts_; }
        ++evaluation_count_;
        total_nonlinear_iterations_ += solved.total_iterations;
        if (used_restart &&
            (!solved.converged || !solved.has_complete_equilibrium())) {
            ++restart_solve_failures_;
            std::cout << "equilibrium_restart_failed=1 initial_stage_restart="
                      << used_initial_restart
                      << " failed_stage=" << solved.failed_stage
                      << " fsqr=" << solved.fsqr << " fsqz=" << solved.fsqz
                      << " fsql=" << solved.fsql << " cold_fallback=1\n";
            request.restart.reset();
            solved = solver_.solve(validated.value(), request);
            ++evaluation_count_;
            total_nonlinear_iterations_ += solved.total_iterations;
        }
        if (!solved.converged || !solved.has_complete_equilibrium()) {
            if (cached_residual_.size() == 0) {
                throw std::runtime_error(
                    "initial cuMES equilibrium failed at objective "
                    "evaluation " +
                    std::to_string(evaluation_count_) + ", stage " +
                    std::to_string(solved.failed_stage) + ", residuals (" +
                    std::to_string(solved.fsqr) + ", " +
                    std::to_string(solved.fsqz) + ", " +
                    std::to_string(solved.fsql) + ")");
            }

            // A failed off-center trust-region trial is outside the feasible
            // equilibrium domain. Return a finite barrier residual so TRF can
            // reject and contract the step instead of terminating the run.
            meow::Vector rejected = meow::Vector::Zero(cached_residual_.size());
            rejected[0] = 100.0 * (1.0 + cached_residual_.norm());
            std::cout << "evaluation=" << evaluation_count_
                      << " rejected_equilibrium=1 failed_stage="
                      << solved.failed_stage << " fsqr=" << solved.fsqr
                      << " fsqz=" << solved.fsqz << " fsql=" << solved.fsql
                      << " penalty_objective=" << rejected.squaredNorm()
                      << '\n';
            return rejected;
        }

        const cumes_meow_example::CompositeQuasisymmetryTarget target =
            calculate_target(solved);

        meow::Vector residual(
            static_cast<Eigen::Index>(target.residuals.size()));
        for (std::size_t index = 0; index < target.residuals.size(); ++index) {
            residual[static_cast<Eigen::Index>(index)] =
                target.residuals[index];
        }
        std::cout << std::setprecision(12) << "evaluation=" << evaluation_count_
                  << " objective=" << target.value << " qs=" << target.qs.value
                  << " aspect=" << target.plasma_size.aspect_ratio
                  << " mean_iota=" << target.iota_integral
                  << " solver_iterations=" << solved.total_iterations << '\n';

        cached_x_ = x;
        cached_residual_ = residual;
        cached_outcome_ = std::move(solved);
        if (warm_start_trials_ && !accepted_equilibrium_.has_value()) {
            accepted_equilibrium_ = cached_outcome_->equilibrium;
        }
        return residual;
    }

    meow::Matrix jacobian(const meow::Vector& x) {
        static_cast<void>((*this)(x));
        if (!cached_outcome_.has_value()) {
            throw std::runtime_error(
                "no equilibrium is available for analytic Jacobian");
        }
        cumes::ProblemSpec problem = trial_problem(x);
        cumes::ValidationResult validated =
            cumes::validate(std::move(problem), validation_options_);
        if (!validated.has_value()) {
            throw std::runtime_error(first_error(
                validated.error(), "Jacobian boundary validation failed"));
        }
        cumes::EquilibriumLinearization linearization(
            validated.value(), cached_outcome_->equilibrium);
        meow::Matrix result(cached_residual_.size(), x.size());
        const auto spec = target_spec(cached_outcome_->report.input_params.nfp);
        cumes::TangentLinearOptions tangent_options;
        tangent_options.max_iterations = 1000;
        tangent_options.restart = 300;
        tangent_options.relative_tolerance = 5.0e-6;
        tangent_options.absolute_tolerance = 1.0e-11;
        std::size_t total_linear_iterations = 0;
        double worst_linear_residual = 0.0;
        double worst_relative_linear_residual = 0.0;
        std::vector<std::size_t> fallback_columns;
        const auto& degrees = boundary_.degrees_of_freedom();
        for (std::size_t column = 0; column < boundary_.size(); ++column) {
            if (cumes_meow_example::requires_blackbox_finite_difference(
                    degrees[column])) {
                fallback_columns.push_back(column);
                continue;
            }
            const cumes::BoundaryTangent boundary_tangent =
                boundary_.tangent(validated.value(), column);
            const cumes::SpectralTangentSolve spectral =
                linearization.solve_boundary_tangent(boundary_tangent,
                                                     tangent_options);
            total_linear_iterations +=
                static_cast<std::size_t>(spectral.iterations);
            worst_linear_residual =
                std::max(worst_linear_residual, spectral.final_residual);
            const double relative_linear_residual =
                spectral.initial_residual == 0.0
                    ? spectral.final_residual
                    : spectral.final_residual / spectral.initial_residual;
            worst_relative_linear_residual = std::max(
                worst_relative_linear_residual, relative_linear_residual);
            if (!spectral.converged ||
                !std::isfinite(relative_linear_residual)) {
                fallback_columns.push_back(column);
                continue;
            }
            const cumes::EquilibriumTangent tangent =
                linearization.materialize_tangent(spectral.state_tangent,
                                                  cached_outcome_->equilibrium,
                                                  cached_outcome_->profiles);
            const std::vector<double> target_tangent =
                is_qa() ? cumes_meow_example::calculate_qa_target_jvp(
                              cached_outcome_->equilibrium,
                              cached_outcome_->profiles, tangent,
                              cached_outcome_->report.input_params, spec,
                              rundown_.target.mean_iota.has_value())
                        : cumes_meow_example::calculate_qh_target_jvp(
                              cached_outcome_->equilibrium,
                              cached_outcome_->profiles, tangent,
                              cached_outcome_->report.input_params, spec);
            if (target_tangent.size() !=
                static_cast<std::size_t>(result.rows())) {
                throw std::runtime_error(
                    "analytic target tangent has the wrong length");
            }
            if (!std::all_of(
                    target_tangent.begin(), target_tangent.end(),
                    [](double value) { return std::isfinite(value); })) {
                fallback_columns.push_back(column);
                continue;
            }
            for (std::size_t row = 0; row < target_tangent.size(); ++row) {
                result(static_cast<Eigen::Index>(row),
                       static_cast<Eigen::Index>(column)) = target_tangent[row];
            }
        }

        const meow::Vector primal_residual = cached_residual_;
        cumes::SolveOutcome primal_outcome = std::move(*cached_outcome_);
        const auto& finite_difference = rundown_.optimizer.jacobian;
        for (const std::size_t column : fallback_columns) {
            const double step =
                std::max(finite_difference.relative_step *
                             std::abs(x[static_cast<Eigen::Index>(column)]),
                         finite_difference.absolute_step);
            meow::Vector perturbed = x;
            perturbed[static_cast<Eigen::Index>(column)] += step;
            const meow::Vector perturbed_residual = (*this)(perturbed);
            if (perturbed_residual.size() != result.rows()) {
                throw std::runtime_error(
                    "black-box fallback residual has the wrong length");
            }
            result.col(static_cast<Eigen::Index>(column)) =
                (perturbed_residual - primal_residual) / step;
        }
        cached_x_ = x;
        cached_residual_ = primal_residual;
        cached_outcome_ = std::move(primal_outcome);

        std::cout << "analytic_jacobian columns=" << result.cols()
                  << " residuals=" << result.rows() << " tangent_columns="
                  << result.cols() -
                         static_cast<Eigen::Index>(fallback_columns.size())
                  << " blackbox_columns=" << fallback_columns.size()
                  << " linear_iterations=" << total_linear_iterations
                  << " worst_linear_residual=" << worst_linear_residual
                  << " worst_relative_linear_residual="
                  << worst_relative_linear_residual << '\n';
        ++analytic_jacobian_evaluations_;
        total_linear_iterations_ += total_linear_iterations;
        return result;
    }

    meow::Matrix restart_finite_difference_jacobian(const meow::Vector& x,
                                                    double absolute_step_floor,
                                                    const char* label,
                                                    bool reset_lambda = false) {
        static_cast<void>((*this)(x));
        if (!cached_outcome_.has_value()) {
            throw std::runtime_error(
                "no equilibrium is available for hot-restart Jacobian");
        }
        const cumes::SolveOutcome& primal = *cached_outcome_;
        const meow::Vector primal_residual = cached_residual_;
        meow::Matrix result(primal_residual.size(), x.size());
        const auto& finite_difference = rundown_.optimizer.jacobian;
        cumes::EquilibriumSnapshot restart_equilibrium = primal.equilibrium;
        if (reset_lambda) {
            std::fill(
                restart_equilibrium.families[cumes::EquilibriumSnapshot::LMNSC]
                    .begin(),
                restart_equilibrium.families[cumes::EquilibriumSnapshot::LMNSC]
                    .end(),
                0.0);
            std::fill(
                restart_equilibrium.families[cumes::EquilibriumSnapshot::LMNCS]
                    .begin(),
                restart_equilibrium.families[cumes::EquilibriumSnapshot::LMNCS]
                    .end(),
                0.0);
        }
        std::size_t jacobian_nonlinear_iterations = 0;
        std::size_t cold_fallbacks = 0;
        std::size_t backward_fallbacks = 0;

        const auto solve_perturbation = [&](const meow::Vector& perturbed_x) {
            cumes::ProblemSpec full_problem = trial_problem(perturbed_x);
            cumes::ProblemSpec hot_problem = full_problem;
            hot_problem.stages = {hot_problem.stages.back()};
            cumes::ValidationResult validated =
                cumes::validate(std::move(hot_problem), validation_options_);
            if (!validated.has_value()) {
                throw std::runtime_error(first_error(
                    validated.error(),
                    "hot-restart Jacobian boundary validation failed"));
            }

            cumes::SolveRequest request;
            request.restart = std::cref(restart_equilibrium);
            if (rundown_.initialization.radial_transfer ==
                meow::config::RadialTransfer::CATMULL_ROM) {
                request.radial_transfer =
                    cumes::RadialTransferPolicy::CATMULL_ROM;
            }
            cumes::SolveOutcome solved =
                solver_.solve(validated.value(), request);
            ++evaluation_count_;
            total_nonlinear_iterations_ += solved.total_iterations;
            jacobian_nonlinear_iterations += solved.total_iterations;
            if (solved.converged && solved.has_complete_equilibrium()) {
                return solved;
            }

            cumes::ValidationResult cold_validated =
                cumes::validate(std::move(full_problem), validation_options_);
            if (!cold_validated.has_value()) {
                throw std::runtime_error(first_error(
                    cold_validated.error(),
                    "cold Jacobian fallback boundary validation failed"));
            }
            request.restart.reset();
            solved = solver_.solve(cold_validated.value(), request);
            ++evaluation_count_;
            ++cold_fallbacks;
            total_nonlinear_iterations_ += solved.total_iterations;
            jacobian_nonlinear_iterations += solved.total_iterations;
            return solved;
        };

        for (std::size_t column = 0; column < boundary_.size(); ++column) {
            const double step = std::max(
                finite_difference.relative_step *
                    std::abs(x[static_cast<Eigen::Index>(column)]),
                std::max(finite_difference.absolute_step, absolute_step_floor));
            meow::Vector perturbed_x = x;
            perturbed_x[static_cast<Eigen::Index>(column)] += step;
            cumes::SolveOutcome solved = solve_perturbation(perturbed_x);
            double difference_direction = 1.0;
            if (!solved.converged || !solved.has_complete_equilibrium()) {
                perturbed_x = x;
                perturbed_x[static_cast<Eigen::Index>(column)] -= step;
                solved = solve_perturbation(perturbed_x);
                difference_direction = -1.0;
                ++backward_fallbacks;
            }
            if (!solved.converged || !solved.has_complete_equilibrium()) {
                throw std::runtime_error(
                    "both finite-difference directions failed for " +
                    boundary_.name(column));
            }

            const cumes_meow_example::CompositeQuasisymmetryTarget target =
                calculate_target(solved);
            if (target.residuals.size() !=
                static_cast<std::size_t>(result.rows())) {
                throw std::runtime_error(
                    "hot-restart target residual has the wrong length");
            }
            for (std::size_t row = 0; row < target.residuals.size(); ++row) {
                result(static_cast<Eigen::Index>(row),
                       static_cast<Eigen::Index>(column)) =
                    difference_direction *
                    (target.residuals[row] -
                     primal_residual[static_cast<Eigen::Index>(row)]) /
                    step;
            }
        }

        std::cout << label << " columns=" << result.cols()
                  << " residuals=" << result.rows()
                  << " nonlinear_iterations=" << jacobian_nonlinear_iterations
                  << " cold_fallbacks=" << cold_fallbacks
                  << " backward_fallbacks=" << backward_fallbacks
                  << " absolute_step_floor=" << absolute_step_floor << '\n';
        ++restart_jacobian_evaluations_;
        return result;
    }

    meow::Matrix hot_restart_jacobian(const meow::Vector& x) {
        return restart_finite_difference_jacobian(x, 1.0e-4,
                                                  "hot_restart_jacobian");
    }

    meow::Matrix warm_restart_jacobian(const meow::Vector& x) {
        return restart_finite_difference_jacobian(x, 0.0,
                                                  "warm_restart_jacobian");
    }

    meow::Matrix geometry_restart_jacobian(const meow::Vector& x) {
        return restart_finite_difference_jacobian(
            x, 0.0, "geometry_restart_jacobian", true);
    }

    meow::Matrix cold_finite_difference_jacobian(const meow::Vector& x) {
        const auto start = std::chrono::steady_clock::now();
        static_cast<void>((*this)(x));
        if (!cached_outcome_.has_value()) {
            throw std::runtime_error(
                "no equilibrium is available for cold Jacobian");
        }
        const meow::Vector primal_residual = cached_residual_;
        cumes::SolveOutcome primal_outcome = std::move(*cached_outcome_);
        meow::Matrix result(primal_residual.size(), x.size());
        const auto& finite_difference = rundown_.optimizer.jacobian;
        for (Eigen::Index column = 0; column < x.size(); ++column) {
            const double step =
                std::max(finite_difference.relative_step * std::abs(x[column]),
                         finite_difference.absolute_step);
            meow::Vector perturbed = x;
            perturbed[column] += step;
            const meow::Vector perturbed_residual = (*this)(perturbed);
            result.col(column) = (perturbed_residual - primal_residual) / step;
        }
        cached_x_ = x;
        cached_residual_ = primal_residual;
        cached_outcome_ = std::move(primal_outcome);
        const double wall_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          start)
                .count();
        std::cout << "cold_finite_difference_jacobian columns=" << result.cols()
                  << " residuals=" << result.rows()
                  << " wall_seconds=" << wall_seconds << '\n';
        return result;
    }

    meow::Matrix parallel_finite_difference_jacobian(
        const meow::Vector& x,
        std::optional<double> relaxed_tolerance = std::nullopt,
        std::size_t worker_count = 2) {
        struct WorkerResult {
            cumes::SolveOutcome solved;
            std::optional<cumes_meow_example::CompositeQuasisymmetryTarget>
                target;
            double validation_seconds = 0.0;
            double solver_construction_seconds = 0.0;
            double solve_seconds = 0.0;
            double target_seconds = 0.0;
            double worker_seconds = 0.0;
        };

        if (worker_count == 0) {
            throw std::invalid_argument(
                "parallel Jacobian worker count must be positive");
        }
        static_cast<void>((*this)(x));
        if (!cached_outcome_.has_value()) {
            throw std::runtime_error(
                "no equilibrium is available for parallel Jacobian");
        }
        const auto start = std::chrono::steady_clock::now();
        const meow::Vector primal_residual = cached_residual_;
        meow::Matrix result(primal_residual.size(), x.size());
        const auto& finite_difference = rundown_.optimizer.jacobian;
        std::size_t jacobian_nonlinear_iterations = 0;
        double problem_construction_seconds = 0.0;
        double validation_seconds = 0.0;
        double solver_construction_seconds = 0.0;
        double solve_seconds = 0.0;
        double device_seconds = 0.0;
        double cumes_setup_seconds = 0.0;
        double cumes_multigrid_seconds = 0.0;
        double cumes_transfer_seconds = 0.0;
        double cumes_total_seconds = 0.0;
        double stage_setup_seconds = 0.0;
        double stage_iteration_seconds = 0.0;
        double stage_output_seconds = 0.0;
        double stage_teardown_seconds = 0.0;
        double multigrid_other_seconds = 0.0;
        double target_seconds = 0.0;
        double assembly_seconds = 0.0;
        double critical_worker_seconds = 0.0;
        double batch_tail_idle_seconds = 0.0;

        const Eigen::Index workers = static_cast<Eigen::Index>(worker_count);
        for (Eigen::Index first = 0; first < x.size(); first += workers) {
            const Eigen::Index batch_size = std::min(workers, x.size() - first);
            std::vector<std::future<WorkerResult>> futures;
            std::vector<double> steps;
            std::vector<double> batch_worker_seconds;
            futures.reserve(static_cast<std::size_t>(batch_size));
            steps.reserve(static_cast<std::size_t>(batch_size));
            batch_worker_seconds.reserve(static_cast<std::size_t>(batch_size));
            for (Eigen::Index offset = 0; offset < batch_size; ++offset) {
                const auto problem_start = std::chrono::steady_clock::now();
                const Eigen::Index column = first + offset;
                const double step = std::max(
                    finite_difference.relative_step * std::abs(x[column]),
                    finite_difference.absolute_step);
                meow::Vector perturbed = x;
                perturbed[column] += step;
                cumes::ProblemSpec problem = trial_problem(perturbed);
                if (relaxed_tolerance.has_value()) {
                    for (cumes::StageRequest& request : problem.stages) {
                        request.tolerance =
                            std::max(request.tolerance, *relaxed_tolerance);
                    }
                }
                problem_construction_seconds +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - problem_start)
                        .count();
                steps.push_back(step);
                futures.emplace_back(std::async(
                    std::launch::async,
                    [this, problem = std::move(problem),
                     validation_options = validation_options_,
                     use_catmull =
                         rundown_.initialization.radial_transfer ==
                         meow::config::RadialTransfer::CATMULL_ROM]() mutable {
                        const auto worker_start =
                            std::chrono::steady_clock::now();
                        const auto validation_start = worker_start;
                        cumes::ValidationResult validated = cumes::validate(
                            std::move(problem), validation_options);
                        if (!validated.has_value()) {
                            throw std::runtime_error(first_error(
                                validated.error(),
                                "parallel Jacobian boundary validation "
                                "failed"));
                        }
                        const auto construction_start =
                            std::chrono::steady_clock::now();
                        cumes::EquilibriumSolver solver;
                        const auto solve_start =
                            std::chrono::steady_clock::now();
                        cumes::SolveRequest request;
                        if (use_catmull) {
                            request.radial_transfer =
                                cumes::RadialTransferPolicy::CATMULL_ROM;
                        }
                        cumes::SolveOutcome solved =
                            solver.solve(validated.value(), request);
                        const auto solve_end = std::chrono::steady_clock::now();
                        std::optional<
                            cumes_meow_example::CompositeQuasisymmetryTarget>
                            target;
                        if (solved.converged &&
                            solved.has_complete_equilibrium()) {
                            target.emplace(calculate_target(solved));
                        }
                        const auto worker_end =
                            std::chrono::steady_clock::now();
                        WorkerResult result;
                        result.solved = std::move(solved);
                        result.target = std::move(target);
                        result.validation_seconds =
                            std::chrono::duration<double>(construction_start -
                                                          validation_start)
                                .count();
                        result.solver_construction_seconds =
                            std::chrono::duration<double>(solve_start -
                                                          construction_start)
                                .count();
                        result.solve_seconds = std::chrono::duration<double>(
                                                   solve_end - solve_start)
                                                   .count();
                        result.target_seconds = std::chrono::duration<double>(
                                                    worker_end - solve_end)
                                                    .count();
                        result.worker_seconds = std::chrono::duration<double>(
                                                    worker_end - worker_start)
                                                    .count();
                        return result;
                    }));
            }
            for (Eigen::Index offset = 0; offset < batch_size; ++offset) {
                WorkerResult worker_result =
                    futures[static_cast<std::size_t>(offset)].get();
                validation_seconds += worker_result.validation_seconds;
                solver_construction_seconds +=
                    worker_result.solver_construction_seconds;
                solve_seconds += worker_result.solve_seconds;
                target_seconds += worker_result.target_seconds;
                batch_worker_seconds.push_back(worker_result.worker_seconds);
                cumes::SolveOutcome solved = std::move(worker_result.solved);
                device_seconds += solved.total_device_time_ms * 1.0e-3;
                cumes_setup_seconds += solved.timings.setup_wall_ms * 1.0e-3;
                cumes_multigrid_seconds +=
                    solved.timings.multigrid_wall_ms * 1.0e-3;
                cumes_transfer_seconds +=
                    solved.timings.final_state_transfer_wall_ms * 1.0e-3;
                cumes_total_seconds += solved.timings.total_wall_ms * 1.0e-3;
                stage_setup_seconds +=
                    solved.timings.stage_setup_wall_ms * 1.0e-3;
                stage_iteration_seconds +=
                    solved.timings.stage_iteration_wall_ms * 1.0e-3;
                stage_output_seconds +=
                    solved.timings.stage_output_wall_ms * 1.0e-3;
                stage_teardown_seconds +=
                    solved.timings.stage_teardown_wall_ms * 1.0e-3;
                multigrid_other_seconds +=
                    solved.timings.multigrid_other_wall_ms * 1.0e-3;
                ++evaluation_count_;
                total_nonlinear_iterations_ += solved.total_iterations;
                jacobian_nonlinear_iterations += solved.total_iterations;
                if (!solved.converged || !solved.has_complete_equilibrium()) {
                    std::ostringstream message;
                    message << "parallel finite-difference equilibrium failed"
                            << " column=" << (first + offset)
                            << " workers=" << worker_count
                            << " failed_stage=" << solved.failed_stage
                            << " iterations=" << solved.total_iterations
                            << " fsqr=" << solved.fsqr
                            << " fsqz=" << solved.fsqz
                            << " fsql=" << solved.fsql;
                    throw std::runtime_error(message.str());
                }
                if (!worker_result.target.has_value()) {
                    throw std::runtime_error(
                        "parallel worker did not produce target residuals");
                }
                const auto& target = *worker_result.target;
                if (target.residuals.size() !=
                    static_cast<std::size_t>(result.rows())) {
                    throw std::runtime_error(
                        "parallel target residual has the wrong length");
                }
                const Eigen::Index column = first + offset;
                const double step = steps[static_cast<std::size_t>(offset)];
                const auto assembly_start = std::chrono::steady_clock::now();
                for (std::size_t row = 0; row < target.residuals.size();
                     ++row) {
                    result(static_cast<Eigen::Index>(row), column) =
                        (target.residuals[row] -
                         primal_residual[static_cast<Eigen::Index>(row)]) /
                        step;
                }
                assembly_seconds +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - assembly_start)
                        .count();
            }
            const double critical_worker = *std::max_element(
                batch_worker_seconds.begin(), batch_worker_seconds.end());
            critical_worker_seconds += critical_worker;
            for (const double worker_seconds : batch_worker_seconds) {
                batch_tail_idle_seconds += critical_worker - worker_seconds;
            }
        }
        const double wall_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          start)
                .count();
        ++parallel_jacobian_evaluations_;
        std::cout << "parallel_finite_difference_jacobian columns="
                  << result.cols() << " residuals=" << result.rows()
                  << " workers=" << worker_count
                  << " nonlinear_iterations=" << jacobian_nonlinear_iterations
                  << " relaxed_tolerance=" << relaxed_tolerance.value_or(0.0)
                  << " wall_seconds=" << wall_seconds
                  << " problem_seconds=" << problem_construction_seconds
                  << " validation_seconds=" << validation_seconds
                  << " solver_construction_seconds="
                  << solver_construction_seconds
                  << " solve_seconds=" << solve_seconds
                  << " device_seconds=" << device_seconds
                  << " cumes_setup_seconds=" << cumes_setup_seconds
                  << " cumes_multigrid_seconds=" << cumes_multigrid_seconds
                  << " cumes_transfer_seconds=" << cumes_transfer_seconds
                  << " cumes_total_seconds=" << cumes_total_seconds
                  << " stage_setup_seconds=" << stage_setup_seconds
                  << " stage_iteration_seconds=" << stage_iteration_seconds
                  << " stage_output_seconds=" << stage_output_seconds
                  << " stage_teardown_seconds=" << stage_teardown_seconds
                  << " multigrid_other_seconds=" << multigrid_other_seconds
                  << " critical_worker_seconds=" << critical_worker_seconds
                  << " batch_tail_idle_seconds=" << batch_tail_idle_seconds
                  << " target_seconds=" << target_seconds
                  << " assembly_seconds=" << assembly_seconds << '\n';
        return result;
    }

    meow::Matrix relaxed_parallel_finite_difference_jacobian(
        const meow::Vector& x,
        std::size_t worker_count) {
        return parallel_finite_difference_jacobian(x, 2.0e-12, worker_count);
    }

    meow::Matrix configured_parallel_finite_difference_jacobian(
        const meow::Vector& x,
        std::size_t worker_count) {
        return parallel_finite_difference_jacobian(x, std::nullopt,
                                                   worker_count);
    }

    meow::Matrix eight_worker_finite_difference_jacobian(
        const meow::Vector& x) {
        return parallel_finite_difference_jacobian(x, std::nullopt, 8);
    }

    std::size_t equilibrium_evaluations() const { return evaluation_count_; }

    std::size_t total_nonlinear_iterations() const {
        return total_nonlinear_iterations_;
    }

    std::size_t analytic_jacobian_evaluations() const {
        return analytic_jacobian_evaluations_;
    }

    std::size_t hot_restart_jacobian_evaluations() const {
        return restart_jacobian_evaluations_;
    }

    std::size_t parallel_jacobian_evaluations() const {
        return parallel_jacobian_evaluations_;
    }

    std::size_t restart_solve_attempts() const {
        return restart_solve_attempts_;
    }

    std::size_t restart_solve_failures() const {
        return restart_solve_failures_;
    }

    std::size_t total_linear_iterations() const {
        return total_linear_iterations_;
    }

    const cumes::EquilibriumSnapshot& equilibrium() const {
        if (!cached_outcome_.has_value()) {
            throw std::runtime_error("no accepted equilibrium is available");
        }
        return cached_outcome_->equilibrium;
    }

    void write_equilibrium(const meow::Vector& x, const std::string& path) {
        static_cast<void>((*this)(x));
        if (!cached_outcome_.has_value()) {
            throw std::runtime_error("no equilibrium is available to write");
        }

        cumes::ProblemSpec problem = trial_problem(x);
        cumes::ValidationResult validated =
            cumes::validate(std::move(problem), validation_options_);
        if (!validated.has_value()) {
            throw std::runtime_error(first_error(
                validated.error(), "equilibrium output validation failed"));
        }
        auto writer = cumes::make_binary_writer();
        if (!writer) {
            throw std::runtime_error("cuMES binary writer is unavailable");
        }
        cumes::OutputSpec output_spec;
        output_spec.format = cumes::OutputFormat::BINARY;
        output_spec.path = path;
        cumes::RunReport report = cached_outcome_->report;
        switch (validation_options_.precision) {
            case cumes::PrecisionPolicy::VERIFY_DOUBLE:
                report.build.scalar_type = "double";
                report.build.precision_policy = "verify-double";
                break;
            case cumes::PrecisionPolicy::FAST_DOUBLE:
                report.build.scalar_type = "double";
                report.build.precision_policy = "fast-double";
                break;
            case cumes::PrecisionPolicy::MIXED_FLOAT:
                report.build.scalar_type = "float";
                report.build.precision_policy = "mixed-float";
                break;
            case cumes::PrecisionPolicy::DEBUG_DOUBLE:
                report.build.scalar_type = "double";
                report.build.precision_policy = "debug-double";
                break;
        }
        const cumes::Status status =
            writer->write_atomic(cached_outcome_->equilibrium, report,
                                 output_spec, validated.value());
        if (!status.has_value()) {
            throw std::runtime_error("equilibrium output failed: " +
                                     status.error());
        }
    }

    cumes::ProblemSpec accept(const meow::Vector& x) {
        static_cast<void>((*this)(x));
        if (!cached_outcome_.has_value()) {
            throw std::runtime_error("no equilibrium is available to accept");
        }
        baseline_ = boundary_.apply(baseline_, x);
        if (rundown_.initialization.track_axis_predictor) {
            cumes_meow_example::refresh_axis_predictor_from_equilibrium(
                baseline_, cached_outcome_->equilibrium);
        }
        if (warm_start_trials_) {
            accepted_equilibrium_ = cached_outcome_->equilibrium;
        }
        return baseline_;
    }

   private:
    cumes::ProblemSpec trial_problem(const meow::Vector& x) const {
        cumes::ProblemSpec problem = boundary_.apply(baseline_, x);
        if (rundown_.initialization.track_axis_predictor) {
            cumes_meow_example::track_axis_predictor_from_accepted_boundary(
                problem, baseline_);
        }
        return problem;
    }

    cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec target_spec(
        int nfp) const {
        cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec spec;
        spec.helicity_m = rundown_.target.helicity_m;
        spec.helicity_n = rundown_.target.helicity_n_per_field_period * nfp;
        spec.target_ntheta = rundown_.target.ntheta;
        spec.target_nzeta = rundown_.target.nzeta;
        spec.normalized_toroidal_flux_surfaces =
            rundown_.target.normalized_toroidal_flux_surfaces;
        spec.surface_weights = rundown_.target.surface_weights;
        return spec;
    }

    bool is_qa() const {
        return rundown_.selected_case == meow::config::QuasisymmetryCase::QA;
    }

    cumes_meow_example::CompositeQuasisymmetryTarget calculate_target(
        const cumes::SolveOutcome& solved) const {
        const auto spec = target_spec(solved.report.input_params.nfp);
        return is_qa() ? cumes_meow_example::calculate_qa_target(
                             solved.equilibrium, solved.profiles,
                             solved.report.input_params, spec,
                             rundown_.target.aspect_ratio,
                             rundown_.target.mean_iota)
                       : cumes_meow_example::calculate_qh_target(
                             solved.equilibrium, solved.profiles,
                             solved.report.input_params, spec,
                             rundown_.target.aspect_ratio);
    }

    cumes::ProblemSpec baseline_;
    cumes_meow_example::StellaratorSymmetricBoundaryParameterization boundary_;
    RelaxationRundown rundown_;
    cumes::SolverOptions validation_options_;
    cumes::EquilibriumSolver solver_;
    const cumes::EquilibriumSnapshot* initial_restart_ = nullptr;
    bool warm_start_trials_ = false;
    std::size_t evaluation_count_ = 0;
    std::size_t total_nonlinear_iterations_ = 0;
    std::size_t analytic_jacobian_evaluations_ = 0;
    std::size_t restart_jacobian_evaluations_ = 0;
    std::size_t parallel_jacobian_evaluations_ = 0;
    std::size_t restart_solve_attempts_ = 0;
    std::size_t restart_solve_failures_ = 0;
    std::size_t total_linear_iterations_ = 0;
    std::optional<meow::Vector> cached_x_;
    meow::Vector cached_residual_;
    std::optional<cumes::SolveOutcome> cached_outcome_;
    std::optional<cumes::EquilibriumSnapshot> accepted_equilibrium_;
};

std::string step_stem(const std::string& directory,
                      meow::config::WorkflowKind workflow,
                      int max_mode,
                      std::size_t iteration,
                      std::string_view phase = {}) {
    std::ostringstream name;
    if (workflow == meow::config::WorkflowKind::CONSTRUCTION) {
        name << meow::config::workflow_kind_name(workflow) << "-";
    }
    name << "mode" << max_mode << "_step_" << std::setw(4) << std::setfill('0')
         << iteration;
    if (!phase.empty()) { name << "-" << phase; }
    return (std::filesystem::path(directory) / name.str()).string();
}

std::string checkpoint_path(const std::string& output_path,
                            meow::config::WorkflowKind workflow,
                            int max_mode) {
    if (workflow == meow::config::WorkflowKind::REFINEMENT) {
        return output_path + ".mode" + std::to_string(max_mode) + ".json";
    }
    return output_path + ".construction.mode" + std::to_string(max_mode) +
           ".json";
}

void set_stage_resolution(cumes::ProblemSpec& problem,
                          const meow::config::RelaxationStepSpec& step) {
    const meow::config::EquilibriumStepSpec& stage = step.equilibrium;
    if (stage.mpol == 0 || stage.ntor == 0) { return; }
    const std::vector<double> source_raxis = problem.raxis_c;
    const std::vector<double> source_zaxis = problem.zaxis_s;
    problem.mpol = stage.mpol;
    problem.ntor = stage.ntor;
    for (cumes::StageRequest& request : problem.stages) {
        request.max_iterations =
            std::max(request.max_iterations,
                     static_cast<std::size_t>(stage.minimum_iterations));
        request.tolerance = std::max(request.tolerance, stage.tolerance_floor);
    }
    const std::size_t axis_size = static_cast<std::size_t>(stage.ntor + 1);
    problem.raxis_c.assign(axis_size, 0.0);
    problem.zaxis_s.assign(axis_size, 0.0);
    for (std::size_t index = 0;
         index < axis_size && index < source_raxis.size(); ++index) {
        problem.raxis_c[index] = source_raxis[index];
    }
    for (std::size_t index = 0;
         index < axis_size && index < source_zaxis.size(); ++index) {
        problem.zaxis_s[index] = source_zaxis[index];
    }
}

std::vector<double> stage_tolerances(const cumes::ProblemSpec& problem) {
    std::vector<double> tolerances;
    tolerances.reserve(problem.stages.size());
    for (const cumes::StageRequest& request : problem.stages) {
        tolerances.push_back(request.tolerance);
    }
    return tolerances;
}

void set_stage_tolerances(cumes::ProblemSpec& problem,
                          const std::vector<double>& qualified,
                          std::optional<double> relaxed_tolerance) {
    if (qualified.size() != problem.stages.size()) {
        throw std::invalid_argument("stage tolerance count changed");
    }
    for (std::size_t index = 0; index < problem.stages.size(); ++index) {
        problem.stages[index].tolerance =
            relaxed_tolerance.has_value()
                ? std::max(qualified[index], *relaxed_tolerance)
                : qualified[index];
    }
}

}  // namespace

int main(int argc, char** argv) {
    struct CliInput {
        std::string rundown_path;
        std::string output_path;
        std::string first_step;
        std::string last_step;
        std::string max_function_evaluations;
        std::string max_accepted_iterations;
        std::string iteration_directory;
        std::string jacobian_method;
        std::string parallel_workers;
        bool dry_run;
    };

    CLAP_BEGIN(CliInput)
    CLAP_ADD_USAGE("[OPTION]... RUNDOWN.json")
    CLAP_ADD_DESCRIPTION(
        "Execute a validated, step-by-step magnetic-equilibrium relaxation "
        "rundown.")
    CLAP_REGISTER_ARG(rundown_path)
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        output_path, "--output", "-o", "override output.path from the rundown")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(first_step, "--first-step",
                                          "begin at the named rundown step")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(last_step, "--last-step",
                                          "finish after the named rundown step")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        max_function_evaluations, "--max-function-evaluations",
        "override the per-phase residual-evaluation limit")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        max_accepted_iterations, "--max-accepted-iterations",
        "override the per-phase accepted-step limit")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        iteration_directory, "--iteration-directory",
        "override the directory for accepted inputs and equilibria")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(jacobian_method, "--jacobian-method",
                                          "override optimizer.jacobian.method")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(parallel_workers,
                                          "--parallel-workers",
                                          "override optimizer.jacobian.workers")
    CLAP_REGISTER_OPTION_WITH_DESCRIPTION(
        dry_run, "--dry-run", "validate and print the rundown without solving")
    CLAP_END(CliInput)

    CliInput cli{};
    try {
        CLAP<CliInput>::parse_input(cli, argc, argv);
    } catch (const std::exception& error) {
        std::cerr << error.what();
        return 2;
    }

    try {
        RelaxationRundown rundown =
            meow::config::read_relaxation_rundown(cli.rundown_path);
        if (!cli.output_path.empty()) { rundown.output.path = cli.output_path; }
        if (!cli.iteration_directory.empty()) {
            rundown.output.iteration_directory = cli.iteration_directory;
        }
        if (!cli.jacobian_method.empty()) {
            static_cast<void>(parse_jacobian_method(cli.jacobian_method));
            rundown.optimizer.jacobian.method = cli.jacobian_method;
        }
        if (!cli.parallel_workers.empty()) {
            rundown.optimizer.jacobian.workers = parse_positive_size(
                cli.parallel_workers, "--parallel-workers", true);
        }
        if (!cli.max_function_evaluations.empty()) {
            rundown.optimizer.max_function_evaluations = parse_positive_size(
                cli.max_function_evaluations, "--max-function-evaluations");
        }
        if (!cli.max_accepted_iterations.empty()) {
            rundown.optimizer.max_accepted_iterations = parse_positive_size(
                cli.max_accepted_iterations, "--max-accepted-iterations");
        }
        const JacobianMethod jacobian_method =
            parse_jacobian_method(rundown.optimizer.jacobian.method);
        const std::size_t parallel_workers = rundown.optimizer.jacobian.workers;

        const auto find_step = [&](const std::string& name,
                                   std::size_t fallback) {
            if (name.empty()) { return fallback; }
            for (std::size_t index = 0; index < rundown.steps.size(); ++index) {
                if (rundown.steps[index].name == name) { return index; }
            }
            throw std::invalid_argument("unknown rundown step '" + name + "'");
        };
        const std::size_t first_step = find_step(cli.first_step, 0);
        const std::size_t last_step =
            find_step(cli.last_step, rundown.steps.size() - 1);
        if (first_step > last_step) {
            throw std::invalid_argument(
                "--first-step must not follow --last-step");
        }
        std::cout << meow::config::summarize_relaxation_rundown(rundown);
        std::cout << "selected_steps=" << rundown.steps[first_step].name << ".."
                  << rundown.steps[last_step].name << '\n';
        if (cli.dry_run) { return 0; }

        cumes::SolverOptions validation_options;
#ifdef CUMES_USE_FLOAT
        validation_options.precision = cumes::PrecisionPolicy::MIXED_FLOAT;
#endif
        cumes::ParsedProblem parsed = cumes::read_problem_spec(
            rundown.equilibrium_input, validation_options);
        if (!parsed.report.ok()) {
            throw std::runtime_error(
                first_error(parsed.report, "input JSON mapping failed"));
        }
        cumes::ProblemSpec current = std::move(parsed.spec);
        if (rundown.initialization.axisymmetric_seed_amplitude.has_value() &&
            first_step == 0 &&
            cumes_meow_example::seed_landreman_qa_construction_boundary(
                current, *rundown.initialization.axisymmetric_seed_amplitude)) {
            std::cout << "seeded_axisymmetric_boundary amplitude="
                      << *rundown.initialization.axisymmetric_seed_amplitude
                      << '\n';
        }
        if (rundown.initialization.track_axis_predictor && first_step == 0) {
            cumes_meow_example::refresh_axis_predictor_from_boundary_centerline(
                current);
        }
        write_problem(rundown.output.path, current, validation_options);
        if (!rundown.output.iteration_directory.empty()) {
            std::filesystem::create_directories(
                rundown.output.iteration_directory);
        }
        std::optional<cumes::EquilibriumSnapshot> continuation_equilibrium;

        for (std::size_t step_index = first_step; step_index <= last_step;
             ++step_index) {
            const meow::config::RelaxationStepSpec& stage =
                rundown.steps[step_index];
            const int max_mode = stage.max_boundary_mode;
            set_stage_resolution(current, stage);
            const std::vector<double> qualified_tolerances =
                stage_tolerances(current);
            cumes_meow_example::StellaratorSymmetricBoundaryParameterization
                boundary(max_mode);
            const int stage_ns =
                static_cast<int>(current.stages.back().radial_surfaces);
            const int stage_mnmax = current.mpol * (current.ntor + 1);
            std::optional<cumes::EquilibriumSnapshot> phase_equilibrium =
                continuation_equilibrium;

            for (const meow::config::RelaxationPhaseSpec& phase :
                 stage.phases) {
                const std::string_view phase_name = phase.name;
                set_stage_tolerances(current, qualified_tolerances,
                                     phase.equilibrium_tolerance_floor);
                const meow::Vector initial = boundary.values(current);
                const cumes::EquilibriumSnapshot* stage_restart = nullptr;
                if (!phase.cold_start && phase_equilibrium.has_value() &&
                    phase_equilibrium->ns == stage_ns &&
                    phase_equilibrium->mnmax == stage_mnmax) {
                    stage_restart = &*phase_equilibrium;
                }
                LandremanResidual residual(
                    current, boundary, rundown, validation_options,
                    stage_restart,
                    rundown.optimizer.jacobian.warm_start_trials);

                meow::TrfOptions options;
                options.finite_difference_step =
                    rundown.optimizer.jacobian.relative_step;
                options.finite_difference_absolute_step =
                    rundown.optimizer.jacobian.absolute_step;
                options.ftol = rundown.optimizer.ftol;
                options.xtol = rundown.optimizer.xtol;
                options.gtol = rundown.optimizer.gtol;
                options.max_function_evaluations =
                    rundown.optimizer.max_function_evaluations;
                options.jacobian_refresh_interval =
                    rundown.optimizer.jacobian.refresh_interval;
                options.broyden_min_reduction_ratio =
                    rundown.optimizer.jacobian.minimum_reduction_ratio;
                options.broyden_max_secant_error =
                    rundown.optimizer.jacobian.maximum_secant_error;
                options.scale_from_jacobian =
                    rundown.optimizer.jacobian.scale_variables;
                if (rundown.optimizer.max_function_evaluations == 0 &&
                    rundown.optimizer.max_accepted_iterations != 0) {
                    const std::size_t evaluations_per_iteration =
                        2 * (boundary.size() + 2);
                    if (rundown.optimizer.max_accepted_iterations >
                        std::numeric_limits<std::size_t>::max() /
                            evaluations_per_iteration) {
                        throw std::overflow_error(
                            "accepted-iteration evaluation budget overflow");
                    }
                    options.max_function_evaluations =
                        rundown.optimizer.max_accepted_iterations *
                        evaluations_per_iteration;
                }
                options.verbose = 1;
                std::vector<double> accepted_objectives;
                options.callback = [&](const meow::Vector& x,
                                       const meow::IterationInfo& info) {
                    current = residual.accept(x);
                    write_problem(rundown.output.path, current,
                                  validation_options);
                    if (!rundown.output.iteration_directory.empty()) {
                        const std::string stem =
                            step_stem(rundown.output.iteration_directory,
                                      rundown.workflow, max_mode,
                                      info.iteration, phase_name);
                        write_problem(stem + "-input.json", current,
                                      validation_options);
                        residual.write_equilibrium(x,
                                                   stem + "-equilibrium.bin");
                    }
                    std::cout << "accepted mode=" << max_mode;
                    if (!phase_name.empty()) {
                        std::cout << " phase=" << phase_name;
                    }
                    std::cout << " iteration=" << info.iteration
                              << " objective=" << 2.0 * info.cost << '\n';
                    accepted_objectives.push_back(2.0 * info.cost);
                    bool phase_stagnated = false;
                    const std::size_t progress_window =
                        phase.stopping.progress_window;
                    if (progress_window != 0 &&
                        accepted_objectives.size() > progress_window) {
                        const double previous =
                            accepted_objectives[accepted_objectives.size() -
                                                progress_window - 1];
                        const double current_objective =
                            accepted_objectives.back();
                        const double relative_progress =
                            (previous - current_objective) /
                            std::max(std::abs(previous),
                                     std::numeric_limits<double>::min());
                        phase_stagnated =
                            info.iteration >=
                                phase.stopping.minimum_iterations &&
                            relative_progress <
                                phase.stopping.minimum_relative_progress;
                    }
                    const bool user_cap_reached =
                        rundown.optimizer.max_accepted_iterations != 0 &&
                        info.iteration >=
                            rundown.optimizer.max_accepted_iterations;
                    return !phase_stagnated && !user_cap_reached;
                };

                std::cout << "beginning max_mode=" << max_mode;
                if (!phase_name.empty()) {
                    std::cout << " phase=" << phase_name;
                }
                std::cout << " variables=" << boundary.size()
                          << " jacobian_method="
                          << jacobian_method_name(jacobian_method)
                          << " parallel_workers=" << parallel_workers
                          << " equilibrium_tolerance="
                          << current.stages.back().tolerance
                          << " stage_restart=" << (stage_restart != nullptr)
                          << '\n';
                if (!rundown.output.iteration_directory.empty()) {
                    const std::string stem =
                        step_stem(rundown.output.iteration_directory,
                                  rundown.workflow, max_mode, 0, phase_name);
                    write_problem(stem + "-input.json", current,
                                  validation_options);
                    residual.write_equilibrium(initial,
                                               stem + "-equilibrium.bin");
                }
                meow::TrfResult result;
                if (jacobian_method == JacobianMethod::GEOMETRY_RESTART_CHECK) {
                    const meow::Matrix cold =
                        residual.cold_finite_difference_jacobian(initial);
                    meow::Matrix restarted =
                        residual.geometry_restart_jacobian(initial);
                    double worst_relative_column_error = 0.0;
                    double mean_relative_column_error = 0.0;
                    for (Eigen::Index column = 0; column < cold.cols();
                         ++column) {
                        const double denominator =
                            std::max(cold.col(column).norm(),
                                     std::numeric_limits<double>::min());
                        const double error =
                            (restarted.col(column) - cold.col(column)).norm() /
                            denominator;
                        worst_relative_column_error =
                            std::max(worst_relative_column_error, error);
                        mean_relative_column_error += error;
                        std::cout << "geometry_restart_column=" << column
                                  << " relative_error=" << error
                                  << " cold_norm=" << denominator << '\n';
                    }
                    mean_relative_column_error /=
                        static_cast<double>(cold.cols());
                    std::cout << "geometry_restart_comparison "
                              << "relative_frobenius_error="
                              << (restarted - cold).norm() / cold.norm()
                              << " worst_relative_column_error="
                              << worst_relative_column_error
                              << " mean_relative_column_error="
                              << mean_relative_column_error << '\n';
                    result.x = initial;
                    result.residual = residual(initial);
                    result.jacobian = std::move(restarted);
                    result.gradient =
                        result.jacobian.transpose() * result.residual;
                    result.cost = 0.5 * result.residual.squaredNorm();
                    result.status = meow::TrfStatus::USER_STOPPED;
                    result.message =
                        "Geometry-restart Jacobian comparison completed.";
                } else if (jacobian_method ==
                           JacobianMethod::PARALLEL_FINITE_DIFFERENCE_CHECK) {
                    const meow::Matrix cold =
                        residual.cold_finite_difference_jacobian(initial);
                    meow::Matrix concurrent =
                        residual.configured_parallel_finite_difference_jacobian(
                            initial, parallel_workers);
                    double worst_relative_column_error = 0.0;
                    double mean_relative_column_error = 0.0;
                    for (Eigen::Index column = 0; column < cold.cols();
                         ++column) {
                        const double denominator =
                            std::max(cold.col(column).norm(),
                                     std::numeric_limits<double>::min());
                        const double error =
                            (concurrent.col(column) - cold.col(column)).norm() /
                            denominator;
                        worst_relative_column_error =
                            std::max(worst_relative_column_error, error);
                        mean_relative_column_error += error;
                        std::cout << "parallel_column=" << column
                                  << " relative_error=" << error
                                  << " cold_norm=" << denominator << '\n';
                    }
                    mean_relative_column_error /=
                        static_cast<double>(cold.cols());
                    std::cout << "parallel_jacobian_comparison "
                              << "relative_frobenius_error="
                              << (concurrent - cold).norm() / cold.norm()
                              << " worst_relative_column_error="
                              << worst_relative_column_error
                              << " mean_relative_column_error="
                              << mean_relative_column_error << '\n';
                    result.x = initial;
                    result.residual = residual(initial);
                    result.jacobian = std::move(concurrent);
                    result.gradient =
                        result.jacobian.transpose() * result.residual;
                    result.cost = 0.5 * result.residual.squaredNorm();
                    result.status = meow::TrfStatus::USER_STOPPED;
                    result.message = "Parallel Jacobian comparison completed.";
                } else if (jacobian_method ==
                           JacobianMethod::
                               RELAXED_PARALLEL_FINITE_DIFFERENCE_CHECK) {
                    const meow::Matrix cold =
                        residual.cold_finite_difference_jacobian(initial);
                    meow::Matrix relaxed =
                        residual.relaxed_parallel_finite_difference_jacobian(
                            initial, parallel_workers);
                    double worst_relative_column_error = 0.0;
                    double mean_relative_column_error = 0.0;
                    for (Eigen::Index column = 0; column < cold.cols();
                         ++column) {
                        const double denominator =
                            std::max(cold.col(column).norm(),
                                     std::numeric_limits<double>::min());
                        const double error =
                            (relaxed.col(column) - cold.col(column)).norm() /
                            denominator;
                        worst_relative_column_error =
                            std::max(worst_relative_column_error, error);
                        mean_relative_column_error += error;
                        std::cout << "relaxed_parallel_column=" << column
                                  << " relative_error=" << error
                                  << " cold_norm=" << denominator << '\n';
                    }
                    mean_relative_column_error /=
                        static_cast<double>(cold.cols());
                    std::cout << "relaxed_parallel_jacobian_comparison "
                              << "relative_frobenius_error="
                              << (relaxed - cold).norm() / cold.norm()
                              << " worst_relative_column_error="
                              << worst_relative_column_error
                              << " mean_relative_column_error="
                              << mean_relative_column_error << '\n';
                    result.x = initial;
                    result.residual = residual(initial);
                    result.jacobian = std::move(relaxed);
                    result.gradient =
                        result.jacobian.transpose() * result.residual;
                    result.cost = 0.5 * result.residual.squaredNorm();
                    result.status = meow::TrfStatus::USER_STOPPED;
                    result.message =
                        "Relaxed parallel Jacobian comparison completed.";
                } else if (jacobian_method ==
                           JacobianMethod::PARALLEL_WORKER_COUNT_CHECK) {
                    const meow::Matrix cold =
                        residual.cold_finite_difference_jacobian(initial);
                    const meow::Matrix two_workers =
                        residual.parallel_finite_difference_jacobian(initial);
                    const meow::Matrix four_workers =
                        residual.configured_parallel_finite_difference_jacobian(
                            initial, 4);
                    meow::Matrix eight_workers =
                        residual.eight_worker_finite_difference_jacobian(
                            initial);
                    std::cout
                        << "parallel_worker_comparison "
                        << "two_relative_frobenius_error="
                        << (two_workers - cold).norm() / cold.norm()
                        << " four_relative_frobenius_error="
                        << (four_workers - cold).norm() / cold.norm()
                        << " two_four_relative_frobenius_difference="
                        << (four_workers - two_workers).norm() / cold.norm()
                        << " eight_relative_frobenius_error="
                        << (eight_workers - cold).norm() / cold.norm()
                        << " four_eight_relative_frobenius_difference="
                        << (eight_workers - four_workers).norm() / cold.norm()
                        << '\n';
                    result.x = initial;
                    result.residual = residual(initial);
                    result.jacobian = std::move(eight_workers);
                    result.gradient =
                        result.jacobian.transpose() * result.residual;
                    result.cost = 0.5 * result.residual.squaredNorm();
                    result.status = meow::TrfStatus::USER_STOPPED;
                    result.message =
                        "Parallel worker-count comparison completed.";
                } else if (jacobian_method == JacobianMethod::ANALYTIC) {
                    result = meow::trf_least_squares(
                        std::ref(residual), initial, options,
                        [&](const meow::Vector& x) {
                            return residual.jacobian(x);
                        });
                } else if (jacobian_method ==
                           JacobianMethod::HOT_FINITE_DIFFERENCE) {
                    result = meow::trf_least_squares(
                        std::ref(residual), initial, options,
                        [&](const meow::Vector& x) {
                            return residual.hot_restart_jacobian(x);
                        });
                } else if (jacobian_method ==
                           JacobianMethod::WARM_FINITE_DIFFERENCE) {
                    result = meow::trf_least_squares(
                        std::ref(residual), initial, options,
                        [&](const meow::Vector& x) {
                            return residual.warm_restart_jacobian(x);
                        });
                } else if (jacobian_method ==
                           JacobianMethod::GEOMETRY_RESTART_FINITE_DIFFERENCE) {
                    result = meow::trf_least_squares(
                        std::ref(residual), initial, options,
                        [&](const meow::Vector& x) {
                            return residual.geometry_restart_jacobian(x);
                        });
                } else if (jacobian_method ==
                               JacobianMethod::PARALLEL_FINITE_DIFFERENCE ||
                           jacobian_method ==
                               JacobianMethod::PARALLEL_AGGRESSIVE_BROYDEN) {
                    result = meow::trf_least_squares(
                        std::ref(residual), initial, options,
                        [&](const meow::Vector& x) {
                            return residual.parallel_finite_difference_jacobian(
                                x, std::nullopt, parallel_workers);
                        });
                } else if (jacobian_method ==
                               JacobianMethod::FOUR_WORKER_AGGRESSIVE_BROYDEN ||
                           jacobian_method ==
                               JacobianMethod::EXTENDED_FOUR_WORKER_BROYDEN) {
                    result = meow::trf_least_squares(
                        std::ref(residual), initial, options,
                        [&](const meow::Vector& x) {
                            return residual
                                .configured_parallel_finite_difference_jacobian(
                                    x, parallel_workers);
                        });
                } else {
                    result = meow::trf_least_squares(std::ref(residual),
                                                     initial, options);
                }
                current = residual.accept(result.x);
                phase_equilibrium = residual.equilibrium();
                write_problem(rundown.output.path, current, validation_options);
                std::cout << "finished max_mode=" << max_mode;
                if (!phase_name.empty()) {
                    std::cout << " phase=" << phase_name;
                }
                std::cout << " status=" << result.message
                          << " objective=" << 2.0 * result.cost
                          << " evaluations=" << result.function_evaluations
                          << " iterations=" << result.iterations
                          << " jacobian_updates=" << result.jacobian_updates
                          << " equilibrium_evaluations="
                          << residual.equilibrium_evaluations()
                          << " nonlinear_iterations="
                          << residual.total_nonlinear_iterations()
                          << " analytic_jacobians="
                          << residual.analytic_jacobian_evaluations()
                          << " hot_restart_jacobians="
                          << residual.hot_restart_jacobian_evaluations()
                          << " parallel_jacobians="
                          << residual.parallel_jacobian_evaluations()
                          << " restart_solve_attempts="
                          << residual.restart_solve_attempts()
                          << " restart_solve_failures="
                          << residual.restart_solve_failures()
                          << " linear_iterations="
                          << residual.total_linear_iterations() << '\n';
            }
            continuation_equilibrium = std::move(phase_equilibrium);
            write_problem(checkpoint_path(rundown.output.path, rundown.workflow,
                                          max_mode),
                          current, validation_options);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cumes-landreman-optimize: " << error.what() << '\n';
        return 1;
    }
}
