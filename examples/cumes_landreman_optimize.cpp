// Reproduce the analytic construction or final refinement from Landreman-Paul.
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <cumes/config/json_reader.hpp>
#include <cumes/config/json_writer.hpp>
#include <cumes/config/validated_problem.hpp>
#include <cumes/io/output_spec.hpp>
#include <cumes/io/writer.hpp>
#include <cumes/solver/equilibrium_linearization.hpp>
#include <cumes/solver/equilibrium_solver.hpp>
#include <meow/cumes/boundary_parameterization.hpp>
#include <meow/cumes/landreman_workflow.hpp>
#include <meow/cumes/quasisymmetry_target.hpp>
#include <meow/cumes/quasisymmetry_target_jvp.hpp>
#include <meow/trf.hpp>

namespace {

using cumes_meow_example::LandremanCase;
using cumes_meow_example::LandremanSelection;
using cumes_meow_example::LandremanWorkflow;

enum class JacobianMethod { ANALYTIC, FINITE_DIFFERENCE };

JacobianMethod parse_jacobian_method(const std::string& name) {
    if (name == "analytic") { return JacobianMethod::ANALYTIC; }
    if (name == "finite-difference") {
        return JacobianMethod::FINITE_DIFFERENCE;
    }
    throw std::invalid_argument(
        "JACOBIAN_METHOD must be analytic or finite-difference");
}

const char* jacobian_method_name(JacobianMethod method) {
    return method == JacobianMethod::ANALYTIC ? "analytic"
                                              : "finite-difference";
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
        LandremanSelection selection,
        cumes::SolverOptions validation_options)
        : baseline_(std::move(baseline)),
          boundary_(std::move(boundary)),
          selection_(selection),
          validation_options_(validation_options) {}

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
        if (selection_.selected_case == LandremanCase::QH) {
            request.radial_transfer = cumes::RadialTransferPolicy::CATMULL_ROM;
        }
        cumes::SolveOutcome solved = solver_.solve(validated.value(), request);
        ++evaluation_count_;
        total_nonlinear_iterations_ += solved.total_iterations;
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

        auto spec = target_spec(solved.report.input_params.nfp);
        const double target_aspect =
            cumes_meow_example::landreman_target_aspect(
                selection_.selected_case);
        const cumes_meow_example::CompositeQuasisymmetryTarget target =
            selection_.selected_case == LandremanCase::QA
                ? cumes_meow_example::calculate_qa_target(
                      solved.equilibrium, solved.profiles,
                      solved.report.input_params, spec, target_aspect,
                      cumes_meow_example::landreman_targets_mean_iota(
                          selection_)
                          ? std::optional<double>(0.42)
                          : std::nullopt)
                : cumes_meow_example::calculate_qh_target(
                      solved.equilibrium, solved.profiles,
                      solved.report.input_params, spec, target_aspect);

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
        std::size_t total_linear_iterations = 0;
        double worst_linear_residual = 0.0;
        for (std::size_t column = 0; column < boundary_.size(); ++column) {
            const cumes::BoundaryTangent boundary_tangent =
                boundary_.tangent(validated.value(), column);
            const cumes::SpectralTangentSolve spectral =
                linearization.solve_boundary_tangent(boundary_tangent,
                                                     tangent_options);
            if (!spectral.converged) {
                throw std::runtime_error(
                    "equilibrium tangent solve failed for " +
                    boundary_.name(column) + " after " +
                    std::to_string(spectral.iterations) +
                    " iterations; initial_residual=" +
                    std::to_string(spectral.initial_residual) +
                    ", final_residual=" +
                    std::to_string(spectral.final_residual));
            }
            total_linear_iterations +=
                static_cast<std::size_t>(spectral.iterations);
            worst_linear_residual =
                std::max(worst_linear_residual, spectral.final_residual);
            const cumes::EquilibriumTangent tangent =
                linearization.materialize_tangent(spectral.state_tangent,
                                                  cached_outcome_->equilibrium,
                                                  cached_outcome_->profiles);
            const std::vector<double> target_tangent =
                selection_.selected_case == LandremanCase::QA
                    ? cumes_meow_example::calculate_qa_target_jvp(
                          cached_outcome_->equilibrium,
                          cached_outcome_->profiles, tangent,
                          cached_outcome_->report.input_params, spec,
                          cumes_meow_example::landreman_targets_mean_iota(
                              selection_))
                    : cumes_meow_example::calculate_qh_target_jvp(
                          cached_outcome_->equilibrium,
                          cached_outcome_->profiles, tangent,
                          cached_outcome_->report.input_params, spec);
            if (target_tangent.size() !=
                static_cast<std::size_t>(result.rows())) {
                throw std::runtime_error(
                    "analytic target tangent has the wrong length");
            }
            for (std::size_t row = 0; row < target_tangent.size(); ++row) {
                result(static_cast<Eigen::Index>(row),
                       static_cast<Eigen::Index>(column)) = target_tangent[row];
            }
        }
        std::cout << "analytic_jacobian columns=" << result.cols()
                  << " residuals=" << result.rows()
                  << " linear_iterations=" << total_linear_iterations
                  << " worst_linear_residual=" << worst_linear_residual << '\n';
        ++analytic_jacobian_evaluations_;
        total_linear_iterations_ += total_linear_iterations;
        return result;
    }

    std::size_t equilibrium_evaluations() const { return evaluation_count_; }

    std::size_t total_nonlinear_iterations() const {
        return total_nonlinear_iterations_;
    }

    std::size_t analytic_jacobian_evaluations() const {
        return analytic_jacobian_evaluations_;
    }

    std::size_t total_linear_iterations() const {
        return total_linear_iterations_;
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
        if (cumes_meow_example::landreman_refreshes_axis_predictor(
                selection_)) {
            cumes_meow_example::refresh_axis_predictor_from_equilibrium(
                baseline_, cached_outcome_->equilibrium);
        }
        return baseline_;
    }

   private:
    cumes::ProblemSpec trial_problem(const meow::Vector& x) const {
        cumes::ProblemSpec problem = boundary_.apply(baseline_, x);
        if (cumes_meow_example::landreman_refreshes_axis_predictor(
                selection_)) {
            cumes_meow_example::track_axis_predictor_from_accepted_boundary(
                problem, baseline_);
        }
        return problem;
    }

    cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec target_spec(
        int nfp) const {
        cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec spec;
        spec.helicity_m = 1;
        spec.helicity_n =
            selection_.selected_case == LandremanCase::QA ? 0 : -nfp;
        const double final_weight =
            cumes_meow_example::landreman_final_surface_weight(selection_);
        for (int index = 0; index <= 10; ++index) {
            spec.normalized_toroidal_flux_surfaces.push_back(index / 10.0);
            spec.surface_weights.push_back(1.0 +
                                           (final_weight - 1.0) * index / 10.0);
        }
        return spec;
    }

    cumes::ProblemSpec baseline_;
    cumes_meow_example::StellaratorSymmetricBoundaryParameterization boundary_;
    LandremanSelection selection_;
    cumes::SolverOptions validation_options_;
    cumes::EquilibriumSolver solver_;
    std::size_t evaluation_count_ = 0;
    std::size_t total_nonlinear_iterations_ = 0;
    std::size_t analytic_jacobian_evaluations_ = 0;
    std::size_t total_linear_iterations_ = 0;
    std::optional<meow::Vector> cached_x_;
    meow::Vector cached_residual_;
    std::optional<cumes::SolveOutcome> cached_outcome_;
};

std::string step_stem(const std::string& directory,
                      LandremanWorkflow workflow,
                      int max_mode,
                      std::size_t iteration) {
    std::ostringstream name;
    if (workflow == LandremanWorkflow::CONSTRUCTION) {
        name << cumes_meow_example::landreman_workflow_name(workflow) << "-";
    }
    name << "mode" << max_mode << "_step_" << std::setw(4) << std::setfill('0')
         << iteration;
    return (std::filesystem::path(directory) / name.str()).string();
}

std::string checkpoint_path(const std::string& output_path,
                            LandremanWorkflow workflow,
                            int max_mode) {
    if (workflow == LandremanWorkflow::REFINEMENT) {
        return output_path + ".mode" + std::to_string(max_mode) + ".json";
    }
    return output_path + ".construction.mode" + std::to_string(max_mode) +
           ".json";
}

void set_stage_resolution(cumes::ProblemSpec& problem,
                          const cumes_meow_example::LandremanStage& stage) {
    if (stage.mpol == 0 || stage.ntor == 0) { return; }
    const std::vector<double> source_raxis = problem.raxis_c;
    const std::vector<double> source_zaxis = problem.zaxis_s;
    problem.mpol = stage.mpol;
    problem.ntor = stage.ntor;
    for (cumes::StageRequest& request : problem.stages) {
        request.max_iterations =
            std::max(request.max_iterations,
                     static_cast<std::size_t>(stage.minimum_niter));
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

void print_usage() {
    std::cerr
        << "usage: cumes_landreman_optimize INPUT.json "
           "qa|qh|qa-construction|qh-construction OUTPUT.json "
           "[MAX_FUNCTION_EVALUATIONS_PER_STAGE [FIRST_MODE [LAST_MODE "
           "[MAX_ACCEPTED_ITERATIONS [ITERATION_DIRECTORY "
           "[JACOBIAN_METHOD]]]]]]\n"
        << "qa/qh run the archived mode-4/mode-5 final refinement. The "
           "*-construction cases start from the analytic boundary and run "
           "modes 1-4 (QA) or 1-5 (QH), with the QA iota target enabled. "
           "Zero or an omitted evaluation limit uses meow's 100*n default. "
           "FIRST_MODE/LAST_MODE select an ordered subset of the chosen "
           "workflow. ITERATION_DIRECTORY stores the "
           "input and native equilibrium for step 0 and every accepted "
           "iteration. JACOBIAN_METHOD is analytic (default) or "
           "finite-difference.\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4 || argc > 10) {
        print_usage();
        return 2;
    }
    try {
        const LandremanSelection selection =
            cumes_meow_example::parse_landreman_selection(argv[2]);
        const auto stages = cumes_meow_example::landreman_stages(selection);
        const std::string output_path = argv[3];
        const std::size_t max_evaluations =
            argc >= 5 ? std::stoull(argv[4]) : 0;
        const int first_mode =
            argc >= 6 ? std::stoi(argv[5]) : stages.front().max_mode;
        const int last_mode =
            argc >= 7 ? std::stoi(argv[6]) : stages.back().max_mode;
        const std::size_t max_accepted_iterations =
            argc >= 8 ? std::stoull(argv[7]) : 0;
        const std::string iteration_directory = argc >= 9 ? argv[8] : "";
        const JacobianMethod jacobian_method =
            parse_jacobian_method(argc >= 10 ? argv[9] : "analytic");
        const auto has_mode = [&](int mode) {
            for (const auto& stage : stages) {
                if (stage.max_mode == mode) { return true; }
            }
            return false;
        };
        if (!has_mode(first_mode) || !has_mode(last_mode) ||
            first_mode > last_mode) {
            throw std::invalid_argument(
                "FIRST_MODE and LAST_MODE must select an ordered subset of "
                "the chosen workflow");
        }

        cumes::SolverOptions validation_options;
#ifdef CUMES_USE_FLOAT
        validation_options.precision = cumes::PrecisionPolicy::MIXED_FLOAT;
#endif
        cumes::ParsedProblem parsed =
            cumes::read_problem_spec(argv[1], validation_options);
        if (!parsed.report.ok()) {
            throw std::runtime_error(
                first_error(parsed.report, "input JSON mapping failed"));
        }
        cumes::ProblemSpec current = std::move(parsed.spec);
        if (selection.workflow == LandremanWorkflow::CONSTRUCTION &&
            selection.selected_case == LandremanCase::QA &&
            first_mode == stages.front().max_mode &&
            cumes_meow_example::seed_landreman_qa_construction_boundary(
                current)) {
            std::cout << "seeded_axisymmetric_qa_boundary amplitude=0.0001\n";
        }
        if (cumes_meow_example::landreman_refreshes_axis_predictor(selection) &&
            first_mode == stages.front().max_mode) {
            cumes_meow_example::refresh_axis_predictor_from_boundary_centerline(
                current);
        }
        write_problem(output_path, current, validation_options);
        if (!iteration_directory.empty()) {
            std::filesystem::create_directories(iteration_directory);
        }

        for (const auto& stage : stages) {
            const int max_mode = stage.max_mode;
            if (max_mode < first_mode || max_mode > last_mode) { continue; }
            set_stage_resolution(current, stage);
            cumes_meow_example::StellaratorSymmetricBoundaryParameterization
                boundary(max_mode);
            const meow::Vector initial = boundary.values(current);
            LandremanResidual residual(current, boundary, selection,
                                       validation_options);

            meow::TrfOptions options;
            const auto finite_difference =
                cumes_meow_example::landreman_finite_difference_policy(
                    selection);
            options.finite_difference_step = finite_difference.relative_step;
            options.finite_difference_absolute_step =
                finite_difference.absolute_step;
            options.max_function_evaluations = max_evaluations;
            if (max_evaluations == 0 && max_accepted_iterations != 0) {
                const std::size_t evaluations_per_iteration =
                    2 * (boundary.size() + 2);
                if (max_accepted_iterations >
                    std::numeric_limits<std::size_t>::max() /
                        evaluations_per_iteration) {
                    throw std::overflow_error(
                        "accepted-iteration evaluation budget overflow");
                }
                options.max_function_evaluations =
                    max_accepted_iterations * evaluations_per_iteration;
            }
            options.verbose = 1;
            options.callback = [&](const meow::Vector& x,
                                   const meow::IterationInfo& info) {
                current = residual.accept(x);
                write_problem(output_path, current, validation_options);
                if (!iteration_directory.empty()) {
                    const std::string stem =
                        step_stem(iteration_directory, selection.workflow,
                                  max_mode, info.iteration);
                    write_problem(stem + "-input.json", current,
                                  validation_options);
                    residual.write_equilibrium(x, stem + "-equilibrium.bin");
                }
                std::cout << "accepted mode=" << max_mode
                          << " iteration=" << info.iteration
                          << " objective=" << 2.0 * info.cost << '\n';
                return max_accepted_iterations == 0 ||
                       info.iteration < max_accepted_iterations;
            };

            std::cout << "beginning max_mode=" << max_mode
                      << " variables=" << boundary.size() << " jacobian_method="
                      << jacobian_method_name(jacobian_method) << '\n';
            if (!iteration_directory.empty()) {
                const std::string stem = step_stem(
                    iteration_directory, selection.workflow, max_mode, 0);
                write_problem(stem + "-input.json", current,
                              validation_options);
                residual.write_equilibrium(initial, stem + "-equilibrium.bin");
            }
            meow::TrfResult result;
            if (jacobian_method == JacobianMethod::ANALYTIC) {
                result = meow::trf_least_squares(
                    std::ref(residual), initial, options,
                    [&](const meow::Vector& x) {
                        return residual.jacobian(x);
                    });
            } else {
                result = meow::trf_least_squares(std::ref(residual), initial,
                                                 options);
            }
            current = residual.accept(result.x);
            write_problem(output_path, current, validation_options);
            write_problem(
                checkpoint_path(output_path, selection.workflow, max_mode),
                current, validation_options);
            std::cout << "finished max_mode=" << max_mode
                      << " status=" << result.message
                      << " objective=" << 2.0 * result.cost
                      << " evaluations=" << result.function_evaluations
                      << " iterations=" << result.iterations
                      << " equilibrium_evaluations="
                      << residual.equilibrium_evaluations()
                      << " nonlinear_iterations="
                      << residual.total_nonlinear_iterations()
                      << " analytic_jacobians="
                      << residual.analytic_jacobian_evaluations()
                      << " linear_iterations="
                      << residual.total_linear_iterations() << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cumes-landreman-optimize: " << error.what() << '\n';
        return 1;
    }
}
