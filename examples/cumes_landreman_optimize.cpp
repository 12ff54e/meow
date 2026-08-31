// Reproduce the final two-stage QA/QH boundary refinement from Landreman-Paul.
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include <cumes/config/json_reader.hpp>
#include <cumes/config/json_writer.hpp>
#include <cumes/config/validated_problem.hpp>
#include <cumes/solver/equilibrium_solver.hpp>
#include <meow/cumes/boundary_parameterization.hpp>
#include <meow/cumes/quasisymmetry_target.hpp>
#include <meow/trf.hpp>

namespace {

enum class LandremanCase { QA, QH };

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
        LandremanCase selected_case,
        cumes::SolverOptions validation_options)
        : baseline_(std::move(baseline)),
          boundary_(std::move(boundary)),
          selected_case_(selected_case),
          validation_options_(validation_options) {}

    meow::Vector operator()(const meow::Vector& x) {
        if (cached_x_.has_value() && cached_x_->size() == x.size() &&
            (cached_x_->array() == x.array()).all()) {
            return cached_residual_;
        }

        cumes::ProblemSpec problem = boundary_.apply(baseline_, x);
        cumes::ValidationResult validated =
            cumes::validate(std::move(problem), validation_options_);
        if (!validated.has_value()) {
            throw std::runtime_error(first_error(
                validated.error(), "cuMES boundary validation failed"));
        }

        cumes::SolveRequest request;
        if (selected_case_ == LandremanCase::QH) {
            request.radial_transfer = cumes::RadialTransferPolicy::CATMULL_ROM;
        }
        cumes::SolveOutcome solved = solver_.solve(validated.value(), request);
        ++evaluation_count_;
        if (!solved.converged || !solved.has_complete_equilibrium()) {
            throw std::runtime_error(
                "cuMES equilibrium failed at objective evaluation " +
                std::to_string(evaluation_count_) + ", stage " +
                std::to_string(solved.failed_stage) + ", residuals (" +
                std::to_string(solved.fsqr) + ", " +
                std::to_string(solved.fsqz) + ", " +
                std::to_string(solved.fsql) + ")");
        }

        auto spec = target_spec(solved.report.input_params.nfp);
        const double target_aspect =
            selected_case_ == LandremanCase::QA ? 6.0 : 8.0;
        const cumes_meow_example::CompositeQuasisymmetryTarget target =
            selected_case_ == LandremanCase::QA
                ? cumes_meow_example::calculate_qa_target(
                      solved.equilibrium, solved.profiles,
                      solved.report.input_params, spec, target_aspect,
                      std::nullopt)
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
        return residual;
    }

   private:
    cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec target_spec(
        int nfp) const {
        cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec spec;
        spec.helicity_m = 1;
        spec.helicity_n = selected_case_ == LandremanCase::QA ? 0 : -nfp;
        const double final_weight =
            selected_case_ == LandremanCase::QA ? 30.0 : 2.0;
        for (int index = 0; index <= 10; ++index) {
            spec.normalized_toroidal_flux_surfaces.push_back(index / 10.0);
            spec.surface_weights.push_back(1.0 +
                                           (final_weight - 1.0) * index / 10.0);
        }
        return spec;
    }

    cumes::ProblemSpec baseline_;
    cumes_meow_example::StellaratorSymmetricBoundaryParameterization boundary_;
    LandremanCase selected_case_;
    cumes::SolverOptions validation_options_;
    cumes::EquilibriumSolver solver_;
    std::size_t evaluation_count_ = 0;
    std::optional<meow::Vector> cached_x_;
    meow::Vector cached_residual_;
};

void print_usage() {
    std::cerr
        << "usage: cumes_landreman_optimize INPUT.json qa|qh OUTPUT.json "
           "[MAX_FUNCTION_EVALUATIONS_PER_STAGE [FIRST_MODE [LAST_MODE]]]\n"
        << "The archived final refinement is run first through boundary mode "
           "4, then through mode 5. Zero or an omitted evaluation limit uses "
           "meow's 100*n default. FIRST_MODE/LAST_MODE can select 4 or 5 for "
           "a checkpointed partial run.\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4 || argc > 7) {
        print_usage();
        return 2;
    }
    try {
        const std::string case_name = argv[2];
        const LandremanCase selected_case =
            case_name == "qa" ? LandremanCase::QA
            : case_name == "qh"
                ? LandremanCase::QH
                : throw std::invalid_argument("case must be qa or qh");
        const std::string output_path = argv[3];
        const std::size_t max_evaluations =
            argc >= 5 ? std::stoull(argv[4]) : 0;
        const int first_mode = argc >= 6 ? std::stoi(argv[5]) : 4;
        const int last_mode = argc >= 7 ? std::stoi(argv[6]) : 5;
        if (first_mode < 4 || last_mode > 5 || first_mode > last_mode) {
            throw std::invalid_argument(
                "FIRST_MODE and LAST_MODE must select an ordered subset of "
                "{4,5}");
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
        write_problem(output_path, current, validation_options);

        for (int max_mode = first_mode; max_mode <= last_mode; ++max_mode) {
            cumes_meow_example::StellaratorSymmetricBoundaryParameterization
                boundary(max_mode);
            const meow::Vector initial = boundary.values(current);
            LandremanResidual residual(current, boundary, selected_case,
                                       validation_options);

            meow::TrfOptions options;
            options.finite_difference_step = 1.0e-5;
            options.finite_difference_absolute_step = 1.0e-9;
            options.max_function_evaluations = max_evaluations;
            options.verbose = 1;
            options.callback = [&](const meow::Vector& x,
                                   const meow::IterationInfo& info) {
                write_problem(output_path, boundary.apply(current, x),
                              validation_options);
                std::cout << "accepted mode=" << max_mode
                          << " iteration=" << info.iteration
                          << " objective=" << 2.0 * info.cost << '\n';
                return true;
            };

            std::cout << "beginning max_mode=" << max_mode
                      << " variables=" << boundary.size() << '\n';
            const meow::TrfResult result =
                meow::trf_least_squares(std::ref(residual), initial, options);
            current = boundary.apply(current, result.x);
            write_problem(output_path, current, validation_options);
            write_problem(
                output_path + ".mode" + std::to_string(max_mode) + ".json",
                current, validation_options);
            std::cout << "finished max_mode=" << max_mode
                      << " status=" << result.message
                      << " objective=" << 2.0 * result.cost
                      << " evaluations=" << result.function_evaluations
                      << " iterations=" << result.iterations << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cumes-landreman-optimize: " << error.what() << '\n';
        return 1;
    }
}
