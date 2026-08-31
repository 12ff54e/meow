#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include <cumes/config/json_reader.hpp>
#include <cumes/config/validated_problem.hpp>
#include <cumes/solver/equilibrium_solver.hpp>
#include <meow/cumes/quasisymmetry_target.hpp>

namespace {

void print_values(const std::vector<double>& values) {
    std::cout << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) std::cout << ',';
        std::cout << values[index];
    }
    std::cout << ']';
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: cumes_landreman_evaluate INPUT.json qa|qh\n";
        return 2;
    }
    try {
        const std::string target_name = argv[2];
        if (target_name != "qa" && target_name != "qh") {
            throw std::invalid_argument("case must be qa or qh");
        }
        cumes::SolverOptions validation_options;
#ifdef CUMES_USE_FLOAT
        validation_options.precision = cumes::PrecisionPolicy::MIXED_FLOAT;
#endif
        cumes::ParsedProblem parsed =
            cumes::read_problem_spec(argv[1], validation_options);
        if (!parsed.report.ok()) {
            const auto errors = parsed.report.errors();
            throw std::runtime_error(
                errors.empty() ? "input JSON mapping failed" : errors.front());
        }
        cumes::ValidationResult validated =
            cumes::validate(std::move(parsed.spec), validation_options);
        if (!validated.has_value()) {
            const auto errors = validated.error().errors();
            throw std::runtime_error(errors.empty() ? "input validation failed"
                                                    : errors.front());
        }

        cumes::EquilibriumSolver solver;
        const cumes::SolveOutcome solved = solver.solve(validated.value());
        if (!solved.converged || !solved.has_complete_equilibrium()) {
            std::cerr << "cuMES equilibrium did not converge: failed_stage="
                      << solved.failed_stage << " fsqr=" << solved.fsqr
                      << " fsqz=" << solved.fsqz << " fsql=" << solved.fsql
                      << '\n';
            return 1;
        }

        cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec spec;
        spec.helicity_m = 1;
        spec.helicity_n =
            target_name == "qa" ? 0 : -solved.report.input_params.nfp;
        for (int index = 0; index <= 10; ++index) {
            spec.normalized_toroidal_flux_surfaces.push_back(index / 10.0);
        }
        const double weight_end = target_name == "qa" ? 30.0 : 2.0;
        for (int index = 0; index <= 10; ++index) {
            spec.surface_weights.push_back(1.0 +
                                           (weight_end - 1.0) * index / 10.0);
        }
        const auto weighted =
            cumes_meow_example::calculate_quasisymmetry_target(
                solved.equilibrium, solved.profiles,
                solved.report.input_params.nfp, spec);
        spec.surface_weights.clear();
        const auto unweighted =
            cumes_meow_example::calculate_quasisymmetry_target(
                solved.equilibrium, solved.profiles,
                solved.report.input_params.nfp, spec);
        const auto plasma_size = cumes_meow_example::calculate_plasma_size(
            solved.equilibrium, solved.report.input_params);
        const double mean_iota =
            cumes_meow_example::integrate_iota(solved.profiles);
        const double target_aspect = target_name == "qa" ? 6.0 : 8.0;
        const double aspect_residual = plasma_size.aspect_ratio - target_aspect;

        std::cout << std::setprecision(17);
        std::cout << "{\n"
                  << "  \"case\": \"" << target_name << "\",\n"
                  << "  \"converged\": true,\n"
                  << "  \"iterations\": " << solved.iterations << ",\n"
                  << "  \"fsqr\": " << solved.fsqr << ",\n"
                  << "  \"fsqz\": " << solved.fsqz << ",\n"
                  << "  \"fsql\": " << solved.fsql << ",\n"
                  << "  \"aspect_ratio\": " << plasma_size.aspect_ratio << ",\n"
                  << "  \"mean_iota\": " << mean_iota << ",\n"
                  << "  \"weighted_qs\": " << weighted.value << ",\n"
                  << "  \"unweighted_qs\": " << unweighted.value << ",\n"
                  << "  \"objective\": "
                  << weighted.value + aspect_residual * aspect_residual << ",\n"
                  << "  \"weighted_qs_profile\": ";
        print_values(weighted.surface_values);
        std::cout << ",\n  \"unweighted_qs_profile\": ";
        print_values(unweighted.surface_values);
        std::cout << "\n}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cumes-landreman-evaluate: " << error.what() << '\n';
        return 1;
    }
}
