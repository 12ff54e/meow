#include "meow/config/relaxation_rundown.hpp"

#include "JsonParser.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace meow::config {
namespace {

using JsonValue = meow::json::Value;

[[noreturn]] void fail(std::string_view path, const std::string& message) {
    throw std::invalid_argument(std::string(path) + ": " + message);
}

const JsonValue::object_container_type& require_object(const JsonValue& value,
                                                       std::string_view path) {
    if (!value.is_object()) { fail(path, "expected an object"); }
    return value.as_object();
}

const JsonValue::array_container_type& require_array(const JsonValue& value,
                                                     std::string_view path) {
    if (!value.is_array()) { fail(path, "expected an array"); }
    return value.as_array();
}

void reject_unknown_keys(const JsonValue& value,
                         std::string_view path,
                         const std::unordered_set<std::string>& allowed) {
    const auto& object = require_object(value, path);
    for (const auto& [key, ignored] : object) {
        static_cast<void>(ignored);
        if (!allowed.contains(key)) {
            fail(std::string(path) + "." + key, "unknown key");
        }
    }
}

const JsonValue& require_key(const JsonValue& object,
                             std::string_view key,
                             std::string_view path) {
    require_object(object, path);
    if (!object.contains(std::string(key))) {
        fail(path, "missing required key '" + std::string(key) + "'");
    }
    return object.at(std::string(key));
}

std::string require_string(const JsonValue& value, std::string_view path) {
    if (!value.is_string()) { fail(path, "expected a string"); }
    const std::string result = value.as_string();
    if (result.empty()) { fail(path, "must not be empty"); }
    return result;
}

bool require_bool(const JsonValue& value, std::string_view path) {
    if (!value.is_boolean()) { fail(path, "expected a boolean"); }
    return value.as_boolean();
}

double require_number(const JsonValue& value, std::string_view path) {
    if (!value.is_number()) { fail(path, "expected a number"); }
    const double result = value.as_number<double>();
    if (!std::isfinite(result)) { fail(path, "must be finite"); }
    return result;
}

long long require_integer(const JsonValue& value, std::string_view path) {
    if (value.value_category() != meow::json::ValueCategory::NumberInt) {
        fail(path, "expected an integer");
    }
    return value.as_number<long long>();
}

int require_int(const JsonValue& value, std::string_view path) {
    const long long result = require_integer(value, path);
    if (result < std::numeric_limits<int>::min() ||
        result > std::numeric_limits<int>::max()) {
        fail(path, "integer is outside the supported range");
    }
    return static_cast<int>(result);
}

std::size_t require_size(const JsonValue& value, std::string_view path) {
    const long long result = require_integer(value, path);
    if (result < 0) { fail(path, "must be nonnegative"); }
    return static_cast<std::size_t>(result);
}

std::optional<double> optional_number(const JsonValue& object,
                                      std::string_view key,
                                      std::string_view path) {
    if (!object.contains(std::string(key))) { return std::nullopt; }
    const JsonValue& value = object.at(std::string(key));
    if (value.is_null()) { return std::nullopt; }
    return require_number(value, std::string(path) + "." + std::string(key));
}

std::vector<double> number_array(const JsonValue& value,
                                 std::string_view path) {
    const auto& array = require_array(value, path);
    std::vector<double> result;
    result.reserve(array.size());
    for (std::size_t index = 0; index < array.size(); ++index) {
        result.push_back(require_number(
            array[index],
            std::string(path) + "[" + std::to_string(index) + "]"));
    }
    return result;
}

QuasisymmetryCase parse_case(const JsonValue& value, std::string_view path) {
    const std::string name = require_string(value, path);
    if (name == "qa") { return QuasisymmetryCase::QA; }
    if (name == "qh") { return QuasisymmetryCase::QH; }
    fail(path, "must be 'qa' or 'qh'");
}

WorkflowKind parse_workflow(const JsonValue& value, std::string_view path) {
    const std::string name = require_string(value, path);
    if (name == "construction") { return WorkflowKind::CONSTRUCTION; }
    if (name == "refinement") { return WorkflowKind::REFINEMENT; }
    fail(path, "must be 'construction' or 'refinement'");
}

RadialTransfer parse_radial_transfer(const JsonValue& value,
                                     std::string_view path) {
    const std::string name = require_string(value, path);
    if (name == "default") { return RadialTransfer::DEFAULT; }
    if (name == "catmull-rom") { return RadialTransfer::CATMULL_ROM; }
    fail(path, "must be 'default' or 'catmull-rom'");
}

TargetSpec parse_target(const JsonValue& value) {
    constexpr std::string_view path = "target";
    reject_unknown_keys(value, path,
                        {"aspect_ratio", "mean_iota", "helicity_m",
                         "helicity_n_per_field_period", "ntheta", "nzeta",
                         "surfaces", "surface_weights"});
    TargetSpec result;
    result.aspect_ratio = require_number(
        require_key(value, "aspect_ratio", path), "target.aspect_ratio");
    result.mean_iota = optional_number(value, "mean_iota", path);
    result.helicity_m = require_int(require_key(value, "helicity_m", path),
                                    "target.helicity_m");
    result.helicity_n_per_field_period =
        require_int(require_key(value, "helicity_n_per_field_period", path),
                    "target.helicity_n_per_field_period");
    result.ntheta =
        require_int(require_key(value, "ntheta", path), "target.ntheta");
    result.nzeta =
        require_int(require_key(value, "nzeta", path), "target.nzeta");
    result.normalized_toroidal_flux_surfaces =
        number_array(require_key(value, "surfaces", path), "target.surfaces");
    result.surface_weights = number_array(
        require_key(value, "surface_weights", path), "target.surface_weights");

    if (!(result.aspect_ratio > 0.0)) {
        fail("target.aspect_ratio", "must be positive");
    }
    if (result.helicity_m == 0) {
        fail("target.helicity_m", "must be nonzero");
    }
    if (result.ntheta < 1 || result.nzeta < 1) {
        fail("target", "ntheta and nzeta must be positive");
    }
    if (result.normalized_toroidal_flux_surfaces.empty()) {
        fail("target.surfaces", "must not be empty");
    }
    if (result.normalized_toroidal_flux_surfaces.size() !=
        result.surface_weights.size()) {
        fail("target.surface_weights", "must have one entry per surface");
    }
    for (std::size_t index = 0;
         index < result.normalized_toroidal_flux_surfaces.size(); ++index) {
        const double surface = result.normalized_toroidal_flux_surfaces[index];
        if (surface < 0.0 || surface > 1.0) {
            fail("target.surfaces[" + std::to_string(index) + "]",
                 "must be in [0, 1]");
        }
        if (index > 0 &&
            surface <= result.normalized_toroidal_flux_surfaces[index - 1]) {
            fail("target.surfaces", "must be strictly increasing");
        }
        if (result.surface_weights[index] < 0.0) {
            fail("target.surface_weights[" + std::to_string(index) + "]",
                 "must be nonnegative");
        }
    }
    return result;
}

InitializationSpec parse_initialization(const JsonValue& value) {
    constexpr std::string_view path = "initialization";
    reject_unknown_keys(value, path,
                        {"axisymmetric_seed_amplitude", "track_axis_predictor",
                         "radial_transfer"});
    InitializationSpec result;
    result.axisymmetric_seed_amplitude =
        optional_number(value, "axisymmetric_seed_amplitude", path);
    result.track_axis_predictor =
        require_bool(require_key(value, "track_axis_predictor", path),
                     "initialization.track_axis_predictor");
    result.radial_transfer =
        parse_radial_transfer(require_key(value, "radial_transfer", path),
                              "initialization.radial_transfer");
    if (result.axisymmetric_seed_amplitude.has_value() &&
        !(*result.axisymmetric_seed_amplitude > 0.0)) {
        fail("initialization.axisymmetric_seed_amplitude",
             "must be positive when present");
    }
    return result;
}

JacobianSpec parse_jacobian(const JsonValue& value) {
    constexpr std::string_view path = "optimizer.jacobian";
    reject_unknown_keys(
        value, path,
        {"method", "workers", "relative_step", "absolute_step",
         "refresh_interval", "minimum_reduction_ratio", "maximum_secant_error",
         "scale_variables", "warm_start_trials"});
    JacobianSpec result;
    result.method = require_string(require_key(value, "method", path),
                                   "optimizer.jacobian.method");
    result.workers = require_size(require_key(value, "workers", path),
                                  "optimizer.jacobian.workers");
    result.relative_step =
        require_number(require_key(value, "relative_step", path),
                       "optimizer.jacobian.relative_step");
    result.absolute_step =
        require_number(require_key(value, "absolute_step", path),
                       "optimizer.jacobian.absolute_step");
    result.refresh_interval =
        require_size(require_key(value, "refresh_interval", path),
                     "optimizer.jacobian.refresh_interval");
    result.minimum_reduction_ratio =
        require_number(require_key(value, "minimum_reduction_ratio", path),
                       "optimizer.jacobian.minimum_reduction_ratio");
    result.maximum_secant_error =
        require_number(require_key(value, "maximum_secant_error", path),
                       "optimizer.jacobian.maximum_secant_error");
    result.scale_variables =
        require_bool(require_key(value, "scale_variables", path),
                     "optimizer.jacobian.scale_variables");
    result.warm_start_trials =
        require_bool(require_key(value, "warm_start_trials", path),
                     "optimizer.jacobian.warm_start_trials");

    const std::set<std::string> methods = {
        "analytic",
        "finite-difference",
        "broyden",
        "aggressive-broyden",
        "parallel-aggressive-broyden",
        "four-worker-aggressive-broyden",
        "extended-four-worker-broyden",
        "jacobian-scaled",
        "two-accuracy",
        "two-accuracy-broyden",
        "geometry-restart-check",
        "geometry-restart-finite-difference",
        "parallel-finite-difference-check",
        "relaxed-parallel-finite-difference-check",
        "parallel-worker-count-check",
        "parallel-finite-difference",
        "hot-finite-difference",
        "warm-finite-difference",
    };
    if (!methods.contains(result.method)) {
        fail("optimizer.jacobian.method",
             "unsupported method '" + result.method + "'");
    }
    if (result.workers == 0) {
        fail("optimizer.jacobian.workers", "must be positive");
    }
    if (result.relative_step < 0.0 || result.absolute_step < 0.0) {
        fail("optimizer.jacobian", "difference steps must be nonnegative");
    }
    if (result.refresh_interval == 0) {
        fail("optimizer.jacobian.refresh_interval", "must be positive");
    }
    if (result.minimum_reduction_ratio < 0.0 ||
        result.maximum_secant_error < 0.0) {
        fail("optimizer.jacobian", "Broyden safeguards must be nonnegative");
    }
    return result;
}

OptimizerSpec parse_optimizer(const JsonValue& value) {
    constexpr std::string_view path = "optimizer";
    reject_unknown_keys(value, path,
                        {"ftol", "xtol", "gtol", "max_function_evaluations",
                         "max_accepted_iterations", "jacobian"});
    OptimizerSpec result;
    result.ftol =
        require_number(require_key(value, "ftol", path), "optimizer.ftol");
    result.xtol =
        require_number(require_key(value, "xtol", path), "optimizer.xtol");
    result.gtol =
        require_number(require_key(value, "gtol", path), "optimizer.gtol");
    result.max_function_evaluations =
        require_size(require_key(value, "max_function_evaluations", path),
                     "optimizer.max_function_evaluations");
    result.max_accepted_iterations =
        require_size(require_key(value, "max_accepted_iterations", path),
                     "optimizer.max_accepted_iterations");
    result.jacobian = parse_jacobian(require_key(value, "jacobian", path));
    if (result.ftol < 0.0 || result.xtol < 0.0 || result.gtol < 0.0 ||
        (result.ftol == 0.0 && result.xtol == 0.0 && result.gtol == 0.0)) {
        fail("optimizer", "tolerances must be nonnegative and not all zero");
    }
    return result;
}

PhaseStoppingSpec parse_phase_stopping(const JsonValue& value,
                                       const std::string& path) {
    reject_unknown_keys(
        value, path,
        {"minimum_iterations", "progress_window", "minimum_relative_progress"});
    PhaseStoppingSpec result;
    result.minimum_iterations =
        require_size(require_key(value, "minimum_iterations", path),
                     path + ".minimum_iterations");
    result.progress_window = require_size(
        require_key(value, "progress_window", path), path + ".progress_window");
    result.minimum_relative_progress =
        require_number(require_key(value, "minimum_relative_progress", path),
                       path + ".minimum_relative_progress");
    if (result.minimum_relative_progress < 0.0) {
        fail(path + ".minimum_relative_progress", "must be nonnegative");
    }
    if ((result.progress_window == 0) !=
        (result.minimum_relative_progress == 0.0)) {
        fail(path,
             "progress_window and minimum_relative_progress must either both "
             "be zero or both be positive");
    }
    return result;
}

RelaxationPhaseSpec parse_phase(const JsonValue& value,
                                const std::string& path) {
    reject_unknown_keys(
        value, path,
        {"name", "equilibrium_tolerance_floor", "cold_start", "stopping"});
    RelaxationPhaseSpec result;
    result.name =
        require_string(require_key(value, "name", path), path + ".name");
    result.equilibrium_tolerance_floor =
        optional_number(value, "equilibrium_tolerance_floor", path);
    result.cold_start = require_bool(require_key(value, "cold_start", path),
                                     path + ".cold_start");
    result.stopping = parse_phase_stopping(require_key(value, "stopping", path),
                                           path + ".stopping");
    if (result.equilibrium_tolerance_floor.has_value() &&
        !(*result.equilibrium_tolerance_floor > 0.0)) {
        fail(path + ".equilibrium_tolerance_floor",
             "must be positive when present");
    }
    return result;
}

RelaxationStepSpec parse_step(const JsonValue& value, std::size_t index) {
    const std::string path = "rundown[" + std::to_string(index) + "]";
    reject_unknown_keys(value, path,
                        {"name", "max_boundary_mode", "equilibrium", "phases"});
    RelaxationStepSpec result;
    result.name =
        require_string(require_key(value, "name", path), path + ".name");
    result.max_boundary_mode =
        require_int(require_key(value, "max_boundary_mode", path),
                    path + ".max_boundary_mode");

    const JsonValue& equilibrium = require_key(value, "equilibrium", path);
    const std::string equilibrium_path = path + ".equilibrium";
    reject_unknown_keys(
        equilibrium, equilibrium_path,
        {"mpol", "ntor", "minimum_iterations", "tolerance_floor"});
    result.equilibrium.mpol =
        require_int(require_key(equilibrium, "mpol", equilibrium_path),
                    equilibrium_path + ".mpol");
    result.equilibrium.ntor =
        require_int(require_key(equilibrium, "ntor", equilibrium_path),
                    equilibrium_path + ".ntor");
    result.equilibrium.minimum_iterations = require_size(
        require_key(equilibrium, "minimum_iterations", equilibrium_path),
        equilibrium_path + ".minimum_iterations");
    result.equilibrium.tolerance_floor = require_number(
        require_key(equilibrium, "tolerance_floor", equilibrium_path),
        equilibrium_path + ".tolerance_floor");

    const auto& phases =
        require_array(require_key(value, "phases", path), path + ".phases");
    for (std::size_t phase = 0; phase < phases.size(); ++phase) {
        result.phases.push_back(parse_phase(
            phases[phase], path + ".phases[" + std::to_string(phase) + "]"));
    }
    if (result.max_boundary_mode < 1) {
        fail(path + ".max_boundary_mode", "must be positive");
    }
    const bool retains_resolution =
        result.equilibrium.mpol == 0 && result.equilibrium.ntor == 0;
    const bool sets_resolution =
        result.equilibrium.mpol > 0 && result.equilibrium.ntor > 0;
    if (!retains_resolution && !sets_resolution) {
        fail(equilibrium_path,
             "mpol and ntor must either both be zero or both be positive");
    }
    if (result.equilibrium.tolerance_floor < 0.0) {
        fail(equilibrium_path + ".tolerance_floor", "must be nonnegative");
    }
    if (result.phases.empty()) { fail(path + ".phases", "must not be empty"); }
    std::set<std::string> phase_names;
    for (const auto& phase : result.phases) {
        if (!phase_names.insert(phase.name).second) {
            fail(path + ".phases", "phase names must be unique within a step");
        }
    }
    return result;
}

OutputSpec parse_output(const JsonValue& value) {
    constexpr std::string_view path = "output";
    reject_unknown_keys(value, path, {"path", "iteration_directory"});
    OutputSpec result;
    result.path =
        require_string(require_key(value, "path", path), "output.path");
    const JsonValue& directory =
        require_key(value, "iteration_directory", path);
    if (!directory.is_string()) {
        fail("output.iteration_directory", "expected a string");
    }
    result.iteration_directory = directory.as_string();
    return result;
}

RelaxationRundown parse_root(JsonValue root, std::string_view source_path) {
    reject_unknown_keys(
        root, "rundown document",
        {"schema", "name", "equilibrium_input", "case", "workflow", "target",
         "initialization", "optimizer", "output", "rundown"});
    RelaxationRundown result;
    result.schema =
        require_string(require_key(root, "schema", "document"), "schema");
    if (result.schema != RELAXATION_RUNDOWN_SCHEMA) {
        fail("schema", "unsupported schema '" + result.schema + "'");
    }
    result.name = require_string(require_key(root, "name", "document"), "name");
    result.equilibrium_input =
        require_string(require_key(root, "equilibrium_input", "document"),
                       "equilibrium_input");
    result.selected_case =
        parse_case(require_key(root, "case", "document"), "case");
    result.workflow =
        parse_workflow(require_key(root, "workflow", "document"), "workflow");
    result.target = parse_target(require_key(root, "target", "document"));
    result.initialization =
        parse_initialization(require_key(root, "initialization", "document"));
    result.optimizer =
        parse_optimizer(require_key(root, "optimizer", "document"));
    result.output = parse_output(require_key(root, "output", "document"));

    if (result.selected_case == QuasisymmetryCase::QA &&
        result.target.helicity_n_per_field_period != 0) {
        fail("target.helicity_n_per_field_period",
             "must be zero for a QA rundown");
    }
    if (result.selected_case == QuasisymmetryCase::QH &&
        result.target.helicity_n_per_field_period == 0) {
        fail("target.helicity_n_per_field_period",
             "must be nonzero for a QH rundown");
    }
    if (result.selected_case == QuasisymmetryCase::QH &&
        result.target.mean_iota.has_value()) {
        fail("target.mean_iota", "is only supported for a QA rundown");
    }

    const auto& steps =
        require_array(require_key(root, "rundown", "document"), "rundown");
    for (std::size_t index = 0; index < steps.size(); ++index) {
        result.steps.push_back(parse_step(steps[index], index));
    }
    if (result.steps.empty()) { fail("rundown", "must not be empty"); }
    std::set<std::string> step_names;
    int previous_mode = 0;
    for (const auto& step : result.steps) {
        if (!step_names.insert(step.name).second) {
            fail("rundown", "step names must be unique");
        }
        if (step.max_boundary_mode <= previous_mode) {
            fail("rundown", "max_boundary_mode must be strictly increasing");
        }
        previous_mode = step.max_boundary_mode;
    }

    if (!source_path.empty()) {
        std::filesystem::path input(result.equilibrium_input);
        if (input.is_relative()) {
            input = std::filesystem::path(source_path).parent_path() / input;
            result.equilibrium_input = input.lexically_normal().string();
        }
    }
    return result;
}

}  // namespace

RelaxationRundown read_relaxation_rundown(const std::string& path) {
    try {
        return parse_root(meow::json::parse_file(path), path);
    } catch (const std::exception& error) {
        throw std::runtime_error(path + ": " + error.what());
    }
}

RelaxationRundown parse_relaxation_rundown(std::string_view document,
                                           std::string_view source_path) {
    try {
        return parse_root(meow::json::parse(std::string(document)),
                          source_path);
    } catch (const std::exception& error) {
        const std::string prefix =
            source_path.empty() ? "rundown" : std::string(source_path);
        throw std::runtime_error(prefix + ": " + error.what());
    }
}

std::string quasisymmetry_case_name(QuasisymmetryCase selected_case) {
    return selected_case == QuasisymmetryCase::QA ? "qa" : "qh";
}

std::string workflow_kind_name(WorkflowKind workflow) {
    return workflow == WorkflowKind::CONSTRUCTION ? "construction"
                                                  : "refinement";
}

std::string radial_transfer_name(RadialTransfer transfer) {
    return transfer == RadialTransfer::CATMULL_ROM ? "catmull-rom" : "default";
}

std::string summarize_relaxation_rundown(const RelaxationRundown& rundown) {
    std::ostringstream output;
    output << "run=" << rundown.name
           << " case=" << quasisymmetry_case_name(rundown.selected_case)
           << " workflow=" << workflow_kind_name(rundown.workflow)
           << " input=" << rundown.equilibrium_input
           << " output=" << rundown.output.path
           << " jacobian=" << rundown.optimizer.jacobian.method
           << " workers=" << rundown.optimizer.jacobian.workers << '\n';
    for (std::size_t index = 0; index < rundown.steps.size(); ++index) {
        const RelaxationStepSpec& step = rundown.steps[index];
        output << "step=" << index << " name=" << step.name
               << " max_mode=" << step.max_boundary_mode
               << " mpol=" << step.equilibrium.mpol
               << " ntor=" << step.equilibrium.ntor
               << " minimum_iterations=" << step.equilibrium.minimum_iterations
               << " tolerance_floor=" << step.equilibrium.tolerance_floor
               << " phases=" << step.phases.size() << '\n';
    }
    return output.str();
}

}  // namespace meow::config
