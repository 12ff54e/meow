#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace meow::config {

inline constexpr std::string_view RELAXATION_RUNDOWN_SCHEMA =
    "meow-relaxation-rundown-v1";

enum class QuasisymmetryCase { QA, QH };
enum class WorkflowKind { CONSTRUCTION, REFINEMENT };
enum class RadialTransfer { DEFAULT, CATMULL_ROM };

struct TargetSpec {
    double aspect_ratio = 0.0;
    std::optional<double> mean_iota;
    int helicity_m = 1;
    int helicity_n_per_field_period = 0;
    int ntheta = 63;
    int nzeta = 64;
    std::vector<double> normalized_toroidal_flux_surfaces;
    std::vector<double> surface_weights;
};

struct InitializationSpec {
    std::optional<double> axisymmetric_seed_amplitude;
    bool track_axis_predictor = false;
    RadialTransfer radial_transfer = RadialTransfer::DEFAULT;
};

struct JacobianSpec {
    std::string method = "finite-difference";
    std::size_t workers = 1;
    double relative_step = 0.0;
    double absolute_step = 0.0;
    std::size_t refresh_interval = 1;
    double minimum_reduction_ratio = 0.1;
    double maximum_secant_error = 0.5;
    bool scale_variables = false;
    bool warm_start_trials = false;
};

struct OptimizerSpec {
    double ftol = 1.0e-8;
    double xtol = 1.0e-8;
    double gtol = 1.0e-8;
    std::size_t max_function_evaluations = 0;
    std::size_t max_accepted_iterations = 0;
    JacobianSpec jacobian;
};

struct EquilibriumStepSpec {
    int mpol = 0;
    int ntor = 0;
    std::size_t minimum_iterations = 0;
    double tolerance_floor = 0.0;
};

struct PhaseStoppingSpec {
    std::size_t minimum_iterations = 0;
    std::size_t progress_window = 0;
    double minimum_relative_progress = 0.0;
};

struct RelaxationPhaseSpec {
    std::string name;
    std::optional<double> equilibrium_tolerance_floor;
    bool cold_start = false;
    PhaseStoppingSpec stopping;
};

struct RelaxationStepSpec {
    std::string name;
    int max_boundary_mode = 0;
    EquilibriumStepSpec equilibrium;
    std::vector<RelaxationPhaseSpec> phases;
};

struct OutputSpec {
    std::string path;
    std::string iteration_directory;
};

struct RelaxationRundown {
    std::string schema;
    std::string name;
    std::string equilibrium_input;
    QuasisymmetryCase selected_case = QuasisymmetryCase::QA;
    WorkflowKind workflow = WorkflowKind::CONSTRUCTION;
    TargetSpec target;
    InitializationSpec initialization;
    OptimizerSpec optimizer;
    OutputSpec output;
    std::vector<RelaxationStepSpec> steps;
};

// Relative equilibrium_input paths are resolved against the rundown file.
RelaxationRundown read_relaxation_rundown(const std::string& path);

// `source_path` is used only to resolve a relative equilibrium_input and to
// prefix diagnostics. It may be empty for an in-memory document.
RelaxationRundown parse_relaxation_rundown(std::string_view document,
                                           std::string_view source_path = {});

std::string quasisymmetry_case_name(QuasisymmetryCase selected_case);
std::string workflow_kind_name(WorkflowKind workflow);
std::string radial_transfer_name(RadialTransfer transfer);
std::string summarize_relaxation_rundown(const RelaxationRundown& rundown);

}  // namespace meow::config
