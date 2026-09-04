#include "test_support.hpp"

#include <functional>
#include <stdexcept>
#include <string>

#include <meow/config/relaxation_rundown.hpp>

namespace {

const char* VALID_RUNDOWN = R"json(
{
  "schema": "meow-relaxation-rundown-v1",
  "name": "qa-construction",
  "equilibrium_input": "equilibria/qa.json",
  "case": "qa",
  "workflow": "construction",
  "target": {
    "aspect_ratio": 6.0,
    "mean_iota": 0.42,
    "helicity_m": 1,
    "helicity_n_per_field_period": 0,
    "ntheta": 63,
    "nzeta": 64,
    "surfaces": [0.0, 0.5, 1.0],
    "surface_weights": [1.0, 1.0, 1.0]
  },
  "initialization": {
    "axisymmetric_seed_amplitude": 0.0001,
    "track_axis_predictor": true,
    "radial_transfer": "default"
  },
  "optimizer": {
    "ftol": 1e-8,
    "xtol": 1e-8,
    "gtol": 1e-8,
    "max_function_evaluations": 0,
    "max_accepted_iterations": 0,
    "jacobian": {
      "method": "parallel-aggressive-broyden",
      "workers": 4,
      "relative_step": 0.0031622776601683794,
      "absolute_step": 1e-7,
      "refresh_interval": 8,
      "minimum_reduction_ratio": 0.05,
      "maximum_secant_error": 0.5,
      "scale_variables": false,
      "warm_start_trials": false
    }
  },
  "output": {
    "path": "qa-result.json",
    "iteration_directory": "qa-steps"
  },
  "rundown": [
    {
      "name": "mode-1",
      "max_boundary_mode": 1,
      "equilibrium": {
        "mpol": 3,
        "ntor": 3,
        "minimum_iterations": 6000,
        "tolerance_floor": 1e-12
      },
      "phases": [
        {
          "name": "qualified",
          "equilibrium_tolerance_floor": null,
          "cold_start": false,
          "stopping": {
            "minimum_iterations": 0,
            "progress_window": 0,
            "minimum_relative_progress": 0.0
          }
        }
      ]
    },
    {
      "name": "mode-2",
      "max_boundary_mode": 2,
      "equilibrium": {
        "mpol": 5,
        "ntor": 5,
        "minimum_iterations": 6000,
        "tolerance_floor": 1e-12
      },
      "phases": [
        {
          "name": "relaxed",
          "equilibrium_tolerance_floor": 1e-9,
          "cold_start": false,
          "stopping": {
            "minimum_iterations": 8,
            "progress_window": 3,
            "minimum_relative_progress": 0.01
          }
        },
        {
          "name": "polish",
          "equilibrium_tolerance_floor": null,
          "cold_start": true,
          "stopping": {
            "minimum_iterations": 0,
            "progress_window": 0,
            "minimum_relative_progress": 0.0
          }
        }
      ]
    }
  ]
}
)json";

void expect_failure(const std::string& document, const std::string& needle) {
    bool failed = false;
    try {
        static_cast<void>(
            meow::config::parse_relaxation_rundown(document, "/tmp/run.json"));
    } catch (const std::runtime_error& error) {
        failed = std::string(error.what()).find(needle) != std::string::npos;
    }
    meow::test::check(failed, "invalid rundown reports " + needle);
}

}  // namespace

int main() {
    using meow::test::check;
    const meow::config::RelaxationRundown rundown =
        meow::config::parse_relaxation_rundown(VALID_RUNDOWN,
                                               "/work/config/run.json");
    check(rundown.schema == meow::config::RELAXATION_RUNDOWN_SCHEMA,
          "schema is retained");
    check(rundown.equilibrium_input == "/work/config/equilibria/qa.json",
          "relative equilibrium input resolves against rundown");
    check(rundown.selected_case == meow::config::QuasisymmetryCase::QA,
          "QA case is parsed");
    check(rundown.target.mean_iota == 0.42,
          "optional mean-iota target is parsed");
    check(rundown.optimizer.jacobian.workers == 4 &&
              rundown.optimizer.jacobian.refresh_interval == 8,
          "parallel Broyden policy is parsed");
    check(rundown.steps.size() == 2 && rundown.steps.back().phases.size() == 2,
          "ordered steps and phases are parsed");
    check(rundown.steps.back().phases.front().stopping.progress_window == 3,
          "phase stopping policy is parsed");
    check(meow::config::summarize_relaxation_rundown(rundown).find(
              "step=1 name=mode-2") != std::string::npos,
          "rundown summary lists steps");

    std::string unknown = VALID_RUNDOWN;
    unknown.replace(unknown.find("\"case\": \"qa\""), 12,
                    "\"case\": \"qa\", \"mystery\": 1");
    expect_failure(unknown, "mystery: unknown key");

    std::string bad_surfaces = VALID_RUNDOWN;
    bad_surfaces.replace(bad_surfaces.find("[0.0, 0.5, 1.0]"), 15,
                         "[0.0, 0.5, 0.4]");
    expect_failure(bad_surfaces, "must be strictly increasing");

    std::string bad_modes = VALID_RUNDOWN;
    const std::size_t second_mode = bad_modes.rfind("\"max_boundary_mode\": 2");
    bad_modes.replace(second_mode, 22, "\"max_boundary_mode\": 1");
    expect_failure(bad_modes, "max_boundary_mode must be strictly increasing");

    std::string bad_stopping = VALID_RUNDOWN;
    const std::size_t window = bad_stopping.find("\"progress_window\": 3");
    bad_stopping.replace(window, 20, "\"progress_window\": 0");
    expect_failure(bad_stopping,
                   "must either both be zero or both be positive");

    return meow::test::summary();
}
