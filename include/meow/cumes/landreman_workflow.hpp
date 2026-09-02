// Archived Landreman-Paul optimizer continuation policies.
#ifndef MEOW_CUMES_LANDREMAN_WORKFLOW_HPP_
#define MEOW_CUMES_LANDREMAN_WORKFLOW_HPP_

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include <meow/cumes/boundary_parameterization.hpp>

namespace cumes_meow_example {

enum class LandremanCase { QA, QH };
enum class LandremanWorkflow { CONSTRUCTION, REFINEMENT };

struct LandremanSelection {
    LandremanCase selected_case;
    LandremanWorkflow workflow;
};

struct LandremanStage {
    int max_mode;
    int mpol;
    int ntor;
    int minimum_niter;
    double tolerance_floor;
};

struct LandremanFiniteDifferencePolicy {
    double relative_step;
    double absolute_step;
};

inline LandremanSelection parse_landreman_selection(std::string_view name) {
    if (name == "qa") {
        return {LandremanCase::QA, LandremanWorkflow::REFINEMENT};
    }
    if (name == "qh") {
        return {LandremanCase::QH, LandremanWorkflow::REFINEMENT};
    }
    if (name == "qa-construction") {
        return {LandremanCase::QA, LandremanWorkflow::CONSTRUCTION};
    }
    if (name == "qh-construction") {
        return {LandremanCase::QH, LandremanWorkflow::CONSTRUCTION};
    }
    throw std::invalid_argument(
        "case must be qa, qh, qa-construction, or qh-construction");
}

inline std::vector<LandremanStage> landreman_stages(
    LandremanSelection selection) {
    if (selection.workflow == LandremanWorkflow::REFINEMENT) {
        // The refinement inputs already carry mpol=8 and ntor=6. Zero means
        // retain that input resolution for the stage.
        return {{4, 0, 0, 0, 0.0}, {5, 0, 0, 0, 0.0}};
    }
    if (selection.selected_case == LandremanCase::QA) {
        return {{1, 3, 3, 6000, 1.0e-12},
                {2, 5, 5, 6000, 1.0e-12},
                {3, 6, 6, 10000, 1.0e-12},
                {4, 6, 6, 30000, 1.0e-12}};
    }
    return {{1, 3, 3, 6000, 1.0e-12},
            {2, 5, 5, 6000, 1.0e-12},
            {3, 6, 6, 10000, 1.0e-12},
            {4, 6, 6, 10000, 1.0e-12},
            {5, 6, 6, 10000, 1.0e-12}};
}

inline LandremanFiniteDifferencePolicy landreman_finite_difference_policy(
    LandremanSelection selection) {
    if (selection.workflow == LandremanWorkflow::REFINEMENT) {
        return {1.0e-5, 1.0e-9};
    }
    if (selection.selected_case == LandremanCase::QA) {
        return {3.1622776601683794e-3, 1.0e-7};
    }
    return {1.0e-3, 1.0e-7};
}

inline double landreman_target_aspect(LandremanCase selected_case) {
    return selected_case == LandremanCase::QA ? 6.0 : 8.0;
}

inline double landreman_final_surface_weight(LandremanSelection selection) {
    if (selection.workflow == LandremanWorkflow::CONSTRUCTION) { return 1.0; }
    return selection.selected_case == LandremanCase::QA ? 30.0 : 2.0;
}

inline bool landreman_targets_mean_iota(LandremanSelection selection) {
    return selection.workflow == LandremanWorkflow::CONSTRUCTION &&
           selection.selected_case == LandremanCase::QA;
}

inline bool landreman_refreshes_axis_predictor(LandremanSelection selection) {
    return selection.workflow == LandremanWorkflow::CONSTRUCTION;
}

// The archived QA construction escaped its exactly axisymmetric stationary
// point through a one-sided finite-difference Jacobian. An analytic Jacobian
// instead needs an explicit 3-D seed because iota and the QA residual have no
// first derivative there. Preserve any user-supplied 3-D boundary unchanged.
inline bool seed_landreman_qa_construction_boundary(cumes::ProblemSpec& problem,
                                                    double amplitude = 1.0e-4) {
    if (!(amplitude > 0.0) || !std::isfinite(amplitude)) {
        throw std::invalid_argument(
            "QA construction seed amplitude must be finite and positive");
    }
    StellaratorSymmetricBoundaryParameterization boundary(1);
    meow::Vector values = boundary.values(problem);
    const auto& degrees = boundary.degrees_of_freedom();
    for (std::size_t index = 0; index < degrees.size(); ++index) {
        if (degrees[index].n != 0 &&
            values[static_cast<Eigen::Index>(index)] != 0.0) {
            return false;
        }
    }
    for (std::size_t index = 0; index < degrees.size(); ++index) {
        const auto& degree = degrees[index];
        if (degree.n == 0) { continue; }
        double sign = 1.0;
        if (degree.family == BoundaryFamily::RBC) {
            sign = degree.m == 0 || degree.n < 0 ? 1.0 : -1.0;
        } else {
            sign = degree.m == 0 ? -1.0 : 1.0;
        }
        values[static_cast<Eigen::Index>(index)] = sign * amplitude;
    }
    problem = boundary.apply(problem, values);
    return true;
}

inline std::string_view landreman_workflow_name(LandremanWorkflow workflow) {
    return workflow == LandremanWorkflow::CONSTRUCTION ? "construction"
                                                       : "refinement";
}

}  // namespace cumes_meow_example

#endif  // MEOW_CUMES_LANDREMAN_WORKFLOW_HPP_
