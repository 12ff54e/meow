// Archived Landreman-Paul optimizer continuation policies.
#ifndef MEOW_CUMES_LANDREMAN_WORKFLOW_HPP_
#define MEOW_CUMES_LANDREMAN_WORKFLOW_HPP_

#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

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
                {4, 6, 6, 30000, 2.0e-12}};
    }
    return {{1, 3, 3, 6000, 1.0e-12},
            {2, 5, 5, 6000, 1.0e-12},
            {3, 6, 6, 10000, 1.0e-12},
            {4, 6, 6, 10000, 1.0e-12},
            {5, 6, 6, 10000, 1.0e-12}};
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

inline std::string_view landreman_workflow_name(LandremanWorkflow workflow) {
    return workflow == LandremanWorkflow::CONSTRUCTION ? "construction"
                                                       : "refinement";
}

}  // namespace cumes_meow_example

#endif  // MEOW_CUMES_LANDREMAN_WORKFLOW_HPP_
