#include "test_support.hpp"

#include <stdexcept>

#include <meow/cumes/landreman_workflow.hpp>

int main() {
    using namespace cumes_meow_example;
    using meow::test::check;

    const LandremanSelection qa_construction =
        parse_landreman_selection("qa-construction");
    const auto qa_stages = landreman_stages(qa_construction);
    check(qa_stages.size() == 4, "QA construction has four mode stages");
    check(qa_stages.front().max_mode == 1 && qa_stages.front().mpol == 3 &&
              qa_stages.front().ntor == 3 &&
              qa_stages.front().minimum_niter == 6000 &&
              qa_stages.front().tolerance_floor == 1.0e-12,
          "QA construction starts at mode 1 and transform resolution 3");
    check(qa_stages.back().max_mode == 4 && qa_stages.back().mpol == 6 &&
              qa_stages.back().minimum_niter == 30000 &&
              qa_stages.back().tolerance_floor == 1.0e-12,
          "QA construction ends at mode 4 with the qualified solve cap");
    check(landreman_targets_mean_iota(qa_construction),
          "QA construction targets mean iota");
    check(landreman_final_surface_weight(qa_construction) == 1.0,
          "construction uses uniform surface weights");
    const auto qa_difference =
        landreman_finite_difference_policy(qa_construction);
    check(qa_difference.relative_step == 3.1622776601683794e-3 &&
              qa_difference.absolute_step == 1.0e-7,
          "QA construction uses run-021 finite differences");

    const LandremanSelection qh_construction =
        parse_landreman_selection("qh-construction");
    const auto qh_stages = landreman_stages(qh_construction);
    check(qh_stages.size() == 5 && qh_stages.back().max_mode == 5,
          "QH construction continues through mode 5");
    check(!landreman_targets_mean_iota(qh_construction),
          "QH construction omits mean iota target");
    const auto qh_difference =
        landreman_finite_difference_policy(qh_construction);
    check(qh_difference.relative_step == 1.0e-3 &&
              qh_difference.absolute_step == 1.0e-7,
          "QH construction uses run-039 finite differences");

    const LandremanSelection qa_refinement = parse_landreman_selection("qa");
    const auto refinement_stages = landreman_stages(qa_refinement);
    check(refinement_stages.size() == 2 &&
              refinement_stages.front().max_mode == 4 &&
              refinement_stages.back().max_mode == 5,
          "existing QA case retains mode-4/mode-5 refinement");
    check(refinement_stages.front().mpol == 0,
          "refinement retains the input equilibrium resolution");
    check(refinement_stages.front().minimum_niter == 0,
          "refinement retains the input iteration caps");
    check(refinement_stages.front().tolerance_floor == 0.0,
          "refinement retains the input tolerances");
    check(landreman_final_surface_weight(qa_refinement) == 30.0,
          "QA refinement retains edge weight 30");
    const auto refinement_difference =
        landreman_finite_difference_policy(qa_refinement);
    check(refinement_difference.relative_step == 1.0e-5 &&
              refinement_difference.absolute_step == 1.0e-9,
          "refinement retains its smaller finite differences");

    bool rejected = false;
    try {
        static_cast<void>(parse_landreman_selection("qa-full-ish"));
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "unknown workflow names are rejected");

    return meow::test::summary();
}
