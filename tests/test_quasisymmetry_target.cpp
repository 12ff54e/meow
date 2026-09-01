#include "test_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>
#include <vector>

#include <meow/cumes/quasisymmetry_target.hpp>
#include <meow/cumes/quasisymmetry_target_jvp.hpp>

int main() {
    using meow::test::check;

    constexpr int ns = 4;
    constexpr int ntheta = 16;
    constexpr int nzeta = 1;
    constexpr double major_radius = 5.0;
    constexpr double minor_radius = 1.0;
    constexpr double b0 = 2.0;
    constexpr double b_amplitude = 0.25;
    constexpr double sqrtg_value = -2.0;
    constexpr double g_covariant = 2.0;
    constexpr double iota = 0.5;

    cumes::EquilibriumSnapshot equilibrium;
    equilibrium.ns = ns;
    equilibrium.mnmax = 2;
    equilibrium.ntheta = ntheta;
    equilibrium.nzeta = nzeta;
    for (auto& family : equilibrium.families) {
        family.assign(equilibrium.family_size(), 0.0);
    }
    for (int surface = 0; surface < ns; ++surface) {
        equilibrium.families[cumes::EquilibriumSnapshot::RMNCC]
                            [static_cast<std::size_t>(surface)] = major_radius;
        equilibrium.families[cumes::EquilibriumSnapshot::RMNCC]
                            [static_cast<std::size_t>(ns + surface)] =
            minor_radius;
        equilibrium.families[cumes::EquilibriumSnapshot::ZMNSC]
                            [static_cast<std::size_t>(ns + surface)] =
            minor_radius;
    }
    for (auto& field : equilibrium.half_fields) {
        field.assign(equilibrium.half_field_size(), 0.0);
    }
    for (auto& field : equilibrium.full_fields) {
        field.assign(equilibrium.full_field_size(), 0.0);
    }

    double i_covariant = 0.0;
    for (int surface = 0; surface < ns - 1; ++surface) {
        for (int theta_index = 0; theta_index < ntheta; ++theta_index) {
            const double theta = 2.0 * std::numbers::pi * theta_index / ntheta;
            const double b = b0 + b_amplitude * std::sin(theta);
            const std::size_t index =
                static_cast<std::size_t>(surface) * ntheta + theta_index;
            equilibrium.half_fields[cumes::EquilibriumSnapshot::SQRTG][index] =
                sqrtg_value;
            equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUPU][index] =
                1.0;
            equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBU][index] =
                b * b;
            equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBV][index] =
                g_covariant;
            if (surface == 0) i_covariant += b * b / ntheta;
        }
    }

    cumes::EquilibriumProfiles profiles;
    profiles.toroidal_flux_derivative.assign(ns - 1, 2.0 * std::numbers::pi);
    // With the negative cuMES Jacobian orientation this gives
    // psi_p'=-4 after division by 2*pi.
    profiles.poloidal_flux_derivative.assign(ns - 1, 8.0 * std::numbers::pi);
    profiles.rotational_transform.assign(ns - 1, iota);
    profiles.poloidal_covariant_field.assign(ns - 1, i_covariant);
    profiles.toroidal_covariant_field.assign(ns - 1, g_covariant);

    cumes::InputParams input;
    input.mpol = 2;
    input.ntor = 0;
    input.nfp = 1;
    input.ntheta = ntheta;
    input.nzeta = nzeta;

    cumes_meow_example::QuasisymmetryTargetSpec qa_spec;
    qa_spec.helicity_m = 1;
    qa_spec.helicity_n = 0;
    qa_spec.half_surface_indices =
        cumes_meow_example::all_half_grid_surfaces(ns);
    qa_spec.flux_gradient =
        cumes_meow_example::QsFluxGradient::NORMALIZED_POLOIDAL;
    const auto qs = cumes_meow_example::calculate_quasisymmetry_target(
        equilibrium, profiles, input.nfp, qa_spec);
    check(qs.residuals.size() ==
              static_cast<std::size_t>((ns - 1) * ntheta * nzeta),
          "QS target: one residual per selected angular point");
    check(std::abs(qs.value) < 1.0e-24,
          "QS target: manufactured QA field has zero residual");

    const auto qa = cumes_meow_example::calculate_qa_target(
        equilibrium, profiles, input, qa_spec, major_radius / minor_radius,
        iota);
    check(qa.residuals.size() == qs.residuals.size() + 2,
          "QA target: aspect and iota residuals appended");
    check(std::abs(qa.value) < 1.0e-24,
          "QA target: matched composite target is zero");

    const auto penalized_qa = cumes_meow_example::calculate_qa_target(
        equilibrium, profiles, input, qa_spec,
        major_radius / minor_radius - 0.25, iota - 0.1);
    check(std::abs(penalized_qa.value - (0.25 * 0.25 + 0.1 * 0.1)) < 1.0e-13,
          "QA target: scalar is squared norm of composite residuals");

    cumes_meow_example::QuasisymmetryTargetSpec qh_spec = qa_spec;
    qh_spec.helicity_n = 1;
    const auto qh = cumes_meow_example::calculate_qh_target(
        equilibrium, profiles, input, qh_spec, major_radius / minor_radius);
    check(qh.residuals.size() == qh.qs.residuals.size() + 1,
          "QH target: aspect residual appended");
    check(std::abs(qh.value - qh.qs.value) < 1.0e-13,
          "QH target: matched aspect leaves QS value");

    cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec paper_grid;
    paper_grid.helicity_m = 1;
    paper_grid.helicity_n = 0;
    paper_grid.normalized_toroidal_flux_surfaces = {0.0, 0.5, 1.0};
    paper_grid.surface_weights = {1.0, 2.0, 3.0};
    paper_grid.target_ntheta = 21;
    paper_grid.target_nzeta = 20;
    paper_grid.flux_gradient =
        cumes_meow_example::QsFluxGradient::NORMALIZED_POLOIDAL;
    const auto resampled = cumes_meow_example::calculate_quasisymmetry_target(
        equilibrium, profiles, input.nfp, paper_grid);
    check(resampled.residuals.size() == 3 * 21 * 20,
          "flux-surface QS target uses the independent angular grid");
    check(resampled.surface_values.size() == 3,
          "flux-surface QS target includes axis, interior, and edge");
    check(std::abs(resampled.value) < 1.0e-23,
          "linear extrapolation and Fourier resampling preserve zero QS");

    const auto paper_qa = cumes_meow_example::calculate_qa_target(
        equilibrium, profiles, input, paper_grid, major_radius / minor_radius,
        std::nullopt);
    check(paper_qa.residuals.size() == resampled.residuals.size() + 1,
          "paper final QA target omits the initial-stage iota residual");
    const auto initial_qa = cumes_meow_example::calculate_qa_target(
        equilibrium, profiles, input, paper_grid, major_radius / minor_radius,
        iota);
    check(initial_qa.residuals.size() == resampled.residuals.size() + 2,
          "paper initial QA target can include the iota residual");

    cumes_meow_example::FluxSurfaceQuasisymmetryTargetSpec archived_grid;
    archived_grid.helicity_m = 1;
    archived_grid.helicity_n = 0;
    archived_grid.flux_gradient =
        cumes_meow_example::QsFluxGradient::NORMALIZED_POLOIDAL;
    for (int index = 0; index <= 10; ++index) {
        archived_grid.normalized_toroidal_flux_surfaces.push_back(index / 10.0);
    }
    const auto archived = cumes_meow_example::calculate_quasisymmetry_target(
        equilibrium, profiles, input.nfp, archived_grid);
    check(archived.residuals.size() == 11 * 63 * 64,
          "Landreman default target has 11*63*64 residuals");

    auto tangent = cumes::EquilibriumTangent::zero_like(equilibrium, profiles);
    for (int surface = 0; surface < ns; ++surface) {
        tangent.equilibrium.families[cumes::EquilibriumSnapshot::RMNCC]
                                    [static_cast<std::size_t>(surface)] = 0.02;
        tangent.equilibrium.families[cumes::EquilibriumSnapshot::RMNCC]
                                    [static_cast<std::size_t>(ns + surface)] =
            0.03;
        tangent.equilibrium.families[cumes::EquilibriumSnapshot::ZMNSC]
                                    [static_cast<std::size_t>(ns + surface)] =
            -0.01;
    }
    for (int surface = 0; surface < ns - 1; ++surface) {
        for (int theta_index = 0; theta_index < ntheta; ++theta_index) {
            const double theta = 2.0 * std::numbers::pi * theta_index / ntheta;
            const std::size_t index =
                static_cast<std::size_t>(surface) * ntheta + theta_index;
            tangent.equilibrium
                .half_fields[cumes::EquilibriumSnapshot::SQRTG][index] =
                0.01 * std::cos(theta);
            tangent.equilibrium
                .half_fields[cumes::EquilibriumSnapshot::BSUPU][index] =
                0.02 * std::sin(theta);
            tangent.equilibrium
                .half_fields[cumes::EquilibriumSnapshot::BSUPV][index] = -0.01;
            tangent.equilibrium
                .half_fields[cumes::EquilibriumSnapshot::BSUBU][index] =
                0.03 * std::cos(2.0 * theta);
            tangent.equilibrium
                .half_fields[cumes::EquilibriumSnapshot::BSUBV][index] =
                0.015 * std::sin(theta);
        }
    }
    tangent.profiles.toroidal_flux_derivative.assign(ns - 1, 0.04);
    tangent.profiles.poloidal_flux_derivative.assign(ns - 1, -0.03);
    tangent.profiles.rotational_transform.assign(ns - 1, 0.02);
    tangent.profiles.poloidal_covariant_field.assign(ns - 1, -0.015);
    tangent.profiles.toroidal_covariant_field.assign(ns - 1, 0.01);

    auto shifted = [&](double scale) {
        auto shifted_equilibrium = equilibrium;
        auto shifted_profiles = profiles;
        for (std::size_t component = 0;
             component < cumes::EquilibriumSnapshot::COUNT; ++component) {
            for (std::size_t index = 0;
                 index < shifted_equilibrium.families[component].size();
                 ++index) {
                shifted_equilibrium.families[component][index] +=
                    scale * tangent.equilibrium.families[component][index];
            }
        }
        for (std::size_t field = 0;
             field < cumes::EquilibriumSnapshot::HALF_FIELD_COUNT; ++field) {
            for (std::size_t index = 0;
                 index < shifted_equilibrium.half_fields[field].size();
                 ++index) {
                shifted_equilibrium.half_fields[field][index] +=
                    scale * tangent.equilibrium.half_fields[field][index];
            }
        }
        auto shift_profile = [scale](std::vector<double>& value,
                                     const std::vector<double>& direction) {
            for (std::size_t index = 0; index < value.size(); ++index) {
                value[index] += scale * direction[index];
            }
        };
        shift_profile(shifted_profiles.toroidal_flux_derivative,
                      tangent.profiles.toroidal_flux_derivative);
        shift_profile(shifted_profiles.poloidal_flux_derivative,
                      tangent.profiles.poloidal_flux_derivative);
        shift_profile(shifted_profiles.rotational_transform,
                      tangent.profiles.rotational_transform);
        shift_profile(shifted_profiles.poloidal_covariant_field,
                      tangent.profiles.poloidal_covariant_field);
        shift_profile(shifted_profiles.toroidal_covariant_field,
                      tangent.profiles.toroidal_covariant_field);
        return std::pair{std::move(shifted_equilibrium),
                         std::move(shifted_profiles)};
    };

    auto paper_qh_grid = paper_grid;
    paper_qh_grid.helicity_n = 1;
    constexpr double epsilon = 1.0e-6;
    const auto plus = shifted(epsilon);
    const auto minus = shifted(-epsilon);
    const auto plus_target = cumes_meow_example::calculate_qh_target(
        plus.first, plus.second, input, paper_qh_grid,
        major_radius / minor_radius);
    const auto minus_target = cumes_meow_example::calculate_qh_target(
        minus.first, minus.second, input, paper_qh_grid,
        major_radius / minor_radius);
    const auto analytic_jvp = cumes_meow_example::calculate_qh_target_jvp(
        equilibrium, profiles, tangent, input, paper_qh_grid);
    check(analytic_jvp.size() == plus_target.residuals.size(),
          "QH target JVP has one value per residual");
    double max_error = 0.0;
    double max_reference = 0.0;
    for (std::size_t index = 0; index < analytic_jvp.size(); ++index) {
        const double reference =
            (plus_target.residuals[index] - minus_target.residuals[index]) /
            (2.0 * epsilon);
        max_error =
            std::max(max_error, std::abs(analytic_jvp[index] - reference));
        max_reference = std::max(max_reference, std::abs(reference));
    }
    check(max_error < 2.0e-8 * std::max(1.0, max_reference),
          "QH target analytic JVP matches centered-difference oracle");

    return meow::test::summary();
}
