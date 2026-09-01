#pragma once

#include "plasma_size_target_jvp.hpp"
#include "quasisymmetry_target.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <vector>

#include <cumes/solver/equilibrium_tangent.hpp>

namespace cumes_meow_example {

namespace qs_jvp_detail {

inline void require_tangent_shape(const cumes::EquilibriumSnapshot& equilibrium,
                                  const cumes::EquilibriumProfiles& profiles,
                                  const cumes::EquilibriumTangent& tangent) {
    if (!tangent.matches(equilibrium, profiles)) {
        throw std::invalid_argument("QS equilibrium tangent shape mismatch");
    }
}

inline double interpolate_profile(const std::vector<double>& values,
                                  int ns,
                                  double surface) {
    return fourier_resampling::interpolate_half_grid_profile(values, ns,
                                                             surface);
}

}  // namespace qs_jvp_detail

// Apply the exact chain rule of the flux-surface QS residual discretization
// to one equilibrium direction. No boundary perturbation or equilibrium solve
// occurs here.
inline std::vector<double> calculate_quasisymmetry_target_jvp(
    const cumes::EquilibriumSnapshot& equilibrium,
    const cumes::EquilibriumProfiles& profiles,
    const cumes::EquilibriumTangent& tangent,
    int nfp,
    const FluxSurfaceQuasisymmetryTargetSpec& spec) {
    qs_jvp_detail::require_tangent_shape(equilibrium, profiles, tangent);
    static_cast<void>(
        calculate_quasisymmetry_target(equilibrium, profiles, nfp, spec));

    namespace resampling = fourier_resampling;
    const auto& dequilibrium = tangent.equilibrium;
    const auto& dprofiles = tangent.profiles;
    const auto& bsups =
        equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUPS];
    const auto& bsupu =
        equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUPU];
    const auto& bsupv =
        equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUPV];
    const auto& bsubs =
        equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBS];
    const auto& bsubu =
        equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBU];
    const auto& bsubv =
        equilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBV];
    const auto& sqrtg =
        equilibrium.half_fields[cumes::EquilibriumSnapshot::SQRTG];
    const auto& dbsups =
        dequilibrium.half_fields[cumes::EquilibriumSnapshot::BSUPS];
    const auto& dbsupu =
        dequilibrium.half_fields[cumes::EquilibriumSnapshot::BSUPU];
    const auto& dbsupv =
        dequilibrium.half_fields[cumes::EquilibriumSnapshot::BSUPV];
    const auto& dbsubs =
        dequilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBS];
    const auto& dbsubu =
        dequilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBU];
    const auto& dbsubv =
        dequilibrium.half_fields[cumes::EquilibriumSnapshot::BSUBV];
    const auto& dsqrtg =
        dequilibrium.half_fields[cumes::EquilibriumSnapshot::SQRTG];

    std::vector<double> native_b(equilibrium.half_field_size());
    std::vector<double> native_db(equilibrium.half_field_size());
    for (std::size_t index = 0; index < native_b.size(); ++index) {
        const double b_squared = bsups[index] * bsubs[index] +
                                 bsupu[index] * bsubu[index] +
                                 bsupv[index] * bsubv[index];
        const double db_squared =
            dbsups[index] * bsubs[index] + bsups[index] * dbsubs[index] +
            dbsupu[index] * bsubu[index] + bsupu[index] * dbsubu[index] +
            dbsupv[index] * bsubv[index] + bsupv[index] * dbsubv[index];
        native_b[index] = std::sqrt(b_squared);
        native_db[index] = 0.5 * db_squared / native_b[index];
    }

    const auto theta_weights = resampling::make_periodic_weights(
        equilibrium.ntheta, spec.target_ntheta);
    const auto zeta_weights =
        resampling::make_periodic_weights(equilibrium.nzeta, spec.target_nzeta);
    const std::size_t target_points =
        static_cast<std::size_t>(spec.target_ntheta * spec.target_nzeta);
    const std::size_t native_points = equilibrium.points_per_surface();
    const double two_pi = 2.0 * std::numbers::pi;
    std::vector<double> result;
    result.reserve(spec.normalized_toroidal_flux_surfaces.size() *
                   target_points);

    for (std::size_t requested = 0;
         requested < spec.normalized_toroidal_flux_surfaces.size();
         ++requested) {
        const double surface =
            spec.normalized_toroidal_flux_surfaces[requested];
        const double surface_weight = spec.surface_weights.empty()
                                          ? 1.0
                                          : spec.surface_weights[requested];
        if (surface_weight == 0.0) {
            result.insert(result.end(), target_points, 0.0);
            continue;
        }
        auto interpolate = [&](const std::vector<double>& field) {
            return resampling::interpolate_half_grid_field(
                field, equilibrium.ns, native_points, surface);
        };
        auto sample = [&](const std::vector<double>& field, bool derivatives) {
            return resampling::resample_2d(
                interpolate(field), equilibrium.ntheta, equilibrium.nzeta,
                theta_weights, zeta_weights, derivatives);
        };

        const auto b = sample(native_b, true);
        const auto db = sample(native_db, true);
        const auto g = sample(sqrtg, false).value;
        const auto dg = sample(dsqrtg, false).value;
        const auto bu = sample(bsupu, false).value;
        const auto dbu = sample(dbsupu, false).value;
        const auto bv = sample(bsupv, false).value;
        const auto dbv = sample(dbsupv, false).value;
        const auto bcu = sample(bsubu, false).value;
        const auto dbcu = sample(dbsubu, false).value;
        const auto bcv = sample(bsubv, false).value;
        const auto dbcv = sample(dbsubv, false).value;

        const double orientation =
            qs_target_detail::orientation_from_jacobian(g, 0, target_points);
        const auto profile = [&](const std::vector<double>& values) {
            return qs_jvp_detail::interpolate_profile(values, equilibrium.ns,
                                                      surface);
        };
        const double public_flux_prime =
            spec.flux_gradient == QsFluxGradient::NORMALIZED_TOROIDAL
                ? profile(profiles.toroidal_flux_derivative)
                : profile(profiles.poloidal_flux_derivative);
        const double dpublic_flux_prime =
            spec.flux_gradient == QsFluxGradient::NORMALIZED_TOROIDAL
                ? profile(dprofiles.toroidal_flux_derivative)
                : profile(dprofiles.poloidal_flux_derivative);
        const double flux_prime = orientation * public_flux_prime / two_pi;
        const double dflux_prime = orientation * dpublic_flux_prime / two_pi;
        const double iota = profile(profiles.rotational_transform);
        const double diota = profile(dprofiles.rotational_transform);
        const double i_cov = profile(profiles.poloidal_covariant_field);
        const double di_cov = profile(dprofiles.poloidal_covariant_field);
        const double g_cov = profile(profiles.toroidal_covariant_field);
        const double dg_cov = profile(dprofiles.toroidal_covariant_field);

        double jacobian_sum = 0.0;
        double djacobian_sum = 0.0;
        for (std::size_t point = 0; point < target_points; ++point) {
            jacobian_sum += std::abs(g[point]);
            djacobian_sum += std::copysign(1.0, g[point]) * dg[point];
        }
        for (std::size_t point = 0; point < target_points; ++point) {
            const double d_b_d_phi = nfp * b.dzeta[point];
            const double dd_b_d_phi = nfp * db.dzeta[point];
            const double b_dot =
                bu[point] * b.dtheta[point] + bv[point] * d_b_d_phi;
            const double db_dot =
                dbu[point] * b.dtheta[point] + bu[point] * db.dtheta[point] +
                dbv[point] * d_b_d_phi + bv[point] * dd_b_d_phi;
            const double cross_numerator =
                bcv[point] * b.dtheta[point] - bcu[point] * d_b_d_phi;
            const double dcross_numerator =
                dbcv[point] * b.dtheta[point] + bcv[point] * db.dtheta[point] -
                dbcu[point] * d_b_d_phi - bcu[point] * dd_b_d_phi;
            const double cross = cross_numerator / g[point];
            const double dcross =
                (dcross_numerator - cross * dg[point]) / g[point];
            const double helicity = static_cast<double>(spec.helicity_n) -
                                    iota * static_cast<double>(spec.helicity_m);
            const double dhelicity =
                -diota * static_cast<double>(spec.helicity_m);
            const double covariant =
                static_cast<double>(spec.helicity_m) * g_cov +
                static_cast<double>(spec.helicity_n) * i_cov;
            const double dcovariant =
                static_cast<double>(spec.helicity_m) * dg_cov +
                static_cast<double>(spec.helicity_n) * di_cov;
            const double numerator =
                helicity * flux_prime * cross + covariant * b_dot;
            const double dnumerator =
                (dhelicity * flux_prime + helicity * dflux_prime) * cross +
                helicity * flux_prime * dcross + dcovariant * b_dot +
                covariant * db_dot;
            const double abs_g = std::abs(g[point]);
            const double dabs_g = std::copysign(1.0, g[point]) * dg[point];
            const double weight = surface_weight * abs_g / jacobian_sum;
            const double dweight =
                weight * (dabs_g / abs_g - djacobian_sum / jacobian_sum);
            const double root_weight = std::sqrt(weight);
            const double field_strength = b.value[point];
            const double inverse_b_cubed =
                1.0 / (field_strength * field_strength * field_strength);
            const double residual = root_weight * numerator * inverse_b_cubed;
            const double dresidual =
                root_weight * dnumerator * inverse_b_cubed +
                0.5 * residual * dweight / weight -
                3.0 * residual * db.value[point] / field_strength;
            if (!std::isfinite(dresidual)) {
                throw std::runtime_error(
                    "flux-surface QS tangent produced a non-finite residual");
            }
            result.push_back(dresidual);
        }
    }
    return result;
}

inline std::vector<double> calculate_qh_target_jvp(
    const cumes::EquilibriumSnapshot& equilibrium,
    const cumes::EquilibriumProfiles& profiles,
    const cumes::EquilibriumTangent& tangent,
    const cumes::InputParams& input,
    const FluxSurfaceQuasisymmetryTargetSpec& spec) {
    std::vector<double> result = calculate_quasisymmetry_target_jvp(
        equilibrium, profiles, tangent, input.nfp, spec);
    result.push_back(
        calculate_plasma_size_jvp(equilibrium, tangent, input).aspect_ratio);
    return result;
}

inline std::vector<double> calculate_qa_target_jvp(
    const cumes::EquilibriumSnapshot& equilibrium,
    const cumes::EquilibriumProfiles& profiles,
    const cumes::EquilibriumTangent& tangent,
    const cumes::InputParams& input,
    const FluxSurfaceQuasisymmetryTargetSpec& spec,
    bool include_iota_residual) {
    std::vector<double> result = calculate_quasisymmetry_target_jvp(
        equilibrium, profiles, tangent, input.nfp, spec);
    result.push_back(
        calculate_plasma_size_jvp(equilibrium, tangent, input).aspect_ratio);
    if (include_iota_residual) {
        result.push_back(integrate_iota(tangent.profiles));
    }
    return result;
}

}  // namespace cumes_meow_example
