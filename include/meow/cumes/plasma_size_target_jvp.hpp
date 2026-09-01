#pragma once

#include "plasma_size_target.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

#include <cumes/solver/equilibrium_tangent.hpp>

namespace cumes_meow_example {

// Directional derivative of the VMEC-compatible size diagnostics. Target
// policy remains in meow; the input tangent is the target-independent public
// equilibrium sensitivity supplied by cuMES.
inline PlasmaSize calculate_plasma_size_jvp(
    const cumes::EquilibriumSnapshot& equilibrium,
    const cumes::EquilibriumTangent& tangent,
    const cumes::InputParams& input) {
    const auto primal = calculate_plasma_size(equilibrium, input);
    const auto& dequilibrium = tangent.equilibrium;
    if (dequilibrium.ns != equilibrium.ns ||
        dequilibrium.mnmax != equilibrium.mnmax ||
        dequilibrium.ntheta != equilibrium.ntheta ||
        dequilibrium.nzeta != equilibrium.nzeta) {
        throw std::invalid_argument("plasma-size tangent shape mismatch");
    }
    for (std::size_t component = 0;
         component < cumes::EquilibriumSnapshot::COUNT; ++component) {
        if (dequilibrium.families[component].size() !=
            equilibrium.families[component].size()) {
            throw std::invalid_argument("plasma-size tangent shape mismatch");
        }
    }

    constexpr double TWO_PI = 2.0 * std::numbers::pi;
    const int boundary = equilibrium.ns - 1;
    double cross_sum = 0.0;
    double volume_sum = 0.0;
    double dcross_sum = 0.0;
    double dvolume_sum = 0.0;
    for (int k = 0; k < equilibrium.nzeta; ++k) {
        const double zeta = TWO_PI * static_cast<double>(k) /
                            static_cast<double>(equilibrium.nzeta);
        for (int l = 0; l < equilibrium.ntheta; ++l) {
            const double theta = TWO_PI * static_cast<double>(l) /
                                 static_cast<double>(equilibrium.ntheta);
            double r = 0.0;
            double z_theta = 0.0;
            double dr = 0.0;
            double dz_theta = 0.0;
            for (int m = 0; m < input.mpol; ++m) {
                const double cos_m = std::cos(static_cast<double>(m) * theta);
                const double sin_m = std::sin(static_cast<double>(m) * theta);
                for (int n = 0; n <= input.ntor; ++n) {
                    const double cos_n =
                        std::cos(static_cast<double>(n) * zeta);
                    const double sin_n =
                        std::sin(static_cast<double>(n) * zeta);
                    const int mode = m * (input.ntor + 1) + n;
                    const std::size_t index =
                        static_cast<std::size_t>(mode) * equilibrium.ns +
                        boundary;
                    const double r_basis_cc = cos_m * cos_n;
                    const double r_basis_ss = sin_m * sin_n;
                    const double z_basis_sc =
                        static_cast<double>(m) * cos_m * cos_n;
                    const double z_basis_cs =
                        -static_cast<double>(m) * sin_m * sin_n;
                    r += equilibrium.families[cumes::EquilibriumSnapshot::RMNCC]
                                             [index] *
                             r_basis_cc +
                         equilibrium.families[cumes::EquilibriumSnapshot::RMNSS]
                                             [index] *
                             r_basis_ss;
                    dr +=
                        dequilibrium.families[cumes::EquilibriumSnapshot::RMNCC]
                                             [index] *
                            r_basis_cc +
                        dequilibrium.families[cumes::EquilibriumSnapshot::RMNSS]
                                             [index] *
                            r_basis_ss;
                    z_theta +=
                        equilibrium.families[cumes::EquilibriumSnapshot::ZMNSC]
                                            [index] *
                            z_basis_sc +
                        equilibrium.families[cumes::EquilibriumSnapshot::ZMNCS]
                                            [index] *
                            z_basis_cs;
                    dz_theta +=
                        dequilibrium.families[cumes::EquilibriumSnapshot::ZMNSC]
                                             [index] *
                            z_basis_sc +
                        dequilibrium.families[cumes::EquilibriumSnapshot::ZMNCS]
                                             [index] *
                            z_basis_cs;
                }
            }
            cross_sum += r * z_theta;
            volume_sum += r * r * z_theta;
            dcross_sum += dr * z_theta + r * dz_theta;
            dvolume_sum += 2.0 * r * dr * z_theta + r * r * dz_theta;
        }
    }

    const double point_count =
        static_cast<double>(equilibrium.ntheta) * equilibrium.nzeta;
    const double cross_raw = TWO_PI * cross_sum / point_count;
    const double volume_raw =
        2.0 * std::numbers::pi * std::numbers::pi * volume_sum / point_count;
    const double dcross =
        std::copysign(1.0, cross_raw) * TWO_PI * dcross_sum / point_count;
    const double dvolume = std::copysign(1.0, volume_raw) * 2.0 *
                           std::numbers::pi * std::numbers::pi * dvolume_sum /
                           point_count;

    PlasmaSize result;
    result.cross_section_area = dcross;
    result.volume = dvolume;
    result.major_radius =
        primal.major_radius *
        (dvolume / primal.volume - dcross / primal.cross_section_area);
    result.minor_radius =
        0.5 * primal.minor_radius * dcross / primal.cross_section_area;
    result.aspect_ratio =
        primal.aspect_ratio * (result.major_radius / primal.major_radius -
                               result.minor_radius / primal.minor_radius);
    return result;
}

}  // namespace cumes_meow_example
