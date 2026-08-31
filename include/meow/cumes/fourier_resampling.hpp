#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace cumes_meow_example::fourier_resampling {

struct PeriodicWeights {
    int source_size = 0;
    int target_size = 0;
    std::vector<double> value;
    std::vector<double> derivative;
};

inline PeriodicWeights make_periodic_weights(int source_size, int target_size) {
    if (source_size < 1 || target_size < 1) {
        throw std::invalid_argument("periodic grid sizes must be positive");
    }
    PeriodicWeights result;
    result.source_size = source_size;
    result.target_size = target_size;
    result.value.resize(static_cast<std::size_t>(source_size * target_size));
    result.derivative.resize(result.value.size());
    for (int target = 0; target < target_size; ++target) {
        const double x = 2.0 * std::numbers::pi * target / target_size;
        for (int source = 0; source < source_size; ++source) {
            const double y = 2.0 * std::numbers::pi * source / source_size;
            const double delta = x - y;
            double value = 1.0;
            double derivative = 0.0;
            const int paired_modes = (source_size - 1) / 2;
            for (int mode = 1; mode <= paired_modes; ++mode) {
                value += 2.0 * std::cos(mode * delta);
                derivative -= 2.0 * mode * std::sin(mode * delta);
            }
            if (source_size % 2 == 0) {
                const int nyquist = source_size / 2;
                value += std::cos(nyquist * delta);
                derivative -= nyquist * std::sin(nyquist * delta);
            }
            const std::size_t index =
                static_cast<std::size_t>(target * source_size + source);
            result.value[index] = value / source_size;
            result.derivative[index] = derivative / source_size;
        }
    }
    return result;
}

struct ResampledField {
    std::vector<double> value;
    std::vector<double> dtheta;
    std::vector<double> dzeta;
};

inline ResampledField resample_2d(const std::vector<double>& source,
                                  int source_ntheta,
                                  int source_nzeta,
                                  const PeriodicWeights& theta,
                                  const PeriodicWeights& zeta,
                                  bool derivatives) {
    if (source.size() !=
            static_cast<std::size_t>(source_ntheta * source_nzeta) ||
        theta.source_size != source_ntheta ||
        zeta.source_size != source_nzeta) {
        throw std::invalid_argument("Fourier resampling shape mismatch");
    }
    const int target_ntheta = theta.target_size;
    const int target_nzeta = zeta.target_size;
    std::vector<double> theta_value(
        static_cast<std::size_t>(target_ntheta * source_nzeta), 0.0);
    std::vector<double> theta_derivative;
    if (derivatives) theta_derivative.assign(theta_value.size(), 0.0);

    for (int source_zeta = 0; source_zeta < source_nzeta; ++source_zeta) {
        for (int target_theta = 0; target_theta < target_ntheta;
             ++target_theta) {
            double value = 0.0;
            double derivative = 0.0;
            for (int source_theta = 0; source_theta < source_ntheta;
                 ++source_theta) {
                const double sample = source[static_cast<std::size_t>(
                    source_theta + source_zeta * source_ntheta)];
                const std::size_t weight = static_cast<std::size_t>(
                    target_theta * source_ntheta + source_theta);
                value += theta.value[weight] * sample;
                if (derivatives) {
                    derivative += theta.derivative[weight] * sample;
                }
            }
            const std::size_t index = static_cast<std::size_t>(
                target_theta + source_zeta * target_ntheta);
            theta_value[index] = value;
            if (derivatives) theta_derivative[index] = derivative;
        }
    }

    ResampledField result;
    const std::size_t target_points =
        static_cast<std::size_t>(target_ntheta * target_nzeta);
    result.value.assign(target_points, 0.0);
    if (derivatives) {
        result.dtheta.assign(target_points, 0.0);
        result.dzeta.assign(target_points, 0.0);
    }
    for (int target_zeta = 0; target_zeta < target_nzeta; ++target_zeta) {
        for (int target_theta = 0; target_theta < target_ntheta;
             ++target_theta) {
            double value = 0.0;
            double dtheta = 0.0;
            double dzeta = 0.0;
            for (int source_zeta = 0; source_zeta < source_nzeta;
                 ++source_zeta) {
                const std::size_t source_index = static_cast<std::size_t>(
                    target_theta + source_zeta * target_ntheta);
                const std::size_t weight = static_cast<std::size_t>(
                    target_zeta * source_nzeta + source_zeta);
                value += zeta.value[weight] * theta_value[source_index];
                if (derivatives) {
                    dtheta +=
                        zeta.value[weight] * theta_derivative[source_index];
                    dzeta +=
                        zeta.derivative[weight] * theta_value[source_index];
                }
            }
            const std::size_t target_index = static_cast<std::size_t>(
                target_theta + target_zeta * target_ntheta);
            result.value[target_index] = value;
            if (derivatives) {
                result.dtheta[target_index] = dtheta;
                result.dzeta[target_index] = dzeta;
            }
        }
    }
    return result;
}

inline std::vector<double> interpolate_half_grid_field(
    const std::vector<double>& field,
    int ns,
    std::size_t points,
    double normalized_toroidal_flux) {
    const int half_surfaces = ns - 1;
    if (half_surfaces < 2 ||
        field.size() != static_cast<std::size_t>(half_surfaces) * points ||
        !std::isfinite(normalized_toroidal_flux) ||
        normalized_toroidal_flux < 0.0 || normalized_toroidal_flux > 1.0) {
        throw std::invalid_argument("invalid half-grid interpolation request");
    }
    const double half_index =
        normalized_toroidal_flux * static_cast<double>(ns - 1) - 0.5;
    int lower = static_cast<int>(std::floor(half_index));
    if (lower < 0) lower = 0;
    if (lower > half_surfaces - 2) lower = half_surfaces - 2;
    const double fraction = half_index - lower;
    std::vector<double> result(points);
    const std::size_t lower_offset = static_cast<std::size_t>(lower) * points;
    const std::size_t upper_offset = lower_offset + points;
    for (std::size_t point = 0; point < points; ++point) {
        result[point] = field[lower_offset + point] +
                        fraction * (field[upper_offset + point] -
                                    field[lower_offset + point]);
    }
    return result;
}

inline double interpolate_half_grid_profile(const std::vector<double>& profile,
                                            int ns,
                                            double normalized_toroidal_flux) {
    const auto value =
        interpolate_half_grid_field(profile, ns, 1, normalized_toroidal_flux);
    return value.front();
}

}  // namespace cumes_meow_example::fourier_resampling
