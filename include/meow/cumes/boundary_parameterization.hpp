// Optimizer-owned parameterization of a stellarator-symmetric VMEC boundary.
#ifndef MEOW_CUMES_BOUNDARY_PARAMETERIZATION_HPP_
#define MEOW_CUMES_BOUNDARY_PARAMETERIZATION_HPP_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cumes/config/problem_spec.hpp>
#include <cumes/io/equilibrium_snapshot.hpp>
#include <meow/trf.hpp>

namespace cumes_meow_example {

enum class BoundaryFamily { RBC, ZBS };

struct BoundaryDegreeOfFreedom {
    BoundaryFamily family;
    int m;
    int n;
};

inline std::vector<double> boundary_centerline_coefficients(
    const std::vector<cumes::BoundaryHarmonic>& harmonics,
    int ntor) {
    const std::size_t axis_size =
        static_cast<std::size_t>(std::max(ntor, 0) + 1);
    std::vector<double> coefficients(axis_size, 0.0);
    for (const auto& harmonic : harmonics) {
        if (harmonic.m == 0 && harmonic.n >= 0 && harmonic.n <= ntor) {
            coefficients[static_cast<std::size_t>(harmonic.n)] = harmonic.value;
        }
    }
    return coefficients;
}

// cuMES uses these arrays only to seed the cold-start interior surfaces.  The
// analytic Landreman inputs ask VMEC to choose this predictor automatically,
// so construction runs approximate that behavior with the current boundary
// centerline instead of retaining a stale iteration-zero axis.
inline void refresh_axis_predictor_from_boundary_centerline(
    cumes::ProblemSpec& problem) {
    const std::size_t axis_size =
        static_cast<std::size_t>(std::max(problem.ntor, 0) + 1);
    problem.raxis_c =
        boundary_centerline_coefficients(problem.rbc, problem.ntor);
    problem.zaxis_s =
        boundary_centerline_coefficients(problem.zbs, problem.ntor);
    problem.raxis_c.resize(axis_size);
    problem.zaxis_s.resize(axis_size);
    problem.has_raxis_c = true;
    problem.has_zaxis_s = true;
}

// Keep all trials in one numerical-Jacobian evaluation tied to the same
// accepted equilibrium axis, while following the trial's centerline change.
inline void track_axis_predictor_from_accepted_boundary(
    cumes::ProblemSpec& problem,
    const cumes::ProblemSpec& accepted) {
    const std::size_t axis_size =
        static_cast<std::size_t>(std::max(problem.ntor, 0) + 1);
    const std::vector<double> trial_r =
        boundary_centerline_coefficients(problem.rbc, problem.ntor);
    const std::vector<double> trial_z =
        boundary_centerline_coefficients(problem.zbs, problem.ntor);
    const std::vector<double> accepted_r =
        boundary_centerline_coefficients(accepted.rbc, problem.ntor);
    const std::vector<double> accepted_z =
        boundary_centerline_coefficients(accepted.zbs, problem.ntor);
    problem.raxis_c = accepted.raxis_c;
    problem.zaxis_s = accepted.zaxis_s;
    problem.raxis_c.resize(axis_size, 0.0);
    problem.zaxis_s.resize(axis_size, 0.0);
    for (std::size_t index = 0; index < axis_size; ++index) {
        problem.raxis_c[index] += trial_r[index] - accepted_r[index];
        problem.zaxis_s[index] += trial_z[index] - accepted_z[index];
    }
    problem.has_raxis_c = true;
    problem.has_zaxis_s = true;
}

inline void refresh_axis_predictor_from_equilibrium(
    cumes::ProblemSpec& problem,
    const cumes::EquilibriumSnapshot& equilibrium) {
    const std::size_t axis_size =
        static_cast<std::size_t>(std::max(problem.ntor, 0) + 1);
    if (equilibrium.ns <= 0 ||
        equilibrium.mnmax < static_cast<int>(axis_size) ||
        equilibrium.component(cumes::EquilibriumSnapshot::RMNCC).size() !=
            equilibrium.family_size() ||
        equilibrium.component(cumes::EquilibriumSnapshot::ZMNCS).size() !=
            equilibrium.family_size()) {
        throw std::invalid_argument(
            "complete equilibrium state is required for the axis predictor");
    }
    problem.raxis_c.assign(axis_size, 0.0);
    problem.zaxis_s.assign(axis_size, 0.0);
    const std::size_t ns = static_cast<std::size_t>(equilibrium.ns);
    for (std::size_t mode = 0; mode < axis_size; ++mode) {
        problem.raxis_c[mode] =
            equilibrium.component(cumes::EquilibriumSnapshot::RMNCC)[mode * ns];
        problem.zaxis_s[mode] = -equilibrium.component(
            cumes::EquilibriumSnapshot::ZMNCS)[mode * ns];
    }
    problem.has_raxis_c = true;
    problem.has_zaxis_s = true;
}

// Matches the SIMSOPT fixed_range(mmin=0, mmax=max_mode,
// nmin=-max_mode, nmax=max_mode) convention used in the Landreman-Paul
// drivers. The redundant negative-n m=0 modes are omitted, as are the fixed
// rbc(0,0) major-radius coefficient and identically-zero zbs(0,0) mode.
class StellaratorSymmetricBoundaryParameterization {
   public:
    explicit StellaratorSymmetricBoundaryParameterization(int max_mode)
        : max_mode_(max_mode), degrees_of_freedom_(make_modes(max_mode)) {}

    int max_mode() const { return max_mode_; }

    std::size_t size() const { return degrees_of_freedom_.size(); }

    const std::vector<BoundaryDegreeOfFreedom>& degrees_of_freedom() const {
        return degrees_of_freedom_;
    }

    meow::Vector values(const cumes::ProblemSpec& problem) const {
        meow::Vector result(static_cast<Eigen::Index>(size()));
        for (std::size_t index = 0; index < size(); ++index) {
            const auto& mode = degrees_of_freedom_[index];
            result[static_cast<Eigen::Index>(index)] =
                coefficient(problem, mode);
        }
        return result;
    }

    // x contains the absolute Fourier coefficients, not coefficient deltas.
    // This is needed for SIMSOPT's relative finite-difference convention.
    cumes::ProblemSpec apply(const cumes::ProblemSpec& baseline,
                             const meow::Vector& x) const {
        if (x.size() != static_cast<Eigen::Index>(size()) ||
            !x.array().isFinite().all()) {
            throw std::invalid_argument(
                "boundary parameterization expects one finite value per "
                "degree of freedom");
        }

        cumes::ProblemSpec result = baseline;
        for (std::size_t index = 0; index < size(); ++index) {
            const auto& mode = degrees_of_freedom_[index];
            set_coefficient(result, mode, x[static_cast<Eigen::Index>(index)]);
        }
        return result;
    }

    std::string name(std::size_t index) const {
        if (index >= size()) {
            throw std::out_of_range("boundary degree-of-freedom index");
        }
        const auto& mode = degrees_of_freedom_[index];
        const char* family = mode.family == BoundaryFamily::RBC ? "rbc" : "zbs";
        return std::string(family) + "(" + std::to_string(mode.m) + "," +
               std::to_string(mode.n) + ")";
    }

   private:
    static std::vector<BoundaryDegreeOfFreedom> make_modes(int max_mode) {
        if (max_mode < 1) {
            throw std::invalid_argument(
                "boundary max_mode must be at least one");
        }
        std::vector<BoundaryDegreeOfFreedom> modes;
        const std::size_t modes_per_family =
            static_cast<std::size_t>(2 * max_mode * (max_mode + 1));
        modes.reserve(2 * modes_per_family);
        for (BoundaryFamily family :
             {BoundaryFamily::RBC, BoundaryFamily::ZBS}) {
            for (int n = 1; n <= max_mode; ++n) {
                modes.push_back({family, 0, n});
            }
            for (int m = 1; m <= max_mode; ++m) {
                for (int n = -max_mode; n <= max_mode; ++n) {
                    modes.push_back({family, m, n});
                }
            }
        }
        return modes;
    }

    static const std::vector<cumes::BoundaryHarmonic>& family(
        const cumes::ProblemSpec& problem,
        BoundaryFamily selected) {
        return selected == BoundaryFamily::RBC ? problem.rbc : problem.zbs;
    }

    static std::vector<cumes::BoundaryHarmonic>& family(
        cumes::ProblemSpec& problem,
        BoundaryFamily selected) {
        return selected == BoundaryFamily::RBC ? problem.rbc : problem.zbs;
    }

    static double coefficient(const cumes::ProblemSpec& problem,
                              const BoundaryDegreeOfFreedom& mode) {
        for (const auto& harmonic : family(problem, mode.family)) {
            if (harmonic.m == mode.m && harmonic.n == mode.n) {
                return harmonic.value;
            }
        }
        return 0.0;
    }

    static void set_coefficient(cumes::ProblemSpec& problem,
                                const BoundaryDegreeOfFreedom& mode,
                                double value) {
        auto& harmonics = family(problem, mode.family);
        for (auto& harmonic : harmonics) {
            if (harmonic.m == mode.m && harmonic.n == mode.n) {
                harmonic.value = value;
                return;
            }
        }
        harmonics.push_back({mode.m, mode.n, value});
    }

    int max_mode_;
    std::vector<BoundaryDegreeOfFreedom> degrees_of_freedom_;
};

}  // namespace cumes_meow_example

#endif  // MEOW_CUMES_BOUNDARY_PARAMETERIZATION_HPP_
