// Optimizer-owned parameterization of a stellarator-symmetric VMEC boundary.
#ifndef MEOW_CUMES_BOUNDARY_PARAMETERIZATION_HPP_
#define MEOW_CUMES_BOUNDARY_PARAMETERIZATION_HPP_

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cumes/config/problem_spec.hpp>
#include <meow/trf.hpp>

namespace cumes_meow_example {

enum class BoundaryFamily { RBC, ZBS };

struct BoundaryDegreeOfFreedom {
    BoundaryFamily family;
    int m;
    int n;
};

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
