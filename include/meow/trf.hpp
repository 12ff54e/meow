#pragma once

#include <cstddef>
#include <functional>
#include <limits>
#include <string>

#include <Eigen/Core>

namespace meow {

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

using ResidualFunction = std::function<Vector(const Vector&)>;
using JacobianFunction = std::function<Matrix(const Vector&)>;

struct Bounds {
    Vector lower;
    Vector upper;

    static Bounds unbounded(Eigen::Index size);
};

enum class TrfStatus : int {
    USER_STOPPED = -2,
    MAX_FUNCTION_EVALUATIONS = 0,
    GRADIENT_TOLERANCE = 1,
    FUNCTION_TOLERANCE = 2,
    STEP_TOLERANCE = 3,
    FUNCTION_AND_STEP_TOLERANCE = 4,
};

struct IterationInfo {
    std::size_t iteration = 0;
    std::size_t function_evaluations = 0;
    std::size_t jacobian_evaluations = 0;
    std::size_t jacobian_updates = 0;
    double cost = 0.0;
    double optimality = 0.0;
    double trust_radius = 0.0;
};

using IterationCallback =
    std::function<bool(const Vector&, const IterationInfo&)>;

struct TrfOptions {
    double ftol = 1e-8;
    double xtol = 1e-8;
    double gtol = 1e-8;

    // Zero selects 100 * number_of_variables. Numerical-Jacobian residual
    // calls are included in this limit.
    std::size_t max_function_evaluations = 0;

    // Forward-difference steps are
    // max(abs(x_j) * finite_difference_step,
    //     finite_difference_absolute_step).
    // If both are zero, the legacy/default rule
    // sqrt(machine epsilon) * max(abs(x_j), 1) is used.
    double finite_difference_step = 0.0;
    double finite_difference_absolute_step = 0.0;

    // Positive characteristic scales for x. An empty vector means all ones.
    Vector x_scale;

    // Derive characteristic variable scales from reciprocal Jacobian column
    // norms, following the standard iterative least-squares scaling rule.
    // This option and an explicit x_scale are mutually exclusive.
    bool scale_from_jacobian = false;

    // One preserves the default behavior of rebuilding the Jacobian after
    // every accepted step. Values above one permit good-Broyden secant
    // updates between exact rebuilds. A rebuild is forced earlier when the
    // trust-region reduction ratio or relative secant defect crosses the
    // safeguards below.
    std::size_t jacobian_refresh_interval = 1;
    double broyden_min_reduction_ratio = 0.1;
    double broyden_max_secant_error = 0.5;

    // Called after each accepted iteration. Return false to stop cleanly.
    IterationCallback callback;
    int verbose = 0;
};

struct TrfResult {
    Vector x;
    Vector residual;
    Matrix jacobian;
    Vector gradient;
    Eigen::VectorXi active_mask;

    double cost = std::numeric_limits<double>::quiet_NaN();
    double optimality = std::numeric_limits<double>::quiet_NaN();
    std::size_t function_evaluations = 0;
    std::size_t jacobian_evaluations = 0;
    std::size_t jacobian_updates = 0;
    std::size_t iterations = 0;
    TrfStatus status = TrfStatus::MAX_FUNCTION_EVALUATIONS;
    bool success = false;
    std::string message;
};

// Minimize 0.5 * ||residual(x)||^2 subject to lower <= x <= upper with a
// trust-region reflective method. If jacobian is empty, a bound-aware forward
// finite-difference Jacobian is used.
TrfResult trf_least_squares(const ResidualFunction& residual,
                            const Vector& initial_x,
                            const Bounds& bounds,
                            const TrfOptions& options = {},
                            const JacobianFunction& jacobian = {});

TrfResult trf_least_squares(const ResidualFunction& residual,
                            const Vector& initial_x,
                            const TrfOptions& options = {},
                            const JacobianFunction& jacobian = {});

}  // namespace meow
