#include "meow/trf.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <Eigen/Eigenvalues>

namespace meow {
namespace {

constexpr double EPSILON = std::numeric_limits<double>::epsilon();

struct Scaling {
    Vector v;
    Vector derivative;
};

struct StepChoice {
    Vector step;
    Vector scaled_step;
    double predicted_reduction = 0.0;
};

struct TrustUpdate {
    double radius = 0.0;
    double ratio = 0.0;
};

bool all_finite(const Vector& value) {
    return value.array().isFinite().all();
}

bool in_bounds(const Vector& x, const Vector& lower, const Vector& upper) {
    return (x.array() >= lower.array()).all() &&
           (x.array() <= upper.array()).all();
}

double infinity_norm(const Vector& value) {
    return value.size() == 0 ? 0.0 : value.cwiseAbs().maxCoeff();
}

Scaling coleman_li_scaling(const Vector& x,
                           const Vector& gradient,
                           const Vector& lower,
                           const Vector& upper) {
    Scaling scaling{Vector::Ones(x.size()), Vector::Zero(x.size())};
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        if (gradient[i] < 0.0 && std::isfinite(upper[i])) {
            scaling.v[i] = upper[i] - x[i];
            scaling.derivative[i] = -1.0;
        } else if (gradient[i] > 0.0 && std::isfinite(lower[i])) {
            scaling.v[i] = x[i] - lower[i];
            scaling.derivative[i] = 1.0;
        }
    }
    return scaling;
}

Vector make_strictly_feasible(Vector x,
                              const Vector& lower,
                              const Vector& upper,
                              double relative_step) {
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        if (std::isfinite(lower[i]) && x[i] <= lower[i]) {
            const double shifted =
                relative_step > 0.0
                    ? lower[i] +
                          relative_step * std::max(1.0, std::abs(lower[i]))
                    : std::nextafter(lower[i], upper[i]);
            x[i] = std::min(shifted, upper[i]);
        }
        if (std::isfinite(upper[i]) && x[i] >= upper[i]) {
            const double shifted =
                relative_step > 0.0
                    ? upper[i] -
                          relative_step * std::max(1.0, std::abs(upper[i]))
                    : std::nextafter(upper[i], lower[i]);
            x[i] = std::max(shifted, lower[i]);
        }
    }
    return x;
}

std::pair<double, Eigen::VectorXi> step_to_bound(const Vector& x,
                                                 const Vector& step,
                                                 const Vector& lower,
                                                 const Vector& upper) {
    Vector distances =
        Vector::Constant(x.size(), std::numeric_limits<double>::infinity());
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        if (step[i] > 0.0 && std::isfinite(upper[i])) {
            distances[i] = (upper[i] - x[i]) / step[i];
        } else if (step[i] < 0.0 && std::isfinite(lower[i])) {
            distances[i] = (lower[i] - x[i]) / step[i];
        }
    }

    const double distance = distances.minCoeff();
    Eigen::VectorXi hits = Eigen::VectorXi::Zero(x.size());
    if (std::isfinite(distance)) {
        const double tolerance =
            16.0 * EPSILON * std::max(1.0, std::abs(distance));
        for (Eigen::Index i = 0; i < x.size(); ++i) {
            if (std::abs(distances[i] - distance) <= tolerance) {
                hits[i] = step[i] > 0.0 ? 1 : -1;
            }
        }
    }
    return {distance, hits};
}

double positive_trust_intersection(const Vector& point,
                                   const Vector& direction,
                                   double radius) {
    const double a = direction.squaredNorm();
    if (a == 0.0) { return 0.0; }
    const double b = point.dot(direction);
    const double c = point.squaredNorm() - radius * radius;
    const double discriminant = std::max(0.0, b * b - a * c);
    return (-b + std::sqrt(discriminant)) / a;
}

double quadratic_value(const Matrix& jacobian,
                       const Vector& gradient,
                       const Vector& diagonal,
                       const Vector& step) {
    return 0.5 * (jacobian * step).squaredNorm() + gradient.dot(step) +
           0.5 * diagonal.dot(step.cwiseProduct(step));
}

std::pair<double, double> minimize_quadratic_line(
    const Matrix& jacobian,
    const Vector& gradient,
    const Vector& diagonal,
    const Vector& direction,
    double lower,
    double upper,
    const Vector* origin = nullptr) {
    Vector zero;
    const Vector& base =
        origin == nullptr ? (zero = Vector::Zero(direction.size())) : *origin;
    const Vector j_direction = jacobian * direction;
    const Vector j_base = jacobian * base;
    const double a = 0.5 * (j_direction.squaredNorm() +
                            diagonal.dot(direction.cwiseProduct(direction)));
    const double b = j_direction.dot(j_base) + gradient.dot(direction) +
                     diagonal.dot(base.cwiseProduct(direction));

    double candidate = lower;
    if (a > 0.0) {
        candidate = std::clamp(-b / (2.0 * a), lower, upper);
    } else if (b < 0.0) {
        candidate = upper;
    }

    const double lower_value =
        quadratic_value(jacobian, gradient, diagonal, base + lower * direction);
    const double upper_value =
        quadratic_value(jacobian, gradient, diagonal, base + upper * direction);
    double value = quadratic_value(jacobian, gradient, diagonal,
                                   base + candidate * direction);
    if (lower_value < value) {
        candidate = lower;
        value = lower_value;
    }
    if (upper_value < value) {
        candidate = upper;
        value = upper_value;
    }
    return {candidate, value};
}

Vector solve_trust_region_quadratic(const Matrix& hessian,
                                    const Vector& gradient,
                                    double radius) {
    if (gradient.norm() == 0.0 || radius == 0.0) {
        return Vector::Zero(gradient.size());
    }

    Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver(hessian);
    if (eigen_solver.info() != Eigen::Success) {
        throw std::runtime_error("TRF trust-region eigensolve failed");
    }
    Vector eigenvalues = eigen_solver.eigenvalues();
    const Matrix eigenvectors = eigen_solver.eigenvectors();
    const Vector transformed_gradient = eigenvectors.transpose() * gradient;
    const double spectral_scale =
        std::max(1.0, eigenvalues.cwiseAbs().maxCoeff());
    const double zero_tolerance = 64.0 * EPSILON * spectral_scale;
    for (Eigen::Index i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues[i] < 0.0 && eigenvalues[i] > -zero_tolerance) {
            eigenvalues[i] = 0.0;
        }
        if (eigenvalues[i] < 0.0) {
            throw std::runtime_error(
                "TRF model Hessian is not positive semidefinite");
        }
    }

    auto transformed_step = [&](double damping) {
        Vector value(transformed_gradient.size());
        for (Eigen::Index i = 0; i < value.size(); ++i) {
            const double denominator = eigenvalues[i] + damping;
            if (denominator <= zero_tolerance) {
                value[i] =
                    transformed_gradient[i] == 0.0
                        ? 0.0
                        : std::copysign(std::numeric_limits<double>::infinity(),
                                        -transformed_gradient[i]);
            } else {
                value[i] = -transformed_gradient[i] / denominator;
            }
        }
        return value;
    };

    Vector spectral_step = transformed_step(0.0);
    if (all_finite(spectral_step) && spectral_step.norm() <= radius) {
        return eigenvectors * spectral_step;
    }

    double lower = 0.0;
    double upper = std::max(1.0, gradient.norm() / radius);
    for (int expansion = 0; expansion < 100; ++expansion) {
        spectral_step = transformed_step(upper);
        if (all_finite(spectral_step) && spectral_step.norm() <= radius) {
            break;
        }
        upper *= 2.0;
    }
    for (int iteration = 0; iteration < 80; ++iteration) {
        const double damping = std::midpoint(lower, upper);
        spectral_step = transformed_step(damping);
        if (!all_finite(spectral_step) || spectral_step.norm() > radius) {
            lower = damping;
        } else {
            upper = damping;
        }
    }
    return eigenvectors * transformed_step(upper);
}

StepChoice choose_reflective_step(const Vector& x,
                                  const Matrix& scaled_jacobian,
                                  const Vector& model_diagonal,
                                  const Vector& scaled_gradient,
                                  Vector step,
                                  Vector scaled_step,
                                  const Vector& transform,
                                  double radius,
                                  const Vector& lower,
                                  const Vector& upper,
                                  double theta) {
    if (in_bounds(x + step, lower, upper)) {
        const double value = quadratic_value(scaled_jacobian, scaled_gradient,
                                             model_diagonal, scaled_step);
        return {std::move(step), std::move(scaled_step), -value};
    }

    const auto [bound_stride, hits] = step_to_bound(x, step, lower, upper);
    Vector reflected_scaled = scaled_step;
    for (Eigen::Index i = 0; i < hits.size(); ++i) {
        if (hits[i] != 0) { reflected_scaled[i] *= -1.0; }
    }
    const Vector reflected = transform.cwiseProduct(reflected_scaled);

    step *= bound_stride;
    scaled_step *= bound_stride;
    const Vector x_on_bound = x + step;

    const double to_trust =
        positive_trust_intersection(scaled_step, reflected_scaled, radius);
    const auto [to_bound, ignored_hits] =
        step_to_bound(x_on_bound, reflected, lower, upper);
    static_cast<void>(ignored_hits);
    const double reflected_limit = std::min(to_bound, to_trust);

    Vector reflection_candidate;
    Vector reflection_scaled_candidate;
    double reflection_value = std::numeric_limits<double>::infinity();
    if (reflected_limit > 0.0 && std::isfinite(reflected_limit)) {
        const double reflected_lower =
            (1.0 - theta) * bound_stride / reflected_limit;
        const double reflected_upper =
            reflected_limit == to_bound ? theta * to_bound : to_trust;
        if (reflected_lower <= reflected_upper) {
            const auto [stride, value] = minimize_quadratic_line(
                scaled_jacobian, scaled_gradient, model_diagonal,
                reflected_scaled, reflected_lower, reflected_upper,
                &scaled_step);
            reflection_scaled_candidate =
                scaled_step + stride * reflected_scaled;
            reflection_candidate =
                transform.cwiseProduct(reflection_scaled_candidate);
            reflection_value = value;
        }
    }

    step *= theta;
    scaled_step *= theta;
    const double truncated_value = quadratic_value(
        scaled_jacobian, scaled_gradient, model_diagonal, scaled_step);

    Vector cauchy_scaled = -scaled_gradient;
    Vector cauchy = transform.cwiseProduct(cauchy_scaled);
    double cauchy_value = std::numeric_limits<double>::infinity();
    if (cauchy_scaled.norm() > 0.0) {
        const double to_cauchy_trust = radius / cauchy_scaled.norm();
        const auto [to_cauchy_bound, cauchy_hits] =
            step_to_bound(x, cauchy, lower, upper);
        static_cast<void>(cauchy_hits);
        const double cauchy_limit = to_cauchy_bound < to_cauchy_trust
                                        ? theta * to_cauchy_bound
                                        : to_cauchy_trust;
        const auto [stride, value] = minimize_quadratic_line(
            scaled_jacobian, scaled_gradient, model_diagonal, cauchy_scaled,
            0.0, cauchy_limit);
        cauchy_scaled *= stride;
        cauchy *= stride;
        cauchy_value = value;
    }

    if (truncated_value <= reflection_value &&
        truncated_value <= cauchy_value) {
        return {std::move(step), std::move(scaled_step), -truncated_value};
    }
    if (reflection_value <= cauchy_value) {
        return {std::move(reflection_candidate),
                std::move(reflection_scaled_candidate), -reflection_value};
    }
    return {std::move(cauchy), std::move(cauchy_scaled), -cauchy_value};
}

TrustUpdate update_trust_radius(double radius,
                                double actual_reduction,
                                double predicted_reduction,
                                double step_norm,
                                bool reached_boundary) {
    double ratio = 0.0;
    if (predicted_reduction > 0.0) {
        ratio = actual_reduction / predicted_reduction;
    } else if (actual_reduction == 0.0 && predicted_reduction == 0.0) {
        ratio = 1.0;
    }

    if (ratio < 0.25) {
        radius = 0.25 * step_norm;
    } else if (ratio > 0.75 && reached_boundary) {
        radius *= 2.0;
    }
    return {radius, ratio};
}

TrfStatus termination_status(double actual_reduction,
                             double old_cost,
                             double step_norm,
                             double x_norm,
                             double ratio,
                             const TrfOptions& options) {
    const bool function_done = options.ftol > 0.0 &&
                               actual_reduction < options.ftol * old_cost &&
                               ratio > 0.25;
    const bool step_done = options.xtol > 0.0 &&
                           step_norm < options.xtol * (options.xtol + x_norm);
    if (function_done && step_done) {
        return TrfStatus::FUNCTION_AND_STEP_TOLERANCE;
    }
    if (function_done) { return TrfStatus::FUNCTION_TOLERANCE; }
    if (step_done) { return TrfStatus::STEP_TOLERANCE; }
    return TrfStatus::MAX_FUNCTION_EVALUATIONS;
}

std::string status_message(TrfStatus status) {
    switch (status) {
        case TrfStatus::USER_STOPPED:
            return "Stopped by the iteration callback.";
        case TrfStatus::MAX_FUNCTION_EVALUATIONS:
            return "Maximum number of residual evaluations reached.";
        case TrfStatus::GRADIENT_TOLERANCE:
            return "The scaled gradient satisfies gtol.";
        case TrfStatus::FUNCTION_TOLERANCE:
            return "The cost reduction satisfies ftol.";
        case TrfStatus::STEP_TOLERANCE:
            return "The step size satisfies xtol.";
        case TrfStatus::FUNCTION_AND_STEP_TOLERANCE:
            return "Both ftol and xtol are satisfied.";
    }
    return "Unknown termination status.";
}

Eigen::VectorXi active_constraints(const Vector& x,
                                   const Vector& lower,
                                   const Vector& upper,
                                   double tolerance) {
    Eigen::VectorXi active = Eigen::VectorXi::Zero(x.size());
    for (Eigen::Index i = 0; i < x.size(); ++i) {
        const double lower_distance =
            std::isfinite(lower[i]) ? x[i] - lower[i]
                                    : std::numeric_limits<double>::infinity();
        const double upper_distance =
            std::isfinite(upper[i]) ? upper[i] - x[i]
                                    : std::numeric_limits<double>::infinity();
        const double lower_threshold =
            tolerance * std::max(1.0, std::abs(lower[i]));
        const double upper_threshold =
            tolerance * std::max(1.0, std::abs(upper[i]));
        if (lower_distance <= lower_threshold &&
            lower_distance <= upper_distance) {
            active[i] = -1;
        } else if (upper_distance <= upper_threshold) {
            active[i] = 1;
        }
    }
    return active;
}

void validate_inputs(const ResidualFunction& residual,
                     const Vector& initial_x,
                     const Bounds& bounds,
                     const TrfOptions& options) {
    if (!residual) {
        throw std::invalid_argument("TRF residual callback is empty");
    }
    if (initial_x.size() == 0 || !all_finite(initial_x)) {
        throw std::invalid_argument(
            "TRF initial_x must be nonempty and finite");
    }
    if (bounds.lower.size() != initial_x.size() ||
        bounds.upper.size() != initial_x.size()) {
        throw std::invalid_argument("TRF bounds must match initial_x size");
    }
    if (!(bounds.lower.array() < bounds.upper.array()).all()) {
        throw std::invalid_argument("TRF requires lower bounds < upper bounds");
    }
    if (!in_bounds(initial_x, bounds.lower, bounds.upper)) {
        throw std::invalid_argument("TRF initial_x is outside the bounds");
    }
    if (options.ftol < 0.0 || options.xtol < 0.0 || options.gtol < 0.0 ||
        (options.ftol == 0.0 && options.xtol == 0.0 && options.gtol == 0.0)) {
        throw std::invalid_argument(
            "TRF tolerances must be nonnegative and not all zero");
    }
    if (options.finite_difference_step < 0.0) {
        throw std::invalid_argument(
            "TRF finite_difference_step must be nonnegative");
    }
    if (options.x_scale.size() != 0 &&
        options.x_scale.size() != initial_x.size()) {
        throw std::invalid_argument(
            "TRF x_scale must be empty or match initial_x size");
    }
    if (options.x_scale.size() != 0 &&
        (!(options.x_scale.array() > 0.0).all() ||
         !all_finite(options.x_scale))) {
        throw std::invalid_argument(
            "TRF x_scale entries must be positive and finite");
    }
}

}  // namespace

Bounds Bounds::unbounded(Eigen::Index size) {
    return {Vector::Constant(size, -std::numeric_limits<double>::infinity()),
            Vector::Constant(size, std::numeric_limits<double>::infinity())};
}

TrfResult trf_least_squares(const ResidualFunction& residual_function,
                            const Vector& initial_x,
                            const Bounds& bounds,
                            const TrfOptions& options,
                            const JacobianFunction& jacobian_function) {
    validate_inputs(residual_function, initial_x, bounds, options);

    const Eigen::Index variable_count = initial_x.size();
    const Vector scale = options.x_scale.size() == 0
                             ? Vector::Ones(variable_count)
                             : options.x_scale;
    const Vector inverse_scale = scale.cwiseInverse();
    const std::size_t max_evaluations =
        options.max_function_evaluations == 0
            ? static_cast<std::size_t>(100 * variable_count)
            : options.max_function_evaluations;
    const std::size_t numerical_jacobian_cost =
        jacobian_function ? 0U : static_cast<std::size_t>(variable_count);
    if (max_evaluations < 1U + numerical_jacobian_cost) {
        throw std::invalid_argument(
            "TRF max_function_evaluations is too small for the initial "
            "Jacobian");
    }

    TrfResult result;
    Vector x =
        make_strictly_feasible(initial_x, bounds.lower, bounds.upper, 1e-10);

    auto evaluate = [&](const Vector& point) {
        Vector value = residual_function(point);
        ++result.function_evaluations;
        if (value.size() == 0) {
            throw std::invalid_argument(
                "TRF residual callback returned an empty vector");
        }
        return value;
    };

    Vector f = evaluate(x);
    if (!all_finite(f)) {
        throw std::invalid_argument(
            "TRF residuals are not finite at initial_x");
    }
    const Eigen::Index residual_count = f.size();

    auto evaluate_jacobian = [&](const Vector& point, const Vector& value) {
        Matrix jacobian;
        if (jacobian_function) {
            jacobian = jacobian_function(point);
        } else {
            jacobian.resize(residual_count, variable_count);
            const double relative_step = options.finite_difference_step > 0.0
                                             ? options.finite_difference_step
                                             : std::sqrt(EPSILON);
            for (Eigen::Index column = 0; column < variable_count; ++column) {
                double step =
                    relative_step * std::max(1.0, std::abs(point[column]));
                if (point[column] + step > bounds.upper[column]) {
                    step = -step;
                }
                if (point[column] + step < bounds.lower[column]) {
                    const double forward_room =
                        bounds.upper[column] - point[column];
                    const double backward_room =
                        point[column] - bounds.lower[column];
                    step = forward_room >= backward_room ? 0.5 * forward_room
                                                         : -0.5 * backward_room;
                }
                if (step == 0.0 || !std::isfinite(step)) {
                    throw std::runtime_error(
                        "TRF could not construct a finite-difference step");
                }
                Vector perturbed = point;
                perturbed[column] += step;
                const Vector perturbed_value = evaluate(perturbed);
                if (perturbed_value.size() != residual_count ||
                    !all_finite(perturbed_value)) {
                    throw std::runtime_error(
                        "TRF numerical Jacobian evaluation returned invalid "
                        "residuals");
                }
                jacobian.col(column) = (perturbed_value - value) / step;
            }
        }
        ++result.jacobian_evaluations;
        if (jacobian.rows() != residual_count ||
            jacobian.cols() != variable_count) {
            throw std::invalid_argument("TRF Jacobian has the wrong shape");
        }
        if (!jacobian.array().isFinite().all()) {
            throw std::invalid_argument(
                "TRF Jacobian contains non-finite values");
        }
        return jacobian;
    };

    Matrix jacobian = evaluate_jacobian(x, f);
    double cost = 0.5 * f.squaredNorm();
    Vector gradient = jacobian.transpose() * f;

    Scaling initial_scaling =
        coleman_li_scaling(x, gradient, bounds.lower, bounds.upper);
    for (Eigen::Index i = 0; i < variable_count; ++i) {
        if (initial_scaling.derivative[i] != 0.0) {
            initial_scaling.v[i] *= inverse_scale[i];
        }
    }
    double trust_radius = (x.cwiseProduct(inverse_scale)
                               .cwiseQuotient(initial_scaling.v.cwiseSqrt()))
                              .norm();
    if (trust_radius == 0.0 || !std::isfinite(trust_radius)) {
        trust_radius = 1.0;
    }

    TrfStatus status = TrfStatus::MAX_FUNCTION_EVALUATIONS;
    bool terminated = false;
    if (options.verbose > 0) {
        std::cout << " iter   nfev              cost        optimality       "
                     "radius\n";
    }

    while (!terminated) {
        Scaling scaling =
            coleman_li_scaling(x, gradient, bounds.lower, bounds.upper);
        const double optimality =
            infinity_norm(gradient.cwiseProduct(scaling.v));

        if (options.verbose > 0) {
            std::cout << std::setw(5) << result.iterations << std::setw(7)
                      << result.function_evaluations << std::setw(19)
                      << std::setprecision(10) << cost << std::setw(18)
                      << optimality << std::setw(13) << trust_radius << '\n';
        }
        if (optimality < options.gtol) {
            status = TrfStatus::GRADIENT_TOLERANCE;
            break;
        }

        // Reserve enough residual calls to form a new numerical Jacobian if
        // the trial point is accepted.
        if (result.function_evaluations + 1U + numerical_jacobian_cost >
            max_evaluations) {
            break;
        }

        for (Eigen::Index i = 0; i < variable_count; ++i) {
            if (scaling.derivative[i] != 0.0) {
                scaling.v[i] *= inverse_scale[i];
            }
        }
        const Vector transform = scaling.v.cwiseSqrt().cwiseProduct(scale);
        const Vector model_diagonal =
            gradient.cwiseProduct(scaling.derivative).cwiseProduct(scale);
        const Vector scaled_gradient = transform.cwiseProduct(gradient);
        const Matrix scaled_jacobian = jacobian * transform.asDiagonal();
        Matrix model_hessian = scaled_jacobian.transpose() * scaled_jacobian;
        model_hessian.diagonal() += model_diagonal;
        const double theta = std::max(0.995, 1.0 - optimality);

        bool accepted = false;
        double actual_reduction = 0.0;
        double step_norm = 0.0;
        double ratio = 0.0;
        Vector accepted_x;
        Vector accepted_f;
        double accepted_cost = cost;

        while (!accepted &&
               result.function_evaluations + 1U + numerical_jacobian_cost <=
                   max_evaluations) {
            const Vector raw_scaled_step = solve_trust_region_quadratic(
                model_hessian, scaled_gradient, trust_radius);
            const Vector raw_step = transform.cwiseProduct(raw_scaled_step);
            StepChoice choice = choose_reflective_step(
                x, scaled_jacobian, model_diagonal, scaled_gradient, raw_step,
                raw_scaled_step, transform, trust_radius, bounds.lower,
                bounds.upper, theta);

            if (!(choice.predicted_reduction > 0.0) ||
                choice.scaled_step.norm() == 0.0) {
                trust_radius *= 0.25;
                if (trust_radius <= EPSILON * (1.0 + x.norm())) {
                    status = TrfStatus::STEP_TOLERANCE;
                    terminated = true;
                    break;
                }
                continue;
            }

            Vector trial_x = make_strictly_feasible(
                x + choice.step, bounds.lower, bounds.upper, 0.0);
            Vector trial_f = evaluate(trial_x);
            const double scaled_step_norm = choice.scaled_step.norm();
            if (trial_f.size() != residual_count || !all_finite(trial_f)) {
                trust_radius = 0.25 * scaled_step_norm;
                continue;
            }

            const double trial_cost = 0.5 * trial_f.squaredNorm();
            actual_reduction = cost - trial_cost;
            const TrustUpdate update = update_trust_radius(
                trust_radius, actual_reduction, choice.predicted_reduction,
                scaled_step_norm, scaled_step_norm > 0.95 * trust_radius);
            trust_radius = update.radius;
            ratio = update.ratio;
            step_norm = (trial_x - x).norm();

            if (actual_reduction > 0.0) {
                accepted = true;
                accepted_x = std::move(trial_x);
                accepted_f = std::move(trial_f);
                accepted_cost = trial_cost;
            } else if (trust_radius <= EPSILON * (1.0 + x.norm())) {
                status = TrfStatus::STEP_TOLERANCE;
                terminated = true;
                break;
            }
        }

        if (terminated) { break; }
        if (!accepted) { break; }

        const TrfStatus tolerance_status = termination_status(
            actual_reduction, cost, step_norm, x.norm(), ratio, options);
        x = std::move(accepted_x);
        f = std::move(accepted_f);
        cost = accepted_cost;
        jacobian = evaluate_jacobian(x, f);
        gradient = jacobian.transpose() * f;
        ++result.iterations;

        if (options.callback) {
            const Scaling callback_scaling =
                coleman_li_scaling(x, gradient, bounds.lower, bounds.upper);
            const IterationInfo info{
                result.iterations,
                result.function_evaluations,
                result.jacobian_evaluations,
                cost,
                infinity_norm(gradient.cwiseProduct(callback_scaling.v)),
                trust_radius,
            };
            if (!options.callback(x, info)) {
                status = TrfStatus::USER_STOPPED;
                break;
            }
        }

        if (tolerance_status != TrfStatus::MAX_FUNCTION_EVALUATIONS) {
            status = tolerance_status;
            break;
        }
    }

    const Scaling final_scaling =
        coleman_li_scaling(x, gradient, bounds.lower, bounds.upper);
    result.x = std::move(x);
    result.residual = std::move(f);
    result.jacobian = std::move(jacobian);
    result.gradient = std::move(gradient);
    result.cost = cost;
    result.optimality =
        infinity_norm(result.gradient.cwiseProduct(final_scaling.v));
    result.active_mask =
        active_constraints(result.x, bounds.lower, bounds.upper, options.xtol);
    result.status = status;
    result.success = static_cast<int>(status) > 0;
    result.message = status_message(status);
    return result;
}

TrfResult trf_least_squares(const ResidualFunction& residual,
                            const Vector& initial_x,
                            const TrfOptions& options,
                            const JacobianFunction& jacobian) {
    return trf_least_squares(residual, initial_x,
                             Bounds::unbounded(initial_x.size()), options,
                             jacobian);
}

}  // namespace meow
