#include "meow/trf.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) { throw std::runtime_error(message); }
}

void require_near(const meow::Vector& actual,
                  const meow::Vector& expected,
                  double tolerance,
                  const std::string& message) {
    if ((actual - expected).cwiseAbs().maxCoeff() > tolerance) {
        throw std::runtime_error(message + ": got " +
                                 std::to_string(actual.transpose()[0]));
    }
}

void test_unconstrained_rosenbrock() {
    const auto residual = [](const meow::Vector& x) {
        meow::Vector value(2);
        value << 10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0];
        return value;
    };
    const auto jacobian = [](const meow::Vector& x) {
        meow::Matrix value(2, 2);
        value << -20.0 * x[0], 10.0, -1.0, 0.0;
        return value;
    };
    meow::Vector x0(2);
    x0 << -2.0, 1.0;

    meow::TrfOptions options;
    options.gtol = 1e-12;
    options.ftol = 1e-12;
    options.xtol = 1e-12;
    const meow::TrfResult result =
        meow::trf_least_squares(residual, x0, options, jacobian);

    require(result.success, "unconstrained Rosenbrock did not converge");
    meow::Vector expected(2);
    expected << 1.0, 1.0;
    require_near(result.x, expected, 2e-7, "wrong Rosenbrock solution");
    require(result.cost < 1e-20, "Rosenbrock cost is too large");
}

void test_active_lower_bound() {
    const auto residual = [](const meow::Vector& x) {
        meow::Vector value(2);
        value << 10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0];
        return value;
    };
    const auto jacobian = [](const meow::Vector& x) {
        meow::Matrix value(2, 2);
        value << -20.0 * x[0], 10.0, -1.0, 0.0;
        return value;
    };

    meow::Vector x0(2);
    x0 << 2.0, 2.0;
    meow::Bounds bounds{
        (meow::Vector(2) << 1.5, -std::numeric_limits<double>::infinity())
            .finished(),
        meow::Vector::Constant(2, std::numeric_limits<double>::infinity())};
    meow::TrfOptions options;
    options.gtol = 1e-11;
    options.ftol = 1e-12;
    options.xtol = 1e-12;

    const meow::TrfResult result =
        meow::trf_least_squares(residual, x0, bounds, options, jacobian);
    meow::Vector expected(2);
    expected << 1.5, 2.25;
    require(result.success, "bounded Rosenbrock did not converge");
    require_near(result.x, expected, 2e-6, "wrong bounded solution");
    require(result.active_mask[0] == -1, "lower bound was not marked active");
}

void test_numerical_jacobian_curve_fit() {
    constexpr double times[] = {0.0, 0.25, 0.5, 1.0, 1.5, 2.0};
    meow::Vector observations(6);
    for (Eigen::Index i = 0; i < observations.size(); ++i) {
        observations[i] = 2.5 * std::exp(-0.7 * times[i]) + 0.15;
    }
    const auto residual = [observations, &times](const meow::Vector& x) {
        meow::Vector value(observations.size());
        for (Eigen::Index i = 0; i < value.size(); ++i) {
            value[i] =
                x[0] * std::exp(-x[1] * times[i]) + x[2] - observations[i];
        }
        return value;
    };

    meow::Vector x0(3);
    x0 << 1.0, 0.2, 0.0;
    meow::Bounds bounds{meow::Vector::Zero(3), meow::Vector::Constant(3, 5.0)};
    meow::TrfOptions options;
    options.gtol = 1e-10;
    options.max_function_evaluations = 1000;
    const meow::TrfResult result =
        meow::trf_least_squares(residual, x0, bounds, options);

    meow::Vector expected(3);
    expected << 2.5, 0.7, 0.15;
    require(result.success, "finite-difference curve fit did not converge");
    require_near(result.x, expected, 2e-5, "wrong curve-fit solution");
    require(result.function_evaluations > result.jacobian_evaluations,
            "finite-difference calls were not counted");
}

void test_callback_stop() {
    const auto residual = [](const meow::Vector& x) {
        return (meow::Vector(1) << x[0] - 10.0).finished();
    };
    meow::Vector x0(1);
    x0 << 0.0;
    meow::TrfOptions options;
    options.callback = [](const meow::Vector&, const meow::IterationInfo&) {
        return false;
    };
    const meow::TrfResult result =
        meow::trf_least_squares(residual, x0, options);
    require(result.status == meow::TrfStatus::USER_STOPPED,
            "callback did not stop the solve");
    require(!result.success, "a callback stop should not report convergence");
}

void test_invalid_bounds() {
    const auto residual = [](const meow::Vector& x) { return x; };
    meow::Vector x0(1);
    x0 << 2.0;
    meow::Bounds bounds{meow::Vector::Zero(1), meow::Vector::Ones(1)};
    bool threw = false;
    try {
        static_cast<void>(meow::trf_least_squares(residual, x0, bounds));
    } catch (const std::invalid_argument&) { threw = true; }
    require(threw, "infeasible initial point was accepted");
}

}  // namespace

int main() {
    try {
        test_unconstrained_rosenbrock();
        test_active_lower_bound();
        test_numerical_jacobian_curve_fit();
        test_callback_stop();
        test_invalid_bounds();
    } catch (const std::exception& error) {
        std::cerr << "test_trf: " << error.what() << '\n';
        return 1;
    }
    std::cout << "All TRF unit tests passed.\n";
    return 0;
}
