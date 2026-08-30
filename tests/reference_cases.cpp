#include "meow/trf.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

void print_result(std::string_view name, const meow::TrfResult& result) {
    std::cout << std::setprecision(17) << name;
    for (Eigen::Index i = 0; i < result.x.size(); ++i) {
        std::cout << ' ' << result.x[i];
    }
    std::cout << ' ' << result.cost << ' ' << result.optimality << '\n';
}

}  // namespace

int main() {
    const auto rosenbrock = [](const meow::Vector& x) {
        return (meow::Vector(2) << 10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0])
            .finished();
    };
    const auto rosenbrock_jacobian = [](const meow::Vector& x) {
        meow::Matrix value(2, 2);
        value << -20.0 * x[0], 10.0, -1.0, 0.0;
        return value;
    };

    meow::TrfOptions options;
    options.ftol = 1e-12;
    options.xtol = 1e-12;
    options.gtol = 1e-12;
    options.max_function_evaluations = 5000;

    meow::Vector x0(2);
    x0 << -2.0, 1.0;
    print_result("rosenbrock", meow::trf_least_squares(rosenbrock, x0, options,
                                                       rosenbrock_jacobian));

    x0 << 2.0, 2.0;
    const meow::Bounds bounds{
        (meow::Vector(2) << 1.5, -std::numeric_limits<double>::infinity())
            .finished(),
        meow::Vector::Constant(2, std::numeric_limits<double>::infinity())};
    print_result("bounded_rosenbrock",
                 meow::trf_least_squares(rosenbrock, x0, bounds, options,
                                         rosenbrock_jacobian));

    const meow::Matrix design =
        (meow::Matrix(4, 2) << 1.0, 2.0, 2.0, -1.0, -1.0, 1.0, 3.0, 1.0)
            .finished();
    const meow::Vector target =
        (meow::Vector(4) << 1.0, 2.0, -1.0, 4.0).finished();
    const auto linear = [design, target](const meow::Vector& x) {
        return design * x - target;
    };
    const auto linear_jacobian = [design](const meow::Vector&) {
        return design;
    };
    x0 << 0.2, 0.2;
    const meow::Bounds box{meow::Vector::Zero(2), meow::Vector::Ones(2)};
    print_result(
        "bounded_linear",
        meow::trf_least_squares(linear, x0, box, options, linear_jacobian));

    constexpr double times[] = {0.0, 0.25, 0.5, 1.0, 1.5, 2.0};
    meow::Vector observations(6);
    for (Eigen::Index i = 0; i < observations.size(); ++i) {
        observations[i] = 2.5 * std::exp(-0.7 * times[i]) + 0.15;
    }
    const auto exponential = [observations, &times](const meow::Vector& x) {
        meow::Vector value(observations.size());
        for (Eigen::Index i = 0; i < value.size(); ++i) {
            value[i] =
                x[0] * std::exp(-x[1] * times[i]) + x[2] - observations[i];
        }
        return value;
    };
    meow::Vector exponential_x0(3);
    exponential_x0 << 1.0, 0.2, 0.0;
    const meow::Bounds positive_box{meow::Vector::Zero(3),
                                    meow::Vector::Constant(3, 5.0)};
    print_result("exponential",
                 meow::trf_least_squares(exponential, exponential_x0,
                                         positive_box, options));
    return 0;
}
