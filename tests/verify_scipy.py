#!/usr/bin/env python3
"""Compare the C++ TRF solutions with scipy.optimize.least_squares."""

from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy.optimize import least_squares


def rosenbrock(x: np.ndarray) -> np.ndarray:
    return np.array([10.0 * (x[1] - x[0] ** 2), 1.0 - x[0]])


def rosenbrock_jacobian(x: np.ndarray) -> np.ndarray:
    return np.array([[-20.0 * x[0], 10.0], [-1.0, 0.0]])


def references() -> dict[str, tuple[np.ndarray, float]]:
    settings = dict(method="trf", ftol=1e-12, xtol=1e-12, gtol=1e-12,
                    max_nfev=5000)
    unconstrained = least_squares(
        rosenbrock, [-2.0, 1.0], jac=rosenbrock_jacobian, **settings)
    bounded = least_squares(
        rosenbrock, [2.0, 2.0], jac=rosenbrock_jacobian,
        bounds=([1.5, -np.inf], [np.inf, np.inf]), **settings)

    design = np.array([[1.0, 2.0], [2.0, -1.0], [-1.0, 1.0], [3.0, 1.0]])
    target = np.array([1.0, 2.0, -1.0, 4.0])
    linear = least_squares(
        lambda x: design @ x - target, [0.2, 0.2], jac=lambda _: design,
        bounds=([0.0, 0.0], [1.0, 1.0]), **settings)
    return {
        "rosenbrock": (unconstrained.x, unconstrained.cost),
        "bounded_rosenbrock": (bounded.x, bounded.cost),
        "bounded_linear": (linear.x, linear.cost),
    }


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_scipy.py TRF_REFERENCE_CASES_BINARY")
    output = subprocess.check_output([sys.argv[1]], text=True)
    expected = references()
    for line in output.splitlines():
        fields = line.split()
        name = fields[0]
        expected_x, expected_cost = expected.pop(name)
        values = np.asarray(fields[1:], dtype=float)
        actual_x = values[: expected_x.size]
        actual_cost = values[expected_x.size]
        np.testing.assert_allclose(actual_x, expected_x, rtol=2e-6, atol=2e-7)
        np.testing.assert_allclose(actual_cost, expected_cost, rtol=2e-6, atol=2e-10)
        print(f"{name}: C++ agrees with SciPy")
    if expected:
        raise AssertionError(f"missing C++ cases: {sorted(expected)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
