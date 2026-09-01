# meow

`meow` currently provides a C++20 trust-region reflective (TRF) solver for
bound-constrained nonlinear least squares:

\[
  \min_x \frac{1}{2}\lVert r(x)\rVert_2^2
  \quad\text{subject to}\quad \ell \le x \le u.
\]

The implementation is independent of SciPy. SciPy is used only by an optional
CTest comparison that checks the final solutions and costs on unconstrained,
active-bound, and box-constrained problems.

## Build and test

```bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

Requirements are CMake 3.20+, a C++20 compiler, and Eigen 3.3+. Python with
NumPy and SciPy enables the comparison test but is not needed by the library.
The build uses strict C++20 and treats warnings as errors by default. Set
`MEOW_WARNINGS_AS_ERRORS=OFF` only when integrating with an unusually noisy
toolchain.

The repository carries the same `.clang-format` and staged-file pre-commit
formatter as cuMES. Activate the versioned hook after cloning with:

```bash
git config core.hooksPath .githooks
```

## C++ API

```cpp
#include <meow/trf.hpp>

meow::Vector x0(2);
x0 << -2.0, 1.0;

auto residual = [](const meow::Vector& x) {
    return (meow::Vector(2) <<
        10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]).finished();
};

meow::TrfResult result = meow::trf_least_squares(residual, x0);
```

Pass `meow::Bounds` for box constraints and a `meow::JacobianFunction` as the
last argument when an analytic Jacobian is available. Otherwise, the solver
uses bound-aware forward differences. The result includes the final residual,
Jacobian, gradient, active-bound mask, termination reason, iteration count, and
all residual/Jacobian evaluation counts.

Algorithmically, the solver uses Coleman–Li interior scaling, an exact dense
trust-region quadratic solve, reflected steps at the first encountered bound,
and a constrained Cauchy-step fallback. It is intended for modest numbers of
optimization variables; residual evaluation may be arbitrarily expensive.

## Boozer symmetry-breaking diagnostic

The optimizer-side postprocessor compares maximum and RMS nonsymmetric
magnetic-field harmonics from two cuMES Boozer-v3 binary results. For QA,
the symmetric family is `n=0`, so all `n != 0` modes contribute:

```bash
python scripts/plot_boozer_symmetry_breaking.py \
  initial-boozer.bin final-boozer.bin \
  --output symmetry-breaking.png
```

The plotted amplitudes are normalized by the local `B_00(s)`. The calculation
accounts for the Boozer-v3 mixed grid through
`alpha_b = alpha + nfp*nu` and its toroidal integration Jacobian. Here
`alpha=nfp*zeta` is the field-period angle and `nu` is stored in physical
toroidal radians. `--helicity H` generalizes the retained symmetry family to
`n=H*m` for quasihelical diagnostics.
