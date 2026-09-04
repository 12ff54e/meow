# cuMES FetchContent integration

## Goal

Allow a fresh meow checkout to obtain and build the cuMES solver through
CMake, while preserving the existing installed-package workflow and keeping a
normal meow-only build independent of CUDA.

## Plan

1. Add an opt-in `MEOW_FETCH_CUMES` CMake option. It has an effect only when
   `MEOW_BUILD_CUMES_INTEGRATION=ON`; otherwise CMake must not download or
   configure cuMES.
2. When fetching is enabled, use `FetchContent` with the public cuMES GitHub
   repository pinned to the audited `interface` revision
   `f61c959334bb62e14c049c66335580b45f63610d`. Expose cache variables for an
   intentional repository or revision override, including local/offline
   development through CMake's standard `FETCHCONTENT_SOURCE_DIR_CUMES`
   override.
3. Configure the embedded dependency as a library: disable its CLI, tests,
   benchmarks, optional output backends, Boozer converter, vacuum-field
   coupling, and verification dumps. Retain its qualified B-spline multigrid
   transfer and fetch only that submodule, whose relative URL resolves over
   HTTPS. Do not fetch the SSH-addressed vacuum-field and magnetic-coordinate
   submodules.
4. Keep `find_package(cuMES CONFIG REQUIRED)` unchanged when fetching is off,
   and accept a pre-existing `cumes::solver` target in either mode so parent
   projects can provide the dependency themselves.
5. Document fresh-clone and installed-package commands, CUDA/toolchain
   requirements, the deliberately disabled optional features, and override
   behavior.
6. Validate three configurations: meow without cuMES (and therefore without a
   fetch), the installed cuMES package path, and FetchContent using both the
   pinned public repository and a local source override. Build and run the
   applicable tests for both integration paths.

## Acceptance criteria

- The default configure performs no network access and still needs only Eigen.
- `MEOW_BUILD_CUMES_INTEGRATION=ON` continues to consume an installed cuMES
  package unless `MEOW_FETCH_CUMES=ON` is explicitly selected.
- The FetchContent build defines `cumes::solver`, builds meow's cuMES-backed
  applications, and does not create cuMES CLI, test, benchmark, or optional
  post-processing targets.
- Reconfiguring as a subproject does not forcibly replace cache choices that a
  parent project supplied before adding meow.

## Implementation and validation

`MEOW_FETCH_CUMES=ON` declares the pinned source with
`OVERRIDE_FIND_PACKAGE`; the integration continues to consume the dependency
through its single `find_package(cuMES CONFIG REQUIRED)` path. The override is
a CMake 3.24 feature, so only fetch mode raises the effective minimum from the
project-wide CMake 3.20 baseline. A parent-provided `cumes::solver` target still
takes precedence.

The public-fetch configuration resolved the exact pinned revision and only
initialized `deps/BSplineInterpolation`; the vacuum-field and
magnetic-coordinate submodules remained uninitialized. Its complete default
build succeeded and all 12 meow tests passed. The same 12 tests passed against
the separately installed cuMES package. A local
`FETCHCONTENT_SOURCE_DIR_CUMES` override and the CUDA-free base configuration
also configured successfully.

cuMES currently declares both its double and float CUDA implementation
libraries as normal build targets. Consequently, a default all-target
FetchContent build compiles both variants even though meow links only the
selected double solver. Avoiding that extra compilation requires an upstream
cuMES target/option change; meow does not patch dependency sources.
