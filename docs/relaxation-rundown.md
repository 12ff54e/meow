# Relaxation rundown schema v1

`cumes_landreman_optimize` executes an ordered magnetic-equilibrium relaxation
from one strict JSON document. The document is optimizer policy: cuMES still
receives an ordinary equilibrium input and remains responsible only for solving
that equilibrium.

The schema identifier is `meow-relaxation-rundown-v1`. Unknown keys, missing
required keys, invalid types, and inconsistent values are errors. Use
`--dry-run` to validate a file and print its resolved input and ordered steps
without starting CUDA work.

## Top-level structure

```json
{
  "schema": "meow-relaxation-rundown-v1",
  "name": "landreman-qa-construction",
  "equilibrium_input": "qa_analytic.json",
  "case": "qa",
  "workflow": "construction",
  "target": {},
  "initialization": {},
  "optimizer": {},
  "output": {},
  "rundown": []
}
```

`equilibrium_input` is resolved relative to the rundown file. Output paths are
interpreted relative to the process working directory, making scheduler output
placement explicit and easy to override.

`case` is `qa` or `qh`. It chooses whether the composite target uses the QA or
QH target evaluator. `workflow` is `construction` or `refinement`; it controls
checkpoint naming only. The numerical policy is always taken from the
remaining JSON fields.

## Target

```json
"target": {
  "aspect_ratio": 6.0,
  "mean_iota": 0.42,
  "helicity_m": 1,
  "helicity_n_per_field_period": 0,
  "ntheta": 63,
  "nzeta": 64,
  "surfaces": [0.0, 0.1, 1.0],
  "surface_weights": [1.0, 1.0, 1.0]
}
```

`mean_iota` may be `null` to omit that residual. It is supported only for QA.
QA requires zero toroidal helicity; QH requires nonzero toroidal helicity. The
runner multiplies `helicity_n_per_field_period` by the equilibrium `nfp` before
constructing the target. Surfaces must be strictly increasing in `[0, 1]`, and
weights must be finite, nonnegative, and have the same extent.

## Initialization and continuation

```json
"initialization": {
  "axisymmetric_seed_amplitude": 1e-4,
  "track_axis_predictor": true,
  "radial_transfer": "default"
}
```

The seed may be `null`. When present, it is applied only at the first selected
step and only if every active 3-D boundary mode is zero. Axis tracking updates
the accepted magnetic-axis predictor and applies the boundary-centerline delta
to trials. Radial transfer is `default` or `catmull-rom`.

## Optimizer and Jacobian

```json
"optimizer": {
  "ftol": 1e-8,
  "xtol": 1e-8,
  "gtol": 1e-8,
  "max_function_evaluations": 0,
  "max_accepted_iterations": 0,
  "jacobian": {
    "method": "parallel-aggressive-broyden",
    "workers": 4,
    "relative_step": 0.0031622776601683794,
    "absolute_step": 1e-7,
    "refresh_interval": 8,
    "minimum_reduction_ratio": 0.05,
    "maximum_secant_error": 0.5,
    "scale_variables": false,
    "warm_start_trials": false
  }
}
```

Zero evaluation/iteration limits retain the optimizer default. The Jacobian
method chooses the existing evaluation backend. Refresh interval and secant
safeguards independently control Broyden reuse, so the complete numerical
policy is visible without relying on a method-name preset. `workers` applies
to parallel Jacobian methods.

## Ordered rundown and phases

```json
"rundown": [
  {
    "name": "mode-1",
    "max_boundary_mode": 1,
    "equilibrium": {
      "mpol": 3,
      "ntor": 3,
      "minimum_iterations": 6000,
      "tolerance_floor": 1e-12
    },
    "phases": [
      {
        "name": "qualified",
        "equilibrium_tolerance_floor": null,
        "cold_start": false,
        "stopping": {
          "minimum_iterations": 0,
          "progress_window": 0,
          "minimum_relative_progress": 0.0
        }
      }
    ]
  }
]
```

Step names are unique and `max_boundary_mode` must increase strictly. Setting
both `mpol` and `ntor` to zero retains the input resolution; otherwise both
must be positive. Minimum iterations and the tolerance floor can only make the
equilibrium solve stricter or give it more work.

Every step has one or more ordered phases. A null phase tolerance restores the
step's qualified equilibrium tolerances. A numeric value relaxes them by
taking the larger tolerance. `cold_start` prevents the previous phase's
equilibrium snapshot from being used as the next phase's initial state.

Phase stagnation is disabled when `progress_window` and
`minimum_relative_progress` are zero. Otherwise, after `minimum_iterations`,
the phase stops when relative objective improvement across the requested
window is smaller than the configured threshold.

## CLI

```bash
build-cumes/cumes_landreman_optimize \
  examples/landreman/qa-construction.rundown.json

build-cumes/cumes_landreman_optimize --dry-run \
  --first-step mode-2 --last-step mode-4 \
  examples/landreman/qa-construction.rundown.json
```

Run `cumes_landreman_optimize --help` for output, step-range, work-limit,
Jacobian-method, worker-count, and iteration-directory overrides. Overrides do
not mutate the rundown file or the equilibrium input.

The checked-in canonical files are:

- `qa-construction.rundown.json`
- `qh-construction.rundown.json`
- `qa-refinement.rundown.json`
- `qh-refinement.rundown.json`
