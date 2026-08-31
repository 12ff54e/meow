# Landreman-Paul QA/QH reproduction plan

## Goal

Reproduce the no-magnetic-well QA and QH optimization cases from Landreman &
Paul, *Magnetic fields with precise quasisymmetry*, using cuMES only as the
fixed-boundary equilibrium black box and `meow` as the least-squares optimizer.
The optimizer owns the boundary degrees of freedom, target definition,
finite-difference policy, and stopping policy. cuMES owns only equilibrium
solution and publication of the quantities required by the target.

The archived final configurations are:

- QA: `configurations/new_QA/input.20210704-01-050_..._nfp2_QA`, target aspect
  ratio 6, two field periods.
- QH: `configurations/new_QH/input.20210704-01-063_..._nfp4_QH_A8_...`, target
  aspect ratio 8, four field periods.

The paths above are relative to the supplemental data directory
`../20211102-01-precise_quasisymmetry_zenodo`.

## Reproduction contract

The target must follow the archived `QuasisymmetryRatioError`, rather than a
nearby discretization:

- `s = 0.0, 0.1, ..., 1.0`, interpreted as normalized toroidal flux.
- Linear interpolation from VMEC/cumes half-grid quantities to requested
  surfaces, with linear extrapolation to the magnetic axis and edge.
- A target-only angular grid of 63 poloidal by 64 toroidal points over one
  field period, independent of the equilibrium grid.
- The normalized toroidal flux derivative
  `d psi / d s = -Phi_edge / (2*pi)` in the archived VMEC sign convention.
- QA helicity `(M,N)=(1,0)`. QH is specified to SIMSOPT as `(m,n)=(1,-1)`;
  SIMSOPT multiplies `n` by `nfp`, so its physical toroidal helicity is `N=-4`.
- Surface weights ramp linearly from 1 to 30 for the final QA refinement and
  from 1 to 2 for the final QH refinement.
- The residual is volume-normalized using the signed VMEC Jacobian, exactly as
  in the archived implementation. cuMES orientation differences must be
  handled once at the coordinate bridge, not hidden in optimizer policy.

The paper's final refinement objective contains aspect ratio and QS residuals.
The earlier QA construction stage additionally targets mean iota 0.42. Hence
an exact workflow reproduction changes objective composition between stages;
the reusable optimizer API will continue to support the user's QA composite
with an iota residual as an explicit option.

The user's `deps/BSplineInterpolation` library is not used for this radial
step. The archived implementation explicitly selects linear interpolation and
linear endpoint extrapolation, while a B-spline would define a different
least-squares problem.

## Why the supplement contains `booz_xform`

The optimization drivers import `QuasisymmetryRatioError` and obtain the QS
residual directly from VMEC-coordinate quantities. Their Boozer imports are
commented out. `booz_xform` is run afterward by the convergence/analysis jobs
to produce an independent Fourier-spectrum diagnostic, make figures, and
support downstream neoclassical analysis. It is therefore part of validating
and presenting the configurations, not part of the optimization target or the
cuMES solver API.

## Progressive implementation

1. **Freeze the reference.** Add a read-only evaluator for the supplemental
   `wout` files and record the archived QA/QH objective components. This first
   verifies target conventions without involving a different equilibrium
   solver.
2. **Accept the vacuum-current inputs.** Permit `ncurr=1`, `curtor=0`, and an
   empty/zero current profile. Keep rejecting a nonzero prescribed total
   current with a zero-normalization profile. Add host and device tests.
3. **Import the boundaries.** Add a narrow VMEC-namelist-to-cuMES JSON utility
   for the fields used by these two inputs, with strict diagnostics for
   unsupported values. Check in generated QA/QH inputs so runs are
   reproducible without a Python/Fortran parser dependency.
4. **Match the target discretization.** Add arbitrary normalized-toroidal-flux
   surfaces, exact linear half-grid interpolation/extrapolation, and 63x64
   Fourier resampling to the optimizer-side QS target. Keep the native-grid
   helper available as a cheaper diagnostic, but make the paper convention
   explicit and tested.
5. **Cross the equilibrium boundary.** Solve the imported final QA and QH
   boundaries with cuMES and compare aspect ratio, iota, and QS profile against
   both the supplied `wout` and the frozen reference metric. VMEC and cuMES have
   different discrete convergence oracles, so comparison tolerances will be
   stated rather than assuming byte-identical equilibria.
6. **Reproduce optimization.** Generalize the meow example from one harmonic
   to the paper's staged stellarator-symmetric `RBC/ZBS` parameter vector:
   modes through 4, then through 5, with `RBC(0,0)` fixed. Use forward finite
   differences with relative step `1e-5` and absolute step `1e-9`, preserving
   failed-solve information as optimizer residual failures.
7. **Validate, then optionally postprocess.** Compare final objective
   components and boundary coefficients. Boozer transformation remains an
   optional, independent validation step and is not linked into either cuMES
   or the target evaluator.

Each numbered implementation slice is committed only after its focused tests
pass. Full GPU regression is run before declaring the reproduction complete.

The implementation executable is `cumes_landreman_optimize`. Its optimizer
variables are absolute Fourier coefficients, matching SIMSOPT's relative-step
definition. It persistently writes strict cuMES input JSON after accepted
iterations, so a long run can be inspected or restarted without translating a
diagnostic normalization document.

Final-grid hot restarts were tested as a finite-difference acceleration and
rejected: at the adapted `1e-12` equilibrium gate, a `1e-9` boundary change can
follow a measurably different near-degenerate equilibrium path, contaminating
the numerical Jacobian. The reproduction therefore uses the same cold
multigrid path for every sample. Mode-range arguments allow the expensive
stages to be checkpointed separately without changing this derivative policy.

## Acceptance evidence

- Reference evaluation reproduces the supplemental SIMSOPT totals to roundoff.
- Imported inputs retain all active boundary coefficients and solver controls.
- Repeated evaluation of an unchanged boundary is deterministic.
- The 11-surface residual contains `11 * 63 * 64` scalar QS residuals.
- QA/QH runs report equilibrium convergence separately from target quality.
- The result report contains aspect, mean iota when requested, weighted and
  unweighted QS totals, and per-surface QS values.

## Frozen archived values

Running `scripts/evaluate_landreman_qs.py` on the terminal optimization wout
files reproduces the original SIMSOPT log values (differences below the last
printed digit):

| case | composite objective | weighted QS | unweighted QS | axis QS | edge QS |
| --- | ---: | ---: | ---: | ---: | ---: |
| QA | 1.31838917076796e-6 | 1.318328430729835e-6 | 9.803160636071010e-8 | 1.692667660112380e-8 | 1.395084878766294e-8 |
| QH | 3.447641837023028e-5 | 3.447629144415016e-5 | 2.110744659368261e-5 | 9.556217390373448e-7 | 5.730803009640856e-6 |

The terminal files are `wout_nfp2_QA_000_000063.nc` in the final QA
optimization directory and `wout_nfp4_QH_000_000090.nc` in the final QH
optimization directory. The `configurations/new_QA` and `configurations/new_QH`
wout files are later high-resolution convergence re-solves. They intentionally
produce different target values, so both baselines must remain distinguishable.

The checked-in `examples/landreman/qa.json` and `qh.json` are generated with
`vmec_namelist_to_cumes_json.py --minimum-ftol 1e-12 --wout-axis <wout>`. The
explicit clamp is required because the archived VMEC convergence study
requests `1e-17` and `2e-17`. Although cuMES accepts tolerances down to
`1e-16`, the imported QA coarse-grid trajectory stalls near `8e-13`; `1e-12`
is the qualified 3-D convergence scale and remains far below the QS target.
The explicit axis import replaces the namelists' zero axis placeholder, for
which VMEC applies an internal automatic initialization that cuMES does not
duplicate. Boundary, flux, current, resolution, and stage-size data are
otherwise retained.

The separate `qa_start.json` and `qh_start.json` files come from iteration zero
of the archived final-refinement directories, not from the terminal or later
high-resolution configurations. They contain modes through 4 and are the
correct inputs to the two-stage `max_mode=4`, then `max_mode=5` reproduction.
Their VMEC caps are raised explicitly with `--minimum-niter` because cuMES's
coarse-grid descent has not reached the adapted `1e-12` gate after the
archived 600 iterations. This changes only the solver work allowance, not the
boundary, target, or convergence threshold.

The QH run explicitly selects the retained Catmull-Rom radial transfer. The
default global B-spline transfer converges through `ns=100` but overshoots near
the LCFS during the `100 -> 150` transition, producing an invalid Jacobian.
Catmull-Rom completes both remaining stages. This is solver-path provenance,
not part of the target definition.

## cuMES final-boundary cross-check

`cumes_landreman_evaluate` solves the checked-in final boundaries and applies
the same 11-surface, 63x64 target in meow. With the solver adaptations above:

| case | composite objective | weighted QS | unweighted QS | aspect | mean iota |
| --- | ---: | ---: | ---: | ---: | ---: |
| QA | 1.780812493182744e-6 | 1.780751753144979e-6 | 1.048284137018133e-7 | 6.000007793589530 | 0.419325650931308 |
| QH | 4.727458337966436e-5 | 4.727445645358373e-5 | 2.784061841868241e-5 | 8.000011266147551 | -1.243256718499347 |

QA ends with `(fsqr,fsqz,fsql) = (9.931e-13,5.281e-13,9.056e-14)`.
QH ends with `(9.966e-13,7.400e-13,1.563e-13)`. These are successful
cuMES equilibria according to the adapted `1e-12` gate, but their QS values
do not equal the archived VMEC terminal metrics. The difference is therefore
reported as a cross-solver/discretization discrepancy, concentrated toward
the heavily weighted edge, rather than hidden by target renormalization.
