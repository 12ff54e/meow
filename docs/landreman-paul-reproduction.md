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
   or the target evaluator. The optimizer-side
   `scripts/plot_boozer_symmetry_breaking.py` utility reproduces the
   supplement's maximum/RMS `B_mn`-versus-radius diagnostic from two
   Boozer-v3 results. The diagnostic applies the field-period relation
   `alpha_b=alpha+nfp*nu`; the optimization target remains independent of this
   post-processing transform.

Each numbered implementation slice is committed only after its focused tests
pass. Full GPU regression is run before declaring the reproduction complete.

## Analytic-start continuation extension

The first implementation covered the archived final mode-4/mode-5
refinements. The complete construction must also expose the earlier analytic
starting boundaries instead of referring to those refinement inputs as the
start of the whole calculation. The extension is split into these reviewable
steps:

1. Import the two sparse analytic boundaries from runs 021 (QA) and 039 (QH),
   using their iteration-zero `wout` files only to supply the magnetic-axis
   predictor that VMEC otherwise initializes internally. Verify that the
   imported boundary contains exactly the documented nonzero harmonics.
2. Add an explicit construction workflow to `cumes_landreman_optimize` while
   retaining the existing refinement workflow. Construction uses uniform
   radial QS weights, includes the QA mean-iota target 0.42, and follows the
   archived mode schedules 1--4 for QA and 1--5 for QH. Refinement continues
   to use modes 4--5, no iota residual, and the final radial weight ramps.
3. Apply the archived construction-stage equilibrium resolutions
   `mpol=ntor={3,5,6,6}` for QA and `{3,5,6,6,6}` for QH. Keep every accepted
   boundary and equilibrium distinguishable by workflow and mode so the
   repeated mode-4 stage cannot overwrite construction output.
4. Add focused tests for the workflow schedules and analytic harmonics, build
   with the cuMES integration enabled, and solve each imported initial state
   once before advertising the commands.

This extension stays entirely in meow. The construction/refinement target
policy and continuation schedule are optimizer concerns; cuMES remains an
unchanged equilibrium solve API.

## QA construction discrepancy investigation

The first cuMES-backed analytic QA run converged to objective `0.044202869362`
with mean iota `0.284125511417`, far from run 021. The archived logs isolate a
policy mismatch before any equilibrium-model comparison is needed:

- Run 021 uses forward differences with relative step
  `3.1622776601683794e-3` and absolute floor `1e-7`.
- The first meow construction run incorrectly reused the final-refinement
  values `1e-5` and `1e-9`.
- At the axisymmetric seed, meow's `1e-9` perturbations changed mean iota only
  at roughly `1e-17`. The numerical Jacobian therefore represented the iota
  target as locally flat and initially optimized almost only aspect ratio.
- The archived mode-1 solve instead reduces the least-squares cost from
  `0.5882` to `0.0048204` in 37 objective evaluations. The incorrect meow run
  stopped at cost `0.0228301` and mean iota about `0.284`.

The corrective work is deliberately staged:

1. Make finite-difference policy part of the Landreman workflow description:
   QA construction uses run 021's `3.1622776601683794e-3` / `1e-7`, QH
   construction uses run 039's `1e-3` / `1e-7`, and refinement retains
   `1e-5` / `1e-9`.
2. Re-run QA mode 1 from the analytic seed and compare its objective, iota,
   accepted step sizes, and boundary coefficients with archived iteration 63
   before permitting modes 2--4.
3. Continue each higher mode only from the corrected checkpoint. Preserve the
   first run under `opt-qa-analytic` as negative evidence rather than
   overwriting it.
4. If the corrected trajectory still diverges materially, compare the
   numerical Jacobian column-by-column against run 021-scale perturbations and
   separate optimizer differences from cuMES/VMEC equilibrium differences.

The corrected-step mode-1 diagnostic did reproduce the archived first accepted
cost (`0.149506` versus approximately `0.14952`), but later stalled at
objective `0.0526046` with 20 accepted iterations.  Its rejected trial solves
exposed a second initialization mismatch: the optimizer was retaining the
iteration-zero numerical magnetic-axis predictor while the optimized
`m=0, n=1` boundary centerline moved.  The original analytic VMEC input instead
sets `RAXIS_CC=0` and `ZAXIS_CS=0`, requesting automatic axis initialization.
For construction only, meow must therefore rebuild the predictor before every
equilibrium solve from the current boundary centerline,
`raxis_c[n]=rbc(0,n)` and `zaxis_s[n]=zbs(0,n)`.  This changes only the cuMES
cold-start guess, not the equilibrium degrees of freedom or target.  The
refinement workflow retains its imported converged-axis predictor.

The next qualification is deliberately limited to QA mode 1.  It must first
reduce the failed-equilibrium barrier trials and track the archived boundary
and cost trajectory more closely; only then are construction modes 2--4 run.

That qualification passed.  With the predictor refreshed, mode 1 terminated
normally after 23 accepted iterations and 236 objective evaluations at
objective `0.0096228845593`, compared with the archived SIMSOPT objective
`0.0096408`.  No equilibrium-barrier trial occurred.  Representative final
coefficients also agree closely:

| coefficient | meow + cuMES | archived run 021 |
| --- | ---: | ---: |
| `rbc(0,1)` | 0.1031565673 | 0.1029584216 |
| `rbc(1,-1)` | 0.0212384768 | 0.0209931493 |
| `rbc(1,0)` | 0.1650055435 | 0.1656929553 |
| `zbs(0,1)` | -0.1323894290 | -0.1326896627 |
| `zbs(1,0)` | 0.2237183921 | 0.2224868566 |

This recovery identifies stale cold-start axis initialization as the cause of
the later trajectory failure.  The smaller remaining differences are
consistent with the independent equilibrium solver and TRF implementations;
they no longer change the mode-1 optimum materially.  Modes 2--4 can therefore
continue from the refreshed-axis checkpoint.

The first higher-mode continuation showed that the boundary centerline alone
is not a sufficiently accurate predictor indefinitely.  Mode 2 stopped at
objective `3.43599596128e-4`, and mode 3 at `1.50742626669e-4`; by then some
finite-difference solves required roughly 6,700 iterations or failed.  The
centerline refresh had removed the original stale-axis error, but discarded
the offset between the centerline and the actual axis after every accepted
solve.

The follow-up policy keeps the mode-1 correction and makes it continuation
aware:

1. Seed the first analytic construction stage from its boundary centerline.
2. After an accepted optimizer step, extract the actual `m=0` axis from the
   converged cuMES equilibrium and store it in the accepted problem.
3. Evaluate every trial from that fixed accepted-axis predictor plus only the
   trial's change in `rbc(0,n)` / `zbs(0,n)`.  Thus every column of one
   finite-difference Jacobian uses the same accepted equilibrium reference and
   remains deterministic.
4. Carry the accepted axis through mode changes, padding only newly introduced
   toroidal modes.  A resumed higher-mode construction likewise keeps the axis
   serialized in its checkpoint.

This is a predictor update only; the converged equilibrium and objective are
still functions of the trial boundary, and no solver state is hot-restarted.

The resulting fresh QA construction completed all four stages within the
100-accepted-step limit:

| max boundary mode | accepted iterations | objective | archived objective |
| ---: | ---: | ---: | ---: |
| 1 | 23 | `9.62215079953e-3` | `9.6408e-3` |
| 2 | 18 | `1.40881439229e-4` | `1.5611e-4` |
| 3 | 13 | `3.74704919717e-6` | `3.5144e-6` |
| 4 | 7 | `4.57544450735e-7` | `2.92794425995e-7` |

There was one rejected oversized mode-2 trust-region proposal and no repeated
failed finite-difference columns.  Finite-difference equilibria remained near
2,100--2,300 solver iterations through modes 3 and 4, instead of climbing to
6,700 and stalling.  The final objective is within a factor 1.56 of run 021,
which is the expected level for this independent equilibrium and least-squares
implementation and is a material reproduction of the archived QA result.  The
65 accepted/initial input-equilibrium pairs, four stage checkpoints, and log
are stored in `../opt-qa-analytic-accepted-axis`.

The corresponding fresh QH construction also completed all five stages within
the 100-accepted-step limit:

| max boundary mode | accepted iterations | objective | archived objective |
| ---: | ---: | ---: | ---: |
| 1 | 21 | `1.40503531818e-1` | approximately `1.39906e-1` |
| 2 | 32 | `1.91355839667e-3` | approximately `2.8614e-3` |
| 3 | 15 | `1.20166888048e-4` | approximately `1.79572e-4` |
| 4 | 8 | `7.30668510539e-5` | approximately `4.1688e-5` |
| 5 | 12 | `5.36939265810e-5` | `3.07486398207e-5` |

Each stage stopped on the optimizer's step tolerance, after 88 accepted steps
in total. There were two isolated rejected oversized proposals, in modes 1
and 2, and no equilibrium failures in modes 3--5. The final aspect ratio is
`8.00001319737` and the final mean rotational transform is
`-1.22245005823`. The mode-1 objective agrees with archived run 039 within
0.5%; the paths then differ across the independent equilibrium and optimizer
implementations, with the final objective 1.746 times the archived construction
value. The 93 accepted/initial input-equilibrium pairs, five stage checkpoints,
latest checkpoint, and complete log are stored in
`../opt-qh-analytic-accepted-axis`.

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

The still-earlier `qa_analytic.json` and `qh_analytic.json` inputs contain the
sparse boundaries from construction runs 021 and 039. Select
`qa-construction` or `qh-construction` to activate their uniform QS weights,
the archived mode-1-up continuation resolution, and (for QA only) the 0.42
mean-iota residual. The archived numerical-Jacobian inputs give initial cuMES
objectives 1.1764 for QA and 10.5719707933 for QH; the QA value is exactly the
sum of its unit aspect residual and squared 0.42 iota residual to the reported
precision. Both sparse boundaries converge through their imported
`ns=12,25,50` equilibrium sequence at the adapted `1e-12` gate.
Construction stages at `mpol=ntor=6` receive at least a 10,000-iteration work
allowance: the first QA mode-3 transition ended its initial 6,000-iteration
attempt at residuals `(1.72e-12,8.44e-13,1.64e-12)`. QA mode 4 receives
30,000 because several valid finite-difference samples on the first, incorrect
trajectory were still near `7e-9` at 10,000 iterations. The corrected
construction retains the common `1e-12` equilibrium gate; only the work cap is
carried forward.

The QH run explicitly selects the retained Catmull-Rom radial transfer. The
default global B-spline transfer converges through `ns=100` but overshoots near
the LCFS during the `100 -> 150` transition, producing an invalid Jacobian.
Catmull-Rom completes both remaining stages. This is solver-path provenance,
not part of the target definition.

## Analytic equilibrium-tangent qualification

The current optimizer evaluates the nonlinear equilibrium once per trial and
assembles the dense target Jacobian from cuMES's retained final-grid
linearization. cuMES solves

```text
F_u du = -F_x dx
```

for each optimizer-owned boundary direction; meow then applies the analytic
QS, aspect-ratio, and mean-iota target JVPs. cuMES remains unaware of the
target definition. At the mode-1 QH start, a strict `1e-6` tangent solve for
`rbc(0,1)` agrees with a centered, fully re-solved nonlinear oracle to 5.2% in
the complete residual vector and 0.14% in the invariant objective derivative.
The production optimizer uses cuMES's `1e-4` linear tolerance, which is more
appropriate for trust-region steps and later 120-column stages.

The sparse QA start needs one additional policy. It is exactly axisymmetric,
so transform and QA error are stationary to first order in the 3-D modes. Run
021 escaped this manifold through its one-sided finite-difference Jacobian.
The analytic workflow instead inserts a deterministic chiral `1e-4` seed only
when the active 3-D coefficients are all exactly zero. With this seed, a fresh
mode-1 run reaches objective `0.0101326414112` and mean iota
`0.386319491437`; the archived mode-1 objective was approximately `0.0096408`.

One accepted analytic step at every construction mode provides the scaling
smoke test below. All tangent columns converged, including all 120 QH mode-5
directions.

| case | mode objectives after one accepted step per stage |
| --- | --- |
| QA | `1.17642 -> 0.23840 -> 0.20441 -> 0.11054 -> 0.09879` |
| QH | `10.57197 -> 4.99554 -> 2.02934 -> 0.64177 -> 0.11575 -> 0.06746` |

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
