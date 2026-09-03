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

### End-to-end tangent benchmark

On 2026-09-02, the analytic and finite-difference Jacobian paths were run
sequentially from the same checked-in sparse analytic boundaries through every
construction stage. Both methods used the same Release, double-precision
executable, target, continuation schedule, equilibrium tolerances, accepted
axis policy, and deterministic QA chiral seed. Iteration-equilibrium output was
disabled so that result serialization did not enter the timing. Every stage in
all four runs stopped because the TRF step satisfied `xtol`; none reached its
function-evaluation safeguard.

The machine used an NVIDIA TITAN Xp with driver 580.173.02. The source states
were cuMES `4bf4d881943a73835a090702689f5fa0e4eb3789` and meow
`8eac01e82b8d1eaa508d53afd8acae2b2f223740`. The raw wall-clock results are:

| case | Jacobian | wall time (s) | speed relative to finite differences | final objective |
| --- | --- | ---: | ---: | ---: |
| QA | forward equilibrium tangent | 473.29 | 3.248x | `6.29151469581e-4` |
| QA | forward finite difference | 1537.26 | 1.000x | `5.94683530877e-7` |
| QH | forward equilibrium tangent | 1387.86 | 2.351x | `1.28435967198e-3` |
| QH | forward finite difference | 3262.78 | 1.000x | `5.36939265810e-5` |

Thus the current tangent implementation reduces wall time by 69.2% for QA and
57.5% for QH. The QH work counters show how it does so even though its less
accurate derivatives require more accepted optimizer iterations:

| QH path | accepted iterations | equilibrium evaluations | nonlinear equilibrium iterations | tangent Jacobians | tangent linear iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| tangent | 218 | 253 | 474,719 | 223 | 635,648 |
| finite difference | 88 | 4,303 | 12,203,424 | 0 | 0 |

This is a performance result, not yet an optimization-equivalence result. The
tangent final objective is 1,058 times the finite-difference objective for QA
and 23.9 times for QH. In particular, the tangent stages can satisfy `xtol`
because their local model predicts a small useful step while finite differences
continue to find strong descent directions. The QH tangent logs also report
worst normalized linear residual diagnostics near 3 in the higher-mode stages.
The forward tangent therefore needs a column-by-column accuracy investigation
and tighter/adaptive linear-solve qualification before its raw speedup can be
treated as acceleration of the same optimization.

The timed commands, differing only in the final Jacobian selector, were:

```bash
./build/meow/cumes_landreman_optimize \
  ../meow/examples/landreman/qa_analytic.json qa-construction RESULT.json \
  0 1 4 0 '' analytic
./build/meow/cumes_landreman_optimize \
  ../meow/examples/landreman/qa_analytic.json qa-construction RESULT.json \
  0 1 4 0 '' finite-difference

./build/meow/cumes_landreman_optimize \
  ../meow/examples/landreman/qh_analytic.json qh-construction RESULT.json \
  0 1 5 0 '' analytic
./build/meow/cumes_landreman_optimize \
  ../meow/examples/landreman/qh_analytic.json qh-construction RESULT.json \
  0 1 5 0 '' finite-difference
```

The result JSON, per-stage checkpoints, full QH logs, and `/usr/bin/time`
records are retained in `../benchmark-forward-tangent-20260902` relative to
the cuMES checkout.

### Tangent-correction plan

The benchmark above is not accepted as an equivalent optimization. The
correction will proceed in the following committed slices, without moving any
target definition into cuMES:

1. Extend the nonlinear Jacobian oracle from one QH mode-1 column to selected
   and worst-case columns at every construction mode. Use the same accepted
   axis predictor and cold nonlinear solve policy as the optimizer. Record
   complete residual-column error, objective directional-derivative error,
   and the error in the optimizer-relevant gradient `J^T r`.
2. Decompose each failing column at the cuMES/meow boundary. First compare the
   solved spectral and published equilibrium-field tangent with the centered
   nonlinear equilibrium difference. Then apply meow's target JVP to that
   nonlinear field difference. This distinguishes an equilibrium tangent
   error from a QS/aspect/iota chain-rule error.
3. Correct the failing layer. For a cuMES linear-solve error, qualify the true
   residual after preconditioning and use an adaptive tolerance tight enough
   for the requested optimizer column; if the near-null equilibrium gauge
   still changes physical fields, constrain it in the same metric as the
   converged nonlinear trajectory. For a meow chain-rule error, correct and
   unit-test only the optimizer-side target JVP.
4. Add an optimizer guard that rejects an analytic Jacobian whose linear
   solves or oracle diagnostics do not meet the qualified accuracy contract.
   A nominal GMRES `converged` flag alone is not sufficient.
5. Re-run QA and QH from their analytic boundaries. Acceptance requires every
   stage to terminate normally, the final objective to agree with the matched
   finite-difference result to a stated tolerance, and the boundary/target
   trajectory to remain physically equivalent. Only then will the wall-clock
   ratio be reported as an optimization speedup.

The diagnostic and target work belongs to meow. Any required change to the
equilibrium residual derivative, gauge, or tangent linear solver belongs to
cuMES and must preserve its existing nonlinear frozen trajectories.

### Tangent correction and hot-restart outcome

The correction plan was executed before accepting any speed claim. The
expanded oracle showed that meow's target chain rule is not the limiting
error: its worst objective-chain discrepancy was about 0.118%. The cuMES
implicit response, however, can satisfy the linearized force equations while
following a different member of a near-null equilibrium branch than the cold
nonlinear black box. Selected complete residual columns differed by roughly
10--14%. Tightening GMRES, reorthogonalizing, removing its preconditioner, and
changing restart lengths did not remove that branch mismatch.

A tight tangent Jacobian with black-box `m=0` columns improved QA to
`4.20074375796e-6` in 2208.79 s. A subsequent 963.35 s cold polish reached
`5.5534985509e-7`, but the combined 3172.14 s was 2.06 times the 1537.26 s
cold control. Thus the corrected forward tangent can recover target quality,
but it does not accelerate this optimization.

The second experiment replaced the implicit response by explicit final-grid
hot-restart differences. Each column uses an absolute step floor of `1e-4`;
a failed hot solve retries the full cold multigrid path, and a perturbation on
the edge of the feasible equilibrium domain may switch from forward to
backward difference. Compatible continuation stages also retain the accepted
equilibrium for their first evaluation. These policies stay in meow; cuMES
still only solves the requested equilibrium.

Fresh analytic-boundary runs on the same TITAN Xp and source states cuMES
`4bf4d881943a73835a090702689f5fa0e4eb3789` and meow
`1e433cd48660aa35f533491f6d9a06251624e700` produced:

| case | Jacobian | wall time (s) | raw speed | final objective | objective / cold |
| --- | --- | ---: | ---: | ---: | ---: |
| QA | hot-restart difference | 1197.67 | 1.284x | `6.64013177288e-6` | 11.17x |
| QA | cold finite difference | 1537.26 | 1.000x | `5.94683530877e-7` | 1.00x |
| QH | hot-restart difference | 1945.78 | 1.677x | `2.77430063630e-4` | 5.17x |
| QH | cold finite difference | 3262.78 | 1.000x | `5.36939265810e-5` | 1.00x |

The raw wall-time reductions are 22.1% for QA and 40.4% for QH, but neither
run is an equivalent-result acceleration. A finer `1e-5` hot step resumed from
the QA endpoint in 49.94 s and stopped at `6.45131354312e-6`, confirming that
the discrepancy is a local-minimum change rather than just the final step
floor. A full cold mode-4 polish from the hot QA endpoint took another 1331.80
s and stopped at `1.30948307823e-6`; the combined 2529.47 s is only 0.608x the
cold-control speed and remains 2.20 times its objective.

Consequently there is **no qualified QA or QH speedup** from the current
forward tangent or hot-restart Jacobian. `finite-difference` remains the CLI
default and reproduction oracle; `analytic` and `hot-finite-difference` are
explicit experimental selectors. The complete hot logs, result JSON,
per-stage checkpoints, and timing records are retained under
`../benchmark-forward-tangent-qualified-20260902/qa-hot-v4` and
`../benchmark-forward-tangent-qualified-20260902/qh-hot-v2` relative to the
cuMES checkout. The QA polish is in the sibling `qa-hot-v4-polish` directory.

### Lambda-independence isolation plan

The QS expression has no explicit lambda term, so lambda must not be assumed
to explain the tangent discrepancy. The isolation is performed in this order:

1. Inject the two stored `(m,n)=(0,0)` lambda coefficients, whose sine basis
   functions vanish identically. Materialize each direction through cuMES and
   require every published magnetic field, flux profile, aspect derivative,
   and meow target JVP to be zero to roundoff. This tests the actual public
   bridge rather than only the spectral input.
2. Do not classify a general nonzero lambda harmonic as a pure gauge: it
   changes the straight-field-line mapping unless accompanied by the matching
   coordinate relabelling of `R` and `Z`. Such harmonics are diagnostics only.
3. For every QH mode-1 boundary column, compare the implicit tangent with a
   converged final-grid hot-restart difference separately for all six state
   families, the target-facing half-grid fields, the five flux profiles,
   aspect ratio, and the complete target residual. Record both reference norms
   and normalized errors so a nearly zero quantity cannot masquerade as a
   large physical discrepancy.
4. Attribute the first material discrepancy to the earliest layer where it
   appears: spectral equilibrium response, published-field materialization,
   or meow target chain rule. Only that owning repository is changed.

Generated logs and fixtures for this investigation belong under `../tmp`, not
the system temporary directory.

The isolation produced the following result on the QH mode-1 analytic start:

- Both stored zero-basis directions, `LMNSC(0,0)` and `LMNCS(0,0)`, give
  exactly zero residual JVP, published field/profile tangent, aspect tangent,
  and target JVP. The actual lambda gauge is therefore absent from the target
  bridge, as required.
- Meow's target chain applied to the nonlinear equilibrium difference agrees
  with the fully re-evaluated target column to at worst 0.116%. The target
  expression and its JVP are not the source of the discrepancy.
- The implicit and nonlinear responses first differ in the spectral state,
  including `R` and `Z`, and consequently in target-facing fields. Worst
  finite-reference errors among the eight columns are 29.4% for `B`, 32.3%
  for `B dot grad(B)`, 26.2% for
  `B cross grad(psi_p) dot grad(B)`, and 39.4% for iota. `G` remains within
  1.48%. Quantities whose nonlinear reference norm is nearly zero, such as
  the QH-start aspect derivative for `rbc(0,1)` and `I`, are not classified by
  relative error.
- The complete target-residual column differs by 3.09--13.42%, while the
  objective directional derivative is much less sensitive: at worst 0.516%
  outside the two optimizer-designated `m=0` black-box columns. This explains
  why a gradient-like check can pass while TRF's `J^T J` model changes.
- Replacing only lambda or only `R/Z` between the two responses increases the
  target-column error to 29--158%. Their separate target contributions are
  almost antiparallel, with cosines from -0.988 to -0.997 and only 3.64--7.63%
  of their summed norms left after cancellation. The mismatch is therefore a
  coupled coordinate/shape response, not an isolated lambda contribution.
- Both responses remain close to the same linearized force nullspace. The
  implicit combined residual is about `5e-6` of the boundary forcing by
  construction; the one-sided nonlinear-difference residual is at most
  `2.36e-4` of that forcing. Their difference is an approximate coupled null
  direction of the discrete equilibrium Jacobian.

The first discrepant owner is consequently cuMES's equilibrium response, not
meow's lambda-independent target. The next correction should characterize and
constrain the coupled null direction using the nonlinear continuation metric;
zeroing lambda or changing the target would destroy the observed physical
cancellation. The complete diagnostic is
`../tmp/lambda-isolation/qh-mode1-fields-final.log`.

### Constraint-tangent wiring correction and repeated benchmark

The branch investigation first exposed a concrete cuMES implementation defect,
not a missing lambda gauge condition. The dual inverse transform produced the
constraint fields `rCon` and `zCon` into tangent-owned buffers, but the
subsequent `ConstraintOperator` read separate internal buffers that had never
received those values. The nonlinear equilibrium path already connected these
views correctly. cuMES commit `8c61395` connects the dual transform directly to
the constraint views and removes the unused buffers; commit `f1a3ed7` also
initializes the dual de-aliasing axis scratch discovered by Compute Sanitizer.
The complete cuMES suite passes 96/96 tests, including sanitizer wrappers, and
meow passes 10/10 tests against the corrected installed library.

This result also clarifies why lambda appears in the equilibrium Jacobian even
though the QS scalar has no explicit lambda term. If the converged equilibrium
is `F(u,x)=0`, for state `u=(R,Z,lambda)` and boundary variables `x`, the
optimizer needs

```text
du/dx = -F_u^-1 F_x
dT/dx = T_x + T_u du/dx.
```

A pure coordinated relabelling must vanish in the physical target after the
second line is applied. A general lambda harmonic with `R` and `Z` held fixed,
however, is not such a relabelling: it changes the straight-field-line mapping
and the field components used by `F`. Removing lambda from `F_u` would alter
the coupled `R/Z` response rather than enforce target invariance. The two
stored zero-basis lambda directions still materialize to exactly zero through
the public equilibrium and target bridge.

After the wiring correction, the centered Solovev restart oracle differs from
the implicit response by 0.228% in aggregate spectral state and 0.514% in the
published fields. At the harder QH mode-1 start, target-residual columns differ
by 2.95--11.68%, while the hybrid objective directional derivative differs by
at most 0.705%. Tightening the relative GMRES tolerance from `5e-6` to `1e-8`
still leaves a 7.77% worst target-column difference. The remaining response is
therefore not qualified by the 1% column gate, but the corrected one-step smoke
trajectories improve substantially:

| case | mode objectives after one accepted step per stage |
| --- | --- |
| QA | `1.17642 -> 0.25383 -> 0.11965 -> 0.06565 -> 0.05646` |
| QH | `10.57197 -> 4.43401 -> 0.69221 -> 0.17152 -> 0.13730 -> 0.02254` |

Fresh full-convergence runs then used cuMES `f1a3ed7`, meow `fc7e66b`, the
same TITAN Xp, analytic starting boundaries, targets, stage schedule, and cold
finite-difference controls recorded above. All stages stopped normally by
`ftol` or `xtol`:

| case | Jacobian | wall time (s) | speed / cold | final objective | objective / cold |
| --- | --- | ---: | ---: | ---: | ---: |
| QA | corrected equilibrium tangent | 2628.35 | 0.585x | `1.34526695087e-6` | 2.26x |
| QA | cold finite difference | 1537.26 | 1.000x | `5.94683530877e-7` | 1.00x |
| QH | corrected equilibrium tangent | 3516.16 | 0.928x | `9.71167722426e-5` | 1.81x |
| QH | cold finite difference | 3262.78 | 1.000x | `5.36939265810e-5` | 1.00x |

The corrected tangent reduces the earlier tangent endpoint by factors of 468
for QA and 13.2 for QH, but it is 1.71 times slower than the QA cold control
and 1.08 times slower than the QH cold control. Since the endpoints also do not
match the finite-difference objectives, there remains **no qualified analytic
Jacobian speedup**. The default remains `finite-difference`; `analytic` remains
experimental. Complete corrected logs, result JSON, checkpoints, and timing
records are retained under
`../benchmark-forward-tangent-corrected-20260903/{qa-analytic,qh-analytic}`
relative to the cuMES checkout. The focused diagnostics are under
`../tmp/lambda-isolation` and the corrected smoke logs under
`../tmp/tangent-fix-smoke`.

### Accepted-restart and secant-reuse experiment plan

The next speed experiment retains the qualified finite-difference target
Jacobian rather than forming one equilibrium tangent solve per boundary
coefficient. It is split into independently reversible steps:

1. Keep a snapshot of the last accepted equilibrium, separate from the most
   recently evaluated trial. Initialize nearby trust-region trials from that
   snapshot and initialize every Jacobian perturbation from the unperturbed
   equilibrium. Use the existing Landreman relative and absolute difference
   steps exactly; do not reuse the earlier `1e-4` hot-difference floor. A failed
   restart retries the complete cold multigrid solve.
2. Compare warm and cold residual/Jacobian columns before timing. The warm path
   is rejected if it changes a materially nonzero target column by more than
   1%, changes an accepted objective trajectory materially, or fails the same
   equilibrium convergence gates.
3. Add an optional good-Broyden update to meow's TRF implementation. The
   default continues to rebuild every Jacobian. The experimental Landreman
   path refreshes the finite-difference Jacobian periodically and immediately
   after poor trust-region agreement, a large secant defect, or any invalid
   update.
4. Run one-accepted-step QA and QH construction smoke tests for accepted
   restart alone, then for accepted restart plus secant reuse. Only a variant
   that preserves the cold objective trajectory proceeds to full convergence.
   Report wall time, equilibrium evaluations and iterations, exact Jacobian
   builds, secant updates, and final objective against the existing cold
   controls.

All optimizer/restart policy remains in meow. cuMES continues to solve the
supplied equilibrium request without knowledge of the target or optimization
method. Generated inputs, equilibria, logs, and timing files belong under the
sibling `tmp/` or benchmark directories, never in either source repository.

### Accepted-restart and secant-reuse experiment outcome

The accepted-equilibrium restart was fast but failed its trajectory gate. A
four-stage, one-accepted-step QA smoke took 98.76 s instead of 201.60 s for the
matching cold control, and every restarted Jacobian perturbation converged
without fallback. Nevertheless, the objectives diverged immediately: mode 1
gave `0.246556030419` warm versus `0.302760535277` cold, and mode 2 gave
`0.138212808326` warm versus `0.119519928704` cold. Exact production
finite-difference steps therefore do not prevent the restart from selecting a
different weak equilibrium branch. `warm-finite-difference` remains an
explicit diagnostic and is not a qualified default.

Safeguarded good-Broyden reuse was then isolated from restarting: every exact
refresh remained a complete cold finite-difference Jacobian. The general TRF
default still refreshes after every accepted step. The experimental Landreman
`broyden` selector uses a five-step maximum age, refreshes when the trust ratio
is below 0.1 or the relative secant defect exceeds 0.1, and reports every
decision. Mode-1 convergence checks gave:

| case | path | wall time (s) | objective | equilibrium evaluations | nonlinear iterations |
| --- | --- | ---: | ---: | ---: | ---: |
| QA | safeguarded Broyden | 34.83 | `9.62278667219e-3` | 136 | 169,745 |
| QA | cold refresh every step | 53.39 | `9.62220283167e-3` | 214 | 264,315 |
| QH | safeguarded Broyden | 40.28 | `0.14050177654` | 161 | 193,894 |
| QH | cold refresh every step | 60.92 | `0.140503531818` | 242 | 295,662 |

Both cases were therefore continued through their complete analytic-boundary
construction on the same TITAN Xp. Every stage stopped normally by `xtol`:

| case | path | wall time (s) | speed / cold | final objective | objective / cold |
| --- | --- | ---: | ---: | ---: | ---: |
| QA | safeguarded Broyden | 1333.45 | 1.153x | `5.09650574804e-7` | 0.857x |
| QA | cold refresh every step | 1537.26 | 1.000x | `5.94683530877e-7` | 1.000x |
| QH | safeguarded Broyden | 3016.60 | 1.082x | `3.84152264277e-5` | 0.715x |
| QH | cold refresh every step | 3262.78 | 1.000x | `5.36939265810e-5` | 1.000x |

QA accepted seven secant updates across all four stages. QH accepted 15 across
five stages; its total equilibrium evaluations fell from 4,303 to 3,979 and
nonlinear iterations from 12,203,424 to 11,052,763. The last QA stage accepted
no secant update, so its endpoint was terminated using exact cold Jacobians.
The sparse earlier updates alter the continuation path and reach lower local
minima in both cases. These are valid improvements in time and achieved
objective, but they are not identical-trajectory speedups. Complete artifacts
are under `../benchmark-broyden-20260904/{qa,qh}` relative to the cuMES
checkout; focused controls are under `../tmp/warm-fd-smoke` and
`../tmp/broyden-smoke`.

### Jacobian-scaling experiment plan

The next independent experiment addresses the unit boundary-variable scaling
used by the current Landreman driver. Meow will add an opt-in TRF policy
equivalent to Jacobian column-norm scaling: after each exact Jacobian or secant
update, the characteristic scale of a variable is the reciprocal of its
column norm, with nondecreasing inverse scales to avoid expanding a direction
as the solve progresses. Fixed user scales and Jacobian-derived scales are
mutually exclusive; the default remains fixed unit scaling.

The first QA/QH tests use cold finite-difference Jacobians after every accepted
step, so scaling is the only changed variable. Mode-1 convergence must retain
a comparable or better objective and reduce wall time or accepted iterations
before a complete construction run is attempted. Only after independent
qualification may Jacobian scaling be combined with safeguarded Broyden reuse.
Artifacts belong under `../tmp/jacobian-scale-smoke` and, if promoted to full
construction, `../benchmark-jacobian-scale-20260904`.

### Jacobian-scaling experiment outcome

The mode-1 checks initially favored Jacobian-derived scaling. QA fell from
18 to 12 accepted steps and took 40.74 s rather than 53.39 s, ending at
`9.62192289293e-3`; QH fell from 21 to 17 steps and took 48.42 s rather than
60.92 s, ending at `0.140500191495`. Both objectives were slightly lower than
their unit-scale controls.

The complete QA construction reversed that result. Scaling took 2,345.81 s,
versus 1,537.26 s for the cold unit-scale control: **0.655x as fast**, or
52.6% more wall time. The final objective, `5.8980706949e-7`, was 0.82% lower
than the control's `5.94683530877e-7`, so this is a valid but slower trajectory,
not a failed equilibrium. The regression arose mainly in the later stages:

| maximum mode | accepted steps | equilibrium evaluations | nonlinear iterations | objective |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 12 | 158 | 190,174 | `9.62192289293e-3` |
| 2 | 18 | 523 | 1,724,222 | `2.41137450101e-4` |
| 3 | 32 | 1,658 | 4,132,021 | `4.24994425983e-6` |
| 4 | 16 | 1,414 | 2,809,064 | `5.8980706949e-7` |

Because standalone scaling failed the complete QA speed gate decisively, it
was not promoted to a multi-hour complete QH construction and it will not be
combined with Broyden as the next primary candidate. The opt-in selector is
retained for controlled experiments. Full QA artifacts are under
`../benchmark-jacobian-scale-20260904/qa`; the QA/QH mode-1 controls are under
`../tmp/jacobian-scale-smoke`.

### Further optimizer-performance experiment plan

The performance investigation will continue after Broyden and scaling. Each
candidate stays opt-in until both QA and QH have been checked, and every timed
comparison records the final full-accuracy objective as well as wall time,
equilibrium evaluations, and nonlinear iterations. The experiments are
ordered to preserve a useful control for each additional change:

1. If standalone Jacobian scaling reaches a competitive full-construction
   endpoint, combine it with the already qualified safeguarded Broyden reuse.
   First run complete mode-1 QA/QH checks, then full constructions only when
   objective and equilibrium-validity gates pass.
2. Measure a two-accuracy equilibrium schedule. Optimize initially with a
   relaxed cuMES force tolerance, then always polish the same stage at the
   qualified `1e-12` tolerance and use only that full-accuracy state for stage
   continuation and final reporting. Tolerance scheduling is an optimizer
   policy in meow; cuMES continues to solve the request it receives.
3. Revisit restart acceleration without carrying the full lambda state. A
   candidate restart may reuse physical geometry while rebuilding or
   canonicalizing the lambda gauge; it must pass direct residual/Jacobian
   comparisons before timing because the full-state restart already selected
   a different weak branch.
4. If the preceding methods leave substantial exact-Jacobian cost, test more
   aggressive secant refresh policies on isolated late stages before full
   runs. Reject policies that save evaluations only by terminating at a worse
   full-accuracy objective.

Generated data remains outside the repositories under dated sibling `tmp/`
or benchmark directories. A negative result is retained in this document so
that it is not silently promoted or repeated.

### Two-accuracy preliminary result

The first implementation restarted the qualified polishing center from the
relaxed equilibrium while forming its finite-difference columns with cold
solves. QH exposed this as the same mixed-branch inconsistency seen in the
earlier restart experiment: polishing accepted no step. The corrected policy
therefore uses cold `1e-12` solves for both the polishing center and all of its
Jacobian columns. The relaxed equilibrium influences only the boundary passed
to polishing, not its equilibrium branch.

With that correction, complete mode-1 runs passed for both cases:

| case | path | wall time (s) | speed / cold | polished objective |
| --- | --- | ---: | ---: | ---: |
| QA | `1e-9` optimize + cold `1e-12` polish | 32.44 | 1.646x | `9.6221676106e-3` |
| QA | cold `1e-12` control | 53.39 | 1.000x | `9.62220283167e-3` |
| QH | `1e-9` optimize + cold `1e-12` polish | 42.58 | 1.431x | `0.140502193535` |
| QH | cold `1e-12` control | 60.92 | 1.000x | `0.140503531818` |

The relaxed QA phase used 67,106 nonlinear equilibrium iterations and its
three-step polish used 76,461; the relaxed QH phase used 114,658 and its
three-step polish used 79,109. Both polished objectives are slightly lower
than their controls, so the method proceeds to complete QA and QH
constructions. Preliminary artifacts are under
`../tmp/two-accuracy-mode1-cold-polish`.

### Complete two-accuracy QA outcome and capped follow-up plan

Fully converging both accuracy phases did not preserve the mode-1 speedup.
The complete QA construction took 2,299.82 s, versus 1,537.26 s for cold
finite differences and 1,333.45 s for safeguarded Broyden. Its final
full-accuracy objective was excellent, `4.59025456054e-7` (22.8% below the
cold control and 9.9% below Broyden), but this is a quality improvement at
0.669x cold speed, not an optimization-performance improvement. The eight
relaxed/polish phase records are:

| mode | phase | steps | equilibrium evaluations | nonlinear iterations | objective |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | relaxed | 13 | 167 | 67,106 | `9.59731261565e-3` |
| 1 | polish | 3 | 62 | 76,461 | `9.6221676106e-3` |
| 2 | relaxed | 17 | 475 | 399,109 | `1.65142242228e-4` |
| 2 | polish | 4 | 149 | 322,052 | `1.59967039239e-4` |
| 3 | relaxed | 12 | 674 | 560,618 | `6.34023072171e-6` |
| 3 | polish | 17 | 929 | 1,915,145 | `3.14248385128e-6` |
| 4 | relaxed | 21 | 1,808 | 1,713,035 | `1.04606370736e-6` |
| 4 | polish | 17 | 1,481 | 3,124,445 | `4.59025456054e-7` |

The late relaxed phases spend many exact-Jacobian builds on diminishing
returns. The next isolated variant therefore caps every relaxed phase at six
accepted steps—the first six steps contain the large reductions in the QA
trace—and uses the already safeguarded good-Broyden policy during the cold
`1e-12` polish. Mode-1 QA/QH controls come first. Only matching-or-better
polished objectives with a timing win proceed to complete constructions.
The complete uncapped QA artifacts are under
`../benchmark-two-accuracy-20260904/qa`.

The fixed six-step cap failed its first mode-1 QA gate. The relaxed phase took
six steps and 26,567 nonlinear iterations, but the cold Broyden polish then
needed 31 steps, 297 equilibrium evaluations, and 367,766 nonlinear
iterations. Total wall time was 81.44 s, slower than the 53.39 s cold control,
even though the objective `9.6219340894e-3` was valid. Only three of the 31
polish steps passed the conservative secant safeguards. No QH run was made
for this rejected fixed-cap rule.

The next cap is based on measured progress instead of stage count: after at
least eight relaxed steps, switch to polishing when the cumulative objective
improvement over the latest three accepted steps is below 1%. On the uncapped
QA trace this retains all useful mode-1 and mode-3 progress, trims the final
two mode-2 steps, and removes nine plateau steps from mode 4. The qualified
polish retains safeguarded Broyden reuse. The rejected fixed-cap artifacts are
under `../tmp/two-accuracy-broyden-mode1`.

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
