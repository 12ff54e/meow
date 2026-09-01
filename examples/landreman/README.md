# Landreman-Paul no-well inputs

`qa_analytic.json` and `qh_analytic.json` are the sparse analytic boundaries
used at the beginning of the original construction runs 021 and 039. Their
boundary formulas are

```text
QA: R = 1 + 0.2 cos(theta),
    Z =     0.2 sin(theta),                         nfp = 2

QH: R = 1 + 0.2 cos(theta) + 0.2 cos(4 zeta),
    Z =     0.2 sin(theta) - 0.2 sin(4 zeta),      nfp = 4
```

Regenerate them from the calculation subtree with:

```bash
CALC="$SUPPLEMENT_ROOT/calculations/20210704-01-simsopt_new_quasisymmetry_metric"
QA_CONSTRUCTION="$CALC/20210704-01-021_QA_multiple_surfaces_relStepScan_A6_iota0.42"
QH_CONSTRUCTION="$CALC/20210704-01-039_QH_multiple_surfaces_relStepScan_A8"

scripts/vmec_namelist_to_cumes_json.py \
  "$QA_CONSTRUCTION/input.nfp2_QA" examples/landreman/qa_analytic.json \
  --minimum-ftol 1e-12 --minimum-niter 6000 \
  --wout-axis "$QA_CONSTRUCTION/abs_step_1.00e-07_rel_step_3.16e-03_forward/wout_nfp2_QA_000_000000.nc"

scripts/vmec_namelist_to_cumes_json.py \
  "$QH_CONSTRUCTION/input.nfp4_QH" examples/landreman/qh_analytic.json \
  --minimum-ftol 1e-12 --minimum-niter 6000 \
  --wout-axis "$QH_CONSTRUCTION/abs_step_1.00e-07_rel_step_1.00e-03_forward/wout_nfp4_QH_000_000000.nc"
```

The imported axis is a numerical initial-guess predictor, not part of the
analytic fixed boundary. The iteration-zero `wout` files were written at
`ntor=3`, so the importer explicitly pads their unused higher axis modes with
zero to match the source input's `ntor=5`.

`qa.json` and `qh.json` are cuMES inputs generated from the final no-magnetic-
well configurations in the paper supplement. Set `SUPPLEMENT_ROOT` to the
downloaded `20211102-01-precise_quasisymmetry_zenodo` directory and regenerate
them with:

```bash
scripts/vmec_namelist_to_cumes_json.py \
  "$SUPPLEMENT_ROOT/configurations/new_QA/input.20210704-01-050_weight_3.00e+01_rel_step_1.00e-05_forward_nfp2_QA" \
  examples/landreman/qa.json --minimum-ftol 1e-12 \
  --wout-axis "$SUPPLEMENT_ROOT/configurations/new_QA/wout_20210704-01-050_weight_3.00e+01_rel_step_1.00e-05_forward_nfp2_QA.nc"

scripts/vmec_namelist_to_cumes_json.py \
  "$SUPPLEMENT_ROOT/configurations/new_QH/input.20210704-01-063_nfp4_QH_A8_weight_2.00e+00_rel_step_1.00e-05_forward_B1" \
  examples/landreman/qh.json --minimum-ftol 1e-12 \
  --wout-axis "$SUPPLEMENT_ROOT/configurations/new_QH/wout_20210704-01-063_nfp4_QH_A8_weight_2.00e+00_rel_step_1.00e-05_forward_B1.nc"
```

`qa_start.json` and `qh_start.json` are the actual initial boundaries of the
paper's final mode-4/mode-5 refinement runs. Regenerate them from the
`calculations/20210704-01-simsopt_new_quasisymmetry_metric` subtree:

```bash
CALC="$SUPPLEMENT_ROOT/calculations/20210704-01-simsopt_new_quasisymmetry_metric"
QA_RUN="$CALC/20210704-01-050_QA_weightRampScan_noIotaTarget_refiningBestFrom021/weight_3.00e+01_rel_step_1.00e-05_forward"
QH_RUN="$CALC/20210704-01-059_QH_A8_refiningBestFrom039_weightRampScan/weight_2.00e+00_rel_step_1.00e-05_forward"

scripts/vmec_namelist_to_cumes_json.py \
  "$QA_RUN/input.nfp2_QA" examples/landreman/qa_start.json \
  --minimum-ftol 1e-12 --minimum-niter 6000 \
  --wout-axis "$QA_RUN/wout_nfp2_QA_000_000000.nc"

scripts/vmec_namelist_to_cumes_json.py \
  "$QH_RUN/input.nfp4_QH" examples/landreman/qh_start.json \
  --minimum-ftol 1e-12 --minimum-niter 10000 \
  --wout-axis "$QH_RUN/wout_nfp4_QH_000_000000.nc"
```

The tolerance clamp and axis import are documented in
`docs/landreman-paul-reproduction.md`. They adapt solver-control and initial-
guess conventions; the fixed boundary coefficients are copied unchanged. The
QH evaluation requests cuMES's retained Catmull-Rom radial transfer because
the default B-spline transfer invalidates the 100-to-150 surface transition
for this configuration.

After configuring meow with `MEOW_BUILD_CUMES_INTEGRATION=ON`, run the final
two-stage refinements with:

```bash
build-cumes/cumes_landreman_optimize \
  examples/landreman/qa_start.json qa qa_optimized.json
build-cumes/cumes_landreman_optimize \
  examples/landreman/qh_start.json qh qh_optimized.json
```

Run the earlier analytic construction stages with:

```bash
build-cumes/cumes_landreman_optimize \
  examples/landreman/qa_analytic.json qa-construction qa_constructed.json
build-cumes/cumes_landreman_optimize \
  examples/landreman/qh_analytic.json qh-construction qh_constructed.json
```

The QA construction enables the archived mean-iota target 0.42 and continues
through maximum boundary modes 1, 2, 3, and 4. The QH construction has no iota
target and continues through modes 1, 2, 3, 4, and 5. Both use uniform radial
QS weights. Their equilibrium transform resolutions follow the archived
drivers: 3, 5, 6, 6 for QA and 3, 5, 6, 6, 6 for QH. The existing `qa` and
`qh` case names remain the separate final-refinement policy, which removes the
QA iota residual and applies the final edge-weight ramps.

At construction resolution 6, cuMES is allowed 10,000 iterations per radial
stage for mode 3 and 30,000 for QA mode 4. The first QA mode-3 transition
reached the previous 6,000-iteration cap at residuals
`(1.72e-12, 8.44e-13, 1.64e-12)`. At mode 4, several valid finite-difference
samples reached the 10,000-iteration cap near `7e-9`; these samples must
converge rather than be mistaken for infeasible trial boundaries. At 30,000
iterations, several of these samples reproducibly stall near
`(1.20e-12, 4.47e-13, 1.27e-12)`, so QA mode 4 uses a qualified `2e-12` gate.
This remains ten orders of magnitude below the construction objective and is
recorded in every emitted stage input.

Construction checkpoints are named
`.construction.modeM.json`, and archived equilibria use
`construction-modeM_step_NNNN-equilibrium.bin`. This keeps construction mode
4 distinct from the later refinement mode 4 when both workflows share an
output directory.

Each run varies modes through 4 and then 5. It updates the requested output
after every accepted optimizer iteration and also writes `.mode4.json` and
`.mode5.json` snapshots. All files use cuMES's strict, read-back-compatible
input schema. The optional fourth argument caps function evaluations per
stage; omit it for the meow default. Numerical Jacobians require 81 equilibrium
solves at mode 4 and 121 at mode 5 before the first optimizer step, so these
are intentionally long GPU runs.

Two further optional arguments select the first and last mode stage. For
example, append `0 4 4` to a refinement run to run mode 4 only, then use its
`.mode4.json` as the input to a separate run ending in `0 5 5`. Appending
`0 1 2` to a construction run selects its first two stages. This is the
supported checkpointed workflow for scheduler time limits.

The next optional argument limits accepted optimizer iterations per stage. A
final directory argument stores `modeM_step_NNNN-input.json` and
`modeM_step_NNNN-equilibrium.bin` for the initial state (`NNNN=0000`) and
every accepted iteration. For example, a 100-step mode-4 QA run is:

```bash
build-cumes/cumes_landreman_optimize \
  examples/landreman/qa_start.json qa ../opt-qa-test/latest.json \
  0 4 4 100 ../opt-qa-test
```
