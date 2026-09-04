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

## Run the reproductions

Configure meow with `MEOW_BUILD_CUMES_INTEGRATION=ON`. The four runnable
workflows are fully described by adjacent `*.rundown.json` files:

```bash
build-cumes/cumes_landreman_optimize \
  examples/landreman/qa-construction.rundown.json
build-cumes/cumes_landreman_optimize \
  examples/landreman/qh-construction.rundown.json

build-cumes/cumes_landreman_optimize \
  examples/landreman/qa-refinement.rundown.json
build-cumes/cumes_landreman_optimize \
  examples/landreman/qh-refinement.rundown.json
```

Validate and inspect a plan without launching an equilibrium solve:

```bash
build-cumes/cumes_landreman_optimize --dry-run \
  examples/landreman/qa-construction.rundown.json
```

For scheduler-sized continuation, select an ordered subset by step name and
override the output location:

```bash
build-cumes/cumes_landreman_optimize \
  --first-step mode-1 --last-step mode-2 \
  --output ../opt-qa/latest.json \
  --iteration-directory ../opt-qa/steps \
  examples/landreman/qa-construction.rundown.json
```

The construction rundowns start at the sparse analytic inputs. QA enables the
archived mean-iota target 0.42 and boundary modes 1--4; QH omits the iota
residual and continues through mode 5. Both use uniform radial QS weights and
the archived transform resolutions and equilibrium work allowances.

The refinement rundowns start at `qa_start.json` and `qh_start.json`, run
boundary modes 4 and 5 at the input equilibrium resolution, omit the iota
residual, and apply the paper's final radial weight ramps (1 to 30 for QA and
1 to 2 for QH).

The qualified reproduction policy remains cold one-sided finite differences:
relative/absolute steps are `3.1622776601683794e-3` / `1e-7` for QA
construction, `1e-3` / `1e-7` for QH construction, and `1e-5` / `1e-9`
for both refinements. QH uses the retained Catmull-Rom radial transfer.

Every accepted iteration updates `output.path`. When
`output.iteration_directory` is nonempty, the run also stores paired
`*-input.json` and `*-equilibrium.bin` artifacts. Construction checkpoints
use `.construction.modeM.json`; refinement checkpoints use `.modeM.json`.

All numerical policy—including experimental parallel, analytic, Broyden, and
multi-accuracy variants—is expressed in the rundown rather than inferred from
a method name. The full field reference and CLI override list are documented
in [the rundown schema](../../docs/relaxation-rundown.md).
