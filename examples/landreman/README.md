# Landreman-Paul no-well inputs

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

Each run varies modes through 4 and then 5. It updates the requested output
after every accepted optimizer iteration and also writes `.mode4.json` and
`.mode5.json` snapshots. All files use cuMES's strict, read-back-compatible
input schema. The optional fourth argument caps function evaluations per
stage; omit it for the meow default. Numerical Jacobians require 81 equilibrium
solves at mode 4 and 121 at mode 5 before the first optimizer step, so these
are intentionally long GPU runs.

Two further optional arguments select the first and last mode stage. For
example, append `0 4 4` to run mode 4 only, then use its `.mode4.json` as the
input to a separate run ending in `0 5 5`. This is the supported checkpointed
workflow for scheduler time limits.
