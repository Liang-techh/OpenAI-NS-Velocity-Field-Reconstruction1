# Resumable explicit-center interval core

The production core for the center family Z in [0.49,0.51] has reached radial degree 15 of 124. This is a freshly generated family of local coefficient series, with the same accepted fourteen-stage analytic pressure datum as the existing candidate. No Z=0.3 tensor is reused.

The new functional_core_step module accepts the center explicitly and advances A, Uz and P together, retaining pressure, swirl and axial derivative couplings. The saved check compares 360 exact interval endpoints against the batch coupled_rows implementation at Z=0.5 and [0.49,0.51]; all comparisons pass.

The candidate_interval_core driver checkpoints exact interval endpoints after every radial order. Resume checks center, construction parameters, source hashes, accepted pressure identity and row depths. The production run resumed from degree 8 and completed 15; the time budget is checked between complete radial steps.

The current maximum normalized mixed C3 analytic Taylor-tail upper bound is 126612833122270085302199092484174206541772425306035599267306885.94725745593822625758709161. The target is 1e-12 and is not yet met. This bound concerns truncation of the analytic candidate relative to accepted construction data. It does not bound interval-family width, errors in the original construction parameters or the full physical Navier-Stokes residual. Physical radial derivatives require Lambda-dependent conversion.

## Continue

1. Run the driver with --seconds 45 --max-steps 8 repeatedly, preserving frozen source identities, until radial degree 124 and a passing tail target are obtained. Record actual completed degree after each run. A time budget may produce fewer than eight updates.
2. Examine interval-family widths independently of the analytic tail. Narrow the axial atlas when broad parameter intervals prevent useful evaluation; do not describe family variation as numerical error.
3. Add a center-aware comparison adapter consuming these exact rows and amplitude data. Preserve the old adapter's Z=0.3 guard.
4. Propagate comparison and transition equations with axial jets, then construct actual functional moment defects and independent repairs.
5. Extend controlled axial coverage before claiming whole-axis five-moment closure.

This is leading-profile radial coefficient generation. Temporal n-dependent background recursion, heat exterior completion, final admissible stress, oscillatory correction and independent Cartesian residual validation remain incomplete.