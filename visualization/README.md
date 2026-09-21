# Visualization hub

Choose the field and dataset first. A newer README does not change the coefficients in an older viewer.

## Latest ST063 parent/child comparison

**Fields:** ST061-P, ST063-G1R and ST063-G2R. The default comparison is parent versus G2R. Use the already delivered **`NS_ST063_MATLAB_Comparison.zip`**, extract it, enter the `NS_ST063_MATLAB_Comparison` directory in MATLAB and run:

```matlab
start_here
```

The package includes the actual `st063_models.mat` data, evaluator and `ns_compare_core.m`. The source branch alone does not include every array or dependency. Its source and study record are pinned at `a3d04d3467361bab7c9fc7c8c6aedf31987ee069`:

- [Comparison source](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/blob/a3d04d3467361bab7c9fc7c8c6aedf31987ee069/visualization/matlab/ns_compare_core.m)
- [Original study record](../docs/research_snapshots/ST063.md)
- [Bundle and candidate hashes](../docs/research_catalog.json)

Both sides share physical axes, camera, time, 200 seed locations, absolute thresholds and color limits. The bottom time slider evaluates the original time polynomials. The axis-rotation plot reveals the longer moderate-strength rotation plateau without stretching the geometry. Streamlines are instantaneous curves, not material particle paths.

**Evidence limit:** ST063's new MATLAB UI has not been run natively in that study. Python MAT roundtrip and reference checks passed; the previews were produced by Python and are not MATLAB screenshots. Its default threshold is 0.25, with 0.15 and 0.35 available. Only the lower threshold gives the reported axis-spanning band; high-threshold strong-core continuity is not achieved.

## Viewer bundled on main: ST054

For the existing checked-in ST054-Q2/M3 dataset:

```matlab
addpath('visualization/matlab');
ns_explorer;
```

[Full controls](matlab/README.md) · [Original native execution record](matlab/VERIFICATION.md)

![Original ST054 native MATLAB screenshot](matlab/tests/output/matlab_explorer.png)

This image is from the ST054 native viewer test, not ST063. The existing ST054 data and viewer remain unchanged. Use zero slice offset to inspect the axis; moving seeds or changing the radial display range changes what is seen, not the physics.

## Which visual assets mean what?

| Asset | Scope |
|---|---|
| ST063 complete comparison ZIP | Latest frozen parent/child geometry comparison and export tests |
| Main `matlab/ns_explorer.m` | ST054 interactive spectral field; original native receipts apply only to that code/data |
| Imported ST052 sampled grids | Interpolated velocity data; do not invent missing pressure or full residual |
| Existing morphology/identity tools | Reuse their candidate-bound protocols; do not copy another field's identity receipt |
| Older optimized_v4 scripts and HTML/PNG exports | Historical demonstrations; not ST063 or evidence of current PDE acceptance |

## Controlled comparison rules

Keep physical axis ratio, camera, seed rules, time and absolute thresholds fixed between fields. Show the unfavorable threshold and L2 results as well as the favorable ones. Local moments depend on the observation window; a spanning band can be truncated by that window. No pixel-match percentage or continuous NS certificate is inferred from a plot.

The source report's G2R aspect improvement is 3.5%–7.3%, not the unachieved 10%–20% design aspiration. The central disk remains, and both original momentum thresholds fail. See [results](../docs/RESEARCH_STATUS.md).
