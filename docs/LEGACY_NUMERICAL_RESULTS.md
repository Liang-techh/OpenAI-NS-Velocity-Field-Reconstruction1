# Research results and limitations

Updated 2026-09-21 from pinned ST061 and ST063 study records. This is a navigation and reporting update, not a new fit or a new scientific validation. The original full momentum maximum and spatial volume L2 targets remain 0.001 and remain unmet.

## Latest geometry experiment: ST063

Source: [unaltered ST063 record](research_snapshots/ST063.md), copied from commit `a3d04d3467361bab7c9fc7c8c6aedf31987ee069`, issue #939.

All selected arrays were frozen before the two holdouts. Each validation used 4,096 Cartesian points and six original times, with separate original spatial, temporal and energy-quadrature refinements. Reported values are worst over the six times at h=0.005 and time step=0.0025.

| Seed | Candidate | Sampled full-vector max | Spatial volume L2 |
|---|---|---:|---:|
| 9216391 | ST061-P | 0.026279255060335002 | 0.034035600656015505 |
| 9216391 | ST063-G1R | 0.024546533709074312 | 0.033767328775136926 |
| 9216391 | ST063-G2R | 0.02438535329193401 | 0.03425749818472021 |
| 9216392 | ST061-P | 0.02566467314396125 | 0.03358854497849373 |
| 9216392 | ST063-G1R | 0.026861392722101148 | 0.033883206585681475 |
| 9216392 | ST063-G2R | 0.02250355909296729 | 0.034206087967209475 |

G2R decreases paired sampled maxima by 7.21%/12.32%, but increases L2 by 0.65%/1.84%. G1R worsens the second maximum and second L2; its independent structure audit also records a radial-pressure sign miss. All six original reports fail both momentum gates. The other original sampled gates passing is not a full-domain structure certificate.

### Geometry, not image similarity

In the declared cylinder r<=0.35, |z|<=0.60, using omega_z squared, centered axial variance, and single-transverse variance r^2/2:

| Time | Parent aspect | G2R aspect |
|---|---:|---:|
| 0.25 | 1.48755304 | 1.59680197 |
| 0.50 | 1.65049839 | 1.73146643 |
| 0.75 | 1.75122740 | 1.81256943 |

At t=0.5, the on-axis angular-velocity factor at z=0.6 changes from 0.02876447 to 0.08350277; at z=0 it decreases from 0.09172868 to 0.08589759. This is redistribution, not stronger rotation everywhere. The low-threshold omega_z>=0.15 near-axis band spans the observation window; its length is observation-limited. Neither field forms the claimed axis-connected strong band at 0.25 or 0.35. The hoped-for 10%–20% aspect improvement is not achieved and the central disk remains.

G2R retains five direction checks on the stated fresh core probes and at least 99.699% signed shear on the fresh midplane probes. Bias remains positive but its minimum decreases. Effective-volume changes are small, not literally zero or certified for every time. Exact windows, tolerances, rejected fits and restart limits remain in the original record.

## Residual-oriented controls: ST061

Source: [unaltered ST061 record](research_snapshots/ST061.md), commit `ad0e6dacf3851a12f4272bb4f6b282cf49506d8e`, issue #900. These are different holdouts from ST063; do not rank candidates by picking a favorable value across experiments.

| Seed | Candidate | Sampled full-vector max | Spatial volume L2 |
|---|---|---:|---:|
| 9206291 | ST060-Q | 0.027559315226219006 | 0.03374187106068435 |
| 9206291 | ST061-D | 0.026788020110924435 | 0.03333465903335943 |
| 9206291 | ST061-P | 0.02580782220423295 | 0.03347432641256566 |
| 9206292 | ST060-Q | 0.027992393733693038 | 0.03457784543286563 |
| 9206292 | ST061-D | 0.027204469154934973 | 0.03426968736612646 |
| 9206292 | ST061-P | 0.02620795221578385 | 0.034440441141806846 |

D is the lower-L2 alternative, P the lower-peak alternative on these paired samples. The finite-budget optimizers were not proved optimal. Unconstrained pressure projection was rejected because it reversed the audited axial pressure directions. This organization update does not re-run those studies or reinterpret internal optimizer failures as convergence.

## Historical compatibility baseline

[ST006](../artifacts/research/ST006/manifest.json) remains the `research_baseline` API default for compatibility. Its original seed 9172801 report records max 0.1082289305, volume L2 0.1075843288 and a failed original divergence maximum gate. It is not the latest scientific result. Its bytes, evidence and API are not changed here.

The old ST006/ST030–ST033 publication discussion is retained at the pre-organization commit and in [the old experiment index](../artifacts/research/experiment_index.json). The newer [catalog](research_catalog.json) is an additional discovery index, not a replacement of historical evidence.

## Acceptance and availability

Spatial volume L2 is sqrt(64*mean(|R|^2)) at each time. It is not RMS, a time average, a color scale or an effective-volume statistic. Finite maxima are not continuous suprema. All source-field identity, PDE-validation and blow-up flags stay false.

ST061/ST063 complete offline bundles contain more than their current GitHub branches. New MATLAB comparison code was exported and checked through Python references but not run natively in ST063. Original ST054 native receipts cannot certify this new interface. See [availability and exact hashes](research_catalog.json) and [visualization](../visualization/README.md).
