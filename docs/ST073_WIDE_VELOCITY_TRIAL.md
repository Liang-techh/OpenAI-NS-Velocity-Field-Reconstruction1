# Endpoint-preserving velocity trial on the wide transition

The pressure-corrected ratio-4 ST073-V bridge still has large angular and
meridional momentum defects. This trial adds eight axisymmetric velocity
modes: four meridional streamfunctions and four swirl profiles, with two
radial and two axial shapes each. The radial factor
`1024 y^5(1-y)^5` preserves the inner-core and outer-heat endpoint jets.
The streamfunction velocity uses the full moving-coordinate `z` derivative,
so its increment is solenoidal. Its axial taper vanishes at `|eta|=0.5`.

The complete nonlinear momentum of the combined field was assembled from
fourth-order Cartesian/time differences of the background and each mode,
then fitted at 15 `k=5.5` bridge points. Direct independent finite-difference
evaluation on disjoint points gave:

| Set | Pressure-only max | With velocity modes max | Pressure-only RMS | With velocity modes RMS |
| --- | ---: | ---: | ---: | ---: |
| Training, `k=5.5` | 5942.90 | 5803.26 | 2316.05 | 1925.12 |
| Spatial holdout, `k=5.5` | 3902.26 | 3853.84 | 1648.34 | **1801.47** |
| Time holdout, `k=5.25` | 3023.96 | 2972.33 | 1272.54 | **1388.00** |

The slight decrease in pointwise maximum is outweighed by higher holdout
RMS. Angular-component maximum nearly doubles on both holdouts: `1519.86`
to `2986.16` in space, `1171.53` to `2301.69` in time. One swirl mode
reaches its coefficient bound. The candidate is **rejected as a robust PDE
improvement**. Its small finite-difference divergence (`<5.3e-8` on the
holdouts) shows this is a momentum-shape failure, not a loss of solenoidality.

The result argues against adding more of these low-degree independently
fitted bubbles at one time. The next velocity correction should be derived
from a coupled meridional vorticity/streamfunction and swirl evolution across
the bridge, with pressure and the radial moments solved consistently; it
must be trained and checked over multiple times and spatial bands. The
paper's [heat exterior and five-moment matching](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
and eventual nonaxisymmetric
stress correction remain separate requirements. This finite-slab trial has
no axial closure or finite global energy.

Run `python experiments/root_st073/wide_velocity_fit.py`; the fitted
coefficients and evidence are in `experiments/root_st073/wide_velocity_fit.json`.
The exact rejected field is callable with
`wide_velocity_fit.load_velocity_candidate().fields(points, tau)` for
comparison and future operator diagnostics, without promoting it as the
project's accepted candidate.
