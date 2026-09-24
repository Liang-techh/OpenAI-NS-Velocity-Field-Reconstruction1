# ST073 two-stage patch evolution

Run `python experiments/root_st073/curl_wave_patch_evolution.py` from the repository root. The output, including the fitted potential-slope and pressure coefficients, is `experiments/root_st073/compact_potential/curl_wave_patch_evolution.json`.

This is the first time-marching version of the compact patch prototype. At each of two time nodes, it evaluates the **full nonlinear momentum residual** of the current frozen potential field. It then fits the time derivative of the two exact-curl nonaxisymmetric vector potentials, an axisymmetric exact-curl mean potential, and their pressures across 25 spatial nodes. An explicit Euler update advances the potentials by `1e-5` in `tau`. The fit is evaluated on 16 interlaced spatial nodes not used for coefficients. The wave amplitude is `0.005`.

| Time node | Held-out frozen max | Held-out max after instantaneous derivative/pressure projection | Harmonic slope norms |
| --- | ---: | ---: | ---: |
| `0.0078125` | `1,071,307` | `345,468` | `4,015`, `2,011` |
| `0.0078225` | `5,705,528` | `743,360` | `43,685`, `19,145` |

The second projection still cuts the local sampled residual sharply, but the slope needed for the first harmonic grows by a factor of about `10.9` after one Euler step. The held-out projected maximum rebounds to roughly the original field's scale. This shows strong coefficient variation on the tested step and explains why a fixed endpoint slope is inadequate; a step-size study would be needed to establish numerical stiffness. It does **not** show a stable evolved pulse. The run has no continuous-in-time pressure interpolation or direct residual certificate between nodes, no small temporal endpoint, and no whole-support or physical-volume L2 estimate. The candidate is rejected.

The next numerical formulation should solve the coupled time-dependent amplitude/pressure system with an implicit or adaptive integrator, include the generated higher harmonics and mean-flow moment constraints, and validate the resulting callable field between nodes. This follows the role of the amplitude-pressure solve, exact-curl correction, and mean-flow correction in Sections 7--9 of the [OpenAI Navier--Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf). The paper's smallness estimates have not been established for this ST073 field.
