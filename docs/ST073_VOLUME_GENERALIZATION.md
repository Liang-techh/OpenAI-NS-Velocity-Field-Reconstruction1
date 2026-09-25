# Whole-support check of the meridional collar trial

The time-scaled meridional modes in `ST073_DYNAMIC_POLOIDAL_TRANSITION.md` reduced full momentum at selected collar points, including an independent set. Before promoting them, `dynamic_poloidal_volume_screen.py` evaluated the entire compact-support cylindrical quadrature at registered `tau=0.0084`. This changed the conclusion: the local fit **does not generalize to spatial volume**.

| Whole-support Gauss order | Dense-swirl baseline max | Local-fit poloidal max | Baseline volume L2 | Local-fit volume L2 |
| ---: | ---: | ---: | ---: | ---: |
| 6 | 9.7092e4 | 2.5201e5 | 158.45 | 318.54 |
| 8 | 1.7610e5 | 2.9096e5 | 128.13 | 146.87 |

The new worst points remain in the axial collar. Even one-tenth of the fitted poloidal amplitude worsens both max and L2 on order 6, while order 8 initially improves. This disagreement shows why a few Gaussian nodes cannot stand in for a continuum maximum.

`dynamic_poloidal_volume_refit.py` then refitted the same six modes jointly to order 6 and 8, combining normalized physical-volume L2, point residual, and a fourth-power peak penalty. It improved **both training grids**:

| Order | Baseline max → refit max | Baseline L2 → refit L2 |
| ---: | ---: | ---: |
| 6 | 9.7092e4 → 6.1814e4 | 158.45 → 141.34 |
| 8 | 1.7610e5 → 6.5444e4 | 128.13 → 46.44 |

But the untouched order-10 grid worsened from max `2.6558e5` to `3.7516e5` and L2 `311.28` to `395.85`. Independent relative collar points at `tau=.0084,.012,.024` also worsened: max values `5.9139e5 → 7.9443e5`, `3.1327e5 → 4.3282e5`, and `9.5934e4 → 1.3863e5`. The refit is therefore **rejected** as an accepted candidate despite its training metrics. The saved coefficients are a diagnostic only.

This is evidence against continuing the same low-degree compact polynomial fitting as the main route. The moving narrow axial cutoff creates strong second derivatives; fit coefficients that cancel selected samples move the peak to unsampled positions. The next construction should derive a time-dependent meridional transition from the pressure-free vorticity/streamfunction equation with inner and outer boundary data, then enforce radial moments and exterior heat/stress matching. The [OpenAI paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf) treats these matching and correction stages explicitly; the present local taper has not replaced them. Any proposed transition must pass adaptive spatial maxima, stable volume quadrature across orders, several time scales, and the smooth compact-force requirement before acceptance.

Reproduce the volume and refit studies with `python experiments/root_st073/dynamic_poloidal_volume_screen.py`, `python experiments/root_st073/dynamic_poloidal_volume_refit.py`, and `python experiments/root_st073/dynamic_poloidal_volume_refit_holdout.py`. Their JSON reports are in `experiments/root_st073/compact_potential/` and all retain `accepted: false`.
