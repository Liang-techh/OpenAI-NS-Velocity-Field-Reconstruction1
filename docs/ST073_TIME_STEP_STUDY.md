# Patch time-step and Hermite interval screen

The first patch-evolution interval was repeated with `--time-step 2.5e-6`, one quarter of the original `1e-5`. The scripts now accept explicit input/output report names so both experiments are reproducible without overwriting each other:

```powershell
python experiments/root_st073/curl_wave_patch_evolution.py --time-step 2.5e-6 --output-name curl_wave_patch_evolution_dt2p5e6.json
python experiments/root_st073/curl_wave_patch_trajectory.py --evolution-name curl_wave_patch_evolution_dt2p5e6.json --output-name curl_wave_patch_trajectory_dt2p5e6.json
```

The trajectory script also constructs a cubic Hermite potential on the first interval. Its potential value starts at zero, ends at the explicit Euler value, and matches the separately fitted slopes at both endpoints. Pressure is interpolated continuously between the two endpoint fits. Since velocity is the curl of the interpolated potential, it is analytically divergence-free at every interior time. The full nonlinear momentum operator is evaluated directly at the interval midpoint on 16 held-out spatial nodes and eight angles per node. The temporal finite-difference stencil is kept inside the interval.

| First interval step | Next-node projected max | Midpoint baseline max | Midpoint linear-potential max | Midpoint Hermite-potential max | Midpoint Hermite RMS |
| --- | ---: | ---: | ---: | ---: | ---: |
| `1e-5` | `743,360` | `726,486` | `2,457,367` | `2,376,202` | `919,513` |
| `2.5e-6` | `412,793` | `728,208` | `696,227` | `687,383` | `321,218` |

With the shorter step, the Hermite field improves the **sampled midpoint** maximum and RMS over the original field (`687,383` versus `728,208`; `321,218` versus `479,995`). The first-harmonic next-node slope norm is `13,445`, compared with `43,685` for the larger step. This shows that time-step control matters and gives a usable small time element for developing the coupled solver. It is not an accepted pulse: only one interior time and a small spatial sample were checked, the next element and temporal endpoint are not matched, and the residual still exceeds `1e-3` by roughly nine orders of magnitude.

The next solver should adapt time steps based on **interior** full-momentum error, continue pressure and potential smoothly across many elements, and expand the spatial basis or handle the generated higher harmonics and mean moments. This remains a prototype of the amplitude/pressure evolution logic, not the estimates or construction of the [OpenAI Navier--Stokes paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
