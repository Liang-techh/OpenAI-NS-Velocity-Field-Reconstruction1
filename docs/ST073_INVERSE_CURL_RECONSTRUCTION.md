# Spatial pulse inverse reconstructed as an exact curl: pressure fixed, viscosity remains

`midplane_inverse_curl_reconstruction.py` converts the five-node
midpoint transverse amplitudes from `ST073_SPATIAL_PULSE_INVERSE_TWO_SCALE.md`
to a compact cylindrical vector potential on `k=11` and `k=19`.
It also fits their local physical-time derivatives and a pressure.
The resulting correction velocity is an analytic curl, hence its
divergence vanishes wherever the field is smooth. This experiment
screens the **complete nonlinear Cartesian momentum** at the pulse
center on eight independent radial/axial points × 16 shifted angles;
it does not construct a temporally supported pulse.

| Potential and pressure reconstruction | `k=11` held-out max / frozen | `k=19` held-out max / frozen |
| :--- | ---: | ---: |
| Five polynomial directions per vector component, normal-pressure values from local inverse | `2.77e5` | `2.81e5` |
| All 27 polynomial directions with column-normalized ridge, same normal-pressure values | `5.11e43` | `8.01e45` |
| All 27 directions with physical coefficient penalty, same pressure values | `15.57` | `15.16` |
| Same regularized potential; compact pressure **gradient** fitted to full momentum on 16 separate spatial nodes | `3.663` | `3.666` |

The unregularized full polynomial fit matches the five node velocity
amplitudes to `~0.3–1.4%` relative, yet explodes between nodes. This
is a severe underdetermination/derivative-control failure; its huge
values should not be interpreted as a physical obstruction. Penalizing
the **physical potential coefficients** prevents that blowup but leaves
`~23%` mode-1 and `~51%` mode-4 velocity-amplitude fit errors.

The normal-pressure scalar values fit to about `1.5e-8` relative at
the five nodes, but their interpolated gradient dominates the momentum
peak: `9.17e10` at `k=11`, versus `1.26e9` from viscosity. Replacing
that scalar interpolation by a compact pressure-gradient fit to the
full sampled momentum lowers the held-out maximum from `9.24e10` to
`2.17e10` at `k=11` and from `3.33e14` to `8.05e13` at `k=19`.
At those new peaks the viscosity term dominates: `1.82e10` and
`6.73e13`, respectively; the pressure term is only `1.69e8` and
`6.75e11`. The same factor-of-`3.66` failure on both scales is a
reproducible scale shape, **not** a contracting recursive update.

The paper's [Lemma 7.7 and residual identity (7.40)](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
keep curl and cutoff remainders in the budget. Our five-node principal
inverse omits their spatial derivatives; the direct reconstruction
shows why they cannot be deferred. The pressure-gradient fit here is
an exploratory PDE projection, not the paper's moving-normal pressure
identity. The next construction should solve for a spatially varying
amplitude with derivative bounds and a wider admissible axial support,
then evolve it along pulse paths and include pressure, mean and moment
repairs in a full residual cycle. Additional sparse-node interpolation
or pressure tuning alone will not address the observed viscous peak.
