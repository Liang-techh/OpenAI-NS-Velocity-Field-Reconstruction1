# Pressure capacity and stress-cone overlap of the adaptive bridge

`adaptive_bridge_pressure_capacity.py` keeps the shared 24-mode bridge
velocity fixed and adds 24 compact scalar pressure modes. Their
Cartesian gradients are computed analytically from the moving
similarity coordinates. The fit minimizes complete momentum on three
axial by five radial points at each of `k=11` and `k=19`, using one
shared coefficient vector. Four pressure coefficients reach the
`±20` bounds.

| Scale | Before extra pressure | After extra pressure | Angular component lower bound |
| ---: | ---: | ---: | ---: |
| Training `k=11` | `390888` | `361740` | `163618` |
| Training `k=19` | `1.5806e9` | `1.4549e9` | `6.8594e8` |

The pressure fit lowers its training maxima by only about `7%–8%`.
At the four off-time holdouts `k=10.5,11.5,18.5,19.5`, it raises the
maxima by roughly `55%–58%`. A scalar pressure gradient has zero
azimuthal component, so the listed angular residual is a hard lower
bound for **any pressure-only change to this fixed velocity**, even
with a larger basis. It remains many orders above `1e-3`.

`adaptive_bridge_cone_overlap.py` then applies the repository's
physical local analogue of the paper's stress-cone test to the
shared 24-mode field. It checks 15 bridge nodes at each scale and
computes the radial stress primitive wherever the necessary local
`lambda_squared` is positive:

| Scale | Positive `lambda_squared` nodes | Strict cone passes | Passes among five largest residuals |
| ---: | ---: | ---: | ---: |
| `k=11` | 8/15 | 0/15 | 0/5 |
| `k=19` | 7/15 | 0/15 | 0/5 |

At the maximum-residual node, `lambda_squared` is about `-1.74e7`
at `k=11` and `-1.47e12` at `k=19`. At other nodes with positive
`lambda_squared`, either the stress target has the wrong sign along
the cone normal or the required cone ratio exceeds one. Thus the
current physical analogue offers no sampled admissible correction
patch covering the bridge defect. This is **not** the paper's
normalized leading-profile cone or a continuum nonexistence proof.

The practical next construction is a coupled mean-profile redesign:
alter meridional shear, swirl slope, and pressure/moment data so that
complete radial/axial momentum improves and a strict cone opens where
the remaining angular defect is concentrated. The radial moments in
the [OpenAI paper's Appendix A and the oscillatory stress mechanism in
Sections 7–9](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
are relevant structural requirements; neither experiment here
implements their full coupled construction. Full-domain momentum max,
spatial-volume L2, finite energy, and the critical-time limit remain
unproved.
