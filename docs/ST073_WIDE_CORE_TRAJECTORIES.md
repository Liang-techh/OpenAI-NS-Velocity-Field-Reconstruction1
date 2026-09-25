# Material transport in the wider low-residual core

The wider ST073-V order-14 recurrence is a local low-residual field, not a
global Navier–Stokes candidate. To test whether its swirl scaling also
corresponds to particle motion, `wide_core_trajectories.py` integrates the
physical-time ODE `dx/dtau=-u(x,tau)` from `tau=0.128` toward `0.0084`.
It stops at the first exit from `X<=3/32, |eta|<=0.5`, and accumulates the
unwrapped azimuthal angle with `dtheta/dtau=-u_theta/r`.

Two midplane starting locations, `X=0.02` and `X=0.05`, both encounter the
**axial** boundary before reaching the final time. The last common sampled
time still inside the domain is `tau=0.04609`:

| Initial `X` | Radius at last sample / initial radius | Accumulated turns through last sample | Angular speed, start → last sample | Last sampled `eta` |
| ---: | ---: | ---: | ---: | ---: |
| 0.02 | 0.1415 | 0.6436 | 31.40 → 73.25 | 0.4094 |
| 0.05 | 0.1396 | 0.6374 | 31.12 → 69.19 | 0.4614 |

These are **material-point** measurements through the last saved in-domain
sample, not the exact exit angle or a vorticity-defined vortex-core radius.
The geometric similarity radius would contract only by
`sqrt(0.04609/0.128)≈0.60` over that interval, whereas these particular
particles contract to about 0.14 of their initial radii. The remaining
local interval is inaccessible to these particles because axial transport
and a shrinking axial support carry them to `|eta|=0.5` first.

This is a specific obstruction to interpreting the current local swirl
profile as sustained material winding near critical time. It does not
invalidate the local Eulerian recurrence or prove all nearby trajectories
exit. The next construction must solve axial transport and outer return flow
with the moving support, then repeat material trajectories across the
inner–transition–outer field. The [paper's](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
slow relative axial scale
`ell_z/ell_r ~ tau^(-h)` alone does not guarantee particle residence.

Run `python experiments/root_st073/wide_core_trajectories.py` from the
repository root. Saved paths and solver counts are in
`experiments/root_st073/wide_core_trajectories.json`. Both code and report
retain the local-domain limitation; no finite-energy global field or
full-domain `1e-3` acceptance is claimed.
