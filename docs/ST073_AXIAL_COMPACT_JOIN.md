# A full-space finite-energy ST073 candidate and its axial cost

The earlier ST073-V-wide14 field and its radial heat join were defined only
for `|eta|<=0.5`; the pure-swirl heat exterior continued without axial decay.
`axial_compact_join.py` now supplies a Cartesian `fields(points, tau)` for
every spatial point and `tau in [0.5/64, 0.5]` (`t=1-tau`), with physical
viscosity `nu=0.01`. It uses the previous pressure-corrected ratio-4 radial
join as its base, leaves that field unchanged for `|eta|<=0.3`, then uses a
smooth cutoff to zero by `|eta|=0.45`.

The meridional velocity is localized at the streamfunction level. If
`u_r=-psi_z/r`, `u_z=psi_r/r`, and `chi=chi(eta(z,tau))`, the new field uses
`u_r=chi*u_r_old-chi_z*psi/r`, `u_z=chi*u_z_old`, and
`u_theta=chi*u_theta_old`. Hence its spatial divergence vanishes by the
streamfunction identity; multiplying all velocity components by `chi`
would not have preserved that identity. Pressure is provisionally `chi*p_old`.
The cutoff is flat at its endpoints, the axis is regular, and field/pressure
are identically zero outside the moving axial support.

This field has finite kinetic energy at each registered positive time.
Inside radius 1 its axial support and spatial region are bounded. Outside
radius 1 the meridional flow is zero and the pure heat swirl obeys
`|u_theta| <= sqrt(nu)*c*(r^2/(2nu))^(-1/2-h)` with `h=0.005>0`.
Thus the exterior energy integrand is bounded by a constant times
`r^(-1-4h)`, which is integrable. Explicit upper bounds on the energy
beyond physical radius 1 are `1.291e-4` at `k=0.4` and `2.244e-5` at
`k=5.5`. These are tail bounds, not total-energy normalization.

The unforced full-momentum diagnostic shows why kinematic closure is not
PDE closure. At `eta=0.375`, sampled physical residual norms at an inner,
bridge, and heat-exterior radial point are:

| Scale | Inner | Bridge | Exterior |
| ---: | ---: | ---: | ---: |
| `k=0.4` | 225.27 | 494.34 | 0.961 |
| `k=5.5` | 43663.79 | **94692.88** | 190.76 |

The nonzero finite-difference divergence at the steep axial collar is
truncation error: for one late inner point it decreases
`6.23e-5 → 3.89e-6 → 2.43e-7` under successive spatial-step halvings,
the expected fourth-order factor of about 16. The complete momentum
residual remains huge. The next task is a dynamically matched axial
return-flow and pressure transition, with its actual time derivative and
viscous terms, rather than a stand-alone cutoff. The [OpenAI paper's](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
outer matching, moment restoration, and oscillatory correction remain
relevant after this kinematic step.

Reproduce with `python experiments/root_st073/axial_compact_join.py`.
The detailed output is `experiments/root_st073/axial_compact_join.json`;
the field is callable through
`axial_compact_join.load_compact_candidate().fields(points, tau)`.
It is nonzero, full-space, solenoidal, and finite-energy on its registered
time interval, but **not** `pde_validated` or a critical-time solution.
