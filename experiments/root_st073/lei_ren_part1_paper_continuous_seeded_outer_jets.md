# Continuous seeded exterior Z-jet adapter

This bounded prototype installs a continuous axial pulse and end-bump value
jet into the actual Rh-seeded outer profile at `Z = .3`. It subclasses the
existing continuous seeded field and keeps the source schedule, angular
correction, pressure adapter, and tail object by identity. The incoming rows
and energy target come from the serialized actual seeded candidate; the
coefficient tangent is obtained by replaying `actual_seeded_input_tangents`
through `ContinuousAxialAlgebra`.

For the pulse region, with `E(Z)` the retained angular amplitude and `g_p(xi)`
the continuous pulse value,

```text
U_z = E a g_p,
U_{z,Z} = E (a_Z g_p + (d log E/dZ) a g_p).
```

For the translated end bumps, the adapter uses the shared continuous bump
basis and applies

```text
U_z = E sum(c_i beta_i),
U_{z,Z} = E sum((c_{i,Z} + (d log E/dZ)c_i) beta_i).
```

The cumulative normalized row is evaluated as

```text
N = base + a I_pulse + sum(c_i I_i),
N_Z = base_Z + a_Z I_pulse + sum(c_{i,Z} I_i),
```

then multiplied by the retained `E exp(-lambda end)` before `R_v` and by
`E_v exp(-end)` after `R_v`. The mean is never assigned zero. At the pulse
sample the radial transport expression uses this analytic axial `U_{z,Z}`;
the radial and log-radius pieces are independently differentiated with two
centered fourth-order stencils.

The receipt was run with the actual source precision `443`, seeded incoming
receipt precision `260`, and continuous algebra / jet precision `200`:

- schedule, angular, tail, and pressure object identities all remain true;
- continuous pulse and mean Z jets are installed, while the angular Z jet is
  explicitly incomplete;
- continuous coefficient tangent row replay is `8.91e-202` and
  `6.66e-201`; energy tangent replay is `0.0` at the stored algebra precision;
- mapped pulse divergence relative cancellation refines from `6.24e-24` at
  `h = 1e-5` to `3.90e-25` at `h = 5e-6`;
- pulse-atom omitted absolute bounds and their signed kernel/mean-jet signs
  are carried into the mean receipt. These are conditional on fixed nominal
  coefficients and quadrature atoms, so they do not enclose incoming,
  coefficient, or quadrature uncertainty;
- the future angular energy derivative remains a fourth-order float-backed
  finite difference supplied by the input tangent provider.

A focused pre-`R_p` smoke evaluation also reaches the explicit inherited
incoming fallback without recursive dispatch; its method is recorded as
`inherited_incoming_float_Z_fallback_before_Rp`.

Run from the worktree with:

```text
python experiments/root_st073/lei_ren_part1_paper_continuous_seeded_outer_jets.py
```

The JSON receipt records signed-log values and the arbitrary-exponent mantissa
for the five sampled radii (`R_p`, pulse, first end bump, `R_v`, and after
`R_v`). It marks full five moments unavailable, keeps quadrature and incoming
uncertainties uncertified, and makes no finite-energy, global analytic Z-jet,
or scale-recursion claim.
