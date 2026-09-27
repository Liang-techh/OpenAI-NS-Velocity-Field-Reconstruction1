# Broad-shear evolution and pressure matching

The positive-energy broad-shear seed failed the actual wave-region stress
cone even though its outer moments and outer cones passed. The next mean
variable is its own amplitude derivative with respect to k. This changes
the physical-time source by dk/dt = 1/(tau log 2), while preserving the
instantaneous velocity and its previously computed energy matrix.

`broad_shear_dynamic_control.py` adds this variable to the nine pressure
and nine outer swirl-slope controls. Its joint problem includes the four
outer moment equations, 22 outer cone locations, and the five wave-region
locations that previously failed. Full problem data are saved before
optimization, so solver failure is not confused with a physical proof.

The completed 19-control solve is feasible: da/dk=-227.6691864863.
Its independent order-96 moment replay maximum is 1.45055e-5. At radial
order 64, all six outer and all five wave-region stress-cone samples pass.
The wave center changes from T dot N = +6141.77 to -3901.40, repairing
the previously identified sign error. The instantaneous velocity and
therefore its energy-growth matrix remain unchanged. No time interval
has yet been evolved with these controls.

The saved problem has 81 inequalities and four equalities in 19 unknowns,
including per-node source vectors, normalized cone maps and physical
locations for future pressure-source updates. Full sampled momentum
still peaks at 7.46467e9 (sample RMS 1.45526e9 on 305 points). Passing
the sampled moment/cone conditions does not satisfy the PDE residual gate.

## Radial pressure repair and its remaining source

`broad_shear_pressure.py` implements a pressure increment for the broad
swirl delta v on a reference swirl V:

```text
g(r,z,t) = (2 V delta_v + delta_v^2)/r
P_inner(r,z,t) = integral from the inner radius to r of g(s,z,t) ds.
```

The velocity is unchanged. This restores the changed radial centrifugal
balance but also changes the axial pressure gradient. The direct trial
uses the previous 18-control mean, not a subsequently evolved field.

| Diagnostic point | Original full momentum norm | With inner-datum pressure |
| --- | ---: | ---: |
| Wave center | 5.81938e7 | 3.15511e7 |
| Previously sampled centrifugal peak | 7.47725e9 | 1.54854e9 |
| Outer radius at eta=0.2 | 75.24 | 4.65556e8 |

At the peak, the radial component falls from -7.31514e9 to 4.27287e4;
the reference without broad swirl has radial component 4.27371e4.
The finite-difference/quadrature difference of about 8.43 is not within
the target tolerance. More fundamentally, axial imbalance is transferred
to the exterior, so the direct primitive is not a matched global repair.

[Section 8, equations (8.3) through (8.7), of the reference paper](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
keeps the radial pressure source together with axial balance and makes a
radial primitive compact by subtracting its total integral with a cutoff.
The resulting cutoff derivative remains in the momentum equation.

The implemented `datum='compact'` option follows that algebra at this
axisymmetric pressure level: P_compact = P_inner - chi*Jg, with chi
transitioning from zero to one on bridge fractions [0.72,0.93]. It is
zero at both radial ends and restores the sampled outer residual to
75.24. However the new cutoff-midpoint residual is 4.71343e9. This is
explicitly retained; it is not discarded or renamed as acceptable forcing.
Both pressure variants remain rejected as complete repairs.

This cannot be fixed merely by tuning the pressure cutoff. For the pure
azimuthal velocity increment, keeping the pressure increment zero at both
radial boundaries gives integral(delta R_r) = -Jg. At eta=0.2 the computed
Jg is 1.96970e6 across an annulus of width 0.00422820. The resulting
lower-bound estimate for the supremum of the ADDED radial imbalance is
4.65847e8. This uses a numerical integral and is not a certified continuum
bound on a future corrected field. It identifies the missing mechanism:
radial/meridional dynamics or actual wave stress must share this balance;
azimuthal time slopes and a compact pressure increment alone cannot remove
the added radial integral.

The next coupled solve must include the induced axial pressure source
and the compact radial remainder. The instantaneous pressure correction
depends on the fixed velocity values, so changing the time-slope controls
does not change this pressure increment at the reference time. That allows
source-vector updates of the saved linear control problem rather than
repeating every velocity/gradient evaluation. The radial budget must be
closed explicitly by a meridional correction and/or actual wave stress.

These are sampled local constructions. Full momentum max and spatial
volume L2 below 1e-3, finite energy, prescribed domain/forcing, actual wave
stress realization, and contraction across scales remain unestablished.
