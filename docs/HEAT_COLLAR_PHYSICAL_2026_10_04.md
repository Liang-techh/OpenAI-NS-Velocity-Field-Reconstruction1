# Actual physical collar stress/remainder decomposition - 2026-10-04

`CompliantCollarPhysicalC2` maps the accepted full collar similarity
stresses into the original physical coordinates and forms the completed
divergence-form tensor. It recovers the leading regional identity
`R_B=-div(T_B)+E_B` from the same angular, energy and pressure moments.
It consumes the accepted main physical velocity/pressure assembly and the
separate `CompliantCollarStressC3` source; the main dispatcher is preserved.

Validation status: the producer's exact identities and source-bound record
are available. The independent Cartesian fixture is still running and its
checker receipt has not been admitted. This stage remains open in the task
list until the current-source receipt passes.

## Physical map and source factors

For constant viscosity nu>0, the original unit-viscosity source is pulled
back with `x_source=x_phys/sqrt(nu)`, `u_phys=sqrt(nu)*u_source` and
`p_phys=nu*p_source`. With tau=T-t>0,

```
lambda^2*(1-Z^2) = tau
r = sqrt(nu)*lambda*sqrt(2R)
z = sqrt(nu)*lambda^(1-delta)*Z
offset = log(R/Rtail) in [0,3]
Ttheta_phys = nu*lambda^(-2-delta)*Ttheta_profile
Tz_phys = nu*lambda^(-2-delta)*Tz_profile
```

Finite physical points have |Z|<1 and r>0. Source endpoints Z=+/-1 are
physical-infinity limits. The API encloses these source limits without
claiming they represent finite points or the critical-time endpoint.

The original positive scale is
`B=Ev0*theta_base*exp(-(1/2+delta/2)*offset)`, with
`Qtheta=sqrt(R/2)*B`, `Qz=sqrt(R/2)*B^2`.
Their logs are combined analytically before interval enclosure. In
particular, the opposite inverse-mu terms in Qtheta cancel exactly.
The completed diagonal uses the separate correlated identity
`sqrt(2R)*Qz=2*Qtheta^2=R*B^2`; summing independently enclosed enormous
radius/amplitude logs would lose this cancellation.

Every physical stress derivative of order N has the viscosity factor
`nu^(1-N/2)`. The divergence and remainder derivatives use
`nu^(1/2-N/2)`. Exact positive scales remain formal source factors;
no radius, amplitude or Ev2 cap endpoint is selected as a field value.

The axial-viscosity operators first combine the exact amplitude rate
`B_y/B=-(1+delta)/2` with the physical coordinate derivatives. Their unit
baseline cancels symbolically before interval evaluation, preserving the
small cutoff/Gamma contributions instead of introducing a rounding-sized
baseline error. The waiting endpoint is an exact zero; the flat Gamma join
also has zero remainder through mixed2.

## Momentum and the completed tensor

The collar is pure axisymmetric swirl: ur=uz=0. The same complete pressure
integral gives `p_r=utheta^2/r`, including the original pressure datum.
Therefore radial momentum and velocity divergence are exactly zero.

The original moment-to-stress equations give

```
Rtheta + nu*partial_zz(utheta) = -(partial_r+2/r)*Ttheta
Rz = partial_z(p) = -(partial_r+1/r)*Tz
E_B = -nu*partial_zz(utheta)*e_theta
```

The angular axial-viscosity remainder is generally nonzero in the collar.
It cannot be discarded as it can in the exact heat exterior. The radial
and axial remainder components vanish in this pure-swirl region.

Use the symmetric tensor with cylindrical entries

```
T_rtheta = Ttheta
T_rz = Tz
T_theta_theta = r*partial_z(Tz)
remaining entries = 0
```

The diagonal cancels the unwanted radial divergence `partial_z(Tz)`.
The implementation derives the complete Cartesian divergence independently,
including the moving cylindrical basis. This tensor is the background
divergence-form tensor; it is not the viscous/pressure Cauchy tensor.

## Available numerical interface and limits

The callable returns signed source enclosures and logarithmic upper bounds:

- Physical cylindrical Ttheta/Tz derivatives through total order 3.
- Cylindrical completed diagonal, stress divergence and angular remainder
  derivatives through total order 2.
- Cartesian tensor, divergence, remainder and momentum-decomposition terms
  at order zero, with the source factors and lambda/viscosity powers explicit.
- The independently accepted physical velocity/pressure spatial4/time1 rows,
  converted to the requested viscosity.

At offset3, the same flat K jets and canonical full-Gamma map make the
physical axial-viscosity remainder zero through mixed2, joining the accepted
zero-stress exterior. The source functional identity holds throughout the
collar and uses the exact original moment recurrences, not sampled residual
values to reset the pressure or moments.

The independent checker differentiates the full velocity/pressure directly
in Cartesian coordinates, solves the implicit physical coordinate map and
uses complete future moment integrals for the stress. Independently obtained
stress first jets differentiate those full integrals under the integral sign
and apply the FTC at their moving lower limits. This avoids repeated nested
high-precision quadrature without using production coefficient helpers.
Those jets form the completed Cartesian tensor; omitted higher tensor jets
do not enter its divergence. Two moderate-parameter fixtures exercise
nu=.01 and .7, six momentum-decomposition components, two radial momenta,
two divergence constraints and nonzero axial-viscosity remainders. These
fixtures check the formulas and units rather than measure the extreme
source family's global NS residual.

The Cartesian fixture independently tests the angular and axial stress
divergence. Its completed radial entry cancels for any supplied Tz_z, so
that fixture is not a standalone check of Tz_z. That derivative uses the
accepted general physical chain rule and differentiation of the same full
moment formula; the direct formula also received a bounded read-only review.

The leading axial-viscosity error has the explicit lambda exponent
`delta-3`; this exponent alone does not make it a flat remainder. Higher
background coefficients remain necessary. No all-region cone, global flat
remainder, finite physical energy, volume L2 norm, genuine recursion or
corrected residual completion is claimed here.

## Reproduction and next work

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage collarphysical
```

Next establish the remaining collar stress-direction cone inequalities and
the directional limits near the zero-stress heat join. Keep the completed
diagonal and axial-viscosity error when assembling other regions. Continue
the finite-width bridge/second-switch closure and leading common-field gates
before coupled n=1, n>=2 recursion, flat remainder and oscillatory correction.
