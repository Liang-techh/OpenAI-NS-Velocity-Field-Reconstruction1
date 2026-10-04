# Actual collar moment-to-stress recovery - 2026-10-04

`CompliantCollarStressC3` recovers both original similarity stress components
through the full heat collar. It uses the admitted velocity, angular repair,
selected energy history and absolute pressure, including the full future
Gamma integrals. It preserves the original datum and all forward histories.
The existing main source dispatcher remains on its accepted pressure provider;
this stress companion consumes that same source without changing its fields.

## Coordinates, units and original moments

Here `offset=log(R/Rtail)` lies in `[0,3]`, `Z` lies in `[-1,1]`,
`a=delta/2`, `k=1-a`, `b=(1-delta)/2`, `d=1-Z^2`, and `L=1-delta*Z^2`.
Axial derivatives are taken at fixed similarity radius R. They are not
fixed-physical-radius derivatives.

```
B = Ev0*theta_base*exp(-(1/2+a)*offset)
Utheta = B*K_full_collar
R = Rtail*exp(offset)
Mz = Mtheta_z = Ur = Uz = 0
Mtheta = sqrt(2)*R^(3/2)*B*A
Mztheta = R*B^2*E/2
P = -B^2*Bp
```

The moment functions have the exact common source definitions

```
A = 1/k + integral_0^infinity exp(k*v)*(1-K(offset+v,Z)) dv
E = integral_0^infinity exp(-delta*v)*K(offset+v,Z)^2 dv
Bp = integral_0^infinity exp(-(1+delta)*v)*K(offset+v,Z)^2/2 dv
```

For offsets beyond 3 the integrands use the complete canonical Gamma heat
function. No finite radial cutoff or finite inverse-radius series defines K.
The original angular-history defect is zero by the admitted waiting/repair
equations; the same full integral transfers this identity throughout the
collar. The selected energy equation and its half normalization, as well as
the cumulative zero meridional histories, are retained.

The current E and Bp multiply the production reference-unit tails by
`exp(delta*offset)` and `exp((1+delta)*offset)`, respectively. Their units
cannot be exchanged. The actual pressure amplitude uses exact
`Ev0^2*theta_base^2`; the runtime Ev2 interval encloses the defining amplitude
and is never selected as its value.
The bridge binds the production inlet U and exact logarithmic scale,
then uses `Ev0=Pstar*U*exp(-13/(2*mu)-13)` explicitly in the unit identities.

## Signed stress functions

Let `A=1/k+Adef`, `E=1/delta+Edef`, `Bp=1/(2*(1+delta))+Pdef`.
The production source has `K=1-epsilon*W-a*S*D`, with exact `S=1/Rtail`;
the numerical S box is only its enclosure. Define positive formal factors
`Qtheta=sqrt(R/2)*B` and `Qz=sqrt(R/2)*B^2`. Original (3.16)-(3.18) give

```
Ttheta/Qtheta = (k*Adef-b*Z*Adef_Z+epsilon*W+a*S*D)/L
               +2*S*exp(-offset)*(K_y-(1+a)*K)
Tz/Qz = (delta*Z*Edef-d*Edef_Z/2
         -2*(1+delta)*Z*Pdef+d*Pdef_Z)/L
```

The code cancels the exact unit baselines before interval evaluation.
This matters because the actual epsilon and delta have exponents of order
minus 4.09e17. Subtracting independently enclosed order-one terms would
bury the intended nonzero stresses under rounding widths.

Defect rows satisfy the same full-function FTC equations

```
Adef_y = Kdef-k*Adef
Edef_y = delta*Edef-(K^2-1)
Pdef_y = (1+delta)*Pdef-(K^2-1)/2
K^2-1 = 2*Kdef+Kdef^2
```

The API supplies every signed stress coefficient with total logR/Z
derivative order at most 3. The positive factors' logarithmic radial rates
are `-a` and `-1/2-delta`; these rates are included in the derivative rows.
Coefficients remain factored, so their small values are not physical
stress magnitudes or norms.

At offset3 the same full future moments, original flat sigma/phi jets and
accepted Gamma moment theorem give exact zero stress and all mixed3
derivatives. This is a functional join; overlapping interval boxes are not
used to establish it. Inward collar stresses are generally nonzero.
The endpoint proof replays the actual inertial and shear code expressions
on the admitted Gamma moments before the interval zero override; the
checker recomputes this binding separately from the serialized zero rows.

## Shear strength and evidence

The same original angular shear gives
`kappa-2=delta-2*K_y/K`. Its full `Z[-1,1] x offset[0,3]` enclosure is
strictly positive. Thus the collar has positive swirl, negative angular
shear and the required shear-strength condition `kappa>2`. Stress-direction
cone inequalities and directional limits near the zero-stress join are
separate open requirements.

The focused producer/checker covers 160 actual finite signed stress bounds
and 40 exact endpoint zero derivatives. Independent moderate fixtures use
full future angular/energy/pressure integrals, 18 axial moment-derivative
checks, 16 original stress/axial/radial derivative checks and four nonzero
stress components. Those fixtures test formulas and normalizations; they
are not global NS residual measurements for the extreme source family.

## Reproduction and next work

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage collarstress
```

The ordered full reproduction runs this stage after the heat companions.
Next map these stress components with `lambda^(-2-delta)` and the constant
viscosity pullback. Construct the completed tensor, retaining its
`theta-theta=r*partial_z(Tz)` entry. For pure swirl the radial momentum
vanishes by the common pressure FTC, while the angular axial-viscosity
remainder must be evaluated rather than discarded. Independently verify
`R_B=-div(T_B)+E_B` using the physical velocity/pressure derivatives.
Then establish collar cone direction/margins and join limits, followed by
remaining regions and flat-remainder bounds/norms.

All-region admissible stress, physical energy, finite-width bridge/second
switch, nonlinear global point selection, coupled n=1 and higher recursion,
oscillatory correction and full corrected NS residual remain open.
