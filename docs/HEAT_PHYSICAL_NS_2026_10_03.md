# Actual heat exterior in physical coordinates - 2026-10-03

`CompliantHeatPhysicalC4` connects the accepted exact-Gamma leading exterior
to its original physical velocity, pressure and completed background stress.
It consumes the current source dispatcher, physical source assembly and
actual terminal-history stress receipts. Its scope is the heat exterior,
not the collar, inner annuli or the complete corrected field.

## Actual velocity and pressure functions

Let tau=T-t_phys>0, r=sqrt(x^2+y^2)>0 and a=delta/2. For constant physical
viscosity nu>0, define

```
c_infinity = Ev0 * theta_base * Rtail^((1+delta)/2)
A = 2^((1+delta)/2) * c_infinity
A_nu = nu^(1+delta/2) * A
xi = 4*nu*tau/r^2
g(r,tau) = A_nu * r^(-1-delta) * H(xi)
[u,v,w] = [-y/r*g, x/r*g, 0]
p = -A_nu^2 * integral_r^infinity q^(-3-2*delta)*H(4*nu*tau/q^2)^2 dq
```

H is the full original positive Gamma expectation. The amplitude is the
same source-defined Ev0, theta_base and exact positive Rtail used by the
current reconstruction; c_infinity is defined from them. No free amplitude,
pressure constant, numerical cap endpoint or fitted tail is selected.
The API retains the amplitude in correlated finite logarithmic source
parts, rather than materializing extreme positive scales.

The source coordinates satisfy

```
lambda^2*(1-Z^2)=tau
r = sqrt(nu)*lambda*sqrt(2R)
z = sqrt(nu)*lambda^(1-delta)*Z
log(R/Rtail)>=3
```

Finite physical points have |Z|<1. The Z=+/-1 source-cover endpoints are
physical-infinity limits, not finite Cartesian points. The formula above is
independent of physical z within this regional domain; the domain itself
depends on z through the source-coordinate condition.

## Physical momentum and stress

The source bridge explicitly consumes the admitted exact radius/Gamma/
constant identities, then binds the actual packet bracket, theta_base and
absolute pressure scale. It proves K=H, xi_source=2*(1-Z^2)/R,
Utheta=c_infinity*R^(-(1+delta)/2)*H and the matching c_infinity-squared
pressure tail. The same bound amplitude is used in the physical proof.

The original lambda powers cancel in the physical exterior. Consequently:

- The full Gamma ODE gives g_t=nu*(g_rr+g_r/r-g/r^2).
- The full pressure integral gives p_r=g^2/r and p_z=0.
- ur=uz=0, and g_z=g_zz=0, so axial viscosity also vanishes.
- Cartesian moving-basis terms give (u dot grad)u=-(g^2/r)e_r.
- All three local physical momentum components and divergence vanish.

The background stress maps with lambda^(-2-delta), and the general
viscosity pullback multiplies the tensor by nu. The completed
divergence-form background tensor T_B is exactly zero in this exterior.
The physical Cauchy tensor has pressure and the generally nonzero shear
nu*(g_r-g/r); it is a separate object. The API names the zero tensor
`completed_background_stress_tensor_components`.

The leading background remainder vanishes locally, including its radial
component. This conclusion uses the centrifugal balance as well as the
angular/axial equations; it does not discard radial momentum.

## Viscosity and derivative units

The source family is normalized to viscosity1. The actual conversion is
x_source=x_physical/sqrt(nu), u_physical=sqrt(nu)*u_source and
p_physical=nu*p_source. An Nth Cartesian spatial derivative therefore
carries nu^((1-N)/2) for velocity and nu^((2-N)/2) for pressure. Fixed-x time
derivatives retain sqrt(nu) and nu respectively. These factors are explicit
in the returned signed source rows and their logarithmic upper bounds.

## Scope of the evidence

The producer supplies source-bound physical mixed4/time1 enclosures and
exact regional stress, momentum and divergence identities. The independent
checker differentiates Cartesian velocity components and the full pressure
integral under its integral sign at moderate parameters with nu=.01 and .7.
These fixtures check the formulas and units; their small numerical residuals
are not the full project NS residual for the extreme source parameters.

The final producer/checker passes after targeted read-only review. It checks
48 exact regional zero bounds,24 finite zeroth physical source rows,6
Cartesian momentum components,2 divergences,10 viscosity derivative-unit
rows and an actual nonunit-viscosity packet. Moderate momentum residuals
are at most approximately4.93e-47; functional source identities establish
the actual regional result, rather than extrapolating from that fixture.

The new regional API owns
`actual_compliant_field_physical_heat_region_certified=True`.
Older source/stress receipts keep their earlier, narrower physical flags.
All-region stress/cone, independent global flat remainder, physical energy,
finite-width bridge/second switch, global nonlinear point evaluation,
genuine coefficient recursion and oscillatory correction remain unfinished.
There is no certification of the critical-time endpoint tau=0.

## Reproduction and next task

```
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_physical_C4.py
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_physical_C4_check.py
```

The ordered workflow provides `--stage heatphysical` after globalphysical.

Next recover the actual collar Ttheta/Tz from its full sigma/phi/Gamma
moment integrals and original (3.16)-(3.17), then map its stress and axial-
viscosity/radial remainder. Do not impose exterior zero stress inside the
collar. Continue the nonzero-stress cone margins and all remaining regions,
while preserving finite-width annular and coefficient-recursion blockers.
