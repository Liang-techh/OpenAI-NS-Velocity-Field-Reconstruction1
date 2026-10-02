# Local rooted velocity, pressure and physical vorticity

The current compliant family now supplies the complete local rooted
(ur,utheta,uz,p) and cylindrical/Cartesian physical vorticity. This uses the
same selected Cstar, analytic preheat pressure datum and coupled nonlinear
core as the previous stages. Outputs are source enclosures in signed
logarithmic factors, not independently chosen midpoint components.

The local domain is rho=Lambda R in[0,4.1], xi in[-4,4], with
Z=a+sigma/sqrt(Lambda)*xi and the shared exact H(a)=0. It is a thin axial
chart, not a completed global physical-point field or all annular values.

## Fresh shared-root coupled rows

The new radial rebuild uses t=Z-a and

```
H(a+t)=Hprime(a)*t+(-12a-j)*t^2-4t^3
```

before interval evaluation. Its constant H coefficient and the scaled
amplitude-gradient seed are exactly zero by the source identity. The
root remains enclosed uniformly; no provisional root center is selected
as the exact parameter. Degree24/depth3 coupled rows retain the original
positive implicit swirl-source Cauchy bound, pressure and all nonlinear
radial recovery terms. Old finite rows are not read.

At the root, chi=H^2/(H^2+sigma^2) has axial valuation2. Thus the Bessel
radial tail n>24 has zero axial coefficients for the requested k<=2.
This removes a broad generic model-tail bound; it does not remove the
nonlinear correction. Admitted correction-norm radial tails are retained.
The Psi model is affine in rho and also has no tail beyond degree24.

Normalize Psi coefficients BEFORE recombining the affine U0=4Z+j.
Subtracting this affine value from a rounded Uz would lose the very small
radial contribution. Transport mixed Phi/Psi/mean-Psi jets from the root
with |b xi| times the FULL analytic derivative bound of order (i,k+1).
The explicit affine offset 4b xi is retained. The radial average always
comes from the same coefficients, divided by n+1, and the same analytic
primitive remainder.

## Velocity and the same moment primitive

Let epsilon=1/Lambda, D=1-Z^2, L=1-delta Z^2, M=M_z/R. Recover

```
Q=[2Z Uz-(1-delta)Z M-D M_Z]/L
ur=lambda^-1*sqrt(epsilon*rho/2)*Q
utheta=lambda^-1-delta*sqrt(2epsilon*rho)*F0*Phi
uz=lambda^-1-delta*Uz
```

Q_rho and Q_Z use the same M, M_Z, M_Zrho and M_ZZ. The actual primitive
identities M+rho M_rho=Uz and M_Z+rho M_Zrho=Uz_Z imply

```
div(u)=lambda^-2*[Q+rho Q_rho+
       (D Uz_Z-(1+delta)Z Uz-2Zrho Uz_rho)/L]=0
```

The checker binds coefficientwise radial/axial primitive identities, not
just interval containment. A separately reported raw interval diagnostic
contains zero with maximum normalized width approximately2.64e-201 in
the initial packets. This finite interval width is distinct from exact
structural divergence and is not a full NS residual or a global time bound.

## Physical vorticity with canceled microscopic width

For lambda^beta g(rho,Z), the physical spatial operators are

```
Dr=lambda^(beta-1)*sqrt(2Lambda*rho)*g_rho
Dz=lambda^(beta+delta-1)/L*[beta Zg+Dg_Z-2Zrho g_rho]
```

Set K=Lambda G, b=sigma*sqrt(epsilon). The exact F0 derivative is
F0_Z/F0=-K_xi/b. Before enclosing the radial vorticity, cancel
sqrt(epsilon)/b=1/sigma. The resulting source formulas are

```
omega_r=lambda^-2*F0/L * [
  sqrt(2epsilon*rho)*((2+delta)Z Phi-D Phi_Z+2Zrho Phi_rho)
  +sqrt(2rho)*D*K_xi*Phi/sigma]
omega_z=2lambda^-2-delta*F0*(Phi+rho Phi_rho)
omega_theta=lambda^-2*sqrt(epsilon*rho/2) * [
  lambda^delta*(D Q_Z-2Z(Q+rho Q_rho))/L
  -2lambda^-delta*Psi_rho]
```

The coefficient 2+delta in omega_r includes the radial prefactor derivative.
No K_xi/b or inverse microscopic radius is formed. The normalized Psi_rho
comes directly from the normalized rows. Transverse velocity and vorticity
vanish structurally on the axis; omega_z has its finite nonzero source
limit. Cartesian vectors use the same cylindrical rotation, including
the signed radial/angular contributions.

## Original compatible pressure

Pressure is the SAME analytic datum plus its original positive increment:

```
P(rho,Z)=Pstar^2*p0(Z)+epsilon*F0(Z)^2*V(rho,Z)
V=int_0^rho Phi(s,Z)^2 ds; V_rho=Phi^2; V(0,Z)=0
p_physical=lambda^-2-2delta*P
```

V is computed from the finite Phi convolution with the directed integrated
error rho*(2Pmax+E)*E, where E includes the nonlinear radial tail and full
axial transport. Positive global Phi bounds further enclose the primitive.
The factor epsilon is retained because dR=epsilon drho. Nothing is added
afterward to fit a residual. Base and increment keep separate scale logs,
including -2logCstar and -2delta loglambda.

## Coordinates and logarithmic factors

Radial position uses r=lambda*sqrt(2epsilon*rho). Physical z is represented
by separate anchor and offset terms with the same lambda^(1-delta) factor;
the microscopic offset is not absorbed into a rounded anchor. Each
velocity/pressure/vorticity component is a signed sum of coefficient
enclosures and separate positive-scale logs. Absolute F0 is never
materialized or replaced with numerical zero. The tiny delta/time term
also stays separate from the leading lambda power.

## Reproduction and checks

```python
from lei_ren_part1_paper_compliant_rooted_core_field import CompliantRootedCoreField
field=CompliantRootedCoreField()
point=field.field(xi='.5',rho='2',log_tau='-10',theta='.7')
axis=field.field(xi='.5',rho=0,log_tau='-10')
```

Focused reproduction:

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage rootedfield
```

The pipeline has115 ordered modules. Initial receipts contain18 packets
at xi=-1,0,1, rho=0,2,4.1 and logtau=-1,-10; six are axis packets.
The independent checker verifies exact root/model-tail/primitive/divergence/
curl identities, the first original Phi radial coefficient,90 actual
coefficientwise moment-jet packets, positive pressure and exact source replay.

A separate finite-scale kinematic fixture inverts the original physical z
map and differentiates the Cartesian velocity directly at three points and
three difference resolutions. Richardson curl errors are approximately
1e-18 or smaller; divergence errors are approximately1e-37. The fixture
is smooth and solenoidal and exercises the physical formulas; its
parameters/profiles are not the paper source and it is not an NS test.

## Remaining work

Local fixed-r Phi-weighted axial maxima and fractional levels are now
resolved by ROOTED_SWIRL_MORPHOLOGY_2026_10_02.md. Whole-core radial swirl
growth makes the original annular values necessary for the global radial
peak. Measure full radial morphology, fit multitime field/vorticity
exponents and integrate true particles/material winding. F0 width and rho=4.1
cutoff reference do not establish a whole-vortex aspect ratio. Numerical
point selection outside this root chart, signed annular integrals, actual
implicit five-bump values and a complete global field remain open.

Required-domain energy/support, admissible divergence stress and independently
flat remainder are not closed. Genuine n-dependent coefficient recursion,
mean/oscillatory correction and full Cartesian forced-NS acceptance remain
later stages. Original unlocalized whole-space energy is still infinite.
