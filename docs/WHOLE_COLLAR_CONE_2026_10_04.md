# Whole original collar cone accepted - 2026-10-04

The original collar now has an admissible two-vector stress cone on the
whole source domain `0 <= offset < 3`, `-1 <= Z <= 1`. The focused
`collar_cone_check` passes with current source hashes and seven positive
source margins. It combines the original sigmoid interval [0,1] with the
accepted heat-part proof [1,3). At offset3 the stress is exactly zero, so
strict inequalities exclude that endpoint.

This closes the local collar cone gate. All-region stress, independent
global flat remainder, physical energy, genuine coefficient recursion,
oscillatory correction and the full corrected NS residual remain open.

## Original sigmoid source

Write t=offset, a=delta/2, k=1-a, b=(1-delta)/2,
L=1-delta*Z^2, eps=.001*delta and S=1/Rtail. The same canonical Gamma
reference is H=H_delta(2*(1-Z^2)*S*exp(-t)). It is an analytic reference;
the actual field and original infinite future moments are retained.

```
K=(1-sigma)*(1-eps)+sigma*(1-eps*phi)*H
phi=exp(-4/(3-t)^2)
W=1-sigma+sigma*phi
f=1-sigma+sigma*eps*phi
D=1-H
G=H-K=eps*W-f*D
G_Z=f*H_Z
```

The production flat sigmoid has
sigma_t=sigma*(1-sigma)*(2/t^3+2/(1-t)^3)>=0 on (0,1).
Its source derivative enclosure supplies the finite absolute cap32;
monotonicity is proved from the exact formula, not from the signed hull.
The AST bridge binds its value, first derivative and reflected evaluation,
the original sigma/phi/Gamma composition and flat endpoint jets at0 and1.

The admitted C3 bridge supplies the actual complete angular/energy/pressure
moment identities, canonical Gamma zero stress and exact zero meridional
histories. This companion consumes that current-hash proof rather than
re-proving its full history transfer.

## Uniform angular and shear margins

The accepted canonical Gamma bounds give
Dcap=2*a*(1+a)*Scap, 0<=D<=Dcap and 0<=H_t<=Dcap.
On t in[0,1], phi in[exp(-1),exp(-4/9)], W>=exp(-1),
W_t<=0, f_t<=0 and abs(phi_t)<=1. Thus

```
Gmin=eps/e-Dcap>0
k*G-b*Z*G_Z >= k*eps/e-Dcap*(k+2*b)>0
-G_t+(1+a)*G
 >=eps*(1+a)/e-Dcap*(sigma_t_cap+eps+2+a)>0
hot_waiting_gap >=eps*(1-exp(-4/9))-Dcap>0
```

The last bound and sigma_t>=0 give K_t>=0. Together with the accepted
whole-collar shear proof, 0<kappa-2<=delta. Subtracting the complete
zero-stress Gamma moments from the actual moments gives

```
ctheta=(k*integral_t^3 exp(k*(q-t))*G dq
       -b*Z*integral_t^3 exp(k*(q-t))*G_Z dq+G)/L
       +2*S*exp(-t)*(-G_t+(1+a)*G)
ctheta>=Gmin>0.
```

Future angular positivity on q>=1 comes from the accepted heat proof.
Compact support of these differences before3 does not truncate or reset
the original infinite moments.

## Common squared-density bound and cone

Let Q=K^2-H^2=-G*(K+H), d=1-Z^2, p=1+delta.
For all future q in[t,3], 0<K<=H<=1, 0<=G<=eps and
d*abs(Q_Z)<=2*Dcap. Regroup the common energy/pressure density first:

```
N=integral_t^3 [
 Z*(delta*wE-p*wP)*Q+d/2*(wP-wE)*Q_Z ] dq
wE=exp(-delta*(q-t)), wP=exp(-p*(q-t))
abs(cz)<=Cz_sigma=3*(2*eps*(1+2*delta)+Dcap)/(1-delta).
```

The actual positive B amplitude is bounded using the same source-log
parts as the physical assembly. Its receipt proves

```
log(delta)+2*log(Bmax)+2*log(Cz_sigma/Gmin)<log(2).
```

Hence ctheta>0 and
2*ctheta^2-(kappa-2)*B^2*cz^2>0 throughout the sigmoid interval.
The source is identical at offset1, with the same full future integrals.
Compose this result with the heat proof to cover the entire original
collar. The accepted heat proof also supplies the uniform positive
e_theta direction as offset approaches3 and the paper's O((3-offset)^6)
stress ratio for fixed source parameters.

## Physical API and reproduction

`CertifiedCollarPhysical` validates the accepted receipt, all input
hashes and the same field/source family before calling the unchanged
`CompliantCollarPhysicalC2.collar` evaluator. It adds local proof metadata
and the receipt digest; it does not project stress or alter velocity,
pressure, scale, viscosity or the retained nonzero axial-viscosity remainder.
Fixed positive viscosity and lambda factors preserve this two-vector cone.

The cone concerns (Ttheta,Tz) against (Stheta,Sz). It does not certify a
cone for the completed tensor diagonal T_theta_theta=r*partial_z(Tz).
The primitive physical companion remains distinct from this admitted
proof adapter. Global flatness cannot be inferred from regional momentum
decomposition.

Run `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py
--stage collarcone` to reproduce producer then checker. The seven checked
positive margins include Gmin, future angular, shear, hot/waiting gap,
Lmin, epsilon and delta. No broad pipeline rerun is needed for this
isolated companion.

Next complete stress/cone and flat-remainder bounds in the remaining
regions, finite-width bridge/second switch and common leading inputs;
then solve the genuine coupled n=1 and n>=2 recursion.
