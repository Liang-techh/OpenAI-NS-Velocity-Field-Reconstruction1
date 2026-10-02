# Actual local physical swirl peaks and axial levels

The same accepted nonlinear core now supplies three fixed-physical-radius
axial swirl maxima and nine true-peak-normalized fractional widths. These
are local rooted-chart source enclosures. Whole-vortex radial width, aspect
ratio, measured blow-up dynamics, energy, stress and temporal recursion
remain incomplete.

## Fixed physical radius, actual Phi and geometry

At each requested time, select an anchor rho_a=1,2,4 and hold
r=lambda(a)*sqrt(2*rho_a/Lambda) fixed during the axial scan. This differs
from holding similarity rho fixed; r is also not fixed across times.

With b=sigma/sqrt(Lambda), Z=a+b*xi, D=1-Z^2 and D_a=1-a^2,

```
rho(Z)=rho_a*D/D_a
u_theta=r*lambda^(-2-delta)*F0*Phi(rho(Z),Z)
T=partial_Z-2Zrho/D*partial_rho
E=Phi_Z-2Zrho*Phi_rho/D
B=E/Phi-(2+delta)Z/D
partial_xi log(u_theta)=-K_xi+b*B
```

The same coupled nonlinear Phi and full mixed derivative tails/transport
enter B and T(B). No independent fitted profile or pressure is introduced.
Original delta, root uncertainty and positive microscopic widths remain.

```
T^2 Phi=Phi_ZZ-4Zrho*Phi_rhoZ/D
        +4Z^2rho^2*Phi_rhorho/D^2-2rho*Phi_rho/D
T(B)=T^2 Phi/Phi-(E/Phi)^2-(2+delta)(1+Z^2)/D^2
kappa=K_xixi-b^2*T(B)>0
```

Strict log concavity on |xi|<=4 and positive B give one local maximum to
the right of the F0 anchor. The peak offset is represented separately as
xi_peak=b*y and Z_peak-a=b^2*y. It is never added to a rounded anchor.
Uniform B/K_xixi bounds and strict endpoint slope signs bracket y.

## Actual peak height and fractional levels

Let B0=B at the anchor, and H=log[u_theta(xi)/u_theta(0)]. Then
H(0)=0, H'(0)=b*B0 and -kappa_hi<=H''<=-kappa_lo. Integration gives

```
b^2*B0_lo^2/(2*kappa_hi) <= max H
                               <= b^2*B0_hi^2/(2*kappa_lo).
```

The lower comparison trial xi=b*B0_lo/kappa_hi is checked to lie inside
the chart. These are bounds using the actual log curvature and anchor
slope, including Phi and lambda derivatives. Positive gain is retained
as separate b,b,gain factors.

The half-maximum, 1/e and 1e-6 levels at each of the three radii are
normalized to this actual local physical peak. Directed original K
remainder, Phi/geometry log-ratio error and peak gain enclose both level
roots, with strict inner/outer signs and uniqueness from log concavity.
Physical axial width uses the original z Jacobian on the whole root
cover. Leading logtau/2 and -delta*logtau/2 factors stay separate.

## Radial result and the next dependency

On |xi|<=4 and 0<=rho<=4.1 the same actual source gives

```
partial_r u_theta=F0*lambda^(-2-delta)*(Phi+2rho*Phi_rho)>0.
```

Thus there is no interior radial swirl maximum in this admitted core.
The core cutoff cannot be a measured radial peak or width. The original
connecting annuli must be resolved to find the global radial swirl peak.
This statement alone does not locate a vorticity core or its width.

Next implement actual coefficientwise core atoms H,M,K,A,B,C from fresh
coupled rows and controlled tails. Feed those into signed prescribed-shear
and first-switch integrals with the same joint histories. Existing bridge
interval covers are not point moments. Preserve the formal comparison
namespace and V100/V110/E identities; caps can bound errors but cannot
define source values. Then recover implicit five-bump values, compose
physical point selection and measure the full radial morphology.

## Reproduction and evidence

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage swirlmorphology
```

The ordered pipeline has117 modules. The focused generator and checker
pass for three source maxima, nine fractional widths and whole-rooted-core
radial monotonicity. Exact fixed-r slope/curvature and second-centered
chart identities are checked independently. A finite-scale synthetic
Phi fixture independently differentiates a scalar physical log velocity
and checks peak location, gain and endpoint signs. It exercises formulas;
it is not the paper source or an NS validation.

Source receipts remain bound to the accepted five-defect family, implicit
source, analytic datum and current file hashes. This local progress does
not close annular values, global fields, required-domain finite energy,
admissible stress/flat remainder, true n-dependent coefficient recursion,
mean/oscillatory correction or full forced-NS residual. Original
unlocalized whole-space kinetic energy remains infinite.
