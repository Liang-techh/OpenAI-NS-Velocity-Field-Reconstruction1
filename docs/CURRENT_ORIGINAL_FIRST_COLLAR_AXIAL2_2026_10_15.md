# Full original n1 first-half collar axial Z2 and q_Z — 2026-10-15

Subsequent progress: [actual selected n1 cutoff and own lower history](CURRENT_ORIGINAL_N1_CUTOFF_NUMERIC_2026_10_15.md)
are now numerically evaluated through Rin and r_minus. The older API below
is the uncut/unrepaired regular-profile interface. Actual cut/own physical
assembly and later modified-chart continuation remain open.

Subsequent progress: [physical first gradient and vorticity](CURRENT_ORIGINAL_FIRST_COLLAR_GRADIENT_2026_10_15.md)
are now evaluated on the admitted first-half collar with original
coordinate/F0/delta factors and a core-boundary first-derivative seam.
Physical second derivatives and complete Cartesian C2 joins remain open.

## Evaluated interface

`OriginalFirstCollarAxial2.evaluate(Z,s)` extends the complete original
n1 core Z2 data through the actual admitted first-half collar:
`|Z|<1, 0<=s<=1/2`, with `rho=4exp(hb*s)`.

```python
from openai_ns_reconstruction.original_first_collar_axial2 import OriginalFirstCollarAxial2
field = OriginalFirstCollarAxial2()
packet = field.evaluate('.3', '.5')
q_value, q_Z = packet['normalized_original_n1']['q']
q_rho, q_rho_Z = packet['native_first_radial_derivatives']['q_rho']
```

The six natural state fields f,u,k,p,fx,ux and own mean have three native
axial Taylor coefficients. q has two. Entry k is `D_Z^k/k!`; multiply
entry 2 and its absolute error by 2 for an ordinary second derivative.
The natural frame is unchanged:
`f=epsilon*F1/F0, u=Uz1/epsilon, k=K1/epsilon, p=P1/(epsilon*Pstar^2)`.
The accepted older Z1 APIs and their receipts are unchanged.

## Actual finite-width transport and errors

The inlet is the complete-error `OriginalN1CoreAxial2` evaluation at rho=4.
Its leading-source, principal-kernel and all-word feedback errors remain
inside every returned collar state enclosure. The actual coupled leading
source is evaluated through Z5 with nonzero implicit swirl and original
pressure data. Weighted second derivatives and Q recovery consume genuine
higher rows, without padding missing derivatives with zero.

The exact original measure `dx=hb*sqrt(rho)*ds/2` first cancels the hb^-1
forcing sector. The frozen singular measure increments remain

```
p:  -8 Z (Psi_rho/Pstar^2) sigma / L
fx: epsilon^2 * 32 Z^2 Phi_rho sigma / L^2
ux: epsilon*Pstar^2 * 64 Z^2 (Psi_rho/Pstar^2) sigma / L^2
L=1-delta*Z^2
```

Their coefficients now retain Taylor orders 0..2. Centered finite-width
variation separately covers rho, sqrt(rho), the actual/comparison Phi
ratio, comparison radial derivatives and the necessary `(1-hb)` factor.
No original width is selected from its bound or exponentiated into an
ordinary number. Zero principal f/u/k increments still carry their actual
finite-width local-integral errors.

The remaining local ODE integral uses the accepted common-inner full
solution bound. Cauchy coefficient bounds `9M/gap^k`, k=0..3, supply the
W_ZZZ row required when differentiating the regular operator through Z2.
This row is a bound used only in the remainder; no numerical W_ZZZ or
q_ZZ is fabricated. The existing width-log proof scales that remainder
before numerical evaluation.

## Own recovery and native mixed derivatives

Writing m=u+k, recovery uses

`q=[2Z*u-(1+delta)Z*m-(1-Z^2)*m_Z]/L`.

The same complete Z2 states provide q_Z. Returned native first radial
rows are

```
f_rho = fx/(2sqrt(rho))       # through Z2
u_rho = ux/(2sqrt(rho))       # through Z2
mean_u_rho = -k/rho          # through Z2
q_rho = [Z*ux/sqrt(rho)+(1+delta)Z*k/rho+(1-Z^2)*k_Z/rho]/L
                              # through Z1, including q_rho_Z
```

Point geometry uses rho=4 as a numerical center; the truth interval
retains the actual original `4exp(hb*s)` geometry. These are fixed-native-
phase Z derivatives. They are not automatically Cartesian derivatives.

## Executed evidence

Evidence is recorded by `scripts/run_original_first_collar_axial2.py` in
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_collar_axial2.json`.
The accepted run exits 0 in **188.525 seconds**, evaluating ten packets:
Z=-0.3,0,0.3 at s=0,0.25,0.5 and the offgrid point Z=0.17,s=0.37.
It passes 18 independent raw-source/singular-measure identities through
Z2, three own-recovery/mixed identities, 54 exact zero-phase increment
coefficients, 69 complete core Z2 seam comparisons, 126 accepted Z1
interface comparisons and 156 precision60/85 complete-error comparisons.
Four unsupported requests are rejected; twenty distinct original width-log
proofs are recorded. Fresh saved-report inspection exits 0 in 2.184 seconds.

At Z=0.3,s=0.5, q_Z is approximately
**94.46400000000000000000086518**, with absolute error below 2.989e-85;
q_rho_Z is approximately **23.61600000000000000000021629**, with error
below 7.472e-86. The pressure increment's Z2 Taylor coefficient is
**5.92587410088389832297**, with error below 4.556e-20. Its corresponding
ordinary second derivative and error are twice these numbers. These are
normalized sample values, not global physical relative-error claims.

At that point the factored p and ux increment coefficients through Z2
have relative error below 7.688e-21; fx is below 2.805e-43. The
1224938-byte report pins 1481 source inputs and 15 implementation files.

The existing read-only Luna/max worker reviews source orders, the W_ZZZ
Cauchy remainder, width variation, recovery formulas and Taylor factorials.
No material mathematical blocker was found.

## Remaining construction

Tasks 14/15 remain IN PROGRESS. The full goal remains ACTIVE / INCOMPLETE.

1. Consume the new q_Z and native mixed rows in physical derivative
   assembly, with original coordinate/amplitude/time factors and complete
   quantitative core/collar joins. Complete Cartesian C2 gluing remains open.
2. Admit and evaluate the remainder of the first, second and lower collars;
   apply source-owned cutoff and own five-moment repair.
3. Complete global physical/chart assembly, explicit tagged-log readers,
   common global N, numerical five-moment closure, heat exterior and energy.
4. Establish admissible stress and flat remainder, genuine temporal n>=2
   coefficient recursion and repairs/summation, oscillatory correction and
   independent corrected Cartesian Navier–Stokes validation.

q_ZZ requires mean Z3. The existing physical collar velocity API is not
silently upgraded to physical derivatives by this separate interface.
Axial Z2 and inverse-width sectors do not establish temporal n=2 recursion,
a global velocity field or a blow-up theorem.
