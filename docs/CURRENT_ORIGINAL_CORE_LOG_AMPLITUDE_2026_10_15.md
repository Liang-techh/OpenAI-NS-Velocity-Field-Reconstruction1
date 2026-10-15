# Original logarithmic amplitude and leading-core velocity terms — 2026-10-15

## Delivered

`OriginalCoreLogAmplitude` computes the original anchored primitive G(Z),
keeps F0 in separate logarithmic terms, resolves the narrow amplitude peak,
and combines the numerical leading core with its physical velocity weights.
It returns physical coordinate and Cartesian u,v,w **signed exponential
sum packets**, with both coefficient and logarithmic uncertainties.

The evaluator loads checked source archives directly. It does not reconstruct
historical owners or select midpoints of source boxes. A new native point is
computed on demand by the original numerical recurrence, so the API is not
restricted to the saved demonstration locations.

```python
from openai_ns_reconstruction.original_core_log_amplitude import OriginalCoreLogAmplitude
field = OriginalCoreLogAmplitude(precision=110)
global_packet = field.global_amplitude('0.3')
peak = field.resolved_peak('1')
velocity = field.leading_velocity_terms('0.17','1.5',
                                       log_tau='-2', azimuth='0.25')
velocity['coordinates']  # physical x,y,z,r and t=1-tau, in factored form
velocity['Cartesian']    # u,v,w: lists of signed exponential terms
```

This is a **leading n0 core** evaluator. It does not include n1, global
annuli, heat/stress/oscillatory layers, a full floating-point vector or a
complete NS solution. The full long-term goal remains ACTIVE / INCOMPLETE.

## Original source and global primitive

The exact source uses j=1/5721363714816000000000, sigma=j/500,
logPstar=exp(40)+11, logLambda=4logPstar+1000, and the physical
`delta=exp(-4logPstar-30)`. The universal tube bound delta<=1e-200 is not
selected as a physical parameter.

`U0=4Z+j`, `L=1-delta Z^2`,

`H=(1-delta)Z/2+(1-Z^2)U0`,

`G(Z)=integral_a^Z L H/(H^2+sigma^2) dZ`,

where a is the unique original H(a)=0 root in [-j,0]. In particular, **a is
not zero**: it is approximately -3.88407787546799802629e-23.

Numerical cubic roots from these source formulas give a nominal six-pole
primitive. The directed calculation separately uses the archived uniformly
certified pole disks, exact tuple endpoints and guarded half-plane logarithms:

`G=sum_poles L(p)/(2Hprime(p)) [log(Z-p)-log(a-p)]`.

The exact G(a)=0 branch preserves root correlation; no near-root subtraction
is used to infer it. The source nonnegativity theorem bounds the lower
endpoint of G by zero, rather than replacing G by its absolute value.

Examples with 110-digit nominal arithmetic and directed original errors:

| Z | G approximation, shortened | Error bound for the full saved decimal |
|---|---:|---:|
| 0 | 1.381024466315153723 | 1.74e-110 |
| 0.3 | 12.589977252638731465 | 3.55e-109 |
| exact anchor | 0 | 0 |

The selected logCstar is the canonical **binary singleton** in the current
physical-norm family. Its two exact endpoint tuples are identical; reading
that singleton does not select a representative from a nondegenerate box.

F0 is retained as `exp(-selected_logCstar - Lambda*G)`. The anchor term and
relative term remain separate, including in physical swirl. Neither exp of
the enormous negative anchor log nor a rounded sum with its small correction
is materialized. Error in G is propagated through Lambda. A small absolute
G error alone does not give usable global relative error in F0; the global
relative-amplitude-control flag remains false.

## Resolved original peak

Use the canonical local coordinate

`Z=exact_shared_anchor + b*xi`, `b=sigma/sqrt(Lambda)`, `|xi|<=4`.

The width b is vastly smaller than the error of any ordinary serialized
anchor approximation. The code never adds b*xi to that approximation. It
uses the exact shared-root identity before forming the original rational
integrand:

`H(a+b*s)=b*s*q(b*s)`,

`q(b*s)=Hprime(a)+(-12a-j)*b*s-4*b^2*s^2`,

`Lambda*G=integral_0^xi s L(a+b*s) q(b*s)/(1+epsilon*s^2*q(b*s)^2) ds`.

All five polynomial numerator terms are retained. The directed omitted-
denominator contribution is at most
`epsilon*Lmax*qmax^3*|xi|^4/4`. Positivity of L and q fixes its sign. This
bounds the actual source correction; the limiting Gaussian is not imposed
as a replacement profile.

For xi=1, the numerical relative amplitude is

`F0(a+b)/F0(a) = 0.105399224561864336783217689240698...`,

with error below 4.59e-111 for the full saved decimal. At xi=0 the ratio is
exactly one. This is a local **axis-amplitude ratio**, not the full physical
swirl maximum, a measured vortex width, material winding or temporal recursion.
The absolute amplitude remains in its tagged logarithmic representation.

## Physical coordinates and velocity semantics

The physical API uses native rho,Z, log(tau) and azimuth phi. Its domain is
0<=rho<=4, |Z|<1, finite log(tau)<=0, all finite real angles in radians;
nu=1 and T=1. It returns the corresponding physical position through

`ell=log(lambda)=[log(tau)-log(1-Z^2)]/2`,

`r=lambda*sqrt(2*epsilon*rho)`, `z=lambda^(1-delta)*Z`,

`x=r*cos(phi)`, `y=r*sin(phi)`, `t=1-tau`.

The original leading profile recovery is

`ur=Q*exp[(log(rho)-logLambda-log(2))/2-ell]`,

`utheta=Phi*exp[-selected_logCstar-Lambda*G
               +(log(rho)-logLambda+log(2))/2-(1+delta)*ell]`,

`uz=Uz*exp[-(1+delta)*ell]`.

The Cartesian rotation is u=ur cos(phi)-utheta sin(phi),
v=ur sin(phi)+utheta cos(phi), w=uz. Q is never mistaken for physical ur.
At rho=0 the radial and swirl terms are structurally zero; the axial term
remains. The same original mean recovery supplies Q, preserving the source
streamfunction relation rather than fitting meridional components separately.

Each term contains a numerical coefficient approximation, original-source
coefficient enclosure, and separate logarithmic scale terms. Its meaning is
`coefficient * exp(sum(log_scale_terms))`; Cartesian components are sums of
these terms. The exact anchor-log tuple is a separate negative term. Such
packets retain positive nonaxis swirl and avoid float underflow. They are
not a precomputed ordinary floating-point velocity vector. Automatic inverse
selection from arbitrary physical Cartesian coordinates remains open.

## Checks and recovery

```powershell
python scripts/run_original_core_log_amplitude.py --write-report
python scripts/run_original_core_log_amplitude.py --inspect-report
```

The final construction exits 0 in 7.390 s. Evidence includes 45 comparisons
between 80- and 110-digit calculations, two independent nominal real
integrals in a sinh coordinate resolving the narrow pole scale, 18 matches
to the accepted root-centered peak source, and 25 physical velocity/coordinate
packets. Nineteen nonaxis packets keep positive swirl coefficients and all
three original scale terms. The extra point Z=0.17,rho=1.5 is computed from
the original recurrence on demand, away from the stored demonstration points.

The real-integral calculations are numerical diagnostics, not a new validated
quadrature certificate. Actual error bounds come from the certified poles,
original source enclosures and local exact-root remainder estimate.

A fresh process checks saved report/source/code pins in 0.406 s with exit 0,
without reconstructing an owner. The artifact is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_core_log_amplitude.json`.
The prior normalized-core report and its implementation remain unchanged.

## Next dependencies

Task 15 remains IN PROGRESS. Local relative amplitude and global logarithmic
amplitude data are now available, while full physical evaluation remains open.

1. Feed the tagged amplitude/source scales into numerical recovery of the
   actual regular n1 equation, retaining all feedback, original axial
   derivatives, nonzero histories and independent moment repair.
2. Assemble leading and coefficient terms through their common streamfunction
   and original n-dependent temporal weights. Keep huge logarithmic terms
   separate and report usable component error, not just G error.
3. Close actual higher inner interfaces, restore inverse physical chart
   selection, common-N/error budgets and full original numerical histories.
4. Continue global five-moment/energy/heat/stress layers, genuine n>=2 repair
   and summation, oscillatory cancellation, measured contraction/material
   winding and independent fixed-forcing corrected Cartesian NS validation.
