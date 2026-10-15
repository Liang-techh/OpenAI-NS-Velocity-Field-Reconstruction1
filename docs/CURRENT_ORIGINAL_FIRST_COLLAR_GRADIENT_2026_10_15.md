# Original physical first gradient and vorticity, first-half collar — 2026-10-15

## Evaluated interface

`OriginalFirstCollarGradient.evaluate(Z,s,...)` now computes the physical
first spatial derivatives and vorticity of the original leading+n1 field
on the actual admitted first-half collar, `|Z|<1, 0<=s<=1/2`. The companion
`core_endpoint(Z,...)` uses original core radial rows at rho=4 directly.
Both use the same original amplitude, time and streamfunction recovery.

```python
from openai_ns_reconstruction.original_first_collar_gradient import OriginalFirstCollarGradient
field = OriginalFirstCollarGradient()
packet = field.evaluate('.3', '.5', log_tau='-2', azimuth='.25')
J = packet['gradient_Cartesian']          # 3 by 3 lists of physical term references
omega = packet['vorticity_Cartesian']    # x, y, z components
comparison = packet['axial_vorticity_order_comparison']
```

Output uses shared `log_atoms` and a `term_pool`. A reference contributes

`factor * term_pool[term].coefficient * exp(sum multiplier*log_atom)`.

Every coefficient, multiplier factor and numeric log atom has an original
source enclosure. The F0 anchor is an exact negative binary tuple; the
positive original hb geometry remains a tagged source. Missing numeric
approximations must not be interpreted as zero. Absolute F0, physical
velocity/gradient exponentials and an ordinary float vector are not
materialized. The existing physical value API retains its interface.

## Physical-coordinate differentiation

At fixed physical time, set

```
d=1-Z^2, L=1-delta*Z^2
lambda=sqrt(tau/d)
r=lambda*sqrt(2*epsilon*rho), z=Z*lambda^(1-delta)
partial_r=(2rho/r)*partial_rho
partial_z=(d*partial_Z-2rho*Z*partial_rho)/(lambda^(1-delta)*L)
```

The source rows are differentiated before physical assembly. In particular
q_Z consumes the full natural state Z2 data delivered in the preceding
milestone. Leading Qpsi_Z uses the same actual own mean through Z2.
The first radial rows use the actual point-phase chi and comparison
radial functions. Differentiating a stored per-hb increment log would miss
the cancellation between ds/drho and hb; the API uses the original native
radial identities instead.

For a physical term `E*c`, where
`E=F0^eta * lambda^(A+B*delta) * Lambda^M * constant*rho^p`,

```
partial_r(E*c) = E/r * (2rho*c_rho+2p*c)
partial_z(E*c) = E/(lambda^(1-delta)*L) *
    [d*c_Z-2rho*Z*c_rho+(A-2p)Z*c+B*delta*Z*c-eta*d*g*c/epsilon]
```

The regular, delta-weight and F0-gauge terms remain separate. The
`1/epsilon=Lambda` gauge factor becomes an extra `log_Lambda` atom; it is
not folded into an enormous numerical coefficient. The delta derivative
has its own `log_delta` atom. Original delta*time logarithms also remain
separate. This prevents the regular derivative from disappearing in a
sum dominated by the original F0 derivative.

The exact L factor is positive. Its logarithm uses the original tiny
a=delta*Z^2 and the proven interval `-2a<=log(1-a)<=-a` for a<1/2.
The actual collar rho is retained in coefficient intervals, with its
physical radial geometry separately tagged by the original hb*s/2.

## Jacobian, vorticity and structural divergence

The orthonormal cylindrical gradient is

```
[ ur_r      -utheta/r    ur_z     ]
[ utheta_r   ur/r        utheta_z ]
[ uz_r       0           uz_z     ]
```

Cartesian assembly uses `R*J*R^T`, including the cylindrical connection
terms. Vorticity is

`(-utheta_z, ur_z-uz_r, utheta_r+utheta/r)`,

then rotated by the same R. The native leading and n1 recoveries still
come from their common streamfunction. The report separately stores the
three normalized meridional divergence coefficients and checks that
their full source enclosures contain zero. These coefficients have their
own physical group scales; they are not ordinary physical divergence
values or a useful combined relative-error claim.

The axial-vorticity order comparison cancels original F0 and radius
exactly:

`log|omega_z,n1/omega_z,leading| = logLambda+2delta*log(lambda)`
`    +log| (rho*f_rho+f)/(rho*Phi_rho+Phi) |`.

It is reported only when both complete coefficient intervals resolve
sign. This compares the two admitted temporal orders; it does not admit
n>=2, the infinite sum, or a measured blow-up law.

## Executed evidence

Driver: `scripts/run_original_first_collar_gradient.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_first_collar_gradient.json`.

An independent finite manufactured field differentiates the physical map
and all three physical velocity components directly, then compares the
production cylindrical/Cartesian gradient and vorticity. This test uses
an explicit Gaussian F0 and finite source scales; it is an algebraic
benchmark, not replacement data for the original construction.

The accepted run exits 0 in **197.753 seconds**. Ten original native
packets cover Z=-0.3,0,0.3 at s=0,0.25,0.5 and Z=0.17,s=0.37. Three
additional packets use source-scaled time S=1e18. It passes:

- 3 independent similarity-chain identities and 27 direct manufactured
  physical velocity/gradient/vorticity comparisons;
- 30 normalized divergence enclosures containing zero;
- 81 accepted full-n1 value/Z/rho overlaps and 18 accepted leading first
  radial overlaps;
- 111 core-boundary physical term coefficient seam overlaps, with equal
  gradient/vorticity reference structure on the two source paths;
- 144 precision60/85 coefficient/log error comparisons and six rejected
  unsupported requests.

Each packet has 37 separately scaled physical terms and 12 shared log
atoms. At Z=0.3,s=0.5, the leading swirl F0-gauge derivative coefficient
is approximately -0.20733638330216882250, with error below 3.942e-45.
This coefficient still multiplies its original physical exponential scale;
it is not an ordinary physical derivative value.

At logtau=-2, sampled omega_z,n1 is larger than omega_z,leading. At the
tested S=1e18, all three Z samples instead have n1 smaller, with log ratios
below -1.058e18. F0 and radial geometry cancel from these comparisons.
These are first-two-order comparisons, not an infinite-sum admission or
an independently measured growth law.

The 4636564-byte saved report pins 1484 source inputs and 23 implementation
files. Fresh saved-report inspection exits 0 in 2.177 seconds.

The existing read-only Luna/max worker independently checks the coordinate
chain, source recoveries, original F0/delta factors, cylindrical connection
terms, vorticity signs and rotations. No material algebraic or scale-factor
error was found.

## Remaining construction

Tasks 14/15 remain IN PROGRESS. Full goal ACTIVE / INCOMPLETE.

1. Admit/evaluate the rest of the original first collar, then the second
   and lower collars. Leading coverage through s=1 alone does not admit
   the n1 solution there. This coverage is the next priority. The existing
   read-only worker traces the concrete admission changes below.

   - `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.py`
     limits collar source bounds to s<=1/2. Its `collar_radial_bounds`
     uses a half-interval sigma hull; integrated forcing uses
     `row.shift(1)*sqrt(rho_max)/4`. A full phase requires a fresh [0,1]
     sigma hull and `sqrt(rho_max)/2`, plus an explicit geometry check
     `4exp(cap)<4.1` instead of only `4exp(cap/2)<4.1`.
   - `..._current_original_n1_weighted_feedback_bound.py` likewise uses
     `algebra.width(v)*a/4`; extend the source measure to `a/2` in a
     separate same-family full-phase receipt. Propagate its enlarged
     M0/C0/C1 through natural-component and regular-series/all-word bounds.
   - The current |rho|<5 geometry and rho_max=4.1 appear to leave room,
     but the complete original complex-floor, contraction and Cauchy-gap
     guards must pass with the regenerated source bounds. Keep accepted
     half-phase receipts unchanged and use a new scope-specific receipt;
     the original function/parameter family must remain the same.
   - Only after that admission may a separate full-first-collar numerical
     transport/physical interface accept s<=1. Relaxing the existing guard
     without new source and feedback bounds is unsupported.
2. Apply source-owned cutoff and own five-moment repair, complete global
   physical/chart assembly, and establish quantitative mixed bounds and
   remaining joins. Physical Hessians and a complete Cartesian C2 join
   remain open; q_ZZ would require mean Z3.
3. Complete common global N, numerical five-moment closure, compatible
   analytic preheat pressure/exact heat exterior and finite energy.
4. Establish admissible stress and flat remainder, genuine temporal n>=2
   coefficients/repairs/summation, oscillatory cancellation and independent
   corrected Cartesian Navier–Stokes validation.

This interface covers nonzero radius on this one collar and the rho=4 core
boundary. It does not establish axis-gradient limits, full core gradients,
remaining collars, a global physical field or a blow-up theorem.
