# Actual selected cut/own physical velocity, first gradient and vorticity — 2026-10-15

Subsequent progress: [actual first-micro modified velocity enclosures](CURRENT_ORIGINAL_FIRST_MICRO_NUMERIC_2026_10_15.md)
continue physical value coverage through native s=1 with the original selected N
and correction errors. This document's physical gradient/vorticity scope still
ends at r_minus; it is not extended by that numerical value adapter.

## Interface and exact selected domain

`OriginalSelectedGradient` assembles the original leading plus unrepaired n=1
field on the actual selected inner cutoff and Rin-to-r_minus lower collar:

```python
from openai_ns_reconstruction.original_selected_gradient import OriginalSelectedGradient
field = OriginalSelectedGradient()
inner = field.evaluate('.3', '.625', log_tau='-2', azimuth='.25')
exit = field.evaluate('.3', '1', log_tau='-2', azimuth='.25')
lower = field.lower_collar('.3', '2', log_tau='-2', azimuth='.25')
later = field.lower_collar('.3', '2', log_tau=None, scaled_time='1e18')
```

The inner fraction is q in [0,1]; the lower fraction lies in [1,2]. In both,
actual phase is fraction*original_sc/4. The exact original sc singleton and
phase_in=sc/4 tuple are retained. These coordinates are not the ordinary raw
uncut first phase. The upper endpoint is r_minus, where actual leading
modifications may begin; the unmodified bridge is not used beyond it.

The physical map remains

```
lambda=sqrt((T-t)/(1-Z^2))
r=lambda*sqrt(2*epsilon*rho)
z=Z*lambda^(1-delta)
rho=4*exp(original_hb*original_phase_in*fraction)
x=r*cos(theta), y=r*sin(theta)
```

The API requires native coordinates. It does not provide a new inverse from
arbitrary Cartesian points to selected charts. The output includes cylindrical
and Cartesian velocity, the first spatial gradient, and vorticity. It is a
finite two-order local field, rather than the complete temporal coefficient sum.

## Sources, own histories and physical scales

The accepted `OriginalN1CutoffNumeric` supplies the actual cutoff f/u, own q
and first rho sectors. Only f/u are multiplied by 1-sigma(4q-2). After q>=3/4,
they vanish exactly, while the retained primitive and own mean produce a
nonzero radial q. The primitive is not restarted at Rin or at r_minus.

Leading Phi/Psi/own mean are enclosed by the existing original first-bridge
source on the selected prefix. No second leading increment is added on top
of those actual leading enclosures. Leading Q uses the (1-delta) recovery;
n=1 q uses its distinct (1+delta) recovery. The positive implicit swirl and
pressure-source contributions already enclosed by the source adapters remain
in the propagated errors. This interface adds no new pressure equation.

Every value, gradient and curl is a sum of references to 43 shared signed
terms. A term has a coefficient/error enclosure and an exponential scale
built from 14 log atoms. The original F0 anchor, time, delta, Lambda, geometry,
original hb and selected phase_in scales remain separate. In particular,
neither hb nor its inverse is assigned a floating-point physical value.

The six extra gradient terms carry the single inverse-width factor

`exp(-original_log_hb-log_phase_in) = 1/(original_hb*original_phase_in)`.

For a native singular radial coefficient c_rho_sing, they have coefficients
`2rho*c_rho_sing` in partial_r and `-2rho*Z*c_rho_sing` in partial_z. The
factor is phase_in, not the current phase. Geometric, delta-weight and F0-gauge
derivatives are added only to the regular sector, not duplicated in this
cutoff sector. All six singular coefficient slots are exactly zero outside
the cutoff transition, including the lower collar.

Point rows remain numerical proxies. Their absolute errors include omitted
finite-width contributions. They must not be interpreted as exact finite-width
centers or combined into an ordinary floating-point physical vector.

## Coordinate derivatives and structural divergence

At fixed physical time:

```
partial_r=(2rho/r)*partial_rho
partial_z=((1-Z^2)*partial_Z-2rho*Z*partial_rho)
          /(lambda^(1-delta)*(1-delta*Z^2))
```

F0's axial derivative and the small delta-weight derivative retain their own
log atoms. The cylindrical orthonormal gradient includes -u_theta/r and
u_r/r connection entries. Cartesian rotation is R*gradient*R^T; curl is
rotated as R*omega, using omega=(-u_theta,z, u_r,z-u_z,r, u_theta,r+u_theta/r).

Normalized meridional divergence is reported separately by temporal group
and radial width sector. The cutoff singular sector cancels by the same
own recovery:

`(1-delta*Z^2)*rho*q_rho_sing - 2rho*Z*u_rho_sing = 0`.

Thus the steep cutoff terms cancel in divergence while remaining present
in the physical gradient and curl.

## Executed evidence

The focused driver is `scripts/run_original_selected_gradient.py`; its saved
report is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_selected_gradient.json`.

The accepted run exits 0 in **146.362 seconds**. It evaluates twelve actual
selected physical packets and one source-scaled-time packet. Thirty
independently differentiated finite-width manufactured physical comparisons
cover velocity, cylindrical/Cartesian gradients and curl. The manufactured
cutoff has an explicitly integrated own primitive, with finite hb and
phase_in, so inverse-width derivatives are observable directly.

Actual-source checks pass 144 accepted native value/Z/radial overlaps,
48 normalized regular/inverse-width divergence enclosures containing zero,
152 exact-zero assertions outside support/transition, 72 inverse-width
scale checks, 80 core/Rin term seam overlaps, and 112 precision60/85
coefficient/log comparisons. Eight unsupported requests are rejected.
Rin Cartesian-gradient reference structures agree exactly, and the same
terminal primitive is used at lower fractions1 and2.

At Z=0.3, transition fraction5/8, the `uz_r_width:u` coefficient is
approximately -811.008 with absolute error below 3.713e-82. Its scale includes
1/(hb*phase_in) and the original physical velocity/radius weights; this is
not an ordinary physical gradient value. At lower fraction2, all n1 f/u
and cutoff-width terms are exactly zero while the n1 radial q value is
approximately -15.8784 with absolute error below 8.857e-85.

The 4779246-byte report pins 1721 source inputs and 27 implementation files.
Fresh saved-report inspection exits 0 in 2.270 seconds. The existing
read-only Luna/max worker independently audits phase transfer, original
source intervals, physical weights, width factors, history recovery and
rotations, and finds no material mathematical/source-fidelity issue.

## Remaining dependency

Tasks 14/15 remain IN PROGRESS and the full goal remains ACTIVE / INCOMPLETE.
The selected inner/lower cut/own physical first derivative layer is now the
local assembly route. Continue from r_minus with the actual modified leading
charts, complete signed five-defect integration and the selected reserved
repair, retain repaired own M/mean/Q/V/P, and assemble the later physical
charts. Do not enlarge the auxiliary uncut raw phase or treat the unmodified
bridge as the later actual field.

Global m13/m15 pressure tails, a common global N, physical Hessians and full
Cartesian C2 joins, axis/full-core gradients, outer chart inverses, global
energy/heat/stress/flat remainder, genuine n>=2 temporal recursion, smooth
summation, oscillatory corrections and corrected Cartesian NS remain open.
The current normalized divergence and local seam checks do not close those
gates or establish a measured blow-up law.
