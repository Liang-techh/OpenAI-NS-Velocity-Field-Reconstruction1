# Actual original n1 cutoff and retained lower-collar history — 2026-10-15

Subsequent progress: [actual selected cut/own physical velocity, gradient and vorticity](CURRENT_ORIGINAL_SELECTED_GRADIENT_2026_10_15.md)
now consumes these states and width sectors. The raw uncut physical APIs retain
their original scope. Read the newer document for current physical coverage.

## Evaluated original field and coordinates

`OriginalN1CutoffNumeric` now consumes the numerical original core and
admitted raw local n1 solution to evaluate the paper-defined flat velocity
cutoff and its own mean/Q/pressure recovery. It has two explicit charts:

```python
from openai_ns_reconstruction.original_n1_cutoff_numeric import OriginalN1CutoffNumeric
field = OriginalN1CutoffNumeric()
inner = field.evaluate('.3', '.625')       # q in [0,1]
exit = field.evaluate('.3', '1')           # Rin
lower = field.lower_collar('.3', '2')     # fraction in [1,2], ending at r_minus
```

Both coordinates represent the actual phase
`phase=fraction*(original_sc/4)`, with `rho=4exp(hb*phase)`.
The original selected sc is an exact binary singleton, retained by its
source tuple, rather than a new value selected from an enclosure. It is
Z-independent and smaller than 1/4. Thus the true selected inner endpoint
is `Rin=Ra*exp(hb*sc/4)`; the lower endpoint is
`r_minus=Ra*exp(hb*sc/2)`. These fractions are not the ordinary raw phase
used by the older first-collar numerical APIs.

The actual cutoff is `chi(q)=1-sigma(4q-2)`: exactly one through q<=1/2,
transitioning on (1/2,3/4), and flat zero from q>=3/4. It multiplies only
natural f and u (physical F1 and Uz1). The current source bindings are
`current_original_n1_inner_cutoff_functions.py` and its numerical bridge;
all accepted files and original parameter/source family remain unchanged.

## Raw local input and own history

Raw data are bounded only on the already admitted selected prefix
`0<=phase<=sc/4`, using the full core Z2 inlet, actual six-row known
sources, singular-source cancellation and the bounded regular operator.
A directed sigma prefix and centered finite-width variation cover the
raw local change. No uncut raw domain is enlarged to ordinary phase1.

Let a=sc/4, m0 be the numerical natural core mean, and

`C(q)=integral_0^q chi(xi) dxi`.

The normalized own primitive is

`J(q)=rho*mean_u=4m0+hb*a*integral_0^q rho(a xi)*chi(xi)*u_raw(a xi) dxi`.

The leading history coefficient is `4u_core*C(q)`; the directed error
retains raw u variation and actual rho-4 variation before multiplying by
the original hb*a. The defining primitive is never reset at a seam.
`C(q)=q` through 1/2 and **C(q)=5/8** from 3/4 onward, by the original
sigma reflection identity. Interior mass uses directed original-cutoff
quadrature with a fourth-derivative remainder. The same terminal history
is reused throughout q>=3/4 and the lower collar.

The physical primitive M1 equals `epsilon^2*J`. Recovered fields are

```
mean_u=J/rho
k=mean_u-u_cut
q=[2Z*u_cut-(1+delta)Z*mean_u-(1-Z^2)*(mean_u)_Z]/(1-delta*Z^2)
```

Mean/k and cut f/u are returned through native Taylor Z2, and q through
Z1. Entry k is D_Z^k/k!: ordinary second derivatives and errors are twice
entry2. q_ZZ remains outside this interface.

Point approximations use the core primitive plus a separately exported
nonzero history coefficient. The unmaterialized hb*a contribution is
covered by the full absolute error. The point value is a numerical core
proxy, not an exact finite-width field value. Explicit factor tags and
source bounds retain the actual positive original history.

## Cutoff derivatives and pressure

First radial rows retain two width sectors:

`key a,b -> original_hb^a * original_phase_in^b * enclosed coefficient`.

The (-1,-1) sector follows the exact product rule

`chi_rho=-4sigma_prime(4q-2)/(hb*phase_in*rho)`.

It is retained separately from chi*f_rho and chi*u_rho. No inverse hb or
inverse selected phase is assigned an ordinary physical value. Own mean
rho recovery uses `(u_cut-mean_u)/rho`; q_rho is derived from that same
history and returns Z1. It includes the cut-u singular sector once.

Pressure integrates the known leading N1p source **once**, together with
the cut F coupling, from the nonzero core pressure inlet. In natural
`p=P1/(epsilon*Pstar^2)` units, the coupling per drho is

`2S*Phi*f_cut/(epsilon^3*Pstar^2)`, `S=epsilon^2*F0^2`.

The known inverse-width source is canceled against dx before enclosure;
the regular known source and actual cut-F history have separate finite-
width error budgets. N1p remains active after F1/Uz1 become zero. The
numerical pressure center omits regular integral contributions only as
an explicitly charged absolute error; it is not an exact terminal
pressure integral and does not close global pressure moments.

## Executed evidence

Driver: `scripts/run_original_n1_cutoff_numeric.py`.
Report: `experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric.json`.

The accepted run exits 0 in **143.999 seconds**. Fifteen selected-inner
packets cover Z=-0.3,0,0.3 at q=0,1/2,5/8,3/4,1; nine lower packets use
fractions 1,3/2,2. It passes four independent identities, 108 exact flat-zero
velocity/radial coefficient slots, 186 core and inner/lower seam error
comparisons, nine identical terminal-primitive source packets, 90 precision
60/85 coefficient comparisons and six unsupported-input rejections.

At Z=0.3 and lower fraction2, cut f/u are exactly zero. The retained mean
is approximately **6.33600000000000000000012759** (error below 3.499e-85),
q is **-15.87839999999999999999975198** (error below 8.857e-85), and q_Z is
**50.68800000000000000000073234** (error below 1.311e-83). Pressure is
approximately -3.95900164979066614032 with error below 1.096e-26.

The own primitive increment coefficient per original hb*phase_in is
**31.68000000000000000000063796**, with error below 3.251e-84. At the
transition q=5/8, the u_rho inverse-hb/inverse-phase coefficient is
**-101.37600000000000000000204147**, with error below 4.641e-83. These
are natural/factored coefficients, not ordinary physical derivatives.

The 1491003-byte saved report pins 1719 source inputs and 16 implementation
files. Fresh saved-report inspection exits 0 in 2.269 seconds.

Independent symbolic checks cover the cutoff radial chain rule, own
q_rho recovery, exact meridional divergence cancellation and pressure
normalization. The existing read-only Luna/max worker audits the selected
domain, source measures, primitive variation, pressure factors and radial
width sectors; it finds no material factor, sign or domain error.

## Correct next dependency

The earlier proposal to extend raw uncut transport to ordinary phase1 is
an auxiliary analytic extension. It is not the next required matched-field
step: the actual two n1 velocities have already been cut before Rin.
The new interface follows that cutoff and the zero-velocity own-history
lower collar instead, without rebuilding or relaxing raw-domain bounds.

Later leading modifications can start at r_minus. Numerical continuation
beyond that point must consume the actual modified leading charts and
retain their pressure source. The unmodified first bridge must not be
substituted for the actual later field merely because its auxiliary
leading values can be computed through ordinary phase1.

Next connect these actual cut/own states and their width sectors to
physical velocity/gradient assembly, then evaluate the actual later
modified charts, complete signed five-defect integrals and reserved
repair/own recovery. Global pressure m13/m15 tails, common global N,
heat exterior/energy, admissible stress/flat remainder, genuine temporal
n>=2 coefficients/repair/summation, oscillatory cancellation and corrected
Cartesian NS validation remain open. Tasks 14/15 remain IN PROGRESS;
full goal ACTIVE / INCOMPLETE.
