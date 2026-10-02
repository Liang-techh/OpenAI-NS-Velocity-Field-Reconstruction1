# All 33 leading source charts in physical Cartesian coordinates

The same admitted analytic core, repaired five-moment family, pressure datum,
annuli and exact heat exterior now have one physical source-bound assembly.
Every chart has spatial derivatives through total order four and the first
time derivative at a fixed physical position. The core and axis use native
nonsingular Cartesian formulas. These are signed factored enclosures of the
source, not selected nonlinear point coefficients or measured blow-up dynamics.

## Coordinates and physical units

The original map is

```
tau = T-t > 0
lambda^2-lambda^(2*delta)*z^2 = tau
Z = z/lambda^(1-delta)
R = (x^2+y^2)/(2*lambda^2)
ur = lambda^-1 * Ur
utheta, uz = lambda^(-1-delta) * (Utheta, Uz)
p = lambda^(-2-2*delta) * P
```

All annular source rows are ordinary derivatives in `(logR,Z)`. The input
provider divides its fully differentiated profile by a FIXED basepoint unit.
This includes angular amplitude and radial prefactor derivatives. The assembly
restores the scalar unit after the linear coordinate operator; it does not
differentiate that normalization again. The actual patch's radial unit uses
`sqrt(Rm/2)`, and is converted to the current-R unit explicitly.

Moving cylindrical basis derivatives are included in Cartesian derivatives.
The time operator is applied at fixed physical `(x,y,z)`, with the basis fixed,
and differentiates the implicit lambda. It is not a time derivative holding
similarity coordinates fixed. For uniform whole-Z bounds, `lambda>=sqrt(tau)`
bounds negative lambda powers; interior Z boxes also expose an exact logarithmic
lambda enclosure. The whole-Z endpoints represent physical infinity.

## Microscopic source factors

The bridge/switch source widths are positive formal sources, not inverted
numerical caps. For the four microscopic charts, `D_logR^k=hb^-k D_phase^k`
is combined with every signed source term before any final bound. The macro
bridge has already differentiated logR rows and receives no extra hb^-k.

The exact source identity

```
2*log(u/Pstar) = log(2) + logR + 2*logF0 + 2*logphi - 2*logPstar
```

converts the original four factor bases into shared physical factors. It is
used algebraically before evaluating enormous source logs. The original
log-radius expressions retain positive hb variations even when a numerical
radius box cannot resolve them. Numerical logR bounds still use the proved
positive-width cap and do not preserve correlations with the formal hb source.
That distinction is explicit in each bridge/switch packet: the radius bound
is conservative, and exact correlated evaluation has not been implemented.
No numerical width cap is inverted. Huge positive or negative source
exponentials are not materialized.

Source precision is explicitly labeled: the four microscopic charts and
macro bridge use `uncapped_factored_rows`; other annuli consume
`provider_prebounded_mixed_rows`. The latter are conservative provider
enclosures and may already lose sharpness. The assembly does not claim to
recover correlations erased by those earlier enclosures.

## API and reproduction

```python
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    CompliantGlobalPhysicalAssembly,
)
field = CompliantGlobalPhysicalAssembly()
packet = field.evaluate('O2_axial', Z='.5', coordinate='.5',
                        log_tau='-10', theta='.7')
axis = field.evaluate('core', Z='.5', coordinate=0,
                      log_tau='-10', axis=True)
```

The API accepts the original explicit chart coordinate and log(tau). It returns
35 Cartesian spatial multiindices for velocity/pressure, fixed-position time
rows, formal scale bases and scope metadata. It does not automatically select
a chart from an arbitrary physical point or return reconstructed point values.

The ordered pipeline has 111 modules, including the subsequent fresh core
coefficient rebuild and anchored amplitude. The focused physical-map stage is

```
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage globalphysical
```

## Independent evidence

The checker binds every current source hash and checks all 33 complete chart
domains plus the nonsingular axis. A separate finite, nonzero, divergence-free
fixture uses an independently solved implicit lambda and direct Cartesian
differentiation: 140 spatial derivatives and four fixed-position time
derivatives are enclosed. Different fixed source units exercise normalization.

Another independent fixture tests 120 signed factor rows at extreme source
scales: 60 microscopic rows with hb conversion and 60 macro rows without it.
The original Mz recovery identity proves divergence symbolically; a direct
Cartesian fixture checks it separately. These fixtures validate operators and
factor algebra, not point coefficients of the admitted nonlinear core.

## Remaining critical work

1. Continue from the fresh compliant core coefficient rebuild described in
   CORE_COEFFICIENT_REBUILD_2026_10_02.md. This supplies directed radial rows
   and controlled infinite remainders, not selected source parameter values.
   The anchored G amplitude now has directed logarithmic values in
   ANCHORED_AXIS_AMPLITUDE_2026_10_02.md. Keep source/pressure uncertainty
   separate from numerical errors; never choose midpoints as exact solutions.
2. Resolve original integral sources and the implicit five-bump branch into
   usable values with error bounds. Preserve shared histories and all five
   functional terminal conditions when composing an actual field evaluator.
3. Provide physical-point chart selection with logarithmic support for extreme
   radii, microscopic boundary separation, explicit axis handling and a clear
   response for points not numerically resolvable at the requested precision.
4. Establish the construction's required physical spacetime energy domain and
   terminal-time bounds. The original unlocalized whole-space energy is
   infinite. The present physical derivative maps do not fix that fact.
5. Construct actual divergence-form admissible stress and an independently
   bounded flat remainder from the assembled background. Report regionwise
   cone margins, maximum norms, volume L2 norms and their scale dependence.
6. Then implement n=1 and n>=2 recovery, per-order moment repair, curl-based
   cutoffs, finite-order remainder and smooth summation. Coordinate scaling
   is not coefficient recursion. Mean/oscillatory corrections and full forced
   Navier-Stokes residual acceptance remain later stages.
