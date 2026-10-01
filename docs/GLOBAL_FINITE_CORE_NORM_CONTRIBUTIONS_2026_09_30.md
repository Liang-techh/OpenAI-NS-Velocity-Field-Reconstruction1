# Whole-axis finite-core norm contributions — 2026-09-30

The accepted coherent preheat pressure now has derivative-error envelopes
through order24 on the entire closed real interval [-1,1]. The q-power
envelope extends to radius1 because 1+Z^2>=1, and the component masses
are radial quantities independent of the chosen axial radius. The existing
radius.8 receipts remain unchanged; the new receipt is
`lei_ren_part1_paper_global_pressure_high_derivatives.json`.

`lei_ren_part1_paper_global_finite_core_bounds.py` encloses the degree18
finite radial reconstruction on x=R/R_a in [0,1], Z in [-1,1], using
analytic axis jets and the same accepted finite preheat datum. All accepted
component masses are checked against the source-alignment receipt before
the recurrence runs. The conditional pressure integral error is included.

It supplies all mixed derivatives with total order<=3 of F, Uz and P.
Each radial coefficient has enough axial entries for these derivatives.
Pressure restores the integral of the complete F-squared polynomial,
including radial powers beyond the pressure-recursion prefix. No pressure
datum is reset or replaced.

The sum norm contributions have logarithmic upper bounds approximately
-5e151 for F, 2.07944154167983718 for Uz and 34.8920759490302556 for P.
These contribute to K but do not by themselves establish the complete K.
In particular, an upper bound for |F| cannot establish a bound for 1/F.

At the actual inlet R_a=4/Lambda, the independently bounded derivative
orders0,1,2 of Uz-4Z and m_z/R_a-4Z give

||Uz-4Z||_C2 + ||m_z/R_a-4Z||_C2 <= 2.00000115626372515e-14

for this finite model on the whole closed axial interval. The affine axis
part is subtracted structurally before interval evaluation; subtracting
two independent [-4,4] intervals would destroy this estimate. The leading
deviation is the supplied j=1e-14 in each of the two quantities.

Fifteen containment checks compare the global derivative/deviation bounds
against the independently stored actual Z=.3 inlet receipt, including its
conditional pressure error. They are regression checks; the uniform
enclosure comes from the directed recurrence and analytic interval data.

## Remaining dependencies

- Certify the infinite radial-series remainder and reconcile it with these
  finite bounds, rather than treating a degree18 polynomial as the theorem core.
- Obtain positive lower bounds for the normalized finite/exact core F,
  then reciprocal-core C3 bounds; retain the supplied logarithmic amplitude.
- Bound the frozen profile D_f,E_f,1/D_f and moments over its required
  radial interval and the entire axial interval.
- Complete the global K and c_star/K1 conditions before claiming the
  edge-factored comparison-error estimate or a finite-width cone certificate.
- Continue five-moment closure, matching, heat exterior, flat remainder,
  time recursion and oscillatory stress correction.

The separate `lei_ren_part1_paper_uniform_axis_log_norms.py` now bounds
axis-only F and 1/F through three actual Z derivatives on [-1,1]. It keeps
the common amplitude outside the normalized derivative recurrence, using
opposite signs for the F and reciprocal equations. Its reciprocal C3 log
upper bound is approximately 5e151, below the assumed log K=1e152 for this
axis contribution. This does not bound 1/F away from the axis in the core.
Forty independent modest-parameter derivative and reciprocal checks pass.
Logarithms and positive norm sums are computed with directed arithmetic;
exact binary interval endpoints are stored separately from display strings.

Original source-parameter derivation, infinite radial error and full K are
explicitly marked uncertified. This is not temporal coefficient recursion.
