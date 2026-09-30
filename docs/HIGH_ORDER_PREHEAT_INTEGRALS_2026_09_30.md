# Taylor-remainder pressure integral enclosures

The first two variable preheat transitions now have a higher-order interval
quadrature path. It evaluates the same stored-parameter continuous formula
and compares the resulting true integral interval with the retained Gauss
mass; it does not modify pressure normalization or replace the source.

## Method

`IntervalTaylor` carries Taylor coefficients f^(k)/k! with directed interval
arithmetic for products, reciprocal and exp. When its constant coefficient
is an entire cell, the resulting derivative coefficients enclose derivatives
at every point in that cell. For a cell centered at c with radius r, integrate
the midpoint polynomial through order N-1 and add the absolute remainder

    2 r^(N+1)/(N+1) * sup_cell |f^(N)/N!|.

The flat switch is evaluated using the stable reciprocal logistic branch.
Interior primitive panels use this Taylor-remainder formula. The two panels
touching0 or1 use monotonicity and sigma symmetry rather than singular
derivative formulas at the flat endpoints. J(1)=1/2 remains exact.

The radial integrator builds log-amplitude derivative coefficients using
J'=sigma, then exponentiates the interval Taylor data. At the center, J comes
from the enclosed primitive prefix grid. Rounded stage lengths are retained
through the affine change of variable. Radial edge panels integrate the
outer exponential branch exactly, with a positive flat correction bound.
This avoids leaving a first-order density bound at the edges.

## Commands and scope

```text
python experiments/root_st073/lei_ren_part1_paper_interval_taylor_fixture.py
python experiments/root_st073/lei_ren_part1_paper_high_order_preheat_fixture.py
python experiments/root_st073/lei_ren_part1_paper_high_order_preheat_integrals_check.py
```

The actual check uses order12 and radial panel counts32 and64. It records
the signed error of the stored192-point Gauss value rather than assuming
that Gauss value is exact. At64 panels the normalized mass error upper bounds
are 2.21850e-13 and 2.97946e-41. Both interval widths shrink from32 panels.
The first error upper bound improves the prior first-order bound0.00240550
by more than ten orders of magnitude. These errors multiply q^-2 exactly,
so the receipt also gives compact axial value/first/second error budgets.

The independent Taylor algebra fixture checks rational/exponential derivative
coefficients and odd/even-order symmetric integration. A separate switch
fixture checks27 derivative coefficients against scalar differentiation,
cell derivative bounds, and three primitive integrals against independent
scalar quadrature. Both pass; the quadrature comparisons are regressions,
while the derivative remainder bounds define the enclosure.

The path now also covers both post-flatten steep transitions; see
STEEP_PREHEAT_INTEGRALS_2026_09_30.md for the signed correction and rounded
origin handling. Variable-beta derivative
errors, original parameter derivations and source/core propagation remain
open. The resulting pressure stage error is not a five-moment closure
certificate or evidence of temporal coefficient recursion.
