# Conditional continuous coefficient enclosure

The shared runtime now exposes `coefficient_enclosure(panels=4096, precision=80)`.
It propagates outward basis and complete pulse bounds through the actual two-row
linear solve and energy quadratic. It does not change the installed coefficients.

The determinant uses `-2*exp(-2+6*mu)*sinh(mu)*B1*B2`, with outward bounds on
`sinh(mu)`, rather than subtracting almost equal matrix products. The inverse
numerators still use interval arithmetic; dependency loss can widen them.
Such width does not prove that a solution is absent.

Let `c0=-inv(M)*base` and `cp=-inv(M)*pulse`, so `c=c0+a*cp`. With
`K_j=exp(-26+{6,2}_j*mu)*energy_gram`, the actual quadratic is

```
A = Kp + mu*sum(K_j*cp_j^2)
B = 2*mu*sum(K_j*c0_j*cp_j)
C = mu*sum(K_j*c0_j^2) - target.
```

The code requires outward `A>0`, `C<0` and discriminant `D>0`, establishing a
unique positive branch for every underlying atom choice. It evaluates
`(-B+sqrt(D))/(2*A)` and, when its denominator is strictly positive,
`-2*C/(B+sqrt(D))`, then intersects these two guaranteed root enclosures.
Interval squares preserve positivity, including intervals crossing zero.

In the installed fixture, nominal amplitude and both coefficients are contained
at 1024 and 4096 panels. Relative coefficient widths decrease from approximately
10.36% to 2.49%. Pulse energy now uses independent positive interval bounds,
including the nested startup primitive. The relative amplitude width is about
`2.05e-9` at 4096 panels. This is still **conditional on fixed incoming rows
and energy target**, not a measurement of total amplitude accuracy. These two
inputs are exact stored dyadic parameters for this conditional calculation;
their originating numerical errors are not bounded. The low-level `solve_atoms`
API also accepts caller-supplied input intervals, including a supplied Kp interval.

Fixed materialized coefficient residual ranges include zero for both rows.
This proves compatibility with the integral bounds, not exact mean cancellation.
The actual nonzero runtime terminal residuals are retained. No finite-energy,
source-uncertainty, full momentum or scale-recursion certificate is issued.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_solve_enclosure.py`.
The JSON saves exact dyadic endpoints, conditional parameter scope, coefficient
widths and fixed-coefficient balance ranges. Decimal displays are approximate.

Next: enclose continuous incoming cutoff integrals and energy target, tighten
correlated pulse/basis bounds, then propagate genuine input intervals before
claiming complete mean and energy closure.
