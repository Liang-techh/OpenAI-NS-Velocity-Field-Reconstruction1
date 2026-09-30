# Continuous incoming axial primitive

`lei_ren_part1_paper_continuous_incoming.py` provides the source incoming
axial factor in source coordinates

```text
y = log(R) - log(R_ref)
Uz(y, Z) = 4 Z c(y).
```

The cutoff is the source schedule's flat transition:

```text
c(y) = 1                                  for y <= 1,
       1 - sigma(log(y) / Md)             for 1 < y < exp(Md),
       0                                  for y >= exp(Md).
```

`sigma` and its derivative are the phase-stable `sigma_pair` definition used
by `continuous_axial_pulse.py`.  The shared candidate uses `Md = .5`; the
provider accepts every finite positive `Md` representable at the selected MP
precision.

## API

Construct `ContinuousIncomingAxial(Md=".5", precision=160,
quadrature_order=48)`.  The selected precision belongs to the instance and
the `Z` argument is kept as an mpmath value; it is never converted to a
binary64 cache key.

Point methods are `c(y)`, `c_y(y)`, `Uz(y, Z)`, `Uz_y(y, Z)`, and
`Uz_Z(y, Z)`.  `values(y, Z)` returns the same point jet as a mapping.

For `y >= 0`, the cumulative definitions are

```text
I_z(y)   = 4 Z [1 + integral_0^y exp(s) c(s) ds],
I_uz2(y)= 16 Z^2 [1 + integral_0^y exp(s) c(s)^2 ds].
```

For `y < 0`, they are `I_z = 4 Z exp(y)` and
`I_uz2 = 16 Z^2 exp(y)`.  The mean and its analytic jets are
`mean = exp(-y) I_z`, `mean_y = Uz - mean`, and `mean_Z` as returned by the
provider.  `weighted_integral(power=1, upper=None, cutoff_power=1)` exposes
the underlying weighted atom.  `upper=None` means `exp(Md)`; a partial
transition endpoint is accepted and values beyond the endpoint saturate.
`full_incoming_rows(Z)` (also `incoming_rows`, `full_rows`, and `moments`)
returns the terminal `I_z` and `I_uz2` rows.

The unit plateau is integrated analytically.  The finite transition is
integrated after the change of variables `y = exp(Md t)`, `0 <= t <= 1`, by
mpmath quadrature.  This makes endpoint splitting explicit and allows an
order refinement through the `order` keyword or instance default.

## Focused receipt

Run:

```text
python experiments/root_st073/lei_ren_part1_paper_continuous_incoming.py
```

The generated JSON receipt uses 120 decimal digits, orders 20 and 36, and
`Z = .375`.  The five-point primitive derivative replay relative errors at
`y = .5`, `1.2`, and the interior midpoint are approximately
`3.33e-22`, `3.33e-17`, and `5.04e-16`.  The order-20/36 changes for the
`exp(y)c`, `exp(y)c^2`, and a half-power diagnostic are recorded explicitly;
they were zero at the printed 60-digit level for this smooth candidate.
Partial queries include the plateau join, two interior transition points,
the exact cutoff endpoint, and a point after the endpoint.  The full rows in
the receipt are approximately `I_z = 5.42218924098` and
`I_uz2 = 7.97958612162` for this diagnostic `Z`.

The receipt records continuity samples at both flat joins and keeps
`quadrature_enclosure_certified`, `global_energy_certified`, and
`recursion_certified` false.  This module does not certify a quadrature
enclosure, a global finite-energy construction, recursive matching, or
installation into the global velocity profile.
