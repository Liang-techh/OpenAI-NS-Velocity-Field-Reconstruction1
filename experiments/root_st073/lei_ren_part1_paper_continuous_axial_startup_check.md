# Continuous axial pulse startup check

This diagnostic exercises the startup branch of `ContinuousAxialPulse.value_jet`.
For `q = 50 xi <= 1/2`, it evaluates the positive primitive

```text
int_0^q sigma(u) du
```

with the endpoint change of variables

```text
u = q / (1 + q^2 v),   v in [0, infinity).
```

Writing `phi(u) = 1/u^2 - 1/(1-u)^2`, the implementation uses

```text
phi(u) - phi(q)
  = 2 v + q^2 v^2
    - dq * (2-u-q) / ((1-u)^2 (1-q)^2),
 dq = -q^3 v / (1 + q^2 v),
```

and integrates the resulting positive scaled integrand on `[0, 1, 4, 16, infinity]`.
The same `sigma_pair` definition is retained, and no bulk-minus-endpoint
subtraction is used. The numerical quadrature remains uncertified.

The focused run checks positive values at `xi = .01, .001, .00001`, compares
the transformed primitive with direct original-coordinate quadrature at `.005`
and `.01`, checks logarithmic finite differences with `h = xi^4`, verifies
continuity at `.01` and `.02`, and evaluates the full pulse energy with the
existing MP energy helper.

The generated receipt is [`lei_ren_part1_paper_continuous_axial_startup_check.json`](lei_ren_part1_paper_continuous_axial_startup_check.json).
The run passed all focused checks. At working precision `60` and Gauss order
`48`, the full energy was approximately

```text
0.2450496202006948234631778106780410702
```

The stored legacy float `K_p` is approximately
`0.24504962020069207`; the nominal relative difference is approximately
`1.1236349501605e-14`. This is a comparison of numerical values, not an
enclosure and not a replacement of source coefficients.
