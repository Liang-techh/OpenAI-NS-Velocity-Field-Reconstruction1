# Continuous axial pulse energy receipt

This adds the public function

```text
continuous_pulse_energy(precision=100, quadrature_order=None) -> mpmath.mpf
```

It evaluates

```text
Kp = integral from 0 to 11 of exp(-2 xi) * gp(xi)^2 dxi
```

using the existing `ContinuousAxialPulse.value_jet` source definition. The
function is importable by a future continuous solver and is not installed into
the global profile.

The support is split at the source transition points:

- `[0, .02]`: MP Gauss-Legendre quadrature of `value_jet`; the shared source
  class performs its startup primitive with nested MP quadrature.
- `[.02, 10]`: analytic antiderivative of
  `exp(-2 xi) (xi - .01)^2`, which is the exact plateau branch of `value_jet`.
- `[10, 11]`: MP Gauss-Legendre quadrature of `value_jet` through the cutoff.

All nodes and weights are `mpmath` values. At precision 100, order 80 gives

```text
Kp = 0.245049620200694480513590319716...
```

At a fixed order of 80, increasing precision from 70 to 100 changes the result
by about `7.54e-72` relative. Increasing the order from 80 to 112 at precision
100 changes the result by about `2.76e-21` relative. The legacy
`pulse_K_p(quadrature_order=192)` float result
is `0.24504962020069207`, about `9.84e-15` relative below the order-112 MP
value.

The startup and cutoff pieces are retained separately in the JSON receipt; the
plateau contribution is the dominant term. These are bounded quadrature
refinements, not rigorous enclosures. The result does not certify a global
finite-energy construction, a continuous incoming solve, profile installation,
or recursive closure.

Run command:

```text
python experiments/root_st073/lei_ren_part1_paper_continuous_pulse_energy.py
```

The machine-readable receipt is
`lei_ren_part1_paper_continuous_pulse_energy.json`.
