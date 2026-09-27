# Outer axial repair across scales: geometry transfers, moments do not interpolate

The frozen-inner/outer-quadratic-axial repair was solved independently
at k=15 and 19, then combined with k=11 into one smooth time-dependent
field using the existing scale partition. All full momentum evaluations
include time derivatives; intermediate scales were not fitted.

| Knot | Max moment, order 96 | Max moment, order 128 | Held-out momentum peak | Inner cone |
| --- | ---: | ---: | ---: | --- |
| 11 | 6.00e-6 | 3.00e-5 | 1.44606e7 | 9/9 |
| 15 | 2.16e-4 | 4.88e-4 | 9.20434e8 | 9/9 |
| 19 | 1.74e-3 | 6.70e-3 | 5.84641e10 | 9/9 |

The finest knot does not pass an absolute 1e-3 moment gate. Training
moments near 1e-9 do not override the independent replay. Quadrature
and finite-difference cancellation require separate investigation at
that scale. No complete momentum gate is close to passing.

| Intermediate k | Max moment, order 48 | Max moment, order 96 | Momentum peak, order 96 | Inner cone |
| --- | ---: | ---: | ---: | --- |
| 13 | 0.238582 | 0.238559 | 7.78824e7 | 9/9 |
| 17 | 3.540215 | 3.540573 | 4.90658e9 | 9/9 |

The retained inner geometry survives both knot transfer and the tested
intermediate times. Moment compatibility does not: smooth interpolation
of independently optimized coefficients leaves substantial intermediate
moment defects. These defects are stable under integration refinement
and cannot be treated as only the finest-knot numerical noise.

Next derive the dependence of each moment on outer coefficient values
and their time derivatives. Determine whether angular moment transport
and axial moment constraints admit a coupled differential-algebraic
continuation with the inner coefficients fixed. Verify the response using
analytic coefficient-time derivatives against direct field replay.
Adding more independent knots alone would not establish a consistent
recursive update or interval bound. Also measure outer stress realizability
and residual contraction before returning to the Section 7 pulse solve.

All candidates remain diagnostic, accepted=false. The complete spatial
maximum, volume L2, finite-energy closure, full pulse interval, and
recursive residual contraction remain unestablished.

Reproduce:
```text
python experiments/root_st073/midplane_outer_axial_repair.py --k 15
python experiments/root_st073/midplane_outer_axial_repair.py --k 19
python experiments/root_st073/midplane_outer_axial_interscale.py
```
The default no-argument repair command still reproduces the k=11 trial.
