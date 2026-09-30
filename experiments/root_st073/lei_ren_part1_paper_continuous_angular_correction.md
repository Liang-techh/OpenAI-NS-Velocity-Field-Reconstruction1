# Continuous angular correction

`install_continuous_angular_correction(correction, precision=None, basis_precision=100)` mutates the supplied `AngularCorrection` in place and returns a `ContinuousAngularCorrectionProvider`. The correction object keeps its identity, while `weights`, `coefficients`, `relative_bump`, `value_jet`, and `weighted_atom` use one shared `ContinuousAxialBump`.

The provider exposes:

```python
provider.coefficients(Z)
provider.relative_bump(t, Z, coefficients=None)
provider.value_jet(t, Z, coefficients=None)
provider.weighted_atom(lam, power=1, lower=None, upper=None)
provider.correction_integral(lam, Z, power=1, lower=None, upper=None)
provider.pressure_correction_integral(Z, lower=None, upper=None)
```

`value_jet` returns MP values with keys `h`, `h_t`, `h_tt`, `h_Z`, and `h_tZ`. It checks both compact supports before obtaining coefficients, so calls outside `[-3-ell,-3+ell]` and `[-1-ell,-1+ell]` return exact zero jets. `weighted_atom` integrates `exp(lam*s) beta(s)**power` directly over the clipped local support. Flat endpoint scaling is used for endpoint intervals; no full-minus-partial subtraction is used.

The coefficient receipt keeps the legacy signed-log fields (`d1`, `d2`, `d1_over_r`, `d2_over_r`, `log_r`, and `log_s`) and adds `Z_exact`, `d1_Z`, `d2_Z`, `d1_over_r_Z`, and `d2_over_r_Z`. Coefficients are cached by the full MP `Z` string. `relative_bump` keeps the signed tiny correction in `relative_correction` and `arbitrary_exponent_value` separately from `log1p_relative_correction`, since `1+h` can round to one.

The continuous weights are evaluated at `integral_working_precision = basis_precision` and are then consumed by the outer MP coefficient solve. This records the actual basis/integral precision instead of implying that a 400 digit outer context certifies a 100 digit bump normalizer. The continuous weights are

```text
A = integral exp((1-mu) s) beta(s) ds
B = integral exp(-(1+2 mu) s) beta(s) ds
D = integral exp(-(1+2 mu) s) beta(s)^2 ds.
```

For the inputs, `C_r` and `C_s` are the central `heat_defects(schedule, 0, order)` first-Taylor constants. The provider uses `r_H=C_r(1-Z^2)`, `s_H=C_s(1-Z^2)` and the analytic derivatives `-2 Z C_r`, `-2 Z C_s`. Its preheat term solves the same baseline `solve_ivp` ODE as `heat_defects` once from `incoming_flattening_X`; fixed quadrature nodes retain the baseline values as floats, while `sigma`, the `Z` weights, and the final integrals use MP arithmetic. This is numerical baseline and heat Taylor input data, with quadrature error and waiting-root error still separate.

Run the focused checks with:

```text
C:/Users/z5242/AppData/Local/Programs/Python/Python311/python.exe experiments/root_st073/lei_ren_part1_paper_continuous_angular_correction.py
```

The generated JSON records derivative checks on normalized coefficients and on `h` jets, even/odd parity, compact support, and the `Z = -1, 0, 1` endpoint cases. The recorded scope keeps `full_outer_closed=false`, `pressure_cancellation_exact=false`, and `quadrature_enclosure_certified=false`.
