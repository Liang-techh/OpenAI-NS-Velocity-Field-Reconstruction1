# Real axis norm bounds for Lei–Ren Part I

This receipt concerns only the real interval `Z ∈ [-1, 1]` and only the axis
factor

\[
G'(Z)=\frac{L(Z)H(Z)}{H(Z)^2+\sigma_0^2},\qquad
H(Z)=\frac{1-\delta}{2}Z+(1-Z^2)(4Z+j),
\]

with `L = 1 - δ Z²` and `σ₀ = j / 500`. It does not bound mixed radial
derivatives, the complex `A_Ω` norm, the pressure amplitude, the nonlinear
core, or a PDE/cone quantity.

The module uses the exact coefficient form

\[
H=-4Z^3-jZ^2+(9/2-\delta/2)Z+j
\]

and bounds each derivative on `[-1,1]` by the sum of the absolute values of
its polynomial coefficients. For `f(h)=h/(h²+σ₀²)`, the quotient identities
are

\[
f' = \frac{\sigma_0^2-h^2}{Q^2},\quad
f'' = \frac{2h(h^2-3\sigma_0^2)}{Q^3},\quad
f''' = -\frac{6(h^4-6h^2\sigma_0^2+\sigma_0^4)}{Q^4},
\quad Q=h^2+\sigma_0^2.
\]

Using `Q ≥ σ₀²` and the real maxima of the normalized rational factors gives

\[
|f|\leq(2\sigma_0)^{-1},\quad |f'|\leq\sigma_0^{-2},\quad
|f''|\leq\frac{3}{2}\sigma_0^{-3},\quad
|f'''|\leq48\sigma_0^{-4}.
\]

Leibniz differentiation of `G' = L f(H)` then yields bounds for
`G'`, `(G')'`, `(G')''`, and `(G')'''`, which are returned by
`axis_norm_bounds(axis)` as `G_derivative_bounds`. The conservative value bound
is `|G| ≤ 2 sup|G'|`; actual `G(-1)` and `G(1)` are retained separately.

For the sample `j = 1e-14`, `Lambda = 1e36`, `delta = 1e-200`, the receipt
reports

```text
G_value_upper = 5.0e16
A_axis_upper  = 1.3476375000000049851e108
```

The pointwise MP derivative replay at `Z = -0.8, -0.2, Z0, 0.2, 0.8`
covered orders `0` through `3` at every sample. This is finite real-axis
evidence only; it is not a theorem certificate for the full construction.
