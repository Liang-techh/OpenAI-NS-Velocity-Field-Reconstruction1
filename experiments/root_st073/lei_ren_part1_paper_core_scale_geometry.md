# Core scale geometry receipt

This bounded experiment uses `build_joined_field().inner.physical_field()` at
the actual regular core point

```text
Z = 0.3
R = 1 / Lambda = 1e-36
Lambda = 1e36
delta = 1e-200
nu = 0.01
```

The point is `R / (4.1 / Lambda) = 0.243902439...`, inside the polynomial
core. For each `alpha` in `0, 2, 4, 6`, the script forms the arbitrary MP chart
scale

```text
logq = -alpha / delta
tau = q * (1 - Z^2)
```

It maps the actual core state through `physical_chart`, calls the actual
`LocalCoreExitField.evaluate` on the returned `(x, y, z, t)`, and measures
coordinates and cylindrical velocities from the roundtrip result. Every
extreme quantity is retained as a sign, `log_abs`, and arbitrary exponent
value in the JSON receipt.

The four rows all return `region = core`. The measured aspect ratio gain has
log gains `0`, `1`, `2`, and `3`, matching the local chart prediction
`alpha / 2`; the largest ordinary ratio is therefore `exp(3) =
20.085536923...`. The radial coordinate has log-ratio gain `-alpha / (2 delta)`; the axial
coordinate has gain `-alpha * (1 - delta) / (2 delta)`. Their difference
`alpha / 2` produces the measured relative elongation. At alpha `6` the
large coordinate log ratios are near `-3e200`. The actual
`utheta` and `uz` each have log gain
`alpha * (1 + delta) / (2 delta)`, and their ratio remains locally constant.
The pointwise winding density
`utheta / (2 pi rho uz)` has log gain near `alpha / (2 delta)`.

The roundtrip relative velocity errors are zero at alpha `0` and have log
absolute sizes around `-572` to `-568` for alpha `2`, `4`, and `6` at the
available 443 digit working precision. Coordinate roundtrip errors have log
absolute sizes around `-573` to `-569` for the extreme rows. Absolute errors
can look large in ordinary notation when a velocity magnitude has a log of order
`1e200`; the signed-log values are the stable evidence.

The axial vorticity component is available from the actual recovered core jets:

```text
omega_z = q^(-(2 + delta) / 2) * (2 F + 2 R F_R)
```

This receipt measures that component and its scale ratio. The full cylindrical
curl is unavailable because the core state does not expose physical `Ur_R` and
`Ur_Z` jets. The winding quantity is a pointwise density, not an integrated
streamline turn count.

The scale values are arbitrary multiprecision chart times, far outside
practical simulation time ranges. This is a local geometry and callable
roundtrip measurement. It does not certify a recursive transition, continuous
time closure, global flow, outer matching, finite energy, or a full vorticity
vector.

Run command:

```text
python experiments/root_st073/lei_ren_part1_paper_core_scale_geometry.py
```

The machine-readable evidence is in
`lei_ren_part1_paper_core_scale_geometry.json`.
