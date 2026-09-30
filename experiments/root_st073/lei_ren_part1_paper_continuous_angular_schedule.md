# Continuous angular schedule installed in the shared field

`install_continuous_angular_schedule` updates the existing schedule object,
preserving its identity for velocity, pressure and cumulative consumers.
The full pre-heat angular primitive uses the paper's five J terms, with
`J(1)=1/2` and exact linear tails. Interior J uses the same continuous MP
primitive as the incoming mixed row. Working precision of this quadrature is
100 digits; source arithmetic uses 443 digits in the current candidate.
Neither count certifies quadrature accuracy.

Individual flat J values remain MP numbers even when their exponents cannot
be represented by Decimal. They are combined into log A before converting
the result for legacy schedule clients. No small primitive is assigned zero
to repair Decimal conversion failures.

The pre-heat flattening uses MP sigma and analytic y/Z slopes:

```
q = sigma((y-yv)/Tf)
log Utheta = log A - log(1+Z^2) + q*(log(1+Z^2)-log(2))
dlog Utheta/dZ = -2Z/(1+Z^2)*(1-q)
dlog Utheta/dy = slope_A + sigma'((y-yv)/Tf)/Tf
                              * (log(1+Z^2)-log(2)).
```

Both public radius entry points preserve MP Z instead of coercing through
binary64. The field adapter now consumes the returned angular Z jet and slope.
The inherited relative angular bump still lacks a complete Z jet and moment
closure; no complete angular-correction certificate is asserted.

Heat normalization is recomputed from the continuous log A and the exact
`1-epsilon` factor. The switch/flat-edge heat interpolation is evaluated in MP
arithmetic to retain epsilon, while the heat factor H and its derivative remain
the existing numerical provider. Its uncertainty remains open. This avoids
introducing a normalization jump while replacing the pre-heat primitive.

Incoming regeneration now recomputes log Ep at Rp before transporting inner
offsets and solving the matching coefficients. Therefore axial pulse values
and means use the same updated schedule amplitude. Inherited swirl and future
energy inputs remain explicitly conditional; pressure/five-moment closure is
not implied by schedule identity or by updating their callable source.

The diagnostic covers reference, mu, Z-flatten, steep and delta transitions.
Independent fourth-order Z differences at steps `1e-4` and `1e-30` pass an
absolute `1e-12` gate. Heat-entry log amplitude differs from its pre-heat limit
by approximately `4.97e-292`. The installed field still passes Rp mass and
mass-Z matching at about `6.92e-83`; both materialized terminal residuals remain
nonzero. All global finite-energy, stress and scale-recursion flags remain false.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_angular_schedule.py`
and `python experiments/root_st073/lei_ren_part1_paper_continuous_incoming_outer.py`.
Next: continuous angular correction jets and cumulative moments, certified heat
factor and energy contributions, inner-offset uncertainty, then complete
coefficient tangent and physical radial-tail closure.
