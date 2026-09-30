# Full incoming axial reference bounds

The actual incoming source has `Uz=4*Z*c(y)` with `y=log(R)-log(R_ref)`.
For `y<=1`, `c=1`; for `1<y<exp(Md)`, `c=sigma(1-log(y)/Md)`;
afterward it is exactly zero. The installed candidate uses `Md=.5`.

The full weighted reference integral is the analytic plateau plus

```
Md * integral_0^1 exp(Md*t + power*exp(Md*t)) * sigma(1-t)^m dt
```

with `y=exp(Md*t)`. Both `m=1` and `m=2` are evaluated on outward interval
rectangles. Monotonicity bounds sigma on each panel using its endpoint nodes.
Exact switch endpoint branches avoid division by zero. All arithmetic uses a
local interval context; panel positions come from integer/rational operations.
No support is omitted and no numerical refinement difference is used as an
error certificate. Positive `power` is required by this initial API.

For `power=1`, the full source rows and analytic Z derivatives are

```
I_z = 4*Z*(1 + weighted_m1)
I_uz2 = 16*Z^2*(1 + weighted_m2)
I_z_Z = 4*(1 + weighted_m1)
I_uz2_Z = 32*Z*(1 + weighted_m2).
```

At `Z=.3`, `Md=.5`, nominal installed rows lie inside both 1024- and
4096-panel bounds. At 4096 panels their relative widths are approximately
`2.221e-4` and `2.156e-4`. Saved endpoints are exact dyadic values; displayed
decimals and relative-width summaries are approximate.

The incoming provider exposes optional `full_incoming_enclosure`. This bounds
the reference axial quadrature for declared real Z and Md. It does not bound
uncertainty of the originating source schedule, inner moment offsets, mixed
angular row, inherited swirl energy, or future heat contribution. These gaps
must remain visible when transporting the intervals into the coefficient solve.
Materialized incoming values and coefficients are unchanged.

The coefficient report now transports these reference errors additively:
`delta_base1=exp(log_scale1-yp-log_Ep)*(Iz_iv-stored_Iz)` and
`delta_target=-mu*exp(-yp-2*log_Ep)*(Iuz2_iv-stored_Iuz2)`.
The minus sign reverses energy-target endpoints as required. This preserves
existing inner offsets exactly once. The mixed second row and non-axial target
contributions remain conditional stored data, rather than being relabeled as
enclosed inputs. Reference Z derivative intervals are available but have not
yet been propagated through the coefficient tangent equations.

Run `python experiments/root_st073/lei_ren_part1_paper_continuous_incoming_enclosure.py`.
