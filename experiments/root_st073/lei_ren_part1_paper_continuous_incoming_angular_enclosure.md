# Continuous incoming angular primitive and mixed-row bounds

On the axial cutoff support `y<=exp(Md)`, all shifted primitives in the
source angular schedule vanish: `yd=11+exp(Md)` and later stages start farther
out. Consequently

```
log A(y) = logPstar + .1*y - .6*J(y),
J(y) = integral_0^y sigma(s) ds.
```

The exact switch symmetry gives `J(1)=1/2`, `J(y)=0` below zero, and
`J(y)=y-1/2` above one. `ContinuousIncomingAngular` implements this definition
with MP quadrature on the nontrivial primitive and exact endpoint/tail branches.
Its point values and analytic angular Z/slope jets are incoming-only; the
complete exterior angular velocity still uses the legacy schedule.

The mixed reference factor multiplying `Z/(1+Z^2)` is exactly

```
4*exp(logPstar) * [1/1.6
 + integral_0^1 exp(1.6*y-.6*J(y))dy
 + exp(.3)*integral_1^exp(Md) exp(y)*c(y)dy].
```

There is no extra sqrt(2) in this dimensionless moment convention. The last
integral is the already enclosed incoming axial transition. The unit integral
uses interval rectangles and a monotone cumulative sigma primitive enclosure.
No nominal nested primitive or quadrature refinement is used as proof. Exact
dyadic endpoints and the analytic mixed-row Z derivative are exported.

At the installed `logPstar=14`, `Md=.5`, `Z=.3` fixture, 1024- and 4096-panel
relative mixed-factor widths are approximately 0.2345% and 0.05858%. Both the
old float-backed nominal input and the regenerated MP input are contained.
The exact half primitive is independently contained by the accumulated bounds.

`regenerate_incoming` now computes its mixed nominal factor with the MP provider,
not `schedule._log_A` or float sigma quadrature. At 100 working digits/order96,
the mixed input changes by about `-9.57e-17` relative to the previous receipt.
Working digits are not certified accuracy. Metadata separately preserves the
legacy complete angular-velocity and inherited swirl-energy limitations.

The coefficient report propagates

```
delta_base2 = exp(log_scale2-1.5*yp-2*log_Ep)
              * (I_theta_z_interval-stored_I_theta_z).
```

This keeps the actual inner theta_z offset once. Both reference matching rows
now have quadrature bounds. The actual inner offsets, complete angular field,
swirl energy, future heat contribution, and full Z-jet closure remain open.
The updated field passes Rp mass and mass-Z matching at about `6.92e-83`;
shared runtime checks retain both nonzero terminal residuals.

Run `continuous_incoming_angular_enclosure.py`, `continuous_incoming_outer.py`,
`continuous_solve_enclosure.py`, and `continuous_axial_runtime.py` in the
`experiments/root_st073` directory (the filenames all have the
`lei_ren_part1_paper_` prefix). No finite-energy, scale-recursion or full
momentum certificate follows from these bounds.
