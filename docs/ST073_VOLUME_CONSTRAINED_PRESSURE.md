# Volume-aware pressure fit with a fixed spacetime cone

The fixed-pressure spacetime patch in `ST073_TIME_PERSISTENT_CONE.md`
passed the sampled local cone but increased global momentum L2. This
experiment fits an **additional constant compact pressure combination**
to the full physical momentum residual, while retaining the local cone
inequalities at three time slices. The 16 radial/axial pressure gradients
enter the residual analytically, and the velocity is unchanged.

The convex quadratic objective uses physical-volume Gauss-6 quadrature
at `tau=0.00825,0.0084,0.00855` and Gauss-8/10 at `tau=0.0084`.
Each grid's L2 is separately constrained not to exceed the incoming
candidate's value. The local cone fit constrains 63 points; a direct
physical finite-difference check at five times confirms **105/105**
points pass, with maximum sampled cone ratio `0.794315`.

| Time | Grid | Incoming L2 | Fitted L2 | Fitted maximum |
| --- | ---: | ---: | ---: | ---: |
| 0.00825 | Gauss-6 | 183.81 | 161.95 | 93,260 |
| 0.0084 | Gauss-6 | 182.21 | 159.90 | 94,068 |
| 0.00855 | Gauss-6 | 180.87 | 158.15 | 94,704 |
| 0.0084 | Gauss-8 | 129.66 | 129.66 | 177,827 |
| 0.0084 | Gauss-10 | 319.09 | 314.44 | 269,257 |

Thus sampled Gauss-6 L2 improves about 12%, Gauss-10 improves 1.5%,
and Gauss-8 is unchanged. This is a genuine step relative to the
immediate incoming candidate, but the grid orders disagree severely,
and no full-domain or spacetime bound has been obtained. Some sampled
pointwise maxima still worsen. The fitted field remains many orders
above the required `1e-3`; it is **not accepted** as a solution.
The calculation also has no pressure-Poisson compatibility certificate.

The next optimization must address both maximum and L2 residual on
adaptive, cross-order spatial grids, while maintaining a cone margin
over a wider spacetime region. A supported nonaxisymmetric wave and
its paper-style transport/viscous control are still missing.

Reproduce:

```powershell
python experiments/root_st073/radial_pressure_volume_constrained.py
python experiments/root_st073/radial_pressure_time_validate.py --volume-constrained
```

Results are in `experiments/root_st073/compact_potential/` as
`radial_pressure_volume_constrained.json` and
`radial_pressure_volume_cone_validate.json`.
