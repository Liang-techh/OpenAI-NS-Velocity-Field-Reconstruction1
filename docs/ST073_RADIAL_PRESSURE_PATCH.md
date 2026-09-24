# Radial pressure patch experiment

At backward time `tau=0.0084`, the registered annular velocity/pressure field
had only 10 of 81 passing nodes in the coarse local cone map and no passing
two-dimensional grid rectangle. A 16-mode compact pressure family was added:
`(1-(r/R)^2)^5 P_i(2(r/R)^2-1) * 256 s^4(1-s)^4 P_j(2s-1)`, with
`i,j=0..3` and `s=(z-zflat)/(zsupport-zflat)` in the positive axial collar.
Pressure leaves the velocity and its divergence unchanged.

Fitting the broad `4 x 5` grid (`r=0.0075..0.013`, `z=0.003..0.0037`)
is linearly infeasible for local cone ratio `<0.8`. Its top `4 x 2`
subrectangle is feasible but failed one of nine midpoint checks. After
including radial and axial midpoints in the fit, the `7 x 3` grid on
`r=0.0075..0.013`, `z=0.00355..0.0037` is feasible. All 21 fitted nodes
pass, with maximum cone ratio `0.7924`; all 12 additional interior half-grid
points pass, with maximum ratio `0.7796`. This is a sampled patch, not a
continuum or time-persistent cone certificate.

At this time, Gauss-6 physical-volume momentum L2 changed from `183.35`
to `182.21`; Gauss-8 changed from `134.02` to `129.66`. Gauss-6 pointwise
maximum grew slightly (`102974` to `103246`), while Gauss-8 maximum fell
(`182243` to `171077`). These values are still many orders above the
required `1e-3`, and the quadrature orders disagree. The pressure fit is
only a local construction step. It has no pressure-Poisson, complete
spacetime momentum, or exact-curl wave certificate.

Reproduce both screens with:

```powershell
python experiments/root_st073/radial_pressure_patch_screen.py
python experiments/root_st073/radial_pressure_patch_screen.py --dense
```

The machine-readable reports are in
`experiments/root_st073/compact_potential/radial_pressure_patch_screen.json`
and `radial_pressure_patch_dense.json`. The registered dense candidate can be
loaded via `radial_pressure_patch_screen.load_dense_candidate()` for the next
spacetime-persistence and wave-scale experiment.
