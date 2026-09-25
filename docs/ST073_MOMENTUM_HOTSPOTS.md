# Momentum hotspot diagnosis and one-mode rejection

The volume-constrained pressure candidate was evaluated at
`tau=0.0084` on Gauss-8 and Gauss-10 physical-volume grids. Its largest
Gauss-8 residual is `177,827` at approximately
`(r,z)=(0.00590,0.00356)`; Gauss-10 reaches `269,257` near
`(0.00704,0.00303)`. Their opposite-z counterparts and nearby
axial-collar nodes are also large. Radial momentum contributes about
`92.6%` of Gauss-8 and `80.6%` of Gauss-10 **L2 squared**. Angular
momentum contributes `4.6%` and `16.4%`, respectively. The remaining
energy is axial. This points to the inner/transition axial collar,
especially its radial transport, as the present dominant defect.

One existing exactly solenoidal `TransitionPoloidalMode` was then added
to the current velocity and screened over signed amplitudes giving
up to one background maximum speed of correction. Full nonlinear
momentum was recombined on Gauss-6/8/10 at the same time. Every
nonzero screened amplitude worsened the worst L2 ratio; the smallest
tested nonzero correction (`2%` relative speed) already slightly
increased L2 on all three grids, and the maximum error rose by at
least `6%` for that sign. The original single collar shape therefore
does not target the hotspot effectively. This is a rejection of one
mode family and amplitude screen, not a general obstruction to
solenoidal poloidal repair.

Next, construct a radial/axial **localized streamfunction basis** with
independent shapes over the hotspot region. Fit its complete nonlinear
momentum jointly across quadrature orders, then recheck the sampled
cone, divergence, time persistence, and finer holdouts. Pressure-only
retuning cannot resolve the angular residual and has so far not
removed the radial peak.

Reproduce:

```powershell
python experiments/root_st073/momentum_hotspot_map.py
python experiments/root_st073/poloidal_hotspot_screen.py
```

Results are in `experiments/root_st073/compact_potential/` as
`momentum_hotspot_map.json` and `poloidal_hotspot_screen.json`.
